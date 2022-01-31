geojsonio
=========



<!-- [![R-CMD-check](https://github.com/ropensci/geojsonio/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/geojsonio/actions?query=workflow%3AR-CMD-check) -->
[![cran checks](https://cranchecks.info/badges/worst/geojsonio)](https://cranchecks.info/pkgs/geojsonio)
[![R-CMD-check-docker](https://github.com/ropensci/geojsonio/workflows/R-CMD-check-docker/badge.svg)](https://github.com/ropensci/geojsonio/actions?query=workflow%3AR-CMD-check-docker)
[![codecov.io](https://codecov.io/github/ropensci/geojsonio/coverage.svg?branch=master)](https://codecov.io/github/ropensci/geojsonio?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/geojsonio)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/geojsonio)](https://cran.r-project.org/package=geojsonio)

__Convert various data formats to GeoJSON or TopoJSON__

This package is a utility to convert geographic data to GeoJSON and TopoJSON formats. Nothing else. We hope to do this one job very well, and handle all reasonable use cases.

Functions in this package are organized first around what you're working with or want to get, GeoJSON or TopoJSON, then convert to or read from various formats:

* `geojson_list()`/`topojson_list()` - convert to GeoJSON/TopoJSON as R list format
* `geojson_json()`/`topojson_json()` - convert to GeoJSON/TopoJSON as JSON
* `geojson_sp()` - convert output of `geojson_list()` or `geojson_json()` to `sp` spatial objects
* `geojson_sf()` - convert output of `geojson_list()` or `geojson_json()` to `sf` objects
* `geojson_read()`/`topojson_read()` - read a GeoJSON/TopoJSON file from file path or URL
* `geojson_write()`/`topojson_write()` - write a GeoJSON/TopoJSON file locally

Each of the above functions have methods for various objects/classes, including `numeric`, `data.frame`, `list`, `SpatialPolygons`, `SpatialLines`, `SpatialPoints`, etc.

Additional functions:

* `map_gist()` - push up a GeoJSON or topojson file as a GitHub gist (renders as an interactive map)
* `map_leaf()` - create a local interactive map using the `leaflet` package

## \*json Info

* GeoJSON - [spec](https://tools.ietf.org/html/rfc7946)
* [GeoJSON lint](https://geojsonlint.com/)
* TopoJSON - [spec](https://github.com/topojson/topojson-specification/blob/master/README.md)


## Install

A note about installing `rgeos` - built on top of C libraries, and installation often causes trouble for Linux users because no binaries are provided on CRAN for those platforms. Other dependencies in `geojsonio` should install easily automatically when you install `geojsonio`.

_Mac_

Install `GDAL` on the command line first, e.g., using `homebrew`

```
brew install gdal
```

Then install `rgeos`


```r
install.packages("rgeos", type = "source")
```

_Linux_

Get deps first

```
sudo apt-get install libgdal1-dev libgdal-dev libgeos-c1 libproj-dev
```

> Note: if you have trouble installing rgeos, try installing `libgeos++-dev`

Then install `rgeos`


```r
install.packages("rgeos", type = "source")
```

__Install geojsonio__

Stable version from CRAN


```r
install.packages("geojsonio")
```

Or development version from GitHub


```r
install.packages("remotes")
remotes::install_github("ropensci/geojsonio")
```


```r
library("geojsonio")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/geojsonio/issues).
* License: MIT
* Get citation information for `geojsonio` in R doing `citation(package = 'geojsonio')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
geojsonio 0.9.4
===============

### BUG FIXES

* fix for `sprintf()` usage within the `projections()` function; only run sprintf on a particular string if it has length > 0 (#172)
* fix for `as.json()` when the input is the output of `topojson_list()` - we weren't constructing the TopoJSON arcs correctly (#160)
* fix to `geojson_read()`: now using package `geojsonsf` to read geojson (#163)

geojsonio 0.9.2
===============

### BUG FIXES

* fix a test for change in `stringsAsFactors` behavior in R v4 (#166) (#167)
* temporarily make `topojson_write()` defunct until we can sort out issues with new sf version (#168)

geojsonio 0.9.0
===============

### NEW FEATURES

* `geojson_sf()` and `geojson_sp()` now accept strings in addition to `json`, `geoson_list` and `geojson_json` types (#164)

### MINOR IMPROVEMENTS

* `topojson_json()` and `topojson_list()` gain params `object_name` and `quantization` to pass through to `geojson_json()` (#158)
* replace httr with crul (#105)
* rgdal replaced with sf throughout the package; all `writeOGR` replaced with `st_write` and `readOGR` with `st_read`; this should not create any user facing changes, but please let us know if you have problems with this version (#41) (#150) (#157)


geojsonio 0.8.0
===============

## NEW FEATURES

* `geojson_read()` gains new S3 method `geojson_read.PqConnection` for connecting to a PostgreSQL database set up with PostGIS. See also `?postgis` for notes on Postgis installation, and setting up some simple data in Postgis from this package (#61) (#155) thanks to @fxi

### MINOR IMPROVEMENTS

* `geojson_read()` instead of going through package `sp` now goes through package `sf` for a significant speed up, see https://github.com/ropensci/geojsonio/issues/136#issuecomment-546123078 (#136)
* `geojson_list()` gains parameter `precision` to adjust number of decimal places used. only applies to classes from packages sp and rgeos (#152) (related to #141) thanks to @ChrisJones687
* improve dependency installation notes in README (#149) (#151) thanks @setgree and @nickto
* move to using markdown docs
* `file_to_geojson()` now using https protocol instead of http for the online ogre service called when using `method = "web"`

### BUG FIXES

* fix `geojson_read()` to fail better when using `method="web"`; and update docs to note that `method="web"` can result if file size issues, but `method="local"` should not have such issues (#153)
* change name of `print.location` method to not conflict with `dplyr` (#154)


geojsonio 0.7.0
===============

## NEW FEATURES

* `geo2topo()` gains a new parameter `quantization` to quantize geometry prior to computing topology. because `topojson_write()` uses `geo2topo()` internally, `topojson_write()` also gains a `quantization` parameter that's passed to `geo2topo()` internally (#138) thanks @pvictor

### MINOR IMPROVEMENTS

* use package `sf` instead of `sp` in `topojson_read()`. note that the return object is now class sf instead of classes from the sp package (#144) (#145)
* the `type` parameter in `topojson_json()` now set to `type="auto"` if the input is an sf/sfc/sfg class object (#139) (#146)
* fix to `geojson_list.sfc()` for changes in sf >= v0.7, which names geometries, but that's not valid geojson (#142)

### DEPRECATED AND DEFUNCT

* The two linting functions in this package, `lint()` and `validate()`
are now defunct. They have been marked as deprecated since `v0.2`. See the package `geojsonlint` on CRAN for linting geojson functionality (#135) (#147)


geojsonio 0.6.0
===============

## NEW FEATURES

* `topojson_write()` gains a new parameter `object_name`. With it you can set the name for the resulting TopoJSON object name created. As part of this `geo2topo()` also gains a new parameter, similarly called `object_name`, that does the same thing as for `topojson_write()`. (#129) (thanks @josiekre) PR (#131)
* As part of PR (#132) we added a new function `geojson_sf()` to convert output of `geojson_list()` or `geojson_json()` to `sf` package classes - as an analog to `geojson_sp()`

### MINOR IMPROVEMENTS

* `geojson_json()` gains option with the `type` parameter to skip a coercion to the `geojson` package class `geoclass`. Using `type = "skip"` you can skip the `geoclass` class coercion, which in some cases with large datasets should have performance improvements (#128) PR (#133)

### BUG FIXES

* A bug arose in `geojson_sp()` with the newest version of `rgdal`. This was resolved by using the `sf` package instead to read GeoJSON. This had a knock-on benefit of speeding up reading GeoJSON. In addition, `sf` is now in `Imports` instead of `Suggests`  (#130) PR (#132)


geojsonio 0.5.0
===============

### NEW FEATURES

* gains new function `geojson_atomize` to "atomize" a FeatureCollection
into its features, or a GeometryCollection into its geometries (#120)
via (#119) thx @SymbolixAU
* gains new functions `topojson_list` and `topojson_json` for converting
many input types with spatial data to TopoJSON, both as lists and
as JSON (#117)
* `geojson_json` uses brief output provided by the `geojson`
package - this makes it less frustrating when you have an especially
large geojson string that prints to console - this instead prints a
brief summary of the GeoJSON object (#86) (#124)

### MINOR IMPROVEMENTS

* doing a much more thorough job of cleaning up temp files that are
necessarily generated due to having to go to disk sometimes (#122)
* @ateucher made improvements to `geojson_json` to make `type`
parameter more flexible (#125)

### BUG FIXES

* Fixe bug in `topojson_write` - we were writing topojson file, but also
a geojson file - we now cleanup the geojson file (#127)


geojsonio 0.4.2
===============

### BUG FIXES

* Fix package so that we load `topojson-server.js` from within the
package instead of from the web. This makes it so that the
package doesn’t make any web requests on load, which prevented package
from loading when no internet connection available. (#118)


geojsonio 0.4.0
===============

### NEW FEATURES

* Gains new functions `geo2topo`, `topo2geo`, `topojson_write`, and `topojson_read` for working with TopoJSON data - associated with this, we
now import `geojson` package (#24) (#100)

### MINOR IMPROVEMENTS

* Updated vignette with details on the GeoJSON specification to the new
specification at <https://tools.ietf.org/html/rfc7946> (#114)



geojsonio 0.3.8
===============

### MINOR IMPROVEMENTS

* `geojson_write` and `geojson_json` now pass `...` argument through to
`rgdal::writeOGR` or `jsonlite::toJSON` depending on the class/method. For
those methods that use the latter, this now allows setting of the `na`
argument to control how `NA` values are represented in json, and the
`pretty` argument to control whether or the resulting json is
pretty-formated or compact (#109) (#111)
* Spelling/grammar fixes, thanks @patperu ! (#106)

### BUG FIXES

* `geojson_json` and `geojson_write` now convert unsupported classes to
their basic class before conversion and/or writing to geojson. This was most
commonly occurring with fields in `sf` objects calculated by `sf::st_area`
and `sf::st_length` which were of class `units`. (#107)
* Fixed a bug occurring with `GDAL` version >= 2.2.0 where the layer name in
a geojson file was not detected properly (#108)


geojsonio 0.3.2
==============

### BUG FIXES

* Fix to tests for internal fxn `convert_wgs84` to do minimal test of
output, and to conditionally test only if `sf` is available (#103)


geojsonio 0.3.0
==============

### NEW FEATURES

* `geojson_json`, `geojson_list`, and `geojson_write` gain new S3 methods:
`sf`, `sfc`, and `sfg` - the three classes in the `sf` package (#95)
* `geojson_json`, `geojson_list`, and `geojson_write` gain two new
parameters each: `convert_wgs84` (boolean) to convert to WGS84 or not (the
projection assumed for GeoJSON)  and `crs` to assign a CRS if known
(#101) (#102)

### MINOR IMPROVEMENTS

* `geojson_json()` for non-sp classes now only keeps seven decimal places
in the coordinates. This follows the default that GDAL uses.
* Now namespacing base package calls for `methods`/`stats`/`utils`
instead of importing them
* Improved documentation for `method` parameter in `geojson_read`
clarifying what the options are for (#93) thanks @bhaskarvk
* Internal fxn `to_json` now defaults to 7 digits, which is used in
`as.json` and `geojson_json` (#96)

### BUG FIXES

* Fix to `geojson_read` to read correctly from a URL - in addition
to file paths (#91) (#92) thanks @lecy
* Fix to `geojson_read` to read non-`.geojson` extensions (#93)
thanks @bhaskarvk


geojsonio 0.2.0
===============

### MINOR IMPROVEMENTS

* Major performance improvement for `geojson_json()` - moved to
reading in json with `readr::read_file()` (#85) thanks @javrucebo !
* Now requiring explicit versions of some package dependencies
* Removed the startup message

### BUG FIXES

* Changed `file_to_geojson()` to use `httr::write_disk()` instead of
`download.file()` (#83) thanks @patperu

### DEPRECATED AND DEFUNCT

* The two linting functions in this package, `lint()` and `validate()`
are now deprecated, and will be defunct in the next version of this
package. See the new package `geojsonlint` on CRAN for linting
geojson functionality (#82)

geojsonio 0.1.8
===============

### NEW FEATURES

* New method `geojson_sp.json()` added to `geojson_sp()` to handle json
class inputs

### MINOR IMPROVEMENTS

* Added `encodin="UTF-8"` to `httr::content()` calls

### BUG FIXES

* `geojson_write()` didn't overwrite existing files despite saying so.
New parameter added to the function `overwrite` to specify whether to
overwrite a function or not, which defaults to `TRUE` (#81)
thanks @Robinlovelace !

geojsonio 0.1.6
===============

### NEW FEATURES

* New function `geojson_sp()` to convert output of `geojson_list()` or
`geojson_json()` to spatial classes (e.g., `SpatialPointsDataFrame`) (#71)

### MINOR IMPROVEMENTS

* Startup message added to notify users to ideally update to `rgdal > v1.1-1`
given fix to make writing multipolygon objects to geojson correct (#69)
* Filled out test suite more (#46)

### BUG FIXES

* Fix to `lint()` function, due to bug in passing data to the Javascript
layer (#73)
* Fixes to `as.json()` (#76)

geojsonio 0.1.4
===============

### NEW FEATURES

* New function `map_leaf()` uses the `leaflet` package to make maps, with
S3 methods for most spatial classes as well as most R classes, including
data.frame's, lists, vectors, file inputs, and more (#48)
* `geojson_read()` now optionally can give back a spatial class object,
just a convenience in case you want to not get back geojson, but a
spatial class (#60)

### MINOR IMPROVEMENTS

* Now that `leaflet` R package is on CRAN, put back in examples using
it to make maps (#49)
* Added a linter for list inputs comined with `geometry="polygon"` to
all `geojson_*()` functions that have `.list` methods. This checks to
make sure inputs have the same first and last coordinate pairs to
close the polygon (#34)

### BUG FIXES

* Importing all non-base R funtions, including from `methods`, `stats` and `utils`
packages (#62)
* Fixed bug in `geojson_write()` in which geojson style names were altered
on accident (#56)

geojsonio 0.1.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.3 Patched
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

I have run R CMD check on the 12 reverse dependencies. Summary at <https://github.com/ropensci/geojsonio/blob/master/revdep/README.md>. No problems were found related to this package.

-------

This version fixes three bugs.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/geojsonio/issues)

### Code contributions

Please keep your text/code/docs within 80 character width in all files. Thanks!

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/geojsonio.git`
* Make sure to track progress upstream (i.e., on our version of `geojsonio` at `ropensci/geojsonio`) by doing `git remote add upstream https://github.com/ropensci/geojsonio.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/geojsonio`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
js scripts/libraries
====================

## geojsonhint

From <https://www.npmjs.com/package/@mapbox/geojsonhint>

## turf extent

From <https://www.npmjs.com/package/turf-extent>

## topojson server

On `v3.0` from <https://unpkg.com/topojson-server@3>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.3 Patched (2020-12-29 r79725) |
|os       |macOS Catalina 10.15.7                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-01-12                                  |

# Dependencies

|package   |old   |new   |Δ  |
|:---------|:-----|:-----|:--|
|geojsonio |0.9.2 |0.9.4 |*  |
|cpp11     |NA    |0.2.5 |*  |

# Revdeps

*Wow, no problems at all. :)*# Check times

|   |package        |version | check_time|
|:--|:--------------|:-------|----------:|
|6  |rmapshaper     |0.3.0   |      135.8|
|7  |rmapzen        |0.3.3   |       48.8|
|5  |repijson       |0.1.0   |         47|
|8  |webglobe       |1.0.2   |       38.8|
|1  |antaresViz     |0.11    |       35.7|
|3  |leaflet.extras |0.2     |       20.4|
|4  |mregions       |0.1.6   |       19.2|
|2  |leaflet.esri   |0.2     |       16.8|


Hi,

This is an automated email to let you know about the release of {{{ my_package }}}, which I'll submit to CRAN on {{{ date }}}.

To check for potential problems, I ran `R CMD check` on your package {{{your_package}}} (v{{{your_version}}}).

I found: {{{your_summary}}}.

{{#you_have_problems}}
{{{your_results}}}

To get the development version of {{{ my_package }}} so you can run the checks yourself, you can run:

    # install.packages("devtools")
    devtools::install_github("{{my_github}}")

To see what's changed visit <https://github.com/{{{my_github}}}/blob/master/NEWS.md>.

{{/you_have_problems}}
{{^you_have_problems}}
It looks like everything is ok, so you don't need to take any action, but you might want to read the NEWS, <https://github.com/{{{my_github}}}/blob/master/NEWS.md>, to see what's changed.
{{/you_have_problems}}


If you have any questions about this email, please feel free to respond directly.

Regards,

{{{ me }}}
*Wow, no problems at all. :)*geojsonio
=========

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
  eval = FALSE
)
```

<!-- [![R-CMD-check](https://github.com/ropensci/geojsonio/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/geojsonio/actions?query=workflow%3AR-CMD-check) -->
[![cran checks](https://cranchecks.info/badges/worst/geojsonio)](https://cranchecks.info/pkgs/geojsonio)
[![R-CMD-check-docker](https://github.com/ropensci/geojsonio/workflows/R-CMD-check-docker/badge.svg)](https://github.com/ropensci/geojsonio/actions?query=workflow%3AR-CMD-check-docker)
[![codecov.io](https://codecov.io/github/ropensci/geojsonio/coverage.svg?branch=master)](https://codecov.io/github/ropensci/geojsonio?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/geojsonio)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/geojsonio)](https://cran.r-project.org/package=geojsonio)

__Convert various data formats to GeoJSON or TopoJSON__

This package is a utility to convert geographic data to GeoJSON and TopoJSON formats. Nothing else. We hope to do this one job very well, and handle all reasonable use cases.

Functions in this package are organized first around what you're working with or want to get, GeoJSON or TopoJSON, then convert to or read from various formats:

* `geojson_list()`/`topojson_list()` - convert to GeoJSON/TopoJSON as R list format
* `geojson_json()`/`topojson_json()` - convert to GeoJSON/TopoJSON as JSON
* `geojson_sp()` - convert output of `geojson_list()` or `geojson_json()` to `sp` spatial objects
* `geojson_sf()` - convert output of `geojson_list()` or `geojson_json()` to `sf` objects
* `geojson_read()`/`topojson_read()` - read a GeoJSON/TopoJSON file from file path or URL
* `geojson_write()`/`topojson_write()` - write a GeoJSON/TopoJSON file locally

Each of the above functions have methods for various objects/classes, including `numeric`, `data.frame`, `list`, `SpatialPolygons`, `SpatialLines`, `SpatialPoints`, etc.

Additional functions:

* `map_gist()` - push up a GeoJSON or topojson file as a GitHub gist (renders as an interactive map)
* `map_leaf()` - create a local interactive map using the `leaflet` package

## \*json Info

* GeoJSON - [spec](https://tools.ietf.org/html/rfc7946)
* [GeoJSON lint](https://geojsonlint.com/)
* TopoJSON - [spec](https://github.com/topojson/topojson-specification/blob/master/README.md)


## Install

A note about installing `rgeos` - built on top of C libraries, and installation often causes trouble for Linux users because no binaries are provided on CRAN for those platforms. Other dependencies in `geojsonio` should install easily automatically when you install `geojsonio`.

_Mac_

Install `GDAL` on the command line first, e.g., using `homebrew`

```
brew install gdal
```

Then install `rgeos`

```{r}
install.packages("rgeos", type = "source")
```

_Linux_

Get deps first

```
sudo apt-get install libgdal1-dev libgdal-dev libgeos-c1 libproj-dev
```

> Note: if you have trouble installing rgeos, try installing `libgeos++-dev`

Then install `rgeos`

```{r}
install.packages("rgeos", type = "source")
```

__Install geojsonio__

Stable version from CRAN

```{r}
install.packages("geojsonio")
```

Or development version from GitHub

```{r}
install.packages("remotes")
remotes::install_github("ropensci/geojsonio")
```

```{r}
library("geojsonio")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/geojsonio/issues).
* License: MIT
* Get citation information for `geojsonio` in R doing `citation(package = 'geojsonio')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "GeoJSON Specification"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{GeoJSON Specification}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE, warning=FALSE, message=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = FALSE
)
```

In `geojsonio` we follow the below guidelines (RFC7946) for GeoJSON, and try
to preserve CRS information, and bbox information when converting formats.

The following are the guidelines for CRS and bounding boxes for geojson,
annotated as needed, get complete guidelines at <https://tools.ietf.org/html/rfc7946>

## The Structure of GeoJSON

See <https://tools.ietf.org/html/rfc7946#section-3> for further information.

__GeoJSON text__: a JSON text and consists of a single GeoJSON object.

__GeoJSON object__: represents a Geometry, Feature, or collection of Features
(i.e., FeatureCollection).

* has a member with the name "type".  The value of the member
MUST be one of the GeoJSON types.
* MAY have a "bbox" member, the value of which MUST be a bounding
box array (see Section 5).
* MAY have other members (see Section 6).

__Geometry object__: represents points, curves, and surfaces in coordinate
space. Every Geometry object is a GeoJSON object no matter where it occurs
in a GeoJSON text.

* The value of a Geometry object's "type" member MUST be one of the
  seven geometry types.
* A GeoJSON Geometry object of any type other than
  "GeometryCollection" has a member with the name "coordinates".
  The value of the "coordinates" member is an array.  The structure
  of the elements in this array is determined by the type of
  geometry.

__Position__: the fundamental geometry construct.  The "coordinates" member of
a Geometry object is composed of either:

* one position in the case of a Point geometry,
* an array of positions in the case of a LineString or MultiPoint geometry,
* an array of LineString or linear ring coordinates in the case of a
Polygon or MultiLineString geometry, or
* an array of Polygon coordinates in the case of a MultiPolygon geometry

__Type of Geometries__:

* Point
* MultiPoint
* LineString
* MultiLineString
* Polygon
* MultiPolygon
* GeometryCollection

__Feature Object__: A Feature object represents a spatially bounded thing.  Every
Feature object is a GeoJSON object no matter where it occurs in a GeoJSON text. A
Feature object has a "type" member with the value "Feature"; has a member with the
name "geometry", the value of which geometry member as defined above or a JSON
null value. A Feature object has a member with the name "properties"; the value
of the properties member is an object (any JSON object or a JSON null value).

__FeatureCollection Object__: A GeoJSON object with the type "FeatureCollection" is a
FeatureCollection object.  A FeatureCollection object has a member with the name
"features".  The value of "features" is a JSON array. Each element of the array
is a Feature object as defined above.  It is possible for this array to be empty.

## CRS (Coordinate Reference System)

See <https://tools.ietf.org/html/rfc7946#page-12> for further information.

* The coordinate reference system for all GeoJSON coordinates is a geographic
coordinate reference system, using the World Geodetic System 1984 (WGS84)
datum, with longitude and latitude units of decimal degrees.  This is equivalent
to the coordinate reference system identified by the Open Geospatial Consortium (OGC)
URN `urn:ogc:def:crs:OGC::CRS84`. An OPTIONAL third-position element SHALL
be the height in meters above or below the WGS 84 reference ellipsoid.  In the
absence of elevation values, applications sensitive to height or depth SHOULD
interpret positions as being at local ground or sea level.
* The crs member has been removed.
* RFC7946 does not that "where all involved parties have a prior arrangement,
alternative coordinate reference systems can be used without risk of data being
misinterpreted."


## Bounding Boxes

See <https://tools.ietf.org/html/rfc7946#page-12> for further information.

To include information on the coordinate range for Geometries, Features, or
FeatureCollections, a GeoJSON object may have a member named `bbox`. The value of the bbox
member must be a 2*n array where n is the number of dimensions represented in the
contained geometries, with the lowest values for all axes followed by the highest
values. The axes order of a bbox follows the axes order of geometries.

Example of a 2D bbox member on a Feature:

```
{
  "type": "Feature",
  "bbox": [-10.0, -10.0, 10.0, 10.0],
  "geometry": {
    "type": "Polygon",
    "coordinates": [[
      [-10.0, -10.0], [10.0, -10.0], [10.0, 10.0], [-10.0, -10.0]
    ]]
  }
  ...
}
```

Example of a 2D bbox member on a FeatureCollection:

```
{
  "type": "FeatureCollection",
  "bbox": [100.0, 0.0, 105.0, 1.0],
  "features": [
    ...
  ]
}
```

Example of a 3D bbox member with a depth of 100 meters on a FeatureCollection:

```
{
  "type": "FeatureCollection",
  "bbox": [100.0, 0.0, -100.0, 105.0, 1.0, 0.0],
  "features": [
    ...
  ]
}
```

## Coordinate Precision

See <https://tools.ietf.org/html/rfc7946#page-18> for further information.

The size of a GeoJSON text in bytes is a major interoperability
consideration, and precision of coordinate values has a large impact
on the size of texts.  A GeoJSON text containing many detailed
Polygons can be inflated almost by a factor of two by increasing
coordinate precision from 6 to 15 decimal places.  For geographic
coordinates with units of degrees, 6 decimal places (a default common
in, e.g., sprintf) amounts to about 10 centimeters, a precision well
within that of current GPS systems.  Implementations should consider
the cost of using a greater precision than necessary.
---
title: "geojsonio vignette"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{geojsonio vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo=FALSE, warning=FALSE, message=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 10,
  fig.path = "../man/figures/"
)
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
```

`geojsonio` converts geographic data to geojson and topojson formats. Nothing else. We hope to do this one job very well, and handle all reasonable use cases.

Functions in this package are organized first around what you're working with or want to get, geojson or topojson, then convert to or read from various formats:

* `geojson_list()`/`topojson_list()` - convert to GeoJSON/TopoJSON as R list format
* `geojson_json()`/`topojson_json()` - convert to GeoJSON/TopoJSON as JSON
* `geojson_sp()` - convert output of `geojson_list()` or `geojson_json()` to spatial objects
* `geojson_read()`/`topojson_read()` - read a GeoJSON/TopoJSON file from file path or URL
* `geojson_write()`/`topojson_write()` - write a GeoJSON/TopoJSON file locally


Each of the above functions have methods for various objects/classes, including `numeric`, `data.frame`, `list`, `SpatialPolygons`, `SpatialLines`, `SpatialPoints`, etc.

Additional functions:

* `map_gist()` - push up a geojson or topojson file as a GitHub gist (renders as an interactive map) - See the _maps with geojsonio_ vignette.
* `map_leaf()` - create a local interactive map with the `leaflet` package - See the _maps with geojsonio_ vignette.

## Install

Install rgdal - in case you can't get it installed from binary , here's what works on a Mac (change to the version of `rgdal` and `GDAL` you have).

```r
install.packages("http://cran.r-project.org/src/contrib/rgdal_1.1-3.tar.gz", repos = NULL, type="source", configure.args = "--with-gdal-config=/Library/Frameworks/GDAL.framework/Versions/1.11/unix/bin/gdal-config --with-proj-include=/Library/Frameworks/PROJ.framework/unix/include --with-proj-lib=/Library/Frameworks/PROJ.framework/unix/lib")
```

Stable version from CRAN

```r
install.packages("geojsonio")
```

Development version from GitHub

```r
remotes::install_github("ropensci/geojsonio")
```

```{r}
library("geojsonio")
```

## GeoJSON

### Convert various formats to geojson

From a `numeric` vector of length 2

as _json_

```{r}
geojson_json(c(32.45, -99.74))
```

as a __list__

```{r output.lines=1:10}
geojson_list(c(32.45, -99.74))
```

From a `data.frame`

as __json__

```{r}
library('maps')
data(us.cities)
geojson_json(us.cities[1:2, ], lat = 'lat', lon = 'long')
```

as a __list__

```{r output.lines=1:10}
geojson_list(us.cities[1:2, ], lat = 'lat', lon = 'long')
```

From `SpatialPolygons` class

```{r}
library('sp')
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
  c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
  c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
```

to __json__

```{r}
geojson_json(sp_poly)
```

to a __list__

```{r output.lines=1:10}
geojson_list(sp_poly)
```

From `SpatialPoints` class

```{r}
x <- c(1, 2, 3, 4, 5)
y <- c(3, 2, 5, 1, 4)
s <- SpatialPoints(cbind(x, y))
```

to __json__

```{r}
geojson_json(s)
```

to a __list__

```{r output.lines=1:10}
geojson_list(s)
```

### Write geojson

```{r}
library('maps')
data(us.cities)
geojson_write(us.cities[1:2, ], lat = 'lat', lon = 'long')
```

### Read geojson

```{r}
library("sp")
file <- system.file("examples", "california.geojson", package = "geojsonio")
out <- geojson_read(file, what = "sp")
plot(out)
```

## Topojson

To JSON

```{r}
topojson_json(c(-99.74,32.45))
```

To a list

```{r output.lines=1:20}
library(sp)
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPoints(cbind(x,y))
topojson_list(s)
```

Read from a file

```{r}
file <- system.file("examples", "us_states.topojson", package = "geojsonio")
out <- topojson_read(file)
summary(out)
```

Read from a URL

```{r}
url <- "https://raw.githubusercontent.com/shawnbot/d3-cartogram/master/data/us-states.topojson"
out <- topojson_read(url)
```

Or use `as.location()` first

```{r}
(loc <- as.location(file))
out <- topojson_read(loc)
```
---
title: "maps with geojsonio"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{maps with geojsonio}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo=FALSE, warning=FALSE, message=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

`geojsonio` creates geojson from various inputs - and can easily feed into tools for making maps with geojson data.

```{r}
library("geojsonio")
```


## Mapping with leaflet

### With geojsonio::map_leaf()

#### From a file

```{r eval=FALSE}
file <- "myfile.geojson"
geojson_write(us_cities[1:20, ], lat='lat', lon='long', file = file)
map_leaf(as.location(file))
```


#### From a SpatialGridDataFrame

```{r eval=FALSE}
sgdim <- c(3, 4)
sg <- SpatialGrid(GridTopology(rep(0, 2), rep(10, 2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
map_leaf(sgdf)
```


### DIY

#### Example 1: Map of California

```{r eval=FALSE}
library("leaflet")
file <- system.file("examples", "california.geojson", package = "geojsonio")
out <- as.json(geojson_read(file))
leaflet() %>% 
  addProviderTiles("Stamen.Toner") %>% 
  setView(lng = -119, lat = 37, zoom = 6) %>%
  addGeoJSON(out)
```


#### Example 2: Map of two polygons

```{r eval=FALSE}
library('sp')
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
    c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
    c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
json <- geojson_json(sp_poly)

leaflet() %>% 
  addProviderTiles("Stamen.Toner") %>% 
  setView(lng = -90, lat = 41, zoom = 4) %>%
  addGeoJSON(json)
```


## Mapping with GitHub gists

### `data.frame`

> Also, can do so from data.frames with polygons, lists, matrices, vectors, and json strings

```{r eval=FALSE}
map_gist(us_cities)
```


### `SpatialPoints` class

```{r eval=FALSE}
library("sp")
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPoints(cbind(x,y))
map_gist(s)
```


### `SpatialPixelsDataFrame` class

```{r eval=FALSE}
library("sp")
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")], data = canada_cities)
)
map_gist(pixelsdf)
```


> Many other spatial classes supported
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bounds.R
\name{bounds}
\alias{bounds}
\title{Get bounds for a list or geo_list}
\usage{
bounds(x, ...)
}
\arguments{
\item{x}{An object of class list or geo_list}

\item{...}{Ignored}
}
\value{
A vector of the form min longitude, min latitude, max longitude,
max latitude
}
\description{
Get bounds for a list or geo_list
}
\examples{
# numeric 
vec <- c(-99.74,32.45)
x <- geojson_list(vec)
bounds(x)

# list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
x <- geojson_list(mylist)
bounds(x)

# data.frame
x <- geojson_list(states[1:20,])
bounds(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_to_geojson.r
\name{file_to_geojson}
\alias{file_to_geojson}
\title{Convert spatial data files to GeoJSON from various formats.}
\usage{
file_to_geojson(
  input,
  method = "web",
  output = ".",
  parse = FALSE,
  encoding = "CP1250",
  verbose = FALSE,
  ...
)
}
\arguments{
\item{input}{The file being uploaded, path to the file on your machine.}

\item{method}{(character) One of "web" (default) or "local". Matches on
partial strings. This parameter determines how the data is
read. "web" means we use the Ogre web service, and "local" means we use
\pkg{sf}. See Details fore more.}

\item{output}{Destination for output geojson file. Defaults to current
working directory, and gives a random alphanumeric file name}

\item{parse}{(logical) To parse geojson to data.frame like structures if
possible. Default: \code{FALSE}}

\item{encoding}{(character) The encoding passed to \code{\link[sf:st_read]{sf::st_read()}}.
Default: CP1250}

\item{verbose}{(logical) Printing of \code{\link[sf:st_read]{sf::st_read()}} progress.
Default: \code{FALSE}}

\item{...}{Additional parameters passed to \code{\link[sf]{st_read}}}
}
\value{
path for the geojson file
}
\description{
You can use a web interface called Ogre, or do conversions locally using the
sf package.
}
\section{Method parameter}{

The web option uses the Ogre web API. Ogre currently has an output size
limit of 15MB. See here \url{http://ogre.adc4gis.com/} for info on the
Ogre web API. The local option uses the function \code{\link[sf]{st_write}}
from the package rgdal.
}

\section{Ogre}{

Note that for Shapefiles, GML, MapInfo, and VRT, you need to send zip files
to Ogre. For other file types (.bna, .csv, .dgn, .dxf, .gxt, .txt, .json,
.geojson, .rss, .georss, .xml, .gmt, .kml, .kmz) you send the actual file
with that file extension.
}

\section{Linting GeoJSON}{

If you're having trouble rendering GeoJSON files, ensure you have a valid
GeoJSON file by running it through the package \pkg{geojsonlint}, which
has a variety of different GeoJSON linters.
}

\section{File size}{

When using \code{method="web"}, be aware of file sizes.
https://ogre.adc4gis.com that we use for this option does not document
what file size is too large, but you should get an error message like
"maximum file length exceeded" when that happens. \code{method="local"}
shouldn't be sensitive to file sizes.
}

\examples{
\dontrun{
file <- system.file("examples", "norway_maple.kml", package = "geojsonio")

# KML type file - using the web method
file_to_geojson(input=file, method='web', output='kml_web')
## read into memory
file_to_geojson(input=file, method='web', output = ":memory:")
file_to_geojson(input=file, method='local', output = ":memory:")

# KML type file - using the local method
file_to_geojson(input=file, method='local', output='kml_local')

# Shp type file - using the web method - input is a zipped shp bundle
file <- system.file("examples", "bison.zip", package = "geojsonio")
file_to_geojson(file, method='web', output='shp_web')

# Shp type file - using the local method - input is the actual .shp file
file <- system.file("examples", "bison.zip", package = "geojsonio")
dir <- tempdir()
unzip(file, exdir = dir)
list.files(dir)
shpfile <- file.path(dir, "bison-Bison_bison-20130704-120856.shp")
file_to_geojson(shpfile, method='local', output='shp_local')

# geojson with .json extension
## this doesn't work anymore, hmmm
# x <- gsub("\n", "", paste0('https://gist.githubusercontent.com/hunterowens/
# 25ea24e198c80c9fbcc7/raw/7fd3efda9009f902b5a991a506cea52db19ba143/
# wards2014.json', collapse = ""))
# res <- file_to_geojson(x)
# jsonlite::fromJSON(res)
# res <- file_to_geojson(x, method = "local")
# jsonlite::fromJSON(res)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojsonio-package.r
\docType{data}
\name{us_cities}
\alias{us_cities}
\title{This is the same data set from the maps library, named differently}
\format{
A list with 6 components, namely "name", "country.etc", "pop",
"lat", "long", and "capital", containing the city name, the state
abbreviation, approximate population (as at January 2006), latitude,
longitude and capital status indication (0 for non-capital, 1 for capital,
2 for state capital.
}
\description{
This database is of us cities of population greater than about 40,000.
Also included are state capitals of any population size.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_atomize.R
\name{geojson_atomize}
\alias{geojson_atomize}
\title{Atomize}
\usage{
geojson_atomize(x, combine = TRUE)
}
\arguments{
\item{x}{(geo_list/geo_json/json/character) input object, either
\code{geo_json}, \code{geo_list}, \code{json}, or \code{character} class.
If \code{character}, must be valid JSON}

\item{combine}{(logical) only applies to \code{geo_json/json} type inputs.
combine valid JSON objects into a single valid JSON object. Default:
\code{TRUE}}
}
\value{
same class as input object, but modified
}
\description{
Atomize
}
\details{
A FeatureCollection is split into many Feature's, and
a GeometryCollection is split into many geometries

Internally we use \pkg{jqr} for JSON parsing
}
\examples{
################# lists 
# featurecollection -> features
mylist <- list(list(latitude=30, longitude=120, marker="red"),
          list(latitude=30, longitude=130, marker="blue"))
(x <- geojson_list(mylist))
geojson_atomize(x)

# geometrycollection -> geometries
mylist <- list(list(latitude=30, longitude=120, marker="red"),
          list(latitude=30, longitude=130, marker="blue"))
(x <- geojson_list(mylist, type = "GeometryCollection"))
geojson_atomize(x)

# sf class
library(sf)
p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
poly <- rbind(c(1,1), c(1,2), c(2,2), c(1,1))
poly_sfg <- st_polygon(list(p1))
(x <- geojson_list(poly_sfg))
geojson_atomize(x)

################# json 
# featurecollection -> features
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
(x <- geojson_json(mylist))
geojson_atomize(x)
geojson_atomize(x, FALSE)

# geometrycollection -> geometries
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
(x <- geojson_json(mylist, type = "GeometryCollection"))
geojson_atomize(x)
geojson_atomize(x, FALSE)

# sf class
library(sf)
nc <- st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE)
(x <- geojson_json(nc))
geojson_atomize(x)
geojson_atomize(x, FALSE)

################# character
# featurecollection -> features
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
(x <- geojson_json(mylist))
geojson_atomize(unclass(x))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_style.R
\name{geojson_style}
\alias{geojson_style}
\title{Style a data.frame or list prior to converting to geojson}
\usage{
geojson_style(
  input,
  var = NULL,
  var_col = NULL,
  var_sym = NULL,
  var_size = NULL,
  var_stroke = NULL,
  var_stroke_width = NULL,
  var_stroke_opacity = NULL,
  var_fill = NULL,
  var_fill_opacity = NULL,
  color = NULL,
  symbol = NULL,
  size = NULL,
  stroke = NULL,
  stroke_width = NULL,
  stroke_opacity = NULL,
  fill = NULL,
  fill_opacity = NULL
)
}
\arguments{
\item{input}{A data.frame or a list}

\item{var}{(character) A single variable to map colors, symbols,
and/or sizes to}

\item{var_col}{(character) A single variable to map colors to.}

\item{var_sym}{(character) A single variable to map symbols to.}

\item{var_size}{(character) A single variable to map size to.}

\item{var_stroke}{(character) A single variable to map stroke to.}

\item{var_stroke_width}{(character) A single variable to map stroke
width to.}

\item{var_stroke_opacity}{(character) A single variable to map stroke
opacity to.}

\item{var_fill}{(character) A single variable to map fill to.}

\item{var_fill_opacity}{(character) A single variable to map fill opacity to}

\item{color}{(character) Valid RGB hex color. Assigned to the variable
\code{marker-color}}

\item{symbol}{(character) An icon ID from the Maki project
https://labs.mapbox.com/maki-icons/
or a single alphanumeric character (a-z or 0-9). Assigned to the variable
\code{marker-symbol}}

\item{size}{(character) One of 'small', 'medium', or 'large'. Assigned
to the variable \code{marker-size}}

\item{stroke}{(character) Color of a polygon edge or line (RGB). Assigned
to the variable \code{stroke}}

\item{stroke_width}{(numeric) Width of a polygon edge or line (number > 0).
Assigned  to the variable \code{stroke-width}}

\item{stroke_opacity}{(numeric) Opacity of a polygon edge or line
(0.0 - 1.0). Assigned to the variable \code{stroke-opacity}}

\item{fill}{(character) The color of the interior of a polygon (GRB).
Assigned to the variable \code{fill}}

\item{fill_opacity}{(character) The opacity of the interior of a polygon
(0.0-1.0). Assigned to the variable \code{fill-opacity}}
}
\description{
This helps you add styling following the Simplestyle Spec. See Details
}
\details{
The parameters color, symbol, size, stroke, stroke_width,
stroke_opacity, fill, and fill_opacity expect a vector of size 1 (recycled),
or exact length of vector being applied to in your input data.

This function helps add styling data to a list or data.frame following the
Simplestyle Spec
(https://github.com/mapbox/simplestyle-spec/tree/master/1.1.0),
used by MapBox and GitHub Gists (that renders geoJSON/topoJSON
as interactive maps).

There are a few other style variables, but deal with polygons

GitHub has a nice help article on geoJSON files
https://help.github.com/articles/mapping-geojson-files-on-github/

Please do get in touch if you think anything should change in this
function.
}
\examples{
\dontrun{
## from data.frames - point data
library("RColorBrewer")
smalluscities <-
   subset(us_cities, country.etc == 'OR' | country.etc == 'NY' | country.etc == 'CA')

### Just color
geojson_style(smalluscities, var = 'country.etc',
   color=brewer.pal(length(unique(smalluscities$country.etc)), "Blues"))
### Just size
geojson_style(smalluscities, var = 'country.etc', size=c('small','medium','large'))
### Color and size
geojson_style(smalluscities, var = 'country.etc',
   color=brewer.pal(length(unique(smalluscities$country.etc)), "Blues"),
   size=c('small','medium','large'))

## from lists - point data
mylist <- list(list(latitude=30, longitude=120, state="US"),
               list(latitude=32, longitude=130, state="OR"),
               list(latitude=38, longitude=125, state="NY"),
               list(latitude=40, longitude=128, state="VT"))
# just color
geojson_style(mylist, var = 'state',
   color=brewer.pal(length(unique(sapply(mylist, '[[', 'state'))), "Blues"))
# color and size
geojson_style(mylist, var = 'state',
   color=brewer.pal(length(unique(sapply(mylist, '[[', 'state'))), "Blues"),
   size=c('small','medium','large','large'))
# color, size, and symbol
geojson_style(mylist, var = 'state',
   color=brewer.pal(length(unique(sapply(mylist, '[[', 'state'))), "Blues"),
   size=c('small','medium','large','large'),
   symbol="zoo")
# stroke, fill
geojson_style(mylist, var = 'state',
   stroke=brewer.pal(length(unique(sapply(mylist, '[[', 'state'))), "Blues"),
   fill=brewer.pal(length(unique(sapply(mylist, '[[', 'state'))), "Greens"))

# from data.frame - polygon data
smallstates <- states[states$group \%in\% 1:3, ]
head(smallstates)
geojson_style(smallstates, var = 'group',
   stroke = brewer.pal(length(unique(smallstates$group)), "Blues"),
   stroke_width = c(1, 2, 3),
   fill = brewer.pal(length(unique(smallstates$group)), "Greens"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_json.R
\name{geojson_json}
\alias{geojson_json}
\title{Convert many input types with spatial data to geojson specified as a json
string}
\usage{
geojson_json(
  input,
  lat = NULL,
  lon = NULL,
  group = NULL,
  geometry = "point",
  type = "FeatureCollection",
  convert_wgs84 = FALSE,
  crs = NULL,
  precision = NULL,
  ...
)
}
\arguments{
\item{input}{Input list, data.frame, spatial class, or sf class. Inputs can
also be dplyr \code{tbl_df} class since it inherits from \code{data.frame}.}

\item{lat}{(character) Latitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{lon}{(character) Longitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{group}{(character) A grouping variable to perform grouping for
polygons - doesn't apply for points}

\item{geometry}{(character) One of point (Default) or polygon.}

\item{type}{(character) The type of collection. One of 'auto' (default
for 'sf' objects), 'FeatureCollection' (default for everything else), or
'GeometryCollection'. "skip" skips the coercion with package \pkg{geojson}
functions; skipping can save significant run time on larger geojson
objects. \code{Spatial} objects can only accept "FeatureCollection" or "skip".
"skip" is not available as an option for \code{numeric}, \code{list},
and \code{data.frame} classes}

\item{convert_wgs84}{Should the input be converted to the
standard CRS system for GeoJSON (https://tools.ietf.org/html/rfc7946)
(geographic coordinate reference system, using
the WGS84 datum, with longitude and latitude units of decimal degrees;
EPSG: 4326). Default is \code{FALSE} though this may change in a future
package version. This will only work for \code{sf} or \code{Spatial}
objects with a CRS already defined. If one is not defined but you know
what it is, you may define it in the \code{crs} argument below.}

\item{crs}{The CRS of the input if it is not already defined. This can be
an epsg code as a four or five digit integer or a valid proj4 string.
This argument will be ignored if \code{convert_wgs84} is \code{FALSE} or
the object already has a CRS.}

\item{precision}{(integer) desired number of decimal places for coordinates.
Using fewer decimal places decreases object sizes (at the
cost of precision). This changes the underlying precision stored in the
data. \verb{options(digits = <some number>)} changes the maximum number of
digits displayed (to find out what yours is set at see
\code{getOption("digits")}); the value of this parameter will change what's
displayed in your console up to the value of \code{getOption("digits")}.
See Precision section for more.}

\item{...}{Further args passed on to internal functions. For Spatial*
classes, it is passed through to
\code{\link[sf:st_write]{sf::st_write()}}. For sf classes, data.frames, lists, numerics,
and geo_lists, it is passed through to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}
}
\value{
An object of class \code{geo_json} (and \code{json})
}
\description{
Convert many input types with spatial data to geojson specified as a json
string
}
\details{
This function creates a geojson structure as a json character
string; it does not write a file - see \code{\link[=geojson_write]{geojson_write()}} for that

Note that all sp class objects will output as \code{FeatureCollection}
objects, while other classes (numeric, list, data.frame) can be output as
\code{FeatureCollection} or \code{GeometryCollection} objects. We're working
on allowing \code{GeometryCollection} option for sp class objects.

Also note that with sp classes we do make a round-trip, using
\code{\link[sf:st_write]{sf::st_write()}} to write GeoJSON to disk, then read it back
in. This is fast and we don't have to think about it too much, but this
disk round-trip is not ideal.

For sf classes (sf, sfc, sfg), the following conversions are made:
\itemize{
\item sfg: the appropriate geometry \verb{Point, LineString, Polygon,  MultiPoint, MultiLineString, MultiPolygon, GeometryCollection}
\item sfc: \code{GeometryCollection}, unless the sfc is length 1, then
the geometry as above
\item sf: \code{FeatureCollection}
}
}
\section{Precision}{

Precision is handled in different ways depending on the class.

The \code{digits} parameter of \code{jsonlite::toJSON} controls precision for classes
\code{numeric}, \code{list}, \code{data.frame}, and \code{geo_list}.

For \code{sp} classes, precision is controlled by \code{sf::st_write}, being passed
down through \code{\link[=geojson_write]{geojson_write()}}, then through internal function
\code{write_geojson()}, then another internal function \code{write_ogr_sf()}

For \code{sf} classes, precision isn't quite working yet.
}

\examples{
\dontrun{
# From a numeric vector of length 2, making a point type
geojson_json(c(-99.74134244,32.451323223))
geojson_json(c(-99.74134244,32.451323223))[[1]]
geojson_json(c(-99.74134244,32.451323223), precision=2)[[1]]
geojson_json(c(-99.74,32.45), type = "GeometryCollection")

## polygon type
### this requires numeric class input, so inputting a list will dispatch
### on the list method
poly <- c(c(-114.345703125,39.436192999314095),
          c(-114.345703125,43.45291889355468),
          c(-106.61132812499999,43.45291889355468),
          c(-106.61132812499999,39.436192999314095),
          c(-114.345703125,39.436192999314095))
geojson_json(poly, geometry = "polygon")

# Lists
## From a list of numeric vectors to a polygon
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0), c(100.0,1.0), 
c(100.0,0.0))
geojson_json(vecs, geometry="polygon")

## from a named list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
geojson_json(mylist, lat='latitude', lon='longitude')

# From a data.frame to points
geojson_json(us_cities[1:2,], lat='lat', lon='long')
geojson_json(us_cities[1:2,], lat='lat', lon='long',
   type="GeometryCollection")

# from data.frame to polygons
head(states)
## make list for input to e.g., rMaps
geojson_json(states[1:351, ], lat='lat', lon='long', geometry="polygon", 
group='group')

# from a geo_list
a <- geojson_list(us_cities[1:2,], lat='lat', lon='long')
geojson_json(a)

# sp classes

## From SpatialPolygons class
library('sp')
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
   c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
   c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
geojson_json(sp_poly)

## Another SpatialPolygons
library("sp")
library("rgeos")
pt <- SpatialPoints(coordinates(list(x = 0, y = 0)), 
 CRS("+proj=longlat +datum=WGS84"))
## transfrom to web mercator becuase geos needs project coords
crs <- gsub("\n", "", 
  paste0("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0
  +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs", collapse = ""))
pt <- spTransform(pt, CRS(crs))
## buffer
pt <- gBuffer(pt, width = 100)
pt <- spTransform(pt, CRS("+proj=longlat +datum=WGS84"))
geojson_json(pt)

## data.frame to geojson
geojson_write(us_cities[1:2,], lat='lat', lon='long') \%>\% as.json

# From SpatialPoints class
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPoints(cbind(x,y))
geojson_json(s)

## From SpatialPointsDataFrame class
s <- SpatialPointsDataFrame(cbind(x,y), mtcars[1:5,])
geojson_json(s)

## From SpatialLines class
library("sp")
c1 <- cbind(c(1,2,3), c(3,2,2))
c2 <- cbind(c1[,1]+.05,c1[,2]+.05)
c3 <- cbind(c(1,2,3),c(1,1.5,1))
L1 <- Line(c1)
L2 <- Line(c2)
L3 <- Line(c3)
Ls1 <- Lines(list(L1), ID = "a")
Ls2 <- Lines(list(L2, L3), ID = "b")
sl1 <- SpatialLines(list(Ls1))
sl12 <- SpatialLines(list(Ls1, Ls2))
geojson_json(sl1)
geojson_json(sl12)

## From SpatialLinesDataFrame class
dat <- data.frame(X = c("Blue", "Green"),
                 Y = c("Train", "Plane"),
                 Z = c("Road", "River"), row.names = c("a", "b"))
sldf <- SpatialLinesDataFrame(sl12, dat)
geojson_json(sldf)
geojson_json(sldf)

## From SpatialGrid
x <- GridTopology(c(0,0), c(1,1), c(5,5))
y <- SpatialGrid(x)
geojson_json(y)

## From SpatialGridDataFrame
sgdim <- c(3,4)
sg <- SpatialGrid(GridTopology(rep(0,2), rep(10,2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
geojson_json(sgdf)

# From SpatialRings
library("rgeos")
r1 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="1")
r2 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="2")
r1r2 <- SpatialRings(list(r1, r2))
geojson_json(r1r2)

# From SpatialRingsDataFrame
dat <- data.frame(id = c(1,2), value = 3:4)
r1r2df <- SpatialRingsDataFrame(r1r2, data = dat)
geojson_json(r1r2df)

# From SpatialPixels
library("sp")
pixels <- suppressWarnings(
 SpatialPixels(SpatialPoints(us_cities[c("long", "lat")])))
summary(pixels)
geojson_json(pixels)

# From SpatialPixelsDataFrame
library("sp")
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")],
 data = canada_cities)
)
geojson_json(pixelsdf)

# From SpatialCollections
library("sp")
library("rgeos")
pts <- SpatialPoints(cbind(c(1,2,3,4,5), c(3,2,5,1,4)))
poly1 <- Polygons(
 list(Polygon(cbind(c(-100,-90,-85,-100), c(40,50,45,40)))), "1")
poly2 <- Polygons(
 list(Polygon(cbind(c(-90,-80,-75,-90), c(30,40,35,30)))), "2")
poly <- SpatialPolygons(list(poly1, poly2), 1:2)
dat <- SpatialCollections(pts, polygons = poly)
geojson_json(dat)

# From sf classes:
if (require(sf)) {
## sfg (a single simple features geometry)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  poly <- rbind(c(1,1), c(1,2), c(2,2), c(1,1))
  poly_sfg <-st_polygon(list(p1))
  geojson_json(poly_sfg)

## sfc (a collection of geometries)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  geojson_json(poly_sfc)

## sf (collection of geometries with attributes)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  poly_sf <- st_sf(foo = c("a", "b"), bar = 1:2, poly_sfc)
  geojson_json(poly_sf)
}

## Pretty print a json string
geojson_json(c(-99.74,32.45))
geojson_json(c(-99.74,32.45)) \%>\% pretty

# skipping the pretty geojson class coercion with the geojson pkg
if (require(sf)) {
  library(sf)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  geojson_json(poly_sfc)
  geojson_json(poly_sfc, type = "skip")
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postgis.R
\name{postgis}
\alias{postgis}
\title{PostGIS setup}
\description{
\code{\link[=geojson_read]{geojson_read()}} allows you to get data out of a PostgreSQL
database set up with PostGIS. Below are steps for setting up data
that we can at the end query with \code{\link[=geojson_read]{geojson_read()}}
}
\details{
If you don't already have PostgreSQL or PostGIS:
\itemize{
\item PostgreSQL installation: https://www.postgresql.org/download/
\item PostGIS installation: https://postgis.net/install/
}

Once you have both of those installed, you can proceed below.
}
\examples{
\dontrun{
if (requireNamespace("DBI") && requireNamespace("RPostgres")) {
library("DBI")
library("RPostgres")

# Create connection
conn <- tryCatch(dbConnect(RPostgres::Postgres()), error = function(e) e)
if (inherits(conn, "PqConnection")) {

# Create database
dbSendQuery(conn, "CREATE DATABASE postgistest")

# New connection to the created database
conn <- dbConnect(RPostgres::Postgres(), dbname = 'postgistest')

# Initialize PostGIS in Postgres
dbSendQuery(conn, "CREATE EXTENSION postgis")
dbSendQuery(conn, "SELECT postgis_full_version()")

# Create table
dbSendQuery(conn, "CREATE TABLE locations(loc_id integer primary key
   , loc_name varchar(70), geog geography(POINT) );")

# Insert data
dbSendQuery(conn, "INSERT INTO locations(loc_id, loc_name, geog)
 VALUES (1, 'Waltham, MA', ST_GeogFromText('POINT(42.40047 -71.2577)') )
   , (2, 'Manchester, NH', ST_GeogFromText('POINT(42.99019 -71.46259)') )
   , (3, 'TI Blvd, TX', ST_GeogFromText('POINT(-96.75724 32.90977)') );")


# Get data (notice warnings of unknown field type for geog)
dbGetQuery(conn, "SELECT * from locations")


# Once you're setup, use geojson_read()
conn <- dbConnect(RPostgres::Postgres(), dbname = 'postgistest')
state <- "SELECT row_to_json(fc)
 FROM (SELECT 'FeatureCollection' As type, array_to_json(array_agg(f)) As features
 FROM (SELECT 'Feature' As type
    , ST_AsGeoJSON(lg.geog)::json As geometry
    , row_to_json((SELECT l FROM (SELECT loc_id, loc_name) As l
      )) As properties
   FROM locations As lg   ) As f )  As fc;"
json <- geojson_read(conn, query = state, what = "json")

## map the geojson with map_leaf()
map_leaf(json)

}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lint.R
\name{lint-defunct}
\alias{lint-defunct}
\alias{lint}
\title{Lint geojson}
\usage{
lint(...)
}
\arguments{
\item{...}{ignored}
}
\description{
Lint geojson
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_sp.R
\name{geojson_sp}
\alias{geojson_sp}
\title{Convert objects to spatial classes}
\usage{
geojson_sp(x, disambiguateFIDs = FALSE, stringsAsFactors = FALSE, ...)
}
\arguments{
\item{x}{Object of class \code{geo_list}, \code{geo_json}, string, or json}

\item{disambiguateFIDs}{Ignored, and will be removed in a future version.
Previously was passed to \code{rgdal::readOGR()}, which is no longer used.}

\item{stringsAsFactors}{Convert strings to Factors? Default \code{FALSE}.}

\item{...}{Further args passed on to \code{\link[sf:st_read]{sf::st_read()}}}
}
\value{
A spatial class object, see Details.
}
\description{
Convert objects to spatial classes
}
\details{
The spatial class object returned will depend on the input GeoJSON.
Sometimes you will get back a \code{SpatialPoints} class, and sometimes a
\code{SpatialPolygonsDataFrame} class, etc., depending on what the
structure of the GeoJSON.

The reading and writing of the CRS to/from geojson is inconsistent. You can
directly set the CRS by passing a valid PROJ4 string or epsg code to the crs
argument in \code{\link[sf:st_read]{sf::st_read()}}
}
\examples{
\dontrun{
library(sp)

# geo_list ------------------
## From a numeric vector of length 2 to a point
vec <- c(-99.74,32.45)
geojson_list(vec) \%>\% geojson_sp

## Lists
## From a list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
geojson_list(mylist) \%>\% geojson_sp
geojson_list(mylist) \%>\% geojson_sp \%>\% plot

## From a list of numeric vectors to a polygon
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0), c(100.0,1.0), c(100.0,0.0))
geojson_list(vecs, geometry="polygon") \%>\% geojson_sp
geojson_list(vecs, geometry="polygon") \%>\% geojson_sp \%>\% plot

# geo_json ------------------
## from point
geojson_json(c(-99.74,32.45)) \%>\% geojson_sp
geojson_json(c(-99.74,32.45)) \%>\% geojson_sp \%>\% plot

# from featurecollectino of points
geojson_json(us_cities[1:2,], lat='lat', lon='long') \%>\% geojson_sp
geojson_json(us_cities[1:2,], lat='lat', lon='long') \%>\% geojson_sp \%>\% plot

# Set the CRS via the crs argument
geojson_json(us_cities[1:2,], lat='lat', lon='long') \%>\%
  geojson_sp(crs = "+init=epsg:4326")

# json ----------------------
x <- geojson_json(us_cities[1:2,], lat='lat', lon='long')
geojson_sp(x)

# character string ----------------------
x <- unclass(geojson_json(c(-99.74,32.45)))
geojson_sp(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_sf.R
\name{geojson_sf}
\alias{geojson_sf}
\title{Convert objects to an sf class}
\usage{
geojson_sf(x, stringsAsFactors = FALSE, ...)
}
\arguments{
\item{x}{Object of class \code{geo_list}, \code{geo_json}, string, or json}

\item{stringsAsFactors}{Convert strings to Factors? Default \code{FALSE}.}

\item{...}{Further args passed on to \code{\link[sf:st_read]{sf::st_read()}}}
}
\value{
An sf class object, see Details.
}
\description{
Convert objects to an sf class
}
\details{
The type of sf object returned will depend on the input GeoJSON.
Sometimes you will get back a \code{POINTS} class, and sometimes a
\code{POLYGON} class, etc., depending on what the structure of the GeoJSON.

The reading and writing of the CRS to/from geojson is inconsistent. You can
directly set the CRS by passing a valid PROJ4 string or epsg code to the crs
argument in \code{\link[sf:st_read]{sf::st_read()}}
}
\examples{
\dontrun{
library(sf)

# geo_list ------------------
## From a numeric vector of length 2 to a point
vec <- c(-99.74,32.45)
geojson_list(vec) \%>\% geojson_sf

## Lists
## From a list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
geojson_list(mylist) \%>\% geojson_sf
geojson_list(mylist) \%>\% geojson_sf \%>\% plot

## From a list of numeric vectors to a polygon
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0), c(100.0,1.0), c(100.0,0.0))
geojson_list(vecs, geometry="polygon") \%>\% geojson_sf
geojson_list(vecs, geometry="polygon") \%>\% geojson_sf \%>\% plot

# geo_json ------------------
## from point
geojson_json(c(-99.74,32.45)) \%>\% geojson_sf
geojson_json(c(-99.74,32.45)) \%>\% geojson_sf \%>\% plot

# from featurecollectino of points
geojson_json(us_cities[1:2,], lat='lat', lon='long') \%>\% geojson_sf
geojson_json(us_cities[1:2,], lat='lat', lon='long') \%>\% geojson_sf \%>\% plot

# Set the CRS via the crs argument
geojson_json(us_cities[1:2,], lat='lat', lon='long') \%>\% geojson_sf(crs = "+init=epsg:4326")

# json ----------------------
x <- geojson_json(us_cities[1:2,], lat='lat', lon='long')
geojson_sf(x)

# character string ----------------------
x <- unclass(geojson_json(c(-99.74,32.45)))
geojson_sf(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topojson_write.R
\name{topojson_write}
\alias{topojson_write}
\title{Write TopoJSON from various inputs}
\usage{
topojson_write(
  input,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "myfile.topojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL,
  object_name = "foo",
  quantization = 0,
  ...
)
}
\arguments{
\item{input}{Input list, data.frame, spatial class, or sf class.
Inputs can  also be dplyr \code{tbl_df} class since it inherits
from \code{data.frame}}

\item{lat}{(character) Latitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{lon}{(character) Longitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{geometry}{(character) One of point (Default) or polygon.}

\item{group}{(character) A grouping variable to perform grouping for
polygons - doesn't apply for points}

\item{file}{(character) A path and file name (e.g., myfile), with the
\code{.geojson} file extension. Default writes to current working
directory.}

\item{overwrite}{(logical) Overwrite the file given in \code{file} with
\code{input}. Default: \code{TRUE}. If this param is \code{FALSE} and
the file already exists, we stop with error message.}

\item{precision}{desired number of decimal places for the coordinates in the
geojson file. Using fewer decimal places can decrease file sizes (at the
cost of precision).}

\item{convert_wgs84}{Should the input be converted to the
standard CRS for GeoJSON (https://tools.ietf.org/html/rfc7946)
(geographic coordinate reference
system, using the WGS84 datum, with longitude and latitude units of decimal
degrees; EPSG: 4326). Default is \code{FALSE} though this may change in a
future package version. This will only work for \code{sf} or \code{Spatial}
objects with a CRS already defined. If one is not defined but you know what
it is, you may define it in the \code{crs} argument below.}

\item{crs}{The CRS of the input if it is not already defined. This can be
an epsg code as a four or five digit integer or a valid proj4 string. This
argument will be ignored if \code{convert_wgs84} is \code{FALSE} or the
object already has a CRS.}

\item{object_name}{(character) name to give to the TopoJSON object created.
Default: "foo"}

\item{quantization}{(numeric) quantization parameter, use this to
quantize geometry prior to computing topology. Typical values are powers of
ten (\code{1e4}, \code{1e5}, ...), default is \code{0} to not perform quantization.
For more information about quantization, see this by Mike Bostock
https://stackoverflow.com/questions/18900022/topojson-quantization-vs-simplification/18921214#18921214}

\item{...}{Further args passed on to internal functions. For Spatial*
classes, data.frames,
regular lists, and numerics, it is passed through to
\code{\link[sf:st_write]{sf::st_write()}}. For sf classes,
geo_lists and json classes, it is passed through to
\code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}.}
}
\value{
A \code{topojson_write} class, with two elements:
\itemize{
\item path: path to the file with the TopoJSON
\item type: type of object the TopoJSON came from, e.g., SpatialPoints
}
}
\description{
\code{topojson_write()} is temporarily defunct; check back later
}
\details{
Under the hood we simply wrap \code{\link[=geojson_write]{geojson_write()}}, then
take the GeoJSON output of that operation, then convert to TopoJSON with
\code{\link[=geo2topo]{geo2topo()}}, then write to disk.

Unfortunately, this process requires a number of round trips to disk, so
speed ups will hopefully come soon.

Any intermediate geojson files are cleaned up (deleted).
}
\seealso{
\code{\link[=geojson_write]{geojson_write()}}, \code{\link[=topojson_read]{topojson_read()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centroid.R
\name{centroid}
\alias{centroid}
\title{Get centroid for a geo_list}
\usage{
centroid(x, ...)
}
\arguments{
\item{x}{An object of class geo_list}

\item{...}{Ignored}
}
\value{
A vector of the form longitude, latitude
}
\description{
Get centroid for a geo_list
}
\examples{
# numeric
vec <- c(-99.74,32.45)
x <- geojson_list(vec)
centroid(x)

# list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
x <- geojson_list(mylist)
centroid(x)

# data.frame
x <- geojson_list(states[1:20,])
centroid(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapleaflet.R
\name{map_leaf}
\alias{map_leaf}
\title{Make an interactive map locally}
\usage{
map_leaf(input, lat = NULL, lon = NULL, basemap = "Stamen.Toner", ...)
}
\arguments{
\item{input}{Input object}

\item{lat}{Name of latitude variable}

\item{lon}{Name of longitude variable}

\item{basemap}{Basemap to use. See \code{leaflet::addProviderTiles}.
Default: \code{Stamen.Toner}}

\item{...}{Further arguments passed on to \code{leaflet::addPolygons},
\code{leaflet::addMarkers}, \code{leaflet::addGeoJSON}, or \code{leaflet::addPolylines}}
}
\description{
Make an interactive map locally
}
\examples{
\dontrun{
# We'll need leaflet below
library("leaflet")

# From file
file <- "myfile.geojson"
geojson_write(us_cities[1:20, ], lat='lat', lon='long', file = file)
map_leaf(as.location(file))

# From SpatialPoints class
library("sp")
x <- c(1,2,3,4,20)
y <- c(3,2,5,3,4)
s <- SpatialPoints(cbind(x,y))
map_leaf(s)

# from SpatialPointsDataFrame class
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPointsDataFrame(cbind(x,y), mtcars[1:5,])
map_leaf(s)

# from SpatialPolygons class
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
   c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
   c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
map_leaf(sp_poly)

# From SpatialPolygonsDataFrame class
sp_polydf <- as(sp_poly, "SpatialPolygonsDataFrame")
map_leaf(sp_poly)

# From SpatialLines class
c1 <- cbind(c(1,2,3), c(3,2,2))
c2 <- cbind(c1[,1]+.05,c1[,2]+.05)
c3 <- cbind(c(1,2,3),c(1,1.5,1))
L1 <- Line(c1)
L2 <- Line(c2)
L3 <- Line(c3)
Ls1 <- Lines(list(L1), ID = "a")
Ls2 <- Lines(list(L2, L3), ID = "b")
sl1 <- SpatialLines(list(Ls1))
sl12 <- SpatialLines(list(Ls1, Ls2))
map_leaf(sl1)
map_leaf(sl12)

# From SpatialLinesDataFrame class
dat <- data.frame(X = c("Blue", "Green"),
                 Y = c("Train", "Plane"),
                 Z = c("Road", "River"), row.names = c("a", "b"))
sldf <- SpatialLinesDataFrame(sl12, dat)
map_leaf(sldf)

# From SpatialGrid
x <- GridTopology(c(0,0), c(1,1), c(5,5))
y <- SpatialGrid(x)
map_leaf(y)

# From SpatialGridDataFrame
sgdim <- c(3,4)
sg <- SpatialGrid(GridTopology(rep(0,2), rep(10,2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
map_leaf(sgdf)

# from data.frame
map_leaf(us_cities)

## another example
head(states)
map_leaf(states[1:351, ])

## From a named list
mylist <- list(list(lat=30, long=120, marker="red"),
               list(lat=30, long=130, marker="blue"))
map_leaf(mylist, lat="lat", lon="long")

## From an unnamed list
poly <- list(c(-114.345703125,39.436192999314095),
             c(-114.345703125,43.45291889355468),
             c(-106.61132812499999,43.45291889355468),
             c(-106.61132812499999,39.436192999314095),
             c(-114.345703125,39.436192999314095))
map_leaf(poly)
## NOTE: Polygons from lists aren't supported yet

# From a json object
map_leaf(geojson_json(c(-99.74, 32.45)))
map_leaf(geojson_json(c(-119, 45)))
map_leaf(geojson_json(c(-99.74, 32.45)))
## another example
map_leaf(geojson_json(us_cities[1:10,], lat='lat', lon='long'))

# From a geo_list object
(res <- geojson_list(us_cities[1:2,], lat='lat', lon='long'))
map_leaf(res)

# From SpatialPixels
pixels <- suppressWarnings(SpatialPixels(SpatialPoints(us_cities[c("long", "lat")])))
summary(pixels)
map_leaf(pixels)

# From SpatialPixelsDataFrame
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")], data = canada_cities)
)
map_leaf(pixelsdf)

# From SpatialRings
library("rgeos")
r1 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="1")
r2 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="2")
r1r2 <- SpatialRings(list(r1, r2))
map_leaf(r1r2)

# From SpatialRingsDataFrame
dat <- data.frame(id = c(1,2), value = 3:4)
r1r2df <- SpatialRingsDataFrame(r1r2, data = dat)
map_leaf(r1r2df)

# basemap toggling ------------------------
map_leaf(us_cities, basemap = "Acetate.terrain")
map_leaf(us_cities, basemap = "CartoDB.Positron")
map_leaf(us_cities, basemap = "OpenTopoMap")

# leaflet options ------------------------
map_leaf(us_cities) \%>\%
   addPopups(-122.327298, 47.597131, "foo bar", options = popupOptions(closeButton = FALSE))

####### not working yet
# From a numeric vector
## of length 2 to a point
## vec <- c(-99.74,32.45)
## map_leaf(vec)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_write.r
\name{geojson_write}
\alias{geojson_write}
\title{Convert many input types with spatial data to a geojson file}
\usage{
geojson_write(
  input,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "myfile.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL,
  ...
)
}
\arguments{
\item{input}{Input list, data.frame, spatial class, or sf class.
Inputs can  also be dplyr \code{tbl_df} class since it inherits
from \code{data.frame}}

\item{lat}{(character) Latitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{lon}{(character) Longitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{geometry}{(character) One of point (Default) or polygon.}

\item{group}{(character) A grouping variable to perform grouping for
polygons - doesn't apply for points}

\item{file}{(character) A path and file name (e.g., myfile), with the
\code{.geojson} file extension. Default writes to current working
directory.}

\item{overwrite}{(logical) Overwrite the file given in \code{file} with
\code{input}. Default: \code{TRUE}. If this param is \code{FALSE} and
the file already exists, we stop with error message.}

\item{precision}{desired number of decimal places for the coordinates in the
geojson file. Using fewer decimal places can decrease file sizes (at the
cost of precision).}

\item{convert_wgs84}{Should the input be converted to the
standard CRS for GeoJSON (https://tools.ietf.org/html/rfc7946)
(geographic coordinate reference
system, using the WGS84 datum, with longitude and latitude units of decimal
degrees; EPSG: 4326). Default is \code{FALSE} though this may change in a
future package version. This will only work for \code{sf} or \code{Spatial}
objects with a CRS already defined. If one is not defined but you know what
it is, you may define it in the \code{crs} argument below.}

\item{crs}{The CRS of the input if it is not already defined. This can be
an epsg code as a four or five digit integer or a valid proj4 string. This
argument will be ignored if \code{convert_wgs84} is \code{FALSE} or the
object already has a CRS.}

\item{...}{Further args passed on to internal functions. For Spatial*
classes, data.frames,
regular lists, and numerics, it is passed through to
\code{\link[sf:st_write]{sf::st_write()}}. For sf classes,
geo_lists and json classes, it is passed through to
\code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}.}
}
\value{
A \code{geojson_write} class, with two elements:
\itemize{
\item path: path to the file with the GeoJSON
\item type: type of object the GeoJSON came from, e.g., SpatialPoints
}
}
\description{
Convert many input types with spatial data to a geojson file
}
\examples{
\dontrun{
# From a data.frame
## to points
geojson_write(us_cities[1:2,], lat='lat', lon='long')

## to polygons
head(states)
geojson_write(input=states, lat='lat', lon='long',
  geometry='polygon', group="group")

## partial states dataset to points (defaults to points)
geojson_write(input=states, lat='lat', lon='long')

## Lists
### list of numeric pairs
poly <- list(c(-114.345703125,39.436192999314095),
          c(-114.345703125,43.45291889355468),
          c(-106.61132812499999,43.45291889355468),
          c(-106.61132812499999,39.436192999314095),
          c(-114.345703125,39.436192999314095))
geojson_write(poly, geometry = "polygon")

### named list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
geojson_write(mylist)

# From a numeric vector of length 2
## Expected order is lon, lat
vec <- c(-99.74, 32.45)
geojson_write(vec)

## polygon from a series of numeric pairs
### this requires numeric class input, so inputting a list will
### dispatch on the list method
poly <- c(c(-114.345703125,39.436192999314095),
          c(-114.345703125,43.45291889355468),
          c(-106.61132812499999,43.45291889355468),
          c(-106.61132812499999,39.436192999314095),
          c(-114.345703125,39.436192999314095))
geojson_write(poly, geometry = "polygon")

# Write output of geojson_list to file
res <- geojson_list(us_cities[1:2,], lat='lat', lon='long')
class(res)
geojson_write(res)

# Write output of geojson_json to file
res <- geojson_json(us_cities[1:2,], lat='lat', lon='long')
class(res)
geojson_write(res)

# From SpatialPolygons class
library('sp')
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
   c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
   c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
geojson_write(sp_poly)

# From SpatialPolygonsDataFrame class
sp_polydf <- as(sp_poly, "SpatialPolygonsDataFrame")
geojson_write(input = sp_polydf)

# From SpatialGrid
x <- GridTopology(c(0,0), c(1,1), c(5,5))
y <- SpatialGrid(x)
geojson_write(y)

# From SpatialGridDataFrame
sgdim <- c(3,4)
sg <- SpatialGrid(GridTopology(rep(0,2), rep(10,2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
geojson_write(sgdf)

# From SpatialRings
library(rgeos)
r1 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="1")
r2 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="2")
r1r2 <- SpatialRings(list(r1, r2))
geojson_write(r1r2)

# From SpatialRingsDataFrame
dat <- data.frame(id = c(1,2), value = 3:4)
r1r2df <- SpatialRingsDataFrame(r1r2, data = dat)
geojson_write(r1r2df)

# From SpatialPixels
library("sp")
pixels <- suppressWarnings(SpatialPixels(SpatialPoints(us_cities[c("long", "lat")])))
summary(pixels)
geojson_write(pixels)

# From SpatialPixelsDataFrame
library("sp")
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")], data = canada_cities)
)
geojson_write(pixelsdf)

# From SpatialCollections
library("sp")
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100), c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90), c(30,40,35,30)))), "2")
poly <- SpatialPolygons(list(poly1, poly2), 1:2)
coordinates(us_cities) <- ~long+lat
dat <- SpatialCollections(points = us_cities, polygons = poly)
geojson_write(dat)

# From sf classes:
if (require(sf)) {
  file <- system.file("examples", "feature_collection.geojson", package = "geojsonio")
  sf_fc <- st_read(file, quiet = TRUE)
  geojson_write(sf_fc)
}
}
}
\seealso{
\code{\link[=geojson_list]{geojson_list()}}, \code{\link[=geojson_json]{geojson_json()}}, \code{\link[=topojson_write]{topojson_write()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojsonio-package.r
\docType{data}
\name{states}
\alias{states}
\title{This is the same data set from the ggplot2 library}
\description{
This is a data.frame with "long", "lat", "group", "order", "region", and
"subregion" columns specifying polygons for each US state.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_list.R
\name{geojson_list}
\alias{geojson_list}
\title{Convert many input types with spatial data to geojson specified as a list}
\usage{
geojson_list(
  input,
  lat = NULL,
  lon = NULL,
  group = NULL,
  geometry = "point",
  type = "FeatureCollection",
  convert_wgs84 = FALSE,
  crs = NULL,
  precision = NULL,
  ...
)
}
\arguments{
\item{input}{Input list, data.frame, spatial class, or sf class. Inputs can
also be dplyr \code{tbl_df} class since it inherits from \code{data.frame}}

\item{lat}{(character) Latitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{lon}{(character) Longitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{group}{(character) A grouping variable to perform grouping for
polygons - doesn't apply for points}

\item{geometry}{(character) One of point (Default) or polygon.}

\item{type}{(character) The type of collection. One of FeatureCollection
(default) or GeometryCollection.}

\item{convert_wgs84}{Should the input be converted to the
standard CRS for GeoJSON (https://tools.ietf.org/html/rfc7946)
(geographic coordinate reference system, using the WGS84 datum, with
longitude and latitude units of decimal degrees; EPSG: 4326).
Default is \code{FALSE} though this may change in a future package version.
This will only work for \code{sf} or \code{Spatial} objects with a CRS
already defined. If one is not defined but you know what it is, you
may define it in the \code{crs} argument below.}

\item{crs}{The CRS of the input if it is not already defined. This can
be an epsg code as a four or five digit integer or a valid proj4 string.
This argument will be ignored if \code{convert_wgs84} is \code{FALSE}
or the object already has a CRS.}

\item{precision}{(integer) desired number of decimal places for coordinates.
Only used with classes from \pkg{sp}\pkg{rgeos} classes; ignored for other
classes. Using fewer decimal places decreases object sizes (at the
cost of precision). This changes the underlying precision stored in the
data. \verb{options(digits = <some number>)} changes the maximum number of
digits displayed (to find out what yours is set at see
\code{getOption("digits")}); the value of this parameter will change what's
displayed in your console up to the value of \code{getOption("digits")}}

\item{...}{Ignored}
}
\description{
Convert many input types with spatial data to geojson specified as a list
}
\details{
This function creates a geojson structure as an R list; it does
not write a file - see \code{\link[=geojson_write]{geojson_write()}} for that.

Note that all sp class objects will output as \code{FeatureCollection} objects,
while other classes (numeric, list, data.frame) can be output as
\code{FeatureCollection} or \code{GeometryCollection} objects. We're working
on allowing \code{GeometryCollection} option for sp class objects.

Also note that with sp classes we do make a round-trip,
using \code{\link[sf:st_write]{sf::st_write()}} to write GeoJSON to disk, then read it back in.
This is fast and we don't have to think
about it too much, but this disk round-trip is not ideal.

For sf classes (sf, sfc, sfg), the following conversions are made:
\itemize{
\item sfg: the appropriate geometry \verb{Point, LineString, Polygon, MultiPoint,  MultiLineString, MultiPolygon, GeometryCollection}
\item sfc: \code{GeometryCollection}, unless the sfc is length 1, then the geometry
as above
\item sf: \code{FeatureCollection}
}

For \code{list} and \code{data.frame} objects, you don't have to pass in \code{lat} and
\code{lon} parameters if they are named appropriately (e.g., lat/latitude,
lon/long/longitude), as they will be auto-detected. If they can not be
found, the function will stop and warn you to specify the parameters
specifically.
}
\examples{
\dontrun{
# From a numeric vector of length 2 to a point
vec <- c(-99.74,32.45)
geojson_list(vec)

# Lists
## From a list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
geojson_list(mylist)

## From a list of numeric vectors to a polygon
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0),
  c(100.0,1.0), c(100.0,0.0))
geojson_list(vecs, geometry="polygon")

# from data.frame to points
(res <- geojson_list(us_cities[1:2,], lat='lat', lon='long'))
as.json(res)
## guess lat/long columns
geojson_list(us_cities[1:2,])
geojson_list(states[1:3,])
geojson_list(states[1:351,], geometry="polygon", group='group')
geojson_list(canada_cities[1:30,])
## a data.frame with columsn not named appropriately, but you can
## specify them
# dat <- data.frame(a = c(31, 41), b = c(-120, -110))
# geojson_list(dat)
# geojson_list(dat, lat="a", lon="b")

# from data.frame to polygons
head(states)
geojson_list(states[1:351, ], lat='lat', lon='long',
  geometry="polygon", group='group')

# From SpatialPolygons class
library('sp')
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
   c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
   c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
geojson_list(sp_poly)

# From SpatialPolygons class with precision agreement
x_coord <- c(-114.345703125, -114.345703125, -106.61132812499999,
  -106.61132812499999, -114.345703125)
y_coord <- c(39.436192999314095, 43.45291889355468, 43.45291889355468,
  39.436192999314095, 39.436192999314095)
coords <- cbind(x_coord, y_coord)
poly <- Polygon(coords)
polys <- Polygons(list(poly), 1)
sp_poly2 <- SpatialPolygons(list(polys))
geojson_list(sp_poly2, geometry = "polygon", precision = 4)
geojson_list(sp_poly2, geometry = "polygon", precision = 3)
geojson_list(sp_poly2, geometry = "polygon", precision = 2)

# From SpatialPoints class with precision
points <- SpatialPoints(cbind(x_coord,y_coord))
geojson_list(points)

# From SpatialPolygonsDataFrame class
sp_polydf <- as(sp_poly, "SpatialPolygonsDataFrame")
geojson_list(input = sp_polydf)

# From SpatialPoints class
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPoints(cbind(x,y))
geojson_list(s)

# From SpatialPointsDataFrame class
s <- SpatialPointsDataFrame(cbind(x,y), mtcars[1:5,])
geojson_list(s)

# From SpatialLines class
library("sp")
c1 <- cbind(c(1,2,3), c(3,2,2))
c2 <- cbind(c1[,1]+.05,c1[,2]+.05)
c3 <- cbind(c(1,2,3),c(1,1.5,1))
L1 <- Line(c1)
L2 <- Line(c2)
L3 <- Line(c3)
Ls1 <- Lines(list(L1), ID = "a")
Ls2 <- Lines(list(L2, L3), ID = "b")
sl1 <- SpatialLines(list(Ls1))
sl12 <- SpatialLines(list(Ls1, Ls2))
geojson_list(sl1)
geojson_list(sl12)
as.json(geojson_list(sl12))
as.json(geojson_list(sl12), pretty=TRUE)

# From SpatialLinesDataFrame class
dat <- data.frame(X = c("Blue", "Green"),
                 Y = c("Train", "Plane"),
                 Z = c("Road", "River"), row.names = c("a", "b"))
sldf <- SpatialLinesDataFrame(sl12, dat)
geojson_list(sldf)
as.json(geojson_list(sldf))
as.json(geojson_list(sldf), pretty=TRUE)

# From SpatialGrid
x <- GridTopology(c(0,0), c(1,1), c(5,5))
y <- SpatialGrid(x)
geojson_list(y)

# From SpatialGridDataFrame
sgdim <- c(3,4)
sg <- SpatialGrid(GridTopology(rep(0,2), rep(10,2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
geojson_list(sgdf)

# From SpatialRings
library("rgeos")
r1 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="1")
r2 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="2")
r1r2 <- SpatialRings(list(r1, r2))
geojson_list(r1r2)

# From SpatialRingsDataFrame
dat <- data.frame(id = c(1,2), value = 3:4)
r1r2df <- SpatialRingsDataFrame(r1r2, data = dat)
geojson_list(r1r2df)

# From SpatialPixels
library("sp")
pixels <- suppressWarnings(
  SpatialPixels(SpatialPoints(us_cities[c("long", "lat")])))
summary(pixels)
geojson_list(pixels)

# From SpatialPixelsDataFrame
library("sp")
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")],
 data = canada_cities)
)
geojson_list(pixelsdf)

# From SpatialCollections
library("sp")
poly1 <- Polygons(
  list(Polygon(cbind(c(-100,-90,-85,-100), c(40,50,45,40)))), "1")
poly2 <- Polygons(
  list(Polygon(cbind(c(-90,-80,-75,-90), c(30,40,35,30)))), "2")
poly <- SpatialPolygons(list(poly1, poly2), 1:2)
coordinates(us_cities) <- ~long+lat
dat <- SpatialCollections(points = us_cities, polygons = poly)
out <- geojson_list(dat)
out$SpatialPoints
out$SpatialPolygons
}

# From sf classes:
if (require(sf)) {
## sfg (a single simple features geometry)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  poly <- rbind(c(1,1), c(1,2), c(2,2), c(1,1))
  poly_sfg <-st_polygon(list(p1))
  geojson_list(poly_sfg)

## sfc (a collection of geometries)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  geojson_list(poly_sfc)

## sf (collection of geometries with attributes)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  poly_sf <- st_sf(foo = c("a", "b"), bar = 1:2, poly_sfc)
  geojson_list(poly_sf)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topojson_read.R
\name{topojson_read}
\alias{topojson_read}
\title{Read topojson from a local file or a URL}
\usage{
topojson_read(x, ...)
}
\arguments{
\item{x}{Path to a local file or a URL.}

\item{...}{Further args passed on to \code{\link[sf:st_read]{sf::st_read()}}. Can use any args
from \code{sf::st_read()} except \code{quiet}, which we have set as \code{quiet = TRUE}
internally already}
}
\value{
an object of class \code{sf}/\code{data.frame}
}
\description{
Read topojson from a local file or a URL
}
\details{
Returns a \code{sf} class, but you can easily and quickly get
this to geojson, see examples.

Note that this does not give you Topojson, but gives you a \code{sf}
class - which you can use then to turn it into geojson as a list or json
}
\examples{
\dontrun{
# From a file
file <- system.file("examples", "us_states.topojson", package = "geojsonio")
topojson_read(file)

# From a URL
url <- "https://raw.githubusercontent.com/shawnbot/d3-cartogram/master/data/us-states.topojson"
topojson_read(url)

# Use as.location first if you want
topojson_read(as.location(file))

# quickly convert to geojson as a list
file <- system.file("examples", "us_states.topojson", package = "geojsonio")
tmp <- topojson_read(file)
geojson_list(tmp)
geojson_json(tmp)

# pass on args
topojson_read(file, quiet = TRUE)
topojson_read(file, stringsAsFactors = FALSE)
}
}
\seealso{
\code{\link[=geojson_read]{geojson_read()}}, \code{\link[=topojson_write]{topojson_write()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topojson_list.R
\name{topojson_list}
\alias{topojson_list}
\title{Convert many input types with spatial data to TopoJSON
as a list}
\usage{
topojson_list(
  input,
  lat = NULL,
  lon = NULL,
  group = NULL,
  geometry = "point",
  type = "FeatureCollection",
  convert_wgs84 = FALSE,
  crs = NULL,
  object_name = "foo",
  quantization = 0,
  ...
)
}
\arguments{
\item{input}{Input list, data.frame, spatial class, or sf class. Inputs can
also be dplyr \code{tbl_df} class since it inherits from \code{data.frame}}

\item{lat}{(character) Latitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{lon}{(character) Longitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{group}{(character) A grouping variable to perform grouping for
polygons - doesn't apply for points}

\item{geometry}{(character) One of point (Default) or polygon.}

\item{type}{(character) The type of collection. One of FeatureCollection
(default) or GeometryCollection.}

\item{convert_wgs84}{Should the input be converted to the
standard CRS for GeoJSON (https://tools.ietf.org/html/rfc7946)
(geographic coordinate reference system, using the WGS84 datum, with
longitude and latitude units of decimal degrees; EPSG: 4326).
Default is \code{FALSE} though this may change in a future package version.
This will only work for \code{sf} or \code{Spatial} objects with a CRS
already defined. If one is not defined but you know what it is, you
may define it in the \code{crs} argument below.}

\item{crs}{The CRS of the input if it is not already defined. This can
be an epsg code as a four or five digit integer or a valid proj4 string.
This argument will be ignored if \code{convert_wgs84} is \code{FALSE}
or the object already has a CRS.}

\item{object_name}{(character) name to give to the TopoJSON object created.
Default: "foo"}

\item{quantization}{(numeric) quantization parameter, use this to
quantize geometry prior to computing topology. Typical values are powers of
ten (\code{1e4}, \code{1e5}, ...), default is \code{0} to not perform quantization.
For more information about quantization, see this by Mike Bostock
https://stackoverflow.com/questions/18900022/topojson-quantization-vs-simplification/18921214#18921214}

\item{...}{args passed down through \code{\link[=topojson_json]{topojson_json()}} to \code{\link[=geojson_json]{geojson_json()}};
see \code{\link[=geojson_json]{geojson_json()}} for help on what's supported here}
}
\value{
a list with TopoJSON
}
\description{
Convert many input types with spatial data to TopoJSON
as a list
}
\details{
Internally, we call \code{\link[=topojson_json]{topojson_json()}}, then use
an internal function to convert that JSON output to a list

The \code{type} parameter is automatically converted to
\code{type="auto"} if a sf, sfc, or sfg class is passed to \code{input}
}
\examples{
\dontrun{
# From a numeric vector of length 2 to a point
vec <- c(-99.74,32.45)
topojson_list(vec)

# Lists
## From a list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
topojson_list(mylist)

## From a list of numeric vectors to a polygon
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0), c(100.0,1.0), c(100.0,0.0))
topojson_list(vecs, geometry="polygon")

# from data.frame to points
(res <- topojson_list(us_cities[1:2,], lat='lat', lon='long'))
as.json(res)
## guess lat/long columns
topojson_list(us_cities[1:2,])
topojson_list(states[1:3,])
topojson_list(states[1:351,], geometry="polygon", group='group')
topojson_list(canada_cities[1:30,])

# from data.frame to polygons
head(states)
topojson_list(states[1:351, ], lat='lat', lon='long', geometry="polygon", group='group')

# From SpatialPolygons class
library('sp')
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
   c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
   c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
topojson_list(sp_poly)

# From SpatialPolygonsDataFrame class
sp_polydf <- as(sp_poly, "SpatialPolygonsDataFrame")
topojson_list(input = sp_polydf)

# From SpatialPoints class
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPoints(cbind(x,y))
topojson_list(s)

# From SpatialPointsDataFrame class
s <- SpatialPointsDataFrame(cbind(x,y), mtcars[1:5,])
topojson_list(s)

# From SpatialLines class
library("sp")
c1 <- cbind(c(1,2,3), c(3,2,2))
c2 <- cbind(c1[,1]+.05,c1[,2]+.05)
c3 <- cbind(c(1,2,3),c(1,1.5,1))
L1 <- Line(c1)
L2 <- Line(c2)
L3 <- Line(c3)
Ls1 <- Lines(list(L1), ID = "a")
Ls2 <- Lines(list(L2, L3), ID = "b")
sl1 <- SpatialLines(list(Ls1))
sl12 <- SpatialLines(list(Ls1, Ls2))
topojson_list(sl1)
topojson_list(sl12)
as.json(topojson_list(sl12))
as.json(topojson_list(sl12), pretty=TRUE)

# From SpatialLinesDataFrame class
dat <- data.frame(X = c("Blue", "Green"),
                 Y = c("Train", "Plane"),
                 Z = c("Road", "River"), row.names = c("a", "b"))
sldf <- SpatialLinesDataFrame(sl12, dat)
topojson_list(sldf)
as.json(topojson_list(sldf))
as.json(topojson_list(sldf), pretty=TRUE)

# From SpatialGrid
x <- GridTopology(c(0,0), c(1,1), c(5,5))
y <- SpatialGrid(x)
topojson_list(y)

# From SpatialGridDataFrame
sgdim <- c(3,4)
sg <- SpatialGrid(GridTopology(rep(0,2), rep(10,2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
topojson_list(sgdf)

# From SpatialRings
library("rgeos")
r1 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="1")
r2 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="2")
r1r2 <- SpatialRings(list(r1, r2))
topojson_list(r1r2)

# From SpatialRingsDataFrame
dat <- data.frame(id = c(1,2), value = 3:4)
r1r2df <- SpatialRingsDataFrame(r1r2, data = dat)
topojson_list(r1r2df)

# From SpatialPixels
library("sp")
pixels <- suppressWarnings(SpatialPixels(SpatialPoints(us_cities[c("long", "lat")])))
summary(pixels)
topojson_list(pixels)

# From SpatialPixelsDataFrame
library("sp")
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")], data = canada_cities)
)
topojson_list(pixelsdf)

# From SpatialCollections
library("sp")
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100), c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90), c(30,40,35,30)))), "2")
poly <- SpatialPolygons(list(poly1, poly2), 1:2)
coordinates(us_cities) <- ~long+lat
dat <- SpatialCollections(points = us_cities, polygons = poly)
out <- topojson_list(dat)
out[[1]]
out[[2]]
}

# From sf classes:
if (require(sf)) {
## sfg (a single simple features geometry)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  poly <- rbind(c(1,1), c(1,2), c(2,2), c(1,1))
  poly_sfg <- st_polygon(list(p1))
  topojson_list(poly_sfg)

## sfc (a collection of geometries)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  topojson_list(poly_sfc)

## sf (collection of geometries with attributes)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  poly_sf <- st_sf(foo = c("a", "b"), bar = 1:2, poly_sfc)
  topojson_list(poly_sf)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pretty.R
\name{pretty}
\alias{pretty}
\title{Convert json input to pretty printed output}
\usage{
pretty(x, indent = 4)
}
\arguments{
\item{x}{Input, character string}

\item{indent}{(integer) Number of spaces to indent}
}
\description{
Convert json input to pretty printed output
}
\details{
Only works with json class input. This is a simple wrapper around
\code{\link[jsonlite:prettify]{jsonlite::prettify()}}, so you can easily use that yourself.
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
% Please edit documentation in R/geo_topo.R
\name{geo2topo}
\alias{geo2topo}
\alias{topo2geo}
\title{GeoJSON to TopoJSON and back}
\usage{
geo2topo(x, object_name = "foo", quantization = 0, ...)

topo2geo(x, ...)
}
\arguments{
\item{x}{GeoJSON or TopoJSON as a character string, json, a file path, or
url}

\item{object_name}{(character) name to give to the TopoJSON object created.
Default: "foo"}

\item{quantization}{(numeric) quantization parameter, use this to
quantize geometry prior to computing topology. Typical values are powers of
ten (\code{1e4}, \code{1e5}, ...), default is \code{0} to not perform quantization.
For more information about quantization, see this by Mike Bostock
https://stackoverflow.com/questions/18900022/topojson-quantization-vs-simplification/18921214#18921214}

\item{...}{for \code{geo2topo} args passed  on to
\code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}}, and for \code{topo2geo} args passed  on to
\code{\link[sf:st_read]{sf::st_read()}}}
}
\value{
An object of class \code{json}, of either GeoJSON or TopoJSON
}
\description{
GeoJSON to TopoJSON and back
}
\examples{
# geojson to topojson
x <- '{"type": "LineString", "coordinates": [ [100.0, 0.0], [101.0, 1.0] ]}'
z <- geo2topo(x)
jsonlite::prettify(z)
\dontrun{
library(leaflet)
leaflet() \%>\%
  addProviderTiles(provider = "Stamen.Terrain") \%>\%
  addTopoJSON(z)
}

# geojson to topojson as a list
x <- list(
 '{"type": "LineString", "coordinates": [ [100, 0], [101, 1] ]}',
 '{"type": "LineString", "coordinates": [ [110, 0], [110, 1] ]}',
 '{"type": "LineString", "coordinates": [ [120, 0], [121, 1] ]}'
)
geo2topo(x)

# change the object name created
x <- '{"type": "LineString", "coordinates": [ [100.0, 0.0], [101.0, 1.0] ]}'
geo2topo(x, object_name = "HelloWorld")
geo2topo(x, object_name = "4")

x <- list(
 '{"type": "LineString", "coordinates": [ [100, 0], [101, 1] ]}',
 '{"type": "LineString", "coordinates": [ [110, 0], [110, 1] ]}',
 '{"type": "LineString", "coordinates": [ [120, 0], [121, 1] ]}'
)
geo2topo(x, "HelloWorld")
geo2topo(x, c("A", "B", "C"))


# topojson to geojson
w <- topo2geo(z)
jsonlite::prettify(w)

## larger examples
file <- system.file("examples", "us_states.topojson", package = "geojsonio")
topo2geo(file)
}
\seealso{
\code{\link[=topojson_write]{topojson_write()}}, \code{\link[=topojson_read]{topojson_read()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.location.R
\name{as.location}
\alias{as.location}
\title{Convert a path or URL to a location object.}
\usage{
as.location(x, ...)
}
\arguments{
\item{x}{Input.}

\item{...}{Ignored.}
}
\description{
Convert a path or URL to a location object.
}
\examples{
\dontrun{
# A file
file <- system.file("examples", "zillow_or.geojson", package = "geojsonio")
as.location(file)

# A URL
url <- "https://raw.githubusercontent.com/glynnbird/usstatesgeojson/master/california.geojson"
as.location(url)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.json.R
\name{as.json}
\alias{as.json}
\title{Convert inputs to JSON}
\usage{
as.json(x, ...)
}
\arguments{
\item{x}{Input}

\item{...}{Further args passed on to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}
}
\description{
Convert inputs to JSON
}
\details{
when the output of \code{\link[=topojson_list]{topojson_list()}} is given to
this function we use a special internal fxn \code{astjl()} to
parse the object - see that fxn and let us know if any
problems you run in to
}
\examples{
\dontrun{
(res <- geojson_list(us_cities[1:2,], lat='lat', lon='long'))
as.json(res)
as.json(res, pretty = TRUE)

vec <- c(-99.74,32.45)
as.json(geojson_list(vec))
as.json(geojson_list(vec), pretty = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojsonio-package.r
\name{geojsonio-defunct}
\alias{geojsonio-defunct}
\title{Defunct functions in geojsonio}
\description{
\itemize{
\item \code{\link[=lint]{lint()}}: See \code{geojsonlint::geojson_hint}
\item \code{\link[=validate]{validate()}}: See \code{geojsonlint::geojson_lint}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/construction.R
\name{geojson-add}
\alias{geojson-add}
\alias{+.geo_list}
\alias{+.json}
\title{Add together geo_list or json objects}
\usage{
\method{+}{geo_list}(x1, x2)

\method{+}{json}(x1, x2)
}
\arguments{
\item{x1}{An object of class \code{geo_list} or \code{json}}

\item{x2}{A component to add to \code{x1}, of class \code{geo_list} or \code{json}}
}
\description{
Add together geo_list or json objects
}
\details{
If the first object is an object of class \code{geo_list}, you
can add another object of class \code{geo_list} or of class \code{json},
and will result  in a \code{geo_list} object.

If the first object is an object of class \code{json}, you can add
another object of class \code{json} or of class \code{geo_list}, and will result
in a \code{json} object.
}
\examples{
\dontrun{
# geo_list + geo_list
## Note: geo_list is the output type from geojson_list, it's just a list with
## a class attached so we know it's geojson :)
vec <- c(-99.74,32.45)
a <- geojson_list(vec)
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0),
  c(100.0,1.0), c(100.0,0.0))
b <- geojson_list(vecs, geometry="polygon")
a + b

# json + json
c <- geojson_json(c(-99.74,32.45))
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0),
  c(100.0,1.0), c(100.0,0.0))
d <- geojson_json(vecs, geometry="polygon")
c + d
(c + d) \%>\% pretty
}
}
\seealso{
\code{\link[=geojson_list]{geojson_list()}}, \code{\link[=geojson_json]{geojson_json()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapgist.R
\name{map_gist}
\alias{map_gist}
\title{Publish an interactive map as a GitHub gist}
\usage{
map_gist(
  input,
  lat = "lat",
  lon = "long",
  geometry = "point",
  group = NULL,
  type = "FeatureCollection",
  file = "myfile.geojson",
  description = "",
  public = TRUE,
  browse = TRUE,
  ...
)
}
\arguments{
\item{input}{Input object}

\item{lat}{Name of latitude variable}

\item{lon}{Name of longitude variable}

\item{geometry}{(character) Are polygons in the object}

\item{group}{(character) A grouping variable to perform grouping for
polygons - doesn't apply for points}

\item{type}{(character) One of FeatureCollection or GeometryCollection}

\item{file}{File name to use to put up as the gist file}

\item{description}{Description for the GitHub gist, or leave to default
(=no description)}

\item{public}{(logical) Want gist to be public or not? Default: \code{TRUE}}

\item{browse}{If \code{TRUE} (default) the map opens in your default browser.}

\item{...}{Further arguments passed on to \code{httr::POST}}
}
\description{
There are two ways to authorize to work with your GitHub
account:
\itemize{
\item PAT - Generate a personal access token (PAT) at
https://help.github.com/articles/creating-an-access-token-for-command-line-use
and record it in the \code{GITHUB_PAT} envar in your \code{.Renviron} file.
\item Interactive - Interactively login into your GitHub account and authorise
with OAuth.
}

Using the PAT method is recommended.

Using the \code{gist_auth()} function you can authenticate separately first, or
if you're not authenticated, this function will run internally with each
function call. If you have a PAT, that will be used, if not, OAuth will
be used.
}
\examples{
\dontrun{
if (!identical(Sys.getenv("GITHUB_PAT"), "")) {

# From file
file <- "myfile.geojson"
geojson_write(us_cities[1:20, ], lat='lat', lon='long', file = file)
map_gist(file=as.location(file))

# From SpatialPoints class
library("sp")
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPoints(cbind(x,y))
map_gist(s)

# from SpatialPointsDataFrame class
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPointsDataFrame(cbind(x,y), mtcars[1:5,])
map_gist(s)

# from SpatialPolygons class
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
   c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
   c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
map_gist(sp_poly)

# From SpatialPolygonsDataFrame class
sp_polydf <- as(sp_poly, "SpatialPolygonsDataFrame")
map_gist(sp_poly)

# From SpatialLines class
c1 <- cbind(c(1,2,3), c(3,2,2))
c2 <- cbind(c1[,1]+.05,c1[,2]+.05)
c3 <- cbind(c(1,2,3),c(1,1.5,1))
L1 <- Line(c1)
L2 <- Line(c2)
L3 <- Line(c3)
Ls1 <- Lines(list(L1), ID = "a")
Ls2 <- Lines(list(L2, L3), ID = "b")
sl1 <- SpatialLines(list(Ls1))
sl12 <- SpatialLines(list(Ls1, Ls2))
map_gist(sl1)

# From SpatialLinesDataFrame class
dat <- data.frame(X = c("Blue", "Green"),
                 Y = c("Train", "Plane"),
                 Z = c("Road", "River"), row.names = c("a", "b"))
sldf <- SpatialLinesDataFrame(sl12, dat)
map_gist(sldf)

# From SpatialGrid
x <- GridTopology(c(0,0), c(1,1), c(5,5))
y <- SpatialGrid(x)
map_gist(y)

# From SpatialGridDataFrame
sgdim <- c(3,4)
sg <- SpatialGrid(GridTopology(rep(0,2), rep(10,2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
map_gist(sgdf)

# from data.frame
## to points
map_gist(us_cities)

## to polygons
head(states)
map_gist(states[1:351, ], lat='lat', lon='long', geometry="polygon", group='group')

## From a list
mylist <- list(list(lat=30, long=120, marker="red"),
               list(lat=30, long=130, marker="blue"))
map_gist(mylist, lat="lat", lon="long")

# From a numeric vector
## of length 2 to a point
vec <- c(-99.74,32.45)
map_gist(vec)

## this requires numeric class input, so inputting a list will dispatch on the list method
poly <- c(c(-114.345703125,39.436192999314095),
          c(-114.345703125,43.45291889355468),
          c(-106.61132812499999,43.45291889355468),
          c(-106.61132812499999,39.436192999314095),
          c(-114.345703125,39.436192999314095))
map_gist(poly, geometry = "polygon")

# From a json object
(x <- geojson_json(c(-99.74,32.45)))
map_gist(x)
## another example
map_gist(geojson_json(us_cities[1:10,], lat='lat', lon='long'))

# From a geo_list object
(res <- geojson_list(us_cities[1:2,], lat='lat', lon='long'))
map_gist(res)

# From SpatialPixels
pixels <- suppressWarnings(SpatialPixels(SpatialPoints(us_cities[c("long", "lat")])))
summary(pixels)
map_gist(pixels)

# From SpatialPixelsDataFrame
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")], data = canada_cities)
)
map_gist(pixelsdf)

# From SpatialRings
library("rgeos")
r1 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="1")
r2 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="2")
r1r2 <- SpatialRings(list(r1, r2))
map_gist(r1r2)

# From SpatialRingsDataFrame
dat <- data.frame(id = c(1,2), value = 3:4)
r1r2df <- SpatialRingsDataFrame(r1r2, data = dat)
map_gist(r1r2df)

}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.R
\name{validate-defunct}
\alias{validate-defunct}
\alias{validate}
\title{Validate a geoJSON file, json object, list, or Spatial class.}
\usage{
validate(...)
}
\arguments{
\item{...}{ignored}
}
\description{
Validate a geoJSON file, json object, list, or Spatial class.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topojson_json.R
\name{topojson_json}
\alias{topojson_json}
\title{Convert many input types with spatial data to TopoJSON
as a JSON string}
\usage{
topojson_json(
  input,
  lat = NULL,
  lon = NULL,
  group = NULL,
  geometry = "point",
  type = "FeatureCollection",
  convert_wgs84 = FALSE,
  crs = NULL,
  object_name = "foo",
  quantization = 0,
  ...
)
}
\arguments{
\item{input}{Input list, data.frame, spatial class, or sf class. Inputs can
also be dplyr \code{tbl_df} class since it inherits from \code{data.frame}.}

\item{lat}{(character) Latitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{lon}{(character) Longitude name. The default is \code{NULL}, and we
attempt to guess.}

\item{group}{(character) A grouping variable to perform grouping for
polygons - doesn't apply for points}

\item{geometry}{(character) One of point (Default) or polygon.}

\item{type}{(character) The type of collection. One of 'auto' (default
for 'sf' objects), 'FeatureCollection' (default for everything else), or
'GeometryCollection'. "skip" skips the coercion with package \pkg{geojson}
functions; skipping can save significant run time on larger geojson
objects. \code{Spatial} objects can only accept "FeatureCollection" or "skip".
"skip" is not available as an option for \code{numeric}, \code{list},
and \code{data.frame} classes}

\item{convert_wgs84}{Should the input be converted to the
standard CRS system for GeoJSON (https://tools.ietf.org/html/rfc7946)
(geographic coordinate reference system, using
the WGS84 datum, with longitude and latitude units of decimal degrees;
EPSG: 4326). Default is \code{FALSE} though this may change in a future
package version. This will only work for \code{sf} or \code{Spatial}
objects with a CRS already defined. If one is not defined but you know
what it is, you may define it in the \code{crs} argument below.}

\item{crs}{The CRS of the input if it is not already defined. This can be
an epsg code as a four or five digit integer or a valid proj4 string.
This argument will be ignored if \code{convert_wgs84} is \code{FALSE} or
the object already has a CRS.}

\item{object_name}{(character) name to give to the TopoJSON object created.
Default: "foo"}

\item{quantization}{(numeric) quantization parameter, use this to
quantize geometry prior to computing topology. Typical values are powers of
ten (\code{1e4}, \code{1e5}, ...), default is \code{0} to not perform quantization.
For more information about quantization, see this by Mike Bostock
https://stackoverflow.com/questions/18900022/topojson-quantization-vs-simplification/18921214#18921214}

\item{...}{args passed down to \code{\link[=geojson_json]{geojson_json()}}; see \code{\link[=geojson_json]{geojson_json()}} for
help on what's supported here}
}
\value{
An object of class \code{geo_json} (and \code{json})
}
\description{
Convert many input types with spatial data to TopoJSON
as a JSON string
}
\details{
The \code{type} parameter is automatically converted to
\code{type="auto"} if a sf, sfc, or sfg class is passed to \code{input}
}
\examples{
\dontrun{
# From a numeric vector of length 2, making a point type
topojson_json(c(-99.74,32.45), pretty=TRUE)
topojson_json(c(-99.74,32.45), type = "GeometryCollection")

## polygon type
### this requires numeric class input, so inputting a list will dispatch on the list method
poly <- c(c(-114.345703125,39.436192999314095),
          c(-114.345703125,43.45291889355468),
          c(-106.61132812499999,43.45291889355468),
          c(-106.61132812499999,39.436192999314095),
          c(-114.345703125,39.436192999314095))
topojson_json(poly, geometry = "polygon", pretty=TRUE)

# Lists
## From a list of numeric vectors to a polygon
vecs <- list(c(100.0,0.0), c(101.0,0.0), c(101.0,1.0), c(100.0,1.0), c(100.0,0.0))
topojson_json(vecs, geometry="polygon", pretty=TRUE)

## from a named list
mylist <- list(list(latitude=30, longitude=120, marker="red"),
               list(latitude=30, longitude=130, marker="blue"))
topojson_json(mylist, lat='latitude', lon='longitude')

# From a data.frame to points
topojson_json(us_cities[1:2,], lat='lat', lon='long', pretty=TRUE)
topojson_json(us_cities[1:2,], lat='lat', lon='long',
   type="GeometryCollection", pretty=TRUE)

# from data.frame to polygons
head(states)
## make list for input to e.g., rMaps
topojson_json(states[1:351, ], lat='lat', lon='long', geometry="polygon", group='group')

# from a geo_list
a <- geojson_list(us_cities[1:2,], lat='lat', lon='long')
topojson_json(a)

# sp classes

## From SpatialPolygons class
library('sp')
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100),
   c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90),
   c(30,40,35,30)))), "2")
sp_poly <- SpatialPolygons(list(poly1, poly2), 1:2)
topojson_json(sp_poly)
topojson_json(sp_poly, pretty=TRUE)

## Another SpatialPolygons
library("sp")
library("rgeos")
pt <- SpatialPoints(coordinates(list(x = 0, y = 0)), CRS("+proj=longlat +datum=WGS84"))
## transfrom to web mercator becuase geos needs project coords
crs <- gsub("\n", "", paste0("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0
   +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs", collapse = ""))
pt <- spTransform(pt, CRS(crs))
## buffer
pt <- gBuffer(pt, width = 100)
pt <- spTransform(pt, CRS("+proj=longlat +datum=WGS84"))
topojson_json(pt)

## data.frame to geojson
geojson_write(us_cities[1:2,], lat='lat', lon='long') \%>\% as.json

# From SpatialPoints class
x <- c(1,2,3,4,5)
y <- c(3,2,5,1,4)
s <- SpatialPoints(cbind(x,y))
topojson_json(s)

## From SpatialPointsDataFrame class
s <- SpatialPointsDataFrame(cbind(x,y), mtcars[1:5,])
topojson_json(s)

## From SpatialLines class
library("sp")
c1 <- cbind(c(1,2,3), c(3,2,2))
c2 <- cbind(c1[,1]+.05,c1[,2]+.05)
c3 <- cbind(c(1,2,3),c(1,1.5,1))
L1 <- Line(c1)
L2 <- Line(c2)
L3 <- Line(c3)
Ls1 <- Lines(list(L1), ID = "a")
Ls2 <- Lines(list(L2, L3), ID = "b")
sl1 <- SpatialLines(list(Ls1))
sl12 <- SpatialLines(list(Ls1, Ls2))
topojson_json(sl1)
topojson_json(sl12)

## From SpatialLinesDataFrame class
dat <- data.frame(X = c("Blue", "Green"),
                 Y = c("Train", "Plane"),
                 Z = c("Road", "River"), row.names = c("a", "b"))
sldf <- SpatialLinesDataFrame(sl12, dat)
topojson_json(sldf)
topojson_json(sldf, pretty=TRUE)

## From SpatialGrid
x <- GridTopology(c(0,0), c(1,1), c(5,5))
y <- SpatialGrid(x)
topojson_json(y)

## From SpatialGridDataFrame
sgdim <- c(3,4)
sg <- SpatialGrid(GridTopology(rep(0,2), rep(10,2), sgdim))
sgdf <- SpatialGridDataFrame(sg, data.frame(val = 1:12))
topojson_json(sgdf)

# From SpatialRings
library("rgeos")
r1 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="1")
r2 <- Ring(cbind(x=c(1,1,2,2,1), y=c(1,2,2,1,1)), ID="2")
r1r2 <- SpatialRings(list(r1, r2))
topojson_json(r1r2)

# From SpatialRingsDataFrame
dat <- data.frame(id = c(1,2), value = 3:4)
r1r2df <- SpatialRingsDataFrame(r1r2, data = dat)
topojson_json(r1r2df)

# From SpatialPixels
library("sp")
pixels <- suppressWarnings(SpatialPixels(SpatialPoints(us_cities[c("long", "lat")])))
summary(pixels)
topojson_json(pixels)

# From SpatialPixelsDataFrame
library("sp")
pixelsdf <- suppressWarnings(
 SpatialPixelsDataFrame(points = canada_cities[c("long", "lat")], data = canada_cities)
)
topojson_json(pixelsdf)

# From SpatialCollections
library("sp")
library("rgeos")
pts <- SpatialPoints(cbind(c(1,2,3,4,5), c(3,2,5,1,4)))
poly1 <- Polygons(list(Polygon(cbind(c(-100,-90,-85,-100), c(40,50,45,40)))), "1")
poly2 <- Polygons(list(Polygon(cbind(c(-90,-80,-75,-90), c(30,40,35,30)))), "2")
poly <- SpatialPolygons(list(poly1, poly2), 1:2)
dat <- SpatialCollections(pts, polygons = poly)
topojson_json(dat)

# From sf classes:
if (require(sf)) {
## sfg (a single simple features geometry)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  poly <- rbind(c(1,1), c(1,2), c(2,2), c(1,1))
  poly_sfg <- st_polygon(list(p1))
  topojson_json(poly_sfg)

## sfc (a collection of geometries)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  topojson_json(poly_sfc)

## sf (collection of geometries with attributes)
  p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
  p2 <- rbind(c(5,5), c(5,6), c(4,5), c(5,5))
  poly_sfc <- st_sfc(st_polygon(list(p1)), st_polygon(list(p2)))
  poly_sf <- st_sf(foo = c("a", "b"), bar = 1:2, poly_sfc)
  topojson_json(poly_sf)
}

## Pretty print a json string
topojson_json(c(-99.74,32.45))
topojson_json(c(-99.74,32.45)) \%>\% pretty
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojsonio-package.r
\docType{data}
\name{canada_cities}
\alias{canada_cities}
\title{This is the same data set from the maps library, named differently}
\format{
A list with 6 components, namely "name", "country.etc", "pop",
"lat", "long", and "capital", containing the city name, the province
abbreviation, approximate population (as at January 2006), latitude,
longitude and capital status indication (0 for non-capital, 1 for capital,
2 for provincial
}
\description{
This database is of Canadian cities of population greater than about 1,000.
Also included are province capitals of any population size.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojson_read.R
\name{geojson_read}
\alias{geojson_read}
\title{Read geojson or other formats from a local file or a URL}
\usage{
geojson_read(
  x,
  parse = FALSE,
  what = "list",
  stringsAsFactors = FALSE,
  query = NULL,
  ...
)
}
\arguments{
\item{x}{(character) Path to a local file or a URL.}

\item{parse}{(logical) To parse geojson to data.frame like structures if
possible. Default: \code{FALSE}}

\item{what}{(character) What to return. One of "list", "sp" (for
Spatial class), or "json". Default: "list". "list" "and" sp run through
package \pkg{sf}. if "json", returns json as character class}

\item{stringsAsFactors}{Convert strings to Factors? Default \code{FALSE}.}

\item{query}{(character) A SQL query, see also \link{postgis}}

\item{...}{Further args passed on to \code{\link[sf:st_read]{sf::st_read()}}}
}
\value{
various, depending on what's chosen in \code{what} parameter
\itemize{
\item list: geojson as a list using \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}}
\item sp: geojson as an sp class object using \code{\link[sf:st_read]{sf::st_read()}}
\item json: geojson as character string, to parse downstream as you wish
}
}
\description{
Read geojson or other formats from a local file or a URL
}
\details{
This function supports various geospatial file formats from a URL,
as well as local kml, shp, and geojson file formats.
}
\section{Linting GeoJSON}{

If you're having trouble rendering GeoJSON files, ensure you have a valid
GeoJSON file by running it through the package \code{geojsonlint}, which
has a variety of different GeoJSON linters.
}

\section{File size}{

We previously used \code{\link[=file_to_geojson]{file_to_geojson()}} in this function, leading to
file size problems; this should no longer be a concern, but let us know
if you run into file size problems
}

\examples{
\dontrun{
# From a file
file <- system.file("examples", "california.geojson", package = "geojsonio")
(out <- geojson_read(file))
geojson_read(file)

# From a URL
url <- "https://raw.githubusercontent.com/glynnbird/usstatesgeojson/master/california.geojson"
geojson_read(url)
geojson_read(url, parse = TRUE)

# Use as.location first if you want
geojson_read(as.location(file))

# output a SpatialClass object
## read kml
file <- system.file("examples", "norway_maple.kml", package = "geojsonio")
geojson_read(as.location(file), what = "sp")
## read geojson
file <- system.file("examples", "california.geojson", package = "geojsonio")
geojson_read(as.location(file), what = "sp")
## read geojson from a url
url <- "https://raw.githubusercontent.com/glynnbird/usstatesgeojson/master/california.geojson"
geojson_read(url, what = "sp")
## read from a shape file
file <- system.file("examples", "bison.zip", package = "geojsonio")
dir <- tempdir()
unzip(file, exdir = dir)
shpfile <- list.files(dir, pattern = ".shp", full.names = TRUE)
geojson_read(shpfile, what = "sp")

x <- "https://raw.githubusercontent.com/johan/world.geo.json/master/countries.geo.json"
geojson_read(x, what = "sp")
geojson_read(x, what = "list")

utils::download.file(x, destfile = basename(x))
geojson_read(basename(x), what = "sp")

# from a Postgres database - your Postgres instance must be running
## MAKE SURE to run the setup in the postgis manual file first!
if (requireNamespace("DBI") && requireNamespace("RPostgres")) {
library(DBI)
conn <- tryCatch(dbConnect(RPostgres::Postgres(), dbname = 'postgistest'), 
 error = function(e) e)
if (inherits(conn, "PqConnection")) {
  state <- "SELECT row_to_json(fc)
   FROM (SELECT 'FeatureCollection' As type, array_to_json(array_agg(f)) As features
   FROM (SELECT 'Feature' As type
     , ST_AsGeoJSON(lg.geog)::json As geometry
     , row_to_json((SELECT l FROM (SELECT loc_id, loc_name) As l
       )) As properties
    FROM locations As lg   ) As f )  As fc;"
  json <- geojson_read(conn, query = state, what = "json")
  map_leaf(json)
 }
}
}
}
\seealso{
\code{\link[=topojson_read]{topojson_read()}}, \code{\link[=geojson_write]{geojson_write()}} \link{postgis}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geojsonio-package.r
\docType{package}
\name{geojsonio}
\alias{geojsonio}
\title{\strong{I/O for GeoJSON}}
\description{
Convert various data formats to/from GeoJSON or TopoJSON. This
package focuses mostly on converting lists, data.frame's, numeric,
SpatialPolygons, SpatialPolygonsDataFrame, and more to GeoJSON with the
help of \pkg{sf}. You can currently read TopoJSON - writing
TopoJSON will come in a future version of this package.
}
\section{Package organization}{

The core functions in this package are organized first around what you're
working with or want to get, GeoJSON or TopoJSON, then convert to or read
from various formats:
\itemize{
\item \code{\link[=geojson_list]{geojson_list()}} / \code{\link[=topojson_list]{topojson_list()}} - convert
to GeoJSON or TopoJSON as R list format
\item \code{\link[=geojson_json]{geojson_json()}} / \code{\link[=topojson_json]{topojson_json()}} - convert
to GeoJSON or TopoJSON as JSON
\item \code{\link[=geojson_sp]{geojson_sp()}} - convert to a spatial object from
\code{geojson_list} or \code{geojson_json}
\item \code{\link[=geojson_sf]{geojson_sf()}} - convert to an sf object from
\code{geojson_list} or \code{geojson_json}
\item \code{\link[=geojson_read]{geojson_read()}} / \code{\link[=topojson_read]{topojson_read()}} - read a
GeoJSON/TopoJSON file from file path or URL
\item \code{\link[=geojson_write]{geojson_write()}} / \code{\link[=topojson_write]{topojson_write()}} - write
a GeoJSON file locally (TopoJSON coming later)
}

Other interesting functions:
\itemize{
\item \code{\link[=map_gist]{map_gist()}} - Create a GitHub gist (renders as an
interactive map)
\item \code{\link[=map_leaf]{map_leaf()}} - Create a local interactive map using the
\code{leaflet} package
\item \code{\link[=geo2topo]{geo2topo()}} - Convert GeoJSON to TopoJSON
\item \code{\link[=topo2geo]{topo2geo()}} - Convert TopoJSON to GeoJSON
}

All of the above functions have methods for various classes, including
\code{numeric} vectors, \code{data.frame}, \code{list}, \code{SpatialPolygons}, \code{SpatialLines},
\code{SpatialPoints}, and many more - which will try to do the right thing
based on the data you give as input.
}

\author{
Scott Chamberlain

Andy Teucher \email{andy.teucher@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projections.r
\name{projections}
\alias{projections}
\title{topojson projections and extensions}
\usage{
projections(
  proj,
  rotate = NULL,
  center = NULL,
  translate = NULL,
  scale = NULL,
  clipAngle = NULL,
  precision = NULL,
  parallels = NULL,
  clipExtent = NULL,
  invert = NULL
)
}
\arguments{
\item{proj}{Map projection name. One of albers, albersUsa, azimuthalEqualArea,
azimuthalEquidistant, conicEqualArea, conicConformal, conicEquidistant, equirectangular,
gnomonic, mercator, orthographic, stereographic, or transverseMercator.}

\item{rotate}{If rotation is specified, sets the projection's three-axis rotation to the
specified angles yaw, pitch and roll (or equivalently longitude, latitude and roll)
in degrees and returns the projection. If rotation is not specified, returns the current
rotation which defaults \verb{[0, 0, 0]}. If the specified rotation has only two values, rather than
three, the roll is assumed to be 0.}

\item{center}{If center is specified, sets the projection's center to the specified location, a
two-element array of longitude and latitude in degrees and returns the projection. If center is
not specified, returns the current center which defaults to (0,0)}

\item{translate}{If point is specified, sets the projection's translation offset to the
specified two-element array \verb{[x, y]} and returns the projection. If point is not specified,
returns the current translation offset which defaults to \verb{[480, 250]}. The translation offset
determines the pixel coordinates of the projection's center. The default translation offset
places (0,0) at the center of a 960x500 area.}

\item{scale}{If scale is specified, sets the projection's scale factor to the specified value
and returns the projection. If scale is not specified, returns the current scale factor which
defaults to 150. The scale factor corresponds linearly to the distance between projected points.
However, scale factors are not consistent across projections.}

\item{clipAngle}{If angle is specified, sets the projection's clipping circle radius to the
specified angle in degrees and returns the projection. If angle is null, switches to
antimeridian cutting rather than small-circle clipping. If angle is not specified, returns the
current clip angle which defaults to null. Small-circle clipping is independent of viewport
clipping via clipExtent.}

\item{precision}{If precision is specified, sets the threshold for the projection's adaptive
resampling to the specified value in pixels and returns the projection. This value corresponds
to the Douglas-Peucker distance. If precision is not specified, returns the projection's current
resampling precision which defaults to Math.SQRT(1/2).}

\item{parallels}{Depends on the projection used! See
https://github.com/mbostock/d3/wiki/Geo-Projections#standard-projections for help}

\item{clipExtent}{If extent is specified, sets the projection's viewport clip extent to the
specified bounds in pixels and returns the projection. The extent bounds are specified as an
array \verb{[[x0, y0], [x1, y1]]}, where x0 is the left-side of the viewport, y0 is the top, x1 is
the right and y1 is the bottom. If extent is null, no viewport clipping is performed. If extent
is not specified, returns the current viewport clip extent which defaults to null. Viewport
clipping is independent of small-circle clipping via clipAngle.}

\item{invert}{Projects backward from Cartesian coordinates (in pixels) to spherical coordinates
(in degrees). Returns an array \verb{[longitude, latitude]} given the input array \verb{[x, y]}.}
}
\description{
topojson projections and extensions
}
\examples{
projections(proj="albers")
projections(proj="albers", rotate='[98 + 00 / 60, -35 - 00 / 60]', scale=5700)
projections(proj="albers", scale=5700)
projections(proj="albers", translate='[55 * width / 100, 52 * height / 100]')
projections(proj="albers", clipAngle=90)
projections(proj="albers", precision=0.1)
projections(proj="albers", parallels='[30, 62]')
projections(proj="albers", clipExtent='[[105 - 87, 40], [105 + 87 + 1e-6, 82 + 1e-6]]')
projections(proj="albers", invert=60)
projections("orthographic")
}
