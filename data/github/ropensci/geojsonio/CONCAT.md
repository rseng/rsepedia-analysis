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
*Wow, no problems at all. :)*