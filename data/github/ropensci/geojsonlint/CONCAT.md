geojsonlint
===========



[![cran checks](https://cranchecks.info/badges/worst/geojsonlint)](https://cranchecks.info/pkgs/geojsonlint)
[![R-check](https://github.com/ropensci/geojsonlint/workflows/R-check/badge.svg)](https://github.com/ropensci/geojsonlint/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/geojsonlint/coverage.svg?branch=master)](https://codecov.io/github/ropensci/geojsonlint?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/geojsonlint)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/geojsonlint)](https://cran.r-project.org/package=geojsonlint)

GeoJSON linters available in `geojsonlint`

* [GeoJSON hint JS library](https://www.npmjs.com/package/geojsonhint) - via `geojson_hint()` - currently using `geojsonhint` `v1.2.1`
* [is-my-json-valid JS library](https://www.npmjs.com/package/is-my-json-valid) - via `geojson_validate()`

Both functions return the same outputs. If the GeoJSON is valid, they return `TRUE`.
If the GeoJSON is invalid, they return `FALSE`, plus reason(s) that the GeoJSON is invalid
in an attribute named _errors_ as a data.frame. The fields in the data.frame's are not
the same across functions unfortunately, but they can be easily coerced to combine via
e.g., `plyr::rbind.fill` or `dplyr::bind_rows` or `data.table::rbindlist(fill = TRUE)`

The parameters for the three functions are similar, though `geojson_validate()` has an
extra parameter `greedy` that's not available in the others.

## Installation

from CRAN


```r
install.packages("geojsonlint")
```

Dev version


```r
remotes::install_github("ropensci/geojsonlint")
```


```r
library("geojsonlint")
```

## Good GeoJSON

geojsonhint JS library


```r
geojson_hint(x = '{"type": "Point", "coordinates": [-100, 80]}')
#> [1] TRUE
```

is-my-json-valid JS library


```r
geojson_validate(x = '{"type": "Point", "coordinates": [-100, 80]}')
#> [1] TRUE
```

## Bad GeoJSON

geojsonhint JS library


```r
geojson_hint('{"type": "FooBar"}')
#> [1] FALSE
```

is-my-json-valid JS library


```r
geojson_validate('{ "type": "FeatureCollection" }')
#> [1] FALSE
```

## Bad GeoJSON - with reason for failure

geojsonhint JS library


```r
geojson_hint('{"type": "FooBar"}', inform = TRUE)
#> [1] FALSE
#> attr(,"errors")
#>   line                    message
#> 1    1 The type FooBar is unknown
```

is-my-json-valid JS library


```r
geojson_validate('{ "type": "FeatureCollection" }', inform = TRUE)
#> [1] FALSE
#> attr(,"errors")
#>   field                             message
#> 1  data no (or more than one) schemas match
```

## Bad GeoJSON - stop on validation failure

geojsonhint JS library


```r
geojson_hint('{"type": "FooBar"}', error = TRUE)
#> Error: Line 1
#>    - The type FooBar is unknown
```

is-my-json-valid JS library


```r
geojson_validate('{ "type": "FeatureCollection" }', error = TRUE)
#> Error: 1 error validating json:
#> 	- data: no (or more than one) schemas match
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/geojsonlint/issues).
* License: MIT
* Get citation information for `geojsonlint` in R doing `citation(package = 'geojsonlint')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.


[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
geojsonlint 0.4.0
=================

### DEFUNCT

* `geojson_lint()` is now defunct; the geojsonlint.com API appears to be gone and there's no sign of it coming back (#20)

## MINOR IMPROVEMENTS

* cran checks identified a failing example - resulting from occassional github bad response - only run these examples if interactive (#19)


geojsonlint 0.3.0
=================

## MINOR IMPROVEMENTS

* update JS library `geojsonhint` to `v2.1.0` (#11)
* replace `httr` with `crul` (#12)
* the parameter `verbose` replaced with `inform` throughout the package. take note if you have `verbose` parameter in use in any R code


geojsonlint 0.2.0
=================

## NEW FEATURES

* Now using a new version of the JS library `geojsonhint` (`v2.0.0-beta2`)
(#9) - see `geojsonhint` Changelog for changes to the JS library:
<https://github.com/mapbox/geojsonhint/blob/master/CHANGELOG.md>


geojsonlint 0.1.0
=================

## NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 3.6.2
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 3 downstream dependencies. (Summary at <https://github.com/ropensci/geojsonlint/blob/master/revdep/README.md>). No problems were found.

---

This version fixes a failing check on CRAN checks for an example; in addition, a function is made defunct as the web service it used has gone away.

I expect there to be failures in revdepchecks for sen2r. Like my recent geojsonio CRAN submission, geojsonio no longer relies on rgdal causing sen2r to fail because it doesn't correctly use rgdal conditionally. I expect the same thing to happen here because this package does not use rgdal either. If there's sen2r failures due to rgdal, that's unrelated to this package.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/geojsonlint/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/geojsonlint.git`
* Make sure to track progress upstream (i.e., on our version of `geojsonlint` at `ropensci/geojsonlint`) by doing `git remote add upstream https://github.com/ropensci/geojsonlint.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/geojsonlint`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
The geojson schema here is taken from the [SchemaStore](https://github.com/SchemaStore/schemastore) repository. It is distributed under the Apache 2.0 license (see LICENSE.md)
# JS dependencies

## geojsonhint

Currently (as of 2019-02-08) using geojsonhint v2.1.0

To recreate `inst/js/geojsohint.js`:

Install `geojsonhint` NPM library

```
npm install -g @mapbox/geojsonhint
```

Browserify

```
echo "global.geojsonhint = require('geojsonhint');" > in.js
browserify in.js -o geojsonhint.js
```

Copy js file into the `inst/js` directory in the `geojsonlint` package

```
cp geojsonhint.js geojsonlint/inst/js
```
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.6.2 Patched (2020-02-04 r77774) |
|os       |macOS Mojave 10.14.6                        |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-02-12                                  |

# Dependencies

|package     |old   |new   |Δ  |
|:-----------|:-----|:-----|:--|
|geojsonlint |0.3.0 |0.4.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*