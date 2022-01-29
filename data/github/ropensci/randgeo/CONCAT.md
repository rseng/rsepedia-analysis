randgeo: random WKT and GeoJSON
===============================



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/randgeo)](https://cranchecks.info/pkgs/randgeo)
[![R-check](https://github.com/ropensci/randgeo/workflows/R-check/badge.svg)](https://github.com/ropensci/randgeo/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/randgeo/coverage.svg?branch=master)](https://codecov.io/github/ropensci/randgeo?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/randgeo?color=C9A115)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/randgeo)](https://cran.r-project.org/package=randgeo)

**randgeo** generates random points and shapes in GeoJSON and WKT formats for use
in examples, teaching, or statistical applications.

Points and shapes are generated in the long/lat coordinate system and with
appropriate spherical geometry; random points are distributed evenly across
the globe, and random shapes are sized according to a maximum great-circle
distance from the center of the shape. 

**randgeo** was adapted from <https://github.com/tmcw/geojson-random> to have a pure R
implementation __without any dependencies__ as well as appropriate geometry. Data generated
by **randgeo** may be processed or displayed of with packages such as
[**sf**](https://cran.r-project.org/package=sf),
[**wicket**](https://cran.r-project.org/package=wicket),
[**geojson**](https://cran.r-project.org/package=geojson),
[**wellknown**](https://cran.r-project.org/package=wellknown),
[**geojsonio**](https://cran.r-project.org/package=geojsonio), or
[**lawn**](https://cran.r-project.org/package=lawn).

Package API:

* `rg_position` - random position (lon, lat)
* `geo_point` - random GeoJSON point
* `geo_linestring` - random GeoJSON linestring
* `geo_polygon` - random GeoJSON polygon
* `wkt_point` - random WKT point
* `wkt_linestring` - random WKT linestring
* `wkt_polygon` - random WKT polygon

## Docs

<https://docs.ropensci.org/randgeo/>

## Install

Stabler CRAN version


```r
install.packages("randgeo")
```

Development version


```r
devtools::install_github("ropensci/randgeo")
```


```r
library("randgeo")
```



## Generate a random position


```r
rg_position()
#> [[1]]
#> [1] 149.33018 -60.94463
```

## Generate random GeoJSON

Random point - evenly distributed across the sphere.  The `bbox` option allows
you to limit points to within long/lat bounds.


```r
geo_point()
#> $type
#> [1] "FeatureCollection"
#> 
#> $features
#> $features[[1]]
#> $features[[1]]$type
#> [1] "Feature"
#> 
#> $features[[1]]$geometry
#> $features[[1]]$geometry$type
#> [1] "Point"
#> 
#> $features[[1]]$geometry$coordinates
#> [1] -76.98977 -41.36819
#> 
#> 
#> $features[[1]]$properties
#> NULL
#> 
#> 
#> 
#> attr(,"class")
#> [1] "geo_list"
```

Random linestring - starting from a random point, with default maximum segment
length and maximum rotation between two segments.


```r
geo_linestring()
#> $type
#> [1] "FeatureCollection"
#> 
#> $features
#> $features[[1]]
#> $features[[1]]$type
#> [1] "Feature"
#> 
#> $features[[1]]$geometry
#> $features[[1]]$geometry$type
#> [1] "LineString"
#> 
#> $features[[1]]$geometry$coordinates
#> $features[[1]]$geometry$coordinates[[1]]
#> $features[[1]]$geometry$coordinates[[1]][[1]]
#> [1] 51.028387 -2.188767
#> 
#> $features[[1]]$geometry$coordinates[[1]][[2]]
#> [1] 51.027774 -2.189005
#> 
#> $features[[1]]$geometry$coordinates[[1]][[3]]
#> [1] 51.028109 -2.189318
#> 
#> $features[[1]]$geometry$coordinates[[1]][[4]]
#> [1] 51.027487 -2.190017
#> 
#> $features[[1]]$geometry$coordinates[[1]][[5]]
#> [1] 51.027774 -2.190379
#> 
#> $features[[1]]$geometry$coordinates[[1]][[6]]
#> [1] 51.027088 -2.191078
#> 
#> $features[[1]]$geometry$coordinates[[1]][[7]]
#> [1] 51.027435 -2.191403
#> 
#> $features[[1]]$geometry$coordinates[[1]][[8]]
#> [1] 51.026922 -2.192147
#> 
#> $features[[1]]$geometry$coordinates[[1]][[9]]
#> [1] 51.027533 -2.192925
#> 
#> $features[[1]]$geometry$coordinates[[1]][[10]]
#> [1] 51.027475 -2.192984
#> 
#> 
#> 
#> 
#> $features[[1]]$properties
#> NULL
#> 
#> 
#> 
#> attr(,"class")
#> [1] "geo_list"
```

Random polygon - centered on a random point, with default maximum size


```r
geo_polygon()
#> $type
#> [1] "FeatureCollection"
#> 
#> $features
#> $features[[1]]
#> $features[[1]]$type
#> [1] "Feature"
#> 
#> $features[[1]]$geometry
#> $features[[1]]$geometry$type
#> [1] "Polygon"
#> 
#> $features[[1]]$geometry$coordinates
#> $features[[1]]$geometry$coordinates[[1]]
#> $features[[1]]$geometry$coordinates[[1]][[1]]
#> [1]  7.375542 21.481953
#> 
#> $features[[1]]$geometry$coordinates[[1]][[2]]
#> [1]  9.557249 13.822460
#> 
#> $features[[1]]$geometry$coordinates[[1]][[3]]
#> [1] 10.882944  6.553386
#> 
#> $features[[1]]$geometry$coordinates[[1]][[4]]
#> [1] 8.192794 5.962077
#> 
#> $features[[1]]$geometry$coordinates[[1]][[5]]
#> [1] 8.305398 5.210405
#> 
#> $features[[1]]$geometry$coordinates[[1]][[6]]
#> [1] 2.573744 9.711262
#> 
#> $features[[1]]$geometry$coordinates[[1]][[7]]
#> [1]  0.4559409 17.8573283
#> 
#> $features[[1]]$geometry$coordinates[[1]][[8]]
#> [1]  5.093829 12.718010
#> 
#> $features[[1]]$geometry$coordinates[[1]][[9]]
#> [1]  2.778335 20.708271
#> 
#> $features[[1]]$geometry$coordinates[[1]][[10]]
#> [1]  5.103798 12.757463
#> 
#> $features[[1]]$geometry$coordinates[[1]][[11]]
#> [1]  7.375542 21.481953
#> 
#> 
#> 
#> 
#> $features[[1]]$properties
#> NULL
#> 
#> 
#> 
#> attr(,"class")
#> [1] "geo_list"
```

Visualize your shapes with **lawn**.


```r
lawn::view(jsonlite::toJSON(geo_polygon(count = 4), auto_unbox = TRUE))
```

![map](https://github.com/ropensci/randgeo/raw/master/tools/plot.png)


## Generate random WKT

Random point


```r
wkt_point()
#> [1] "POINT (50.3923570 -70.3787919)"
```

Random linestring


```r
wkt_linestring()
#> [1] "LINESTRING (42.7817546 19.4598109, 42.7818265 19.4597713, 42.7819039 19.4597543, 42.7819206 19.4597476, 42.7819732 19.4597341, 42.7820674 19.4596920, 42.7821133 19.4596554, 42.7821233 19.4596390, 42.7821761 19.4595728, 42.7821867 19.4595509)"
```

Random polygon


```r
wkt_polygon()
#> [1] "POLYGON ((-164.5191674 48.1392788, -160.7019535 50.0321237, -162.9547961 47.6318361, -154.1477728 46.2428280, -164.4173904 45.9485672, -160.5131717 43.5274697, -167.1833639 41.2046996, -164.5419419 45.9658491, -171.4209181 42.8764619, -166.7761058 46.2084577, -164.5191674 48.1392788))"
```

## Contributors

* Scott Chamberlain
* Noam Ross
* Samuel Bosch

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/randgeo/issues).
* License: MIT
* Get citation information for `randgeo` in R doing `citation(package = 'randgeo')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
randgeo 0.3.0
=============

### NEW FEATURES

* Gains `geo_linestring` and `wkt_linestring` ([#13](https://github.com/ropensci/randgeo/pull/13)) thanks to @samuelbosch
* All `geo_*` functions output lists - they all gain a S3 class `geo_list` to be in line with `geojsonio` outputs

### MINOR IMPROVEMENTS

* fix placement of image used in readme that's CRAN compliant


randgeo 0.2.0
=============

### NEW FEATURES

All changes thanks to @noamross, and come from (#6) (#7)

Internals of the package are re-worked to have the features:

* Random points are now evenly distributed on the sphere
* Random shape vertex distances from the center are now based on equal-scale 
great-circle arcs
* Everything is in units of degrees latitude (e.g., ~69 miles or 111 km)

### MINOR IMPROVEMENTS

* More tests were added


randgeo 0.1.0
=============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 3.5.0
* ubuntu 12.04 (on travis-ci), R 3.5.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

  License components with restrictions and base license permitting such:
    MIT + file LICENSE
  File 'LICENSE':
    YEAR: 2018
    COPYRIGHT HOLDER: Scott Chamberlain

## Reverse dependencies

There are no reverse dependencies.

---

This version adds two new functions `geo_linestring` and `wkt_linestring` 
to make randgeom linestrings as geojson or wkt, respectively. Also fixes 
placement of images used in readme.

Thanks!
Scott Chamberlain
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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/randgeo/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/randgeo.git`
* Make sure to track progress upstream (i.e., on our version of `randgeo` at `ropensci/randgeo`) by doing `git remote add upstream https://github.com/ropensci/randgeo.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/randgeo`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
