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
mregions
========

[![cran checks](https://cranchecks.info/badges/worst/mregions)](https://cranchecks.info/pkgs/mregions)
[![Build Status](https://travis-ci.org/ropensci/mregions.svg)](https://travis-ci.org/ropensci/mregions)
[![codecov.io](https://codecov.io/github/ropensci/mregions/coverage.svg?branch=master)](https://codecov.io/github/ropensci/mregions?branch=master)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/mregions?color=FAB657)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/mregions)](https://cran.r-project.org/package=mregions)
[![](https://badges.ropensci.org/53_status.svg)](https://github.com/ropensci/onboarding/issues/53)

`mregions` - Get data from <http://www.marineregions.org>

Some data comes from the [Flanders Marine Institute (VLIZ) geoserver](http://geo.vliz.be/geoserver/web/)

`mregions` is useful to a wide diversity of R users because you get access to all of the
data MarineRegions has, which can help in a variety of use cases:

* Visualize marine regions alone
* Visualize marine regions with associated data paired with analysis
* Use marine region geospatial boundaries to query data providers (e.g., OBIS (<http://www.iobis.org>))
* Geocode - get geolocation data from place names
* Reverse Geocode - get place names from geolocation data

## Install


```r
install.packages("mregions")
```

Development version


```r
devtools::install_github("ropensci/mregions")
```


```r
library("mregions")
```

## GeoJSON

Get region


```r
res1 <- mr_geojson(key = "Morocco:dam")
```

Get helper library


```r
install.packages("leaflet")
```

Plot data


```r
library('leaflet')
leaflet() %>%
  addProviderTiles(provider = 'OpenStreetMap') %>%
  addGeoJSON(geojson = res1$features) %>%
  setView(-3.9, 35, zoom = 10)
```

![map](tools/img/leaf1.png)

## Shape

Get region


```r
res2 <- mr_shp(key = "MarineRegions:eez_iho_union_v2", maxFeatures = 5)
```

Get helper library


```r
install.packages("leaflet")
```

Plot data


```r
library('leaflet')
leaflet() %>%
  addProviderTiles(provider = 'OpenStreetMap') %>%
  addPolygons(data = res2)
```

![map2](tools/img/leaf2.png)

## Convert to WKT

From geojson


```r
res3 <- mr_geojson(key = "Morocco:dam")
mr_as_wkt(res3, fmt = 5)

#> [1] "MULTIPOLYGON (((41.573732 -1.659444, 45.891882 ... cutoff
```

From shp object (`SpatialPolygonsDataFrame`) or file, both work


```r
mr_as_wkt(mr_shp(key = "MarineRegions:eez_iho_union_v2"))

#> [1] "GEOMETRYCOLLECTION (POLYGON ((-7.25 ... cutoff
```

## Get OBIS EEZ ID


```r
mr_obis_eez_id("bulgarian exclusive economic zone")
```

```
## [1] 71
```

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Francois Michonneau](https://github.com/fmichonneau)
* [Pieter Provoost](https://github.com/pieterprovoost)
* [Michael Sumner](https://github.com/mdsumner)
* [Lennert Schepers](https://github.com/LennertSchepers)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/mregions/issues).
* License: MIT
* Get citation information for `mregions` in R doing `citation(package = 'mregions')`
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

[![rofooter](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
mregions 0.1.6
==============

### BUG FIXES

* Bug fixes for slight changes in the Marineregions web services


mregions 0.1.4
==============

<http://marineregions.org/> changed some of their services. Thus,
we had to change the way some functions work, remove functionality
of some parameters, and add new functions. We are still in the process
of making it all work smoothly.

### NEW FEATURES

* new function `mr_features_get()` to fetch features of many different
data types, including geojson, shp, kml, and more.
* new function `mr_layers()` to list layers

### MINOR IMPROVEMENTS

* Tidying man pages throughout package to reduce down to 80 line width
* `mr_names()` changed behavior. You used to be able to get all names/regions,
but now you have to specify a layer you want information for.
* `mr_names_search()` internals changed to account for changes in
`mr_names()`.



mregions 0.1.0
==============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 3.4.2 patched
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

  License components with restrictions and base license permitting such:
     MIT + file LICENSE
   File 'LICENSE':
     YEAR: 2017
     COPYRIGHT HOLDER: Scott Chamberlain

## Reverse dependencies

There are no reverse dependencies.

--------

This version includes some bug fixes for changes in the webservices this
package works with.

Thanks!
Scott Chamberlain
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{mregions introduction}
%\VignetteEncoding{UTF-8}
-->



mregions introduction
=====================

`mregions` is useful to a wide diversity of R users because you get access to all of the
data MarineRegions has, which can help in a variety of use cases:

* Visualize marine regions alone
* Visualize marine regions with associated data paired with analysis
* Use marine region geospatial boundaries to query data providers (e.g., OBIS (<http://www.iobis.org>))
* Geocode - get geolocation data from place names
* Reverse Geocode - get place names from geolocation data

## Install

Stable version


```r
install.packages("mregions")
```

Dev version


```r
devtools::install_github("ropensci/mregions")
install.packages("leaflet")
```


```r
library("mregions")
```

## Get list of place types


```r
res <- mr_place_types()
head(res$type)
#> [1] "Town"                      "Arrondissement"           
#> [3] "Department"                "Province (administrative)"
#> [5] "Country"                   "Continent"
```

## Get Marineregions records by place type


```r
res1 <- mr_records_by_type(type = "EEZ")
head(res1)
#>   MRGID
#> 1  3293
#> 2  5668
#> 3  5669
#> 4  5670
#> 5  5672
#> 6  5673
#>                                                                                                                                                                                                             gazetteerSource
#> 1 Flanders Marine Institute (2016). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 9. Available online at http://www.marineregions.org/. http://dx.doi.org/10.14284/242
#> 2 Flanders Marine Institute (2016). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 9. Available online at http://www.marineregions.org/. http://dx.doi.org/10.14284/242
#> 3 Flanders Marine Institute (2016). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 9. Available online at http://www.marineregions.org/. http://dx.doi.org/10.14284/242
#> 4 Flanders Marine Institute (2016). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 9. Available online at http://www.marineregions.org/. http://dx.doi.org/10.14284/242
#> 5 Flanders Marine Institute (2016). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 9. Available online at http://www.marineregions.org/. http://dx.doi.org/10.14284/242
#> 6 Flanders Marine Institute (2016). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 9. Available online at http://www.marineregions.org/. http://dx.doi.org/10.14284/242
#>   placeType latitude longitude minLatitude minLongitude maxLatitude
#> 1       EEZ 51.46483  2.704458    51.09111     2.238118    51.87000
#> 2       EEZ 53.61508  4.190675    51.26203     2.539443    55.76500
#> 3       EEZ 54.55970  8.389231    53.24281     3.349999    55.91928
#> 4       EEZ 40.87030 19.147094    39.63863    18.461940    41.86124
#> 5       EEZ 42.94272 29.219062    41.97820    27.449580    43.74779
#> 6       EEZ 43.42847 15.650844    41.62201    13.001390    45.59079
#>   maxLongitude precision            preferredGazetteerName
#> 1     3.364907  58302.49   Belgian Exclusive Economic Zone
#> 2     7.208364 294046.10     Dutch Exclusive Economic Zone
#> 3    14.750000 395845.50    German Exclusive Economic Zone
#> 4    20.010030 139751.70  Albanian Exclusive Economic Zone
#> 5    31.345280 186792.50 Bulgarian Exclusive Economic Zone
#> 6    18.552360 313990.30  Croatian Exclusive Economic Zone
#>   preferredGazetteerNameLang   status accepted
#> 1                    English standard     3293
#> 2                    English standard     5668
#> 3                    English standard     5669
#> 4                    English standard     5670
#> 5                    English standard     5672
#> 6                    English standard     5673
```

## Get a data.frame of region names


```r
rnames <- mr_names("MarineRegions:iho")
```

## Search region names

Either pass output of `mr_names()`


```r
mr_names_search(rnames, "IHO")
#> # A tibble: 7 x 6
#>               layer    name_first name_second     id
#>               <chr>         <chr>       <chr>  <chr>
#> 1 MarineRegions:iho MarineRegions         iho  iho.1
#> 2 MarineRegions:iho MarineRegions         iho  iho.7
#> 3 MarineRegions:iho MarineRegions         iho iho.18
#> 4 MarineRegions:iho MarineRegions         iho iho.40
#> 5 MarineRegions:iho MarineRegions         iho iho.53
#> 6 MarineRegions:iho MarineRegions         iho iho.76
#> 7 MarineRegions:iho MarineRegions         iho iho.94
#> # ... with 2 more variables: name <chr>, mrgid <chr>
```

or don't (but then `mr_names_search()` call takes longer)


```r
mr_names_search("iho", q = "Sea")
#> # A tibble: 73 x 6
#>                layer    name_first name_second     id
#>                <chr>         <chr>       <chr>  <chr>
#>  1 MarineRegions:iho MarineRegions         iho  iho.3
#>  2 MarineRegions:iho MarineRegions         iho  iho.4
#>  3 MarineRegions:iho MarineRegions         iho  iho.6
#>  4 MarineRegions:iho MarineRegions         iho  iho.7
#>  5 MarineRegions:iho MarineRegions         iho  iho.8
#>  6 MarineRegions:iho MarineRegions         iho iho.10
#>  7 MarineRegions:iho MarineRegions         iho iho.15
#>  8 MarineRegions:iho MarineRegions         iho iho.16
#>  9 MarineRegions:iho MarineRegions         iho iho.17
#> 10 MarineRegions:iho MarineRegions         iho iho.27
#> # ... with 63 more rows, and 2 more variables: name <chr>, mrgid <chr>
```

## Get a region - geojson


```r
res3 <- mr_geojson(key = "Morocco:dam")
class(res3)
#> [1] "mr_geojson"
names(res3)
#> [1] "type"          "totalFeatures" "features"      "crs"
```

## Get a region - shp


```r
res4 <- mr_shp(key = "Morocco:dam")
class(res4)
#> [1] "SpatialPolygonsDataFrame"
#> attr(,"package")
#> [1] "sp"
```

## Get OBIS EEZ ID


```r
mr_obis_eez_id("Bulgarian Exclusive Economic Zone")
#> [1] 71
```

## Convert to WKT

From geojson or shp. Here, geojson


```r
res7 <- mr_geojson(key = "Morocco:dam")
mr_as_wkt(res7, fmt = 5)
#> [1] "MULTIPOLYGON (((41.573732 -1.659444, 45.891882 ... cutoff
```

## Dealing with bigger WKT

What if you're WKT string is super long?  It's often a problem because some online species occurrence databases that accept WKT to search by geometry bork due to
limitations on length of URLs if your WKT string is too long (about 8000 characters,
including remainder of URL). One way to deal with it is to reduce detail - simplify.


```r
install.packages("rmapshaper")
```

Using `rmapshaper` we can simplify a spatial object, then search with that.


```r
shp <- mr_shp(key = "MarineRegions:eez_iho_union_v2", maxFeatures = 5)
```

Visualize


```r
library(leaflet)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = shp)
```

![map2](figure/complex.png)

Simplify


```r
library("rmapshaper")
shp <- ms_simplify(shp)
```

It's simplified:


```r
library(leaflet)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = shp)
```

![map3](figure/simple.png)

[mr]: https://github.com/ropensci/mregions
