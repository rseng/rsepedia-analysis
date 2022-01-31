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

```{r eval=FALSE}
install.packages("mregions")
```

Development version

```{r eval=FALSE}
devtools::install_github("ropensci/mregions")
```

```{r}
library("mregions")
```

## GeoJSON

Get region

```{r}
res1 <- mr_geojson(key = "Morocco:dam")
```

Get helper library

```{r eval=FALSE}
install.packages("leaflet")
```

Plot data

```{r eval=FALSE}
library('leaflet')
leaflet() %>%
  addProviderTiles(provider = 'OpenStreetMap') %>%
  addGeoJSON(geojson = res1$features) %>%
  setView(-3.9, 35, zoom = 10)
```

![map](tools/img/leaf1.png)

## Shape

Get region

```{r}
res2 <- mr_shp(key = "MarineRegions:eez_iho_union_v2", maxFeatures = 5)
```

Get helper library

```{r eval=FALSE}
install.packages("leaflet")
```

Plot data

```{r eval=FALSE}
library('leaflet')
leaflet() %>%
  addProviderTiles(provider = 'OpenStreetMap') %>%
  addPolygons(data = res2)
```

![map2](tools/img/leaf2.png)

## Convert to WKT

From geojson

```{r eval=FALSE}
res3 <- mr_geojson(key = "Morocco:dam")
mr_as_wkt(res3, fmt = 5)

#> [1] "MULTIPOLYGON (((41.573732 -1.659444, 45.891882 ... cutoff
```

From shp object (`SpatialPolygonsDataFrame`) or file, both work

```{r eval=FALSE}
mr_as_wkt(mr_shp(key = "MarineRegions:eez_iho_union_v2"))

#> [1] "GEOMETRYCOLLECTION (POLYGON ((-7.25 ... cutoff
```

## Get OBIS EEZ ID

```{r}
mr_obis_eez_id("bulgarian exclusive economic zone")
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
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{mregions introduction}
%\VignetteEncoding{UTF-8}
-->

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("mregions")
```

Dev version

```{r eval=FALSE}
devtools::install_github("ropensci/mregions")
install.packages("leaflet")
```

```{r}
library("mregions")
```

## Get list of place types

```{r}
res <- mr_place_types()
head(res$type)
```

## Get Marineregions records by place type

```{r}
res1 <- mr_records_by_type(type = "EEZ")
head(res1)
```

## Get a data.frame of region names

```{r}
rnames <- mr_names("MarineRegions:iho")
```

## Search region names

Either pass output of `mr_names()`

```{r}
mr_names_search(rnames, "IHO")
```

or don't (but then `mr_names_search()` call takes longer)

```{r}
mr_names_search("iho", q = "Sea")
```

## Get a region - geojson

```{r}
res3 <- mr_geojson(key = "Morocco:dam")
class(res3)
names(res3)
```

## Get a region - shp

```{r}
res4 <- mr_shp(key = "Morocco:dam")
class(res4)
```

## Get OBIS EEZ ID

```{r}
mr_obis_eez_id("Bulgarian Exclusive Economic Zone")
```

## Convert to WKT

From geojson or shp. Here, geojson

```{r eval=FALSE}
res7 <- mr_geojson(key = "Morocco:dam")
mr_as_wkt(res7, fmt = 5)
#> [1] "MULTIPOLYGON (((41.573732 -1.659444, 45.891882 ... cutoff
```

## Dealing with bigger WKT

What if you're WKT string is super long?  It's often a problem because some online species occurrence databases that accept WKT to search by geometry bork due to
limitations on length of URLs if your WKT string is too long (about 8000 characters,
including remainder of URL). One way to deal with it is to reduce detail - simplify.

```{r eval=FALSE}
install.packages("rmapshaper")
```

Using `rmapshaper` we can simplify a spatial object, then search with that.

```{r}
shp <- mr_shp(key = "MarineRegions:eez_iho_union_v2", maxFeatures = 5)
```

Visualize

```{r eval=FALSE}
library(leaflet)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = shp)
```

![map2](figure/complex.png)

Simplify

```{r}
library("rmapshaper")
shp <- ms_simplify(shp)
```

It's simplified:

```{r eval=FALSE}
library(leaflet)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = shp)
```

![map3](figure/simple.png)

[mr]: https://github.com/ropensci/mregions
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/region_names_search.R
\name{mr_names_search}
\alias{mr_names_search}
\title{Search for region names}
\usage{
mr_names_search(x, q = NULL, ...)
}
\arguments{
\item{x, q}{Either a \code{tbl_df}, returned from \code{\link[=mr_names]{mr_names()}}, or
a query as a character string. If a \code{tbl_df}, you must pass a query
string to \code{q}. If a query string (character) is passed to \code{x},
leave \code{q} as \code{NULL}}

\item{...}{Parameters passed on to \code{\link[=agrep]{agrep()}}}
}
\value{
\code{NULL} if no matches found, or a data.frame, or tibble, of class
\code{tbl_df}, with slots:
\itemize{
\item name (character) - name of the region, which is a combination of the
name_first and name_second, e.g., Morocco:elevation_10m
\item title (character) - title for the region
\item name_first (character) - first part of the name, e.g., Morocco
\item name_second (character) - second part of the name, e.g., elevation_10m
}
}
\description{
Search for region names
}
\examples{
\dontrun{
# Get region names with mr_names() function
(res <- mr_names("MarineRegions:eez"))

# to save time, pass in the result from mr_names()
mr_names_search(res, q = "Amer")

# if you don't pass in the result from mr_names(), we have to
# call mr_names() internally, adding some time
mr_names_search(x = "iho", q = "Black")
mr_names_search(x = "iho", q = "Sea")

# more examples
mr_names_search("iho", "Sea")
(res <- mr_names("MarineRegions:iho"))
mr_names_search(res, q = "Sea")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_wkt.R
\name{mr_as_wkt}
\alias{mr_as_wkt}
\title{Convert data to WKT}
\usage{
mr_as_wkt(x, fmt = 16, ...)
}
\arguments{
\item{x}{Output from \code{\link[=mr_geojson]{mr_geojson()}}, \code{\link[=mr_shp]{mr_shp()}},
or a \code{SpatialPolygonsDataFrame}}

\item{fmt}{(integer) The number of digits to display after the decimal
point when formatting coordinates. Ignored when shp files or
\code{SpatialPolygonsDataFrame} passed in}

\item{...}{Further args passed on to \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} only
in the event of json passed as a character string. Ignored when shp files or
\code{SpatialPolygonsDataFrame} passed in}
}
\value{
a character string of WKT data
}
\description{
Convert data to WKT
}
\details{
WKT, or Well Known Text, is a way to encode spatial data. It's
somewhat similar to GeoJSON, but instead of being in JSON format, it's a
character string (though can also be encoded in binary format). WKT is
often used in SQL databases, and many species occurrence APIs allow only
WKT. You could do the conversion to WKT yourself, but we provide
\code{as_wkt} as a convenience
}
\examples{
\dontrun{
res <- mr_geojson(key = "Morocco:dam")
mr_as_wkt(res, fmt = 5)

# shp files
## path to wkt
mr_as_wkt(mr_shp(key = "Morocco:dam", read = FALSE))

## spatial object to wkt
mr_as_wkt(mr_shp(key = "Morocco:dam", read = TRUE))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mregions-package.R
\docType{package}
\name{mregions-package}
\alias{mregions-package}
\alias{mregions}
\title{Marine regions data from Marineregions}
\description{
Tools to get marine regions data from
\url{https://www.marineregions.org/}. Includes tools to get region metadata,
as well as data in 'GeoJSON' format, as well as Shape files. Use cases
include using data downstream to visualize 'geospatial' data by marine
region, mapping variation among different regions, and more.
}
\details{
mregions gets data from \url{http://www.marineregions.org/}
}
\section{Use-cases}{

\pkg{mregions} is useful to a wide diversity of R users because you get
access to all of the data MarineRegions has, which can help in a
variety of use cases:
\itemize{
\item Visualize marine regions alone
\item Visualize marine regions with associated data paired with analysis
\item Use marine region geospatial boundaries to query data providers
(e.g., OBIS (\url{http://www.iobis.org}))
\item Geocode - get geolocation data from place names
\item Reverse Geocode - get place names from geolocation data
}
}

\examples{
\dontrun{
## GeoJSON
### Get region
res <- mr_geojson(key = "Morocco:dam")

### Plot data
if (!requireNamespace("leaflet")) {
 install.packages("leaflet")
}
library('leaflet')
leaflet() \%>\%
  addProviderTiles(provider = 'OpenStreetMap') \%>\%
  addGeoJSON(geojson = res$features) \%>\%
  setView(-3.98, 35.1, zoom = 11)

## Shape
### Get region
res <- mr_shp(key = "MarineRegions:eez_iho_union_v2")
library('leaflet')
leaflet() \%>\%
  addProviderTiles(provider = 'OpenStreetMap') \%>\%
  addPolygons(data = res)

## Convert to WKT
### From geojson
res <- mr_geojson(key = "Morocco:dam")
mr_as_wkt(res, fmt = 5)

### From shp object (`SpatialPolygonsDataFrame`) or file, both work
mr_as_wkt(mr_shp(key = "Morocco:dam", read = FALSE))
## spatial object to wkt
mr_as_wkt(mr_shp(key = "Morocco:dam", read = TRUE))
}
}
\author{
Scott Chamberlain

Francois Michonneau

Pieter Provoost

Michael Sumner

Lennert Schepers \email{lennert.schepers@vliz.be}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/region_geojson.R
\name{mr_geojson}
\alias{mr_geojson}
\title{Get a Marineregions geojson file}
\usage{
mr_geojson(key = NULL, name = NULL, maxFeatures = 50, ...)
}
\arguments{
\item{key}{(character) Region key, of the form \code{x:y}, where
\code{x} is a namespace (e.g., \code{MarineRegions}), and \code{y} is
a region (e.g., \code{eez_33176})}

\item{name}{(character) Region name, if you supply this, we search
against titles via \code{\link[=mr_names]{mr_names()}} function}

\item{maxFeatures}{(integer) Number of features to return. Default: \code{50}}

\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\value{
an S3 class of type \code{mr_geojson}, just a thin wrapper around
a list. The list has names:
\itemize{
\item type (character) - the geojson type (e.g., FeatureCollection)
\item totalFeatures (integer) - the
\item features (list) - the features, with slots for each feature: type,
id, geometry, geometry_name, and properties
\item crs (list) - the coordinate reference system
\item bbox (list) - the bounding box that encapsulates the object
}
}
\description{
Get a Marineregions geojson file
}
\examples{
\dontrun{
# by key
res1 <- mr_geojson(key = "Morocco:dam")

# by name -- not working right now

if (requireNamespace("geojsonio")) {
  library("geojsonio")
  as.json(unclass(res1)) \%>\% map_leaf

  # MEOW - marine ecoregions
  as.json(unclass(mr_geojson("Ecoregions:ecoregions"))) \%>\% map_leaf()
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_features_get.R
\name{mr_features_get}
\alias{mr_features_get}
\title{Get features}
\usage{
mr_features_get(
  type,
  featureID,
  maxFeatures = 100,
  format = "json",
  path = NULL,
  version = "2.0.0",
  ...
)
}
\arguments{
\item{type}{(character) a region type, e.g., "MarineRegions:eez". required}

\item{featureID}{(character) a feature ID. required}

\item{maxFeatures}{(integer) Number of features. Default: 100}

\item{format}{(character) output format, see Details for allowed options.
Default: json}

\item{path}{(character) required when \code{format="SHAPE-ZIP"},
otherwise, ignored}

\item{version}{(character) either 1.0.0 or 2.0.0 (default). In v1.0.0, the
coordinates are in format y,x (long,lat), while in 2.0.0 the coordinates
are in format x,y (lat,long)}

\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\value{
depends on the \code{format} option used, usually a text string
}
\description{
Get features
}
\details{
Allowed options for the \code{format} parameter:
\itemize{
\item \verb{text/xml; subtype=gml/3.2}
\item \code{GML2}
\item \code{KML}
\item \code{SHAPE-ZIP}
\item \verb{application/gml+xml; version=3.2}
\item \code{application/json}
\item \verb{application/vnd.google-earth.kml xml}
\item \code{application/vnd.google-earth.kml+xml}
\item \code{csv}
\item \code{gml3}
\item \code{gml32}
\item \code{json}
\item \verb{text/xml; subtype=gml/2.1.2}
\item \verb{text/xml; subtype=gml/3.1.1}
}
}
\examples{
\dontrun{
# json by default
mr_features_get(type = "MarineRegions:eez", featureID = "eez.3")
# csv
mr_features_get(type = "MarineRegions:eez", featureID = "eez.3",
  format = "csv")
# KML
mr_features_get(type = "MarineRegions:eez", featureID = "eez.3",
  format = "KML")

# if you want SHAPE-ZIP, give a file path
# FIXME - shape files not working right now
# file <- tempfile(fileext = ".zip")
# mr_features_get(type = "MarineRegions:eez", featureID = "eez.3",
#   format = "SHAPE-ZIP", path = file)
# file.exists(file)
# unlink(file)

# glm32
mr_features_get(type = "MarineRegions:eez", featureID = "eez.3",
  format = "gml32")

# version parameter
## notice the reversed coordinates
mr_features_get(type = "MarineRegions:eez", featureID = "eez.3")
mr_features_get(type = "MarineRegions:eez", featureID = "eez.3",
  version = "1.0.0")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/records_by_type.R
\name{mr_records_by_type}
\alias{mr_records_by_type}
\title{Get Marineregions records by place type}
\usage{
mr_records_by_type(type, offset = 0, ...)
}
\arguments{
\item{type}{(character) One place type name. See
\code{\link[=mr_place_types]{mr_place_types()}} for place type names}

\item{offset}{(numeric) Offset to start at. Each request can return up to
100 results. e.g., an offset of 200 will give records 200 to 299.}

\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\value{
If no results, an empty list. If results found, a data.frame with the columns:
\itemize{
\item MRGID (integer)
\item gazetteerSource (character)
\item placeType (character)
\item latitude (numeric)
\item longitude (numeric)
\item minLatitude (numeric)
\item minLongitude (numeric)
\item maxLatitude (numeric)
\item maxLongitude (numeric)
\item precision (numeric)
\item preferredGazetteerName (character)
\item preferredGazetteerNameLang (character)
\item status (character)
\item accepted (integer)
}
}
\description{
Get Marineregions records by place type
}
\details{
Internally we use the \code{getGazetteerRecordsByType.json} API
method, which searches for Marineregions records by user supplied place type
}
\examples{
\dontrun{
# Get records of type 'EEZ', then inspect data.frame
res <- mr_records_by_type(type="EEZ")
head(res)

# You can use mr_place_types() function to get types
## then pass those into this function
types <- mr_place_types()
mr_records_by_type(types$type[1])
mr_records_by_type(types$type[10])

# use regex to find a type name matching a pattern
x <- grep("MEOW", types$type, value = TRUE)

# then pass to the function
mr_records_by_type(x)
mr_records_by_type(x, offset = 100)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rev_geo_code.R
\name{mr_rev_geo_code}
\alias{mr_rev_geo_code}
\title{Reverse Geocode with Marineregions}
\usage{
mr_rev_geo_code(lat, lon, lat_radius = 1, lon_radius = 1, ...)
}
\arguments{
\item{lat}{(numeric) Latitude for the coordinates (decimal format)}

\item{lon}{(numeric) Longitude for the coordinates (decimal format)}

\item{lat_radius}{(numeric) Extends search to include the range from
\code{lat}-\code{lat_radius} to \code{lat}+\code{lat_radius}}

\item{lon_radius}{(numeric) Extends search to include the range from
\code{lon}-\code{lon_radius} to \code{lon}+\code{lon_radius}}

\item{...}{curl options to be passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\value{
If no results, an empty list. If results found, a data.frame with the columns:
\itemize{
\item MRGID (integer)
\item gazetteerSource (character)
\item placeType (character)
\item latitude (numeric)
\item longitude (numeric)
\item minLatitude (numeric)
\item minLongitude (numeric)
\item maxLatitude (numeric)
\item maxLongitude (numeric)
\item precision (numeric)
\item preferredGazetteerName (character)
\item preferredGazetteerNameLang (character)
\item status (character)
\item accepted (integer)
}
}
\description{
Retrieve the names of geographic objects from coordinates (and
optionally a radius around them).
}
\examples{
\dontrun{
# Setting radius to 0.5
mr_rev_geo_code(-21.5, 55.5, lat_radius=0.5, lon_radius=0.5)

# radius to 3
mr_rev_geo_code(-21.5, 55.5, lat_radius=3, lon_radius=3)

# radius to 1
mr_rev_geo_code(-15, 45, lat_radius=1, lon_radius=1)
}
}
\author{
Francois Michonneau \href{mailto:francois.michonneau@gmail.com}{francois.michonneau@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/place_types.R
\name{mr_place_types}
\alias{mr_place_types}
\title{Get Marineregions place types}
\usage{
mr_place_types(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\value{
A data.frame with the columns:
\itemize{
\item type (character) the place type
\item description (character) description of the place type
}
}
\description{
Get Marineregions place types
}
\examples{
\dontrun{
res <- mr_place_types()
head(res)
res$type
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_names2.R
\name{mr_names}
\alias{mr_names}
\title{Get region names - v2}
\usage{
mr_names(layer, ...)
}
\arguments{
\item{layer}{A layer name, one of MarineRegions:eez,
MarineRegions:eez_boundaries, MarineRegions:iho, MarineRegions:fao,
or MarineRegions:lme}

\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\value{
a data.frame, or tibble, of class tbl_df (basically, a compact
data.frame), with slots:
\itemize{
\item layer (character) - name of the layer (e.g. MarineRegions:eez)
\item name_first (character) - first part of the name, e.g., MarineRegions
\item name_second (character) - second part of the name, e.g., eez
\item id (character) - the feature ID
}

additional columns vary by layer
}
\description{
Get region names - v2
}
\examples{
\dontrun{
# mr_names gives a tidy data.frame
(res <- mr_names("MarineRegions:eez"))
(res <- mr_names('MarineRegions:eez_boundaries'))
(res <- mr_names('MarineRegions:iho'))
(res <- mr_names('MarineRegions:fao'))
(res <- mr_names('MarineRegions:lme'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/region_shp.R
\name{mr_shp}
\alias{mr_shp}
\title{Get a region shp file}
\usage{
mr_shp(
  key = NULL,
  name = NULL,
  maxFeatures = 500,
  overwrite = TRUE,
  read = TRUE,
  filter = NULL,
  ...
)
}
\arguments{
\item{key}{(character) Region key, of the form \code{x:y}, where
\code{x} is a namespace (e.g., \code{MarineRegions}), and \code{y} is
a region (e.g., \code{eez_33176})}

\item{name}{(character) Region name, if you supply this, we search
against titles via \code{\link[=mr_names]{mr_names()}} function}

\item{maxFeatures}{(integer) Number of features}

\item{overwrite}{(logical) Overwrite file if already exists.
Default: \code{FALSE}}

\item{read}{(logical) To read in as spatial object. If \code{FALSE} a path
given back. if \code{TRUE}, you need the \code{rgdal} package installed.
Default: \code{FALSE}}

\item{filter}{(character) String to filter features on}

\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}. since we
use caching, note that if you've made the exact same request before and the
file is still in cache, we grab the cached file and don't make an HTTP
request, so any curl options passed would be ignored.}
}
\value{
A \code{SpatialPolygonsDataFrame} if \code{read = TRUE}, or a path to
a SHP file on disk if \code{read = FALSE}.
}
\description{
Get a region shp file
}
\details{
We use \pkg{rappdirs} to determine where to cache data depening on
your operating system. See \code{rappdirs::user_cache_dir("mregions")} for
location on your machine

We cache based on the name of the region plus the \code{maxFeatures}
parameter. That is to say, you can query the same region name, but
with different \code{maxFeatures} parameter values, and they will get
cached separately. You can clear the cache by going to the directory at
\code{rappdirs::user_cache_dir("mregions")} and deleting the files.

We use \code{stringsAsFactors = FALSE} inside of \code{rgdal::readOGR()}
so that character variables aren't converted to factors.
}
\note{
the parameter \code{name} is temporarily not useable. MarineRegions
updated their web services, and we haven't sorted out yet how to make
this feature work. We may bring it back in future version of this pacakge.
}
\examples{
\dontrun{
## just get path
mr_shp(key = "MarineRegions:eez_iho_union_v2", read = FALSE)
## read shp file into spatial object
res <- mr_shp(key = "MarineRegions:eez_iho_union_v2", read = TRUE)

mr_shp(key = "SAIL:w_marinehabitatd")

# maxFeatures
library(sp)
plot(mr_shp(key = "MarineRegions:eez_iho_union_v2"))
plot(mr_shp(key = "MarineRegions:eez_iho_union_v2", maxFeatures = 5))

# vizualize with package leaflet
if (requireNamespace("leaflet")) {
  library('leaflet')
  leaflet() \%>\%
    addTiles() \%>\%
    addPolygons(data = res)
}

# use `filter` param to get a subset of a region
library(sp)
pp <- mr_shp(key = "MarineRegions:eez_iho_union_v2")
plot(pp)
rr <- mr_shp(key = "MarineRegions:eez_iho_union_v2",
  filter = "North Atlantic Ocean")
plot(rr)

# get Samoan Exclusive Economic Zone
res <- mr_shp(
  key = "MarineRegions:eez",
  filter = "Samoan Exclusive Economic Zone"
)
sp::plot(res)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/place_relations.R
\name{mr_place_relations}
\alias{mr_place_relations}
\title{Related records}
\usage{
mr_place_relations(
  mrgid,
  direction = c("upper", "lower", "both"),
  type = c("partof", "partlypartof", "adjacentto", "similarto", "administrativepartof",
    "influencedby", "all"),
  ...
)
}
\arguments{
\item{mrgid}{(numeric) the MRGID (Marineregions Global Identifier) for the
record of interest}

\item{direction}{(character) in which direction of the geographical hierarchy
should the records be retrieved? Default: \code{upper}}

\item{type}{(character) what kind of relations should the records retrieve
have with the place? Default: \code{partof}}

\item{...}{curl options to be passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\description{
Get related records based on their MRGID.
}
\examples{
\dontrun{
## geocode to get geospatial data for a place name
(tikehau <- mr_geo_code("tikehau"))

## then pass in in an MRGID as the first parameter
mr_place_relations(tikehau$MRGID)

## Set direction='both'
mr_place_relations(tikehau$MRGID, direction = "both")

## Set type to various other options
mr_place_relations(307, type = "adjacentto")
mr_place_relations(414, type = "similarto")
mr_place_relations(4177, type = "all")
}
}
\author{
Francois Michonneau \href{mailto:francois.michonneau@gmail.com}{francois.michonneau@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_layers.R
\name{mr_layers}
\alias{mr_layers}
\title{list layers}
\usage{
mr_layers(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\description{
list layers
}
\examples{
\dontrun{
res <- mr_layers()
vapply(res, '[[', '', 'Name')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obis_eez_id.R
\name{mr_obis_eez_id}
\alias{mr_obis_eez_id}
\title{Get OBIS EEZ id}
\usage{
mr_obis_eez_id(x)
}
\arguments{
\item{x}{(character) An Exclusive Economic Zone name}
}
\value{
An integer EEZ ID if a match found in list of EEZ's, or
\code{NULL} if no match found.
}
\description{
Get OBIS EEZ id
}
\details{
internally we use the OBIS API to retrieve an EEZ id.

Matching internally is case insensitive, as we convert your input and match
against EEZ names that are all lower case.
}
\examples{
\dontrun{
# You can get EEZ names via the mr_names() function
(res <- mr_names('MarineRegions:eez_boundaries'))
mr_obis_eez_id(res$eez1[19])

# Or pass in a name
mr_obis_eez_id("Bulgarian Exclusive Economic Zone")

# case doesn't matter
mr_obis_eez_id("bulgarian exclusive economic zone")

# No match, gives NULL
mr_obis_eez_id("stuff things")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_code.R
\name{mr_geo_code}
\alias{mr_geo_code}
\title{Geocode with Marineregions}
\usage{
mr_geo_code(place, like = TRUE, fuzzy = FALSE, ...)
}
\arguments{
\item{place}{(character) a place name}

\item{like}{(logical) adds a percent-sign before and after place value
(a SQL LIKE function). Default: \code{TRUE}}

\item{fuzzy}{(logical) Uses Levenshtein query to find nearest matches.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\value{
If no results, an empty list. If results found, a data.frame with the columns:
\itemize{
\item MRGID (integer)
\item gazetteerSource (character)
\item placeType (character)
\item latitude (numeric)
\item longitude (numeric)
\item minLatitude (numeric)
\item minLongitude (numeric)
\item maxLatitude (numeric)
\item maxLongitude (numeric)
\item precision (numeric)
\item preferredGazetteerName (character)
\item preferredGazetteerNameLang (character)
\item status (character)
\item accepted (integer)
}
}
\description{
Geocode with Marineregions
}
\examples{
\dontrun{
# search for 'oost', like=TRUE, and not fuzzy
mr_geo_code(place = "oost", like = TRUE, fuzzy = FALSE)

# search for 'oost', like=FALSE, and not fuzzy
mr_geo_code(place = "oost", like = FALSE, fuzzy = FALSE)

# search for 'oost', like=FALSE, and fuzzy
mr_geo_code(place = "oost", like = FALSE, fuzzy = TRUE)

# search for 'oost', like=TRUE, and fuzzy
mr_geo_code(place = "oost", like = TRUE, fuzzy = TRUE)

# search for 'ast', like=TRUE, and fuzzy
mr_geo_code(place = "ast", like = TRUE, fuzzy = TRUE)
}
}
