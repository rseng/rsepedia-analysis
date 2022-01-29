---
title: osmdata
tags:
    - openstreetmap
    - spatial
    - R
    - Simple Features
authors:
    - name: Mark Padgham
      affiliation: 1
    - name: Robin Lovelace
      affiliation: 3
    - name: Maëlle Salmon
      affiliation: 2
    - name: Bob Rudis
      affiliation: 4
affiliations:
    - name: Department of Geoinformatics, University of Salzburg, Austria
      index: 1
    - name: ISGlobal, Centre for Research in Environmental Epidemiology,Universitat Pompeu Fabra, CIBER Epidemiología y Salud Pública, Barcelona, Spain.
      index: 2
    - name: Institute of Transport Studies, University of Leeds, U.K.
      index: 3
    - name: Rapid7
      index: 4
date: 8 March 2017
bibliography: vignettes/osmdata-refs.bib
nocite: |
  @*
---

# Summary

`osmdata` imports OpenStreetMap (OSM) data into R as either Simple Features or
`R` Spatial objects, respectively able to be processed with the R packages `sf`
and `sp`.  OSM data are extracted from the Overpass API and processed with very
fast C++ routines for return to R.  The package enables simple Overpass queries
to be constructed without the user necessarily understanding the syntax of the
Overpass query language, while retaining the ability to handle arbitrarily
complex queries. Functions are also provided to enable recursive searching
between different kinds of OSM data (for example, to find all lines which
intersect a given point). The package is faster than current alternatives for importing 
OSM data into R and is the only one compatible with `sf`.

# References
<!-- README.md is generated from README.Rmd. Please edit that file -->

# osmdata <a href='https://docs.ropensci.org/osmdata/'><img src='man/figures/logo.png' align="right" height=210 width=182/></a>

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/osmdata/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/osmdata/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/osmdata/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/osmdata)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/osmdata)](https://cran.r-project.org/package=osmdata/)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/osmdata?color=orange)](https://cran.r-project.org/package=osmdata)

<!--![](./man/figures/title.png)-->

[![](https://badges.ropensci.org/103_status.svg)](https://github.com/ropensci/software-review/issues/103)
[![status](https://joss.theoj.org/papers/10.21105/joss.00305/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00305)

<!-- badges: end -->

`osmdata` is an R package for accessing the data underlying
OpenStreetMap (OSM), delivered via the [Overpass
API](https://wiki.openstreetmap.org/wiki/Overpass_API). (Other packages
such as
[`OpenStreetMap`](https://cran.r-project.org/package=OpenStreetMap) can
be used to download raster tiles based on OSM data.)
[Overpass](https://overpass-turbo.eu) is a read-only API that extracts
custom selected parts of OSM data. Data can be returned in a variety of
formats, including as [Simple Features
(`sf`)](https://cran.r-project.org/package=sf), [Spatial
(`sp`)](https://cran.r-project.org/package=sp), or [Silicate
(`sc`)](https://github.com/hypertidy/silicate) objects. The package is
designed to allow access to small-to-medium-sized OSM datasets (see
[`osmextract`](https://github.com/ropensci/osmextract) for an approach
for reading-in bulk OSM data extracts).

## Installation

To install latest CRAN version:

``` r
install.packages("osmdata")
```

Alternatively, install the development version with any one of the
following options:

``` r
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~mpadge/osmdata")
remotes::install_bitbucket("mpadge/osmdata")
remotes::install_gitlab("mpadge/osmdata")
remotes::install_github("ropensci/osmdata")
```

To load the package and check the version:

``` r
library(osmdata)
#> Data (c) OpenStreetMap contributors, ODbL 1.0. https://www.openstreetmap.org/copyright
packageVersion("osmdata")
#> [1] '0.1.8.13'
```

## Usage

[Overpass API](https://wiki.openstreetmap.org/wiki/Overpass_API) queries
can be built from a base query constructed with `opq` followed by
`add_osm_feature`. The corresponding OSM objects are then downloaded and
converted to [Simple Feature
(`sf`)](https://cran.r-project.org/package=sf) objects with
`osmdata_sf()`, [Spatial (`sp`)](https://cran.r-project.org/package=sp)
objects with `osmdata_sp()` or [Silicate
(`sc`)](https://github.com/hypertidy/silicate) objects with
`osmdata_sc()`. For example,

``` r
x <- opq(bbox = c(-0.27, 51.47, -0.20, 51.50)) %>% # Chiswick Eyot in London, U.K.
    add_osm_feature(key = 'name', value = "Thames", value_exact = FALSE) %>%
    osmdata_sf()
x
```

    #> Object of class 'osmdata' with:
    #>                  $bbox : 51.47,-0.27,51.5,-0.2
    #>         $overpass_call : The call submitted to the overpass API
    #>                  $meta : metadata including timestamp and version numbers
    #>            $osm_points : 'sf' Simple Features Collection with 24548 points
    #>             $osm_lines : 'sf' Simple Features Collection with 2219 linestrings
    #>          $osm_polygons : 'sf' Simple Features Collection with 33 polygons
    #>        $osm_multilines : 'sf' Simple Features Collection with 6 multilinestrings
    #>     $osm_multipolygons : 'sf' Simple Features Collection with 3 multipolygons

OSM data can also be downloaded in OSM XML format with `osmdata_xml()`
and saved for use with other software.

``` r
osmdata_xml(q1, "data.osm")
```

### Bounding Boxes

All `osmdata` queries begin with a bounding box defining the area of the
query. The [`getbb()`
function](https://docs.ropensci.org/osmdata/reference/getbb.html) can be
used to extract bounding boxes for specified place names.

``` r
getbb ("astana kazakhstan")
#>        min      max
#> x 71.22444 71.78519
#> y 51.00068 51.35111
```

The next step is to convert that to an overpass query object with the
[`opq()`
function](https://docs.ropensci.org/osmdata/reference/opq.html):

``` r
q <- opq (getbb ("astana kazakhstan"))
q <- opq ("astana kazakhstan") # identical result
```

It is also possible to use bounding polygons rather than rectangular
boxes:

``` r
b <- getbb ("bangalore", format_out = "polygon")
class (b); head (b [[1]])
#> [1] "matrix" "array"
#> [1] 77.4601
```

### Features

The next step is to define features of interest using the
[`add_osm_feature()`
function](https://docs.ropensci.org/osmdata/reference/add_osm_feature.html).
This function accepts `key` and `value` parameters specifying desired
features in the [OSM key-vale
schema](https://wiki.openstreetmap.org/wiki/Map_Features). Multiple
`add_osm_feature()` calls may be combined as illustrated below, with the
result being a logical AND operation, thus returning all amenities that
are labelled both as restaurants and also as pubs:

``` r
q <- opq ("portsmouth usa") %>%
    add_osm_feature(key = "amenity", value = "restaurant") %>%
    add_osm_feature(key = "amenity", value = "pub") # There are none of these
```

Negation can also be specified by pre-pending an exclamation mark so
that the following requests all amenities that are NOT labelled as
restaurants and that are not labelled as pubs:

``` r
q <- opq ("portsmouth usa") %>%
    add_osm_feature(key = "amenity", value = "!restaurant") %>%
    add_osm_feature(key = "amenity", value = "!pub") # There are a lot of these
```

Additional arguments allow for more refined matching, such as the
following request for all pubs with “irish” in the name:

``` r
q <- opq ("washington dc") %>%
    add_osm_feature(key = "amenity", value = "pub") %>%
    add_osm_feature(key = "name", value = "irish",
                    value_exact = FALSE, match_case = FALSE)
```

Logical OR combinations can be constructed using the separate
[`add_osm_features()`
function](https://docs.ropensci.org/osmdata/reference/add_osm_features.html).
The first of the above examples requests all features that are both
restaurants AND pubs. The following query will request data on
restaurants OR pubs:

``` r
q <- opq ("portsmouth usa") %>%
    add_osm_features(features = c ("\"amenity\"=\"restaurant\"",
                                   "\"amenity\"=\"pub\""))
```

The vector of `features` contains key-value pairs separated by an
[overpass “filter”
symbol](https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#By_tag_.28has-kv.29)
such as `=`, `!=`, or `~`. Each key and value must be enclosed in
escape-delimited quotations as shown above.

Full lists of available features and corresponding tags are available in
the functions
[`?available_features`](https://docs.ropensci.org/osmdata/reference/available_features.html)
and
[`?available_tags`](https://docs.ropensci.org/osmdata/reference/available_tags.html).

### Data Formats

An overpass query constructed with the `opq()` and `add_osm_feature()`
functions is then sent to the [overpass
server](https://overpass-turbo.eu) to request data. These data may be
returned in a variety of formats, currently including:

1.  XML data (downloaded locally) via
    [`osmdata_xml()`](https://docs.ropensci.org/osmdata/reference/osmdata_xml.html);
2.  [Simple Features (sf)](https://cran.r-project.org/package=sf) format
    via
    [`osmdata_sf()`](https://docs.ropensci.org/osmdata/reference/osmdata_sf.html);
3.  [R Spatial (sp)](https://cran.r-project.org/package=sp) format via
    [`osmdata_sp()`](https://docs.ropensci.org/osmdata/reference/osmdata_sp.html);
    and
4.  [Silicate (SC)](https://github.com/hypertidy/silicate) format via
    [`osmdata_sc()`](https://docs.ropensci.org/osmdata/reference/osmdata_sc.html).

### Additional Functionality

Data may also be trimmed to within a defined polygonal shape with the
[`trim_osmdata()`](https://docs.ropensci.org/osmdata/reference/trim_osmdata.html)
function. Full package functionality is described on the
[website](https://docs.ropensci.org/osmdata/)

## Citation

``` r
citation ("osmdata")
#> 
#> To cite osmdata in publications use:
#> 
#>   Mark Padgham, Bob Rudis, Robin Lovelace, Maëlle Salmon (2017).
#>   osmdata Journal of Open Source Software, 2(14). URL
#>   https://doi.org/10.21105/joss.00305
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {osmdata},
#>     author = {Mark Padgham and Bob Rudis and Robin Lovelace and Maëlle Salmon},
#>     journal = {The Journal of Open Source Software},
#>     year = {2017},
#>     volume = {2},
#>     number = {14},
#>     month = {jun},
#>     publisher = {The Open Journal},
#>     url = {https://doi.org/10.21105/joss.00305},
#>     doi = {10.21105/joss.00305},
#>   }
```

## Data licensing

All data that you access using `osmdata` is licensed under
[OpenStreetMap’s license, the Open Database
Licence](https://wiki.osmfoundation.org/wiki/Licence). Any derived data
and products must also carry the same licence. You should make sure you
understand that licence before publishing any derived datasets.

## Other approaches

<!-- todo: add links to other packages -->

-   [osmextract](https://docs.ropensci.org/osmextract/) is an R package
    for downloading and importing compressed ‘extracts’ of OSM data
    covering large areas (e.g. all roads in a country). The package
    represents data in [`sf`](https://github.com/r-spatial/sf) format
    only, and only allows a single “layer” (such as points, lines, or
    polygons) to be read at one time. It is nevertheless recommended
    over osmdata for large queries of single layers, or where
    relationships between layers are not important.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Contributors


<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the [`allcontributors` package](https://github.com/ropenscilabs/allcontributors) following the [all-contributors](https://allcontributors.org) specification. Contributions of any kind are welcome!

### Code

<table>

<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/Robinlovelace">
<img src="https://avatars.githubusercontent.com/u/1825120?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=Robinlovelace">Robinlovelace</a>
</td>
<td align="center">
<a href="https://github.com/hrbrmstr">
<img src="https://avatars.githubusercontent.com/u/509878?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=hrbrmstr">hrbrmstr</a>
</td>
<td align="center">
<a href="https://github.com/virgesmith">
<img src="https://avatars.githubusercontent.com/u/19323577?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=virgesmith">virgesmith</a>
</td>
<td align="center">
<a href="https://github.com/maelle">
<img src="https://avatars.githubusercontent.com/u/8360597?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=maelle">maelle</a>
</td>
<td align="center">
<a href="https://github.com/agila5">
<img src="https://avatars.githubusercontent.com/u/22221146?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=agila5">agila5</a>
</td>
<td align="center">
<a href="https://github.com/espinielli">
<img src="https://avatars.githubusercontent.com/u/891692?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=espinielli">espinielli</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/idshklein">
<img src="https://avatars.githubusercontent.com/u/12258810?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=idshklein">idshklein</a>
</td>
<td align="center">
<a href="https://github.com/anthonynorth">
<img src="https://avatars.githubusercontent.com/u/391385?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=anthonynorth">anthonynorth</a>
</td>
<td align="center">
<a href="https://github.com/jeroen">
<img src="https://avatars.githubusercontent.com/u/216319?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=jeroen">jeroen</a>
</td>
<td align="center">
<a href="https://github.com/neogeomat">
<img src="https://avatars.githubusercontent.com/u/2562658?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=neogeomat">neogeomat</a>
</td>
<td align="center">
<a href="https://github.com/angela-li">
<img src="https://avatars.githubusercontent.com/u/15808896?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=angela-li">angela-li</a>
</td>
<td align="center">
<a href="https://github.com/Mashin6">
<img src="https://avatars.githubusercontent.com/u/5265707?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=Mashin6">Mashin6</a>
</td>
<td align="center">
<a href="https://github.com/odeleongt">
<img src="https://avatars.githubusercontent.com/u/1044835?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=odeleongt">odeleongt</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/Tazinho">
<img src="https://avatars.githubusercontent.com/u/11295192?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=Tazinho">Tazinho</a>
</td>
<td align="center">
<a href="https://github.com/ec-nebi">
<img src="https://avatars.githubusercontent.com/u/48711241?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=ec-nebi">ec-nebi</a>
</td>
<td align="center">
<a href="https://github.com/karpfen">
<img src="https://avatars.githubusercontent.com/u/11758039?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=karpfen">karpfen</a>
</td>
<td align="center">
<a href="https://github.com/arfon">
<img src="https://avatars.githubusercontent.com/u/4483?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=arfon">arfon</a>
</td>
<td align="center">
<a href="https://github.com/brry">
<img src="https://avatars.githubusercontent.com/u/8860095?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=brry">brry</a>
</td>
<td align="center">
<a href="https://github.com/ccamara">
<img src="https://avatars.githubusercontent.com/u/706549?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=ccamara">ccamara</a>
</td>
<td align="center">
<a href="https://github.com/danstowell">
<img src="https://avatars.githubusercontent.com/u/202965?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=danstowell">danstowell</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/dpprdan">
<img src="https://avatars.githubusercontent.com/u/1423562?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=dpprdan">dpprdan</a>
</td>
<td align="center">
<a href="https://github.com/JimShady">
<img src="https://avatars.githubusercontent.com/u/2901470?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=JimShady">JimShady</a>
</td>
<td align="center">
<a href="https://github.com/karthik">
<img src="https://avatars.githubusercontent.com/u/138494?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=karthik">karthik</a>
</td>
<td align="center">
<a href="https://github.com/MHenderson">
<img src="https://avatars.githubusercontent.com/u/23988?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=MHenderson">MHenderson</a>
</td>
<td align="center">
<a href="https://github.com/patperu">
<img src="https://avatars.githubusercontent.com/u/82020?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=patperu">patperu</a>
</td>
<td align="center">
<a href="https://github.com/stragu">
<img src="https://avatars.githubusercontent.com/u/1747497?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=stragu">stragu</a>
</td>
<td align="center">
<a href="https://github.com/fzenoni">
<img src="https://avatars.githubusercontent.com/u/6040873?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=fzenoni">fzenoni</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/rgzn">
<img src="https://avatars.githubusercontent.com/u/1675905?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/commits?author=rgzn">rgzn</a>
</td>
</tr>

</table>


### Issue Authors

<table>

<tr>
<td align="center">
<a href="https://github.com/sytpp">
<img src="https://avatars.githubusercontent.com/u/8035937?u=8efe7a4f4c3088bb35974e7488950c25658693ae&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Asytpp">sytpp</a>
</td>
<td align="center">
<a href="https://github.com/niklaas">
<img src="https://avatars.githubusercontent.com/u/705637?u=7f54fc15b926d15b2e990c93c2ab3c1bca5f271f&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Aniklaas">niklaas</a>
</td>
<td align="center">
<a href="https://github.com/RoyalTS">
<img src="https://avatars.githubusercontent.com/u/702580?u=e7d21835a6f7ba3a2f1ea7a573266708d62b1af7&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3ARoyalTS">RoyalTS</a>
</td>
<td align="center">
<a href="https://github.com/lrob">
<img src="https://avatars.githubusercontent.com/u/1830221?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Alrob">lrob</a>
</td>
<td align="center">
<a href="https://github.com/mem48">
<img src="https://avatars.githubusercontent.com/u/15819577?u=0c128db4e7567656c23e83e4314111fcea424526&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amem48">mem48</a>
</td>
<td align="center">
<a href="https://github.com/beingalink">
<img src="https://avatars.githubusercontent.com/u/871741?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Abeingalink">beingalink</a>
</td>
<td align="center">
<a href="https://github.com/yaakovfeldman">
<img src="https://avatars.githubusercontent.com/u/17687145?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Ayaakovfeldman">yaakovfeldman</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/gregor-d">
<img src="https://avatars.githubusercontent.com/u/33283245?u=3d70f9d18b0be2c20cf08a9c7d51353797d61208&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Agregor-d">gregor-d</a>
</td>
<td align="center">
<a href="https://github.com/gregmacfarlane">
<img src="https://avatars.githubusercontent.com/u/2234830?u=954f7029df0417634df181e7a27c5e163ebc8c6d&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Agregmacfarlane">gregmacfarlane</a>
</td>
<td align="center">
<a href="https://github.com/legengliu">
<img src="https://avatars.githubusercontent.com/u/7606454?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Alegengliu">legengliu</a>
</td>
<td align="center">
<a href="https://github.com/mtennekes">
<img src="https://avatars.githubusercontent.com/u/2444081?u=918a9a672af9b784480c8a6a9134ab6e3a6c3276&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amtennekes">mtennekes</a>
</td>
<td align="center">
<a href="https://github.com/lbuk">
<img src="https://avatars.githubusercontent.com/u/7860160?u=82d4376c97dbee9ec8bebd1ca0de3da8e5ddb300&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Albuk">lbuk</a>
</td>
<td align="center">
<a href="https://github.com/prokulski">
<img src="https://avatars.githubusercontent.com/u/19608488?u=cf3c1f9249688cd14fd04004efa67f4d9c67cf1e&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Aprokulski">prokulski</a>
</td>
<td align="center">
<a href="https://github.com/waholulu">
<img src="https://avatars.githubusercontent.com/u/2868000?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Awaholulu">waholulu</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/ibarraespinosa">
<img src="https://avatars.githubusercontent.com/u/27447280?u=047e026d41aa370f04fdb97b7c487eb600a012a9&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Aibarraespinosa">ibarraespinosa</a>
</td>
<td align="center">
<a href="https://github.com/tbuckl">
<img src="https://avatars.githubusercontent.com/u/98956?u=9580c2ee3c03cbbe44ac8180b0f6a6725b0415f0&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Atbuckl">tbuckl</a>
</td>
<td align="center">
<a href="https://github.com/morellek">
<img src="https://avatars.githubusercontent.com/u/38642291?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amorellek">morellek</a>
</td>
<td align="center">
<a href="https://github.com/mdsumner">
<img src="https://avatars.githubusercontent.com/u/4107631?u=c04c3e58dcca3b8c7a49ef4a3ccc6552df195e1b&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amdsumner">mdsumner</a>
</td>
<td align="center">
<a href="https://github.com/michielvandijk">
<img src="https://avatars.githubusercontent.com/u/5227806?u=956e61310e9c7ee08749ddb95458c571eafa76e3&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amichielvandijk">michielvandijk</a>
</td>
<td align="center">
<a href="https://github.com/loreabad6">
<img src="https://avatars.githubusercontent.com/u/10034237?u=53193bed2fad4f0808b55a227f99897a8d63ebc2&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Aloreabad6">loreabad6</a>
</td>
<td align="center">
<a href="https://github.com/slow-data">
<img src="https://avatars.githubusercontent.com/u/20839947?u=cd0522e56560daff7a7ed3bfedaa0ca6c85699f2&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Aslow-data">slow-data</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/mroorda">
<img src="https://avatars.githubusercontent.com/u/41475296?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amroorda">mroorda</a>
</td>
<td align="center">
<a href="https://github.com/MiKatt">
<img src="https://avatars.githubusercontent.com/u/19970683?u=1d21f231f6c2b14ce65c740014612d5e1e2ff080&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3AMiKatt">MiKatt</a>
</td>
<td align="center">
<a href="https://github.com/alanlzl">
<img src="https://avatars.githubusercontent.com/u/15748113?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Aalanlzl">alanlzl</a>
</td>
<td align="center">
<a href="https://github.com/PublicHealthDataGeek">
<img src="https://avatars.githubusercontent.com/u/43342160?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3APublicHealthDataGeek">PublicHealthDataGeek</a>
</td>
<td align="center">
<a href="https://github.com/mgageo">
<img src="https://avatars.githubusercontent.com/u/2681495?u=a98e4f2bcb64aa79f87f9e16029c8a0d3cd69768&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amgageo">mgageo</a>
</td>
<td align="center">
<a href="https://github.com/polettif">
<img src="https://avatars.githubusercontent.com/u/17431069?u=757eac2821736acbb02e7c90b456411d256d5780&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Apolettif">polettif</a>
</td>
<td align="center">
<a href="https://github.com/marcusyoung">
<img src="https://avatars.githubusercontent.com/u/10391966?u=0a04c8fedb59cd34404dabc66979bf91c4a1978c&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amarcusyoung">marcusyoung</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/barryrowlingson">
<img src="https://avatars.githubusercontent.com/u/888980?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Abarryrowlingson">barryrowlingson</a>
</td>
<td align="center">
<a href="https://github.com/ChrisWoodsSays">
<img src="https://avatars.githubusercontent.com/u/42043980?u=023bdaa73d20b313355286fec61a9f7401be0e5e&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3AChrisWoodsSays">ChrisWoodsSays</a>
</td>
<td align="center">
<a href="https://github.com/daluna1">
<img src="https://avatars.githubusercontent.com/u/60740817?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Adaluna1">daluna1</a>
</td>
<td align="center">
<a href="https://github.com/khzannat26">
<img src="https://avatars.githubusercontent.com/u/63047666?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Akhzannat26">khzannat26</a>
</td>
<td align="center">
<a href="https://github.com/gdkrmr">
<img src="https://avatars.githubusercontent.com/u/12512930?u=75e643ebcbe5e613fe9eeff8e2cf749d43ead9ea&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Agdkrmr">gdkrmr</a>
</td>
<td align="center">
<a href="https://github.com/dipenpatel235">
<img src="https://avatars.githubusercontent.com/u/8135097?u=57ce3616c4b1eb8928d0eb049d58866f7990e43c&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Adipenpatel235">dipenpatel235</a>
</td>
<td align="center">
<a href="https://github.com/robitalec">
<img src="https://avatars.githubusercontent.com/u/16324625?u=a7a98d4e17a14bf97383a5059ef4a079e15438d7&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Arobitalec">robitalec</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/nfruehADA">
<img src="https://avatars.githubusercontent.com/u/69671715?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3AnfruehADA">nfruehADA</a>
</td>
<td align="center">
<a href="https://github.com/orlandoandradeb">
<img src="https://avatars.githubusercontent.com/u/48104481?u=66d48bb0e7efb664a94eace3472aa6a06960a7f4&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Aorlandoandradeb">orlandoandradeb</a>
</td>
<td align="center">
<a href="https://github.com/changwoo-lee">
<img src="https://avatars.githubusercontent.com/u/45101999?u=2c054abd53e520d846f654b19daa3606c3e598e0&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Achangwoo-lee">changwoo-lee</a>
</td>
<td align="center">
<a href="https://github.com/maellecoursonnais">
<img src="https://avatars.githubusercontent.com/u/64737131?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Amaellecoursonnais">maellecoursonnais</a>
</td>
<td align="center">
<a href="https://github.com/Suspicis">
<img src="https://avatars.githubusercontent.com/u/78321010?u=0b4fbe51ef6fed8d90b4d4d1dabd5608f64bfc66&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3ASuspicis">Suspicis</a>
</td>
<td align="center">
<a href="https://github.com/AlbertRapp">
<img src="https://avatars.githubusercontent.com/u/65388595?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3AAlbertRapp">AlbertRapp</a>
</td>
<td align="center">
<a href="https://github.com/dmag-ir">
<img src="https://avatars.githubusercontent.com/u/89243490?u=8f64a3cd937d87a5de9d1484f25b789c960c6947&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3Admag-ir">dmag-ir</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/FlxPo">
<img src="https://avatars.githubusercontent.com/u/5145583?u=cbd02ee0a0fa0447429f38bd7e3a1da57c841239&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+author%3AFlxPo">FlxPo</a>
</td>
</tr>

</table>


### Issue Contributors

<table>

<tr>
<td align="center">
<a href="https://github.com/sckott">
<img src="https://avatars.githubusercontent.com/u/577668?u=c54eb1ce08ff22365e094559a109a12437bdca40&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Asckott">sckott</a>
</td>
<td align="center">
<a href="https://github.com/nsfinkelstein">
<img src="https://avatars.githubusercontent.com/u/2919482?u=eb162d42c4563f2cef29a6eef1d8e9e28862242d&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Ansfinkelstein">nsfinkelstein</a>
</td>
<td align="center">
<a href="https://github.com/gawbul">
<img src="https://avatars.githubusercontent.com/u/321291?u=56cf9ff94bf27ed9a5c1e6734629ce9da7969e3e&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Agawbul">gawbul</a>
</td>
<td align="center">
<a href="https://github.com/edzer">
<img src="https://avatars.githubusercontent.com/u/520851?u=9bc892c3523be428dc211f2ccbcf04e8e0e564ff&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Aedzer">edzer</a>
</td>
<td align="center">
<a href="https://github.com/MAnalytics">
<img src="https://avatars.githubusercontent.com/u/27354347?u=47f4c742c95c72b88a07ac1cb6406c9e1d186a54&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3AMAnalytics">MAnalytics</a>
</td>
<td align="center">
<a href="https://github.com/richardellison">
<img src="https://avatars.githubusercontent.com/u/10625733?u=8d7cd55a61f1a1b3f9973ddff5adbb45e0b193c6&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Arichardellison">richardellison</a>
</td>
<td align="center">
<a href="https://github.com/cboettig">
<img src="https://avatars.githubusercontent.com/u/222586?u=dfbe54d3b4d538dc2a8c276bb5545fdf4684752f&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Acboettig">cboettig</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/prise6">
<img src="https://avatars.githubusercontent.com/u/6558161?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Aprise6">prise6</a>
</td>
<td align="center">
<a href="https://github.com/PaoloFrac">
<img src="https://avatars.githubusercontent.com/u/38490683?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3APaoloFrac">PaoloFrac</a>
</td>
<td align="center">
<a href="https://github.com/Dris101">
<img src="https://avatars.githubusercontent.com/u/11404162?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3ADris101">Dris101</a>
</td>
<td align="center">
<a href="https://github.com/TomBor">
<img src="https://avatars.githubusercontent.com/u/8322713?u=bf72198850753d4eb709b2b17d89b4afa68936a1&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3ATomBor">TomBor</a>
</td>
<td align="center">
<a href="https://github.com/matkoniecz">
<img src="https://avatars.githubusercontent.com/u/899988?u=1a682cd39f51bb0224a52c7640a040c849b73ae8&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Amatkoniecz">matkoniecz</a>
</td>
<td align="center">
<a href="https://github.com/urswilke">
<img src="https://avatars.githubusercontent.com/u/13970666?u=0c6b83fb03792d052736768a8832300661c84370&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/osmdata/issues?q=is%3Aissue+commenter%3Aurswilke">urswilke</a>
</td>
</tr>

</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->


[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
0.1.9
===================

Major changes:

- New function `opq_around` to query features within a specified radius
  *around* a defined location; thanks to @barryrowlingson via #199 and
  @maellecoursonnais via #238
- New vignette on splitting large queries thanks to @Machin6 (via #262)

Minor changes:

- New dependency on `reproj` package, so that `trim_osmdata()` can be applied
  to re-projected coordinates.

0.1.8
===================

Minor changes:

- Fix some failing CRAN checks (no change to functionality)


0.1.7
===================

Minor changes:

- `add_osm_feature` bug fix to revert AND behaviour (#240 thanks to @anthonynorth)

0.1.6
===================

Major changes:

- New function `add_osm_features` to enable OR-combinations of features in
  single queries.

0.1.5
===================

Minor changes:

- Bug fix in `getbb()` via #232, thanks to @changwoo-lee
- hard-code WKT string for EPSG:4326, to avoid obsolete proj4strings (#218)
- bug fix in `print` method via #236; thanks to @odeleongt 

0.1.4
===================

Major changes:

- New `osm_enclosing()` function; thanks to @barryrowlingson via #199
- `opq()` now has additional `datetime` and `datetime2` parameters which can be
  used to extract historical data prior to `datetime`, or differences between
  two datetimes by specifying `datetime2`; thanks to @neogeomat for the idea in
  issue#179.
- opq() also has additional `nodes_only` parameter to return nodes as points
  only, for efficient extraction of strictly point-based OSM data; thanks to
  @gdkrmr for the idea in issue#221.

Minor changes:

- New contributor Enrico Spinielli (@espinielli), via #207, #210, #211, #212 - Thanks!


0.1.3
===================

Major changes:

- `osmdata_pbf` function removed as the overpass server no longer provides the
  experimental API for pbf-format data.
- Remove deprecated `add_feature()` function; entirely replaced by
  `add_osm_feature()`.
- `get_bb()` with polygon output formats now returns ALL polygon and
  multipolygon objects by default (issue#195)

Minor changes:

- New Contributors: Andrea Gilardi (@agila5)
- Bug fix for issue#205

0.1.2
===================

Major changes:

- New function `unname_osmdata_sf`, to remove row names from `sf`-format
  geometry objects that may cause issues with some plotting routines such as
  leaflet.

Minor changes:

- `getbb` now allows arbitrary `featuretype` specification, no longer just
  those pertaining to settlement forms.
- available_tags returns tags with underscore precisely as required for
  `add_osm_feature` - previous version returned text values with spaces instead
  of underscore.
- Fix bug in `osmdata_sf` for data with no names and/or no key-val pairs
- Fix bug in `trim_osmdata` for multi\* objects; thanks to @stragu
- Implement `trim_osmdata.sc` method
- retry httr calls to nominatim, which has lately been timing out quite often

0.1.1
===================

Minor changes:

- bug fix in `trim_osmdata` function

0.1.0
===================

Major changes:

- New function, `osm_elevation` to insert elevation data into `SC`-format data
  returned by `osmdata_sc` function.
- New vignette on `osmdata_sc` function and elevation data.
- `opq()` function now accepts polygonal bounding boxes generated with
  `getbb(..., format_out = "polygon")`.

0.0.10
===================

Minor changes:

- Bug fix for vectorized lists of values in `add_osm_feature`, so only listed
  items are returns (see #139; thanks @loreabad6)
- But fix to ensure all `sf` `data.frame` objects have `stringsAsFactors =
  FALSE`

0.0.9
===================

Major changes:

- New function `osmdata_sc` to return data in `silicate::SC` format (see
  github.com/hypertidy/silicate; this also requires additional dependency on
  `tibble`)
- Structure of `osmdata` object modified to replace former `$timestamp` field
  with `$meta` field containing a list of `$timestamp`, `$OSM_version`
  (currently 0.6), and `$overpass_version`.
- add_osm_feature() now accepts vectors of multiple values (see #139).
- osmdata_sf() objects default to character vectors, not factors (see #44).

Minor changes:

- vignette updated
- Overpass URL now randomly selected from the four primary servers (see
  https://wiki.openstreetmap.org/wiki/Overpass_API#Public_Overpass_API_instances),
  thanks to @JimShady.
- bug fix for osmdata_sp() (see #56)
- osmdata_sp() fixed to return osm_id values (see #131; thanks @JimShady).

0.0.8
===================
- Fix bug in `trim_osmdata` so that all sf attributes are reinstated, and also
  issue message that sf-preload is necessary for this function
- Fix bug with opq (key_exact = FALSE) so value_exact is always also set to
  FALSE

0.0.7
===================
- Fix bug in `c` method so it works when `sf` not loaded
- Fix bug in overpass query syntax to match new QL requirements

0.0.6
===================
- Add new function 'osm_poly2line()' to coerce the 'osmdata$odm_polygons' object
  for 'osmdata_sf' objects to lines, and append to 'osmdata$osm_lnes'. This is
  important for street networks ('add_osm_objects (key = "highway")'), which are
  otherwise separated between these two components. 
- Add new function `opq_osm_id` to query by OSM identifier alone
- Add `timeout` and `memsize` options to `opq()` to improve handling large
  queries.
- Return useful information from overpass server when it returns neither error
  nor useful data
- Make C++ code interruptible so long processing can be cancelled
- Fix minor yet important C++ code lines that prevented package being used as
  dependency by other packages on some systems

0.0.5
===================
- Add extraction of bounding polygons with `getbb (..., format_out = "polygon")`
- Add `trim_osmdata` function to trim an `osmdata` object to within a bounding
  polygon (thanks @sytpp)
- Add `unique_osmdata` function which reduces each component of an `osmdata`
  object only to unique elements (so `$osm_points`, for example, only contains
  points that are not represented in other - line, polygon, whatever -
  objects).
- Rename `add_feature` to `add_osm_feature` (and deprecate old version)


0.0.4
===================
- Enable alternative overpass API services through `get_overpass_url()` and
  `set_overpass_url()` functions
- Extend and improve vignette

0.0.3
===================
- Change tests only, no functional difference

0.0.2
===================
- Rename function `opq_to_string()` to `opq_string()`

0.0.1 (19 May 2017)
===================
- Remove configure and Makevars files
- Fix tests

0.0.0 (18 May 2017)
===================
- Initial CRAN release
# Contributing to osmdata

## Opening issues

The easiest way to note any behavioural curiosities or to request any new
features is by opening a [github issue](https://github.com/ropensci/osmdata/issues).


## Development guidelines

If you'd like to contribute changes to `osmdata`, we use [the GitHub
flow](https://guides.github.com/introduction/flow/index.html) for proposing,
submitting, reviewing, and accepting changes. If you haven't done this before,
there's a nice overview of git [here](http://r-pkgs.had.co.nz/git.html), as well
as best practices for submitting pull requests
[here](http://r-pkgs.had.co.nz/git.html#pr-make).

The `osmdata` coding style diverges somewhat from [this commonly used R style
guide](http://adv-r.had.co.nz/Style.html), primarily in the following two ways,
both of which improve code readability: (1) All curly braces are vertically aligned:
```r
this <- function ()
{
    x <- 1
}
```
and **not**
```r
this <- function(){
    x <- 1
}
```
and (2) Also highlighted in that code is the additional whitespace which
permeates `osmdata` code. Words of text are separated by whitespace, and so
code words should be too:
```r
this <- function1 (function2 (x))
```
and **not**
```r
this <- function1(function2(x))
```
with the natural result that one ends up writing
```r
this <- function ()
```
with a space between `function` and `()`. That's it.


## Code of Conduct

We want to encourage a warm, welcoming, and safe environment for contributing to
this project. See the [code of
conduct](https://github.com/ropensci/osmdata/blob/master/CODE_OF_CONDUCT.md) for
more information.
# CRAN notes for osmdata_0.1.9 submission

## Test environments

This submission generates NO notes on:

* Linux (via github actions): R-release, R-oldrelease
* Windows (via github actions): R-release, R-oldrelease, R-devel
* win-builder: R-oldrelease, R-release, R-devel

Package also checked using `Clang++ -Weverything and local memory sanitzer with clean results.
---
name: Bug report
about: Create a report to help us improve
title: "[BUG] <description of bug>"
labels: ''
assignees: ''

---

**Steps to follow in reporting a bug with `osmdata`**

Please follow all of the following steps, and only submit once you've checked all boxes (where appropriate).

- [ ] Download data locally via `osmdata_xml(q, filename = "myfile.xml")` (where `q` is your query)
- [ ] Use the [`reprex` package](https://reprex.tidyverse.org/) to reproduce your bug, including a `setwd()` command at the start to ensure you are in the directory where you've downloaded your data
- [ ] Include comparison with equivalent results from the `sf` package using the code shown below
- [ ] Paste the result in the indicated section below
- [ ] Delete the `setwd()` line from the pasted result
- [ ] Include the output of both `packageVersion("osmdata")` and `R.Version()$version.string`
- [ ] Alternatively, include the full output of `sessionInfo()`


### Paste `reprex` output here

Please include the following lines:
``` r
# <your reprex code here>

library(sf)
st_layers("myfile.xml") # give information on available layers
st_read("myfile.xml", layer = <desired_layer>)

packageVersion("osmdata")
R.Version()$version.string
#sessionInfo()
```



***If you are constructing or using a specialized `overpass` query***

- [ ] I have tried my query on [overpass-turbo.eu](https://overpass-turbo.eu), and it works
- [ ] I confirm that the data returned on [overpass-turbo.eu](https://overpass-turbo.eu) are identical to those returned from the equivalent call to `osmdata_xml(q, filename = "myfile.xml`)`.

Thanks! :smile:
---
name: Feature Request
about: Request a new or modified feature
title: "[FEATURE] <feature request>"
labels: ''
assignees: ''

---

**Steps to follow in requesting a feature with `osmdata`**

1. Please ensure that what you are requesting is compatible with the underlying structure of Open Street Map (OSM) itself, for example by searching the [OSM wiki](https://wiki.openstreetmap.org/wiki/Main_Page) or the overview of the [overpass Query Language](https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL) for all relevant information. Please include links to any relevant pages in your issue.
2. Please provide reproducible code generated by the [`reprex` package](https://reprex.tidyverse.org) to demonstrate current behaviour of `osmdata` which you would like changed or improved, including the following lines in your `reprex`:
``` r
# <your reprex code here>

packageVersion("osmdata")
R.Version()$version.string
#sessionInfo()
```
3. Please use the results of the preceding step to explain precisely what you would like modified from or added to the results as given.


Thanks! :smile:
