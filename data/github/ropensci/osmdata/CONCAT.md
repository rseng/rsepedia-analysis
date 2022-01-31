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
---
title: "osmdata, an R package for OpenStreetMap data"
keywords: "open street map, openstreetmap, overpass API, OSM"
output:
  rmarkdown::html_vignette:
    self_contained: no

  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r opts, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  message = TRUE,
  width = 120,
  comment = "#>",
  fig.retina = 2,
  fig.path = "README-"
)
```


# osmdata <a href='https://docs.ropensci.org/osmdata/'><img src='man/figures/logo.png' align="right" height=210 width=182/></a>



<!-- badges: start -->

[![R build
status](https://github.com/ropensci/osmdata/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/osmdata/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/osmdata/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/osmdata)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/osmdata)](https://cran.r-project.org/package=osmdata/) 
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/osmdata?color=orange)](https://cran.r-project.org/package=osmdata)

<!--![](./man/figures/title.png)-->

[![](https://badges.ropensci.org/103_status.svg)](https://github.com/ropensci/software-review/issues/103)
[![status](https://joss.theoj.org/papers/10.21105/joss.00305/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00305)

<!-- badges: end -->


`osmdata` is an R package for accessing the data underlying OpenStreetMap
(OSM), delivered via the [Overpass
API](https://wiki.openstreetmap.org/wiki/Overpass_API).  (Other packages such
as
[`OpenStreetMap`](https://cran.r-project.org/package=OpenStreetMap)
can be used to download raster tiles based on OSM data.)
[Overpass](https://overpass-turbo.eu) is a read-only API that extracts custom
selected parts of OSM data. Data can be returned in a variety of formats,
including as [Simple Features (`sf`)](https://cran.r-project.org/package=sf),
[Spatial (`sp`)](https://cran.r-project.org/package=sp), or [Silicate
(`sc`)](https://github.com/hypertidy/silicate) objects. The package is designed
to allow access to small-to-medium-sized OSM datasets (see
[`osmextract`](https://github.com/ropensci/osmextract) for an approach for
reading-in bulk OSM data extracts).


## Installation

To install latest CRAN version:
```{r cran-install, eval=FALSE}
install.packages("osmdata")
```
Alternatively, install the development version with any one of the following
options:
```{r remotes, eval = FALSE}
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~mpadge/osmdata")
remotes::install_bitbucket("mpadge/osmdata")
remotes::install_gitlab("mpadge/osmdata")
remotes::install_github("ropensci/osmdata")
```

To load the package and check the version:
```{r, eval=TRUE}
library(osmdata)
packageVersion("osmdata")
```

## Usage

[Overpass API](https://wiki.openstreetmap.org/wiki/Overpass_API) queries can be
built from a base query constructed with `opq` followed by `add_osm_feature`. The
corresponding OSM objects are then downloaded and converted to [Simple
Feature (`sf`)](https://cran.r-project.org/package=sf) objects with
`osmdata_sf()`, [Spatial (`sp`)](https://cran.r-project.org/package=sp)
objects with `osmdata_sp()` or [Silicate (`sc`)](https://github.com/hypertidy/silicate)
objects with `osmdata_sc()`.  For example, 

```{r query-thames, eval=FALSE}
x <- opq(bbox = c(-0.27, 51.47, -0.20, 51.50)) %>% # Chiswick Eyot in London, U.K.
    add_osm_feature(key = 'name', value = "Thames", value_exact = FALSE) %>%
    osmdata_sf()
x
```
```{r, echo=FALSE}
msg <- c ("Object of class 'osmdata' with:\n",
"                 $bbox : 51.47,-0.27,51.5,-0.2\n",
"        $overpass_call : The call submitted to the overpass API\n",
"                 $meta : metadata including timestamp and version numbers\n",
"           $osm_points : 'sf' Simple Features Collection with 24548 points\n",
"            $osm_lines : 'sf' Simple Features Collection with 2219 linestrings\n",
"         $osm_polygons : 'sf' Simple Features Collection with 33 polygons\n",
"       $osm_multilines : 'sf' Simple Features Collection with 6 multilinestrings\n",
"    $osm_multipolygons : 'sf' Simple Features Collection with 3 multipolygons")
message (msg)
```

OSM data can also be downloaded in OSM XML format with `osmdata_xml()` and saved
for use with other software.

```r
osmdata_xml(q1, "data.osm")
```

### Bounding Boxes

All `osmdata` queries begin with a bounding box defining the area of the query.
The [`getbb()`
function](https://docs.ropensci.org/osmdata/reference/getbb.html) can be used
to extract bounding boxes for specified place names.
```{r getbb-astana}
getbb ("astana kazakhstan")
```
The next step is to convert that to an overpass query object with the [`opq()`
function](https://docs.ropensci.org/osmdata/reference/opq.html):
```{r opq}
q <- opq (getbb ("astana kazakhstan"))
q <- opq ("astana kazakhstan") # identical result
```
It is also possible to use bounding polygons rather than rectangular boxes:
```{r getbb-haines}
b <- getbb ("bangalore", format_out = "polygon")
class (b); head (b [[1]])
```

### Features

The next step is to define features of interest using the [`add_osm_feature()`
function](https://docs.ropensci.org/osmdata/reference/add_osm_feature.html).
This function accepts `key` and `value` parameters specifying desired features
in the [OSM key-vale schema](https://wiki.openstreetmap.org/wiki/Map_Features).
Multiple `add_osm_feature()` calls may be combined as illustrated below, with
the result being a logical AND operation, thus returning all amenities that
are labelled both as restaurants and also as pubs:
```{r key-val1}
q <- opq ("portsmouth usa") %>%
    add_osm_feature(key = "amenity", value = "restaurant") %>%
    add_osm_feature(key = "amenity", value = "pub") # There are none of these
```
Negation can also be specified by pre-pending an exclamation mark so that the
following requests all amenities that are NOT labelled as restaurants and that
are not labelled as pubs:
```{r key-val2}
q <- opq ("portsmouth usa") %>%
    add_osm_feature(key = "amenity", value = "!restaurant") %>%
    add_osm_feature(key = "amenity", value = "!pub") # There are a lot of these
```

Additional arguments allow for more refined matching, such as the following
request for all pubs with "irish" in the name:
```{r key-val3}
q <- opq ("washington dc") %>%
    add_osm_feature(key = "amenity", value = "pub") %>%
    add_osm_feature(key = "name", value = "irish",
                    value_exact = FALSE, match_case = FALSE)
```

Logical OR combinations can be constructed using the separate
[`add_osm_features()`
function](https://docs.ropensci.org/osmdata/reference/add_osm_features.html).
The first of the above examples requests all features that are both restaurants
AND pubs. The following query will request data on restaurants OR pubs:

```{r features}
q <- opq ("portsmouth usa") %>%
    add_osm_features(features = c ("\"amenity\"=\"restaurant\"",
                                   "\"amenity\"=\"pub\""))
```

The vector of `features` contains key-value pairs separated by an [overpass
"filter"
symbol](https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#By_tag_.28has-kv.29)
such as `=`, `!=`, or `~`. Each key and value must be enclosed in
escape-delimited quotations as shown above.

Full lists of available features and corresponding tags are available in the functions
[`?available_features`](https://docs.ropensci.org/osmdata/reference/available_features.html)
and
[`?available_tags`](https://docs.ropensci.org/osmdata/reference/available_tags.html).


### Data Formats

An overpass query constructed with the `opq()` and `add_osm_feature()`
functions is then sent to the [overpass server](https://overpass-turbo.eu) to
request data. These data may be returned in a variety of formats, currently
including:

1. XML data (downloaded locally) via
   [`osmdata_xml()`](https://docs.ropensci.org/osmdata/reference/osmdata_xml.html);
2. [Simple Features (sf)](https://cran.r-project.org/package=sf) format via
   [`osmdata_sf()`](https://docs.ropensci.org/osmdata/reference/osmdata_sf.html);
3. [R Spatial (sp)](https://cran.r-project.org/package=sp) format via
   [`osmdata_sp()`](https://docs.ropensci.org/osmdata/reference/osmdata_sp.html);
   and
4. [Silicate (SC)](https://github.com/hypertidy/silicate) format via
   [`osmdata_sc()`](https://docs.ropensci.org/osmdata/reference/osmdata_sc.html).


### Additional Functionality {#additional}


Data may also be trimmed to within a defined polygonal shape with the
[`trim_osmdata()`](https://docs.ropensci.org/osmdata/reference/trim_osmdata.html)
function.  Full package functionality is described on the
[website](https://docs.ropensci.org/osmdata/)


## Citation

```{r}
citation ("osmdata")
```

## Data licensing

All data that you access using `osmdata` is licensed under
[OpenStreetMap's license, the Open Database Licence](https://wiki.osmfoundation.org/wiki/Licence).
Any derived data and products must also carry the same licence. You should make
sure you understand that licence before publishing any derived datasets.

## Other approaches

<!-- todo: add links to other packages -->
- [osmextract](https://docs.ropensci.org/osmextract/) is an R package for downloading and importing compressed 'extracts' of OSM data covering large areas (e.g. all roads in a country).
The package represents data in [`sf`](https://github.com/r-spatial/sf) format only, and only allows a single "layer" (such as points, lines, or polygons) to be read at one time.
It is nevertheless recommended over osmdata for large queries of single layers, or where relationships between layers are not important.


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








[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

# hex sticker script (issue 185)

## Get a background map, here of Kichijoji in Tokyo

```{r}
library (osmdata)
library (osmplotr)
library (magrittr)
bb <- getbb ("kichijoji japan")
bb [2, 1] <- bb [2, 1] - 0.2 * diff (bb [2, ])
bb [2, 2] <- mean (bb [2, ])
hw <- opq (bb) %>%
    add_osm_feature (key = "highway") %>%
    osmdata_sf (quiet = FALSE) %>%
    osm_poly2line () %>%
    extract2 ("osm_lines")
b <- opq (bb) %>%
    add_osm_feature (key = "building") %>%
    osmdata_sf (quiet = FALSE) %>%
    extract2 ("osm_polygons")
g1 <- opq (bb) %>%
    add_osm_feature (key = "landuse", value = "grass") %>%
    osmdata_sf (quiet = FALSE) %>%
    extract2 ("osm_polygons")
g2 <- opq (bb) %>%
    add_osm_feature (key = "leisure", value = "park") %>%
    osmdata_sf (quiet = FALSE) %>%
    extract2 ("osm_polygons")
g <- sf::st_sf (osm_id = c (g1$osm_id, g2$osm_id),
                 geometry = c (g1$geometry, g2$geometry))
w <- opq (bb) %>%
    add_osm_feature (key = "natural", value = "water") %>%
    osmdata_sf (quiet = FALSE) %>%
    extract2 ("osm_polygons")

#osm_structures (col_scheme = "light")
map <- osm_basemap (bbox = bb, bg = "gray95") %>%
    add_osm_objects (hw, col = "#969696FF") %>%
    add_osm_objects (b, col = "#C8C8C8FF") %>%
    add_osm_objects (g, col = "#C8FFC8FF") %>%
    add_osm_objects (w, col = "#C8C8DCFF")
#print_osm_map (map, filename = "kichijoji.png")
saveRDS (map, file = "map.Rds")
```

Reduce image to square:
```{r}
library (magrittr)
map <- readRDS ("map.Rds")
f <- "kichijoji.png"
osmplotr::print_osm_map (map, filename = f)
magick::image_read (f) %>%
    magick::image_trim () %>%
    magick::image_write (f)
i <- magick::image_read (f) %>%
    magick::image_info ()
chop <- floor (i$width - i$height) / 2
chop <- paste0 (i$height, "x", i$height, "+", chop, "+0")
magick::image_read (f) %>%
    magick::image_crop (chop) %>%
    magick::image_write (f)
```



# crop image to hex

image is 1260 h X 2100 w

```{r}
f <- "kichijoji.png"
library (magrittr)
img <- png::readPNG (f)
# define hexagon and set all outer pixels to 1
s3 <- sqrt (3) / 2
border <- data.frame (x = 1 + c (rep (-s3, 2), 0, rep (s3, 2), 0, -s3),
                      y = 1 + c (0.5, -0.5, -1, -0.5, 0.5, 1, 0.5))
border$x <- round (dim (img) [2] * border$x / 2)
border$y <- round (dim (img) [1] * border$y / 2)

#h <- 7
#w <- h * diff (range (border$x)) / diff (range (border$y))
#x11 (width = w, height = h)
#plot (border, type = "l")
border

library (sp)
library (raster)
p <- Polygon (border)
p <- SpatialPolygons (list (Polygons (list (p), "p")))

n1 <- dim (img) [1]
n2 <- dim (img) [2]
r <- raster  (nrows = n1, ncols = n2, xmn = 1, xmx = n2, ymn = 1, ymx = n1,
              vals = TRUE)
r [mask (r, p)] <- FALSE
r <- !r
#plot (r)
r <- as.matrix (r)
index <- which (!r)
index_not <- which (r)

for (i in 1:3) {
    i1 <- img [, , i]
    i1 [index] <- 0
    img [, , i] <- i1
}
mmand::display (img)
# Then add a 4th channel for alpha values
img4 <- array (dim = c (dim (img) [1:2], 4))
for (i in 1:3) {
    i1 <- img [, , i]
    img4 [, , i] <- i1
    i1 [index] <- 0
    i1 [index_not] <- 1
    img4 [, , 4] <- i1
}
png::writePNG (img4, f)
```

## make the hex logo

```{r}
library (ggplot2)
# trace outline of hexagon from centre bottom point in anti-clockwise direction
s3 <- sqrt (3) / 2
border <- data.frame (x = 1 + c (rep (-s3, 2), 0, rep (s3, 2), 0, -s3),
                      y = 1 + c (0.5, -0.5, -1, -0.5, 0.5, 1, 0.5))
asp <- diff (range (border$x)) / diff (range (border$y)) # aspect ratio for image

f <- "kichijoji.png"
d <- data.frame(x = 1, y = 1, image = f)
size <- 1.0
hex <- ggplot() +
    ggimage::geom_image (aes_ (x = ~x, y = ~y, image = ~image), d,
                         size = 1.05, asp = asp) +
    geom_polygon (aes_ (x = ~x, y = ~y), data = border,
                 size = 5, fill = NA, color = "#555555")

#extrafont::loadfonts ()
lab_dat <- data.frame (x = 1 - 0.0001,
                       y = 1 + 0.0001,
                       lab = 'osmdata')
aes <- ggplot2::aes (x, y, label = lab)
fs <- 30 # font size
hex <- hex + ggplot2::geom_text (dat = lab_dat,
                                 mapping = aes,
                                 size = fs,
                                 colour = 'gray80',
                                 family = 'SF Alien Encounters', 
                                 fontface = 1,
                                 nudge_y = 0.0001,
                                 nudge_x = 0.0001)
hex <- hex + ggplot2::geom_text (dat = lab_dat,
                                 mapping = aes,
                                 size = fs,
                                 colour = 'black',
                                 fontface = 1,
                                 family = 'SF Alien Encounters')

th <- theme_minimal ()
th$panel.background <- element_rect (fill = "transparent", size = 0)
th$line <- element_blank ()
th$axis.text <- element_blank ()
th$axis.title <- element_blank ()
th$plot.margin <- margin (rep (unit (0, 'null'), 4))
#th$plot.margin <- margin (rep (unit (-0.5, 'line'), 4))
th$legend.position <- 'none'
th$axis.ticks.length <- unit (0, 'null')

hex <- hex + th

print (hex)
```
```{r}
asp <- 1
fname <- file.path (here::here (), "man", "figures", "logo.png")
ggsave (hex, filename = fname, width = 7, height = 7 * asp)
```

it is then necessary to read the png back in and re-convert the border pixels
to alpha = 0
```{r}
fname <- file.path (here::here (), "man", "figures", "logo.png")
img <- png::readPNG (fname)
img4 <- array (1, dim = c (dim (img) [1:2], 4))
index_out <- which (img [, , 1] == 1 & img [, , 2] == 1 & img [, , 3] == 1)
index_in <- which (!seq (img [, , 1]) %in% index_out)
for (i in 1:3) {
    i1 <- img [, , i]
    img4 [, , i] <- i1
    i1 [index_in] <- 1
    i1 [index_out] <- 0
    img4 [, , 4] <- i1
}
fname <- file.path (here::here (), "man", "figures", "logo.png")
png::writePNG (img4, fname)
```
---
title: "4. Splitting large queries"
author: 
  - "Mark Padgham"
  - "Martin Machyna"
date: "`r Sys.Date()`"
bibliography: osmdata-refs.bib
output: 
    html_document:
        toc: true
        toc_float: true
        number_sections: false
        theme: flatly
vignette: >
  %\VignetteIndexEntry{4. query-split}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## 1. Introduction

The `osmdata` package retrieves data from the [`overpass`
server](https://overpass-api.de) which is primarily designed to deliver small
subsets of the full Open Street Map (OSM) data set, determined both by specific
bounding coordinates and specific OSM key-value pairs. The server has internal
routines to limit delivery rates on queries for excessively large data sets,
and may ultimately fail for large queries. This vignette describes one approach
for breaking overly large queries into a set of smaller queries, and for
re-combining the resulting data sets into a single `osmdata` object reflecting
the desired, large query.


## 2. Query splitting

Complex or data-heavy queries may exhaust the time or memory limits of the
`overpass` server. One way to get around this problem is to split the bounding
box (bbox) of a query into several smaller fragments, and then to re-combine
the data and remove duplicate objects. This section demonstrates how that may
be done, starting with a large bounding box.

```{r get-bbox, eval = FALSE}
library(osmdata)

bb <- getbb("Southeastern Connecticut COG", featuretype = "boundary")
bb
```
```{r out1, eval = FALSE}
        min       max
x -72.46677 -71.79315
y  41.27591  41.75617
```

The following lines then divide that bounding box into two smaller areas:

```{r bbox-split, eval = FALSE}
dx <- (bb["x", "max"] - bb["x", "min"]) / 2

bbs <- list(bb, bb)

bbs[[1]]["x", "max"] <- bb["x", "max"] - dx
bbs[[2]]["x", "min"] <- bb["x", "min"] + dx

bbs
```
```{r out2, eval = FALSE}
[[1]]
        min       max
x -72.46677 -72.12996
y  41.27591  41.75617

[[2]]
        min       max
x -72.12996 -71.79315
y  41.27591  41.75617
```

These two bounding boxes can then be used to submit two separate overpass
queries:

```{r opq-2x, eval = FALSE}
res <- list()

res[[1]] <- opq(bbox = bbs[[1]]) |>
        add_osm_feature(key="admin_level", value="8") |>
        osmdata_sf()
res[[2]] <- opq(bbox = bbs[[2]]) |>
        add_osm_feature(key="admin_level", value="8") |>
        osmdata_sf()
```

The retrieved `osmdata` objects can then be merged using the`c(...)` function,
which automatically removes duplicate objects.

```{r opq-merge, eval = FALSE}
res <- c(res[[1]], res[[2]])
```


## 3. Automatic bbox splitting

The previous code demonstrated how to divide a bounding box into two, smaller
regions. It will generally not be possible to know in advance how small a
bounding box should be for a query for work, and so we need a more general
version of that functionality to divide a bounding box into a arbitrary number
of sub-regions.

We can automate this process by monitoring the exit status of `opq() |>
osmdata_sf()` and in case of a failed query we can keep recursively splitting
the current bounding box into increasingly smaller fragments until the overpass
server returns a result. The following function demonstrates splitting a
bounding box into a list of four equal-sized bounding boxes in a 2-by-2 grid.

```{r bbox-auto-split, eval = FALSE}
split_bbox <- function(bbox, grid = 2) {
    xmin <- bbox["x", "min"]
    ymin <- bbox["y", "min"]
    dx <- (bbox["x", "max"] - bbox["x", "min"]) / grid
    dy <- (bbox["y", "max"] - bbox["y", "min"]) / grid

    bboxl <- list()

    for (i in 1:grid) {
        for (j in 1:grid) {
            b <- matrix(c(xmin + ((i-1) * dx),
                          ymin + ((j-1) * dy),
                          xmin + (i * dx),
                          ymin + (j * dy)),
                        nrow = 2,
                        dimnames = dimnames(bbox))

            bboxl <- append(bboxl, list(b))
        }
    }
    bboxl
}
```

We pre-split our area and create a queue of bounding boxes that we will use for 
submitting queries.

```{r bbox-pre-split, eval = FALSE}
bb <- getbb("Connecticut", featuretype = NULL)  
queue <- split_bbox(bb)
result <- list()
```

Now we can create a loop that will monitor the exit status of our query and in 
case of success remove the bounding box from the queue. If our query fails for 
some reason, we split the failed bounding box into four smaller fragments and
add them to our queue, repeating until all results have been successfully
delivered.

```{r auto-query, eval = FALSE}
while (length(queue) > 0) {

    print(queue[[1]])

    opres <- NULL
    opres <- try({
                opq(bbox = queue[[1]], timeout = 25) |>
                    add_osm_feature(key="natural", value="tree") |>
                    osmdata_sf()
              })

    if (class(opres)[1] != "try-error") {
        result <- append(result, list(opres))
        queue <- queue[-1]
    } else {
        bboxnew <- split_bbox(queue[[1]])
        queue <- append(bboxnew, queue[-1])
    }
}
```

All retrieved `osmdata` objects stored in the `result` list can then be
combined using the `c(...)` operator. Note that for large datasets this process
can be quite time consuming.

```{r merge-result-list, eval = FALSE}
final <- do.call(c, result)
```
---
title: "3. Translation of OSM to Simple Features"
author: "Mark Padgham"
date: "`r Sys.Date()`"
output:
    html_document:
        toc: true
        toc_float: true
        theme: flatly
vignette: >
  %\VignetteIndexEntry{3. OSM to Simple Features}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


## 1. OpenStreetMap Data Structure

[OpenStreetMap (OSM)](https://www.openstreetmap.org/) data has a unique structure that
is not directly reconcilable with other modes of representing spatial data,
notably including the widely-adopted 
[Simple Features (SF)](https://www.ogc.org/standards/sfa) 
scheme of the 
[Open Geospatial Consortium (OGC)](https://www.ogc.org). 
The three primary spatial objects of OSM are:

1. `nodes`, which are directly translatable to spatial points

2. `ways`, which may be closed, in which case they form polygons, or unclosed,
   in which case they are (non-polygonal) lines.

3. `relations` which are higher-level objects used to specify relationships
   between collections of ways and nodes. While there are several recognised
   categories of `relations`, in spatial terms these may be reduced to a binary
   distinction between:

   - `multipolygon` relations, which specify relationships between an
   exterior polygon (through designating `role='outer'`) and possible inner
   polygons (`role='inner'`). These may or may not be designated with
   `type=multipolygon`. Political boundaries, for example, often have
   `type=boundary` rather than explicit `type-multipolygon`. `osmdata`
   identifies multipolygons as those `relation` objects having at least one
   member with `role=outer` or `role=inner`.

   - In the absence of `inner` and `outer` roles, an OSM relation is assumed to
   be non-polygonal, and to instead form a collection of non-enclosing lines.

--------

## 2. Simple Features Data Structure

The representation of spatial objects as Simple Features is 
[described at length by the OGC](https://www.ogc.org/standards/sfo),
with this document merely reviewing relevant aspects. The SF system assumes that
spatial features can be represented in one of seven distinct primary classes,
which by convention are referred to in all capital letters.  Relevant classes
for OSM data are:

1. POINT
2. MULTIPOINT
3. LINESTRING
4. MULTILINESTRING
5. POLYGON
6. MULTIPOLYGON

(The seventh primary class is `GEOMETRYCOLLECTION`, which contains several
objects with different geometries.) An SF (where that acronym may connote both
singular and plural) consists of a sequence of spatial coordinates, which for
OSM data are only ever `XY` coordinates represented as strings enclosed within
brackets.  In addition to coordinate data and associated coordinate reference
systems, an SF may include any number of additional data which quantify or
qualify the feature of interest.  In the 
[`sf` extension to `R`](https://cran.r-project.org/package=sf), for example, a
single SF is represented by one row of a `data.frame`, with the geometry stored
in a single column, and any number of other columns containing these additional
data.

Simple Feature geometries are referred to in this vignette using all capital
letters (such as `POLYGON`), while OSM geometries use lower case (such as
`polygon`). Similarly, the Simple Features standard of the OGC is referred to as
`SF`, while the `R` package of the same name is referred to as `R::sf`--upper
case `R` followed by lower case `sf`. Much functionality of `R::sf` is
determined by the underlying 
[Geospatial Data Abstraction Library](https://gdal.org) (`GDAL`; described
below). Representations of data are often discussed here with reference to
`GDAL/sf`, in which case it may always be assumed that the translation and
representation of data are determined by `GDAL` and not directly by the creators
of `R::sf`.


--------

## 3. How `osmdata` translates OSM into Simple Features

### 3.1. OSM Nodes

OSM nodes translate directly into `SF::POINT` objects, with all OSM `key-value`
pairs stored in additional `data.frame` columns.

### 3.2. OSM Ways

OSM ways may be either polygons or (non-polygonal) lines. `osmdata` translates
these into `SF::LINESTRING` and `SF::POLYGON` objects, respectively. Although
polygonal and non-polygonal ways may have systematically different `key` fields,
they are conflated here to the single set of `key` values common to all `way`
objects regardless of shape. This enables direct comparison and uniform
operation on both `SF::LINESTRING` and `SF::POLYGON` objects.

### 3.3 OSM Relations

OSM relations comprising members with `role=outer` or `role=inner` are
translated into `SF::MULTIPOLYGON` objects; otherwise they form
`SF::MULTILINESTRING` objects. As in the preceding case of OSM ways, potentially
systematic differences between OSM `key` fields for `multipolygon` and other
`relation` objects are ignored in favour of returning identical `key` fields in
both cases, whether or not `value` fields for those `key`s exist.

#### 3.3(a) Multipolygon Relations

An OSM multipolygon is translated by `osmdata` into a single `SF::MULTIPOLYGON`
object which has an additional column specifying `num_members`. The `SF`
geometry thus consists of a list (an `R::List` object) of this number of
polygons, the first of which is the `outer` polygon, with all subsequent members
forming closed inner rings (either individually or in combination).

Each of these inner polygons are also represented as one or more OSM objects, which
will generally include detailed data on the individual components **not** able
to be represented in the single multipolygon representation.  Each inner polygon
is therefore additionally stored in the `sf::MULTIPOLYGON` `data.frame` along
with all associated data.  Thus the row containing a multipolygon of
`num_polygon` polygons is followed by `num_polygon - 1` rows containing the data
for each `inner` polygon.

Note that OSM `relation` objects generally have fewer (or different) `key-value`
pairs than do OSM `way` objects. In the OSM system, data describing the detailed
properties of the constituent ways of a given OSM relation are stored with those
`ways` rather than with the `relation`. `osmdata` follows this general
principle, and stored the geometry of all `ways` of a `relation` with the
relation itself (that is, as part of the `MULTIPOLYGON` or `MULTILINESTRING`
object), while those ways are also stored themselves as `LINESTRING` (or
potentially `POLYGON`) objects, from where their additional `key-value` data may
be accessed.


#### 3.3(b) Multilinestring Relations

OSM `relations` that are not `multipolygons` are translated into
`SF::MULTILINESTRING` objects. Each member of any OSM relation is attributed a
`role`, which may be empty. `osmdata` collates all ways within a relation
according to their `role` attributes. Thus, unlike multipolygon relations which
are always translated into a single `sf::MULTIPOLYGON` object, multilinestring
relations are translated by `omsdata` into potentially several
`sf::MULTILINESTRING` objects, one for each unique role.

This is particularly useful because `relations` are often used to designated
extended `highways` (for example, designated bicycle routes or motorways), yet
these often exist in `primary` and `alternative` forms, with these categories
specified in roles. Separating these roles enables ready access to any desired
role.

These multilinestring objects also have a column specifying `num_members`, as
for multipolygons, with the primary member followed by `num_members` rows, one
for each member of the multilinestring.

---------

## 4. `GDAL` Translation of OSM into Simple Features

The `R` package [`sf`](https://cran.r-project.org/package=sf) provide an `R`
implementation of Spatial Features, and provides a wrapper around GDAL for
reading geospatial data. `GDAL` provides a ['driver' to read OSM
data](https://gdal.org/drv_osm.html), and thus `sf` can also be used to read
`OSM` data in `R`, 
[as detailed in the main `osmdata` vignette](https://docs.ropensci.org/osmdata/articles/osmdata.html).
However, the `GDAL` translation of OSM data differs in several important ways
from the `osmdata` translation.

The primary difference is that GDAL only returns *unique* objects of each
spatial (SF) type. Thus `sf::POINT` objects consist of only those points that
are not otherwise members of some 'higher' object (`line`, `polygon`, or
`relation` objects).  Although a given set of OSM data may actually contain a
great many points, attempting to load these with
```{r, eval=FALSE}
sf::st_read (file, layer = 'points')
```
will generally return surprisingly few points.

### 4.1. OSM Nodes

Apart from the numerical difference arising through `osmdata` returning an
`sf::POINTS` structure containing **all** nodes within a given set of OSM data,
while `sf::st_read (file, layer='points')` returns only those points not
represented in other structure, the representation of points remains otherwise
broadly similar. The only other major difference is that `osmdata` retains all
`key-value` pairs present in a given set of OSM data, whereas `GDAL/sf` only
retains a select few of these. Moreover, the `keys` returned by `GDAL/sf` are
pre-defined and invariant, meaning that data returned from `sf::st_read (...)`
may often contain `key` columns in the resultant `data.frame` which contain no
(non-`NA`) data. This difference is illustrated in an example repeated here from
the  
[main `osmdata` vignette](https://docs.ropensci.org/osmdata/articles/osmdata.html),
with the same principles applying to all of the following classes of OSM data.

The following three lines define a query and download the resultant data to an
`XML` file.
```{r trentham, eval=FALSE}
q <- opq (bbox = 'Trentham, Australia')
q <- add_osm_feature (q, key = 'name') # any named objects
osmdata_xml (q, 'trentham.osm')
```
These data may then be converted into SF representations using either `R::sf` or
`osmdata`, with OSM `keys` being the column names of the resultant `data.frame`
objects.
```{r, eval=FALSE}
names (sf::st_read ('trentham.osm', layer = 'points', quiet = TRUE))
```
```{r, echo=FALSE}
c ("osm_id",     "name",       "barrier",    "highway",    "ref",
"address",    "is_in",      "place",      "man_made",   "other_tags",
"geometry")
```
```{r, eval=FALSE}
names (osmdata_sf (q, 'trentham.osm')$osm_points)
```
```{r, echo=FALSE}
c ("osm_id",           "name",             "X_description_",   "X_waypoint_",
"addr.city",        "addr.housenumber", "addr.postcode",    "addr.street",
"amenity",          "barrier",          "denomination",     "foot",
"ford",             "highway",          "leisure",          "note_1",
"phone",            "place",            "railway",  "railway.historic",
"ref",              "religion",         "shop",             "source",
"tourism",          "waterway",         "geometry")
```
`osmdata` returns far more `key` fields than does `GDAL/sf`. More importantly,
however, `GDAL/sf` returns pre-defined `key` fields regardless of whether they
contain any data:
```{r, eval=FALSE}
addr <- sf::st_read ('trentham.osm', layer = 'points', quiet = TRUE)$address
all (is.na (addr))
```
```{r, echo=FALSE}
TRUE
```
In contrast, `osmdata` returns only those `key` fields which contain data (and
so excludes `address` in the above example).


### 4.2. OSM Ways

As for points, `GDAL/sf` only returns those ways that are not represented or
contained in 'higher' objects (OSM relations interpreted as `SF::MULTIPOLYGON`
or `SF::MULTILINESTRING` objects). `osmdata` returns all ways, and thus enables,
for example, examination of the full attributes of any member of a multigeometry
object. This is not possible with the `GDAL/sf` translation.  As for points, the
only additional difference between `osmdata` and `GDAL/sf` is that `osmdata`
retains all `key-value` pairs, whereas `GDAL` retains only a select few.

### 4.3 OSM Relations

Translation of OSM relations into Simple Features differs more significantly
between `osmdata` and `GDAL/sf`.

#### 4.3(a) Multipolygon Relations

As indicated above, multipolygon relations are translated in broadly comparable
ways by both `osmdata` and `sf/GDAL`. Note, however, the `way` members of an OSM
relation may be specified in arbitrary order, and the multipolygonal way may not
necessarily be traced through simply following the segments in the order
returned by `sf/GDAL`.

#### 4.3(b) Multilinestring Relations

Linestring relations are simply read by GDAL directly in terms of the their
constituent ways, resulting in a single `SF::MULTILINESTRING` object that
contains exactly the same number of lines as the ways in the OSM relation,
regardless of their `role` attributes. Note that `roles` are frequently used to
specify `alternative` multi-way routes through a single OSM relation. Such
distinctions between primary and alternative are erased with `GDAL/sf` reading.

## 5 Examples

### 5.1 Routing

Navigable paths, routes, and ways are all tagged within OSM as `highway`,
readily enabling an `overpass` query to return only `ways` that can be used for
routing purposes. Routes are nevertheless commonly assembled within OSM
relations, particularly where they form major, designated transport ways such as
long-distance foot or bicycle paths or major motorways.

#### 5.1(a) Routing with `sf/GDAL`

A query for `key=highway` translated through `GDAL/sf` will return those ways
not part of any 'higher' structure as `SF::LINESTRING` objects, but components
of an entire transport network might also be returned as:

1. `SF::MULTIPOLYGON` objects, holding all single ways which form simple
   polygons (that is, in which start and end points are the same); 
2. `SF::MULTIPOLYGON` objects holding all single (non-polygonal) ways which
   combine to form an `OSM multipolygon` relation (that is, in which the
   collection of ways ultimately forms a closed `role=outer` polygon).
3. `SF::MULTILINESTRING` objects holding all single (non-polygonal) ways which
   combine to form an OSM relation that is not a multipolygon.

Translating these data into a single form usable for routing purposes is not
simple. A particular problem that is extremely difficult to resolve is
reconciling the `SF::MULTIPOLYGON` objects with the geometry of the
`SF::LINESTRING` objects. Highway components contained in `SF::MULTIPOLYGON`
objects need to be re-connected with the network represented by the
`SF::LINESTRING` objects, yet the OSM identifiers of the `MULTIPOLYGON`
components are removed by `sf/GDAL`, preventing these components from being
directly re-connected. The only way to ensure connection would be to re-connect
those geographic points sharing identical coordinates. This would require code
too long and complicated to be worthwhile demonstrating here.

#### 5.1(b) Routing with `osmdata`

`osmdata` retains all of the underlying ways of 'higher' structures
(`SF::MULTIPOLYGON` or `SF::MULTILINESTRING` objects) as 
`SF::LINESTRING` or `SF::POLYGON` objects. The geometries of the latter objects
duplicate those of the 'higher' relations, yet contain additional `key-value`
pairs corresponding to each way. Most importantly, the OSM ID values for all
members of a `relation` are stored within that relation, readily enabling the
individual ways (`LINESTRING` or `POLYGON` objects) to be identified from the
`relation` (`MULTIPOLYGON` or `MULTILINESTRING` object).

The `osmdata` translation thus readily enables a singularly complete network to
be reconstructed by simply combining the `SF::LINESTRING` layer with the
`SF::POLYGON` layer. These layers will always contain entirely independent
members, and so will always be able to be directly combined without duplicating
any objects.
---
title: "2. Elevation data and OSM: The osmdata_sc function"
author: 
  - "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        number_sections: false
        theme: flatly
vignette: >
  %\VignetteIndexEntry{2. osmdata_sc}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## 1. Introduction


[`silicate`](https://github.com/hypertidy/silicate) is a new
form for representing spatial data in **R**. In contrast to all
other forms (such as [`sp`](https://cran.r-project.org/package=sp) or
[`sf`](https://cran.r-project.org/package=sf)),
[`silicate`](https://github.com/hypertidy/silicate) is *multi-tabular*, and
primarily consists of one table for point entities; one table for binary
relationships between these point entities -- spatial ''edges'' -- and
additional tables for higher-order inter-relationships between these objects.
The new `osmdata` function `osmdata_sc()` returns Open Street Map (OSM) data in
[`silicate`](https://github.com/hypertidy/silicate) form. This form also closely
resembles the data storage scheme of Open Street Map itself, and in this case
consists of the following tables:

1. A `vertex` table holding the coordinates of all OSM nodes;
2. An `edge` table mapping all edge connections between vertices;
3. An `object_link_edge` table linking all `edge` entities to the OSM objects of
   which they are part;
4. An `object` table holding all 'key'--'value' pairs for each OSM *way* object;
5. A `relation_properties` table holding all 'key'--'value' pairs for each OSM
   *relation* object;
6. A `relation_members` table holding all members of each OSM relation; and
7. A `nodes` table holding all 'key'--'value' pairs for each OSM *node* object.

The translation of the underlying OSM data structure -- consisting of nodes,
way, and relations -- into 
[Simple Features (SF)](https://www.ogc.org/standards/sfa) via the
[`osmdata_sf()`](https://docs.ropensci.org/osmdata/reference/osmdata_sf.html)
function is less than 100% faithful, and results in some representational loss
compared with the original OSM structure (for details, see the [vignette on
translation of OSM into
SF](https://docs.ropensci.org/osmdata/articles/osm-sf-translation.html)). In
contrast, the [`osmdata_sc()`
function](https://docs.ropensci.org/osmdata/reference/osmdata_sc.html) delivers
a representation that is entirely faithful to the underlying OSM representation.

One of the advantages of [`silicate`](https://github.com/hypertidy/silicate)
format offered by the `osmdata` package is enabling elevation data to be
combined with OSM data. The result is
a [`silicate`](https://github.com/hypertidy/silicate)-format object which is
able to be submitted directly to the [`dodgr`
package](https://github.com/ATFutures/dodgr) to enable routing on street
networks that accounts for elevation changes.


## 2. Elevation Data

Incorporating elevation data with OSM data currently requires local storage of
desired elevation data. These must be downloaded for the desired region from 
[http://srtm.csi.cgiar.org/srtmdata](https://srtm.csi.cgiar.org/srtmdata/) in Geo
TIFF format. Elevation data may then be incorporated with 
[`silicate`](https://github.com/hypertidy/silicate)-format data generated by 
[`x <- osmdata_sc()`](https://docs.ropensci.org/osmdata/reference/osmdata_sc.html) 
through the `osm_elevation()` function. The entire procedure is demonstrated
with the following lines:
```{r omaha, eval = FALSE}
dat <- opq ("omaha nebraska") %>%
    add_osm_feature (key = "highway") %>%
    osmdata_sc ()
```
This object has a `vertex` table like this:
```{r dat_vertex, eval = FALSE}
dat$vertex
```
```{r dat_vertex_dat, echo = FALSE}
n <- 345239
x_ <- c (-95.9, -95.9, -95.9, -95.9, -95.9, -95.9, -96.2, -96.2, -96.3, -96.3)
y_ <- c (41.2, 41.2, 41.2, 41.2, 41.2, 41.2, 41.3, 41.3, 41.3, 41.3)
z_ <- c (291.0, 295.0, 297.0, 301.0, 295.0, 300.0, 359.0, 359.0, 358.0, 358.0)
vertex_ <- paste0 (c (31536366, 31536367, 31536368, 31536370, 31536378, 31536379,
                      133898322, 133898328, 133898340, 133898342))

tibble::tibble (x_ = c (x_, rep (NA, n - 10)),
                y_ = c (y_, rep (NA, n - 10)),
                vertex_ = c (vertex_, rep (NA, n - 10)))
```

Incorporating elevation data is then as simple as
```{r osm_elevation, eval = FALSE}
dat <- osm_elevation (dat, elev_file = "/path/to/elevation/data/filename.tiff")
```
```{r osm_elevation2, echo = FALSE}
message ("Loading required namespace: raster\n",
         "Elevation data from Consortium for Spatial Information; see ",
         "https://srtm.csi.cgiar.org/srtmdata/")
```
This function then simply appends the elevation values to the `vertex_` table,
so that it now looks like this:
```{r dat_vertex2, eval = FALSE}
dat$vertex_
```
```{r, dat_vertex_dat2, echo = FALSE}
tibble::tibble (x_ = c (x_, rep (NA, n - 10)),
                y_ = c (y_, rep (NA, n - 10)),
                z_ = c (z_, rep (NA, n - 10)),
                vertex_ = c (vertex_, rep (NA, n - 10)))
```

### Example usage of elevation data

The [`silicate`](https://github.com/hypertidy/silicate) format is very easy to
manipulate using standard [`dplyr`](https://dplyr.tidyverse.org) verbs. The
following code uses the [`mapdeck`
package](https://github.com/SymbolixAU/mapdeck) package to colour the street
network and elevation data downloaded and processed in the preceding lines by
the elevation of each network edge. We first join the vertex elevation data on
to the edges, and calculate the mean elevation of each edge.

```{r edges, eval = FALSE}
edges <- dplyr::left_join (dat$edge, dat$vertex, by = c (".vx0" = "vertex_")) %>%
    dplyr::rename (".vx0_x" = x_, ".vx0_y" = y_, ".vx0_z" = z_) %>%
    dplyr::left_join (dat$vertex, by = c (".vx1" = "vertex_")) %>%
    dplyr::rename (".vx1_x" = x_, ".vx1_y" = y_, ".vx1_z" = z_) %>%
    dplyr::mutate ("zmn" = (.vx0_z + .vx1_z) / 2) %>%
    dplyr::select (-c (.vx0_z, .vx1_z))
edges
```
```{r edges-dat, echo = FALSE}
n <- 376370
x <- paste0 (c (1903265686, 1903265664, 1903265638, 1903265710, 1903265636,
                 1903265685, 1903265678, 1903265646, 1903265714, 1903265659))
y <- paste0 (c (1903265664, 1903265638, 1903265710, 1903265636, 1903265685,
                1903265678, 1903265646, 1903265714, 1903265659, 1903265702))
edge <- c ("V6kgqvWjtM", "mX4HQkykiD", "26e5NHT8nI", "9TOmVAvGH4", "hYbpf832vX",
           "ctvd1FWGEw", "mvaAOdSOKA", "dSVFPNDFty", "uc8L3jGR87", "MpjXnvIvcF")
x0_x <- c (-96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2)
x0_y <- c (41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3)
x1_x <- c (-96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2, -96.2)
x1_y <- c (41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3, 41.3)
z <- c (351.0, 352.0, 352.0, 352.0, 352.0, 352.0, 352.0, 352.0, 352.0, 352.0)

tibble::tibble (".vx0" = c (x, rep (NA, n - 10)),
                ".vx1" = c (y, rep (NA, n - 10)),
                "edge_" = c (edge, rep (NA, n - 10)),
                ".vx0_x" = c (x0_x, rep (NA, n - 10)),
                ".vx0_y" = c (x0_y, rep (NA, n - 10)),
                ".vx1_x" = c (x1_x, rep (NA, n - 10)),
                ".vx1_y" = c (x1_y, rep (NA, n - 10)),
                "zmn" = c (z, rep (NA, n - 10)))
```

Those data can then be submitted directly to
[`mapdeck`](https://github.com/SymbolixAU/mapdeck) to generate an interactive
plot with the following code:
```{r, eval = FALSE}
library (mapdeck)
set_token (Sys.getenv ("MAPBOX_TOKEN")) # load local token for MapBox
mapdeck (style = mapdeck_style ("dark")) %>%
    add_line (edges,
              origin = c (".vx0_x", ".vx0_y"),
              destination = c (".vx1_x", ".vx1_y"),
              stroke_colour = "z",
              legend = TRUE)
```
(The result is not shown here, but can be directly inspected by simply running
the above lines.)
---
title: "1. osmdata"
author: 
  - "Mark Padgham"
  - "Robin Lovelace"
date: "`r Sys.Date()`"
bibliography: osmdata-refs.bib
output: 
    html_document:
        toc: true
        toc_float: true
        number_sections: false
        theme: flatly
vignette: >
  %\VignetteIndexEntry{1. osmdata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## 1. Introduction

`osmdata` is an R package for downloading and using data from OpenStreetMap
([OSM](https://www.openstreetmap.org/)).  OSM is a global open access mapping project,
which is free and open under the 
[ODbL licence](https://www.openstreetmap.org/copyright) [@OpenStreetMap].  This
has many benefits, ensuring transparent data provenance and ownership, enabling
real-time evolution of the database and, by allowing anyone to contribute,
encouraging democratic decision making and citizen science
[@johnson_models_2017].  See the 
[OSM wiki](https://wiki.openstreetmap.org/wiki/Contribute_map_data) to find out
how to contribute to the world's open geographical data commons. 

Unlike the [`OpenStreetMap`](https://cran.r-project.org/package=OpenStreetMap)
package, which facilitates the download of raster tiles, `osmdata` provides
access to the vector data underlying OSM.

`osmdata` can be installed from CRAN with
```{r, eval = FALSE}
install.packages("osmdata")
```
and then loaded in the usual way:
```{r}
library(osmdata)
```
The development version of `osmdata` can be installed with the `remotes`
package using the following command:

```{r install, eval = FALSE}
remotes::install_github('ropensci/osmdata')
```

`osmdata` uses the [`overpass` API](https://overpass-api.de) to download

OpenStreetMap (OSM) data and can convert the results to a variety of formats,
including both Simple Features (typically of class `sf`) and Spatial objects
(e.g. `SpatialPointsDataFrame`), as defined by the packages
[`sf`](https://cran.r-project.org/package=sf) and
[`sp`](https://cran.r-project.org/package=sp) packages respectively.

`overpass` is a C++ library that serves OSM data over the web.  All `overpass`
queries begin with a bounding box, defined in `osmdata` with the function
`opq()`:

```{r opq1}
q <- opq(bbox = c(51.1, 0.1, 51.2, 0.2))
```

The following sub-section provides more detail on bounding boxes.  Following the
initial `opq()` call, `osmdata` queries are built by adding one or more
'features', which are specified in terms of `key-value` pairs.  For example, all
paths, ways, and roads are designated in OSM with `key=highway`, so that a query
all motorways in greater London (UK) can be constructed as follows:

```{r opq3, eval = FALSE}
q <- opq(bbox = 'greater london uk') %>%
    add_osm_feature(key = 'highway', value = 'motorway')
```
```{r, echo = FALSE}
q <- opq (bbox = c (51.2867602, -0.510375, 51.6918741, 0.3340155)) %>%
    add_osm_feature(key = 'highway', value = 'motorway')
```

A detailed description of features is provided at the 
[OSM wiki](https://wiki.openstreetmap.org/wiki/Map_Features), or the 
`osmdata` function `available_features()` can be used to retrieve the
comprehensive list of feature keys currently used in OSM.

```{r available-features, eval=FALSE}
head (available_features ())
```

```{r available-features-results, echo=FALSE}
c ("4wd only", "abandoned", "abutters", "access", "addr", "addr:city")
```

There are two primary `osmdata` functions for obtaining data from a query:
`osmdata_sf()` and `osmdata_sp()`, which return data in 
[Simple Features (`sf`)](https://cran.r-project.org/package=sf)
and [Spatial (`sp`)](https://cran.r-project.org/package=sp) formats,
respectively. The typical workflow for extracting OSM data with `osmdata` thus
consists of the three lines:
```{r workflow, eval = FALSE}
x <- opq(bbox = 'greater london uk') %>%
    add_osm_feature(key = 'highway', value = 'motorway') %>%
    osmdata_sf ()
```
The return object (`x`) is described in the third section below.


### 1.1 Bounding boxes: the `getbb()` function

While bounding boxes may be explicitly specified for the `opq()` function, they
are more commonly obtained from the `getbb()` function, which accepts character
strings. As illustrated in the above example, the `opq()` function also accepts
character strings, which are simply passed directly to `getbb()` to convert them
to rectangular bounding boxes.
```{r opq2, eval = FALSE}
bb <- getbb('Greater London, U.K.')
q <- opq(bbox = bb)
```

Note that the text string is not case sensitive, as illustrated in the following
code:
```{r eval=FALSE}
identical(q, opq(bbox = 'greater london uk'))
## TRUE
```

Note also that `getbb()` can return a data frame reporting multiple matches or
matrices representing bounding polygons of matches:
```{r, eval=FALSE}
bb_df <- getbb(place_name = "london", format_out = "data.frame")
bb_poly <- getbb(place_name = "london", format_out = "polygon")
```
The [`overpass API`](https://www.overpass-api.de) only accepts simple rectangular
bounding boxes, and so data requested with a bounding polygon will actually be
all data within the corresponding rectangular bounding box, but such data may be
subsequently trimmed to within the polygon with the `trim_osmdata()` function,
demonstrated in the code immediately below.

All highways from within the polygonal boundary of Greater London can be
extracted with,
```{r trim-osmdata, eval = FALSE}
bb <- getbb ('london uk', format_out = 'polygon')
x <- opq(bbox = bb) %>%
    add_osm_feature(key = 'highway', value = 'motorway') %>%
    osmdata_sf () %>%
    trim_osmdata (bb)
```
See `?trim_osmdata()` for further ways to obtain polygonally bounded sets of OSM
data.

The `getbb()` function also allows specification of an explicit `featuretype`,
such as street, city, county, state, or country. The default value of
`settlement` combines all results below country and above streets. See `?getbb`
for more details.


## 2. The overpass API

As mentioned, `osmdata` obtains OSM data from the 
[`overpass API`](https://www.overpass-api.de), 
[which is](https://wiki.openstreetmap.org/wiki/Overpass_API)

> a read-only API that serves up custom selected parts of the OSM map data.

The syntax of `overpass` queries is powerful yet hard to learn.  This
section briefly introduces the structure of `overpass` queries in order to help
construct more efficient and powerful queries.  Those wanting to skip straight
onto query construction in `osmdata` may safely jump ahead to the [query example
below](#query-example).

`osmdata` simplifies queries so that OSM data can be extracted with very little
understanding of the `overpass` query syntax, although it is still possible to
submit arbitrarily complex `overpass`  queries via `osmdata`.  An excellent
place to explore `overpass` queries specifically and OSM data in general is the 
online interactive query builder at [overpass-turbo](https://overpass-turbo.eu/),
which includes a helpful corrector function for incorrectly formatted queries.
Examples of its functionality in action can be found on the 
[OpenStreetMap wiki](https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_API_by_Example),
with full details of the `overpass`
query language given in the 
[Query Language Guide](https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL) 
as well as the 
[overpass API Language Guide](https://wiki.openstreetmap.org/wiki/Overpass_API/Language_Guide).

By default, `osmdata` sends queries to one of the four main [`overpass` server
instances](https://wiki.openstreetmap.org/wiki/Overpass_API#Public_Overpass_API_instances),
such as `https://overpass-api.de/api/interpreter` but other servers listed on the
page linked to above can be used, thanks to functions that *get* and *set* the
base url:

```{r}
get_overpass_url()
new_url <- "https://overpass.openstreetmap.ie/api/interpreter"
```

```{r, eval=FALSE}
set_overpass_url(new_url) # reset the base url (not run)
```

`osmdata` queries are lists of class `overpass_query`. The actual query passed
to the `overpass API` with a query can be obtained with the function
`opq_string()`.  Applied to the preceding query, this function gives:

```{r, eval=FALSE}
opq_string(q)
## [out:xml][timeout:25];
## (
##   node
##     ["highway"="motorway"]
##     (51.2867602,-0.510375,51.6918741,0.3340155);
##   way
##     ["highway"="motorway"]
##     (51.2867602,-0.510375,51.6918741,0.3340155);
##   relation
##     ["highway"="motorway"]
##     (51.2867602,-0.510375,51.6918741,0.3340155);
## );
## (._;>);out body;
```

The resultant output may be pasted directly into the 
[overpass-turbo](https://overpass-turbo.eu/) online interactive query builder.
(The output of `opq_string` has been somewhat reformatted here to reflect
the format typically used in `overpass-turbo`.)

### 2.1. osmdata queries

As demonstrated above, an `osmdata` query begins by specifying a bounding box
with the function `opq()`, followed by specifying desired OSM features with
`add_osm_feature()`.  

```{r kunming1, eval = FALSE}
q <- opq(bbox = 'Kunming, China') %>%
    add_osm_feature(key = 'natural', value = 'water')
```
This query will request all natural water water bodies in Kunming, China. A
particular water body may be requested through appending a further call to
`add_osm_feature()`:

```{r kunming2, eval = FALSE}
q <- opq(bbox = 'Kunming, China') %>%
    add_osm_feature(key = 'natural', value = 'water') %>%
    add_osm_feature(key = 'name:en', value = 'Dian', value_exact = FALSE)
```
```{r, echo = FALSE}
q <- opq(bbox = c(102.5417638, 24.8915153, 102.8617638, 25.2115153)) %>%
    add_osm_feature(key = 'natural', value = 'water') %>%
    add_osm_feature(key = 'name:en', value = 'Dian', value_exact = FALSE)
```
Each successive call to `add_osm_feature()` **adds** features to a query. This query
is thus a request for all bodies of natural water **and** those with English
names that include 'Dian'. The requested data may be extracted through calling
one of the `osmdata_xml/sp/sf()` functions.

Single queries are always constructed through **adding** features, and therefore
correspond to logical **AND** operations: natural water bodies **AND** those
whose names include 'Dian'.  The equivalent **OR** combination can be extracted
with the [`add_osm_features()`
function](https://docs.ropensci.org/osmdata/reference/add_osm_features.html).
The following query represents the OR-equivalent of the above query, requesting
data on both all natural features with the value of `"water"` OR all features
whose English name is `"Dian"`.

```{r add_osm_features-fakey, eval = FALSE}
q <- opq(bbox = 'Kunming, China') %>%
    add_osm_features(features = c ("\"natural\"=\"water\"",
                                   "\"name:en\"=\"Dian\""))
```
```{r add_osm_features, echo = FALSE}
q <- list (
    bbox = "24.388848,102.1697441,26.548485,103.6683522",
    prefix = "[out:xml][timeout:25];\n(\n",
    suffix = ");\n(._;>;);\nout body;",
    features = c ("[\"natural\"=\"water\"]", "[\"name:en\"=\"Dian\"]")
    )
attr(q, "class") <- c ("list", "overpass_query")
attr(q, "nodes_only") <- FALSE
```

Note that the `"="` symbols here requests features whose values exactly match
the given values. Other "filter" symbols are possible, as described in the
[overpass query language
definition](https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#By_tag_.28has-kv.29),
including symbols for negation (`!=`), or approximate matching (`~`).

Passing this query to
[`osmdata_sf()`](https://docs.ropensci.org/osmdata/reference/osmdata_sf.html)
will return identical data to the following way to explicitly construct an OR
query through using the inbuilt `c` operator
of `osmdata`.

```{r kunming3, eval = FALSE}
dat1 <- opq(bbox = 'Kunming, China') %>%
    add_osm_feature(key = 'natural', value = 'water') %>%
    osmdata_sf ()
dat2 <- opq(bbox = 'Kunming, China') %>%
    add_osm_feature(key = 'name:en', value = 'Dian', value_exact = FALSE) %>%
    osmdata_sf ()
dat <- c (dat1, dat2)
```

While the "filter" symbols may be explicitly specified in 
[the `add_osm_features()`
function](https://docs.ropensci.org/osmdata/reference/add_osm_features.html),
the single-feature version of 
[`add_osm_feature()`
function](https://docs.ropensci.org/osmdata/reference/add_osm_feature.html) has
several logical parameters to control matching without needing to remember
precise overpass syntax:

- `key_exact` can be set to `FALSE` to approximately match given keys;
- `value_exact` can be set to `FALSE` to approximately match given values; and
- `match_case` can be set to `FALSE` to match keys and values in both lower and
  upper case forms.

The previous query with `key = 'name:end'` and `value = 'Dian'` could thus be
replaced by the following:

```{r kunming4, eval = FALSE}
add_osm_feature(key = 'name', value = 'dian',
                key_exact = FALSE,
                value_exact = FALSE,
                match_case = FALSE)
```

### 2.2 Extracting `OSM` data from a query

The primary `osmdata` functions `osmdata_sf()` or `osmdata_sp()` pass these
queries to `overpass` and return OSM data in corresponding `sf` or `sp` format,
respectively.  Both of these functions also accept direct `overpass` queries,
such as those produced by the `osmdata` function `opq_string()`, or copied
directly from the [`overpass-turbo` query builder](https://overpass-turbo.eu).

```{r, eval=FALSE}
osmdata_sf(opq_string(q))
## Object of class 'osmdata' with:
##                  $bbox :
##         $overpass_call : The call submitted to the overpass API
##             $timestamp : [ Thurs 5 May 2017 14:33:54 ]
##            $osm_points : 'sf' Simple Features Collection with 360582 points
##            ...
```

Note that the result contains no value for `bbox`, because that information is
lost when the full `osmdata_query`, `q`, is converted to a string.
Nevertheless, the results of the two calls `osmdata_sf (opq_string (q))` and
`osmdata_sf (q)` differ only in the values of `bbox` and `timestamp`, while
returning otherwise identical data.

In summary, `osmdata` queries are generally simplified versions of potentially
more complex `overpass` queries, although arbitrarily complex `overpass` queries
may be passed directly to the primary `osmdata` functions.  As illustrated
above, `osmdata` queries are generally constructed through initiating a query
with `opq()`, and then specifying OSM features in terms of `key-value` pairs
with `add_osm_feature()`, along with judicious usage of the `key_exact`,
`value_exact`, and `match_case` parameters.

The simplest way to use `osmdata` is to simply request all data within a given
bounding box (warning - not intended to run):

```{r, eval=FALSE}
q <- opq(bbox = 'London City, U.K.')
lots_of_data <- osmdata_sf(q)
```

Queries are, however, usually more useful when refined through using
`add_osm_feature()`, which minimally requires a single `key` and returns all objects
specifying any value for that `key`: 

```{r opq-london, eval = FALSE}
not_so_much_data <- opq(bbox = 'city of london uk') %>%
    add_osm_feature(key = 'highway') %>%
    add_osm_feature(key = 'name') %>%
    osmdata_sf()
```

`osmdata` will use that query to return all named highways within the requested
bounding box. Note that `key` specifications are requests for features which
must include those keys, yet most features will also include many other keys,
and thus `osmdata` objects generally list a large number of distinct keys, as
demonstrated below.

### 2.3. Query example

To appreciate query building in more concrete terms, let's imagine that we
wanted to find all cycle paths in Seville, Spain:

```{r opq-seville-plot, eval = FALSE}
q1 <- opq('Sevilla') %>%
    add_osm_feature(key = 'highway', value = 'cycleway')
cway_sev <- osmdata_sp(q1)
sp::plot(cway_sev$osm_lines)
```

![](https://cloud.githubusercontent.com/assets/1825120/23708182/fba4ed96-040c-11e7-90cd-3274d394030a.png)

Now imagine we want to make a more specific query that only extracts designated
cycleways or those which are bridges. Combining these into one query will return
only those that are designated cycleways **AND** that are bridges:

```{r des-bike1, eval=FALSE}
des_bike <- osmdata_sf(q1)
q2 <- add_osm_feature(q1, key = 'bridge', value = 'yes')
des_bike_and_bridge <- osmdata_sf(q2)
nrow(des_bike_and_bridge$osm_points); nrow(des_bike_and_bridge$osm_lines)
## [1] 99
## [1] 32
```

That query returns only 99 points and 32 lines.  Designed cycleways **OR** bridges
can be obtained through simply combining multiple `osmdata` objects with the `c`
operator:

```{r, des-bike2, eval=FALSE}
q2 <- opq('Sevilla') %>%
    add_osm_feature(key = 'bridge', value = 'yes')
bridge <- osmdata_sf(q2)
des_bike_or_bridge <- c(des_bike, bridge)
nrow(des_bike_or_bridge$osm_points); nrow(des_bike_or_bridge$osm_lines)
## [1] 9757
## [1] 1061
```

And as expected, the `OR` operation produces more data than the equivalent
`AND`, showing the utility of combining `osmdata` objects with the generic
function `c()`. 

## 3. The `osmdata` object

The `osmdata` extraction functions (`osmdata_sf()` and `osmdata_sp()`), both
return objects of class `osmdata`.  The structure of `osmdata` objects are clear
from their default print method, illustrated using the `bridge` example from the
previous section:

```{r, eval=FALSE}
bridge
##  Object of class 'osmdata' with:
##                   $bbox : 37.3002036,-6.0329182,37.4529579,-5.819157
##          $overpass_call : The call submitted to the overpass API
##              $timestamp : [ Thurs 5 May 2017 14:41:19 ]
##             $osm_points : 'sf' Simple Features Collection with 69 points
##              $osm_lines : 'sf' Simple Features Collection with 25 linestrings
##           $osm_polygons : 'sf' Simple Features Collection with 0 polygons
##         $osm_multilines : 'sf' Simple Features Collection with 0 multilinestrings
##      $osm_multipolygons : 'sf' Simple Features Collection with 0 multipolygons
```

As the results show, all `osmdata` objects should contain:

- A bounding box (which can be accessed with `bridge$bbox`)
- A time-stamp of the query (`bridge$timestamp`, useful for checking data is
  up-to-date)
- The spatial data, consisting of `osm_points`, `osm_lines`, `osm_polygons`,
  `osm_multilines` and `osm_multipolygons`.

 Some or all of these can be empty: the example printed above contains only
 points and lines. The more complex features of `osm_multilines` and
 `osm_multipolygons` refer to OSM relations than contain multiple lines and
 polygons.

The actual spatial data contained in an `osmdata` object are of either `sp`
format when extracted with `osmdata_sp()` or `sf` format when extracted with
`osmdata_sf()`. 

```{r osmdata_with_files3a, eval=FALSE}
class(osmdata_sf(q)$osm_lines)
## [1] "sf"         "data.frame"
```

```{r osmdata_with_files3b, eval=FALSE}
class(osmdata_sp(q)$osm_lines)
## [1] "SpatialLinesDataFrame"
## attr(,"package")
## [1] "sp"
```

In addition to these two functions, `osmdata` provides a third function,
`osmdata_xml()`, which allows raw OSM data to be returned and optionally saved
to disk in XML format. The following code demonstrates this function, beginning
with a new query.

```{r osmdata_xml-london-buildings, eval = FALSE}
dat <- opq(bbox = c(-0.12, 51.51, -0.11, 51.52)) %>%
    add_osm_feature(key = 'building') %>%
    osmdata_xml(file = 'buildings.osm')
class(dat)
## [1] "xml_document" "xml_node"
```

This call both returns the same data as the object `dat` and saves them to the
file `buildings.osm`. Downloaded XML data can be converted to `sf` or `sp`
formats by simply passing the data to the respective `osmdata` functions, either
as the name of a file or an XML object:

```{r osmdata_with_files, eval = FALSE}
q <- opq(bbox = c(-0.12, 51.51, -0.11, 51.52)) %>%
    add_osm_feature(key = 'building')
doc <- osmdata_xml(q, 'buildings.osm')
dat1 <- osmdata_sf(q, doc)
dat2 <- osmdata_sf(q, 'buildings.osm')
identical(dat1, dat2)
## [1] TRUE
```

The following sub-sections now explore these three functions in more detail,
beginning with `osmdata_xml()`.

### 3.1. The `osmdata_xml()` function

`osmdata_xml()` returns OSM data in native XML format, and also allows these
data to be saved directly to disk (conventionally using the file suffix `.osm`,
although any suffix may be used). The `XML` data are formatting using the `R`
package `xml2`, and may be processed within `R` using any methods compatible
with such data, or may be processed by any other software able to load the `XML`
data directly from disk.

The first few lines of the XML data downloaded above look like this:

```{r, eval=FALSE}
readLines('buildings.osm')[1:6]
## [1] "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
## [2] "<osm version=\"0.6\" generator=\"Overpass API\">"
## [3] "  <note>The data included in this document is from www.openstreetmap.org. The data is made available under ODbL.</note>"
## [4] "  <meta osm_base=\"2017-03-07T09:28:03Z\"/>"
## [5] "  <node id=\"21593231\" lat=\"51.5149566\" lon=\"-0.1134203\"/>"
## [6] "  <node id=\"25378129\" lat=\"51.5135870\" lon=\"-0.1115193\"/>"
```

These data can be used in any other programs able to read and process XML data,
such as the open source GIS [QGIS](https://qgis.org/en/site/) or the OSM data
editor [JOSM](https://wiki.openstreetmap.org/wiki/JOSM).  The remainder of this
vignette assumes that not only do you want to get OSM data using R, you also
want to import and eventually process it, using R. For that you'll need to
import the data into a native R class.

As demonstrated above, downloaded data can be directly processed by passing
either filenames or the `R` objects containing those data to the
`osmdata_sf/sp()` functions:

```{r, eval=FALSE}
dat_sp <- osmdata_sp(q, 'buildings.osm')
dat_sf <- osmdata_sf(q, 'buildings.osm')
```

### 3.2. The `osmdata_sf()` function

`osmdata_sf()` returns OSM data in 
[Simple Features (SF)](https://www.ogc.org/standards/sfo) 
format, defined by the
[Open Geospatial Consortium](https://www.ogc.org), and implemented in
the `R` package [`sf`](https://cran.r-project.org/package=sf). This package
provides a direct interface to the `C++` 
[Graphical Data Abstraction Library (GDAL)](https://gdal.org) which also includes
a so-called ['driver' for OSM data](https://gdal.org/drv_osm.html). This
means that OSM data may also be read directly with `sf`, rather than using
`osmdata`. In this case, data must first be saved to disk, which can 
be readily achieved using `osmdata_xml()` described above, or through
downloading directly from the [overpass interactive query
builder](https://overpass-turbo.eu). 

The following example is based on this query:

```{r trentham, eval = FALSE}
opq(bbox = 'Trentham, Australia') %>%
    add_osm_feature(key = 'name') %>%
    osmdata_xml(filename = 'trentham.osm')
```

`sf` can then read such data independent of `osmdata` though:

```{r sf1, eval=FALSE}
sf::st_read('trentham.osm', layer = 'points')
## Reading layer `points' from data source `trentham.osm' using driver `OSM'
## Simple feature collection with 38 features and 10 fields
## geometry type:  POINT
## dimension:      XY
## bbox:           xmin: 144.2894 ymin: -37.4846 xmax: 144.3893 ymax: -37.36012
## epsg (SRID):    4326
## proj4string:    +proj=longlat +datum=WGS84 +no_defs
```

The `GDAL` drivers used by `sf` can only load single 'layers' of features, for
example, `points`, `lines`, or `polygons`. In contrast, `osmdata` loads all
features simultaneously:

```{r osmdata_sf2, eval=FALSE}
osmdata_sf(q, 'trentham.osm')
## Object of class 'osmdata' with:
##                  $bbox : -37.4300874,144.2863388,-37.3500874,144.3663388
##         $overpass_call : The call submitted to the overpass API
##             $timestamp : [ Thus 5 May 2017 14:42:19 ]
##            $osm_points : 'sf' Simple Features Collection with 7106 points
##             $osm_lines : 'sf' Simple Features Collection with 263 linestrings
##          $osm_polygons : 'sf' Simple Features Collection with 38 polygons
##        $osm_multilines : 'sf' Simple Features Collection with 1 multilinestrings
##     $osm_multipolygons : 'sf' Simple Features Collection with 6 multipolygons
```

Even for spatial objects of the same type (the same 'layers' in `sf`
terminology), `osmdata` returns considerably more objects--7,166 points compared
.with just 38. The raw sizes of data returned can be compared with:

```{r object-sizes, eval=FALSE}
s1 <- object.size(osmdata_sf(q, 'trentham.osm')$osm_points)
s2 <- object.size(sf::st_read('trentham.osm', layer = 'points', quiet = TRUE))
as.numeric(s1 / s2)
## [1] 511.4193
```

And the `osmdata points` contain over 500 times as much data.  The primary
difference between `sf/GDAL` and `osmdata` is that the former returns only those
objects unique to each category of spatial object.  Thus OSM nodes (`points` in
`sf/osmdata` representations) include, in `sf/GDAL` representation, only those
points which are not part of any other objects (such as lines or polygons). In
contrast, the `osm_points` object returned by `osmdata` includes all points
regardless of whether or not these are represented in other spatial objects.
Similarly, `line` objects in `sf/GDAL` exclude any lines that are part of other
objects such as `multipolygon` or `multiline` objects.

This processing of data by `sf/GDAL` has two important implications:

1. An implicit hierarchy of spatial objects is enforced through including
elements of objects only at their 'highest' level of representation, where
`multipolygon` and `multiline` objects are assumed to be at 'higher' levels
than `polyon` or `line` objects, and these in turn are at 'higher' levels than
`point` objects. `osmdata` makes no such hierarchical assumptions.

2. All OSM are structured by giving each object a unique identifier so that the
components of any given object (the nodes of a line, for example, or the lines
of a multipolygon) can be described simply by giving these identifiers.  This
enables the components of any OSM object to be examined in detail.
The `sf/GDAL` representation obviates this ability
through removing these IDs and reducing everything to geometries alone (which
is, after all, why it is called '*Simple* Features'). This means, for example,
that the `key-value` pairs of the `line` or `polygon` components of
`multipolygon` can never be extracted from an `sf/GDAL` representation. In
contrast, `osmdata` retains all unique identifiers for all OSM objects, and so
readily enables, for example, the properties of all `point` objects of a `line`
to be extracted.

Another reason why `osmdata` returns more data than `GDAL/sf` is that the latter
extracts only a restricted list of OSM `keys`, whereas `osmdata` returns all
`key` fields present in the requested data:

```{r, eval=FALSE}
names(sf::st_read('trentham.osm', layer = 'points', quiet = TRUE)) # the keys
## [1] "osm_id"     "name"       "barrier"    "highway"
## [5] "ref"        "address"    "is_in"      "place"
## [9] "man_made"   "other_tags" "geometry"
```

```{r, eval=FALSE}
names(osmdata_sf(q, 'trentham.osm')$osm_points)
## [1] "osm_id"           "name"             "X_description_"   "X_waypoint_"
## [5] "addr.city"        "addr.housenumber" "addr.postcode"    "addr.street"
## [9] "amenity"          "barrier"          "denomination"     "foot"
## [13] "ford"             "highway"          "leisure"          "note_1"
## [17] "phone"            "place"            "railway"          "railway.historic"
## [21] "ref"              "religion"         "shop"             "source"
## [25] "tourism"          "waterway"         "geometry"
```

`key` fields which are not specified in a given set of OSM data are not returned
by `osmdata`, while `GDAL/sf` returns the same `key` fields regardless of
whether any values are specified.

```{r, eval=FALSE}
addr <- sf::st_read('trentham.osm', layer = 'points', quiet = TRUE)$address
all(is.na(addr))
## TRUE
```

and `key=address` contains no data yet is still returned by `GDAL/sf`.

Finally, note that `osmdata` will generally extract OSM data considerably faster
than equivalent `sf/GDAL` routines (as detailed
[here](https://github.com/ropensci/osmdata/wiki/Timing-benchmarks)).

### 3.3. The `osmdata_sp()` function

As with `osmdata_sf()` described above, OSM data may be converted to `sp`
format without using `osmdata` via the `sf` functions demonstrated below: 

```{r sf_sp, eval=FALSE}
dat <- sf::st_read('buildings.osm', layer = 'multipolygons', quiet = TRUE)
dat_sp <- as(dat, 'Spatial')
class(dat_sp)
## [1] "SpatialPolygonsDataFrame"\nattr(,"package")\n[1] "sp"
```

These data are extracted using the GDAL, and so suffer all of the same
shortcomings mentioned above. Note differences in the amount of data returned:

```{r, eval=FALSE}
dim(dat_sp)
## [1] 560  25
```

```{r, eval=FALSE}
dim(osmdata_sp(q, doc = 'buildings.osm')$osm_polygons)
## [1] 566 114
```

```{r, eval=FALSE}
dim(osmdata_sp(q, doc = 'buildings.osm')$osm_multipolygons)
## [1] 15 52
```

## 4. Recursive searching

As described above, `osmdata` returns all data of each type and so allows the
components of any given spatial object to be examined in their own right.  This
ability to extract, for example, all points of a line, or all polygons which
include a given set of points, is referred to as *recursive searching*.

Recursive searching is not possible with `GDAL/sf`, because OSM identifiers are
removed, and only the unique data of each type of object are retained. To
understand both recursive searching and why it is useful, note that OSM data are
structured in three hierarchical levels:

1. `nodes` representing spatial points 

2. `ways` representing lines, both as `polygons` (with connected ends) and
   non-polygonal `lines`

3. `relations` representing more complex objects generally comprising
   collections of `ways` and/or `nodes`. Examples include 
   `multipolygon relations` comprising an outer polygon (which may itself be
   made of several distinct `ways` which ultimately connect to form a single
   circle), and several inner polygons.

Recursive searching allows for objects within any one of these hierarchical
levels to be extracted based on components in any other level. Recursive
searching is performed in `osmdata` with the following functions:

1. `osm_points()`, which extracts all `point` or `node` objects 

2. `osm_lines()`, which extracts all `way` objects that are `lines` (that are, that are
   not `polygons`) 

3. `osm_polygons()`, which extracts all `polygon` objects

4. `osm_multilines()`, which extracts all `multiline` objects; and

5. `osm_multipolygons()`, which extracts all `multipolygon` objects.

Each of these functions accepts as an argument a vector of OSM identifiers. To
demonstrate these functions, we first re-create the example above of named
objects from Trentham, Australia:

```{r, eval = FALSE}
tr <- opq(bbox = 'Trentham, Australia') %>%
    add_osm_feature(key = 'name') %>%
    osmdata_sf()
```

### 4.1. Example

Then imagine we are interested in the `osm_line` object describing the 'Coliban
River':

```{r, eval=FALSE}
i <- which(tr$osm_lines$name == 'Coliban River')
coliban <- tr$osm_lines[i, ]
coliban[which(!is.na(coliban))]
## Simple feature collection with 1 feature and 3 fields
## geometry type:  LINESTRING
## dimension:      XY
## bbox:           xmin: 144.3235 ymin: -37.37162 xmax: 144.3335 ymax: 37.36366
## epsg (SRID):    4326
## proj4string:    +proj=longlat +datum=WGS84 +no_defs
##            osm_id          name waterway                       geometry
## 87104907 87104907 Coliban River    river LINESTRING(144.323471069336...
```

The locations of the points of this line can be extracted directly from the `sf`
object with:

```{r, eval=FALSE}
coliban$geometry[[1]]
## LINESTRING(144.323471069336 -37.3716201782227, 144.323944091797 -37.3714790344238, 144.324356079102 -37.3709754943848, 144.324493408203 -37.3704833984375, 144.324600219727 -37.370174407959, 144.324981689453 -37.3697204589844, 144.325149536133 -37.369441986084, 144.325393676758 -37.3690567016602, 144.325714111328 -37.3686943054199, 144.326080322266 -37.3682441711426)
```

The output contains nothing other than geometries (because, to reiterate, these
are '*Simple* Features'), and no further information regarding those points can
be extracted. The Coliban River has a waterfall in Trentham, and one of the
`osm_points` objects describes this waterfall.  The information necessary to
locate this waterfall is removed from the `GDAL/sf` representation, but can be
extracted with `osmdata` with the following lines, noting that the 
OSM ID of the line `coliban` is given by `rownames(coliban)`.

```{r, eval=FALSE}
pts <- osm_points(tr, rownames(coliban))
wf <- pts[which(pts$waterway == 'waterfall'), ]
wf[which(!is.na(wf))]
## Simple feature collection with 1 feature and 4 fields
## geometry type:  POINT
## dimension:      XY
## bbox:           xmin: 144.3246 ymin: -37.37017 xmax: 144.3246 ymax: -37.37017
## epsg (SRID):    4326
## proj4string:    +proj=longlat +datum=WGS84 +no_defs
##                osm_id           name    tourism  waterway
## 1013064837 1013064837 Trentham Falls attraction waterfall
##                                  geometry
## 1013064837 POINT(144.324600219727 -37....
```

This point could be used as the basis for further recursive searches. For
example, all `multipolygon` objects which include Trentham Falls could be
extracted with:

```{r, eval=FALSE}
mp <- osm_multipolygons(tr, rownames(wf))
```

Although this returns no data in this case, it nevertheless demonstrates the
usefulness and ease of recursive searching with `osmdata`.

```{r, echo=FALSE}
for (f in list.files(pattern = "\\.osm"))
    if (file.exists(f)) file.remove(f)
```

### 4.2 Relation example

A special type of OSM object is a relation. These can be defined by their name,
which can join many divers features into a single object.
The following example extracts the London Route Network Route 9,
which is composed of many (over 100) separate lines:

```{r, eval = FALSE}
lcnr9 <- opq ('greater london uk') %>%
    add_osm_feature (key = "name", value = "LCN 9",
                 value_exact = FALSE) %>%
    osmdata_sp()
sp::plot(lcnr9$osm_lines)
```

![](https://cloud.githubusercontent.com/assets/1825120/23709879/c98e2c2c-0412-11e7-86b8-1ffc95aab5a1.png)

## 5. Additional Functionality

This section briefly describes a few of additional functions, with additional
detail provided in the help files for each of these function.

1. The `trim_osmdata()` function, as described above in the sub-section on
   bounding boxes, trims an `osmdata` object to within a defined bounding
   *polygon*, rather than bounding box.
2. The `opq_osm_id()` function allows queries for particular OSM objects by
   their OSM-allocated ID values.
3. The `osm_poly2line()` function converts all `$osm_polygons` items of an
   `osmdata` object to `$osm_lines`. These objects remain polygonal in form,
   through sharing identical start and end points, but can then be treated as
   simple lines. This is important for polygonal highways, which are
   automatically classified as `$osm_polygons` simply because they form closed
   loops. The function enables all highways to be grouped together (as
   `$osm_lines`) regardless of the form.
4. The `unique_osmdata()` function removes redundant items from the different
   components of an `osmdata` object. A multilinestring, for example, is
   composed of multiple lines, and each line is composed of multiple points. For
   a multilinestring, an `osmdata` object will thus contain several
   `$osm_lines`, and for each of these several `$osm_points`. This function
   removes all of these redundant objects, so that `$osm_lines` only contains
   lines which are not part of any higher-level objects, and `$osm_points` only
   contains points which are not part of any higher-level objects.

A further additional function is the ability to extract data as represented in
the OSM database prior to a specified date, or within a specified range of
dates. This is achieved by passing one or both values to the [`opq()`
function](https://docs.ropensci.org/osmdata/reference/opq.html) of `datetime`
and `datetime2`. The resultant data extracted with one or more
`add_osm_feature()` calls and an extraction function (`osmdata_sf/sp/sc/xml`)
will then contain only those data present prior to the specified date (when
`datetime` only given), or between the two specified dates (when both
`datetime` and `datetime2` given).


## 6. Related Packages

@eugster_osmar:_2012 describe `osmar`, an R package for handling OSM data
that enables visualisation, search and even rudimentary routing operations.
`osmar` is not user friendly or able to download OSM data flexibly,
as reported in an [early tutorial](https://eprints.whiterose.ac.uk/77643/)
comparing R and QGIS for handling OSM data [@lovelace_harnessing_2014]. Note
also that the `osmar` package does not work at present, and can not be used for
accessing OSM data.

`osmdata` builds on two previous R packages:
`osmplotr`, a package [available from CRAN](https://cran.r-project.org/package=osmplotr)
for accessing and plotting OSM data [@osmplotr]
and `overpass`, a [GitHub package](https://github.com/hrbrmstr/overpass)
by Bob Rudis that provides an R interface to the 
[overpass](https://overpass-api.de/) API.

## 7. References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-osmdata.R
\name{osmdata_sp}
\alias{osmdata_sp}
\title{Return an OSM Overpass query as an \link{osmdata} object in \pkg{sp}
format.}
\usage{
osmdata_sp(q, doc, quiet = TRUE)
}
\arguments{
\item{q}{An object of class \code{overpass_query} constructed with
\link{opq} and \link{add_osm_feature}. May be be omitted,
in which case the \link{osmdata} object will not include the
query.}

\item{doc}{If missing, \code{doc} is obtained by issuing the overpass query,
\code{q}, otherwise either the name of a file from which to read data,
or an object of class \pkg{XML} returned from
\link{osmdata_xml}.}

\item{quiet}{suppress status messages.}
}
\value{
An object of class \code{osmdata} with the OSM components (points, lines,
and polygons) represented in \pkg{sp} format.
}
\description{
Return an OSM Overpass query as an \link{osmdata} object in \pkg{sp}
format.
}
\examples{
\dontrun{
hampi_sp <- opq ("hampi india") \%>\%
            add_osm_feature (key="historic", value="ruins") \%>\%
            osmdata_sp ()
}
}
\seealso{
Other extract: 
\code{\link{osmdata_sc}()},
\code{\link{osmdata_sf}()},
\code{\link{osmdata_xml}()}
}
\concept{extract}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poly2line.R
\name{osm_poly2line}
\alias{osm_poly2line}
\title{Convert osmdata polygons into lines}
\usage{
osm_poly2line(osmdat)
}
\arguments{
\item{osmdat}{An \link{osmdata} object.}
}
\value{
Modified version of same object with all \code{osm_polygons}
objects merged into \code{osm_lines}.
}
\description{
Street networks downloaded with \code{add_osm_object(key = "highway")} will
store any circular highways in \code{osm_polygons}. this function combines
those with the \code{osm_lines} component to yield a single \pkg{sf}
\code{data.frame} of all highways, whether polygonal or not.
}
\note{
The \code{osm_polygons} field is retained, with those features also
repeated as \code{LINESTRING} objects in \code{osm_lines}.
}
\examples{
\dontrun{
dat <- opq ("colchester uk") \%>\%
            add_osm_feature (key="highway") \%>\%
            osmdata_sf ()
# colchester has lots of roundabouts, and these are stored in 'osm_polygons'
# rather than 'osm_lines'. The former can be merged with the latter by:
dat2 <- osm_poly2line (dat)
# 'dat2' will have more lines than 'dat', but the same number of polygons
# (they are left unchanged.)
}
}
\seealso{
Other transform: 
\code{\link{osm_elevation}()},
\code{\link{trim_osmdata}()},
\code{\link{unique_osmdata}()},
\code{\link{unname_osmdata_sf}()}
}
\concept{transform}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unique-osmdata.R
\name{unique_osmdata}
\alias{unique_osmdata}
\title{unique_osmdata}
\usage{
unique_osmdata(dat)
}
\arguments{
\item{dat}{An \link{osmdata} object}
}
\value{
Equivalent object reduced to only unique objects of each type
}
\description{
Reduce the components of an \link{osmdata} object to only unique items of
each type. That is, reduce \verb{$osm_points} to only those points not
present in other objects (lines, polygons, etc.); reduce \verb{$osm_lines} to
only those lines not present in multiline objects; and reduce
\verb{$osm_polygons} to only those polygons not present in multipolygon
objects. This renders an \link{osmdata} object more directly compatible with
typical output of \pkg{sf}.
}
\seealso{
Other transform: 
\code{\link{osm_elevation}()},
\code{\link{osm_poly2line}()},
\code{\link{trim_osmdata}()},
\code{\link{unname_osmdata_sf}()}
}
\concept{transform}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osm-extract.R
\name{osm_multipolygons}
\alias{osm_multipolygons}
\title{Extract all \code{osm_multipolygons} from an osmdata object}
\usage{
osm_multipolygons(dat, id)
}
\arguments{
\item{dat}{An object of class \link{osmdata}}

\item{id}{OSM identification of one or more objects for which multipolygons
are to be extracted}
}
\value{
An \pkg{sf} Simple Features Collection of multipolygons
}
\description{
\code{id} must be of an \code{osm_points}, \code{osm_lines}, or
\code{osm_polygons} object. \code{osm_multipolygons} returns any multipolygon
object(s) which contain the object specified by \code{id}.
}
\examples{
\dontrun{
# find all multipolygons which contain the single polygon called
# "Chiswick Eyot" (which is an island).
dat <- opq ("London UK") \%>\%
    add_osm_feature (key="name", value="Thames", exact=FALSE) \%>\%
    osmdata_sf ()
index <- which (dat$osm_multipolygons$name == "Chiswick Eyot")
id <- rownames (dat$osm_polygons [id, ])
osm_multipolygons (dat, id)
# That multipolygon is the Thames itself, but note that
nrow (dat$osm_multipolygons) # = 14 multipolygon objects
nrow (osm_multipolygons (dat, id)) # = 1 - the main Thames multipolygon
}
}
\seealso{
Other search: 
\code{\link{osm_lines}()},
\code{\link{osm_multilines}()},
\code{\link{osm_points}()},
\code{\link{osm_polygons}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opq.R
\name{opq}
\alias{opq}
\title{Build an Overpass query}
\usage{
opq(
  bbox = NULL,
  nodes_only = FALSE,
  datetime = NULL,
  datetime2 = NULL,
  timeout = 25,
  memsize
)
}
\arguments{
\item{bbox}{Either (i) four numeric values specifying the maximal and minimal
longitudes and latitudes, in the form \code{c(xmin, ymin, xmax, ymax)}
or (ii) a character string in the form \code{xmin,ymin,xmax,ymax}. These
will be passed to \link{getbb} to be converted to a numerical bounding
box. Can also be (iii) a matrix representing a bounding polygon as
returned from \code{getbb(..., format_out = "polygon")}.}

\item{nodes_only}{If \code{TRUE}, query OSM nodes only. Some OSM structures such
as \code{place = "city"} or \code{highway = "traffic_signals"} are represented by
nodes only. Queries are built by default to return all nodes, ways, and
relation, but this can be very inefficient for node-only queries.
Setting this value to \code{TRUE} for such cases makes queries more
efficient, with data returned in the \code{osm_points} list item.}

\item{datetime}{If specified, a date and time to extract data from the OSM
database as it was up to the specified date and time, as described at
\url{https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#date}.
This \emph{must} be in ISO8601 format ("YYYY-MM-DDThh:mm:ssZ"), where
both the "T" and "Z" characters must be present.}

\item{datetime2}{If specified, return the \emph{difference} in the OSM
database between \code{datetime} and \code{datetime2}, where
\code{datetime2 > datetime}. See
\url{https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#Delta_between_two_dates_.28.22diff.22.29}.}

\item{timeout}{It may be necessary to increase this value for large queries,
because the server may time out before all data are delivered.}

\item{memsize}{The default memory size for the 'overpass' server in \emph{bytes};
may need to be increased in order to handle large queries.}
}
\value{
An \code{overpass_query} object
}
\description{
Build an Overpass query
}
\note{
See
\url{https://wiki.openstreetmap.org/wiki/Overpass_API#Resource_management_options_.28osm-script.29}
for explanation of \code{timeout} and \code{memsize} (or \code{maxsize} in overpass terms).
Note in particular the comment that queries with arbitrarily large \code{memsize}
are likely to be rejected.
}
\examples{
\dontrun{
q <- getbb ("portsmouth", display_name_contains = "United States") \%>\%
            opq () \%>\%
            add_osm_feature("amenity", "restaurant") \%>\%
            add_osm_feature("amenity", "pub")
osmdata_sf (q) # all objects that are restaurants AND pubs (there are none!)
q1 <- getbb ("portsmouth", display_name_contains = "United States") \%>\%
                opq () \%>\%
                add_osm_feature("amenity", "restaurant")
q2 <- getbb ("portsmouth", display_name_contains = "United States") \%>\%
                opq () \%>\%
                add_osm_feature("amenity", "pub")
c (osmdata_sf (q1), osmdata_sf (q2)) # all restaurants OR pubs

# Use nodes_only to retrieve single point data only, such as for central
# locations of cities.
opq <- opq (bbox, nodes_only = TRUE) \%>\%
    add_osm_feature (key = "place", value = "city") \%>\%
    osmdata_sf (quiet = FALSE)
}
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osm-extract.R
\name{osm_multilines}
\alias{osm_multilines}
\title{Extract all \code{osm_multilines} from an osmdata object}
\usage{
osm_multilines(dat, id)
}
\arguments{
\item{dat}{An object of class \link{osmdata}}

\item{id}{OSM identification of one of more objects for which multilines are
to be extracted}
}
\value{
An \pkg{sf} Simple Features Collection of multilines
}
\description{
\code{id} must be of an \code{osm_points} or \code{osm_lines} object (and can
not be the \code{id} of an \code{osm_polygons} object because multilines by
definition contain no polygons.  \code{osm_multilines} returns any multiline
object(s) which contain the object specified by \code{id}.
}
\examples{
\dontrun{
dat <- opq ("London UK") \%>\%
    add_osm_feature (key="name", value="Thames", exact=FALSE) \%>\%
    osmdata_sf ()
# Get ids of lines called "The Thames":
id <- rownames (dat$osm_lines [which (dat$osm_lines$name == "The Thames"),])
# and find all multilinestring objects which include those lines:
osm_multilines (dat, id)
# Now note that
nrow (dat$osm_multilines) # = 24 multiline objects
nrow (osm_multilines (dat, id)) # = 1 - the recursive search selects the
                                # single multiline containing "The Thames"
}
}
\seealso{
Other search: 
\code{\link{osm_lines}()},
\code{\link{osm_multipolygons}()},
\code{\link{osm_points}()},
\code{\link{osm_polygons}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opq.R
\name{opq_around}
\alias{opq_around}
\title{opq_around}
\usage{
opq_around(lon, lat, radius = 15, key = NULL, value = NULL, timeout = 25)
}
\arguments{
\item{lon}{Longitude of desired point}

\item{lat}{Latitude of desired point}

\item{radius}{Radius in metres around the point for which data should be
extracted. Queries with  large values for this parameter may fail.}

\item{key}{(Optional) OSM key of enclosing data}

\item{value}{(Optional) OSM value matching 'key' of enclosing data}

\item{timeout}{It may be necessary to increase this value for large queries,
because the server may time out before all data are delivered.}
}
\description{
Find all features around a given point, and optionally match specific
'key'-'value' pairs. This function is \emph{not} intended to be combined with
\link{add_osm_feature}, rather is only to be used in the sequence
\link{opq_around} -> \link{osmdata_xml} (or other extraction function). See
examples for how to use.
}
\examples{
\dontrun{
# Get all benches ("amenity=bench") within 100m of a particular point
lat <- 53.94542
lon <- -2.52017
key <- "amenity"
value <- "bench"
radius <- 100
x <- opq_around (lon, lat, radius, key, value) \%>\%
    osmdata_sf ()
}
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osm-extract.R
\name{osm_points}
\alias{osm_points}
\title{Extract all \code{osm_points} from an osmdata object}
\usage{
osm_points(dat, id)
}
\arguments{
\item{dat}{An object of class \link{osmdata}}

\item{id}{OSM identification of one or more objects for which points are to
be extracted}
}
\value{
An \pkg{sf} Simple Features Collection of points
}
\description{
Extract all \code{osm_points} from an osmdata object
}
\examples{
\dontrun{
tr <- opq ("trentham australia") \%>\% osmdata_sf ()
coliban <- tr$osm_lines [which (tr$osm_lines$name == "Coliban River"),]
pts <- osm_points (tr, rownames (coliban)) # all points of river
waterfall <- pts [which (pts$waterway == "waterfall"),] # the waterfall point
}
}
\seealso{
Other search: 
\code{\link{osm_lines}()},
\code{\link{osm_multilines}()},
\code{\link{osm_multipolygons}()},
\code{\link{osm_polygons}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trim-osmdata.R
\name{trim_osmdata}
\alias{trim_osmdata}
\title{trim_osmdata}
\usage{
trim_osmdata(dat, bb_poly, exclude = TRUE)
}
\arguments{
\item{dat}{An \link{osmdata} object returned from \link{osmdata_sf} or
\link{osmdata_sp}.}

\item{bb_poly}{A matrix representing a bounding polygon obtained with
\code{getbb (..., format_out = "polygon")} (and possibly selected from
resultant list where multiple polygons are returned).}

\item{exclude}{If TRUE, objects are trimmed exclusively, only retaining those
strictly within the bounding polygon; otherwise all objects which partly
extend within the bounding polygon are retained.}
}
\value{
A trimmed version of \code{dat}, reduced only to those components
lying within the bounding polygon.
}
\description{
Trim an \link{osmdata} object to within a bounding polygon
}
\note{
It will generally be necessary to pre-load the \pkg{sf} package for
this function to work correctly.

Caution is advised when using polygons obtained from Nominatim via
\code{getbb(..., format_out = "polygon"|"sf_polygon")}. These shapes can be
outdated and thus could cause the trimming operation to not give results
expected based on the current state of the OSM data.
}
\examples{
\dontrun{
dat <- opq ("colchester uk") \%>\%
            add_osm_feature (key="highway") \%>\%
            osmdata_sf (quiet = FALSE)
bb <- getbb ("colchester uk", format_out = "polygon")
library (sf) # required for this function to work
dat_tr <- trim_osmdata (dat, bb)
bb <- getbb ("colchester uk", format_out = "sf_polygon")
class (bb) # sf data.frame
dat_tr <- trim_osmdata (dat, bb)
bb <- as (bb, "Spatial")
class (bb) # SpatialPolygonsDataFrame
dat_tr <- trim_osmdata (dat, bb)
}
}
\seealso{
Other transform: 
\code{\link{osm_elevation}()},
\code{\link{osm_poly2line}()},
\code{\link{unique_osmdata}()},
\code{\link{unname_osmdata_sf}()}
}
\concept{transform}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-osmdata.R
\name{osmdata_sc}
\alias{osmdata_sc}
\title{Return an OSM Overpass query as an \link{osmdata} object in
\code{silicate} (\code{SC}) format.}
\usage{
osmdata_sc(q, doc, quiet = TRUE)
}
\arguments{
\item{q}{An object of class \code{overpass_query} constructed with
\link{opq} and \link{add_osm_feature}. May be be omitted,
in which case the \link{osmdata} object will not include the
query.}

\item{doc}{If missing, \code{doc} is obtained by issuing the overpass query,
\code{q}, otherwise either the name of a file from which to read data,
or an object of class \pkg{XML} returned from
\link{osmdata_xml}.}

\item{quiet}{suppress status messages.}
}
\value{
An object of class \code{osmdata} representing the original OSM hierarchy
of nodes, ways, and relations.
}
\description{
Return an OSM Overpass query as an \link{osmdata} object in
\code{silicate} (\code{SC}) format.
}
\note{
The \code{silicate} format is currently highly experimental, and
recommended for use only if you really know what you're doing.
}
\examples{
\dontrun{
hampi_sf <- opq ("hampi india") \%>\%
            add_osm_feature (key="historic", value="ruins") \%>\%
            osmdata_sc ()
}
}
\seealso{
Other extract: 
\code{\link{osmdata_sf}()},
\code{\link{osmdata_sp}()},
\code{\link{osmdata_xml}()}
}
\concept{extract}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osm-extract.R
\name{osm_lines}
\alias{osm_lines}
\title{Extract all \code{osm_lines} from an osmdata object}
\usage{
osm_lines(dat, id)
}
\arguments{
\item{dat}{An object of class \link{osmdata}}

\item{id}{OSM identification of one or more objects for which lines are to be
extracted}
}
\value{
An \pkg{sf} Simple Features Collection of linestrings
}
\description{
If \code{id} is of a point object, \code{osm_lines} will return all lines
containing that point. If \code{id} is of a line or polygon object,
\code{osm_lines} will return all lines which intersect the given line or
polygon.
}
\examples{
\dontrun{
dat <- opq ("hengelo nl") \%>\% add_osm_feature (key="highway") \%>\%
     osmdata_sf ()
bus <- dat$osm_points [which (dat$osm_points$highway == "bus_stop"),] \%>\%
        rownames () # all OSM IDs of bus stops
osm_lines (dat, bus) # all highways containing bus stops

# All lines which intersect with Piccadilly Circus in London, UK
dat <- opq ("Fitzrovia London") \%>\% add_osm_feature (key="highway") \%>\%
    osmdata_sf ()
i <- which (dat$osm_polygons$name == "Piccadilly Circus")
id <- rownames (dat$osm_polygons [i,])
osm_lines (dat, id)
}
}
\seealso{
Other search: 
\code{\link{osm_multilines}()},
\code{\link{osm_multipolygons}()},
\code{\link{osm_points}()},
\code{\link{osm_polygons}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getbb.R
\name{bbox_to_string}
\alias{bbox_to_string}
\title{Convert a named matrix or a named or unnamed vector to a string}
\usage{
bbox_to_string(bbox)
}
\arguments{
\item{bbox}{bounding box as character, matrix or vector.
If character, the bbox will be found (geocoded) and extracted with
\link{getbb}. Unnamed vectors will be sorted appropriately and must merely be
in the order (x, y, x, y).}
}
\value{
A character string representing min x, min y, max x, and max y
bounds. For example: \code{"15.3152361,76.4406446,15.3552361,76.4806446"} is
the bounding box for Hampi, India.
}
\description{
This function converts a bounding box into a string for use in web apis
}
\examples{
\dontrun{
bbox_to_string (getbb ("hampi india"))
}
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{get_overpass_url}
\alias{get_overpass_url}
\title{get_overpass_url}
\usage{
get_overpass_url()
}
\value{
The overpass API URL
}
\description{
Return the URL of the specified overpass API. Default is
\url{https://overpass-api.de/api/interpreter/}.
}
\seealso{
\code{\link[=set_overpass_url]{set_overpass_url()}}

Other overpass: 
\code{\link{set_overpass_url}()}
}
\concept{overpass}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-osmdata.R
\name{osmdata_sf}
\alias{osmdata_sf}
\title{Return an OSM Overpass query as an \link{osmdata} object in \pkg{sf}
format.}
\usage{
osmdata_sf(q, doc, quiet = TRUE, stringsAsFactors = FALSE)
}
\arguments{
\item{q}{An object of class \code{overpass_query} constructed with
\link{opq} and \link{add_osm_feature}. May be be omitted,
in which case the \link{osmdata} object will not include the
query.}

\item{doc}{If missing, \code{doc} is obtained by issuing the overpass query,
\code{q}, otherwise either the name of a file from which to read data,
or an object of class \pkg{XML} returned from
\link{osmdata_xml}.}

\item{quiet}{suppress status messages.}

\item{stringsAsFactors}{Should character strings in 'sf' 'data.frame' be
coerced to factors?}
}
\value{
An object of class \code{osmdata} with the OSM components (points, lines,
and polygons) represented in \pkg{sf} format.
}
\description{
Return an OSM Overpass query as an \link{osmdata} object in \pkg{sf}
format.
}
\examples{
\dontrun{
hampi_sf <- opq ("hampi india") \%>\%
            add_osm_feature (key="historic", value="ruins") \%>\%
            osmdata_sf ()
}
}
\seealso{
Other extract: 
\code{\link{osmdata_sc}()},
\code{\link{osmdata_sp}()},
\code{\link{osmdata_xml}()}
}
\concept{extract}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osmdata-class.R, R/osmdata-package.R
\docType{package}
\name{osmdata}
\alias{osmdata}
\title{osmdata class def}
\usage{
osmdata(
  bbox = NULL,
  overpass_call = NULL,
  meta = NULL,
  osm_points = NULL,
  osm_lines = NULL,
  osm_polygons = NULL,
  osm_multilines = NULL,
  osm_multipolygons = NULL
)
}
\arguments{
\item{bbox}{bounding box}

\item{overpass_call}{overpass_call}

\item{meta}{metadata of overpass query, including timestamps and version
numbers}

\item{osm_points}{OSM nodes as \pkg{sf} Simple Features Collection of points
or \pkg{sp} SpatialPointsDataFrame}

\item{osm_lines}{OSM ways \pkg{sf} Simple Features Collection of linestrings
or \pkg{sp} SpatialLinesDataFrame}

\item{osm_polygons}{OSM ways as \pkg{sf} Simple Features Collection of
polygons or \pkg{sp} SpatialPolygonsDataFrame}

\item{osm_multilines}{OSM relations as \pkg{sf} Simple Features Collection
of multilinestrings or \pkg{sp} SpatialLinesDataFrame}

\item{osm_multipolygons}{OSM relations as \pkg{sf} Simple Features
Collection of multipolygons or \pkg{sp}
SpatialPolygonsDataFrame}
}
\description{
Imports OpenStreetMap (OSM) data into R as either 'sf' or 'sp' objects.  OSM
data are extracted from the overpass API and processed with very fast C++
routines for return to R.  The package enables simple overpass queries to be
constructed without the user necessarily understanding the syntax of the
overpass query language, while retaining the ability to handle arbitrarily
complex queries. Functions are also provided to enable recursive searching
between different kinds of OSM data (for example, to find all lines which
intersect a given point).
}
\note{
Class constructor should never be used directly, and is only exported
to provide access to the print method
}
\section{Functions to Prepare Queries}{

\itemize{
\item \link{getbb}: Get bounding box for a given place name
\item \link{bbox_to_string}: Convert a named matrix or a named vector
(or an unnamed vector) return a string
\item \link{overpass_status}: Retrieve status of the overpass API
\item \link{opq}: Build an overpass query
\item \link{add_osm_feature}: Add a feature to an overpass query
\item \link{opq_string}: Convert an osmdata query to overpass API
string
}
}

\section{Functions to Get Additional OSM Information}{

\itemize{
\item \link{available_features}: List recognised features in OSM
\item \link{available_tags}: List tags associated with a feature
}
}

\section{Functions to Extract OSM Data}{

\itemize{
\item \link{osmdata_sf}: Return OSM data in \pkg{sf} format
\item \link{osmdata_sp}: Return OSM data in \pkg{sp} format
\item \link{osmdata_xml}: Return OSM data in \pkg{XML} format
}
}

\section{Functions to Search Data}{

\itemize{
\item \code{osm_points}: Extract all \code{osm_points} objects
\item \code{osm_lines}: Extract all \code{osm_lines} objects
\item \code{osm_polygons}: Extract all \code{osm_polygons} objects
\item \code{osm_multilines}: Extract all \code{osm_multilines} objects
\item \code{osm_multipolygons}: Extract all \code{osm_multipolygons} objects
}
}

\author{
Mark Padgham, Bob Rudis, Robin Lovelace, Maëlle Salmon
}
\concept{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-osmdata.R
\name{osmdata_xml}
\alias{osmdata_xml}
\title{Return an OSM Overpass query in XML format
Read an (XML format) OSM Overpass response from a string, a connection,
or a raw vector.}
\usage{
osmdata_xml(q, filename, quiet = TRUE, encoding)
}
\arguments{
\item{q}{An object of class \code{overpass_query} constructed with
\link{opq} and \link{add_osm_feature}.}

\item{filename}{If given, OSM data are saved to the named file}

\item{quiet}{suppress status messages.}

\item{encoding}{Unless otherwise specified XML documents are assumed to be
encoded as UTF-8 or UTF-16. If the document is not UTF-8/16, and lacks
an explicit encoding directive, this allows you to supply a default.}
}
\value{
An object of class \code{XML::xml_document} containing the result of the
overpass API query.
}
\description{
Return an OSM Overpass query in XML format
Read an (XML format) OSM Overpass response from a string, a connection,
or a raw vector.
}
\note{
Objects of class \code{xml_document} can be saved as \code{.xml} or
\code{.osm} files with \code{xml2::write_xml}.
}
\examples{
\dontrun{
q <- opq ("hampi india")
q <- add_osm_feature (q, key="historic", value="ruins")
osmdata_xml (q, filename="hampi.osm")
}
}
\seealso{
Other extract: 
\code{\link{osmdata_sc}()},
\code{\link{osmdata_sf}()},
\code{\link{osmdata_sp}()}
}
\concept{extract}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elevation.R
\name{osm_elevation}
\alias{osm_elevation}
\title{osm_elevation}
\usage{
osm_elevation(dat, elev_file)
}
\arguments{
\item{dat}{An \code{SC} object produced by \link{osmdata_sc}.}

\item{elev_file}{A vector of one or more character strings specifying paths
to \code{.tif} files containing global elevation data.}
}
\value{
A modified version of the input \code{dat} with an additional \code{z_} column
appended to the vertices.
}
\description{
Add elevation data to a previously-extracted OSM data set, using a
pre-downloaded global elevation file from
\url{https://srtm.csi.cgiar.org/srtmdata/}. Currently only works for
\code{SC}-class objects returned from \link{osmdata_sc}.
}
\seealso{
Other transform: 
\code{\link{osm_poly2line}()},
\code{\link{trim_osmdata}()},
\code{\link{unique_osmdata}()},
\code{\link{unname_osmdata_sf}()}
}
\concept{transform}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osm-extract.R
\name{osm_polygons}
\alias{osm_polygons}
\title{Extract all \code{osm_polygons} from an osmdata object}
\usage{
osm_polygons(dat, id)
}
\arguments{
\item{dat}{An object of class \link{osmdata}}

\item{id}{OSM identification of one or more objects for which polygons are to
be extracted}
}
\value{
An \pkg{sf} Simple Features Collection of polygons
}
\description{
If \code{id} is of a point object, \code{osm_polygons} will return all
polygons containing that point. If \code{id} is of a line or polygon object,
\code{osm_polygons} will return all polygons which intersect the given line
or polygon.
}
\examples{
\dontrun{
Extract polygons which intersect Conway Street in London
dat <- opq ("Marylebone London") \%>\% add_osm_feature (key="highway") \%>\%
    osmdata_sf ()
conway <- which (dat$osm_lines$name == "Conway Street")
id <- rownames (dat$osm_lines [conway,])
osm_polygons (dat, id)
}
}
\seealso{
Other search: 
\code{\link{osm_lines}()},
\code{\link{osm_multilines}()},
\code{\link{osm_multipolygons}()},
\code{\link{osm_points}()}
}
\concept{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osmdata-package.R
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
% Please edit documentation in R/opq.R
\name{opq_enclosing}
\alias{opq_enclosing}
\title{opq_enclosing}
\usage{
opq_enclosing(
  lon = NULL,
  lat = NULL,
  key = NULL,
  value = NULL,
  enclosing = "relation",
  timeout = 25
)
}
\arguments{
\item{lon}{Longitude of desired point}

\item{lat}{Latitude of desired point}

\item{key}{(Optional) OSM key of enclosing data}

\item{value}{(Optional) OSM value matching 'key' of enclosing data}

\item{enclosing}{Either 'relation' or 'way' for whether to return enclosing
objects of those respective types (where generally 'relation' will correspond
to multipolygon objects, and 'way' to polygon objects).}

\item{timeout}{It may be necessary to increase this value for large queries,
because the server may time out before all data are delivered.}
}
\description{
Find all features which enclose a given point, and optionally match specific
'key'-'value' pairs. This function is \emph{not} intended to be combined with
\link{add_osm_feature}, rather is only to be used in the sequence
\link{opq_enclosing} -> \link{opq_string} -> \link{osmdata_xml} (or other
extraction function). See examples for how to use.
}
\examples{
\dontrun{
# Get water body surrounding a particular point:
lat <- 54.33601
lon <- -3.07677
key <- "natural"
value <- "water"
x <- opq_enclosing (lon, lat, key, value) \%>\%
    opq_string () \%>\%
    osmdata_sf ()
}
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opq.R
\name{opq_string}
\alias{opq_string}
\alias{opq_to_string}
\title{Convert an overpass query into a text string}
\usage{
opq_string(opq)
}
\arguments{
\item{opq}{An \code{overpass_query} object}
}
\value{
Character string to be submitted to the overpass API
}
\description{
Convert an osmdata query of class opq to a character string query to
be submitted to the overpass API.
}
\examples{
\dontrun{
q <- opq ("hampi india")
opq_string (q)
}
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R
\name{available_tags}
\alias{available_tags}
\title{List tags associated with a feature}
\usage{
available_tags(feature)
}
\arguments{
\item{feature}{feature to retrieve}
}
\value{
character vector of all known tags for a feature
}
\description{
List tags associated with a feature
}
\note{
requires internet access
}
\examples{
\dontrun{
available_tags("aerialway")
}
}
\references{
\url{https://wiki.openstreetmap.org/wiki/Map_Features}
}
\seealso{
Other osminfo: 
\code{\link{available_features}()}
}
\concept{osminfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R
\name{available_features}
\alias{available_features}
\title{List recognized features in OSM}
\usage{
available_features()
}
\value{
character vector of all known features
}
\description{
List recognized features in OSM
}
\note{
requires internet access
}
\examples{
\dontrun{
available_features()
}
}
\references{
\url{https://wiki.openstreetmap.org/wiki/Map_Features}
}
\seealso{
Other osminfo: 
\code{\link{available_tags}()}
}
\concept{osminfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opq.R
\name{opq_osm_id}
\alias{opq_osm_id}
\title{Add a feature specified by OSM ID to an Overpass query}
\usage{
opq_osm_id(id = NULL, type = NULL, open_url = FALSE)
}
\arguments{
\item{id}{One or more official OSM identifiers (long-form integers), which
must be entered as either a character or \emph{numeric} value (because R does not
support long-form integers).}

\item{type}{Type of object; must be either \code{node}, \code{way}, or \code{relation}}

\item{open_url}{If \code{TRUE}, open the OSM page of the specified object in web
browser. Multiple objects (\code{id} values) will be opened in multiple pages.}
}
\value{
\link{opq} object
}
\description{
Add a feature specified by OSM ID to an Overpass query
}
\note{
Extracting elements by ID requires explicitly specifying the type of
element. Only elements of one of the three given types can be extracted in a
single query, but the results of multiple types can nevertheless be combined
with the \link{c} operation of \link{osmdata}.
}
\examples{
\dontrun{
id <- c (1489221200, 1489221321, 1489221491)
dat1 <- opq_osm_id (type = "node", id = id) \%>\%
    opq_string () \%>\%
    osmdata_sf ()
dat1$osm_points # the desired nodes
id <- c (136190595, 136190596)
dat2 <- opq_osm_id (type = "way", id = id) \%>\%
    opq_string () \%>\%
    osmdata_sf ()
dat2$osm_lines # the desired ways
dat <- c (dat1, dat2) # The node and way data combined
}
}
\references{
\url{https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#By_element_id}
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_string}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opq.R
\name{add_osm_feature}
\alias{add_osm_feature}
\title{Add a feature to an Overpass query}
\usage{
add_osm_feature(
  opq,
  key,
  value,
  key_exact = TRUE,
  value_exact = TRUE,
  match_case = TRUE,
  bbox = NULL
)
}
\arguments{
\item{opq}{An \code{overpass_query} object}

\item{key}{feature key}

\item{value}{value for feature key; can be negated with an initial
exclamation mark, \code{value = "!this"}, and can also be a vector,
\code{value = c ("this", "that")}.}

\item{key_exact}{If FALSE, \code{key} is not interpreted exactly; see
\url{https://wiki.openstreetmap.org/wiki/Overpass_API}}

\item{value_exact}{If FALSE, \code{value} is not interpreted exactly}

\item{match_case}{If FALSE, matching for both \code{key} and \code{value} is
not sensitive to case}

\item{bbox}{optional bounding box for the feature query; must be set if no
opq query bbox has been set}
}
\value{
\link{opq} object
}
\description{
Add a feature to an Overpass query
}
\note{
\code{key_exact} should generally be \code{TRUE}, because OSM uses a
reasonably well defined set of possible keys, as returned by
\link{available_features}. Setting \code{key_exact = FALSE} allows matching
of regular expressions on OSM keys, as described in Section 6.1.5 of
\url{https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL}. The actual
query submitted to the overpass API can be obtained from
\link{opq_string}.
}
\section{\code{add_osm_feature} vs \code{add_osm_features}}{

Features defined within an \link{add_osm_features} call are combined with a
logical OR.

Chained calls to either \link{add_osm_feature} or \code{add_osm_features()} combines
features from these calls in a logical AND; this is analagous to chaining
\code{dplyr::filter()} on a data frame.

\code{add_osm_features()} with only one feature is logically equivalent to
\code{add_osm_feature()}.
}

\examples{
\dontrun{
q <- opq ("portsmouth usa") \%>\%
                add_osm_feature(key = "amenity",
                                value = "restaurant") \%>\%
                add_osm_feature(key = "amenity", value = "pub")
osmdata_sf (q) # all objects that are restaurants AND pubs (there are none!)
q1 <- opq ("portsmouth usa") \%>\%
                add_osm_feature(key = "amenity",
                                value = "restaurant")
q2 <- opq ("portsmouth usa") \%>\%
                add_osm_feature(key = "amenity", value = "pub")
c (osmdata_sf (q1), osmdata_sf (q2)) # all restaurants OR pubs
# Use of negation to extract all non-primary highways
q <- opq ("portsmouth uk") \%>\%
        add_osm_feature (key = "highway", value = "!primary")
}
}
\references{
\url{https://wiki.openstreetmap.org/wiki/Map_Features}
}
\seealso{
\link{add_osm_features}

Other queries: 
\code{\link{add_osm_features}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getbb.R
\name{getbb}
\alias{getbb}
\title{Get bounding box for a given place name}
\usage{
getbb(
  place_name,
  display_name_contains = NULL,
  viewbox = NULL,
  format_out = "matrix",
  base_url = "https://nominatim.openstreetmap.org",
  featuretype = "settlement",
  limit = 10,
  key = NULL,
  silent = TRUE
)
}
\arguments{
\item{place_name}{The name of the place you're searching for}

\item{display_name_contains}{Text string to match with display_name field
returned by \url{https://wiki.openstreetmap.org/wiki/Nominatim}}

\item{viewbox}{The bounds in which you're searching}

\item{format_out}{Character string indicating output format: matrix
(default), string (see \code{\link[=bbox_to_string]{bbox_to_string()}}), data.frame (all 'hits' returned
by Nominatim), sf_polygon (for polygons that work with the sf package) or
polygon (full polygonal bounding boxes for each match).}

\item{base_url}{Base website from where data is queried}

\item{featuretype}{The type of OSM feature (settlement is default; see Note)}

\item{limit}{How many results should the API return?}

\item{key}{The API key to use for services that require it}

\item{silent}{Should the API be printed to screen? TRUE by default}
}
\value{
Defaults to a matrix in the form:
\code{
  min   max
x ...   ...
y ...   ...
}

If \code{format_out = "polygon"}, one or more two-columns matrices of polygonal
longitude-latitude points. Where multiple \code{place_name} occurrences are found
within \code{nominatim}, each item of the list of coordinates may itself contain
multiple coordinate matrices where multiple exact matches exist. If one one
exact match exists with potentially multiple polygonal boundaries (for
example, "london uk" is an exact match, but can mean either greater London or
the City of London), only the first is returned. See examples below for
illustration.
}
\description{
This function uses the free Nominatim API provided by OpenStreetMap to find
the bounding box (bb) associated with place names.
}
\details{
It was inspired by the functions
\code{bbox} from the \pkg{sp} package,
\code{bb} from the \pkg{tmaptools} package and
\code{bb_lookup} from the github package \pkg{nominatim} package,
which can be found at \url{https://github.com/hrbrmstr/nominatim}.

See \url{https://wiki.openstreetmap.org/wiki/Nominatim} for details.
}
\note{
Specific values of \code{featuretype} include "street", "city",
\url{https://wiki.openstreetmap.org/wiki/Nominatim} for details). The default
\code{featuretype = "settlement"} combines results from all intermediate
levels below "country" and above "streets". If the bounding box or polygon of
a city is desired, better results will usually be obtained with
\code{featuretype = "city"}.
}
\examples{
\dontrun{
getbb("Salzburg")
# select based on display_name, print query url
getbb("Hereford", display_name_contains = "USA", silent = FALSE)
# top 3 matches as data frame
getbb("Hereford", format_out = "data.frame", limit = 3)
# Examples of polygonal boundaries
bb <- getbb ("london uk", format_out = "polygon") # single match
dim(bb[[1]][[1]]) # matrix of longitude/latitude pairs
bb_sf = getbb("kathmandu", format_out = "sf_polygon")
# sf:::plot.sf(bb_sf) # can be plotted if sf is installed
getbb("london", format_out = "sf_polygon")
getbb("accra", format_out = "sf_polygon") # rectangular bb
# Using an alternative service (locationiq requires an API key)
# add LOCATIONIQ=type_your_api_key_here to .Renviron:
key <- Sys.getenv("LOCATIONIQ")
if(nchar(key) ==  32) {
  getbb(place_name,
        base_url = "https://locationiq.org/v1/search.php",
        key = key)
}
}
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{set_overpass_url}
\alias{set_overpass_url}
\title{set_overpass_url}
\usage{
set_overpass_url(overpass_url)
}
\arguments{
\item{overpass_url}{The desired overpass API URL}
}
\value{
The overpass API URL
}
\description{
Set the URL of the specified overpass API. Possible APIs with global coverage
are:
\itemize{
\item "https://overpass-api.de/api/interpreter" (default)
\item "https://overpass.kumi.systems/api/interpreter"
\item "https://overpass.osm.rambler.ru/cgi/interpreter"
\item "https://api.openstreetmap.fr/oapi/interpreter"
\item "https://overpass.osm.vi-di.fr/api/interpreter"
}
Additional APIs with limited local coverage include:
\itemize{
\item "https://overpass.osm.ch/api/interpreter" (Switzerland)
\item "https://overpass.openstreetmap.ie/api/interpreter" (Ireland)
}
}
\details{
For further details, see
\url{https://wiki.openstreetmap.org/wiki/Overpass_API}
}
\seealso{
\code{\link[=get_overpass_url]{get_overpass_url()}}

Other overpass: 
\code{\link{get_overpass_url}()}
}
\concept{overpass}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overpass-query.R
\name{overpass_status}
\alias{overpass_status}
\title{Retrieve status of the Overpass API}
\usage{
overpass_status(quiet = FALSE)
}
\arguments{
\item{quiet}{if \code{FALSE} display a status message}
}
\value{
an invisible list of whether the API is available along with the
text of the message from Overpass and the timestamp of the
next available slot
}
\description{
Retrieve status of the Overpass API
}
\seealso{
Other queries: 
\code{\link{add_osm_features}()},
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{opq}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opq.R
\name{add_osm_features}
\alias{add_osm_features}
\title{Add multiple features to an Overpass query}
\usage{
add_osm_features(opq, features, bbox = NULL)
}
\arguments{
\item{opq}{An \code{overpass_query} object}

\item{features}{Character vector of key-value pairs with keys and values
enclosed in escape-formatted quotations (see examples).}

\item{bbox}{optional bounding box for the feature query; must be set if no
opq query bbox has been set}
}
\value{
\link{opq} object
}
\description{
Alternative version of \link{add_osm_feature} for creating single queries
with multiple features. Key-value matching may be controlled by using
the filter symbols described in
\url{https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#By_tag_.28has-kv.29}.
}
\section{\code{add_osm_feature} vs \code{add_osm_features}}{

Features defined within an \link{add_osm_features} call are combined with a
logical OR.

Chained calls to either \link{add_osm_feature} or \code{add_osm_features()} combines
features from these calls in a logical AND; this is analagous to chaining
\code{dplyr::filter()} on a data frame.

\code{add_osm_features()} with only one feature is logically equivalent to
\code{add_osm_feature()}.
}

\examples{
\dontrun{
q <- opq ("portsmouth usa") \%>\%
     add_osm_features (features = c ("\"amenity\"=\"restaurant\"",
                                     "\"amenity\"=\"pub\""))
# This extracts in a single query the same result as the following:
q1 <- opq ("portsmouth usa") \%>\%
                add_osm_feature(key = "amenity",
                                value = "restaurant")
q2 <- opq ("portsmouth usa") \%>\%
                add_osm_feature(key = "amenity", value = "pub")
c (osmdata_sf (q1), osmdata_sf (q2)) # all restaurants OR pubs
}
}
\references{
\url{https://wiki.openstreetmap.org/wiki/Map_Features}
}
\seealso{
\link{add_osm_feature}

Other queries: 
\code{\link{add_osm_feature}()},
\code{\link{bbox_to_string}()},
\code{\link{getbb}()},
\code{\link{opq_around}()},
\code{\link{opq_enclosing}()},
\code{\link{opq_osm_id}()},
\code{\link{opq_string}()},
\code{\link{opq}()},
\code{\link{overpass_status}()}
}
\concept{queries}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{unname_osmdata_sf}
\alias{unname_osmdata_sf}
\title{unname_osmdata_sf}
\usage{
unname_osmdata_sf(x)
}
\arguments{
\item{x}{An 'osmdata_sf' object returned from function of same name}
}
\value{
Same object, yet with no row names on geometry objects.
}
\description{
Remove names from \code{osmdata} geometry objects, for cases in which these cause
issues, particularly with plotting, such as
\url{https://github.com/rstudio/leaflet/issues/631}, or
\url{https://github.com/r-spatial/sf/issues/1177}. Note that removing these
names also removes any ability to inter-relate the different components of an
\code{osmdata} object, so use of this function is only recommended to resolve
issues such as those linked to above.
}
\examples{
\dontrun{
hampi_sf <- opq ("hampi india") \%>\%
            add_osm_feature (key="historic", value="ruins") \%>\%
            osmdata_sf ()
hampi_clean <- unname_osmdata_sf (hampi_sf)

# All coordinate matrices include rownames with OSM ID values:
head (as.matrix (hampi_sf$osm_lines$geometry [[1]]))
# But 'unname_osmdata_sf' removes both row and column names:
head (as.matrix (hampi_clean$osm_lines$geometry [[1]]))
}
}
\seealso{
Other transform: 
\code{\link{osm_elevation}()},
\code{\link{osm_poly2line}()},
\code{\link{trim_osmdata}()},
\code{\link{unique_osmdata}()}
}
\concept{transform}
