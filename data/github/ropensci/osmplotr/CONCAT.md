<!-- README.md is generated from README.Rmd. Please edit that file -->

[![R build
status](https://github.com/ropensci/osmplotr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/osmplotr/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/osmplotr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/osmplotr)
[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)

![](man/figures/map1.png)

[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/osmplotr?color=orange)](http://cran.r-project.org/package=osmplotr/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/osmplotr)](http://cran.r-project.org/package=osmplotr/)
[![](https://badges.ropensci.org/27_status.svg)](https://github.com/ropensci/software-review/issues/27)

R package to produce visually impressive customisable images of
OpenStreetMap (OSM) data downloaded internally from the [overpass
api](http://overpass-api.de/). The above map was produced directly from
`osmplotr` with no further modification. This `README` briefly
demonstrates the following functionality:

[1. Quick Introduction](#1%20intro)

[2. Installation](#2%20installation)

[3. A Simple Map](#3%20simple%20map)

[4. Highlighting Selected Areas](#4%20highlighting%20areas)

[5. Highlighting Clusters](#5%20highlighting%20clusters)

[6. Highlighting Areas Bounded by Named
Highways](#6%20highlighting%20with%20highways)

[7. Data Surfaces](#7%20data%20surfaces)

[8. Gallery](#8%20gallery)

------------------------------------------------------------------------

## <a name="1 intro"></a>1. Quick Introduction

But first the easy steps to map making:

1.  Specify the bounding box for the desired region

    ``` r
    bbox <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
    ```

2.  Download the desired data—in this case, all building perimeters.

    ``` r
    dat_B <- extract_osm_objects (key = "building", bbox = bbox)
    ```

3.  Initiate an `osm_basemap` with desired background (`bg`) colour

    ``` r
    map <- osm_basemap (bbox = bbox, bg = "gray20")
    ```

4.  Overlay objects on plot in the desired colour.

    ``` r
    map <- add_osm_objects (map, dat_B, col = "gray40")
    ```

5.  Print the map to graphics device of choice

    ``` r
    print_osm_map (map)
    ```

------------------------------------------------------------------------

## <a name="2 installation"></a>2. Installation

First install the package

``` r
install.packages ("osmplotr")
```

or the development version

``` r
devtools::install_github ("ropensci/osmplotr")
```

And then load it in the usual way

``` r
library (osmplotr)
```

------------------------------------------------------------------------

## <a name="3 simple map"></a>3. A Simple Map

Simple maps can be made by overlaying different kinds of OSM data in
different colours:

``` r
dat_H <- extract_osm_objects (key = "highway", bbox = bbox)
dat_P <- extract_osm_objects (key = "park", bbox = bbox)
dat_G <- extract_osm_objects (key = "landuse", value = "grass", bbox = bbox)
```

``` r
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_objects (map, dat_B, col = "gray40")
map <- add_osm_objects (map, dat_H, col = "gray80")
map <- add_osm_objects (map, dat_P, col = "darkseagreen")
map <- add_osm_objects (map, dat_G, col = "darkseagreen1")
print_osm_map (map)
```

<!--
![](./man/figures/map2.png)
-->

<img src="man/figures/map2.png" width = "80%"/>

------------------------------------------------------------------------

## <a name="4 highlighting areas"></a>4. Highlighting Selected Areas

`osmplotr` is primarily intended as a data visualisation tool,
particularly through enabling selected regions to be highlighted.
Regions can be defined according to simple point boundaries:

``` r
pts <- sp::SpatialPoints (cbind (c (-0.115, -0.13, -0.13, -0.115),
                             c (51.505, 51.505, 51.515, 51.515)))
```

OSM objects within the defined regions can then be highlighted with
different colour schemes. `cols` defines colours for each group (with
only one here), while `bg` defines the colour of the remaining,
background area.

``` r
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_groups (map, dat_B, groups = pts, cols = "orange", bg = "gray40")
map <- add_osm_objects (map, london$dat_P, col = "darkseagreen1")
map <- add_osm_groups (map, london$dat_P, groups = pts, cols = "darkseagreen1",
                   bg = "darkseagreen", boundary = 0)
print_osm_map (map)
```

<!--
![](./man/figures/map3.png)
-->

<img src="man/figures/map3.png" width = "80%"/>

Note the `border = 0` argument on the last call divides the park
polygons precisely along the border. The same map highlighted in
dark-on-light:

``` r
map <- osm_basemap (bbox = bbox, bg = "gray95")
map <- add_osm_groups (map, dat_B, groups = pts, cols = "gray40", bg = "gray85")
map <- add_osm_groups (map, dat_H, groups = pts, cols = "gray20", bg = "gray70")
print_osm_map (map)
```

<!--
![](./man/figures/map4.png)
-->

<img src="man/figures/map4.png" width = "80%"/>

------------------------------------------------------------------------

## <a name="5 highlighting clusters"></a>5. Highlighting Clusters

`add_osm_groups` also enables plotting an entire region as a group of
spatially distinct clusters of defined colours. Groups can be defined by
simple spatial points denoting their centres:

``` r
set.seed (2)
ngroups <- 12
x <- bbox [1, 1] + runif (ngroups) * diff (bbox [1, ])
y <- bbox [2, 1] + runif (ngroups) * diff (bbox [2, ])
groups <- cbind (x, y)
groups <- apply (groups, 1, function (i)
              sp::SpatialPoints (matrix (i, nrow = 1, ncol = 2)))
```

Calling `add_osm_groups` with no `bg` argument forces all points lying
outside those defined groups to be allocated to the nearest groups, and
thus produces an inclusive grouping extending across an entire region.

``` r
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_groups (map, dat_B, groups = groups,
                       cols = rainbow (length (groups)), border_width = 2)
print_osm_map (map)
```

<!--
![](./man/figures/map5.png)
-->

<img src="man/figures/map5.png" width = "80%"/>

------------------------------------------------------------------------

## <a name="6 highlighting with highways"></a>6. Highlighting Areas Bounded by Named Highways

An alternative way of defining highlighted groups is by naming the
highways encircling desired regions.

``` r
# These highways extend beyond the previous, smaller bbox
bbox_big <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
highways <- c ("Davies.St", "Berkeley.Sq", "Berkeley.St", "Piccadilly",
               "Regent.St", "Oxford.St")
highways1 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Regent.St", "Oxford.St", "Shaftesbury")
highways2 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Piccadilly", "Shaftesbury.Ave", "Charing.Cross.R",
               "Saint.Martin", "Trafalgar.Sq", "Cockspur.St",
               "Pall.Mall", "St.James")
highways3 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Charing.Cross", "Duncannon.St", "Strand", "Aldwych",
               "Kingsway", "High.Holborn", "Shaftesbury.Ave")
highways4 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Kingsway", "Holborn", "Farringdon.St", "Strand",
               "Fleet.St", "Aldwych")
highways5 <- connect_highways (highways = highways, bbox = bbox_big)
groups <- list (highways1, highways2, highways3, highways4, highways5)
```

And then passing these lists of groups returned by `connect_highways` to
`add_osm_groups`, this time with some Wes Anderson flair.

``` r
map <- osm_basemap (bbox = bbox, bg = "gray20")
library (wesanderson)
cols <- wes_palette ("Darjeeling", 5)
map <- add_osm_groups (map, dat_B, groups = groups, boundary = 1,
                       cols = cols, bg = "gray40", colmat = FALSE)
map <- add_osm_groups (map, dat_H, groups = groups, boundary = 0,
                       cols = cols, bg = "gray70", colmat = FALSE)
print_osm_map (map)
```

<!--
![](./man/figures/map6.png)
-->

<img src="man/figures/map6.png" width = "80%"/>

------------------------------------------------------------------------

## <a name="7 data surfaces"></a>7. Data Surfaces

Finally, `osmplotr` contains a function `add_osm_surface` that spatially
interpolates a given set of spatial data points and colours OSM objects
according to a specified colour gradient. This is illustrated here with
the `volcano` data projected onto the `bbox`.

``` r
x <- seq (bbox [1, 1], bbox [1, 2], length.out = dim (volcano)[1])
y <- seq (bbox [2, 1], bbox [2, 2], length.out = dim (volcano)[2])
xy <- cbind (rep (x, dim (volcano) [2]), rep (y, each = dim (volcano) [1]))
z <- as.numeric (volcano)
dat <- data.frame (x = xy [, 1], y = xy [, 2], z = z)
```

``` r
map <- osm_basemap (bbox = bbox, bg = "gray20")
cols <- gray (0:50 / 50)
map <- add_osm_surface (map, dat_B, dat = dat, cols = cols)
# Darken cols by ~20%
map <- add_osm_surface (map, dat_H, dat = dat,
                        cols = adjust_colours (cols, -0.2))
map <- add_colourbar (map, cols = cols, zlims = range (volcano))
map <- add_axes (map)
print_osm_map (map)
```

<!--
![](./man/figures/map7.png)
-->

<img src="man/figures/map7.png" width = "80%"/>

------------------------------------------------------------------------

## <a name="8 gallery"></a>8. Gallery

Got a nice `osmplotr` map? Please contribute in one of the following
ways:

1.  Fork repo, add link to `README.md/.Rmd`, and send pull request; or

2.  Open issue with details; or

3.  Send email to address in
    [`DESCRIPTION`](https://github.com/ropensci/osmplotr/blob/master/DESCRIPTION).

------------------------------------------------------------------------

See package vignettes ([basic
maps](https://docs.ropensci.org/osmplotr/articles/basic-maps.html) and
[data maps](https://docs.ropensci.org/osmplotr/articles/data-maps.html))
for a lot more detail and further capabilities of `osmplotr`. Please
note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree
to abide by its terms.

------------------------------------------------------------------------

[![ropensci\_footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
osmplotr v0.3.3
===============

Minor changes
-------
- Changes in response to `spatstat` v2 updates

Minor changes
-------
- 'add_osm_surface' functions changed to directly calculate and plot colours of
  objects, rather than rely on `ggplot2::scale_fill_gradientn`.

osmplotr v0.3.2
===============

Minor changes
-------
* 'verbose' parameter of 'extract_osm_objects' renamed to 'quiet'

osmplotr v0.3.1
===============

Major changes
-------
* New function 'osm_line2poly' enables plotting polygonal shapes delineated only
  by lines, through tracing around the bounding box to form full polygons.
* New vignette to describe this functionality, "maps-with-ocean".

Minor changes
-------
* 'osm_basemap' now accepts an 'sf' object instead of explicit 'bbox' values,
  and extracts the corresponding 'bbox' directly from that object.


osmplotr v0.3.0
===============

Major changes
-------
* Major re-structure to use 'osmdata' package instead of 'osmar', with
  concomitantly enormous increase in speed of 'extract_osm_objects'
* Package is now also 'sf'-compatible: objects to be plotted can be either 'sp'
  or 'sf' format, with all 'osmplotr' functions defaulting to 'sf'

Minor changes
-------
* Title in DESCRIPTION changed from "Customisable Images of OpenStreetMap Data"
  to "Bespoke Images of 'OpenStreetMap' Data"
* Better control of timeout errors when calling the overpass API
* Git host transferred from ropenscilabs to ropensci
* Acknowledge OSM contributors on startup
* Rename 'borderWidth' parameter of 'add_osm_groups' to 'border_width'
* 'connect_highways' also entirely re-coded to be much more efficient, but
  this should not affect functionality at all.

osmplotr v0.2.3
===============

* add 'return_type' argument to 'extract_osm_objects' to enable explicit
  specification of return type (points, lines, polygons)
* fix tests so they pass even if download fails

osmplotr v0.2.2
===============

* 'add_osm_surface' did not previously work properly for different bboxes
    (and so zooming was not possible). Now fixed.
* both 'add_osm_surface' and 'add_osm_groups' now enable maps to be zoomed
* fix make_osm_map to produce maps even when not all requested data exists

osmplotr v0.2.1
===============

* vignette 'making-maps' renamed 'basic-maps' and tidied
* vignette 'making-maps-with-data' renamed 'data-maps' and tidied
* 'plot_osm_basemap' renamed 'osm_basemap', and now uses
    'ggplot2::coord_equal()' to ensure maps are scaled to bounding boxes.
* 'print_osm_map' added to enable device proportions to be automatically scaled
    to bounding boxes.
* manual entries cleaned up to remove non-exported functions

osmplotr v0.2.0
===============

Major update with (almost) all plotting routines shifted from 'graphics::plot'
to 'ggplot2'. All previous parameters specifying graphics devices (such as
heights and widths) no longer apply.

Changes:

* vignette 'downloading-data' removed (incorporate in 'making-maps')
* vignette 'making-maps' extended
* vignette 'making-maps-with-data' added
* Extensive examples added to most functions
* 'click_map' removed
* 'connect_highways' renamed 'get_highway_cycle'
* 'highways2polygon' renamed 'connect_highways'
* 'extract_highway', 'extract_highways', 'order_lines' no longer exported
* 'extract_osm_objects' now just returns objects (instead of '$obj' and
    '$warn'), and dumps warnings direct to screen.
* 'add_osm_groups' now accepts lists of simple spatial points as groups
* Coordinate reference system properly attributed to all objects
* many tests added
* Change to 'ggplot2' has considerably changed structure of many functions. For
  details see function examples and vignettes


osmplotr v0.1-3
===============

Changes:

* added 'add_axes' to plot lat-lon axes
* added 'add_osm_surface' to spatially interpolate continuous surfaces from
user-defined data
* added 'add_colourbar' to plot a colourbar legend for 'add_osm_surface'
* renamed 'group_osm_objects' to 'add_osm_groups' 
* added 'adjust_colours' to allow colours to be lightened or darkened
* all usages of 'xylims' (vectors of four components) and 'get_xylims' changed
to 'bbox' (2-by-2 matrices) for consistency with sp and tmap
* reduce size of 'london' data (through smaller bbox), with corresponding
changes in vignettes

osmplotr v0.1-1, 0.1-2
=================

Published on CRAN
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
# Contributing to osmplotr

## Opening issues

The easiest way to note any behavioural curiosities or to request any new
features is by opening a [github issue](https://github.com/ropensci/osmplotr/issues).


## Development guidelines

If you'd like to contribute changes to `osmplotr`, we use [the GitHub
flow](https://guides.github.com/introduction/flow/index.html) for proposing,
submitting, reviewing, and accepting changes. If you haven't done this before,
there's a nice overview of git [here](http://r-pkgs.had.co.nz/git.html), as well
as best practices for submitting pull requests
[here](http://r-pkgs.had.co.nz/git.html#pr-make).

The `osmplotr` coding style diverges somewhat from [this commonly used R style
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
permeates `osmplotr` code. Words of text are separated by whitespace, and so
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
conduct](https://github.com/ropensci/osmplotr/blob/master/CODE_OF_CONDUCT.md)
for more information.
# CRAN notes for osmplotr_v0.3.3 submission

This submission rectifies the previously failing checks due to the breakup of the `spatstat` package into new sub-packages.

The single note regarding installed size of ~6MB is due to the vignettes. These produce many graphical files illustrating the package's functionality. Halving the current resolution of these images (from 72 to 36 dpi) only decreases the final package size by around 200 kB.

## Test environments

Other than the above, this submission generates NO notes on:

- Linux (via github actions): R-release, R-oldrelease
- Windows (via github actions): R-release, R-oldrelease, R-devel
- win-builder: R-oldrelease, R-release, R-devel
Compiled vignettes from the cran version of 'osmplotr' are here:

* [making maps](https://cran.r-project.org/web/packages/osmplotr/vignettes/making-maps.html)

* [making maps with data](https://cran.r-project.org/web/packages/osmplotr/vignettes/making-maps-with-data.html)

---
title: "osmplotr, an R package for making maps with OpenStreetMap data"
keywords: "open street map, openstreetmap, OSM, map, visualisation, visualization"
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

[![R build
status](https://github.com/ropensci/osmplotr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/osmplotr/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/osmplotr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/osmplotr)
[![Project Status: Active](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)

![](man/figures/map1.png)

[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/osmplotr?color=orange)](http://cran.r-project.org/package=osmplotr/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/osmplotr)](http://cran.r-project.org/package=osmplotr/)
[![](https://badges.ropensci.org/27_status.svg)](https://github.com/ropensci/software-review/issues/27)

R package to produce visually impressive customisable images of OpenStreetMap
(OSM) data downloaded internally from the 
[overpass api](http://overpass-api.de/). The above map was produced directly
from `osmplotr` with no further modification. This `README` briefly demonstrates
the following functionality:

[1. Quick Introduction](#1 intro)

[2. Installation](#2 installation)

[3. A Simple Map](#3 simple map)

[4. Highlighting Selected Areas](#4 highlighting areas)

[5. Highlighting Clusters](#5 highlighting clusters)

[6. Highlighting Areas Bounded by Named Highways](#6 highlighting with highways)

[7. Data Surfaces](#7 data surfaces)

[8. Gallery](#8 gallery)

---------------

## <a name="1 intro"></a>1. Quick Introduction

But first the easy steps to map making:
```{r, echo = FALSE, message = FALSE, eval = TRUE}
library (osmplotr)
```


1. Specify the bounding box for the desired region
    ```{r}
    bbox <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
    ```
2. Download the desired data---in this case, all building perimeters.
    ```{r, eval = FALSE}
    dat_B <- extract_osm_objects (key = "building", bbox = bbox)
    ```
3. Initiate an `osm_basemap` with desired background (`bg`) colour
    ```{r map1, eval = FALSE}
    map <- osm_basemap (bbox = bbox, bg = "gray20")
    ```
4. Overlay objects on plot in the desired colour.
    ```{r, eval = FALSE}
    map <- add_osm_objects (map, dat_B, col = "gray40")
    ```
5. Print the map to graphics device of choice
    ```{r, eval = FALSE}
    print_osm_map (map)
    ```

```{r london2, echo = FALSE, eval = FALSE}
library (osmdata)
bbox <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
q0 <- opq (bbox)
q1 <- add_osm_feature (q0, key = "building")
dat_B <- osmdata_sf (q1, quiet = FALSE)$osm_polygons
q1 <- add_osm_feature (q0, key = "highway")
dat_H <- osmdata_sf (q1, quiet = FALSE)$osm_lines
q1 <- add_osm_feature (q0, key = "leisure", value = "park")
dat_P <- osmdata_sf (q1, quiet = FALSE)$osm_polygons
q1 <- add_osm_feature (q0, key = "landuse", value = "grass")
dat_G <- osmdata_sf (q1, quiet = FALSE)$osm_polygons
london2 <- list (dat_B = dat_B, dat_H = dat_H, dat_P = dat_P, dat_G = dat_G)
save (london2, file = "london2.rda")
```

---------------

## <a name="2 installation"></a>2. Installation

First install the package
```{r, eval = FALSE}
install.packages ("osmplotr")
```
or the development version
```{r, eval = FALSE}
devtools::install_github ("ropensci/osmplotr")
```
And then load it in the usual way
```{r, eval = FALSE}
library (osmplotr)
```

---------------

## <a name="3 simple map"></a>3. A Simple Map

Simple maps can be made by overlaying different kinds of OSM data in different
colours:
```{r, eval = FALSE}
dat_H <- extract_osm_objects (key = "highway", bbox = bbox)
dat_P <- extract_osm_objects (key = "park", bbox = bbox)
dat_G <- extract_osm_objects (key = "landuse", value = "grass", bbox = bbox)
```
```{r, echo = FALSE, eval = FALSE}
load ("london2.rda")
dat_B <- london2$dat_B
dat_H <- london2$dat_H
dat_G <- london2$dat_G
dat_P <- london2$dat_P
```
```{r map2, eval = FALSE}
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_objects (map, dat_B, col = "gray40")
map <- add_osm_objects (map, dat_H, col = "gray80")
map <- add_osm_objects (map, dat_P, col = "darkseagreen")
map <- add_osm_objects (map, dat_G, col = "darkseagreen1")
print_osm_map (map)
```
```{r map2-print, eval = FALSE, echo = FALSE}
print_osm_map (map, file = "map2.png", width = 600, units = "px", dpi = 72)
```

<!--
![](./man/figures/map2.png)
-->
<img src="man/figures/map2.png" width = "80%"/>

---------------

## <a name="4 highlighting areas"></a>4. Highlighting Selected Areas

`osmplotr` is primarily intended as a data visualisation tool, particularly
through enabling selected regions to be highlighted. Regions can be defined
according to simple point boundaries:

```{r}
pts <- sp::SpatialPoints (cbind (c (-0.115, -0.13, -0.13, -0.115),
                             c (51.505, 51.505, 51.515, 51.515)))
```
OSM objects within the defined regions can then be highlighted with different
colour schemes. `cols` defines colours for each group (with only one here),
while `bg` defines the colour of the remaining, background area.
```{r map3, eval = FALSE}
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_groups (map, dat_B, groups = pts, cols = "orange", bg = "gray40")
map <- add_osm_objects (map, london$dat_P, col = "darkseagreen1")
map <- add_osm_groups (map, london$dat_P, groups = pts, cols = "darkseagreen1",
                   bg = "darkseagreen", boundary = 0)
print_osm_map (map)
```
```{r map3-print, eval = FALSE, echo = FALSE}
print_osm_map (map, filename = "map3.png", width = 600, units = "px", dpi = 72)
```
<!--
![](./man/figures/map3.png)
-->
<img src="man/figures/map3.png" width = "80%"/>

Note the `border = 0` argument on the last call divides the park polygons
precisely along the border. The same map highlighted in dark-on-light:
```{r map4, eval = FALSE}
map <- osm_basemap (bbox = bbox, bg = "gray95")
map <- add_osm_groups (map, dat_B, groups = pts, cols = "gray40", bg = "gray85")
map <- add_osm_groups (map, dat_H, groups = pts, cols = "gray20", bg = "gray70")
print_osm_map (map)
```
```{r map4-print, eval = FALSE, echo = FALSE}
print_osm_map (map, filename = "map4.png", width = 600, units = "px", dpi = 72)
```
<!--
![](./man/figures/map4.png)
-->
<img src="man/figures/map4.png" width = "80%"/>

---------------

## <a name="5 highlighting clusters"></a>5. Highlighting Clusters

`add_osm_groups` also enables plotting an entire region as a group of
spatially distinct clusters of defined colours. Groups can be defined by simple
spatial points denoting their centres:

```{r, echo = TRUE}
set.seed (2)
ngroups <- 12
x <- bbox [1, 1] + runif (ngroups) * diff (bbox [1, ])
y <- bbox [2, 1] + runif (ngroups) * diff (bbox [2, ])
groups <- cbind (x, y)
groups <- apply (groups, 1, function (i)
              sp::SpatialPoints (matrix (i, nrow = 1, ncol = 2)))
```
Calling `add_osm_groups` with no `bg` argument forces all points lying outside
those defined groups to be allocated to the nearest groups, and thus produces an
inclusive grouping extending across an entire region.

```{r map5, eval = FALSE}
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_groups (map, dat_B, groups = groups,
                       cols = rainbow (length (groups)), border_width = 2)
print_osm_map (map)
```
```{r map5-print, eval = FALSE, echo = FALSE}
print_osm_map (map, filename = "map5.png", width = 600, units = "px", dpi = 72)
```
<!--
![](./man/figures/map5.png)
-->
<img src="man/figures/map5.png" width = "80%"/>

---------------

## <a name="6 highlighting with highways"></a>6. Highlighting Areas Bounded by Named Highways

An alternative way of defining highlighted groups is by naming the highways
encircling desired regions.
```{r, eval = FALSE}
# These highways extend beyond the previous, smaller bbox
bbox_big <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
highways <- c ("Davies.St", "Berkeley.Sq", "Berkeley.St", "Piccadilly",
               "Regent.St", "Oxford.St")
highways1 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Regent.St", "Oxford.St", "Shaftesbury")
highways2 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Piccadilly", "Shaftesbury.Ave", "Charing.Cross.R",
               "Saint.Martin", "Trafalgar.Sq", "Cockspur.St",
               "Pall.Mall", "St.James")
highways3 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Charing.Cross", "Duncannon.St", "Strand", "Aldwych",
               "Kingsway", "High.Holborn", "Shaftesbury.Ave")
highways4 <- connect_highways (highways = highways, bbox = bbox_big)
highways <- c ("Kingsway", "Holborn", "Farringdon.St", "Strand",
               "Fleet.St", "Aldwych")
highways5 <- connect_highways (highways = highways, bbox = bbox_big)
groups <- list (highways1, highways2, highways3, highways4, highways5)
```
And then passing these lists of groups returned by `connect_highways` to
`add_osm_groups`, this time with some Wes Anderson flair.
```{r map 6, eval = FALSE}
map <- osm_basemap (bbox = bbox, bg = "gray20")
library (wesanderson)
cols <- wes_palette ("Darjeeling", 5)
map <- add_osm_groups (map, dat_B, groups = groups, boundary = 1,
                       cols = cols, bg = "gray40", colmat = FALSE)
map <- add_osm_groups (map, dat_H, groups = groups, boundary = 0,
                       cols = cols, bg = "gray70", colmat = FALSE)
print_osm_map (map)
```
```{r map6-print, eval = FALSE, echo = FALSE}
print_osm_map (map, filename = "map6.png", width = 600, units = "px", dpi = 72)
```
<!--
![](./man/figures/map6.png)
-->
<img src="man/figures/map6.png" width = "80%"/>


---------------

## <a name="7 data surfaces"></a>7. Data Surfaces

Finally, `osmplotr` contains a function `add_osm_surface` that spatially
interpolates a given set of spatial data points and colours OSM objects
according to a specified colour gradient. This is illustrated here with the
`volcano` data projected onto the `bbox`.
```{r}
x <- seq (bbox [1, 1], bbox [1, 2], length.out = dim (volcano)[1])
y <- seq (bbox [2, 1], bbox [2, 2], length.out = dim (volcano)[2])
xy <- cbind (rep (x, dim (volcano) [2]), rep (y, each = dim (volcano) [1]))
z <- as.numeric (volcano)
dat <- data.frame (x = xy [, 1], y = xy [, 2], z = z)
```
```{r map7, eval = FALSE}
map <- osm_basemap (bbox = bbox, bg = "gray20")
cols <- gray (0:50 / 50)
map <- add_osm_surface (map, dat_B, dat = dat, cols = cols)
# Darken cols by ~20%
map <- add_osm_surface (map, dat_H, dat = dat,
                        cols = adjust_colours (cols, -0.2))
map <- add_colourbar (map, cols = cols, zlims = range (volcano))
map <- add_axes (map)
print_osm_map (map)
```
```{r map7-print, eval = FALSE, echo = FALSE}
print_osm_map (map, filename = "map7.png", width = 600, units = "px", dpi = 72)
```

```{r map1-print, eval = FALSE, echo = FALSE}
# This is map1 used as the title
# extrafont::loadfonts ()
lab_dat <- data.frame (x = mean (bbox[1, ]), y = mean (bbox [2, ]),
                       lab = "osmplotr")
aes <- ggplot2::aes (x, y, label = lab)

bbox <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
map <- osm_basemap (bbox = bbox, bg = "gray20")
cols <- gray (0:50 / 50)
map <- add_osm_surface (map, dat_B, dat = dat, cols = cols)
map <- add_osm_surface (map, dat_H, dat = dat,
                        cols = adjust_colours (cols, -0.2))

#map2 <- map + ggplot2::geom_text (dat = dat, mapping = aes, size = 60,
#                                  colour = "white",
#                                 family = "Lato Light", nudge_y = 0.0015)
map2 <- map + ggplot2::geom_text (dat = lab_dat, mapping = aes, size = 45,
                                  colour = "black",
                                 family = "Purisa", fontface = 2,
                                 nudge_y = 0.0005, nudge_x = 0.0005)
map2 <- map2 + ggplot2::geom_text (dat = lab_dat, mapping = aes, size = 45,
                                   colour = "white", family = "Purisa",
                                   nudge_y = 0.001, fontface = 2)
print_osm_map (map2, filename = "map1.png", width = 800, units = "px",
               dpi = 72)
```
<!--
![](./man/figures/map7.png)
-->
<img src="man/figures/map7.png" width = "80%"/>

---------------

## <a name="8 gallery"></a>8. Gallery

Got a nice `osmplotr` map? Please contribute in one of the following ways:

1. Fork repo, add link to `README.md/.Rmd`, and send pull request; or

2. Open issue with details; or

3. Send email to address in
   [`DESCRIPTION`](https://github.com/ropensci/osmplotr/blob/master/DESCRIPTION).

---------------

See package vignettes 
([basic maps](https://docs.ropensci.org/osmplotr/articles/basic-maps.html) and
[data maps](https://docs.ropensci.org/osmplotr/articles/data-maps.html)) for a
lot more detail and further capabilities of `osmplotr`.  Please note that this
project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By
participating in this project you agree to abide by its terms.

--------------

[![ropensci\_footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
# scrips to generate both `data/london.rda` and `inst/extdata/hwys.rda`

```{r load, echo = FALSE}
devtools::load_all (".", export_all = FALSE)
library (magrittr)
```


### test data for highway cycles

These are the `/inst/extdata/hwys.rda` data

```{r}
bbox <- get_bbox (c(-0.15, 51.50, -0.10, 51.52))
highways1 <- c ('Monmouth.St', 'Short.?s.Gardens', 'Endell.St', 'Long.Acre',
               'Upper.Saint.Martin') %>%
                osmplotr:::extract_highways (bbox = bbox)
highways2 <- c ('Endell.St', 'High.Holborn', 'Drury.Lane', 'Long.Acre') %>%
                osmplotr:::extract_highways (bbox = bbox)
highways3 <- c ('Drury.Lane', 'High.Holborn', 'Kingsway', 'Great.Queen.St') %>%
                osmplotr:::extract_highways (bbox = bbox)
highways4 <- c ('Kingsway', 'Holborn', 'Farringdon.St', 'Strand',
               'Fleet.St', 'Aldwych') %>%
                osmplotr:::extract_highways (bbox = bbox)
hwys <- list (highways1 = highways1, highways2 = highways2,
              highways3 = highways3, highways4 = highways4)
fname <- system.file ('extdata', 'hwys.rda', package = 'osmplotr')
save (hwys, file = fname)
format (file.size (fname), big.mark = ',')
```


### The main London data


The `london` data are stripped of all columns except the 2 primary ones. The
names can't be stored because they fail `R CMD check` due to non-ASCII strings.

```{r}
col_names <- c ('osm_id', 'geometry')
bbox <- get_bbox (c (-0.13, 51.51, -0.11, 51.52))
dat_H <- extract_osm_objects (key = 'highway', value = '!primary',
                              bbox = bbox)
indx <- which (names (dat_H) %in% col_names)
dat_H <- dat_H [, indx]
dat_HP <- extract_osm_objects (key = 'highway', value = 'primary',
                              bbox = bbox)
indx <- which (names (dat_HP) %in% col_names)
dat_HP <- dat_HP [, indx]
dat_BNR <- extract_osm_objects (key = 'building', value = '!residential',
                              bbox = bbox)
indx <- which (names (dat_BNR) %in% col_names)
dat_BNR <- dat_BNR [, indx]
dat_BR <- extract_osm_objects (key = 'building', value = 'residential',
                              bbox = bbox)
indx <- which (names (dat_BR) %in% col_names)
dat_BR <- dat_BR [, indx]
dat_BC <- extract_osm_objects (key = 'building', value = 'commercial',
                              bbox = bbox)
indx <- which (names (dat_BC) %in% col_names)
dat_BC <- dat_BC [, indx]
dat_A <- extract_osm_objects (key = 'amenity', bbox = bbox,
                              return_type = 'polygon')
indx <- which (names (dat_A) %in% col_names)
dat_A <- dat_A [, indx]
dat_P <- extract_osm_objects (key = 'park', bbox = bbox)
indx <- which (names (dat_P) %in% col_names)
dat_P <- dat_P [, indx]
dat_T <- extract_osm_objects (key = 'tree', bbox = bbox)
indx <- which (names (dat_T) %in% col_names)
dat_T <- dat_T [, indx]
bbox <- get_bbox (c (-0.13, 51.50, -0.11, 51.52))
dat_RFH <- extract_osm_objects (key = 'building', bbox = bbox,
                                extra_pairs = c ('name',
                                                 'Royal.Festival.Hall'))
extra_pairs <- list (c ('addr:street', 'Stamford.St'),
                     c ('addr:housenumber', '150'))
dat_ST <- extract_osm_objects (key = 'building', extra_pairs = extra_pairs,
                            bbox = bbox)
```
```{r}
london <- list (dat_H = dat_H, dat_HP = dat_HP, dat_BNR = dat_BNR,
                dat_BR = dat_BR, dat_BC = dat_BC, dat_A = dat_A, dat_P = dat_P,
                dat_T = dat_T, dat_RFH = dat_RFH, dat_ST = dat_ST)
devtools::use_data (london, overwrite = TRUE, compress = 'xz')
format (file.size ('./data/london.rda'), big.mark = ',') # 189,984
```

---
title: "Basic Maps"
author: "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        #number_sections: true
        theme: flatly
vignette: >
  %\VignetteIndexEntry{Basic Maps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


The **R** package `osmplotr` uses OpenStreetMap (OSM) data to produce highly
customisable maps. Data are downloaded via the 
[`osmdata` package](https://cran.r-project.org/package=osmdata), and different
aspects of map data - such as roads, buildings, parks, or water bodies - are
able to be visually customised.  This vignette demonstrates both data
downloading and the creation of simple maps. The subsequent vignette
(['data-maps'](https://cran.r-project.org/package=osmplotr))
demonstrates how `osmplotr` enables user-defined data to be visualised using OSM
data.  The maps in this vignette represent a small portion of central London,
U.K. 


# 1. Introduction

A map can be generated using the following simple steps:
```{r, echo = FALSE, message = FALSE}
library (osmplotr)
map_dpi <- 72 # dpi res for all maps
```

```{r, echo = FALSE, message = FALSE}
dat_B <- rbind (london$dat_BR, london$dat_BNR)
dat_H <- rbind (london$dat_H, london$dat_HP)
dat_T <- london$dat_T
```
1. Specify the bounding box for the desired region
```{r get-bbox}
bbox <- get_bbox (c (-0.13, 51.51, -0.11, 51.52))
```
2. Download the desired data---in this case, all building perimeters.
```{r extract buildings, eval = FALSE}
dat_B <- extract_osm_objects (key = "building", bbox = bbox)
```
3. Initiate an `osm_basemap` with desired background (`bg`) colour
```{r basemap1}
map <- osm_basemap (bbox = bbox, bg = "gray20")
```
4. Add desired plotting objects in the desired colour.
```{r add objects1}
map <- add_osm_objects (map, dat_B, col = "gray40")
```
5. Print the map
```{r, eval = FALSE}
print_osm_map (map)
```
```{r map1-print, echo = FALSE}
print_osm_map (map, filename = "map_a1.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a1.png)

The function `print_osm_map` creates a graphics device that is scaled to the
bounding box of the map.  Note also that `osmplotr` maps contain no margins and
fill the entire plot area, reflecting the general layout of most printed maps.
Additional capabilities of `osmplotr` are described in the following sections,
beginning with downloading and extraction of data.

# 2. Downloading Data

The package [`osmdata`](https://cran.r-project.org/package=osmdata) is used to
download data from 'OpenStreetMap' using the 'overpass' API [overpass
API](https://overpass-api.de). Data may be returned in either
['Simple Features' (`sf`)](https://cran.r-project.org/package=sf) or
['R Spatial' (`sp`)](https://cran.r-project.org/package=sp) form. `osmplotr` has
a convenience function, `extract_osm_objects`, to allow direct import, or the
functions of [`osmdata`](https://cran.r-project.org/package=osmdata) can also be
used directly. 

Data of a particular type can be extracted by specifying the appropriate OSM
`key`, as in the above example:
```{r buildings-highways, eval = FALSE}
bbox <- get_bbox (c(-0.13, 51.51, -0.11, 51.52))
dat_B <- extract_osm_objects (key = "building", bbox = bbox)
dat_H <- extract_osm_objects (key = "highway", bbox = bbox)
```
These objects are of appropriate `Spatial` classes:
```{r}
class (dat_B); class (dat_H)
class (dat_B$geometry); class (dat_H$geometry)
```
`Spatial` ([`sp`](https://cran.r-project.org/package=sf)) objects may be
returned with,
```{r, eval = FALSE}
dat_B <- extract_osm_objects (key = "building", bbox = bbox, sf = FALSE)
```
otherwise `sf` is used as the default format.  The Simple Features (`sf`)
objects with polygons of London buildings and linestrings of highways
respectively contain 
```{r}
nrow (dat_B); nrow (dat_H)
```
... 1,759 building polygons and 1,133 highway lines.  `extract_osm_objects` also
accepts `key-value` pairs which are passed to the
[overpass API](https://overpass-api.de) :
```{r trees, eval = FALSE}
dat_T <- extract_osm_objects (key = "natural", value = "tree", bbox = bbox)
```
Trees are located by single coordinates and are thus point objects:
```{r}
class (dat_T$geometry); nrow (dat_T)
```


## 2.1 `osmdata`

The [`osmdata`](https://cran.r-project.org/package=osmdata) package provides a
more powerful interface for downloading OSM data, and may be used directly with
`osmplotr`. The `osmplotr` function `extract_osm_objects` is effectively just a
convenience wrapper around `omsdata` functionality. The primary differences
between the two are:

1. `osmdata` returns *all* spatial data for a given query; that is, *all*
   points, lines, polygons, multilines, and multipolygons, while `osmplotr`
   returns a single specified geometric type.
2. `osmplotr` accepts multiple `key-value` pairs in a single call to
   `extract_osm_objects`, which the equivalent `osmdata` function,
   `add_feature`, accepts only a single `key-value` pair, with queries
   successively build through multiple calls to `add_feature`.

These differences are illustrated in the following code which generates
identical results in both cases (with namespaces explicitly given to aid
clarity),
```{r, eval = FALSE}
dat1 <- osmplotr::extract_osm_objects (key = "highway", value = "!primary",
                                       bbox = bbox)
dat2 <- osmdata::opq (bbox = bbox) %>%
    add_feature (key = "highway") %>%
    add_feature (key = "highway", value = "!primary") %>%
    osmdata_sf ()
dat2 <- dat2$osm_lines
```
The `osmdata` function `opq()` constructs an overpass query, with successive
calls to `add_feature` extending the query until it is finally submitted to
overpass by `osmdata_sf()` (or the `sp` version `osmdata_sp()`).

Note that `add_feature()` has to be called twice in this case, because a single
call to `add_feature (key = 'highway", value = "!primary")` would request *all*
features that are not primary highways. The initial query for `key = "highway"`
ensures that only npn-primary highways are returned.


## 2.2 Negation

As demonstrated above, negation can be specified by pre-pending `!` to the
`value` argument so that, for example, all `natural` objects that are **not**
trees can be extracted with
```{r not trees, eval = FALSE}
dat_NT <- extract_osm_objects (bbox = bbox, key = "natural", value = "!tree")
```
```{r, echo = FALSE}
message ("Cannot determine return type; maybe specify explicitly?")
```
The message is generated because of course a request for anything that is not a
tree could be for any kind of spatial object. `osmplotr` makes several educated
guesses in the absence of specified return types, but these can always be forced
with the `return_type` parameter:
```{r not tree points, eval = FALSE}
pts_NT <- extract_osm_objects (bbox = bbox, key = "natural", value = "!tree",
                               return_type = "points")
```

`london$dat_H` contains all non-primary highways, and was extracted with the
call demonstrated above, while `london$dat_HP` contains the corresponding set of
exclusively primary highways.  An `osmplotr` request for `key = "highway"`
automatically returns line objects (although, again, other kinds of objects may
be forced through specifying `return_type`).

## 2.3 Additional `key-value` pairs

Any number of `key-value` pairs may be passed to `extract_osm_objects`. For
example, a named building can be extracted with
```{r royal festival hall, eval = FALSE}
bbox <- get_bbox (c(-0.13, 51.50, -0.11, 51.52))
extra_pairs <- c ("name", "Royal.Festival.Hall")
dat <- extract_osm_objects (key = "building", extra_pairs = extra_pairs,
                                       bbox = bbox)
```
These data are stored in `london$dat_RFH`. Note that periods or dots are used for
white space, and in fact symbolise (in `grep` terms) any character whatsoever.
The polygon of a building at a particular street address can be extracted with
```{r stamford st 150, eval = FALSE}
extra_pairs <- list (c ("addr:street", "Stamford.St"),
                     c ("addr:housenumber", "150"))
dat <- extract_osm_objects (key = "building", extra_pairs = extra_pairs,
                                      bbox = bbox)
```
These data are stored as `london$dat_ST`.  Note that addresses generally require
combining both `addr:street` with `addr:housenumber`.


## 2.4 Downloading with `osm_structures` and `make_osm_map`

The functions `osm_structures` and `make_osm_map` aid both downloading multiple
OSM data types and plotting (with the latter described below).  `osm_structures`
returns a `data.frame` of OSM structure types, associated `key-value` pairs,
unique suffices which may be appended to data structures for storage purposes,
and suggested colours. Passing this list to `make_osm_map` will return a list of
the requested OSM data items, named through combining the `dat_prefix` specified
in `make_osm_map` and the suffices specified in `osm_structures`.

```{r}
osm_structures ()
```
Many structures are identified by keys only, in which cases the values are empty
strings.
```{r}
osm_structures()$value [1:4]
```
The last row of `osm_structures` exists only to define the background colour of
the map, as explained below 
([4.3 Automating map production](#4.3 Automating map production)).

The suffices include as many letters as are necessary to represent all unique
structure names.  `make_osm_map` returns a list of two components:

1. `osm_data` containing the data objects passed in the `osm_structures`
   argument. Any existing `osm_data` may also be submitted to `make_osm_map`, in
   which case any objects not present in the submitted data will be
   appended to the returned version. If `osm_data` is not submitted, all objects in
   `osm_structures` will be downloaded and returned.
2. `map` containing the `ggplot2` map objects with layers overlaid according to
   the sequence and colour schemes specified in `osm_structures`

The data specified in `osm_structures` can then be downloaded simply by calling:
```{r, eval = FALSE}
dat <- make_osm_map (structures = osm_structures (), bbox = bbox)
```
```{r, echo = FALSE}
dat1 <- list (dat_BU = NULL, dat_A = NULL, dat_W = NULL, dat_G = NULL,
              dat_N = NULL, dat_P = NULL, dat_H = NULL, dat_BO = NULL,
              dat_T = NULL)
dat <- list (osm_data = dat1, map = ggplot2::ggplot ())
```
```{r}
names (dat); sapply (dat, class); names (dat$osm_data)
```
The requested data are contained in `dat$osm_data`.  A list of desired structures
can also be passed to this function, for example,
```{r}
osm_structures (structures = c("building", "highway"))
```
Passing this to `make_osm_map` will download only these two structures.
Finally, note that the example of,
```{r}
osm_structures (structures = "grass")
```
demonstrates that `osm_structures` converts a number of common `keys` to
OSM-appropriate `key-value` pairs.

### 2.4.1 The `london` data of `osmplotr`

To illustrate the use of `osm_structures` to download data, this section
reproduces the code that was used to generate the `london` data object which
forms part of the `osmplotr` package.
```{r}
structures <- c ("highway", "highway", "building", "building", "building",
                 "amenity", "park", "natural", "tree")
structs <- osm_structures (structures = structures, col_scheme = "dark")
structs$value [1] <- "!primary"
structs$value [2] <- "primary"
structs$suffix [2] <- "HP"
structs$value [3] <- "!residential"
structs$value [4] <- "residential"
structs$value [5] <- "commercial"
structs$suffix [3] <- "BNR"
structs$suffix [4] <- "BR"
structs$suffix [5] <- "BC"
```
Suffices are generated automatically from structure names only, not values, and
the suffices for negated forms must therefore be specified manually.  The
`london` data can then be downloaded by simply calling `make_osm_map`:
```{r, eval = FALSE}
london <- make_osm_map (structures = structs, bbox = bbox)$osm_data
```
The requested data are contained in the `$osm_data` list item. `make_osm_map`
also returns a `$map` item which is described below 
(see [4.3 Automating map production](#4.3 make-osm-map)).

# 3. Downloading connected highways

The visualisation functions described in the second `osmplotr` vignette
([Data maps](https://cran.r-project.org/package=osmplotr))
enable particular regions of maps
to be highlighted. While it may often be desirable to highlight regions
according to a user's own data, `osmplotr` also enables regions to be defined by
providing a list of the names of encircling highways. The function which
achieves this is `connect_highways`, which returns a sequential matrix of
coordinates from those segments of the named highways which connected
continuously and sequentially to form a single enclosed space. An example is,
```{r, echo = FALSE}
load (system.file ("extdata", "hwys.rda", package = "osmplotr"))
highways1 <- hwys [[1]]
highways2 <- hwys [[2]]
highways3 <- hwys [[3]]
```
```{r, eval = FALSE}
highways <- c ("Monmouth.St", "Short.?s.Gardens", "Endell.St", "Long.Acre",
               "Upper.Saint.Martin")
highways1 <- connect_highways (highways = highways, bbox = bbox)
```
Note the use of the [regex](https://en.wikipedia.org/wiki/Regular_expression)
character `?` which declares that the previous character is optional. This
matches both "Shorts Gardens" and "Short's Gardens", both of which appear in OSM
data.
```{r}
class (highways1); length (highways1); highways1 [[1]] [[1]]
```

The extraction of bounding polygons from named highways is not fail-safe, and may
generate various warning messages.  To understand the kinds of conditions under
which it may not work, it is useful to examine `connect_highways` in more
detail.

## 3.1 `connect_highways` in detail

`connect_highways` takes a list of OpenStreetMap highways and sequentially
connects closest nodes of adjacent highways until the set of named highways
connects to form a cycle.  Cases where no circular connection is possible
generate an error message. The routine proceeds through the three stages of,

1. Adding intersection nodes to junctions of ways where these don't already
   exist

2. Filling a connectivity matrix between the listed highways and extracting the
   **longest** cycle connecting all of them 

3. Inserting extra connections between highways until the length of the longest
   cycle is equal to `length (highways)`.

This procedure can not be guaranteed fail-safe owing both to the inherently
unpredictable nature of OpenStreetMap, as well as to the unknown relationships
between named highways. To enable problematic cases to be examined and hopefully
resolved, `connect_highways` has a `plot` option:
```{r connect_highways, fig.width = 4, message = FALSE, eval = FALSE}
bbox_big <- get_bbox (c(-0.15, 51.5, -0.10, 51.52))
highways <- c ("Kingsway", "Holborn", "Farringdon.St", "Strand",
               "Fleet.St", "Aldwych")
highway_list <- connect_highways (highways = highways, bbox = bbox_big,
                                  plot = TRUE)
```
```{r connect-highways-manual-plot, fig.width = 4, echo = FALSE}
load (system.file ("extdata", "hwys.rda", package = "osmplotr"))
ways <- hwys$highways4
osmplotr:::plot_highways (ways)
ways <- osmplotr:::connect_single_ways (ways)
ways <- osmplotr:::get_highway_cycle (ways)

conmat <- osmplotr:::get_conmat (ways)
cycles <- try (ggm::fundCycles (conmat), TRUE)
cyc <- cycles [[which.max (sapply (cycles, nrow))]]

path <- osmplotr:::sps_through_cycle (ways, cyc)
lines (path [, 1], path [, 2], lwd = 2, lty = 2)
```

The plot depicts each highway in a different colour, along with numbers at start
and end points of each segment. This plot reveals in this case that highway#6
("Aldwych") is actually nested within two components of highway#4 ("Strand").
`connect_highways` searches for the shortest path connecting all named highways,
and since "Strand" connects to both highways#1 and #5, the shortest path
excludes #6. This exclusion of one of the named components generates the
warning message. 

These connected polygons returned from `connect_highways` can then be used to
highlight the enclosed regions within maps, as demonstrated in the second
vignette,
['Data Maps'](https://cran.r-project.org/package=osmplotr).


# 4. Producing maps

Maps will generally contain multiple kinds of OSM data, for example,
```{r, eval = FALSE}
dat_B <- extract_osm_objects (key = "building", bbox = bbox)
dat_H <- extract_osm_objects (key = "highway", bbox = bbox)
dat_T <- extract_osm_objects (key = "natural", value = "tree", bbox = bbox)
```

As illustrated above, plotting maps requires first making a basemap with a
specified background colour. Portions of maps can also be plotted by creating a
`basemap` with a smaller bounding box.
```{r map2}
bbox_small <- get_bbox (c(-0.13, 51.51, -0.11, 51.52))
map <- osm_basemap (bbox = bbox_small, bg = "gray20")
map <- add_osm_objects (map, dat_H, col = "gray70")
map <- add_osm_objects (map, dat_B, col = "gray40")
```
`map` is then a `ggplot2` which may be viewed simply by passing it to
`print_osm_map`:
```{r, eval = FALSE}
print_osm_map (map)
```
```{r map2-print, echo = FALSE}
print_osm_map (map, filename = "map_a2.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a2.png)


Other graphical parameters can also be passed to `add_osm_objects`, such as
border colours or line widths and types. For example,
```{r map3, eval = FALSE}
map <- osm_basemap (bbox = bbox_small, bg = "gray20")
map <- add_osm_objects (map, dat_B, col = "gray40", border = "orange",
                        size = 0.2)
print_osm_map (map)
```
```{r map3-print, echo = FALSE}
map <- osm_basemap (bbox = bbox_small, bg = "gray20")
map <- add_osm_objects (map, dat_B, col = "gray40", border = "orange",
                        size = 0.2)
print_osm_map (map, filename = "map_a3.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a3.png)

The `size` argument is passed to the corresponding `ggplot2` routine for
plotting polygons, lines, or points, and respectively determines widths of lines
(for polygon outlines and for lines), and sizes of points.  The `col` argument
determines the fill colour of polygons, or the colour of lines or points.

```{r map4, eval = FALSE}
map <- add_osm_objects (map, dat_H, col = "gray70", size = 0.7)
map <- add_osm_objects (map, dat_T, col = "green", size = 2, shape = 1)
print_osm_map (map)
```
```{r map4-print, echo = FALSE}
map <- add_osm_objects (map, dat_H, col = "gray70", size = 0.7)
map <- add_osm_objects (map, dat_T, col = "green", size = 2, shape = 1)
print_osm_map (map, filename = "map_a4.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a4.png)

Note also that the `shape` parameter determines the point shape, for details of
which see `?ggplot2::shape`. Also note that plot order affects the final
outcome, because components are sequentially overlaid and thus the same map
components plotted in a different order will generally produce a different
result.

## 4.1 Saving Maps

The function `print_osm_map()` can be used to print either to on-screen
graphical devices or to graphics files (see, for example, `?png` for a list of
possible graphics devices). Sizes and resolutions of devices may be
specified with the appropriate parameters. Device dimensions are scaled by
default to the proportions of the bounding box (although this can be
over-ridden).

A screen-based device simply requires
```{r, eval = FALSE}
print_osm_map (map)
```
while examples of writing higher resolution versions to files include:
```{r, eval = FALSE}
print_osm_map (map, filename = "map.png", width = 10,
               units = "in", dpi = map_dpi)
print_osm_map (map, filename = "map.eps", width = 1000,
               units = "px", dpi = map_dpi)
print_osm_map (map, filename = "map", device = "jpeg", width = 10, units = "cm")
```

## 4.2 Plotting different OSM Structures

The ability demonstrated above to use negation in `extract-osm-objects` allows
different kinds of the same object to be visually contrasted, for example
primary and non-primary highways:

```{r, eval = FALSE}
dat_HP <- extract_osm_objects (key = "highway", value = "primary", bbox = bbox)
dat_H <- extract_osm_objects (key = "highway", value = "!primary", bbox = bbox)
```
```{r, echo = FALSE}
dat_HP <- london$dat_HP
dat_H <- london$dat_H
```
```{r map5, eval = FALSE}
map <- osm_basemap (bbox = bbox_small, bg = "gray20")
map <- add_osm_objects (map, dat_H, col = "gray50")
map <- add_osm_objects (map, dat_HP, col = "gray80", size = 2)
print_osm_map (map)
```
```{r map5-print, echo = FALSE}
map <- osm_basemap (bbox = bbox_small, bg = "gray20")
map <- add_osm_objects (map, dat_H, col = "gray50")
map <- add_osm_objects (map, dat_HP, col = "gray80", size = 2)
print_osm_map (map, filename = "map_a5.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a5.png)

The additional `key-value` pairs demonstrated above (for Royal Festival Hall,
`dat_RFH` and 150 Stamford Street, `dat_ST`) also demonstrated above allow for
highly customised maps in which distinct objects are plotting with different
colour schemes.
```{r, echo = FALSE}
dat_RFH <- london$dat_RFH
dat_ST <- london$dat_ST
```
```{r map7, eval = FALSE}
bbox_small2 <- get_bbox (c (-0.118, 51.504, -0.110, 51.507))
map <- osm_basemap (bbox = bbox_small2, bg = "gray95")
map <- add_osm_objects (map, dat_H, col = "gray80")
map <- add_osm_objects (map, dat_HP, col = "gray20", size = 2)
map <- add_osm_objects (map, dat_RFH, col = "orange", border = "red", size = 2)
map <- add_osm_objects (map, dat_ST, col = "skyblue", border = "blue", size = 2)
print_osm_map (map)
```
```{r map7-print, echo = FALSE}
bbox_small2 <- get_bbox (c (-0.118, 51.504, -0.110, 51.507))
map <- osm_basemap (bbox = bbox_small2, bg = "gray95")
map <- add_osm_objects (map, dat_H, col = "gray80")
map <- add_osm_objects (map, dat_HP, col = "gray60", size = 2)
map <- add_osm_objects (map, dat_RFH, col = "orange", border = "red", size = 2)
map <- add_osm_objects (map, dat_ST, col = "skyblue", border = "blue", size = 2)
print_osm_map (map, filename = "map_a7.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a7.png)

## 4.3 Filling within boundary lines

Different portions of a map may sometimes be delineated by lines, for example
with coastlines which are always represented in OpenStreetMap as lines. Plotting
the water or land either side of a coastline in a single block of colour
requires the regions to be polygons, not lines. `osmplotr` has a function
`osm_line2poly()` which converts boundary lines extending beyond a given
bounding box into polygons encircling the perimeter of the bounding box. An
example is given in `?osm_line2poly`, using both the `osmdata` package to obtain
the bounding box of a named region, and the `magrittr` pipe operator.
```{r, eval = FALSE}
library (osmdata)
bb <- osmdata::getbb ("melbourne, australia")
coast <- extract_osm_objects (bbox = bb, key = "natural", value = "coastline",
                              return_type = "line")
coast <- osm_line2poly (coast, bbox = bb)
map <- osm_basemap (bbox = bb) %>%
        add_osm_objects (coast [[1]], col = "lightsteelblue") %>%
        print_osm_map ()
```
The `osm_line2poly()` function returns a list of two `sf` polygons. For
coastline, one of these will correspond to water, one to land. In the preceding
example, the first polygon is the ocean, which is coloured in
`"lightsteelblue"`. Users must determine for themselves which polygon is to be
plotted in which colour. Note that `osm_line2poly()` only accepts `sf`-formatted
data, and not `sp`.


## 4.4 Automating map production

As indicated above 
([2.4 Downloading with `osm_structures` and `make_osm_map`](#2.4 downloading2)),
the production of maps overlaying various type of OSM objects is facilitated
with `make_osm_map`.  The structure of a map is defined by `osm_structures` as
described above.

Producing a map with customised data is as simple as,
```{r map8, eval = FALSE}
structs <- c ("highway", "building", "park", "tree")
structures <- osm_structures (structures = structs, col_scheme = "light")
dat <- make_osm_map (structures = structures, bbox = bbox_small)
print_osm_map (dat$map)
```
```{r, echo = FALSE}
structs <- c ("highway", "building", "park", "tree")
structures <- osm_structures (structures = structs, col_scheme = "light")
osm_dat <- list (dat_B = dat_B, dat_H = dat_H, dat_P = london$dat_P,
                 dat_A = london$dat_A, dat_P = london$dat_P,
                 dat_T = london$dat_T)
dat <- make_osm_map (structures = structures, osm_data = osm_dat,
                     bbox = bbox)
print_osm_map (dat$map, filename = "map_a8.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a8.png)

Calling `make_osm_map()` downloads the requested structures within the given
`bbox` and returns a list of two components, the first of which contains the
downloaded data:
```{r}
names (dat); names (dat$osm_data)
```

Pre-downloaded data may also be passed to `make_osm_map()`
```{r map9, eval = FALSE}
dat <- make_osm_map (structures = structures, osm_data = dat$osm_data,
                     bbox = bbox)
print_osm_map (dat$map)
```
```{r map9-print, echo = FALSE}
dat <- make_osm_map (structures = structures, osm_data = osm_dat, bbox = bbox)
print_osm_map (dat$map, filename = "map_a9.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a9.png)

Note that omitting the bounding box argument (`bbox`) produces a map with a bounding
box is extracted as the **largest** box spanning all objects in `osm_data`. This
may be considerably larger than the desired boundaries, particularly because
highways are returned by `overpass` in their entirety, and will generally extend
well beyond the specified bounding box.

Finally, objects in maps are overlaid on the plot according to the order of rows
in `osm_structures`, with the single exception that `background` is plotted
first.  This order can be readily changed or restricted simply by submitting
structures in a desired order.

```{r}
structs <- c ("amenity", "building", "highway", "park")
osm_structures (structs, col_scheme = "light")
```


## 4.5 Axes

Axes may be added to maps using the `add_axes` function. In contrast to many `R`
packages for producing maps, maps in `osmplotr` fill the entire plotting space,
and axes are added *internal* to this space. The separate function for adding
axes allows them to be overlaid on top of all previous layers.

Axes added to a dark version of the previous map look like this:
```{r map10}
structures <- osm_structures (structures = structs, col_scheme = "dark")
dat <- make_osm_map (structures = structures, osm_data = dat$osm_dat,
                     bbox = bbox_small)
map <- add_axes (dat$map, colour = "black")
```
Note that, as described above, `make_osm_map` returns a list of two items: (i)
potentially modified data (in `$osm_data`) and (ii) the map object (in `$map`).
All other `add_` functions take a map object as one argument and return the
single value of the modified map object.
```{r, eval = FALSE}
print_osm_map (map)
```
```{r map10-print, echo = FALSE}
print_osm_map (map, filename = "map_a10.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a10.png)

This map reveals that the axes and labels are printed above semi-transparent
background rectangles, with transparency controlled by the `alpha` parameter.
Axes are always plotted on the left and lower side, but positions can be
adjusted with the `pos` parameter which specifies the positions of axes and
labels relative to entire plot device

```{r map11, eval = FALSE}
map <- add_axes (map, colour = "blue", pos = c(0.1, 0.2),
                      fontsize = 5, fontface = 3, fontfamily = "Times")
print_osm_map (map)
```
```{r map11-print, echo = FALSE}
map <- add_axes (map, colour = "blue", pos = c(0.1, 0.2),
                      fontsize = 5, fontface = 3, fontfamily = "Times")
print_osm_map (map, filename = "map_a11.png", width = 600,
               units = "px", dpi = map_dpi)
```
![](map_a11.png)

The second call to `add_axes` overlaid additional axes on a map that already had
axes from the previous call.  This call also demonstrates how sizes and other
font characteristics of text labels can be specified.  

Finally, the current version of `osmplotr` does not allow text labels of axes to
be rotated. (This is because the semi-transparent underlays are generated with
`ggplot2::geom_label` which currently prevents rotation.)

Click on the following link to proceed to the second `osmplotr` vignette:
[Data maps](https://cran.r-project.org/package=osmplotr)
---
title: "Data Maps"
author: "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        theme: flatly
vignette: >
  %\VignetteIndexEntry{Data Maps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This vignette extends from the vignette 
([Basic maps](https://cran.r-project.org/package=osmplotr))
to demonstrate how `osmplotr` enables the graphical properties of OpenStreetMap
objects to be modified according to user-provided data.  Categorical data can be
plotted by highlighting defined regions with different colours using
`add_osm_groups`, while continuous data can be plotted with `add_osm_surface`.


```{r load, message = FALSE, eval = TRUE}
library (osmplotr)
```
```{r, echo = FALSE, message = FALSE}
map_dpi <- 72 # dpi res for all maps
```
As in the first vignette, maps produced in this vignette contain data for a
small portion of central London, U.K.
```{r}
bbox <- get_bbox (c(-0.13, 51.51, -0.11, 51.52))
```

```{r, echo = FALSE, message= FALSE}
dat_B <- rbind (london$dat_BR, london$dat_BNR)
dat_H <- rbind (london$dat_H, london$dat_HP)
dat_HP <- london$dat_HP
```

# 1. Categorical data: `add_osm_groups`

The function `add_osm_groups` enables spatially-defined groups to be plotted in
different colours.  The two primary arguments are `obj`, which defines the OSM
structure to be used for plotting the regions, and `groups` which is a list of
geometric coordinates defining the desired regions. An example of an `obj` is
the Simple Features ([`sf`](https://cran.r-project.org/package=sf)) `data.frame`
of building polygons downloaded in the first vignette with the following line
```{r, eval = FALSE}
dat_B <- extract_osm_objects (key = "building", bbox = bbox)
```
These data may be obtained by simply combining the data provided with the
package of residential and non-residential buildings to give all buildings as
```{r}
dat_B <- rbind (london$dat_BNR, london$dat_BR)
```

The most direct way to define `groups` is through specifying coordinates of
boundary points:
```{r map1, eval = FALSE}
pts <- cbind (c (-0.115, -0.125, -0.125, -0.115),
              c (51.513, 51.513, 51.517, 51.517))

map <- osm_basemap (bbox = bbox,
                    bg = "gray20")

map <- add_osm_groups (map,
                       dat_B,
                       groups = pts,
                       cols = "orange",
                       bg = "gray40")

print_osm_map (map)
```
```{r map1-print, echo = FALSE}

pts <- cbind (c (-0.115, -0.125, -0.125, -0.115),
              c (51.513, 51.513, 51.517, 51.517))

map <- osm_basemap (bbox = bbox,
                    bg = "gray20")

map <- add_osm_groups (map,
                       dat_B,
                       groups = pts,
                       cols = "orange",
                       bg = "gray40")

print_osm_map (map,
               filename = "map_b1.png",
               width = 600,
               units = "px", dpi = map_dpi)
```
![](map_b1.png)

Multiple groups can be defined by passing a list of multiple sets of point
coordinates to the `groups` argument of `add_osm_groups`, and specifying
corresponding colours.

```{r map2, eval = FALSE}
pts2 <- cbind (c (-0.111, -0.1145, -0.1145, -0.111),
               c (51.517, 51.517, 51.519, 51.519))

map <- osm_basemap (bbox = bbox,
                    bg = "gray20")

map <- add_osm_groups (map,
                       dat_B,
                       groups = list (pts, pts2),
                       cols = c ("orange", "tomato"),
                       bg = "gray40")

print_osm_map (map)
```
```{r map2-print, echo = FALSE}

pts2 <- cbind (c (-0.111, -0.1145, -0.1145, -0.111),
               c (51.517, 51.517, 51.519, 51.519))

map <- osm_basemap (bbox = bbox,
                    bg = "gray20")

map <- add_osm_groups (map,
                       dat_B,
                       groups = list (pts, pts2),
                       cols = c ("orange", "tomato"),
                       bg = "gray40")

print_osm_map (map,
               filename = "map_b2.png",
               width = 600,
               units = "px",
               dpi = map_dpi)
```
![](map_b2.png)

The `bg` argument specifies the colour of any objects lying outside the
boundaries of the specified groups. If this argument is not given, then all
objects are assigned to the nearest group, so that the groups fill the entire
map.

```{r map3, eval = FALSE}

map <- osm_basemap (bbox = bbox,
                    bg = "gray20")

map <- add_osm_groups (map,
                       dat_B,
                       groups = list (pts, pts2),
                       cols = c ("orange", "tomato"))

print_osm_map (map)
```
```{r map3-print, echo = FALSE}

map <- osm_basemap (bbox = bbox,
                    bg = "gray20")

map <- add_osm_groups (map,
                       dat_B,
                       groups = list (pts, pts2),
                       cols = c ("orange", "tomato"))

print_osm_map (map,
               filename = "map_b3.png",
               width = 600,
               units = "px",
               dpi = map_dpi)
```
![](map_b3.png)

Now that you've seen the general workflow of `osmplotr`, let's repeat the previous code, but streamline it with `magrittr`'s `%>%` function. This allows us to pipe the functions together instead of re-assigning the `map` variable. 

```{r map_b3_pipe, eval = FALSE}

library(magrittr)

osm_basemap(bbox = bbox,
            bg = "gray20") %>%
  add_osm_groups(dat_B,
                 groups = list(pts, pts2),
                 cols = c("orange", "tomato")) %>%
  print_osm_map()
```
```{r map_b3_pipe-print, echo = FALSE}

library(magrittr)

map <- osm_basemap (bbox = bbox,
                    bg = "gray20")

map <- add_osm_groups (map,
                       dat_B,
                       groups = list (pts, pts2),
                       cols = c ("orange", "tomato"))

print_osm_map (map,
               filename = "map_b3_pipe.png",
               width = 600,
               units = "px",
               dpi = map_dpi)
```

![](map_b3_pipe.png)


## 1.1 Hulls around groups

`add_osm_groups` includes the argument `make_hull` which specifies whether
convex hulls should be fitted around the points defining the provided `groups`,
or whether the `groups` already define their own boundaries (the default
behaviour). If a point is added internal to the four points defining the first
of the above groups, then the group boundary will connect to that point and
create a concave shape.

```{r map5, eval = FALSE}
pts <- rbind (pts, c (-0.12, 51.515))

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = pts,
                  cols = "orange",
                  bg = "gray40") %>%
  print_osm_map ()
```
```{r map5-print, echo = FALSE}

pts <- rbind (pts, c (-0.12, 51.515))

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = pts,
                  cols = "orange",
                  bg = "gray40") %>%
  print_osm_map (filename = "map_b5.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b5.png)

The previous points started in the south-east and ended in the north-east, and
thus the concave boundary extends in between the two easterly points. Setting
`make_hull = TRUE` defines groups by the convex hulls surrounding them, which in
this case would revert this map to the initial map with the group defined by a
regular, convex perimeter.

## 1.2 Inclusive, exclusive, and bisected polygons

The highlighted regions of the previous maps are irregular because the default
behaviour of `add_osm_groups` is to include within a group only those OSM
objects which lie entirely within a group boundary.  `add_osm_groups` has a
`boundary` argument which defines whether objects should be assigned to groups
inclusively (`boundary > 0`) or exclusively (`boundary < 0`), or whether they
should be precisely bisected by a group boundary (`boundary = 0`).  The previous
maps illustrate the default option (`boundary = -1`), while the two other
options produce the following maps. 

```{r map6, eval = FALSE}

osm_basemap (bbox = bbox, bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c("orange", "tomato"),
                  bg = "gray40",
                  boundary = 1) %>%
  print_osm_map ()
```
```{r map6-print, echo = FALSE}

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c("orange", "tomato"),
                  bg = "gray40",
                  boundary = 1) %>%
  print_osm_map (filename = "map_b6.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b6.png)

The inclusive option (`boundary>0`) includes all objects which have any points
lying within a boundary, meaning more objects are included resulting in larger
regions than the previous default exclusive option. Precisely
bisecting boundaries produces the following map.

```{r map7, eval = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c ("orange", "tomato"),
                  bg = "gray40",
                  boundary = 0) %>%
  print_osm_map ()
```
```{r map7-print, echo = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c ("orange", "tomato"),
                  bg = "gray40",
                  boundary = 0) %>%
  print_osm_map (filename = "map_b7.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b7.png)


The ability to combine different kinds of boundaries is particularly useful when
highlighting areas which partially contain large polygons such as parks. The
parks within the following maps were downloaded with
```{r, eval = FALSE}
dat_P <- extract_osm_objects (key = "park", bbox = bbox)
```
(Noting that, as described in the first vignette,
[Basic maps](https://cran.r-project.org/package=osmplotr),
both `extract_osm_objects` and `make_osm_map` convert several common keys to
appropriate `key-value` pairs, so
```{r}
osm_structures (structure = "park")
```
reveals that this `key` is actually converted to `key = "leisure"` and
`value = "park"`.) These data are also provided with the package as
`london$dat_P`.
```{r, echo = FALSE}
dat_P <- london$dat_P
```
Plotting buildings inclusively within each group and overlaying parks bisected
by the group boundaries produces the following map:
```{r map8, eval = FALSE}
col_park_in <- rgb (50, 255, 50, maxColorValue = 255)
col_park_out <- rgb (50, 155, 50, maxColorValue = 255)

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c("orange", "tomato"),
                  bg = "gray40",
                  boundary = 0) %>%
  add_osm_groups (dat_P,
                  groups = list (pts, pts2),
                  cols = rep (col_park_in, 2),
                  bg = col_park_out,
                  boundary = 0) %>%
  print_osm_map ()
```
```{r map8-print, echo = FALSE}
col_park_in <- rgb (50, 255, 50, maxColorValue = 255)
col_park_out <- rgb (50, 155, 50, maxColorValue = 255)

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c("orange", "tomato"),
                  bg = "gray40",
                  boundary = 0) %>%
  add_osm_groups (dat_P,
                  groups = list (pts, pts2),
                  cols = rep (col_park_in, 2),
                  bg = col_park_out,
                  boundary = 0) %>%
  print_osm_map (filename = "map_b8.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b8.png)

Bisection divides single polygons to form one polygon of points lying within a
given boundary and one polygon of points lying outside the boundary. The two
resultant polygons are often separated by visible gaps between locations at
which they are defined.  Because the layers of a plot are progressively
overlaid, such gaps can be avoided by initially plotting underlying layers using
`add_osm_objects` prior to grouping objects:

```{r map9, eval = FALSE}

map <- osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_objects (dat_P,
                   col = col_park_out) %>%
  add_osm_groups (dat_P,
                  groups = list (pts, pts2),
                  cols = rep (col_park_in, 2),
                  bg = col_park_out,
                  boundary = 0) %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c ("orange", "tomato"),
                  bg = "gray40",
                  boundary = 0)

map %>%
  print_osm_map ()
```
```{r map9-print, echo = FALSE}
map <- osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_objects (dat_P,
                   col = col_park_out) %>%
  add_osm_groups (dat_P,
                  groups = list (pts, pts2),
                  cols = rep (col_park_in, 2),
                  bg = col_park_out,
                  boundary = 0) %>%
  add_osm_groups (dat_B,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = c ("orange", "tomato"),
                  bg = "gray40",
                  boundary = 0)

map %>%
  print_osm_map (filename = "map_b9.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b9.png)

Bisections with `boundary = 0` will only be as accurate as the underlying OSM
data. This example was chosen to highlight that bisection may be inaccurate if
actual OSM points do not lie near to a desired bisection line. The larger a
map, the less visually evident are likely to be any such inaccuracies. Finally,
note that the plot order was changed to allow the building within the park to be
overlaid upon the grass surfaces. Plot order, whether controlled manually or
with `make_osm_map`, may often have to be tweaked to appropriately visualise all
objects.

The `boundary` argument has no effect if `bg` is not given, because in this case
all objects will be assigned to a group and there will be no boundaries between
groups and other, non-grouped objects.

## 1.3 Adjusting colours with `adjust_colours`

The `adjust_colours` function allows different groups to be highlighted with
slightly different colours for different kinds of OSM objects. For example, the
following code adds highways to the above map in slightly darkened versions of
the highlight colours (using `boundary = 1`, so any highways with any points lying
within the bounding box are included in the groups):
```{r map10, eval = FALSE}
#create separate data for all highways and primary highways
dat_H <- rbind (london$dat_H, london$dat_HP)
dat_HP <- london$dat_HP

# darken colours by aboud 20%
cols_adj <- adjust_colours (c ("orange", "tomato"),
                            adj = -0.2)

map %>%
  add_osm_groups (dat_HP,
                groups = list (pts, pts2),
                make_hull = TRUE,
                cols = cols_adj,
                bg = adjust_colours("gray40",
                                    adj = -0.4),
                boundary = 1, size = 2) %>%
  add_osm_groups (dat_H,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = cols_adj,
                  bg = adjust_colours ("gray40",
                                       adj = -0.2),
                  boundary = 1,
                  size = 1) %>%
  print_osm_map ()
```
```{r map10-print, echo = FALSE}
# darken colours by aboud 20%
cols_adj <- adjust_colours (c ("orange", "tomato"),
                            adj = -0.2)

map %>%
  add_osm_groups (dat_HP,
                groups = list (pts, pts2),
                make_hull = TRUE,
                cols = cols_adj,
                bg = adjust_colours("gray40",
                                    adj = -0.4),
                boundary = 1, size = 2) %>%
  add_osm_groups (dat_H,
                  groups = list (pts, pts2),
                  make_hull = TRUE,
                  cols = cols_adj,
                  bg = adjust_colours ("gray40",
                                       adj = -0.2),
                  boundary = 1,
                  size = 1) %>%
  print_osm_map (filename = "map_b10.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b10.png)

And of course `adjust_colours ("gray40", adj = -0.2)` is nothing other than
"gray32", and `adj = -0.4` gives "gray24".

## 1.4 Dark-on-Light Highlights 

A particularly effective way to highlight single regions within a map is through
using dark colours upon otherwise light coloured maps.

```{r map11, eval = FALSE}
osm_basemap (bbox = bbox, bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = pts,
                  cols = "gray40",
                  bg = "gray85",
                  boundary = 1) %>%
  add_osm_groups (dat_H,
                  groups = pts,
                  cols = "gray20",
                  bg = "gray70",
                  boundary = 0) %>%
  add_osm_groups (dat_HP,
                  groups = pts,
                  cols = "gray10",
                  bg = "white",
                  boundary = 0,
                  size = 1) %>%
  print_osm_map ()
```
```{r map11-print, echo = FALSE}
osm_basemap (bbox = bbox, bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = pts,
                  cols = "gray40",
                  bg = "gray85",
                  boundary = 1) %>%
  add_osm_groups (dat_H,
                  groups = pts,
                  cols = "gray20",
                  bg = "gray70",
                  boundary = 0) %>%
  add_osm_groups (dat_HP,
                  groups = pts,
                  cols = "gray10",
                  bg = "white",
                  boundary = 0,
                  size = 1) %>%
  print_osm_map (filename = "map_b11.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b11.png)


## 1.5 Visualising clustering data

One of the most likely uses of `add_osm_groups` is to visualise statistical
clusters. Clustering algorithms will generally produce membership lists which
may be mapped onto spatial locations. Each cluster can be defined as a matrix of
points in a single list of `groups`.  A general approach is illustrated here
with `groups` defined by single, randomly generated points.
```{r}
set.seed (2)
ngroups <- 12
x <- bbox [1, 1] + runif (ngroups) * diff (bbox [1, ])
y <- bbox [2, 1] + runif (ngroups) * diff (bbox [2, ])
groups <- as.list (data.frame (t (cbind (x, y))))
```
(The last line just transforms each row of the matrix into a list item.) Having
generated the points, a map of corresponding clusters can be generated by the
following simple code.

```{r map12, eval = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  cols = rainbow (length (groups))) %>%
  print_osm_map ()
```
```{r map12-print, echo = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  cols = rainbow (length (groups))) %>%
  print_osm_map (filename = "map_b12.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b12.png)

Although individual groups will generally be defined by collections of multiple
points, this example illustrates that they can also be defined by single points.
In such cases, the `bg` option should of course be absent, so that all remaining
points are allocated to the nearest groups.

This map also illustrates the kind of visual mess that may arise in attempts to
specify colours, particularly because the sequence of colours passed to
`add_osm_groups` will generally not map on to any particular spatial order, so
even if a pleasing colour scheme is submitted, the results may still be less
than desirable. Although it may be possible to devise pleasing schemes for small
numbers of groups, manually defined colour schemes are likely to become
impractical for larger numbers of groups.

```{r map13, eval = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  border_width = 2,
                  cols = heat.colors (length (groups))) %>%
  print_osm_map ()
```
```{r map13-print, echo = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  border_width = 2,
                  cols = heat.colors (length (groups))) %>%
  print_osm_map (filename = "map_b13.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b13.png)

Note the submitting any positive values to the additional `border_width` argument
causes `add_osm_groups` to drawn convex hull borders around the different
groups. Even this is not sufficient, however, to render the result particularly
visually pleasing or intelligible. To overcome this, `add_osm_groups` includes
an option described in the following section to generate spatially sensible
colour schemes for colouring distinct groups.

## 1.6 The Colour Matrix: Colouring Several Regions

An additional argument which may be passed to `add_osm_groups` is `colmat`,
an abbreviation of 'colour matrix'. If set to true (the default is `FALSE`),
group colours are specified by the function `colour_mat`. This function takes
a vector of four or more colours as input, wraps them around the four corners of a
rectangular grid, and spatially interpolates a chromatically regular grid between
these corners. To visual different schemes, it has a `plot` argument:

```{r, eval = FALSE}
cmat <- colour_mat (rainbow (4), plot = TRUE)
```
```{r, fig.width = 4, echo = FALSE}
plot.new ()
cmat <- colour_mat (rainbow (4), plot = TRUE)
```

This grid illustrates the default colours, `rainbow (4)`. The two-dimensional
colour field produced by `colour_mat` may also be rotated by a specified number
of degrees using the `rotate` argument.

```{r, eval = FALSE}
cmat <- colour_mat (rainbow (4), n = c(4, 8), rotate = 90, plot = TRUE)
```
```{r, fig.width = 4, echo = FALSE}
plot.new ()
cmat <- colour_mat (rainbow (4), n = c(4, 8), rotate = 90, plot = TRUE)
```

This example also illustrates that the size of colour matrices may also be
arbitrarily specified. Using the `colmat` option in `add_osm_groups` enables the
previous maps to be redrawn like this:
```{r map14, eval = FALSE}

osm_basemap (bbox = bbox,
             bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  border_width = 2,
                  colmat = TRUE,
                  cols = c("red", "green", "yellow", "blue"),
                  rotate = 180) %>%
  print_osm_map ()
```
```{r map14-print, echo = FALSE}
osm_basemap (bbox = bbox, bg = "gray95") %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  border_width = 2,
                  colmat = TRUE,
                  cols = c("red", "green", "yellow", "blue"),
                  rotate = 180) %>%
  print_osm_map (filename = "map_b14.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b14.png)

Note both that when `add_osm_groups` is called with `colmat = TRUE`, then `cols`
need only be of length 4, to specify the four corners of the colour matrix, and
also that the `rotate` argument can be submitted to `add_osm_groups` and passed
on to `colour_mat`.


## 1.7 Bounding areas within named highways

As explained in the first vignette,
[Basic maps](https://cran.r-project.org/package=osmplotr),
the function `connect_highways` takes a list of OSM highway names and a bounding
box, and returns the boundary of a polygon encircling the named highways. This
can be used to highlight selected regions simply by naming the highways which
encircle them, producing maps which look like this:
```{r connect-highways, eval = FALSE}
highways <- c ("Monmouth.St", "Short.?s.Gardens", "Endell.St", "Long.Acre",
               "Upper.Saint.Martin")
highways1 <- connect_highways (highways = highways, bbox = bbox)
highways <- c ("Endell.St", "High.Holborn", "Drury.Lane", "Long.Acre")
highways2 <- connect_highways (highways = highways, bbox = bbox)
highways <- c ("Drury.Lane", "High.Holborn", "Kingsway", "Great.Queen.St")
highways3 <- connect_highways (highways = highways, bbox = bbox)
```
Note the use of the [regex](https://en.wikipedia.org/wiki/Regular_expression)
character `?` in the first list of highway names, denoting the previous
character as optional.  This is necessary here because there are OSM sections
named both "Shorts Gardens" and "Short's Gardens".
```{r, echo = FALSE}
load (system.file ("extdata", "hwys.rda", package = "osmplotr"))
highways1 <- hwys [[1]]
highways2 <- hwys [[2]]
highways3 <- hwys [[3]]
con_h <- function (hwy) {

    hwy <- osmplotr:::connect_single_ways (hwy)
    hwy <- osmplotr:::get_highway_cycle (hwy)
    conmat <- osmplotr:::get_conmat (hwy)
    cycles <- try (ggm::fundCycles (conmat), TRUE)
    cyc <- cycles [[which.max (sapply (cycles, nrow))]]
    osmplotr:::sps_through_cycle (hwy, cyc)
}
highways1 <- con_h (highways1)
highways2 <- con_h (highways2)
highways3 <- con_h (highways3)
```
```{r}
class (highways1); nrow (highways1); nrow (highways2); nrow (highways3)
```
`connect_highways` returns a list of `SpatialPoints` representing the shortest
path that sequentially connects all of the listed highways. (Connecting all
listed highways may not necessarily be possible, in which case warnings will be
issued. As described in the first vignette, 
[Basic maps](https://cran.r-project.org/package=osmplotr),
`connect_highways` also has a `plot` option allowing problematic cases to be
visually inspected and hopefully corrected.)

These lists of highway coordinates can then be used to highlight the areas they
encircle. First group the highways and establish a colour scheme for the map:
```{r, echo = TRUE}
groups <- list (highways1, highways2, highways3)
cols_B <- c ("red", "orange", "tomato") # for the 3 groups
cols_H <- adjust_colours (cols_B, -0.2)
bg_B <- "gray40"
bg_H <- "gray60"
```
And then plot the map.
```{r map15, eval = FALSE}
osm_basemap (bbox = bbox, bg = "gray20") %>%
  add_osm_objects (dat_P,
                   col = col_park_out) %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  boundary = 1,
                  bg = bg_B,
                  cols = cols_B) %>%
  add_osm_groups (dat_H,
                  groups = groups,
                  boundary = 1,
                  bg = bg_H,
                  cols = cols_H) %>%
  add_osm_groups (dat_HP,
                  groups = groups,
                  boundary = 0,
                  cols = cols_H,
                  bg = bg_H,
                  size = 1) %>%
  print_osm_map ()

```
```{r map15-print, echo = FALSE}
osm_basemap (bbox = bbox, bg = "gray20") %>%
  add_osm_objects (dat_P,
                   col = col_park_out) %>%
  add_osm_groups (dat_B,
                  groups = groups,
                  boundary = 1,
                  bg = bg_B,
                  cols = cols_B) %>%
  add_osm_groups (dat_H,
                  groups = groups,
                  boundary = 1,
                  bg = bg_H,
                  cols = cols_H) %>%
  add_osm_groups (dat_HP,
                  groups = groups,
                  boundary = 0,
                  cols = cols_H,
                  bg = bg_H,
                  size = 1) %>%
  print_osm_map (filename = "map_b15.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b15.png)

These encircling highways are included in the `london` data provided with
`osmplotr`. 

# 2. Continuous data: `add_osm_surface`

The `add_osm_surface` function enables a continuous data surface to be overlaid
on a map. User-provided data is spatially interpolated across a map region and
OSM items coloured according to a specified continuous colour gradient. The data
must be provided as a data frame with three columns, '(x,y,z)', where '(x,y)'
are the coordinates of points at which data are given, and 'z' are the values to
be spatially interpolated across the map.

A simple data frame can be constructed as
```{r}
n <- 5
x <- seq (bbox [1, 1], bbox [1, 2], length.out = n)
y <- seq (bbox [2, 1], bbox [2, 2], length.out = n)
dat <- data.frame (
    x = as.vector (array (x, dim = c(n, n))),
    y = as.vector (t (array (y, dim = c(n, n)))),
    z = x * y
    )
head (dat)
```
And then passed to `add_osm_surface`
```{r map16, eval = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30)) %>%
  print_osm_map ()
```
```{r map16-print, echo = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30)) %>%
  print_osm_map (filename = "map_b16.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b16.png)

At present, `add_osm_surface` generates an warning if it is applied more than
once to any one kind of `Spatial` object (polygons or lines), as illustrated in
the following code (in which both `dat_H` and `dat_HP` are of class
`SpatialLinesDataFrame`:
```{r, eval = TRUE}
osm_basemap (bbox = bbox, bg = "gray20") %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = heat.colors (30)) %>%
  add_osm_surface (dat_H,
                   dat = dat,
                   cols = heat.colors (30))
```
This is because `add_osm_surface` creates new `ggplot2` aesthetic schemes for
each kind of object, and these schemes are not intended to be modified or
replaced within a single plot. The above map may still be printed, but the
warning means that the last provided colour scheme will be applied to all
objects of that class. This means that `osmplotr` can only overlay two distinct
colour schemes: one for all objects of class `SpatialLines`, and a potentially
different one for all objects of class `SpatialPolygons`.

Of course, any number of additional objects may be overlaid with
`add_osm_objects`, for example,
```{r map17, eval = FALSE}
cols_adj <- adjust_colours (heat.colors (30), -0.2)

map <- osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30)) %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = cols_adj,
                   size = 1.5) %>%
  add_osm_objects (dat_P,
                   col = rgb (0.1, 0.3, 0.1)) %>%
  add_osm_objects (dat_H,
                   col = "gray60")

map %>%
  print_osm_map ()
```
```{r map17-print, echo = FALSE}
cols_adj <- adjust_colours (heat.colors (30), -0.2)

map <- osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30)) %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = cols_adj,
                   size = 1.5) %>%
  add_osm_objects (dat_P,
                   col = rgb (0.1, 0.3, 0.1)) %>%
  add_osm_objects (dat_H,
                   col = "gray60")

map %>%
  print_osm_map (filename = "map_b17.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b17.png)


## 2.1 Colourbar legends for data surfaces

A colourbar legend for the surface may be added with `add_colourbar`. As with
`add_axes`, this function is provided separately to allow colourbars to be
overlaid only after all desired map items have been added. The only parameters
required for `add_colourbar` are the limits of the data (`zlims`) and the
colours (along with the `map`, a modified version of which is returned).
```{r map18, eval = FALSE}
map %>%
  add_colourbar (cols = terrain.colors (100),
                 zlims = range (dat$z)) %>%
  print_osm_map ()
```
```{r map18-print, echo = FALSE}
map %>%
  add_colourbar (cols = terrain.colors (100),
                 zlims = range (dat$z)) %>%
  print_osm_map (filename = "map_b18.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b18.png)

Note that the colours submitted to `add_colourbar` need not be the same as those
used to plot the surface. (Although using different colours is rarely likely to
be useful.) As for `add_axes`, and explained in the first vignette,
[Basic maps](https://cran.r-project.org/package=osmplotr),
the transparency of the boxes surrounding the elements of the colourbar may be
controlled by specifying the value of `alpha`. Both alignment and position may
also be adjusted, as illustrated in this example.
```{r map19, eval = FALSE}
cols_adj <- adjust_colours (heat.colors (30), -0.2)

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30)) %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = cols_adj,
                   size = 1.5) %>%
add_colourbar (cols = heat.colors (100),
               zlims = range (dat$z),
               alpha = 0.9,
               vertical = FALSE,
               barwidth = c(0.1, 0.12),
               barlength = c(0.5, 0.9),
               text_col = "blue",
               fontsize = 5,
               fontface = 3,
               fontfamily = "Times") %>%
  add_axes (colour = "blue",
            fontsize = 5,
            fontface = 3,
            fontfamily = "Times") %>%
  print_osm_map ()
```
```{r map19-print, echo = FALSE}
cols_adj <- adjust_colours (heat.colors (30), -0.2)

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30)) %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = cols_adj,
                   size = 1.5) %>%
add_colourbar (cols = heat.colors (100),
               zlims = range (dat$z),
               alpha = 0.9,
               vertical = FALSE,
               barwidth = c(0.1, 0.12),
               barlength = c(0.5, 0.9),
               text_col = "blue",
               fontsize = 5,
               fontface = 3,
               fontfamily = "Times") %>%
  add_axes (colour = "blue",
            fontsize = 5,
            fontface = 3,
            fontfamily = "Times") %>%
  print_osm_map (filename = "map_b19.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b19.png)

Both `barwidth` and `barlength` can be specified in terms of one or two numbers.
A single value for `barwidth` determines its relative width (0-1) from the
border (the right side if `vertical = TRUE` or the top if `vertical = FALSE`), while
two values determine the relative start and end positions of the sides of the
bar. A single value for `barlength` produces a bar of the given length centred
in the middle of the map, while two values determine its respective upper and
lower points (for `vertical = TRUE`) or left and right points (for
`vertical = FALSE`).

This example also demonstrates how colours, sizes, and other font
characteristics of text labels can be specified (with `text_col` determining the
colour of all elements of the colourbar other than the gradient itself).
Finally, as for `add_axes`, the text labels of colourbars are not currently able
to be rotated because `ggplot2` does not permit rotation for the `geom_label`
function used to produce these labels.

## 2.1 Surfaces and data perimeters

It may often be that user-provided data only extend across a portion of a map,
leaving a perimeter beyond the data boundary for which interpolation should not
be applied. `add_osm_surface` has a `bg` parameter specifying a background
colour for objects beyond the perimeter of the data surface. Passing this
parameter to `add_osm_surface` causes objects beyond the data perimeter to be
coloured within this 'background' colour.

To illustrate, trim the above data to within a circular range of the centre of
the map.
```{r}
d <- sqrt ((dat$x - mean (dat$x)) ^ 2 + (dat$y - mean (dat$y)) ^ 2)
range (d)
```
Remove from `dat` all rows translating to `d>0.01`:
```{r map20, eval = FALSE}
dat <- dat [which (d < 0.01), ]
cols_adj <- adjust_colours (heat.colors (30), -0.2)

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30),
                   bg = "gray40") %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = cols_adj,
                   size = c (1.5, 0.5),
                   bg = "gray70") %>%
  print_osm_map ()
```
```{r map20-print, echo = FALSE}
dat <- dat [which (d < 0.01), ]
cols_adj <- adjust_colours (heat.colors (30), -0.2)

osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_B,
                   dat = dat,
                   cols = heat.colors (30),
                   bg = "gray40") %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = cols_adj,
                   size = c (1.5, 0.5),
                   bg = "gray70") %>%
  print_osm_map (filename = "map_b20.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b20.png)

(The perimeter is irregular because of the positions of the points in `dat`.)

## 2.3 Further control of surface appearance

The final `add_osm_surface` call in the above code (for `dat_HP`) illustrates
additional parameters that may be passed for further control of map appearance.
In this case, the two `size` parameters control the size of the lines within the
data surface and beyond its perimeter. Single values may also be passed, in
which case they determine the width of lines in both cases. One or two `shape`
parameters may also be passed, with these also determining the shapes of
`SpatialPoints`, as illustrated in the next example, which overlays trees on the
map.

Both lines and points use the same `ggplot2` colour gradient, and so adding the
second of these again generates an error and means that the actual colour scheme
will be determined by the final call to add either lines or points.

```{r, eval = FALSE}
dat_T <- extract_osm_objects (key = "tree", bbox = bbox)
```
```{r, echo = FALSE}
dat_T <- london$dat_T
```

```{r map21, eval = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = terrain.colors (30),
                   size = c (1.5, 0.5),
                   bg = "gray70") %>%
  add_osm_surface (dat_H,
                   dat = dat,
                   cols = terrain.colors (30),
                   size = c (1, 0.5),
                   bg = "gray70") %>%
  add_osm_surface (dat_T,
                   dat = dat,
                   cols = heat.colors (30),
                   bg = "lawngreen",
                   size = c(3, 2),
                   shape = c(8, 1)) %>%
  print_osm_map ()
```
```{r map21-print, echo = FALSE}
osm_basemap (bbox = bbox,
             bg = "gray20") %>%
  add_osm_surface (dat_HP,
                   dat = dat,
                   cols = terrain.colors (30),
                   size = c (1.5, 0.5),
                   bg = "gray70") %>%
  add_osm_surface (dat_H,
                   dat = dat,
                   cols = terrain.colors (30),
                   size = c (1, 0.5),
                   bg = "gray70") %>%
  add_osm_surface (dat_T,
                   dat = dat,
                   cols = heat.colors (30),
                   bg = "lawngreen",
                   size = c(3, 2),
                   shape = c(8, 1)) %>%
  print_osm_map (filename = "map_b21.png",
                 width = 600,
                 units = "px",
                 dpi = map_dpi)
```
![](map_b21.png)

The first two colour specifications (`terrain.colors`) have been ignored, and
all added items are coloured according to the final value of `heat.colors (30)`.
Other aspects such as line sizes and point shapes are nevertheless respected.
---
title: "Rendering Ocean"
author: "Richard Beare"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        #number_sections: true
        theme: flatly
vignette: >
  %\VignetteIndexEntry{Maps with Ocean}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# 1 Introduction

It may often be desirable to separately colour two or more portions of a map
that are separated by a line object. This is not directly possible because only
polygons can be filled with colour, not lines. The `osm_line2poly()` function
comes to the rescue here by converting a line to a polygon surrounding a given
plotting region. The classic example where this arises is with coastlines. These
are always represented in OpenStreetMap as line objects, preventing any ability
to simply colour land and ocean separately.

This vignette illustrates the general principles of the `osm_line2poly()`
function, along with several ancilliary issues such as plotting coastal islands.
Although this functionality has been primarily developed with coastlines in
mind, the `osm_line2poly()` function has been designed in a suffciiently general
manner to be readily adaptible to any other cases where such line-to-polygon
conversion may be desirable.

This vignette explores the example of the coastline around Greater Melbourne,
Australia, first demonstrating how to extract coastline and convert it to land
and sea polygons, and then demonstrating how these areas may be delinated on a
plot.

# 2 Data Extraction and Conversion

We use [`osmdata`](https://cran.r-project.org/package=osmdata) to extract the
coastline within the bounding box of Greater Melbourne.
```{r, echo = FALSE}
map_dpi <- 72 # dpi res for all maps
fetch_osm <- FALSE
```
```{r GMFuncs, message=FALSE, eval = fetch_osm}
library (osmplotr)
library (osmdata)
library (magrittr)

bbox <- osmdata::getbb ("greater melbourne, australia")
coast <- opq (bbox = bbox) %>%
    add_osm_feature (key = "natural", value = "coastline") %>%
    osmdata_sf (quiet = FALSE)
```
This coastline object consists of several types of structure
```{r, eval = FALSE}
coast
```
```{r, echo = FALSE}
message (paste0 ("Object of class 'osmdata' with:\n",
"                 $bbox : -38.49937,144.44405,-37.40175,146.1925\n",
"        $overpass_call : The call submitted to the overpass API\n",
"            $timestamp : [ Thurs 5 Oct 2017 10:23:18 ]\n",
"           $osm_points : 'sf' Simple Features Collection with 13635 points\n",
"            $osm_lines : 'sf' Simple Features Collection with 73 linestrings\n",
"         $osm_polygons : 'sf' Simple Features Collection with 12 polygons\n",
"       $osm_multilines : 'sf' Simple Features Collection with 0 multilinestrings\n",
"    $osm_multipolygons : 'sf' Simple Features Collection with 0 multipolygons"))
```
Because OpenStreetMap represents coastline as line objects, all coastline data
is contained within the `$osm_lines` object. The `osm_line2poly()` function can
then converte these lines to polygons which can be used to plot filled areas.
```{r, eval = fetch_osm}
coast_poly <- osm_line2poly (coast$osm_lines, bbox)
names(coast_poly)
```
```{r, echo = FALSE}
c ("sea", "land", "islands")
```
Note that reflecting its envisioned primary usage, the function always returns
objects named `"sea"`, `"land"`, and `"islands"`. For usages other than
coastline, these names will of course reflect other kinds of object.  The
`"islands"` item contains any polygons which are separate to those originally
queried. Each item of this list is an `sf::data.frame` object:
```{r, eval = FALSE}
class (coast_poly$sea)
```
```{r, echo = FALSE}
c ("sf", "data.frame")
```

# 3 Plotting

The `list` items returned by `osm_line2poly()` may then be used to provide a map
background which distinguishes ocean from land. Here we first colour the entire
map using the background colour for the ocean, and overlay the land and island
polygons on top of that.

```{r, eval = fetch_osm}
map <- osm_basemap (bbox = bbox, bg = "cadetblue2") %>%
  add_osm_objects (coast_poly$land, col = "lightyellow1") %>%
  add_osm_objects (coast_poly$islands, col="orange") %>%
  add_osm_objects (coast$osm_polygons, col="purple", border = "black") %>%
  add_osm_objects (coast$osm_lines, col="black") %>%
  print_osm_map ()
```
```{r, echo=FALSE, eval = fetch_osm}
print_osm_map (map, filename = 'melb_a1.png', width = 600,
               units = 'px', dpi = map_dpi)
```

![](melb_a1.png)

The gaudy colours differentiate the source of polygons. Purple islands were
returned by the original osm query, while the orange ones were constructed from
fragments by `osm_line2poly`.

# Further Demonstrations

The `osm_line2poly()` function works by identifying lines which extend at at
least two points beyond a given bounding box. For coastline, OpenStreetMap is
designed so that land always lies to the left side in the direction of the
line, enabling water and land to be systematically distinguished. The following
test cases demonstrate the reliability of this distinction.

```{r, echo=FALSE}
  getCoast <- function(bbox)
  {
    qry <- opq(bbox)
    qry <- add_osm_feature(qry, key = "natural", value = "coastline")
    return(osmdata_sf(qry))
  }
  testPlot <- function(coast, bbox)
  {
    if (!dev.cur()) dev.off()
    map <- osm_basemap(bbox=bbox)
    map <- add_osm_objects(map, coast$osm_lines)
    print_osm_map(map)
  }
  testPlotPoly <- function(coast, bbox, fname)
  {
    ## trouble doing this check properly on Travis
    if (nrow(coast$osm_lines) > 0) {
      coastp <- osm_line2poly(coast$osm_lines, bbox=bbox)
      map <- osm_basemap(bbox=bbox)
      map <- add_osm_objects(map, coastp$sea, col='cadetblue2')
      map <- add_osm_objects(map, coastp$land, col='sienna2')
      print_osm_map(map,filename = fname, width = 200,
                    units = 'px', dpi = map_dpi)
    } else {
      warning("osm query probably failed - not plotting")
      invisible(NULL)
    }
    
  }
  
```
```{r, eval = fetch_osm}
test_plot <- function (bbox)
{
    dat <- opq (bbox) %>%
        add_osm_feature (key = "natural", value = "coastline") %>%
        osmdata_sf (quiet = FALSE)
    coast <- osm_line2poly (dat$osm_lines, bbox)
    osm_basemap (bbox = bbox) %>%
        add_osm_objects(coast$sea, col = 'cadetblue2') %>%
        add_osm_objects(coast$land, col = 'sienna2')
}
```
```{r, eval = fetch_osm, echo = FALSE}
test_plot <- function (bbox, filename, map_dpi)
{
    dat <- opq (bbox) %>%
        add_osm_feature (key = "natural", value = "coastline") %>%
        osmdata_sf (quiet = FALSE)
    coast <- osm_line2poly (dat$osm_lines, bbox)
    osm_basemap (bbox = bbox) %>%
        add_osm_objects(coast$sea, col = 'cadetblue2') %>%
        add_osm_objects(coast$land, col = 'sienna2') %>%
        print_osm_map (file = filename, width = 200,
                       units = "px", dpi = map_dpi)
}
```
Fetch the test data. A variable name with `WE` indicates that the coast enters the bounding box
on the western side and exits on the east. The land is on the left when following that path.
```{r, eval = fetch_osm}
    bbWE <- get_bbox (c(142.116906, -38.352713, 142.205162, -38.409661))
    coastWE <- getCoast(bbWE)

    bbEW <- get_bbox(c(144.603127, -38.104003, 144.685557, -38.135596))
    coastEW <- getCoast(bbEW)

    bbNS <- get_bbox(c(143.807998, -39.770986, 143.906494, -39.918643))
    coastNS <- getCoast(bbNS)

    bbSN <- get_bbox(c(144.073544, -39.854586, 144.149318, -39.960047))
    coastSN <- getCoast(bbSN)

    bbWW <- get_bbox(c(144.904865, -37.858295, 144.923679, -37.874367))
    coastWW <- getCoast(bbWW)

    bbEE <- get_bbox(c(144.643383, -38.294671, 144.692197, -38.336022))
    coastEE <- getCoast(bbEE)

    bbNN <- get_bbox(c(145.856321, -38.831642, 146.050920, -38.914031))
    coastNN <- getCoast(bbNN)

    bbSS <- get_bbox(c(146.363768, -38.770345, 146.486389, -38.837287))
    coastSS <- getCoast(bbSS)

    bbEN <- get_bbox(c(144.738212, -38.337690, 144.758053, -38.346966))
    coastEN <- getCoast(bbEN)

    bbEWWS <- get_bbox(c(144.693077, -38.307526, 144.729113, -38.343997 ))
    coastEWWS <- getCoast(bbEWWS)

    bbWS <- get_bbox(c(143.164906 ,-38.704885, 143.2075563, -38.7462058 ))
    coastWS <- getCoast(bbWS)

```

```{r, eval = fetch_osm}
testPlotPoly(coastWE, bbWE, "testWE.png")
```
![](testWE.png)
```{r, eval = fetch_osm}
testPlotPoly(coastEW, bbEW, "testEW.png")
```
![](testEW.png)

```{r, eval = fetch_osm}
testPlotPoly(coastNS, bbNS, "testNS.png")
```
![](testNS.png)

```{r, eval = fetch_osm}
testPlotPoly(coastSN, bbSN, "testSN.png")
```
![](testSN.png)

```{r, eval = fetch_osm}
testPlotPoly(coastWW, bbWW, "testWW.png")
```
![](testWW.png)

```{r, eval = fetch_osm}
testPlotPoly(coastEE, bbEE, "testEE.png")
```
![](testEE.png)

```{r, eval = fetch_osm}
testPlotPoly(coastNN, bbNN, "testNN.png")
```
![](testNN.png)

```{r, eval = fetch_osm}
testPlotPoly(coastSS, bbSS, "testSS.png")
```
![](testSS.png)

```{r, eval = fetch_osm}
testPlotPoly(coastEN, bbEN, "testEN.png")
```
![](testEN.png)

```{r, eval = fetch_osm}
testPlotPoly(coastEWWS, bbEWWS, "testEWWS.png")
```
![](testEWWS.png)

```{r, eval = fetch_osm}
testPlotPoly(coastWS, bbWS, "testWS.png")
```
![](testWS.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-osm-objects.R
\name{add_osm_objects}
\alias{add_osm_objects}
\title{add_osm_objects}
\usage{
add_osm_objects(map, obj, col = "gray40", border = NA, hcol, size, shape)
}
\arguments{
\item{map}{A \code{ggplot2} object to which the objects are to be added.}

\item{obj}{A spatial (\code{sp}) data frame of polygons, lines, or points,
typically as returned by \code{\link{extract_osm_objects}}.}

\item{col}{Colour of lines or points; fill colour of polygons.}

\item{border}{Border colour of polygons.}

\item{hcol}{(Multipolygons only) Vector of fill colours for holes}

\item{size}{Size argument passed to \code{ggplot2} (polygon, path, point)
functions: determines width of lines for (polygon, line), and sizes of
points.  Respective defaults are (0, 0.5, 0.5).}

\item{shape}{Shape of points or lines (the latter passed as \code{linetype});
see \code{\link[ggplot2]{shape}}.}
}
\value{
modified version of \code{map} to which objects have been added.
}
\description{
Adds layers of spatial objects (polygons, lines, or points generated by
\code{\link{extract_osm_objects}}) to a graphics object initialised with
\code{\link{osm_basemap}}.
}
\examples{
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
map <- osm_basemap (bbox = bbox, bg = "gray20")

\dontrun{
# The 'london' data used below were downloaded as:
dat_BNR <- extract_osm_objects (bbox = bbox, key = 'building',
                                value = '!residential')
dat_HP <- extract_osm_objects (bbox = bbox, key = 'highway',
                               value = 'primary')
dat_T <- extract_osm_objects (bbox = bbox, key = 'tree')
}
map <- add_osm_objects (map, obj = london$dat_BNR,
                        col = "gray40", border = "yellow")
map <- add_osm_objects (map, obj = london$dat_HP, col = "gray80",
                        size = 1, shape = 2)
map <- add_osm_objects (map, london$dat_T, col = "green",
                        size = 2, shape = 1)
print_osm_map (map)

# Polygons with different coloured borders
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_objects (map, obj = london$dat_HP, col = "gray80")
map <- add_osm_objects (map, london$dat_T, col = "green")
map <- add_osm_objects (map, obj = london$dat_BNR, col = "gray40",
                        border = "yellow", size = 0.5)
print_osm_map (map)
}
\seealso{
\code{\link{osm_basemap}}, \code{\link{extract_osm_objects}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osm-structures.R
\name{osm_structures}
\alias{osm_structures}
\title{osm_structures}
\usage{
osm_structures(
  structures = c("building", "amenity", "waterway", "grass", "natural", "park",
    "highway", "boundary", "tree"),
  col_scheme = "dark"
)
}
\arguments{
\item{structures}{The vector of types of structures (defaults listed in
\code{\link{extract_osm_objects}}).}

\item{col_scheme}{Colour scheme for the plot (current options include
\code{dark} and \code{light}).}
}
\value{
\code{data.frame} of structures, \code{key-value} pairs,
corresponding prefixes, and colours.
}
\description{
For the given vector of structure types returns a \code{data.frame}
containing two columns of corresponding OpenStreetMap \code{key-value} pairs,
one column of unambiguous suffixes to be appended to the objects returned by
\code{\link{extract_osm_objects}}, and one column specifying colours. This
\code{data.frame} may be subsequently modified as desired, and ultimately
passed to \code{\link{make_osm_map}} to automate map production.
}
\examples{
# Default structures:
osm_structures ()
# user-defined structures:
structures <- c ("highway", "park", "ameniiy", "tree")
structs <- osm_structures (structures = structures, col_scheme = "light")
# make_osm_map returns potentially modified list of data
\dontrun{
dat <- make_osm_map (osm_data = london, structures = structs)
# map contains updated $osm_data and actual map in $map
print_osm_map (dat$map)
}
}
\seealso{
\code{\link{make_osm_map}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make-osm-map.R
\name{make_osm_map}
\alias{make_osm_map}
\title{make_osm_map}
\usage{
make_osm_map(
  bbox,
  osm_data,
  structures = osm_structures(),
  dat_prefix = "dat_"
)
}
\arguments{
\item{bbox}{The bounding box for the map.  A 2-by-2 matrix of 4 elements with
columns of min and max values, and rows of x and y values.  If \code{NULL},
\code{bbox} is taken from the largest extent of OSM objects in
\code{osm_data}.}

\item{osm_data}{A list of OSM objects as returned from
\code{\link{extract_osm_objects}}.  These objects may be included in the plot
without downloading. These should all be named with the stated
\code{dat_prefix} and have suffixes as given in \code{structures}.}

\item{structures}{A \code{data.frame} specifying types of OSM structures as
returned from \code{\link{osm_structures}}, and potentially modified to alter
lists of structures to be plotted, and their associated colours. Objects are
overlaid on plot according to the order given in \code{structures}.}

\item{dat_prefix}{Prefix for data structures (default \code{dat_}). Final
data structures are created by appending the suffixes from
\code{\link{osm_structures}}.}
}
\value{
List of two components:
\enumerate{
  \item List of OSM structures each as
     \code{Spatial(Points/Lines/Polygons)DataFrame} and appended to
     \code{osm_data} (which is \code{NULL} by default), and
  \item The \code{map} as a \code{ggplot2} object
}
}
\description{
Makes an entire OSM map for the given bbox using the submitted data, or by
downloading data if none submitted. This is a convenience function enabling
an entire map to be produced according to the graphical format specified with
the \code{structures} argument.
}
\section{Note}{

If \code{osm_data} is not given, then data will be downloaded, which can take
some time.  Progress is dumped to screen.
}

\examples{
structures <- c ("highway", "park")
structs <- osm_structures (structures = structures, col_scheme = "light")
# make_osm_map returns potentially modified list of data using the provided
# 'london' data for highways and parks.
dat <- make_osm_map (osm_data = london, structures = structs)
# or download data automatically using a defined bounding boox
bbox <- get_bbox (c(-0.15,51.5,-0.10,51.52))
\dontrun{
dat <- make_osm_map (bbox = bbox, structures = structs)
print_osm_map (dat$map)
}
}
\seealso{
\code{\link{osm_basemap}}, \code{\link{add_osm_objects}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osmplotr.R
\docType{package}
\name{osmplotr}
\alias{osmplotr}
\title{osmplotr.}
\description{
Produces customisable images of OpenStreetMap (OSM) data and enables data
visualisation using OSM objects.  Extracts data using the overpass API.
Contains the following functions, data, and vignettes.
}
\section{Data Functions}{

\itemize{
\item \code{\link{extract_osm_objects}}: Download arbitrary OSM objects
\item \code{\link{connect_highways}}: Returns points sequentially connecting
list of named highways
}
}

\section{Basic Plotting Functions (without data)}{

\itemize{
\item \code{\link{add_axes}}: Overlay longitudinal and latitudinal axes on
plot
\item \code{\link{add_osm_objects}}: Overlay arbitrary OSM objects
\item \code{\link{make_osm_map}}: Automate map production with structures
defined in \code{\link{osm_structures}}
\item \code{\link{osm_structures}}: Define structures and graphics schemes
for automating map production
\item \code{\link{osm_basemap}}: Initiate a \code{ggplot2} object for an OSM
map
\item \code{\link{print_osm_map}}: Print a map to specified graphics
device
}
}

\section{Advanced Plotting Functions (with data)}{

\itemize{
\item \code{\link{add_osm_groups}}: Overlay groups of objects using specified
colour scheme
\item \code{\link{add_osm_surface}}: Overlay data surface by interpolating
given data
\item \code{\link{add_colourbar}}: Overlay a scaled colourbar for data added
with \code{\link{add_osm_surface}}
}
}

\section{Colour Manipulation Functions}{

\itemize{
\item \code{\link{adjust_colours}}: Lighted or darken given colours by
specified amount
\item \code{\link{colour_mat}}: Generate continuous 2D spatial matrix of
colours
}
}

\section{Other Functions}{

\itemize{
\item \code{\link{get_bbox}}: return bounding box from input vector
}
}

\section{Data}{

\itemize{
\item \code{\link{london}}: OSM Data from a small portion of central London
}
}

\section{Vignettes}{

\itemize{
\item \code{basic-maps}: Describes basics of downloading data and making
custom maps
\item \code{data-maps}: Describes how map elements can be coloured according
to user-provided data, whether categorical or continuous
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/line2poly.R
\name{osm_line2poly}
\alias{osm_line2poly}
\title{osm_line2poly}
\usage{
osm_line2poly(obj, bbox)
}
\arguments{
\item{obj}{A Simple Features (\code{sf}) data frame of lines, typically as
returned by \code{\link{extract_osm_objects}}, or by
\code{osmdata::osmdata_sf}.}

\item{bbox}{bounding box (Latitude-longitude range) to be plotted.  A 2-by-2
matrix of 4 elements with columns of min and max values, and rows of x and y
values. Can also be an object of class \code{sf}, for example as returned
from \code{extract_osm_objects} or the \code{osmdata} package, in which case
the bounding box will be extracted from the object coordinates.}
}
\value{
A list of three Simple Features (\code{sf}) data frames, labelled sea
islands and land.
}
\description{
Converts \code{sf::sfc_LINSTRING} objects to polygons by connecting end
points around the given bounding box. This is particularly useful for
plotting water and land delineated by coastlines. Coastlines in OpenStreetMap
are lines, not polygons, and so there is no directly way to plot ocean water
distinct from land. This function enables that by connecting the end points
of coastline \code{LINESTRING} objects to form closed polygons.
}
\details{
This is a tricky problem for a number of reasons, and the current
implementation may not be correct, although it does successfully deal with a
few tough situations. Some of the issues are: an osm coastline query returns
a mixture of "ways" and polygons.

Polygons correspond to islands, but not all islands are polygons. A "way" is
a connected set of points with the land on the left. A piece of coastline in
a bounding box may consist of multiple ways, which need to be connected
together to create a polygon. Also, ways extend outside the query bounding
box, and may join other ways that enter the bounding box (e.g ends of a
peninsula). The degree to which this happens depends on the scale of the
bounding box. Coastlines may enter at any bounding box edge and exit at any
other, including the one they entered from.
}
\examples{
# This example uses the \code{osmdata} package to extract data from
# a named bounding box
\dontrun{
library (magrittr)
library (osmdata)
bb <- osmdata::getbb ("melbourne, australia")
coast <- extract_osm_objects (bbox = bb,
                              key = "natural",
                              value = "coastline",
                              return_type = "line")
coast <- osm_line2poly (coast, bbox = bb)
# The following map then colours in just the ocean:
map <- osm_basemap (bbox = bb) \%>\%
    add_osm_objects (coast$sea, col = "lightsteelblue") \%>\%
    print_osm_map ()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-bbox.R
\name{get_bbox}
\alias{get_bbox}
\title{get_bbox}
\usage{
get_bbox(latlon)
}
\arguments{
\item{latlon}{A vector of (longitude, latitude, longitude, latitude) values.}
}
\value{
A 2-by-2 matrix of 4 elements with columns of min and max values, and
rows of x and y values.
}
\description{
Converts a string of latitudes and longitudes into a square matrix to be
passed as a \code{bbox} argument (to \code{\link{extract_osm_objects}},
\code{\link{osm_basemap}}, or \code{\link{make_osm_map}}).
}
\examples{
bbox <- get_bbox (c (-0.15, 51.5, -0.1, 51.52))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-axes.R
\name{add_axes}
\alias{add_axes}
\title{add_axes}
\usage{
add_axes(
  map,
  colour = "black",
  pos = c(0.02, 0.03),
  alpha = 0.4,
  fontsize = 3,
  fontface,
  fontfamily,
  ...
)
}
\arguments{
\item{map}{A \code{ggplot2} object to which the axes are to be added.}

\item{colour}{Colour of axis (determines colour of all elements: lines,
ticks, and labels).}

\item{pos}{Positions of axes and labels relative to entire plot device.}

\item{alpha}{alpha value for semi-transparent background surrounding axes and
labels (lower values increase transparency).}

\item{fontsize}{Size of axis font (in \code{ggplot2} terms; default=3).}

\item{fontface}{Fontface for axis labels (1:4=plain,bold,italic,bold-italic).}

\item{fontfamily}{Family of axis font (for example, `\code{Times}').}

\item{...}{Mechanism to allow many parameters to be passed with alternative
names (\code{color} for \code{colour} and \code{xyz} for \code{fontxyz}.}
}
\value{
Modified version of \code{map} with axes added.
}
\description{
Adds axes to the internal region of an OSM plot.
}
\examples{
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_objects (map, london$dat_BNR, col = "gray40")
map <- add_axes (map)
print (map)

# Map items are added sequentially, so adding axes prior to objects will
# produce a different result.
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_axes (map)
map <- add_osm_objects (map, london$dat_BNR, col = "gray40")
print_osm_map (map)
}
\seealso{
\code{\link{osm_basemap}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connect-highways.R
\name{connect_highways}
\alias{connect_highways}
\title{connect_highways}
\usage{
connect_highways(highways, bbox, plot = FALSE)
}
\arguments{
\item{highways}{A vector of highway names passed directly to the Overpass
API. Wildcards and whitespaces are `.'; for other options see online help for
the overpass API.}

\item{bbox}{the bounding box for the map.  A 2-by-2 matrix of 4 elements with
columns of min and max values, and rows of x and y values.}

\item{plot}{If \code{TRUE}, then all OSM data for each highway is plotted and
the final cycle overlaid.}
}
\value{
A single set of \code{SpatialPoints} containing the lat-lon
coordinates of the cyclic line connecting all given streets.
}
\description{
Takes a list of highways names which must enclose an internal area, and
returns a \code{SpatialLines} object containing a sequence of OSM nodes which
cyclically connect all highways. Will fail if the streets do not form a
cycle.
}
\note{
\enumerate{
\item \code{connect_highways} is primarily intended to provide a means to
define boundaries of groups which can then be highlighted using
\code{\link{add_osm_groups}}.
\item This function can not be guaranteed failsafe owing both to the
inherently unpredictable nature of OpenStreetMap, as well as to the unknown
relationships between named highways. The \code{plot} option enables
problematic cases to be examined and hopefully resolved.  The function is
still experimental, so please help further improvements by reporting any
problems!
}
}
\examples{
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
\dontrun{
highways <- c ("Monmouth.St", "Short.?s.Gardens", "Endell.St", "Long.Acre",
               "Upper.Saint.Martin")
# Note that dots signify "anything", including whitespace and apostrophes,
# and that '?' denotes optional previous character and so here matches
# both "Shorts Gardens" and "Short's Gardens"
highways1 <- connect_highways (highways = highways, bbox = bbox, plot = TRUE)
highways <- c ("Endell.St", "High.Holborn", "Drury.Lane", "Long.Acre")
highways2 <- connect_highways (highways = highways, bbox = bbox, plot = TRUE)

# Use of 'connect_highways' to highlight a region on a map
map <- osm_basemap (bbox = bbox, bg = "gray20")
# dat_B <- extract_osm_data (key = "building",
#                            value = "!residential",
#                            bbox = bbox)
# Those data are part of 'osmplotr':
dat_BNR <- london$dat_BNR # Non-residential buildings
groups <- list (highways1, highways2)
map <- add_osm_groups (map, obj = dat_BNR, groups = groups,
                       cols = c("red", "blue"), bg = "gray40")
print_osm_map (map)
}
}
\seealso{
\code{\link{add_osm_groups}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-colourbar.R
\name{add_colourbar}
\alias{add_colourbar}
\title{add_colorbar}
\usage{
add_colourbar(
  map,
  barwidth = 0.02,
  barlength = 0.7,
  zlims,
  cols,
  vertical = TRUE,
  alpha = 0.4,
  text_col = "black",
  fontsize = 3,
  fontface,
  fontfamily,
  ...
)
}
\arguments{
\item{map}{A \code{ggplot2} object to which the colourbar is to be added.}

\item{barwidth}{Relative width of the bar (perpendicular to its direction),
either a single number giving distance from right or upper margin, or two
numbers giving left/right or lower/upper limits.}

\item{barlength}{Relative length of the bar (parallel to its direction),
either a single number giving total length of centred bar, or two numbers
giving lower/upper or left/right limits.}

\item{zlims}{Vector of (min,max) values for scale of colourbar. These should
be the values returned from \code{\link{add_osm_surface}}.}

\item{cols}{Vector of colours.}

\item{vertical}{If \code{FALSE}, colourbar is aligned horizontally instead of
default vertical alignment.}

\item{alpha}{Transparency level of region immediately surrounding colourbar,
including behind text. Lower values are more transparent.}

\item{text_col}{Colour of text, tick marks, and lines on colourbar.}

\item{fontsize}{Size of text labels (in \code{ggplot2} terms; default=3).}

\item{fontface}{Fontface for colourbar labels
(1:4=plain,bold,italic,bold-italic).}

\item{fontfamily}{Family of colourbar font (for example, `\code{Times}').}

\item{...}{Mechanism to allow many parameters to be passed with alternative
names (such as \code{xyz} for \code{fontxyz}).}
}
\value{
Modified version of \code{map} with colourbar added.
}
\description{
Adds a colourbar to an existing map. Intended to be used in combination with
\code{\link{add_osm_surface}}. At present, only plots on right side of map.
}
\examples{
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
map <- osm_basemap (bbox = bbox, bg = "gray20")
# Align volcano data to lat-lon range of bbox
dv <- dim (volcano)
x <- seq (bbox [1,1], bbox [1,2], length.out = dv [1])
y <- seq (bbox [2,1], bbox [2,2], length.out = dv [2])
dat <- data.frame (
                  x = rep (x, dv [2]),
                  y = rep (y, each = dv [1]),
                  z = as.numeric (volcano)
                  )
map <- add_osm_surface (map, obj = london$dat_BNR, dat = dat,
                        cols = heat.colors (30))
map <- add_axes (map)
# Note colours of colourbar can be artibrarily set, and need not equal those
# passed to 'add_osm_surface'
map <- add_colourbar (map, zlims = range (volcano), cols = heat.colors(100),
                      text_col = "black")
print_osm_map (map)

# Horizontal colourbar shifted away from margins:
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_surface (map, obj = london$dat_BNR, dat = dat,
                        cols = heat.colors (30))
map <- add_colourbar (map, zlims = range (volcano), cols = heat.colors(100),
                      barwidth = c(0.1,0.15), barlength = c(0.5, 0.9),
                      vertical = FALSE)
print_osm_map (map)
}
\seealso{
\code{\link{osm_basemap}}, \code{\link{add_osm_surface}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract-osm-objects.R
\name{extract_osm_objects}
\alias{extract_osm_objects}
\title{extract_osm_objects}
\usage{
extract_osm_objects(
  bbox,
  key,
  value,
  extra_pairs,
  return_type,
  sf = TRUE,
  geom_only = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{bbox}{the bounding box within which all key-value objects should be
downloaded.  A 2-by-2 matrix of 4 elements with columns of min and
max values, and rows of x and y values.}

\item{key}{OSM key to search for. Useful keys include \code{building},
\code{waterway}, \code{natural}, \code{grass}, \code{park}, \code{amenity},
\code{shop}, \code{boundary}, and \code{highway}. Others will be passed
directly to the overpass API and may not necessarily return results.}

\item{value}{OSM value to match to key. If \code{NULL}, all keys will be
returned.  Negation is specified by \code{!value}.}

\item{extra_pairs}{A list of additional \code{key-value} pairs to be passed
to the overpass API.}

\item{return_type}{If specified, force return of spatial (\code{point},
\code{line}, \code{polygon}, \code{multiline}, \code{multipolygon}) objects.
\code{return_type = 'line'} will, for example, always return a
SpatialLinesDataFrame. If not specified, defaults to 'sensible' values (for
example, \code{lines} for highways, \code{points} for trees, \code{polygons}
for buildings).}

\item{sf}{If \code{TRUE}, return Simple Features (\code{sf}) objects;
otherwise Spatial (\code{sp}) objects.}

\item{geom_only}{If \code{TRUE}, return only those OSM data describing the
geometric object; otherwise return all data describing each object.}

\item{quiet}{If \code{FALSE}, provides notification of progress.}
}
\value{
Either a \code{SpatialPointsDataFrame}, \code{SpatialLinesDataFrame},
or \code{SpatialPolygonsDataFrame}.
}
\description{
Downloads OSM XML objects and converts to \code{sp} objects
(\code{SpatialPointsDataFrame}, \code{SpatialLinesDataFrame}, or
\code{SpatialPolygonsDataFrame}).
}
\examples{
\dontrun{
bbox <- get_bbox (c(-0.13,51.50,-0.11,51.52))
dat_B <- extract_osm_objects (key = 'building', bbox = bbox)
dat_H <- extract_osm_objects (key = 'highway', bbox = bbox)
dat_BR <- extract_osm_objects (key = 'building',
                               value = 'residential',
                               bbox = bbox)
dat_HP <- extract_osm_objects (key = 'highway',
                               value = 'primary',
                               bbox = bbox)
dat_HNP <- extract_osm_objects (key = 'highway',
                                value = '!primary',
                                bbox = bbox)
extra_pairs <- c ('name', 'Royal.Festival.Hall')
dat <- extract_osm_objects (key = 'building', extra_pairs = extra_pairs,
                            bbox = bbox)
}
}
\seealso{
\code{\link{add_osm_objects}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colour-mat.R
\name{colour_mat}
\alias{colour_mat}
\title{colour_mat}
\usage{
colour_mat(cols, n = c(10, 10), rotate, plot = FALSE)
}
\arguments{
\item{cols}{vector of length >= 4 of colors (example, default = \code{rainbow
(4)}, or \code{RColorBrewer::brewer.pal (4, 'Set1')}).
\code{cols} are wrapped clockwise around the corners from top left to bottom
left.}

\item{n}{number of rows and columns of colour matrix (default = 10; if length
2, then dimensions of rectangle).}

\item{rotate}{rotates the entire colour matrix by the specified angle (in
degrees).}

\item{plot}{plots the colour matrix.}
}
\value{
\code{Matrix} of colours.
}
\description{
Generates a 2D matrix of graduated colours by interpolating between the given
colours specifying the four corners.
}
\examples{
cm <- colour_mat (n = 5, cols = rainbow(4), rotate = 90, plot = TRUE)

# 'colour_mat' is intended primarily for use in colouring groups added with
# 'add_osm_groups' using the 'colmat = TRUE' option:
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
# Generate random points to serve as group centres
set.seed (2)
ngroups <- 6
x <- bbox [1,1] + runif (ngroups) * diff (bbox [1,])
y <- bbox [2,1] + runif (ngroups) * diff (bbox [2,])
groups <- cbind (x, y)
groups <- apply (groups, 1, function (i)
                 sp::SpatialPoints (matrix (i, nrow = 1, ncol = 2)))
# plot a basemap and add groups
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_groups (map, obj = london$dat_BNR, group = groups,
                       cols = rainbow (4), colmat = TRUE, rotate = 90)
print_osm_map (map)
}
\seealso{
\code{\link{add_osm_groups}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osmplotr.R
\docType{data}
\name{london}
\alias{london}
\title{london}
\format{
A list of spatial objects
}
\description{
A list of \code{Simple Features} (\code{sf}) \code{data.frame} objects
containing OpenStreetMap polygons, lines, and points for various
OpenStreetMap structures in a small part of central London, U.K.  (\code{bbox
= -0.13, 51.51, -0.11, 51.52}). The list includes:
\enumerate{
 \item \code{dat_H}: 974 non-primary highways as linestrings
 \item \code{dat_HP}: 159 primary highways as linestrings
 \item \code{dat_BNR}: 1,716 non-residential buildings as polygons
 \item \code{dat_BR}: 43 residential buildings as polygons
 \item \code{dat_BC}: 67 commerical buildings as polygons
 \item \code{dat_A}: 372 amenities as polygons
 \item \code{dat_P}: 13 parks as polygons
 \item \code{dat_T}: 688 trees as points
 \item \code{dat_RFH}: 1 polygon representing Royal Festival Hall
 \item \code{dat_ST}: 1 polygon representing 150 Stamford Street
}
}
\details{
The vignette \code{basic-maps} details how these data were downloaded. Note
that these internal versions have had all descriptive data removed other than
their names, geometries, and their OSM identification numbers.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-osm-surface.R
\name{add_osm_surface}
\alias{add_osm_surface}
\title{add_osm_surface}
\usage{
add_osm_surface(
  map,
  obj,
  dat,
  method = "idw",
  grid_size = 100,
  cols = heat.colors(30),
  bg,
  size,
  shape
)
}
\arguments{
\item{map}{A \code{ggplot2} object to which the surface are to be added}

\item{obj}{An \code{sp} \code{SpatialPolygonsDataFrame} or
\code{SpatialLinesDataFrame} (list of polygons or lines) returned by
\code{\link{extract_osm_objects}}}

\item{dat}{A matrix or data frame of 3 columns (x, y, z), where (x, y) are
(longitude, latitude), and z are the values to be interpolated}

\item{method}{Either \code{idw} (Inverse Distance Weighting as
\code{spatstat.core::idw}; default), \code{Gaussian} for kernel smoothing (as
\code{spatstat.core::Smooth.ppp}), or any other value to avoid interpolation.
In this case, \code{dat} must be regularly spaced in \code{x} and \code{y}.}

\item{grid_size}{size of interpolation grid}

\item{cols}{Vector of colours for shading z-values (for example,
\code{terrain.colors (30)})}

\item{bg}{If specified, OSM objects outside the convex hull surrounding
\code{dat} are plotted in this colour, otherwise they are included in the
interpolation (which will generally be inaccurate for peripheral values)}

\item{size}{Size argument passed to \code{ggplot2} (polygon, path, point)
functions: determines width of lines for (polygon, line), and sizes of
points.  Respective defaults are (0, 0.5, 0.5). If \code{bg} is provided and
\code{size} has 2 elements, the second determines the \code{size} of the
background objects.}

\item{shape}{Shape of lines or points, for details of which see
\code{?ggplot::shape}. If \code{bg} is provided and \code{shape} has 2
elements, the second determines the \code{shape} of the background objects.}
}
\value{
modified version of \code{map} to which surface has been added
}
\description{
Adds a colour-coded surface of spatial objects (polygons, lines, or points
generated by \code{\link{extract_osm_objects}} to a graphics object
initialised with \code{\link{osm_basemap}}. The surface is spatially
interpolated between the values given in \code{dat}, which has to be a matrix
of \code{data.frame} of 3 columns (x, y, z), where (x,y) are (longitude,
latitude), and z are the values to be interpolated. Interpolation uses
\code{spatstat.core::Smooth.ppp}, which applies a Gaussian kernel smoother
optimised to the given data, and is effectively non-parametric.
}
\note{
Points beyond the spatial boundary of \code{dat} are included in the surface
if \code{bg} is not given. In such cases, values for these points may exceed
the range of provided data because the surface will be extrapolated beyond
its domain.  Actual plotted values are therefore restricted to the range of
given values, so any extrapolated points greater or less than the range of
\code{dat} are simply set to the respective maximum or minimum values. This
allows the limits of \code{dat} to be used precisely when adding colourbars
with \code{\link{add_colourbar}}.
}
\examples{
# Get some data
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
# dat_B <- extract_osm_objects (key = 'building', bbox = bbox)
# These data are also provided in
dat_B <- london$dat_BNR # actuall non-residential buildings
# Make a data surface across the map coordinates, and remove periphery
n <- 5
x <- seq (bbox [1,1], bbox [1,2], length.out = n)
y <- seq (bbox [2,1], bbox [2,2], length.out = n)
dat <- data.frame (
    x = as.vector (array (x, dim = c(n, n))),
    y = as.vector (t (array (y, dim = c(n, n)))),
    z = x * y
    )
\dontrun{
map <- osm_basemap (bbox = bbox, bg = 'gray20')
map <- add_osm_surface (map, dat_B, dat = dat, cols = heat.colors (30))
print_osm_map (map)
}

# If data do not cover the entire map region, then the peripheral remainder
# can be plotted by specifying the 'bg' colour. First remove periphery from
# 'dat':
d <- sqrt ((dat$x - mean (dat$x)) ^ 2 + (dat$y - mean (dat$y)) ^ 2)
dat <- dat [which (d < 0.01),]
\dontrun{
map <- osm_basemap (bbox = bbox, bg = 'gray20')
map <- add_osm_surface (map, dat_B, dat = dat,
                        cols = heat.colors (30), bg = 'gray40')
print_osm_map (map)
}

# Polygons and (lines/points) can be overlaid as data surfaces with different
# colour schemes.
# dat_HP <- extract_osm_objects (key = 'highway',
#                                value = 'primary',
#                                bbox = bbox)
# These data are also provided in
dat_HP <- london$dat_HP
cols <- adjust_colours (heat.colors (30), adj = -0.2) # darken by 20\%
\dontrun{
map <- add_osm_surface (map, dat_HP, dat, cols = cols,
                        bg = 'gray60', size = c(1.5,0.5))
print_osm_map (map)
}

# Adding multiple surfaces of either polygons or (lines/points) produces a
# 'ggplot2' warning, and forces the colour gradient to revert to the last
# given value.
dat_T <- london$dat_T # trees
\dontrun{
map <- osm_basemap (bbox = bbox, bg = 'gray20')
map <- add_osm_surface (map, dat_B, dat = dat,
                        cols = heat.colors (30), bg = 'gray40')
map <- add_osm_surface (map, dat_HP, dat,
                        cols = heat.colors (30), bg = 'gray60',
                        size = c(1.5,0.5))
map <- add_osm_surface (map, dat_T, dat, cols = topo.colors (30),
                        bg = 'gray70', size = c(5,2), shape = c(8, 1))
print_osm_map (map) # 'dat_HP' is in 'topo.colors' not 'heat.colors'
}

# Add axes and colourbar
\dontrun{
map <- add_axes (map)
map <- add_colourbar (map, cols = heat.colors (100), zlims = range (dat$z),
                      barwidth = c(0.02), barlength = c(0.6,0.99),
                      vertical = TRUE)
print_osm_map (map)
}
}
\seealso{
\code{\link{osm_basemap}}, \code{\link{add_colourbar}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-osm-groups.R
\name{add_osm_groups}
\alias{add_osm_groups}
\title{add_osm_groups}
\usage{
add_osm_groups(
  map,
  obj,
  groups,
  cols,
  bg,
  make_hull = FALSE,
  boundary = -1,
  size,
  shape,
  border_width = 1,
  colmat,
  rotate
)
}
\arguments{
\item{map}{A \code{ggplot2} object to which the grouped objects are to be
added.}

\item{obj}{An \code{sp} \code{SpatialPointsDataFrame},
\code{SpatialPolygonsDataFrame}, or \code{SpatialLinesDataFrame} (list of
polygons or lines) returned by \code{\link{extract_osm_objects}}.}

\item{groups}{A list of spatial points objects, each of which contains the
coordinates of points defining one group.}

\item{cols}{Either a vector of >= 4 colours passed to \code{colour_mat} (if
\code{colmat = TRUE}) to arrange as a 2-D map of visually distinct colours
(default uses \code{rainbow} colours), or (if \code{colmat = FALSE}), a
vector of the same length as groups specifying individual colours for each.}

\item{bg}{If given, then any objects not within groups are coloured this
colour, otherwise (if not given) they are assigned to nearest group and
coloured accordingly (\code{boundary} has no effect in this latter case).}

\item{make_hull}{Either a single boolean value or a vector of same length as
groups specifying whether convex hulls should be constructed around all
groups (\code{TRUE}), or whether the group already defines a hull (convex or
otherwise; \code{FALSE}).}

\item{boundary}{(negative, 0, positive) values define whether the boundary of
groups should (exclude, bisect, include) objects which straddle the precise
boundary. (Has no effect if \code{bg} is given).}

\item{size}{Size argument passed to \code{ggplot2} (polygon, path, point)
functions: determines width of lines for (polygon, line), and sizes of
points.  Respective defaults are (0, 0.5, 0.5).}

\item{shape}{Shape of points or lines (the latter passed as \code{linetype});
see \code{\link[ggplot2]{shape}}.}

\item{border_width}{If given, draws convex hull borders around entire groups
in same colours as groups (try values around 1-2).}

\item{colmat}{If \code{TRUE} generates colours according to
\code{colour_mat}, otherwise the colours of groups are specified directly by
the vector of \code{cols}.}

\item{rotate}{Passed to \code{colour_mat} to rotate colours by the specified
number of degrees clockwise.}
}
\value{
Modified version of \code{map} with groups added.
}
\description{
Plots spatially distinct groups of OSM objects in different colours.
}
\section{Note}{

Any group that is entirely contained within any other group is assumed to
represent a hole, such that points internal to the smaller contained group
are *excluded* from the group, while those outside the smaller yet inside the
bigger group are included.
}

\examples{
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
# Download data using 'extract_osm_objects'
\dontrun{
dat_HP <- extract_osm_objects (key = 'highway',
                               value = 'primary',
                               bbox = bbox)
dat_T <- extract_osm_objects (key = 'tree', bbox = bbox)
dat_BNR <- extract_osm_objects (key = 'building', value = '!residential',
bbox = bbox)
}
# These data are also provided in
dat_HP <- london$dat_HP
dat_T <- london$dat_T
dat_BNR <- london$dat_BNR

# Define a function to easily generate a basemap
bmap <- function ()
{
    map <- osm_basemap (bbox = bbox, bg = "gray20")
    map <- add_osm_objects (map, dat_HP, col = "gray70", size = 1)
    add_osm_objects (map, dat_T, col = "green")
}

# Highlight a single region using all objects lying partially inside the
# boundary (via the boundary = 1 argument)
pts <- sp::SpatialPoints (cbind (c (-0.115, -0.125, -0.125, -0.115),
                                 c (51.505, 51.505, 51.515, 51.515)))
\dontrun{
dat_H <- extract_osm_objects (key = 'highway', bbox = bbox) # all highways
map <- bmap ()
map <- add_osm_groups (map, dat_BNR, groups = pts, cols = "gray90",
                       bg = "gray40", boundary = 1)
map <- add_osm_groups (map, dat_H, groups = pts, cols = "gray80",
                       bg = "gray30", boundary = 1)
print_osm_map (map)
}

# Generate random points to serve as group centres
set.seed (2)
ngroups <- 6
x <- bbox [1,1] + runif (ngroups) * diff (bbox [1,])
y <- bbox [2,1] + runif (ngroups) * diff (bbox [2,])
groups <- cbind (x, y)
groups <- apply (groups, 1, function (i)
                 sp::SpatialPoints (
                     matrix (i, nrow = 1, ncol = 2)))
# plot a basemap and add groups
map <- bmap ()
cols <- rainbow (length (groups))
\dontrun{
map <- add_osm_groups (map,
                       obj = london$dat_BNR,
                       group = groups,
                       cols = cols)
cols <- adjust_colours (cols, -0.2)
map <- add_osm_groups (map, obj = london$dat_H, groups = groups, cols = cols)
print_osm_map (map)

# Highlight convex hulls containing groups:
map <- bmap ()
map <- add_osm_groups (map,
                       obj = london$dat_BNR,
                       group = groups,
                       cols = cols,
                       border_width = 2)
print_osm_map (map)
}
}
\seealso{
\code{\link{colour_mat}}, \code{\link{add_osm_objects}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osm-basemap.R
\name{osm_basemap}
\alias{osm_basemap}
\title{osm_basemap}
\usage{
osm_basemap(bbox, structures, bg = "gray20")
}
\arguments{
\item{bbox}{bounding box (Latitude-longitude range) to be plotted.  A 2-by-2
matrix of 4 elements with columns of min and max values, and rows of x and y
values. Can also be an object of class \code{sf}, for example as returned
from \code{extract_osm_objects} or the \code{osmdata} package, in which case
the bounding box will be extracted from the object coordinates.}

\item{structures}{Data frame returned by \code{\link{osm_structures}} used
here to specify background colour of plot; if missing, the colour is
specified by \code{bg}.}

\item{bg}{Background colour of map (default = \code{gray20}) only if
\code{structs} not given).}
}
\value{
A \code{ggplot2} object containing the base \code{map}.
}
\description{
Generates a base OSM plot ready for polygon, line, and point objects to be
overlain with \code{\link{add_osm_objects}}.
}
\examples{
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_objects (map, london$dat_BNR, col = "gray40")
print_osm_map (map)
}
\seealso{
\code{\link{add_osm_objects}}, \code{\link{make_osm_map}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print-osm-map.R
\name{print_osm_map}
\alias{print_osm_map}
\title{print_osm_map}
\usage{
print_osm_map(
  map,
  width,
  height,
  filename,
  device,
  units = c("in", "cm", "mm", "px"),
  dpi = 300
)
}
\arguments{
\item{map}{The map to be printed; a \code{ggplot2} object produced by
\code{osmplotr}.}

\item{width}{Desired width of graphics device.}

\item{height}{Desired height of graphics device. Ignored if width specified.}

\item{filename}{Name of file to which map is to be printed.}

\item{device}{Type of graphics device (extracted from filename extension if
not explicitly provided).}

\item{units}{Units for height and width of graphics device.}

\item{dpi}{Resolution of graphics device (dots-per-inch).}
}
\description{
Prints an OSM map produced with \code{osmplotr} to a specified graphics
device.
}
\examples{
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
map <- osm_basemap (bbox = bbox, bg = "gray20")
map <- add_osm_objects (map, london$dat_BNR, col = "gray40")
print_osm_map (map, width = 7) # prints to screen device
\dontrun{
print_osm_map (map, file = "map.png", width = 500, units = "px")
}
}
\seealso{
\code{\link{osm_basemap}}, \code{\link{add_osm_objects}},
\code{\link{make_osm_map}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust-colours.R
\name{adjust_colours}
\alias{adjust_colours}
\title{adjust_colours}
\usage{
adjust_colours(cols, adj = 0, plot = FALSE)
}
\arguments{
\item{cols}{A vector of \code{R} colours (for allowable formats of which, see
\code{?col2rgb}).}

\item{adj}{A number between -1 and 1 determining how much to lighten
(positive values) or darken (negative values) the colours.}

\item{plot}{If \code{TRUE}, generates a plot to allow visual comparison of
original and adjusted colours.}
}
\value{
Corresponding vector of adjusted colours (as hexadecimal strings).
}
\description{
Adjusts a given colour by lightening or darkening it by the specified amount
(relative scale of -1 to 1).  Adjustments are made in RGB space, for
limitations of which see \code{?convertColor}
}
\examples{
cols <- adjust_colours (cols = heat.colors (10), adj = -0.2, plot = TRUE)

# 'adjust_colours' also offers an easy way to adjust the default colour
# schemes provided by 'osm_structures'. The following lines darken the
# highway colour of the 'light' colour scheme by 20\%
structures <- osm_structures (structures = c("building", "highway", "park"),
                              col_scheme = "light")
structures$cols [2] <- adjust_colours (structures$cols [2], adj = -0.2)
# Plot these structures:
bbox <- get_bbox (c (-0.13, 51.5, -0.11, 51.52))
\dontrun{
dat_B <- extract_osm_objects (key = "building", bbox = bbox)
dat_H <- extract_osm_objects (key = "highway", bbox = bbox)
dat_P <- extract_osm_objects (key = "park", bbox = bbox)
}
# These data are also included in the 'london' data of 'osmplotr'
osm_data <- list (dat_B = london$dat_BNR,
                  dat_H = london$dat_HP,
                  dat_P = london$dat_P)
dat <- make_osm_map (structures = structures,
                     osm_data = osm_data,
                     bbox = bbox)
print_osm_map (dat$map)
}
\seealso{
\code{\link{osm_structures}}, \code{?col2rgb}.
}
