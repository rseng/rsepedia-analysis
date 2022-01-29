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

