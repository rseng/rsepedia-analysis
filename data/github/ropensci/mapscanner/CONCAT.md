<!-- README.md is generated from README.Rmd. Please edit that file -->

# mapscanner <a href='https://docs.ropensci.org/mapscanner/'><img src='man/figures/logo.png' align="right" height=210 width=182/></a>

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/mapscanner/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/mapscanner/actions?query=workflow%3AR-CMD-check)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/mapscanner)](https://cran.r-project.org/package=mapscanner/)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mapscanner?color=orange)](https://cran.r-project.org/package=mapscanner)
[![codecov](https://codecov.io/gh/ropensci/mapscanner/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/mapscanner)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[![](https://badges.ropensci.org/330_status.svg)](https://github.com/ropensci/software-review/issues/330)

<!-- badges: end -->

## What does this package do for me?

`mapscanner` is an **R** package that enables lines drawn by hand on
maps to be converted to spatial objects. The package has two primary
functions: one for producing maps, and one for rectifying hand-drawn
lines to the coordinate system of the original map. The package is
intended for use in social surveys and similar endeavours in which
hand-drawn markings on maps need to be converted to spatial objects.
Maps can be either paper- or screen-based. Markings on paper maps need
to be scanned, photographed, or otherwise digitised, while maps with
screen-based markings need to be saved as `.png`-format images.

## Installation

The current stable version can be installed from CRAN with:

``` r
install.packages ("mapscanner")
```

Alternatively, the development version can be installed via [rOpenSci’s
r-universe](https://ropensci.r-universe.dev/) by running the following
prior to calling `install.packages()`:

``` r
options(repos = c(
                  ropensci = 'https://ropensci.r-universe.dev',
                  CRAN = 'https://cloud.r-project.org'))
```

The package can then be loaded for usage in a R session with

``` r
library (mapscanner)
```

## Usage

The package is designed to enable the following workflow:

1.  Generate a map with the
    [`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
    function, which automatically produces both `.pdf` and `.png`
    versions;

2.  Either print the `.pdf` version to use as desired in any kind of
    survey environment, or use either the `.pdf` or `.png` versions in
    digital form for screen-based surveys.

3.  Draw on the map;

4.  For paper maps, digitise the drawn-on (from here on, “modified”)
    map, converting it to either `.pdf` or `.png` format; and

5.  Rectify the modified version against the original via the
    [`ms_rectify_map()`](https://docs.ropensci.org/mapscanner/reference/ms_rectify_map.html)
    function, which distinguishes individual annotations, and converts
    each one to a spatial object able to be analysed in any desired
    manner.

### Practical tips

The `mapscanner` package is intended to aid a *practical* workflow, and
so a few practical tips may be recommended here to ensure best results:

1.  The original digital files generated with
    [`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
    are necessary to rectify subsequently drawn-on and scanned maps, and
    so must be retained at all times.
2.  Marks drawn on maps should be *coloured* – any black or grey
    markings will be ignored. This has the advantage that individual
    annotations *not* intended to be converted to spatial objects (such
    as unique identification or participant codes) may be made on maps
    in black or grey.
3.  For drawings of areas, best results will be achieved through
    ensuring that all lines form closed polygons. While the default
    `type = "hulls"` argument should work even when lines are not
    closed, the `type = "polygons"` argument will generally produce more
    accurate results, yet should only be used when all lines form closed
    polygons (see below for details on how these two differ).
4.  Digitised versions of paper maps should contain *white* borders, so
    do not, for example, photograph modified maps lying on dark
    surfaces. If maps are to be photographed, then best results can be
    achieved by simply placing them on a larger, enclosing sheet of
    white paper.

The following two sections describe the two primary functions of the
`mapscanner` package, corresponding to the two primary steps of
producing maps to be used in surveys (or other activities), and
rectifying modified maps against these originals in order to extract
spatial objects. The second of these sections also describes the kinds
of markings able to be recognised, and the kinds of spatial objects to
which these may be converted.

### Mapbox API tokens

Map generation with `mapscanner` requires a personal token or key from
[`mapbox`](https://www.mapbox.com/), which can be obtained by following
the links from
[https://docs.mapbox.com/api](https://docs.mapbox.com/api/#access-tokens-and-token-scopes/).
If you already have a token, the easiest way to use it with `mapscanner`
is to create (or edit) a file `~/.Renviron`, and insert a line,

``` bash
MAPBOX_TOKEN=<my_mapbox_token>
```

This will then be available every time you start R, without any need to
explicitly set the token each time you want to use the package. The
token may be given any unique name that includes “mapbox” (case
insensitive). Alternatively, if you wish to keep your token truly
private, and only use it for your current R session, you may load
`mapscanner`, and then run `set_mapbox_token(<my_mapbox_token>)`.

### Map generation

Having obtained and set a [`mapbox`](https://www.mapbox.com/) token as
described above, `mapscanner` may then be used to generate maps. The
package comes with a sample map of Omaha, Nebraska, USA, and one with
some red lines drawn on it: ![](man/figures/omaha-polygons.png)

That’s just a standard `png` image with no notion of geographical
coordinates. The original map was generated with

``` r
bbox <- rbind (c (-96.12923, -96.01011),
               c (41.26145, 41.32220)) # portion of omaha
ms_generate_map (bbox, max_tiles = 16L, mapname = "omaha")
```

    #> Successfully generated 'omaha.pdf' and 'omaha.png'

As indicated, the function generates a map in both `.pdf` and `.png`
formats. These files must be retained as the “master” maps against which
subsequently modified – drawn-over and scanned-in – versions will be
rectified.

### Map rectification

The magic within the `mapscanner` package happens via the [`RNiftyReg`
package](https://github.com/jonclayden/RNiftyReg), itself primarily
intended to align brain scans and other medical images, but which is
precisely the tool needed here. The package comes with two sample `.png`
images which can be used to demonstrate map rectification. In the
following code, `f_modified` is the image shown above, modified from the
original by drawing a red line around a particular region of Omaha.

``` r
f_orig <- system.file ("extdata", "omaha.png", package = "mapscanner")
f_mod <- system.file ("extdata", "omaha-polygons.png", package = "mapscanner")
res <- ms_rectify_map (f_orig, f_mod, type = "polygons")
#> ══ mapscanner ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
#> ✔ Image [/usr/lib/R/library/mapscanner/extdata/omaha.png] reduced in size by factor of 2
#> ❯ Rectifying the two maps ✔ Rectified the two maps  
#> ❯ Estimating optimal signal-to-noise threshold✔ Estimated optimal signal-to-noise threshold
#> ✔ Identified 2 objects
#> ❯ Converting to spatial format ✔ Converted to spatial format
res
#> Simple feature collection with 2 features and 0 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -96.11814 ymin: 41.26638 xmax: -96.02722 ymax: 41.30109
#> Geodetic CRS:  WGS 84
#>                         geometry
#> 1 POLYGON ((-96.11589 41.2663...
#> 2 POLYGON ((-96.03544 41.2927...
```

The rectification can take quite some time, during which [`RNiftyReg`
package](https://github.com/jonclayden/RNiftyReg) is constructing the
best transformation of the modified image back on to the original. The
result of `ms_rectify_map()` is a spatial object in
[`sf`](https://cran.r-project.org/package=sf)-format in which each drawn
component is represented as a separate polygon. Finally, we can plot the
result as an interactive map using packages like
[`mapview`](https://github.com/r-spatial/mapview) with the following
commands:

``` r
library (mapview)
mapview (res)
```

or [`mapdeck`](https://github.com/symbolixAU/mapdeck), which similarly
requires a mapbox token:

``` r
library (mapdeck)
set_token (Sys.getenv ("<my_mapbox_token>"))
mapdeck () %>%
    add_polygon (res, fill_colour = "#ffff00cc",
                 stroke_colour = "#ff0000", stroke_width = 20)
```

![](man/figures/leaflet-1.png)

And our hand-drawn lines shown above have been converted to standard
spatial objects able to be analysed in any desired way. See the [package
vignette](https://docs.ropensci.org/mapscanner/articles/mapscanner.html)
for more detail of what the `mapscanner` package can do.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the
[`allcontributors`
package](https://github.com/ropenscilabs/allcontributors) following the
[all-contributors](https://allcontributors.org/) specification.
Contributions of any kind are welcome!

### Code

<table>
<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars1.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/mdsumner">
<img src="https://avatars2.githubusercontent.com/u/4107631?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/commits?author=mdsumner">mdsumner</a>
</td>
<td align="center">
<a href="https://github.com/potterzot">
<img src="https://avatars0.githubusercontent.com/u/477294?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/commits?author=potterzot">potterzot</a>
</td>
</tr>
</table>

### Issues

<table>
<tr>
<td align="center">
<a href="https://github.com/ThomasDier">
<img src="https://avatars0.githubusercontent.com/u/42271539?u=4194553fedeadf5bf3a23efcc829dfd51615b549&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3AThomasDier">ThomasDier</a>
</td>
<td align="center">
<a href="https://github.com/dcooley">
<img src="https://avatars0.githubusercontent.com/u/8093396?u=2c8d9162f246d90d433034d212b29a19e0f245c1&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Adcooley">dcooley</a>
</td>
<td align="center">
<a href="https://github.com/SymbolixAU">
<img src="https://avatars2.githubusercontent.com/u/18344164?u=022e0d3bdcca3e224021bae842672bda12b599df&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3ASymbolixAU">SymbolixAU</a>
</td>
<td align="center">
<a href="https://github.com/khondula">
<img src="https://avatars3.githubusercontent.com/u/6106733?u=91bac57101b4e8047b2a96b8ad67437cc32e6144&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Akhondula">khondula</a>
</td>
<td align="center">
<a href="https://github.com/sckott">
<img src="https://avatars0.githubusercontent.com/u/577668?u=c54eb1ce08ff22365e094559a109a12437bdca40&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Asckott">sckott</a>
</td>
<td align="center">
<a href="https://github.com/tomroh">
<img src="https://avatars1.githubusercontent.com/u/6668593?u=9e585864a75453fb972c50d461f58aae35d65fac&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Atomroh">tomroh</a>
</td>
<td align="center">
<a href="https://github.com/stefaniebutland">
<img src="https://avatars1.githubusercontent.com/u/11927811?u=ea2b36cbdc6c1d4b5cd9231b397c03998d730626&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Astefaniebutland">stefaniebutland</a>
</td>
</tr>
<tr>
<td align="center">
<a href="https://github.com/SantoshSrinivas79">
<img src="https://avatars0.githubusercontent.com/u/1036163?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3ASantoshSrinivas79">SantoshSrinivas79</a>
</td>
<td align="center">
<a href="https://github.com/asitemade4u">
<img src="https://avatars0.githubusercontent.com/u/11460106?u=2d81ef9fbaa1ebcd7b9bcf65cd4894e85850d68e&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Aasitemade4u">asitemade4u</a>
</td>
</tr>
</table>
<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
0.0.6
===================

Minor changes:

- Fix initial submission which failed on solaris due to s2 geometry issues
  (unrelated to this pkg.


0.0.5
===================

First CRAN release
# Contributing to mapscanner

## Opening issues

The easiest way to note any behavioural curiosities or to request any new
features is by opening a [github issue](https://github.com/ropensci/mapscanner/issues).


## Development guidelines

If you'd like to contribute changes to `mapscanner`, we use [the GitHub
flow](https://guides.github.com/introduction/flow/index.html) for proposing,
submitting, reviewing, and accepting changes. If you haven't done this before,
there's a nice overview of git [here](http://r-pkgs.had.co.nz/git.html), as well
as best practices for submitting pull requests
[here](http://r-pkgs.had.co.nz/git.html#pr-make).

The `mapscanner` coding style diverges somewhat from [this commonly used R style
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
permeates `mapscanner` code. Words of text are separated by whitespace, and so
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

We want to encourage a warm, welcoming, and safe environment for contributing
to this project. See the [code of
conduct](https://github.com/ropensci/mapscanner/blob/master/CODE_OF_CONDUCT.md)
for more information.
# CRAN notes for mapscanner_0.0.6 submission

This submission fixes the Solaris failure on previous submission, along with failing example and test on previous submission attempt. These failures are due to dependencies in external packages which rely on external system libraries. Those packages themselves skip tests of the functions used on Solaris machines, which is what this submission now also does. The submission was checked on Solaris via r-hub with no errors or notes.

The submission still generates NOTEs on some systems regarding large installed size, which slightly exceeds 5MB (6.4MB at most on apple-darwin17.0). The docs are slightly under half of this, yet are necessary as the package enables digital rectification of hand-drawn marks on maps. Including example maps as image files is therefore an essential part of the documentation. Every effort has been made to ensure numbers and sizes of included files are as small as possible. Any further reduction in file sizes renders images markedly less useful, and would decrease the ability of people to understand how the package is intended to be used.

The package has been checked on all environments listed below, and generates only the one note regarding the multiple licenses used here - one for the main package, and another + LICENSE file for internally bundled code.

## Test environments

GitHub actions:
* Linux: R-release, R-devel, R-oldrelease
* OSX: R-release
* Windows: R3.6, R4.0, R-devel

CRAN win-builder:
* R-oldrelease, R-release, R-devel

Package also checked using `Clang++ -Weverything`, and both local memory sanitzer and `rocker/r-devel-san` with clean results.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
```

# mapscanner <a href='https://docs.ropensci.org/mapscanner/'><img src='man/figures/logo.png' align="right" height=210 width=182/></a>

<!-- badges: start -->
[![R build
status](https://github.com/ropensci/mapscanner/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/mapscanner/actions?query=workflow%3AR-CMD-check)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/mapscanner)](https://cran.r-project.org/package=mapscanner/) 
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mapscanner?color=orange)](https://cran.r-project.org/package=mapscanner)
[![codecov](https://codecov.io/gh/ropensci/mapscanner/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/mapscanner)
[![Project Status: Concept](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[![](https://badges.ropensci.org/330_status.svg)](https://github.com/ropensci/software-review/issues/330)

<!-- badges: end -->


## What does this package do for me?

`mapscanner` is an **R** package that enables lines drawn by hand on maps to be
converted to spatial objects. The package has two primary functions: one for
producing maps, and one for rectifying hand-drawn lines to the coordinate
system of the original map. The package is intended for use in social surveys
and similar endeavours in which hand-drawn markings on maps need to be
converted to spatial objects. Maps can be either paper- or screen-based.
Markings on paper maps need to be scanned, photographed, or otherwise
digitised, while maps with screen-based markings need to be saved as
`.png`-format images.

## Installation

The current stable version can be installed from CRAN with:
```{r cran-inst, eval = FALSE}
install.packages ("mapscanner")
```

Alternatively, the development version can be installed via [rOpenSci's
r-universe](https://ropensci.r-universe.dev/) by running the following prior to
calling `install.packages()`:
```{r r-univ, eval = FALSE}
options(repos = c(
                  ropensci = 'https://ropensci.r-universe.dev',
                  CRAN = 'https://cloud.r-project.org'))
```

The package can then be loaded for usage in a R session with
```{r library}
library (mapscanner)
```

## Usage

The package is designed to enable the following workflow:

1. Generate a map with the
   [`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
   function, which automatically produces both `.pdf` and `.png` versions;

2. Either print the `.pdf` version to use as desired in any kind of survey
   environment, or use either the `.pdf` or `.png` versions in digital form for
   screen-based surveys.

3. Draw on the map;

4. For paper maps, digitise the drawn-on (from here on, "modified") map,
   converting it to either `.pdf` or `.png` format; and

5. Rectify the modified version against the original via the
   [`ms_rectify_map()`](https://docs.ropensci.org/mapscanner/reference/ms_rectify_map.html)
   function, which distinguishes individual annotations,
   and converts each one to a spatial object able to be analysed in any desired
   manner.

### Practical tips

The `mapscanner` package is intended to aid a *practical* workflow, and so
a few practical tips may be recommended here to ensure best results:

1. The original digital files generated with
   [`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
   are necessary to rectify subsequently drawn-on and scanned maps, and so must
   be retained at all times.
2. Marks drawn on maps should be *coloured* -- any black or grey markings will
   be ignored. This has the advantage that individual annotations *not*
   intended to be converted to spatial objects (such as unique identification
   or participant codes) may be made on maps in black or grey.
3. For drawings of areas, best results will be achieved through ensuring that
   all lines form closed polygons. While the default `type = "hulls"` argument
   should work even when lines are not closed, the `type = "polygons"` argument
   will generally produce more accurate results, yet should only be used when
   all lines form closed polygons (see below for details on how these two differ).
4. Digitised versions of paper maps should contain *white* borders, so do not,
   for example, photograph modified maps lying on dark surfaces. If maps are to
   be photographed, then best results can be achieved by simply placing them on
   a larger, enclosing sheet of white paper.

The following two sections describe the two primary functions of the
`mapscanner` package, corresponding to the two primary steps of producing maps
to be used in surveys (or other activities), and rectifying modified maps
against these originals in order to extract spatial objects. The second of
these sections also describes the kinds of markings able to be recognised, and
the kinds of spatial objects to which these may be converted.


### Mapbox API tokens

Map generation with `mapscanner` requires a personal token or key from
[`mapbox`](https://www.mapbox.com/), which can be obtained by following the
links from
[https://docs.mapbox.com/api](https://docs.mapbox.com/api/#access-tokens-and-token-scopes/).
If you already have a token, the easiest way to use it with `mapscanner` is to
create (or edit) a file `~/.Renviron`, and insert a line,

``` bash
MAPBOX_TOKEN=<my_mapbox_token>
```
This will then be available every time you start R, without any need to
explicitly set the token each time you want to use the package. The token may
be given any unique name that includes "mapbox" (case insensitive).
Alternatively, if you wish to keep your token truly private, and only use it
for your current R session, you may load `mapscanner`, and then run
`set_mapbox_token(<my_mapbox_token>)`.

### Map generation

Having obtained and set a [`mapbox`](https://www.mapbox.com/) token as
described above, `mapscanner` may then be used to generate maps. The package
comes with a sample map of Omaha, Nebraska, USA, and one with some red lines
drawn on it: ![](man/figures/omaha-polygons.png)

That's just a standard `png` image with no notion of geographical coordinates.
The original map was generated with

```{r omaha-fakey, eval = FALSE}
bbox <- rbind (c (-96.12923, -96.01011),
               c (41.26145, 41.32220)) # portion of omaha
ms_generate_map (bbox, max_tiles = 16L, mapname = "omaha")
```
```{r omaha, echo = FALSE}
bbox <- rbind (c (-96.12923, -96.01011),
               c (41.26145, 41.32220))
message ("Successfully generated 'omaha.pdf' and 'omaha.png'")
```
As indicated, the function generates a map in both `.pdf` and `.png` formats.
These files must be retained as the "master" maps against which subsequently
modified -- drawn-over and scanned-in -- versions will be rectified. 


### Map rectification

The magic within the `mapscanner` package happens via the [`RNiftyReg`
package](https://github.com/jonclayden/RNiftyReg), itself primarily intended to
align brain scans and other medical images, but which is precisely the tool
needed here. The package comes with two sample `.png` images which can be used
to demonstrate map rectification. In the following code, `f_modified` is the
image shown above, modified from the original by drawing a red line around
a particular region of Omaha.
```{r scan_maps}
f_orig <- system.file ("extdata", "omaha.png", package = "mapscanner")
f_mod <- system.file ("extdata", "omaha-polygons.png", package = "mapscanner")
res <- ms_rectify_map (f_orig, f_mod, type = "polygons")
res
```
The rectification can take quite some time, during which [`RNiftyReg`
package](https://github.com/jonclayden/RNiftyReg) is constructing the best
transformation of the modified image back on to the original. The result of
`ms_rectify_map()` is a spatial object in
[`sf`](https://cran.r-project.org/package=sf)-format in which each drawn
component is represented as a separate polygon. Finally, we can plot the result
as an interactive map using packages like
[`mapview`](https://github.com/r-spatial/mapview) with the following commands:
```{r mapview, eval = FALSE}
library (mapview)
mapview (res)
```
or [`mapdeck`](https://github.com/symbolixAU/mapdeck), which similarly requires
a mapbox token:
```{r mapdeck, eval = FALSE}
library (mapdeck)
set_token (Sys.getenv ("<my_mapbox_token>"))
mapdeck () %>%
    add_polygon (res, fill_colour = "#ffff00cc",
                 stroke_colour = "#ff0000", stroke_width = 20)
```


```{r leaflet, echo = FALSE, eval = FALSE}
library(leaflet)

leaflet (res) %>%
    addPolygons (color = "#FF1111", weight = 1, opacity = 1.0,
                 fillOpacity = 0.5) %>%
    addProviderTiles ("CartoDB.Positron") %>%
    setView (lng = mean (bbox [1, ]),
             lat = mean (bbox [2, ]),
             zoom = 12)
```

![](man/figures/leaflet-1.png)

And our hand-drawn lines shown above have been converted to standard spatial
objects able to be analysed in any desired way. See the [package
vignette](https://docs.ropensci.org/mapscanner/articles/mapscanner.html) for
more detail of what the `mapscanner` package can do.


## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


## Contributors




<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the [`allcontributors` package](https://github.com/ropenscilabs/allcontributors) following the [all-contributors](https://allcontributors.org/) specification. Contributions of any kind are welcome!

### Code

<table>

<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars1.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/mdsumner">
<img src="https://avatars2.githubusercontent.com/u/4107631?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/commits?author=mdsumner">mdsumner</a>
</td>
<td align="center">
<a href="https://github.com/potterzot">
<img src="https://avatars0.githubusercontent.com/u/477294?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/commits?author=potterzot">potterzot</a>
</td>
</tr>

</table>


### Issues

<table>

<tr>
<td align="center">
<a href="https://github.com/ThomasDier">
<img src="https://avatars0.githubusercontent.com/u/42271539?u=4194553fedeadf5bf3a23efcc829dfd51615b549&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3AThomasDier">ThomasDier</a>
</td>
<td align="center">
<a href="https://github.com/dcooley">
<img src="https://avatars0.githubusercontent.com/u/8093396?u=2c8d9162f246d90d433034d212b29a19e0f245c1&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Adcooley">dcooley</a>
</td>
<td align="center">
<a href="https://github.com/SymbolixAU">
<img src="https://avatars2.githubusercontent.com/u/18344164?u=022e0d3bdcca3e224021bae842672bda12b599df&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3ASymbolixAU">SymbolixAU</a>
</td>
<td align="center">
<a href="https://github.com/khondula">
<img src="https://avatars3.githubusercontent.com/u/6106733?u=91bac57101b4e8047b2a96b8ad67437cc32e6144&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Akhondula">khondula</a>
</td>
<td align="center">
<a href="https://github.com/sckott">
<img src="https://avatars0.githubusercontent.com/u/577668?u=c54eb1ce08ff22365e094559a109a12437bdca40&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Asckott">sckott</a>
</td>
<td align="center">
<a href="https://github.com/tomroh">
<img src="https://avatars1.githubusercontent.com/u/6668593?u=9e585864a75453fb972c50d461f58aae35d65fac&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Atomroh">tomroh</a>
</td>
<td align="center">
<a href="https://github.com/stefaniebutland">
<img src="https://avatars1.githubusercontent.com/u/11927811?u=ea2b36cbdc6c1d4b5cd9231b397c03998d730626&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Astefaniebutland">stefaniebutland</a>
</td>
</tr>


<tr>
<td align="center">
<a href="https://github.com/SantoshSrinivas79">
<img src="https://avatars0.githubusercontent.com/u/1036163?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3ASantoshSrinivas79">SantoshSrinivas79</a>
</td>
<td align="center">
<a href="https://github.com/asitemade4u">
<img src="https://avatars0.githubusercontent.com/u/11460106?u=2d81ef9fbaa1ebcd7b9bcf65cd4894e85850d68e&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropensci/mapscanner/issues?q=is%3Aissue+author%3Aasitemade4u">asitemade4u</a>
</td>
</tr>

</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
---
title: "mapscanner styles"
author: "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document
vignette: >
  %\VignetteIndexEntry{mapscanner styles}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r pkg-load, echo = FALSE, message = FALSE}
library (mapscanner)
```

The `mapscanner` package is able to generate maps in the [three styles provided
by mapbox](https://docs.mapbox.com/api/maps/#styles). [The `ms_generate_map()`
function](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
has a `bw` parameter with a default of `TRUE` to generate black-and-white maps,
including conversion of coloured mapbox styles to black-and-white. The use of
coloured maps, generated with `bw = FALSE`, is not generally recommended in
`mapscanner`, and overdrawn markings are identified exclusively by colour, and
any markings in a colour which is also present in the underlying maps will
generally not be accurately identified.

Because illustrations of these map styles requires files too large to be placed
inside an R package, they are not included with this vignette. The styles can
nevertheless be easily examined from [the main mapbox styles
page](https://docs.mapbox.com/api/maps/styles/), or by clicking the hyperlinked
text below, and are:

1. `style = "light"` - the generally recommended style for `mapscanner`,
   ["designed to provide geographic context while highlighting the data on your
   analytics dashboard, data visualization, or data
   overlay](https://www.mapbox.com/maps/light/).
2. `style = "streets"`, which is,
   ["a comprehensive, general-purpose map that emphasizes accurate, legible
   styling of road and transit networks](https://www.mapbox.com/maps/streets/).
3. `style = "outdoors"`, which is, ["a general-purpose map with curated
   tilesets and specialized styling tailored to hiking, biking, and the most
   adventurous use cases"](https://www.mapbox.com/maps/outdoors/).

Alternatively, these different styles can be examined for any desired area with
the following code.

```{r map-gen, eval = FALSE}
bb <- osmdata::getbb ("<name of location>")
shrink <- 0.3 # shrink that bounding box to 30% size
bb <- t (apply (bb, 1, function (i)
                 mean (i) + c (-shrink, shrink) * diff (i) / 2))

ms_generate_map (bbox = bb,
                 max_tiles = 16L,
                 mapname = "map-light",
                 style = "light")
ms_generate_map (bbox = bb,
                 max_tiles = 16L,
                 mapname = "omaha-streets",
                 style = "streets",
                 bw = FALSE)
ms_generate_map (bbox = bb,
                 max_tiles = 16L,
                 mapname = "omaha-outdoors",
                 style = "outdoors",
                 bw = FALSE)
```
---
title: "mapscanner"
author: "Mark Padgham and Michael D. Sumner"
date: "`r Sys.Date()`"
output: 
    html_document
vignette: >
  %\VignetteIndexEntry{mapscanner}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r pkg-load, echo = FALSE, message = FALSE}
library (mapscanner)
library (ggplot2)
requireNamespace ("sf")
# vwidth/height are passed via mapview::mapshot to webshot::webshot,
# and default to (992,744), which bloats package size (issue #27).
red <- 2
vwidth <- 992 / red
vheight <- 744 / red
```

<!--
A NOTE regarding the size of this vignette, for which it was difficult to
produce a sufficiently small sized end result, largely because of the mapshot
images. The above values can be passed to mapshot to reduce those. The
following table shows sizes of resultant installed versions following various
options:

ggplot | mapshot } inst size | vignettes size
---  | --- | --- | ---
 -   |  -  | 2.2MB | 1.2MB
 yes |  -  | 2.4MB | 1.5MB
 yes |  1  | 2.9MB | 2.0MB (original mapshot image only)
 yes |  3  | 3.4MB | 2.5MB

 At present, the vignette is rendered as in the final row, but can potentially
 be reduced in size in the future according to this table. -->


## What does this package do for me?

`mapscanner` is an **R** package that enables lines drawn by hand on maps to be
converted to spatial objects. The package has two primary functions: one for
producing maps, and one for rectifying hand-drawn lines to the coordinate
system of the original map. The package is intended for use in social surveys
and similar endeavours in which hand-drawn markings on maps need to be
converted to spatial objects. Maps can be either paper- or screen-based.
Markings on paper maps need to be scanned, photographed, or otherwise
digitised, while maps with screen-based markings need to be saved as
`.png`- or `.pdf`-format images.


## How do I use it?

The package is designed to enable the following workflow:

1. Generate a map with the
   [`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
   function, which automatically produces both `.pdf` and `.png` versions;

2. Either print the `.pdf` version to use as desired in any kind of survey
   environment, or use either the `.pdf` or `.png` versions in digital form for
   screen-based surveys.

3. Draw on the map;

4. For paper maps, digitise the drawn-on (from here on, "modified") map,
   converting it to either `.pdf` or `.png` format; and

5. Rectify the modified version against the original via the
   [`ms_rectify_map()`](https://docs.ropensci.org/mapscanner/reference/ms_rectify_map.html)
   function, which distinguishes individual annotations,
   and converts each one to a spatial object able to be analysed in any desired
   manner.

### Practical tips

The `mapscanner` package is intended to aid a *practical* workflow, and so
a few practical tips may be recommended here to ensure best results:

1. The original digital files generated with
   [`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
   are necessary to rectify subsequently drawn-on and scanned maps, and so must
   be retained at all times.
2. Marks drawn on maps should be *coloured* -- any black or grey markings will
   be ignored. This has the advantage that individual annotations *not*
   intended to be converted to spatial objects (such as unique identification
   or participant codes) may be made on maps in black or grey.
3. For drawings of areas, best results will be achieved through ensuring that
   all lines form closed polygons. While the default `type = "hulls"` argument
   should work even when lines are not closed, the `type = "polygons"` argument
   will generally produce more accurate results, yet should only be used when
   all lines form closed polygons (see below for details on how these two differ).
4. Digitised versions of paper maps should contain *white* borders, so do not,
   for example, photograph modified maps lying on dark surfaces. If maps are to
   be photographed, then best results can be achieved by simply placing them on
   a larger, enclosing sheet of white paper.

The following two sections describe the two primary functions of the
`mapscanner` package, corresponding to the two primary steps of producing maps
to be used in surveys (or other activities), and rectifying modified maps
against these originals in order to extract spatial objects. The second of
these sections also describes the kinds of markings able to be recognised, and
the kinds of spatial objects to which these may be converted.

## Mapbox API tokens

Map generation with `mapscanner` requires a personal token or key from
[`mapbox`](https://www.mapbox.com/), which can be obtained by following the
links from
[https://docs.mapbox.com/api/](https://docs.mapbox.com/api/#access-tokens-and-token-scopes/).
If you already have a token, the easiest way to use it with `mapscanner` is to
create (or edit) a file `~/.Renviron`, and insert a line,

``` bash
MAPBOX_TOKEN=<my_mapbox_token>
```
This will then be available every time you start R, without any need to
explicitly set the token each time you want to use the package. The token may
be given any unique name that includes "mapbox" (case insensitive).
Alternatively, if you wish to keep your token truly private, and only use it
for your current R session, you may load `mapscanner`, and then run
`set_mapbox_token(<my_mapbox_token>)`.



## Map generation

Having obtained and set a [`mapbox`](https://www.mapbox.com/) token as
described above, the
[`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
function can be used to generate printable maps for a specified bounding box in
both `.pdf` and `.png` formats. Usage is a simple as,
```{r generate-fakey, eval = FALSE}
ms_generate_map ("chennai india", mapname = "chennai")
```
```{r generate, echo = FALSE}
message ("Successfully generated 'chennai.pdf' and 'chennai.png'")
```

The two generated maps are saved in the current working directory (`getwd()`).
To save maps in alternative locations, the `mapname` parameter can optionally
specify paths. To provide finer control over the scales of maps, precise
bounding boxes can also be submitted. To determine desired bounding boxes, we
recommend using the ['openstreetmap.org'
website](https://www.openstreetmap.org/), zooming to a desired area, then
clicking the "Export" button. A window will appear which includes the bounding
coordinates of the current screen. Even finer control can be gained by clicking
beneath this coordinate window on the line which says, "Manually select
a different area," which brings a drag-able rectangle onto the current screen.
The coordinates in the bounding box then simply need to be entered in to the
`bbox` parameter of
[`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
in the order (`xmin`, `ymin`, `xmax`, `ymax`) -- or anti-clockwise from the
left-hand coordinate.

The amount of detail in resultant maps is controlled by the `max_tiles`
argument, with larger values producing more detail, and resulting in larger
file sizes. The default value of `max_tiles = 16L` (where the `L` symbol tells
`R`to treat the value as an integer) should produce acceptable results for maps
extending across hundreds of metres to a few kilometres. Smaller-scale maps may
require higher values, and vice-versa. Map generation is relatively fast, and
so different values can be readily trialled.

Maps are generated in two formats, because the `.pdf` version will generally
be the most convenient for printing, while the `png` version should be retained
as the "master" copy against which to rectify subsequently scanned-in version.
Behind the scenes, the function downloads a series of vector map tiles from
[mapbox](https://www.mapbox.com/), and converts them to a `rasterBrick` object
from the [`raster` package](https://cran.r-project.org/package=raster). This
`rasterBrick` object is invisibly returned from the function:

```{r chennaiBrick-fakey, eval = FALSE}
x <- ms_generate_map ("chennai india", mapname = "chennai")
```
```{r, echo = FALSE}
message ("Successfully generated 'chennai.pdf' and 'chennai.png'")
```
```{r x, eval = FALSE}
x
```
```{r chennaiBrick, echo = FALSE}
x <- paste0 ("class      : RasterBrick\n",
             "dimensions : 1147, 562, 644614, 3  (nrow, ncol, ncell, nlayers)\n",
             "resolution : 38.21851, 38.21851  (x, y)\n",
             "extent     : 8921157, 8942635, 1442787, 1486624  (xmin, xmax, ymin, ymax)\n",
             "crs        : +proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs\n",
             "source     : memory\n",
             "names      : index_1.1, index_1.2, index_1.3\n",
             "min values :       107,       107,       107\n",
             "max values :       254,       254,       254")
message (x)
```
This `rasterBrick` object contains raster information for the three colour
channels of the image, and so may also be used for immediate viewing within
**R** with `raster::plotRGB(x)`.

Standard uses of the package should not need to explicitly access or modify
these data, but it is nevertheless possible to do so, and then use a
custom-modified object to produce the external `.pdf` and `.png` files by
submitting the `rasterBrick` object to
[`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html):
```{r chennaiBrick2, eval = FALSE}
ms_generate_map (raster_brick = x, mapname = "chennai")
```


## Map rectification

Having produced digital maps using the
[`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
function as described above, and having printed, variously drawn-on, and, for
paper maps, scanned the result back in to digital form, the package can then be
used to rectify the hand-drawn markings against the original map with the
[`ms_rectify_map()`](https://docs.ropensci.org/mapscanner/reference/ms_rectify_map.html)
function, which returns the drawn-on objects as spatial objects in [Simple
Features (`sf`)](https://cran.r-project.org/package=sf) format. The only
requirement is that the drawn-on objects are coloured; black or grey objects
will be ignored. As described above, this has the advantage that maps may be
annotated in ways not intended to be converted to spatial objects (such as
adding unique identification or participant codes), through simply providing
such annotations in grey or black.

The
[`ms_rectify_map()`](https://docs.ropensci.org/mapscanner/reference/ms_rectify_map.html)
function has two primary arguments, specifying the names (and locations) of the
original and modified map files -- in that order: `ms_rectify_map(original,
modified)`. These files should ideally be in `.png` formats, but will be
auto-converted from `.pdf` if needed. The package comes with two sample maps,
both in `.png` format. The first is the reference version needed for
rectification, while the second has two red lines drawn upon it:

<!--
<img src="omaha-polygons.png" width="600" height="400"/>
-->

```{r omaha-poly-png, echo = FALSE, out.width = "75%", eval = TRUE}
op <- system.file ("extdata", "omaha-polygons.png", package = "mapscanner")
f <- file.path (tempdir (), "omama-polygons.png")
chk <- file.copy (op, f)
knitr::include_graphics (f)
```


Converting the lines on this scanned image file is then as simple as:

```{r rectify, message = FALSE, eval = FALSE}
f_orig <- system.file ("extdata", "omaha.png", package = "mapscanner")
f_mod <- system.file ("extdata", "omaha-polygons.png", package = "mapscanner")
xy <- ms_rectify_map (f_orig, f_mod, nitems = 2)
xy
```
```{r rectify-output, echo = FALSE}
m <- c (
"Simple feature collection with 2 features and 0 fields",
"Geometry type: POLYGON",
"Dimension:     XY",
"Bounding box:  xmin: -96.11801 ymin: 41.26638 xmax: -96.02722 ymax: 41.30108",
"Geodetic CRS:  WGS 84",
"                        geometry",
"1 POLYGON ((-96.10692 41.2684...",
"2 POLYGON ((-96.0276 41.2964,...")
m <- paste0 (m, collapse = "\n")
message (m)
```

The result of
[`ms_rectify_map()`](https://docs.ropensci.org/mapscanner/reference/ms_rectify_map.html)
can be plotted using any standard option for plotting spatial data, such as
through online mapping packages such as `mapview`:
```{r mapshot-hulls, eval = FALSE, echo = FALSE}
xy$id <- seq (nrow (xy))
if (!file.exists ("mapshot-polys.png"))
{
    x <- mapview::mapview (xy)
    mapview::mapshot (x, file = "mapshot-hulls.png",
                      remove_controls = c ("layersControl", "zoomControl"),
                      vwidth = vwidth, vheight = vheight)
}
```
```{r mapview-fakey, eval = FALSE}
xy$id <- seq (nrow (xy))
mapview::mapview (xy)
```

<!--
![](mapshot-polys.png)
-->
```{r mapshot-hulls-png, echo = FALSE, out.width = "75%"}
knitr::include_graphics ("./mapshot-hulls.png")
#plot (xy)
```



### Types of map markings and types of spatial objects

The
[`ms_rectify_map()`](https://docs.ropensci.org/mapscanner/reference/ms_rectify_map.html)
function has an additional argument, `type`, which takes the following values:

1. `type = "hulls"` (the default), which returns convex or concave hulls around
   distinct sets of contiguously marked lines, regardless of whether those
   lines form closed polygons or not (see function help for details).
2. `type = "polygons"`, which returns the outlines traced around each
   individual drawn object. This tracing is pixel-based, resulting in polygons
   with one spatial point for each scanned pixel. This may generate spatial
   objects that are both overly large as well as visually pixillated. The
   function includes an additional `downsample` parameter which down-samples
   and smooths the resultant polygons by the specified multiple.
3. `type = "points"`, which returns single points (as geometric centroids) for
   each object. This is useful for identification of individual point locations
   regardless of the kinds of marks actually drawn on a map (dots, circles,
   crosses, or any shape, should all give equivalent results).

The type of `polygons` assumes -- and indeed requires -- that the drawn objects
are *closed* polygons (as illustrated in the first of the above figures), so
care must be taken to ensure this is in fact the case. Any lines that do not
form closed circles will not be appropriately translated. Algorithms for
extracting objects with `type = "polygons"`are fundamentally different from
`type = "hulls"`. The latter applies convex or concave-hull tracing algorithms,
while the former explicitly traces every individual pixel of a contiguous
object, and returns the external boundary comprised of the coordinates of all
pixels lying on that boundary. This will thus often produce more accurate and
detailed results, yet as mentioned should only be applied where markings form
strictly closed polygons. All other cases in which areal rather than
point-based results are desired should use the default `type = "hulls"`.
Examples include participants being asked to colour particular areas using any
desired kind of marks, enabling areal-filling scribbles can be converted to
polygons representing the outer boundaries.

The `"polygon"` and `"point"` types are illustrated in the following maps:
```{r type-polys, eval = FALSE}
f_orig <- system.file ("extdata", "omaha.png", package = "mapscanner")
f_mod <- system.file ("extdata", "omaha-polygons.png", package = "mapscanner")
xy <- ms_rectify_map (f_orig, f_mod, type = "polygons", quiet = TRUE)
```
```{r mapview-polys, eval = FALSE, echo = FALSE}
xy$id <- seq (nrow (xy))
if (!file.exists ("mapshot-polys.png"))
{
    x <- mapview::mapview (xy)
    mapview::mapshot (x, file = "mapshot-polys.png",
                      remove_controls = c ("layersControl", "zoomControl"),
                      vwidth = vwidth, vheight = vheight)
}
```
```{r mapview-polys-fakey, eval = FALSE}
xy$id <- seq (nrow (xy))
mapview::mapview (xy)
```

<!--
![](mapshot-polys.png)
-->
```{r mapshot-polys-png, echo = FALSE, out.width = "75%"}
knitr::include_graphics ("./mapshot-polys.png")
#plot (xy)
```



```{r type-points, eval = FALSE}
f_orig <- system.file ("extdata", "omaha.png", package = "mapscanner")
f_mod <- system.file ("extdata", "omaha-polygons.png", package = "mapscanner")
xy <- ms_rectify_map (f_orig, f_mod, type = "points", quiet = TRUE)
```
```{r mapview-points, eval = FALSE, echo = FALSE}
xy$id <- seq (nrow (xy))
if (!file.exists ("mapshot-points.png"))
{
    x <- mapview::mapview (xy)
    mapview::mapshot (x, file = "mapshot-points.png",
                      remove_controls = c ("layersControl", "zoomControl"),
                      vwidth = vwidth, vheight = vheight)
}
```
```{r mapview-points-fakey, eval = FALSE}
xy$id <- seq (nrow (xy))
mapview::mapview (xy)
```

<!--
![](mapshot-points.png)
-->
```{r mapshot-points-png, echo = FALSE, out.width = "75%"}
knitr::include_graphics ("./mapshot-points.png")
#plot (xy)
```




## Bonus Feature: Polygon Aggregation

Maps are typically used in social surveys to delineate participants'
understanding or perception of particular regions or areas. In such contexts,
surveys often result in numerous polygonal shapes representing different
perceptions of a particular region. The `mapscanner` package provides an
additional function, `ms_aggregate_polys()`, to aggregate polygons into
a single "heat map" containing vector outlines of aggregated polygons. Each
component of these aggregated polygons defines the region within which `n`
polygons overlap.

The function is now illustrated with a slightly more complicated version of the
example provided for `ms_aggregate_polys()`, starting by generating a series of
polygons as convex hulls surrounding random points.
```{r random-polys-fakey, eval = FALSE}
n <- 5 # number of polygons
polys <- lapply (seq (n), function (i) {
                     nxy <- 20 # number of points used to generate hull
                     xy <- matrix (runif (2 * nxy), ncol = 2)
                     h <- chull (xy)
                     sf::st_polygon (list (xy [c (h, h [1]), ]))
            })
polys <- sf::st_sf (n = seq (n), geometry = polys)
```
The `polys` object is then a Simple Features
([`sf`](https://cran.r-project.org/package=sf))
`data.frame` with `n` overlapping polygons, and an additional row, `n`, to
identify each polygon. The following lines then convert these to aggregated,
overlapping polygons, and plot the result:
```{r random-polys-fakey2, eval = FALSE}
aggr <- ms_aggregate_polys (polys)
polys$type <- "raw polygons"
aggr$type <- "aggregated polygons"
polys <- rbind (polys, aggr)
# Convert type to factor to order facets:
polys$type <- factor (polys$type, levels = c ("raw polygons", "aggregated polygons"))
library (ggplot2)
ggplot (polys, aes (fill = n)) + geom_sf () + facet_wrap (~type)
```
```{r random-polys, echo = FALSE, eval = TRUE}
set.seed (1)
n <- 5 # number of polygons
polys <- lapply (seq (n), function (i) {
                     nxy <- 20 # number of points used to generate hull
                     xy <- matrix (runif (2 * nxy), ncol = 2)
                     h <- chull (xy)
                     sf::st_polygon (list (xy [c (h, h [1]), ]))
            })
polys <- sf::st_sf (n = seq (n), geometry = polys)
aggr <- ms_aggregate_polys (polys)
polys$type <- "raw polygons"
aggr$type <- "aggregated polygons"
polys <- rbind (polys, aggr)
polys$type <- factor (polys$type, levels = c ("raw polygons", "aggregated polygons"))
theme_set (theme_minimal ())
ggplot (polys, aes (fill = n)) + geom_sf () + facet_wrap (~type)
```

The left panel of that figure shows the random polygons in raw form
successively overlaid upon one another. The right panel shows the
aggregated contours of successive overlap from 1 to 5. The object returned from
[`ms_aggregate_polys()`](https://docs.ropensci.org/mapscanner/reference/ms_aggregate_polys.html)
contains polygons ordered by level of aggregation (`n`), so the first entirely
encloses the second; the second encloses the third; and so on. Particular
contours can then be directly selected by filtering for desired values of `n`:
```{r aggr13, eval = TRUE}
ggplot (aggr [aggr$n %in% c (1, 3, 5), ], aes (fill = n)) + geom_sf ()
```

Polygon aggregation enables many interesting analyses to be performed, such as
relationships between aggregation level and area:
```{r aggr-area, eval = TRUE}
aggr$area <- sf::st_area (aggr)
ggplot (aggr, aes (x = n, y = area)) + geom_line (size = 2)
```

That result is of course (roughly) linear, because it was derived from random
data. In actual usage, results such as that are likely to generate direct
insight into consensus of opinion regarding how people understand particular
areas.



## How it works

(This section is not necessary for package usage, and merely provides detail
for those interested in how the process actually works.) `mapscanner` primarily
relies on the [`RNiftyReg` package](https://github.com/jonclayden/RNiftyReg) to
rectify the images. This package is itself primarily aimed at rectifying
medical scans, but also happens to be the perfect tool for the present
purposes. Being an image analysis software, the library requires image
objects, and not `pdf` files, which is why the
[`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html)
function produces both kinds of files - the `.pdf` for printing and the `.png`
for rectifying with the `RNiftyReg` package.

Rectification re-projects a scanned image back on to the coordinate system of
an original image. This coordinate system translates here in to a defined
bounding box (which may differ slightly from the values input into the
function, due to the cutting and stitching of the vector tiles). This bounding
box is embedded as meta-information in both the files produced by
[`ms_generate_map()`](https://docs.ropensci.org/mapscanner/reference/ms_generate_map.html);
in the `.pdf` as standard meta information accessible in **R** like this:

```{r pdfinfo-fakey, eval = FALSE}
pdftools::pdf_info ("chennai.pdf")$keys
```
```{r pdfinfo, echo = FALSE}
x <- list ("Title" = "EX8921118.44521949+1442748.9088827+8942673.68719591+1486623.76311839",
           "Producer" = "R 4.0.0",
           "Creator" = "R")
x
```

or in a terminal via `pdfinfo` (or non-linux equivalent), and embedded in the
`.png` file as comment, accessible like this:

```{r img-comment-fakey, eval = FALSE}
img <- magick::image_read ("chennai.png")
magick::image_comment (img)
```
```{r img-comment, echo = FALSE}
"EX8921118.44521949+1442748.9088827+8942673.68719591+1486623.76311839"
```
or in a terminal via `identify -verbose` command (itself part of the
`imagemagick` library which drives the [`magick` **R**
package](https://github.com/ropensci/magick)).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapscanner-package.R
\docType{package}
\name{mapscanner-package}
\alias{mapscanner}
\alias{mapscanner-package}
\title{mapscanner: Print Maps, Draw on Them, Scan Them Back in}
\description{
Enables preparation of maps to be printed and drawn on. Modified maps can then be scanned back in, and hand-drawn marks converted to spatial objects.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/mapscanner}
  \item Report bugs at \url{https://github.com/ropensci/mapscanner/issues}
}

}
\author{
\strong{Maintainer}: Mark Padgham \email{mark.padgham@email.com}

Authors:
\itemize{
  \item Michael D. Sumner \email{mdsumner@gmail.com}
}

Other contributors:
\itemize{
  \item Kelly Hondula (Kelly reviewed the package for rOpenSci, see https://github.com/ropensci/software-review/issues/330) [reviewer]
  \item Nicholas Potter (Nichola reviewed the package for rOpenSci, see https://github.com/ropensci/software-review/issues/330) [reviewer]
  \item Stanislaw Adaszewski (author of include concaveman-cpp code) [copyright holder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tokens.R
\name{set_mapbox_token}
\alias{set_mapbox_token}
\title{Set 'mapbox' token}
\usage{
set_mapbox_token(token)
}
\arguments{
\item{token}{Personal mapbox API token, obtained from
\url{https://docs.mapbox.com/api/#access-tokens-and-token-scopes}.}
}
\value{
\code{TRUE} if the token was able to be set; otherwise \code{FALSE}.
}
\description{
Set a mapbox token for use with the \link{ms_generate_map} function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_polys.R
\name{ms_aggregate_polys}
\alias{ms_aggregate_polys}
\title{Aggregate disparate polygons}
\usage{
ms_aggregate_polys(p)
}
\arguments{
\item{p}{input (multi-)polygons (assumed to be overlapping)}
}
\value{
Set of \pkg{sf}-format polygons with additional column, \code{n}, denoting
number of overlaps contributing to each of the resultant polygons.
}
\description{
Planar partition from disparate polygon inputs. Overlaps aggregate to \code{n}.
}
\details{
Input is a single simple features polygon data frame. No attribute data is
considered.
}
\examples{
g <- sf::st_sfc(list(sf::st_point(cbind(0, 0)),
                     sf::st_point(cbind(0, 1)),
                     sf::st_point(cbind(1, 0))))
pts <- sf::st_sf(a = 1:3,  geometry = g)
overlapping_polys <- sf::st_buffer(pts, 0.75)

## decompose and count space-filling from overlapping polygons
x <- ms_aggregate_polys(overlapping_polys)
plot(x)
\dontrun{
library(ggplot2)
ggplot(x) + geom_sf() + facet_wrap(~n)
}

library(sf)
set.seed(6)
pts <- expand.grid (x = 1:8, y = 1:10) \%>\% st_as_sf (coords = c("x", "y"))
xsf <- sf::st_buffer (pts, runif (nrow (pts), 0.2, 1.5))
\dontrun{
out <- ms_aggregate_polys (xsf)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate-maps.R
\name{ms_generate_map}
\alias{ms_generate_map}
\title{Generate maps for 'mapscanner' use}
\usage{
ms_generate_map(
  bbox,
  max_tiles = 16L,
  mapname = NULL,
  bw = TRUE,
  style = "light",
  raster_brick = NULL
)
}
\arguments{
\item{bbox}{Either a string specifying the location, or a numeric bounding
box as a single vector of (xmin, ymin, xmax, ymax), or a 2-by-2 matrix with
columns of (min, max) and rows of (x, y), respectively.}

\item{max_tiles}{Maximum number of tiles to use to create map}

\item{mapname}{Name of map to be produced, optionally including full path.
Extension will be ignored.}

\item{bw}{If \code{FALSE}, print maps in colour, otherwise black-and-white. Note
that the default \code{style = "light"} is monochrome, and that this parameter
only has effect for \code{style} values of \code{"streets"} or \code{"outdoors"}.}

\item{style}{The style of the map to be generated; one of 'light', 'streets',
or 'outdoors', rendered in black and white. See
\url{https://docs.mapbox.com/api/maps/#styles/} for examples.}

\item{raster_brick}{Instead of automatically downloading tiles within a given
\code{bbox}, a pre-downloaded \code{raster::rasterBrick} object may be submitted and
used to generate the \code{.pdf} and \code{.png} equivalents.}
}
\value{
Invisibly returns a \code{rasterBrick} object from the \pkg{raster}
package containing all data used to generate the map.
}
\description{
Generate a map image for a specified area or bounding box, by downloading
tiles from \url{https://www.mapbox.com/}. Map is automatically saved in both
\code{.pdf} and \code{.png} formats, by default in current working directory, or
alternative location when \code{mapname} includes the full path.
}
\examples{
\dontrun{
# code used to generate internal files for a portion of Omaha:
bb <- osmdata::getbb ("omaha nebraska")
shrink <- 0.3 # shrink that bb to 30\% size
bb <- t (apply (bb, 1, function (i)
                mean (i) + c (-shrink, shrink) * diff (i) / 2))
ms_generate_map (bb, max_tiles = 16L, mapname = "omaha")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rectify-maps.R
\name{ms_rectify_map}
\alias{ms_rectify_map}
\title{Rectify one map to another}
\usage{
ms_rectify_map(
  map_original,
  map_modified,
  nitems = NULL,
  non_linear = 1,
  type = "hulls",
  downsample = 10,
  concavity = 0,
  length_threshold = 10,
  quiet = FALSE
)
}
\arguments{
\item{map_original}{File name of the original map without anything drawn over
it (either a \code{.pdf} or \code{.png}; extension will be ignored).}

\item{map_modified}{File name of the modified version with drawings (either a
\code{.pdf} or \code{.png}; extension will be ignored).}

\item{nitems}{Optional parameter to explicitly specify number of distinct
items to be extracted from map; if possible, specifying this parameter may
improve results.}

\item{non_linear}{Integer value of 0, 1, or 2 representing degree of
non-linearity in modified image - see Note.}

\item{type}{Currently either "points", "polygons", or "hulls", where
"points" simply reduces each distinct object to a single, central point;
"polygons" identifies individual groups and returns the polygon representing
the outer boundary of each; and "hulls" constructs convex or concave polygons
around each group.}

\item{downsample}{Factor by which to downsample \code{type = "polygons"}, noting
that polygons initially include every outer pixel of image, so can generally
be downsampled by at least an order or magnitude (\code{n = 10}). Higher values
may be used for higher-resolution images; lower values will generally only be
necessary for very low lower resolution images.}

\item{concavity}{For \code{type = "hulls"}, a value between 0 and 1, with 0 giving
convex hulls and 1 giving highly concave hulls.}

\item{length_threshold}{For \code{type = "hulls"}, the minimal length of a segment
to be made more convex. Low values will produce highly detailed hulls which
may cause problems; if in doubt, or if odd results appear, increase this
value.}

\item{quiet}{If \code{FALSE}, display progress information on screen}
}
\value{
An \pkg{sf} object representing the drawn additions to map_modified.
}
\description{
Rectify two previously scanned-in pdf or png maps with \code{RNiftyReg}, and
return the modifications in \code{map_modified} as spatial objects in \pkg{sf}
format.
}
\note{
The \code{non-linear} parameter should generally set according to how the
modified maps were digitised. A value of 0 will give fastest results, and
should be used for directly scanned or photocopied images. A value of 1 (the
default) still presumes modified images have been linearly translated, and
will apply affine transformations (rotations, contractions, dilations). This
value should be used when modified images have been photographed (potentially
from an oblique angle). A value of 2 should only be used when modified maps
have somehow been non-linearly distorted, for example through having been
crumpled or screwed up. Rectification with \code{non-linear = 2} will likely take
considerably longer than with lower values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapscanner-package.R
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
% Please edit documentation in R/file-manipulation.R
\name{ms_rotate_map}
\alias{ms_rotate_map}
\title{Rotate maps}
\usage{
ms_rotate_map(map_original, map_modified, rotation = 0, apply_rotation = FALSE)
}
\arguments{
\item{map_original}{File name of the original map without anything drawn over
it (either a \code{.pdf} or \code{.png}; extension will be ignored).}

\item{map_modified}{File name of the modified version with drawings (either a
\code{.pdf} or \code{.png}; extension will be ignored).}

\item{rotation}{Rotation value to be applied, generally +/- 90}

\item{apply_rotation}{If \code{FALSE}, display results of rotation without
actually applying it; otherwise transform the specified \code{map_modified} image
according to the specified rotation.}
}
\value{
No return value. Function either modifies files on disk by rotating
images by the specified amount (if \code{apply_rotation = TRUE}), or displays a
rotated version of \code{map_original} (if \code{apply_rotation = FALSE}).
}
\description{
Display original and modified maps to determine necessary rotation
}
\note{
If a call to \link{ms_rectify_map} detects potential image rotation,
that function will stop and suggest that rotation be applied using this
function in order to determine the required degree of image rotation. Values
for \code{rotation} can be trialled in order to determine the correct value,
following which that value can be entered with \code{apply_rotation = TRUE} in
order to actually apply that rotation to the modified image.
}
