
<!-- README.md is generated from README.Rmd. Please edit that file -->

# terrainr: Landscape Visualization in R and Unity <a href='https://docs.ropensci.org/terrainr/'><img src="man/figures/logo.png" align="right" height="138.5"/></a>

<!-- badges: start -->

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04060/status.svg)](https://doi.org/10.21105/joss.04060)
[![License:
MIT](https://img.shields.io/badge/license-MIT-green)](https://choosealicense.com/licenses/mit/)
[![CRAN
status](https://www.r-pkg.org/badges/version/terrainr)](https://cran.r-project.org/package=terrainr)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![codecov](https://codecov.io/gh/ropensci/terrainr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/terrainr)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build
status](https://github.com/ropensci/terrainr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/terrainr/actions)
[![rOpenSci Review
Status](https://badges.ropensci.org/416_status.svg)](https://github.com/ropensci/software-review/issues/416)

<!-- badges: end -->

## Overview

terrainr makes it easy to retrieve elevation and base map image tiles
for areas of interest within the United States from the [National Map
family of APIs](https://apps.nationalmap.gov/services), and then process
that data into larger, joined images or crop it into tiles that can be
imported into the Unity 3D rendering engine.

There are three main utilities provided by terrainr. First, users are
able to download data from the National Map via the `get_tiles`
function, downloading data tiles for the area represented by an `sf` or
`Raster` object:

``` r
library(terrainr)
library(sf)

# Optional way to display a progress bar while your tiles download:
library(progressr)
handlers("progress")
handlers(global = TRUE)

location_of_interest <- tmaptools::geocode_OSM("Hyampom California")$coords

location_of_interest <- data.frame(
  x = location_of_interest[["x"]],
  y = location_of_interest[["y"]]
)

location_of_interest <- st_as_sf(
  location_of_interest, 
  coords = c("x", "y"), 
  crs = 4326
)

location_of_interest <- set_bbox_side_length(location_of_interest, 8000)

output_tiles <- get_tiles(location_of_interest,
                          services = c("elevation", "ortho"),
                          resolution = 30 # pixel side length in meters
                          )
```

Once downloaded, these images are in standard GeoTIFF or PNG formats and
can be used as expected with other utilities:

``` r
raster::plot(raster::raster(output_tiles[["elevation"]][[1]]))
```

<img src="man/figures/20210728elevation.jpg" width="100%" />

``` r
raster::plotRGB(raster::brick(output_tiles[["ortho"]][[1]]), scale = 1)
```

<img src="man/figures/20210728naip.jpg" width="100%" />

Finally, terrainr helps you visualize this data, both natively in R via
the new `geom_spatial_rgb` geom:

``` r
library(ggplot2)
ggplot() + 
  geom_spatial_rgb(data = output_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  coord_sf(crs = 4326) + 
  theme_void()
```

<img src="man/figures/20210728ggplot.jpg" width="100%" />

As well as with the Unity 3D rendering engine, allowing you to fly or
walk through your downloaded data sets in 3D and VR:

``` r
with_progress( # When not specifying resolution, default is 1m pixels
  output_tiles <- get_tiles(location_of_interest,
                            services = c("elevation", "ortho"))
)

merged_dem <- merge_rasters(output_tiles[["elevation"]], 
                            tempfile(fileext = ".tif"))
merged_ortho <- merge_rasters(output_tiles[["ortho"]], 
                              tempfile(fileext = ".tif"))

make_manifest(output_tiles$elevation,
              output_tiles$ortho)
```

We can then import these tiles to Unity (following the [Import
Vignette](https://docs.ropensci.org/terrainr/articles/unity_instructions.html))
to create:

<img src="man/figures/20210728unity.jpg" width="100%" />

The more time intensive processing steps can all be monitored via the
[progressr](https://github.com/HenrikBengtsson/progressr) package, so
you’ll be more confident that your computer is still churning along and
not just stalled out. For more information, check out [the introductory
vignette](https://docs.ropensci.org/terrainr//articles/overview.html)
and [the guide to importing your data into
Unity\!](https://docs.ropensci.org/terrainr//articles/unity_instructions.html)

## Citing terrainr

The United States Geological Survey provides guidelines for citing USGS
data products (as downloaded from `get_tiles`) at
<https://www.usgs.gov/faqs/how-should-i-cite-datasets-and-services-national-map>
.

To cite terrainr in publications please use:

> Mahoney et al., (2022). terrainr: An R package for creating immersive
> virtual environments. Journal of Open Source Software, 7(69), 4060,
> <https://doi.org/10.21105/joss.04060>

A BibTeX entry for LaTeX users is:

``` bibtex
  @Article{,
    year = {2022},
    publisher = {The Open Journal},
    volume = {7},
    number = {69},
    pages = {4060},
    author = {Michael J. Mahoney and Colin M. Beier and Aidan C. Ackerman},
    title = {{terrainr}: An R package for creating immersive virtual environments},
    journal = {Journal of Open Source Software},
    doi = {10.21105/joss.04060},
    url = {https://doi.org/10.21105/joss.04060},
  }
```

## Available Datasets

The following datasets can currently be downloaded using `get_tiles` or
`hit_national_map_api`:

  - [3DEPElevation](https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer):
    The USGS 3D Elevation Program (3DEP) Bare Earth DEM.
  - [USGSNAIPPlus](https://services.nationalmap.gov/arcgis/rest/services/USGSNAIPPlus/MapServer):
    National Agriculture Imagery Program (NAIP) and high resolution
    orthoimagery (HRO).
  - [nhd](https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer):
    A comprehensive set of digital spatial data that encodes information
    about naturally occurring and constructed bodies of surface water
    (lakes, ponds, and reservoirs), paths through which water flows
    (canals, ditches, streams, and rivers), and related entities such as
    point features (springs, wells, stream gauges, and dams).
  - [govunits](https://carto.nationalmap.gov/arcgis/rest/services/govunits/MapServer):
    Major civil areas for the Nation, including States or Territories,
    counties (or equivalents), Federal and Native American areas,
    congressional districts, minor civil divisions, incorporated places
    (such as cities and towns), and unincorporated places.
  - [contours](https://carto.nationalmap.gov/arcgis/rest/services/contours/MapServer):
    The USGS Elevation Contours service.
  - [geonames](https://carto.nationalmap.gov/arcgis/rest/services/geonames/MapServer):
    Information about physical and cultural geographic features,
    geographic areas, and locational entities that are generally
    recognizable and locatable by name.
  - [NHDPlus\_HR](https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer):
    A comprehensive set of digital spatial data comprising a nationally
    seamless network of stream reaches, elevation-based catchment areas,
    flow surfaces, and value-added attributes.
  - [structures](https://carto.nationalmap.gov/arcgis/rest/services/structures/MapServer):
    The name, function, location, and other core information and
    characteristics of selected manmade facilities.
  - [transportation](https://carto.nationalmap.gov/arcgis/rest/services/transportation/MapServer):
    Roads, railroads, trails, airports, and other features associated
    with the transport of people or commerce.
  - [wbd](https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer):
    Hydrologic Unit (HU) polygon boundaries for the United States,
    Puerto Rico, and the U.S. Virgin Islands.

(All descriptions above taken from the [National Map API
descriptions](https://apps.nationalmap.gov/services).)

## Installation

You can install terrainr from CRAN via:

``` r
install.packages("terrainr")
```

Or, if you want the newest patches and features, you can install the
development version of terrainr from
[GitHub](https://github.com/ropensci/terrainr) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/terrainr")
```

Be aware that the development version is not stable, and features that
haven’t been published on CRAN may change at any time\!

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# terrainr 0.5.1
* New features:
    * A new endpoint, `ecosystems`, has been added to `get_tiles` and 
      `hit_national_map_api`. 
* Improvements and bug fixes:
    * `merge_rasters` gains an argument, `overwrite`, which allows you to 
      specify whether or not to overwrite `output_raster` if it exists. Previous
      versions expected you to pass "-overwrite" to `options`. If a file exists
      at `output_raster` and `overwrite` is FALSE, `merge_rasters` will throw an
      error.
* Dependency changes:
    * `sf` now has a minimum dependency of 1.0-5, to take advantage of an 
      upstream bug fix (relating to `merge_rasters` overwrite)

# terrainr 0.5.0.
* New features:
    * A new function, `make_manifest`, now helps automate the import of terrain
      and imagery to Unity. It fully replaces `raster_to_raw_tiles` (see 
      Deprecations below). Documentation updates are forthcoming.
* Deprecations:
    * `raster_to_raw_tiles` is now deprecated and will be removed in a future
      release (no earlier than 2022). Use `make_manifest` instead.
    * The method `get_tiles.list` is now deprecated and will be removed in a 
      future release (unexported in Fall 2021, removed no earlier than 2022).
      Convert your list to an `sf` object instead.
    * The `bbox` argument to `hit_national_map_api` is now documented as 
      "An object from [sf::st_bbox]." This is a change from the earlier options 
      of a length 2 list or terrainr_bounding_box object. Those methods are 
      currently still supported, but undocumented; they will be removed in a 
      future release (no earlier than 2022).
* Improvements and bug fixes:
    * `get_tiles` no longer mangles data with projected coordinates (via a 
      fix to the internal function `split_bbox`). If for some reason you want 
      the old behavior back, set the new argument `projected` to `FALSE` while
      providing projected data.
    * The documentation for `add_bbox_buffer` and `set_bbox_side_length` now 
      specifies that they should only be used with geographic coordinate 
      systems. If you use these functions with projected data, they will warn;
      this may be upgraded to an error in future versions.
    * The README images are now of beautiful Hyampom, California, a somewhat 
      more appealing vista than the original Mt Marcy scene.
    * The "Import to Unity" vignette has been rewritten to use make_manifest, as
      has the overview vignette and other documentation.
    * Typos in the message `merge_rasters` gives when using the fallback method
      have been fixed.
    * `merge_rasters` gains an argument `force_fallback` which, if TRUE, will 
      use the older, slower method for merging tiles. This is not recommended, 
      but is useful for testing.
* Internal changes:
    * The slow removal of all `terrainr_*` custom classes marches on! These 
      classes should no longer be present in any user-facing, non-deprecated 
      code; the only functions still relying on custom classes are internal 
      utilities and the `split_bbox` function responsible for tiling `get_tiles`
      requests.
    * `split_bbox` should now run faster, particularly for large tile sets, as
      some nested loops have been vectorized.
    * Improvements to test coverage and CI.

# terrainr 0.4.1
* Improvements and bug fixes:
    * `get_tiles` now displays a bulleted list of endpoints (again?), rather 
      than a jumble of raw markdown
    * `add_bbox_buffer` properly sets the CRS of the output when attempting to 
      buffer geodesic coordinates.
    * Typo fixes to an error message in `combine_overlays`
* Internal changes:
    * Added `importFrom` tag to `terrainr-package.R` to silence R CMD CHECK NOTE.

# terrainr 0.4.0
* Breaking changes:
    * Three changes in how `vector_to_overlay` deals with missing CRS in 
      `vector_data`:
        * A new argument, `error_crs`, behaves just like `error_crs` in 
          `add_bbox`: if `NULL`, the function will give a warning when assuming
          CRS; if `FALSE`, the function will assume a CRS silently, and if 
          `TRUE`, the function will error if `vector_data` is missing a CRS.
        * `target_crs` has been removed. `vector_data` will be given the CRS of 
          `reference_raster` if it doesn't have its own CRS, and will always be
          projected to the CRS of `reference_raster`.
        * `error_crs` has been added to mirror `add_bbox_buffer`: if `NULL` and
          your input data has no CRS, `vector_to_overlay` will warn about 
          assuming the raster CRS. Set to `TRUE` to error or `FALSE` to ignore
          the warning.
    * NAIP imagery is now downloaded with `transparent = "false"` to
      minimize the number of times the backup method to `merge_rasters` (see 
      below) is called. To restore the old behavior, set `transparent = "true"` 
      in either `get_tiles` or `hit_national_map_api`.
    * `get_tiles` will now infer `bboxSR` and `imageSR` from provided `sf` or
      `Raster` objects if not otherwise specified. To restore the old behavior,
      set `bboxSR` and `imageSR` to `4326` in `get_tiles` (or set your data's 
      CRS to 4326 before calling `get_tiles`).
* Improvements and bug fixes:
    * `merge_rasters` can once again handle merging mixed-band rasters (such as
      NAIP images with and without alpha bands). At the moment this is using the
      older, slower implementation and will raise a warning. (#30, #32).
* Internal changes:
    * Removed code to check for `ggplot2` from `vector_to_overlay` now that
      `ggplot2` is required
    * `calc_haversine_distance` (not exported) now assumes it's been provided
      with degrees. `coord_units` has been removed as an argument.
    * `get_tiles.terrainr_bounding_box` has been removed; it should no longer be
      possible for users to have `terrainr_bounding_box` objects unless they 
      were using non-exported functionality.

# terrainr 0.3.1
* First CRAN submission!
* This is the smallest of patch releases, with almost no user-facing changes.
* Internal changes:
    * Added rOpenSci reviewers to DESCRIPTION.
    * Changed USGS API link to new website.
    * Added rOpenSci badge to README.
    * Changed most PNG images to _slightly_ smaller JPGs.
    * Edited URLs for new rOpenSci website.
    * Moved lifecycle badge href to new site.
    * Some small spelling issues have been fixed.
    * Added \value tags to non-exported point_to_distance and 
      terrainr_bounding_box functions
    * Added single quotes around Unity in the DESCRIPTION

# terrainr 0.3.0
* Breaking changes:
    * `terrainr_*` classes have been effectively removed and are no longer 
      exported. Functions which previously expected these objects now generally 
      accept `sf` and `Raster` class objects instead. Functions which previously
      returned these objects now generally return `sf` objects instead (#24).
    * The list returned by `get_tiles` now uses the service names provided by
      the user, not the endpoint names. This means that 
      `get_tiles(..., services = "elevation")` will now use the name `elevation`
      instead of `3DEPElevation`, and remain standard across versions (#12).
    * `get_bbox` and `get_coordinate_bbox` have been removed. Functions that 
      used to expect `terrainr_bounding_box` objects now accept objects of class
      `sf` or `raster` (#24).
    * `add_bbox_buffer` loses the `divisible` argument. For precise control over
      side length, use `set_bbox_side_length` (which should be more accurate, if
      slightly more conservative, than the `divisible` system ever was) (#17).
    * `convert_distance` has been removed (internally replaced by the 
      `units` package) (#7).
    * `merge_rasters` loses the `input_images` and `output_image` function, as
      most downloaded files are now already georeferenced. To recreate this 
      functionality, georeference image tiles directly via 
      `output <- georeference_overlay(img_tiles, ref_tiles, tempfile(fileext = ".tif"))`
      and then provide `output` to `merge_rasters`.
    * A handful of utility functions are no longer exported:
        * `calc_haversine_distance`
        * `point_from_distance`
        * `rad_to_deg`
        * `deg_to_rad`
* New features:
    * Two new functions, `geom_spatial_rgb` and `stat_spatial_rgb`, allow you to 
      use RGB map tiles as backgrounds for further plotting.
    * `calc_haversine_distance` gains an argument `coord_units` allowing it to 
      handle coordinates in radians as well as degrees.
* Improvements and bug fixes:
    * `georeference_overlay` provides `tempfile(fileext = ".tif")` as a default
      output location if no `output_file` is provided.
    * `get_tiles` now tells you what tiles it's retrieving, not retriving.
* Internal changes:
    * `calc_haversine_distance` has been internally simplified somewhat to 
      reduce code duplication.
    * All `services` arguments to `hit_national_map_api` and `get_tiles` can 
      now handle both base64 and binary returns, removing the need to manually
      categorize endpoints (54ad9fb).
        * `hit_national_map_api` auto-detects whether API endpoints are 
          returning base64 or binary and handles them appropriately
        * `get_tiles` now auto-detects whether `hit_national_map_api` is 
          returning base64 or binary and writes to file appropriately.
    * `hit_national_map_api` is now more likely to fail with a human-friendly
      error message if API endpoints return a non-200 status (54ad9fb).
    * `hit_national_map_api` (and by extension `get_tiles`) now register a user
      agent.
* Changes in dependencies:
    * `gdalUtilities` has been removed, with functionality replaced by `sf`.
    * `rlang` has been removed, with functionality removed.
    * `units` has been added.
    * `ggplot2` has been moved to Imports (was previously in Suggests) due to 
      the new `geom_spatial_rgb` and `stat_spatial_rgb` functions.

# terrainr 0.2.1
* Improvements and bug fixes:
    * The `transportation` endpoint has moved servers, and is now handled by the
      same function that handles DEMs and orthoimages
* Internal changes:
    * The main branch of `terrainr` is now `main`
    * Tests run on a schedule on Monday/Wednesday/Friday mornings, to alert to 
      endpoint changes
    * Restyled code

# terrainr 0.2.0
* Breaking changes:
    * `merge_rasters` loses the argument `merge_raster`. For the "georeference 
      a single image" use case, see the new `georeference_overlay` function.
    * `get_tiles` gains an argument `resolution` (details below) between 
      `side_length` and `services`. No functionality should be changed, but code 
      with unnamed arguments to `services`, `verbose`, or `georeference` may be 
      impacted. 
* New features:
    * A new family of functions for dealing with overlay creation:
        * `vector_to_overlay` lets users quickly produce image overlays from 
          vector data.
        * `georeference_overlay` replaces the use of merge_raster for creating 
          single-file georeferenced overlay files.
        * `combine_overlays` lets users, well, combine overlays into a single 
          image
    * `get_tiles` gains an argument, `resolution`, specifying the number of 
      meters each pixel should represent (so higher images result in smaller 
      downloads). 
    * `get_bbox` provides an S3 generic to create `terrainr_bounding_box` 
      objects. In this version, that means users can use `get_bbox` to get 
      bounding boxes from `sf` and `RasterLayer` objects, and it means adding 
      methods will be easier going forward. The generic `get_bbox` method 
      is equivalent to `get_coord_bbox`
    * `raster_to_raw_tiles` handles rectangles appropriately
* Improvements and bug fixes:
    * `get_tiles`, `raster_to_raw_tiles`, and `merge_rasters` are now much more 
      conscientious about deleting tempfiles when they're done with them.
    * `merge_rasters` no longer fails when handed a mix of 3- and 4-band raster
      files. The current implementation will cast all 4 band rasters to 3 band
      images and then return a 3 band raster image.
    * The `output_image` argument to `merge_rasters` now has a default value of 
      `tempfile(fileext = ".tif")` to be a little more friendly to users.
    * Arguments `lat` and `lng` to `get_bbox` (and `get_coord_bbox`) no longer 
      need to be quoted -- either the tidyverse-feeling NSE approach or the 
      more standard quoted argument approach will work.
* Internal changes:
    * All terrainr-provided functions now explicitly use the terrainr:: 
      namespace.
* Changes in dependencies:
    * `sf` has been added as an explicit import due to `vector_to_overlay`. `sf`
      is required by `gdalUtilities`, also imported by this package, so this 
      change should have no impact on users.
    * `rlang` is added as a dependency to allow `lat` and `lng` be unquoted in 
      `get_bbox`.
    * `ggplot2` has been added to `Suggests` due to `vector_to_overlay`. 
    * `jpeg` and `tiff` have been added to `Suggests` due to 
      `georeference_overlay`. I'd expect more image libraries to join this list
      over time. 

# terrainr 0.1.0
* New features:
    * set_bbox_side_length wraps add_bbox_buffer to set each side of the 
      bounding box to an equal length (within ~1% accuracy)
* First version released on GitHub 

# terrainr 0.0.0.9001

* First development version
* Supports retrieval from 3DEP and NAIP data sources
* Supports export to Unity-friendly format
* Functions in this version:
    * Utility functions: 
        * add_bbox_buffer
        * calc_haversine_distance
        * convert_distance
        * deg_to_rad
        * get_bbox_centroid
        * get_coord_bbox
        * point_from_distance
        * rad_to_deg
    * Data retrieval functions:
      * get_tiles 
      * hit_national_map_api
    * Data processing functions:
      * merge_rasters
      * raster_to_raw_tiles
    * Classes and class utility functions:
      * terrainr_bounding_box (class)
      * terrainr_coordinate_pair (class)
      * terrainr_bounding_box (creation utility)
      * terrainr_coordinate_pair (creation utility)
      * export_bounding_box
      * export_coord_pair
## Resubmission

This is a resubmission. In this version, I fixed the incorrect DOI in the CITATION file and edited the URL format in the description. Thank you for pointing that out!

## Test environments
* local R installation, R 4.1.2
* ubuntu 20.04 (on GitHub Actions), R devel, 4.1.2, 4.0.5
* MacOS X 10.15.7 (on GitHub Actions), R 4.1.2
* Windows Server 2019 (on GitHub Actions), R 4.1.2
* win-builder (devel, 4.1.2, 4.0.5)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies
There are no downstream dependencies on CRAN for terrainr at the time of 
submission.
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

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the GitHub build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the terrainr project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# terrainr: Landscape Visualization in R and Unity <a href='https://docs.ropensci.org/terrainr/'><img src="man/figures/logo.png" align="right" height="138.5"/></a>

<!-- badges: start -->
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04060/status.svg)](https://doi.org/10.21105/joss.04060)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](https://choosealicense.com/licenses/mit/) [![CRAN status](https://www.r-pkg.org/badges/version/terrainr)](https://cran.r-project.org/package=terrainr) [![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing) [![codecov](https://codecov.io/gh/ropensci/terrainr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/terrainr) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![R build status](https://github.com/ropensci/terrainr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/terrainr/actions) [![rOpenSci Review Status](https://badges.ropensci.org/416_status.svg)](https://github.com/ropensci/software-review/issues/416)

<!-- badges: end -->

## Overview

terrainr makes it easy to retrieve elevation and base map image tiles for areas 
of interest within the United States from the 
[National Map family of APIs](https://apps.nationalmap.gov/services), and 
then process that data into larger, joined images or crop it into tiles that can 
be imported into the Unity 3D rendering engine.

There are three main utilities provided by terrainr. First, users are able to 
download data from the National Map via the `get_tiles` function, downloading 
data tiles for the area represented by an `sf` or `Raster` object:

```{r eval=FALSE}
library(terrainr)
library(sf)

# Optional way to display a progress bar while your tiles download:
library(progressr)
handlers("progress")
handlers(global = TRUE)

location_of_interest <- tmaptools::geocode_OSM("Hyampom California")$coords

location_of_interest <- data.frame(
  x = location_of_interest[["x"]],
  y = location_of_interest[["y"]]
)

location_of_interest <- st_as_sf(
  location_of_interest, 
  coords = c("x", "y"), 
  crs = 4326
)

location_of_interest <- set_bbox_side_length(location_of_interest, 8000)

output_tiles <- get_tiles(location_of_interest,
                          services = c("elevation", "ortho"),
                          resolution = 30 # pixel side length in meters
                          )
```

Once downloaded, these images are in standard GeoTIFF or PNG formats and can be 
used as expected with other utilities:

```{r eval=FALSE}
raster::plot(raster::raster(output_tiles[["elevation"]][[1]]))
```

```{r, echo = FALSE}
knitr::include_graphics("man/figures/20210728elevation.jpg")
```

```{r eval=FALSE}
raster::plotRGB(raster::brick(output_tiles[["ortho"]][[1]]), scale = 1)
```

```{r echo=FALSE}
knitr::include_graphics("man/figures/20210728naip.jpg")
```

Finally, terrainr helps you visualize this data, both natively in R via the new
`geom_spatial_rgb` geom:

```{r eval = FALSE}
library(ggplot2)
ggplot() + 
  geom_spatial_rgb(data = output_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  coord_sf(crs = 4326) + 
  theme_void()
```

```{r echo = FALSE}
knitr::include_graphics("man/figures/20210728ggplot.jpg")
```

As well as with the Unity 3D rendering engine, allowing you to fly or walk 
through your downloaded data sets in 3D and VR:

```{r, results = FALSE, eval=FALSE}
with_progress( # When not specifying resolution, default is 1m pixels
  output_tiles <- get_tiles(location_of_interest,
                            services = c("elevation", "ortho"))
)

merged_dem <- merge_rasters(output_tiles[["elevation"]], 
                            tempfile(fileext = ".tif"))
merged_ortho <- merge_rasters(output_tiles[["ortho"]], 
                              tempfile(fileext = ".tif"))

make_manifest(output_tiles$elevation,
              output_tiles$ortho)
```

We can then import these tiles to Unity (following the 
[Import Vignette](https://docs.ropensci.org/terrainr/articles/unity_instructions.html))
to create:

```{r echo=FALSE}
knitr::include_graphics("man/figures/20210728unity.jpg")
```

The more time intensive processing steps can all be monitored via the [progressr](https://github.com/HenrikBengtsson/progressr) package, so you'll be 
more confident that your computer is still churning along and not just stalled 
out. For more information, check out [the introductory vignette](https://docs.ropensci.org/terrainr//articles/overview.html) and 
[the guide to importing your data into Unity!](https://docs.ropensci.org/terrainr//articles/unity_instructions.html)

## Citing terrainr

The United States Geological Survey provides guidelines for citing USGS data products (as
downloaded from `get_tiles`) at 
https://www.usgs.gov/faqs/how-should-i-cite-datasets-and-services-national-map .

To cite terrainr in publications please use:

> Mahoney et al., (2022). terrainr: An R package for creating immersive virtual environments. Journal of Open Source Software, 7(69), 4060, https://doi.org/10.21105/joss.04060

A BibTeX entry for LaTeX users is:

```{bibtex}
  @Article{,
    year = {2022},
    publisher = {The Open Journal},
    volume = {7},
    number = {69},
    pages = {4060},
    author = {Michael J. Mahoney and Colin M. Beier and Aidan C. Ackerman},
    title = {{terrainr}: An R package for creating immersive virtual environments},
    journal = {Journal of Open Source Software},
    doi = {10.21105/joss.04060},
    url = {https://doi.org/10.21105/joss.04060},
  }
```


## Available Datasets

The following datasets can currently be downloaded using `get_tiles` or `hit_national_map_api`:

-   [3DEPElevation](https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer): The USGS 3D Elevation Program (3DEP) Bare Earth DEM.
-   [USGSNAIPPlus](https://services.nationalmap.gov/arcgis/rest/services/USGSNAIPPlus/MapServer): National Agriculture Imagery Program (NAIP) and high resolution orthoimagery (HRO).
-   [nhd](https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer): A comprehensive set of digital spatial data that encodes information about naturally occurring and constructed bodies of surface water (lakes, ponds, and reservoirs), paths through which water flows (canals, ditches, streams, and rivers), and related entities such as point features (springs, wells, stream gauges, and dams).
-   [govunits](https://carto.nationalmap.gov/arcgis/rest/services/govunits/MapServer): Major civil areas for the Nation, including States or Territories, counties (or equivalents), Federal and Native American areas, congressional districts, minor civil divisions, incorporated places (such as cities and towns), and unincorporated places.
-   [contours](https://carto.nationalmap.gov/arcgis/rest/services/contours/MapServer): The USGS Elevation Contours service.
-   [geonames](https://carto.nationalmap.gov/arcgis/rest/services/geonames/MapServer): Information about physical and cultural geographic features, geographic areas, and locational entities that are generally recognizable and locatable by name.
-   [NHDPlus\_HR](https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer): A comprehensive set of digital spatial data comprising a nationally seamless network of stream reaches, elevation-based catchment areas, flow surfaces, and value-added attributes.
-   [structures](https://carto.nationalmap.gov/arcgis/rest/services/structures/MapServer): The name, function, location, and other core information and characteristics of selected manmade facilities.
-   [transportation](https://carto.nationalmap.gov/arcgis/rest/services/transportation/MapServer): Roads, railroads, trails, airports, and other features associated with the transport of people or commerce.
-   [wbd](https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer): Hydrologic Unit (HU) polygon boundaries for the United States, Puerto Rico, and the U.S. Virgin Islands.

(All descriptions above taken from the [National Map API descriptions](https://apps.nationalmap.gov/services).)

## Installation

You can install terrainr from CRAN via:

``` r
install.packages("terrainr")
```

Or, if you want the newest patches and features, you can install the development 
version of terrainr from [GitHub](https://github.com/ropensci/terrainr) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/terrainr")
```

Be aware that the development version is not stable, and features that haven't
been published on CRAN may change at any time!

## Code of Conduct

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "A gentle introduction to terrainr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{A gentle introduction to terrainr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The goal of the terrainr package is to make it easier to visualize landscapes,
both by providing functions to download elevation data and base maps from the 
USGS National Map and by adding utilities to manipulate base maps for use in 
visualizations using [ggplot2](https://ggplot2.tidyverse.org/) or the 
freely available [Unity 3D rendering engine](https://unity.com/).
This vignette will walk through the core functions available in terrainr and 
how they interact.

Let's load the terrainr package to get started:

```{r setup}
library(terrainr)
```

We're going to work with data for Mount Elbert today, the highest point in the 
Rocky Mountain range. I'm just choosing this location for the dramatic scenery;
the National Map can be used to retrieve data for the entire United States and 
much of Canada. Let's simulate some data for the area right around Mt. Elbert,
such as the point data we might get from some field collection:

```{r}
mt_elbert_points <- data.frame(
  lat = runif(100, min = 39.11144, max = 39.12416),
  lng = runif(100, min = -106.4534, max = -106.437)
)
```

terrainr is built to play nicely with functions from the 
[sf](https://r-spatial.github.io/sf/) and [raster](https://rspatial.org/raster/) 
packages. In order to get our simulated points into the right format, we need to 
use the `st_as_sf` function from the `sf` package:

```{r}
mt_elbert_points <- sf::st_as_sf(mt_elbert_points,
                                 coords = c("lng", "lat"))
mt_elbert_points <- sf::st_set_crs(mt_elbert_points, 4326)
```

Now that we've got our data in the right format, it's time to retrieve our data. 
terrainr currently supports downloading DEMs from the USGS 3D Elevation Program 
as well as orthoimages from the National Agricultural Imagery Program, in 
addition to other base map images from the National Map. These programs each 
have slightly different APIs and different restrictions on file types and the 
size of image you can download at once. 
Rather than make you think about this, terrainr handles all the edges of making 
API requests for you, including splitting your request into tiles and 
formatting the query. 

For this vignette, we'll retrieve both elevation and orthoimagery using the 
`get_tiles` function. We can either use the generic "elevation" and "ortho"
shorthands to get our data, or we can specify "3DEPElevation" and "USGSNAIPPlus"
to make sure we're using the same specific service -- the short codes aren't
guaranteed to download data from the same service between releases!

One last note -- all the longer-running terrainr functions can print out 
progress bars, if the user requests them via the 
[progressr](https://CRAN.R-project.org/package=progressr) 
package. We'll demonstrate that syntax here:

```{r eval = FALSE}
library(progressr)
handlers("progress")
with_progress(
  output_files <- get_tiles(mt_elbert_points,
                            output_prefix = tempfile(),
                            services = c("elevation", "ortho"))
  )
```

And just like that, we have our data tiles! To make multi-step processing 
easier, terrainr functions which deal with these tiles typically return lists of 
the file paths they saved your data to.

```{r eval = FALSE}
output_files
```

```{r echo = FALSE}
output_files <- list(
  elevation = "/tmp/RtmphTFQvZ/file65e5d859e628_3DEPElevation_1_1.tif",
  ortho = "/tmp/RtmphTFQvZ/file65e5d859e628_USGSNAIPPlus_1_1.tif"
)
output_files
```

If we were requesting more data than we can download at once, each element of 
the list would be a character vector containing the file paths for all of our 
downloaded tiles. Since we're sticking with a relatively small area for this 
example, we only have one tile for each service. 

As a quick aside, note that 
you can control where these files save to via the `output_prefix` argument 
(which appends the suffix `servicename_xindex_yindex.tif` to each tile it 
downloads) -- you don't need to save them to a temporary directory (and 
redownload every time you launch R) as we're doing here!

If all you want is to access these endpoints to download data, this is probably
the only terrainr function you'll need -- the files produced by this function
can be processed just like any other spatial data:

```{r eval=FALSE}
raster::plot(raster::raster(output_files[[1]]))
```

```{r echo = FALSE}
knitr::include_graphics("example_dem.jpg")
```

```{r eval=FALSE}
raster::plotRGB(raster::brick(output_files[[2]]), scale = 1)
```

```{r echo = FALSE}
knitr::include_graphics("example_ortho.jpg")
```

In addition to the regular methods for plotting rasters in R, terrainr makes it 
a bit easier to use ggplot2 for plotting the data returned by `get_tiles`. 
Plotting single-band rasters, like our elevation file, is already well-supported
in base ggplot2:

```{r eval=FALSE}
library(ggplot2)

elevation_raster <- raster::raster(output_files[[1]])
elevation_df <- as.data.frame(elevation_raster, xy = TRUE)
elevation_df <- setNames(elevation_df, c("x", "y", "elevation"))

ggplot() + 
  geom_raster(data = elevation_df, aes(x = x, y = y, fill = elevation)) + 
  scale_fill_distiller(palette = "BrBG") + 
  coord_sf(crs = 4326)
```

```{r echo = FALSE}
knitr::include_graphics("elevation_ggplot.jpg")
```

terrainr adds the ability to plot using multi-band RGB rasters, like the tiles
downloaded for non-elevation endpoints, using the new `geom_spatial_rgb`
function (or its partner, `stat_spatial_rgb`):

```{r eval = FALSE}
ortho_raster <- raster::stack(output_files[[2]])
ortho_df <- as.data.frame(ortho_raster, xy = TRUE)
ortho_df <- setNames(ortho_df, c("x", "y", "red", "green", "blue", "alpha"))

ggplot() + 
  geom_spatial_rgb(data = ortho_df,
                   # Required aesthetics r/g/b specify color bands:
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  coord_sf(crs = 4326)
```

```{r}
knitr::include_graphics("ortho_ggplot.jpg")
```


Note that `geom_spatial_rgb` is a little different from other ggplot2 geoms in 
that it can also accept RasterStack objects directly:

```{r eval = FALSE}
ggplot() + 
  geom_spatial_rgb(data = ortho_raster,
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  coord_sf(crs = 4326)
```

```{r echo = FALSE}
knitr::include_graphics("ortho_ggplot.jpg")
```

Or length 1 character vectors with a path to a file that can be read by 
`raster::stack`:

```{r eval = FALSE}
ggplot() + 
  geom_spatial_rgb(data = output_files[[2]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  coord_sf(crs = 4326)
```

```{r echo = FALSE}
knitr::include_graphics("ortho_ggplot.jpg")
```

You can then use these multi-band rasters as base maps for further plotting as 
desired.

```{r eval = FALSE}
ggplot() + 
  geom_spatial_rgb(data = output_files[[2]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  geom_sf(data = mt_elbert_points)
```

```{r echo = FALSE}
knitr::include_graphics("with_points.jpg")
```

In case you find this visualization falls a little bit flat, terrainr also 
provides the ability to bring your landscapes into Unity to visualize in 3D. 
Our first step in this process is going to be replicating our image of 
base-map-plus-field-sites (above) in a format that we can import into Unity 
directly. First, we'll need to use the `vector_to_overlay` function to create
an image overlay from our point data:

```{r, eval = FALSE}
mt_elbert_overlay <- vector_to_overlay(mt_elbert_points,
                                       output_files[[2]],
                                       size = 15,
                                       color = "red")
knitr::include_graphics(mt_elbert_overlay)
```

```{r, echo = FALSE}
knitr::include_graphics("mt_elbert_overlay.jpg")
```

Note that `vector_to_overlay` can be used with any sf object, not just point 
data.

These overlays may be stacked on top of one another or downloaded imagery using
the `combine_overlays` function:

```{r, eval=FALSE}
ortho_with_points <- combine_overlays(
  # Overlays are stacked in order, with the first file specified on the bottom
  output_files[[2]],
  mt_elbert_overlay,
  output_file = tempfile(fileext = ".png")
  )
knitr::include_graphics(ortho_with_points)
```

```{r, echo = FALSE}
knitr::include_graphics("combined_overlay.jpg")
```

Unfortunately, this image processing strips the georeferencing on the image. 
We can restore the original georeferencing via the `georeference_overlay`
function:

```{r eval = FALSE}
georef_overlay <- georeference_overlay(
  ortho_with_points,
  output_files[[2]]
)
```

We've been working so far with a single tile, but Unity is able to handle much,
much larger rasters than we would normally work with in R. In order to create
overlays for these larger rasters, it's usually best to create an overlay for
smaller image tiles which can then be joined back together with `merge_rasters`:

```{r, eval = FALSE}
tile_overlays <- lapply(output_files[[2]],
                        function(x) vector_to_overlay(mt_elbert_points, 
                                                      x, 
                                                      size = 15, 
                                                      color = "red", 
                                                      na.rm = TRUE))

combined_tiles <- mapply(function(x, y) {
  combine_overlays(x, y, output_file = tempfile(fileext = ".png"))
  },                            
  output_files[[2]],
  tile_overlays)

georef_tiles <- mapply(georeference_overlay, combined_tiles, output_files[[2]])

merged_tiles <- merge_rasters(georef_tiles)
```

Of course, since we're only working with a single tile, `georef_tiles` is 
identical to `merged_tiles`. But when working with larger areas, `merged_tiles`
is particularly useful for joining the separate tiles downloaded by `get_tiles`
into a single raster file.

In particular, having a single joined raster is necessary for the function
`make_manifest`, which is designed to turn these larger rasters into tiles 
in a format that can be imported into the Unity 3D rendering engine. You can 
find more information about that process in 
[the Unity vignette](unity_instructions.html).

```{r eval = FALSE}
elevation_tile <- output_files[[1]]
make_manifest(elevation_tile, georef_tiles)
```

After that function runs, it's a matter of minutes to create beautiful -- 
and fully physically-simulated -- landscape visualizations of your area of 
interest in Unity:

```{r echo = FALSE}
knitr::include_graphics("ebert_unity.jpg")
```

For more instructions on how to create these 3D simulations in Unity, check out
the [Unity vignette](https://docs.ropensci.org/terrainr/articles/unity_instructions.html).
---
title: "Importing terrainr tiles into Unity"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Importing terrainr tiles into Unity}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

terrainr advertises itself as a package for landscape visualization in R and 
[Unity](https://unity.com/). This vignette focuses on the Unity half of the 
equation -- specifically, on how to mostly-automatically import tiles into 
Unity. If you're interested in the R half of the package, refer to 
(the overview vignette)[overview.html]. Note that this vignette will assume you 
already have Unity installed on your computer.

In order to import terrain tiles into Unity, we will first need to have some 
data worth importing! For the purposes of this vignette, we'll be using data 
from the USGS National Map downloaded using `get_tiles`, but note that you can
use _any_ raster data for this process.

First things first, I'm going to use the `geocode_OSM` function from `tmaptools` 
to get the latitude and longitude of Zion National Park, out in Utah. We'll use
this area for our visualization today:

```{r, eval = FALSE}
zion <- tmaptools::geocode_OSM("Zion National Park")$coords
```

Our `zion` object now contains the x and y coordinates for a spot near the 
middle of Zion National Park. Let's go ahead and turn that into an `sf` object,
then use `set_bbox_side_length` to add a buffer around the point coordinates -- 
we'll download data for the entire 8 kilometer square around the central point:

```{r, eval = FALSE}
library(terrainr)
library(sf)
library(magrittr)

zion <- data.frame(x = zion[["x"]], y = zion[["y"]]) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  set_bbox_side_length(8000)
```

And now we can go ahead and download data for our area of interest, then merge 
the downloaded tiles into individual rasters. For more on this process or what
these functions do, check out (the overview vignette)[overview.html]. 

Fair warning -- downloading this amount of data can take a bit of time! You can
optionally add a progress bar to the `get_tiles` download by calling 
`library(progressr)` and then `handlers(global = TRUE)` before running this 
code.

```{r, eval = FALSE}
merged_tiles <- zion %>%
  get_tiles(services = c("elevation", "ortho")) %>% 
  lapply(merge_rasters)
```

We've now got our data downloaded! Our next step is to turn it into a data 
format we can import into Unity.

As of terrainr 0.5.0, the way to do this is via the function `make_manifest`. 
The first argument to this function is the elevation raster you want to use as
a heightmap -- in our case, `merged_tiles$elevation`. The second argument 
optionally takes the image overlay you want to put on top of that heightmap. 
In our case, that means everything we need to provide is in the `merged_tiles`
list:

```{r, eval = FALSE}
make_manifest(merged_tiles$elevation, 
              merged_tiles$ortho)
```

After a moment, this function will spit out a number of files: our heightmap and
overlay tiles, all prefixed with `import_`, a C# file named `import_terrain.cs`,
and a final file named `terrainr.manifest` (note that all of these names can be 
changed via arguments to `make_manifest`, but for simplicity's sake I'm using 
the default names now). 

```{r, echo = FALSE}
knitr::include_graphics("generated_files.png")
```

Now go ahead and open Unity. From the main "Hub" menu, click "New" to create a
new project. Set the project name to whatever you want, then click "Create".

```{r, echo = FALSE}
knitr::include_graphics("new_unity.jpg")
```

Go ahead and move all the files from `make_manifest` into the root directory of
your new Unity project. Then move the `import_terrain.cs` file into the `Assets`
directory inside that folder (but leave everything else in the root directory!).

Go back to Unity now. A second after you click into the window, you should 
notice a "terrainr" menu appear in the top bar. Click that menu, then the only
option in the drop-down. A menu should appear; click "Import" to import your 
tiles into Unity.

```{r, echo = FALSE}
knitr::include_graphics("manifest_import.png")
```

The importer menu will disappear, then Unity will take a minute or two to import
all your tiles. Depending on your data, you may see something that looks like
this:

```{r, echo = FALSE}
knitr::include_graphics("ominous.png")
```

That's perfectly fine! Right click on the "Scene" window in the middle, and then
press and hold "S" on your keyboard to move the camera back. After a 
second, you should see your terrain surface!

```{r, echo = FALSE}
knitr::include_graphics("terrain_surface.jpg")
```

You can now move around your surface by right clicking on the image and moving
around with the W-A-S-D keys on your keyboard. Note that your movement speed 
starts off very slow, and then accelerates over time.

And ta-da, you have a surface in Unity! You can go ahead and customize the scene
further (I'll usually then click on "Directional Light" and change "Render Mode" 
to "Not Important" and "Shadow Type" to "No Shadows"), fly around and across
the scene, or do whatever else you want!


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{deg_to_rad}
\alias{deg_to_rad}
\title{Convert decimal degrees to radians}
\usage{
deg_to_rad(deg)
}
\arguments{
\item{deg}{A vector of values, in decimal degrees, to convert to radians}
}
\value{
A vector of the same length in radians
}
\description{
Convert decimal degrees to radians
}
\seealso{
Other utilities: 
\code{\link{addbuff}},
\code{\link{calc_haversine_distance}()},
\code{\link{get_centroid}()},
\code{\link{rad_to_deg}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rad_to_deg}
\alias{rad_to_deg}
\title{Convert radians to degrees}
\usage{
rad_to_deg(rad)
}
\arguments{
\item{rad}{A vector of values, in radians, to convert to decimal degrees}
}
\value{
A vector of the same length in decimal degrees
}
\description{
Convert radians to degrees
}
\seealso{
Other utilities: 
\code{\link{addbuff}},
\code{\link{calc_haversine_distance}()},
\code{\link{deg_to_rad}()},
\code{\link{get_centroid}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{point_from_distance}
\alias{point_from_distance}
\title{Find latitude and longitude for a certain distance and azimuth from a point.}
\usage{
point_from_distance(
  coord_pair,
  distance,
  azimuth,
  distance_unit = "meters",
  azimuth_unit = c("degrees", "radians")
)
}
\arguments{
\item{coord_pair}{A numeric vector of length 2 with names "lat" and "lng"}

\item{distance}{A distance (in meters) representing the distance away from
the original point to apply}

\item{azimuth}{A azimuth (in units specified in \code{azimuth_unit})
representing the direction to apply the distance from the original point in}

\item{distance_unit}{A string passed to [convert_distance]
indicating the units of the provided distance.}

\item{azimuth_unit}{A string (either \code{degrees} or \code{radians})
indicating the units of the \code{azimuth} argument}
}
\value{
An object of class [terrainr_coordinate_pair].
}
\description{
Find latitude and longitude for a certain distance and azimuth from a point.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hit_api.R
\name{hit_national_map_api}
\alias{hit_national_map_api}
\title{Hit the USGS 3DEP API and retrieve an elevation heightmap}
\usage{
hit_national_map_api(
  bbox,
  img_width,
  img_height,
  service,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{bbox}{An object from \link[sf:st_bbox]{sf::st_bbox}.}

\item{img_width}{The number of pixels in the x direction to retrieve}

\item{img_height}{The number of pixels in the y direction to retrieve}

\item{service}{A string indicating what service API to use. For a full list
of available services, see \link{get_tiles}. Short codes are not accepted by this
function.}

\item{verbose}{Logical: Print out the number of tries required to pull each
tile? Default \code{FALSE}.}

\item{...}{Additional arguments passed to the National Map API.
These can be used to change default query parameters or as additional options
for the National Map services. See below for more information.}
}
\value{
A raw vector.
}
\description{
This function retrieves a single tile of data from a single National Map
service and returns the raw response. End users are recommended to use
\link{get_tiles} instead, as it does much more validation and provides
a more friendly interface. For a description of the datasets provided by the
National Map, see \url{https://apps.nationalmap.gov/services}
}
\section{Additional Arguments}{

The \code{...} argument can be used to pass additional arguments to the
National Map API or to edit the hard-coded defaults used by this function.
Some of the most useful options that can be changed include:
\itemize{
\item \code{bboxSR}: The spatial reference of the bounding box given to this function.
If not specified, assumed to be
\href{https://spatialreference.org/ref/epsg/wgs-84/}{4326}.
\item \code{imageSR}: The spatial reference of the image downloaded.
If not specified, assumed to be
\href{https://spatialreference.org/ref/epsg/wgs-84/}{4326}.
\item layers: Which data layers to download. If the National Map API returns data
without specifying layers, this argument isn't used by default. When the
National Map requires this argument, the default value is 0.
\item format: The image format to be downloaded. Defaults depend on the service
being used and are set to be compatible with \link{get_tiles}.
}

Pass these arguments to \code{hit_national_map_api} like you would any other
argument to substitute new values. Note that \code{...} values are never
validated before being used; passing invalid parameters to \code{...} will
cause data retrieval to fail.
}

\examples{
\dontrun{
hit_national_map_api(
  bbox = list(
    c(lat = 44.10438, lng = -74.01231),
    c(lat = 44.17633, lng = -73.91224)
  ),
  img_width = 8000,
  img_height = 8000,
  service = "3DEPElevation"
)
}

}
\seealso{
\link{get_tiles} for a friendlier interface to the National
Map API.

Other data retrieval functions: 
\code{\link{get_tiles}()}
}
\concept{data retrieval functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{terrainr_coordinate_pair}
\alias{terrainr_coordinate_pair}
\title{Construct a terrainr_coordinate_pair object.}
\usage{
terrainr_coordinate_pair(coords, coord_units = c("degrees", "radians"))
}
\arguments{
\item{coords}{A vector of length 2 containing a latitude and longitude. If
unnamed, coordinates are assumed to be in (latitude, longitude) format; if
named, the function will attempt to figure out which value represents which
coordinate. Currently this function understands "lat", "latitude", and "y" as
names for latitude and "lng", "long", "longitude", and "x" for longitude.}

\item{coord_units}{String indicating whether coordinates are in degrees or
radians. Degrees stored in radians will be converted to degrees.}
}
\value{
\code{terrainr_coordinate_pair} object
}
\description{
In order to simplify code, most \code{terrainr} functions expect a set S4
class representation of coordinate pairs and bounding boxes. If the provided
data isn't in the expected S4 format, these functions are used to cast the
data into the target class.
}
\seealso{
Other classes and related functions: 
\code{\link{terrainr_coordinate_pair-class}}
}
\concept{classes and related functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terrainr-package.R
\docType{package}
\name{terrainr-package}
\alias{terrainr}
\alias{terrainr-package}
\title{terrainr: Landscape Visualizations in R and 'Unity'}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Functions for the retrieval, manipulation, and visualization of 'geospatial' data, with an aim towards producing '3D' landscape visualizations in the 'Unity' '3D' rendering engine. Functions are also provided for retrieving elevation data and base map tiles from the 'USGS' National Map <https://apps.nationalmap.gov/services/>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/terrainr/}
  \item \url{https://github.com/ropensci/terrainr}
  \item Report bugs at \url{https://github.com/ropensci/terrainr/issues}
}

}
\author{
\strong{Maintainer}: Michael Mahoney \email{mike.mahoney.218@gmail.com} (\href{https://orcid.org/0000-0003-2402-304X}{ORCID})

Other contributors:
\itemize{
  \item Mike Johnson (Mike reviewed the package (v. 0.2.1) for rOpenSci, see <https://github.com/ropensci/software-review/issues/416>) [reviewer]
  \item Sydney Foks (Sydney reviewed the package (v. 0.2.1) for rOpenSci, see <https://github.com/ropensci/software-review/issues/416>) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_spatial_rgb.R
\docType{data}
\name{geom_spatial_rgb}
\alias{geom_spatial_rgb}
\alias{StatSpatialRGB}
\alias{stat_spatial_rgb}
\title{Plot RGB rasters in ggplot2}
\usage{
geom_spatial_rgb(
  mapping = NULL,
  data = NULL,
  stat = "spatialRGB",
  position = "identity",
  ...,
  hjust = 0.5,
  vjust = 0.5,
  interpolate = FALSE,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  scale = NULL
)

stat_spatial_rgb(
  mapping = NULL,
  data = NULL,
  geom = "raster",
  position = "identity",
  na.rm = FALSE,
  show.legend = FALSE,
  inherit.aes = TRUE,
  scale = NULL,
  ...
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[ggplot2:aes]{aes()}} or
\code{\link[ggplot2:aes_]{aes_()}}. If specified and \code{inherit.aes = TRUE} (the
default), it is combined with the default mapping at the top level of the
plot. You must supply \code{mapping} if there is no plot mapping.}

\item{data}{The data to be displayed in this layer. In addition to the three
options described in [ggplot2::geom_raster], there are two additional
methods:

If a `RasterStack` object (see [raster::stack]), this function will coerce
the stack to a data frame and assume the raster bands are in RGB order
(while allowing for, but ignoring, a fourth alpha band).

If a length-1 character vector, this function will attempt to load the object
via [raster::stack].}

\item{stat}{The statistical transformation to use on the data for this
layer, as a string.}

\item{position}{Position adjustment, either as a string, or the result of
a call to a position adjustment function.}

\item{...}{Other arguments passed on to \code{\link[ggplot2:layer]{layer()}}. These are
often aesthetics, used to set an aesthetic to a fixed value, like
\code{colour = "red"} or \code{size = 3}. They may also be parameters
to the paired geom/stat.}

\item{hjust}{horizontal and vertical justification of the grob.  Each
justification value should be a number between 0 and 1.  Defaults to 0.5
for both, centering each pixel over its data location.}

\item{vjust}{horizontal and vertical justification of the grob.  Each
justification value should be a number between 0 and 1.  Defaults to 0.5
for both, centering each pixel over its data location.}

\item{interpolate}{If \code{TRUE} interpolate linearly, if \code{FALSE}
(the default) don't interpolate.}

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with
a warning. If \code{TRUE}, missing values are silently removed.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[ggplot2:borders]{borders()}}.}

\item{scale}{Integer. Maximum (possible) value in the three channels.
If `NULL`, attempts to infer proper values from data -- if all RGB values
are <= 1 then 1, <= 255 then 255, and otherwise 65535.}

\item{geom}{The geometric object to use display the data}
}
\description{
`geom_spatial_rgb` and `stat_spatial_rgb` allow users to plot three-band RGB
rasters in [ggplot2], using these layers as background base maps for other
spatial plotting. Note that unlike [ggplot2::geom_sf], this function does
_not_ force [ggplot2::coord_sf]; for accurate mapping, add
[ggplot2::coord_sf] with a `crs` value matching your input raster as a layer.
}
\examples{
\dontrun{

simulated_data <- data.frame(
  id = seq(1, 100, 1),
  lat = runif(100, 44.04905, 44.17609),
  lng = runif(100, -74.01188, -73.83493)
)

simulated_data <- sf::st_as_sf(simulated_data, coords = c("lng", "lat"))
simulated_data <- sf::st_set_crs(simulated_data, 4326)

output_tiles <- get_tiles(simulated_data,
  services = c("ortho"),
  resolution = 120
)

merged_ortho <- tempfile(fileext = ".tif")
merge_rasters(output_tiles[["ortho"]], merged_ortho)

merged_stack <- raster::stack(merged_ortho)

library(ggplot2)

ggplot() +
  geom_spatial_rgb(
    data = merged_ortho,
    mapping = aes(
      x = x,
      y = y,
      r = red,
      g = green,
      b = blue
    )
  ) +
  geom_sf(data = simulated_data) +
  coord_sf(crs = 4326)

ggplot() +
  geom_spatial_rgb(
    data = merged_stack,
    mapping = aes(
      x = x,
      y = y,
      r = red,
      g = green,
      b = blue
    )
  ) +
  geom_sf(data = simulated_data) +
  coord_sf(crs = 4326)
}

}
\seealso{
Other visualization functions: 
\code{\link{combine_overlays}()},
\code{\link{raster_to_raw_tiles}()},
\code{\link{vector_to_overlay}()}
}
\concept{visualization functions}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_centroid}
\alias{get_centroid}
\title{Get the great-circle centroid for latitude/longitude data}
\usage{
get_centroid(lat, lng)
}
\arguments{
\item{lat}{A vector of latitudes in degrees.}

\item{lng}{A vector of longitudes in degrees.}
}
\value{
A latitude/longitude
}
\description{
Get the great-circle centroid for latitude/longitude data
}
\seealso{
Other utilities: 
\code{\link{addbuff}},
\code{\link{calc_haversine_distance}()},
\code{\link{deg_to_rad}()},
\code{\link{rad_to_deg}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raster_to_raw_tiles.R
\name{raster_to_raw_tiles}
\alias{raster_to_raw_tiles}
\title{Crop a raster and convert the output tiles into new formats.}
\usage{
raster_to_raw_tiles(input_file, output_prefix, side_length = 4097, raw = TRUE)
}
\arguments{
\item{input_file}{File path to the input TIFF file to convert.}

\item{output_prefix}{The file path to prefix output tiles with.}

\item{side_length}{The side length, in pixels, for the .raw tiles.}

\item{raw}{Logical: Convert the cropped tiles to .raw? When \code{FALSE}
returns a .png.}
}
\value{
Invisibly, a character vector containing the file paths that were
written to.
}
\description{
This function has been deprecated as of terrainr 0.5.0 in favor of the new
function, [make_manifest]. While it will be continued to be exported until
at least 2022, improvements and bug fixes will only be made to the new
function. Please open an issue if any features you relied upon is
missing from the new function!
}
\details{
This function crops input raster files into smaller square tiles and then
converts them into either .png or .raw files which are ready to be imported
into the Unity game engine.
}
\examples{
\dontrun{
if (!isTRUE(as.logical(Sys.getenv("CI")))) {
  simulated_data <- data.frame(
    id = seq(1, 100, 1),
    lat = runif(100, 44.04905, 44.17609),
    lng = runif(100, -74.01188, -73.83493)
  )
  simulated_data <- sf::st_as_sf(simulated_data, coords = c("lng", "lat"))
  output_files <- get_tiles(simulated_data)
  temptiff <- tempfile(fileext = ".tif")
  merge_rasters(output_files["elevation"][[1]], temptiff)
  raster_to_raw_tiles(temptiff, tempfile())
}
}

}
\seealso{
Other data manipulation functions: 
\code{\link{combine_overlays}()},
\code{\link{georeference_overlay}()},
\code{\link{merge_rasters}()},
\code{\link{vector_to_overlay}()}

Other visualization functions: 
\code{\link{combine_overlays}()},
\code{\link{geom_spatial_rgb}()},
\code{\link{vector_to_overlay}()}
}
\concept{data manipulation functions}
\concept{visualization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_bbox_buffer.R
\name{addbuff}
\alias{addbuff}
\alias{add_bbox_buffer}
\alias{add_bbox_buffer.sf}
\alias{add_bbox_buffer.Raster}
\alias{set_bbox_side_length}
\alias{set_bbox_side_length.sf}
\alias{set_bbox_side_length.Raster}
\title{Add a uniform buffer around a bounding box for geographic coordinates}
\usage{
add_bbox_buffer(data, distance, distance_unit = "meters", error_crs = NULL)

\method{add_bbox_buffer}{sf}(data, distance, distance_unit = "meters", error_crs = NULL)

\method{add_bbox_buffer}{Raster}(data, distance, distance_unit = "meters", error_crs = NULL)

set_bbox_side_length(
  data,
  distance,
  distance_unit = "meters",
  error_crs = NULL
)

\method{set_bbox_side_length}{sf}(
  data,
  distance,
  distance_unit = "meters",
  error_crs = NULL
)

\method{set_bbox_side_length}{Raster}(
  data,
  distance,
  distance_unit = "meters",
  error_crs = NULL
)
}
\arguments{
\item{data}{The original data to add a buffer around. Must be either an `sf`
or `Raster` object.}

\item{distance}{The distance to add or to set side lengths equal to.}

\item{distance_unit}{The units of the distance to add to the buffer, passed
to [units::as_units].}

\item{error_crs}{Logical: Should this function error if `data` has no CRS?
If `TRUE`, function errors; if `FALSE`, function quietly assumes EPSG:4326.
If `NULL`, the default, function assumes EPSG:4326 with a warning.}
}
\value{
An `sfc` object (from [sf::st_as_sfc]).
}
\description{
[add_bbox_buffer] calculates the great circle distance both corners of
your bounding box are from the centroid and extends those by a set distance.
Due to using Haversine/great circle distance, latitude/longitude calculations
will not be exact.

[set_bbox_side_length] is a thin wrapper around [add_bbox_buffer] which sets
all sides of the bounding box to (approximately) a specified length.

Both of these functions are intended to be used with geographic coordinate
systems (data using longitude and latitude for position). For projected
coordinate systems, a more sane approach is to use [sf::st_buffer] to add a
buffer, or combine [sf::st_centroid] with the buffer to set a specific side
length.
}
\examples{

df <- data.frame(
  lat = c(44.04905, 44.17609),
  lng = c(-74.01188, -73.83493)
)

df_sf <- sf::st_as_sf(df, coords = c("lng", "lat"))
df_sf <- sf::st_set_crs(df_sf, 4326)

add_bbox_buffer(df_sf, 10)

df <- data.frame(
  lat = c(44.04905, 44.17609),
  lng = c(-74.01188, -73.83493)
)

df_sf <- sf::st_as_sf(df, coords = c("lng", "lat"))
df_sf <- sf::st_set_crs(df_sf, 4326)

set_bbox_side_length(df_sf, 4000)
}
\seealso{
Other utilities: 
\code{\link{calc_haversine_distance}()},
\code{\link{deg_to_rad}()},
\code{\link{get_centroid}()},
\code{\link{rad_to_deg}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tiles.R
\name{get_tiles}
\alias{get_tiles}
\alias{get_tiles.sf}
\alias{get_tiles.sfc}
\alias{get_tiles.Raster}
\alias{get_tiles.list}
\title{A user-friendly way to get USGS National Map data tiles for an area}
\usage{
get_tiles(
  data,
  output_prefix = tempfile(),
  side_length = NULL,
  resolution = 1,
  services = "elevation",
  verbose = FALSE,
  georeference = TRUE,
  projected = NULL,
  ...
)

\method{get_tiles}{sf}(
  data,
  output_prefix = tempfile(),
  side_length = NULL,
  resolution = 1,
  services = "elevation",
  verbose = FALSE,
  georeference = TRUE,
  projected = NULL,
  ...
)

\method{get_tiles}{sfc}(
  data,
  output_prefix = tempfile(),
  side_length = NULL,
  resolution = 1,
  services = "elevation",
  verbose = FALSE,
  georeference = TRUE,
  projected = NULL,
  ...
)

\method{get_tiles}{Raster}(
  data,
  output_prefix = tempfile(),
  side_length = NULL,
  resolution = 1,
  services = "elevation",
  verbose = FALSE,
  georeference = TRUE,
  projected = NULL,
  ...
)

\method{get_tiles}{list}(
  data,
  output_prefix = tempfile(),
  side_length = NULL,
  resolution = 1,
  services = "elevation",
  verbose = FALSE,
  georeference = TRUE,
  projected = NULL,
  ...
)
}
\arguments{
\item{data}{An \code{sf} or \code{Raster} object; tiles will be downloaded for the full
extent of the provided object.}

\item{output_prefix}{The file prefix to use when saving tiles.}

\item{side_length}{The length, in meters, of each side of tiles to download.
If \code{NULL}, defaults to the maximum side length permitted by the least
permissive service requested.}

\item{resolution}{How many meters are represented by each pixel? The default
value of 1 means that 1 pixel = 1 meter, while a value of 2 means that
1 pixel = 2 meters, and so on.}

\item{services}{A character vector of services to download data from. Current
options include "3DEPElevation", "USGSNAIPPlus", and "nhd". Users can also
use short codes to download a specific type of data without specifying the
source; current options for short codes include "elevation" (equivalent to
"3DEPElevation"), "ortho" (equivalent to "USGSNAIPPlus), and "hydro" ("nhd").
Short codes are
not guaranteed to refer to the same source across releases. Short codes are
converted to their service name and then duplicates are removed, so any given
source will only be queried once per tile.}

\item{verbose}{Logical: should tile retrieval functions run in verbose mode?}

\item{georeference}{Logical: should tiles be downloaded as PNGs without
georeferencing, or should they be downloaded as georeferenced TIFF files?
This option does nothing when only elevation data is being downloaded.}

\item{projected}{Logical: is \code{data} in a projected coordinate reference
system? If \code{NULL}, the default, inferred from \link[sf:st_is_longlat]{sf::st_is_longlat}.}

\item{...}{Additional arguments passed to \link{hit_national_map_api}.
These can be used to change default query parameters or as additional options
for the National Map services. See below for more details.}
}
\value{
A list of the same length as the number of unique services requested,
containing named vectors of where data files were saved to. Returned
invisibly.
}
\description{
This function splits the area contained within a bounding box into a set of
tiles, and retrieves data from the USGS National map for each tile. As of
version 0.5.0, the method for lists has been deprecated.
}
\section{Available Datasources}{

The following services are currently available
(with short codes in parentheses where applicable). See links for API
documentation.
\itemize{
\item \href{https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer}{3DEPElevation}
(short code: elevation)
\item \href{https://services.nationalmap.gov/arcgis/rest/services/USGSNAIPPlus/MapServer}{USGSNAIPPlus}
(short code: ortho)
\item \href{https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer}{nhd}
(short code: hydro)
\item \href{https://carto.nationalmap.gov/arcgis/rest/services/govunits/MapServer}{govunits}
\item \href{https://carto.nationalmap.gov/arcgis/rest/services/contours/MapServer}{contours}
\item \href{https://carto.nationalmap.gov/arcgis/rest/services/geonames/MapServer}{geonames}
\item \href{https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer}{NHDPlus_HR}
\item \href{https://carto.nationalmap.gov/arcgis/rest/services/structures/MapServer}{structures}
\item \href{https://carto.nationalmap.gov/arcgis/rest/services/transportation/MapServer}{transportation}
\item \href{https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer}{wbd}
("short code": watersheds)
\item \href{https://www.usgs.gov/centers/geosciences-and-environmental-change-science-center/science/global-ecosystems}{ecosystems}
}
}

\section{Additional Arguments}{

The \code{...} argument can be used to pass additional arguments to the
National Map API or to edit the hard-coded defaults used by this function.
More information on common arguments to change can be found in
\link{hit_national_map_api}. Note that \code{...} can also be used to change
the formats returned by the server, but that doing so while using this
function will likely cause the function to error (or corrupt the output
data). To download files in different formats, use \link{hit_national_map_api}.
}

\examples{
\dontrun{
simulated_data <- data.frame(
  id = seq(1, 100, 1),
  lat = runif(100, 44.04905, 44.17609),
  lng = runif(100, -74.01188, -73.83493)
)

simulated_data <- sf::st_as_sf(simulated_data, coords = c("lng", "lat"))

get_tiles(simulated_data, tempfile())
}

}
\seealso{
Other data retrieval functions: 
\code{\link{hit_national_map_api}()}
}
\concept{data retrieval functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vector_to_overlay.R
\name{vector_to_overlay}
\alias{vector_to_overlay}
\title{Turn spatial vector data into an image overlay}
\usage{
vector_to_overlay(
  vector_data,
  reference_raster,
  output_file = NULL,
  transparent = "#ffffff",
  ...,
  error_crs = NULL
)
}
\arguments{
\item{vector_data}{The spatial vector data set to be transformed into an
overlay image. Users may provide either an \code{sf} object or a length 1
character vector containing a path to a file readable by \link[sf:st_read]{sf::read_sf}.}

\item{reference_raster}{The raster file to produce an overlay for. The output
overlay will have the same extent and resolution as the input raster. Users
may provide either a Raster* object or a length 1 character
vector containing a path to a file readable by \link[raster:raster]{raster::raster}.}

\item{output_file}{The path to save the image overlay to. If \code{NULL}, saves to
a tempfile.}

\item{transparent}{The hex code for a color to be made transparent in the
final image. Set to \code{FALSE} to not set any colors to transparent.}

\item{...}{Arguments passed to \code{...} in either \link[ggplot2:geom_point]{ggplot2::geom_point} (for
point vector data), \link[ggplot2:geom_path]{ggplot2::geom_line} (for line data),
or \link[ggplot2:geom_polygon]{ggplot2::geom_polygon} (for all other data types).}

\item{error_crs}{Logical: Should this function error if \code{data} has no CRS?
If \code{TRUE}, function errors; if \code{FALSE}, function quietly assumes EPSG:4326.
If \code{NULL}, the default, function assumes EPSG:4326 with a warning.}
}
\value{
\code{output_file}, invisibly.
}
\description{
This function allows users to quickly transform any vector data into an
image overlay, which may then be imported as a texture into Unity.
}
\examples{
\dontrun{

# Generate points to download raster tiles for
set.seed(123)
simulated_data <- data.frame(
  id = seq(1, 100, 1),
  lat = runif(100, 44.1114, 44.1123),
  lng = runif(100, -73.92273, -73.92147)
)

# Create an sf object from our original simulated data

simulated_data_sf <- sf::st_as_sf(simulated_data, coords = c("lng", "lat"))
sf::st_crs(simulated_data_sf) <- sf::st_crs(4326)

# Download data!

downloaded_tiles <- get_tiles(simulated_data_sf, tempfile())

merged_file <- merge_rasters(
  downloaded_tiles[[1]],
  tempfile(fileext = ".tif")
)


# Create an overlay image
vector_to_overlay(simulated_data_sf, merged_file[[1]], na.rm = TRUE)
}

}
\seealso{
Other data manipulation functions: 
\code{\link{combine_overlays}()},
\code{\link{georeference_overlay}()},
\code{\link{merge_rasters}()},
\code{\link{raster_to_raw_tiles}()}

Other overlay creation functions: 
\code{\link{combine_overlays}()},
\code{\link{georeference_overlay}()}

Other visualization functions: 
\code{\link{combine_overlays}()},
\code{\link{geom_spatial_rgb}()},
\code{\link{raster_to_raw_tiles}()}
}
\concept{data manipulation functions}
\concept{overlay creation functions}
\concept{visualization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{terrainr_bounding_box}
\alias{terrainr_bounding_box}
\title{Construct a terrainr_bounding_box object}
\usage{
terrainr_bounding_box(bl, tr, coord_units = "degrees")
}
\arguments{
\item{bl, tr}{The bottom left (\code{bl}) and top right (\code{tr}) corners of
the bounding box, either as a [terrainr_coordinate_pair] object
or a coordinate pair. If the coordinate pair is not named, it is assumed to
be in (lat, lng) format; if it is named, the function will attempt to
properly identify coordinates.}

\item{coord_units}{Arguments passed to [terrainr_coordinate_pair].
If \code{bl} and \code{tr} are already [terrainr_coordinate_pair]
objects, these arguments are not used.}
}
\value{
An object of class [terrainr_bounding_box].
}
\description{
In order to simplify code, most \code{terrainr} functions expect a set S4
class representation of coordinate pairs and bounding boxes. If the provided
data is not in the expected S4 format, these functions are used to cast the
data into the target class.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{terrainr_bounding_box-class}
\alias{terrainr_bounding_box-class}
\title{S4 class for bounding boxes in the format expected by \code{terrainr}
functions.}
\description{
S4 class for bounding boxes in the format expected by \code{terrainr}
functions.
}
\section{Slots}{

\describe{
\item{\code{bl}}{A \code{\link{terrainr_coordinate_pair}} representing the bottom
left corner of the bounding box}

\item{\code{tr}}{A \code{\link{terrainr_coordinate_pair}} representing the top right
corner of the bounding box}
}}

\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_overlays.R
\name{combine_overlays}
\alias{combine_overlays}
\title{Combine multiple image overlays into a single file}
\usage{
combine_overlays(
  ...,
  output_file = tempfile(fileext = ".png"),
  transparency = 0
)
}
\arguments{
\item{...}{File paths for images to be combined. Note that combining TIFF
images requires the \code{tiff} package be installed.}

\item{output_file}{The path to save the resulting image to. Can
be any format accepted by \link[magick:editing]{magick::image_read}. Optionally, can be set to
\code{NULL}, in which case this function will return the image as a \code{magick}
object instead of writing to disk.}

\item{transparency}{A value indicating how much transparency should be added
to each image. If less than 1, interpreted as a proportion (so a value of
0.1 results in each image becoming 10\% more transparent); if between 1 and
100, interpreted as a percentage (so a value of 10 results in each image
becoming 10\% more transparent.) A value of 0 is equivalent to no
additional transparency.}
}
\value{
If \code{output_file} is not null, \code{output_file}, invisibly. If
\code{output_file} is null, a \code{magick} image object.
}
\description{
This function combines any number of images into a single file, which may
then be further processed as an image or transformed into an image overlay.
}
\examples{
\dontrun{
# Generate points and download orthoimagery
mt_elbert_points <- data.frame(
  lat = runif(100, min = 39.11144, max = 39.12416),
  lng = runif(100, min = -106.4534, max = -106.437)
)

mt_elbert_sf <- sf::st_as_sf(mt_elbert_points, coords = c("lng", "lat"))
sf::st_crs(mt_elbert_sf) <- sf::st_crs(4326)

output_files <- get_tiles(
  mt_elbert_sf,
  output_prefix = tempfile(),
  services = c("ortho")
)

# Merge orthoimagery into a single file
ortho_merged <- merge_rasters(
  input_rasters = output_files[1],
  output_raster = tempfile(fileext = ".tif")
)

# Convert our points into an overlay
mt_elbert_overlay <- vector_to_overlay(mt_elbert_sf,
  ortho_merged[[1]],
  size = 15,
  color = "red",
  na.rm = TRUE
)

# Combine the overlay with our orthoimage
ortho_with_points <- combine_overlays(
  ortho_merged[[1]],
  mt_elbert_overlay
)
}

}
\seealso{
Other data manipulation functions: 
\code{\link{georeference_overlay}()},
\code{\link{merge_rasters}()},
\code{\link{raster_to_raw_tiles}()},
\code{\link{vector_to_overlay}()}

Other overlay creation functions: 
\code{\link{georeference_overlay}()},
\code{\link{vector_to_overlay}()}

Other visualization functions: 
\code{\link{geom_spatial_rgb}()},
\code{\link{raster_to_raw_tiles}()},
\code{\link{vector_to_overlay}()}
}
\concept{data manipulation functions}
\concept{overlay creation functions}
\concept{visualization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_rasters.R
\name{merge_rasters}
\alias{merge_rasters}
\title{Merge multiple raster files into a single raster}
\usage{
merge_rasters(
  input_rasters,
  output_raster = tempfile(fileext = ".tif"),
  options = character(0),
  overwrite = FALSE,
  force_fallback = FALSE
)
}
\arguments{
\item{input_rasters}{A character vector containing the file paths to the
georeferenced rasters you want to use.}

\item{output_raster}{The file path to save the merged georeferenced raster
to.}

\item{options}{Optionally, a character vector of options to be passed
directly to [sf::gdal_utils]. If the fallback is used and any options (other
than "-overwrite") are specified, this will issue a warning.}

\item{overwrite}{Logical: overwrite `output_raster` if it exists? If FALSE
and the file exists, this function will fail with an error. The behavior if
this argument is TRUE and "-overwrite" is passed to `options` directly is
not stable.}

\item{force_fallback}{Logical: if TRUE, uses the much slower fallback method
by default. This is used for testing purposes and is not recommended for use
by end users.}
}
\value{
`output_raster`, invisibly.
}
\description{
Some functions like [get_tiles] return multiple separate files
when it can be useful to have a single larger raster instead. This function
is a thin wrapper over [sf::gdal_utils(util = "warp")], making it easy to
collapse those multiple raster files into a single TIFF.
}
\examples{
\dontrun{
simulated_data <- data.frame(
  lat = c(44.10379, 44.17573),
  lng = c(-74.01177, -73.91171)
)

simulated_data <- sf::st_as_sf(simulated_data, coords = c("lng", "lat"))

img_files <- get_tiles(simulated_data)
merge_rasters(img_files[[1]])
}

}
\seealso{
Other data manipulation functions: 
\code{\link{combine_overlays}()},
\code{\link{georeference_overlay}()},
\code{\link{raster_to_raw_tiles}()},
\code{\link{vector_to_overlay}()}
}
\concept{data manipulation functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{terrainr_coordinate_pair-class}
\alias{terrainr_coordinate_pair-class}
\title{S4 class for coordinate points in the format expected by
\code{terrainr} functions.}
\description{
S4 class for coordinate points in the format expected by
\code{terrainr} functions.
}
\section{Slots}{

\describe{
\item{\code{lat}}{Numeric latitude, in decimal degrees}

\item{\code{lng}}{Numeric longitude, in decimal degrees}
}}

\seealso{
Other classes and related functions: 
\code{\link{terrainr_coordinate_pair}}
}
\concept{classes and related functions}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_haversine_distance.R
\name{calc_haversine_distance}
\alias{calc_haversine_distance}
\title{Extract latitude and longitude from a provided object}
\usage{
calc_haversine_distance(point_1, point_2)
}
\arguments{
\item{point_1, point_2}{Coordinate pairs (as length-2 numeric vectors with the
names "lat" and "lng") to calculate distance between.}
}
\value{
A vector of length 1 containing distance between points
}
\description{
This is an internal utility function to convert bounding boxes into
coordinate pairs.
}
\seealso{
Other utilities: 
\code{\link{addbuff}},
\code{\link{deg_to_rad}()},
\code{\link{get_centroid}()},
\code{\link{rad_to_deg}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_manifest.R
\name{make_manifest}
\alias{make_manifest}
\title{Transform rasters and write manifest file for import into Unity}
\usage{
make_manifest(
  heightmap,
  overlay = NULL,
  output_prefix = "import",
  manifest_path = "terrainr.manifest",
  importer_path = "import_terrain.cs"
)
}
\arguments{
\item{heightmap}{File path to the heightmap to transform.}

\item{overlay}{Optionally, file path to the image overlay to transform.}

\item{output_prefix}{The file path to prefix output tiles with.}

\item{manifest_path}{File path to write the manifest file to.}

\item{importer_path}{File name to write the importer script to. Set to NULL
to not copy the importer script. Will overwrite any file at the same path.}
}
\value{
`manifest_path`, invisibly.
}
\description{
This function crops input raster files into smaller square tiles and then
converts them into either .png or .raw files which are ready to be imported
into the Unity game engine. It also writes a "manifest" file and importer
script which may be used to automatically import the tiles into Unity.
}
\examples{
\dontrun{
if (!isTRUE(as.logical(Sys.getenv("CI")))) {
  simulated_data <- data.frame(
    id = seq(1, 100, 1),
    lat = runif(100, 44.04905, 44.17609),
    lng = runif(100, -74.01188, -73.83493)
  )
  simulated_data <- sf::st_as_sf(simulated_data, coords = c("lng", "lat"))
  output_files <- get_tiles(simulated_data)
  temptiff <- tempfile(fileext = ".tif")
  merge_rasters(output_files["elevation"][[1]], temptiff)
  make_manifest(temptiff, output_prefix = tempfile(), importer_path = NULL)
}
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/georeference_overlay.R
\name{georeference_overlay}
\alias{georeference_overlay}
\title{Georeference image overlays based on a reference raster}
\usage{
georeference_overlay(
  overlay_file,
  reference_raster,
  output_file = tempfile(fileext = ".tif")
)
}
\arguments{
\item{overlay_file}{The image overlay to georeference. File format will be
detected automatically from file extension; options include `jpeg/jpg`,
`png`, and `tif/tiff`.}

\item{reference_raster}{The raster file to base georeferencing on. The output
image will have the same extent and CRS as the reference raster. Accepts both
Raster* objects from the `raster` package or a file readable by
[raster::raster].}

\item{output_file}{The path to write the georeferenced image file to. Must
be a TIFF.}
}
\value{
The file path written to, invisibly.
}
\description{
This function georeferences an image overlay based on a reference raster,
setting the extent and CRS of the image to those of the raster file. To
georeference multiple images and merge them into a single file, see
[merge_rasters].
}
\examples{
\dontrun{
simulated_data <- data.frame(
  id = seq(1, 100, 1),
  lat = runif(100, 44.1114, 44.1123),
  lng = runif(100, -73.92273, -73.92147)
)

simulated_data <- sf::st_as_sf(simulated_data, coords = c("lng", "lat"))

downloaded_tiles <- get_tiles(simulated_data,
  services = c("elevation", "ortho"),
  georeference = FALSE
)

georeference_overlay(
  overlay_file = downloaded_tiles[[2]],
  reference_raster = downloaded_tiles[[1]],
  output_file = tempfile(fileext = ".tif")
)
}

}
\seealso{
Other data manipulation functions: 
\code{\link{combine_overlays}()},
\code{\link{merge_rasters}()},
\code{\link{raster_to_raw_tiles}()},
\code{\link{vector_to_overlay}()}

Other overlay creation functions: 
\code{\link{combine_overlays}()},
\code{\link{vector_to_overlay}()}
}
\concept{data manipulation functions}
\concept{overlay creation functions}
