
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
