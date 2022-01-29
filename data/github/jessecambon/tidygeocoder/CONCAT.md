
<!-- README.md is generated from README.Rmd. Please edit that file directly and reknit -->

# tidygeocoder<a href='https://jessecambon.github.io/tidygeocoder/'><img src="man/figures/tidygeocoder_hex.png" align="right" height="139px"/></a>

<!-- badges: start -->

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03544/status.svg)](https://doi.org/10.21105/joss.03544)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/jessecambon/tidygeocoder/blob/master/LICENSE.md)
[![CRAN](https://www.r-pkg.org/badges/version/tidygeocoder)](https://cran.r-project.org/package=tidygeocoder)
[![CRAN Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/tidygeocoder)](https://CRAN.R-project.org/package=tidygeocoder)
[![CRAN Downloads Per
Month](http://cranlogs.r-pkg.org/badges/tidygeocoder)](https://cran.r-project.org/package=tidygeocoder)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R Build
Status](https://github.com/jessecambon/tidygeocoder/workflows/R-CMD-check/badge.svg)](https://github.com/jessecambon/tidygeocoder/actions?workflow=R-CMD-check)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.5627341.svg)](https://doi.org/10.5281/zenodo.5627341)
<!-- badges: end -->

Tidygeocoder makes getting data from geocoding services easy. A unified
high-level interface is provided for a selection of [supported geocoding
services](https://jessecambon.github.io/tidygeocoder/articles/geocoder_services.html)
and results are returned in [tibble](https://tibble.tidyverse.org/)
(dataframe) format.

**Features:**

-   Forward geocoding (addresses ⮕ coordinates)
-   Reverse geocoding (coordinates ⮕ addresses)
-   Batch geocoding (geocoding multiple addresses or coordinates in a
    single query) is automatically used if applicable.
-   Duplicate, NA, and blank input data is handled elegantly; only
    unique inputs are submitted in queries, but the rows in the original
    data are preserved by default.
-   The maximum rate of querying is automatically set according to the
    usage policies of the selected geocoding service.

In addition to the usage examples below, see the [Getting Started
Vignette](https://jessecambon.github.io/tidygeocoder/articles/tidygeocoder.html)
and [blog posts on
tidygeocoder](https://jessecambon.github.io/tag/tidygeocoder).

## Installation

To install the stable version from CRAN (the official R package
servers):

``` r
install.packages('tidygeocoder')
```

Alternatively, you can install the latest development version from
GitHub:

``` r
devtools::install_github("jessecambon/tidygeocoder")
```

## Usage

In this first example we will geocode a few addresses using the
`geocode()` function and plot them on a map with ggplot.

``` r
library(dplyr, warn.conflicts = FALSE)
library(tidygeocoder)

# create a dataframe with addresses
some_addresses <- tibble::tribble(
~name,                  ~addr,
"White House",          "1600 Pennsylvania Ave NW, Washington, DC",
"Transamerica Pyramid", "600 Montgomery St, San Francisco, CA 94111",     
"Willis Tower",         "233 S Wacker Dr, Chicago, IL 60606"                                  
)

# geocode the addresses
lat_longs <- some_addresses %>%
  geocode(addr, method = 'osm', lat = latitude , long = longitude)
#> Passing 3 addresses to the Nominatim single address geocoder
#> Query completed in: 3 seconds
```

The `geocode()` function geocodes addresses contained in a dataframe.
The [Nominatim (“osm”)](https://nominatim.org/) geocoding service is
used here, but other services can be specified with the `method`
argument. Only latitude and longitude are returned from the geocoding
service in this example, but `full_results = TRUE` can be used to return
all of the data from the geocoding service. See the `geo()` function
documentation for details.

| name                 | addr                                       | latitude |  longitude |
|:---------------------|:-------------------------------------------|---------:|-----------:|
| White House          | 1600 Pennsylvania Ave NW, Washington, DC   | 38.89770 |  -77.03655 |
| Transamerica Pyramid | 600 Montgomery St, San Francisco, CA 94111 | 37.79520 | -122.40279 |
| Willis Tower         | 233 S Wacker Dr, Chicago, IL 60606         | 41.87535 |  -87.63576 |

Now that we have the longitude and latitude coordinates, we can use
ggplot to plot our addresses on a map.

``` r
library(ggplot2)

ggplot(lat_longs, aes(longitude, latitude), color = "grey99") +
  borders("state") + geom_point() +
  ggrepel::geom_label_repel(aes(label = name)) +
  theme_void()
```

<img src="man/figures/README-usamap-1.png" style="display: block; margin: auto;" />

To perform reverse geocoding (obtaining addresses from geographic
coordinates), we can use the `reverse_geocode()` function. The arguments
are similar to the `geocode()` function, but now we specify the input
data columns with the `lat` and `long` arguments. The input dataset used
here is the results of the geocoding query above.

The single line address is returned in a column named by the `address`
argument and all columns from the geocoding service results are returned
because `full_results = TRUE`. See the `reverse_geo()` function
documentation for more details.

<!-- 
Removing the licence column is done just to prevent a note from 
occurring in automated CRAN checks for an improper/old link.
-->

``` r
reverse <- lat_longs %>%
  reverse_geocode(lat = latitude, long = longitude, method = 'osm',
                  address = address_found, full_results = TRUE) %>%
  select(-addr, -licence)
#> Passing 3 coordinates to the Nominatim single coordinate geocoder
#> Query completed in: 3 seconds
```

| name                 | latitude |  longitude | address\_found                                                                                                                                         | place\_id | osm\_type |   osm\_id | osm\_lat           | osm\_lon            | office      | house\_number | road                          | city          | state                | postcode | country       | country\_code | boundingbox                                          | tourism              | neighbourhood | county        | suburb |
|:---------------------|---------:|-----------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|:----------|----------:|:-------------------|:--------------------|:------------|:--------------|:------------------------------|:--------------|:---------------------|:---------|:--------------|:--------------|:-----------------------------------------------------|:---------------------|:--------------|:--------------|:-------|
| White House          | 38.89770 |  -77.03655 | White House, 1600, Pennsylvania Avenue Northwest, Washington, District of Columbia, 20500, United States                                               | 159983331 | way       | 238241022 | 38.897699700000004 | -77.03655315        | White House | 1600          | Pennsylvania Avenue Northwest | Washington    | District of Columbia | 20500    | United States | us            | 38.8974908 , 38.897911 , -77.0368537, -77.0362519    | NA                   | NA            | NA            | NA     |
| Transamerica Pyramid | 37.79520 | -122.40279 | Transamerica Pyramid, 600, Montgomery Street, Chinatown, San Francisco, San Francisco City and County, San Francisco, California, 94111, United States | 106008002 | way       |  24222973 | 37.795200550000004 | -122.40279267840137 | NA          | 600           | Montgomery Street             | San Francisco | California           | 94111    | United States | us            | 37.7948854 , 37.7954472 , -122.4031399, -122.4024317 | Transamerica Pyramid | Chinatown     | San Francisco | NA     |
| Willis Tower         | 41.87535 |  -87.63576 | South Wacker Drive, Printer’s Row, Loop, Chicago, Cook County, Illinois, 60606, United States                                                          | 182238096 | way       | 337681342 | 41.8753503         | -87.6357587         | NA          | NA            | South Wacker Drive            | Chicago       | Illinois             | 60606    | United States | us            | 41.8749718 , 41.8757997 , -87.6361005, -87.6354602   | NA                   | Printer’s Row | Cook County   | Loop   |

## In the Wild

For inspiration, here are a few articles (with code) that leverage
tidygeocoder:

-   [Exercises: Spatial Data Wrangling with
    sf](http://www2.stat.duke.edu/courses/Spring21/sta323.001/exercises/lec_12.html) -
    part of a [statistical computing
    course](http://www2.stat.duke.edu/courses/Spring21/sta323.001/) at
    Duke
-   [Geocoding the Minard
    Map](https://www.jla-data.net/eng/minard-map-tidygeocoder/) -
    recreating a famous infographic with geocoding
-   [Mapping a network of women in
    demography](https://www.monicaalexander.com/posts/2021-21-02-mapping/) -
    using rvest and tidygeocoder to map Google Scholar data
-   [Mapping
    Routes](https://bensstats.wordpress.com/2021/10/21/robservations-15-i-reverse-engineered-atlas-co-well-some-of-it/) -
    mapping routes with tidygeocoder and osrm
-   [Road Routing in
    R](https://www.jla-data.net/eng/routing-in-r-context/) -
    demonstration of three different routing APIs
-   [Mapping Texas Ports With
    R](https://www.sharpsightlabs.com/blog/mapping-texas-ports-with-r-part1/) -
    mapping the Texas coast with rnaturalearth and sf

## Contributing

Contributions to the tidygeocoder package are welcome. File [an
issue](https://github.com/jessecambon/tidygeocoder/issues) for bug fixes
or suggested features. If you would like to contribute code such as
adding support for a new geocoding service, reference the [developer
notes](https://jessecambon.github.io/tidygeocoder/articles/developer_notes.html)
for instructions and documentation.

## Citing tidygeocoder

Use the `citation()` function:

``` r
citation('tidygeocoder')
```

</br>

<blockquote>


    To cite tidygeocoder use:

      Cambon J, Hernangómez D, Belanger C, Possenriede D (2021).
      tidygeocoder: An R package for geocoding. Journal of Open Source
      Software, 6(65), 3544, https://doi.org/10.21105/joss.03544 (R package
      version 1.0.5)

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {tidygeocoder: An R package for geocoding},
        author = {Jesse Cambon and Diego Hernangómez and Christopher Belanger and Daniel Possenriede},
        year = {2021},
        journal = {Journal of Open Source Software},
        publisher = {The Open Journal},
        doi = {10.21105/joss.03544},
        url = {https://doi.org/10.21105/joss.03544},
        volume = {6},
        number = {65},
        pages = {3544},
        note = {R package version 1.0.5},
      }

</blockquote>

Or refer to the [citation
page](https://jessecambon.github.io/tidygeocoder/authors.html).
# tidygeocoder 1.0.5

- Corrected documentation for the `quiet` parameter in `geo()` and `reverse_geo()`.
- To fix an issue that occurred on Mac CRAN checks ([#152](https://github.com/jessecambon/tidygeocoder/issues/152)) the vignette is now precomputed so that it does not run during `R CMD check` (ie. `devtools::check()`).
- Changed some default parameter values to facilitate the use of the [memoise package](https://memoise.r-lib.org/) ([#154](https://github.com/jessecambon/tidygeocoder/pull/154), [@dpprdan](https://github.com/dpprdan)).

# tidygeocoder 1.0.4

### New Features

- Added support for the [Geoapify](https://www.geoapify.com/) service (thanks [@dpprdan](https://github.com/dpprdan)). This service supports batch geocoding, but this capability is currently not implemented in tidygeocoder (see [#119](https://github.com/jessecambon/tidygeocoder/issues/119)).
- Added the functions `geocode_combine()` and `geo_combine()` to combine the results of multiple geocoding queries. These functions are meant to replace `geo(method = "cascade")`.
- Deprecated `method = "cascade"` and the `cascade_order`, `param_error`, and `batch_limit_error` arguments for `geo()`.
- Deprecated the `return_type`, `geocodio_v`, `mapbox_permanent`, `mapquest_open`, `iq_region`, and `here_request_id` arguments in favor of the new `api_options` parameter for the `geo()` and `reverse_geo()` functions.
- Added a `api_options = list(geocodio_hipaa = TRUE/FALSE))` option to `geo()` and `reverse_geo()` functions to allow the toggling of the HIPAA-compliant Geocodio API endpoint ([#137](https://github.com/jessecambon/tidygeocoder/issues/137)).

### Console Output

- Added a progress bar for single input geocoding (ie. not batch geocoding) ([#38](https://github.com/jessecambon/tidygeocoder/issues/38)). The progress bar can be disabled with `progress_bar = FALSE` (a new parameter for the `geo()` and `reverse_geo()` functions). Similar to the [readr](https://readr.tidyverse.org/reference/show_progress.html) package, progress bars are shown by default if the session is interactive and the code being run is not part of an RStudio Notebook chunk or R Markdown knitting process.
- Some console messages related to geocoding queries are now shown by default. These messages show the number of inputs (addresses or coordinates) submitted, the geocoding service used, and how long the query took to execute. To suppress these messages you can set `quiet = TRUE` (a new parameter for the `geo()` and `reverse_geo()` functions).
- Default arguments with `options()` for `verbose`, `quiet`, and `progress_bar`. For instance `options(tidygeocoder.verbose = TRUE)` changes the default value of `verbose` from FALSE to TRUE.

### Bugfixes

- Fixed a bug for Bing forward geocoding `geo()` when no results are found ([#112](https://github.com/jessecambon/tidygeocoder/issues/112)).
- Fixed a bug that occurred in reverse geocoding when passing a set of exclusively duplicate coordinates (ie. 1 unique coordinate) ([#129](https://github.com/jessecambon/tidygeocoder/issues/129)).


# tidygeocoder 1.0.3

### New Features

- Added support for reverse geocoding with the new `reverse_geo()` and `reverse_geocode()` functions. 
- Added support for the [OpenCage](https://opencagedata.com/) geocoding service ([#67](https://github.com/jessecambon/tidygeocoder/issues/67)) (thanks [@dpprdan](https://github.com/dpprdan)).
- Added support for the [HERE](https://developer.here.com/products/geocoding-and-search) ([#74](https://github.com/jessecambon/tidygeocoder/issues/74)), [Mapbox](https://docs.mapbox.com/api/search/) ([#71](https://github.com/jessecambon/tidygeocoder/issues/71)), [MapQuest](https://developer.mapquest.com/documentation/geocoding-api/) ([#85](https://github.com/jessecambon/tidygeocoder/issues/85)),  [TomTom](https://developer.tomtom.com/search-api/search-api-documentation/geocoding) ([#76](https://github.com/jessecambon/tidygeocoder/issues/76)), [Bing](https://docs.microsoft.com/en-us/bingmaps/rest-services/locations/) ([#92](https://github.com/jessecambon/tidygeocoder/issues/92)), and [ArcGIS](https://developers.arcgis.com/rest/geocode/api-reference/overview-world-geocoding-service.htm) ([#98](https://github.com/jessecambon/tidygeocoder/issues/98)) geocoding services (thanks [@dieghernan](https://github.com/dieghernan)). Note that the batch geocoding capabilities for the Mapbox and ArcGIS services are not currently implemented (see [#73](https://github.com/jessecambon/tidygeocoder/issues/73) and [#102](https://github.com/jessecambon/tidygeocoder/issues/102)).
- Added the ` mapbox_permanent`, `here_request_id`, and `mapquest_open` parameters to the `geo()` and `reverse_geo()` functions.
- The `limit` argument can now be used with the "google" and "census" methods to control the number of results returned. These two services do not have limit arguments in their APIs so the limit is applied after the results are returned.
- `batch_limit` is now automatically set according to the specified geocoding service unless otherwise specified.

### Other Changes

- Changed default `method` to `"osm"` (Nominatim) for the `geo()` function (it was previously `"census"`).
- The `geo_<method>` functions are now deprecated.
- Added the `return_input` argument to `geocode()` and `reverse_geocode()` to provide more flexibility when using dataframes as inputs in geocoder queries.
- `limit = NULL` can now be passed to use the default `limit` value for the geocoding service.
- Added the `min_time_reference`, `batch_limit_reference`, `api_key_reference`, and `api_info_reference` datasets to more accessibly store values for `min_time`, `batch_limit`, the names of environmental variables for API keys, and information for documentation (such as API documentation links).
- `geocode()` and `reverse_geocode()` now require `limit = 1` (default) unless `return_input = FALSE`. This fixes a bug where geocoding results could be misaligned with the input dataset when `limit > 1`.  ([#88](https://github.com/jessecambon/tidygeocoder/issues/88)).
- An error is now thrown by default if the number of inputs exceeds the batch query limit. For forward geocoding, this behavior can be toggled with the new `batch_limit_error` argument in the `geo()` function and `batch_limit_error` is set to FALSE if `method = "cascade"`. When `batch_limit_error` is FALSE then the batch query size is limited to the batch limit and executed (which was the default behavior in the previous version).
- The `address_list` argument of `query_api()` has been renamed to `input_list` to reflect that it is used for both forward and reverse queries when using the Geocodio service for batch geocoding.
- The `query_api()` function now returns a named list which contains the response content and the HTTP status code. The `geo()` and `reverse_geo()` functions now use the HTTP status code directly to determine if a response is valid.
- Added [external tests](https://github.com/jessecambon/tidygeocoder/blob/main/external/online_tests.R) to more thoroughly test the package with live queries (internal package tests don't run queries).
- Added functions to generate package documentation from built-in datasets (ex. the methods documentation in `geo()` and `reverse_geo()`).
- Converted package documentation from standard roxygen syntax to Markdown.

# tidygeocoder 1.0.2

- Added support for the [Google](https://developers.google.com/maps/documentation/geocoding/overview) geocoding service ([#34](https://github.com/jessecambon/tidygeocoder/issues/34)) (thanks [@chris31415926535](https://github.com/chris31415926535)).
- An error is now thrown if invalid parameters are passed to geocoding services (the parameters checked are limit, address, street, city, county, state, postalcode, and country) ([#53](https://github.com/jessecambon/tidygeocoder/issues/53)). This behavior can be toggled with the new `param_error` parameter in `geo()` (or `geocode()`).
- Leading zeros on Census FIPs geography columns are now preserved ([\#47](https://github.com/jessecambon/tidygeocoder/issues/47)).
- Bug fix for `custom_query` argument with Geocodio batch geocoding ([\#48](https://github.com/jessecambon/tidygeocoder/issues/48)).
- Bug fix for vctrs datatype error with cascade method ([\#49](https://github.com/jessecambon/tidygeocoder/issues/49)).
- Added more comprehensive testing for internal package functions such as package_addresses, unpackage_addresses, and get_api_query ([#58](https://github.com/jessecambon/tidygeocoder/issues/58)).
- Per CRAN request, `order()` is no longer called on data frames ([\#57](https://github.com/jessecambon/tidygeocoder/issues/57)).

# tidygeocoder 1.0.1

-   Fixed an issue that prevented installation on R \< 4.0. ([\#35](https://github.com/jessecambon/tidygeocoder/issues/35)).
-   Updated package documentation. Added examples to utility functions `query_api()` and `get_api_query()`.

# tidygeocoder 1.0.0

### New Functionality

- **New geocoding services**: Support for the [Geocodio](https://www.geocod.io/) and [Location IQ](https://locationiq.com/) services has been added.
- **Batch geocoding** (geocoding multiple addresses per query) is now available for both the Census and Geocodio services.
- **Full results** from the geocoding services can now be returned by using `full_results = TRUE`. This will return all data provided by the geocoding service instead of just latitude and longitude coordinates. Additionally, the `return_type = 'geographies'` argument for the Census geocoder will return geography columns.
- **Address component arguments**: As an alternative to specifying a single-line address, address component arguments are now supported (`street`, `city`, `county`, `postalcode`, `country`).
- **Customizable queries**: Geocoding queries can now be customized using the `limit` and `custom_query` arguments (see the `geo()` function for details).
- **Smart address handling**: Only unique addresses are passed to geocoding services, but the rows in the original data are preserved.
- **Usage limits**: The OSM and IQ services by default are now limited to submitting one query per second (per the `min_time` argument in `geo()`) to respect usage limits. This should fix the past issue of users being locked out of the OSM service due to usage limit violations.
- The `cascade` method can now be customized by using the `cascade_order` argument (see `geo()` documentation)
- **Custom API URLs** can now be specified. This will allow users to specify their own local Nominatim server, for instance.
- The parameters passed to the geocoding service are now displayed to the console when `verbose = TRUE`

### Under the Hood Improvements

- **Reduced dependencies**: The package has been overhauled so that the only remaining dependencies are `tibble`, `dplyr`, `jsonlite`, and `httr`. The package no longer has direct dependencies on `tmaptools`, `stringr`, `purrr`, `tidyr`, and `rlang`.
- All geocoding queries are now directly executed by `httr`. The inbuilt `api_parameter_reference` dataset is used to map standard "generic" parameter names to the parameter names used by each specific geocoding service.
- All geocoding functionality has been centralized in the `geo()` function. Users can still use `geocode()`, `geo_osm()`, and `geo_census()` as before. However, `geo_osm()` and `geo_census()` are now just convenience functions that call `geo()` and `geocode()` passes all addresses to `geo()` for geocoding.

# tidygeocoder 0.2.5

Per CRAN request, fixed an issue where the example for `R/geocode.R` failed. The only change required was to add a `library(dplyr)` statement.

# tidygeocoder 0.2.4

Initial CRAN release. Per CRAN request:

- Replaced `print()` with `warning()` to make suppressing console output possible.
- Replaced `\dontrun` with `\donttest` in .R files## Test environments
* local Ubuntu install, R 4.1.0
* GitHub Actions [R-CMD-check](https://github.com/jessecambon/tidygeocoder/blob/main/.github/workflows/check-full.yaml)
* winbuilder r-old devel: `devtools::check_win_oldrelease()`
* winbuilder r-devel : `devtools::check_win_devel()`
* Other environments checked via `rhub::check_for_cran()`

## R CMD check results

0 errors | 0 warnings | 1 notes

## Notes 

- This release was required to prevent the package from being removed from CRAN due to CRAN check issues that occured without network connectivity.
- RE the inactive Zenodo URL, this DOI is reserved and will become active once released.---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

Describe the bug and include a small reproducible example: https://www.tidyverse.org/help/ 

The datapasta package can be useful for including data in the reproducible example (see the `tribble_paste()` function): https://milesmcbain.github.io/datapasta/

It can also be helpful to post the results of `devtools::session_info()` to show information about your R session including what package versions you are using.
