
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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file directly and reknit -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  fig.width = 8,
  fig.height = 5,
  fig.align = 'center'
)
options(tibble.print_min = 5, tibble.print_max = 5)
```

# tidygeocoder<a href='https://jessecambon.github.io/tidygeocoder/'><img src="man/figures/tidygeocoder_hex.png" align="right" height="139px"/></a>

<!-- badges: start -->
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03544/status.svg)](https://doi.org/10.21105/joss.03544)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/jessecambon/tidygeocoder/blob/master/LICENSE.md)
[![CRAN](https://www.r-pkg.org/badges/version/tidygeocoder)](https://cran.r-project.org/package=tidygeocoder)
[![CRAN Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/tidygeocoder)](https://CRAN.R-project.org/package=tidygeocoder)
[![CRAN Downloads Per Month](http://cranlogs.r-pkg.org/badges/tidygeocoder)](https://cran.r-project.org/package=tidygeocoder)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R Build Status](https://github.com/jessecambon/tidygeocoder/workflows/R-CMD-check/badge.svg)](https://github.com/jessecambon/tidygeocoder/actions?workflow=R-CMD-check)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.5627341.svg)](https://doi.org/10.5281/zenodo.5627341)
<!-- badges: end -->

Tidygeocoder makes getting data from geocoding services easy. A unified high-level interface is provided for a selection of [supported geocoding services](https://jessecambon.github.io/tidygeocoder/articles/geocoder_services.html) and results are returned in [tibble](https://tibble.tidyverse.org/) (dataframe) format.

**Features:**

- Forward geocoding (addresses ⮕ coordinates)
- Reverse geocoding (coordinates ⮕ addresses)
- Batch geocoding (geocoding multiple addresses or coordinates in a single query) is automatically used if applicable.
- Duplicate, NA, and blank input data is handled elegantly; only unique inputs are submitted in queries, but the rows in the original data are preserved by default.
- The maximum rate of querying is automatically set according to the usage policies of the selected geocoding service.

In addition to the usage examples below, see the [Getting Started Vignette](https://jessecambon.github.io/tidygeocoder/articles/tidygeocoder.html) and [blog posts on tidygeocoder](https://jessecambon.github.io/tag/tidygeocoder).

## Installation

To install the stable version from CRAN (the official R package servers):

```{r, eval = FALSE}
install.packages('tidygeocoder')
```

Alternatively, you can install the latest development version from GitHub:

```{r, eval = FALSE}
devtools::install_github("jessecambon/tidygeocoder")
```

## Usage

In this first example we will geocode a few addresses using the `geocode()` function and plot them on a map with ggplot.

```{r}
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
```

The `geocode()` function geocodes addresses contained in a dataframe. The [Nominatim ("osm")](https://nominatim.org/) geocoding service is used here, but other services can be specified with the `method` argument. Only latitude and longitude are returned from the geocoding service in this example, but `full_results = TRUE` can be used to return all of the data from the geocoding service. See the `geo()` function documentation for details.

```{r, echo = FALSE}
knitr::kable(lat_longs)
```

Now that we have the longitude and latitude coordinates, we can use ggplot to plot our addresses on a map.

```{r usamap}
library(ggplot2)

ggplot(lat_longs, aes(longitude, latitude), color = "grey99") +
  borders("state") + geom_point() +
  ggrepel::geom_label_repel(aes(label = name)) +
  theme_void()
```

To perform reverse geocoding (obtaining addresses from geographic coordinates), we can use the `reverse_geocode()` function. The arguments are similar to the `geocode()` function, but now we specify the input data columns with the `lat` and `long` arguments. The  input dataset used here is the results of the geocoding query above. 

The single line address is returned in a column named by the `address` argument and all columns from the geocoding service results are returned because `full_results = TRUE`. See the `reverse_geo()` function documentation for more details.

<!-- 
Removing the licence column is done just to prevent a note from 
occurring in automated CRAN checks for an improper/old link.
-->
```{r}
reverse <- lat_longs %>%
  reverse_geocode(lat = latitude, long = longitude, method = 'osm',
                  address = address_found, full_results = TRUE) %>%
  select(-addr, -licence)
```

```{r, echo = FALSE}
knitr::kable(reverse)
```

## In the Wild

For inspiration, here are a few articles (with code) that leverage tidygeocoder:

- [Exercises: Spatial Data Wrangling with sf](http://www2.stat.duke.edu/courses/Spring21/sta323.001/exercises/lec_12.html) - part of a [statistical computing course](http://www2.stat.duke.edu/courses/Spring21/sta323.001/) at Duke
- [Geocoding the Minard Map](https://www.jla-data.net/eng/minard-map-tidygeocoder/) - recreating a famous infographic with geocoding
- [Mapping a network of women in demography](https://www.monicaalexander.com/posts/2021-21-02-mapping/) - using rvest and tidygeocoder to map Google Scholar data
- [Mapping Routes](https://bensstats.wordpress.com/2021/10/21/robservations-15-i-reverse-engineered-atlas-co-well-some-of-it/) - mapping routes with tidygeocoder and osrm
- [Road Routing in R](https://www.jla-data.net/eng/routing-in-r-context/) - demonstration of three different routing APIs
- [Mapping Texas Ports With R](https://www.sharpsightlabs.com/blog/mapping-texas-ports-with-r-part1/) - mapping the Texas coast with rnaturalearth and sf

## Contributing

Contributions to the tidygeocoder package are welcome. File [an issue](https://github.com/jessecambon/tidygeocoder/issues) for bug fixes or suggested features. If you would like to contribute code such as adding support for a new geocoding service, reference the [developer notes](https://jessecambon.github.io/tidygeocoder/articles/developer_notes.html) for instructions and documentation.

## Citing tidygeocoder

Use the `citation()` function:

``` r
citation('tidygeocoder')
```

</br>

<blockquote>
```{r, comment = '', echo = FALSE}
citation('tidygeocoder')
```
</blockquote>

Or refer to the [citation page](https://jessecambon.github.io/tidygeocoder/authors.html).
---
title: "Geocoding Services"
output: rmarkdown::html_vignette
description: >
  Documentation on the supported geocoding services and their query parameters and settings
vignette: >
  %\VignetteIndexEntry{Geocoding Services}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
set.seed(42)

library(DT)
library(tidygeocoder)
library(gt)
library(dplyr)
```

## Overview

The supported geocoding services are shown in the table below. The `method` is used to select the geocoding service in tidygeocoder functions such as `geo()` and `reverse_geo()`. The usage rate limitations are listed for the free tier of the service when applicable and many services have faster rates available with paid plans.

Also note that there are many other considerations when selecting a geocoding service such as if the service uses open source data with permissive licensing or if there are restrictions on how you can use the data provided by the service. Refer to each service's documentation for details.


```{r, echo = FALSE, message = FALSE, warning = FALSE, output = 'asis'}
library(dplyr)
check_mark <- "\U2705" #unicode character for heavy white check mark

geocoder_summary_table <-
  tidygeocoder::api_info_reference %>%
    mutate(
      service = paste0(
        '[', method_display_name, '](', site_url, ')'
      ),
      batch_geocoding = ifelse(method %in% names(tidygeocoder:::batch_func_map), check_mark, ''),
      api_key_required = ifelse(method %in% tidygeocoder::api_key_reference[['method']], check_mark, ''),
      api_documentation = paste0(
        '[docs](', api_documentation_url, ')'
      )
    ) %>%
    left_join(tidygeocoder::min_time_reference %>% select(method, description), by = 'method') %>%
    select(service, method, api_key_required, batch_geocoding, usage_limitations = description, api_documentation) %>%
    mutate(across(method, function(x) stringr::str_c('`', x, '`'))) %>% # format method column
    tidyr::replace_na(list(usage_limitations = ''))
  
# Format column names
colnames(geocoder_summary_table) <- colnames(geocoder_summary_table) %>%
  stringr::str_replace_all('_', ' ') %>%
  stringr::str_to_title() %>%
  stringr::str_replace_all('Api', 'API')

geocoder_summary_table %>%
  knitr::kable()
```

**Highlights:**

-   The US Census service does not support reverse geocoding and is limited to the United States.
-   Geocodio is limited to the United States and Canada.
-   The Census and Nominatim ("osm") services are free and do not require an API key.
-   ArcGIS can be used for free without an API key, but some features require a paid subscription and API key.
-   Geoapify, Geocodio, Location IQ, OpenCage, Mapbox, HERE, TomTom, MapQuest, and Bing are commercial services that offer free tiers.
-   The Google service has no free tier and [bills per query](https://developers.google.com/maps/documentation/geocoding/usage-and-billing).

**Other Notes:**

-   The US Census service supports street-level addresses only (ie. "11 Wall St New York, NY" is OK but "New York, NY" is not).
-   The Mapbox service is capable of performing batch geocoding when using the [permanent endpoint](https://docs.mapbox.com/api/search/geocoding/#batch-geocoding), but this capability is not currently implemented in tidygeocoder. If you'd like to add this capability to the package see [issue \#73](https://github.com/jessecambon/tidygeocoder/issues/73).
-   The ArcGIS service is capable of performing [batch geocoding](https://developers.arcgis.com/rest/geocode/api-reference/geocoding-geocode-addresses.htm) but this capability is not currently implemented in tidygeocoder. If you'd like to add this capability see [#102](https://github.com/jessecambon/tidygeocoder/issues/102).
-   For the ArcGIS service, an API Key is not strictly required if the service is used for search capabilities only (see [Free vs. paid operations](https://developers.arcgis.com/rest/geocode/api-reference/geocoding-free-vs-paid.htm)). It is possible to include an API Key on the request via the `custom_query` parameter:

``` r
tidygeocoder::geo(address = "New York, USA", method = "arcgis",
  custom_query = list(token = "<API_KEY>"))
```

## Usage Notes

- When used in batch mode, the **US Census** geocoder will return NA data when there are multiple results available for an address. The expectation is that you would see that a "Tie" is indicated and use single address geocoding to return the results for these addresses. See [#87](https://github.com/jessecambon/tidygeocoder/issues/87) for details.
- When performing reverse geocoding, the **Mapbox** service requires a `types` parameter to be set if `limit > 1`. See [#104](https://github.com/jessecambon/tidygeocoder/issues/104).
- The **Bing** batch geocoder does not use the `limit` parameter ([#106](https://github.com/jessecambon/tidygeocoder/issues/106)).

## API Parameters

The `api_parameter_reference` maps the API parameters for each geocoding service to a common set of "generic" parameters. The `generic_name` below is the generic parameter name while the `api_name` is the parameter name for the specified geocoding service (`method`). Refer to `?api_parameter_reference` for more details.

```{r, echo = FALSE}
api_parameter_reference %>% 
  mutate(across(c(method, generic_name, api_name), as.factor)) %>%
  datatable(filter = 'top', rownames = FALSE, 
  options = list(
    lengthMenu = c(5, 10, 15, 20, nrow(.)),
    pageLength = 10,
    autoWidth = TRUE)
  )
```

## API Key Retrieval

API keys are retrieved from environmental variables. The name of the environmental variable used for each service is stored in the `api_key_reference` dataset. See `?api_key_reference`.

```{r, echo = FALSE}
api_key_reference %>%
  gt() %>%
  opt_table_outline() %>%
  opt_table_lines() %>%
  tab_options(column_labels.font.weight = 'bold')
```

## Minimum Time Per Query

The minimum time (in seconds) required per query to comply with the usage limitations policies of each geocoding service is stored in the `min_time_reference` dataset. See `?min_time_reference`.

```{r, echo = FALSE}
min_time_reference %>%
  gt() %>%
  opt_table_outline() %>%
  opt_table_lines() %>%
  tab_options(column_labels.font.weight = 'bold')
```

<!-- https://bookdown.org/yihui/rmarkdown-cookbook/results-asis.html -->
Links to the usage policies for each geocoding service:

```{r, echo = FALSE, results = 'asis'}
cat(tidygeocoder:::get_api_usage_bullets(), sep = '\n')
```

## Batch Query Size Limits

The maximum number of inputs (geographic coordinates or addresses) per batch query for each geocoding service is stored in the `batch_limit_reference` dataset. See `?batch_limit_reference`.

```{r, echo = FALSE}
batch_limit_reference %>%
  gt() %>%
  fmt_number(columns = 'batch_limit', decimals = 0) %>%
  opt_table_outline() %>%
  opt_table_lines() %>%
  tab_options(column_labels.font.weight = 'bold')
```---
title: "Getting Started"
output: rmarkdown::html_vignette
description: >
  Start here if this is your first time using tidygeocder.
vignette: >
  %\VignetteIndexEntry{Getting Started}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---





## Introduction

Tidygeocoder provides a unified interface for performing both forward and reverse geocoding queries with a variety of geocoding services. In forward geocoding you provide an address to the geocoding service and you get latitude and longitude coordinates in return. In reverse geocoding you provide the latitude and longitude and the geocoding service will return that location's address. In both cases, other data about the location can be provided by the geocoding service.

The `geocode()` and `geo()` functions are for forward geocoding while the `reverse_geocode()` and `reverse_geo()` functions perform reverse geocoding. The `geocode()` and `reverse_geocode()` functions extract either addresses (forward geocoding) or coordinates (reverse geocoding) from the input dataframe and pass this data to the `geo()` and `reverse_geo()` functions respectively which execute the geocoding queries. All extra arguments (`...`) given to `geocode()` are passed to `geo()` and extra arguments given to `reverse_geocode()` are passed to `reverse_geo()`.

## Forward Geocoding


```r
library(tibble)
library(dplyr)
library(tidygeocoder)

address_single <- tibble(singlelineaddress = c(
  "11 Wall St, NY, NY",
  "600 Peachtree Street NE, Atlanta, Georgia"
))
address_components <- tribble(
  ~street, ~cty, ~st,
  "11 Wall St", "NY", "NY",
  "600 Peachtree Street NE", "Atlanta", "GA"
)
```

You can use the `address` argument to specify single-line addresses. Note that when multiple addresses are provided, the batch geocoding functionality of the Census geocoding service is used. Additionally, `verbose = TRUE` displays logs to the console.


```r
census_s1 <- address_single %>%
  geocode(address = singlelineaddress, method = "census", verbose = TRUE)
#> 
#> Number of Unique Addresses: 2
#> Executing batch geocoding...
#> Batch limit: 10,000
#> Passing 2 addresses to the US Census batch geocoder
#> Querying API URL: https://geocoding.geo.census.gov/geocoder/locations/addressbatch
#> Passing the following parameters to the API:
#> format : "json"
#> benchmark : "Public_AR_Current"
#> vintage : "Current_Current"
#> Query completed in: 1 seconds
```


|singlelineaddress                         |      lat|      long|
|:-----------------------------------------|--------:|---------:|
|11 Wall St, NY, NY                        | 40.70747| -74.01122|
|600 Peachtree Street NE, Atlanta, Georgia | 33.77085| -84.38505|

Alternatively you can run the same query with the `geo()` function by passing the address values from the dataframe directly. In either `geo()` or `geocode()`, the `lat` and `long` arguments are used to name the resulting latitude and longitude fields. Here the `method` argument is used to specify the "osm" (Nominatim) geocoding service. Refer to the `geo()` function documentation for the possible values of the `method` argument.


```r
osm_s1 <- geo(
  address = address_single$singlelineaddress, method = "osm",
  lat = latitude, long = longitude
)
#> Passing 2 addresses to the Nominatim single address geocoder
#> Query completed in: 2 seconds
```

|address                                   | latitude| longitude|
|:-----------------------------------------|--------:|---------:|
|11 Wall St, NY, NY                        | 40.70707| -74.01117|
|600 Peachtree Street NE, Atlanta, Georgia | 33.77086| -84.38614|

Instead of single-line addresses, you can use any combination of the following arguments to specify your addresses: `street`, `city`, `state`, `county`, `postalcode`, and `country`. 


```r
census_c1 <- address_components %>%
  geocode(street = street, city = cty, state = st, method = "census")
#> Passing 2 addresses to the US Census batch geocoder
#> Query completed in: 2.5 seconds
```

|street                  |cty     |st |      lat|      long|
|:-----------------------|:-------|:--|--------:|---------:|
|11 Wall St              |NY      |NY | 40.70747| -74.01122|
|600 Peachtree Street NE |Atlanta |GA | 33.77085| -84.38505|

To return the full geocoding service results (not just latitude and longitude), specify `full_results = TRUE`. Additionally, for the Census geocoder you can get fields for geographies such as Census tracts by specifying `api_options = list(census_return_type = 'geographies')`. Be sure to use `full_results = TRUE` with the "geographies" return type in order to allow the Census geography columns to be returned.


```r
census_full1 <- address_single %>% geocode(
  address = singlelineaddress,
  method = "census", full_results = TRUE, api_options = list(census_return_type = 'geographies')
)
#> Passing 2 addresses to the US Census batch geocoder
#> Query completed in: 1.2 seconds
```


|singlelineaddress                         |      lat|      long| id|input_address                                  |match_indicator |match_type |matched_address                      |tiger_line_id |tiger_side |state_fips |county_fips |census_tract |census_block |
|:-----------------------------------------|--------:|---------:|--:|:----------------------------------------------|:---------------|:----------|:------------------------------------|:-------------|:----------|:----------|:-----------|:------------|:------------|
|11 Wall St, NY, NY                        | 40.70747| -74.01122|  1|11 Wall St, NY, NY, , ,                        |Match           |Exact      |11 WALL ST, NEW YORK, NY, 10005      |59659656      |R          |36         |061         |000700       |1004         |
|600 Peachtree Street NE, Atlanta, Georgia | 33.77085| -84.38505|  2|600 Peachtree Street NE, Atlanta, Georgia, , , |Match           |Non_Exact  |600 PEACHTREE ST, ATLANTA, GA, 30308 |17343689      |L          |13         |121         |001902       |2003         |

As mentioned earlier, the `geocode()` function passes addresses in dataframes to the `geo()` function for geocoding so we can also directly use the `geo()` function in a similar way:

<!-- 
Removing the licence column is done just to prevent a note from 
occurring in automated CRAN checks for an improper/old link.
-->

```r
salz <- geo("Salzburg, Austria", method = "osm", full_results = TRUE) %>%
  select(-licence)
#> Passing 1 address to the Nominatim single address geocoder
#> Query completed in: 1 seconds
```


|address           |      lat|     long| place_id|osm_type |   osm_id|boundingbox                                    |display_name               |class |type | importance|icon                                                                     |
|:-----------------|--------:|--------:|--------:|:--------|--------:|:----------------------------------------------|:--------------------------|:-----|:----|----------:|:------------------------------------------------------------------------|
|Salzburg, Austria | 47.79813| 13.04648|   147539|node     | 34964314|47.6381346, 47.9581346, 12.8864806, 13.2064806 |Salzburg, 5020, Österreich |place |city |  0.6854709|https://nominatim.openstreetmap.org/ui/mapicons//poi_place_city.p.20.png |

## Reverse Geocoding

For reverse geocoding you'll use `reverse_geocode()` instead of `geocode()` and `reverse_geo()` instead of `geo()`. Note that the reverse geocoding functions are structured very similarly to the forward geocoding functions and share many of the same arguments (`method`, `limit`, `full_results`, etc.). For reverse geocoding you will provide latitude and longitude coordinates as inputs and the location's address will be returned by the geocoding service. 

Below, the `reverse_geocode()` function is used to geocode coordinates in a dataframe. The `lat` and `long` arguments specify the columns that contain the latitude and longitude data. The `address` argument can be used to specify the single line address column name that is returned from the geocoder. Just as with forward geocoding, the `method` argument is used to specify the geocoding service.


```r
lat_longs1 <- tibble(
  latitude = c(38.895865, 43.6534817),
  longitude = c(-77.0307713, -79.3839347)
)

rev1 <- lat_longs1 %>%
  reverse_geocode(lat = latitude, long = longitude, address = addr, method = "osm")
#> Passing 2 coordinates to the Nominatim single coordinate geocoder
#> Query completed in: 2 seconds
```


| latitude| longitude|addr                                                                                                                                               |
|--------:|---------:|:--------------------------------------------------------------------------------------------------------------------------------------------------|
| 38.89587| -77.03077|Freedom Plaza, 1455, Pennsylvania Avenue Northwest, Washington, District of Columbia, 20004, United States                                         |
| 43.65348| -79.38393|Toronto City Hall, 100, Queen Street West, Financial District, Spadina—Fort York, Old Toronto, Toronto, Golden Horseshoe, Ontario, M5H 2N2, Canada |

The same query can also be performed by passing the latitude and longitudes directly to the `reverse_geo()` function. Here we will use `full_results = TRUE` so that the full results are returned (not just the single line address column).


```r
rev2 <- reverse_geo(
  lat = lat_longs1$latitude,
  long = lat_longs1$longitude,
  method = "osm",
  full_results = TRUE
)
#> Passing 2 coordinates to the Nominatim single coordinate geocoder
#> Query completed in: 2 seconds
glimpse(rev2)
#> Rows: 2
#> Columns: 22
#> $ lat            <dbl> 38.89587, 43.65348
#> $ long           <dbl> -77.03077, -79.38393
#> $ address        <chr> "Freedom Plaza, 1455, Pennsylvania Avenue Northwest, Washington, District of Columbia, 20004, United States", "Toronto City Hall, 1…
#> $ place_id       <int> 284009208, 148364261
#> $ licence        <chr> "Data © OpenStreetMap contributors, ODbL 1.0. https://osm.org/copyright", "Data © OpenStreetMap contributors, ODbL 1.0. https://osm…
#> $ osm_type       <chr> "relation", "way"
#> $ osm_id         <int> 8060882, 198500761
#> $ osm_lat        <chr> "38.895849999999996", "43.6536032"
#> $ osm_lon        <chr> "-77.03077367444483", "-79.38400546703345"
#> $ tourism        <chr> "Freedom Plaza", NA
#> $ house_number   <chr> "1455", "100"
#> $ road           <chr> "Pennsylvania Avenue Northwest", "Queen Street West"
#> $ city           <chr> "Washington", "Old Toronto"
#> $ state          <chr> "District of Columbia", "Ontario"
#> $ postcode       <chr> "20004", "M5H 2N2"
#> $ country        <chr> "United States", "Canada"
#> $ country_code   <chr> "us", "ca"
#> $ boundingbox    <list> <"38.8956276", "38.896068", "-77.03182", "-77.0297273">, <"43.6529946", "43.6541458", "-79.3848438", "-79.3830415">
#> $ amenity        <chr> NA, "Toronto City Hall"
#> $ neighbourhood  <chr> NA, "Financial District"
#> $ quarter        <chr> NA, "Spadina—Fort York"
#> $ state_district <chr> NA, "Golden Horseshoe"
```

## Working With Messy Data

Only unique input data (either addresses or coordinates) is passed to geocoding services even if your data contains duplicates. NA and blank inputs are excluded from queries. Input latitudes and longitudes are also limited to the range of possible values. 

Below is an example of how duplicate and missing data is handled by tidygeocoder. As the console messages shows, only the two unique addresses are passed to the geocoding service.


```r
# create a dataset with duplicate and NA addresses
duplicate_addrs <- address_single %>%
  bind_rows(address_single) %>%
  bind_rows(tibble(singlelineaddress = rep(NA, 3)))

duplicates_geocoded <- duplicate_addrs %>%
  geocode(singlelineaddress, verbose = TRUE)
#> 
#> Number of Unique Addresses: 2
#> Passing 2 addresses to the Nominatim single address geocoder
#> 
#> Number of Unique Addresses: 1
#> Querying API URL: https://nominatim.openstreetmap.org/search
#> Passing the following parameters to the API:
#> q : "11 Wall St, NY, NY"
#> limit : "1"
#> format : "json"
#> HTTP Status Code: 200
#> Query completed in: 0.2 seconds
#> Total query time (including sleep): 1 seconds
#> 
#> 
#> Number of Unique Addresses: 1
#> Querying API URL: https://nominatim.openstreetmap.org/search
#> Passing the following parameters to the API:
#> q : "600 Peachtree Street NE, Atlanta, Georgia"
#> limit : "1"
#> format : "json"
#> HTTP Status Code: 200
#> Query completed in: 0.2 seconds
#> Total query time (including sleep): 1 seconds
#> 
#> Query completed in: 2 seconds
```


|singlelineaddress                         |      lat|      long|
|:-----------------------------------------|--------:|---------:|
|11 Wall St, NY, NY                        | 40.70707| -74.01117|
|600 Peachtree Street NE, Atlanta, Georgia | 33.77086| -84.38614|
|11 Wall St, NY, NY                        | 40.70707| -74.01117|
|600 Peachtree Street NE, Atlanta, Georgia | 33.77086| -84.38614|
|NA                                        |       NA|        NA|
|NA                                        |       NA|        NA|
|NA                                        |       NA|        NA|

As shown above, duplicates will not be removed from your results by default. However, you can return only unique results by using `unique_only = TRUE`. Note that passing `unique_only = TRUE` to `geocode()` or `reverse_geocode()` will result in the original dataframe format (including column names) to be discarded in favor of the standard field names (ie. "address", 'lat, 'long', etc.).


```r
uniqueonly1 <- duplicate_addrs %>%
  geocode(singlelineaddress, unique_only = TRUE)
#> Passing 2 addresses to the Nominatim single address geocoder
#> Query completed in: 2 seconds
```

|address                                   |      lat|      long|
|:-----------------------------------------|--------:|---------:|
|11 Wall St, NY, NY                        | 40.70707| -74.01117|
|600 Peachtree Street NE, Atlanta, Georgia | 33.77086| -84.38614|


## Combining Multiple Queries

The `geocode_combine()` function allows you to execute and combine the results of multiple geocoding queries. The queries are specified as a list of lists with the `queries` parameter and are executed in the order provided. By default only addresses that are not found are passed to the next query, but this behavior can be toggled with the `cascade` argument. 

In the first example below, the US Census service is used for the first query while the Nominatim ("osm") service is used for the second query. The `global_params` argument passes the `address` column from the input dataset to both queries.


```r
addresses_combine <- tibble(
  address = c('100 Wall Street NY, NY', 'Paris', 'Not An Address')
)

cascade_results1 <- addresses_combine %>%
  geocode_combine(
    queries = list(
      list(method = 'census'),
      list(method = 'osm')
    ),
    global_params = list(address = 'address')
  )
#> 
#> Passing 3 addresses to the US Census batch geocoder
#> Query completed in: 0.2 seconds
#> Passing 2 addresses to the Nominatim single address geocoder
#> Query completed in: 2 seconds
```


|address                |      lat|       long|query  |
|:----------------------|--------:|----------:|:------|
|100 Wall Street NY, NY | 40.70516| -74.007350|census |
|Paris                  | 48.85889|   2.320041|osm    |
|Not An Address         |       NA|         NA|       |

If `cascade` is set to FALSE then all addresses are attempted by each query regardless of if the address was found initially or not.


```r
no_cascade_results1 <- addresses_combine %>%
  geocode_combine(
    queries = list(
      list(method = 'census'),
      list(method = 'osm')
    ),
    global_params = list(address = 'address'),
    cascade = FALSE
  )
#> 
#> Passing 3 addresses to the US Census batch geocoder
#> Query completed in: 0.3 seconds
#> Passing 3 addresses to the Nominatim single address geocoder
#> Query completed in: 3 seconds
```


|address                |      lat|       long|query  |
|:----------------------|--------:|----------:|:------|
|100 Wall Street NY, NY | 40.70516| -74.007350|census |
|100 Wall Street NY, NY | 40.70522| -74.006800|osm    |
|Paris                  |       NA|         NA|census |
|Paris                  | 48.85889|   2.320041|osm    |
|Not An Address         |       NA|         NA|census |
|Not An Address         |       NA|         NA|osm    |

Additionally, the results from each query can be returned in separate list items by setting `return_list = TRUE`.

## Customizing Queries

The `limit` argument can be specified to allow multiple results (rows) per input if available. The maximum value for the `limit` argument is often 100 for geocoding services. To use the default `limit` value for the selected geocoding service you can use `limit = NULL` which will prevent the limit parameter from being included in the query.


```r
geo_limit <- geo(
  c("Lima, Peru", "Cairo, Egypt"),
  method = "osm",
  limit = 3, full_results = TRUE
)
#> Passing 2 addresses to the Nominatim single address geocoder
#> Query completed in: 2 seconds
glimpse(geo_limit)
#> Rows: 6
#> Columns: 13
#> $ address      <chr> "Lima, Peru", "Lima, Peru", "Lima, Peru", "Cairo, Egypt", "Cairo, Egypt", "Cairo, Egypt"
#> $ lat          <dbl> -11.96784, -12.03089, -11.95785, 30.04439, 30.08695, 30.22503
#> $ long         <dbl> -77.01094, -77.09072, -77.04139, 31.23573, 31.96162, 31.69733
#> $ place_id     <int> 128601425, 118029294, 128061402, 283020077, 216381765, 272362266
#> $ licence      <chr> "Data © OpenStreetMap contributors, ODbL 1.0. https://osm.org/copyright", "Data © OpenStreetMap contributors, ODbL 1.0. https://osm.o…
#> $ osm_type     <chr> "way", "way", "way", "relation", "way", "way"
#> $ osm_id       <int> 116948976, 71187508, 115715296, 5466227, 544272017, 914722135
#> $ boundingbox  <list> <"-11.9678367", "-11.9672995", "-77.0117952", "-77.0109387">, <"-12.0308994", "-12.0304024", "-77.0911106", "-77.090276">, <"-11.9608…
#> $ display_name <chr> "Peru, San Juan de Lurigancho, Huascar, Lima, Lima Metropolitana, Lima, 15423, Perú", "Peru, San Martín de Porres, Lima, Lima Metropo…
#> $ class        <chr> "highway", "highway", "highway", "place", "highway", "highway"
#> $ type         <chr> "residential", "residential", "residential", "city", "motorway", "trunk"
#> $ importance   <dbl> 0.3200000, 0.3200000, 0.3200000, 0.6960286, 0.1000000, 0.1000000
#> $ icon         <chr> NA, NA, NA, "https://nominatim.openstreetmap.org/ui/mapicons//poi_place_city.p.20.png", NA, NA
```

To directly specify specific API parameters for a given `method` you can use the `custom_query` parameter. For example, [the Nominatim (OSM) geocoder has a 'polygon_geojson' argument](https://nominatim.org/release-docs/develop/api/Details/#parameters) that can be used to return GeoJSON geometry content. To pass this parameter you can insert it with a named list using the `custom_query` argument:


```r
cairo_geo <- geo("Cairo, Egypt",
  method = "osm", full_results = TRUE,
  custom_query = list(polygon_geojson = 1), verbose = TRUE
)
#> 
#> Number of Unique Addresses: 1
#> Passing 1 address to the Nominatim single address geocoder
#> 
#> Number of Unique Addresses: 1
#> Querying API URL: https://nominatim.openstreetmap.org/search
#> Passing the following parameters to the API:
#> q : "Cairo, Egypt"
#> limit : "1"
#> polygon_geojson : "1"
#> format : "json"
#> HTTP Status Code: 200
#> Query completed in: 0.2 seconds
#> Total query time (including sleep): 1 seconds
#> 
#> Query completed in: 1 seconds
glimpse(cairo_geo)
#> Rows: 1
#> Columns: 15
#> $ address             <chr> "Cairo, Egypt"
#> $ lat                 <dbl> 30.04439
#> $ long                <dbl> 31.23573
#> $ place_id            <int> 283020077
#> $ licence             <chr> "Data © OpenStreetMap contributors, ODbL 1.0. https://osm.org/copyright"
#> $ osm_type            <chr> "relation"
#> $ osm_id              <int> 5466227
#> $ boundingbox         <list> <"29.7483062", "30.3209168", "31.2200331", "31.9090054">
#> $ display_name        <chr> "القاهرة, محافظة القاهرة, مصر"
#> $ class               <chr> "place"
#> $ type                <chr> "city"
#> $ importance          <dbl> 0.6960286
#> $ icon                <chr> "https://nominatim.openstreetmap.org/ui/mapicons//poi_place_city.p.20.png"
#> $ geojson.type        <chr> "Polygon"
#> $ geojson.coordinates <list> <<array[1 x 119 x 2]>>
```

To test a query without sending any data to a geocoding service, you can use `no_query = TRUE` (NA results are returned).


```r
noquery1 <- geo(c("Vancouver, Canada", "Las Vegas, NV"),
  no_query = TRUE,
  method = "arcgis"
)
#> 
#> Number of Unique Addresses: 2
#> Passing 2 addresses to the ArcGIS single address geocoder
#> 
#> Number of Unique Addresses: 1
#> Querying API URL: https://geocode.arcgis.com/arcgis/rest/services/World/GeocodeServer/findAddressCandidates
#> Passing the following parameters to the API:
#> SingleLine : "Vancouver, Canada"
#> maxLocations : "1"
#> f : "json"
#> 
#> Number of Unique Addresses: 1
#> Querying API URL: https://geocode.arcgis.com/arcgis/rest/services/World/GeocodeServer/findAddressCandidates
#> Passing the following parameters to the API:
#> SingleLine : "Las Vegas, NV"
#> maxLocations : "1"
#> f : "json"
#> Query completed in: 0 seconds
```

|address           | lat| long|
|:-----------------|---:|----:|
|Vancouver, Canada |  NA|   NA|
|Las Vegas, NV     |  NA|   NA|

Additional usage notes for the `geocode()`, `geo()`, `reverse_geocode()`, and `reverse_geo()` functions:

- Set `quiet = TRUE` to silence console logs displayed by default (how many inputs were submitted, to what geocoding service, and the elapsed time). 
- Use the `progress_bar` argument to control if a progress bar is displayed.
- The `verbose`, `quiet`, and `progress_bar` arguments can be set globally with `options`. For instance `options(tidygeocoder.verbose = TRUE)` will set verbose to `TRUE` for all queries by default.
- To customize the API endpoint, use the `api_options` or `api_url` arguments. See `?geo` or `?reverse_geo` for details.
- If not specified manually, the `min_time` argument will default to a value based on the maximum query rate of the given geocoding service. If you are using a local Nominatim server or have a commercial geocoder plan that has less restrictive usage limits then you can manually set `min_time` to a lower value (such as 0).
- The `mode` argument can be used to specify whether the batch geocoding or single address/coordinate geocoding should be used. By default batch geocoding will be used if available when more than one address or coordinate is provided (with some noted exceptions for slower batch geocoding services).
- The `return_addresses` and `return_coords` parameters (for forward and reverse geocoding respectively) can be used to toggle whether the input addresses or coordinates are returned. Setting these parameters to `FALSE` is necessary to use batch geocoding if `limit` is greater than 1 or NULL. 
- For the `reverse_geocode()` and `geocode()` functions, the `return_input` argument can be used to toggle if the input dataset is included in the returned dataframe.
- For use with the [memoise package](https://memoise.r-lib.org/), you may need to use character values when specifying input data columns in the `geocode()` and `reverse_geocode()` functions. See [#154](https://github.com/jessecambon/tidygeocoder/pull/154) for details.
---
title: "Developer Notes"
output: rmarkdown::html_vignette
description: >
  Documentation for developers including how to add support for an additional geocoding service
vignette: >
  %\VignetteIndexEntry{Developer Notes}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tibble.print_min = 5, tibble.print_max = 5)

library(dplyr)
library(tidygeocoder)
```

This page contains documentation relevant for those wishing to contribute to the package and specific instructions for how to add support for a new geocoding service.

## Introduction

The two core functions to focus on in the package are [geo()](https://github.com/jessecambon/tidygeocoder/blob/main/R/geo.R) and [reverse_geo()](https://github.com/jessecambon/tidygeocoder/blob/main/R/reverse_geo.R). These functions have very similar layouts, but `geo()` is for forward geocoding while `reverse_geo()` is for reverse geocoding. The `geocode()` and `reverse_geocode()` functions only extract input data from a dataframe and pass it to the `geo()` and `reverse_geo()` functions respectively for geocoding.

Both the `geo()` and `reverse_geo()` functions take inputs (either addresses or coordinates) and call other functions as needed to deduplicate the inputs, pause to comply with API usage rate policies, and execute queries. Key parameters and settings for geocoding are stored for easy access and display in built-in datasets. 

Consider this query:

```{r}
library(dplyr)
library(tidygeocoder)

df <- tibble(
  id = c(1, 2, 1),
  locations = c('tokyo', 'madrid', 'tokyo')
  )

df %>%
  geocode(address = locations, method = 'osm', full_results = TRUE, verbose = TRUE)
```

Here is what is going on behind the scenes:

- The `geocode()` function extracts the address data from the input dataframe and passes it to the `geo()` function.
- The `geo()` function looks for unique inputs and prepares them for geocoding. In this case, there is one duplicate input so we only have two unique inputs.
- The `geo()` function must figure out whether to use *single address geocoding* (1 address per query) or *batch geocoding* (multiple addresses per query). In this case the specified Nominatim ("osm") geocoding service does not have a batch geocoding function so single address geocoding is used.
- Because single address geocoding is used, the `geo()` function is called once for each input to geocode all addresses (twice in this case) and the results are combined. If batch geocoding was used then the appropriate batch geocoding function would be called based on the geocoding service specified.
- Because the specified geocoding service has a usage limit, the rate of querying is limited accordingly. By default this is based on the `min_time_reference` dataset. This behavior can be modified with the `min_time` argument.
- Since the input data was deduplicated, the results must be aligned to the original inputs (which contained duplicates) so that the original data structure is preserved. Alternatively, if you only want to return unique results, you can specify `unique_only = TRUE`.
- This combined data is returned by `geo()` to the `geocode()` function. The `geocode()` function then combines the returned data with the original dataset. 

Refer to the notes below on adding a geocoding service for more specific documentation on the code structure.

## Adding a New Geocoding Service

This section documents how to add support for a new geocoding service to the package. Required changes are organized by file. If anything isn't clear, feel free to [file an issue](https://github.com/jessecambon/tidygeocoder/issues).

**Base all changes on the main branch**.

### Files to Update

* **[R/api_url.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/api_url.R)**
    * Add a standalone function for obtaining the API URL and update the `get_api_url()` function accordingly. If arguments need to be added to the `get_api_url()` function, make sure to adjust the calls to this function in the `geo()` and `reverse_geo()` functions accordingly.
* **[data-raw/api_parameter_reference.R](https://github.com/jessecambon/tidygeocoder/blob/main/data-raw/api_parameter_reference.R)**
    * Add rows to the [api_parameter_reference](https://jessecambon.github.io/tidygeocoder/articles/tidygeocoder.html#api-reference) dataset to include the geocoding service. Each service is referred to by a short name in the `method` column (which is how the service is specified in the `geo()` and `geocode()` functions). The `generic_name` column has the universal parameter name that is used across geocoding services (ie. "address", "limit", etc.) while the `api_name` column stores the parameter names that are specific to the geocoding service. 
    * Note that there is no need to include parameters that are only used for reverse geocoding or parameters that have no equivalent in other geocoding services (ie. there is no `generic_name`) **unless** the parameters are required. Parameters can always be passed to services directly with the `custom_query` argument in `geo()` or `reverse_geo()`. 
* **[data-raw/api_references.R](https://github.com/jessecambon/tidygeocoder/blob/main/data-raw/api_references.R)**
   * Add a row to `min_time_reference` with the minimum time each query should take (in seconds) according to the geocoding service's free tier usage restrictions.
   * Add a row to `api_key_reference` if the service requires an API key.
   * If the service you are adding has batch geocoding capabilities, add the maximum batch size (as a row) to `batch_limit_reference`.
   * Add a row to `api_info_reference` with links to the service's website, documentation, and usage policy.
* **[R/geo.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/geo.R)**
    * If the service supports batch geocoding then add a new function in **[R/batch_geocoding.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/batch_geocoding.R)** and add it to the `batch_func_map` named list.
* **[R/reverse_geo.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/reverse_geo.R)**
    * Update the `get_coord_parameters()` function based on how the service passed latitude and longitude coordinates for reverse geocoding.
    * If the service supports reverse batch geocoding then add a new function in **[R/reverse_batch_geocoding.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/reverse_batch_geocoding.R)** and add it to the `reverse_batch_func_map` named list.
* **[R/results_processing.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/results_processing.R)**
    * Update the `extract_results()` function which is used for parsing single addresses (ie. not batch geocoding). You can see examples of how I've tested out parsing the results of geocoding services [here](https://github.com/jessecambon/tidygeocoder/tree/main/sandbox/query_debugging).
    * In a similar fashion, update the `extract_reverse_results()` function for reverse geocoding.
    * Update the `extract_errors_from_results()` function to extract error messages for invalid queries.
* If applicable, add new tests to the scripts in the [tests directory](https://github.com/jessecambon/tidygeocoder/tree/main/tests/testthat) for the method. Note that tests should avoid making a HTTP query (ie. use `no_query = TRUE` in the `geo()` and `geocode()` functions).
* **[R/global_variables.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/global_variables.R)**
    * If applicable, add your service to one of the global variables.

### Other Files

These files don't necessarily need to be updated. However, you might need to make changes to these files if the service you are implementing requires some non-standard workarounds. 

* **[R/query_factory.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/query_factory.R)**
    * Houses the functions used to create and execute API queries.
* **[R/documentation.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/documentation.R)**
    * Functions for producing rmarkdown package documentation.
* **[R/data.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/data.R)**
    * Documentation for in-built datasets.
* **[R/utils.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/utils.R)**
    * Common utility functions.
* **[R/input_handling.R](https://github.com/jessecambon/tidygeocoder/blob/main/R/input_handling.R)**
    * Handles the deduplication of input data.
* **[external/create_logo.Rmd](https://github.com/jessecambon/tidygeocoder/blob/main/external/create_logo.Rmd)**
    * Creates the package logo.
* **[vignettes/tidygeocoder.Rmd.orig](vignettes/tidygeocoder.Rmd.orig)**
    * This file produces the vignette. See the knit command and comments at the top of this file.

### Testing

- Test out the new service to make sure it behaves as expected. You can reference tests and example code in [the 'sandbox' folder](https://github.com/jessecambon/tidygeocoder/tree/main/sandbox).
- Run **`devtools::check()`** to make sure the package still passes all tests and checks, but note that these tests are designed to work offline so they do not make queries to geocoding services.
- As a final check, run **[external/online_tests.R](https://github.com/jessecambon/tidygeocoder/blob/main/external/online_tests.R)** to test making queries to the geocoding services. These tests are not included in the internal package tests (`devtools::test()`) because they require API keys which would not exist on all systems and are dependent on the geocoding services being online at that the time of the test.
- Run the commands detailed in **[cran-comments.md](https://github.com/jessecambon/tidygeocoder/blob/main/cran-comments.md)** to test the package on other environments. Note that these tests should also be included in the automated GitHub actions tests for pull requests.

### Releasing a New Version

To release a new version of tidygeocoder:

- Run the tests detailed above
- Update the package version in [DESCRIPTION](DESCRIPTION)
- Update the package version for the citation note in [inst/CITATION](inst/CITATION)
- Reserve a DOI for a new package version with [Zenodo](https://zenodo.org/)
- Update the Zenodo DOI in [README.Rmd](README.Rmd) and reknit to update [README.md](README.md)
- Update the site with `pkgdown::build_site()`
- Use `urlchecker::url_check()` to check all package URLs

Lastly, run `devtools::release()` to release the new version---
title: "Regression Testing"
---

Regression testing. Make sure to cover services with API keys since those methods aren't covered in the vignette

```{r setup}
library(tibble)
library(dplyr)
library(knitr)
library(tidygeocoder)

sample1 <- tibble(address = c('11 Wall St New York, NY', NA, '',
    '1600 Pennsylvania Ave NW Washington, DC', '11 Wall St New York, NY', 
    'Toronto, Canada'))

sample2 <- tibble(street = c('1600 Pennsylvania Ave NW', '11 Wall Street', ''), 
  city = c('Washington', 'New York', 'Nashville'), state = c('DC', 'NY', 'TN'))
```



test batch limit

```{r, error = TRUE}
sample1 %>%
  geocode(address, method = 'census', verbose = FALSE, batch_limit = 1)
```


```{r}
louisville_check <- louisville %>% 
  geocode(street = street, city = city, state = state, postalcode = zip, verbose = TRUE, full_results = TRUE) %>%
  mutate(lat_diff = lat - latitude, long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_check
```

```{r}
louisville_geocodio <- louisville %>% 
  geocode(street = street, city = city, state = state, postalcode = zip, verbose = TRUE, method = 'geocodio', full_results = TRUE) %>%
   mutate(lat_diff = lat - latitude, long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_geocodio
```



```{r}
geo(street='1449 ST JAMES CT', city =  'Louisville', state = 'KY', postalcode =  '40208', full_results = TRUE, verbose = T)	
```


```{r}
geocode_google1 <- sample1 %>% 
  geocode(address = address, method = 'google', full_results = TRUE, verbose = TRUE)
geocode_google1
```

```{r}
geocode_google1_notflat <- sample1 %>% slice(1,2,4,5) %>%
  geocode(address = address, method = 'google', full_results = TRUE, verbose = TRUE, flatten = FALSE, unique_only = TRUE)
geocode_google1_notflat
```


```{r}
cascade_r1 <- geo(c('', NA, sample_addresses$addr, NA, ''), method = 'cascade', cascade_order = c('census', 'geocodio'), verbose = T)
cascade_r1
```


```{r}
gc1 <- geocode_geocodio1 <- sample1 %>% 
  geocode(address = address, method = 'geocodio', verbose = TRUE)
gc1
```

```{r}
iq1 <- geocode_geocodio1 <- sample1 %>% 
  geocode(address = address, method = 'iq', full_results = TRUE, verbose = TRUE, unique_only = TRUE)
iq1
```



```{r}
gc_c1 <- geo(method = 'geocodio', street = c('1600 Pennsylvania Ave NW', '11 Wall Street', ''), 
  city = c('Washington', 'New York', 'Nashville'), state = c('DC', 'NY', 'TN'), verbose = TRUE, full_results = TRUE)
gc_c1
```

```{r}
iq_c1 <- sample2 %>% 
  geocode(street = street, city = city, state = state, method = 'iq', verbose = TRUE, full_results = TRUE)
iq_c1
```


```{r, error = TRUE}
geo('blank address', method = 'makeanerror')
```


Check geocodio error catching

```{r, error = TRUE}
geo(' -----', method ='geocodio')
```

# Check mapbox

```{r}

louisville_mapbox <- louisville %>%
  slice(1:10) %>%
  geocode(
    address = street,
    method = 'mapbox',
    full_results = TRUE,
    custom_query = list(
      proximity = paste0(louisville[1, c("longitude", "latitude")], 
                         collapse = ","))) %>%
  mutate(lat_diff = lat - latitude,
         long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_mapbox
```


```{r}
geocode_mapbox1 <- sample1 %>%
  mutate(addr = address) %>%
  select(-address) %>%
  geocode(
    address = addr,
    method = 'mapbox',
    full_results = TRUE,
    verbose = TRUE
  )
geocode_mapbox1
```


```{r}

geocode_mapbox2 <- sample1 %>%
  geocode(
    address = address,
    method = 'mapbox',
    full_results = FALSE,
    verbose = TRUE
  )
geocode_mapbox2
```

## Check mapbox error catching

```{r}
geo('Testing',
    method = 'mapbox',
    custom_query = list(country = "ERROR"))
```
# Check HERE

```{r}
# single geocoding

louisville_here <- louisville %>%
  geocode(
    address = street,
    method = 'here',
    full_results = TRUE,
    custom_query = list(at = paste0(louisville[1, c("latitude", "longitude")],
                                    collapse = ","))
  ) %>%
  mutate(lat_diff = lat - latitude,
         long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_here
```

```{r}
# Force batch geocoding

louisville_here2 <- louisville %>%
  mutate(addr = paste(street, city, state, "USA", dlm=", ")) %>%
  select(addr, 
         longitude, latitude) %>%
  geocode(
    address = addr,
    method = 'here',
    mode = 'batch',
    verbose = TRUE,
    full_results = TRUE) %>%
  mutate(lat_diff = lat - latitude,
         long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_here2
```


```{r}
geocode_here1 <- sample1 %>%
  geocode(
    address = address,
    method = 'here',
    full_results = TRUE
  )
geocode_here1
```

# Check tomtom

```{r}

louisville_tomtom <- louisville %>%
  geocode(
    address = street,
    method = 'tomtom',
    full_results = TRUE,
    custom_query = list(
      lat = louisville[1, c("latitude")],
      lon = louisville[1, c("longitude")]
      )) %>%
  mutate(lat_diff = lat - latitude,
         long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_tomtom
```


```{r}
geocode_tomtom1 <- sample1 %>%
  mutate(addr = address) %>%
  select(-address) %>%
  geocode(
    address = addr,
    method = 'tomtom',
    mode = 'single',
    full_results = TRUE,
    return_input = FALSE,
    verbose = TRUE,
    limit = 3
  )
geocode_tomtom1
```


```{r}
geocode_tomtom2 <- sample1 %>%
  geocode(
    address = address,
    method = 'tomtom',
    full_results = FALSE,
    verbose = TRUE
  )
geocode_tomtom2
```

# Check mapquest

```{r}
louisville_mapquest <- louisville %>%
  mutate(st = street) %>%
  select(st, latitude, longitude) %>%
  geocode(
    address = st,
    method = 'mapquest',
    full_results = TRUE,
    custom_query = list(
      lat = louisville[1, c("latitude")],
      lon = louisville[1, c("longitude")]
      )) %>%
  mutate(lat_diff = lat - latitude,
         long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_mapquest
```


```{r}
geocode_mapquest1 <- sample1 %>%
  mutate(addr = address) %>%
  select(-address) %>%
  geocode(
    address = addr,
    method = 'mapquest',
    mode = 'single',
    return_input = FALSE,
    full_results = TRUE,
    verbose = TRUE,
    limit = 3
  )
geocode_mapquest1
```


```{r}
geocode_mapquest2 <- sample1 %>%
  geocode(
    address = address,
    method = 'mapquest',
    full_results = FALSE,
    verbose = TRUE
  )
geocode_mapquest2
```


# Check bing

```{r}

louisville_bing <- louisville %>%
  slice(1:20) %>%
  geocode(
    address = street,
    method = 'bing',
    full_results = TRUE) %>%
  mutate(lat_diff = lat - latitude,
         long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_bing
```


```{r}
geocode_bing1 <- sample1 %>%
  mutate(addr = address) %>%
  select(-address) %>%
  geocode(
    address = addr,
    method = 'bing',
    mode = 'single',
    return_input = FALSE,
    full_results = TRUE,
    verbose = TRUE,
    limit = 3
  )
geocode_bing1
```


```{r}

geocode_bing2 <- sample1 %>%
  geocode(
    address = address,
    method = 'bing',
    full_results = FALSE,
    verbose = TRUE
  )
geocode_bing2
```

# Check arcgis

```{r}
louisville_arcgis <- louisville %>%
  sample_n(10) %>%
  mutate(full = paste0(street, ", Louisville")) %>%
  geocode(
    address = full,
    method = 'arcgis',
    full_results = TRUE,
    verbose = TRUE,
    custom_query = list(
      outFields = "*"
      )) %>%
  mutate(lat_diff = lat - latitude,
         long_diff = long - longitude) %>%
  select(lat_diff, long_diff, everything())

louisville_arcgis
```


```{r}
geocode_arcgis1 <- sample1 %>%
  mutate(addr = address) %>%
  select(-address) %>%
  geocode(
    address = addr,
    method = 'arcgis',
    mode = 'single',
    full_results = TRUE,
    verbose = TRUE,
    limit = 1
  )
geocode_arcgis1
```---
title: "Operation Logo"
---

Create the package hex sticker

https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
https://stackoverflow.com/questions/43207947/whole-earth-polygon-for-world-map-in-ggplot2-and-sf
https://github.com/GuangchuangYu/hexSticker


```{r setup}
library(tidyverse)
library(showtext)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(hexSticker)
library(here)
```

```{r}
world_coastlines <- ne_coastline(scale = 'medium', returnclass = 'sf')

#crs <- "+proj=robin +ellps=WGS84 +lat_0=20 +lon_0=30"
#crs <- "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

sphere <- st_graticule(ndiscr = 10000, margin = 10e-6) %>%
  st_transform(crs = 3035) %>%
  st_convex_hull() %>%
  summarise(geometry = st_union(geometry))
```


Make the globeplot

```{r}
globe_background <- ggplot()  +
  geom_sf(data = sphere, color = 'black', size = 0.23, fill = 'white') +
  geom_sf(data = world_coastlines, color = "black", size = 0.15) +
  theme_void() + theme_transparent()
globe_background + theme(panel.background = element_rect(fill = 'black'))
```



## Create Hex

Load google font http://www.google.com/fonts

```{r}
# Load google font http://www.google.com/fonts

# candidate fonts: Crete Round; Encode Sans Semi Expanded

font_name <- "Encode Sans Semi Expanded"
font_add_google(font_name)

# Save and plot hex 

s <- sticker(globe_background,
    package = "tidygeocoder", p_size = 18, s_x = 1, s_y = .7, s_width = 1.1, s_height = 1.1, p_family = font_name,
    h_color = 'dimgrey', h_fill = 'black', p_color = 'white', filename = here("man/figures/tidygeocoder_hex.png"))

plot(s)
```% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{api_parameter_reference}
\alias{api_parameter_reference}
\title{Geocoding service API parameter reference}
\format{
A tibble dataframe
\describe{
\item{method}{Geocoding service name}
\item{generic_name}{Universal parameter name}
\item{api_name}{Name of the parameter for the specified geocoding service}
\item{default_value}{Default value of the parameter}
\item{required}{Is the parameter required by the specified geocoding service?}
}
}
\usage{
api_parameter_reference
}
\description{
This dataset contains the mapping that allows this package to use a
universal syntax to specify parameters for different geocoding services.
Note that latitude and longitude input parameters for reverse geocoding
are not in this dataset and are instead handled directly by the \link{reverse_geo} function.

The \code{generic_name} column is a universal parameter name that is shared between services.
The \code{api_name} column is the parameter name for the given geocoding service specified by the
\code{method} column. When \code{generic_name} is missing
this means the parameter is specific to that geocoding service.

While the "census" and "google" services do not have a \code{limit}
argument in their APIs, tidygeocoder provides a passthrough so you can still
use the \code{limit} argument in \link{geo} and \link{reverse_geo} to limit the
number of results per input.

Note that some geocoding services only use the \code{limit} argument for forward geocoding.
Refer to API documentation of each service for more information.

Reference the documentation for \link{geo} and \link{reverse_geo} for more information.
Also reference \code{vignette("tidygeocoder")} for more details on constructing API queries.
}
\details{
The API documentation for each service is linked to below:
\itemize{
\item \href{https://nominatim.org/release-docs/develop/api/Search/}{Nominatim}
\item \href{https://www.census.gov/programs-surveys/geography/technical-documentation/complete-technical-documentation/census-geocoder.html}{US Census}
\item \href{https://developers.arcgis.com/rest/geocode/api-reference/overview-world-geocoding-service.htm}{ArcGIS}
\item \href{https://www.geocod.io/docs/}{Geocodio}
\item \href{https://locationiq.com/docs}{Location IQ}
\item \href{https://developers.google.com/maps/documentation/geocoding/overview}{Google}
\item \href{https://opencagedata.com/api}{OpenCage}
\item \href{https://docs.mapbox.com/api/search/geocoding/}{Mapbox}
\item \href{https://developer.here.com/documentation/geocoding-search-api/dev_guide/index.html}{HERE}
\item \href{https://developer.tomtom.com/search-api/search-api-documentation-geocoding/geocode}{TomTom}
\item \href{https://developer.mapquest.com/documentation/geocoding-api/}{MapQuest}
\item \href{https://docs.microsoft.com/en-us/bingmaps/rest-services/locations/}{Bing}
\item \href{https://apidocs.geoapify.com/docs/geocoding/api/}{Geoapify}
}
}
\seealso{
\link{geo} \link{reverse_geo} \link{get_api_query} \link{query_api} \link{min_time_reference} \link{batch_limit_reference}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geocode_combine.R
\name{geo_combine}
\alias{geo_combine}
\title{Combine multiple geocoding queries}
\usage{
geo_combine(
  queries,
  global_params = list(),
  address = NULL,
  street = NULL,
  city = NULL,
  county = NULL,
  state = NULL,
  postalcode = NULL,
  country = NULL,
  lat = lat,
  long = long,
  ...
)
}
\arguments{
\item{queries}{a list of queries, each provided as a list of parameters. The queries are
executed by the \link{geocode} function in the order provided.
(ex. \code{list(list(method = 'osm'), list(method = 'census'), ...)})}

\item{global_params}{a list of parameters to be used for all queries
(ex. \code{list(address = 'address', full_results = TRUE)})}

\item{address}{single line address (ie. '1600 Pennsylvania Ave NW, Washington, DC').
Do not combine with the address component arguments below
(\code{street}, \code{city}, \code{county}, \code{state}, \code{postalcode}, \code{country}).}

\item{street}{street address (ie. '1600 Pennsylvania Ave NW')}

\item{city}{city (ie. 'Tokyo')}

\item{county}{county (ie. 'Jefferson')}

\item{state}{state (ie. 'Kentucky')}

\item{postalcode}{postalcode (ie. zip code if in the United States)}

\item{country}{country (ie. 'Japan')}

\item{lat}{latitude column name. Can be quoted or unquoted (ie. \code{lat} or \code{"lat"}).}

\item{long}{longitude column name. Can be quoted or unquoted (ie. \code{long} or `"long"``).}

\item{...}{arguments passed to the \link{geocode_combine} function}
}
\value{
tibble (dataframe)
}
\description{
Passes address inputs in character vector form to the
\link{geocode_combine} function for geocoding.

Note that address inputs must be specified for queries either with the \code{queries} parameter (for each query)
or the \code{global_params} parameter (for all queries). For example \code{global_params = list(address = 'address')}
passes addresses provided in the \code{address} parameter to all queries.
}
\examples{
\donttest{

options(tidygeocoder.progress_bar = FALSE)
example_addresses <- c("100 Main St New York, NY", "Paris", "Not a Real Address")

geo_combine(
    queries = list(
        list(method = 'census'),
        list(method = 'osm')
    ),
    address = example_addresses,
    global_params = list(address = 'address')
  )
  
geo_combine(
  queries = list(
      list(method = 'arcgis'), 
      list(method = 'census', mode = 'single'),
      list(method = 'census', mode = 'batch')
  ),
  global_params = list(address = 'address'),
  address = example_addresses,
  cascade = FALSE,
  return_list = TRUE
)

geo_combine(
   queries = list(
      list(method = 'arcgis', address = 'city'),
      list(method = 'osm', city = 'city', country = 'country')
   ),
   city = c('Tokyo', 'New York'),
   country = c('Japan', 'United States'),
   cascade = FALSE
)
}
}
\seealso{
\link{geocode_combine} \link{geo} \link{geocode}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results_processing.R
\name{extract_reverse_results}
\alias{extract_reverse_results}
\title{Extract reverse geocoding results}
\usage{
extract_reverse_results(
  method,
  response,
  full_results = TRUE,
  flatten = TRUE,
  limit = 1
)
}
\arguments{
\item{method}{method name}

\item{response}{content from the geocoding service (returned by the \link{query_api} function)}

\item{full_results}{if TRUE then the full results (not just an address column)
will be returned.}

\item{flatten}{if TRUE then flatten any nested dataframe content}

\item{limit}{only used for the "google"
method(s). Limits number of results per coordinate.}
}
\value{
geocoding results in tibble format
}
\description{
Parses the output of the \link{query_api} function for reverse geoocoding.
The address is extracted into the first column
of the returned dataframe. This function is not used for batch
geocoded results. Refer to \link{query_api} for example
usage.
}
\seealso{
\link{get_api_query} \link{query_api} \link{reverse_geo}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geocode.R
\name{geocode}
\alias{geocode}
\title{Geocode addresses in a dataframe}
\usage{
geocode(
  .tbl,
  address = NULL,
  street = NULL,
  city = NULL,
  county = NULL,
  state = NULL,
  postalcode = NULL,
  country = NULL,
  lat = "lat",
  long = "long",
  return_input = TRUE,
  limit = 1,
  return_addresses = NULL,
  unique_only = FALSE,
  ...
)
}
\arguments{
\item{.tbl}{dataframe containing addresses}

\item{address}{single line street address column name. Do not combine with
address component arguments (\code{street}, \code{city}, \code{county}, \code{state}, \code{postalcode}, \code{country})}

\item{street}{street address column name}

\item{city}{city column name}

\item{county}{county column name}

\item{state}{state column name}

\item{postalcode}{postalcode column name (zip code if in the United States)}

\item{country}{country column name}

\item{lat}{latitude column name. Can be quoted or unquoted (ie. lat or 'lat').}

\item{long}{longitude column name. Can be quoted or unquoted (ie. long or 'long').}

\item{return_input}{if TRUE then the input dataset will be combined with the geocoder query results
and returned. If FALSE only the geocoder results will be returned.}

\item{limit}{maximum number of results to return per input address. For many geocoding services
the maximum value of the limit parameter is 100. Pass \code{limit = NULL} to use
the default \code{limit} value of the selected geocoding service.
For batch geocoding, limit must be set to 1 (default) if \code{return_addresses = TRUE}.
To use \code{limit > 1} or \code{limit = NULL} set return_input to FALSE.
Refer to \link{api_parameter_reference} for more details.}

\item{return_addresses}{if TRUE return input addresses. Defaults to TRUE if \code{return_input} is FALSE
and FALSE if \code{return_input} is TRUE. This argument is passed to the \code{geo()} function.}

\item{unique_only}{if TRUE then only unique results will be returned and
return_input will be set to FALSE.}

\item{...}{arguments passed to the \link{geo} function}
}
\value{
tibble (dataframe)
}
\description{
Takes a dataframe containing addresses as an input and returns
the results from a specified geocoding service in a dataframe format using the
\link{geo} function. See example usage in \code{vignette("tidygeocoder")}.

This function passes all additional parameters (\code{...}) to the
\link{geo} function, so you can refer to its documentation for more details
on possible arguments.

Note that the arguments used for specifying address columns (\code{address},
\code{street}, \code{city}, \code{county}, \code{state}, \code{postalcode}, and \code{country}) accept either
quoted or unquoted column names (ie. \code{"address_col"} and \code{address_col} are
both acceptable).
}
\examples{
\donttest{
library(dplyr, warn.conflicts = FALSE)
sample_addresses \%>\% slice(1:2) \%>\%
 geocode(addr, method = 'arcgis')

louisville \%>\% head(2) \%>\%
 geocode(street = street, city = city, state = state,
  postalcode = zip, method = 'census', full_results = TRUE)

sample_addresses \%>\% slice(8:9) \%>\%
 geocode(addr, method = 'osm', limit = 2,
  return_input = FALSE, full_results = TRUE)

sample_addresses \%>\% slice(4:5) \%>\%
 geocode(addr, method = 'arcgis',
  lat = latitude, long = longitude,
  full_results = TRUE)
}
}
\seealso{
\link{geo}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reverse_geo.R
\name{reverse_geo}
\alias{reverse_geo}
\title{Reverse geocode coordinates}
\usage{
reverse_geo(
  lat,
  long,
  method = "osm",
  address = "address",
  limit = 1,
  full_results = FALSE,
  mode = "",
  unique_only = FALSE,
  return_coords = TRUE,
  min_time = NULL,
  progress_bar = show_progress_bar(),
  quiet = getOption("tidygeocoder.quiet", FALSE),
  api_url = NULL,
  timeout = 20,
  flatten = TRUE,
  batch_limit = NULL,
  verbose = getOption("tidygeocoder.verbose", FALSE),
  no_query = FALSE,
  custom_query = list(),
  api_options = list(),
  iq_region = "us",
  geocodio_v = 1.6,
  mapbox_permanent = FALSE,
  here_request_id = NULL,
  mapquest_open = FALSE
)
}
\arguments{
\item{lat}{latitude values (input data)}

\item{long}{longitude values (input data)}

\item{method}{the geocoding service to be used. API keys are loaded from environmental variables. Run \code{usethis::edit_r_environ()} to open your .Renviron file and add an API key as an environmental variable. For example, add the line \code{GEOCODIO_API_KEY="YourAPIKeyHere"}.
\itemize{
\item \code{"osm"}: \href{https://nominatim.org}{Nominatim}.
\item \code{"arcgis"}: \href{https://developers.arcgis.com/rest/geocode/api-reference/overview-world-geocoding-service.htm}{ArcGIS}.
\item \code{"geocodio"}: \href{https://www.geocod.io/}{Geocodio}. Geographic coverage is limited to the United States and Canada. An API key must be stored in the environmental variable "GEOCODIO_API_KEY". Batch geocoding is supported.
\item \code{"iq"}: \href{https://locationiq.com/}{Location IQ}.  An API key must be stored in the environmental variable "LOCATIONIQ_API_KEY".
\item \code{"google"}: \href{https://developers.google.com/maps/documentation/geocoding/overview}{Google}.  An API key must be stored in the environmental variable "GOOGLEGEOCODE_API_KEY".
\item \code{"opencage"}: \href{https://opencagedata.com}{OpenCage}.  An API key must be stored in the environmental variable "OPENCAGE_KEY".
\item \code{"mapbox"}: \href{https://docs.mapbox.com/api/search/}{Mapbox}.  An API key must be stored in the environmental variable "MAPBOX_API_KEY".
\item \code{"here"}: \href{https://developer.here.com/products/geocoding-and-search}{HERE}.  An API key must be stored in the environmental variable "HERE_API_KEY". Batch geocoding is supported, but must be explicitly called with \code{mode = "batch"}.
\item \code{"tomtom"}: \href{https://developer.tomtom.com/search-api/search-api-documentation/geocoding}{TomTom}.  An API key must be stored in the environmental variable "TOMTOM_API_KEY". Batch geocoding is supported.
\item \code{"mapquest"}: \href{https://developer.mapquest.com/documentation/geocoding-api/}{MapQuest}.  An API key must be stored in the environmental variable "MAPQUEST_API_KEY". Batch geocoding is supported.
\item \code{"bing"}: \href{https://docs.microsoft.com/en-us/bingmaps/rest-services/locations/}{Bing}.  An API key must be stored in the environmental variable "BINGMAPS_API_KEY". Batch geocoding is supported, but must be explicitly called with \code{mode = "batch"}.
\item \code{"geoapify"}: \href{https://www.geoapify.com/geocoding-api}{Geoapify}.  An API key must be stored in the environmental variable "GEOAPIFY_KEY".
}}

\item{address}{name of the address column (in the output data)}

\item{limit}{maximum number of results to return per input coordinate. For many geocoding services
the maximum value of the limit parameter is 100. Pass \code{limit = NULL} to use
the default \code{limit} value of the selected geocoding service.
For batch geocoding, limit must be set to 1 (default) if \code{return_coords = TRUE}.
Refer to \link{api_parameter_reference} for more details.}

\item{full_results}{returns all available data from the geocoding service if TRUE.
If FALSE (default) then only a single address column is returned from the geocoding service.}

\item{mode}{set to 'batch' to force batch geocoding or 'single' to force single coordinate
geocoding (one coordinate per query). If not specified then batch geocoding will
be used if available (given method selected) when multiple coordinates are
provided; otherwise single address geocoding will be used. For the "here" and "bing" methods the
batch mode should be explicitly specified with \code{mode = 'batch'}.}

\item{unique_only}{only return results for unique inputs if TRUE}

\item{return_coords}{return input coordinates with results if TRUE. Note that
most services return the input coordinates with \code{full_results = TRUE} and setting
\code{return_coords} to FALSE does not prevent this.}

\item{min_time}{minimum amount of time for a query to take (in seconds). If NULL
then min_time will be set to the default value specified in \link{min_time_reference}.}

\item{progress_bar}{if TRUE then a progress bar will be displayed
for single input geocoding (1 input per query). By default the progress bar
will not be shown for code executed when knitting R Markdown files or code within
an RStudio notebook chunk. Can be set permanently with \code{options(tidygeocoder.progress_bar = FALSE)}.}

\item{quiet}{if TRUE then console messages that are displayed by default
regarding queries will be suppressed. FALSE is default.
Can be set permanently with \code{options(tidygeocoder.quiet = TRUE)}.}

\item{api_url}{custom API URL. If specified, the default API URL will be overridden.
This parameter can be used to specify a local Nominatim server, for instance.}

\item{timeout}{query timeout (in minutes)}

\item{flatten}{if TRUE (default) then any nested dataframes in results are flattened if possible.
Note that in some cases results are flattened regardless such as for Geocodio batch geocoding.}

\item{batch_limit}{limit to the number of coordinates in a batch geocoding query.
Defaults to the value in \link{batch_limit_reference} if not specified.}

\item{verbose}{if TRUE then detailed logs are output to the console. FALSE is default. Can be set
permanently with \code{options(tidygeocoder.verbose = TRUE)}}

\item{no_query}{if TRUE then no queries are sent to the geocoding service and verbose is set to TRUE.
Used for testing.}

\item{custom_query}{API-specific parameters to be used, passed as a named list
(ex. \code{list(extratags = 1)}.}

\item{api_options}{a named list of parameters specific to individual services.
(ex. \code{list(geocodio_v = 1.6, geocodio_hipaa = TRUE)}). Each parameter begins
with the name of the \code{method} (service) it applies to. The possible parameters
are shown below with their default values.
\itemize{
\item \code{census_return_type} (default: \code{"locations"}): set to \code{"geographies"} to return
additional geography columns. Make sure to use \code{full_results = TRUE} if using
the "geographies" setting.
\item \code{iq_region} (default: \code{"us"}): set to "eu" to use the European Union API endpoint
\item \code{geocodio_v} (default: \code{1.6}): the version number of the Geocodio API to be used
\item \code{geocodio_hipaa} (default: \code{FALSE}): set to \code{TRUE} to use the HIPAA compliant
Geocodio API endpoint
\item \code{mapbox_permanent} (default: \code{FALSE}): set to \code{TRUE} to use the \code{mapbox.places-permanent}
endpoint. Note that this option should be used only if you have applied for a permanent
account. Unsuccessful requests made by an account that does not have access to the
endpoint may be billable.
\item \code{mapbox_open} (default: \code{FALSE}): set to \code{TRUE} to use the Open Geocoding endpoint which
relies solely on OpenStreetMap data
\item \code{here_request_id} (default: \code{NULL}): this parameter would return a previous HERE batch job,
identified by its RequestID. The RequestID of a batch job is displayed
when \code{verbose} is TRUE. Note that this option would ignore the
current \code{address} parameter on the request, so the \code{return_addresses} or \code{return_coords}
parameters need to be FALSE.
}}

\item{iq_region}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{geocodio_v}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{mapbox_permanent}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{here_request_id}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{mapquest_open}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}
}
\value{
tibble (dataframe)
}
\description{
Reverse geocodes geographic coordinates (latitude and longitude) given as numeric values.
Latitude and longitude inputs are limited to possible values. Latitudes must be between -90 and 90 and
longitudes must be between -180 and 180. Invalid values will not be sent to the geocoding service.
The \link{reverse_geocode} function utilizes this function on coordinates contained in dataframes.
See example usage in \code{vignette("tidygeocoder")}.

Refer to \link{api_parameter_reference},
\link{min_time_reference}, and \link{batch_limit_reference} for more details on
geocoding service parameters and usage.

This function uses the \link{get_api_query}, \link{query_api}, and
\link{extract_reverse_results} functions to create, execute, and parse geocoder
API queries.
}
\examples{
\donttest{
options(tidygeocoder.progress_bar = FALSE)

 reverse_geo(lat = 38.895865, long = -77.0307713, method = 'osm')
 
 reverse_geo(
   lat = c(38.895865, 43.6534817, 300), 
   long = c(-77.0307713, -79.3839347, 600),
   method = 'osm', full_results = TRUE
 )
 
}
}
\seealso{
\link{reverse_geocode} \link{api_parameter_reference} \link{min_time_reference} \link{batch_limit_reference}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{batch_limit_reference}
\alias{batch_limit_reference}
\title{Geocoding batch size limits}
\format{
A tibble dataframe
\describe{
\item{method}{Geocoding service name}
\item{batch_limit}{The maximum number of addresses or coordinates allowed per batch}
}
}
\usage{
batch_limit_reference
}
\description{
The \link{geo} and \link{reverse_geo} functions use this dataset to set the
maximum batch query size for each service.
}
\seealso{
\link{geo} \link{reverse_geo}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geocode_combine.R
\name{geocode_combine}
\alias{geocode_combine}
\title{Combine multiple geocoding queries}
\usage{
geocode_combine(
  .tbl,
  queries,
  global_params = list(),
  return_list = FALSE,
  cascade = TRUE,
  query_names = NULL,
  lat = "lat",
  long = "long"
)
}
\arguments{
\item{.tbl}{dataframe containing addresses}

\item{queries}{a list of queries, each provided as a list of parameters. The queries are
executed by the \link{geocode} function in the order provided.
(ex. \code{list(list(method = 'osm'), list(method = 'census'), ...)})}

\item{global_params}{a list of parameters to be used for all queries
(ex. \code{list(address = 'address', full_results = TRUE)})}

\item{return_list}{if TRUE then results from each service will be returned as separate
dataframes. If FALSE (default) then all results will be combined into a single dataframe.}

\item{cascade}{if TRUE (default) then only addresses that are not found by a geocoding
service will be attempted by subsequent queries. If FALSE then all queries will
attempt to geocode all addresses.}

\item{query_names}{optional vector with one label for each query provided
(ex. \code{c('geocodio batch', 'geocodio single')}).}

\item{lat}{latitude column name. Can be quoted or unquoted (ie. lat or 'lat').}

\item{long}{longitude column name. Can be quoted or unquoted (ie. long or 'long').}
}
\value{
tibble (dataframe)
}
\description{
Executes multiple geocoding queries on a dataframe input and combines
the results. To use a character vector input instead, see the \link{geo_combine} function.
Queries are executed by the \link{geocode} function. See example usage
in \code{vignette("tidygeocoder")}.

Query results are by default labelled to show which query produced each result. Labels are either
placed in a \code{query} column (if \code{return_list = FALSE}) or used as the names of the returned list
(if \code{return_list = TRUE}). By default the \code{method} parameter value of each query is used as a query label.
If the same \code{method} is used in multiple queries then a number is added according
to the order of the queries (ie. \code{osm1}, \code{osm2}, ...). To provide your own custom query labels
use the \code{query_names} parameter.
}
\examples{
\donttest{

library(dplyr, warn.conflicts = FALSE)

sample_addresses \%>\%
  geocode_combine(
    queries = list(list(method = 'census'), list(method = 'osm')),
    global_params = list(address = 'addr'), cascade = TRUE)

more_addresses <- tibble::tribble(
     ~street_address, ~city, ~state,        ~zip_cd,
     "624 W DAVIS ST #1D",   "BURLINGTON", "NC", 27215,
     "201 E CENTER ST #268", "MEBANE",     "NC", 27302,
     "100 Wall Street",      "New York",   "NY", 10005,
     "Bucharest",            NA,           NA,   NA
     )
 
 more_addresses \%>\%        
   geocode_combine( 
     queries = list(
         list(method = 'census', mode = 'batch'),
         list(method = 'census', mode = 'single'),
         list(method = 'osm')
      ),
     global_params = list(street = 'street_address', 
       city = 'city', state = 'state', postalcode = 'zip_cd'),
     query_names = c('census batch', 'census single', 'osm')
   )
   
 more_addresses \%>\%
   geocode_combine( 
     queries = list(
         list(method = 'census', mode = 'batch', street = 'street_address', 
       city = 'city', state = 'state', postalcode = 'zip_cd'),
         list(method = 'arcgis', address = 'street_address')
      ),
     cascade = FALSE,
     return_list = TRUE
   )
}
}
\seealso{
\link{geo_combine} \link{geo} \link{geocode}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidygeocoder.R
\docType{package}
\name{tidygeocoder-package}
\alias{tidygeocoder}
\alias{tidygeocoder-package}
\title{The tidygeocoder package makes getting data from geocoder services easy.}
\description{
The \link{geocode} and \link{geo} functions are for forward geocoding while
the \link{reverse_geocode} and \link{reverse_geo} functions are for reverse geocoding.
Refer to the documentation on the \code{method} argument in the \link{geo} and \link{reverse_geo} functions
for more details on the available geocoding services. Also see \code{vignette("tidygeocoder")}
for example usage.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://jessecambon.github.io/tidygeocoder/}
  \item \url{https://github.com/jessecambon/tidygeocoder}
  \item Report bugs at \url{https://github.com/jessecambon/tidygeocoder/issues}
}

}
\author{
\strong{Maintainer}: Jesse Cambon \email{jesse.cambon@gmail.com} (\href{https://orcid.org/0000-0001-6854-1514}{ORCID})

Authors:
\itemize{
  \item Diego Hernangómez \email{diego.hernangomezherrero@gmail.com} (\href{https://orcid.org/0000-0001-8457-4658}{ORCID})
  \item Christopher Belanger \email{christopher.a.belanger@gmail.com} (\href{https://orcid.org/0000-0003-2070-5721}{ORCID})
  \item Daniel Possenriede \email{possenriede+r@gmail.com} (\href{https://orcid.org/0000-0002-6738-9845}{ORCID})
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_methods.R
\name{geo_census}
\alias{geo_census}
\alias{geo_osm}
\alias{geo_geocodio}
\alias{geo_iq}
\alias{geo_google}
\alias{geo_opencage}
\alias{geo_mapbox}
\alias{geo_here}
\alias{geo_tomtom}
\alias{geo_mapquest}
\alias{geo_bing}
\alias{geo_arcgis}
\alias{geo_cascade}
\title{Convenience functions for calling \code{geo()}}
\usage{
geo_census(...)

geo_osm(...)

geo_geocodio(...)

geo_iq(...)

geo_google(...)

geo_opencage(...)

geo_mapbox(...)

geo_here(...)

geo_tomtom(...)

geo_mapquest(...)

geo_bing(...)

geo_arcgis(...)

geo_cascade(...)
}
\arguments{
\item{...}{arguments to be passed to the \code{geo} function}
}
\description{
The \code{method} for \code{geo()} is specified in the function name.

\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

Use the \link{geo} function directly instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{louisville}
\alias{louisville}
\title{Louisville, Kentucky street addresses}
\format{
A tibble dataframe with component street addresses
\describe{
\item{street}{Description of the address}
\item{city}{Single line address}
\item{state}{state}
\item{zip}{zip code}
}
}
\source{
Downloaded from \href{https://results.openaddresses.io/sources/us/ky/jefferson}{OpenAddresses.io}
on June 1st 2020
}
\usage{
louisville
}
\description{
Louisville, Kentucky street addresses
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{api_info_reference}
\alias{api_info_reference}
\title{Geocoding service links and information}
\format{
A tibble dataframe
\describe{
\item{method}{Geocoding service name}
\item{method_display_name}{Geocoding service display name}
\item{site_url}{Link to the main site of the geocoding service}
\item{api_documentation_url}{Link to API documentation}
\item{api_usage_policy_url}{Link to the usage policy}
}
}
\usage{
api_info_reference
}
\description{
This dataset is used for generating package documentation.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results_processing.R
\name{extract_results}
\alias{extract_results}
\title{Extract forward geocoding results}
\usage{
extract_results(
  method,
  response,
  full_results = TRUE,
  flatten = TRUE,
  limit = 1
)
}
\arguments{
\item{method}{method name}

\item{response}{content from the geocoding service (returned by the \link{query_api} function)}

\item{full_results}{if TRUE then the full results (not just latitude and longitude)
will be returned.}

\item{flatten}{if TRUE then flatten any nested dataframe content}

\item{limit}{only used for "census" and "google" methods. Limits number of results per address.}
}
\value{
geocoding results in tibble format
}
\description{
Parses the output of the \link{query_api} function for single
address geocoding (ie. not batch geocoding).
Latitude and longitude are extracted into the first two columns
of the returned dataframe.  Refer to  \link{query_api} for example
usage.
}
\seealso{
\link{get_api_query} \link{query_api} \link{geo}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{api_key_reference}
\alias{api_key_reference}
\title{API key environmental variables}
\format{
A tibble dataframe
\describe{
\item{method}{Geocoding service name}
\item{env_var}{Environmental variable name}
}
}
\usage{
api_key_reference
}
\description{
API keys are obtained from environmental variables.
The \link{geo} and \link{reverse_geo} functions use this dataset
to know which environmental variable to use for
each geocoding service.
}
\seealso{
\link{geo} \link{reverse_geo}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reverse_geocode.R
\name{reverse_geocode}
\alias{reverse_geocode}
\title{Reverse geocode coordinates in a dataframe}
\usage{
reverse_geocode(
  .tbl,
  lat,
  long,
  address = "address",
  return_input = TRUE,
  limit = 1,
  return_coords = NULL,
  unique_only = FALSE,
  ...
)
}
\arguments{
\item{.tbl}{dataframe containing coordinates}

\item{lat}{latitude column name (input data). Can be quoted or unquoted (ie. lat or 'lat').}

\item{long}{longitude column name (input data). Can be quoted or unquoted (ie. long or 'long').}

\item{address}{address column name (output data). Can be quoted or unquoted (ie. addr or 'addr').}

\item{return_input}{if TRUE then the input dataset will be combined with the geocoder query results
and returned. If FALSE only the geocoder results will be returned.}

\item{limit}{maximum number of results to return per input coordinate. For many geocoding services
the maximum value of the limit parameter is 100. Pass \code{limit = NULL} to use
the default \code{limit} value of the selected geocoding service.
For batch geocoding, limit must be set to 1 (default) if \code{return_coords = TRUE}.
To use \code{limit > 1} or \code{limit = NULL} set return_input to FALSE.
Refer to \link{api_parameter_reference} for more details.}

\item{return_coords}{if TRUE return input coordinates. Defaults to TRUE if \code{return_input} is
FALSE and FALSE if \code{return_input} is TRUE. This argument is passed to the \code{reverse_geo()} function.}

\item{unique_only}{if TRUE then only unique results will be returned and
return_input will be set to FALSE.}

\item{...}{arguments passed to the \link{reverse_geo} function}
}
\value{
tibble (dataframe)
}
\description{
Takes a dataframe containing coordinates (latitude and longitude) and returns
the reverse geocoding query results from a specified service by using the
\link{reverse_geo} function. See example usage in \code{vignette("tidygeocoder")}.

This function passes all additional parameters (\code{...}) to the
\link{reverse_geo} function, so you can refer to its documentation for more details
on possible arguments.
}
\examples{
\donttest{
library(tibble)
library(dplyr, warn.conflicts = FALSE)

tibble(
    latitude = c(38.895865, 43.6534817),
    longitude = c(-77.0307713,-79.3839347)
  ) \%>\%
  reverse_geocode(
    lat = latitude,
    long = longitude,
    method = 'osm',
    full_results = TRUE
  )

louisville \%>\% head(3) \%>\% 
  reverse_geocode(lat = latitude, long = longitude, 
  method = 'arcgis')

louisville \%>\% head(2) \%>\% 
  reverse_geocode(lat = latitude, long = longitude,  
  method = 'osm',
  limit = 2, return_input = FALSE)

}
}
\seealso{
\link{reverse_geo}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{min_time_reference}
\alias{min_time_reference}
\title{Minimum time required per query}
\format{
A tibble dataframe
\describe{
\item{method}{Geocoding service name}
\item{min_time}{The minimum number of seconds required per query to comply with usage restrictions}
\item{description}{A description of the usage rate restriction}
}
}
\usage{
min_time_reference
}
\description{
The \link{geo} and \link{reverse_geo} functions use this dataset
to set the maximum query rate for each geocoding service.
This rate is based on the usage restriction policies for
each geocoding service.
}
\details{
Links to the usage policies of each geocoding service are below:
\itemize{
\item \href{https://operations.osmfoundation.org/policies/nominatim/}{Nominatim}
\item \href{https://www.census.gov/programs-surveys/geography/technical-documentation/complete-technical-documentation/census-geocoder.html}{US Census}
\item \href{https://developers.arcgis.com/rest/geocode/api-reference/geocoding-free-vs-paid.htm}{ArcGIS}
\item \href{https://www.geocod.io/pricing/}{Geocodio}
\item \href{https://locationiq.com/pricing}{Location IQ}
\item \href{https://developers.google.com/maps/documentation/geocoding/usage-and-billing}{Google}
\item \href{https://opencagedata.com/pricing}{OpenCage}
\item \href{https://www.mapbox.com/pricing/}{Mapbox}
\item \href{https://developer.here.com/pricing}{HERE}
\item \href{https://developer.tomtom.com/store/maps-api}{TomTom}
\item \href{https://developer.mapquest.com/plans}{MapQuest}
\item \href{https://docs.microsoft.com/en-us/bingmaps/spatial-data-services/geocode-and-data-source-limits}{Bing}
\item \href{https://www.geoapify.com/term-and-conditions}{Geoapify}
}
}
\seealso{
\link{geo} \link{reverse_geo}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_factory.R
\name{get_api_query}
\alias{get_api_query}
\title{Construct a geocoder API query}
\usage{
get_api_query(method, generic_parameters = list(), custom_parameters = list())
}
\arguments{
\item{method}{method name (ie. 'census')}

\item{generic_parameters}{universal 'generic' parameters}

\item{custom_parameters}{custom api-specific parameters}
}
\value{
API parameters as a named list
}
\description{
The geocoder API query is created using universal "generic" parameters
and optional api-specific "custom" parameters. Generic parameters
are converted into api parameters using the  \link{api_parameter_reference}
dataset.

The \link{query_api} function executes the queries created
by this function.
}
\examples{
get_api_query("osm", list(address = 'Hanoi, Vietnam'))

get_api_query("census", list(street = '11 Wall St', city = "NY", state = 'NY'),
  list(benchmark = "Public_AR_Census2010"))

}
\seealso{
\link{query_api} \link{api_parameter_reference} \link{geo} \link{reverse_geo}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sample_addresses}
\alias{sample_addresses}
\title{Sample addresses for testing}
\format{
A tibble dataframe with single line addresses
\describe{
\item{name}{Description of the address}
\item{addr}{Single line address}
}
}
\usage{
sample_addresses
}
\description{
Sample addresses for testing
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_factory.R
\name{query_api}
\alias{query_api}
\title{Execute a geocoder API query}
\usage{
query_api(
  api_url,
  query_parameters,
  mode = "single",
  batch_file = NULL,
  input_list = NULL,
  content_encoding = "UTF-8",
  timeout = 20,
  method = ""
)
}
\arguments{
\item{api_url}{Base URL of the API. query parameters are appended to this}

\item{query_parameters}{api query parameters in the form of a named list}

\item{mode}{determines the type of query to execute\preformatted{- "single": geocode a single input (all methods)
- "list": batch geocode a list of inputs (ex. geocodio)
- "file": batch geocode a file of inputs (ex. census)
}}

\item{batch_file}{a csv file of input data to upload (for \code{mode = 'file'})}

\item{input_list}{a list of input data (for \code{mode = 'list'})}

\item{content_encoding}{Encoding to be used for parsing content}

\item{timeout}{timeout in minutes}

\item{method}{if 'mapquest' or 'arcgis' then the query status code is changed appropriately}
}
\value{
a named list containing the response content (\code{content}) and the HTTP request status (\code{status})
}
\description{
The \link{get_api_query} function can create queries for this
function to execute.
}
\examples{
\donttest{
raw1 <- query_api("http://nominatim.openstreetmap.org/search", 
   get_api_query("osm", list(address = 'Hanoi, Vietnam')))
   
raw1$status
   
extract_results('osm', jsonlite::fromJSON(raw1$content))

raw2 <- query_api("http://nominatim.openstreetmap.org/reverse", 
   get_api_query("osm", custom_parameters = list(lat = 38.895865, lon = -77.0307713)))
   
extract_reverse_results('osm', jsonlite::fromJSON(raw2$content))
}

}
\seealso{
\link{get_api_query} \link{extract_results} \link{extract_reverse_results} \link{geo} \link{reverse_geo}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo.R
\name{geo}
\alias{geo}
\title{Geocode addresses}
\usage{
geo(
  address = NULL,
  street = NULL,
  city = NULL,
  county = NULL,
  state = NULL,
  postalcode = NULL,
  country = NULL,
  method = "osm",
  cascade_order = c("census", "osm"),
  lat = "lat",
  long = "long",
  limit = 1,
  full_results = FALSE,
  mode = "",
  unique_only = FALSE,
  return_addresses = TRUE,
  min_time = NULL,
  progress_bar = show_progress_bar(),
  quiet = getOption("tidygeocoder.quiet", FALSE),
  api_url = NULL,
  timeout = 20,
  flatten = TRUE,
  batch_limit = NULL,
  batch_limit_error = TRUE,
  verbose = getOption("tidygeocoder.verbose", FALSE),
  no_query = FALSE,
  custom_query = list(),
  api_options = list(),
  return_type = "locations",
  iq_region = "us",
  geocodio_v = 1.6,
  param_error = TRUE,
  mapbox_permanent = FALSE,
  here_request_id = NULL,
  mapquest_open = FALSE
)
}
\arguments{
\item{address}{single line address (ie. '1600 Pennsylvania Ave NW, Washington, DC').
Do not combine with the address component arguments below
(\code{street}, \code{city}, \code{county}, \code{state}, \code{postalcode}, \code{country}).}

\item{street}{street address (ie. '1600 Pennsylvania Ave NW')}

\item{city}{city (ie. 'Tokyo')}

\item{county}{county (ie. 'Jefferson')}

\item{state}{state (ie. 'Kentucky')}

\item{postalcode}{postalcode (ie. zip code if in the United States)}

\item{country}{country (ie. 'Japan')}

\item{method}{the geocoding service to be used. API keys are loaded from environmental variables. Run \code{usethis::edit_r_environ()} to open your .Renviron file and add an API key as an environmental variable. For example, add the line \code{GEOCODIO_API_KEY="YourAPIKeyHere"}.
\itemize{
\item \code{"osm"}: \href{https://nominatim.org}{Nominatim}.
\item \code{"census"}: \href{https://geocoding.geo.census.gov/}{US Census}. Geographic coverage is limited to the United States.  Batch geocoding is supported.
\item \code{"arcgis"}: \href{https://developers.arcgis.com/rest/geocode/api-reference/overview-world-geocoding-service.htm}{ArcGIS}.
\item \code{"geocodio"}: \href{https://www.geocod.io/}{Geocodio}. Geographic coverage is limited to the United States and Canada. An API key must be stored in the environmental variable "GEOCODIO_API_KEY". Batch geocoding is supported.
\item \code{"iq"}: \href{https://locationiq.com/}{Location IQ}.  An API key must be stored in the environmental variable "LOCATIONIQ_API_KEY".
\item \code{"google"}: \href{https://developers.google.com/maps/documentation/geocoding/overview}{Google}.  An API key must be stored in the environmental variable "GOOGLEGEOCODE_API_KEY".
\item \code{"opencage"}: \href{https://opencagedata.com}{OpenCage}.  An API key must be stored in the environmental variable "OPENCAGE_KEY".
\item \code{"mapbox"}: \href{https://docs.mapbox.com/api/search/}{Mapbox}.  An API key must be stored in the environmental variable "MAPBOX_API_KEY".
\item \code{"here"}: \href{https://developer.here.com/products/geocoding-and-search}{HERE}.  An API key must be stored in the environmental variable "HERE_API_KEY". Batch geocoding is supported, but must be explicitly called with \code{mode = "batch"}.
\item \code{"tomtom"}: \href{https://developer.tomtom.com/search-api/search-api-documentation/geocoding}{TomTom}.  An API key must be stored in the environmental variable "TOMTOM_API_KEY". Batch geocoding is supported.
\item \code{"mapquest"}: \href{https://developer.mapquest.com/documentation/geocoding-api/}{MapQuest}.  An API key must be stored in the environmental variable "MAPQUEST_API_KEY". Batch geocoding is supported.
\item \code{"bing"}: \href{https://docs.microsoft.com/en-us/bingmaps/rest-services/locations/}{Bing}.  An API key must be stored in the environmental variable "BINGMAPS_API_KEY". Batch geocoding is supported, but must be explicitly called with \code{mode = "batch"}.
\item \code{"geoapify"}: \href{https://www.geoapify.com/geocoding-api}{Geoapify}.  An API key must be stored in the environmental variable "GEOAPIFY_KEY".
\item \code{"cascade"} \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use \link{geocode_combine} or \link{geo_combine} instead.
The "cascade" method first uses one geocoding service and then uses
a second geocoding service if the first service didn't return results.
The services and order is specified by the cascade_order argument.
Note that this is not compatible with \code{full_results = TRUE} as geocoding
services have different columns that they return.
}}

\item{cascade_order}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} a vector with two character values for the
method argument in the order in which the geocoding services will be attempted for \code{method = "cascade"}
(ie. \code{c("census", "geocodio")})}

\item{lat}{latitude column name. Can be quoted or unquoted (ie. \code{lat} or \code{"lat"}).}

\item{long}{longitude column name. Can be quoted or unquoted (ie. \code{long} or `"long"``).}

\item{limit}{maximum number of results to return per input address. For many geocoding services
the maximum value of the limit parameter is 100. Pass \code{limit = NULL} to use
the default \code{limit} value of the selected geocoding service.
For batch geocoding, limit must be set to 1 (default) if \code{return_addresses = TRUE}.
Refer to \link{api_parameter_reference} for more details.}

\item{full_results}{returns all available data from the geocoding service if TRUE.
If FALSE (default) then only latitude and longitude columns are returned from the geocoding service.}

\item{mode}{set to 'batch' to force batch geocoding or 'single' to force single address
geocoding (one address per query). If not specified then batch geocoding will
be used if available (given method selected) when multiple addresses are
provided; otherwise single address geocoding will be used. For the "here" and "bing" methods the
batch mode should be explicitly specified with \code{mode = 'batch'}.}

\item{unique_only}{only return results for unique inputs if TRUE}

\item{return_addresses}{return input addresses with results if TRUE. Note that
most services return the input addresses with \code{full_results = TRUE} and setting
return_addresses to FALSE does not prevent this.}

\item{min_time}{minimum amount of time for a query to take (in seconds). If NULL
then min_time will be set to the default value specified in \link{min_time_reference}.}

\item{progress_bar}{if TRUE then a progress bar will be displayed
for single input geocoding (1 input per query). By default the progress bar
will not be shown for code executed when knitting R Markdown files or code within
an RStudio notebook chunk. Can be set permanently with \code{options(tidygeocoder.progress_bar = FALSE)}.}

\item{quiet}{if TRUE then console messages that are displayed by default
regarding queries will be suppressed. FALSE is default.
Can be set permanently with \code{options(tidygeocoder.quiet = TRUE)}.}

\item{api_url}{custom API URL. If specified, the default API URL will be overridden.
This parameter can be used to specify a local Nominatim server, for instance.}

\item{timeout}{query timeout (in minutes)}

\item{flatten}{if TRUE (default) then any nested dataframes in results are flattened if possible.
Note that in some cases results are flattened regardless such as for Geocodio batch geocoding.}

\item{batch_limit}{limit to the number of addresses in a batch geocoding query.
Defaults to the value in \link{batch_limit_reference} if not specified.}

\item{batch_limit_error}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} if TRUE then an error is thrown if the number of addresses exceeds the batch limit.
(if executing a batch query). This is reverted to FALSE when using the cascade method.}

\item{verbose}{if TRUE then detailed logs are output to the console. FALSE is default. Can be set
permanently with \code{options(tidygeocoder.verbose = TRUE)}}

\item{no_query}{if TRUE then no queries are sent to the geocoding service and verbose is set to TRUE.
Used for testing.}

\item{custom_query}{API-specific parameters to be used, passed as a named list
(ex. \code{list(extratags = 1)}.}

\item{api_options}{a named list of parameters specific to individual services.
(ex. \code{list(geocodio_v = 1.6, geocodio_hipaa = TRUE)}). Each parameter begins
with the name of the \code{method} (service) it applies to. The possible parameters
are shown below with their default values.
\itemize{
\item \code{census_return_type} (default: \code{"locations"}): set to \code{"geographies"} to return
additional geography columns. Make sure to use \code{full_results = TRUE} if using
the "geographies" setting.
\item \code{iq_region} (default: \code{"us"}): set to "eu" to use the European Union API endpoint
\item \code{geocodio_v} (default: \code{1.6}): the version number of the Geocodio API to be used
\item \code{geocodio_hipaa} (default: \code{FALSE}): set to \code{TRUE} to use the HIPAA compliant
Geocodio API endpoint
\item \code{mapbox_permanent} (default: \code{FALSE}): set to \code{TRUE} to use the \code{mapbox.places-permanent}
endpoint. Note that this option should be used only if you have applied for a permanent
account. Unsuccessful requests made by an account that does not have access to the
endpoint may be billable.
\item \code{mapbox_open} (default: \code{FALSE}): set to \code{TRUE} to use the Open Geocoding endpoint which
relies solely on OpenStreetMap data
\item \code{here_request_id} (default: \code{NULL}): this parameter would return a previous HERE batch job,
identified by its RequestID. The RequestID of a batch job is displayed
when \code{verbose} is TRUE. Note that this option would ignore the
current \code{address} parameter on the request, so the \code{return_addresses} or \code{return_coords}
parameters need to be FALSE.
}}

\item{return_type}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{iq_region}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{geocodio_v}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{param_error}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} if TRUE then an error will be thrown if any address
parameters are used that are invalid for the selected service (\code{method}).
If \code{method = "cascade"} then no errors will be thrown.}

\item{mapbox_permanent}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{here_request_id}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}

\item{mapquest_open}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} use the \code{api_options} parameter instead}
}
\value{
tibble (dataframe)
}
\description{
Geocodes addresses given as character values. The \link{geocode}
function utilizes this function on addresses contained in dataframes.
See example usage in \code{vignette("tidygeocoder")}.

Note that not all geocoding services support certain address component
parameters. For example, the Census geocoder only covers the United States
and does not have a "country" parameter.

Refer to \link{api_parameter_reference},
\link{min_time_reference}, and \link{batch_limit_reference} for more details on
geocoding service parameters and usage.

This function uses the \link{get_api_query}, \link{query_api}, and
\link{extract_results} functions to create, execute, and parse geocoder
API queries.
}
\examples{
\donttest{
options(tidygeocoder.progress_bar = FALSE)

geo(street = "600 Peachtree Street NE", city = "Atlanta",
 state = "Georgia", method = "census")

geo(address = c("Tokyo, Japan", "Lima, Peru", "Nairobi, Kenya"),
 method = 'osm')
 
geo("100 Main St New York, NY",  full_results = TRUE,
 method = "census", api_options = list(census_return_type = 'geographies'))

geo(county = 'Jefferson', state = "Kentucky", country = "US",
     method = 'osm')
}
}
\seealso{
\link{geocode} \link{api_parameter_reference} \link{min_time_reference} \link{batch_limit_reference}
}
