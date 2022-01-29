
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

# osmextract <a href='https://docs.ropensci.org/osmextract/'><img src='man/figures/logo.svg' align="right" height=275/></a>

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/osmextract/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/osmextract/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/osmextract/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/osmextract?branch=master)
[![peer-review](https://badges.ropensci.org/395_status.svg)](https://github.com/ropensci/software-review/issues/395)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN
status](https://www.r-pkg.org/badges/version/osmextract)](https://CRAN.R-project.org/package=osmextract)
<!-- badges: end -->

The goal of `osmextract` is to make it easier for people to access
OpenStreetMap (OSM) data for reproducible research. OSM data is the
premier source of freely available, community created geographic data
worldwide. We aim to enable you to extract it for data-driven work in
the public interest.

`osmextract` matches, downloads, converts and imports bulk OSM data
hosted by providers such as [Geofabrik
GmbH](http://download.geofabrik.de) and
[bbbike](https://download.bbbike.org/osm/). For information on
alternative providers and how to add them see the [providers
vignette](https://docs.ropensci.org/osmextract/articles/providers.html).

## Why osmextract?

The package answers a common question for researchers who use OSM data:
how to get it into a statistical environment, in an appropriate format,
as part of a computationally efficient and reproducible workflow? Other
packages answer parts of this question.
[`osmdata`](https://github.com/ropensci/osmdata), for example, is an R
package that provides an R interface to the [Overpass
API](https://wiki.openstreetmap.org/wiki/Overpass_API), which is ideal
for downloading small OSM datasets. However, the API is rate limited,
making it hard to download large datasets. As a case study, try to
download all cycleways in England using `osmdata`:

``` r
library(osmdata)
cycleways_england = opq("England") %>% 
  add_osm_feature(key = "highway", value = "cycleway") %>% 
  osmdata_sf()
# Error in check_for_error(doc) : General overpass server error; returned:
# The data included in this document is from www.openstreetmap.org. The data is made available under ODbL. runtime error: Query timed out in "query" at line 4 after 26 seconds. 
```

The query stops with an error message after around 30 seconds. The same
query can be made with `osmextract` as follows, which reads-in almost
100k linestrings in less than 10 seconds, after the data has been
downloaded in the compressed `.pbf` format and converted to the open
standard `.gpkg` format. The download-and-conversion operation of the
OSM extract associated to England takes approximately a few minutes, but
this operation must be executed only once. The following code chunk is
not evaluated.

``` r
library(osmextract)

cycleways_england = oe_get(
  "England",
  quiet = FALSE,
  query = "SELECT * FROM 'lines' WHERE highway = 'cycleway'"
)
par(mar = rep(0.1, 4))
plot(sf::st_geometry(cycleways_england))
```

<img src="man/figures/89990770-22bf2480-dc83-11ea-9092-764594534959.png" width="80%" style="display: block; margin: auto;" />

The package is designed to complement `osmdata`, which has advantages
over `osmextract` for small datasets: `osmdata` is likely to be quicker
for datasets less than a few MB in size, provides up-to-date data and
has an intuitive interface. `osmdata` can provide data in a range of
formats, while `osmextract` only returns
[`sf`](https://github.com/r-spatial/sf) objects. `osmextract`’s niche is
that it provides a fast way to download large OSM datasets in the highly
compressed `pbf` format and read them in via the fast C library
[GDAL](https://gdal.org/drivers/vector/osm.html) and the popular R
package for working with geographic data
[`sf`](https://github.com/r-spatial/sf).

## Installation

You can install the released version of `osmextract` from
[CRAN](https://cran.r-project.org/package=osmextract) with:

``` r
install.packages("osmextract")
```

You can install the development version from
[GitHub](https://github.com/ropensci/osmextract) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/osmextract")
```

Load the package with:

``` r
library(osmextract)
#> Data (c) OpenStreetMap contributors, ODbL 1.0. https://www.openstreetmap.org/copyright.
#> Check the package website, https://docs.ropensci.org/osmextract/, for more details.
```

To use alongside functionality in the `sf` package, we also recommend
attaching this geographic data package as follows:

``` r
library(sf)
#> Linking to GEOS 3.9.1, GDAL 3.2.1, PROJ 7.2.1
```

### Warnings:

The functions defined in this package may return a warning message like

    st_crs<- : replacing crs does not reproject data; use st_transform for that 

if the user is running an old version of GDAL (&lt;= 3.0.0) or PROJ
(&lt;= 6.0.0). See [here](https://github.com/r-spatial/sf/issues/1419)
for more details. Nevertheless, every function should still work
correctly. Please, raise [a new
issue](https://github.com/ropensci/osmextract/issues) if you find any
odd behaviour.

## Basic usage

Give `osmextract` a place name and it will try to find it in a list of
names in the specified provider
([Geofabrik](https://www.geofabrik.de/data/download.html) by default).
If the name you give it matches a place, it will download and import the
associated data into R. The function `oe_get()` downloads (if not
already downloaded) and reads-in data from OSM providers as `sf`
objects. By default `oe_get()` imports the `lines` layer, but any layer
can be read-in by changing the `layer` argument:

``` r
osm_lines = oe_get("Isle of Wight", stringsAsFactors = FALSE, quiet = TRUE)
osm_points = oe_get("Isle of Wight", layer = "points", stringsAsFactors = FALSE, quiet = TRUE)
nrow(osm_lines)
#> [1] 48281
nrow(osm_points)
#> [1] 61329
par(mar = rep(0, 4))
plot(st_geometry(osm_lines), xlim = c(-1.59, -1.1), ylim = c(50.5, 50.8))
plot(st_geometry(osm_points), xlim = c(-1.59, -1.1), ylim = c(50.5, 50.8))
```

<img src="man/figures/README-points-lines-iow-1.png" width="50%" /><img src="man/figures/README-points-lines-iow-2.png" width="50%" />

The figures above give an insight into the volume and richness of data
contained in OSM extracts. Even for a small island such as the Isle of
Wight, it contains over 50k features including ferry routes, shops and
roads. The column names in the `osm_lines` object are as follows:

``` r
names(osm_lines) # default variable names
#>  [1] "osm_id"     "name"       "highway"    "waterway"   "aerialway" 
#>  [6] "barrier"    "man_made"   "z_order"    "other_tags" "geometry"
```

Once imported, you can use all functions for data frames in base R and
other packages. You can also use functions from the `sf` package for
spatial analysis and visualisation. Let’s plot all the major, secondary
and residential roads, for example:

``` r
ht = c("primary", "secondary", "tertiary", "unclassified") # highway types of interest
osm_major_roads = osm_lines[osm_lines$highway %in% ht, ]
plot(osm_major_roads["highway"], key.pos = 1)
```

<img src="man/figures/README-iow1-1.png" width="100%" />

The same steps can be used to get other OSM datasets (examples not run):

``` r
malta = oe_get("Malta", quiet = TRUE)
andorra = oe_get("Andorra", extra_tags = "ref")
leeds = oe_get("Leeds")
goa = oe_get("Goa", query = "SELECT highway, geometry FROM 'lines'")
```

If the input place does not match any of the existing names in the
supported providers, then `oe_get()` will try to geocode it via
[Nominatim
API](https://nominatim.org/release-docs/develop/api/Overview/), and it
will select the smallest OSM extract intersecting the area. For example
(not run):

``` r
oe_get("Milan") # Warning: It will download more than 400MB of data
#> No exact match found for place = Milan and provider = geofabrik. Best match is Iran.
#> Checking the other providers.
#> No exact match found in any OSM provider data. Searching for the location online.
#> ... (extra messages here)
```

For further details on using the package, see the [Introducing
osmextract
vignette](https://docs.ropensci.org/osmextract/articles/osmextract.html).

## Persistent download directory

The default behaviour of `oe_get()` is to save all the files in a
temporary directory, which is erased every time you restart your R
session. If you want to set a directory that will persist, you can add
`OSMEXT_DOWNLOAD_DIRECTORY=/path/for/osm/data` in your `.Renviron` file,
e.g. with:

``` r
usethis::edit_r_environ()
# Add a line containing: OSMEXT_DOWNLOAD_DIRECTORY=/path/to/save/files
```

We strongly advise you setting a persistent directory since working with
`.pbf` files is an expensive operation, that is skipped by `oe_*()`
functions if they detect that the input `.pbf` file was already
downloaded.

You can always check the default `download_directory` used by `oe_get()`
with:

``` r
oe_download_directory()
```

<!-- The following section was removed since now oe_download sets the timeout value. See https://github.com/ropensci/osmextract/issues/222 -->
<!-- ## Troubleshooting -->
<!-- Depending on the `.pbf` file selected and your connection speed, you may experience an error stating `Timeout of 60 seconds was reached`.  -->
<!-- If so, before calling `oe_get()`, you can adjust the timeout using `options(timeout = 300)`, choosing an appropriate value.  -->
<!-- This setting affects all calls to [download.file()](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/download.file), so you may need to reset it for the rest of your script. -->
<!-- If you need to update an existing `.pbf` file or replace an incomplete extract, you can use the argument `force_download`, i.e `oe_get("some-place", force_download = TRUE)`.    -->
<!-- Check `?oe_get` and `?oe_download` for more details.  -->

## Next steps

We would love to see more providers added (see the [Add new
OpenStreetMap
providers](https://docs.ropensci.org/osmextract/articles/providers.html)
for details) and see what people can do with OSM datasets of the type
provided by this package in a reproducible and open statistical
programming environment for the greater good. Any contributions to
support this or any other improvements to the package are very welcome
via our issue tracker.

## Licence

We hope this package will provide easy access to OSM data for
reproducible research in the public interest, adhering to the condition
of the [OdBL licence](https://opendatacommons.org/licenses/odbl/) which
states that

> Any Derivative Database that You Publicly Use must be only under the
> terms of:

-   1.  This License;

-   2.  A later version of this License similar in spirit to this

See the [Introducing osmextract
vignette](https://docs.ropensci.org/osmextract/articles/osmextract.html)
for more details.

## Other approaches

<!-- todo: add links to other packages -->

-   [osmdata](https://github.com/ropensci/osmdata) is an R package for
    importing small datasets directly from OSM servers
-   [geofabrik](https://cran.r-project.org/package=geofabrik) is an R
    package to download OSM data from
    [Geofabrik](https://download.geofabrik.de/)
-   [pyrosm](https://pyrosm.readthedocs.io/en/latest/) is a Python
    package for reading .pbf files
-   [pydriosm](https://pypi.org/project/pydriosm/) is a Python package
    to download, read and import OSM extracts
-   [osmium](https://pypi.org/project/osmium/) provides python bindings
    for the Libosmium C++ library
-   [OpenStreetMapX.jl](https://github.com/pszufe/OpenStreetMapX.jl) is
    a Julia package for reading and analysing .osm files
-   [PostGIS](https://www.bostongis.com/PrinterFriendly.aspx?content_name=loading_osm_postgis)
    is an established spatial database that works well with large OSM
    datasets
-   Any others? Let us know!

## Contribution

We very much look forward to comments, questions and contributions. If
you have any question or if you want to suggest a new approach, feel
free to create a new discussion in the [github
repository](https://github.com/ropensci/osmextract/discussions). If you
found a bug, or if you want to add a new OSM extracts provider, create a
new issue in the [issue
tracker](https://github.com/ropensci/osmextract/issues) or a new [pull
request](https://github.com/ropensci/osmextract/pulls). We always try to
build the most intuitive user interface and write the most informative
error messages, but if you think that something is not clear and could
have been explained better, please let us know.

## Contributor Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

<!-- :) -->
<!-- :) -->
<!-- :) -->
# osmextract (development version)

### MINOR CHANGES
* The `boundary` argument can be specified using `bbox` objects. The `bbox` object is converted to `sfc` object with `sf::st_as_sfc` and preserves the same CRS. 
* Added a more informative error message when `oe_get()` or `oe_read()` are run with empty or unnamed arguments in `...` (#234).

### DOCUMENTATION FIXES
* Update description for `boundary` and `boundary_type` arguments. 

# osmextract 0.4.0 

### MAJOR CHANGES

* Import two new packages: [*httr*](https://cran.r-project.org/package=httr) and [jsonlite](https://cran.r-project.org/package=jsonlite) (#231, #232). 
* Improved the approach adopted to download files from the web. In particular, the functions `oe_download()` and `oe_search()` now take advantage of `httr` functionalities. They return informative messages in case of errors (#231, #232). 
* Vignettes and examples do not require internet connection. 

### BUG FIXES

* Fixed a bug in `oe_vectortranslate()` that occurred when reading `multilinestrings` or `other_relations` layers with one or more extra tags (#229). 
* Fixed a bug in `oe_get()`/`oe_read()` that could return a warning message when reading an existing GPKG file with a `query` argument. 

### MINOR CHANGES

* The duplicated fields in `extra_tags` are now removed before modifying the `osmconf.ini` file. Duplicated tags means something like `extra_tags = c("A", "A")` or even fields that are included by default (i.e. `extra_tags = "highway"` for the `lines` layer). See discussion in #229. 

# osmextract 0.3.1

### MAJOR CHANGES

* Added a new (but still experimental) function named `oe_get_network()` to import a road network used by a specific mode of transport. For the moment, we support the following modes of transport: cycling (default), walking, and driving. Check `?oe_get_network` for more details and examples (#218). 

### MINOR CHANGES

* The `layer` argument is now converted to lower case before checking if the required layer is admissible. 
* Adjusted the code behind `oe_get()` and `oe_vectortranslate()` for `sf` v1.0.2.
* Remove the call to `suppressMessages()` in `oe_match()` (#217).

### DOCUMENTATION FIXES

* Slightly changed the description of the package. 
* Added a `.Rd` file documenting the whole package. 
* Slightly changed the description of parameter `place`. 

# osmextract 0.3.0

### MAJOR CHANGES

* The `oe_get_keys()` function can be used to extract the values associated with all or some keys. We also defined an ad-hoc printing method and fixed several bugs. The examples were improved. Moreover, the function tries to match an input `zone` with one of the OSM extracts previously downloaded (#201 and #196). 
* If the parameter `place` represents an `sf`/`sfc`/`bbox` object with missing CRS, then `oe_match()` raises a warning message and sets `CRS = 4326`. This has relevant consequences on other functions (like `oe_get()`) that wrap `oe_match()`. 
* Starting from `sf` > 0.9.8, the function `oe_vectortranslate()` stops with an error when there is a problem in the argument `vectortranslate_options` and `quiet = FALSE` (instead of raising a warning or crashing the `R` session). See [here](https://github.com/r-spatial/sf/issues/1680) for more details. 
* The options `c("-f", "GPKG", "-overwrite", "-oo", "CONFIG_FILE=", path-to-config-file, "-lco", "GEOMETRY_NAME=geometry", layer)` are always appended at the end of `vectortranslate_options` argument unless the user explicitly sets different default parameters for the arguments `-f`, `-oo` and `-lco` (#200). We believe those are sensible defaults and can help users creating less verbose specifications for `ogr2ogr` utility. 
* We create two new arguments in `oe_vectortranslate()` (therefore also in `oe_get()` and `oe_read()`) named `boundary` and `boundary_type`. They can be used to create an ad-hoc spatial filter during the vectortranslate operations (and create even less verbose specifications in `vectortranslate_options` argument). See docs and introductory vignette for more details. 
* The argument `provider` was removed from `oe_match_pattern()` since the function automatically checks all available providers (#208). 

### BUG FIXES

* The parameter `force_vectortranslate` is checked before reading the layers of an existing `gpkg` file. If `force_vectortranslate` is `TRUE`, then `oe_vectortranslate()` doesn't check the existing layers. This is important for user that run `oe_vectortranslate()` after stopping the vectortranslate process.  
* The arguments `extra_tags` and `osmconf_ini` are not ignored when `vectortranslate_options` is not `NULL` (#182). 
* Fix the provider's data objects for `sf` v1.0 (#194). 

### MINOR IMPROVEMENTS

* The arguments passed to `oe_read()` via `...` are compared with the formals of `st_read.character`, `st_as_sf.data.frame`, and `read_sf`.  
* Added a new method to `oe_match` for `bbox` objects (#185).
* The `oe_get_keys()` function can be applied to `.osm.pbf` objects (#188). 

### DOCUMENTATION FIXES

* Improved several examples and fixed a small bug in the documentation of `oe_match()`.
* Fix several typos in the vignettes and docs. 

### OTHERS

* Created a new space in the github repo named _Discussion_ to have conversations, ask questions and post answers without opening issues. Link: https://github.com/ropensci/osmextract/discussions.
* Tests that require an internet connection are now skipped on CRAN (#189). 

# osmextract 0.2.1

This is a minor release. 

* We modified several examples and tests to fix several errors noticed during CRAN tests (#175). 

# osmextract 0.2.0

Published on CRAN! 

### NEW FEATURES

* Add a `level` parameter to `oe_match()`. It is used to choose between multiple hierarchically nested OSM extracts. The default behaviour is to select the smallest administrative unit (#160).
* Modify the behaviour of `oe_match()`. The function checks all implemented providers in case the input `place` is not matched with any geographical zone for the chosen provider (#155).
* Add a simple interface to Nominatim API that enables `oe_match()` to geolocate text strings that cannot be found in the providers (#155). 

### MINOR IMPROVEMENTS

* Normalise the paths managed by `oe_download_directory` and `oe_download` (#150 and #161) 
* `oe_get_keys` returns an informative error when there is no other_tags field in the input file (#158)

### BUG FIXES

* Fix the structure of `geofabrik_zones` object (#167)
* Fix warning messages related to ... in `oe_get()` (#152)

### DOCUMENTATION FIXES

* Simplify several warning messages in case of spatial matching
* Simplify startup message (#156)
* Add more details related to download timeouts (#145)
* Documented values returned by `oe_find()` and `oe_search()`

# osmextract 0.1.0

* Finish development of the main functions
* Submit to rOpenSci for peer-review
## Release summary

- We decided two import two new packages (httr and jsonlite) and changed the approach for downloading resources from the web. We also fixed some bugs. Details are listed in the NEWS file.  

## Test environments

- local Windows, R 4.0.5
- Github Actions: ubuntu-latest (r-release)
- debian-gcc-devel, fedora-gcc-devel, and macos-highsierra-release-cran via rhub

## R CMD check results

0 errors √ | 0 warnings √ | 0 note √
# oe_get_keys + values: printing method

    Found 38 unique keys, printed in ascending order of % NA values. The first 10 keys are: 
    surface (91% NAs) = {#asphalt = 12; #paved = 3; #cobblestone = 1; #paving_sto...}
    lanes (92% NAs) = {#2 = 9; #1 = 7}
    bicycle (92% NAs) = {#yes = 10; #designated = 5}
    lit (92% NAs) = {#yes = 15}
    access (93% NAs) = {#permissive = 12; #yes = 2}
    oneway (93% NAs) = {#yes = 13}
    maxspeed (94% NAs) = {#30 mph = 12}
    ref (95% NAs) = {#A660 = 9; #4184 = 1}
    foot (95% NAs) = {#yes = 5; #designated = 4}
    natural (96% NAs) = {#tree_row = 7}
    [Truncated output...]

---

    Found 38 unique keys, printed in ascending order of % NA values. The first 10 keys are: 
    surface (91% NAs) = {#asphalt = 12; #paved = 3; #cobblestone = 1; #paving_sto...}
    lanes (92% NAs) = {#2 = 9; #1 = 7}
    bicycle (92% NAs) = {#yes = 10; #designated = 5}
    lit (92% NAs) = {#yes = 15}
    access (93% NAs) = {#permissive = 12; #yes = 2}
    oneway (93% NAs) = {#yes = 13}
    maxspeed (94% NAs) = {#30 mph = 12}
    ref (95% NAs) = {#A660 = 9; #4184 = 1}
    foot (95% NAs) = {#yes = 5; #designated = 4}
    natural (96% NAs) = {#tree_row = 7}
    [Truncated output...]

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
*  Look at the Gihub Actions build status before and after making changes.
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR. The only difference is that we adopt the `=` operator for assignment 
instead of `<-`. 
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the osmextract project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project you agree to abide by its terms.

See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html) for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Describe the steps used to reproduce the bug. Please note that including a reproducible example (consider using a ["reprex"](https://cran.rstudio.com/web/packages/reprex/)) can be extremely helpful for diagnosing the problem and solving the bug. 

**Expected behaviour**
A clear and concise description of what you expected to happen.

**Additional context**
Add any other context about the problem here. In particular, along with your query, please paste your `devtools::session_info()` or `sessionInfo()` into the code block below. 

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
name: Feature request
about: Suggest an idea for this project
title: "[FEATURE]"
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
