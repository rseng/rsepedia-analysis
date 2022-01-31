
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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# osmextract <a href='https://docs.ropensci.org/osmextract/'><img src='man/figures/logo.svg' align="right" height=275/></a>

<!-- badges: start -->

[![R build status](https://github.com/ropensci/osmextract/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/osmextract/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/osmextract/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/osmextract?branch=master)
[![peer-review](https://badges.ropensci.org/395_status.svg)](https://github.com/ropensci/software-review/issues/395)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN status](https://www.r-pkg.org/badges/version/osmextract)](https://CRAN.R-project.org/package=osmextract)
<!-- badges: end -->

The goal of `osmextract` is to make it easier for people to access OpenStreetMap (OSM) data for reproducible research.
OSM data is the premier source of freely available, community created geographic data worldwide.
We aim to enable you to extract it for data-driven work in the public interest.

`osmextract` matches, downloads, converts and imports bulk OSM data hosted by providers such as [Geofabrik GmbH](http://download.geofabrik.de) and [bbbike](https://download.bbbike.org/osm/).
For information on alternative providers and how to add them see the [providers vignette](https://docs.ropensci.org/osmextract/articles/providers.html).

## Why osmextract?

The package answers a common question for researchers who use OSM data:
how to get it into a statistical environment, in an appropriate format, as part of a computationally efficient and reproducible workflow?
Other packages answer parts of this question.
[`osmdata`](https://github.com/ropensci/osmdata), for example, is an R package that provides an R interface to the [Overpass API](https://wiki.openstreetmap.org/wiki/Overpass_API), which is ideal for downloading small OSM datasets.
However, the API is rate limited, making it hard to download large datasets.
As a case study, try to download all cycleways in England using `osmdata`:

```r
library(osmdata)
cycleways_england = opq("England") %>% 
  add_osm_feature(key = "highway", value = "cycleway") %>% 
  osmdata_sf()
# Error in check_for_error(doc) : General overpass server error; returned:
# The data included in this document is from www.openstreetmap.org. The data is made available under ODbL. runtime error: Query timed out in "query" at line 4 after 26 seconds. 
```

The query stops with an error message after around 30 seconds.
The same query can be made with `osmextract` as follows, which reads-in almost 100k linestrings in less than 10 seconds, after the data has been downloaded in the compressed `.pbf` format and converted to the open standard `.gpkg` format. 
The download-and-conversion operation of the OSM extract associated to England takes approximately a few minutes, but this operation must be executed only once. 
The following code chunk is not evaluated.

```{r, eval = FALSE}
library(osmextract)

cycleways_england = oe_get(
  "England",
  quiet = FALSE,
  query = "SELECT * FROM 'lines' WHERE highway = 'cycleway'"
)
par(mar = rep(0.1, 4))
plot(sf::st_geometry(cycleways_england))
```

```{r, echo = FALSE, out.width="80%", fig.align='center'}
knitr::include_graphics("man/figures/89990770-22bf2480-dc83-11ea-9092-764594534959.png")
```

The package is designed to complement `osmdata`, which has advantages over `osmextract` for small datasets: `osmdata` is likely to be quicker for datasets less than a few MB in size, provides up-to-date data and has an intuitive interface. 
`osmdata` can provide data in a range of formats, while `osmextract` only returns [`sf`](https://github.com/r-spatial/sf) objects. 
`osmextract`'s niche is that it provides a fast way to download large OSM datasets in the highly compressed `pbf` format and read them in via the fast C library [GDAL](https://gdal.org/drivers/vector/osm.html) and the popular R package for working with geographic data [`sf`](https://github.com/r-spatial/sf).

## Installation

You can install the released version of `osmextract` from [CRAN](https://cran.r-project.org/package=osmextract) with:

``` r
install.packages("osmextract")
```

You can install the development version from [GitHub](https://github.com/ropensci/osmextract) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/osmextract")
```

Load the package with:

```{r}
library(osmextract)
```

To use alongside functionality in the `sf` package, we also recommend attaching this geographic data package as follows:

```{r}
library(sf)
```

### Warnings: 

The functions defined in this package may return a warning message like 

```
st_crs<- : replacing crs does not reproject data; use st_transform for that 
```

if the user is running an old version of GDAL (<= 3.0.0) or PROJ (<= 6.0.0). 
See [here](https://github.com/r-spatial/sf/issues/1419) for more details. 
Nevertheless, every function should still work correctly. 
Please, raise [a new issue](https://github.com/ropensci/osmextract/issues) if you find any odd behaviour. 

## Basic usage

Give `osmextract` a place name and it will try to find it in a list of names in the specified provider ([Geofabrik](https://www.geofabrik.de/data/download.html) by default).
If the name you give it matches a place, it will download and import the associated data into R.
The function `oe_get()` downloads (if not already downloaded) and reads-in data from OSM providers as `sf` objects.
By default `oe_get()` imports the `lines` layer, but any layer can be read-in by changing the `layer` argument:

```{r points-lines-iow, fig.show = 'hold', out.width = "50%"}
osm_lines = oe_get("Isle of Wight", stringsAsFactors = FALSE, quiet = TRUE)
osm_points = oe_get("Isle of Wight", layer = "points", stringsAsFactors = FALSE, quiet = TRUE)
nrow(osm_lines)
nrow(osm_points)
par(mar = rep(0, 4))
plot(st_geometry(osm_lines), xlim = c(-1.59, -1.1), ylim = c(50.5, 50.8))
plot(st_geometry(osm_points), xlim = c(-1.59, -1.1), ylim = c(50.5, 50.8))
```

The figures above give an insight into the volume and richness of data contained in OSM extracts.
Even for a small island such as the Isle of Wight, it contains over 50k features including ferry routes, shops and roads.
The column names in the `osm_lines` object are as follows:

```{r}
names(osm_lines) # default variable names
```

Once imported, you can use all functions for data frames in base R and other packages.
You can also use functions from the `sf` package for spatial analysis and visualisation.
Let's plot all the major, secondary and residential roads, for example:

```{r iow1}
ht = c("primary", "secondary", "tertiary", "unclassified") # highway types of interest
osm_major_roads = osm_lines[osm_lines$highway %in% ht, ]
plot(osm_major_roads["highway"], key.pos = 1)
```

The same steps can be used to get other OSM datasets (examples not run):

```{r, eval = FALSE}
malta = oe_get("Malta", quiet = TRUE)
andorra = oe_get("Andorra", extra_tags = "ref")
leeds = oe_get("Leeds")
goa = oe_get("Goa", query = "SELECT highway, geometry FROM 'lines'")
```

If the input place does not match any of the existing names in the supported providers, then `oe_get()` will try to geocode it via [Nominatim API](https://nominatim.org/release-docs/develop/api/Overview/), and it will select the smallest OSM extract intersecting the area. 
For example (not run): 

```{r, eval = FALSE}
oe_get("Milan") # Warning: It will download more than 400MB of data
#> No exact match found for place = Milan and provider = geofabrik. Best match is Iran.
#> Checking the other providers.
#> No exact match found in any OSM provider data. Searching for the location online.
#> ... (extra messages here)
```

For further details on using the package, see the [Introducing osmextract vignette](https://docs.ropensci.org/osmextract/articles/osmextract.html).

## Persistent download directory

The default behaviour of `oe_get()` is to save all the files in a temporary directory, which is erased every time you restart your R session. 
If you want to set a directory that will persist, you can add `OSMEXT_DOWNLOAD_DIRECTORY=/path/for/osm/data` in your `.Renviron` file, e.g. with:

```{r, eval = FALSE}
usethis::edit_r_environ()
# Add a line containing: OSMEXT_DOWNLOAD_DIRECTORY=/path/to/save/files
```

We strongly advise you setting a persistent directory since working with `.pbf` files is an expensive operation, that is skipped by `oe_*()` functions if they detect that the input `.pbf` file was already downloaded.

You can always check the default `download_directory` used by `oe_get()` with: 

```{r, eval = FALSE}
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

We would love to see more providers added (see the [Add new OpenStreetMap providers](https://docs.ropensci.org/osmextract/articles/providers.html) for details) and see what people can do with OSM datasets of the type provided by this package in a reproducible and open statistical programming environment for the greater good.
Any contributions to support this or any other improvements to the package are very welcome via our issue tracker.

## Licence

We hope this package will provide easy access to OSM data for reproducible research in the public interest, adhering to the condition of the [OdBL licence](https://opendatacommons.org/licenses/odbl/) which states that

> Any Derivative Database that You Publicly Use must be only under the terms of:

- i. This License;
- ii. A later version of this License similar in spirit to this

See the [Introducing osmextract vignette](https://docs.ropensci.org/osmextract/articles/osmextract.html) for more details.

## Other approaches

<!-- todo: add links to other packages -->
- [osmdata](https://github.com/ropensci/osmdata) is an R package for importing small datasets directly from OSM servers
- [geofabrik](https://cran.r-project.org/package=geofabrik) is an R package to download OSM data from [Geofabrik](https://download.geofabrik.de/)
- [pyrosm](https://pyrosm.readthedocs.io/en/latest/) is a Python package for reading .pbf files
- [pydriosm](https://pypi.org/project/pydriosm/) is a Python package to download, read and import OSM extracts
- [osmium](https://pypi.org/project/osmium/) provides python bindings for the Libosmium C++ library
- [OpenStreetMapX.jl](https://github.com/pszufe/OpenStreetMapX.jl) is a Julia package for reading and analysing .osm files
- [PostGIS](https://www.bostongis.com/PrinterFriendly.aspx?content_name=loading_osm_postgis) is an established spatial database that works well with large OSM datasets
- Any others? Let us know!

## Contribution

We very much look forward to comments, questions and contributions. 
If you have any question or if you want to suggest a new approach, feel free to create a new discussion in the [github repository](https://github.com/ropensci/osmextract/discussions). 
If you found a bug, or if you want to add a new OSM extracts provider, create a new issue in the [issue tracker](https://github.com/ropensci/osmextract/issues) or a new [pull request](https://github.com/ropensci/osmextract/pulls). 
We always try to build the most intuitive user interface and write the most informative error messages, but if you think that something is not clear and could have been explained better, please let us know. 

## Contributor Code of Conduct
Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.

<!-- :) -->
<!-- :) -->
<!-- :) -->
<!-- :) -->
---
title: "Comparing the supported OSM providers"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Comparing the supported OSM providers}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>", 
  fig.align = "center"
)

# save user's pars
user_par = par(no.readonly = TRUE)
```

This vignette presents a simple comparison between the OSM providers supported by `osmextract`, explaining their pros and cons. 
We decided to write this vignette since, as you will see in the following examples, even if you always start from the same pre-defined `place`, you can get significantly different OSM extracts according to the chosen `provider`. 
Hence, we want to help you choose the best suitable provider for a given situation. 

We assume that you are already familiar with the basic functions in `osmextract`, otherwise please check the "Get Started" vignette for a more detailed introduction. 
Now, let's start with an example, but, first of all, we have to load the package: 

```{r}
library(osmextract)
library(sf)
```

We geocode the coordinates of Lima, the Capital of Peru, 

```{r, eval = FALSE}
lima = tmaptools::geocode_OSM("Lima, Peru")$coords
```

```{r, echo = FALSE}
lima = c(-77.0365256, -12.0621065)
```

and look for a match in the OSM extracts using `oe_match()`: 

```{r, warning = FALSE, message = FALSE}
oe_match(lima, provider = "geofabrik")
oe_match(lima, provider = "bbbike")
oe_match(lima, provider = "openstreetmap_fr")
```

We can see that: 

* when we used `geofabrik` provider (which is also the default provider), then the input `place` was matched with an OSM extract corresponding to Peru region; 
* when we used the `bbbike` provider, then the input `place` was matched with an OSM extract corresponding to the city of Lima; 
* when we used `openstreetmap_fr` provider, then the input data was matched with an OSM extract covering the whole of South America. 

The reason behind these differences is that each OSM provider divides the geographical space into different discrete chunks, and, in the following paragraphs, we will show the tessellation used by each provider. 

## Geofabrik

`geofabrik` is a society that provides map-based services and free downloads of OSM extracts that are updated daily. 
These extracts are based on a division of the world into different regions, covering a whole continent (plus Russian Federation): 

```{r}
par(mar = rep(0, 4))
plot(geofabrik_zones[geofabrik_zones$level == 1, "name"], key.pos = NULL, main = NULL)
```

or several countries all around the world: 

```{r}
plot(geofabrik_zones[geofabrik_zones$level == 2, "name"], key.pos = NULL, main = NULL)
```

Geofabrik also defines several special zones, such as Alps, Britain and Ireland, Germany, Austria and Switzerland, US Midwest, US Northeast, US Pacific, US South and US West. 
Moreover, it contains extracts relative to some administrative subregions, mainly in Europe, Russia, Canada and South America: 

```{r}
plot(geofabrik_zones[geofabrik_zones$level == 3, "name"], key.pos = NULL, main = NULL)
```

Check `?geofabrik_zones` and the [provider's webpage](http://download.geofabrik.de/) for more details. 

## Openstreetmap.fr

`openstreetmap_fr` extracts are taken from http://download.openstreetmap.fr/, a web-service that provides OSM data updated every few minutes. 
The extracts are based on several regions, such as the continents: 

```{r}
# Russian federation is considered as a level 1 zone
plot(openstreetmap_fr_zones[openstreetmap_fr_zones$level == 1, "name"], key.pos = NULL, main = NULL)
```

or some countries around the world (less than `geofabrik`): 

```{r}
plot(openstreetmap_fr_zones[openstreetmap_fr_zones$level == 2, "name"], key.pos = NULL, main = NULL)
```

It can be noticed that there are several holes (such as Peru, which is the reason why, in the first example, Lima was matched with South America data), implying that `openstreetmap_fr` cannot always be used for geographical matching of a `place`. 
Nevertheless, it provides extremely detailed extracts for some regions of the world, like China, 

```{r}
plot(openstreetmap_fr_zones[openstreetmap_fr_zones$parent == "china", "name"], key.pos = NULL, main = NULL)
```

India, 

```{r}
plot(openstreetmap_fr_zones[openstreetmap_fr_zones$parent == "india", "name"], key.pos = NULL, main = NULL)
```

France, 

```{r}
ids_2 = openstreetmap_fr_zones$parent %in% "france"
ids_3 = openstreetmap_fr_zones$parent %in% openstreetmap_fr_zones$id[ids_2]

plot(openstreetmap_fr_zones[ids_2 | ids_3, "name"], key.pos = NULL, main = NULL)
```

and Brazil

```{r}
ids_2 = openstreetmap_fr_zones$parent %in% "brazil"
ids_3 = openstreetmap_fr_zones$parent %in% openstreetmap_fr_zones$id[ids_2]

plot(openstreetmap_fr_zones[ids_2 | ids_3, "name"], key.pos = NULL, main = NULL)
```

## BBBike

`bbbike` provider is based on https://download.bbbike.org/osm/bbbike/. 
It is quite different from any other provider supported in `osmextract` since it contains OSM data for more than 200 cities worldwide. 

```{r, eval = FALSE}
par(mar = rep(0, 4))
plot(sf::st_geometry(spData::world))
plot(sf::st_geometry(bbbike_zones), border = "darkred", add = TRUE, lwd = 3)
```

```{r, echo = FALSE, out.width="100%"}
knitr::include_graphics(
  "../man/figures/96640949-3f7f7480-1324-11eb-9dca-a971c8103a4e.png"
)
```

`bbbike` provider is the safest choice if you are looking for OSM data relative to a particular city in the world.

```{r, include=FALSE}
par(user_par)
```

---
title: "Introducing osmextract"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introducing osmextract}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.align = "center"
)

# save user's options and pars
user_options = options()
user_par = par(no.readonly = TRUE)

# set new options
options(width = 100)
```

This vignette provides an introduction to using the package, building on the README which covers installation and our motivations for creating it.

Loading the package generates important messages about the license associated with OSM data.

```{r}
library(osmextract)
```

The first thing to say is: **do not ignore this message**!
The Open Street Map (OSM) extracts are stored by external providers such as [Geofabrik](https://download.geofabrik.de/), [Bbbike](https://download.bbbike.org/osm/), or [OpenStreetMap.fr](http://download.openstreetmap.fr/). 
There are important legal considerations that you should be aware of before using OSM data, especially if you are working in a for-profit capacity.

# Legal considerations

Anyone using OSM data is bound by law to adhere to the [ODbL](https://opendatacommons.org/licenses/odbl/summary/), which means that you must:

- **Attribute**: You must attribute any public use of the database, or works produced from the database, in the manner specified in the ODbL. For any use or redistribution of the database, or works produced from it, you must make clear to others the license of the database and keep intact any notices on the original database.
- **Share-Alike**: If you publicly use any adapted version of this database, or works produced from an adapted database, you must also offer that adapted database under the ODbL.
- **Keep open**: If you redistribute the database or an adapted version of it, then you may use technological measures that restrict the work (such as DRM) as long as you also redistribute a version without such measures.

In short, publicly using OSM data without attribution or selling datasets derived from it is illegal.
See the [License/Use Cases page on the OSM wiki](https://wiki.openstreetmap.org/wiki/License/Use_Cases) for detailed use cases.

# Main package functions

The package is composed of the following main functions: 

1. `oe_providers()`: Show which OSM providers are available;
1. `oe_match()`: Match an input place with one of the files stored by the OSM providers;
1. `oe_download()`: Download the chosen file;
1. `oe_vectortranslate()`: Convert between `.pbf` and `.gpkg` formats;
1. `oe_read()`: Read `.pbf` and `.gpkg` files;
1. `oe_get()`: Match, download, (vector)translate, and import data, all in one step.

For many users who just want to get OSM data quickly, `oe_get()` may be sufficient, as covered in the README.
We will demonstrate each function in turn, following the same order in which they are typically used. 
As you can see, the name of the most important functions in this package start with `oe_*` prefix, which means that you can easily use auto-completion features (with Rstudio or similar IDE(s)). 

## `oe_providers()`: List providers

`oe_providers()` lists the providers that are currently available with the version of `osmextract` you have installed.

```{r}
oe_providers()
```

Each element in the column `database_name` is a data object that is packaged with `osmextract`. 
You can read a detailed description of each provider data running, for example, `?geofabrik_zones` or `?bbbike_zones`. 

Perhaps, the best known bulk OSM data provider is [Geofabrik](https://www.geofabrik.de/), and its extracts are summarised as a `data.frame` in the packaged object `geofabrik_zones`.

```{r}
class(geofabrik_zones)
```

Note that in addition to being a data frame with rows and columns, `geofabrik_zones` is also an `sf` object, as defined in the [package](https://r-spatial.github.io/sf/) of the same name.
When working with `sf` objects, it makes sense to have the package loaded:

```{r}
library(sf)
```

That gives you access to many functions for working with geographic vector data of the type provided by `osmextract`.
Each row of data in an `sf` object contains a geometry, representing the area covered by a provider zone, meaning you can plot the data as follows:

```{r}
par(mar = rep(0.1, 4))
plot(st_geometry(geofabrik_zones))
```

The plot above shows how the provider divides geographic space into discrete chunks.
Different providers have other zoning systems. 
For example: 

```{r, eval = FALSE}
par(mar = rep(0.1, 4))
plot(st_geometry(spData::world), xlim = c(-2, 10), ylim = c(35, 60))
plot(st_geometry(bbbike_zones), xlim = c(-2, 10), ylim = c(35, 60), col = "darkred", add = TRUE)
```

```{r, echo = FALSE, out.width="80%"}
knitr::include_graphics(
  path = "../man/figures/94461461-772e4d00-01ba-11eb-950c-804ad177729f.png"
)
```

As shown in the visualisation above of [BBBike.org](https://download.bbbike.org/osm/) zones in Europe, that provider offers rectangular extracts of the major cities.
We are working on adding support for manually selected regions from the BBBike website (see https://github.com/ropensci/osmextract/issues/100).

Check the ["Comparing the supported OSM providers"](https://docs.ropensci.org/osmextract/articles/providers_comparisons.html) vignette for some simple guidelines on how to choose the best provider. 

## `oe_match()`: Match an input place with an OSM extract

The function `oe_match()` takes in input a string through the parameter `place`, and it returns a named list of length two with the URL and the size (in bytes) of a `.osm.pbf`^[The `.pbf` format is a highly optimised binary format used by OSM providers to store and share OSM extracts.] file representing a geographical zone stored by one of the supported providers. 
For example: 

```{r}
oe_match("Italy")
oe_match("Leeds", provider = "bbbike")
```

The geographical zone is chosen by calculating the Approximate String Distance (`?adist()`) between the input `place` and one of the fields in the provider’s dataset.
Then, the function selects the closest match. 
By default, `oe_match()` uses the `name` field and `Geofabrik` provider, but you can select a different field via the argument `match_by`. 
We refer to the providers' help pages for a detailed description of all available fields. 
If you are using Geofabrik provider, a useful and interesting alternative field is represented by the (unique and unambiguous) [`iso3166-1 alpha2` codes](https://it.wikipedia.org/wiki/ISO_3166-1_alpha-2): 

```{r}
oe_match("RU", match_by = "iso3166_1_alpha2")
oe_match("US", match_by = "iso3166_1_alpha2")
```

There are a few scenarios where the `iso3166-1 alpha2` codes in `geofabrik_data` cannot be used since there are no per-country extracts (e.g. Israel and Palestine):

```{r, error = TRUE}
oe_match("PS", match_by = "iso3166_1_alpha2", quiet = TRUE)
oe_match("IL", match_by = "iso3166_1_alpha2", quiet = TRUE)
```

For this reason, we coded a function named `oe_match_pattern()` to explore the matching operations for all available providers according to a pre-defined pattern. 
It returns a named list where the names are the id(s) of the supported OSM providers and the values are the matched names. 
For example:

```{r}
oe_match_pattern("London")
oe_match_pattern("Yorkshire")
oe_match_pattern("Russia")
oe_match_pattern("Palestine")
```

The default field is `name`, but we can change that as follows: 

```{r}
oe_match_pattern("US", match_by = "iso3166_2")
```

If we set `full_row = TRUE`, then `oe_match_pattern()` will return the complete row(s) from each provider's data: 

```{r}
lapply(oe_match_pattern("Israel", full_row = TRUE), function(x) x[, 1:3])
```

We can combine the two functions as follows: 

```{r}
oe_match_pattern("Valencia")
oe_match("Comunitat Valenciana", provider = "openstreetmap_fr")
```

The parameter `max_string_dist` (default value is 1) represents the maximum tolerable distance between the input place and the closest match in `match_by` column. 
This value can always be increased to help the matching operations, but that can lead to false matches:

```{r, error = TRUE}
# erroneous match
oe_match("Milan", max_string_dist = 2)
```

The parameter `max_string_dist` is always set to 0 if `match_by` argument is equal to `iso3166_1_alpha2` or `iso3166_2`. 

If the approximate string distance between the closest match and the input `place` is greater than `max_string_dist`, then `oe_match()` will also check the other supported providers. 
For example: 

```{r}
oe_match("Leeds")
oe_match("London")
oe_match("Vatican City")
```

Finally, if there is no exact match with any of the supported providers and `match_by` argument is equal to `"name"`, then `oe_match()` will use the [Nominatim API](https://nominatim.org/release-docs/develop/api/Overview/) to geolocate the input place and perform a spatial matching operation (explained below): 

```{r, eval = FALSE}
oe_match("Milan")
#> No exact match found for place = Milan and provider = geofabrik. Best match is Iran. 
#> Checking the other providers.
#> No exact match found in any OSM provider data. Searching for the location online.
#> The input place was matched with Nord-Ovest. 
#> $url
#> [1] "https://download.geofabrik.de/europe/italy/nord-ovest-latest.osm.pbf"
#> $file_size
#> [1] 416306623
```

### Finding zones based on geographic inputs

The input `place` can also be specified using an `sf`, `sfc`, or `bbox` object with arbitrary CRS^[If the input spatial object has no CRS, then `oe_match()` raises a warning message and sets `CRS = 4326`.], as documented in the following example. 
`oe_match()` will return a named list of length two with the URL and the size of a `.pbf` file representing a zone that geographically contains the `sf` or `sfc` object  (or an error if the input is not contained into any geographical area).

```{r}
milan_duomo = sf::st_sfc(sf::st_point(c(1514924, 5034552)), crs = 3003)
oe_match(milan_duomo)
```

If the input `place` intersects multiple geographically nested areas and the argument `level` is equal to `NULL` (the default value), then `oe_match()` automatically returns the extract with the highest `level`. 
In particular, we could roughly say that smaller geographical areas are associated with higher `level`(s). 
For example, `level = 1` may correspond to continent-size extracts, `2` is for countries, `3` represents the regions and `4` the subregions:

```{r, error = TRUE}
yak = c(-120.51084, 46.60156)
oe_match(yak, level = 1, quiet = TRUE)
oe_match(yak, level = 2, quiet = TRUE) # the default
oe_match(yak, level = 3, quiet = TRUE) # error
```

If there are multiple OSM extract intersecting the input `place` at the same `level`, then `oe_match()` will return the area whose centroid is closest to the input `place`.

If you specify more than one geometry into the `sf` or `sfc` object, then `oe_match()` will select an area that contains all of them. 

```{r}
milan_leeds = st_sfc(
  st_point(c(9.190544, 45.46416)), # Milan
  st_point(c(-1.543789, 53.7974)), # Leeds
  crs = 4326
)
oe_match(milan_leeds)
```

The same operations work with `LINESTRING` or `POLYGON` objects: 

```{r}
milan_leeds_linestring = st_sfc(
  st_linestring(
    rbind(c(9.190544, 45.46416), c(-1.543789, 53.7974))
  ), 
  crs = 4326
)
oe_match(milan_leeds_linestring)
```

The input `place` can also be specified using a numeric vector of coordinates. 
In that case, the CRS is assumed to be [EPSG:4326](https://spatialreference.org/ref/epsg/4326/):

```{r}
oe_match(c(9.1916, 45.4650)) # Duomo di Milano using EPSG: 4326
```

Finally, to reduce unnecessary computational resources and save bandwidth/electricity, we will use a small OSM extract in subsequent sections that can be matched as follows:

```{r}
# ITS stands for Institute for Transport Studies: https://environment.leeds.ac.uk/transport
(its_details = oe_match("ITS Leeds"))
```

## `oe_download()`: Download OSM extracts 

The `oe_download()` function is used to download `.pbf` files representing OSM extracts.
It takes in input a URL, through the parameter `file_url`, and it downloads the requested data in a directory (specified by the parameter `download_directory`):

```{r, eval = FALSE}
oe_download(
  file_url = its_details$url, 
  file_size = its_details$file_size,
  provider = "test",
  download_directory = # path-to-a-directory
)
```

The argument `provider` can be omitted if the input `file_url` is associated with one of the supported providers.
The default value for `download_directory` is `tempdir()` (see `?tempdir`), but, if you want to point to a directory that will persist, you can add `OSMEXT_DOWNLOAD_DIRECTORY=/path/for/osm/data` in your `.Renviron` file, e.g. with:

```{r, eval = FALSE}
usethis::edit_r_environ()
# Add a line containing: OSMEXT_DOWNLOAD_DIRECTORY=/path/for/osm/data
```

You can always check the default `download_directory` used by `oe_download()` with: 

```{r}
oe_download_directory()
```

We strongly advise you setting a persistent directory since downloading and converting (see the next sub-section) `.pbf` files are expensive operations, that can be skipped if the functions detect that the requested extract was already downloaded and/or converted.

More precisely, `oe_download()` runs several checks before actually downloading a new file, to avoid overloading the OSM providers. 
The first step is the definition of the path associated with the input `file_url`. 
The path is created by pasting together the `download_directory`, the name of the chosen provider (specified by `provider` argument or inferred from the input URL), and the `basename()` of the URL. 
For example, if `file_url` is equal to `"https://download.geofabrik.de/europe/italy-latest.osm.pbf"`, and `download_directory = "/tmp/`, then the path is built as `/tmp/geofabrik_italy-latest.osm.pbf`. 
In the second step, the function checks if the new path/file already exists (using `file.exists()`) and, in that case, it returns the path, without downloading anything^[The parameter `force_download` can be used to override this behaviour in case you need to update an old OSM extract.]. 
Otherwise, it downloads a new file (using `download.file()` with `mode = "wb"`) and then it returns the path.

## `oe_vectortranslate()`: Convert to gpkg format

The `oe_vectortranslate()` function translates a `.pbf` file into `.gpkg` format^[The GeoPackage (`.gpkg`) is an *open, standards-based, platform-independent, portable, self-descripting, compact format for transferring geospatial information*. See [here](http://www.geopackage.org/).]. 
It takes in input a string representing the path to an existing `.pbf` file, and it returns the path to the newly generated `.gpkg` file. 
The `.gpkg` file is created in the same directory as the input `.pbf` file and with the same name. 
The conversion is performed using [ogr2ogr](https://gdal.org/programs/ogr2ogr.html#ogr2ogr) through `vectortranslate` utility in `sf::gdal_utils()`.

We decided to adopt this approach following [the suggestions](https://github.com/OSGeo/gdal/issues/2100#issuecomment-565707053) of the maintainers of GDAL.
Moreover, GeoPackage files have database capabilities like random access and querying that are extremely important for OSM data (see below). 

Let's start with an example.
First, we download the `.pbf` file associated with ITS example: 

```{r, include=FALSE}
its_pbf = file.path(oe_download_directory(), "test_its-example.osm.pbf")
file.copy(
  from = system.file("its-example.osm.pbf", package = "osmextract"), 
  to = its_pbf, 
  overwrite = TRUE
)
```

```{r, eval = 2}
its_pbf = oe_download(its_details$url, provider = "test", quiet = TRUE) # skipped online, run it locally
list.files(oe_download_directory(), pattern = "pbf|gpkg")
```

and then we convert it to `.gpkg` format: 

```{r}
its_gpkg = oe_vectortranslate(its_pbf)
list.files(oe_download_directory(), pattern = "pbf|gpkg")
```

The vectortranslate operation can be customised in several ways modifying the parameters `layer`, `extra_tags`, `osmconf_ini`, `vectortranslate_options`, `boundary` and `boundary_type`.  

### `layer` argument

The `.pbf` files processed by GDAL are usually categorized into 5 layers, named `points`, `lines`, `multilinestrings`, `multipolygons` and `other_relations`^[Check the first paragraphs [here](https://gdal.org/drivers/vector/osm.html) for more details.]. 
The `oe_vectortranslate()` function can covert only one layer at a time. 
Nevertheless, several layers with different names can be stored in the same `.gpkg` file. 
By default, the function will convert the `lines` layer (which is the most common one according to our experience), but you can change that using the parameter `layer`. 

The `.pbf` files always contain all five layers: 

```{r}
st_layers(its_pbf, do_count = TRUE)
```

while, by default, `oe_vectortranslate` convert only the `lines` layer: 

```{r}
st_layers(its_gpkg, do_count = TRUE)
```

We can add another layer as follows: 

```{r}
its_gpkg = oe_vectortranslate(its_pbf, layer = "points")
st_layers(its_gpkg, do_count = TRUE)
```

### `osmconf_ini` and `extra_tags`

The arguments `osmconf_ini` and `extra_tags` are used to modify how GDAL reads and processes a `.pbf` file. 
More precisely, several operations that GDAL performs on a `.pbf` file are governed by a `CONFIG` file, that you can check at the following [link](https://github.com/OSGeo/gdal/blob/master/data/osmconf.ini). 
The package stores a local copy which is used as the standard `CONFIG` file.

The basic components of OSM data are called [*elements*](https://wiki.openstreetmap.org/wiki/Elements) and they are divided into *nodes*, *ways* or *relations*. 
Hence, for example, the code at line 7 of that `CONFIG` file is used to determine which *ways* are assumed to be *polygons* if they are closed.

The parameter `osmconf_ini` can be used to specify the path to a different `CONFIG` file in case you need more control over GDAL operations. 
See the next sub-sections for an example. 
If `osmconf_ini` is equal to `NULL` (the default), then `oe_vectortranslate()` function uses the standard `CONFIG` file.

Another example can be presented as follows. 
OSM data is usually described using several [*tags*](https://wiki.openstreetmap.org/wiki/Tags), i.e. pairs of two items: a *key* and a *value*.
The code at lines 33, 53, 85, 103, and 121 of the default `CONFIG` file determines, for each layer, which tags are explicitly reported as fields, while all the other tags are stored in the `other_tags` column (see [here](https://gdal.org/drivers/vector/osm.html#other-tags-field) for more details). 
The parameter `extra_tags` (default value: `NULL`) governs which tags are explicitly reported in the `.gpkg` file and are omitted from the `other_tags` field. 
The default tags are always included (unless you modify the `CONFIG` file or the `vectortranslate_options`). 
Please note that the argument `extra_tags` is ignored if `osmconf_ini` is not `NULL` (since we do not know how you generated the new `.ini` file). 

Lastly, the `oe_get_keys()` function can be used to check all `keys` that are stored in the `other_tags` field for a given `.gpkg` or `.pbf` file. 
For example, 

```{r}
oe_get_keys(its_gpkg, layer = "lines")
```

Starting from version `0.3.0`, if you set `values = TRUE`, then `oe_get_keys` returns the values associated to each key (we also defined an ad-hoc printing method): 

```{r}
oe_get_keys(its_gpkg, layer = "lines", values = TRUE)
```

Check `?oe_get_keys` for more details. 

We can always re-create the `.gpkg` file adding one or more new `tags`: 

```{r}
its_gpkg = oe_vectortranslate(its_pbf, extra_tags = c("bicycle", "foot"))
```

Check the next sections for more complex, useful, and realistic use-cases. 

### `vectortranslate_options` argument

The parameter `vectortranslate_options` is used to control the arguments that are passed to `ogr2ogr` via `sf::gdal_utils()` when converting between `.pbf` and `.gpkg` formats. 
The utility `ogr2ogr` can perform various operations during the translation process, such as spatial filters or queries. 
These operations can be tuned using the `vectortranslate_options` argument. 
If `NULL` (default value), then `vectortranslate_options` is set equal to `c("-f", "GPKG", "-overwrite", "-oo", paste0("CONFIG_FILE=", osmconf_ini),  "-lco", "GEOMETRY_NAME=geometry", layer)`. 
Explanation:

* `"-f", "GPKG"` says that the output format is `GPKG`. This is mandatory for GDAL < 2.3;
* `"-overwrite` is used to delete an existing layer and recreate it empty;
* `"-oo", paste0("CONFIG_FILE=", osmconf_ini)` is used to modify the [open options](https://gdal.org/drivers/vector/osm.html#open-options) for the `.osm.pbf` file and set the path of the `CONFIG` file;
* `"-lco", "GEOMETRY_NAME=geometry"` adjust the [layer creation options](https://gdal.org/drivers/vector/gpkg.html?highlight=gpkg#layer-creation-options) for the `.gpkg` file, modifying the name of the geometry column; 
* `layer` indicates which layer should be converted.

Starting from version 0.3.0, the options `c("-f", "GPKG", "-overwrite", "-oo", "CONFIG_FILE=", paste0("CONFIG_FILE=", osmconf_ini), "-lco", "GEOMETRY_NAME=geometry", layer)` are always appended at the end of `vectortranslate_options` unless you explicitly set different default parameters for the arguments `-f`, `-oo` and `-lco`. 

### `boundary` and `boundary_type` arguments 

According to our experience, spatial filters are the most common operations added to the (default) vectortranslate process (usually to select a smaller area lying in a larger OSM extract). 
Hence, starting from version 0.3.0, we defined two new arguments named `boundary` and `boundary_type` that can be used to easily apply a spatial filter directly when converting the compressed OSM extract. 
These new arguments are exemplified in the next sections and can help all users creating less verbose `vectortranslate_options`. 

### Other notes

By default, the vectortranslate operations are skipped if `oe_vectortranslate()` function detects a file having the same path as the input file, `.gpkg` extension and a layer with the same name as the parameter `layer` with all `extra_tags`. 
In that case, the function will return the path of the `.gpkg` file. 
This behaviour can be overwritten by setting `force_vectortranslate = TRUE`. 
If the arguments `osmconf_ini`, `vectortranslate_options` or `boundary` parameters are not `NULL`, the vectortranslate operations are never skipped.

Starting from `sf` version [0.9.6](https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13), if `quiet` argument is equal to `FALSE`, then `oe_vectortranslate()` will display a progress bar during he vectortranslate process.

## `oe_read()`: Read-in OSM data 

The `oe_read()` function is a wrapper around `oe_download()`, `oe_vectortranslate()`, and `sf::st_read()`. 
It is used for reading-in a `.pbf` or `.gpkg` file that is specified using its path or URL. 

So, for example, the following code can be used for reading-in the `its-gpkg` file: 

```{r}
oe_read(its_gpkg)
```

If the input `file_path` points to a `.osm.pbf` file, the vectortranslate operations can be skipped using the parameter `skip_vectortranslate`. 
In that case, `oe_read()` will ignore the conversion step. 

```{r}
oe_read(its_pbf, skip_vectortranslate = TRUE, quiet = FALSE)
```

We can see that the output data includes nine fields (i.e. the default tags), while the previous example had 11 fields (i.e. the default tags + `bicycle` and `foot` tags, that were added to the `.gpkg` file a few chunks above). 

We can also read an object starting from a URL (not evaluated here): 

```{r, eval = FALSE}
my_url = "https://github.com/ropensci/osmextract/raw/master/inst/its-example.osm.pbf"
oe_read(my_url, provider = "test", quiet = TRUE, force_download = TRUE, force_vectortranslate = TRUE)
```

Please note that if you are reading from a URL which is not linked with any of the supported providers, you need to specify the `provider` parameter. 
The `test_its-example.osm.pbf` file already exists in the `download_directory`, but we forced the download and vectortranslate operations. 

## `oe_get()`: Do it all in one step 

To simplify the steps outlined above, while enabling modularity if needs be, we packaged them all into a single function that works as follows:

```{r}
its_lines = oe_get("ITS Leeds")
par(mar = rep(0.1, 4))
plot(its_lines["highway"], lwd = 2, key.pos = NULL)
```

The function `oe_get()` is a wrapper around `oe_match()` and `oe_read()`, and it summarizes the algorithm that we use for importing OSM extracts: 

1. Match the input `place` with the URL of a `.pbf` file through `oe_match()`;
2. If necessary, download the corresponding `.pbf` file using `oe_download()`; 
3. Convert it into `.gpkg` format using `oe_vectortranslate()`; 
4. Read-in one layer of the `.gpkg` file using `sf::st_read()`. 

The following commands (not evaluated here) show how `oe_get()` can be used to import the OSM extracts associated with the desired input `place`, after downloading the `.pbf` file and performing the vectortranslate operations. 
We suggest you run the commands and check the output. 

```{r, eval = FALSE}
oe_get("Andorra")
oe_get("Leeds")
oe_get("Goa")
oe_get("Malta", layer = "points", quiet = FALSE)
oe_match("RU", match_by = "iso3166_1_alpha2", quiet = FALSE)

oe_get("Andorra", download_only = TRUE)
oe_get_keys("Andorra")
oe_get_keys("Andorra", values = TRUE)
oe_get_keys("Andorra", values = TRUE, which_keys = c("oneway", "surface", "maxspeed"))

oe_get("Andorra", extra_tags = c("maxspeed", "oneway", "ref", "junction"), quiet = FALSE)
oe_get("Andora", stringsAsFactors = FALSE, quiet = TRUE, as_tibble = TRUE) # like read_sf

# Geocode the capital of Goa, India
(geocode_panaji = tmaptools::geocode_OSM("Panaji, India"))
oe_get(geocode_panaji$coords, quiet = FALSE) # Large file
oe_get(geocode_panaji$coords, provider = "bbbike", quiet = FALSE)
oe_get(geocode_panaji$coords, provider = "openstreetmap_fr", quiet = FALSE)

# Spatial match starting from the coordinates of Arequipa, Peru
geocode_arequipa = c(-71.537005, -16.398874)
oe_get(geocode_arequipa, quiet = FALSE)
oe_get(geocode_arequipa, provider = "bbbike", quiet = FALSE) # Error
oe_get(geocode_arequipa, provider = "openstreetmap_fr", quiet = FALSE) # No country-specific extract
```

The arguments `osmconf_ini`, `vectortranslate_options`, `boundary`, `boundary_type`, `query` and `wkt_filter` (the last two arguments are defined in `sf::st_read()`) can be used to further optimize the process of getting OSM extracts into R.  

### `osmconf_ini`

The following example shows how to create an ad-hoc `CONFIG` file, which is used by GDAL to read a `.pbf` file in a customised way. 
First, we load a local copy of the default `osmconf.ini` file, taken from the following [link](https://github.com/OSGeo/gdal/blob/master/data/osmconf.ini). 

```{r}
custom_osmconf_ini = readLines(system.file("osmconf.ini", package = "osmextract"))
```

Then, we modify the code at lines 18 and 21 asking GDAL to report all nodes and ways (even without any significant tag).

```{r}
custom_osmconf_ini[[18]] = "report_all_nodes=yes"
custom_osmconf_ini[[21]] = "report_all_ways=yes"
```

We change also the code at lines 45 and 53, removing the `osm_id` field and changing the default attributes: 

```{r}
custom_osmconf_ini[[45]] = "osm_id=no"
custom_osmconf_ini[[53]] = "attributes=highway,lanes"
```

Another relevant parameter that could be customised during the creating of an ad-hoc `osmconf.ini` file is `closed_ways_area_polygons` (see lines 5-7 of the default CONFIG file). 
We can now write a local copy of the `custom_osmconf_ini` file: 

```{r}
temp_ini = tempfile(fileext = ".ini")
writeLines(custom_osmconf_ini, temp_ini)
```

and read the ITS Leeds file with the new `osmconf.ini` file: 

```{r}
oe_get("ITS Leeds", provider = "test", osmconf_ini = temp_ini, quiet = FALSE)
```

If we compare it with the default output: 

```{r}
oe_get("ITS Leeds", provider = "test", quiet = FALSE, force_vectortranslate = TRUE)
```

we can see that there are 2 extra features in the `sf` object that was read-in using the customized `CONFIG` file (i.e. 191 features instead of 189 since we set `"report_all_nodes=yes"` and `"report_all_ways=yes"`) and just 4 field: `highway`, `lanes` (see the code a few chunks above), `z_order` (check the code [here](https://github.com/OSGeo/gdal/blob/9f31018839b32aeeafad7663a8de662153a956c3/gdal/data/osmconf.ini#L65-L71)), and `other_tags`.

Please note that the argument `extra_tags` is always ignored (with a warning message), if you are using an ad-hoc `osmconf.ini` file: 

```{r}
oe_get("ITS Leeds", provider = "test", osmconf_ini = temp_ini, quiet = FALSE, extra_tags = "foot")
```

### `vectortranslate_options` + `boundary` and `boundary_type`

The parameter `vectortranslate_options` is used to modify the options that are passed to [ogr2ogr](https://gdal.org/programs/ogr2ogr.html#ogr2ogr). 
This is extremely important because if we tune the `vectortranslate_options` parameter, then we can analyse small parts of an enormous `.pbf` files without fully reading it in memory. 

The first example, reported in the following chunk, shows how to use the argument `-t_srs` to modify the CRS of the output `.gpkg` object (i.e. transform from [`EPSG:4326`](https://epsg.io/4326) to [`EPSG:27700`](https://epsg.io/27700)) while performing vectortranslate operations: 

```{r}
# Check the CRS
oe_get("ITS Leeds", vectortranslate_options = c("-t_srs", "EPSG:27700"), quiet = FALSE)
```

The default CRS of all OSM extracts obtained by Geofabrik and several other providers is `EPSG:4326`, i.e. latitude and longitude coordinates expressed via WGS84 ellipsoid, while the code `EPSG:27700` indicates the British National Grid. 
Hence, the parameter `-t_srs` can be used to transform geographical data into projected coordinates, which may be essential for some statistical software like `spatstat`. 
The same operation can also be performed in `R` with the `sf` package (e.g. `?st_transform()`), but the conversion can be slow for large spatial objects.
Please note that the default options (i.e. `c("-f", "GPKG", "-overwrite", "-oo", "CONFIG_FILE=", paste0("CONFIG_FILE=", osmconf_ini), "-lco", "GEOMETRY_NAME=geometry", layer)`) are internally appended to the `vectortranslate_options` argument. 

The next example illustrates how to apply an SQL-like query during the vectortranslate process. 
More precisely, we can use the arguments `-select` and `-where` to create an SQL-like query that is run during the vectortranslate process. 
Check [here](https://gdal.org/user/ogr_sql_dialect.html) for more details on the OGR SQL dialect. 

First of all, we need to build a character vector with the options that will be passed to `ogr2ogr: `

```{r}
my_vectortranslate = c(
  "-t_srs", "EPSG:27700", 
  # SQL-like query where we select only the following fields
  "-select", "osm_id, highway", 
  # SQL-like query where we filter only the features where highway is equal to footway or cycleway
  "-where", "highway IN ('footway', 'cycleway')"
)
```

and then we can process the file: 

```{r}
its_leeds = oe_get("ITS Leeds", vectortranslate_options = my_vectortranslate, quiet = FALSE)
```

The same procedure can be repeated using an ad-hoc `osmconf.ini` file.

These arguments are fundamental if you need to work with a small portion of a bigger `.pbf` file. 
For example, the following code (not run in the vignette) is used to extract all `primary`, `secondary` and `tertiary` roads from the `.pbf` file of Portugal stored by Geofabrik servers.
After downloading the data, it takes approximately 35 seconds to run the code using an HP ENVY Notebook with Intel i7-7500U processor and 8GB of RAM using Windows 10:

```{r, eval = FALSE}
# 1. Download the data and skip gpkg conversion
oe_get("Portugal", download_only = TRUE, skip_vectortranslate = TRUE)

# 2. Define the vectortranslate options
my_vectortranslate = c(
  # SQL-like query where we select only the features where highway in (primary, secondary, tertiary)
  "-select", "osm_id, highway",
  "-where", "highway IN ('primary', 'secondary', 'tertiary')"
)

# 3. Convert and read-in
system.time({
  portugal1 = oe_get("Portugal", vectortranslate_options = my_vectortranslate)
})
#  user  system elapsed 
# 17.39    9.93   25.53 
```

while the classical approach (also not run in the vignette) is slower and provides identical results: 

```{r, eval = FALSE}
system.time({
  portugal2 = oe_get("Portugal", quiet = FALSE, force_vectortranslate = TRUE)
  portugal2 = portugal2 %>% 
    dplyr::select(osm_id, highway) %>% 
    dplyr::filter(highway %in% c('primary', 'secondary', 'tertiary'))
})
#   user  system elapsed 
# 131.05   28.70  177.03

nrow(portugal1) == nrow(portugal2)
#> TRUE
```

Starting from version 0.3.0, the arguments `boundary` and `boundary_type` can be used to perform spatial filter operations during the vectortranslate process. 
In particular, a spatial boundary can be created using an `sf` or `sfc` object (with `POLYGON` or `MULTIPOLYGON` geometry) via the argument `boundary`: 

```{r}
its_bbox = st_bbox(c(xmin = -1.559184 , ymin = 53.807739 , xmax = -1.557375 , ymax = 53.808094), crs = 4326) %>% 
  st_as_sfc()

its_small = oe_get ("ITS Leeds", boundary = its_bbox)
```

This is the output, where the bounding box was highlighted in black, the intersecting streets in red and all the other roads in grey. 

```{r, echo = FALSE, out.width="85%"}
its_leeds = oe_get("ITS Leeds", force_vectortranslate = TRUE, quiet = TRUE)

par(mar = rep(0.1, 4))
plot(st_geometry(its_leeds), reset = FALSE, col = "grey")
plot(st_geometry(its_small), lwd = 3, col = "darkred", add = TRUE)
plot(st_as_sfc(st_bbox(c(xmin = -1.559184 , ymin = 53.807739 , xmax = -1.557375 , ymax = 53.808094), crs = 4326)), add = TRUE, lwd = 3)
```

Finally, the argument `boundary_type` can be used to select among different types of spatial filters. 
For the moment we support only two types of filters: `"spat"` (default value) and `"clipsrc"`.
The former option implies that the spatial filter selects all features that intersect a given area (as shown above), while the latter option implies that the features are also cropped. 
In both cases, the polygonal boundary must be specified as an `sf` or `sfc` object.   

The following example shows how to download from Geofabrik servers the `.pbf` extract associated with Malta and apply a spatial filter while performing vectortranslate operations. 
We select and clip only the road segments that intersect a 5 kilometres circular buffer centred in La Valletta, the capital. 

```{r, eval = FALSE}
# 1. Define the polygonal boundary
la_valletta = st_sfc(st_point(c(456113.1, 3972853)), crs = 32633) %>%
  st_buffer(5000)

# 2. Define the vectortranslate options
my_vectortranslate = c(
  "-t_srs", "EPSG:32633",
  "-select", "highway",
  "-where", "highway IN ('primary', 'secondary', 'tertiary', 'unclassified')",
  "-nlt", "PROMOTE_TO_MULTI"
)

# 3. Download data
oe_get("Malta", skip_vectortranslate = TRUE, download_only = TRUE)

# 4. Read-in data
system.time({
  oe_get("Malta", vectortranslate_options = my_vectortranslate, boundary = la_valletta, boundary_type = "clipsrc")
})
# The input place was matched with: Malta
# The chosen file was already detected in the download directory. Skip downloading.
# Start with the vectortranslate operations on the input file!
# 0...10...20...30...40...50...60...70...80...90...100 - done.
# Finished the vectortranslate operations on the input file!
# Reading layer `lines' from data source `C:\Users\Utente\AppData\Local\Temp\RtmpYVijx8\geofabrik_malta-latest.gpkg' using driver `GPKG'
# Simple feature collection with 1205 features and 1 field
# Geometry type: MULTILINESTRING
# Dimension:     XY
# Bounding box:  xmin: 451113.7 ymin: 3967858 xmax: 460364.8 ymax: 3976642
# Projected CRS: WGS 84 / UTM zone 33N
#    user  system elapsed 
#    0.55    0.11    0.61 
```

The options `-t_srs`, `-select` and `-where` have the same interpretation as before. 
The spatial filter may return invalid `LINESTRING` geometries (due to the cropping operation). 
For this reason, the `-nlt` and `PROMOTE_TO_MULTI` options are used to override the default geometry type and promote the `LINESTRING`(s) into `MULTILINESTRING`(s). 
You can use `st_cast()` to convert the `MULTILINESTRING` into `LINESTRING` (which may be important for some packages or functions). 

The following map represent the result, where we highlight the bounding circle and the road segments within using a dark-red colour, while all the other road segments are coloured in grey.

```{r, echo = FALSE, eval = FALSE}
malta_regular = oe_get("Malta", force_vectortranslate = TRUE) %>% 
  dplyr::filter(highway %in% c('primary', 'secondary', 'tertiary', 'unclassified')) %>% 
  st_transform(32633)
malta_small = oe_get("Malta", vectortranslate_options = my_vectortranslate, boundary = la_valletta, boundary_type = "clipsrc")

par(mar = rep(0.1, 4))
plot(st_geometry(malta_regular), col = "grey", reset = FALSE)
plot(st_boundary(la_valletta), add = TRUE, lwd = 2)
plot(st_geometry(malta_small), add = TRUE, col = "darkred", lwd = 2)
```

```{r, echo=FALSE, out.width="80%", fig.align="center"}
knitr::include_graphics("../man/figures/104240598-9d6fb400-545c-11eb-93b5-3563908ff4af.png")
```

The process takes approximately 1 or 2 seconds, while the equivalent R code, reported below, is slower: 

```{r, eval = FALSE}
system.time({
  malta_crop = oe_get("Malta", force_vectortranslate = TRUE) %>% 
    dplyr::filter(highway %in% c('primary', 'secondary', 'tertiary', 'unclassified')) %>% 
    st_transform(32633) %>% 
    st_crop(la_valletta)
})
#> user  system elapsed 
#> 4.61    1.67    7.69
```

The time difference gets more and more relevant for larger OSM extracts. 
Moreover, the R code crops the road segments using a rectangular boundary instead of the proper circular polygon: 

```{r, echo = FALSE, eval = FALSE}
malta_regular = oe_get("Malta", force_vectortranslate = TRUE) %>% 
  dplyr::filter(highway %in% c('primary', 'secondary', 'tertiary', 'unclassified')) %>% 
  st_transform(32633)

par(mar = rep(0.1, 4))
plot(st_geometry(malta_regular), col = "grey", reset = FALSE)
plot(st_boundary(la_valletta), add = TRUE, lwd = 2)
plot(st_geometry(malta_crop), add = TRUE, col = "darkred", lwd = 2)
```

```{r, echo=FALSE, out.width="80%", fig.align='center'}
knitr::include_graphics("../man/figures/104241581-32bf7800-545e-11eb-896b-3f535dd1af5e.png")
```

### `query` and `wkt_filter` arguments

The last two options that we introduce are `query` and `wkt_filter`. 
They are defined in the `R` package `sf` and represent a useful compromise between the `GDAL` and the `R` approaches explained above, especially when a user needs to apply different queries to the same (typically small or medium-size) OSM extract. 
In fact, the two parameters create regular queries and spatial filters, respectively, that are applied immediately before reading-in the `.gpkg` file. 
The following code, for example, mimics the operations illustrated above, reading-in the road segments that intersect the circular buffer defined around La Valletta:

```{r, eval = FALSE}
malta_small = oe_get(
  "Malta", 
  query = "
  SELECT highway, geometry 
  FROM 'lines' 
  WHERE highway IN ('primary', 'secondary', 'tertiary', 'unclassified')", 
  wkt_filter = st_as_text(st_transform(la_valletta, 4326)),
  force_vectortranslate = TRUE
)
```

This is the output and we can see that it applies a circular spatial filter but it doesn't crop the features: 

```{r, echo = FALSE, eval = FALSE}
malta_regular = oe_get("Malta", force_vectortranslate = TRUE) %>% 
  dplyr::filter(highway %in% c('primary', 'secondary', 'tertiary', 'unclassified'))

par(mar = rep(0.1, 4))
plot(st_geometry(malta_regular), col = "grey", reset = FALSE)
plot(st_boundary(la_valletta) %>% st_transform(4326), add = TRUE, lwd = 2)
plot(st_geometry(malta_small), col = "darkred", add = TRUE, lwd = 2)
```

```{r, echo = FALSE, fig.align="center", out.width="80%"}
knitr::include_graphics("../man/figures/104243054-4966ce80-5460-11eb-951b-ca1ce9d09f33.png")
```

This approach has its pros and cons. 
First of all, it’s slightly slower than the GDAL routines, mainly because several unnecessary features are being converted to the `.gpkg` format. 
Hence, it may become unfeasible for larger `.pbf` files, probably starting from 300/500MB. 
We will test more cases and add more benchmarks in the near future. 
On the other side, it does not require a new time-consuming `ogr2ogr` conversion every time a user defines a new query. 
For these reasons, this is the suggested approach for querying a small OSM extract.

Last but not least, we can use the function `hstore_get_value` to extract one of the tags saved in the `other_tags` column without using `ogr2ogr` and rerunning the `oe_vectortranslate()` function:: 

```{r}
# No extra tag
colnames(oe_get("ITS Leeds", quiet = TRUE))

# Check extra tags
oe_get_keys("ITS Leeds")

# Add extra tag
colnames(oe_get(
  "ITS Leeds", 
  provider = "test", 
  query = "SELECT *, hstore_get_value(other_tags, 'bicycle') AS bicycle FROM lines"
))
```

# Other providers

The package supports downloading, reading and extracting OpenStreetMap data from various providers.
A list of providers can be found at [wiki.openstreetmap.org](https://wiki.openstreetmap.org/wiki/Processed_data_providers).
The first provider supported was [Geofabrik](http://download.geofabrik.de/).
The second was [bbbike](https://download.bbbike.org/osm/bbbike/).
The package can be extended to support additional providers, as seen in the following [commit](https://github.com/ropensci/osmextract/commit/dbf131667a80e5a6837a6c8eb3b967075e1aba16) that adds a working provider.

For information on adding new providers to the package, see the [providers vignette](https://docs.ropensci.org/osmextract/articles/providers.html).

# More on OpenStreetMap

There is a world of knowledge, convention and wisdom contained in OSM data that we hope this package helps you discover and use this knowledge for public benefit.
To learn more about the structure of OSM data and the various tagging systems and conventions, the [Elements page on the OSM wiki](https://wiki.openstreetmap.org/wiki/Elements) is an ideal place to start.
You will find much more excellent content on the OSM wiki pages.

# Contributing to OSM

The final thing to say in this introductory vignette is that as a citizen-led project like Wikipedia, OSM relies on a participatory culture, where people not only consume but contribute data, to survive.
On that note, we urge anyone reading this to at least sign-up to get an OSM account at [osm.org](https://www.openstreetmap.org).

We highly recommend contributing to the world's geographic commons.
The step from being a user to being a contributor to OSM data is a small one and can be highly rewarding.
If you find any issues with OSM data, people in the OpenStreetMap will be very happy for you to correct the data.
Once logged-in, you can contribute by using editors such as the excellent ID editor, which you can get to by zooming into anywhere you want at [www.openstreetmap.org](https://www.openstreetmap.org/) and clicking "Edit".

To learn more about contributing to the amazing OSM community, we recommend checking out the [OSM Beginners Guide](https://wiki.openstreetmap.org/wiki/Beginners_Guide_1.3).

```{r, include=FALSE}
# reset par and options
options(user_options)
par(user_par)
```
---
title: "Add new OpenStreetMap providers"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Add new OpenStreetMap providers}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This vignette aims to provide a simple guide on adding a new provider to `osmextract`. 
Let' start loading the package: 

```{r setup}
library(osmextract)
```

As of summer 2020, there are several services providing bulk OSM datasets listed  [here](https://wiki.openstreetmap.org/wiki/Processed_data_providers) and [here](https://wiki.openstreetmap.org/wiki/Planet.osm#Country_and_area_extracts). 
At the moment, we support the following providers: 

```{r}
oe_providers()
```

Check the "Comparing the supported OSM providers" for more details on the existing providers. 

This package is designed to make it easy to add new providers. 
There are three main steps to add a new provider: creating the zones, adding the provider and documenting it. 
They are outlined below. 

# Adding a `provider_zones` object to the package

The first and hardest step is to create an `sf` object analogous to the `test_zones` object shown below:

```{r}
names(test_zones)
str(test_zones[, c(2, 6, 7)])
```

The output shows the three most important column names:

1. The zone `name` (that is used for matching the input `place`, see `oe_match()`); 
1. The URL endpoint where `.pbf` files associated with each zone can be downloaded; 
1. The geometry, representing the spatial extent of the dataset.

The object must also include the fields `level` and `id`, which are used, respectively, for spatial matching and updating. 
See `oe_match()` and `oe_update()`. 

The best way to start creating a new `_zones` object for a new provider is probably by looking at the code we wrote for the first supported provider in [`data-raw/geofabrik_zones.R`](https://github.com/ropensci/osmextract/blob/master/data-raw/geofabrik_zones.R). 
The following commands will clone this repo and open the relevant file:

```bash
git clone git@github.com:ropensci/osmextract
rstudio osmextract/osmextract.Rproj
```
Then in RStudio:

```{r, eval = FALSE}
file.edit("data-raw/geofabrik_zones.R")
```

Create a new script to document the code that generates the new object, e.g. for `bbbike`:

```{r, eval = FALSE}
file.edit("data-raw/bbbike_zones.R")
# or, even better, use
usethis::use_data_raw("bbbike_zones")
```

After you have created the new provider `_zones` file, it's time to add the provider to the package.

# Adding the new provider to the package

Once you have created your overview `_zones` file as outlined in the previous step, you need to modify the following files for the provider to be available for others:

- [data.R](https://github.com/ropensci/osmextract/blob/master/R/data.R), where you'll need to document the new dataset;
- [globals.R](https://github.com/ropensci/osmextract/blob/master/R/globals.R), where you'll need to add the new object name;
- [providers.R](https://github.com/ropensci/osmextract/blob/master/R/providers.R), where you'll need to add the new object name in `oe_available_providers()` and `load_provider_data()`. 

# Documenting the provider

The final step is also the most fun: documenting and using the provider.
Add an example, mention it in the README and tell others about what this new provider can do!
If you want to ask for help on adding a new provider, feel free to open in a new issue in the [github repository](https://github.com/ropensci/osmextract)! 

# Conclusion

This vignette talks through the main steps needed to extend `osmextract` by adding new OSM data providers.
To see the same information in code form, see the PR that implemented the `openstreetmap_fr` provider here: https://github.com/ropensci/osmextract/commit/dbf131667a80e5a6837a6c8eb3b967075e1aba16
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-network.R
\name{oe_get_network}
\alias{oe_get_network}
\title{Import transport network used by a specific mode of transport}
\usage{
oe_get_network(place, mode = c("cycling", "driving", "walking"), ...)
}
\arguments{
\item{place}{Description of the geographical area that should be matched with
a \code{.osm.pbf} file. Can be either a length-1 character vector, an
\code{sf}/\code{sfc}/\code{bbox} object, or a numeric vector of coordinates with length 2.
In the last case, it is assumed that the EPSG code is 4326 specified as
c(LON, LAT), while you can use any CRS with \code{sf}/\code{sfc}/\code{bbox} objects. See
Details and Examples in \code{\link[=oe_match]{oe_match()}}.}

\item{mode}{A character string denoting the desired mode of transport. Can be
abbreviated. Currently cycling (the default), driving and walking are supported.}

\item{...}{Additional arguments passed to \code{oe_get()} such as \code{boundary} or
\code{force_download}.}
}
\value{
An \code{sf} object.
}
\description{
This function is a wrapper around \code{oe_get()} and can be used to import a road
network given a \code{place} and a mode of transport. Check the Details for a
precise description of the procedures used to filter each mode of transport.
}
\details{
The definition of usable transport network was taken from the Python
packages \href{https://github.com/gboeing/osmnx/blob/main/osmnx/downloader.py}{osmnx} and
\href{https://pyrosm.readthedocs.io/en/latest/}{pyrosm} and several other
documents found online, i.e.
\url{https://wiki.openstreetmap.org/wiki/OSM_tags_for_routing/Access_restrictions},
\url{https://wiki.openstreetmap.org/wiki/Key:access}. See also the discussion
in \url{https://github.com/ropensci/osmextract/issues/153}.

The \code{cycling} mode of transport (i.e. the default value for \code{mode}
parameter) selects the OSM ways that meet the following conditions:
\itemize{
\item The \code{highway} tag is not missing and is not equal to \code{abandonded},
\code{bus_guideway}, \code{byway}, \code{construction}, \code{corridor}, \code{elevator}, \code{fixme},
\code{escalator}, \code{gallop}, \code{historic}, \code{no}, \code{planned}, \code{platform}, \code{proposed},
\code{raceway} or \code{steps};
\item The \code{highway} tag is not equal to \code{motorway}, \code{motorway_link}, \code{footway},
\code{bridleway} or \code{pedestrian} unless the tag \code{bicycle} is equal to \code{yes} (see
\href{https://wiki.openstreetmap.org/wiki/Bicycle#Bicycle_Restrictions}{here})
for more details;
\item The \code{access} tag is not equal to \code{private} or \code{no} unless \code{bicycle} tag
is equal to \code{yes};
\item The \code{bicycle} tag is not equal to \code{no}, \code{use_sidepath}, \code{private}, pr
\code{restricted};
\item The \code{service} tag does not contain the string \code{private} (i.e. \code{private};
\code{private_access} and similar);
}

The \code{walking} mode of transport selects the OSM ways that meet the
following conditions:
\itemize{
\item The \code{highway} tag is not missing and is not equal to \code{abandonded},
\code{bus_guideway}, \code{byway}, \code{construction}, \code{corridor}, \code{elevator}, \code{fixme},
\code{escalator}, \code{gallop}, \code{historic}, \code{no}, \code{planned}, \code{platform}, \code{proposed},
\code{raceway}, \code{motorway} or \code{motorway_link};
\item The \code{highway} tag is not equal to \code{cycleway} unless the \code{foot} tag is
equal to \code{yes};
\item The \code{access} tag is not equal to \code{private} or \code{no} unless \code{foot} tag
is equal to \code{yes};
\item The \code{foot} tag is not equal to \code{no}, \code{use_sidepath}, \code{private}, pr
\code{restricted};
\item The \code{service} tag does not contain the string \code{private} (i.e. \code{private};
\code{private_access} and similar).
}

The \code{driving} mode of transport selects the OSM ways that meet the
following conditions:
\itemize{
\item The \code{highway} tag is not missing and is not equal to \code{abandonded},
\code{bus_guideway}, \code{byway}, \code{construction}, \code{corridor}, \code{elevator}, \code{fixme},
\code{escalator}, \code{gallop}, \code{historic}, \code{no}, \code{planned}, \code{platform}, \code{proposed},
\code{cycleway}, \code{pedestrian}, \code{bridleway}, \code{path}, or \code{footway};
\item The \code{access} tag is not equal to \code{private} or \code{no};
\item The \code{service} tag does not contain the string \code{private} (i.e. \code{private};
\code{private_access} and similar).
}

Feel free to start a new issue in the \href{https://github.com/ropensci/osmextract}{github repo} if you want to suggest
modifications to the current filters or propose new values for alternative
modes of transport.
}
\examples{
\dontshow{
  its = file.copy(
    from = system.file("its-example.osm.pbf", package = "osmextract"),
    to = file.path(tempdir(), "test_its-example.osm.pbf"),
    overwrite = TRUE
)}
# default value returned by OSM
its = oe_get("ITS Leeds", quiet = TRUE, download_directory = tempdir())
plot(its["highway"], lwd = 2, key.pos = 4, key.width = lcm(2.75))
# walking mode of transport
its_walking = oe_get_network("ITS Leeds", mode = "walking", quiet = TRUE)
plot(its_walking["highway"], lwd = 2, key.pos = 4, key.width = lcm(2.75))
# driving mode of transport
its_driving = oe_get_network("ITS Leeds", mode = "driving", quiet = TRUE)
plot(its_driving["highway"], lwd = 2, key.pos = 4, key.width = lcm(2.75))

# Remove .pbf and .gpkg files in tempdir
# (since they may interact with other examples)
file.remove(list.files(path = tempdir(), pattern = "(pbf|gpkg)", full.names = TRUE))
}
\seealso{
\code{\link[=oe_get]{oe_get()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{oe_download_directory}
\alias{oe_download_directory}
\title{Return the download directory used by the package}
\usage{
oe_download_directory()
}
\value{
A character vector representing the path for the download directory
used by the package.
}
\description{
By default, the download directory is equal to \code{tempdir()}. You can set a
persistent download directory by adding the following command to your
\code{.Renviron} file (e.g. with \code{edit_r_environ} function in \code{usethis} package):
\verb{OSMEXT_DOWNLOAD_DIRECTORY=/path/to/osm/data}.
}
\examples{
oe_download_directory()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vectortranslate.R
\name{oe_vectortranslate}
\alias{oe_vectortranslate}
\title{Translate a .osm.pbf file into .gpkg format}
\usage{
oe_vectortranslate(
  file_path,
  layer = "lines",
  vectortranslate_options = NULL,
  osmconf_ini = NULL,
  extra_tags = NULL,
  force_vectortranslate = FALSE,
  never_skip_vectortranslate = FALSE,
  boundary = NULL,
  boundary_type = c("spat", "clipsrc"),
  quiet = FALSE
)
}
\arguments{
\item{file_path}{Character string representing the path of the input
\code{.pbf} or \code{.osm.pbf} file.}

\item{layer}{Which \code{layer} should be read in? Typically \code{points}, \code{lines}
(the default), \code{multilinestrings}, \code{multipolygons} or \code{other_relations}. If
you specify an ad-hoc query using the argument \code{query} (see introductory
vignette and examples), then \code{oe_get()} and \code{oe_read()} will read the layer
specified in the query and ignore \code{layer}. See also
\href{https://github.com/ropensci/osmextract/issues/122}{#122}.}

\item{vectortranslate_options}{Options passed to the \code{\link[sf:gdal_utils]{sf::gdal_utils()}}
argument \code{options}. Set by default. Check details in the introductory
vignette and the help page of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}.}

\item{osmconf_ini}{The configuration file. See documentation at
\href{https://gdal.org/drivers/vector/osm.html}{gdal.org}. Check details in the
introductory vignette and the help page of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}. Set by
default.}

\item{extra_tags}{Which additional columns, corresponding to OSM tags, should
be in the resulting dataset? \code{NULL} by default. Check the introductory
vignette and the help pages of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}} and \code{\link[=oe_get_keys]{oe_get_keys()}}.
Ignored when \code{osmconf_ini} is not \code{NULL}.}

\item{force_vectortranslate}{Boolean. Force the original \code{.pbf} file to be
translated into a \code{.gpkg} file, even if a \code{.gpkg} with the same name
already exists? \code{FALSE} by default. If tags in \code{extra_tags} match data in
previously translated \code{.gpkg} files no translation occurs (see
\href{https://github.com/ropensci/osmextract/issues/173}{#173} for details).
Check the introductory vignette and the help page of
\code{\link[=oe_vectortranslate]{oe_vectortranslate()}}.}

\item{never_skip_vectortranslate}{Boolean. This is used in case the user
passed its own \code{.ini} file or vectortranslate options (since, in those
case, it's too difficult to determine if an existing \code{.gpkg} file was
generated following the same options.)}

\item{boundary}{An \code{sf}/\code{sfc}/\code{bbox} object that will be used to create a
spatial filter during the vectortranslate operations. The type of filter
can be chosen using the argument \code{boundary_type}.}

\item{boundary_type}{A character vector of length 1 specifying the type of
spatial filter. The \code{spat} filter selects only those features that
intersect a given area, while \code{clipsrc} also clips the geometries. Check
the examples and also \href{https://gdal.org/programs/ogr2ogr.html}{here} for
more details.}

\item{quiet}{Boolean. If \code{FALSE}, the function prints informative messages.
Starting from \code{sf} version
\href{https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13}{0.9.6},
if \code{quiet} is equal to \code{FALSE}, then vectortranslate operations will
display a progress bar.}
}
\value{
Character string representing the path of the \code{.gpkg} file.
}
\description{
This function is used to translate a \code{.osm.pbf} file into \code{.gpkg} format.
The conversion is performed using
\href{https://gdal.org/programs/ogr2ogr.html#ogr2ogr}{ogr2ogr} via the
\code{vectortranslate} utility in \code{\link[sf:gdal_utils]{sf::gdal_utils()}} . It was created following
\href{https://github.com/OSGeo/gdal/issues/2100#issuecomment-565707053}{the suggestions}
of the maintainers of GDAL. See Details and Examples to understand the basic
usage, and check the introductory vignette for more complex use-cases.
}
\details{
The new \code{.gpkg} file is created in the same directory as the input
\code{.osm.pbf} file. The translation process is performed using the
\code{vectortranslate} utility in \code{\link[sf:gdal_utils]{sf::gdal_utils()}}. This operation can be
customized in several ways modifying the parameters \code{layer}, \code{extra_tags},
\code{osmconf_ini}, \code{vectortranslate_options}, \code{boundary} and \code{boundary_type}.

The \code{.osm.pbf} files processed by GDAL are usually categorized into 5
layers, named \code{points}, \code{lines}, \code{multilinestrings}, \code{multipolygons} and
\code{other_relations}. Check the first paragraphs
\href{https://gdal.org/drivers/vector/osm.html}{here} for more details. This
function can covert only one layer at a time, and the parameter \code{layer} is
used to specify which layer of the \code{.osm.pbf} file should be converted.
Several layers with different names can be stored in the same \code{.gpkg} file.
By default, the function will convert the \code{lines} layer (which is the most
common one according to our experience).

The arguments \code{osmconf_ini} and \code{extra_tags} are used to modify how GDAL
reads and processes a \code{.osm.pbf} file. More precisely, several operations
that GDAL performs on the input \code{.osm.pbf} file are governed by a \code{CONFIG}
file, that can be checked at the following
\href{https://github.com/OSGeo/gdal/blob/master/data/osmconf.ini}{link}.
The basic components of OSM data are called
\href{https://wiki.openstreetmap.org/wiki/Elements}{\emph{elements}} and they are
divided into \emph{nodes}, \emph{ways} or \emph{relations}, so, for example, the code at
line 7 of that file is used to determine which \emph{ways} are assumed to be
polygons (according to the simple-feature definition of polygon) if they
are closed. Moreover, OSM data is usually described using several
\href{https://wiki.openstreetmap.org/wiki/Tags}{\emph{tags}}, i.e pairs of two items:
a key and a value. The code at lines 33, 53, 85, 103, and 121 is used to
determine, for each layer, which tags should be explicitly reported as
fields (while all the other tags are stored in the \code{other_tags} column).
The parameter \code{extra_tags} is used to determine which extra tags (i.e.
key/value pairs) should be added to the \code{.gpkg} file (other than the
default ones).

By default, the vectortranslate operations are skipped if the function
detects a file having the same path as the input file, \code{.gpkg} extension, a
layer with the same name as the parameter \code{layer} and all \code{extra_tags}. In
that case the function will simply return the path of the \code{.gpkg} file.
This behaviour can be overwritten setting \code{force_vectortranslate = TRUE}.
The vectortranslate operations are never skipped if \code{osmconf_ini},
\code{vectortranslate_options}, \code{boundary} or \code{boundary_type} arguments are not
\code{NULL}.

The parameter \code{osmconf_ini} is used to pass your own \code{CONFIG} file in case
you need more control over the GDAL operations. Check the package
introductory vignette for an example. If \code{osmconf_ini} is equal to \code{NULL}
(the default value), then the function uses the standard \code{osmconf.ini} file
defined by GDAL (but for the extra tags).

The parameter \code{vectortranslate_options} is used to control the options that
are passed to \code{ogr2ogr} via \code{\link[sf:gdal_utils]{sf::gdal_utils()}} when converting between
\code{.osm.pbf} and \code{.gpkg} formats. \code{ogr2ogr} can perform various operations
during the conversion process, such as spatial filters or SQL queries.
These operations can be tuned using the \code{vectortranslate_options} argument.
If \code{NULL} (the default value), then \code{vectortranslate_options} is set equal
to

\code{c("-f", "GPKG", "-overwrite", "-oo", paste0("CONFIG_FILE=", osmconf_ini), "-lco", "GEOMETRY_NAME=geometry", layer)}.

Explanation:
\itemize{
\item \verb{"-f", "GPKG"} says that the output format is \code{GPKG};
\item \verb{"-overwrite} is used to delete an existing layer and recreate
it empty;
\item \verb{"-oo", paste0("CONFIG_FILE=", osmconf_ini)} is used to set the
\href{https://gdal.org/drivers/vector/osm.html#open-options}{Open Options}
for the \code{.osm.pbf} file and change the \code{CONFIG} file (in case the user
asks for any extra tag or a totally different CONFIG file);
\item \verb{"-lco", "GEOMETRY_NAME=geometry"} is used to change the
\href{https://gdal.org/drivers/vector/gpkg.html?highlight=gpkg#layer-creation-options}{layer creation options}
for the \code{.gpkg} file and modify the name of the geometry column;
\item \code{layer} indicates which layer should be converted.
}

If \code{vectortranslate_options} is not \code{NULL}, then the options \code{c("-f", "GPKG", "-overwrite", "-oo", "CONFIG_FILE=", path-to-config-file, "-lco", "GEOMETRY_NAME=geometry", layer)} are always appended unless the user
explicitly sets different default parameters for the arguments \code{-f}, \code{-oo},
\code{-lco}, and \code{layer}.

The arguments \code{boundary} and \code{boundary_type} can be used to set up a
spatial filter during the vectortranslate operations (and speed up the
process) using an \code{sf} or \code{sfc} object (\code{POLYGON} or \code{MULTIPOLYGON}). The
default arguments create a rectangular spatial filter which selects all
features that intersect the area. Setting \code{boundary_type = "clipsrc"} clips
the geometries. In both cases, the appropriate options are automatically
added to the \code{vectortranslate_options} (unless a user explicitly sets
different default options). Check Examples in \code{oe_get()} and the
introductory vignette.

See also the help page of \code{\link[sf:gdal_utils]{sf::gdal_utils()}} and
\href{https://gdal.org/programs/ogr2ogr.html}{ogr2ogr} for more examples and
extensive documentation on all available options that can be tuned during
the vectortranslate process.
}
\examples{
# First we need to match an input zone with a .osm.pbf file
its_match = oe_match("ITS Leeds")
\dontshow{file.copy(
  from = system.file("its-example.osm.pbf", package = "osmextract"),
  to = file.path(tempdir(), "test_its-example.osm.pbf"),
  overwrite = TRUE
)}
# The we can download the .osm.pbf file (is it was not already downloaded)
its_pbf = oe_download(
  file_url = its_match$url,
  file_size = its_match$file_size,
  download_directory = tempdir(),
  provider = "test"
)

# Check that the file was downloaded
list.files(tempdir(), pattern = "pbf|gpkg", full.names = TRUE)

# Convert to gpkg format
its_gpkg = oe_vectortranslate(its_pbf)

# Now there is an extra .gpkg file
list.files(tempdir(), pattern = "pbf|gpkg", full.names = TRUE)

# Check the layers of the .gpkg file
sf::st_layers(its_gpkg, do_count = TRUE)

# Add points layer
its_gpkg = oe_vectortranslate(its_pbf, layer = "points")
sf::st_layers(its_gpkg, do_count = TRUE)

# Add extra tags to the lines layer
names(sf::st_read(its_gpkg, layer = "lines", quiet = TRUE))
its_gpkg = oe_vectortranslate(
  its_pbf,
  extra_tags = c("oneway", "maxspeed")
)
names(sf::st_read(its_gpkg, layer = "lines", quiet = TRUE))

# Adjust vectortranslate options and convert only 10 features
# for the lines layer
oe_vectortranslate(
  its_pbf,
  vectortranslate_options = c("-limit", 10)
)
sf::st_layers(its_gpkg, do_count = TRUE)


# Remove .pbf and .gpkg files in tempdir
# (since they may interact with other examples)
file.remove(list.files(path = tempdir(), pattern = "(pbf|gpkg)", full.names = TRUE))
}
\seealso{
\code{\link[=oe_get_keys]{oe_get_keys()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match.R
\name{oe_match_pattern}
\alias{oe_match_pattern}
\title{Check patterns in the provider's databases}
\usage{
oe_match_pattern(pattern, match_by = "name", full_row = FALSE)
}
\arguments{
\item{pattern}{Character string representing the pattern that should be
explored.}

\item{match_by}{Column name of the provider's database that will be used to
find the match.}

\item{full_row}{Boolean. Return all columns for the matching rows? \code{FALSE} by
default.}
}
\value{
A list of character vectors or \code{sf} objects (according to the value
of the parameter \code{full_row}). If no OSM zone can be matched with the input
string, then the function returns an empty list.
}
\description{
This function is used to explore the provider's databases and look for
patterns. This function can be useful in combination with \code{\link[=oe_match]{oe_match()}} and
\code{\link[=oe_get]{oe_get()}} for an easy match. See Examples.
}
\examples{
oe_match_pattern("Yorkshire")

res = oe_match_pattern("Yorkshire", full_row = TRUE)
lapply(res, function(x) sf::st_drop_geometry(x)[, 1:3])

oe_match_pattern("ABC")
oe_match_pattern("Yorkshire", match_by = "ABC")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.R
\name{oe_get}
\alias{oe_get}
\title{Find, download, translate and read OSM extracts from several providers}
\usage{
oe_get(
  place,
  layer = "lines",
  ...,
  provider = "geofabrik",
  match_by = "name",
  max_string_dist = 1,
  level = NULL,
  download_directory = oe_download_directory(),
  force_download = FALSE,
  max_file_size = 5e+08,
  vectortranslate_options = NULL,
  osmconf_ini = NULL,
  extra_tags = NULL,
  force_vectortranslate = FALSE,
  boundary = NULL,
  boundary_type = c("spat", "clipsrc"),
  download_only = FALSE,
  skip_vectortranslate = FALSE,
  never_skip_vectortranslate = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{place}{Description of the geographical area that should be matched with
a \code{.osm.pbf} file. Can be either a length-1 character vector, an
\code{sf}/\code{sfc}/\code{bbox} object, or a numeric vector of coordinates with length 2.
In the last case, it is assumed that the EPSG code is 4326 specified as
c(LON, LAT), while you can use any CRS with \code{sf}/\code{sfc}/\code{bbox} objects. See
Details and Examples in \code{\link[=oe_match]{oe_match()}}.}

\item{layer}{Which \code{layer} should be read in? Typically \code{points}, \code{lines}
(the default), \code{multilinestrings}, \code{multipolygons} or \code{other_relations}. If
you specify an ad-hoc query using the argument \code{query} (see introductory
vignette and examples), then \code{oe_get()} and \code{oe_read()} will read the layer
specified in the query and ignore \code{layer}. See also
\href{https://github.com/ropensci/osmextract/issues/122}{#122}.}

\item{...}{(Named) arguments that will be passed to \code{\link[sf:st_read]{sf::st_read()}}, like
\code{query}, \code{wkt_filter} or \code{stringsAsFactors}.  Check the introductory
vignette to understand how to create your own (SQL-like) queries.}

\item{provider}{Which provider should be used to download the data? Available
providers can be found with the following command: \code{\link[=oe_providers]{oe_providers()}}. For
\code{\link[=oe_get]{oe_get()}} and \code{\link[=oe_match]{oe_match()}}, if \code{place} is equal to \verb{ITS Leeds}, then
\code{provider} is set equal to \code{test}. This is just for simple examples and
internal tests.}

\item{match_by}{Which column of the provider's database should be used for
matching the input \code{place} with a \code{.osm.pbf} file? The default is \code{"name"}.
Check Details and Examples in \code{\link[=oe_match]{oe_match()}} to understand how this parameter
works. Ignored if \code{place} is not a character vector since the matching is
performed through a spatial operation.}

\item{max_string_dist}{Numerical value greater or equal than 0. What is the
maximum distance in fuzzy matching (i.e. Approximate String Distance, see
\code{\link[=adist]{adist()}}) between input \code{place} and \code{match_by} column to tolerate before
testing alternative providers or looking for geographical matching with
Nominatim API? This parameter is set equal to 0 if \code{match_by} is equal to
\code{iso3166_1_alpha2} or \code{iso3166_2}. Check Details and Examples in
\code{\link[=oe_match]{oe_match()}} to understand why this parameter is important. Ignored if
\code{place} is not a character vector since the matching is performed through a
spatial operation.}

\item{level}{An integer representing the desired hierarchical level in case
of spatial matching. For the \code{geofabrik} provider, for example, \code{1}
corresponds with continent-level datasets, \code{2} for countries, \code{3}
corresponds to regions and \code{4} to subregions. Hence, we could approximately
say that smaller administrative units correspond to bigger levels. If
\code{NULL}, the default, the \verb{oe_*} functions will select the highest available
level. See Details and Examples in \code{oe_match()}.}

\item{download_directory}{Where to download the file containing the OSM data?
By default this is equal to \code{\link[=oe_download_directory]{oe_download_directory()}}, which is equal to
\code{\link[=tempdir]{tempdir()}} and it changes each time you restart R. You can set a
persistent \code{download_directory} by adding the following to your \code{.Renviron}
file (e.g. with \code{edit_r_environ} function in \code{usethis} package):
\verb{OSMEXT_DOWNLOAD_DIRECTORY=/path/to/osm/data}.}

\item{force_download}{Should the \code{.osm.pbf} file be updated if it has already
been downloaded? \code{FALSE} by default. This parameter is used to update old
\code{.osm.pbf} files.}

\item{max_file_size}{The maximum file size to download without asking in
interactive mode. Default: \code{5e+8}, half a gigabyte.}

\item{vectortranslate_options}{Options passed to the \code{\link[sf:gdal_utils]{sf::gdal_utils()}}
argument \code{options}. Set by default. Check details in the introductory
vignette and the help page of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}.}

\item{osmconf_ini}{The configuration file. See documentation at
\href{https://gdal.org/drivers/vector/osm.html}{gdal.org}. Check details in the
introductory vignette and the help page of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}. Set by
default.}

\item{extra_tags}{Which additional columns, corresponding to OSM tags, should
be in the resulting dataset? \code{NULL} by default. Check the introductory
vignette and the help pages of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}} and \code{\link[=oe_get_keys]{oe_get_keys()}}.
Ignored when \code{osmconf_ini} is not \code{NULL}.}

\item{force_vectortranslate}{Boolean. Force the original \code{.pbf} file to be
translated into a \code{.gpkg} file, even if a \code{.gpkg} with the same name
already exists? \code{FALSE} by default. If tags in \code{extra_tags} match data in
previously translated \code{.gpkg} files no translation occurs (see
\href{https://github.com/ropensci/osmextract/issues/173}{#173} for details).
Check the introductory vignette and the help page of
\code{\link[=oe_vectortranslate]{oe_vectortranslate()}}.}

\item{boundary}{An \code{sf}/\code{sfc}/\code{bbox} object that will be used to create a
spatial filter during the vectortranslate operations. The type of filter
can be chosen using the argument \code{boundary_type}.}

\item{boundary_type}{A character vector of length 1 specifying the type of
spatial filter. The \code{spat} filter selects only those features that
intersect a given area, while \code{clipsrc} also clips the geometries. Check
the examples and also \href{https://gdal.org/programs/ogr2ogr.html}{here} for
more details.}

\item{download_only}{Boolean. If \code{TRUE}, then the function only returns the
path where the matched file is stored, instead of reading it. \code{FALSE} by
default.}

\item{skip_vectortranslate}{Boolean. If \code{TRUE}, then the function skips all
vectortranslate operations and it reads (or simply returns the path) of the
\code{.osm.pbf} file. \code{FALSE} by default.}

\item{never_skip_vectortranslate}{Boolean. This is used in case the user
passed its own \code{.ini} file or vectortranslate options (since, in those
case, it's too difficult to determine if an existing \code{.gpkg} file was
generated following the same options.)}

\item{quiet}{Boolean. If \code{FALSE}, the function prints informative messages.
Starting from \code{sf} version
\href{https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13}{0.9.6},
if \code{quiet} is equal to \code{FALSE}, then vectortranslate operations will
display a progress bar.}
}
\value{
An \code{sf} object.
}
\description{
This function is used to find, download, translate and read OSM extracts
obtained from several providers. It is a wrapper around \code{\link[=oe_match]{oe_match()}} and
\code{\link[=oe_read]{oe_read()}}. Check the introductory vignette, the examples and the help pages
of the wrapped functions to understand the details behind all parameters.
}
\details{
The algorithm that we use for importing an OSM extract data into R
is divided into 4 steps: 1) match the input \code{place} with the url of a
\code{.pbf} file; 2) download the \code{.pbf} file; 3) convert it into \code{.gpkg} format
and 4) read-in the \code{.gpkg} file. The function \code{oe_match()} is used to
perform the first operation and the function \code{oe_read()} (which is a
wrapper around \code{oe_download()}, \code{oe_vectortranslate()} and \code{sf::st_read()})
performs the other three operations.
}
\examples{
# Copy ITS file to tempdir so that the examples do not require internet
# connection. You can skip the next few lines when running the examples
# locally.
its_pbf = file.path(tempdir(), "test_its-example.osm.pbf")
file.copy(
  from = system.file("its-example.osm.pbf", package = "osmextract"),
  to = its_pbf,
  overwrite = TRUE
)

# Match, download (not really) and convert OSM extracts associated to a simple test.
its = oe_get("ITS Leeds", quiet = FALSE, download_directory = tempdir())
class(its)
unique(sf::st_geometry_type(its))

# Get another layer from ITS Leeds extract
its_points = oe_get("ITS Leeds", layer = "points")
unique(sf::st_geometry_type(its_points))

# Get the .osm.pbf and .gpkg files paths
oe_get("ITS Leeds", download_only = TRUE, quiet = TRUE)
oe_get("ITS Leeds", download_only = TRUE, skip_vectortranslate = TRUE, quiet = TRUE)
# See also ?oe_find()

# Add additional tags
its_with_oneway = oe_get("ITS Leeds", extra_tags = "oneway")
names(its_with_oneway)
table(its_with_oneway$oneway, useNA = "ifany")

# Use the query argument to get only oneway streets:
q = "SELECT * FROM 'lines' WHERE oneway == 'yes'"
its_oneway = oe_get("ITS Leeds", query = q)
its_oneway[, c(1, 3, 9)]

# Apply a spatial filter during the vectortranslate operations
its_poly = sf::st_sfc(
  sf::st_polygon(
    list(rbind(
      c(-1.55577, 53.80850),
      c(-1.55787, 53.80926),
      c(-1.56096, 53.80891),
      c(-1.56096, 53.80736),
      c(-1.55675, 53.80658),
      c(-1.55495, 53.80749),
      c(-1.55577, 53.80850)
    ))
  ),
  crs = 4326
)
its_spat = oe_get("ITS Leeds", boundary = its_poly)
its_clipped = oe_get("ITS Leeds", boundary = its_poly, boundary_type = "clipsrc", quiet = TRUE)

plot(sf::st_geometry(its), reset = FALSE, col = "lightgrey")
plot(sf::st_boundary(its_poly), col = "black", add = TRUE)
plot(sf::st_boundary(sf::st_as_sfc(sf::st_bbox(its_poly))), col = "black", add = TRUE)
plot(sf::st_geometry(its_spat), add = TRUE, col = "darkred")
plot(sf::st_geometry(its_clipped), add = TRUE, col = "orange")

# More complex examples
\dontrun{
west_yorkshire = oe_get("West Yorkshire")
# If you run it again, the function will not download the file
# or convert it again
west_yorkshire = oe_get("West Yorkshire")
# Match with place name
oe_get("Milan") # Warning: the .pbf file is 400MB
oe_get("Vatican City") # Check all providers
oe_get("Zurich") # Use Nominatim API for geolocating places

# Match with coordinates (any EPSG)
milan_duomo = sf::st_sfc(sf::st_point(c(1514924, 5034552)), crs = 3003)
oe_get(milan_duomo, quiet = FALSE) # Warning: the .pbf file is 400MB
# Match with numeric coordinates (EPSG = 4326)
oe_match(c(9.1916, 45.4650), quiet = FALSE)

# Check also alternative providers
baku = oe_get(place = "Baku")

# Other examples:
oe_get("RU", match_by = "iso3166_1_alpha2", quiet = FALSE)
# The following example mimics read_sf
oe_get("Andora", stringsAsFactors = FALSE, quiet = TRUE, as_tibble = TRUE)}

# Remove .pbf and .gpkg files in tempdir
# (since they may interact with other examples)
file.remove(list.files(path = tempdir(), pattern = "(pbf|gpkg)", full.names = TRUE))
}
\seealso{
\code{\link[=oe_match]{oe_match()}}, \code{\link[=oe_download]{oe_download()}}, \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}, and
\code{\link[=oe_read]{oe_read()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{oe_download}
\alias{oe_download}
\title{Download a file given a url}
\usage{
oe_download(
  file_url,
  provider = NULL,
  file_basename = basename(file_url),
  download_directory = oe_download_directory(),
  file_size = NA,
  force_download = FALSE,
  max_file_size = 5e+08,
  quiet = FALSE
)
}
\arguments{
\item{file_url}{A URL pointing to a \code{.osm.pbf} file that should be
downloaded.}

\item{provider}{Which provider stores the file? If \code{NULL} (the default), it
may be inferred from the URL, but it must be specified for non-standard
cases. See details and examples.}

\item{file_basename}{The basename of the file. The default behaviour is to
auto-generate it from the URL using \code{basename()}.}

\item{download_directory}{Where to download the file containing the OSM data?
By default this is equal to \code{\link[=oe_download_directory]{oe_download_directory()}}, which is equal to
\code{\link[=tempdir]{tempdir()}} and it changes each time you restart R. You can set a
persistent \code{download_directory} by adding the following to your \code{.Renviron}
file (e.g. with \code{edit_r_environ} function in \code{usethis} package):
\verb{OSMEXT_DOWNLOAD_DIRECTORY=/path/to/osm/data}.}

\item{file_size}{How big is the file? Optional. \code{NA} by default. If it's
bigger than \code{max_file_size} and the function is run in interactive mode,
then an interactive menu is displayed, asking for permission for
downloading the file.}

\item{force_download}{Should the \code{.osm.pbf} file be updated if it has already
been downloaded? \code{FALSE} by default. This parameter is used to update old
\code{.osm.pbf} files.}

\item{max_file_size}{The maximum file size to download without asking in
interactive mode. Default: \code{5e+8}, half a gigabyte.}

\item{quiet}{Boolean. If \code{FALSE}, the function prints informative messages.
Starting from \code{sf} version
\href{https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13}{0.9.6},
if \code{quiet} is equal to \code{FALSE}, then vectortranslate operations will
display a progress bar.}
}
\value{
A character string representing the file's path.
}
\description{
This function is used to download a file given a URL. It focuses on OSM
extracts with \code{.osm.pbf} format stored by one of the providers implemented in
the package. The URL is specified through the parameter \code{file_url}.
}
\details{
This function runs several checks before actually downloading a new
file to avoid overloading the OSM providers. The first step is the
definition of the file's path associated to the input \code{file_url}. The path
is created by pasting together the \code{download_directory}, the name of chosen
provider (which may be inferred from the URL) and the \code{basename()} of the
URL. For example, if \code{file_url} is equal to
\code{"https://download.geofabrik.de/europe/italy-latest.osm.pbf"}, and
\code{download_directory = "/tmp"}, then the path is built as
\code{"/tmp/geofabrik_italy-latest.osm.pbf"}. Thereafter, the function checks
the existence of that file and, if it founds it, then it returns the path.
The parameter \code{force_download} is used to modify this behaviour. If there
is no file associated with the new path, then the function downloads a new
file using \code{\link[=download.file]{download.file()}} with \code{mode = "wb"}, and, again, it returns the
path.
}
\examples{
its_match = oe_match("ITS Leeds", quiet = TRUE)

\dontrun{
oe_download(
  file_url = its_match$url,
  file_size = its_match$file_size,
  provider = "test",
  download_directory = tempdir()
)
iow_url = oe_match("Isle of Wight")
oe_download(
  file_url = iow_url$url,
  file_size = iow_url$file_size,
  download_directory = tempdir()
)
Sucre_url = oe_match("Sucre", provider = "bbbike")
oe_download(
  file_url = Sucre_url$url,
  file_size = Sucre_url$file_size,
  download_directory = tempdir()
)}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{oe_update}
\alias{oe_update}
\title{Update all the .osm.pbf files saved in a directory}
\usage{
oe_update(
  download_directory = oe_download_directory(),
  quiet = FALSE,
  delete_gpkg = TRUE,
  max_file_size = 5e+08,
  ...
)
}
\arguments{
\item{download_directory}{Character string of the path of the directory
where the \code{.osm.pbf} files are saved.}

\item{quiet}{Boolean. If \code{FALSE} the function prints informative
messages. See Details.}

\item{delete_gpkg}{Boolean. if \code{TRUE} the function deletes the old \code{.gpkg}
files. We added this parameter to minimize the probability of accidentally
reading-in old and not-synchronized \code{.gpkg} files. See Details. Defaults to
\code{TRUE}.}

\item{max_file_size}{The maximum file size to download without asking in
interactive mode. Default: \code{5e+8}, half a gigabyte.}

\item{...}{Additional parameter that will be passed to \code{\link[=oe_get]{oe_get()}} (such as
\code{stringsAsFactors} or \code{query}).}
}
\value{
The path(s) of the \code{.osm.pbf} file(s) that were updated.
}
\description{
This function is used to re-download all \code{.osm.pbf} files stored in
\code{download_directory} that were firstly downloaded through \code{\link[=oe_get]{oe_get()}}. See
Details.
}
\details{
This function is used to re-download \code{.osm.pbf} files that are
stored in a directory (specified by \code{download_directory} param) and that
were firstly downloaded through \code{\link[=oe_get]{oe_get()}} . The name of the files must
begin with the name of one of the supported providers (see
\code{\link[=oe_providers]{oe_providers()}}) and it must end with \code{.osm.pbf}. All other
files in the directory that do not match this format are ignored.

The process for re-downloading the \code{.osm.pbf} files is performed using the
function \code{\link[=oe_get]{oe_get()}} . The appropriate provider is determined by looking at
the first word in the path of the \code{.osm.pbf} file. The place is determined
by looking at the second word in the file path and the matching is
performed through the \code{id} column in the provider's database. So, for
example, the path \code{geofabrik_italy-latest-update.osm.pbf} will be matched
with the provider \code{"geofabrik"} and the geographical zone \code{italy} through
the column \code{id} in \code{geofabrik_zones}.

The parameter \code{delete_gpkg} is used to delete all \code{.gpkg} files in
\code{download_directory}. We decided to set its default value to \code{TRUE} to
minimize the possibility of reading-in old and non-synchronized \code{.gpkg}
files. If you set \code{delete_gpkg = FALSE}, then you need to manually
reconvert all files using \code{\link[=oe_get]{oe_get()}} or \code{\link[=oe_vectortranslate]{oe_vectortranslate()}} .

If you set the parameter \code{quiet} to \code{FALSE}, then the function will print
some useful messages regarding the characteristics of the files before and
after updating them. More precisely, it will print the output of the
columns \code{size}, \code{mtime} and \code{ctime} from \code{\link[=file.info]{file.info()}}. Please note that
the meaning of \code{mtime} and \code{ctime} depends on the OS and the file system.
Check \code{\link[=file.info]{file.info()}}.
}
\examples{
\dontrun{
# Set up a fake directory with .pbf and .gpkg files
fake_dir = tempdir()
# Fill the directory
oe_get("Andorra", download_directory = fake_dir, download_only = TRUE)
# Check the directory
list.files(fake_dir, pattern = "gpkg|pbf")
# Update all .pbf files and delete all .gpkg files
oe_update(fake_dir)
list.files(fake_dir, pattern = "gpkg|pbf")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/providers.R
\name{oe_providers}
\alias{oe_providers}
\title{Summary of available providers}
\usage{
oe_providers(quiet = FALSE)
}
\arguments{
\item{quiet}{Boolean. If \code{FALSE}, the function prints informative messages.
Starting from \code{sf} version
\href{https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13}{0.9.6},
if \code{quiet} is equal to \code{FALSE}, then vectortranslate operations will
display a progress bar.}
}
\value{
A \code{data.frame} with 4 columns representing the name of each available
provider, the name of the corresponding database and the number of features
and fields.
}
\description{
This function is used to display a short summary of the major characteristics
of the databases associated to all available providers.
}
\examples{
oe_providers()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match.R
\name{oe_match}
\alias{oe_match}
\alias{oe_match.default}
\alias{oe_match.bbox}
\alias{oe_match.sf}
\alias{oe_match.sfc}
\alias{oe_match.numeric}
\alias{oe_match.character}
\title{Match input place with a url}
\usage{
oe_match(place, ...)

\method{oe_match}{default}(place, ...)

\method{oe_match}{bbox}(place, ...)

\method{oe_match}{sf}(place, ...)

\method{oe_match}{sfc}(place, provider = "geofabrik", level = NULL, quiet = FALSE, ...)

\method{oe_match}{numeric}(place, provider = "geofabrik", quiet = FALSE, ...)

\method{oe_match}{character}(
  place,
  provider = "geofabrik",
  quiet = FALSE,
  match_by = "name",
  max_string_dist = 1,
  ...
)
}
\arguments{
\item{place}{Description of the geographical area that should be matched with
a \code{.osm.pbf} file. Can be either a length-1 character vector, an
\code{sf}/\code{sfc}/\code{bbox} object, or a numeric vector of coordinates with length 2.
In the last case, it is assumed that the EPSG code is 4326 specified as
c(LON, LAT), while you can use any CRS with \code{sf}/\code{sfc}/\code{bbox} objects. See
Details and Examples in \code{\link[=oe_match]{oe_match()}}.}

\item{...}{arguments passed to other methods}

\item{provider}{Which provider should be used to download the data? Available
providers can be found with the following command: \code{\link[=oe_providers]{oe_providers()}}. For
\code{\link[=oe_get]{oe_get()}} and \code{\link[=oe_match]{oe_match()}}, if \code{place} is equal to \verb{ITS Leeds}, then
\code{provider} is set equal to \code{test}. This is just for simple examples and
internal tests.}

\item{level}{An integer representing the desired hierarchical level in case
of spatial matching. For the \code{geofabrik} provider, for example, \code{1}
corresponds with continent-level datasets, \code{2} for countries, \code{3}
corresponds to regions and \code{4} to subregions. Hence, we could approximately
say that smaller administrative units correspond to bigger levels. If
\code{NULL}, the default, the \verb{oe_*} functions will select the highest available
level. See Details and Examples in \code{oe_match()}.}

\item{quiet}{Boolean. If \code{FALSE}, the function prints informative messages.
Starting from \code{sf} version
\href{https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13}{0.9.6},
if \code{quiet} is equal to \code{FALSE}, then vectortranslate operations will
display a progress bar.}

\item{match_by}{Which column of the provider's database should be used for
matching the input \code{place} with a \code{.osm.pbf} file? The default is \code{"name"}.
Check Details and Examples in \code{\link[=oe_match]{oe_match()}} to understand how this parameter
works. Ignored if \code{place} is not a character vector since the matching is
performed through a spatial operation.}

\item{max_string_dist}{Numerical value greater or equal than 0. What is the
maximum distance in fuzzy matching (i.e. Approximate String Distance, see
\code{\link[=adist]{adist()}}) between input \code{place} and \code{match_by} column to tolerate before
testing alternative providers or looking for geographical matching with
Nominatim API? This parameter is set equal to 0 if \code{match_by} is equal to
\code{iso3166_1_alpha2} or \code{iso3166_2}. Check Details and Examples in
\code{\link[=oe_match]{oe_match()}} to understand why this parameter is important. Ignored if
\code{place} is not a character vector since the matching is performed through a
spatial operation.}
}
\value{
A list with two elements, named \code{url} and \code{file_size}. The first
element is the URL of the \code{.osm.pbf} file associated with the input
\code{place}, while the second element is the size of the file in bytes (which
may be \code{NULL} or \code{NA})
}
\description{
This function is used to match an input \code{place} with the URL of a \code{.osm.pbf}
file (and its file-size, if present). The URLs are stored in several
provider's databases. See \code{\link[=oe_providers]{oe_providers()}} and examples.
}
\details{
If the input place is specified as a spatial object (either \code{sf} or \code{sfc}),
then the function will return a geographical area that completely contains
the object (or an error). The argument \code{level} (which must be specified as an
integer between 1 and 4, extreme values included) is used to select between
multiple geographically nested areas. We could roughly say that smaller
administrative units correspond to higher levels. Check the help page of the
chosen provider for more details on \code{level} field. By default, \code{level = NULL}, which means that \code{oe_match()} will return the area corresponding to
the highest available level. If there is no geographical area at the desired
level, then the function will return an error. If there are multiple areas at
the same \code{level} intersecting the input place, then the function will return
the area whose centroid is closest to the input place.

If the input place is specified as a character vector and there are multiple
plausible matches between the input place and the \code{match_by} column, then the
function will return a warning and it will select the first match. See
Examples. On the other hand, if the approximate string distance between the
input \code{place} and the best match in \code{match_by} column is greater than
\code{max_string_dist}, then the function will look for exact matches (i.e.
\code{max_string_dist = 0}) in the other supported providers. If it finds an exact
match, then it will return the corresponding URL. Otherwise, if \code{match_by} is
equal to \code{"name"}, then it will try to geolocate the input \code{place} using the
\href{https://nominatim.org/release-docs/develop/api/Overview/}{Nominatim API},
and then it will perform a spatial matching operation (see Examples and
introductory vignette), while, if \code{match_by != "name"}, then it will return
an error.

The fields \code{iso3166_1_alpha2} and \code{iso3166_2} are used by Geofabrik provider
to perform matching operations using \href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}{ISO 3166-1 alpha-2} and \href{https://en.wikipedia.org/wiki/ISO_3166-2}{ISO 3166-2} codes. See
\link{geofabrik_zones} for more details.
}
\examples{
# The simplest example:
oe_match("Italy")

# The default provider is "geofabrik", but we can change that:
oe_match("Leeds", provider = "bbbike")

# By default, the matching operations are performed through the column
# "name" in the provider's database but this can be a problem. Hence,
# you can perform the matching operations using other columns:
oe_match("RU", match_by = "iso3166_1_alpha2")
# Run oe_providers() for reading a short description of all providers and
# check the help pages of the corresponding databases to learn which fields
# are present.

# You can always increase the max_string_dist argument, but it can be
# dangerous:
oe_match("London", max_string_dist = 3, quiet = FALSE)

# Match the input zone using an sfc object:
milan_duomo = sf::st_sfc(sf::st_point(c(1514924, 5034552)), crs = 3003)
oe_match(milan_duomo, quiet = FALSE)
leeds = sf::st_sfc(sf::st_point(c(430147.8, 433551.5)), crs = 27700)
oe_match(leeds, provider = "bbbike")

# If you specify more than one sfg object, then oe_match will select the OSM
# extract that covers all areas
milan_leeds = sf::st_sfc(
  sf::st_point(c(9.190544, 45.46416)), # Milan
  sf::st_point(c(-1.543789, 53.7974)), # Leeds
  crs = 4326
)
oe_match(milan_leeds)

# Match the input zone using a numeric vector of coordinates
# (in which case crs = 4326 is assumed)
oe_match(c(9.1916, 45.4650)) # Milan, Duomo using CRS = 4326

# The following returns a warning since Berin is matched both
# with Benin and Berlin
oe_match("Berin", quiet = FALSE)

# If the input place does not match any zone in the chosen provider, then the
# function will test the other providers:
oe_match("Leeds")

# If the input place cannot be exactly matched with any zone in any provider,
# then the function will try to geolocate the input and then it will perform a
# spatial match:
\dontrun{
oe_match("Milan")}

# The level parameter can be used to select smaller or bigger geographical
# areas during spatial matching
yak = c(-120.51084, 46.60156)
\dontrun{
oe_match(yak, level = 3) # error
oe_match(yak, level = 2) # by default, level is equal to the maximum value
oe_match(yak, level = 1)}
}
\seealso{
\code{\link[=oe_providers]{oe_providers()}} and \code{\link[=oe_match_pattern]{oe_match_pattern()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osmextract-package.R
\docType{package}
\name{osmextract-package}
\alias{osmextract}
\alias{osmextract-package}
\title{osmextract: Download and Import Open Street Map Data Extracts}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Match, download, convert and import Open Street Map data extracts obtained from several providers.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/osmextract/}
  \item \url{https://github.com/ropensci/osmextract}
  \item Report bugs at \url{https://github.com/ropensci/osmextract/issues}
}

}
\author{
\strong{Maintainer}: Andrea Gilardi \email{andrea.gilardi@unimib.it} (\href{https://orcid.org/0000-0002-9424-7439}{ORCID})

Authors:
\itemize{
  \item Robin Lovelace (\href{https://orcid.org/0000-0001-5679-6536}{ORCID})
}

Other contributors:
\itemize{
  \item Barry Rowlingson (\href{https://orcid.org/0000-0002-8586-6625}{ORCID}) [contributor]
  \item Salva Fernández (Salva reviewed the package (v. 0.1) for rOpenSci, see <https://github.com/ropensci/software-review/issues/395>) [reviewer]
  \item Nicholas Potter (Nicholas reviewed the package (v. 0.1) for rOpenSci, see <https://github.com/ropensci/software-review/issues/395>) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{oe_search}
\alias{oe_search}
\title{Search for a place and return an sf data frame locating it}
\usage{
oe_search(
  place,
  base_url = "https://nominatim.openstreetmap.org",
  destfile = tempfile(fileext = ".geojson"),
  ...
)
}
\arguments{
\item{place}{Text string containing the name of a place the location of
which is to be found, such as \code{"Leeds"} or \code{"Milan"}.}

\item{base_url}{The URL of the nominatim server to use. The main
open server hosted by OpenStreetMap is the default.}

\item{destfile}{The name of the destination file where the output
of the search query, a \code{.geojson} file, should be saved.}

\item{...}{Extra arguments that are passed to \code{sf::st_read}.}
}
\value{
An \code{sf} object corresponding to the input place. The \code{sf} object is
read by \code{sf::st_read()} and it is based on a \code{geojson} file returned by
Nominatim API.
}
\description{
This (only internal and experimental) function provides a simple
interface to the \href{https://nominatim.openstreetmap.org}{nominatim} service for
finding the geographical location of place names.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find.R
\name{oe_find}
\alias{oe_find}
\title{Get the location of files}
\usage{
oe_find(
  place,
  provider = "geofabrik",
  download_directory = oe_download_directory(),
  download_if_missing = FALSE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{place}{Description of the geographical area that should be matched with
a \code{.osm.pbf} file. Can be either a length-1 character vector, an
\code{sf}/\code{sfc}/\code{bbox} object, or a numeric vector of coordinates with length 2.
In the last case, it is assumed that the EPSG code is 4326 specified as
c(LON, LAT), while you can use any CRS with \code{sf}/\code{sfc}/\code{bbox} objects. See
Details and Examples in \code{\link[=oe_match]{oe_match()}}.}

\item{provider}{Which provider should be used to download the data? Available
providers can be found with the following command: \code{\link[=oe_providers]{oe_providers()}}. For
\code{\link[=oe_get]{oe_get()}} and \code{\link[=oe_match]{oe_match()}}, if \code{place} is equal to \verb{ITS Leeds}, then
\code{provider} is set equal to \code{test}. This is just for simple examples and
internal tests.}

\item{download_directory}{Directory where the files downloaded by osmextract
are stored. By default it is equal to \code{\link[=oe_download_directory]{oe_download_directory()}}.}

\item{download_if_missing}{Attempt to download the file if it cannot be
found? \code{FALSE} by default.}

\item{quiet}{Boolean. If \code{FALSE}, the function prints informative messages.
Starting from \code{sf} version
\href{https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13}{0.9.6},
if \code{quiet} is equal to \code{FALSE}, then vectortranslate operations will
display a progress bar.}

\item{...}{Extra arguments that are passed to \code{\link[=oe_match]{oe_match()}} and \code{\link[=oe_get]{oe_get()}}.
Please note that you cannot modify the argument \code{download_only}.}
}
\value{
A character vector of length one (or two) representing the path(s) of the
corresponding \code{.pbf} (and \code{.gpkg}) files.
}
\description{
This function takes a \code{place} name and it returns the path of \code{.pbf} and
\code{.gpkg} files associated with it.
}
\details{
The matching between the existing files (saved in a directory
specified by \code{download_directory} parameter) and the input \code{place} is
performed using \code{list.files}, setting a pattern equal to the basename of
the URL associated to the input \code{place}. For example, if you specify
\code{place = "Isle of Wight"}, then the input \code{place} is matched with a URL of
a \code{.osm.pbf} file (via \code{\link[=oe_match]{oe_match()}}) and the matching is performed setting a
pattern equal to the basename of that URL.

If there is no file in \code{download_directory} that can be matched with the
basename and \code{download_if_missing} parameter is equal to \code{TRUE}, then the
function tries to download and translate a new file from the chosen
provider (\code{geofabrik} is the default provider). If \code{download_if_missing}
parameter is equal to \code{FALSE} (default value), then the function stops with
an error.
}
\examples{
\dontshow{
res = file.copy(
  from = system.file("its-example.osm.pbf", package = "osmextract"),
  to = file.path(tempdir(), "test_its-example.osm.pbf"),
  overwrite = TRUE
)}
res = oe_get("ITS Leeds", quiet = TRUE, download_directory = tempdir())
oe_find("ITS Leeds", provider = "test", download_directory = tempdir())

\dontrun{
oe_find("Isle of Wight", download_directory = tempdir())
oe_find("Malta", download_if_missing = TRUE, download_directory = tempdir())
oe_find(
  "Leeds",
  provider = "bbbike",
  download_if_missing = TRUE,
  download_directory = tempdir()
)}

# Remove .pbf and .gpkg files in tempdir
# (since they may interact with other examples)
file.remove(list.files(path = tempdir(), pattern = "(pbf|gpkg)", full.names = TRUE))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{openstreetmap_fr_zones}
\alias{openstreetmap_fr_zones}
\title{An sf object of geographical zones taken from download.openstreetmap.fr}
\format{
An \code{sf} object with 903 rows and
7 columns:
\describe{
\item{id}{A unique ID for each area. It is used by \code{oe_update()}.}
\item{name}{The, usually English, long-form name of the city.}
\item{parent}{The identifier of the next larger excerpts that contains
this one, if present.}
\item{level}{An integer code between 1 and 4. Check
\url{http://download.openstreetmap.fr/polygons/} to see the hierarchical
structure of the zones. 1L correspond to the biggest areas. This is used
only for matching operations in case of spatial input.}
\item{pbf}{Link to the latest \code{.osm.pbf} file for this region.}
\item{pbf_file_size}{Size of the pbf file in bytes.}
\item{geometry}{The \code{sfg} for that geographical region, rectangular.}
}
}
\source{
\url{https://download.bbbike.org/osm/}
}
\usage{
openstreetmap_fr_zones
}
\description{
An \code{sf} object containing the URLs, names and file-sizes of the OSM
extracts stored at \url{http://download.openstreetmap.fr/}.
}
\seealso{
Other provider's-database: 
\code{\link{bbbike_zones}},
\code{\link{geofabrik_zones}}
}
\concept{provider's-database}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bbbike_zones}
\alias{bbbike_zones}
\title{An sf object of geographical zones taken from bbbike.org}
\format{
An \code{sf} object with 235 rows and
11 columns:
\describe{
\item{name}{The, usually English, long-form name of the city.}
\item{last_modified}{When was it last modified?}
\item{type}{empty}
\item{pbf_file_size}{Size of the pbf file in bytes.}
\item{base_url}{The base URL for the city.}
\item{poly_url}{The \code{.poly} file location.}
\item{pbf}{Link to the latest \code{.osm.pbf} file for this region.}
\item{level}{An integer code always equal to 3 (since the bbbike
data represent non-hierarchical geographical zones). This is used only for
matching operations in case of spatial input. The oe_* functions will
select the geographical area closest to the input place with the highest
"level". See \link{geofabrik_zones} for an example of a hierarchical structure.}
\item{geometry}{The \code{sfg} for that geographical region, rectangular.}
}
}
\source{
\url{https://download.bbbike.org/osm/}
}
\usage{
bbbike_zones
}
\description{
Start bicycle routing for... everywhere!
}
\details{
An \code{sf} object containing the URLs, names and file_size of the OSM extracts.
}
\seealso{
Other provider's-database: 
\code{\link{geofabrik_zones}},
\code{\link{openstreetmap_fr_zones}}
}
\concept{provider's-database}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{test_zones}
\alias{test_zones}
\title{An sf object of geographical zones taken from download.openstreetmap.fr}
\format{
An object of class \code{sf} (inherits from \code{data.frame}) with 2 rows and 7 columns.
}
\usage{
test_zones
}
\description{
This object represent a minimal provider's database and it should be used
only for examples and tests.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{oe_read}
\alias{oe_read}
\title{Read a .pbf or .gpkg object from file or url}
\usage{
oe_read(
  file_path,
  layer = "lines",
  ...,
  provider = NULL,
  download_directory = oe_download_directory(),
  file_size = NULL,
  force_download = FALSE,
  max_file_size = 5e+08,
  download_only = FALSE,
  skip_vectortranslate = FALSE,
  vectortranslate_options = NULL,
  osmconf_ini = NULL,
  extra_tags = NULL,
  force_vectortranslate = FALSE,
  never_skip_vectortranslate = FALSE,
  boundary = NULL,
  boundary_type = c("spat", "clipsrc"),
  quiet = FALSE
)
}
\arguments{
\item{file_path}{A URL or the path of a \code{.pbf} or \code{.gpkg} file. If a URL,
then it must be specified using HTTP/HTTPS protocol.}

\item{layer}{Which \code{layer} should be read in? Typically \code{points}, \code{lines}
(the default), \code{multilinestrings}, \code{multipolygons} or \code{other_relations}. If
you specify an ad-hoc query using the argument \code{query} (see introductory
vignette and examples), then \code{oe_get()} and \code{oe_read()} will read the layer
specified in the query and ignore \code{layer}. See also
\href{https://github.com/ropensci/osmextract/issues/122}{#122}.}

\item{...}{(Named) arguments that will be passed to \code{\link[sf:st_read]{sf::st_read()}}, like
\code{query}, \code{wkt_filter} or \code{stringsAsFactors}.  Check the introductory
vignette to understand how to create your own (SQL-like) queries.}

\item{provider}{Which provider should be used to download the data? Available
providers can be found with the following command: \code{\link[=oe_providers]{oe_providers()}}. For
\code{\link[=oe_get]{oe_get()}} and \code{\link[=oe_match]{oe_match()}}, if \code{place} is equal to \verb{ITS Leeds}, then
\code{provider} is set equal to \code{test}. This is just for simple examples and
internal tests.}

\item{download_directory}{Where to download the file containing the OSM data?
By default this is equal to \code{\link[=oe_download_directory]{oe_download_directory()}}, which is equal to
\code{\link[=tempdir]{tempdir()}} and it changes each time you restart R. You can set a
persistent \code{download_directory} by adding the following to your \code{.Renviron}
file (e.g. with \code{edit_r_environ} function in \code{usethis} package):
\verb{OSMEXT_DOWNLOAD_DIRECTORY=/path/to/osm/data}.}

\item{file_size}{How big is the file? Optional. \code{NA} by default. If it's
bigger than \code{max_file_size} and the function is run in interactive mode,
then an interactive menu is displayed, asking for permission to download
the file.}

\item{force_download}{Should the \code{.osm.pbf} file be updated if it has already
been downloaded? \code{FALSE} by default. This parameter is used to update old
\code{.osm.pbf} files.}

\item{max_file_size}{The maximum file size to download without asking in
interactive mode. Default: \code{5e+8}, half a gigabyte.}

\item{download_only}{Boolean. If \code{TRUE}, then the function only returns the
path where the matched file is stored, instead of reading it. \code{FALSE} by
default.}

\item{skip_vectortranslate}{Boolean. If \code{TRUE}, then the function skips all
vectortranslate operations and it reads (or simply returns the path) of the
\code{.osm.pbf} file. \code{FALSE} by default.}

\item{vectortranslate_options}{Options passed to the \code{\link[sf:gdal_utils]{sf::gdal_utils()}}
argument \code{options}. Set by default. Check details in the introductory
vignette and the help page of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}.}

\item{osmconf_ini}{The configuration file. See documentation at
\href{https://gdal.org/drivers/vector/osm.html}{gdal.org}. Check details in the
introductory vignette and the help page of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}. Set by
default.}

\item{extra_tags}{Which additional columns, corresponding to OSM tags, should
be in the resulting dataset? \code{NULL} by default. Check the introductory
vignette and the help pages of \code{\link[=oe_vectortranslate]{oe_vectortranslate()}} and \code{\link[=oe_get_keys]{oe_get_keys()}}.
Ignored when \code{osmconf_ini} is not \code{NULL}.}

\item{force_vectortranslate}{Boolean. Force the original \code{.pbf} file to be
translated into a \code{.gpkg} file, even if a \code{.gpkg} with the same name
already exists? \code{FALSE} by default. If tags in \code{extra_tags} match data in
previously translated \code{.gpkg} files no translation occurs (see
\href{https://github.com/ropensci/osmextract/issues/173}{#173} for details).
Check the introductory vignette and the help page of
\code{\link[=oe_vectortranslate]{oe_vectortranslate()}}.}

\item{never_skip_vectortranslate}{Boolean. This is used in case the user
passed its own \code{.ini} file or vectortranslate options (since, in those
case, it's too difficult to determine if an existing \code{.gpkg} file was
generated following the same options.)}

\item{boundary}{An \code{sf}/\code{sfc}/\code{bbox} object that will be used to create a
spatial filter during the vectortranslate operations. The type of filter
can be chosen using the argument \code{boundary_type}.}

\item{boundary_type}{A character vector of length 1 specifying the type of
spatial filter. The \code{spat} filter selects only those features that
intersect a given area, while \code{clipsrc} also clips the geometries. Check
the examples and also \href{https://gdal.org/programs/ogr2ogr.html}{here} for
more details.}

\item{quiet}{Boolean. If \code{FALSE}, the function prints informative messages.
Starting from \code{sf} version
\href{https://r-spatial.github.io/sf/news/index.html#version-0-9-6-2020-09-13}{0.9.6},
if \code{quiet} is equal to \code{FALSE}, then vectortranslate operations will
display a progress bar.}
}
\value{
An \code{sf} object or a character vector when the \code{download_only}
argument is \code{TRUE}.
}
\description{
This function is used to read a \code{.pbf} or \code{.gpkg} object from file or URL. It
is a wrapper around \code{\link[=oe_download]{oe_download()}}, \code{\link[=oe_vectortranslate]{oe_vectortranslate()}}, and
\code{\link[sf:st_read]{sf::st_read()}}, creating an easy way to download, convert, and read a \code{.pbf}
or \code{.gpkg} file. Check the introductory vignette and the help pages of the
wrapped function for more details.
}
\details{
The arguments \code{provider}, \code{download_directory}, \code{file_size},
\code{force_download}, and \code{max_file_size} are ignored if \code{file_path} points to
an existing \code{.pbf} or \code{.gpkg} file.

Please note that you cannot add any field to an existing \code{.gpkg} file using
the argument \code{extra_tags} without rerunning the vectortranslate process on
the corresponding \code{.pbf} file. On the other hand, you can extract some of
the tags in \code{other_tags} field as new columns. See examples and
\code{\link[=oe_get_keys]{oe_get_keys()}} for more details.
}
\examples{
# Read an existing .pbf file. First we need to copy a .pbf file into a
# temporary directory
file.copy(
  from = system.file("its-example.osm.pbf", package = "osmextract"),
  to = file.path(tempdir(), "its-example.osm.pbf")
)

my_pbf = file.path(tempdir(), "its-example.osm.pbf")
oe_read(my_pbf)

# Read a new layer
oe_read(my_pbf, layer = "points")

# The following example shows how to add new tags
names(oe_read(my_pbf, extra_tags = c("oneway", "ref"), quiet = TRUE))

# Read an existing .gpkg file. This file was created by oe_read
my_gpkg = file.path(tempdir(), "its-example.gpkg")
oe_read(my_gpkg)

# You cannot add any layer to an existing .gpkg file but you can extract some
# of the tags in other_tags. Check oe_get_keys() for more details.
names(oe_read(my_gpkg, extra_tags = c("maxspeed"))) # doesn't work
# Instead, use the query argument
names(oe_read(
  my_gpkg,
  quiet = TRUE,
  query =
  "SELECT *,
  hstore_get_value(other_tags, 'maxspeed') AS maxspeed
  FROM lines
  "
))

# Read from a URL
my_url = "https://github.com/ropensci/osmextract/raw/master/inst/its-example.osm.pbf"
# Please note that if you read from a URL which is not linked to one of the
# supported providers, you need to specify the provider parameter:
\dontrun{
oe_read(my_url, provider = "test", quiet = FALSE)}

# Remove .pbf and .gpkg files in tempdir
# (since they may interact with other examples)
file.remove(list.files(path = tempdir(), pattern = "(pbf|gpkg)", full.names = TRUE))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{geofabrik_zones}
\alias{geofabrik_zones}
\title{An sf object of geographical zones taken from Geofabrik}
\format{
An sf object with 430 rows and
15 columns:
\describe{
\item{id}{A unique identifier. It contains letters, numbers and potentially
the characters "-" and "/".}
\item{name}{The, usually English, long-form name of the area.}
\item{parent}{The identifier of the next larger excerpts that contains this
one, if present.}
\item{level}{An integer code between 1 and 4. If level = 1 then the zone
corresponds to one of the continents plus the Russian Federation: Africa,
Antarctica, Asia, Australia and Oceania, Central America, Europe, North
America, Russian Federation and South America. If level = 2 then the zone
corresponds to the continent's subregions (i.e. the countries, such as
Italy, Great Britain, Spain, USA, Mexico, Belize, Morocco, Peru and so on).
There are also some exceptions that correspond to the Special Sub Regions
(according to their Geofabrik definition), which are: South Africa
(includes Lesotho), Alps, Britain and Ireland, Germany + Austria +
Switzerland, US Midwest, US Northeast, US Pacific, US South, US West and
all US states. Level = 3L correspond to the subregions of each state (or
each level 2 zone). For example the West Yorkshire, which is a subregion of
England, is a level 3 zone. Finally, level = 4L correspond to the
subregions of the third level and it is mainly related to some small areas
in Germany. This field is used only for matching operations in case of
spatial input.}
\item{iso3166-1_alpha2}{A character vector of two-letter \href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}{ISO3166-1 codes}. This will be set
on the smallest extract that still fully (or mostly) contains the entity
with that code; e.g. the code "DE" will be given for the Germany extract
and not for Europe even though Europe contains Germany. If an extract
covers several countries and no per-country extracts are available (e.g.
Israel and Palestine), then several ISO codes will be given (such as "PS
IL" for "Palestine and Israel").}
\item{iso3166_2}{A character vector of usually five-character \href{https://en.wikipedia.org/wiki/ISO_3166-2}{ISO3166-2 codes}. The same rules as above
apply. Some entities have both an \emph{iso3166-1} and \emph{iso3166-2} code. For
example, the \emph{iso3166_2} code of each US State is "US - " plus the code of
the state.}
\item{pbf}{Link to the latest \code{.osm.pbf} file for this region.}
\item{bz2}{Link to the latest \code{.osm.bz2} file for this region.}
\item{shp}{Link to the latest shape file for this region.}
\item{pbf.internal}{Link to the latest \code{.osm.pbf} file with user data for
this region (requires OSM login).}
\item{history}{Link to the latest history file for this region (requires
OSM login).}
\item{taginfo}{Link to the Geofabrik taginfo instance for this region.}
\item{updates}{Link to the updates directory (append /state.txt for status
file).}
\item{geometry}{The sfc for that geographical region. These are not the
country boundaries but a buffer around countries.}
\item{pbf_file_size}{Size of the \code{.pbf} file in bytes.}
}
}
\source{
\url{https://download.geofabrik.de/}
}
\usage{
geofabrik_zones
}
\description{
An \code{sf} object containing the URLs, names and file-sizes of the OSM
extracts stored at \url{https://download.geofabrik.de/}. You can read more
details about these data at the following link:
\url{https://download.geofabrik.de/technical.html}.
}
\seealso{
Other provider's-database: 
\code{\link{bbbike_zones}},
\code{\link{openstreetmap_fr_zones}}
}
\concept{provider's-database}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-key-values.R
\name{oe_get_keys}
\alias{oe_get_keys}
\alias{oe_get_keys.default}
\alias{oe_get_keys.character}
\alias{oe_get_keys.sf}
\alias{print.oe_key_values_list}
\title{Return keys and (optionally) values stored in "other_tags" column}
\usage{
oe_get_keys(zone, layer = "lines", values = FALSE, which_keys = NULL)

\method{oe_get_keys}{default}(zone, layer = "lines", values = FALSE, which_keys = NULL)

\method{oe_get_keys}{character}(zone, layer = "lines", values = FALSE, which_keys = NULL)

\method{oe_get_keys}{sf}(zone, layer = "lines", values = FALSE, which_keys = NULL)

\method{print}{oe_key_values_list}(x, n = getOption("oe_max_print_keys", 10L), ...)
}
\arguments{
\item{zone}{An \code{sf} object with an \code{other_tags} field or a character vector
(of length 1) that can be linked to or pointing to a \code{.osm.pbf} or \code{.gpkg}
file with an \code{other_tags} field. Character vectors are linked to \code{.osm.pbf}
files using \code{oe_find()}.}

\item{layer}{Which \code{layer} should be read in? Typically \code{points}, \code{lines}
(the default), \code{multilinestrings}, \code{multipolygons} or \code{other_relations}. If
you specify an ad-hoc query using the argument \code{query} (see introductory
vignette and examples), then \code{oe_get()} and \code{oe_read()} will read the layer
specified in the query and ignore \code{layer}. See also
\href{https://github.com/ropensci/osmextract/issues/122}{#122}.}

\item{values}{Logical. If \code{TRUE}, then function returns the keys and the
corresponding values, otherwise only the keys. Defaults to \code{FALSE. }}

\item{which_keys}{Character vector used to subset only some keys and
corresponding values. Ignored if \code{values} is \code{FALSE}. See examples.}

\item{x}{object of class \code{oe_key_values_list}}

\item{n}{Maximum number of keys (and corresponding values) to print; can be
set globally by \code{options(oe_max_print_keys=...)}. Default value is 10.}

\item{...}{Ignored.}
}
\value{
If the argument \code{values} is \code{FALSE} (the default), then the function
returns a character vector with the names of all keys stored in the
\code{other_tags} field. If \code{values} is \code{TRUE}, then the function returns named
list which stores all keys and the corresponding values. In the latter
case, the returned object has class \code{oe_key_values_list} and we defined an
ad-hoc printing method. See Details.
}
\description{
This function returns the OSM keys and (optionally) the values stored in the
\code{other_tags} field. See Details. In both cases, the keys are sorted according
to the number of occurrences, which means that the most common keys are
stored first.
}
\details{
OSM data are typically documented using several
\href{https://wiki.openstreetmap.org/wiki/Tags}{\code{tags}}, i.e. pairs of two
items, namely a \code{key} and a \code{value}. The conversion between \code{.osm.pbf} and
\code{.gpkg} formats is governed by a \code{CONFIG} file that lists which tags must
be explicitly added to the \code{.gpkg} file. All the other keys are
automatically stored using an \code{other_tags} field with a syntax compatible
with the PostgreSQL HSTORE type. See
\href{https://gdal.org/drivers/vector/osm.html#driver-capabilities}{here} for
more details.

When the argument \code{values} is \code{TRUE}, then the function returns a named list
of class \code{oe_key_values_list} that, for each key, summarises the
corresponding values. The key-value pairs are stored using the following
format:
\verb{list(key1 = c("value1", "value1", "value2", ...), key2 = c("value1", ...) ...)}.
We decided to implement an ad-hoc method for printing objects of class
\code{oe_key_values_list} using the following structure:\preformatted{
key1 = {#value1 = n1; #value2 = n2; #value3 = n3, ...}
key2 = {#value1 = n1; #value2 = n2; ...}
key3 = {#value1 = n1}
...
}
where \code{n1} denotes the number of times that value1 is repeated, \code{n2} denotes
the number of times that value2 is repeated and so on. Also the values are
listed according to the number of occurrences in decreasing order. By
default, the function prints only the ten most common keys, but the number
can be adjusted using the option \code{oe_max_print_keys}.

Finally, the \code{hstore_get_value()} function can be used inside the \code{query}
argument in \code{oe_get()} to extract one particular tag from an existing file.
Check the introductory vignette and see examples.
}
\examples{
# Copy ITS file to tempdir so that the examples do not require internet
# connection. You can skip the next few lines (and start directly with
# oe_get_keys) when running the examples ocally.
its_pbf = file.path(tempdir(), "test_its-example.osm.pbf")
file.copy(
  from = system.file("its-example.osm.pbf", package = "osmextract"),
  to = its_pbf,
  overwrite = TRUE
)

# Get keys
oe_get_keys("ITS Leeds")

# Get keys and values
oe_get_keys("ITS Leeds", values = TRUE)

# Subset some keys
oe_get_keys("ITS Leeds", values = TRUE, which_keys = c("surface", "lanes"))

# Print all (non-NA) values for a given set of keys
oe_get_keys("ITS Leeds", values = TRUE)["surface"]

# Get keys from an existing sf object
\dontrun{
its = oe_get("ITS Leeds")
oe_get_keys(its, values = TRUE)}

# Get keys from a character vector pointing to a file (might be faster than
# reading the complete file)
its_path = oe_get("ITS Leeds", download_only = TRUE, download_directory = tempdir())
oe_get_keys(its_path, values = TRUE)

# Add a key to an existing .gpkg file without repeating the
# vectortranslate operations
\dontrun{
colnames(its)
colnames(oe_read(
  its_path,
  query = "SELECT *, hstore_get_value(other_tags, 'oneway') AS oneway FROM lines",
  quiet = TRUE
))}

# Remove .pbf and .gpkg files in tempdir
# (since they may interact with other examples)
file.remove(list.files(path = tempdir(), pattern = "(pbf|gpkg)", full.names = TRUE))
}
\seealso{
\code{oe_vectortranslate()}
}
