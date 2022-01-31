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

<!-- README.md is generated from README.Rmd. Please edit that file -->
antanym
=======

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/antanym/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/antanym/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/antanym/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/antanym?branch=master)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://badges.ropensci.org/198_status.svg)](https://github.com/ropensci/onboarding/issues/198)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/antanym)](http://cran.r-project.org/web/packages/antanym) 
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/antanym)
<!-- badges: end -->

Overview
--------

This R package provides easy access to Antarctic geographic place name information, and tools for working with those names.

The authoritative source of place names in Antarctica is the Composite Gazetteer of Antarctica (CGA), which is produced by the Scientific Committee on Antarctic Research (SCAR). The CGA consists of approximately 37,000 names corresponding to 19,000 distinct features. It covers features south of 60 °S, including terrestrial and undersea or under-ice.

There is no single naming authority responsible for place names in Antarctica because it does not fall under the sovereignty of any one nation. In general, individual countries have administrative bodies that are responsible for their national policy on, and authorisation and use of, Antarctic names. The CGA is a compilation of place names that have been submitted by representatives of national names committees from 22 countries.

The composite nature of the CGA means that there are often multiple names associated with a given feature. Consider using the `an_preferred()` function for resolving a single name per feature.

For more information, see the [CGA home page](http://data.aad.gov.au/aadc/gaz/scar/). The CGA was begun in 1992. Since 2008, Italy and Australia have jointly managed the CGA, the former taking care of the editing, the latter maintaining the database and website. The SCAR [Standing Committee on Antarctic Geographic Information (SCAGI)](http://www.scar.org/data-products/scagi) coordinates the project. This R package is a product of the SCAR [Expert Group on Antarctic Biodiversity Informatics](http://www.scar.org/ssg/life-sciences/eg-abi) and SCAGI.

### Citing

The SCAR Composite Gazetteer of Antarctica is made available under a CC-BY license. If you use it, please cite it:

> Composite Gazetteer of Antarctica, Scientific Committee on Antarctic Research. GCMD Metadata (<http://gcmd.nasa.gov/records/SCAR_Gazetteer.html>)

Installing
----------

``` r
install.packages("remotes")
remotes::install_github("ropensci/antanym")
```

Example usage
-------------

Start by fetching the names data from the host server. Here we use a temporary cache so that we can re-load it later in the session without needing to re-download it:

``` r
library(antanym)
g <- an_read(cache = "session")
```

How many names do we have in total?

``` r
nrow(g)
#> [1] 37661
```

Corresponding to how many distinct features?

``` r
length(unique(g$scar_common_id))
#> [1] 19601
```

Find names starting with "Slom":

``` r
an_filter(g, query = "^Slom")[, c("place_name", "longitude", "latitude")]
#> # A tibble: 3 x 3
#>   place_name     longitude latitude
#>   <chr>              <dbl>    <dbl>
#> 1 Sloman Glacier     -68.6    -67.7
#> 2 Sloman Glacier     -68.6    -67.7
#> 3 Slomer Cove        -59.4    -63.8
```

Find islands within 20km of 100 °E, 66 °S:

``` r
nms <- an_near(an_filter(g, feature_type = "Island"), loc = c(100, -66), max_distance = 20)

## or equivalently, using the pipe operator
nms <- g %>% an_filter(feature_type = "Island") %>% an_near(loc = c(100, -66), max_distance = 20)

nms[, c("place_name", "longitude", "latitude")]
#> # A tibble: 3 x 3
#>   place_name    longitude latitude
#>   <chr>             <dbl>    <dbl>
#> 1 Foster Island       100    -66.1
#> 2 Severnyj holm       100    -66.1
#> 3 Foster Island       100    -66.1
```

Resolving multiple names per feature
------------------------------------

As noted above, the CGA is a composite gazetteer and so there are often multiple names associated with a given feature. For example, we can see all names associated with feature 1589 (Booth Island) and the country of origin of each name:

``` r
an_filter(g, feature_ids = 1589)[, c("place_name", "origin")]
#> # A tibble: 7 x 2
#>   place_name   origin                  
#>   <chr>        <chr>                   
#> 1 Booth, isla  Argentina               
#> 2 Wandel, Ile  Belgium                 
#> 3 Booth, Isla  Chile                   
#> 4 Boothinsel   Germany                 
#> 5 Booth Island United Kingdom          
#> 6 Booth Island Russia                  
#> 7 Booth Island United States of America
```

The `an_preferred` function can help with finding one name per feature. It takes an `origin` parameter that specifies one or more preferred naming authorities (countries or organisations). For features that have multiple names (e.g. have been named by multiple countries) a single name will be chosen, preferring names from the specified naming authorities where possible.

We start with 37661 names in the full CGA, corresponding to 19601 distinct features. Choose one name per feature, preferring the Polish name where there is one, and the German name as a second preference:

``` r
g <- an_preferred(g, origin = c("Poland", "Germany"))
```

Now we have 19601 names in our data frame, corresponding to the same 19601 distinct features.

Name suggestions
----------------

Antanym includes an experimental function that will suggest which features might be best to to name on a given map. These suggestions are based on maps prepared by expert cartographers, and the features that were explicitly named on those maps.

See the package vignette and the `an_suggest` function for more information.

Recent changes
--------------

### antanym 0.4.0

General revisions following rOpenSci review. Note several breaking changes:

-   `an_read` now takes a `cache` parameter instead of `cache_directory` (and now can have special values "session" and "persistent")
-   `an_filter` and `an_suggest` now take an `origin` parameter that replaces the previous `origin_country` and `cga_source` parameters
-   the default data structure (returned by `an_read(..., simplified = TRUE)` no longer contains the "country\_name" or "cga\_source\_gazetteer columns, but if needed these are available via `an_read(..., simplified = FALSE)`

Other map examples
------------------

A [leaflet app](https://australianantarcticdatacentre.github.io/antanym-demo/leaflet.html) using Mercator projection and clustered markers for place names.

<a href="https://australianantarcticdatacentre.github.io/antanym-demo/leaflet.html"><img src="vignettes/README-leaflet.png" width="40%" /></a>

And a similar example using a [polar stereographic projection](https://australianantarcticdatacentre.github.io/antanym-demo/leafletps.html).

<a href="https://australianantarcticdatacentre.github.io/antanym-demo/leafletps.html"><img src="vignettes/README-leafletps.png" width="40%" /></a>

See the [antanym-demo](https://github.com/AustralianAntarcticDataCentre/antanym-demo) repository for the source code of these examples.

Other packages
--------------

The [geonames package](https://cran.r-project.org/package=geonames) also provides access to geographic place names, including from the SCAR Composite Gazetteer. If you need *global* place name coverage, geonames may be a better option. However, the composite nature of the CGA is not particularly well suited to geonames, and at the time of writing the geonames database did not include the most current version of the CGA. The geonames package requires a login for some functionality, and because it makes calls to api.geonames.org it isn't easily used while offline.

[![ropensci\_footer](https://ropensci.org/public_images/scar_footer.png)](https://ropensci.org)
# antanym 0.4.0

General revisions following rOpenSci review. Note several breaking changes:

  - `an_read` now takes a `cache` parameter instead of `cache_directory` (and now can have special values "session" and "persistent")
  - `an_filter` and `an_suggest` now take an `origin` parameter that replaces the previous `origin_country` and `cga_source` parameters
  - the default data structure (returned by `an_read(..., simplified = TRUE)` no longer contains the "country_name" or "cga_source_gazetteer columns, but if needed these are available via `an_read(..., simplified = FALSE)`
  

# antanym 0.3.0

* Added a `NEWS.md` file to track changes to the package.

# Contributions

Suggestions, requests, bug reports, and code are welcome. Please open an [issue](https://github.com/SCAR/antanym/issues) or [fork this repository](https://help.github.com/articles/fork-a-repo/) and submit a pull request.

# Code of conduct

Maintainers and contributors must follow this repository's [code of conduct](CONDUCT.md).
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "vignettes/README-"
)
```

# antanym

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/antanym/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/antanym/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/antanym/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/antanym?branch=master)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://badges.ropensci.org/198_status.svg)](https://github.com/ropensci/onboarding/issues/198)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/antanym)](http://cran.r-project.org/web/packages/antanym) 
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/antanym)
<!-- badges: end -->

## Overview

This R package provides easy access to Antarctic geographic place name information, and tools for working with those names.

The authoritative source of place names in Antarctica is the Composite Gazetteer of Antarctica (CGA), which is produced by the Scientific Committee on Antarctic Research (SCAR). The CGA consists of approximately 37,000 names corresponding to 19,000 distinct features. It covers features south of 60 &deg;S, including terrestrial and undersea or under-ice.

There is no single naming authority responsible for place names in Antarctica because it does not fall under the sovereignty of any one nation. In general, individual countries have administrative bodies that are responsible for their national policy on, and authorisation and use of, Antarctic names. The CGA is a compilation of place names that have been submitted by representatives of national names committees from 22 countries.

The composite nature of the CGA means that there are often multiple names associated with a given feature. Consider using the `an_preferred()` function for resolving a single name per feature.

For more information, see the [CGA home page](http://data.aad.gov.au/aadc/gaz/scar/). The CGA was begun in 1992. Since 2008, Italy and Australia have jointly managed the CGA, the former taking care of the editing, the latter maintaining the database and website. The SCAR [Standing Committee on Antarctic Geographic Information (SCAGI)](http://www.scar.org/data-products/scagi) coordinates the project. This R package is a product of the SCAR [Expert Group on Antarctic Biodiversity Informatics](http://www.scar.org/ssg/life-sciences/eg-abi) and SCAGI.


### Citing

The SCAR Composite Gazetteer of Antarctica is made available under a CC-BY license. If you use it, please cite it:

> Composite Gazetteer of Antarctica, Scientific Committee on Antarctic Research. GCMD Metadata (http://gcmd.nasa.gov/records/SCAR_Gazetteer.html)

## Installing


```{r, eval = FALSE}
install.packages("remotes")
remotes::install_github("ropensci/antanym")
```

## Example usage

Start by fetching the names data from the host server. Here we use a temporary cache so that we can re-load it later in the session without needing to re-download it:

```{r message = FALSE, warning = FALSE}
library(antanym)
g <- an_read(cache = "session")
```

How many names do we have in total?
```{r message = FALSE, warning = FALSE}
nrow(g)
```

Corresponding to how many distinct features?
```{r message = FALSE, warning = FALSE}
length(unique(g$scar_common_id))
```

Find names starting with "Slom":
```{r message = FALSE, warning = FALSE}
an_filter(g, query = "^Slom")[, c("place_name", "longitude", "latitude")]
```

Find islands within 20km of 100 &deg;E, 66 &deg;S:
```{r message = FALSE, warning = FALSE}
nms <- an_near(an_filter(g, feature_type = "Island"), loc = c(100, -66), max_distance = 20)

## or equivalently, using the pipe operator
nms <- g %>% an_filter(feature_type = "Island") %>% an_near(loc = c(100, -66), max_distance = 20)

nms[, c("place_name", "longitude", "latitude")]
```

## Resolving multiple names per feature

As noted above, the CGA is a composite gazetteer and so there are often multiple names associated with a given feature. For example, we can see all names associated with feature 1589 (Booth Island) and the country of origin of each name:

```{r message = FALSE, warning = FALSE}
an_filter(g, feature_ids = 1589)[, c("place_name", "origin")]
```

The `an_preferred` function can help with finding one name per feature. It takes an `origin` parameter that specifies one or more preferred naming authorities (countries or organisations). For features that have multiple names (e.g. have been named by multiple countries) a single name will be chosen, preferring names from the specified \code{origin} naming authorities where possible.

We start with `r nrow(g)` names in the full CGA, corresponding to `r length(unique(g$scar_common_id))` distinct features. Choose one name per feature, preferring the Polish name where there is one, and the German name as a second preference:
```{r}
g <- an_preferred(g, origin = c("Poland", "Germany"))
```

Now we have `r nrow(g)` names in our data frame, corresponding to the same `r length(unique(g$scar_common_id))` distinct features.

## Name suggestions

Antanym includes an experimental function that will suggest which features might be best to to name on a given map. These suggestions are based on maps prepared by expert cartographers, and the features that were explicitly named on those maps.

See the package vignette and the `an_suggest` function for more information.

## Recent changes

### antanym 0.4.0

General revisions following rOpenSci review. Note several breaking changes:

  - `an_read` now takes a `cache` parameter instead of `cache_directory` (and now can have special values "session" and "persistent")
  - `an_filter` and `an_suggest` now take an `origin` parameter that replaces the previous `origin_country` and `cga_source` parameters
  - the default data structure (returned by `an_read(..., simplified = TRUE)` no longer contains the "country_name" or "cga_source_gazetteer columns, but if needed these are available via `an_read(..., simplified = FALSE)`


## Other map examples

A [leaflet app](https://australianantarcticdatacentre.github.io/antanym-demo/leaflet.html) using Mercator projection and clustered markers for place names.

<a href="https://australianantarcticdatacentre.github.io/antanym-demo/leaflet.html"><img src="vignettes/README-leaflet.png" width="40%" /></a>

And a similar example using a [polar stereographic projection](https://australianantarcticdatacentre.github.io/antanym-demo/leafletps.html).

<a href="https://australianantarcticdatacentre.github.io/antanym-demo/leafletps.html"><img src="vignettes/README-leafletps.png" width="40%" /></a>

See the [antanym-demo](https://github.com/AustralianAntarcticDataCentre/antanym-demo) repository for the source code of these examples.


## Other packages

The [geonames package](https://cran.r-project.org/package=geonames) also provides access to geographic place names, including from the SCAR Composite Gazetteer. If you need *global* place name coverage, geonames may be a better option. However, the composite nature of the CGA is not particularly well suited to geonames, and at the time of writing the geonames database did not include the most current version of the CGA. The geonames package requires a login for some functionality, and because it makes calls to api.geonames.org it isn't easily used while offline.

[![ropensci_footer](https://ropensci.org/public_images/scar_footer.png)](https://ropensci.org)
---
title: "antanym"
author: "Ben Raymond, Michael Sumner"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{antanym}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r chunkopts, eval=TRUE, echo=FALSE}
knitr::opts_chunk$set(eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, tidy=FALSE, cache=FALSE, include=TRUE, dpi=72, fig.width=9, fig.height=6, fig.align="center", results="markup")
```

## Overview

This R package provides easy access to Antarctic geographic place name information, and tools for working with those names.

The authoritative source of place names in Antarctica is the Composite Gazetteer of Antarctica (CGA), which is produced by the Scientific Committee on Antarctic Research (SCAR). The CGA consists of approximately 37,000 names corresponding to 19,000 distinct features. It covers features south of 60 &deg;S, including terrestrial and undersea or under-ice features.

There is no single naming authority responsible for place names in Antarctica because it does not fall under the sovereignty of any one nation. In general, individual countries have administrative bodies that are responsible for their national policy on, and authorisation and use of, Antarctic names. The CGA is a compilation of place names that have been submitted by representatives of national names committees from 22 countries.

The composite nature of the CGA means that there may be multiple names associated with a given feature. Consider using the `an_preferred()` function for resolving a single name per feature.

For more information, see the [CGA home page](http://data.aad.gov.au/aadc/gaz/scar/). The CGA was begun in 1992. Since 2008, Italy and Australia have jointly managed the CGA, the former taking care of the editing, the latter maintaining the database and website. The SCAR [Standing Committee on Antarctic Geographic Information (SCAGI)](http://www.scar.org/data-products/scagi) coordinates the project. This R package is a product of the SCAR [Expert Group on Antarctic Biodiversity Informatics](http://www.scar.org/ssg/life-sciences/eg-abi) and SCAGI.


### Citing

The SCAR Composite Gazetteer of Antarctica is made available under a CC-BY license. If you use it, please cite it:

> Composite Gazetteer of Antarctica, Scientific Committee on Antarctic Research. GCMD Metadata (http://gcmd.nasa.gov/records/SCAR_Gazetteer.html)

## Installing


```{r eval = FALSE}
install.packages("remotes")
remotes::install_github("ropensci/antanym")
```

## Usage

Start by fetching the names data from the host server. Here we use a temporary cache so that we can re-load it later in the session without needing to re-download it:

```{r}
library(antanym)
g <- an_read(cache = "session")
```

If you want to be able to work offline, you can cache the data to a directory that will persist between R sessions. The location of this directory is determined by `rappdirs::user_cache_dir()` (see `an_cache_directory("persistent")` if you want to know what it is):

```{r eval = FALSE}
g <- an_read(cache = "persistent")
```

If you prefer working with `sp` objects, then this will return a `SpatialPointsDataFrame`:
```{r as_sp, eval = FALSE}
gsp <- an_read(sp = TRUE, cache = "session")
```

### Data summary

How many names do we have in total?
```{r}
nrow(g)
```

Corresponding to how many distinct features?
```{r}
length(unique(g$scar_common_id))
```

We can get a list of the feature types in the data:
```{r}
head(an_feature_types(g), 10)
```

### Data structure

A description of the data structure can be found in the `an_load` function help, (the same information is available via `an_cga_metadata`). By default, the gazetteer data consists of these columns:
```{r echo = FALSE}
knitr::kable(an_cga_metadata())
```

A few additional (but probably less useful) columns can be included by calling `an_read(..., simplified = FALSE)`. Consult the `an_read` documentation for more information.

### Finding names

The `an_filter` function provides a range of options to find the names you are interested in.


#### Basic searching

A simple search for any place name containing the word 'William':
```{r}
an_filter(g, query = "William")
```

We can filter according to the country or organisation that issued the name. Which bodies (countries or organisations) provided the names in our data?
```{r}
an_origins(g)
```

Find names containing "William" and originating from Australia or the USA:
```{r}
an_filter(g, query = "William", origin = "Australia|United States of America")
```

Compound names can be slightly trickier. This search will return no matches, because the actual place name is 'William Scoresby Archipelago':
```{r}
an_filter(g, query = "William Archipelago")
```

To get around this, we can split the search terms so that each is matched separately (only names matching both "William" and "Archipelago" will be returned):
```{r}
an_filter(g, query = c("William", "Archipelago"))
```

Or a simple text query but additionally constrained by feature type:
```{r}
an_filter(g, query = "William", feature_type = "Archipelago")
```

Or we can use a regular expression:
```{r}
an_filter(g, query = "William .* Archipelago")
```

### Regular expressions

Regular expressions also allow more complex searches. For example, names matching "West" or "East":
```{r}
an_filter(g, query = "West|East")
```

Names **starting** with "West" or "East":
```{r}
an_filter(g, query = "^(West|East)")
```

Names with "West" or "East" appearing as complete words in the name ("\\b" matches a word boundary; see `help("regex")`):
```{r}
an_filter(g, query = "\\b(West|East)\\b")
```

#### Spatial searching

We can filter by spatial extent:
```{r}
an_filter(g, extent = c(100, 120, -70, -65))
```

The extent can also be passed as Spatial or Raster object, in which case its extent will be used:
```{r}
my_sp <- sp::SpatialPoints(cbind(c(100, 120), c(-70, -65)))
an_filter(g, extent = my_sp)
```

Searching by proximity to a given point, for example islands within 20km of 100 &deg;E, 66 &deg;S:
```{r}
an_near(an_filter(g, feature_type = "Island"), loc = c(100, -66), max_distance = 20)
```

### Resolving multiple names

As noted above, the CGA is a composite gazetteer and so there may be multiple names associated with a given feature. 

Find all names associated with feature 1589 (Booth Island) and show the country of origin of each name:
```{r}
an_filter(g, feature_ids = 1589)[, c("place_name", "origin")]
```

The `an_preferred` function can help with finding one name per feature. It takes an `origin` parameter that specifies one or more preferred naming authorities (countries or organisations). For features that have multiple names (e.g. have been named by multiple countries) a single name will be chosen, preferring names from the specified \code{origin} naming authorities where possible.

We start with `r nrow(g)` names in the full CGA, corresponding to `r length(unique(g$scar_common_id))` distinct features. Choose one name per feature, preferring the Polish name where there is one, and the German name as a second preference:
```{r}
g <- an_preferred(g, origin = c("Poland", "Germany"))
```

Now we have `r nrow(g)` names in our data frame, corresponding to the same `r length(unique(g$scar_common_id))` distinct features.

Features that have not been named by either of our preferred countries will have a name chosen from another country, with those countries in random order of preference.


### Using the pipe operator

All of the above functionality can be achieved in a piped workflow, if that's your preference, e.g.:

```{r}
nms <- g %>% an_filter(feature_type = "Island") %>% an_near(loc = c(100, -66), max_distance = 20)
nms[, c("place_name", "longitude", "latitude")]
```

### Working with sp

If you prefer to work with Spatial objects, the gazetteer data can be converted to a SpatialPointsDataFrame when loaded:

```{r}
gsp <- an_read(cache = "session", sp = TRUE)
```

And the above functions work in the same way, for example:
```{r}
my_sp <- sp::SpatialPoints(cbind(c(100, 120), c(-70, -65)))
an_filter(gsp, extent = my_sp)
```

### Selecting names for plotting

Let's say we are preparing a figure of the greater Prydz Bay region (60-90 &deg;E, 65-70 &deg;S), to be shown at 80mm x 80mm in size (this is approximately a 1:10M scale map). Let's plot all of the place names in this region:

```{r map0}
my_longitude <- c(60, 90)
my_latitude <- c(-70, -65)

this_names <- an_filter(g, extent = c(my_longitude, my_latitude))

if (!requireNamespace("rworldmap", quietly = TRUE)) {
  message("Skipping map figure - install the rworldmap package to see it.")
} else {
  library(rworldmap)
  map <- getMap(resolution = "low")
  plot(map, xlim = my_longitude + c(-7, 4), ylim = my_latitude, col = "grey50")
  ## allow extra xlim space for labels
  points(this_names$longitude, this_names$latitude, pch = 21, bg = "green", cex = 2)
  ## alternate the positions of labels to reduce overlap
  pos <- rep(c(1, 2, 3, 4), ceiling(nrow(this_names)/4))
  pos[order(this_names$longitude)] <- pos[1:nrow(this_names)]
  text(this_names$longitude, this_names$latitude, labels = this_names$place_name, pos = pos)
}
```

Oooooo-kay. That's not ideal.

Antanym includes an experimental function that will suggest which features might be best to add names to on a given map. These suggestions are based on maps prepared by expert cartographers, and the features that were explicitly named on those maps. We can ask for suggested names to show on our example map:

```{r}
suggested <- an_suggest(g, map_extent = c(my_longitude, my_latitude), map_dimensions = c(80, 80))
```

Plot the top ten suggested names purely by score:

```{r map1}
this_names <- head(suggested, 10)

if (!requireNamespace("rworldmap", quietly = TRUE)) {
  message("Skipping map figure - install the rworldmap package to see it.")
} else {
  plot(map, xlim = my_longitude + c(-7, 4), ylim = my_latitude, col = "grey50")
  points(this_names$longitude, this_names$latitude, pch = 21, bg = "green", cex = 2)
  pos <- rep(c(1, 2, 3, 4), ceiling(nrow(this_names)/4))
  pos[order(this_names$longitude)] <- pos[1:nrow(this_names)]
  text(this_names$longitude, this_names$latitude, labels = this_names$place_name, pos = pos)
}
```

Or the ten best suggested names considering both score and spatial coverage:

```{r map2}
this_names <- an_thin(suggested, n = 10)

if (!requireNamespace("rworldmap", quietly = TRUE)) {
  message("Skipping map figure - install the rworldmap package to see it.")
} else {
  plot(map, xlim = my_longitude + c(-7, 4), ylim = my_latitude, col = "grey50")
  points(this_names$longitude, this_names$latitude, pch = 21, bg = "green", cex = 2)
  pos <- rep(c(1, 2, 3, 4), ceiling(nrow(this_names)/4))
  pos[order(this_names$longitude)] <- pos[1:nrow(this_names)]
  text(this_names$longitude, this_names$latitude, labels = this_names$place_name, pos = pos)
}
```

## Other map examples

A [leaflet app](https://australianantarcticdatacentre.github.io/antanym-demo/leaflet.html) using Mercator projection and clustered markers for place names.

<a href="https://australianantarcticdatacentre.github.io/antanym-demo/leaflet.html"><img src="README-leaflet.png" width="40%" /></a>

And a similar example using a [polar stereographic projection](https://australianantarcticdatacentre.github.io/antanym-demo/leafletps.html).

<a href="https://australianantarcticdatacentre.github.io/antanym-demo/leafletps.html"><img src="README-leafletps.png" width="40%" /></a>

See the [antanym-demo](https://github.com/AustralianAntarcticDataCentre/antanym-demo) repository for the source code of these examples.


## Future directions

Antanym currently only provides information from the SCAR CGA. This does not cover other features that may be of interest to Antarctic researchers, such as those on subantarctic islands or features that have informal names not registered in the CGA. Antanym may be expanded to cover extra gazetteers containing such information, at a later date.

## Other packages

The [geonames package](https://cran.r-project.org/package=geonames) also provides access to geographic place names, including from the SCAR Composite Gazetteer. If you need *global* place name coverage, geonames may be a better option. However, the composite nature of the CGA is not particularly well suited to geonames, and at the time of writing the geonames database did not include the most current version of the CGA. The geonames package requires a login for some functionality, and because it makes calls to api.geonames.org it isn't easily used while offline.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gazetteers.R
\name{an_gazetteers}
\alias{an_gazetteers}
\title{The place name gazetteers available}
\usage{
an_gazetteers(gaz)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: (optional) as returned by \code{\link{an_read}}, \code{\link{an_preferred}}, or \code{\link{an_filter}}}
}
\value{
character vector. If \code{gaz} was provided, this will be a list of all gazetteers present in \code{gaz}. Otherwise, it will be a list of all gazetteers available through the antanym package
}
\description{
Return a character vector that lists all of the gazetteers present in the \code{gaz} data, or (if \code{gaz} was not provided) all of the gazetteers available through the antanym package. Currently only one gazetteer is available: the Composite Gazetteer of Antarctica.
}
\examples{

an_gazetteers()

\dontrun{
 g <- an_read(cache = "session")
 an_gazetteers(g)
}

}
\seealso{
\code{\link{an_read}}, \code{\link{an_filter}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/suggest.R
\name{an_suggest}
\alias{an_suggest}
\title{Suggest names for a map (experimental)}
\usage{
an_suggest(gaz, map_scale, map_extent, map_dimensions)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: as returned by \code{\link{an_read}}, \code{\link{an_preferred}}, or \code{\link{an_filter}}}

\item{map_scale}{numeric: the scale of the map (e.g. 20e6 for a 1:20M map). If \code{map_scale} is not provided, it will be estimated from \code{map_extent} and \code{map_dimensions}}

\item{map_extent}{vector of c(longitude_min, longitude_max, latitude_min, latitude_max): the extent of the area for which name suggestions are sought. This is required if \code{map_scale} is not provided, and optional if \code{map_scale} is provided (if \code{map_extent} is provided in this situation then the \code{gaz} data frame will be filtered to this extent before the suggestion algorithm is applied; otherwise all names in \code{gaz} will be considered). \code{map_extent} can also be passed as a raster Extent object, a Raster object (in which case its extent will be used), a Spatial object (in which case the bounding box of the object will be used as the extent), or a matrix (in which case it will be assumed to be the output of \code{sp::bbox})}

\item{map_dimensions}{numeric: 2-element numeric giving width and height of the map, in mm. Not required if \code{map_scale} is provided}
}
\value{
data.frame of names with a "score" column added. Score values range from 0 to 1. The data frame will be sorted in descending score order. Names with higher scores are those that are suggested as the most suitable for display.
}
\description{
Features are given a suitability score based on maps prepared by expert cartographers. Data were tabulated from a collection of such maps, indicating for each feature whether it was named on a given map, along with details (such as scale) of the map. These data are used as the basis of a recommendation algorithm, which suggests the best features to name on a map given its properties (extent and scale). This is an experimental function and currently only implemented for \code{map_scale} values of 10 million or larger.
}
\examples{
\dontrun{
 g <- an_read(cache = "session")

 ## get a single name per feature, preferring the
 ##  Australian name where there is one
 g <- an_preferred(g, origin = "Australia")

 ## suggested names for a 100x100 mm map covering 60-90E, 70-60S
 ##  (this is about a 1:12M scale map)
 suggested <- an_suggest(g, map_extent = c(60, 90, -70, -60), map_dimensions = c(100, 100))
 head(suggested, 20) ## top 20 names

 ## an equivalent result can be achieved by supplying map scale and extent
 suggested <- an_suggest(g, map_scale = 12e6, map_extent = c(60, 90, -70, -60))
}
}
\seealso{
\code{\link{an_read}} \code{\link{an_thin}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{an_mapscale}
\alias{an_mapscale}
\title{Calculate approximate map scale}
\usage{
an_mapscale(map_dimensions, map_extent)
}
\arguments{
\item{map_dimensions}{numeric: 2-element numeric giving width and height of the map, in mm}

\item{map_extent}{vector of c(longitude_min, longitude_max, latitude_min, latitude_max): the geographic extent of the map. \code{map_extent} can also be passed as a raster Extent object, a Raster object (in which case its extent will be used), a Spatial object (in which case the bounding box of the object will be used as the extent), or a matrix (in which case it will be assumed to be the output of \code{sp::bbox})}
}
\value{
numeric
}
\description{
Calculate approximate map scale
}
\examples{
## an A3-sized map of the Southern Ocean (1:20M)
an_mapscale(map_dimensions = c(400, 570), map_extent = c(-180, 180, -90, -40))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cga_preferred_names.R
\name{an_preferred}
\alias{an_preferred}
\title{Find one name per feature in the Composite Gazetteer}
\usage{
an_preferred(gaz, origin, unmatched = "random")
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: as returned by \code{\link{an_read}} or \code{\link{an_filter}}}

\item{origin}{character: vector of preferred name origins (countries or organisations), in order of preference. If a given feature has been named by one of these bodies, this place name will be chosen. If the feature in question has not been given a name by any of these bodies, a place name given by another body will be chosen, with preference according to the \code{unmatched} parameter. For valid \code{origin} values, see \code{\link{an_origins}}}

\item{unmatched}{string: how should names be chosen for features that have not been been named by one of the preferred \code{origin} bodies? Valid values are "random" (the non-preferred originating bodies will be randomly ordered) or "count" (the non-preferred originating bodies will be ordered by their number of entries, with the largest first)}
}
\value{
data.frame of results
}
\description{
The Composite Gazetteer of Antarctica is a compilation of place names provided by different countries and organisations. The composite nature of the CGA means that there may be multiple names associated with a single feature. The \code{an_preferred} function can be used to resolve a single name per feature. Provide one or more \code{origin} entries and the input \code{gaz} will be filtered to a single name per feature. For features that have multiple names (e.g. have been named by multiple countries) a single name will be chosen, preferring names from the specified \code{origin} bodies where possible.
}
\examples{
\dontrun{
 g <- an_read(cache = "session")

 ## get a single name per feature, preferring the
 ##  Polish name where there is one
 pnames <- an_preferred(g, origin = "Poland")

 ## names starting with "Sm", preferring US names then
 ##  Australian ones if available
 g \%>\% an_filter("^Sm") \%>\%
       an_preferred(origin = c("United States of America", "Australia"))
}

}
\references{
\url{https://data.aad.gov.au/aadc/gaz/scar/}, \url{https://www.scar.org/data-products/place-names/}
}
\seealso{
\code{\link{an_read}}, \code{\link{an_origins}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{an_read}
\alias{an_read}
\title{Load Antarctic place name data}
\usage{
an_read(
  gazetteers = "all",
  sp = FALSE,
  cache,
  refresh_cache = FALSE,
  simplified = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{gazetteers}{character: vector of gazetteers to load. For the list of available gazetteers, see \code{\link{an_gazetteers}}. Use \code{gazetteers = "all"} to load all available gazetteers. Currently only one gazetteer is available: the Composite Gazetteer of Antarctica}

\item{sp}{logical: if FALSE return a data.frame; if TRUE return a SpatialPointsDataFrame}

\item{cache}{string: the gazetteer data can be cached locally, so that it can be used offline later. Valid values are \code{"session"}, \code{"persistent"}, or a directory name. Specifying \code{cache = "session"} will use a temporary directory that persists only for the current session. \code{cache = "persistent"} will use \code{rappdirs::user_cache_dir()} to determine the appropriate directory to use. Otherwise, if a string is provided it will be assumed to be the path to the directory to use. In this case, an attempt will be made to create the cache directory if it does not exist. A warning will be given if a cached copy of the data exists and is more than 30 days old}

\item{refresh_cache}{logical: if TRUE, and a data file already exists in the cache, it will be refreshed. If FALSE, the cached copy will be used}

\item{simplified}{logical: if TRUE, only return a simplified set of columns (see details in "Value", below)}

\item{verbose}{logical: show progress messages?}
}
\value{
a data.frame or SpatialPointsDataFrame, with the following columns (note that not all information is populated for all place names):
\itemize{
  \item gaz_id - the unique identifier of each gazetteer entry. Note that the same feature (e.g. "Browns Glacier") might have multiple gazetteer entries, each with their own \code{gaz_id}, because the feature has been named multiple times by different naming authorities. The \code{scar_common_id} for these entries will be identical, because \code{scar_common_id} identifies the feature itself
  \item scar_common_id - the unique identifier (in the Composite Gazetteer of Antarctica) of the feature. A single feature may have multiple names, given by different naming authorities
  \item place_name - the name of the feature
  \item place_name_transliterated - the name of the feature transliterated to simple ASCII characters (e.g. with diacritical marks removed)
  \item longitude and latitude - the longitude and latitude of the feature (negative values indicate degrees west or south). Note that many features are not point features (e.g. mountains, lakes), in which case the \code{longitude} and \code{latitude} values are indicative only, generally of the centroid of the feature
  \item altitude - the altitude of the feature, in metres relative to sea level. Negative values indicate features below sea level
  \item feature_type_name - the feature type (e.g. "Archipelago", "Channel", "Mountain")
  \item date_named - the date on which the feature was named
  \item narrative - a text description of the feature; may include a synopsis of the history of its name
  \item named_for - the person after whom the feature was named, or other reason for its naming. For historical reasons the distinction between "narrative" and "named for" is not always obvious
  \item origin - the naming authority that provided the name. This is a country name, or organisation name for names that did not come from a national source
  \item relic - if \code{TRUE}, this name is associated with a feature that no longer exists (e.g. an ice shelf feature that has disappeared)
  \item gazetteer - the gazetteer from which this information came (currently only "CGA")
}
If \code{simplified} is FALSE, these additional columns will also be included:
\itemize{
  \item meeting_date - the date on which the name was formally approved by the associated naming authority. This is not available for many names: see the \code{date_named} column
  \item meeting_paper - references to papers or documents associated with the naming of the feature
  \item remote_sensor_info - text describing the remote sensing information (e.g. satellite platform name and image details) used to define the feature, if applicable
  \item coordinate_accuracy - an indicator of the accuracy of the coordinates, in metres
  \item altitude_accuracy - an indicator of the accuracy of the altitude value, in metres
  \item cga_source_gazetteer - for the Composite Gazetteer, this entry gives the source gazetteer from which this entry was taken. This is currently either a three-letter country code (e.g. "ESP", "USA") or "GEBCO" (for the GEBCO gazetteer of undersea features)
  \item country_name - the full name of the country where \code{cga_source_gazetteer} is a country
  \item source_name - the cartographic/GIS/remote sensing source from which the coordinates were derived
  \item source_publisher - where coordinates were derived from a map, the publisher of that map
  \item source_scale - the scale of the map from which the coordinates were derived
  \item source_institution - the institution from which the coordinate information came
  \item source_person - the contact person at the source institution, if applicable
  \item source_country_code - the country from which the coordinate information came
  \item source_identifier - where a coordinate or elevation was derived from a map, the identifier of that map
  \item comments - comments about the name or naming process
}
}
\description{
Place name data will be downloaded and optionally cached locally. If you wish to be able to use \code{antanym} offline, consider using \code{cache = "persistent"} so that the cached data will persist from one R session to the next. See \code{\link{an_cache_directory}} to get the path to the cache directory.
}
\examples{
\dontrun{
 ## download without caching
 g <- an_read()

 ## download to session cache, in sp format
 g <- an_read(cache = "session", sp = TRUE)

 ## download and cache to a persistent directory for later, offline use
 g <- an_read(cache = "persistent")

 ## refresh the cached copy
 g <- an_read(cache = "persistent", refresh_cache = TRUE)

 ## download and cache to a persistent directory of our choice
 g <- an_read(cache = "c:/my/cache/directory")
}

}
\references{
\url{https://data.aad.gov.au/aadc/gaz/scar/}, \url{https://www.scar.org/data-products/place-names/}
}
\seealso{
\code{\link{an_cache_directory}}, \code{\link{an_gazetteers}}, \code{\link{an_cga_metadata}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{an_filter}
\alias{an_filter}
\title{Filter a collection of place names by various criteria}
\usage{
an_filter(
  gaz,
  query,
  feature_ids,
  extent,
  feature_type,
  origin,
  origin_gazetteer,
  ignore_case = TRUE,
  as_regex = TRUE
)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: as returned by \code{\link{an_read}} or \code{\link{an_preferred}}}

\item{query}{character: vector of place name terms to search for. Returned place names will be those that match all entries in \code{query}}

\item{feature_ids}{numeric: return only place names associated with the features identified by these identifiers. Currently these values can only be \code{scar_common_id} values}

\item{extent}{vector of c(longitude_min, longitude_max, latitude_min, latitude_max): if provided, search only for names within this bounding box. \code{extent} can also be passed as a raster Extent object, a Raster object (in which case its extent will be used), a Spatial object (in which case the bounding box of the object will be used as the extent), or a matrix (in which case it will be assumed to be the output of \code{sp::bbox})}

\item{feature_type}{string: return only place names corresponding to feature types matching this pattern. For valid feature type names see \code{\link{an_feature_types}}}

\item{origin}{string: return only place names originating from bodies (countries or organisations) matching this pattern. For valid \code{origin} values see \code{link{an_origins}}}

\item{origin_gazetteer}{string: return only place names originating from gazetteers matching this pattern. For valid gazetteer names see \code{\link{an_gazetteers}}}

\item{ignore_case}{logical: if \code{TRUE}, use case-insensitive text matching}

\item{as_regex}{logical: if \code{TRUE}, treat \code{query} and other string input parameters as regular expressions. If \code{FALSE}, they will be treated as fixed strings to match against}
}
\value{
data.frame of results
}
\description{
A data frame of place names can be filtered according to name, geographic location, feature type, or other criteria. All text-related matches are by default treated as regular expressions and are case-insensitive: you can change this behaviour via the \code{ignore_case} and \code{as_regex} parameters.
}
\examples{
\dontrun{
 g <- an_read(cache = "session")

 ## simple search for any place name containing the word 'William'
 an_filter(g, query = "William")

 ## which bodies (countries or organisations) provided the names in our data?
 an_origins(g)

 ## find names containing "William" and originating from Australia or the USA
 an_filter(g, query = "William", origin = "Australia|United States of America")

 ## this search will return no matches
 ## because the actual place name is 'William Scoresby Archipelago'
 an_filter(g, query = "William Archipelago")

 ## we can split the search terms so that each is matched separately
 an_filter(g, query = c("William", "Archipelago"))

 ## or use a regular expression
 an_filter(g, query = "William .* Archipelago")

 ## or refine the search using feature type
 an_filter(g, query = "William", feature_type = "Archipelago")

 ## what feature types do we have in our data?
 an_feature_types(g)

 ## for more complex text searching, use regular expressions
 ## e.g. names matching "West" or "East"
 an_filter(g, query = "West|East")

 ## names starting with "West" or "East"
 an_filter(g, query = "^(West|East)")

 ## names with "West" or "East" appearing as complete words in the name
 ## ["\\b" matches a word boundary: see help("regex") ]
 an_filter(g, query = "\\\\b(West|East)\\\\b")

 ## filtering by spatial extent
 nms <- an_filter(g, extent = c(100, 120, -70, -65), origin = "Australia")
 with(nms, plot(longitude, latitude))
 with(nms, text(longitude, latitude, place_name))

 ## searching within the extent of an sp object
 my_sp <- sp::SpatialPoints(cbind(c(100, 120), c(-70, -65)))
 an_filter(g, extent = my_sp)

 ## or equivalently
 an_filter(g, extent = bbox(my_sp))

 ## or using the sp form of the gazetteer data
 gsp <- an_read(cache = "session", sp = TRUE)
 an_filter(gsp, extent = my_sp)

 ## using the pipe operator
 g \%>\% an_filter(query = "Ross", feature_type = "Ice shelf|Mountain")

 g \%>\% an_near(loc = c(100, -66), max_distance = 20) \%>\%
       an_filter(feature_type = "Island")

 ## find all names for feature 1589 and the naming
 ##  authority for each name
 an_filter(g, feature_ids = 1589)[, c("place_name", "origin")]
}
}
\references{
\url{https://data.aad.gov.au/aadc/gaz/scar/}, \url{https://www.scar.org/data-products/place-names/}
}
\seealso{
\code{\link{an_read}}, \code{\link{an_gazetteers}}, \code{\link{an_origins}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{an_feature_types}
\alias{an_feature_types}
\title{List feature types present in gazetteer data}
\usage{
an_feature_types(gaz)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: as returned by \code{\link{an_read}}, \code{\link{an_preferred}}, or \code{\link{an_filter}}}
}
\value{
character vector of country names
}
\description{
The gazetteer place names are associated with different feature types (e.g. "Hill", "Mountain", "Water body"). This function lists the feature types that are present in a given data frame.
}
\examples{
\dontrun{
 g <- an_read(cache = "session")

 ## what feature types do we have in our data?
 an_feature_types(g)
}
}
\seealso{
\code{\link{an_filter}} for filtering data according to feature type
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{an_get_url}
\alias{an_get_url}
\title{Get links to gazetteer entries}
\usage{
an_get_url(gaz)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: as returned by \code{\link{an_read}}, \code{\link{an_preferred}}, or \code{\link{an_filter}}}
}
\value{
character vector, where each component is a URL to a web page giving more information about the associated gazetteer entry
}
\description{
Each entry in the Composite Gazetteer of Antarctica has its own web page. The \code{an_url} function will return the URL of the page associated with a given gazetteer entry.
}
\examples{
\dontrun{
 g <- an_read(cache = "session")
 my_url <- an_get_url(an_filter(g, query = "Ufs Island")[1, ])
 browseURL(my_url)
}
}
\references{
\url{https://data.aad.gov.au/aadc/gaz/scar/}, \url{https://www.scar.org/data-products/place-names/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/re_exports.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{an_cache_directory}
\alias{an_cache_directory}
\title{The cache directory used by antanym}
\usage{
an_cache_directory(cache)
}
\arguments{
\item{cache}{string: the gazetteer data can be cached locally, so that it can be used offline later. Valid values are \code{"session"}, \code{"persistent"}, or a directory name. Specifying \code{cache="session"} will use a temporary directory that persists only for the current session. \code{cache="persistent"} will use \code{rappdirs::user_cache_dir()} to determine the appropriate directory to use. Otherwise, the input string will be assumed to be the path to the directory to use}
}
\value{
directory path
}
\description{
The cache directory used by antanym
}
\examples{
## per-session caching
an_cache_directory(cache = "session")

## persistent caching that will keep the data from one R session to the next
an_cache_directory(cache = "persistent")

}
\seealso{
\code{\link{an_read}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cga_metadata.R
\name{an_cga_metadata}
\alias{an_cga_metadata}
\title{Information about the Composite Gazetteer of Antarctica data structure}
\usage{
an_cga_metadata(simplified = TRUE)
}
\arguments{
\item{simplified}{logical: if TRUE, only describe the simplified set of columns (see the equivalent parameter in \code{\link{an_read}})}
}
\value{
a data frame with columns "field" and "description"
}
\description{
The Composite Gazetteer of Antarctica data structure (as returned by \code{\link{an_read}}):
}
\examples{

an_cga_metadata()

}
\references{
\url{https://data.aad.gov.au/aadc/gaz/scar/}, \url{https://www.scar.org/data-products/place-names/}
}
\seealso{
\code{\link{an_read}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thin.R
\name{an_thin}
\alias{an_thin}
\title{Thin names to give approximately uniform spatial coverage}
\usage{
an_thin(gaz, n, score_col = "score", score_weighting = 5, row_limit = 2000)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: typically as returned by \code{\link{an_suggest}}}

\item{n}{numeric: number of names to return}

\item{score_col}{string: the name of the column that gives the relative score of each name (e.g. as returned by \code{an_suggest}). Names with higher scores will be
preferred by the thinning process. If the specified \code{score_col} column is not present in \code{gaz}, or if all values within that column are equal, then the
thinning will be based entirely on the spatial distribution of the features}

\item{score_weighting}{numeric: weighting of scores relative to spatial distribution. A lower \code{score_weighting} value will tend to choose lower-scored names
in order to achieve better spatial uniformity. A higher \code{score_weighting} value will trade spatial uniformity in favour of selecting
higher-scored names}

\item{row_limit}{integer: the maximum number of rows allowed in \code{gaz}; see Details. Data frames larger than this will not be processed (with an error).}
}
\value{
data.frame
}
\description{
The provided data.frame of names will be thinned down to a smaller number of names. The thinning process attempts to select a subset of names that are uniformly spatially distributed, while simultaneously choosing the most important names (according to their relative score in the \code{score_col} column.
}
\details{
Note that the algorithm calculates all pairwise distances between the rows of \code{gaz}. This is memory-intensive, and so if \code{gaz} has many rows the algorithm will fail or on some platforms might crash. Input \code{gaz} data.frames with more than \code{row_limit} rows will not be processed for this reason. You can try increasing \code{row_limit} from its default value if necessary.
}
\examples{
\dontrun{
 g <- an_read(cache = "session")

 ## get a single name per feature, preferring the
 ##  Japanese name where there is one
 g <- an_preferred(g, origin = "Japan")

 ## suggested names for a 100x100 mm map covering 60-90E, 70-60S
 ##  (this is about a 1:12M scale map)
 suggested <- an_suggest(g, map_extent = c(60, 90, -70, -60), map_dimensions = c(100, 100))

 ## find the top 20 names by score
 head(suggested, 20)

 ## find the top 20 names chosen for spatial coverage and score
 an_thin(suggested, 20)
}

}
\seealso{
\code{\link{an_read}}, \code{\link{an_suggest}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{an_near}
\alias{an_near}
\title{Find placenames near a given location}
\usage{
an_near(gaz, loc, max_distance)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: as returned by \code{\link{an_read}}, \code{\link{an_preferred}}, or \code{\link{an_filter}}}

\item{loc}{numeric: target location (a two-element numeric vector giving longitude and latitude, or a SpatialPoints object)}

\item{max_distance}{numeric: maximum search distance in kilometres}
}
\value{
data.frame of results
}
\description{
Find placenames near a given location
}
\examples{
\dontrun{
 g <- an_read(cache = "session")

 ## named features within 10km of 110E, 66S
 an_near(g, loc = c(110, -66), max_distance = 10)

 ## using pipe operator
 g \%>\% an_near(loc = c(100, -66), max_distance = 10)

 ## with sp objects
 gsp <- an_read(cache = "session", sp = TRUE)
 loc <- sp::SpatialPoints(matrix(c(110, -66), nrow = 1),
   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
 an_near(gsp, loc = loc, max_distance = 10)
}

}
\references{
\url{https://data.aad.gov.au/aadc/gaz/scar/}, \url{https://www.scar.org/data-products/place-names/}
}
\seealso{
\code{\link{an_read}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antanym.R
\docType{package}
\name{antanym}
\alias{antanym}
\alias{antanym-package}
\title{\pkg{antanym}}
\description{
Antarctic geographic place names from the Composite Gazetteer of Antarctica, and functions for working with those place names.
}
\references{
\url{http://data.aad.gov.au/aadc/gaz/scar}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{an_origins}
\alias{an_origins}
\title{List the origins of place names present in gazetteer data}
\usage{
an_origins(gaz)
}
\arguments{
\item{gaz}{data.frame or SpatialPointsDataFrame: as returned by \code{\link{an_read}}, \code{\link{an_preferred}}, or \code{\link{an_filter}}}
}
\value{
character vector of origin names (countries or organisations)
}
\description{
The Composite Gazetteer of Antarctica is a compilation of place names provided by different countries and organisations. This function lists the originating bodies that provided the names in a given data frame.
}
\examples{
\dontrun{
 g <- an_read(cache = "session")

 ## which bodies (countries or organisations) provided the names in our data?
 an_origins(g)
}
}
\seealso{
\code{\link{an_filter}} for filtering data according to origin
}
