mapr
====



[![R-check](https://github.com/ropensci/mapr/workflows/R-check/badge.svg)](https://github.com/ropensci/mapr/actions/)
[![cran checks](https://cranchecks.info/badges/worst/mapr)](https://cranchecks.info/pkgs/mapr)
[![codecov](https://codecov.io/gh/ropensci/mapr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/mapr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/mapr?color=FAB657)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/mapr)](https://cran.r-project.org/package=mapr)


Helper for making maps of species occurrence data, including for:

* spocc (https://github.com/ropensci/spocc)
* rgbif (https://github.com/ropensci/rgbif)
* and some `sp` classes

This package has utilities for making maps with:

* base R
* ggplot2
* Leaflet - via `leaflet` pkg
* GitHub Gists - via `gistr` package

Get started with the docs: https://docs.ropensci.org/mapr/

## Installation

Install `mapr`


```r
install.packages("mapr")
```

Or the development version from GitHub


```r
remotes::install_github("ropensci/mapr")
```


```r
library("mapr")
library("spocc")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/mapr/issues).
* License: MIT
* Get citation information for `mapr` in R doing `citation(package = 'mapr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
mapr 0.5.2
==========

### MINOR IMPROVEMENTS

* fix vignette issue (#40)


mapr 0.5.0
==========

### DEFUNCT

* `map_ggmap()` is defunct. authentication has changed, more trouble than its worth (#34)

### NEW FEATURES

* all mapping functions gain the `name` parameter to specify the column that holds the taxon name - if not given, we look for a column "name" (#32)

### MINOR IMPROVEMENTS

* fix vignette missing title (#37)

### BUG FIXES

* fix non-ascii strings in the two package datasets - and script added to make those datasets reproducible, including fixing non-ascii strings (#39)
* remove linked references to pkgs in docs that are not imported/suggested (#38)
* `map_plot()` speed up: using `maps::map()` instead of `rworldmap::getMap()`, faster (#35)
* improve internal handling of `name` parameter users can pass down through mapping functions (#36)
* `rgbif` added to Suggests - was used in examples but wasn't in Suggests - used conditionally in examples now


mapr 0.4.0
==========

### MINOR IMPROVEMENTS

* All `map_*()` functions now support `gbif_data` class from the `rgbif` package, which is the output of `rgbif::occ_data()` (#29)
* Added `.github` files for contributing, issue and PR templates (#30)


mapr 0.3.4
==========

### MINOR IMPROVEMENTS

* Now using markdown for docs (#27)
* Replaced `httr` with `crul` as http client (#26)

### Problem with ggmap

* Note that there is a problem with `map_ggmap` due to a bug in 
`ggmap`. It is fixed in the `ggmap` dev version, so should be fixed
in the CRAN version soon, hopefully.


mapr 0.3.0
==========

### NEW FEATURES

* Now in all functions, when there's more than 1 taxon, we'll do a separate
color for each taxon and draw a legend if applicable (#21) (#22)
* Added support for adding convex hulls to some of the plot types (#23)
thanks to @rossmounce for the feature request
* `map_leaflet()` now adds metadata as a popup to each marker (#18) (#25)


mapr 0.2.0
==========

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.3 RC
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes
     
## Reverse dependencies

There are no reverse dependencies.

---

This version includes a fix to the vignette that was causing a warning on two CRAN platforms.

Thanks! 
Scott Chamberlain
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/mapr/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/mapr.git`
* Make sure to track progress upstream (i.e., on our version of `mapr` at `ropensci/mapr`) by doing `git remote add upstream https://github.com/ropensci/mapr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/mapr`

### Vignette changes

If you want to contribute to the vignette, or add a new one, those are kept in the `inst/vign/` directory. The main vignette is in `inst/vign/mapr_vignette.Rmd`.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```

</details>
mapr
====

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![R-check](https://github.com/ropensci/mapr/workflows/R-check/badge.svg)](https://github.com/ropensci/mapr/actions/)
[![cran checks](https://cranchecks.info/badges/worst/mapr)](https://cranchecks.info/pkgs/mapr)
[![codecov](https://codecov.io/gh/ropensci/mapr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/mapr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/mapr?color=FAB657)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/mapr)](https://cran.r-project.org/package=mapr)


Helper for making maps of species occurrence data, including for:

* spocc (https://github.com/ropensci/spocc)
* rgbif (https://github.com/ropensci/rgbif)
* and some `sp` classes

This package has utilities for making maps with:

* base R
* ggplot2
* Leaflet - via `leaflet` pkg
* GitHub Gists - via `gistr` package

Get started with the docs: https://docs.ropensci.org/mapr/

## Installation

Install `mapr`

```{r eval=FALSE}
install.packages("mapr")
```

Or the development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/mapr")
```

```{r}
library("mapr")
library("spocc")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/mapr/issues).
* License: MIT
* Get citation information for `mapr` in R doing `citation(package = 'mapr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{An R Markdown Vignette made with knitr}
-->

rgbif vignette - Seach and retrieve data from the Global Biodiverity Information Facilty (GBIF)
======

### About the package

`rgbif` is an R package to search and retrieve data from the Global Biodiverity Information Facilty (GBIF). `rgbif` wraps R code around the [GBIF API][gbifapi] to allow you to talk to the BISON database from R. 

********************

#### Install rgbif

```{r install, comment=NA, warning=FALSE}
# install.packages("devtools"); library(devtools); install_github("rbison", "ropensci")
library(rgbif); library(XML); library(RCurl); library(plyr); library(ggplot2); library(maps)
```

********************

#### Get a list of the data networks in GBIF - and you can use the networkkey number to seach for occurrences for the specific provider in other functions

```{r networks, comment=NA, warning=FALSE}
# Test the function for a few networks
networks(maxresults=5)

# By name
networks('ORNIS')
```

********************

#### Get a list of the data providers in GBIF - and you can use the dataproviderkey number to seach for occurrences for the specific provider in other functions

```{r providers, comment=NA, warning=FALSE}
# Test the function for a few providers
providers(maxresults=5)

# By data provider name
providers('University of Texas-Austin')
```

********************

#### Get a list of the data resources in GBIF - and you can use the resourcekey number to seach for occurrences for the specific resource in other functions

```{r resources, comment=NA, warning=FALSE}
# Test the function for a few resources
resources(maxresults=5)

# By name
head(resources('Flora'))
```

********************

#### Get number of occurrences for a set of search parameters
```{r occurrencecount, comment=NA, warning=FALSE}
occurrencecount(scientificname = 'Accipiter erythronemius', coordinatestatus = TRUE)
occurrencecount(scientificname = 'Helianthus annuus', coordinatestatus = TRUE, year=2005, maxlatitude=20)
```

********************

#### Get possible values to be used in taxonomic rank arguments in functions

```{r taxrank, comment=NA, warning=FALSE}
taxrank()
```

********************

#### Seach by taxon to retrieve number of records per taxon found in GBIF
```{r taxoncount, comment=NA, warning=FALSE}
taxoncount(scientificname = 'Puma concolor')
taxoncount(scientificname = 'Helianthus annuus')
taxoncount(rank = 'family')
```

********************

#### Get taxonomic information on a specific taxon or taxa in GBIF by their taxon concept keys

```{r taxonget, comment=NA, warning=FALSE}
(out <- taxonsearch(scientificname = 'Puma concolor'))
```

********************

#### Search for taxa in GBIF

```{r taxonsearch, comment=NA, warning=FALSE}
taxonsearch(scientificname = 'Puma concolor', rank="species", maxresults=10)
taxonsearch(scientificname = 'Puma concolor', rank="species", dataproviderkey=1)
```

********************

#### Get data for a single occurrence. Note that data is returned as a list, so you have to convert to a data.frame, etc. as you wish
```{r occurrenceget, comment=NA, warning=FALSE}
occurrenceget(key = 13749100)$dataProvider$dataResources$dataResource$occurrenceRecords$TaxonOccurrence[1:10]
```

********************

```{r occurrencelist, comment=NA, warning=FALSE}
out <- occurrencelist(scientificname = 'Puma concolor', coordinatestatus=TRUE, maxresults=20)
```

Note that the default object printed from a call to `occurrencelist` is a list that contains:

+ NumberFound: number of occurrences found in search results.
+ TaxonNames: Unique list of taxonomic names in search results.
+ Coordinates: Min and max latitude and longitude of all occurrences.
+ Countries: Countries contained in results set.

```{r defaultoccurencelist, comment=NA, warning=FALSE}
out
```

Where do you get data after a call to the `occurrencelist` function? This is where `gbifdata` comes in. By default a call to `gbifdata` prints a minimal data.frame with just rows *name*, *latitude*, and *longitude*.

```{r minimaltrue, comment=NA, warning=FALSE}
gbifdata(out)
```

Though you can get more detailed data by calling *minimal=FALSE*.

```{r minimalfalse, comment=NA, warning=FALSE}
head( gbifdata(out, minimal=FALSE)[,1:6] )
```

And you can get all possible data by specifying *format=darwin*.

```{r occurrencelistdarwin, comment=NA, warning=FALSE}
out <- occurrencelist(scientificname = 'Puma concolor', coordinatestatus=TRUE, format="darwin", maxresults=20)
head( gbifdata(out, minimal=FALSE)[,1:6] )
```

********************

#### Maps

```{r gbifmap1, comment=NA, warning=FALSE}
splist <- c('Accipiter erythronemius', 'Junco hyemalis', 'Aix sponsa', 'Ceyx fallax', 'Picoides lignarius', 'Campephilus leucopogon')
out <- occurrencelist_many(splist, coordinatestatus = TRUE, maxresults = 20)
gbifmap_list(out)
```

Another example, setting scientificname="*" so we just get any species, and then mapping points only within the state of Texas in the US.

```{r gbifmap2, comment=NA, warning=FALSE}
out <- occurrencelist(scientificname="*", minlatitude=30, maxlatitude=35, minlongitude=-100, maxlongitude=-95, coordinatestatus = TRUE, maxresults = 200)
gbifmap_list(input=out, mapdatabase="state", region="texas", geom=geom_jitter, jitter=position_jitter(width = 0.3, height = 0.3))
```

********************

[gbifapi]: http://data.gbif.org/tutorial/services---
title: "mapr introduction"
author: "Scott Chamberlain"
date: "2020-07-29"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{mapr introduction}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



## Load spocc and mapr


```r
library("spocc")
library("mapr")
```

## Interactive maps

### Leaflet.js

[Leaflet JS](http://leafletjs.com/) is an open source mapping library that can leverage various layers from multiple sources. Using the [`leaflet`](https://cran.r-project.org/package=leaflet) library, we can generate a local interactive map of species occurrence data.

An example:


```r
spp <- c('Danaus plexippus','Accipiter striatus','Pinus contorta')
dat <- occ(query = spp, from = 'gbif', has_coords = TRUE, limit = 100)
map_leaflet(dat)
```

![leaflet](../man/figures/leaflet.png)

### Geojson map as a Github gist

You can also create interactive maps via the `mapgist` function. You have to have a Github account to use this function. Github accounts are free though, and great for versioning and collaborating on code or papers. When you run the `map_gist` function it will ask for your Github username and password. You can alternatively store those in your `.Rprofile` file by adding entries for username (`options(github.username = 'username')`) and password (`options(github.password = 'password')`).


```r
spp <- c('Danaus plexippus', 'Accipiter striatus', 'Pinus contorta')
dat <- occ(query = spp, from = 'gbif', has_coords = TRUE, limit = 100)
dat <- fixnames(dat)
map_gist(dat, color = c("#976AAE", "#6B944D", "#BD5945"))
```

![gist](../man/figures/gist.png)

## Static maps

### base plots

Base plots, or the built in plotting facility in R accessed via `plot()`, is quite fast, but not easy or efficient to use, but are good for a quick glance at some data.


```r
spnames <- c('Accipiter striatus', 'Setophaga caerulescens', 'Spinus tristis')
out <- occ(query = spnames, from = 'gbif', has_coords = TRUE, limit = 100)
map_plot(out, size = 1, pch = 10)
```

![plot of chunk unnamed-chunk-5](../man/figures/unnamed-chunk-5-1.png)

### ggplot2

`ggplot2` is a powerful package for making visualizations in R. Read more about it [here](https://cran.r-project.org/package=ggplot2).


```r
dat <- occ(query = 'Lynx rufus californicus', from = 'gbif', has_coords = TRUE, limit = 200)
map_ggplot(dat, map = "usa")
```

![plot of chunk unnamed-chunk-6](../man/figures/unnamed-chunk-6-1.png)

### via dismo

if that's your jam, though you might find `rgbif` easier


```r
library("dismo")
g <- gbif('Batrachoseps', 'luciae', geo = TRUE, end = 300)
map_leaflet(g, lon = "lon", lat = "lat", name = "scientificName")
```

![dismomap](http://f.cl.ly/items/2u2V0n0B3Y2y0p1d0f1M/Screen%20Shot%202016-01-22%20at%204.46.12%20PM.png)

## The supported inputs

All functions take the following kinds of inputs:

* An object of class `occdat`, from the package `spocc`. An object of
this class is composed of many objects of class `occdatind`
* An object of class `occdatind`, from the package `spocc`
* An object of class `gbif`, from the package `rgbif`
* An object of class `data.frame`. This data.frame can have any columns, but
must include a column for taxonomic names (e.g., `name`), and for latitude
and longitude (we guess your lat/long columns, starting with the default
`latitude` and `longitude`).
* An object of class `SpatialPoints`
* An object of class `SpatialPointsDatFrame`
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapr-package.R
\docType{data}
\name{occdat_eg1}
\alias{occdat_eg1}
\title{Example dataset: output from call to \code{\link[spocc:occ]{spocc::occ()}}}
\format{
A data frame with 25 rows and 62 variables
}
\description{
A dataset with 25 rows, and 62 columns, from the query:
\code{occ(query='Accipiter striatus', from='gbif', limit=25, has_coords=TRUE)}
}
\details{
See \code{inst/ignore/datasets.R} for the code to prepare this dataaset
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_ggmap.R
\name{map_ggmap}
\alias{map_ggmap}
\title{ggpmap visualization of species occurences}
\usage{
map_ggmap(...)
}
\arguments{
\item{...}{Ignored}
}
\description{
THIS FUNCTION IS DEFUNCT
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ2sp.R
\name{occ2sp}
\alias{occ2sp}
\title{Create a spatial points dataframe from a spocc search}
\usage{
occ2sp(x, coord_string = "+proj=longlat +datum=WGS84", just_coords = FALSE)
}
\arguments{
\item{x}{The resuslts of a spocc search called by \code{\link[spocc:occ]{spocc::occ()}}}

\item{coord_string}{A valid EPGS cooridate string from the sp package,
the default is WSGS 84}

\item{just_coords}{Return data frame with specios names and provenance or
just a spatial points object, which is the default.}
}
\description{
Create a spatial points dataframe from a spocc search
}
\details{
This function will return either a spatial points dataframe or
spatial points object. Conversion to spatial points objects allows spocc
searches to interact with other spatial data sources. More coordinate system
codes can be found at the EPGS registry
}
\examples{
\dontrun{
### See points on a map
library("maptools")
library("spocc")
data(wrld_simpl)
plot(wrld_simpl[wrld_simpl$NAME == "United States", ], xlim = c(-70, -60))
out <- occ(query = "Accipiter striatus", from = c("vertnet", "gbif"),
  limit = 50)
xx <- occ2sp(out, just_coords = TRUE)
points(xx, col = 2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_github.r
\name{style_geojson}
\alias{style_geojson}
\title{Style a data.frame prior to converting to geojson.}
\usage{
style_geojson(
  input,
  var = NULL,
  var_col = NULL,
  var_sym = NULL,
  var_size = NULL,
  color = NULL,
  symbol = NULL,
  size = NULL
)
}
\arguments{
\item{input}{A data.frame}

\item{var}{A single variable to map colors, symbols, and/or sizes to.}

\item{var_col}{The variable to map colors to.}

\item{var_sym}{The variable to map symbols to.}

\item{var_size}{The variable to map size to.}

\item{color}{Valid RGB hex color}

\item{symbol}{An icon ID from the Maki project
https://labs.mapbox.com/maki-icons/ or a single alphanumeric character
(a-z or 0-9).}

\item{size}{One of 'small', 'medium', or 'large'}
}
\description{
Style a data.frame prior to converting to geojson.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_ggplot.R
\name{map_ggplot}
\alias{map_ggplot}
\title{ggplot2 mapping}
\usage{
map_ggplot(
  x,
  map = "world",
  point_color = "#86161f",
  color = NULL,
  size = 3,
  lon = "longitude",
  lat = "latitude",
  name = NULL,
  ...
)
}
\arguments{
\item{x}{The data. An object of class \code{occdat}, \code{occdatind},
\code{gbif}, \code{gbif_data}, \code{SpatialPoints},
\code{SpatialPointsDataFrame}, or \code{data.frame}. The package
\pkg{spocc} needed for
the first two, and \pkg{rgbif} needed for the third. When \code{data.frame}
input, any number of columns allowed, but with at least the following:
name (the taxonomic name), latitude (in dec. deg.), longitude (in dec. deg.)}

\item{map}{(character) One of world, world2, state, usa, county, france,
italy, or nz}

\item{point_color}{Default color of your points. Deprecated, use
\code{color}}

\item{color}{Default color of your points.}

\item{size}{point size, Default: 3}

\item{lon, lat}{(character) Longitude and latitude variable names. Ignored
unless \code{data.frame} input to \code{x} parameter. We attempt to guess,
but if nothing close, we stop. Default: \code{longitude} and
\code{latitude}}

\item{name}{(character) the column name that contains the name to use in
creating the map. If left \code{NULL} we look for a "name" column.
default: \code{NULL}}

\item{...}{Ignored}
}
\value{
A ggplot2 map, of class \code{gg/ggplot}
}
\description{
ggplot2 mapping
}
\examples{
# map spocc output, here using a built in object
data(occdat_eg1)
map_ggplot(occdat_eg1)

# map rgbif output, here using a built in object
data(gbif_eg1)
map_ggplot(gbif_eg1)

\dontrun{
## spocc
library("spocc")
ddat <- occ(query = 'Lynx rufus californicus', from = 'gbif', limit=100)
map_ggplot(ddat)
map_ggplot(ddat$gbif)
map_ggplot(ddat$gbif, "usa")
map_ggplot(ddat, "county")

### usage of occ2sp()
#### SpatialPoints
spdat <- occ2sp(ddat)
map_ggplot(spdat)
#### SpatialPointsDataFrame
spdatdf <- as(spdat, "SpatialPointsDataFrame")
map_ggplot(spdatdf)

## rgbif
if (requireNamespace("rgbif")) {
library("rgbif")
library("ggplot2")
### occ_search() output
res <- occ_search(scientificName = "Puma concolor", limit = 100)
map_ggplot(res)

### occ_data() output
res <- occ_data(scientificName = "Puma concolor", limit = 100)
map_ggplot(res)

#### many taxa
res <- occ_data(scientificName = c("Puma concolor", "Quercus lobata"), 
   limit = 30)
map_ggplot(res)

### add a convex hull
hull(map_ggplot(res))
}

## data.frame
df <- data.frame(name = c('Poa annua', 'Puma concolor', 'Foo bar'),
                 longitude = c(-120, -121, -121),
                 latitude = c(41, 42, 45), stringsAsFactors = FALSE)
map_ggplot(df)

# many species, each gets a different color
library("spocc")
spp <- c('Danaus plexippus', 'Accipiter striatus', 'Pinus contorta')
dat <- occ(spp, from = 'gbif', limit = 30, has_coords = TRUE)
map_ggplot(dat, color = c('#976AAE', '#6B944D', '#BD5945'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapr-package.R
\docType{package}
\name{mapr-package}
\alias{mapr-package}
\alias{mapr}
\title{mapr}
\description{
Visualize species occurrence data
}
\section{Many inputs}{

All functions take the following kinds of inputs:
\itemize{
\item An object of class \code{occdat}, from the package \pkg{spocc}.
An object of this class is composed of many objects of class
\code{occdatind}
\item An object of class \code{occdatind}, from the package \pkg{spocc}
\item An object of class \code{gbif}, from the package \pkg{rgbif}
\item An object of class \code{data.frame}. This data.frame can have any
columns, but must include a column for taxonomic names (e.g., \code{name}),
and for latitude and longitude (we guess your lat/long columns, starting
with the default \code{latitude} and \code{longitude})
\item An object of class \code{SpatialPoints}
\item An object of class \code{SpatialPointsDatFrame}
}
}

\section{Package API}{

\itemize{
\item \code{\link[=map_plot]{map_plot()}} - static Base R plots
\item \code{\link[=map_ggplot]{map_ggplot()}} - static ggplot2 plots
\item \code{\link[=map_leaflet]{map_leaflet()}} - interactive Leaflet.js interactive maps
\item \code{\link[=map_gist]{map_gist()}} - ineractive, shareable maps on GitHub Gists
}
}

\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapr-package.R
\docType{data}
\name{gbif_eg1}
\alias{gbif_eg1}
\title{Example dataset: output from call to \code{rgbif::occ_search()}}
\format{
A data frame with 50 rows and 101 variables
}
\description{
A dataset with 50 rows, and 101 columns, from the query:
\code{rgbif::occ_search(scientificName = "Puma concolor", limit = 100)}
}
\details{
See \code{inst/ignore/datasets.R} for the code to prepare this dataaset
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hull.R
\name{hull}
\alias{hull}
\title{Add a convex hull to a map}
\usage{
hull(x, ...)
}
\arguments{
\item{x}{input}

\item{...}{ignored}
}
\value{
Adds a convex hull to the plot. See \code{\link[grDevices:chull]{grDevices::chull()}}
for info.
}
\description{
Add a convex hull to a map
}
\details{
Can be used with \code{\link[=map_leaflet]{map_leaflet()}}, \code{\link[=map_plot]{map_plot()}},
and \code{\link[=map_ggplot]{map_ggplot()}}. Other methods in this package may be supported
in the future.
}
\examples{
# map spocc output, here using a built in object
data(occdat_eg1)
map_plot(occdat_eg1, hull = TRUE)

# map rgbif output, here using a built in object
hull(map_ggplot(occdat_eg1))

\dontrun{
# leaflet
library("spocc")
spp <- c('Danaus plexippus', 'Accipiter striatus', 'Pinus contorta')
dat <- occ(spp, from = 'gbif', limit = 30, has_coords = TRUE)
hull(map_leaflet(dat))

# ggplot
if (requireNamespace("rgbif")) {
library("rgbif")
res <- occ_search(scientificName = "Puma concolor", limit = 100)
hull(map_ggplot(res))
}

# base plots
library("spocc")
out <- occ(query='Accipiter striatus', from='gbif', limit=25,
  has_coords=TRUE)
map_plot(out, hull = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_gist.R
\name{map_gist}
\alias{map_gist}
\title{Make an interactive map to view in the browser as a GitHub gist}
\usage{
map_gist(
  x,
  description = "",
  public = TRUE,
  browse = TRUE,
  lon = "longitude",
  lat = "latitude",
  name = NULL,
  ...
)
}
\arguments{
\item{x}{The data. An object of class \code{occdat}, \code{occdatind},
\code{gbif}, \code{gbif_data}, \code{SpatialPoints},
\code{SpatialPointsDataFrame}, or \code{data.frame}. The package
\pkg{spocc} needed for
the first two, and \pkg{rgbif} needed for the third. When \code{data.frame}
input, any number of columns allowed, but with at least the following:
name (the taxonomic name), latitude (in dec. deg.), longitude (in dec. deg.)}

\item{description}{Description for the Github gist, or leave to
default (=no description)}

\item{public}{(logical) Whether gist is public (default: \code{TRUE})}

\item{browse}{If \code{TRUE} (default) the map opens in your default browser.}

\item{lon, lat}{(character) Longitude and latitude variable names. Ignored
unless \code{data.frame} input to \code{x} parameter. We attempt to guess,
but if nothing close, we stop. Default: \code{longitude} and
\code{latitude}}

\item{name}{(character) the column name that contains the name to use in
creating the map. If left \code{NULL} we look for a "name" column.}

\item{...}{Further arguments passed on to \code{\link[=style_geojson]{style_geojson()}}}
}
\description{
Make an interactive map to view in the browser as a GitHub gist
}
\details{
See \code{\link[gistr:gist_auth]{gistr::gist_auth()}} for help on authentication

Does not support adding a convex hull via \code{\link[=hull]{hull()}}
}
\examples{
\dontrun{
## spocc
library("spocc")
spp <- c('Danaus plexippus', 'Accipiter striatus', 'Pinus contorta')
dat <- spocc::occ(spp, from=c('gbif','ecoengine'), limit=30,
  gbifopts=list(hasCoordinate=TRUE))

# Define colors
map_gist(dat, color=c('#976AAE','#6B944D','#BD5945'))

# Define colors and marker size
map_gist(dat, color=c('#976AAE','#6B944D','#BD5945'),
  size=c('small','medium','large'))

# Define symbols
map_gist(dat, symbol=c('park','zoo','garden'))

## rgbif
if (requireNamespace("rgbif")) {
library("rgbif")
### occ_search() output
res <- occ_search(scientificName = "Puma concolor", limit = 100)
map_gist(res)

### occ_data() output
res <- occ_data(scientificName = "Puma concolor", limit = 100)
map_gist(res)

#### many taxa
res <- occ_data(scientificName = c("Puma concolor", "Quercus lobata"), 
   limit = 30)
res
map_gist(res)
}

## data.frame
df <- data.frame(name = c('Poa annua', 'Puma concolor', 'Foo bar'),
                 longitude = c(-120, -121, -121),
                 latitude = c(41, 42, 45), stringsAsFactors = FALSE)
map_gist(df)

### usage of occ2sp()
#### SpatialPoints
spdat <- occ2sp(dat)
map_gist(spdat)
#### SpatialPointsDataFrame
spdatdf <- as(spdat, "SpatialPointsDataFrame")
map_gist(spdatdf)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_plot.R
\name{map_plot}
\alias{map_plot}
\title{Base R mapping}
\usage{
map_plot(
  x,
  lon = "longitude",
  lat = "latitude",
  color = NULL,
  size = 1,
  pch = 16,
  hull = FALSE,
  name = NULL,
  ...
)
}
\arguments{
\item{x}{The data. An object of class \code{occdat}, \code{occdatind},
\code{gbif}, \code{gbif_data}, \code{SpatialPoints},
\code{SpatialPointsDataFrame}, or \code{data.frame}. The package
\pkg{spocc} needed for
the first two, and \pkg{rgbif} needed for the third. When \code{data.frame}
input, any number of columns allowed, but with at least the following:
name (the taxonomic name), latitude (in dec. deg.), longitude (in dec. deg.)}

\item{lon, lat}{(character) Longitude and latitude variable names. Ignored
unless \code{data.frame} input to \code{x} parameter. We attempt to guess,
but if nothing close, we stop. Default: \code{longitude} and
\code{latitude}}

\item{color}{Default color of your points.}

\item{size}{point size, passed to \code{cex} Default: 1}

\item{pch}{point symbol shape, Default: 16}

\item{hull}{(logical) whether to add a convex hull. Default: \code{FALSE}}

\item{name}{(character) the column name that contains the name to use in
creating the map. If left \code{NULL} we look for a "name" column.}

\item{...}{Further args to \code{\link[graphics:points]{graphics::points()}}}
}
\value{
Plots a world scale map
}
\description{
Base R mapping
}
\examples{
# map spocc output, here using a built in object
data(occdat_eg1)
map_plot(occdat_eg1)

# map rgbif output, here using a built in object
data(gbif_eg1)
map_plot(gbif_eg1)

\dontrun{
## spocc
library("spocc")
(out <- occ(query='Accipiter striatus', from='gbif', limit=25,
  has_coords=TRUE))
### class occdat
map_plot(out)
map_plot(out, hull = TRUE)
### class occdatind
map_plot(out$gbif)
map_plot(out$gbif, hull = TRUE)

## rgbif
if (requireNamespace("rgbif")) {
library("rgbif")
### occ_search() output
res <- occ_search(scientificName = "Puma concolor", limit = 100)
map_plot(res)
map_plot(res, hull = TRUE)

### occ_data() output
res <- occ_data(scientificName = "Puma concolor", limit = 100)
map_plot(res)
#### many taxa
res <- occ_data(scientificName = c("Puma concolor", "Quercus lobata"), 
   limit = 30)
res
map_plot(res)
}


## data.frame
df <- data.frame(
  name = c('Poa annua', 'Puma concolor', 'Foo bar', 'Stuff things'),
  longitude = c(-125, -123, -121, -110),
  latitude = c(41, 42, 45, 30), stringsAsFactors = FALSE)
map_plot(df)
map_plot(df, hull = TRUE)

### usage of occ2sp()
#### SpatialPoints
spdat <- occ2sp(out)
map_plot(spdat)
map_plot(spdat, hull = TRUE)

#### SpatialPointsDataFrame
spdatdf <- as(spdat, "SpatialPointsDataFrame")
map_plot(spdatdf)
map_plot(spdatdf, hull = TRUE)

# many species, each gets a different color
library("spocc")
spp <- c('Danaus plexippus', 'Accipiter striatus', 'Pinus contorta',
  'Ursus americanus')
dat <- occ(spp, from = 'gbif', limit = 30, has_coords = TRUE,
  gbifopts = list(country = 'US'))
map_plot(dat)
map_plot(dat, hull = TRUE)
## diff. color for each taxon
map_plot(dat, color = c('#976AAE', '#6B944D', '#BD5945', 'red'))
map_plot(dat, color = c('#976AAE', '#6B944D', '#BD5945', 'red'), hull = TRUE)

# add a convex hull
if (requireNamespace("rgbif")) {
library("rgbif")
res <- occ_search(scientificName = "Puma concolor", limit = 100)
map_plot(res, hull = FALSE)
map_plot(res, hull = TRUE)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_leaflet.R
\name{map_leaflet}
\alias{map_leaflet}
\title{Make interactive maps with Leaflet.js}
\usage{
map_leaflet(
  x,
  lon = "longitude",
  lat = "latitude",
  color = NULL,
  size = 13,
  name = NULL,
  ...
)
}
\arguments{
\item{x}{The data. An object of class \code{occdat}, \code{occdatind},
\code{gbif}, \code{gbif_data}, \code{SpatialPoints},
\code{SpatialPointsDataFrame}, or \code{data.frame}. The package
\pkg{spocc} needed for
the first two, and \pkg{rgbif} needed for the third. When \code{data.frame}
input, any number of columns allowed, but with at least the following:
name (the taxonomic name), latitude (in dec. deg.), longitude (in dec. deg.)}

\item{lon, lat}{(character) Longitude and latitude variable names. Ignored
unless \code{data.frame} input to \code{x} parameter. We attempt to guess,
but if nothing close, we stop. Default: \code{longitude} and
\code{latitude}}

\item{color}{Default color of your points.}

\item{size}{point size, Default: 13}

\item{name}{(character) the column name that contains the name to use in
creating the map. If left \code{NULL} we look for a "name" column.}

\item{...}{Ignored}
}
\value{
a Leaflet map in Viewer in Rstudio, or in your default browser
otherwise
}
\description{
Make interactive maps with Leaflet.js
}
\details{
We add popups by default, and add all columns to the popup. The
html is escaped with \code{htmltools::htmlEscape()}
}
\examples{
\dontrun{
## spocc
library("spocc")
(out <- occ(query='Accipiter striatus', from='gbif', limit=50,
  has_coords=TRUE))
### with class occdat
map_leaflet(out)
### with class occdatind
map_leaflet(out$gbif)
### use occ2sp
map_leaflet(occ2sp(out))

## rgbif
if (requireNamespace("rgbif")) {
library("rgbif")
res <- occ_search(scientificName = "Puma concolor", limit = 100)
map_leaflet(res)
}

## SpatialPoints class
library("sp")
df <- data.frame(longitude = c(-120,-121),
                 latitude = c(41, 42), stringsAsFactors = FALSE)
x <- SpatialPoints(df)
map_leaflet(x)

## SpatialPointsDataFrame class
if (requireNamespace("rgbif")) {
library("rgbif")
### occ_search() output
res <- occ_search(scientificName = "Puma concolor", limit = 100)
x <- res$data
library("sp")
x <- x[stats::complete.cases(x$decimalLatitude, x$decimalLongitude), ]
coordinates(x) <- ~decimalLongitude+decimalLatitude
map_leaflet(x)

### occ_data() output
res <- occ_data(scientificName = "Puma concolor", limit = 100)
map_leaflet(res)
}

#### many taxa
res <- occ_data(scientificName = c("Puma concolor", "Quercus lobata"), 
   limit = 30)
res
map_leaflet(res)


## data.frame
df <- data.frame(name = c('Poa annua', 'Puma concolor'),
                 longitude = c(-120,-121),
                 latitude = c(41, 42), stringsAsFactors = FALSE)
map_leaflet(df)

# many species
library("spocc")
spp <- c('Danaus plexippus', 'Accipiter striatus', 'Pinus contorta')
dat <- occ(spp, from = 'gbif', limit = 50, has_coords = TRUE)
map_leaflet(dat)
map_leaflet(dat, color = c('#AFFF71', '#AFFF71', '#AFFF71'))
map_leaflet(dat, color = c('#976AAE', '#6B944D', '#BD5945'))

# add a convex hull
## map_leaflet(dat) \%>\% hull()  # using pipes
hull(map_leaflet(dat))
}
}
