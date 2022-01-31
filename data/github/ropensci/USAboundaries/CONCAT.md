---
title: "USAboundaries: Historical and Contemporary Boundaries of the United States of America"
tags: 
  - R
  - spatial
  - history
  - digital-history
authors:
  - name: Lincoln A. Mullen
    orcid: 0000-0001-5103-6917
    affiliation: 1
  - name: Jordan Bratt
    orcid: 0000-0001-9051-7203
    affiliation: 1
affiliations:
  - name: Department of History and Art History, George Mason University
    index: 1
date: 20 March 2018
bibliography: paper.bib
---

The USAboundaries package for R provides contemporary and historical boundaries of the United States of America [@USAboundaries]. (The package is available on [GitHub](https://github.com/ropensci/USAboundaries/) and archived on [Zenodo](https://doi.org/10.5281/zenodo.825218).) Historical data in the package includes state and county boundaries from 1629 to 2000 from the Newberry Library's "Atlas of Historical County Boundaries" [@ahcb]. Also included is historical city population data from Erik Steiner's "United States Historical City Populations, 1790-2010" [@steiner-cities]. Contemporary data in the package includes state, county, and Congressional district boundaries, as well as Zip Code Tabulation Area centroids. These data are all drawn from the U.S. Census Bureau [@census].

These historical and contemporary boundaries are provided at different resolutions suitable for national and state-level mapping. A consistent interface provides a way to easily select historical boundaries for any specific date. The package includes helper functions and datasets, including tables of state names, abbreviations and FIPS codes for joining to attribute data, as well as functions and data to get State Plane Coordinate System projections as EPSG codes or PROJ.4 strings [@stateplane]. A first step in many spatial analyses is joining data of interest to spatial data, which the datasets in this package enable.

This package underlies the [*Mapping Early American Elections*](http://earlyamericanelections.org/) project created by a team at the Roy Rosenzweig Center for History and New Media [@meae]. That project maps Congressional elections and state legislative elections from 1787 to 1825. The USAboundaries package provides access to the frequently changing state and county boundaries during that time period.  

*Development of this package was funded in part by a Humanities Collections and Reference Resources grant from the Division of Preservation and Access at the National Endowment for the Humanities (grant number PW-234776-16). The package is part of [rOpenSci](https://ropensci.org/).*

# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

# USAboundaries

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/USAboundaries)](https://cran.r-project.org/package=USAboundaries)
[![JOSS
Status](https://joss.theoj.org/papers/3458a33133aa6c069ab4dd8df0b5f3b5/status.svg)](https://doi.org/10.21105/joss.00314)
[![R-CMD-check](https://github.com/ropensci/USAboundaries/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/USAboundaries/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/USAboundaries/master.svg)](https://codecov.io/github/ropensci/USAboundaries?branch=master)

## Overview

This R package includes contemporary state, county, and Congressional
district boundaries, as well as zip code tabulation area centroids. It
also includes historical boundaries from 1629 to 2000 for states and
counties from the Newberry Library’s [Atlas of Historical County
Boundaries](https://publications.newberry.org/ahcbp/), as well as
historical city population data from Erik Steiner’s “[United States
Historical City Populations,
1790-2010](https://github.com/cestastanford/historical-us-city-populations).”
The package has some helper data, including a table of state names,
abbreviations, and FIPS codes, and functions and data to get [State
Plane Coordinate
System](https://en.wikipedia.org/wiki/State_Plane_Coordinate_System)
projections as EPSG codes or PROJ.4 strings.

This package can serve a number of purposes. The spatial data can be
joined to any other kind of data in order to make thematic maps. Unlike
other R packages, this package also contains historical data for use in
analyses of the recent or more distant past. See the [“A sample analysis
using
USAboundaries”](http://lincolnmullen.com/software/usaboundaries/articles/usaboundaries-sample-analysis.html)
vignette for an example of how the package can be used for both
historical and contemporary maps.

## Citation

If you use this package in your research, we would appreciate a
citation.

``` r
citation("USAboundaries")
#> 
#> To cite the USAboundaries package in publications, please cite the
#> paper in the Journal of Open Source Software:
#> 
#>   Lincoln A. Mullen and Jordan Bratt, "USAboundaries: Historical and
#>   Contemporary Boundaries of the United States of America," Journal of
#>   Open Source Software 3, no. 23 (2018): 314,
#>   https://doi.org/10.21105/joss.00314.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {{USAboundaries}: Historical and Contemporary Boundaries
#> of the United States of America},
#>     author = {Lincoln A. Mullen and Jordan Bratt},
#>     journal = {Journal of Open Source Software},
#>     year = {2018},
#>     volume = {3},
#>     issue = {23},
#>     pages = {314},
#>     url = {https://doi.org/10.21105/joss.00314},
#>     doi = {10.21105/joss.00314},
#>   }
```

## Installation

You can install this package from CRAN.

    install.packages("USAboundaries")

Almost all of the data for this package is provided by the
[USAboundariesData
package](https://github.com/ropensci/USAboundariesData). That package
will be automatically installed (with your permission) from the
[rOpenSci package repository](https://ropensci.r-universe.dev) the first
time that you need it.

Or you can install the development versions from GitHub using
[remotes](https://remotes.r-lib.org).

    # install.packages("remotes")
    remotes::install_github("ropensci/USAboundaries")
    remotes::install_github("ropensci/USAboundariesData")

## Use

This package provides a set of functions, one for each of the types of
boundaries that are available. These functions have a consistent
interface.

Passing a date to `us_states()`, `us_counties()`, and `us_cities()`
returns the historical boundaries for that date. If no date argument is
passed, then contemporary boundaries are returned. The functions
`us_congressional()` and `us_zipcodes()` only offer contemporary
boundaries.

For almost all functions, pass a character vector of state names or
abbreviations to the `states =` argument to return only those states or
territories.

For certain functions, more or less detailed boundary information is
available by passing an argument to the `resolution =` argument.

See the examples below to see how the interface works, and see the
documentation for each function for more details.

``` r
library(USAboundaries) 
library(sf) # for plotting and projection methods
#> Linking to GEOS 3.9.1, GDAL 3.3.2, PROJ 8.1.1

states_1840 <- us_states("1840-03-12")
plot(st_geometry(states_1840))
title("U.S. state boundaries on March 3, 1840")
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

``` r
states_contemporary <- us_states()
plot(st_geometry(states_contemporary))
title("Contemporary U.S. state boundaries")
```

![](man/figures/README-unnamed-chunk-3-2.png)<!-- -->

``` r
counties_va_1787 <- us_counties("1787-09-17", states = "Virginia")
plot(st_geometry(counties_va_1787))
title("County boundaries in Virginia in 1787")
```

![](man/figures/README-unnamed-chunk-3-3.png)<!-- -->

``` r
counties_va <- us_counties(states = "Virginia")
plot(st_geometry(counties_va))
title("Contemporary county boundaries in Virginia")
```

![](man/figures/README-unnamed-chunk-3-4.png)<!-- -->

``` r
counties_va_highres <- us_counties(states = "Virginia", resolution = "high")
plot(st_geometry(counties_va_highres))
title("Higher resolution contemporary county boundaries in Virginia")
```

![](man/figures/README-unnamed-chunk-3-5.png)<!-- -->

``` r
congress <- us_congressional(states = "California")
plot(st_geometry(congress))
title("Congressional district boundaries in California")
```

![](man/figures/README-unnamed-chunk-3-6.png)<!-- -->

## State plane projections

The `state_plane()` function returns EPSG codes and PROJ.4 strings for
the State Plane Coordinate System. You can use these to use suitable
projections for specific states.

``` r
va <- us_states(states = "VA", resolution = "high")
plot(st_geometry(va), graticule = TRUE)
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

``` r
va_projection <- state_plane("VA")
va <- st_transform(va, va_projection)
plot(st_geometry(va), graticule = TRUE)
```

![](man/figures/README-unnamed-chunk-4-2.png)<!-- -->

## Related packages

Each function returns an `sf` object from the
[sf](https://cran.r-project.org/package=sf) package, which can be mapped
using the [leaflet](https://cran.r-project.org/package=leaflet) or
[ggplot2](https://cran.r-project.org/package=ggplot2) packages.

If you need U.S. Census Bureau boundary files which are not provided by
this package, consider using the
[tigris](https://cran.r-project.org/package=tigris) package, which
downloads those shapefiles.

## License

The historical boundary data provided in this package is available under
the CC BY-NC-SA 2.5 license from John H. Long, et al., [Atlas of
Historical County Boundaries](https://publications.newberry.org/ahcbp/),
Dr. William M. Scholl Center for American History and Culture, The
Newberry Library, Chicago (2010). Please cite that project if you use
this package in your research and abide by the terms of their license if
you use the historical information.

The historical population data for cities is provided by U.S. Census
Bureau and Erik Steiner, Spatial History Project, Center for Spatial and
Textual Analysis, Stanford University. See the data in [this
repository](https://github.com/cestastanford/historical-us-city-populations).

The contemporary data is provided by the U.S. Census Bureau and is in
the public domain.

All code in this package is copyright [Lincoln
Mullen](http://lincolnmullen.com) and is released under the MIT license.

------------------------------------------------------------------------

[![rOpenSci
footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
# USAboundaries 0.3.1

- New vignette demonstrating the package's functionality (#40).
- Additions and clarifications to documentation following @AndySouth's suggestions for JOSS peer review (#38).
- `us_cities()` now returns an `sf` object rather than a data frame (#36).
- `us_cities()` gains a `states` argument to match other functions in the package (#35).
- Citation to JOSS paper.

# USAboundaries 0.4.0

- Update all data files to use current version of `sf` package.
- Update Census Bureau data from 2016 to 2020.
- Remove the `us_boundaries()` function was a needless wrapper around other functions.
- Prompt user to install data package rather than installing it for them.

# USAboundaries 0.3.0

- Moved most data to USAboundariesData. This improves loading time and permits more frequent updates to the user-facing package.
- Added state plane projections table and functions (@jfbratt). Now users can get an appropriate projection for a state.
- Converted all boundary objects to `sf` objects.
- Updated all contemporary census boundaries from the 2014 to the 2016 versions.
- Added zipcode tabulation area centroids.
- Added historical city populations compiled by Erik Steiner at CESTA/Stanford University.

# USAboundaries 0.2.0

-  Added contemporary boundaries for states, counties, and congressional districts.
-  Import many fewer packages. The `us_boundaries()` function no longer has an option to return a fortified data frame. It is assumed that users will convert the `SpatialPolygonsDataFrame` objects to whatever format they need.
-  High resolution data is now available in the USAboundariesData package.

# USAboundaries 0.1.1

-  Fix to README.md as requested by CRAN.

# USAboundaries 0.1

-   Initial release.
-   `us_boundaries()` returns an sp object or a data frame which can be
    plotted.
This is an update release of the 'USAboundaries' package, updating Census Bureau
data. It also fixes problems brought about by a new version of the `sf` package.

This is a resubmission following CRAN guidance. It updates the date field and corrects changed URLs as requested.

## Test environments

* local OS X install, R-release
* GitHub Actions, R-devel, R-release, R-oldrel 
* win-builder (R-devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There are several NOTEs about suggested package in an additional repository which adds extra data. This 'USAboundariesData' package is provided in the rOpenSci CRAN-style repository, and it is installed on user request.
---
output: github_document
pagetitle: Historical and Contemporary Boundaries of the United States of America
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)
set.seed(87737)
```

# USAboundaries

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/USAboundaries)](https://cran.r-project.org/package=USAboundaries)
[![JOSS Status](https://joss.theoj.org/papers/3458a33133aa6c069ab4dd8df0b5f3b5/status.svg)](https://doi.org/10.21105/joss.00314)
[![R-CMD-check](https://github.com/ropensci/USAboundaries/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/USAboundaries/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/USAboundaries/master.svg)](https://codecov.io/github/ropensci/USAboundaries?branch=master)

## Overview

This R package includes contemporary state, county, and Congressional district boundaries, as well as zip code tabulation area centroids. It also includes historical boundaries from 1629 to 2000 for states and counties from the Newberry Library's [Atlas of Historical County Boundaries](https://publications.newberry.org/ahcbp/), as well as historical city population data from Erik Steiner's "[United States Historical City Populations, 1790-2010](https://github.com/cestastanford/historical-us-city-populations)." The package has some helper data, including a table of state names, abbreviations, and FIPS codes, and functions and data to get [State Plane Coordinate System](https://en.wikipedia.org/wiki/State_Plane_Coordinate_System) projections as EPSG codes or PROJ.4 strings.

This package can serve a number of purposes. The spatial data can be joined to any other kind of data in order to make thematic maps. Unlike other R packages, this package also contains historical data for use in analyses of the recent or more distant past. See the ["A sample analysis using USAboundaries"](http://lincolnmullen.com/software/usaboundaries/articles/usaboundaries-sample-analysis.html) vignette for an example of how the package can be used for both historical and contemporary maps.

## Citation

If you use this package in your research, we would appreciate a citation.

```{r}
citation("USAboundaries")
```

## Installation

You can install this package from CRAN.

```
install.packages("USAboundaries")
```

Almost all of the data for this package is provided by the [USAboundariesData package](https://github.com/ropensci/USAboundariesData). That package will be automatically installed (with your permission) from the [rOpenSci package repository](https://ropensci.r-universe.dev) the first time that you need it.

Or you can install the development versions from GitHub using [remotes](https://remotes.r-lib.org). 

```
# install.packages("remotes")
remotes::install_github("ropensci/USAboundaries")
remotes::install_github("ropensci/USAboundariesData")
```

## Use

This package provides a set of functions, one for each of the types of boundaries that are available. These functions have a consistent interface. 

Passing a date to `us_states()`, `us_counties()`, and `us_cities()` returns the historical boundaries for that date. If no date argument is passed, then contemporary boundaries are returned. The functions `us_congressional()` and `us_zipcodes()` only offer contemporary boundaries.

For almost all functions, pass a character vector of state names or abbreviations to the `states =` argument to return only those states or territories.

For certain functions, more or less detailed boundary information is available by passing an argument to the `resolution =` argument.

See the examples below to see how the interface works, and see the documentation for each function for more details.

```{r}
library(USAboundaries) 
library(sf) # for plotting and projection methods

states_1840 <- us_states("1840-03-12")
plot(st_geometry(states_1840))
title("U.S. state boundaries on March 3, 1840")

states_contemporary <- us_states()
plot(st_geometry(states_contemporary))
title("Contemporary U.S. state boundaries")

counties_va_1787 <- us_counties("1787-09-17", states = "Virginia")
plot(st_geometry(counties_va_1787))
title("County boundaries in Virginia in 1787")

counties_va <- us_counties(states = "Virginia")
plot(st_geometry(counties_va))
title("Contemporary county boundaries in Virginia")

counties_va_highres <- us_counties(states = "Virginia", resolution = "high")
plot(st_geometry(counties_va_highres))
title("Higher resolution contemporary county boundaries in Virginia")

congress <- us_congressional(states = "California")
plot(st_geometry(congress))
title("Congressional district boundaries in California")
```

## State plane projections

The `state_plane()` function returns EPSG codes and PROJ.4 strings for the State Plane Coordinate System. You can use these to use suitable projections for specific states.

```{r}
va <- us_states(states = "VA", resolution = "high")
plot(st_geometry(va), graticule = TRUE)

va_projection <- state_plane("VA")
va <- st_transform(va, va_projection)
plot(st_geometry(va), graticule = TRUE)
```

## Related packages

Each function returns an `sf` object from the [sf](https://cran.r-project.org/package=sf) package, which can be mapped using the [leaflet](https://cran.r-project.org/package=leaflet) or [ggplot2](https://cran.r-project.org/package=ggplot2) packages.

If you need U.S. Census Bureau boundary files which are not provided by this package, consider using the [tigris](https://cran.r-project.org/package=tigris) package, which downloads those shapefiles.

## License

The historical boundary data provided in this package is available under the CC BY-NC-SA 2.5 license from John H. Long, et al., [Atlas of Historical County Boundaries](https://publications.newberry.org/ahcbp/), Dr. William M. Scholl Center for American History and Culture, The Newberry Library, Chicago (2010).  Please cite that project if you use this package in your research and abide by the terms of their license if you use the historical information.

The historical population data for cities is provided by U.S. Census Bureau and Erik Steiner, Spatial History Project, Center for Spatial and Textual Analysis, Stanford University. See the data in [this repository](https://github.com/cestastanford/historical-us-city-populations).

The contemporary data is provided by the U.S. Census Bureau and is in the public domain.

All code in this package is copyright [Lincoln Mullen](http://lincolnmullen.com) and is released under the MIT license.

---
[![rOpenSci footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
---
title: "A sample analysis using USAboundaries"
author: "Lincoln Mullen"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{A sample analysis using USAboundaries}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The USAboundaries package provides both historical and contemporary boundaries of the United States. This vignette gives an example of using the package to create maps of Maryland Congressional elections (to the 115th Congress) and in 1792 (to the 3rd Congress). We will use the USAboundaries package data for the spatial data. See the documentation for the functions used for citations to the underlying datasets. For the 2016 Congressional votes, we will use [data available on the website](http://history.house.gov/Institution/Election-Statistics/Election-Statistics/) of the U.S. House of Representatives. For the 1792 elections we will use data from the [*Mapping Early American Elections*](https://earlyamericanelections.org/data/) project.

We begin by loading the packages that we need.

```{r message=FALSE}
library(USAboundaries)
library(dplyr)
library(sf)
library(leaflet)

# Check that USAboundariesData is available, since it is not on CRAN
avail <- requireNamespace("USAboundariesData", quietly = TRUE)
```

We can load the voting data for the 115th Congress into `md_votes`. We use this package to get the contemporary (2016) boundaries for Congressional districts in Maryland. They are returned as an `sf` object from the [sf](https://github.com/r-spatial/sf/) package.

```{r}
if (avail) {
  md_votes <- read.csv(system.file("extdata", "md-115.csv",
                                   package = "USAboundaries"),
                       colClasses = c("character", "integer", "character", "integer",
                                      "character", "integer", "integer", "integer",
                                      "integer", "integer", "integer",  "numeric",
                                      "numeric", "numeric", "numeric", "numeric"))
  md_districts <- us_congressional(states = "MD")
}
```

In order to map the data, we have to join the voting data to the spatial data.

```{r}
if (avail) {
  md_115 <- md_districts %>% 
   left_join(md_votes, by = c("cd116fp" = "district")) 
}
```

The data has the percentage of votes for Democrats and Republicans. We are going to map the margin that Republicans had over an even vote, so that we can see both the Democrat and Republican results. For a more careful analysis we would have to take account of third parties and write-in votes, but this will serve our purpose.

```{r}
if (avail) {
  md_115$margin_percent <- md_115$percentage_republican - 0.5
}
```

We create a color palette to turn the margin of victory into red and blue colors.

```{r}
if (avail) {
  palette <- colorBin("RdBu", domain = c(-0.3, 0.3), bins = 7, reverse = TRUE)
}
```

We can use the [leaflet](https://rstudio.github.io/leaflet/) package to make the map.

```{r}
if (avail) {
  leaflet(md_115) %>% 
    addTiles() %>% 
    addPolygons(color = "black",
                opacity = 1,
                weight = 1,
                fillColor = ~palette(margin_percent),
                fillOpacity = 1,
                label = ~paste0("District: ", cd116fp))
}
```

We can do the same analysis for the 3rd Congress, thanks to the historical data included in the package. In this case, we will have to map counties rather than Congressional districts. The parties of interest are the Federalists and the Democratic-Republicans.

We start by loading and joining the data. 

```{r}
if (avail) {
  md_03_votes <- read.csv(system.file("extdata", "md-003.csv",
                                      package = "USAboundaries"),
                          stringsAsFactors = FALSE)
  md_03_districts <- us_counties(map_date = "1792-10-01",
                                 resolution = "high", states = "MD")
  md_003 <- md_03_districts %>% 
    left_join(md_03_votes, by = c("id" = "county_ahcb"))
  md_003$margin_percent <- md_003$federalist_percentage - 0.5
}
```

Then we create a new palette function and make the map.

```{r}
if (avail) {
  palette_03 <- colorBin("PRGn", domain = c(-0.5, 0.5), bins = 9)
  leaflet(md_003) %>% 
    addTiles() %>% 
    addPolygons(color = "black",
                opacity = 1,
                weight = 1,
                fillColor = ~palette_03(margin_percent),
                fillOpacity = 1,
                label = ~paste0("County: ", name))
}  
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/us_cities.R
\name{us_cities}
\alias{us_cities}
\title{City locations and populations (historical and contemporary)}
\usage{
us_cities(map_date = NULL, states = NULL)
}
\arguments{
\item{map_date}{If \code{NULL}, then city populations from the 2010 census
(the most recent census) are returned. This parameter accepts a \code{Date}
object or a character string coercible to a \code{Date} object, as well as
numeric values representing a year. If a year or date is used, then city
populations from the decennial census from 1790 to 2010 \emph{prior} to
that year is returned. For example, \code{1805} or \code{"1805-07-04"}
would return city populations from the 1800 census.}

\item{states}{A character vector of state or territory names or
abbreviations. Only boundaries for those states/territories will be
returned. If \code{NULL}, all boundaries will be returned.}
}
\description{
This function returns an \code{sf} object of cities (or populated places)
with their populations and latitudes and longitudes. Population data is taken
from the U.S. Census.
}
\examples{
if (require(USAboundariesData)) {
  us_cities(1805)
  us_cities("1828-05-08")
  us_cities()
}

}
\references{
The data was compiled by Erik Steiner and Jason Heppler at the
  Center for Spatial and Textual Analysis, Stanford University. See their
  \href{https://github.com/cestastanford/historical-us-city-populations}{description
   of the data} for a fuller accounting of how the data was gathered.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state_plane.R
\name{state_plane}
\alias{state_plane}
\title{Projections from the State Plane Coordinate System}
\usage{
state_plane(state, plane_id = NULL, type = c("epsg", "proj4"))
}
\arguments{
\item{state}{The postal code for the state.}

\item{plane_id}{The state plane geographic zone identifier. A \code{NULL}
value will return a projection for the entire state.}

\item{type}{The type of output to return: either an EPSG code or PROJ4
string.}
}
\value{
Either a PROJ4 string as a character vector or an EPSG code as an
  integer.
}
\description{
Get EPSG codes or PROJ.4 codes for projections from the
\href{https://en.wikipedia.org/wiki/State_Plane_Coordinate_System}{State
Plane Coordinate System}.
}
\examples{
if (require(USAboundariesData)) {
  state_plane(state = "MA", type = "epsg")
  state_plane(state = "MA", type = "proj4")
  state_plane(state = "MA", plane_id = "island", type = "epsg")
  state_plane(state = "MA", plane_id = "island", type = "proj4")

  # Show the difference made by a state plane projection
  if (require(sf)) {
    va <- us_states(states = "VA", resolution = "high")
    plot(st_geometry(va), graticule = TRUE)
    va <- st_transform(va, state_plane("VA"))
    plot(st_geometry(va), graticule = TRUE)
  }
}

}
\seealso{
For documentation of the underlying State Plane Coordinate System
  projection data frame, see \code{\link{state_proj}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/us_counties.R
\name{us_counties}
\alias{us_counties}
\title{County boundaries (contemporary and historical)}
\usage{
us_counties(map_date = NULL, resolution = c("low", "high"), states = NULL)
}
\arguments{
\item{map_date}{The date of the boundaries as some object coercible to a date
with \code{as.Date()}; the easiest option is a character vector following
the \href{https://en.wikipedia.org/wiki/ISO_8601}{ISO 8601} data format. If
\code{NULL} (the default) the contemporary boundaries will be returned.}

\item{resolution}{The resolution of the map.}

\item{states}{A character vector of state or territory names or
abbreviations. Only boundaries for those states/territories will be
returned. If \code{NULL}, all boundaries will be returned.}
}
\value{
An \code{sf} object.
}
\description{
Get the current (2020) boundaries for U.S counties from the U.S. Census
Bureau, or get historical county boundaries for dates between 30 December
1636 and 31 December 2000.
}
\examples{
if (require(USAboundariesData) && require(sf)) {
  contemporary_us  <- us_counties()
  historical_us    <- us_counties("1820-07-04")
  contemporary_ne  <- us_counties(states = c("Massachusetts", "Vermont", "Maine",
                                             "New Hampshire", "Rhode Island",
                                             "Connecticut"))
  historical_ne    <- us_counties("1803-04-28",
                                  states = c("Massachusetts", "Vermont", "Maine",
                                             "New Hampshire", "Rhode Island",
                                             "Connecticut"),
                                  resolution = "high")

  plot(st_geometry(contemporary_us))
  plot(st_geometry(historical_us))
  plot(st_geometry(contemporary_ne))
  plot(st_geometry(historical_ne))
}

}
\seealso{
For documentation of and citation to the underlying shapefiles for
  contemporary data from the U.S. Census Bureau, see the
  \code{census_boundaries} help file in the USAboundariesData package. For
  documentation of and citation to the underlying shapefiles for historical
  data from the Atlas of Historical County Boundaries, see the
  \code{ahcb_boundaries} help file in the USAboundariesData package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/usboundaries-package.r
\docType{package}
\name{USAboundaries}
\alias{USAboundaries}
\title{USAboundaries: Historical and Contemporary Boundaries of the United States of
America}
\description{
This package provides contemporary (2014) boundaries for states, counties,
zip code tabulation areas, and congressional districts in the United States
of America. This data is provided by the U.S. Census Bureau.
}
\details{
This package also provides spatial objects with historical boundaries of
states or counties in the United States of America from 1629 to 2000. It
provides data from the \href{https://publications.newberry.org/ahcbp/}{Atlas
of Historical County Boundaries}. The copyright to the historical data used
in this package is owned by the Newberry Library, and it is included in the
\code{USAboundariesData} package under the terms of the
\href{https://creativecommons.org/licenses/by-nc-sa/2.5/}{Creative Commons
Attribution-NonCommercial-ShareAlike 2.5 Generic} (CC BY-NC-SA 2.5) license.

The code in this package is copyrighted by
\href{https://lincolnmullen.com}{Lincoln Mullen}, and is released under the
terms of the \href{https://opensource.org/licenses/MIT}{MIT License}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/us_states.R
\name{us_states}
\alias{us_states}
\title{State boundaries (contemporary and historical)}
\usage{
us_states(map_date = NULL, resolution = c("low", "high"), states = NULL)
}
\arguments{
\item{map_date}{The date of the boundaries as some object coercible to a date
with \code{as.Date()}; the easiest option is a character vector following
the \href{https://en.wikipedia.org/wiki/ISO_8601}{ISO 8601} data format. If
\code{NULL} (the default) the contemporary boundaries will be returned.}

\item{resolution}{The resolution of the map.}

\item{states}{A character vector of state or territory names or
abbreviations. Only boundaries for those states/territories will be
returned. If \code{NULL}, all boundaries will be returned.}
}
\value{
An \code{sf} object.
}
\description{
Get the current (2020) boundaries for U.S states from the U.S. Census Bureau,
or get historical state boundaries for dates between 3 September 1783 and 31
December 2000.
}
\examples{
contemporary_us <- us_states()

if (require(USAboundariesData) && require(sf)) {
  historical_us   <- us_states("1820-07-04")
  contemporary_ne <- us_states(states = c("Massachusetts", "Vermont", "Maine",
                                          "New Hampshire", "Rhode Island",
                                          "Connecticut"))
  historical_ne   <- us_states(as.Date("1805-03-12"),
                               states = c("Massachusetts", "Vermont", "Maine",
                                          "New Hampshire", "Rhode Island",
                                          "Connecticut"),
                               resolution = "high")
   plot(st_geometry(contemporary_us))
   plot(st_geometry(historical_us))
   plot(st_geometry(contemporary_ne))
   plot(st_geometry(historical_ne))
}

}
\seealso{
For documentation of and citation to the underlying shapefiles for
  contemporary data from the U.S. Census Bureau, see \code{census_boundaries}
  documentation in the USAboundariesData package. For documentation of and
  citation to the underlying shapefiles for historical data from the Atlas of
  Historical County Boundaries, see the \code{ahcb_boundaries} documentation
  in the USAboundariesData package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install-data-pkg.R
\name{install_data_package}
\alias{install_data_package}
\title{Install the \code{USAboundariesData} package after checking with the user}
\usage{
install_data_package()
}
\description{
Install the \code{USAboundariesData} package after checking with the user
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/us_congressional.R
\name{us_congressional}
\alias{us_congressional}
\title{Congressional district boundaries (contemporary)}
\usage{
us_congressional(resolution = c("low", "high"), states = NULL)
}
\arguments{
\item{resolution}{The resolution of the boundaries.}

\item{states}{A character vector of state or territory names. Only boundaries
inside these states/territories will be returned. If \code{NULL}, all
boundaries will be returned.}
}
\value{
An \code{sf} object.
}
\description{
Get the current (2020) boundaries for U.S. Congressional districts.
}
\examples{
if (require(USAboundariesData) && require(sf)) {
  us_congressional <- us_congressional()
  va_congressional <- us_congressional(states = "Virginia", resolution = "high")
  plot(st_geometry(us_congressional))
  plot(st_geometry(va_congressional))
}

}
\seealso{
For documentation of and citation to the underlying shapefiles for
  contemporary data from the U.S. Census Bureau, see the
  \code{census_boundaries} help file in the USAboundariesData package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state_proj.R
\docType{data}
\name{state_proj}
\alias{state_proj}
\title{Data for projections from the State Plane Coordinate System}
\format{
A data frame with 123 rows and 5 variables:
\describe{
\item{state}{The state or territory abbreviation.}
\item{zone}{Name of the state plane zone.}
\item{epsg}{The EPSG code for each state plane zone.}
\item{proj4_string}{The PROJ4 string for the state plane projection.}
\item{statewide_proj}{State plane zone for projecting the entire state.}
}
}
\usage{
state_proj
}
\description{
This data frame includes state abbreviations, EPSG codes, and proj4 strings for projections from the State Plane Coordinate System.
}
\references{
\href{https://en.wikipedia.org/wiki/State_Plane_Coordinate_System}{State Plane Coordinate System}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state_codes.R
\docType{data}
\name{state_codes}
\alias{state_codes}
\title{State codes and abbreviations for U.S. states and territories}
\format{
A data.frame with 69 rows and 4 variables:
\describe{
\item{state_name}{The state or territory name}
\item{state_abbr}{The two character abbreviation for the state or territory.}
\item{state_code}{A three digit numeric FIPS code for the state or territory.}
\item{jurisdiction_type}{One of \code{state}, \code{territory}, or \code{district}.}
}
}
\usage{
state_codes
}
\description{
This data frame includes abbreviations and codes for states and territories
in the United States. It is intended as a lookup table.
}
\references{
U.S. Census Bureau,
  \href{https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.2014.html}{U.S. Gazeteer Files} (2014).

  \href{https://en.wikipedia.org/wiki/Federal_Information_Processing_Standard_state_code}{Federal
   Information Processing Standard state code}, Wikipedia (accessed July 23,
   2015).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-doc.R
\docType{data}
\name{states_contemporary_lores}
\alias{states_contemporary_lores}
\title{U.S. state boundaries}
\format{
An object of class \code{sf} (inherits from \code{data.frame}) with 52 rows and 13 columns.
}
\usage{
states_contemporary_lores
}
\description{
The U.S. Census Bureau provides
\href{https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html}{cartographic
boundary files} for current U.S. boundaries. This package has only the
low-resolution contemporary state boundaries. Other census boundary files are
provided by and documented in the USAboundariesData package.
}
\seealso{
For citations for the other Census boundary files provided by the
  USAboundariesData package, see the \code{census_boundaries} documentation
  in that package.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/us_zipcodes.R
\name{us_zipcodes}
\alias{us_zipcodes}
\title{Zip Code Tabulation Areas (contemporary)}
\usage{
us_zipcodes()
}
\value{
An \code{sf} object.
}
\description{
Get the current (2019) centroids for U.S Zipcode Tabulation Areas from the
U.S. Census Bureau. The centroids were calculated from the ZCTA boundary
files available on the U.S. Census Bureau website.
}
\examples{
if (require(USAboundariesData)) {
  us_zipcodes()
}

}
\seealso{
For documentation of and citation to the underlying shapefiles for
  contemporary data from the U.S. Census Bureau, see the
  \code{census_boundaries} documentation in the USAboundariesData package.
}
