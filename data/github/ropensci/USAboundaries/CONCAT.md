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
