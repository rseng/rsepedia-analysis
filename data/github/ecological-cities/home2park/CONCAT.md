
<!-- README.md is generated from README.Rmd. Please edit that file -->

# home2park: Spatial Provision of Urban Parks

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/ecological-cities/home2park/workflows/R-CMD-check/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![test-coverage](https://github.com/ecological-cities/home2park/workflows/test-coverage/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03609/status.svg)](https://doi.org/10.21105/joss.03609)
<!-- badges: end -->

<a href='https://ecological-cities.github.io/home2park/'><img src='man/figures/logo.png' align="right" height="175" /></a>

`home2park` is an R package for assessing the spatial provision of urban
parks to residential buildings city-wide. Refer to the [package
website](https://ecological-cities.github.io/home2park/) for
demonstrations of how the package may be used.

## Installation

Install the development version of `home2park` from GitHub:

``` r
devtools::install_github("ecological-cities/home2park", ref = "main")
```

## Setup

Load the package:

``` r
library(home2park)
```

## Citation

To cite `home2park` or acknowledge its use, please cite the following:

*Song, X. P., Chong, K. Y. (2021). home2park: An R package to assess the
spatial provision of urban parks. Journal of Open Source Software,
6(65), 3609. <https://doi.org/10.21105/joss.03609>*

<br>

The get a BibTex entry, run `citation("home2park")`.

<br>

## Background

Parks are important spaces for recreation and leisure in cities.
Conventional measures of park provision tend to rely on summaries of
park area within a given region (e.g. per capita park area). However,
there is a need to characterise the wide variety of parks (e.g. nature
areas, gardens, waterfront parks, outdoor playgrounds, etc.) that serve
different groups of people. When planning at fine spatial scales, such
current metrics are also limited by their coarse spatial resolution and
the presence of artificial boundaries.

<br>

## Using home2park

The package `home2park` provides a way to measure *multiple aspects* of
park provision to homes, at the resolution of *individual buildings*.
The key features include the ability to:

-   Download relevant data from OpenStreetMap (OSM) such as buildings,
    parks and features related to recreation. The user may also supply
    data for new buildings and parks, for the purpose of future scenario
    planning.
-   Redistribute coarse-scale population data per census unit into
    residential buildings, also known as ‘dasymetric mapping’, which
    helps highlight specific areas where more people will benefit from
    the presence of parks.
-   Summarise at each park multiple attributes that are important for
    recreation (e.g. dense vegetation, length of waterfronts, open
    spaces, trails, etc.).
-   Calculate the supply (provision) of the park attributes to each
    residential building, while accounting for ‘distance decay’, or the
    fact that supply from parks further away are reduced.

The following sections provide a high-level overview of the various
steps required to measure the spatial provision of parks. Further
details and code examples can be found in the package vignette ‘[Get
started](articles/home2park.html)’.

### 1. Process city population

Residential buildings (homes) are an important component of the
analysis. These may be obtained, for example, by downloading building
polygons from OpenStreetMap (OSM), and subsetting the dataset to areas
within ‘residential’ land use zones.

In addition, having the population count per residential building allows
us to calculate the total spatial provision of parks to all residents,
and can help highlight important areas where more people will benefit
from the presence of parks. Coarse-scale population census data can be
redistributed into the residential buildings, via a technique known as
‘dasymetric mapping’. The number of building ‘levels’ from OSM can be
used as a proxy for population density (i.e. more residents per unit
area). Here’s an example screenshot showing an overlay of multiple
example datasets in the package (for the city of Singapore), which were
used to redistribute population data per census unit (subzones) across
residential buildings.

<div class="figure" style="text-align: center">

<img src="man/figures/README-dasymetric-mapping.png" alt="Example screenshot showing an overlay of multiple datasets used to redistribute the population across buildings within residential land use zones. The legends are ordered (top to bottom) by increasing spatial resolution." width="80%" />
<p class="caption">
Example screenshot showing an overlay of multiple datasets used to
redistribute the population across buildings within residential land use
zones. The legends are ordered (top to bottom) by increasing spatial
resolution.
</p>

</div>

<br>

Residential building polygons in Singapore each with a population count
can be found in the following example dataset:

``` r
data(buildings_pop_sgp)
head(buildings_pop_sgp)
#> Simple feature collection with 6 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 103.8412 ymin: 1.297955 xmax: 103.8968 ymax: 1.426561
#> Geodetic CRS:  WGS 84
#>    popcount                       geometry
#> 1 32.367111 POLYGON ((103.8412 1.420676...
#> 2 20.880897 POLYGON ((103.8808 1.320019...
#> 3  2.311594 POLYGON ((103.8643 1.320283...
#> 4 56.747406 POLYGON ((103.8504 1.303723...
#> 5  6.368404 POLYGON ((103.8516 1.426561...
#> 6  9.502521 POLYGON ((103.8966 1.298226...
```

<br>

### 2. Process parks

Parks are the other important component of the analysis. These may be
downloaded from OSM and processed using this package. The following
example dataset contains parks in Singapore with selected attributes
related to recreation/leisure:

``` r
data(parks_sgp)
head(parks_sgp[, 28:33]) # subset to relevant columns
#> Simple feature collection with 6 features and 6 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 103.7809 ymin: 1.248586 xmax: 103.8704 ymax: 1.28586
#> Geodetic CRS:  WGS 84
#>               area     perimeter playground_count playground_ptdensity
#> 1   2454.365 [m^2]  346.5265 [m]                0            0 [1/m^2]
#> 2   1964.286 [m^2]  305.3844 [m]                0            0 [1/m^2]
#> 3 219319.203 [m^2] 3019.2241 [m]                0            0 [1/m^2]
#> 4  22513.834 [m^2]  615.1039 [m]                0            0 [1/m^2]
#> 5 571259.766 [m^2] 3973.5567 [m]                0            0 [1/m^2]
#> 6  67533.345 [m^2] 1092.2758 [m]                0            0 [1/m^2]
#>    trails_length trails_length_perim_ratio                       geometry
#> 1     0.0000 [m]             0.0000000 [1] MULTIPOLYGON (((103.8471 1....
#> 2     0.0000 [m]             0.0000000 [1] MULTIPOLYGON (((103.8446 1....
#> 3  4072.5481 [m]             1.3488724 [1] MULTIPOLYGON (((103.8059 1....
#> 4   837.9212 [m]             1.3622434 [1] MULTIPOLYGON (((103.8233 1....
#> 5 20947.6579 [m]             5.2717652 [1] MULTIPOLYGON (((103.8613 1....
#> 6   905.5080 [m]             0.8290105 [1] MULTIPOLYGON (((103.7812 1....
```

<br>

### 3. Recreation supply

With the processed building and park polygons, the provision of park
attributes per residential building can be calculated. The total supply
*S* of each park attribute to a building is calculated based on the
following equation. Its value depends on the distances between that
particular building and all parks; attributes from parks further away
are reduced as a result of the negative exponential function
*e<sup>-cd</sup>*, an effect also known as the ‘distance decay’ ([Rossi
et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008)).

<img src="man/figures/README-equation.png" width="18%" style="display: block; margin: auto;" />

where

-   *S* = Total supply of a specific park attribute to the building from
    parks *i*; *i* = 1,2,3, … *n* where *n* = total number of parks
    citywide.

-   *s*<sub>*i*</sub> = Supply of a specific park attribute from park
    *i*. A perfect positive linear association is assumed, since the
    focus is on supply metrics.

-   *d*<sub>*i*</sub> = Distance in kilometres from the building to park
    *i* (e.g. Euclidean, Manhattan, etc.).

-   *c* = Coefficient determining rate of decay in supply
    *s*<sub>*i*</sub> with increasing distance.

<br>

Note that the value of Coefficient *c* depends on both park and park
visitors’ attributes, such as socio-demographic factors and preferences
for activities that may impel shorter or longer travel ([Rossi et al.,
2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008); [Tu et
al.](https://doi.org/10.1016/j.ufug.2020.126689)). A lower value implies
that parks further away are accessible or frequently visited by
residents (i.e. still contributes to the ‘recreation supply’ of a
particular building).

<div class="figure" style="text-align: center">

<img src="man/figures/README-c-sensitivity.png" alt="Figure: The value of Coefficient c and its effect on the distance decay between a building and park." width="50%" />
<p class="caption">
Figure: The value of Coefficient c and its effect on the distance decay
between a building and park.
</p>

</div>

<br>

<div class="figure" style="text-align: center">

<img src="man/figures/README-c-sensitivity-map.png" alt="Screenshot: Examples showing the supply of OSM park area to residential buildings in Singapore for the year 2020 when the value of Coefficient c is 0.1 (left panel) and 1 (right panel). Each building is denoted as a point (a random subset is shown); the color palette is binned according to quantile values." width="100%" />
<p class="caption">
Screenshot: Examples showing the supply of OSM park area to residential
buildings in Singapore for the year 2020 when the value of Coefficient c
is 0.1 (left panel) and 1 (right panel). Each building is denoted as a
point (a random subset is shown); the color palette is binned according
to quantile values.
</p>

</div>

<br>

To calculate the supply of each park attribute, we first calculate the
pairwise distances between all buildings and parks (a distance matrix).
This output is supplied to the function `recre_supply()`. For example,
we can calculate the supply of park *area* to each building. This supply
value can then be multiplied by the population count per building, to
obtain the total supply to all residents.

``` r
# transform buildings & parks to projected crs
buildings_pop_sgp <- sf::st_transform(buildings_pop_sgp, sf::st_crs(32648))
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648))


# convert buildings to points (centroids), then calculate distances to every park
m_dist <- buildings_pop_sgp %>%
  sf::st_centroid() %>%
  sf::st_distance(parks_sgp) %>% # euclidean distance
    units::set_units(NULL)

m_dist <- m_dist / 1000 # convert distances to km


# new column for the supply of park area
buildings_pop_sgp$area_supply <- recre_supply(park_attribute = parks_sgp$area, 
                                              dist_matrix = m_dist, 
                                              c = 0.302) # e.g. from Tu et al. (2020)

# supply to all residents per building
buildings_pop_sgp$area_supplytopop <- buildings_pop_sgp$area_supply * buildings_pop_sgp$popcount
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-supply-parkarea-to-building-residents.png" alt="Screenshot: Supply of park area to building residents in Singapore based on OSM data (2020). Each building is denoted as a point (a random subset is shown). The value for Coefficient c was set at 0.302. The color palette is binned according to quantile values." width="100%" />
<p class="caption">
Screenshot: Supply of park area to building residents in Singapore based
on OSM data (2020). Each building is denoted as a point (a random subset
is shown). The value for Coefficient c was set at 0.302. The color
palette is binned according to quantile values.
</p>

</div>

<br>

## Data sources

-   Singapore census data from the [Department of Statistics
    Singapore](https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data).
    Released under the terms of the [Singapore Open Data Licence version
    1.0](https://data.gov.sg/open-data-licence).

-   Singapore subzone polygons from the [Singapore Master Plan
    Subzones](https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea).
    Released under the terms of the [Singapore Open Data Licence version
    1.0](https://data.gov.sg/open-data-licence).

-   Singapore Master Plan Land Use Zones for the years
    [2014](https://data.gov.sg/dataset/master-plan-2014-land-use) and
    [2019](https://data.gov.sg/dataset/master-plan-2019-land-use-layer).
    Released under the terms of the [Singapore Open Data
    License](https://data.gov.sg/open-data-licence).

-   Building polygons derived from map data
    [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap
    contributors and available from <https://www.openstreetmap.org>.
    Released under the terms of the [ODbL
    License](https://opendatacommons.org/licenses/odbl/summary/).

-   Park polygons and summarised attributes (trails, playgrounds)
    derived from map data
    [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap
    contributors and available from <https://www.openstreetmap.org>.
    Released under the terms of the [ODbL
    License](https://opendatacommons.org/licenses/odbl/summary/).

<br>

## References

Rossi, S. D., Byrne, J. A., & Pickering, C. M. (2015). The role of
distance in peri-urban national park use: Who visits them and how far do
they travel?. Applied Geography, 63, 77-88.

Tu, X., Huang, G., Wu, J., & Guo, X. (2020). How do travel distance and
park size influence urban park visits?. Urban Forestry & Urban Greening,
52, 126689.

<!-- README.md is generated from README.Rmd. Please edit that file -->

# home2park <a href='https://ecological-cities.github.io/home2park/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/ecological-cities/home2park/workflows/R-CMD-check/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![test-coverage](https://github.com/ecological-cities/home2park/workflows/test-coverage/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03609/status.svg)](https://doi.org/10.21105/joss.03609)
<!-- badges: end -->

## Spatial Provision of Urban Parks

Assess the spatial provision of urban parks to residential buildings
city-wide. Refer to [package
website](https://ecological-cities.github.io/home2park/) for
demonstrations of how the package may be used.

## Installation

Install the development version of `home2park` from GitHub:

``` r
devtools::install_github("ecological-cities/home2park", ref = "main")
```

## Setup

Load the package:

``` r
library(home2park)
```

## Citation

To cite `home2park` or acknowledge its use, please cite the following:

*Song, X. P., Chong, K. Y. (2021). home2park: An R package to assess the
spatial provision of urban parks. Journal of Open Source Software,
6(65), 3609. <https://doi.org/10.21105/joss.03609>*

<br>

To get a BibTex entry, run `citation("home2park")`.

<br>

## Background

Parks are important spaces for recreation and leisure in cities.
Conventional measures of park provision tend to rely on summaries of
park area within a given region (e.g. per capita park area). However,
there is a need to characterise the wide variety of parks (e.g. nature
areas, gardens, waterfront parks, outdoor playgrounds, etc.) that serve
different groups of people. When planning at fine spatial scales, such
current metrics are also limited by their coarse spatial resolution and
the presence of artificial boundaries.

<br>

## Using home2park

The package `home2park` provides a way to measure *multiple aspects* of
park provision to homes, at the resolution of *individual buildings*.
The key features include the ability to:

-   Download relevant data from OpenStreetMap (OSM) such as buildings,
    parks and features related to recreation. The user may also supply
    data for new buildings and parks, for the purpose of future scenario
    planning.
-   Redistribute coarse-scale population data per census unit into
    residential buildings, also known as ‘dasymetric mapping’, which
    helps highlight specific areas where more people will benefit from
    the presence of parks.
-   Summarise at each park multiple attributes that are important for
    recreation (e.g. dense vegetation, length of waterfronts, open
    spaces, trails, etc.).
-   Calculate the supply (provision) of the park attributes to each
    residential building, while accounting for ‘distance decay’, or the
    fact that supply from parks further away are reduced.

The following sections provide a high-level overview of the various
steps required to measure the spatial provision of parks. Further
details and code examples can be found in the package vignette ‘[Get
started](https://ecological-cities.github.io/home2park/articles/home2park.html)’.

### 1. Process city population

Residential buildings (homes) are an important component of the
analysis. These may be obtained, for example, by downloading building
polygons from OpenStreetMap (OSM), and subsetting the dataset to areas
within ‘residential’ land use zones.

In addition, having the population count per residential building allows
us to calculate the total spatial provision of parks to all residents,
and can help highlight important areas where more people will benefit
from the presence of parks. Coarse-scale population census data can be
redistributed into the residential buildings, via a technique known as
‘dasymetric mapping’. The number of building ‘levels’ from OSM can be
used as a proxy for population density (i.e. more residents per unit
area). Here’s an example screenshot showing an overlay of multiple
example datasets in the package (for the city of Singapore), which were
used to redistribute population data per census unit (subzones) across
residential buildings.

<div class="figure" style="text-align: center">

<img src="man/figures/README-dasymetric-mapping.png" alt="Example screenshot showing an overlay of multiple datasets used to redistribute the population across buildings within residential land use zones. The legends are ordered (top to bottom) by increasing spatial resolution." width="80%" />
<p class="caption">
Example screenshot showing an overlay of multiple datasets used to
redistribute the population across buildings within residential land use
zones. The legends are ordered (top to bottom) by increasing spatial
resolution.
</p>

</div>

<br>

Residential building polygons in Singapore each with a population count
can be found in the following example dataset:

``` r
data(buildings_pop_sgp)
head(buildings_pop_sgp)
#> Simple feature collection with 6 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 103.8412 ymin: 1.297955 xmax: 103.8968 ymax: 1.426561
#> Geodetic CRS:  WGS 84
#>    popcount                       geometry
#> 1 32.367111 POLYGON ((103.8412 1.420676...
#> 2 20.880897 POLYGON ((103.8808 1.320019...
#> 3  2.311594 POLYGON ((103.8643 1.320283...
#> 4 56.747406 POLYGON ((103.8504 1.303723...
#> 5  6.368404 POLYGON ((103.8516 1.426561...
#> 6  9.502521 POLYGON ((103.8966 1.298226...
```

<br>

### 2. Process parks

Parks are the other important component of the analysis. These may be
downloaded from OSM and processed using this package. The following
example dataset contains parks in Singapore with selected attributes
related to recreation/leisure:

``` r
data(parks_sgp)
head(parks_sgp[, 28:33]) # subset to relevant columns
#> Simple feature collection with 6 features and 6 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 103.7809 ymin: 1.248586 xmax: 103.8704 ymax: 1.28586
#> Geodetic CRS:  WGS 84
#>               area     perimeter playground_count playground_ptdensity
#> 1   2454.365 [m^2]  346.5265 [m]                0            0 [1/m^2]
#> 2   1964.286 [m^2]  305.3844 [m]                0            0 [1/m^2]
#> 3 219319.203 [m^2] 3019.2241 [m]                0            0 [1/m^2]
#> 4  22513.834 [m^2]  615.1039 [m]                0            0 [1/m^2]
#> 5 571259.766 [m^2] 3973.5567 [m]                0            0 [1/m^2]
#> 6  67533.345 [m^2] 1092.2758 [m]                0            0 [1/m^2]
#>    trails_length trails_length_perim_ratio                       geometry
#> 1     0.0000 [m]             0.0000000 [1] MULTIPOLYGON (((103.8471 1....
#> 2     0.0000 [m]             0.0000000 [1] MULTIPOLYGON (((103.8446 1....
#> 3  4072.5481 [m]             1.3488724 [1] MULTIPOLYGON (((103.8059 1....
#> 4   837.9212 [m]             1.3622434 [1] MULTIPOLYGON (((103.8233 1....
#> 5 20947.6579 [m]             5.2717652 [1] MULTIPOLYGON (((103.8613 1....
#> 6   905.5080 [m]             0.8290105 [1] MULTIPOLYGON (((103.7812 1....
```

<br>

### 3. Recreation supply

With the processed building and park polygons, the provision of park
attributes per residential building can be calculated. The total supply
*S* of each park attribute to a building is calculated based on the
following equation. Its value depends on the distances between that
particular building and all parks; attributes from parks further away
are reduced as a result of the negative exponential function
*e<sup>-cd</sup>*, an effect also known as the ‘distance decay’ ([Rossi
et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008)).

<img src="man/figures/README-equation.png" width="18%" style="display: block; margin: auto;" />

where

-   *S* = Total supply of a specific park attribute to the building from
    parks *i*; *i* = 1,2,3, … *n* where *n* = total number of parks
    citywide.

-   *s*<sub>*i*</sub> = Supply of a specific park attribute from park
    *i*. A perfect positive linear association is assumed, since the
    focus is on supply metrics.

-   *d*<sub>*i*</sub> = Distance in kilometres from the building to park
    *i* (e.g. Euclidean, Manhattan, etc.).

-   *c* = Coefficient determining rate of decay in supply
    *s*<sub>*i*</sub> with increasing distance.

<br>

Note that the value of Coefficient *c* depends on both park and park
visitors’ attributes, such as socio-demographic factors and preferences
for activities that may impel shorter or longer travel ([Rossi et al.,
2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008); [Tu et
al.](https://doi.org/10.1016/j.ufug.2020.126689)). A lower value implies
that parks further away are accessible or frequently visited by
residents (i.e. still contributes to the ‘recreation supply’ of a
particular building).

<div class="figure" style="text-align: center">

<img src="man/figures/README-c-sensitivity.png" alt="Figure: The value of Coefficient c and its effect on the distance decay between a building and park." width="50%" />
<p class="caption">
Figure: The value of Coefficient c and its effect on the distance decay
between a building and park.
</p>

</div>

<br>

<div class="figure" style="text-align: center">

<img src="man/figures/README-c-sensitivity-map.png" alt="Screenshot: Examples showing the supply of OSM park area to residential buildings in Singapore for the year 2020 when the value of Coefficient c is 0.1 (left panel) and 1 (right panel). Each building is denoted as a pount (a random subset is shown). The color palette for the buildings (points) is binned according to quantile values." width="100%" />
<p class="caption">
Screenshot: Examples showing the supply of OSM park area to residential
buildings in Singapore for the year 2020 when the value of Coefficient c
is 0.1 (left panel) and 1 (right panel). Each building is denoted as a
pount (a random subset is shown). The color palette for the buildings
(points) is binned according to quantile values.
</p>

</div>

<br>

To calculate the supply of each park attribute, we first calculate the
pairwise distances between all buildings and parks (a distance matrix).
This output is supplied to the function `recre_supply()`. For example,
we can calculate the supply of park *area* to each building. This supply
value can then be multiplied by the population count per building, to
obtain the total supply to all residents.

``` r
# transform buildings & parks to projected crs
buildings_pop_sgp <- sf::st_transform(buildings_pop_sgp, sf::st_crs(32648))
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648))


# convert buildings to points (centroids), then calculate distances to every park
m_dist <- buildings_pop_sgp %>%
  sf::st_centroid() %>%
  sf::st_distance(parks_sgp) %>% # euclidean distance
    units::set_units(NULL)

m_dist <- m_dist / 1000 # convert distances to km


# new column for the supply of park area
buildings_pop_sgp$area_supply <- recre_supply(park_attribute = parks_sgp$area, 
                                              dist_matrix = m_dist, 
                                              c = 0.302) # e.g. from Tu et al. (2020)

# supply to all residents per building
buildings_pop_sgp$area_supplytopop <- buildings_pop_sgp$area_supply * buildings_pop_sgp$popcount
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-supply-parkarea-to-building-residents.png" alt="Screenshot: Supply of park area to building residents in Singapore based on OSM data (2020). Each building is denoted as a point (a random subset is shown). The value for Coefficient c was set at 0.302. The color palette is binned according to quantile values." width="100%" />
<p class="caption">
Screenshot: Supply of park area to building residents in Singapore based
on OSM data (2020). Each building is denoted as a point (a random subset
is shown). The value for Coefficient c was set at 0.302. The color
palette is binned according to quantile values.
</p>

</div>

<br>

## Data sources

-   Singapore census data from the [Department of Statistics
    Singapore](https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data).
    Released under the terms of the [Singapore Open Data Licence version
    1.0](https://data.gov.sg/open-data-licence).

-   Singapore subzone polygons from the [Singapore Master Plan
    Subzones](https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea).
    Released under the terms of the [Singapore Open Data Licence version
    1.0](https://data.gov.sg/open-data-licence).

-   Singapore Master Plan Land Use Zones for the years
    [2014](https://data.gov.sg/dataset/master-plan-2014-land-use) and
    [2019](https://data.gov.sg/dataset/master-plan-2019-land-use-layer).
    Released under the terms of the [Singapore Open Data
    License](https://data.gov.sg/open-data-licence).

-   Building polygons derived from map data
    [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap
    contributors and available from <https://www.openstreetmap.org>.
    Released under the terms of the [ODbL
    License](https://opendatacommons.org/licenses/odbl/summary/).

-   Park polygons and summarised attributes (trails, playgrounds)
    derived from map data
    [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap
    contributors and available from <https://www.openstreetmap.org>.
    Released under the terms of the [ODbL
    License](https://opendatacommons.org/licenses/odbl/summary/).

<br>

## References

Rossi, S. D., Byrne, J. A., & Pickering, C. M. (2015). The role of
distance in peri-urban national park use: Who visits them and how far do
they travel?. *Applied Geography*, 63, 77-88.

Tu, X., Huang, G., Wu, J., & Guo, X. (2020). How do travel distance and
park size influence urban park visits?. *Urban Forestry & Urban
Greening*, 52, 126689.

# home2park (development version)

## Bug fixes 
* Remove lines in `rasterise_buildings()` that clear temporary files used by the `terra` package, as these may cause the function output to go missing later on (if export directory is not specified).


# home2park 0.1.1

## Bug fixes 
* Helper function `raster_class_area()`. The use of `terra::extract()` may output an additional third `area` column that is unused. Only keep first 2 columns of output.

* Fix error in `rasterise_buildings()`. Getting the files for land use rasters resulted in an error (e.g. multiple files). Syntax for re-naming for the output using the `glue::glue()` function also did not work within package environment.

* In `get_buildings_osm()`, add removal of invalid building polygons which remained in spite of `st_make_valid()`. Related to https://github.com/r-spatial/sf/issues/1649. Also, add arguments for `sf::st_write()` (e.g. `driver`, `delete_dsn`, `append`).

* In `rasterise_buildings_osm()`

    - Change default argument for 'year' to `NULL`, as multiple years may not necessarily be analysed. Add section to rasterise buildings without information on the 'year'. Assumes that only one population/landuse rasters were generated in pre-processing step too; if there are multiple, uses the first object for analysis.
    - Allow raster template (`terra::rast()` object) to be input directly as argument, or alternatively supply file path to the raster file.
    - Check if CRS are similar between supplied objects (e.g. similar crs between `sf_buildings` and `sf_pop` if supplied; similar crs between `sf_buildings` and `rastertemplate`).

* In `pop_dasymap()`

    - Remove argument `raster_template`, and replaced with raster from argument `land_relative_density`.

<br>

# home2park 0.1.0

* Finish development of main functions
* Prepare for journal submission
# Contributing to home2park

This outlines how to propose a change to home2park. 
For more detailed info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib). 

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("ecological-cities/home2park", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  
---
title: 'home2park: An R package to assess the spatial provision of urban parks'
tags:
  - R
  - Distance decay
  - Dasymetric map
  - Outdoor recreation
  - Park provision
  - Recreation supply
  - Service area
  - Spatial distribution
  - Urban green space
  - Urban parks
authors:
  - name: Xiao Ping Song
    orcid: 0000-0002-8825-195X
    affiliation: 1
  - name: Kwek Yan Chong
    orcid: 0000-0003-4754-8957
    affiliation: 1
affiliations:
 - name: Department of Biological Sciences, National University of Singapore
   index: 1
citation_author: Song and Chong
date: 11 Jul 2021
year: 2021
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

In the wake of the COVID-19 pandemic, significant changes to human mobility and working environments have prompted a re-think of land use distribution in the city, particularly for office, residential and public spaces such as parks. One prominent trend is a rising demand for the outdoors, as urban dwellers turn to local parks for recreation and leisure [@Venter2020]. However, current metrics of park provision tend to focus on summaries of park area within a given region [@Tan2017]. There is a need to characterise the wide variety of parks that serve different groups of people [@Song2020a], and to deal with the limitations of artificial boundaries when planning at fine spatial scales.

The package ``home2park`` provides a way to measure multiple aspects of park provision to homes, at the resolution of individual buildings. The key features include the ability to:  

1.	Download relevant data from OpenStreetMap (OSM) from any city worldwide (e.g. buildings, parks and features related to recreation), using the R package ``osmextract`` [@Lovelace2019]. The user may also supply data for new buildings and parks, for the purpose of future scenario planning. 
2.	Summarise at each park multiple attributes that are important for recreation (e.g. dense vegetation, length of waterfronts, open spaces, trails, etc.).
3.	Calculate the supply (provision) of the park attributes to residential buildings while accounting for ‘distance decay’ [@Rossi2015; @Tu2020b], or the fact that supply from parks further away are reduced.
4.	Redistribute coarse-scale population data per census unit into residential buildings, also known as ‘dasymetric mapping’ [@Dmowska], which helps highlight specific areas where more people will benefit from the presence of parks.

See the [package website](https://ecological-cities.github.io/home2park/) for detailed documentation and examples of how regional summaries derived from building-level metrics vary widely from conventional metrics of park area provision.

# Statement of Need

Numerous metrics have been used to measure the spatial provision of parks in cities. These include summaries of the per capita park area (i.e. park provision ratio) or ‘service radius’ (a distance buffer) around parks within geographical areal units such as census blocks or administrative boundaries [@Tan2017]. Such metrics have been used to investigate park use and accessibility, as well as social-cultural issues related to spatial equity and environmental justice [@Sister2010; @Tan2017]. However, the presence of artificial boundaries gives rise to the ‘modifiable areal unit problem’ , where differences in spatial scale can result in considerable variation in the value of summarised metrics [@Sister2010]. Since urban development typically occurs at fine spatial scales, metrics at the scale of individual buildings provided by ``home2park`` can help overcome such limitations [@Gao2017]. Finally, most metrics of park provision do not differentiate between the wide variety of parks that would likely affect usage by different groups of people. Other than park area, features such as waterfronts, open lawns, natural vegetation, trails and fitness amenities each relate to different types of park use and user groups [@Song2020a; @Song2020b]. A more holistic measure of park provision should thus include a variety of park attributes that are important for recreation and leisure. ``home2park`` enables the user to summarise a variety of spatial attributes such as points (e.g. playgrounds, sports and fitness amenities), lines (e.g. cycling and walking trails) and rasters (e.g. land cover types) at each park. Customisable parameters provide the user with flexibility when summarising specific attributes, such as specifying a minimum patch size for a land cover class, or including spatial features within a certain distance beyond park boundaries. By providing fine-grained measures for a variety of park attributes, we hope that ``home2park`` will contribute toward a more nuanced understanding of recreation that will improve the planning and design of parks and homes across the city.

## State of the Field

To our knowledge, there are currently no software (e.g. R packages) that measure the spatial provision of parks and recreation. 


# Acknowledgements

We thank Edwin Y. W. Tan and Justin K. H. Nai for useful discussions on the methodology. 

# References
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

# home2park <a href='https://ecological-cities.github.io/home2park/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/ecological-cities/home2park/workflows/R-CMD-check/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![test-coverage](https://github.com/ecological-cities/home2park/workflows/test-coverage/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03609/status.svg)](https://doi.org/10.21105/joss.03609)
<!-- badges: end -->

## Spatial Provision of Urban Parks

Assess the spatial provision of urban parks to residential buildings city-wide. Refer to [package website](https://ecological-cities.github.io/home2park/) for demonstrations of how the package may be used.

## Installation

Install the development version of `home2park` from GitHub:

```{r, eval = FALSE}
devtools::install_github("ecological-cities/home2park", ref = "main")
```

## Setup

Load the package:

```{r eval = FALSE}
library(home2park)
```

```{r include = FALSE, eval = TRUE}
devtools::load_all() # to knit manually w latest changes (not built/installed yet)
```

## Citation

To cite `home2park` or acknowledge its use, please cite the following:

*Song, X. P., Chong, K. Y. (2021). home2park: An R package to assess the spatial provision of urban parks. Journal of Open Source Software, 6(65), 3609. https://doi.org/10.21105/joss.03609*  


<br>

To get a BibTex entry, run `citation("home2park")`.

<br>

## Background

Parks are important spaces for recreation and leisure in cities. Conventional measures of park provision tend to rely on summaries of park area within a given region (e.g. per capita park area). However, there is a need to characterise the wide variety of parks (e.g. nature areas, gardens, waterfront parks, outdoor playgrounds, etc.) that serve different groups of people. When planning at fine spatial scales, such current metrics are also limited by their coarse spatial resolution and the presence of artificial boundaries.

<br>

## Using home2park

The package `home2park` provides a way to measure _multiple aspects_ of park provision to homes, at the resolution of _individual buildings_. The key features include the ability to:

- Download relevant data from OpenStreetMap (OSM) such as buildings, parks and features related to recreation. The user may also supply data for new buildings and parks, for the purpose of future scenario planning.
- Redistribute coarse-scale population data per census unit into residential buildings, also known as ‘dasymetric mapping’, which helps highlight specific areas where more people will benefit from the presence of parks.
- Summarise at each park multiple attributes that are important for recreation (e.g. dense vegetation, length of waterfronts, open spaces, trails, etc.).
- Calculate the supply (provision) of the park attributes to each residential building, while accounting for 'distance decay', or the fact that supply from parks further away are reduced.

The following sections provide a high-level overview of the various steps required to measure the spatial provision of parks. Further details and code examples can be found in the package vignette '[Get started](https://ecological-cities.github.io/home2park/articles/home2park.html)'.

### 1. Process city population

Residential buildings (homes) are an important component of the analysis. These may be obtained, for example, by downloading building polygons from OpenStreetMap (OSM), and subsetting the dataset to areas within 'residential' land use zones. 

In addition, having the population count per residential building allows us to calculate the total spatial provision of parks to all residents, and can help highlight important areas where more people will benefit from the presence of parks. Coarse-scale population census data can be redistributed into the residential buildings, via a technique known as 'dasymetric mapping'. The number of building 'levels' from OSM can be used as a proxy for population density (i.e. more residents per unit area). Here's an example screenshot showing an overlay of multiple example datasets in the package (for the city of Singapore), which were used to redistribute population data per census unit (subzones) across residential buildings. 

```{r, echo=FALSE, fig.align = "center", out.width="80%", fig.cap = paste0("Example screenshot showing an overlay of multiple datasets used to redistribute the population across buildings within residential land use zones. The legends are ordered (top to bottom) by increasing spatial resolution.")}
knitr::include_graphics("man/figures/README-dasymetric-mapping.png")
```

<br>

Residential building polygons in Singapore each with a population count can be found in the following example dataset:

```{r}
data(buildings_pop_sgp)
head(buildings_pop_sgp)
```


<br>

### 2. Process parks

Parks are the other important component of the analysis. These may be downloaded from OSM and processed using this package. The following example dataset contains parks in Singapore with selected attributes related to recreation/leisure:

```{r}
data(parks_sgp)
head(parks_sgp[, 28:33]) # subset to relevant columns
```

<br>

### 3. Recreation supply

With the processed building and park polygons, the provision of park attributes per residential building can be calculated. The total supply $S$ of each park attribute to a building is calculated based on the following equation. Its value depends on the distances between that particular building and all parks; attributes from parks further away are reduced as a result of the negative exponential function _e<sup>-cd</sup>_, an effect also known as the 'distance decay' ([Rossi et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008)).

```{r, echo=FALSE, fig.align = "center", out.width="18%"}
knitr::include_graphics("man/figures/README-equation.png")
```

where  

- $S$ = Total supply of a specific park attribute to the building from parks $i$; $i$ = 1,2,3, ... $n$ where $n$ = total number of parks citywide.

- $s_{i}$ = Supply of a specific park attribute from park $i$. A perfect positive linear association is assumed, since the focus is on supply metrics.  

- $d_{i}$ = Distance in kilometres from the building to park $i$ (e.g. Euclidean, Manhattan, etc.).

- $c$ = Coefficient determining rate of decay in supply $s_{i}$ with increasing distance. 

<br>

Note that the value of Coefficient $c$ depends on both park and park visitors' attributes, such as socio-demographic factors and preferences for activities that may impel shorter or longer travel ([Rossi et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008); [Tu et al.](https://doi.org/10.1016/j.ufug.2020.126689)). A lower value implies that parks further away are accessible or frequently visited by residents (i.e. still contributes to the 'recreation supply' of a particular building).


```{r, echo=FALSE, fig.align = "center", out.width="50%", fig.cap = paste0("Figure: The value of Coefficient c and its effect on the distance decay between a building and park.")}
knitr::include_graphics("man/figures/README-c-sensitivity.png")
```
<br>

```{r, echo=FALSE, fig.align = "center", out.width="100%", fig.cap = paste0("Screenshot: Examples showing the supply of OSM park area to residential buildings in Singapore for the year 2020 when the value of Coefficient c is 0.1 (left panel) and 1 (right panel). Each building is denoted as a pount (a random subset is shown). The color palette for the buildings (points) is binned according to quantile values.")}
knitr::include_graphics("man/figures/README-c-sensitivity-map.png")
```


<br>

To calculate the supply of each park attribute, we first calculate the pairwise distances between all buildings and parks (a distance matrix). This output is supplied to the function `recre_supply()`. For example, we can calculate the supply of park _area_ to each building. This supply value can then be multiplied by the population count per building, to obtain the total supply to all residents.

```{r calc distance matrix and supply of park area, message = FALSE, warning = FALSE}

# transform buildings & parks to projected crs
buildings_pop_sgp <- sf::st_transform(buildings_pop_sgp, sf::st_crs(32648))
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648))


# convert buildings to points (centroids), then calculate distances to every park
m_dist <- buildings_pop_sgp %>%
  sf::st_centroid() %>%
  sf::st_distance(parks_sgp) %>% # euclidean distance
    units::set_units(NULL)

m_dist <- m_dist / 1000 # convert distances to km


# new column for the supply of park area
buildings_pop_sgp$area_supply <- recre_supply(park_attribute = parks_sgp$area, 
                                              dist_matrix = m_dist, 
                                              c = 0.302) # e.g. from Tu et al. (2020)

# supply to all residents per building
buildings_pop_sgp$area_supplytopop <- buildings_pop_sgp$area_supply * buildings_pop_sgp$popcount

```

```{r, echo=FALSE, fig.align = "center", out.width="100%", fig.cap = paste0("Screenshot: Supply of park area to building residents in Singapore based on OSM data (2020). Each building is denoted as a point (a random subset is shown). The value for Coefficient c was set at 0.302. The color palette is binned according to quantile values.")}
knitr::include_graphics("man/figures/README-supply-parkarea-to-building-residents.png")
```


```{r include = FALSE}
rm(buildings_pop_sgp, parks_sgp, pop_sgp,
   m_dist)
```


<br>

## Data sources

- Singapore census data from the [Department of Statistics Singapore](https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore subzone polygons from the [Singapore Master Plan Subzones](https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore Master Plan Land Use Zones for the years [2014](https://data.gov.sg/dataset/master-plan-2014-land-use) and [2019](https://data.gov.sg/dataset/master-plan-2019-land-use-layer). Released under the terms of the [Singapore Open Data License](https://data.gov.sg/open-data-licence).

- Building polygons derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).

- Park polygons and summarised attributes (trails, playgrounds) derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).

<br>

## References

Rossi, S. D., Byrne, J. A., & Pickering, C. M. (2015). The role of distance in peri-urban national park use: Who visits them and how far do they travel?. _Applied Geography_, 63, 77-88.

Tu, X., Huang, G., Wu, J., & Guo, X. (2020). How do travel distance and park size influence urban park visits?. _Urban Forestry & Urban Greening_, 52, 126689.
---
output: github_document
always_allow_html: true
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/"
)
```

# home2park: Spatial Provision of Urban Parks
<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/ecological-cities/home2park/workflows/R-CMD-check/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![test-coverage](https://github.com/ecological-cities/home2park/workflows/test-coverage/badge.svg)](https://github.com/ecological-cities/home2park/actions)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03609/status.svg)](https://doi.org/10.21105/joss.03609)
<!-- badges: end -->

<a href='https://ecological-cities.github.io/home2park/'><img src='man/figures/logo.png' align="right" height="175" /></a>

`home2park` is an R package for assessing the spatial provision of urban parks to residential buildings city-wide. Refer to the [package website](https://ecological-cities.github.io/home2park/) for demonstrations of how the package may be used. 

## Installation

Install the development version of `home2park` from GitHub:

```{r, eval = FALSE}
devtools::install_github("ecological-cities/home2park", ref = "main")
```

## Setup

Load the package:

```{r eval = FALSE}
library(home2park)
```

```{r include = FALSE, eval = TRUE}
devtools::load_all() # to knit manually w latest changes (not built/installed yet)
```

## Citation

To cite `home2park` or acknowledge its use, please cite the following:

*Song, X. P., Chong, K. Y. (2021). home2park: An R package to assess the spatial provision of urban parks. Journal of Open Source Software, 6(65), 3609. https://doi.org/10.21105/joss.03609*  

<br>

The get a BibTex entry, run `citation("home2park")`.

<br>

## Background

Parks are important spaces for recreation and leisure in cities. Conventional measures of park provision tend to rely on summaries of park area within a given region (e.g. per capita park area). However, there is a need to characterise the wide variety of parks (e.g. nature areas, gardens, waterfront parks, outdoor playgrounds, etc.) that serve different groups of people. When planning at fine spatial scales, such current metrics are also limited by their coarse spatial resolution and the presence of artificial boundaries.

<br>

## Using home2park

The package `home2park` provides a way to measure _multiple aspects_ of park provision to homes, at the resolution of _individual buildings_. The key features include the ability to:

- Download relevant data from OpenStreetMap (OSM) such as buildings, parks and features related to recreation. The user may also supply data for new buildings and parks, for the purpose of future scenario planning.
- Redistribute coarse-scale population data per census unit into residential buildings, also known as ‘dasymetric mapping’, which helps highlight specific areas where more people will benefit from the presence of parks.
- Summarise at each park multiple attributes that are important for recreation (e.g. dense vegetation, length of waterfronts, open spaces, trails, etc.).
- Calculate the supply (provision) of the park attributes to each residential building, while accounting for 'distance decay', or the fact that supply from parks further away are reduced.

The following sections provide a high-level overview of the various steps required to measure the spatial provision of parks. Further details and code examples can be found in the package vignette '[Get started](articles/home2park.html)'.

### 1. Process city population

Residential buildings (homes) are an important component of the analysis. These may be obtained, for example, by downloading building polygons from OpenStreetMap (OSM), and subsetting the dataset to areas within 'residential' land use zones. 

In addition, having the population count per residential building allows us to calculate the total spatial provision of parks to all residents, and can help highlight important areas where more people will benefit from the presence of parks. Coarse-scale population census data can be redistributed into the residential buildings, via a technique known as 'dasymetric mapping'. The number of building 'levels' from OSM can be used as a proxy for population density (i.e. more residents per unit area). Here's an example screenshot showing an overlay of multiple example datasets in the package (for the city of Singapore), which were used to redistribute population data per census unit (subzones) across residential buildings. 

```{r, echo=FALSE, fig.align = "center", out.width="80%", fig.cap = paste0("Example screenshot showing an overlay of multiple datasets used to redistribute the population across buildings within residential land use zones. The legends are ordered (top to bottom) by increasing spatial resolution.")}
knitr::include_graphics("man/figures/README-dasymetric-mapping.png")
```

<br>

Residential building polygons in Singapore each with a population count can be found in the following example dataset:

```{r}
data(buildings_pop_sgp)
head(buildings_pop_sgp)
```


<br>

### 2. Process parks

Parks are the other important component of the analysis. These may be downloaded from OSM and processed using this package. The following example dataset contains parks in Singapore with selected attributes related to recreation/leisure:

```{r}
data(parks_sgp)
head(parks_sgp[, 28:33]) # subset to relevant columns
```

<br>

### 3. Recreation supply

With the processed building and park polygons, the provision of park attributes per residential building can be calculated. The total supply $S$ of each park attribute to a building is calculated based on the following equation. Its value depends on the distances between that particular building and all parks; attributes from parks further away are reduced as a result of the negative exponential function _e<sup>-cd</sup>_, an effect also known as the 'distance decay' ([Rossi et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008)).

```{r, echo=FALSE, fig.align = "center", out.width="18%"}
knitr::include_graphics("man/figures/README-equation.png")
```

where  

- $S$ = Total supply of a specific park attribute to the building from parks $i$; $i$ = 1,2,3, ... $n$ where $n$ = total number of parks citywide.

- $s_{i}$ = Supply of a specific park attribute from park $i$. A perfect positive linear association is assumed, since the focus is on supply metrics.  

- $d_{i}$ = Distance in kilometres from the building to park $i$ (e.g. Euclidean, Manhattan, etc.).

- $c$ = Coefficient determining rate of decay in supply $s_{i}$ with increasing distance. 

<br>

Note that the value of Coefficient $c$ depends on both park and park visitors' attributes, such as socio-demographic factors and preferences for activities that may impel shorter or longer travel ([Rossi et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008); [Tu et al.](https://doi.org/10.1016/j.ufug.2020.126689)). A lower value implies that parks further away are accessible or frequently visited by residents (i.e. still contributes to the 'recreation supply' of a particular building).


```{r, echo=FALSE, fig.align = "center", out.width="50%", fig.cap = paste0("Figure: The value of Coefficient c and its effect on the distance decay between a building and park.")}
knitr::include_graphics("man/figures/README-c-sensitivity.png")
```
<br>

```{r, echo=FALSE, fig.align = "center", out.width="100%", fig.cap = paste0("Screenshot: Examples showing the supply of OSM park area to residential buildings in Singapore for the year 2020 when the value of Coefficient c is 0.1 (left panel) and 1 (right panel). Each building is denoted as a point (a random subset is shown); the color palette is binned according to quantile values.")}
knitr::include_graphics("man/figures/README-c-sensitivity-map.png")
```


<br>

To calculate the supply of each park attribute, we first calculate the pairwise distances between all buildings and parks (a distance matrix). This output is supplied to the function `recre_supply()`. For example, we can calculate the supply of park _area_ to each building. This supply value can then be multiplied by the population count per building, to obtain the total supply to all residents.

```{r calc distance matrix and supply of park area, message = FALSE, warning = FALSE}

# transform buildings & parks to projected crs
buildings_pop_sgp <- sf::st_transform(buildings_pop_sgp, sf::st_crs(32648))
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648))


# convert buildings to points (centroids), then calculate distances to every park
m_dist <- buildings_pop_sgp %>%
  sf::st_centroid() %>%
  sf::st_distance(parks_sgp) %>% # euclidean distance
    units::set_units(NULL)

m_dist <- m_dist / 1000 # convert distances to km


# new column for the supply of park area
buildings_pop_sgp$area_supply <- recre_supply(park_attribute = parks_sgp$area, 
                                              dist_matrix = m_dist, 
                                              c = 0.302) # e.g. from Tu et al. (2020)

# supply to all residents per building
buildings_pop_sgp$area_supplytopop <- buildings_pop_sgp$area_supply * buildings_pop_sgp$popcount

```

```{r, echo=FALSE, fig.align = "center", out.width="100%", fig.cap = paste0("Screenshot: Supply of park area to building residents in Singapore based on OSM data (2020). Each building is denoted as a point (a random subset is shown). The value for Coefficient c was set at 0.302. The color palette is binned according to quantile values.")}
knitr::include_graphics("man/figures/README-supply-parkarea-to-building-residents.png")
```


```{r include = FALSE}
rm(buildings_pop_sgp, parks_sgp, pop_sgp,
   m_dist)
```


<br>

## Data sources

- Singapore census data from the [Department of Statistics Singapore](https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore subzone polygons from the [Singapore Master Plan Subzones](https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore Master Plan Land Use Zones for the years [2014](https://data.gov.sg/dataset/master-plan-2014-land-use) and [2019](https://data.gov.sg/dataset/master-plan-2019-land-use-layer). Released under the terms of the [Singapore Open Data License](https://data.gov.sg/open-data-licence).

- Building polygons derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).

- Park polygons and summarised attributes (trails, playgrounds) derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).

<br>

## References

Rossi, S. D., Byrne, J. A., & Pickering, C. M. (2015). The role of distance in peri-urban national park use: Who visits them and how far do they travel?. Applied Geography, 63, 77-88.

Tu, X., Huang, G., Wu, J., & Guo, X. (2020). How do travel distance and park size influence urban park visits?. Urban Forestry & Urban Greening, 52, 126689.
---
title: "Get started with home2park"
author: "Song, Xiao Ping"
date: "`r Sys.Date()`"
opengraph:
  image: 
    src: "man/figures/logo.png"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get started with home2park}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Overview

Parks are important spaces for recreation and leisure in cities. Conventional metrics that measure the provision of parks to city residents usually summarise the park area within a given region (e.g. per capita park area). However, 
there is a need to characterise the wide variety of parks (e.g. nature areas, gardens, waterfront parks, outdoor playgrounds, etc.) that serve different groups of people. When planning at fine spatial scales, such metrics are also limited by their coarse spatial resolution and the presence of artificial boundaries.

The package `home2park` provides a way to measure _multiple aspects_ of park provision to homes, at the resolution of _individual buildings_. The key features include the ability to:

- Download relevant data from OpenStreetMap (OSM) such as buildings, parks and features related to recreation. The user may also supply data for new buildings and parks, for the purpose of future scenario planning.
- Redistribute coarse-scale population data per census unit into residential buildings, also known as ‘dasymetric mapping’, which helps highlight specific areas where more people will benefit from the presence of parks.
- Summarise at each park multiple attributes that are important for recreation (e.g. dense vegetation, length of waterfronts, open spaces, trails, etc.).
- Calculate the supply (provision) of the park attributes to each residential building, while accounting for 'distance decay', or the fact that supply from parks further away are reduced.

## Installation

Install the development version of `home2park` from GitHub:

```{r, eval = FALSE}
devtools::install_github("ecological-cities/home2park", ref = "main")
```

## Setup

Load the package:

```{r eval = TRUE}
library(home2park)
```

```{r include = FALSE, eval = FALSE}
devtools::load_all() # to knit manually w latest changes (not built/installed yet)
# vignette will be checked in R-CMD check & package build (should use library(home2park))
```

<br>

## 1. Process city population

The first step is to process data related to the residential population (e.g. population census counts, residential buildings), in order to get the population count per building. The package comes with example data from the city of Singapore. This includes [population counts](https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data) per [census unit](https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea) in the years 2018 and 2020, and land use zones based on the Master Plan released in the years [2014](https://data.gov.sg/dataset/master-plan-2014-land-use) and [2019](https://data.gov.sg/dataset/master-plan-2019-land-use-layer). To get the residential buildings for the population census conducted in 2020, we can download building polygons from OpenStreetMap on 2021-01-01, and subset them to areas within 'residential' land use zones in the 2019 Master Plan. The processed polygons can then be rasterised in preparation for dasymetric mapping.

```{r load data, eval = FALSE}

data(pop_sgp) # population count per census unit (years 2018 & 2020)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs
  

# merge all census units for chosen year (2020) into single multi-polygon 
# next function requires that polygons are merged
city_boundaries <- pop_sgp %>%
    dplyr::filter(year == 2020) %>%
    sf::st_union() %>% 
    sf::st_as_sf() %>%
    smoothr::fill_holes(threshold = units::set_units(1, "km^2"))  %>% # clean up
    smoothr::drop_crumbs(threshold = units::set_units(1, "km^2"))  %>%
    sf::st_make_valid()


# get buildings from OpenStreetMap (use closest matching archive date)
# raw data saved in tempdir() unless otherwise specified
buildings <- get_buildings_osm(place = city_boundaries,
                               date = as.Date(c("2021-01-01"))) %>%
  dplyr::mutate(year = 2020)



# rasterise population - specify relevant column names
pop_rasters <- rasterise_pop(pop_sgp,
                             census_block = "subzone_n",
                             pop_count = "pop_count",
                             year = "year")


# rasterise land use - specify relevant column names and residential land use classes
data(landuse_sgp) 
landuse_sgp <- sf::st_transform(landuse_sgp, sf::st_crs(32648)) # transform to projected crs

landuse_rasters <- rasterise_landuse(landuse_sgp,
                                     land_use = 'lu_desc',
                                     subset = c('1' = 'RESIDENTIAL',
                                                '2' = 'COMMERCIAL & RESIDENTIAL',
                                                '3' = 'RESIDENTIAL WITH COMMERCIAL AT 1ST STOREY',
                                                '4' = 'RESIDENTIAL / INSTITUTION'),
                                     sf_pop = pop_sgp)


# rasterise buildings - specify relevant column names
# the column 'levels' is a proxy for how densely populated the building can be (i.e. more residents per unit area)
buildings_rasters <- rasterise_buildings(buildings,
                                         proxy_pop_density = 'levels',
                                         year = 'year',
                                         sf_pop = pop_sgp,
                                         sf_landuse = landuse_sgp)

```

Next, perform dasymetric mapping to redistribute the population counts per census unit region (year 2020) into the residential buildings. The number of building 'levels' from OSM can be used as a proxy for population density (i.e. more residents per unit area).  

```{r dasymetric mapping, eval = FALSE}

# list element 2 is for the year 2020
popdens_raster <- pop_dasymap(pop_polygons = pop_rasters$pop_polygons[[2]],
                               pop_perblock_count = pop_rasters$pop_count[[2]],
                               pop_perblock_density = pop_rasters$pop_density[[2]],
                               land_relative_density = buildings_rasters[[2]])

```

<br>

Here's a screenshot showing an overlay of the multiple datasets used to redistribute the population data per census unit (subzones) across residential buildings. 

```{r, echo=FALSE, fig.align = "center", out.width="80%", fig.cap = paste0("Example screenshot showing an overlay of multiple datasets used to redistribute the population across buildings within residential land use zones. The legends are ordered (top to bottom) by increasing spatial resolution.")}
knitr::include_graphics("../man/figures/README-dasymetric-mapping.png")
```

<br>

Finally, we can convert this raster to building polygons each with a population count. Note that this is provided as an example dataset `data(buildings_pop_sgp)` in this package.

```{r eval = FALSE}
buildings_pop_sgp <- pop_density_polygonise(input_raster = popdens_raster)

head(buildings_pop_sgp)
```

```{r echo = FALSE}
data(buildings_pop_sgp)
buildings_pop_sgp <- sf::st_transform(buildings_pop_sgp, sf::st_crs(32648)) # transform to projected crs
head(buildings_pop_sgp)
```

<br>

## 2. Process parks

The example dataset `data(parks_sgp)` contains parks in Singapore and selected attributes related to recreation/leisure (downloaded from OSM on 2021-01-01). This was downloaded and processed as follows. Parks within other cities (provided polygon boundaries) can be downloaded and processed in a similar manner.

```{r get park and their attributes, eval = FALSE}

# get park polygons
parks_sgp <- get_parks_osm(city_boundaries,
                           date = as.Date(c("2021-01-01")))


# get playground points
playgrounds <- get_playgrounds_osm(place = city_boundaries,
                                   date = as.Date(c("2021-01-01")))

# get trail lines
trails <- get_trails_osm(place = city_boundaries,
                         date = as.Date(c("2021-01-01")))


# convert to lists
point_list <- list(playgrounds) # calculation of attributes requires input to be named list 
names(point_list) <- c("playground")

line_list <- list(trails)
names(line_list) <- c("trails")


# per park, calculate point count/density & line length/ratio 
parks_sgp <- parks_calc_attributes(parks = parks_sgp, 
                                   data_points = point_list,
                                   data_lines = line_list,
                                   relative = TRUE)

head(parks_sgp[, 28:33]) # view park attributes
```

```{r echo = FALSE}
data(parks_sgp)
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648)) # transform to projected crs
head(parks_sgp[, 28:33])
```

<br>

## 3. Recreation supply

The processed data can now be used to calculate the provision of park attributes per residential building. The total supply $S$ of each park attribute to a building is calculated based on the following equation. Its value depends on the distances between that particular building and all parks; supply from parks further away are generally reduced as a result of the negative exponential function _e<sup>-cd</sup>_, an effect also known as the 'distance decay' ([Rossi et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008)).

```{r, echo=FALSE, fig.align = "center", out.width="18%"}
knitr::include_graphics("../man/figures/README-equation.png")
```

where  

- $S$ = Total supply of a specific park attribute to the building from parks $i$; $i$ = 1,2,3, ... $n$ where $n$ = total number of parks citywide.

- $s_{i}$ = Supply of a specific park attribute from park $i$. A perfect positive linear association is assumed, since the focus is on supply metrics.  

- $d_{i}$ = Distance in kilometres from the building to park $i$ (e.g. Euclidean, Manhattan, etc.).

- $c$ = Coefficient determining rate of decay in supply $s_{i}$ with increasing distance. 

<br>

The value of coefficient $c$ depends on both park and park visitors' attributes, such as socio-demographic factors and preferences for activities that may impel shorter or longer travel ([Rossi et al., 2015](http://dx.doi.org/10.1016/j.apgeog.2015.06.008)). For example, (non-Euclidean) distance decay in the ratio of visitors at urban parks in Beijing in [Tu et al. (2020)](https://doi.org/10.1016/j.ufug.2020.126689) fit the negative exponential curve _e<sup>-cd</sup>_, with coefficient $c$ empirically determined as _0.302_. Empirical cut-off distances were 1–2km where most (> 50%) people would visit a park and 5–10km beyond which few (< 5%) people would visit a park. 


```{r, echo=FALSE, fig.align = "center", out.width="50%", fig.cap = paste0("Figure: The value of Coefficient c and its effect on the distance decay between a building and park.")}
knitr::include_graphics("../man/figures/README-c-sensitivity.png")
```

<br> 

```{r, echo=FALSE, fig.align = "center", out.width="100%", fig.cap = paste0("Screenshot: Examples showing the supply of OSM park area to residential buildings in Singapore for the year 2020 when the value of Coefficient c is 0.1 (left panel) and 1 (right panel). Each building is denoted as a pount (a random subset is shown). The color palette for the buildings (points) is binned according to the quantile values.")}
knitr::include_graphics("../man/figures/README-c-sensitivity-map.png")
```


<br>

To calculate the supply of each park attribute, we first calculate the pairwise distances between all buildings and parks (a distance matrix). This output is supplied to the function `recre_supply()`. For example, we can calculate the supply of park _area_ to each building. This supply value can then be multiplied by the population count per building, to obtain the total supply to all residents. This can help highlight important areas where more people would benefit from the presence of parks. 

```{r calc distance matrix and supply of park area, message = FALSE, warning = FALSE}

# convert buildings to points (centroids), then calculate distances to every park
m_dist <- buildings_pop_sgp %>%
  sf::st_centroid() %>%
  sf::st_distance(parks_sgp) %>% # euclidean distance
    units::set_units(NULL)

m_dist <- m_dist / 1000 # convert distances to km


# new column for the supply of park area
buildings_pop_sgp$area_supply <- recre_supply(park_attribute = parks_sgp$area, 
                                              dist_matrix = m_dist, 
                                              c = 0.302) # e.g. from Tu et al. (2020)


# supply to all residents per building
buildings_pop_sgp$area_supplytopop <- buildings_pop_sgp$area_supply * buildings_pop_sgp$popcount
```

```{r include = FALSE}
rm(m_dist)
```

We can visualise the supply of park area to building residents using `library(tmap)`:  

```{r code to plot using building-level metrics using tmap, eval = FALSE}

# convert buildings to points (centroids) for plotting
buildings_pop_sgp <- buildings_pop_sgp %>% 
  sf::st_centroid() %>%
  dplyr::mutate(across(everything(), as.vector)) %>% # remove attributes from columns for plotting
  dplyr::mutate(across(.cols = contains("area"), # convert from m2 to km2
              .fns = function(x) x*1e-6))


tmap::tmap_mode("view")

tmap::tm_basemap("Esri.WorldGrayCanvas") +
  tmap::tm_shape(parks_sgp) +
    tmap::tm_polygons(col = "#33a02c",
                      alpha = 0.6,
                      border.col = "grey50", 
                      border.alpha = 0.5) +
  tmap::tm_shape(buildings_pop_sgp) +
    tmap::tm_dots(title = "Supply of park area (km<sup>2</sup>)",
                  col = "area_supplytopop", 
                  border.col = "transparent",
                  palette = viridis::viridis(5),
                  style = "quantile",
                  size = 0.01,
                  alpha = 0.7,
                  showNA = FALSE)

```


```{r, echo=FALSE, fig.align = "center", out.width="100%", fig.cap = paste0("Screenshot: Supply of park area to building residents in Singapore based on OSM data (2020). Each building is denoted as a point (a random subset is shown). The color palette is binned according to the quantile values.")}
knitr::include_graphics("../man/figures/README-supply-parkarea-to-building-residents.png") 
```

<br>

Finally, it is also possible to summarise the supply values of buildings within coarser spatial regions. This allows us to compare our results with conventional park provision metrics (per capita park area) for the selected year 2020.

```{r echo = FALSE}
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs
```

```{r}
# calculate for year 2020
pop_2020 <- pop_sgp %>%
  dplyr::filter(year == 2020)
```

Calculate summary statistics (e.g. sum, median, mean) for the supply of park area within each region:

```{r summarise supply of park area per subzone, message = FALSE, warning = FALSE}

# append subzone info to buildings if building intersects with subzone  
buildings_subzones <- buildings_pop_sgp  %>%
  sf::st_join(pop_2020,
              join = sf::st_intersects,
              left = TRUE) %>%
  sf::st_set_geometry(NULL)


# summarise total/median/average supply value for all buildings within subzone
buildings_subzones <- buildings_subzones %>%
  dplyr::group_by(subzone_n) %>%
  dplyr::summarise(across(.cols = c(area_supply, area_supplytopop), 
                          .fns = sum, .names = "{.col}_sum"),
                   across(.cols = c(area_supply, area_supplytopop), 
                          .fns = median, .names = "{.col}_median"),
                   across(.cols = c(area_supply, area_supplytopop), 
                          .fns = mean, .names = "{.col}_mean")) 
  

# join information to pop_2020, calculate supply per capita
pop_2020 <- pop_2020 %>%
  dplyr::left_join(buildings_subzones) %>%
  dplyr::mutate(area_supplyperpop = area_supplytopop_sum / pop_count * 1e-6) %>% # convert to km2
  dplyr::mutate(area_supplyperpop = ifelse(is.infinite(area_supplyperpop), NA, area_supplyperpop)) # make infinite NA

```

Calculate the (conventional) per capita park area for each census unit:

```{r include = FALSE}
sf::st_agr(pop_2020) = "constant"
sf::st_agr(parks_sgp) = "constant"
# https://github.com/r-spatial/sf/issues/406
```


```{r conventional park provision metric - per capita park area, message = FALSE, warning = FALSE}

subzones_parks <- sf::st_intersection(pop_2020, parks_sgp) # subset census unit polygons to those that intersect parks, append park info
subzones_parks$parkarea_m2 <- sf::st_area(subzones_parks) # calc area of each polygon


# calculate total park area per census unit
subzones_parks <- subzones_parks %>%
  dplyr::group_by(subzone_n) %>%
  dplyr::summarise(parkarea_m2 = as.numeric(sum(parkarea_m2))) %>%
  sf::st_drop_geometry()


# join information to pop_2020, calculate park area per capita
pop_2020 <- pop_2020 %>%
  dplyr::left_join(subzones_parks) %>%
  dplyr::mutate(parkarea_m2 = ifelse(is.na(parkarea_m2), 0, parkarea_m2)) %>% # make NA 0
  dplyr::mutate(parkperpop_m2 = parkarea_m2 / pop_count * 1e-6) %>% # convert to km2
  dplyr::mutate(parkperpop_m2 = ifelse(is.infinite(parkperpop_m2), # make infinite NA
                                NA, parkperpop_m2)) %>%
  dplyr::mutate(parkperpop_m2 = ifelse(is.nan(parkperpop_m2), # make NaN NA
                                NA, parkperpop_m2))

```

```{r include = FALSE}
rm(buildings_subzones, subzones_parks)
```

Similarly, we can visualise the park metrics summarised per region:  

```{r code to plot subzone-level metrics using tmap, eval = FALSE}

tmap::tm_basemap("Esri.WorldGrayCanvas") +
  tmap::tm_shape(pop_2020) +
    tmap::tm_polygons(title = "Park area per capita (km<sup>2</sup>)",
                      group = "Subzones: Per capita park area (conventional)",
                      col = "parkperpop_m2",
                      palette = "Greens", alpha = 0.7,
                      style = "quantile",
                      border.col = "white", border.alpha = 0.5, lwd = 1) +
  tmap::tm_shape(pop_2020) +
    tmap::tm_polygons(title = "Supply of park area per capita (km<sup>2</sup>)",
                      group = "Subzones: Per capita supply of park area",
                      col = "area_supplyperpop", 
                      palette = "Greens", alpha = 0.7,
                      style = "quantile",
                      border.col = "white", border.alpha = 0.5, lwd = 1)

```

<br>

The following figure provides a visual comparison between the provision of park area per capita summarised per region (subzone), based on conventional metrics (left panel) and derived from each residential building (right panel). The interactive map is viewable online at https://ecological-cities.github.io/home2park/articles/online/metric-comparisons.html.


```{r, echo=FALSE, fig.align = "center", out.width="100%", fig.cap = paste0("Screenshot: Provision of park area per capita within each subzone in Singapore, as measured by the total park area per subzone (left panel); and by the supply of park area to residents within residential buildings (right panel). Based on OSM data (2020). For the calculation of the supply value per building, the value for Coefficient c was set at 0.302. The color palettes are binned according to quantile values.")}
knitr::include_graphics("../man/figures/vignette-subzone-park-area.png")
```


```{r include = FALSE}
rm(tm, buildings_pop_sgp, parks_sgp, pop_sgp, pop_2020)
```


<br>

## Data sources

- Singapore census data from the [Department of Statistics Singapore](https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore subzone polygons from the [Singapore Master Plan Subzones](https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore Master Plan Land Use Zones for the years [2014](https://data.gov.sg/dataset/master-plan-2014-land-use) and [2019](https://data.gov.sg/dataset/master-plan-2019-land-use-layer). Released under the terms of the [Singapore Open Data License](https://data.gov.sg/open-data-licence).

- Building polygons derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).

- Park polygons and summarised attributes (trails, playgrounds) derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).

<br>

## References

Rossi, S. D., Byrne, J. A., & Pickering, C. M. (2015). The role of distance in peri-urban national park use: Who visits them and how far do they travel?. Applied Geography, 63, 77-88.

Tu, X., Huang, G., Wu, J., & Guo, X. (2020). How do travel distance and park size influence urban park visits?. Urban Forestry & Urban Greening, 52, 126689.


---
title: "Park area provision: Metric comparisons"
author: "Song, Xiao Ping"
date: "`r Sys.Date()`"
opengraph:
  image: 
    src: "man/figures/logo.png"
output: rmarkdown::html_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r eval = TRUE, include = FALSE}
library(home2park)
```

```{r include = FALSE, eval = FALSE}
devtools::load_all() # to knit manually w latest changes (not built/installed yet)
# vignette will be checked in R-CMD check & package build (should use library(home2park))
```

```{r load data, include = FALSE}
data(buildings_pop_sgp)
buildings_pop_sgp <- sf::st_transform(buildings_pop_sgp, sf::st_crs(32648)) # transform to projected crs

data(parks_sgp)
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648)) # transform to projected crs

data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs
```

```{r calc distance matrix and supply of park area, message = FALSE, warning = FALSE, include = FALSE}

# convert buildings to points (centroids), then calculate distances to every park
m_dist <- buildings_pop_sgp %>%
  sf::st_centroid() %>%
  sf::st_distance(parks_sgp) %>% # euclidean distance
    units::set_units(NULL)

m_dist <- m_dist / 1000 # convert distances to km


# new column for the supply of park area
buildings_pop_sgp$area_supply <- recre_supply(park_attribute = parks_sgp$area, 
                                              dist_matrix = m_dist, 
                                              c = 0.302) # e.g. from Tu et al. (2020)

# supply to all residents per building
buildings_pop_sgp$area_supplytopop <- buildings_pop_sgp$area_supply * buildings_pop_sgp$popcount

rm(m_dist)
```

```{r include = FALSE}
# calculate for year 2020
pop_2020 <- pop_sgp %>%
  dplyr::filter(year == 2020)
```

```{r summarise supply of park area per census unit, message = FALSE, warning = FALSE, include = FALSE}

# append subzone info to buildings if building intersects with subzone  
buildings_subzones <- buildings_pop_sgp  %>%
  sf::st_join(pop_2020,
              join = sf::st_intersects,
              left = TRUE) %>%
  sf::st_set_geometry(NULL)


# summarise total/median/average supply value for all buildings within subzone
buildings_subzones <- buildings_subzones %>%
  dplyr::group_by(subzone_n) %>%
  dplyr::summarise(across(.cols = c(area_supply, area_supplytopop), 
                          .fns = sum, .names = "{.col}_sum"),
                   across(.cols = c(area_supply, area_supplytopop), 
                          .fns = median, .names = "{.col}_median"),
                   across(.cols = c(area_supply, area_supplytopop), 
                          .fns = mean, .names = "{.col}_mean")) 
  

# join information to pop_2020, calculate supply per capita
pop_2020 <- pop_2020 %>%
  dplyr::left_join(buildings_subzones) %>%
  dplyr::mutate(area_supplyperpop = area_supplytopop_sum / pop_count * 1e-6) %>% # convert to km2
  dplyr::mutate(area_supplyperpop = ifelse(is.infinite(area_supplyperpop), NA, area_supplyperpop)) # make infinite NA

```

```{r include = FALSE}
sf::st_agr(pop_2020) = "constant"
sf::st_agr(parks_sgp) = "constant"
# https://github.com/r-spatial/sf/issues/406
```


```{r conventional park provision metric - per capita park area, message = FALSE, warning = FALSE, include = FALSE}

subzones_parks <- sf::st_intersection(pop_2020, parks_sgp) # subset census unit polygons to those that intersect parks, append park info
subzones_parks$parkarea_m2 <- sf::st_area(subzones_parks) # calc area of each polygon


# calculate total park area per census unit
subzones_parks <- subzones_parks %>%
  dplyr::group_by(subzone_n) %>%
  dplyr::summarise(parkarea_m2 = as.numeric(sum(parkarea_m2))) %>%
  sf::st_drop_geometry()


# join information to pop_2020, calculate park area per capita
pop_2020 <- pop_2020 %>%
  dplyr::left_join(subzones_parks) %>%
  dplyr::mutate(parkarea_m2 = ifelse(is.na(parkarea_m2), 0, parkarea_m2)) %>% # make NA 0
  dplyr::mutate(parkperpop_m2 = parkarea_m2 / pop_count * 1e-6) %>% # convert to km2
  dplyr::mutate(parkperpop_m2 = ifelse(is.infinite(parkperpop_m2), # make infinite NA
                                NA, parkperpop_m2)) %>%
  dplyr::mutate(parkperpop_m2 = ifelse(is.nan(parkperpop_m2), # make NaN NA
                                NA, parkperpop_m2))

```

```{r include = FALSE}
rm(buildings_subzones, subzones_parks)
```

<br>

The following interactive map provides a visual comparison between metrics for park area provision calculated in the '[Get started](https://ecological-cities.github.io/home2park/articles/home2park.html)' vignette. Toggle the map layers to view the per capita provision of park area summarised per region (subzone).

```{r  plot park metrics on map, echo = FALSE, message = FALSE, warning = FALSE, dpi = 300, fig.height = 2.0, fig.width = 2.35, fig.cap=paste0("**Map: Supply of park area in Singapore based on OSM data (2020).** Each building is denoted as a point (a random subset is shown). All color palettes are binned according to quantile values.")}

# convert buildings to points (centroids) for plotting
buildings_pop_sgp <- buildings_pop_sgp %>% 
  sf::st_centroid() %>%
  dplyr::mutate(across(everything(), as.vector)) # remove attributes from columns for plotting


# for area related vars, convert m2 to km2
buildings_pop_sgp <- buildings_pop_sgp %>% 
  dplyr::mutate(across(.cols = contains("area"), 
              .fns = function(x) x*1e-6))


# random sampling of buildings
# set.seed(123) 
# # some processing to get proper range for color scale later
# # get max values per building (for particular round)
# buildings_max <- buildings_pop_sgp %>%
#   tidyr::pivot_longer(cols = c("area_supply", "area_supplytopop"),
#                       names_to = "supply", 
#                       values_to = "value") %>%
#      dplyr::group_by(supply) %>%
#      dplyr::slice(which.max(value)) %>%
#   tidyr::pivot_wider(names_from = "supply", values_from = "value")
# 
# # subset random sample (already done with example dataset)
# buildings_pop_sgp <-
#   dplyr::slice_sample(buildings_pop_sgp,
#                       n = nrow(buildings_pop_sgp)/20)
# 
# # append max buildings to the random subset 
# buildings_pop_sgp <- buildings_pop_sgp %>%
#   dplyr::bind_rows(buildings_max)





tmap::tmap_mode("view")
tmap::tmap_options(check.and.fix = TRUE)

tm <- 
  tmap::tm_basemap(c("Esri.WorldGrayCanvas", "CartoDB.DarkMatter", "OpenStreetMap")) +

  # parks  
  tmap::tm_shape(parks_sgp %>% dplyr::select(id, name, area)) +
    tmap::tm_polygons(group = "Parks",
                      col = "#33a02c",
                      alpha = 0.6,
                      border.col = "grey50",  border.alpha = 0.5) +

  # buildings
  tmap::tm_shape(buildings_pop_sgp) +
    tmap::tm_dots(title = "Buildings: Supply of park area (km<sup>2</sup>)",
                  group = "Buildings: Supply of park area",
                  col = "area_supply", border.col = "transparent",
                  palette = viridis::viridis(5),
                  style = "quantile",
                  size = 0.01,
                  alpha = 0.7,
                  interactive = FALSE, # some glitch in values of hover text/popup
                  popup.vars = NULL,
                  showNA = FALSE) +
  
  tmap::tm_shape(buildings_pop_sgp) +
    tmap::tm_dots(title = "All building residents: Supply of park area (km<sup>2</sup>)",
                  group = "Buildings: Supply of park area to residents",
                  col = "area_supplytopop", border.col = "transparent",
                  palette = viridis::viridis(5),
                  style = "quantile",
                  size = 0.01,
                  alpha = 0.7,
                  interactive = FALSE, # some glitch in values of hover text/popup
                  popup.vars = NULL,
                  showNA = FALSE) +
    
  # census units
  tmap::tm_shape(pop_2020) +
    tmap::tm_polygons(title = "Supply of park area per capita (km<sup>2</sup>)",
                      group = "Subzones: Per capita supply of park area",
                      col = "area_supplyperpop", 
                      palette = "Greens", alpha = 0.7,
                      style = "quantile",
                      border.col = "white", border.alpha = 0.5, lwd = 1) +
  tmap::tm_shape(pop_2020) +
    tmap::tm_polygons(title = "Park area per capita (km<sup>2</sup>)",
                      group = "Subzones: Per capita park area (conventional)",
                      col = "parkperpop_m2",
                      palette = "Greens", alpha = 0.7,
                      style = "quantile",
                      border.col = "white", border.alpha = 0.5, lwd = 1)
  
  
# Pipe the tmap object into tmap_leaflet() to create a leaflet widget,
# so that we can use leaflet::hideGroup().
tm %>% 
  tmap::tmap_leaflet() %>%
  leaflet::hideGroup("Buildings: Supply of park area") %>%
  leaflet::hideGroup("Subzones: Per capita supply of park area") %>% 
  leaflet::hideGroup("Subzones: Per capita park area (conventional)") 

```


```{r include = FALSE}
rm(tm, buildings_pop_sgp, parks_sgp, pop_sgp, pop_2020)
```


<br>

## Data sources

- Singapore census data from the [Department of Statistics Singapore](https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore subzone polygons from the [Singapore Master Plan Subzones](https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea). Released under the terms of the [Singapore Open Data Licence version 1.0](https://data.gov.sg/open-data-licence).

- Singapore Master Plan Land Use Zones for the years [2014](https://data.gov.sg/dataset/master-plan-2014-land-use) and [2019](https://data.gov.sg/dataset/master-plan-2019-land-use-layer). Released under the terms of the [Singapore Open Data License](https://data.gov.sg/open-data-licence).

- Building polygons derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).

- Park polygons and summarised attributes (trails, playgrounds) derived from map data [copyrighted](https://www.openstreetmap.org/copyright) OpenStreetMap contributors and available from https://www.openstreetmap.org. Released under the terms of the [ODbL License](https://opendatacommons.org/licenses/odbl/summary/).


---
title: 'home2park: An R package to assess the spatial provision of urban parks'
tags:
  - R
  - Distance decay
  - Dasymetric map
  - Outdoor recreation
  - Park provision
  - Recreation supply
  - Service area
  - Spatial distribution
  - Urban green space
  - Urban parks
authors:
  - name: Xiao Ping Song
    orcid: 0000-0002-8825-195X
    affiliation: 1
  - name: Kwek Yan Chong
    orcid: 0000-0003-4754-8957
    affiliation: 1
affiliations:
 - name: Department of Biological Sciences, National University of Singapore
   index: 1
citation_author: Song and Chong
date: 11 Jul 2021
year: 2021
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

In the wake of the COVID-19 pandemic, significant changes to human mobility and working environments have prompted a re-think of land use distribution in the city, particularly for office, residential and public spaces such as parks. One prominent trend is a rising demand for the outdoors, as urban dwellers turn to local parks for recreation and leisure [@Venter2020]. However, current metrics of park provision tend to focus on summaries of park area within a given region [@Tan2017]. There is a need to characterise the wide variety of parks that serve different groups of people [@Song2020a], and to deal with the limitations of artificial boundaries when planning at fine spatial scales.

The package ``home2park`` provides a way to measure multiple aspects of park provision to homes, at the resolution of individual buildings. The key features include the ability to:  

1.	Download relevant data from OpenStreetMap (OSM) from any city worldwide (e.g. buildings, parks and features related to recreation), using the R package ``osmextract`` [@Lovelace2019]. The user may also supply data for new buildings and parks, for the purpose of future scenario planning. 
2.	Summarise at each park multiple attributes that are important for recreation (e.g. dense vegetation, length of waterfronts, open spaces, trails, etc.).
3.	Calculate the supply (provision) of the park attributes to residential buildings while accounting for ‘distance decay’ [@Rossi2015; @Tu2020b], or the fact that supply from parks further away are reduced.
4.	Redistribute coarse-scale population data per census unit into residential buildings, also known as ‘dasymetric mapping’ [@Dmowska], which helps highlight specific areas where more people will benefit from the presence of parks.

See the [package website](https://ecological-cities.github.io/home2park/) for detailed documentation and examples of how regional summaries derived from building-level metrics vary widely from conventional metrics of park area provision.

# Statement of Need

Numerous metrics have been used to measure the spatial provision of parks in cities. These include summaries of the per capita park area (i.e. park provision ratio) or ‘service radius’ (a distance buffer) around parks within geographical areal units such as census blocks or administrative boundaries [@Tan2017]. Such metrics have been used to investigate park use and accessibility, as well as social-cultural issues related to spatial equity and environmental justice [@Sister2010; @Tan2017]. However, the presence of artificial boundaries gives rise to the ‘modifiable areal unit problem’, where differences in spatial scale can result in considerable variation in the value of summarised metrics [@Sister2010]. Since urban development typically occurs at fine spatial scales, metrics at the scale of individual buildings provided by ``home2park`` can help overcome such limitations [@Gao2017]. Finally, most metrics of park provision do not differentiate between the wide variety of parks that would likely affect usage by different groups of people. Other than park area, features such as waterfronts, open lawns, natural vegetation, trails and fitness amenities each relate to different types of park use and user groups [@Song2020a; @Song2020b]. A more holistic measure of park provision should thus include a variety of park attributes that are important for recreation and leisure. ``home2park`` enables the user to summarise a variety of spatial attributes such as points (e.g. playgrounds, sports and fitness amenities), lines (e.g. cycling and walking trails) and rasters (e.g. land cover types) at each park. Customisable parameters provide the user with flexibility when summarising specific attributes, such as specifying a minimum patch size for a land cover class, or including spatial features within a certain distance beyond park boundaries. By providing fine-grained measures for a variety of park attributes, we hope that ``home2park`` will contribute toward a more nuanced understanding of recreation that will improve the planning and design of parks and homes across the city.

## State of the Field

To our knowledge, there are currently no software (e.g. R packages) that measure the spatial provision of parks and recreation. 


# Acknowledgements

We thank Edwin Y. W. Tan and Justin K. H. Nai for useful discussions on the methodology. 

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_dasymap.R
\name{pop_dasymap}
\alias{pop_dasymap}
\title{Perform dasymetric mapping of (rasterised) population data}
\usage{
pop_dasymap(
  pop_polygons,
  pop_perblock_count,
  pop_perblock_density,
  land_relative_density,
  filename = NULL,
  overwrite = TRUE,
  ...
)
}
\arguments{
\item{pop_polygons}{\code{sf} polygons of the population data for a single time period.
May be extracted from the \code{pop_polygons} list element generated by \code{rasterise_pop()}.}

\item{pop_perblock_count}{Raster of population count per census block for a single time period.
May be extracted from the \code{pop_count_rasters} list element generated by \code{rasterise_pop()}.}

\item{pop_perblock_density}{Raster of population density per census block for a single time period.
May be extracted from the \code{pop_density_rasters} list element generated by \code{rasterise_pop()}.}

\item{land_relative_density}{Raster of the relative density (e.g. suitability, habitability) of land for a single time period.
May be extracted from results generated by \code{rasterise_buildings()} or \code{rasterise_landuse()}.}

\item{filename}{character (optional). Export output raster to disk.}

\item{overwrite}{logical. Argument passed to \code{terra::writeRaster()}. If \code{TRUE}, \code{filename} is overwritten.}

\item{...}{Other arguments passed to \code{terra::writeRaster()}.}
}
\value{
Raster of population density, based on relative density values of the land defined in \code{land_relative_density}.
}
\description{
Distribute population data per census block across inhabitable land as defined by \code{land_relative_density}.
Raster data provided should be a snapshot of the population for a single time period (e.g. a year),
created with other package functions (e.g. \code{rasterise_pop()}, \code{rasterise_buildings()}).
If there are multiple years present in the intermediate datasets, the user will have to extract the relevant list element for
the single year of interest (see examples). All input data should have a similar projected coordinate reference system
specific to the target area.
}
\examples{
\dontrun{
data(pop_sgp) # population census block polygons
data(landuse_sgp) # land use polygons


# transform to projected crs
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648))
landuse_sgp <- sf::st_transform(landuse_sgp, sf::st_crs(32648))


# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

buildings <- get_buildings_osm(place = city_boundaries,
                               date = as.Date('2021-01-01')) \%>\%
   mutate(year = 2020)


# rasterise population, landuse & buildings
pop_rasters <- rasterise_pop(pop_sgp,
                             census_block = "subzone_n",
                             pop_count = "pop_count")

landuse_rasters <- rasterise_landuse(landuse_sgp,
                                     land_use = 'lu_desc',
                                     subset = c('1' = 'RESIDENTIAL',
                                                '2' = 'COMMERCIAL & RESIDENTIAL',
                                                '3' = 'RESIDENTIAL WITH COMMERCIAL AT 1ST STOREY',
                                                '4' = 'RESIDENTIAL / INSTITUTION'),
                                     sf_pop = pop_sgp,
                                     match_landuse_pop = 'recent')

buildings_rasters <- rasterise_buildings(buildings,
                                         proxy_pop_density = 'levels',
                                         year = 'year',
                                         sf_pop = pop_sgp,
                                         sf_landuse = landuse_sgp,
                                         match_buildings_pop = 'closest')


# finally, perform dasymetric mapping on selected year (2020)
popdens_raster <- pop_dasymap(pop_polygons = pop_rasters$pop_polygons[[2]],
                              pop_perblock_count = pop_rasters$pop_count[[2]],
                              pop_perblock_density = pop_rasters$pop_density[[2]],
                              land_relative_density = buildings_rasters[[2]],
                              filename = 'buildings_popdensity.tif',
                              wopt = list(gdal=c('COMPRESS=LZW')))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchyear.R
\name{matchyear}
\alias{matchyear}
\title{Match the year between the dataset of interest and another data source}
\usage{
matchyear(
  data,
  data_tomatch,
  match = c("exact", "closest", "recent", "soonest"),
  year = NULL
)
}
\arguments{
\item{data}{Tabular object (e.g. \code{data.frame}, \code{sf}, \code{tibble}) containing the data.}

\item{data_tomatch}{Tabular object (e.g. \code{data.frame}, \code{sf}, \code{tibble}) to match the data to.}

\item{match}{Type of matching; either \code{'exact'}, \code{'closest'}, \code{'recent'} or \code{'soonest'}.}

\item{year}{Specify column name for the year within both datasets.
Columns in both datasets should be numeric.}
}
\value{
\code{data} with additional column \code{year_match},
representing the matching year of the other dataset.
}
\description{
Helper function to match the year in the data with either the (1) exact;
(2) most recent (past); (3) closest (past or future); or (4) soonest (future) year in another dataset.
The column name within both datasets that represent the year must be specified.
Dataset containing multiple years must be in the 'long' format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-singapore.R
\docType{data}
\name{parks_sgp}
\alias{parks_sgp}
\title{Public parks in Singapore}
\format{
\code{sf} polygons.
}
\source{
Map data \href{https://www.openstreetmap.org/copyright}{copyrighted} OpenStreetMap contributors and available from https://www.openstreetmap.org.
Made available via the \href{https://opendatacommons.org/licenses/odbl/summary/}{ODbL License}.
}
\usage{
parks_sgp
}
\description{
Example dataset of park polygons in Singapore downloaded from OpenStreetMap
(data snapshot on \code{2021-01-01} from the \href{https://download.geofabrik.de}{Geofabrik database}),
using the function \code{get_parks_osm()}. Includes summaries (columns) of selected attributes related
to outdoor recreation, i.e., playgrounds (see \code{get_playgrounds_osm()}) and trails (see \code{get_trails_osm()})
calculated using the function \code{parks_calc_attributes()}.
}
\examples{
data(parks_sgp)
head(parks_sgp)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-singapore.R
\docType{data}
\name{pop_sgp}
\alias{pop_sgp}
\title{Population counts within census block polygons in Singapore}
\format{
\code{sf} polygons with census data for years 2018 and 2020.
Data is in the 'long' format (rows are repeated for each census year).
The following columns are used in this package:
\describe{
\item{subzone_n}{Census block name}
\item{year}{Census year}
\item{pop_count}{Population count within census block}
}
}
\source{
Contains census data from the \href{https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data}{Department of Statistics Singapore}
and polygons from the \href{https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea}{Singapore Master Plan Subzones},
both of which are made available under the terms of the \href{https://data.gov.sg/open-data-licence}{Singapore Open Data Licence version 1.0}.
}
\usage{
pop_sgp
}
\description{
Example dataset containing Singapore census data for years 2018 and 2020 from the \href{https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data}{Department of Statistics Singapore},
joined by name to the Master Plan Subzone polygons from \href{https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea}{data.gov.sg}.
Both datasets are released under the \href{https://data.gov.sg/open-data-licence}{Singapore Open Data License}.
}
\examples{
data(pop_sgp)
head(pop_sgp)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rasterise_buildings.R
\name{rasterise_buildings}
\alias{rasterise_buildings}
\title{Rasterise building polygons}
\usage{
rasterise_buildings(
  sf_buildings,
  proxy_pop_density = NULL,
  year = NULL,
  sf_pop = NULL,
  sf_landuse = NULL,
  match_buildings_pop = "closest",
  match_buildings_landuse = "closest",
  raster_template = NULL,
  dir_processing = tempdir(),
  dir_export = NULL,
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW")),
  ...
)
}
\arguments{
\item{sf_buildings}{\code{sf }polygons of the buildings.
Data should be in a projected coordinate reference system.}

\item{proxy_pop_density}{character. Specify column name of the building attribute that is used
as an indicator of population density (e.g. building height, no. of levels).
This column will be used to assign values to the output raster. If not provided, population density
across all building pixels are assumed to be similar.}

\item{year}{character. Specify column name for the year within \code{sf_buildings}, \code{sf_landuse}
and \code{sf_pop} (if provided). Defaults to \code{NULL}, assuming that there is no 'year' column. Column data should be numeric.}

\item{sf_pop}{(optional) \code{sf} polygons containing the population census data with column containing the census year.
If absent, the output (building rasters) will not be associated with specific population census year(s).}

\item{sf_landuse}{(optional) \code{sf} polygons of the land use zones.}

\item{match_buildings_pop}{character. Type of matching between \code{sf_buildings} and \code{sf_pop} (if provided), passed on to \code{match} argument in function \code{matchyear()}.
Either \code{'exact'}, \code{'closest'}, \code{'recent'} or \code{'soonest'}. Defaults to \code{'closest'}.}

\item{match_buildings_landuse}{character. Type of matching between \code{sf_buildings} and \code{sf_landuse} (if provided), passed on to \code{match} argument in function \code{matchyear()}.
Either \code{'exact'}, \code{'closest'}, \code{'recent'} or \code{'soonest'}. Defaults to \code{'closest'}.
This argument is used only if \code{sf_pop} is not provided.}

\item{raster_template}{Either a \code{terra::rast} object, or a filepath to the raster used to define the pixel resolution, extent, nrow, ncol of
the output raster. Defaults to raster template in \code{tempdir()} processed in previous steps.
Object is passed to the '\code{y}' argument in \code{terra::rasterize()}.
Defaults to the template raster generated by the function \code{rasterise_pop()} within \code{dir_processing} (see next argument).}

\item{dir_processing}{character. Directory to get intermediate files generated in previous steps (e.g. rasterised population/land use data).
Defaults to \code{tempdir()}.}

\item{dir_export}{character (optional). File path to directory to export output raster(s) to.}

\item{overwrite}{logical. Argument passed to \code{terra::writeRaster()}. Defaults to \code{TRUE}.}

\item{wopt}{list. Argument passed to \code{terra::writeRaster()}.}

\item{...}{Other arguments passed to \code{terra::writeRaster()}.}
}
\value{
List of raster files. Zero values are converted to \code{NA}.
}
\description{
Rasterise building polygons (\code{sf_buildings}) with reference to population (\code{sf_pop}) and/or
land use data (\code{sf_landuse}), if supplied. Multiple output (a list of) rasters will be generated if
building data for multiple years are present.
}
\details{
If population (\code{sf_pop}) and/or land use data (\code{sf_landuse}) are supplied,
their raster(s) must have been previously generated in the \code{tempdir()} using
the functions \code{rasterise_pop()} and \code{rasterise_landuse()}, respectively.
Rasterised buildings will be masked away (i.e., convert pixels to \code{NA})
at areas with no (zero) population and/or land use data, according to the respective matching year(s).
Custom ways to match the years between datasets can be set via the helper function \code{matchyear()}.
The argument \code{match_buildings_pop} provides ways to match each year in the \code{sf_pop} (if provided)
to a specific year in \code{sf_buildings}. Building raster(s) will then be associated
with a specific population census year, and removed if there are no matching population census years.

If necessary, the argument \code{match_buildings_landuse} provides ways to match each year in \code{sf_buildings}
with a specific year in \code{sf_landuse}, but only if \code{sf_pop} is not provided. This is because the main
reference point for all the datasets is \code{sf_pop} (primary focus is the analysis of population data).
}
\examples{
\dontrun{
# load data
data(pop_sgp) # population census block polygons
data(landuse_sgp) # land use polygons


# transform to projected crs
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648))
landuse_sgp <- sf::st_transform(landuse_sgp, sf::st_crs(32648))


# get osm buildings based on census block polygons (year 2020)
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

buildings <- get_buildings_osm(place = city_boundaries,
                               date = as.Date('2021-01-01')) \%>\%
   mutate(year = 2020)


# run function
buildings_rasters <- rasterise_buildings(buildings,
                                         proxy_pop_density = 'levels',
                                         year = 'year',
                                         sf_pop = pop_sgp,
                                         sf_landuse = landuse_sgp,
                                         match_buildings_pop = 'closest')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_parks_osm.R
\name{get_beaches_osm}
\alias{get_beaches_osm}
\title{Get public beaches from OpenStreetMap}
\usage{
get_beaches_osm(
  place,
  date = NULL,
  mutually_exclusive_with = list(),
  snap_tolerance = 5,
  min_area = units::set_units(0, "m^2"),
  aggregate_polygons = 15,
  dir_raw = osmextract::oe_download_directory(),
  filename = NULL,
  ...
)
}
\arguments{
\item{place}{\code{sf} object (with projected coordinate reference system). Geographical area to match with the (\code{.osm.pbf}) file in the data archive.
Argument passed to \code{osmextract::oe_match()}.}

\item{date}{Date of OSM data snapshot to download. Object of class "Date" in format \verb{\%Y-\%m-\%d}. Refer to https://download.geofabrik.de
for the specific dates available. Defaults to \code{NULL} (download the latest available data).}

\item{mutually_exclusive_with}{list of \code{sf} object(s). This may be used to
ensure that polygons (e.g. parks, beaches, informal nature areas) are mutually-exclusive (i.e. non-overlapping).
Remove output polygons \href{https://postgis.net/docs/ST_Contains.html}{contained} within, as well as \href{https://postgis.net/docs/ST_Intersection.html}{intersections} between,
each element of this list. Should have the same coordinate reference system as \code{place}.}

\item{snap_tolerance}{numeric. Argument for \code{tolerance} level passed to \code{sf::st_snap()},
used to rectify nearly coincident edges between polygons before processing (e.g. \code{sf::st_contains()}, \code{sf::st_covers()}).
Provided either as a units object (see \code{units::set_units()}), or a number in the units of the coordinate reference system.
Defaults to \code{5}. Set to \code{0} if you do not wish to rectify minor overlaps.}

\item{min_area}{numeric. Specify minimum area of each polygon to be retained in the output,
passed to argument \code{threshold} in \code{smoothr::drop_crumbs()}.
Provided either as a units object (see \code{units::set_units()}), or a number in the units of
the coordinate reference system. Defaults to \code{0} m^2.}

\item{aggregate_polygons}{numeric. Argument for \code{dist} passed to \code{sf::st_buffer()}.
Buffered polygons that overlap will be aggregated into multipolygons.
Set to \code{NULL} if you do not wish to aggregate to multipolygons.}

\item{dir_raw}{character. Directory to download the raw unprocessed OSM data. Passed to
argument \code{download_directory} in \code{osmextract::oe_read()}.}

\item{filename}{character (optional). File path to export output data (GeoJSON format).}

\item{...}{Other arguments passed to \code{osmextract::oe_read()}.}
}
\value{
The processed beach polygons (\code{sf} object).
}
\description{
Download and process OpenStreetMap (OSM) beach polygons within a specified geographical \code{place},
from the \href{https://download.geofabrik.de}{Geofabrik database}. It is a wrapper around
functions in the package \href{https://docs.ropensci.org/osmextract/index.html}{\code{osmextract}}, and
processes the downloaded files for subsequent analyses. Refer to package \code{osmextract} for
more details and options for input arguments when downloading the data.
}
\details{
OSM polygons are filtered by key-value attributes, where \code{natural:beach}, and \verb{access:} is not \code{no} or \code{private}.
If \code{mutually_exclusive_with} is provided, \href{https://postgis.net/docs/ST_Intersection.html}{intersections} between
the output and these polygon(s) will be excluded using \code{polygons_mutually_exclude()}.
Polygons are then cleaned up using \code{polygons_clean()}.
}
\examples{
\dontrun{
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs

# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

data(parks_sgp) # to exclude beaches within/intersecting these polygons
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648)) # transform to projected crs

# run function
get_beaches_osm(place = city_boundaries,
                date = as.Date('2021-01-01'),
                mutually_exclusive_with = list(parks_sgp),
                snap_tolerance = 5,
                aggregate_polygons = 15,
                filename = 'public-beaches_osm-polygons_2021-01-01.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rasterise_pop.R
\name{rasterise_pop}
\alias{rasterise_pop}
\title{Rasterise population counts within census block polygons}
\usage{
rasterise_pop(
  sf,
  res = 10,
  census_block = NULL,
  pop_count = NULL,
  year = "year",
  dir_processing = tempdir()
)
}
\arguments{
\item{sf}{\code{sf}polygons containing the population census data.
Data should be in a projected coordinate reference system (for calculation of pixel resolution).}

\item{res}{number. Specify the pixel resolution for rasterisation (in metres). Defaults to \code{10}.}

\item{census_block}{character. Specify column name of the census block unique identifier (e.g. name).}

\item{pop_count}{character. Specify column name of the population counts per census block.}

\item{year}{character. Specify column name of the census year. Defaults to \code{'year'}.}

\item{dir_processing}{character. Directory to store intermediate files. Defaults to \code{tempdir()}.
Set to \code{NULL} if you do not wish to export intermediate (temporary) files for subsequent processing.}
}
\value{
List containing (1) rasterised population count per census block,
(2) rasterised population density per census block, and (3) processed polygons used to generate
these rasters. Named \code{pop_count_rasters}, \code{pop_density_rasters} and \code{pop_polygons}, respectively.
List is nested if multiple years are present (one sub-list for each census year).
Zero values for population counts/density are converted to \code{NA}. Intermediate (raster) files are exported to \code{tempdir()} for further processing.
}
\description{
Convert population counts per census block (polygons)
into population density grid (raster) for subsequent analyses.
Each polygon's count is divided by its area, such that the integrated
density over each census block should equal the original count.
Census data for multiple years can be processed to output a list of rasters (one raster per year).
The column name for the \code{'year'} in \code{sf} must be specified, even if data does not contain multiple years
(combined data across multiple years must be in the 'long' format).
}
\examples{
\dontrun{
data(pop_sgp)


# transform to projected crs
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648))


# run function
pop_rasters <- rasterise_pop(pop_sgp,
                             res = 10,
                             census_block = "subzone_n",
                             pop_count = "pop_count",
                             year = "year")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parks_calc_attributes_helpers.R
\name{raster_class_area}
\alias{raster_class_area}
\title{Calculate the total class area of a (classified) raster for each park polygon.}
\usage{
raster_class_area(
  polygons,
  raster,
  raster_min_patch_size = units::set_units(0, "m^2"),
  relative = TRUE,
  ...
)
}
\arguments{
\item{polygons}{\code{sf} (with projected coordinate reference system).}

\item{raster}{\code{SpatRaster} object from \code{terra::rast()}. Should have a (projected) coordinate reference system similar to \code{polygons}.}

\item{raster_min_patch_size}{Minimum patch size to be included in results.
Provided either as a units object (see \code{units::set_units()}),
or a number in the units of the coordinate reference system. Defaults to \code{0} m^2.}

\item{relative}{logical. Whether or not to calculate relative amounts
(i.e. percentage area). Defaults to \code{TRUE}.}

\item{...}{Other arguments passed to \code{terra::extract()}.}
}
\value{
\code{polygons} with added column(s) \verb{< class value >_area}, and \verb{< class value >_area_pct} if \code{relative} is set to \code{TRUE}.
Note that the value \code{0} will be summarised; convert pixels to \code{NA} if you wish to exclude them.
}
\description{
Helper function within \code{parks_calc_attributes()}.
For each unique pixel value (each representing a specific raster class/category),
the total area per class is calculated.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_amenities_osm.R
\name{get_sportfitness_osm}
\alias{get_sportfitness_osm}
\title{Get public sports/fitness amenities from OpenStreetMap}
\usage{
get_sportfitness_osm(
  place,
  date = NULL,
  dir_raw = osmextract::oe_download_directory(),
  filename = NULL,
  ...
)
}
\arguments{
\item{place}{\code{sf} object (with projected coordinate reference system). Geographical area to match with the (\code{.osm.pbf}) file in the data archive.
Argument passed to \code{osmextract::oe_match()}.}

\item{date}{Date of OSM data snapshot to download. Object of class "Date" in format \verb{\%Y-\%m-\%d}. Refer to https://download.geofabrik.de
for the specific dates available. Defaults to \code{NULL} (download the latest available data).}

\item{dir_raw}{character. Directory to download the raw unprocessed OSM data. Passed to
argument \code{download_directory} in \code{osmextract::oe_read()}.}

\item{filename}{character (optional). File path to export output data (GeoJSON format).}

\item{...}{Other arguments passed to \code{osmextract::oe_read()}.}
}
\value{
The processed sport/fitness amenities (\code{sf} object).
}
\description{
Download and process OpenStreetMap (OSM) public sports/fitness amenities (points) within \code{place},
from the \href{https://download.geofabrik.de}{Geofabrik database}. It is a wrapper around
functions in the package \href{https://docs.ropensci.org/osmextract/index.html}{\code{osmextract}}, and
processes the downloaded files for subsequent analyses. Refer to package \code{osmextract} for
more details and options for input arguments when downloading the data.
}
\details{
OSM points filtered by key-value attributes, where \code{leisure:fitness_station} or \verb{sport:*},
and \verb{access:} is not \code{no} or \code{private}.
}
\examples{
\dontrun{
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs

# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

# run function
get_sportfitness_osm(place = city_boundaries,
                    date = as.Date('2021-01-01'),
                    filename = 'sport-fitness_osm-points_2021-01-01.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-singapore.R
\docType{data}
\name{buildings_pop_sgp}
\alias{buildings_pop_sgp}
\title{Population count per residential building in Singapore}
\format{
\code{sf} polygons.
}
\source{
Building polygons \href{https://www.openstreetmap.org/copyright}{copyrighted} OpenStreetMap contributors and available from https://www.openstreetmap.org.
Made available via the \href{https://opendatacommons.org/licenses/odbl/summary/}{ODbL License}.

Population data and census block polygons \code{data(pop_sgp)} are from the \href{https://www.singstat.gov.sg/find-data/search-by-theme/population/geographic-distribution/latest-data}{Department of Statistics Singapore};
and \href{https://data.gov.sg/dataset/master-plan-2019-subzone-boundary-no-sea}{Singapore Master Plan Subzones}, respectively;
land use polygons \code{data(landuse_sgp)} are from the Singapore Land Use Master Plan released in \href{https://data.gov.sg/dataset/master-plan-2019-land-use-layer}{2019}.
All are made available under the terms of the \href{https://data.gov.sg/open-data-licence}{Singapore Open Data Licence version 1.0}.
}
\usage{
buildings_pop_sgp
}
\description{
Example (random 5\% subset) dataset of residential building polygons in Singapore,
each with a population count (column \code{popcount})
estimated via dasymetric mapping.
}
\details{
Building polygons were downloaded from OpenStreetMap
(data snapshot on \code{2021-01-01} from the \href{https://download.geofabrik.de}{Geofabrik database}),
using the function \code{get_buildings_osm()}. The population count per census block in the year 2020
was re-distributed across the buildings located within residential land use zones,
by performing dasymetric mapping using the functions \code{pop_dasymap()} and \code{pop_density_polygonise()}.
See vignette and examples in \code{pop_density_polygonise()} for more details.
The dataset is a random 5\% subset of the resulting polygons.
}
\examples{
data(buildings_pop_sgp)
head(buildings_pop_sgp)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_parks_osm.R
\name{get_informalnature_osm}
\alias{get_informalnature_osm}
\title{Get informal nature areas from OpenStreetMap}
\usage{
get_informalnature_osm(
  place,
  date = NULL,
  mutually_exclusive_with = list(),
  snap_tolerance = 5,
  min_area = units::set_units(0, "m^2"),
  min_trails = units::set_units(0, "m"),
  aggregate_polygons = 15,
  dir_raw = osmextract::oe_download_directory(),
  filename = NULL,
  ...
)
}
\arguments{
\item{place}{\code{sf} object (with projected coordinate reference system). Geographical area to match with the (\code{.osm.pbf}) file in the data archive.
Argument passed to \code{osmextract::oe_match()}.}

\item{date}{Date of OSM data snapshot to download. Object of class "Date" in format \verb{\%Y-\%m-\%d}. Refer to https://download.geofabrik.de
for the specific dates available. Defaults to \code{NULL} (download the latest available data).}

\item{mutually_exclusive_with}{list of \code{sf} object(s). This may be used to
ensure that polygons (e.g. parks, beaches, informal nature areas) are mutually-exclusive (i.e. non-overlapping).
Remove output polygons \href{https://postgis.net/docs/ST_Contains.html}{contained} within, as well as \href{https://postgis.net/docs/ST_Intersection.html}{intersections} between,
each element of this list. Should have the same coordinate reference system as \code{place}.}

\item{snap_tolerance}{numeric. Argument for \code{tolerance} level passed to \code{sf::st_snap()},
used to rectify nearly coincident edges between polygons before processing (e.g. \code{sf::st_contains()}, \code{sf::st_covers()}).
Provided either as a units object (see \code{units::set_units()}), or a number in the units of the coordinate reference system.
Defaults to \code{5}. Set to \code{0} if you do not wish to rectify minor overlaps.}

\item{min_area}{numeric. Specify minimum area of each polygon to be retained in the output,
passed to argument \code{threshold} in \code{smoothr::drop_crumbs()}.
Provided either as a units object (see \code{units::set_units()}), or a number in the units of
the coordinate reference system. Defaults to \code{0} m^2.}

\item{min_trails}{numeric. Specify minimum length of OSM trail lines that has to be within
nature area polygons, for the polygons to be retained in the output.
Provided either as a units object (see \code{units::set_units()}), or a number in the units of
the coordinate reference system. Defaults to \code{0} m.}

\item{aggregate_polygons}{numeric. Argument for \code{dist} passed to \code{sf::st_buffer()}.
Buffered polygons that overlap will be aggregated into multipolygons.
Set to \code{NULL} if you do not wish to aggregate to multipolygons.}

\item{dir_raw}{character. Directory to download the raw unprocessed OSM data. Passed to
argument \code{download_directory} in \code{osmextract::oe_read()}.}

\item{filename}{character (optional). File path to export output data (GeoJSON format).}

\item{...}{Other arguments passed to \code{osmextract::oe_read()}.}
}
\value{
The processed 'informal nature areas' (polygons) that are publicly accessible (\code{sf} object).
}
\description{
Download and process OpenStreetMap (OSM) 'informal nature' polygons within a specified geographical \code{place},
from the \href{https://download.geofabrik.de}{Geofabrik database}. It is a wrapper around
functions in the package \href{https://docs.ropensci.org/osmextract/index.html}{\code{osmextract}}, and
processes the downloaded files for subsequent analyses. Refer to package \code{osmextract} for
more details and options for input arguments when downloading the data.
}
\details{
OSM polygons are filtered by key-value attributes, where \verb{landuse:} is \code{forest} or \code{meadow},
or \verb{natural:} is \code{wood}, \code{scrub}, \code{heath}, \code{grassland}, \code{wetland}, \code{marsh}, \code{fell} or \code{tundra},
and \verb{access:} is not \code{no}, \code{private} or \code{restricted}.
Polygons \href{https://postgis.net/docs/ST_Contains.html}{contained} within \code{leisure:golf_course} are excluded after mutually-\href{https://postgis.net/docs/ST_Snap.html}{snapping} vertices between the two
(see argument \code{snap_tolerance}). \href{https://postgis.net/docs/ST_Intersection.html}{Intersections} with \code{landuse:military} polygons are also removed.
To exclude nature areas that are inaccessible, OSM trail lines are extracted (\verb{highway:} is \code{track}, \code{path}, \code{footway} or \code{cycleway},
and \verb{access:} is not \code{no} or \code{private}). The remaining polygons without such trails \href{https://postgis.net/docs/ST_Within.html}{within} their boundaries are removed.
If \code{mutually_exclusive_with} is provided, \href{https://postgis.net/docs/ST_Intersection.html}{intersections} between
the output and these polygon(s) will be excluded using \code{polygons_mutually_exclude()}.
Polygons are then cleaned up using \code{polygons_clean()}.
}
\examples{
\dontrun{
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs

# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

data(parks_sgp) # to exclude nature areas within/intersecting these polygons
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648)) # transform to projected crs

# run function
get_informalnature_osm(place = city_boundaries,
                       date = as.Date('2021-01-01'),
                       mutually_exclusive_with = list(parks_sgp),
                       snap_tolerance = 5,
                       min_area = units::set_units(2500, 'm^2'),
                       min_trails = units::set_units(0, 'm'),
                       aggregate_polygons = 15,
                       filename = 'nature-areas_osm-polygons_2021-01-01.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rasterise_landuse.R
\name{rasterise_landuse}
\alias{rasterise_landuse}
\title{Rasterise residential land use polygons}
\usage{
rasterise_landuse(
  sf_landuse,
  subset = NULL,
  land_use = NULL,
  year = "year",
  sf_pop = NULL,
  match_landuse_pop = "recent",
  dir_rastertemplate = NULL,
  dir_processing = tempdir(),
  dir_export = tempdir(),
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW")),
  ...
)
}
\arguments{
\item{sf_landuse}{\code{sf} polygons of the land use zones.
Data should be in a projected coordinate reference system.}

\item{subset}{Named vector containing the residential land use classes of interest to retain in column \code{land_use}.
Each element is named according to the integer that will be used
to represent this class in the output raster. Defaults to \code{NULL} (integers will be alphabetically
assigned to land use classes).}

\item{land_use}{character. Specify column name of the land use classes.}

\item{year}{character. Specify column name for the year within \code{sf_landuse}
and \code{sf_pop} (if provided). Defaults to \code{'year'}. Column data should be numeric.}

\item{sf_pop}{(optional) \code{sf} polygons containing the population census data with column containing the census year.
If absent, the output (land use rasters) will not be associated with specific population census year(s).}

\item{match_landuse_pop}{character. Type of matching between \code{sf_landuse} and \code{sf_pop}, passed on to \code{match} argument in function \code{matchyear()}.
Either \code{'exact'}, \code{'closest'}, \code{'recent'} or \code{'soonest'}. Defaults to \code{'recent'}, i.e.
assumes that land use is pre-planned (and followed by population changes), rather than a result of post-monitoring (e.g. remotely-sensed).}

\item{dir_rastertemplate}{character. Filepath to the raster used to define the pixel resolution, extent, nrow, ncol of
the output  raster; object is passed to the '\code{y}' argument in \code{terra::rasterize()}.
Defaults to the template raster generated by the function \code{rasterise_pop()} within \code{dir_processing} (see next argument).}

\item{dir_processing}{character. Directory to get intermediate files generated in previous steps (e.g. rasterised population data).
Defaults to \code{tempdir()}.}

\item{dir_export}{character. File path to directory to export output raster(s) to. Defaults to \code{tempdir()}.
Set to \code{NULL} if you do not wish to export output for subsequent processing.}

\item{overwrite}{logical. Argument passed to \code{terra::writeRaster()}. Defaults to \code{TRUE}.}

\item{wopt}{list. Argument passed to \code{terra::writeRaster()}.}

\item{...}{Other arguments passed to \code{terra::writeRaster()}.}
}
\value{
List of raster files. Zero values are converted to \code{NA}.
Each raster is also an intermediate file exported to the directory as defined by \code{dir_export}
(defaults to \code{tempdir()}) for further processing.
}
\description{
Rasterise land use zones (\code{sf_landuse}) with reference to population data (\code{sf_pop}), if supplied.
Multiple (a list of) output rasters will be generated if land use data for multiple years are present.
If the land use data includes unnecessary land use classes, the argument \code{subset} and \code{landuse}
allows the user to first subset the data to polygons defined as 'residential' land use.
}
\details{
If population data (\code{sf_pop}) is supplied, its raster(s) must have been previously generated in the \code{tempdir()} using
the function \code{rasterise_pop()}. Each year in the population data is matched to a specific year in the land use data
(using the helper function \code{matchyear()}). Land use raster(s) will then
be associated with a specific population census year, and removed if there are no matching population census years.
Rasterised land use zones will be masked away (i.e., convert pixels to \code{NA})
at areas with no (zero) population data (for the respective matching year).
}
\examples{
\dontrun{
data(pop_sgp)
data(landuse_sgp)


# transform to projected crs
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648))
landuse_sgp <- sf::st_transform(landuse_sgp, sf::st_crs(32648))


# run function
landuse_rasters <- rasterise_landuse(landuse_sgp,
                                     land_use = 'lu_desc',
                                     subset = c('1' = 'RESIDENTIAL',
                                                '2' = 'COMMERCIAL & RESIDENTIAL',
                                                '3' = 'RESIDENTIAL WITH COMMERCIAL AT 1ST STOREY',
                                                '4' = 'RESIDENTIAL / INSTITUTION'),
                                     year = 'year',
                                     sf_pop = pop_sgp,
                                     match_landuse_pop = 'recent')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_parks_osm.R
\name{get_parks_osm}
\alias{get_parks_osm}
\title{Get public parks from OpenStreetMap}
\usage{
get_parks_osm(
  place,
  date = NULL,
  mutually_exclusive_with = list(),
  snap_tolerance = 5,
  min_area = units::set_units(0, "m^2"),
  aggregate_polygons = 15,
  dir_raw = osmextract::oe_download_directory(),
  filename = NULL,
  ...
)
}
\arguments{
\item{place}{\code{sf} object (with projected coordinate reference system). Geographical area to match with the (\code{.osm.pbf}) file in the data archive.
Argument passed to \code{osmextract::oe_match()}.}

\item{date}{Date of OSM data snapshot to download. Object of class "Date" in format \verb{\%Y-\%m-\%d}. Refer to https://download.geofabrik.de
for the specific dates available. Defaults to \code{NULL} (download the latest available data).}

\item{mutually_exclusive_with}{list of \code{sf} object(s). This may be used to
ensure that polygons (e.g. parks, beaches, informal nature areas) are mutually-exclusive (i.e. non-overlapping).
Remove output polygons \href{https://postgis.net/docs/ST_Contains.html}{contained} within, as well as \href{https://postgis.net/docs/ST_Intersection.html}{intersections} between,
each element of this list. Should have the same coordinate reference system as \code{place}.}

\item{snap_tolerance}{numeric. Argument for \code{tolerance} level passed to \code{sf::st_snap()},
used to rectify nearly coincident edges between polygons before processing (e.g. \code{sf::st_contains()}, \code{sf::st_covers()}).
Provided either as a units object (see \code{units::set_units()}), or a number in the units of the coordinate reference system.
Defaults to \code{5}. Set to \code{0} if you do not wish to rectify minor overlaps.}

\item{min_area}{numeric. Specify minimum area of each polygon to be retained in the output,
passed to argument \code{threshold} in \code{smoothr::drop_crumbs()}.
Provided either as a units object (see \code{units::set_units()}), or a number in the units of
the coordinate reference system. Defaults to \code{0} m^2.}

\item{aggregate_polygons}{numeric. Argument for \code{dist} passed to \code{sf::st_buffer()}.
Buffered polygons that overlap will be aggregated into multipolygons.
Set to \code{NULL} if you do not wish to aggregate to multipolygons.}

\item{dir_raw}{character. Directory to download the raw unprocessed OSM data. Passed to
argument \code{download_directory} in \code{osmextract::oe_read()}.}

\item{filename}{character (optional). File path to export output data (GeoJSON format).}

\item{...}{Other arguments passed to \code{osmextract::oe_read()}.}
}
\value{
The processed park polygons (\code{sf} object).
}
\description{
Download and process OpenStreetMap (OSM) park polygons within a specified geographical \code{place},
from the \href{https://download.geofabrik.de}{Geofabrik database}. It is a wrapper around
functions in the package \href{https://docs.ropensci.org/osmextract/index.html}{\code{osmextract}}, and
processes the downloaded files for subsequent analyses. Refer to package \code{osmextract} for
more details and options for input arguments when downloading the data.
}
\details{
OSM polygons are filtered by key-value attributes, where \verb{leisure:} is \code{park}, \code{garden} or \code{nature_reserve},
and \verb{access:} is not \code{no} or \code{private}, as well as the key-value pair \code{tourism:zoo}.
If \code{mutually_exclusive_with} is provided, \href{https://postgis.net/docs/ST_Intersection.html}{intersections} between
the output and these polygon(s) will be excluded using \code{polygons_mutually_exclude()}.
Polygons are then cleaned up using \code{polygons_clean()}.
}
\examples{
\dontrun{
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs

# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

# run function
get_parks_osm(place = city_boundaries,
              date = as.Date('2021-01-01'),
              snap_tolerance = 5,
              aggregate_polygons = 15,
              filename = 'public-parks_osm-polygons_2021-01-01.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_amenities_osm.R
\name{get_playgrounds_osm}
\alias{get_playgrounds_osm}
\title{Get public playgrounds from OpenStreetMap}
\usage{
get_playgrounds_osm(
  place,
  date = NULL,
  dir_raw = osmextract::oe_download_directory(),
  filename = NULL,
  ...
)
}
\arguments{
\item{place}{\code{sf} object (with projected coordinate reference system). Geographical area to match with the (\code{.osm.pbf}) file in the data archive.
Argument passed to \code{osmextract::oe_match()}.}

\item{date}{Date of OSM data snapshot to download. Object of class "Date" in format \verb{\%Y-\%m-\%d}. Refer to https://download.geofabrik.de
for the specific dates available. Defaults to \code{NULL} (download the latest available data).}

\item{dir_raw}{character. Directory to download the raw unprocessed OSM data. Passed to
argument \code{download_directory} in \code{osmextract::oe_read()}.}

\item{filename}{character (optional). File path to export output data (GeoJSON format).}

\item{...}{Other arguments passed to \code{osmextract::oe_read()}.}
}
\value{
The processed playgrounds (\code{sf} object).
}
\description{
Download and process OpenStreetMap (OSM) public playgrounds (points) within a specified geographical \code{place},
from the \href{https://download.geofabrik.de}{Geofabrik database}. It is a wrapper around
functions in the package \href{https://docs.ropensci.org/osmextract/index.html}{\code{osmextract}}, and
processes the downloaded files for subsequent analyses. Refer to package \code{osmextract} for
more details and options for input arguments when downloading the data.
}
\details{
OSM points are filtered by key-value attributes, where \code{leisure:playground}, and \verb{access:} is not \code{no} or \code{private}.
}
\examples{
\dontrun{
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs

# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

# run function
get_playgrounds_osm(place = city_boundaries,
                    date = as.Date('2021-01-01'),
                    filename = 'public-playgrounds_osm-points_2021-01-01.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygons_mutually_exclude.R
\name{polygons_mutually_exclude}
\alias{polygons_mutually_exclude}
\title{Make polygon input mutually exclusive (no overlaps) with another}
\usage{
polygons_mutually_exclude(input, mutually_exclude, snap_tolerance = 0)
}
\arguments{
\item{input}{\code{sf} object to process (with projected coordinate reference system).}

\item{mutually_exclude}{\code{sf} object. \code{input} polygons \href{https://postgis.net/docs/ST_Contains.html}{contained}
within this object, and \href{https://postgis.net/docs/ST_Intersection.html}{intersections} between the two, will be removed.
Should have the same coordinate reference system as \code{input}.}

\item{snap_tolerance}{numeric. Argument for \code{tolerance} level passed to \code{sf::st_snap()},
used to rectify nearly coincident edges between polygons before running \code{sf::st_contains()}.
Provided either as a units object (see \code{units::set_units()}),
or a number in the units of the coordinate reference system. Defaults to \code{0}.}
}
\value{
The processed \code{input} polygons (\code{sf} object).
}
\description{
Helper function to ensure that polygons (e.g. parks, beaches, informal nature areas)
are mutually-exclusive (i.e. non-overlapping).
}
\details{
Polygons vertices first are \href{https://postgis.net/docs/ST_Snap.html}{snapped} at a tolerance level
based on \code{snap_tolerance} (if specified) to rectify nearly coincident edges,
before removing \code{input} polygons \href{https://postgis.net/docs/ST_Contains.html}{contained} in \code{mutually_exclude}.
\href{https://postgis.net/docs/ST_Intersection.html}{Intersections} between the two are then removed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_amenities_osm.R
\name{get_trails_osm}
\alias{get_trails_osm}
\title{Get accessible trails from OpenStreetMap}
\usage{
get_trails_osm(
  place,
  date = NULL,
  dir_raw = osmextract::oe_download_directory(),
  filename = NULL,
  ...
)
}
\arguments{
\item{place}{\code{sf} object (with projected coordinate reference system). Geographical area to match with the (\code{.osm.pbf}) file in the data archive.
Argument passed to \code{osmextract::oe_match()}.}

\item{date}{Date of OSM data snapshot to download. Object of class "Date" in format \verb{\%Y-\%m-\%d}. Refer to https://download.geofabrik.de
for the specific dates available. Defaults to \code{NULL} (download the latest available data).}

\item{dir_raw}{character. Directory to download the raw unprocessed OSM data. Passed to
argument \code{download_directory} in \code{osmextract::oe_read()}.}

\item{filename}{character (optional). File path to export output data (GeoJSON format).}

\item{...}{Other arguments passed to \code{osmextract::oe_read()}.}
}
\value{
The processed trail lines (\code{sf} object).
}
\description{
Download and process OpenStreetMap (OSM) accessible trail (lines) within \code{place},
from the \href{https://download.geofabrik.de}{Geofabrik database}. It is a wrapper around
functions in the package \href{https://docs.ropensci.org/osmextract/index.html}{\code{osmextract}}, and
processes the downloaded files for subsequent analyses. Refer to package \code{osmextract} for
more details and options for input arguments when downloading the data.
}
\details{
OSM lines filtered by key-value attributes, where \verb{highway:} is \code{track}, \code{path}, \code{footway} or \code{cycleway}, and \verb{access:} is not \code{no} or \code{private}.
}
\examples{
\dontrun{
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs

# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

# run function
get_trails_osm(place = city_boundaries,
               date = as.Date('2021-01-01'),
               filename = 'accessible-trails_osm-lines_2021-01-01.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_density_polygonise.R
\name{pop_density_polygonise}
\alias{pop_density_polygonise}
\title{Polygonise population density raster to population count per (building) polygon}
\usage{
pop_density_polygonise(
  input_raster,
  write = FALSE,
  dsn,
  driver = "GeoJSON",
  overwrite = TRUE,
  delete_dsn = TRUE,
  ...
)
}
\arguments{
\item{input_raster}{Population density raster as a SpatRaster object (\code{terra::rast()}).
Output from \code{pop_dasymap()} may be used.}

\item{write}{Whether or not to export the output. Defaults to \code{FALSE}.}

\item{dsn}{Argument passed to \code{sf::st_write()}.}

\item{driver}{character. Argument passed to \code{sf::st_write()}. Defaults to \code{'GeoJSON'}.}

\item{overwrite}{logical. Argument passed to \code{sf::st_write()}. Defaults to \code{TRUE}.}

\item{delete_dsn}{Argument passed to \code{sf::st_write()}. Defaults to \code{TRUE}.}

\item{...}{Other arguments passed to \code{sf::st_write()}.}
}
\value{
Building \code{sf }polygons with column \code{popcount}.
}
\description{
Convert raster of population density (e.g. from \code{pop_dasymap()} output) to population count per polygon
(adjacent polygons with similar pixel values are merged).
Thus, if rasters supplied to \code{pop_dasymap()} are at a spatial resolution small enough to delineate individual buildings,
conversion to polygons using this function would reflect the population count per building.
Input data should have a projected coordinate reference system specific to the target area.
}
\examples{
\dontrun{
data(pop_sgp) # population census block polygons
data(landuse_sgp) # land use polygons


# transform to projected crs
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648))
landuse_sgp <- sf::st_transform(landuse_sgp, sf::st_crs(32648))


# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

buildings <- get_buildings_osm(place = city_boundaries,
                               date = as.Date('2021-01-01')) \%>\%
   dplyr::mutate(year = 2020)


# rasterise population, landuse & buildings
pop_rasters <- rasterise_pop(pop_sgp,
                             census_block = "subzone_n",
                             pop_count = "pop_count")

landuse_rasters <- rasterise_landuse(landuse_sgp,
                                     land_use = 'lu_desc',
                                     subset = c('1' = 'RESIDENTIAL',
                                                '2' = 'COMMERCIAL & RESIDENTIAL',
                                                '3' = 'RESIDENTIAL WITH COMMERCIAL AT 1ST STOREY',
                                                '4' = 'RESIDENTIAL / INSTITUTION'),
                                     sf_pop = pop_sgp,
                                     match_landuse_pop = 'recent')

buildings_rasters <- rasterise_buildings(buildings,
                                         proxy_pop_density = 'levels',
                                         year = 'year',
                                         sf_pop = pop_sgp,
                                         sf_landuse = landuse_sgp,
                                         match_buildings_pop = 'closest')


# perform dasymetric mapping on selected year (2020)
popdens_raster <- pop_dasymap(pop_polygons = pop_rasters$pop_polygons[[2]],
                              pop_perblock_count = pop_rasters$pop_count[[2]],
                              pop_perblock_density = pop_rasters$pop_density[[2]],
                              land_relative_density = buildings_rasters[[2]],
                              filename = 'buildings_popdensity.tif',
                              wopt = list(gdal=c('COMPRESS=LZW')))


# finally, convert to population count per building polygon
pop_density_polygonise(input_raster = popdens_raster,
                       write = TRUE,
                       dsn = 'buildings_popcount.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recre_supply.R
\name{recre_supply}
\alias{recre_supply}
\title{Calculate cumulative supply of an attribute from all parks to each building}
\source{
Rossi, S. D., Byrne, J. A., & Pickering, C. M. (2015). The role of distance in peri-urban national park use:
Who visits them and how far do they travel?. Applied Geography, 63, 77-88.

Tu, X., Huang, G., Wu, J., & Guo, X. (2020). How do travel distance and park size influence urban park visits?.
Urban Forestry & Urban Greening, 52, 126689.
}
\usage{
recre_supply(park_attribute, dist_matrix, c = 1)
}
\arguments{
\item{park_attribute}{numeric vector. Amount of a specific attribute per park.
Length of vector is equal to the number of parks considered.}

\item{dist_matrix}{Matrix containing buildings (rows) and their pairwise distances to each park (columns).
Order of parks (columns) should be identical to the order of parks (elements) in \code{park_attribute}.}

\item{c}{Coefficient determining rate of decay in recreation supply with increasing distance.}
}
\value{
A numeric vector of the cumulative supply value per building (row) in the input matrix \code{dist_matrix}.
The length of the vector equals to the number of buildings considered.
}
\description{
Sum up the total amount of a specific park attribute supplied to each building polygon.
The amount per building depends on the distances between that particular building and all parks;
attributes from parks further away are generally reduced, an effect also known as
the 'distance decay' (Rossi et al., 2015; Tu et al., 2020).
}
\details{
The supply \eqn{S} of the park attribute is calculated based on the following equation:

\eqn{S = \sum\limits_{i=1}^{n} s_{i} \cdot e^{-cd_{i}}}

where

\eqn{S} = Total supply of a specific park attribute to the building from parks \eqn{i};
\eqn{i = 1,2,3}... \eqn{n}, \eqn{n} = total number of parks

\eqn{s_{i}} = Supply of a specific park attribute from park \eqn{i}.
A perfect positive linear association is assumed, since the focus is on supply metrics.

\eqn{d_{i}} = Distance in kilometres from the building to park \eqn{i} (e.g. Euclidean, Manhattan, etc.).

\eqn{c} = Coefficient determining rate of decay in supply \eqn{i} with increasing distance.
}
\examples{
\dontrun{
data(parks_sgp) # load park polygons
data(buildings_pop_sgp) # load building polygons w population counts

# transform to projected crs
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648))
buildings_pop_sgp <- sf::st_transform(buildings_pop_sgp, sf::st_crs(32648))


# Calculate pairwise distances between (the centroid of) each building & all parks
d_matrix <- buildings_pop_sgp \%>\%
  st_centroid() \%>\%
  st_distance(parks_sgp) # euclidean distance

m_dist <- m_dist / 1000 # convert distances to km


# run function for a specific park attribute (e.g. area)
recre_supply(park_attribute = parks_sgp$area,
             dist_matrix = d_matrix,
             c = 0.3) # example value for distance decay coefficient c
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_buildings_osm.R
\name{get_buildings_osm}
\alias{get_buildings_osm}
\title{Get building polygons from OpenStreetMap}
\usage{
get_buildings_osm(
  place,
  date = NULL,
  dir_raw = osmextract::oe_download_directory(),
  filename = NULL,
  driver = "GeoJSON",
  delete_dsn = TRUE,
  append = NA,
  ...
)
}
\arguments{
\item{place}{\code{sf} object (with projected coordinate reference system). Geographical area to match with the (\code{.osm.pbf}) file in the data archive.
Argument passed to \code{osmextract::oe_match()}.}

\item{date}{Date of OSM data snapshot to download. Object of class "Date" in format \verb{\%Y-\%m-\%d}. Refer to https://download.geofabrik.de
for the specific dates available. Defaults to \code{NULL} (download the latest available data).}

\item{dir_raw}{character. Directory to download the raw unprocessed OSM data. Passed to
argument \code{download_directory} in \code{osmextract::oe_read()}.}

\item{filename}{character (optional). File path to export output data.}

\item{driver}{character (optional). Name of driver used to export output data, passed to \code{sf::st_write()}.
Defaults to "GeoJSON".}

\item{delete_dsn}{logical (optional). Passed to \code{sf::st_write()}.}

\item{append}{defaults to \code{NA}, which raises an error if a layer exists. Passed to \code{sf::st_write()}.}

\item{...}{Other arguments passed to \code{osmextract::oe_read()}.}
}
\value{
The processed building polygons (\code{sf} object).
}
\description{
Download and process OpenStreetMap (OSM) building polygons within a specified geographical \code{place},
from the \href{https://download.geofabrik.de}{Geofabrik database}. It is a wrapper around
functions in the package \href{https://docs.ropensci.org/osmextract/index.html}{\code{osmextract}}, and
processes the downloaded files for subsequent analyses. Refer to package \code{osmextract} for
more details and options for input arguments when downloading the data.
}
\details{
Data is filtered by key-value attributes, where \verb{building:} is not \code{NULL}.
The column \code{levels} is derived from \code{building:levels}; values were set to \code{1} if the
extracted value is empty or \code{NA}, and set to \code{NA} if \verb{≤ 0} (i.e. underground);
values were then rounded up to the nearest whole number.
The column \code{area_m2} represents the building footprint area,
and \code{floorarea_m2} is calculated by multiplying the \code{area_m2} by the number of \code{levels}.
}
\examples{
\dontrun{
data(pop_sgp)
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648)) # transform to projected crs

# merge all census blocks for chosen year (2020) into single multi-polygon
# function requires that polygons are merged
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

# run function
get_buildings_osm(place = city_boundaries,
                  date = as.Date('2021-01-01'),
                  filename = 'buildings_osm-polygons_2021-01-01.geojson')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parks_calc_attributes.R
\name{parks_calc_attributes}
\alias{parks_calc_attributes}
\title{Calculate attributes per park polygon based on supplied data}
\usage{
parks_calc_attributes(
  parks,
  data_points = list(),
  data_lines = list(),
  data_rasters = list(),
  rasters_summarise_fun = NULL,
  raster_min_patch_size = units::set_units(0, "m^2"),
  raster_edge = NULL,
  relative = TRUE,
  filename = NULL,
  ...
)
}
\arguments{
\item{parks}{\code{sf} polygons (with projected coordinate reference system).}

\item{data_points}{Named list of \code{sf} object(s) containing data of geometry type \code{POINT} or \code{MULTIPOINT}.
List names are used to name the park attributes (columns) in the output.}

\item{data_lines}{Named list of \code{sf} object(s) containing data of geometry type \code{LINESTRING} or \code{MULTILINESTRING}.
List names are used to name the park attributes (columns) in the output.}

\item{data_rasters}{Named list of single \code{SpatRaster} object(s) from \code{terra::rast()}.
List names are used to name the park attributes (columns) in the output.}

\item{rasters_summarise_fun}{Function to summarise the raster data passed to \code{terra::extract()}.
Defaults to \code{NULL}, which tabulates the sum of each unique value in the raster (i.e. classified raster,
with numbers each representing a specific class/category) using the function \code{raster_class_area()}.}

\item{raster_min_patch_size}{Minimum patch size to be included when tabulating the sum of each raster class.
Only relevant if \code{rasters_summarise_fun = NULL}. Provided either as a units object (see \code{units::set_units()}),
or a number in the units of the coordinate reference system. Defaults to \code{0} m^2.}

\item{raster_edge}{numeric. Option to calculate total edge length of (classified) raster(s)
associated with each of the \code{parks} polygons, using the function \code{raster_edge_length()}. The total edge length of raster patches contained \href{https://postgis.net/docs/ST_Within.html}{within}
\code{parks} will be calculated, and patches in close proximity can be included by increasing the value of this argument;
the length of \code{parks} borders that \href{https://postgis.net/docs/ST_Intersection.html}{intersect} the buffered raster will be calculated.
This provides a way include patches in close proximity to the \code{parks} (e.g. total waterfront length close to parks).
Note that each class (unique value in the raster) will be summarised, including \code{0}; convert pixels to \code{NA} if you wish to exclude them.
Provided either as a units object (see \code{units::set_units()}),
or a number in the units of the coordinate reference system. Defaults to \code{NULL} (not calculated). Set to \code{0} to include patches
\href{https://postgis.net/docs/ST_Within.html}{within} \code{parks} only.}

\item{relative}{logical. Whether or not to calculate relative amounts
(e.g. point density, ratio of line-to-perimeter length, proportional area). Defaults to \code{TRUE}.}

\item{filename}{character (optional). File path to export output data (GeoJSON format).}

\item{...}{Other arguments passed to \code{terra::extract()}.}
}
\value{
\code{parks} with added columns containing the summaries of basic park attributes
(area and perimeter), as well as summaries of each of the supplied datasets.
The summary method depends on geometry type (i.e. points, lines or rasters) and supplied arguments.
Some examples:
\describe{
\item{area}{Area of park polygon.}
\item{perimeter}{Perimeter of park polygon.}
\item{< object name in \code{data_points} >_count}{Total count of points within park polygon.}
\item{< object name in \code{data_points} >_ptdensity}{Total count divided by area of the park polygon (point density).
Included if argument \code{relative} set to \code{TRUE}.}
\item{< object name in \code{data_lines} >_length}{Sum of line lengths within park polygon.}
\item{< object name in \code{data_lines} >_length_perim_ratio}{Ratio of line-to-perimeter length of the park polygon.
Included if argument \code{relative} set to \code{TRUE}.}
\item{< object name in \code{data_rasters} >}{Summarised values of a specific raster class
(depends on function provided in the \code{rasters_summarise_fun} argument).}
\item{< object name in \code{data_rasters} >< class value >_area}{Total area of a specific raster class
(if \code{rasters_summarise_fun = NULL}).}
\item{< object name in \code{data_rasters} >< class value >_area_pct}{Percentage area of a specific raster class
(if \code{rasters_summarise_fun = NULL}). Included if argument \code{relative} set to \code{TRUE}.}
\item{< object name in \code{data_rasters} >< class value >_length}{Total edge length of a specific raster class.}
\item{< object name in \code{data_rasters} >< class value >_length_perim_ratio}{Edge-to-perimeter length of a specific raster class.
Included if argument \code{relative} set to \code{TRUE}.}
}
}
\description{
Summaries will be calculated for each of the supplied datasets (\code{data_points}, \code{data_lines}, or \code{data_rasters})
and appended to the \code{parks} data as additional columns.
Ensure that all have a (projected) coordinate reference system similar to \code{parks}.
}
\examples{
\dontrun{
data(pop_sgp) # city census blocks
data(parks_sgp) # park polygons


# transform to projected crs
pop_sgp <- sf::st_transform(pop_sgp, sf::st_crs(32648))
parks_sgp <- sf::st_transform(parks_sgp, sf::st_crs(32648))


# get playground points (example attribute)
city_boundaries <- pop_sgp \%>\%
   dplyr::filter(year == 2020) \%>\%
   sf::st_union() \%>\%
   sf::st_as_sf() \%>\%
   smoothr::fill_holes(threshold = units::set_units(1, 'km^2'))  \%>\% # clean up
   smoothr::drop_crumbs(threshold = units::set_units(1, 'km^2'))  \%>\%
   sf::st_make_valid()

playgrounds <- get_playgrounds_osm(place = city_boundaries,
                                   date = as.Date('2021-01-01'))
point_list <- list(playgrounds) # convert to list (can add other point data too)
names(point_list) <- c("playground") # name each element in list


# calculate playground point count & density per park
parks_calc_attributes(parks = parks_sgp,
                      data_points = point_list,
                      relative = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parks_calc_attributes_helpers.R
\name{raster_edge_length}
\alias{raster_edge_length}
\title{Calculate the total edge length of a raster associated with each park polygon}
\usage{
raster_edge_length(
  polygons,
  raster,
  raster_min_patch_size = units::set_units(0, "m^2"),
  raster_edge_buffer = units::set_units(0, "m^2"),
  relative = TRUE
)
}
\arguments{
\item{polygons}{\code{sf} (with projected coordinate reference system).}

\item{raster}{\code{SpatRaster} object from \code{terra::rast()}. Should have a (projected) coordinate reference system similar to \code{polygons}.}

\item{raster_min_patch_size}{Minimum patch size to be included in results.
Provided either as a units object (see \code{units::set_units()}),
or a number in the units of the coordinate reference system. Defaults to \code{0} m^2.}

\item{raster_edge_buffer}{numeric. Specify buffer distance to add to polygonised raster; the total edge length of \code{polygons} that
\href{https://postgis.net/docs/ST_Intersection.html}{intersect} the buffered \code{raster}
will be summed up together with the total edge length contained \href{https://postgis.net/docs/ST_Within.html}{within} the \code{polygons}.
Defaults to \code{0} (only patches fully contained \href{https://postgis.net/docs/ST_Within.html}{within} polygons will be considered).
Provided either as a units object (see \code{units::set_units()}), or a number in the units of the coordinate reference system.}

\item{relative}{logical. Whether or not to calculate relative amounts
(i.e. ratio of edge-to-perimeter length). Defaults to \code{TRUE}.}
}
\value{
\code{polygons} with added column(s) \verb{< class value >_length}, and \verb{< class value >_length_perim_ratio} if \code{relative} is set to \code{TRUE}.
Note that the value \code{0} will be summarised; convert pixels to \code{NA} if you wish to exclude them.
}
\description{
Helper function within \code{parks_calc_attributes()}.
The total edge length of (classified) raster patches contained \href{https://postgis.net/docs/ST_Within.html}{within}
the \code{polygons} will be calculated. Additionally, the argument \code{raster_edge_buffer} provides a way include patches
in close proximity to the \code{polygons}. The total edge lengths are summed together and appended to the \code{polygons}
data as additional columns (or one column, if there is only one raster class).
Note that this operation may take a while to complete, as it involves the conversion of rasters to polygons (and vice versa).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygons_clean.R
\name{polygons_clean}
\alias{polygons_clean}
\title{Clean up polygons after download}
\usage{
polygons_clean(
  input,
  snap_tolerance = 5,
  min_area = units::set_units(0, "m^2"),
  aggregate_polygons = 15
)
}
\arguments{
\item{input}{\code{sf} object (with projected coordinate reference system).}

\item{snap_tolerance}{numeric. Argument for \code{tolerance} level passed to \code{sf::st_snap()},
used to rectify nearly coincident edges between polygons. Provided either as a units object (see \code{units::set_units()}),
or a number in the units of the coordinate reference system.
After \code{sf::st_snap()} is run, a polygon that is \href{https://postgis.net/docs/ST_Covers.html}{covered} by another will be removed.
Defaults to \code{5}. Set to \code{0} if you do not wish to rectify minor overlaps.}

\item{min_area}{numeric. Specify minimum area of each polygon to be retained in the output,
passed to argument \code{threshold} in \code{smoothr::drop_crumbs()}.
Provided either as a units object (see \code{units::set_units()}), or a number in the units of
the coordinate reference system. Defaults to \code{0} m^2.}

\item{aggregate_polygons}{numeric. Argument for \code{dist} passed to \code{sf::st_buffer()}.
Buffered polygons that overlap will be aggregated into multipolygons.
Set to \code{NULL} if you do not wish to aggregate to multipolygons.}
}
\value{
The processed polygons (\code{sf} object).
}
\description{
Helper function to process and clean up polygons downloaded from OpenStreetMap.
}
\details{
Polygons \href{https://postgis.net/docs/ST_Contains.html}{contained} in or \href{https://postgis.net/docs/ST_Overlaps.html}{overlapping} others are removed.
Polygons vertices are then \href{https://postgis.net/docs/ST_Snap.html}{snapped}
at a tolerance level based on \code{snap_tolerance} (if specified) to rectify nearly coincident edges,
and polygons \href{https://postgis.net/docs/ST_Covers.html}{covered} by others are subsequently removed. Polygons are then aggregated
to multipolygons based on distance set in \code{aggregate_polygons} (if specified).
All spatial relations (italicised) follow the \href{https://en.wikipedia.org/wiki/DE-9IM}{DE-9IM} standard.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-singapore.R
\docType{data}
\name{landuse_sgp}
\alias{landuse_sgp}
\title{Land use master plan for the city of Singapore}
\format{
\code{sf} polygons of land use released in years 2014 and 2019.
Data is in the 'long' format (rows are repeated for each year).
The following columns are used in this package:
\describe{
\item{lu_desc}{Land use type (description)}
\item{year}{Year of Master Plan}
}
}
\source{
Polygons from the Singapore Land Use Master Plans released in the years
\href{https://data.gov.sg/dataset/master-plan-2014-land-use}{2014} and
\href{https://data.gov.sg/dataset/master-plan-2019-land-use-layer}{2019}.
Made available under the terms of the \href{https://data.gov.sg/open-data-licence}{Singapore Open Data Licence version 1.0}.
}
\usage{
landuse_sgp
}
\description{
Example dataset containing Singapore Master Plan Land Use Zones for the years
\href{https://data.gov.sg/dataset/master-plan-2014-land-use}{2014} and
\href{https://data.gov.sg/dataset/master-plan-2019-land-use-layer}{2019}.
Released under the \href{https://data.gov.sg/open-data-licence}{Singapore Open Data License}.
}
\examples{
data(landuse_sgp)
head(landuse_sgp)
}
\keyword{datasets}
