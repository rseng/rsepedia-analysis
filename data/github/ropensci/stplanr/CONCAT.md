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

# stplanr <a href='https://docs.ropensci.org/stplanr/'><img src='man/figures/logo.png' align="right" height=215/></a>

<!-- [![Build Status](https://travis-ci.org/ropensci/stplanr.svg?branch=master)](https://travis-ci.org/ropensci/stplanr) -->

[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/stplanr)](https://github.com/r-hub/cranlogs.app)
[![](https://cranlogs.r-pkg.org/badges/grand-total/stplanr)](https://cran.r-project.org/package=stplanr)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/stplanr)](https://cran.r-project.org/package=stplanr)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/10_status.svg)](https://github.com/ropensci/software-review/issues/10)
[![R-CMD-check](https://github.com/ropensci/stplanr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/stplanr/actions)

**stplanr** is a package for sustainable transport planning with R.

It provides functions for solving common problems in transport planning
and modelling, such as how to best get from point A to point B. The
overall aim is to provide a reproducible, transparent and accessible
toolkit to help people better understand transport systems and inform
policy, as outlined in a
[paper](https://journal.r-project.org/archive/2018/RJ-2018-053/index.html)
about the package, and the potential for open source software in
transport planning in general, published in the [R
Journal](https://journal.r-project.org/).

The initial work on the project was funded by the Department of
Transport
([DfT](https://www.gov.uk/government/organisations/department-for-transport))
as part of the development of the Propensity to Cycle Tool (PCT), a web
application to explore current travel patterns and cycling potential at
zone, desire line, route and route network levels (see
[www.pct.bike](https://www.pct.bike/) and click on a region to try it
out). The basis of the methods underlying the PCT is origin-destination
data, which are used to highlight where many short distance trips are
being made, and estimate how many could switch to cycling. The results
help identify where cycleways are most needed, an important component of
sustainable transport planning infrastructure engineering and policy
[design](https://www.icevirtuallibrary.com/doi/abs/10.1680/dfct.63495.001).

See the package vignette (e.g. via `vignette("introducing-stplanr")`) or
an [academic paper on the Propensity to Cycle Tool
(PCT)](https://dx.doi.org/10.5198/jtlu.2016.862) for more information on
how it can be used. This README provides some basics.

Much of the work supports research undertaken at the Leeds’ Institute
for Transport Studies ([ITS](https://environment.leeds.ac.uk/transport))
but **stplanr** should be useful to transport researchers and
practitioners needing free, open and reproducible methods for working
with geographic data everywhere.

## Key functions

Data frames representing flows between origins and destinations must be
combined with geo-referenced zones or points to generate meaningful
analyses and visualisations of ‘flows’ or origin-destination (OD) data.
**stplanr** facilitates this with `od2line()`, which takes flow and
geographical data as inputs and outputs spatial data. Some example data
is provided in the package:

``` r
library(stplanr)
```

Let’s take a look at this data:

``` r
od_data_sample[1:3, 1:3] # typical form of flow data
#> # A tibble: 3 x 3
#>   geo_code1 geo_code2   all
#>   <chr>     <chr>     <dbl>
#> 1 E02002361 E02002361   109
#> 2 E02002361 E02002363    38
#> 3 E02002361 E02002367    10
cents_sf[1:3,] # points representing origins and destinations
#>       geo_code  MSOA11NM percent_fem  avslope             geometry
#> 1708 E02002384 Leeds 055    0.458721 2.856563 -1.546463, 53.809517
#> 1712 E02002382 Leeds 053    0.438144 2.284782 -1.511861, 53.811611
#> 1805 E02002393 Leeds 064    0.408759 2.361707 -1.524205, 53.804098
```

These datasets can be combined as follows:

``` r
travel_network <- od2line(flow = od_data_sample, zones = cents_sf)
w <- flow$all / max(flow$all) *10
plot(travel_network, lwd = w)
```

<img src="man/figures/README-plot1-1.png" width="100%" />

**stplanr** has many functions for working with OD data. See the
[`stplanr-od`](https://docs.ropensci.org/stplanr/articles/stplanr-od.html)
vignette for details.

The package can also allocate flows to the road network, e.g. with
[CycleStreets.net](https://www.cyclestreets.net/api/) and the
OpenStreetMap Routing Machine
([OSRM](https://github.com/Project-OSRM/osrm-backend)) API interfaces.
These are supported in `route_*()` functions such as
`route_cyclestreets` and `route_osrm()`:

Routing can be done using a range of back-ends and using lat/lon or
desire line inputs with the `route()` function, as illustrated by the
following commands which calculates the route between Fleet Street and
Southwark Street over the River Thames on Blackfriars Bridge in London:

``` r
library(osrm)
#> Data: (c) OpenStreetMap contributors, ODbL 1.0 - http://www.openstreetmap.org/copyright
#> Routing: OSRM - http://project-osrm.org/
trip <- route(
  from = c(-0.11, 51.514),
  to = c(-0.10, 51.506),
  route_fun = osrmRoute,
  returnclass = "sf"
  )
#> Most common output is sf
mapview::mapview(trip)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

You can also use and place names, found using the Google Map API:

``` r
trip2 <- route(
  from = "Leeds",
  to = "Bradford",
  route_fun = osrmRoute,
  returnclass = "sf"
  )
#> Most common output is sf
mapview::mapview(trip2)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

We can replicate this call multiple times with the `l` argument in
`route()`:

``` r
desire_lines <- travel_network[2:6, ]
```

Next, we’ll calculate the routes:

``` r
routes <- route(
  l = desire_lines,
  route_fun = osrmRoute,
  returnclass = "sf"
  )
mapview::mapview(routes) +
  mapview::mapview(desire_lines, color = "red")
```

<img src="man/figures/README-plot2-1.png" width="100%" />

<!-- The resulting routes will look something like this: -->

For more examples, `example("route")`.

`overline()` takes a series of route-allocated lines, splits them into
unique segments and aggregates the values of overlapping lines. This can
represent where there will be most traffic on the transport system, as
demonstrated in the following code chunk.

``` r
routes$foot <- desire_lines$foot
rnet <- overline(routes, attrib = "foot")
#> 2020-09-03 22:24:08 constructing segments
#> 2020-09-03 22:24:08 building geometry
#> 2020-09-03 22:24:08 simplifying geometry
#> 2020-09-03 22:24:08 aggregating flows
#> 2020-09-03 22:24:08 rejoining segments into linestrings
```

The resulting route network, with segment totals calculated from
overlapping parts for the routes for walking, can be visualised as
follows:

``` r
plot(rnet["foot"], lwd = rnet$foot)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

The above plot represents the number walking trips made (the ‘flow’)
along particular segments of a transport network.

<!-- (results not shown): -->

## Policy applications

The examples shown above, based on tiny demonstration datasets, may not
seem particularly revolutionary. At the city scale, however, this type
of analysis can be used to inform sustainable transport policies, as
described in papers [describing the Propensity to Cycle
Tool](https://www.jtlu.org/index.php/jtlu/article/view/862/859) (PCT),
and its [application to calculate cycling to school
potential](https://www.sciencedirect.com/science/article/pii/S2214140518301257)
across England.

Results generated by **stplanr** are now part of national government
policy: the PCT is the recommended tool for local and regional
authorities developing strategic cycle network under the Cycling and
Walking Infrastructure Strategy
([CWIS](https://www.gov.uk/government/publications/cycling-and-walking-investment-strategy)),
which is part of the Infrastructure Act
[2015](https://www.legislation.gov.uk/ukpga/2015/7/contents/enacted).
**stplanr** is helping dozens of local authorities across the UK to
answer the question: where to prioritise investment in cycling? In
essence, stplanr was designed to support sustainable transport policies.

There are many other research and policy questions that functions in
**stplanr**, and other open source software libraries and packages, can
help answer. At a time of climate, health and social crises, it is
important that technology is not only sustainable itself (e.g. as
enabled by open source communities and licenses) but that it contributes
to a sustainable future.

## Installation

To install the stable version, use:

``` r
install.packages("stplanr")
```

The development version can be installed using **devtools**:

``` r
# install.packages("devtools") # if not already installed
devtools::install_github("ropensci/stplanr")
library(stplanr)
```

stplanr depends on rgdal, which can be tricky to install.

### Installing stplanr on Linux and Mac

**stplanr** depends on **sf**. Installation instructions for Mac, Ubuntu
and other Linux distros can be found here:
<https://github.com/r-spatial/sf#installing>

## Functions, help and contributing

The current list of available functions can be seen on the package’s
website at
[docs.ropensci.org/stplanr/](https://docs.ropensci.org/stplanr/), or
with the following command:

``` r
lsf.str("package:stplanr", all = TRUE)
```

To get internal help on a specific function, use the standard way.

``` r
?od2line
```

To contribute, report bugs or request features, see the [issue
tracker](https://github.com/ropensci/stplanr/issues).

## Further resources / tutorials

Want to learn how to use open source software for reproducible
sustainable transport planning work? Now is a great time to learn.
Transport planning is a relatively new field of application in R.
However, there are already some good resources on the topic, including
(any further suggestions: welcome):

  - The Transport chapter of *Geocomputation with R*, which provides a
    broad introduction from a geographic data perspective:
    <https://geocompr.robinlovelace.net/transport.html>
  - The **stplanr** paper, which describes the context in which the
    package was developed:
    <https://journal.r-project.org/archive/2018/RJ-2018-053/index.html>
    (please cite this if you use **stplanr** in your work)
  - The `dodgr` vignette, which provides an introduction to routing in
    R: <https://atfutures.github.io/dodgr/>

## Meta

  - Please report issues, feature requests and questions to the [github
    issue tracker](https://github.com/ropensci/stplanr/issues)
  - License: MIT
  - Get citation information for **stplanr** in R doing
    `citation(package = 'stplanr')`
  - This project is released with a [Contributor Code of
    Conduct](https://github.com/ropensci/stplanr/blob/master/CONDUCT.md).
    By participating in this project you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# stplanr 0.8.6 (November 2021)

- `route()` checks the CRS and gives an appropriate warning if it is projected
- The package now Suggests igraph rather than depending on it making it easier to install

# stplanr 0.8.5

- No longer Suggests stats19 package

# stplanr 0.8.4

- No longer Suggests bench package
- Tests pass when internet unavailable (#469)

# stplanr 0.8.3

- Outputs of `line_via()` now have correct CRS
- `calc_catchment()` have been deprecated because the generate warnings
- Minor fixes and improvements in the package's documentation

# stplanr 0.8.2

- Bug fixed related to the `route()` function: it silently failed when `pbapply` not installed. The package was promoted from Suggests to Imports (#460)

# stplanr 0.8.1

- Thanks to the `styler` package, code in `stplanr` now adheres to a consistent style (using arrow `<-` assignment despite the maintainer's default of equals `=` assignment - many `=` had been introduced accidentally!)
- New function `rnet_group()` groups route network elements (#403)
- `overline()` now converts `MULTILINESTRING` geometries to `LINESTRINGS` automatically
- Routing on the network using `sum_network_routes()` now returns a linestring for routes that start where they end (i.e. no travel) (#444)
- Routing using `route_local()` fixed: the `l` argument now works (#448)
- New `route_osrm()` function (#449)
- `rnet_breakup_vertices()` is now way faster, thanks to Andrea Gilardi and others (#416)
- `rnet_group()` is now a generic function that works with `sfNetwork` objects (#455)
- `route_osrm()` provides easy access to multi-modal routing (#449)
- Bug in `route()` function's calculation of bbox attributes fixed (#452)

# stplanr 0.8.0

- New function `line_breakup()` breaks-up lines (#434)
- Minor documentation fixes, including (#431)

# stplanr 0.7.2

- Support `sf` objects for `toptail_buff()`
- Updated examples for tests and README to enable tests to pass on legacy versions of PROJ (#423)
- Reduce vulnerability to upcoming upstream changes (#426)

# stplanr 0.7.1

- Tweaks to the documentation and examples for CRAN tests

# stplanr 0.7.0

- Issue with `line2points()` on some set-ups fixed (#418)
- Old `mapshape()`, `line_match()` and `gclip()` functions deprecated, maintained alternatives can now be found in other packages.
- `sp` objects updated for latest version of `sp` (#364)
- `sf` objects updated to support more recent CRS encoding (#393)
- Deprecated functions including `od_aggregate()`, `onewayid()`, `gtfs2sldf()`, and `od_radiation()` have been removed

# stplanr 0.6.2

## New features

- New interface to Google Directions API via `mapsapi` package (#410)
- New `quiet` argument in `overline()` makes it less chatty

# stplanr 0.6.1

## BUG FIXES

- `route()` failed when `data.table` not installed (#408)

# stplanr 0.6.0

## BUG FIXES

- Bug in `SpatialLinesNetwork()` fixed thanks to Andrea Gilardi (#394)
- Updated documentation for finding shortest paths
- Check `start` and `end` arguments in short path calculations are numeric
- `dodgr` removed as Suggests until it's back on CRAN (#398)
- Updates to `dplyr` code to prevent warnings when using the dev version (#395)

## NEW FEATURES

- Improvements to `route()` allowing you to save a list of raw outputs and use `data.table` for faster performance if available
- A fleet of new `route_rolling_*()` functions have been added:
  - `route_rolling_gradient()` calculates a rolling gradient from elevation and distance data at the segment level
  - `route_rolling_average()` calculates the rolling average of values
  - `route_rolling_diff()` calculates the rolling difference between each value and the next
- `geo_toptail()` function now has `sf` implementation (#390)
- New `cl` argument in `route()` for parallel routing (#388)
- New and updated datasets representing `od_data_sample` in different ways: `od_data_lines` and `od_data_routes`
- `route_graphhopper()` deprecated (#389)
- Old functions that use legacy `sp` and `dplyr` code, `sp_aggregate` and `od_aggregate`, have been deprecated
- New work-in-progress `route_slope*()` functions

# stplanr 0.5.2

## BUG FIXES

- Documentation fixes (#384)
- A bug affecting `route()` function calls when `pbapply` was not installed has been fixed (#386)
- `oneway()` has been deprecated in favour of faster and easier-to-maintain function `od_oneway()` (also see the in-development `od` package) (#387)
- Various other changes have been made to accomodate the dev version of `dplyr` (#383)

# stplanr 0.5.1

- Changes for compatibility with R 4.0.0

# stplanr 0.5.0

- `route_graphhopper()` should now work with a local graphhopper instance. See https://github.com/ropensci/stplanr/pull/369
- The old `line2route()` function now works with routing functions that return `sf` objects
- The new `route()` function is now more resilient, providing a progress bar if you have `pbapply` package installed, returning a result even if some routes failed, and telling you which routes failed
- The package has fewer dependencies, `readr`, `openxlsx` and `lubridate` removed
- Deprecated function `buff_geo()` removed

## BUG FIXES

- `destination` now works again as an argument in `line2route()` (#368)
- `overline()` now accepts `sf` objects regardless of the name of the geometry column
- `line2points()` now works with `sfc` objects (#361)

# stplanr 0.4.1

## NEW FEATURES

- Better error messages if `od2line()` fails due to non-matching ids
- Improved documentation of `od2line()` in the vignette: https://docs.ropensci.org/stplanr/articles/stplanr-od.html#non-matching-ids

# stplanr 0.4.0

## NEW FEATURES

- A family of new functions, including `route_split()`, `rnet_add_node()` and `sln_add_node()` for adding new nodes to routes, route networks and `sfNetwork` objects, closing [#342](https://github.com/ropensci/stplanr/issues/342)
- Updated vignette on route networks, solving [#237](https://github.com/ropensci/stplanr/issues/237), which can be found here: https://docs.ropensci.org/stplanr/articles/stplanr-route-nets.html
- [Fix](https://github.com/ropensci/stplanr/commit/592fba2a6d191135d036af73e7c902c3ef4f758c) in `line2points()`
- `line_to_points()` function removed in favour of `line2point()`.
- New function `sln_clean_graph()` removes unconnected elements of `sfNetwork` objects. Credit to Andrea Gilardi. See (#344).
- New functions `rnet_breakup_vertices()` and `line2vertices()` for breaking up linestrings representing route networks into smaller segments, in preparation for routing. See (#282) (which these functions address) and PR (#347) for details.

## BUG FIXES

- Bugs in `route_dodgr()` and associated examples fixed (#348)
- Annoying message printed on load removed (#355)

## OTHER

- Andrea Gilardi added as author.
- Deprecated functions related to road crash (STATS19) data removed

# stplanr 0.3.1

- stplanr now has a logo! See [#334](https://github.com/ropensci/stplanr/issues/334)
- `line_to_points()` depreciated in favour of `od2line()`, the latter function name being more consistent with the package's other functions
- `line2pointsn()` now works with `sf` objects
- Documentation fixes - see [#329](https://github.com/ropensci/stplanr/issues/329)

## OTHER

- Various improvements made to the `stplanr-od` vignette, thanks to Edward Leigh
- URLs updated to link to stplanr's new, official website: https://docs.ropensci.org/stplanr/

# stplanr 0.3.0

## NEW FEATURES

- New functions `od_to_odmatrix()` and `odmatrix_to_od()` to convert between matrix forms of origin-destination data
- New function `od_oneway()` replaces `onewayid()`, works better and is twice as fast
- New `od_id*()` functions provide a range of ways to convert origin-destination pair IDs into a single ID. See [Stackoverflow](https://stackoverflow.com/questions/57235601/how-to-identify-duplicated-ordered-pairs-efficiently/57236658#57236658) and the [issue tracker](https://github.com/ropensci/stplanr/issues/321)
- New vignette [`stplanr-od`](https://docs.ropensci.org/stplanr/articles/stplanr-od.html) provides detailed documentation on the package's OD data handling capabilities

# stplanr 0.2.10

- Fix in documentation. See [#311](https://github.com/ropensci/stplanr/issues/311)

# stplanr 0.2.9

## NEW FEATURES

- New functions `od_aggregate_from()` and `od_aggregate_to()` provide easy ways to aggregate origin-destination pairs. See [#303](https://github.com/ropensci/stplanr/pull/303).
- Updated `overline2()` is now faster and better documented (#307)
- Updates to `route_dodgr()` function, which provides an interface to the [dodgr](https://github.com/ATFutures/dodgr) package, accepts wider range of inputs
- Better website and updated function list. See https://ropensci.github.io/stplanr/index.html
- The `sf` method for `overline()` has been updated so it calls the much faster `overline2()` function
- Updated documentation for `route_local()`

## BUG FIXES

- Bug in `sum_network_routes()` fixed. See [#267](https://github.com/ropensci/stplanr/issues/267)

# stplanr 0.2.8

## NEW FEATURES

- The stplanr paper has been published! See it here: https://journal.r-project.org/archive/2018/RJ-2018-053/index.html
- STATS19 functions such as `dl_stats19()` are depreciated. They have been split-out into the new package [`stats19`](https://github.com/ropensci/stats19)
- `route_dodgr()` has now been implemented
- A new function `overline2()` has been added, thanks to Malcolm Morgan. This is faster than `overline()`.
- A substantial refactoring operation has begun. This has resulted in fewer lines of code in the `od` functions, a new `stplanr::od_coords2line()` function, and more support of `sf`
- `route_dodgr()` has been added
- A new example dataset, `osm_net_example`, has been added for local routing purposes.
- A citation to the package has been added. Try `citation("stplanr")`
- The package has a shiny new website thanks to `@maelle`: https://ropensci.github.io/stplanr/
- The package looses its Imports dependency on rgdal, which has been demoted to a Suggests

## BUG FIXES

- An issue with `route_graphhopper()` has been fixed, see https://github.com/ropensci/stplanr/pull/297

# stplanr 0.2.7

## NEW FEATURES

* Various changes to support `dplyr` [0.8.0]: https://github.com/ropensci/stplanr/pull/275


## BUG FIXES

* Fixed [#272](https://github.com/ropensci/stplanr/issues/272) by removing byvars argument of overline in preparation for overdue overhaul of overline function.

## OTHER

* No longer suggests **tmap** to reduce install times: `install.packages()` installs suggested packages by default

# stplanr 0.2.6

## NEW FEATURES

* New function `route_local()`
* New argument in `line2route()`: `time_sleep` waits a period between each route request

## BUG FIXES

* Issue with `dl_stats19()`, see [#270](https://github.com/ropensci/stplanr/issues/270)
* Make style consistent, see [commit](https://github.com/ropensci/stplanr/commit/521598c329c15ff2dfeb159b0da4ba7b1507b060)
* Various small fixes to documentation and style

# stplanr 0.2.5

## NEW FEATURES

* New function `line_via()` for identifying intermediary points on a transport network

## BUG FIXES

* Bug associated with `SpatialLinesNetwork()` fixed (see [#249](https://github.com/ropensci/stplanr/issues/249))

# stplanr 0.2.4

## NEW FEATURES

* New function `geo_length()` returns numeric vector of line lengths from **sp** or **sf** objects.

## DOCUMENTATION

- `?route_graphhopper` no longer mentions the depreciated 'bike2' profile - see [#246](https://github.com/ropensci/stplanr/issues/246)
- `?route_osrm` mentions that the public API only routes for cars - see [#246](https://github.com/ropensci/stplanr/issues/246)
- Updated `introducing-stplanr` vignette to show new function and make more robust

# stplanr 0.2.3

## NEW FEATURES

* **stplanr** now imports **lwgeom**, needed for `sf::st_length()`, used in `SpatialLinesNetwork()`.
* Plotting behaviour updated for `sfNetwork` objects: now only plots the geometry by default.
* Improved documentation for `SpatialLinesNetwork()` and `plot()` for spatial networks.

## BUG FIXES

* Bug in `sum_network_routes()` fixed (see [#240](https://github.com/ropensci/stplanr/issues/240)).

# stplanr 0.2.2

## NEW FEATURES

* In this release **sp** is demoted from a Depends to an Imports, meaning that all its functions will not be attached to your namespace (it will not be loaded) when you run `library(stplanr)`, making it less tied to **sp**. This is a continuation of the work to support **sf** and will make it easier for the package to work with alternative representations of geographic data.

## BUG FIXES

* Bug in `geo_select_aeq.sf()` was fixed by Jakub Nowosad in pull [#238](https://github.com/ropensci/stplanr/pull/238)
* An issue with `od_aggregate.sf()` was fixed making it much faster

# stplanr 0.2.0

## NEW FEATURES

* This is the largest release since the package was created, with dozens of changes to support simple features - see https://github.com/ropensci/stplanr/pull/198 for details.
* Support for **sf**. The package now support the new spatial class system for most functions.
* New function `geo_bb()` supercedes `bb2poly()`. The new function can return polygons, points and matrix objects determined by the `output` argument. It also allows bounding boxes to be extended in metres, and scaled in x and y dimensions.
* `geo_code()` now uses nominatim by default to find locations on the maps.
* New function `od_coords()` takes a wide range of input data types to return a consistent output representing OD data as a data frame of origin and destination coordinates. This is used behind the scenes to make other functions more modular.

## WORK IN PROGRESS

Plans for the next release

* New generic `route()` function for routing. This is more flexible and user-friendly than the existing `line2route()` and `route_*()` functions it enhances.
* Updated function names to make using **stplanr** easier and more intuitive.


# stplanr 0.1.9

## NEW FEATURES

* Dependency cull: we have removed dependencies on foreach and doParallel
* `route_cyclestreet()` now also called (correctly) `route_cyclestreets()`
* New `geo_code()` function replaces dependency on RGoogleMaps

## BUG FIXES

* See issues closed after the last release with this search term: https://github.com/ropensci/stplanr/issues?utf8=%E2%9C%93&q=is%3Aissue%20closed%3A%3E2017-06-01%20
* Bug with `google_dist()` fixed
* Fixed fails due to breaking changes in dplyr

# stplanr 0.1.8

## NEW FEATURES

* New argument `combinations` added to `sum_network_routes()` so it runs quicker - see [pull/177](https://github.com/ropensci/stplanr/pull/177).
* New examples added to `sum_network_routes()`, `weightfield()` and `find_network_nodes()` - see e.g. `example(sum_network_routes)` for details.
* New dataset `l_poly` [added](https://github.com/ropensci/stplanr/commit/7641760fbd6718352ed74142e5c339f6216afea4).
* **stplanr** now has a website! See [ropensci.github.io/stplanr/](https://ropensci.github.io/stplanr/).

## BUG FIXES

* Serious bug with `SpatialLinesNetwork()` [fixed](https://github.com/ropensci/stplanr/pull/186).
* Depreciated `_each()` **dplyr** functions replaced with equivalent `_at` or `_all` functions.

# stplanr 0.1.7

## NEW FEATURES

* There is a new vignette! See [vignettes/stplanr-paper.Rmd](https://github.com/ropensci/stplanr/blob/master/vignettes/stplanr-paper.Rmd) and `vignette("stplanr-paper")` for details.
* The original [`introducing-stplanr`](https://github.com/ropensci/stplanr/blob/master/vignettes/stplanr.Rmd) vignette has been updated. It now provides a more basic introduction for people new to R for spatial and transport data.
* `line2route()` has been refactored to improve error detection and allow `n_processes` arguments. Thanks @nikolai-b. See [pull/151](https://github.com/ropensci/stplanr/pull/151) for details.
* `line_match()` function added, a wrapper around `rgeos::gDistance()`, to find similar routes.
* **RCurl** and **data.table** dependencies have been [removed](https://github.com/ropensci/stplanr/pull/169)
* **leaflet** has been demoted from an import to a suggest. This should reduce install times.
* New functions `od_aggregate()` and `sp_aggregate()` have been [added](https://github.com/ropensci/stplanr/pull/165), to enable OD data to be aggregated to new geographic levels.


## BUG FIXES

* `#141` fixed: `viaroute()` works again.
* [#153](https://github.com/ropensci/stplanr/issues/153) fixed: `bidirectional = TRUE` returns a different result in `line_bearing()` now.

## FUTURE PLANS

* A new branch that uses **sf** is being [tested](https://github.com/ropensci/stplanr/pull/164). We may eventually transition to using simple features classes instead of **sp** classes.

# stplanr 0.1.6

## NEW FEATURES

* `onewayid()` is now a generic function, meaning it can handle spatial and non-spatial data
* New arguments provided for `line2route()` allow you to specify variables to join-by - also has updated and more sensible defaults
* New function `od_id_order()` to put origin-destination ids in order, to identify 2 way duplicates (split out from `onewayid()`)

## BUG FIXES

* See the [issue tracker](https://github.com/ropensci/stplanr/issues?q=is%3Aissue+is%3Aclosed)
* Bug in `route_cyclestreet()` leading `change_elev` and `av_incline` being wrong now fixed
* Bug making variable names with spaces in the id columns failed - now fixed [#138](https://github.com/ropensci/stplanr/issues/138)

stplanr 0.1.5
----------------------------------------------------------------

NEW FEATURES

* New argument destinations added to `od2line()`. See `example(od2line)` for an example.
* New dataset `destinations` for showing how OD matrix with destinations can be converted to spatial data
* New argument `list_output` allows the route information to be saved as a list, allowing `save_raw = TRUE` (which does not return a `Spatial` object) to be passed to the `route_` function.
* tmap dependency removed for faster installs

BUG FIXES

* Bug with `line2route()` (#124) fixed
* Various improvements to documentation

stplanr 0.1.4
----------------------------------------------------------------

NEW FEATURES

* New function `reproject()` is a simple wrapper around `spTransform()` that uses
  `crs_select_aeq()` to convert a spatial object in geographic (lat/lon) coordinates
  into on with projected coordinates, with units of 1 m. This is useful for various
  spatial operations, such as finding the length and area of an object.

* Implement `gprojected()`, a function for performing GIS operations on a temporary, projected, version
  of spatial objects.

* Addition of `line_bearing()` to return the bearing of lines based on start and end points.

* Addition of `angle_diff()` for finding the angular difference between lines: are they roughly parallel or perpendicular?

BUG FIXES

* `line2df()` now works on lines with multiple vertices and is faster.

* Fixes in the examples used to illustrate how `od_dist()` works.

stplanr 0.1.3
----------------------------------------------------------------

NEW FEATURES

* Update to OSRM functions to support API v5.

* New parameter `byvars` in the `overline()` function, to allow disaggregation of results by a grouping variable (see `example(overline)`).

* Faster implementation of `od2line()`: `od2line2()`. Plan is to replace the original if no issues are found with new implementation.

* New function `od2odf()` which converts OD data into a dataframe of origins and destinations (feeds `od2line2()` but also useful as self-standing function).

* New argument `new_proj` in `buff_geo()` allows the results to be exported to any coordinate reference system (CRS).

* New function `gprojected()` generalises concept of `buff_geo()`, building on `crs_select_aeq()` to allow any GIS query to be conducted on a temporary projected version of spatial objects with geographical CRSs.

* New function `od_dist()` can quickly calculate Euclidean distances of OD pairs without converting to spatial objects.

BUG FIXES

* Bug fix in `onewayid()` so it captures all lines.

* Various improvements to documentation and code.

stplanr 0.1.2
----------------------------------------------------------------

NEW FEATURES

* Interface to the Google Distance Matrix `API with dist_google`.

* New transport planning API added, with `route_transportapi_public` (for testing).

* Update to `line2route`, allowing it to accept different routing funtions via the new argument `route_fun` (for testing - tested with `route_fun = route_cyclestreet`).

* New functions for creating origin-destination data frames (`point2odf`) and SpatialLinesDataFrames (`points2flow`).

* Addition of `n_vertices` and `is_linepoint` for identifying the number of vertices in spatial objects and whether the 'line' is really a point.

BUG FIXES

* `line2route` refactored, with 10 fold speed increases on large (1000+) batches of lines.

stplanr 0.1.0
----------------------------------------------------------------

NEW FEATURES

* Addition of new class definition `SpatialLinesNetwork`, methods for `plot`
  and `summary` and functions `calc_network_routes` and `find_network_nodes`
  allowing fast route calculations via igraph and other network analysis
  functions.

* Functions for removing beginning and end of lines: `toptail` and
  `toptailgs`. Helper functions `buff_geo`,
  `crs_select_aeq` and `line2points` added.

* Functionality for reading in the UK's stats19 data: `read_stats19_*`
  functions download, unzip and re-categorise the data.

* `read_table` functions added for reading Australian OD data.

* `decode_gl` added to decode Google polylines and other functions for
  querying and reading data from OSRM services.

* `gtfs2sldf` added to import GTFS routes as SpatialLinesDataFrames.

stplanr 0.0.2
----------------------------------------------------------------

* Published on CRAN
Update to reduce dependencies on igraph

It seems the binary versions of stplanr 0.8.6 did not build successfully.

See https://github.com/ropensci/stplanr/issues/473 for details.

## Test environments

* local R installation, R 4.1.1
* various OSs (GitHub Actions), R 4.1.1
* win-builder (devel): https://win-builder.r-project.org/D13453liE2AE

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
Leaflet-providers
=================
An extension to [Leaflet](https://leafletjs.com/) that contains configurations for various free<sup>[1](#what-is-free)</sup> tile providers.

# Usage
Leaflet-providers [providers](#providers) are refered to with a `provider[.<variant>]`-string. Let's say you want to add the nice [Watercolor](https://maps.stamen.com/#watercolor/) style from Stamen to your map, you pass `Stamen.Watercolor` to the `L.tileLayer.provider`-constructor, which will return a [L.TileLayer](https://leafletjs.com/reference.html#tilelayer) instance for Stamens Watercolor tile layer.

```Javascript
// add Stamen Watercolor to map.
L.tileLayer.provider('Stamen.Watercolor').addTo(map);
```

## Protocol relativity (`https://`-urls)

Leaflet-providers tries to use `https://` if the page uses `https://` and the provider supports it.
You can force the use of `https://` by passing `force_http: true` in the options argument.

## Retina tiles

Some providers have retina tiles for which the URL only needs to be slightly adjusted, e.g. `-----@2x.png`. For this, add the retina option in the URL, e.g. `-----{retina}.png`, and set a retina value in the options, e.g. `retina: '@2x'`. If Leaflet detects a retina screen (`L.Browser.retina`), the retina option passed to the tileLayer is set to the value supplied, otherwise it's replaced by an empty string.

# Providers

Leaflet-providers provides tile layers from different providers, including *OpenStreetMap*, *Stamen*, *Esri* and *OpenWeatherMap*. The full listing of free to use layers can be [previewed](https://leaflet-extras.github.io/leaflet-providers/preview/index.html). The page will show you the name to use with `leaflet-providers.js` and the code to use it without dependencies.

## Providers requiring registration

In addition to the providers you are free<b id="what-is-free">1</b> to use, we support some layers which require registration.

### HERE (formerly Nokia).

In order to use HERE layers, you must [register](https://developer.here.com/). Once registered, you can create an `app_id` and `app_code` which you have to pass to `L.tileLayer.provider` in the options:

```Javascript
L.tileLayer.provider('HERE.terrainDay', {
    app_id: '<insert ID here>',
    app_code: '<insert ID here>'
}).addTo(map);
```

[Available HERE layers](https://leaflet-extras.github.io/leaflet-providers/preview/#filter=HERE)

### Mapbox

In order to use Mapbox maps, you must [register](https://tiles.mapbox.com/signup). You can get map ID and ACCESS_TOKEN from [Mapbox projects](https://www.mapbox.com/projects):
```JavaScript
L.tileLayer.provider('MapBox', {id: 'ID', accessToken: 'ACCESS_TOKEN'}).addTo(map);
```

### Esri/ArcGIS

In order to use ArcGIS maps, you must [register](https://developers.arcgis.com/en/sign-up/) and abide by the [terms of service](https://developers.arcgis.com/en/terms/). No special syntax is required.

[Available Esri layers](https://leaflet-extras.github.io/leaflet-providers/preview/#filter=Esri)

# Attribution

This work was inspired from <https://gist.github.com/1804938>, and originally created by [Stefan Seelmann](https://github.com/seelmann).

### What do we mean by *free*?
<b id="what-is-free">1</b>
We try to maintain leaflet-providers in such a way that you'll be able to use the layers we include without paying money.
This doesn't mean no limits apply, you should always check before using these layers for anything serious.
So you want to add a layer?
=======

Yay! go add it to the leaflet-providers.js as long as it follows the following 
rules:

- Don't violate a providers TOS (if it exists, include a link to it)
- Don't pre-populate api keys with working keys.
- It should be a basic tile source, no exteral libraries etc.
- The owner hasn't asked us to remove it (hasn't happened yet)