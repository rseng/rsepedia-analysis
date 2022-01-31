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
- The owner hasn't asked us to remove it (hasn't happened yet)---
output: github_document
---


<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# stplanr <a href='https://docs.ropensci.org/stplanr/'><img src='man/figures/logo.png' align="right" height=215/></a>

<!-- [![Build Status](https://travis-ci.org/ropensci/stplanr.svg?branch=master)](https://travis-ci.org/ropensci/stplanr) -->
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/stplanr)](https://github.com/r-hub/cranlogs.app)
[![](https://cranlogs.r-pkg.org/badges/grand-total/stplanr)](https://cran.r-project.org/package=stplanr)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/stplanr)](https://cran.r-project.org/package=stplanr)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/10_status.svg)](https://github.com/ropensci/software-review/issues/10)
[![R-CMD-check](https://github.com/ropensci/stplanr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/stplanr/actions)

```{r, echo=FALSE, message=FALSE, warning=FALSE}
library(stplanr)
```

**stplanr** is a package for sustainable transport planning with R.

It provides functions for solving common problems in transport planning and modelling, such as how to best get from point A to point B.
The overall aim is to provide a reproducible, transparent and accessible toolkit to help people better understand transport systems and inform policy, as outlined in a [paper](https://journal.r-project.org/archive/2018/RJ-2018-053/index.html) about the package, and the potential for open source software in transport planning in general, published in the [R Journal](https://journal.r-project.org/).

The initial work on the project was funded by the Department of Transport
([DfT](https://www.gov.uk/government/organisations/department-for-transport))
as part of the development of the Propensity to Cycle Tool
(PCT), a web application to explore current travel patterns and cycling potential at zone, desire line, route and route network levels (see [www.pct.bike](https://www.pct.bike/) and click on a region to try it out).
The basis of the methods underlying the PCT is origin-destination data, which are used to highlight where many short distance trips are being made, and estimate how many could switch to cycling.
The results help identify where cycleways are most needed, an important component of sustainable transport planning infrastructure engineering and policy [design](https://www.icevirtuallibrary.com/doi/abs/10.1680/dfct.63495.001).

See the package vignette (e.g. via `vignette("introducing-stplanr")`) 
or an [academic paper on the Propensity to Cycle Tool (PCT)](https://dx.doi.org/10.5198/jtlu.2016.862)
for more information on how it can be used.
This README provides some basics.

Much of the work supports research undertaken at the Leeds' Institute for Transport Studies ([ITS](https://environment.leeds.ac.uk/transport)) but  **stplanr** should be useful to transport researchers and practitioners needing free, open and reproducible methods for working with geographic data everywhere.

## Key functions

Data frames representing flows between origins and destinations
must be combined with geo-referenced zones or points to generate meaningful
analyses and visualisations of 'flows' or origin-destination (OD) data.
**stplanr** facilitates this with
`od2line()`, which takes flow and geographical data as inputs and
outputs spatial data. Some example data is provided in the package:

```{r, results='hide', message=FALSE}
library(stplanr)
```

Let's take a look at this data:

```{r}
od_data_sample[1:3, 1:3] # typical form of flow data
cents_sf[1:3,] # points representing origins and destinations
```

These datasets can be combined as follows:

```{r plot1, warning=FALSE}
travel_network <- od2line(flow = od_data_sample, zones = cents_sf)
w <- flow$all / max(flow$all) *10
plot(travel_network, lwd = w)
```


**stplanr** has many functions for working with OD data.
See the [`stplanr-od`](https://docs.ropensci.org/stplanr/articles/stplanr-od.html) vignette for details.

The package can also allocate flows to the road network, e.g. with [CycleStreets.net](https://www.cyclestreets.net/api/) and the OpenStreetMap Routing Machine ([OSRM](https://github.com/Project-OSRM/osrm-backend)) API interfaces.
These are supported in `route_*()` functions such as `route_cyclestreets` and `route_osrm()`:

Routing can be done using a range of back-ends and using lat/lon or desire line inputs with the `route()` function, as illustrated by the following commands which calculates the route between Fleet Street and Southwark Street over the River Thames on Blackfriars Bridge in London:

```{r, eval=FALSE, echo=FALSE}
tmaptools::geocode_OSM("fleet street london")
tmaptools::geocode_OSM("southwark street london")
```


```{r}
library(osrm)
trip <- route(
  from = c(-0.11, 51.514),
  to = c(-0.10, 51.506),
  route_fun = osrmRoute,
  returnclass = "sf"
  )
mapview::mapview(trip)
```

You can also use and place names, found using the Google Map API:

```{r cycle-trip, message=FALSE, warning=FALSE, eval=FALSE, echo=FALSE}
if(!Sys.getenv("CYCLESTREETS") == ""){
  trip <- route_cyclestreets("Bradford, UK", "Leeds, Yorkshire", plan = "balanced")
  plot(trip)
}
```

```{r}
trip2 <- route(
  from = "Leeds",
  to = "Bradford",
  route_fun = osrmRoute,
  returnclass = "sf"
  )
mapview::mapview(trip2)
```

We can replicate this call multiple times with the `l` argument in `route()`:

```{r}
desire_lines <- travel_network[2:6, ]
```

```{r, echo=FALSE}
# Sys.sleep(2) # wait a moment 
```


Next, we'll calculate the routes:

```{r plot2, results='hide', message=FALSE}
routes <- route(
  l = desire_lines,
  route_fun = osrmRoute,
  returnclass = "sf"
  )
mapview::mapview(routes) +
  mapview::mapview(desire_lines, color = "red")
```

<!-- The resulting routes will look something like this: -->

```{r routes, echo=FALSE, eval=FALSE}
lwd <- desire_lines$foot
routes <- routes_fast_sf[2:6, ]
plot(routes$geometry, lwd = lwd)
plot(desire_lines$geometry, col = "green", lwd = lwd, add = TRUE)
```

For more examples, `example("route")`.

`overline()` takes a series of route-allocated lines,
splits them into unique segments and aggregates
the values of overlapping lines. This can represent where there will be
most traffic on the transport system, as demonstrated in the following code chunk. 

```{r rnet, warning=FALSE}
routes$foot <- desire_lines$foot
rnet <- overline(routes, attrib = "foot")
```

The resulting route network, with segment totals calculated from overlapping parts for the routes for walking, can be visualised as follows:

```{r}
plot(rnet["foot"], lwd = rnet$foot)
```

The above plot represents the number walking trips made (the 'flow') along particular segments of a transport network.

<!-- (results not shown): -->

```{r, routes-leaf, eval=FALSE, echo=FALSE}
library(leaflet)
pal = leaflet::colorNumeric(palette = "YlGnBu", domain = rnet$all)
leaflet(data = rnet) %>%
  addProviderTiles(providers$OpenStreetMap.BlackAndWhite) %>%
  addPolylines(weight = rnet$all / 3, color = ~pal(all), opacity = 0.9) %>% 
  addLegend(pal = pal, values = ~all)
```

## Policy applications

The examples shown above, based on tiny demonstration datasets, may not seem particularly revolutionary.
At the city scale, however, this type of analysis can be used to inform sustainable transport policies, as described in papers [describing the Propensity to Cycle Tool](https://www.jtlu.org/index.php/jtlu/article/view/862/859) (PCT), and its [application to calculate cycling to school potential](https://www.sciencedirect.com/science/article/pii/S2214140518301257) across England.

Results generated by **stplanr** are now part of national government policy: the PCT is the recommended tool for local and regional authorities developing strategic cycle network under the Cycling and Walking Infrastructure Strategy ([CWIS](https://www.gov.uk/government/publications/cycling-and-walking-investment-strategy)), which is part of the Infrastructure Act [2015](https://www.legislation.gov.uk/ukpga/2015/7/contents/enacted).
**stplanr** is helping dozens of local authorities across the UK to answer the question: where to prioritise investment in cycling?
In essence, stplanr was designed to support sustainable transport policies.

There are many other research and policy questions that functions in **stplanr**, and other open source software libraries and packages, can help answer.
At a time of climate, health and social crises, it is important that technology is not only sustainable itself (e.g. as enabled by open source communities and licenses) but that it contributes to a sustainable future.

## Installation

To install the stable version, use:

```{r, eval=FALSE}
install.packages("stplanr")
```

The development version can be installed using **devtools**:

```{r, eval=FALSE}
# install.packages("devtools") # if not already installed
devtools::install_github("ropensci/stplanr")
library(stplanr)
```

stplanr depends on rgdal, which can be tricky to install.

### Installing stplanr on Linux and Mac

**stplanr** depends on **sf**. Installation instructions for Mac, Ubuntu and other Linux distros can be found here: https://github.com/r-spatial/sf#installing

## Funtions, help and contributing

The current list of available functions can be seen on the package's website at [docs.ropensci.org/stplanr/](https://docs.ropensci.org/stplanr/), or with the following command:

```{r, eval=FALSE}
lsf.str("package:stplanr", all = TRUE)
```

To get internal help on a specific function, use the standard way.

```{r, eval=FALSE}
?od2line
```

To contribute, report bugs or request features, see the [issue tracker](https://github.com/ropensci/stplanr/issues).

```{r, eval=FALSE, echo=FALSE}
# Aim: explore dependencies
desc = read.dcf("DESCRIPTION")
headings = dimnames(desc)[[2]]
fields = which(headings %in% c("Depends", "Imports", "Suggests"))
pkgs = paste(desc[fields], collapse = ", ")
pkgs = gsub("\n", " ", pkgs)
strsplit(pkgs, ",")[[1]]
install.packages("miniCRAN")
library(miniCRAN)
tags <- "stplanr"
pkgDep(tags)
dg <- makeDepGraph(tags, enhances = TRUE)
set.seed(1)
plot(dg, legendPosition = c(-1, 1), vertex.size = 20)
library(DiagrammeR)
DiagrammeR::visnetwork(dg)
visNetwork::visIgraph(dg)
```

## Further resources / tutorials

Want to learn how to use open source software for reproducible sustainable transport planning work?
Now is a great time to learn.
Transport planning is a relatively new field of application in R.
However, there are already some good resources on the topic, including (any further suggestions: welcome):

- The Transport chapter of *Geocomputation with R*, which provides a broad introduction from a geographic data perspective: https://geocompr.robinlovelace.net/transport.html
- The **stplanr** paper, which describes the context in which the package was developed: https://journal.r-project.org/archive/2018/RJ-2018-053/index.html (please cite this if you use **stplanr** in your work)
- The `dodgr` vignette, which provides an introduction to routing in R: https://atfutures.github.io/dodgr/

## Meta

* Please report issues, feature requests and questions to the [github issue tracker](https://github.com/ropensci/stplanr/issues)
* License: MIT
* Get citation information for **stplanr** in R doing `citation(package = 'stplanr')`
* This project is released with a [Contributor Code of Conduct](https://github.com/ropensci/stplanr/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

## Create cents dataset

```{r}
cents <- rgdal::readOGR(dsn = "/home/robin/npct/pct-bigdata/cents.geojson", layer = "OGRGeoJSON")
# library(geojsonio) # load with the ropensci package geojsonio if rgdal fails
# cents <- geojsonio::geojson_read(x = "~/repos/pct/pct-data/national/cents.geojson")
crs <- sp::CRS("+init=epsg:4326")
crsuk <- sp::CRS("+init=epsg:27700")
cents <- sp::spTransform(x = cents, CRSobj = crsuk)
home <- geo_code("LS7 3HB")
home <- sp::SpatialPoints(matrix(home, ncol = 2), proj4string = crs)
home <- sp::spTransform(x = home, CRSobj = crsuk)
buf <- rgeos::gBuffer(home, width = 2000)
# Check it saved the points OK
cents <- cents[buf, ]
plot(buf)
points(cents)
cents <- sp::spTransform(x = cents, CRSobj = crs)
cents$geo_code <- as.character(cents$geo_code)
library(devtools)
# use_data(cents, overwrite = TRUE)
cents_sf <- sf::st_as_sf(cents)
devtools::use_data(cents_sf)

cents <- rgdal::readOGR(dsn = "/home/robin/npct/pct-bigdata/cents.geojson", layer = "OGRGeoJSON")
# library(geojsonio) # load with the ropensci package geojsonio if rgdal fails
# cents <- geojsonio::geojson_read(x = "~/repos/pct/pct-data/national/cents.geojson")
crs <- sp::CRS("+init=epsg:4326")
crsuk <- sp::CRS("+init=epsg:27700")
cents <- sp::spTransform(x = cents, CRSobj = crsuk)
home <- geo_code("LS7 3HB")
home <- sp::SpatialPoints(matrix(home, ncol = 2), proj4string = crs)
home <- sp::spTransform(x = home, CRSobj = crsuk)
buf <- rgeos::gBuffer(home, width = 2000)
# Check it saved the points OK
cents <- cents[buf, ]
plot(buf)
points(cents)
cents <- sp::spTransform(x = cents, CRSobj = crs)
cents$geo_code <- as.character(cents$geo_code)
library(devtools)
# use_data(cents, overwrite = TRUE)
cents_sf <- sf::st_as_sf(cents)
devtools::use_data(cents_sf)
```

```{r, eval=FALSE, echo=FALSE}
# sample od data
remotes::install_github("ITSLeeds/pct")
library(pct)
od_data_all = pct::get_od()
sel_local = 
  od_data_all$geo_code1 %in% cents_sf$geo_code &
  od_data_all$geo_code2 %in% cents_sf$geo_code 
od_data_sample = od_data_all[sel_local, ]
od_data_lines = od2line(od_data_sample, cents_sf)
plot(od_data_lines)
od_data_routes_fast = route(l = od_data_lines, route_fun = cyclestreets::journey, plan = "fastest")
pryr::object_size(od_data_routes_fast) # 500+ kB
names(od_data_routes_fast)
od_data_routes = od_data_routes_fast[c(19:ncol(od_data_routes_fast))]
plot(od_data_routes)
pryr::object_size(od_data_routes) # 300+ kB
usethis::use_data(od_data_routes, overwrite = TRUE)
usethis::use_data(od_data_lines)
file.size("data/od_data_routes.rda") # 24 kB - it's fine
system('ls -hal data/od_data_routes.rda')
```

The code below shows how the `od_lnd` and `c_lnd` datasets were created.

```{r}
library(dplyr)

# get nationwide OD data
od_all = pct::get_od()
nrow(od_all)
#> 2402201
od_all$Active = (od_all$bicycle + od_all$foot) /
    od_all$all * 100
centroids_all = pct::get_centroids_ew() %>% sf::st_transform(4326)
z_london = pct::get_pct_zones(region = "london") %>% 
  select(geo_code, all, foot)
nrow(centroids_all)
#> 7201
london = pct::pct_regions %>% filter(region_name == "london")
centroids_london = centroids_all[london, ]
od_london = od_all %>%
  filter(geo_code1 %in% centroids_london$msoa11cd) %>% 
  filter(geo_code2 %in% centroids_london$msoa11cd)
od_london = od_all[
  od_all$geo_code1 %in% centroids_london$msoa11cd &
  od_all$geo_code2 %in% centroids_london$msoa11cd , 
]
```

```{r}
# aim: create a reproducible OD dataset
od_lnd = od_london %>% 
  select(-matches("rail|name|moto|car|tax|home|la_|Active")) %>% 
  filter(geo_code2 == "E02000001") %>% 
  top_n(4, wt = all)
c_lnd = centroids_london %>% 
  filter(msoa11cd %in% c(od$geo_code1, od$geo_code2))
z_lnd = z_london %>% 
  filter(geo_code %in% c_lnd$msoa11cd)
usethis::use_data(od_lnd)
usethis::use_data(c_lnd)
usethis::use_data(z_lnd)
```

---
title: "Route networks with stplanr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Route networks with stplanr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## This vignette is work in progress - watch this space!

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = all(c(curl::has_internet(), requireNamespace("igraph", quietly = TRUE))) 
)
```

```{r setup, message=FALSE}
library(stplanr)
library(sf)
```

# Introduction

Route networks represent the network of highways, cycleways, footways and other ways along which transport happens.
You can get route network data from OpenStreetMap (e.g. via the `osmdata` R package) and other providers or transport network data.


# Creating route networks from overlapping routes

Unlike routes, each segment geometry in a route network can only appear once.

**stplanr** can be used to convert a series of routes into a route network, using the function `overline()`, as illustrated below:

```{r, out.width="40%", fig.show='hold', fig.width=5, message=FALSE}
library(stplanr)
library(sf)
sample_routes <- routes_fast_sf[2:6, 1]
sample_routes$value <- rep(1:3, length.out = 5)
rnet <- overline(sample_routes, attrib = "value")
plot(sample_routes["value"], lwd = sample_routes$value, main = "Routes")
plot(rnet["value"], lwd = rnet$value, main = "Route network")
```

The above figure shows how `overline()` breaks the routes into segments with the same values and removes overlapping segments.
It is a form of geographic aggregation.

<!-- The figure below shows in more detail how the function works with 2, 3 and then 6 lines (see the vignette's source code to reproduce the plot): -->

```{r rnets1, message=FALSE, warning=FALSE, out.width="100%", fig.width=6, fig.height=6, echo=FALSE}
# knitr::include_graphics("route-networks.png")
```

# Identifying route network groups

Route networks can be represented as a graph.
Usually all segments are connected together, meaning the graph is [connected](https://en.wikipedia.org/wiki/Connectivity_(graph_theory)#Connected_graph).
We can show that very simple network above is connected as follows:

```{r}
touching_list = st_intersects(sample_routes)
g = igraph::graph.adjlist(touching_list)
igraph::is_connected(g)
```

A more complex network may not be connected in this way, as shown in the example below:

```{r, eval=FALSE}
# piggyback::pb_download_url("r_key_roads_test.Rds")
u = "https://github.com/ropensci/stplanr/releases/download/0.6.0/r_key_roads_test.Rds"
rnet_disconnected = readRDS(url(u))
touching_list = sf::st_intersects(rnet_disconnected)
g = igraph::graph.adjlist(touching_list)
igraph::is_connected(g)
#> [1] FALSE
sf:::plot.sfc_LINESTRING(rnet_disconnected$geometry)
```

```{r, echo=FALSE}
knitr::include_graphics("https://i.imgur.com/827f761.png")
```

The elements of the network are clearly divided into groups.
We can identify these groups as follows:

```{r, eval=FALSE}
rnet_disconnected$group = rnet_igroup(rnet_disconnected)
```

# SpatialLineNetworks

An important feature of route networks is that they are simultaneously spatial and graph entities.
This duality is captured in `sfNetwork` objects, which can be created by the function `SpatialLinesNetwork()`:

```{r rnet-routing1}
sln <- SpatialLinesNetwork(rnet)
class(sln)
```

`sln` has both spatial and graph components, with the number of lines equal to the number graph edges:

```{r}
class(sln@sl)
nrow(sln@sl)
class(sln@g)
length(igraph::edge.attributes(sln@g)[["weight"]])
class(sln@nb)
length(unique(unlist(sln@nb)))
identical(sln@sl$geometry, rnet$geometry)
```

```{r}
sln_nodes <- sln2points(sln)
nrow(sln_nodes)
length(sln@nb)
```

```{r}
rnet_coordinates <- sf::st_coordinates(rnet)
set.seed(85)
x <- runif(n = 2, min = min(rnet_coordinates[, 1]), max = max(rnet_coordinates[, 1]))
y <- runif(n = 2, min = min(rnet_coordinates[, 2]), max = max(rnet_coordinates[, 2]))
crs <- sf::st_crs(rnet)
xy_sf <- sf::st_as_sf(data.frame(n = 1:2, x, y), coords = c("x", "y"), crs = crs)
xy_nodes <- stplanr::find_network_nodes(sln = sln, x = x, y = y)
```

# Routing on route networks

Currently not running due to issues with dev version of `dplyr`:

https://github.com/ropensci/stplanr/issues/383

```{r, out.width="49%", fig.show='hide'}
# plot(rnet$geometry)
# plot(sln_nodes, add = TRUE)
# xy_path <- sum_network_routes(sln = sln, start = xy_nodes[1], end = xy_nodes[2], sumvars = "length")
# # xy_path = sum_network_links(sln = sln, start = xy_nodes[1], end = xy_nodes[2])
# plot(rnet$geometry)
# plot(xy_sf$geometry, add = TRUE)
# plot(xy_path$geometry, add = TRUE, lwd = 5)
```

# Adding new nodes

New nodes can be added to the network, although this should be done before the graph representation is created.
Imagine we want to create a point half way along the the most westerly route segment in the network, near the coordinates -1.540, 53.826:

```{r netpoint}
new_point_coordinates <- c(-1.540, 53.826)
p <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(new_point_coordinates)), crs = crs)
```

We can identify the nearest point on the network at this point and use that to split the associated linestring:

```{r, fig.show='hold', out.width="49%"}
sln_new <- sln_add_node(sln = sln, p = p)
route_new <- route_local(sln = sln_new, from = p, to = xy_sf[1, ])
plot(sln_new)
plot(p, add = TRUE)
plot(route_new, lwd = 5, add = TRUE)
```

# Other approaches

Other approaches to working with route networks include:

- [sDNA](https://github.com/fiftysevendegreesofrad/sdna_open), an open source C++ library for analysing route networks and estimating flows at segments across network segments
- [sfnetworks](https://github.com/luukvdmeer/sfnetworks), an R package that provides an alternative igraph/sf spatial network class
- [dodgr](https://github.com/ATFutures/dodgr), an R package providing functions for calculating distances on directed graphs
- [cppRouting](https://cran.r-project.org/package=cppRouting), a package for routing in C++
- Chapter [10 of Geocomputation with R](https://geocompr.robinlovelace.net/transport.html), which provides context and demonstrates a transport planning workflow in R.

```{r, message=FALSE, warning=FALSE, out.width="100%", fig.width=6, fig.height=6, echo=FALSE, eval=FALSE}
# Show solutions to https://github.com/ropensci/stplanr/issues/237

library(tmap)
# sample_routes = routes_fast_sf
# sample_routes$id <- 1:nrow(sample_routes)
# sample_routes <- sample_routes[!is.na(sample_routes$plan),]
sample_routes <- routes_fast_sf[c(22, 38, 39, 46, 47), NULL]
sample_routes$type <- " Routes"
v <- 1:5
n <- c(2, 3, 5)
# route_list = purrr::map(n, ~sample_routes[1:., ])
route_list <- lapply(n, function(x) {
  l <- sample_routes[1:x, ]
  l$n <- x
  l$value <- rep(v, length.out = x)
  l
})
routes_all <- do.call(rbind, route_list)
# rnet_list = purrr::map(route_list, function(x) {
rnet_list <- lapply(route_list, function(x) {
  l <- overline2(x, "value")
  l$n <- mean(x$n)
  l
})
rnet_all <- do.call(rbind, rnet_list)
rnet_all$type <- "Route network"
all_routes <- rbind(routes_all, rnet_all)
p <- sf::st_centroid(all_routes)
tmap_mode("plot")
m <- tm_shape(all_routes, bbox = tmaptools::bb(all_routes)) +
  tm_lines(
    col = "value", lwd = 1, palette = "OrRd", scale = 8, alpha = 0.8, breaks = 0:6,
    legend.lwd.show = FALSE, labels = as.character(1:6),
    legend.col.show = FALSE
  ) +
  tm_text("value") +
  tm_facets(by = c("n", "type")) +
  tm_layout(scale = 1.5)
# m
tmap_save(m, "vignettes/route-networks.png")
```

```{r, eval=FALSE, echo=FALSE}
# test code:

sample_routes2 <- sample_routes5[2:3, ]
sample_routes3 <- sample_routes5[2:4, ]
rnet2 <- overline2(sample_routes2, attrib = "value")
rnet3 <- overline2(sample_routes3, attrib = "value")
rnet5 <- overline2(sample_routes5, attrib = "value")

b <- 0:6
bb <- tmaptools::bb(rnet, ext = 1.1)

rnet5$n <- 5
rnet5$type <- "Route network"
sample_routes5$n <- 5
sample_routes5$type <- " Routes"

all_routes <- rbind(rnet5, sample_routes5)

m2 <- tm_shape(sample_routes[1:2, ], bbox = bb) +
  tm_lines(col = "value", lwd = "value", palette = "magma", scale = 8, alpha = 0.5, breaks = b) +
  tm_layout(title = "2 Routes", legend.show = FALSE)
r2 <- tm_shape(rnet2, bbox = bb) +
  tm_lines(
    col = "value", palette = "viridis", scale = 10, alpha = 0.5, breaks = b,
    legend.lwd.show = FALSE, labels = as.character(1:6)
  )
m3 <- tm_shape(sample_routes[1:3, ], bbox = bb) +
  tm_lines(col = "value", lwd = "value", palette = "viridis", scale = 15, alpha = 0.5, breaks = b) +
  tm_layout(title = "3 Routes", legend.show = FALSE)
r3 <- tm_shape(rnet3, bbox = bb) +
  tm_lines(col = "value", palette = "viridis", scale = 8, alpha = 0.5, breaks = b) +
  tm_layout(legend.show = FALSE)
m6 <- tm_shape(sample_routes, bbox = bb) +
  tm_lines(col = "value", lwd = "value", palette = "viridis", scale = 10, alpha = 0.5, breaks = b) +
  tm_layout(title = "6 Routes", legend.show = FALSE)
r6 <- tm_shape(rnet, bbox = bb) +
  tm_lines(col = "value", palette = "viridis", scale = 8, alpha = 0.5, breaks = b) +
  tm_layout(legend.show = FALSE)



tmap_arrange(m2, r2, m3, r3, m6, r6, nrow = 3)
```

---
title: "Parallel routing and performance with stplanr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Parallel routing and performance with stplanr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

The code in this vignette provides a basic introduction to parallel routing.
Routing is something that is highly parallelisable because each route can be calculated independently of the others.
The code should be fairly self-explanatory.
No results are shown and the code is not run to reduce package build times.

```{r setup}
library(stplanr)
library(sf)
library(dplyr)
library(tmap)
library(parallel)
library(cyclestreets)
```

# With old route_cyclestreets function

```{r}
# ?route
l = flowlines_sf %>% 
  dplyr::filter()
t1 = Sys.time()
routes_route_cyclestreet = line2route(l)
Sys.time() - t1
ncol(routes_route_cyclestreet)
nrow(routes_route_cyclestreet)
names(routes_route_cyclestreet)
routes_route_cyclestreet_joined = dplyr::inner_join(routes_route_cyclestreet, sf::st_drop_geometry(l))
Sys.time() - t1
rnet_go_dutch = overline(routes_route_cyclestreet_joined, "All")
Sys.time() - t1
tm_shape(rnet_go_dutch) +
  tm_lines(lwd = 5, col = "All", breaks = c(0, 10, 100, 500, 1000), palette = "viridis")
```

# With new route function

```{r}
# ?route
t1 = Sys.time()
routes_journey = route(l = l, route_fun = cyclestreets::journey)
ncol(routes_journey)
nrow(routes_journey)

Sys.time() - t1
names(routes_journey)
rnet_go_dutch_journey = overline(routes_journey, "All")
Sys.time() - t1
rnet_go_dutch_agg = overline(routes_journey, "All")
Sys.time() - t1
tm_shape(rnet_go_dutch_agg) +
  tm_lines(lwd = 5, col = "All", breaks = c(0, 10, 100, 500, 1000), palette = "viridis")
```

# With new route function in parallel

```{r}
# ?route
t1 = Sys.time()


# load parallel stuff
cl <- makeCluster(detectCores())
clusterExport(cl, c("journey"))
Sys.time() - t1
routes_journey_par = route(l = l, route_fun = cyclestreets::journey, cl = cl) # multi-core
stopCluster(cl) # kill cluster

Sys.time() - t1
Sys.time() - t1
names(routes_journey_par)
rnet_go_dutch_journey = overline(routes_journey_par, "All")
Sys.time() - t1
rnet_go_dutch_agg = overline(routes_journey_par, "All")
Sys.time() - t1
tm_shape(rnet_go_dutch_agg) +
  tm_lines(lwd = 5, col = "All", breaks = c(0, 10, 100, 500, 1000), palette = "viridis")
```

# In parallel with quietness plan


```{r}
# ?route
t1 = Sys.time()


# load parallel stuff
library(parallel)
library(cyclestreets)
cl <- makeCluster(detectCores())
clusterExport(cl, c("journey"))
Sys.time() - t1
routes_journey_par = route(l = l, route_fun = cyclestreets::journey, cl = cl, plan = "quietest") # multi-core
stopCluster(cl) # kill cluster

Sys.time() - t1
Sys.time() - t1
names(routes_journey_par)
rnet_go_dutch_journey = overline(routes_journey_par, "All")
Sys.time() - t1
rnet_go_dutch_agg = overline(routes_journey_par, "All")
Sys.time() - t1
tm_shape(rnet_go_dutch_agg) +
  tm_lines(lwd = 5, col = "All", breaks = c(0, 10, 100, 500, 1000), palette = "viridis")
```

## Tests

```{r}
routes_journey_aggregated = routes_journey %>% # already has data from data frame in there!
  group_by(id) %>% 
  summarise(All = median(All)) %>% 
  sf::st_cast("LINESTRING")


rnet_journey_dplyr = routes_journey %>% # already has data from data frame in there!
  group_by(name, distances) %>% 
  summarise(All = sum(All)) 
Sys.time() - t1
tm_shape(rnet_journey_dplyr) +
  tm_lines(lwd = 5, col = "All", breaks = c(0, 10, 100, 500, 1000), palette = "viridis") # quite different...


rnet_journey_go_dutch = routes_journey %>% 
  group_by(start_longitude, start_latitude, finish_longitude, finish_latitude) %>% 
  summarise(All = sum(All))


```

---
title: "Introducing stplanr"
author: "Robin Lovelace"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introducing stplanr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography:
  - references.bib
  - stplanr-citation.bib
---

```{r, include=FALSE}
library(stplanr)
knitr::opts_chunk$set(eval = curl::has_internet())
```

# Introduction

The purpose of this vignette is to get you up-to-speed with the basics and provide useful links for doing transport research with R.

**stplanr** was initially developed to answer a practical question: how to convert official data on travel behaviour into geographic objects that can be plotted on a map and analysed using methods from geographical information systems (GIS)?
Specifically, how can origin-destination (OD) data, such as the open datasets provided by the UK Data Services WICID portal (see [wicid.ukdataservice.ac.uk/](https://wicid.ukdataservice.ac.uk/)), be used to estimate cycling potential down to the street levels at city and national levels?
The project was initially developed to support the Propensity to Cycle Tool (PCT), which has now been deployed as a national web application hosted at www.pct.bike and written-up as an academic paper [@lovelace_propensity_2017].

**stplanr** has since grown to include a wide range of functions for transport planning.
The package was [reviewed](https://github.com/ropensci/software-review/issues/10) through the rOpenSci package review process and the package is now hosted on their site. See the website at [docs.ropensci.org/stplanr](https://docs.ropensci.org/stplanr/).
A more detailed overview of the package's aims and capabilities is contained in a [longer vignette](https://github.com/ropensci/stplanr/blob/master/vignettes/stplanr-paper.Rmd), which has since been published in the R Journal [@lovelace_stplanr_2018].

# Installing stplanr

If you're new to programming and transport data, we recommend using **stplanr** interactively in an Integrated Development Environment (IDE), [RStudio](https://www.rstudio.com/products/rstudio/download/).
Broader guidance on R set-up can be found in [Efficient R Programming](https://csgillespie.github.io/efficientR/set-up.html) [@gillespie_efficient_2016], RStudio's [Education pages](https://education.rstudio.com/learn/beginner/) and on [CRAN](https://cran.r-project.org/).

Once you have an R set-up you are happy with, the latest version can be installed as follows:

```{r, eval=FALSE}
install.packages("stplanr")
```

To install the development version, which may have new features, can be installed as follows:

```{r, eval=FALSE}
remotes::install_github("ropensci/stplanr")
```

Load the package as follows:

```{r}
library(stplanr)
```

**stplanr** contains many datasets for testing and demonstrating how R can be used for transport planning.
The names of these datasets (which are loaded 'lazily' into your namespace when you attach **stplanr**) are listed below:

```{r}
data(package = "stplanr")$result[, "Item"]
```

A more complete list of functions in the package can be found here: https://docs.ropensci.org/stplanr/reference/index.html.

# OD data to desire lines and routes

Transport data can take many forms.
R is an appropriate language for handling transport data, as it can read-in data in such a wide range of formats, e.g. with packages such as **haven** and **foreign**.
This section focusses on OD datasets, and their conversion to *desire lines* and *routes* because these are foundational data types for many transport research applications. (**stplanr** also contains functions for: the analysis of road traffic casualty data, interfacing with various routing APIs, 'travel watershed' analyis and access to Google's Travel Matrix API.)

Origin-destination (OD) data is simply data in the following form:

```{r}
od_eg <- read.csv(
  text =
  "origin, destination, V1, V2
  1, 2, 100, 3
  1, 3, 50, 5"
)
knitr::kable(od_eg)
```

What this example OD table means is that 100 units of 'V1' and 3 units of V2 travel between zone 1 and zone 2. There is also movement represented between Zone 1 and 3. 

This dataset can also be represent as an 'od matrix', where rows represent the origins and columns destinations. However, for multiple variables (e.g. modes of transport) and to prevent giant and unwieldy sparse matrices, the 'long' form represented above is much more common.

Now, imagine that V1 represents the total number of people travelling between the origin and destination and that V2 represents the number who regularly cycle. From this we can get a good indication of where people cycle at the desire line level. (Note: a good source of open OD data has been made available from the [wicid.ukdataservice.ac.uk](https://wicid.ukdataservice.ac.uk/) website).

To extract useful information from this OD dataset, we need to be able to place the lines on the map. What kind of place does a desire line originate from? What about the destination? What is the environment like that it passes through? To answer all these questions we need a geographic representation of the OD table illustrated above.

# Converting OD data to desire lines with R

One problem with OD data is that the rows do not tend to have geography inherently built in. They could contain a variables called `lat_origin`, `lon_origin`, `lat_destination` and `lon_destination`. But generally they only contain the IDs of geographic zones.

Work is needed to convert the OD data into 'desire lines'. Desire lines are straight lines between the origin and destination and represent where people would go if they were not constrained by the route network (see Figure 3 from [@lovelace_propensity_2017]).

To show how these desire lines are created, we'll switch to using real OD data provided by **stplanr**. The first three of these is shown below:

```{r}
head(flow[c(1:3, 12)])
```

This shows that, between zone E02002361 and E02002361 (i.e. intrazonal flow) there were 109 people travelling to work by all modes in the 2011 census. 2 of them cycled. The equivalent numbers for the OD pair E02002361 to E02002371 were 44 and 3. But how to make this data geographical?

For that we need another dataset, also provided by **stplanr**:

```{r}
head(cents_sf)
```

The `cents_sf` dataset is *spatial*, as defined in the `sf` package.
The default `plot()` method for `sf` objects creates a map, as illustrated below:

```{r}
library(sf)
class(cents_sf)
plot(cents_sf)
```

**stplanr** creates desire lines using the `od2line()` function, which links geographical and non-geographical datasets together. 
Note: this functionality has been superseded by functions in the `od` package.
In this case, it will join the non-geographical `flow` data with the geographical `cents_sf` data plotted above. Let's take a single OD pair, E02002361 to E02002371, the fourth row represented in the table above, to see how this works:

```{r}
flow_single_line <- od_data_sample[2:3, ] # select only the first line
desire_line_single <- od2line(flow = flow_single_line, zones = cents_sf)
```

This can be plotted as follows:

```{r}
plot(desire_line_single$geometry, lwd = 5)
plot(cents_sf, add = TRUE, cex = 5)
```

The following command creates desire lines longer than than 2km in distance via the `geo_length()` function --- omitting 'internal flows' via the `sel` object below --- represented in the dataset `flowlines`:

```{r}
l <- od2line(flow = flow, zones = cents_sf)
# identify 'intrazone flows'
sel_intra <- l$Area.of.residence == l$Area.of.workplace
# find distances
l_distances <- geo_length(l)
summary(l_distances)
sel_dist <- l_distances > 2000
sel <- !sel_intra & sel_dist
l <- l[sel, ]
```

This creates the geographic data object `l`, which can be visualised as follows:

```{r, eval=FALSE}
plot(l)
```

Now the data is set-up, we can change the visual appearance of the desire lines with a single extra argument passed to the plotting function. Let's make width depend on the total number of people travelling along the desire line:

```{r, echo=FALSE}
l_bb <- sf::st_bbox(l)
# l_bb[1] <- NA
no_na_in_bb <- !any(is.na(as.numeric(l_bb)))
knitr::opts_chunk$set(eval = no_na_in_bb)
```

```{r}
lwd <- l$All / mean(l$All)
plot(st_geometry(l), lwd = lwd)
```

Another useful visualisation involves setting the colour relative to the number of people cycling:

```{r}
plot(l["Bicycle"], lwd = lwd)
```

Finally, we can convert these desire lines into routes as follows (other routing functions can be used, but may require API keys to work - see the [`cyclestreets`](https://cran.r-project.org/package=cyclestreets) package documentation for example):

```{r, message=FALSE, warning=FALSE}
# if the next line returns FALSE the code will not run
(has_internet <- curl::has_internet())
(cs_key <- nchar(Sys.getenv("CYCLESTREETS")))
if (has_internet & cs_key == 16) {
  r <- route(l = l, route_fun = cyclestreets::journey)
  r <- aggregate(r[c(3, 12)], by = list(r[[1]], r[[2]]), FUN = mean)
} else {
  r <- routes_fast_sf[sel, ]
}
```

These routes contain the same information on origin and destination, but have additional spatial information about the route network.
The routes can be plotted in the same way as the desire lines were plotted:

```{r, out.width="500", out.height="500", eval=FALSE}
plot(r$geometry, lwd = lwd * 3, reset = FALSE)
```

```{r, out.width="500", out.height="500", echo=FALSE, eval=FALSE}
# alternative showing buildings:
r_sf <- st_sf(l, geometry = st_as_sfc(r))
if (require(osmdata)) {
  buildings <- opq(st_bbox(l)) %>%
    add_osm_feature(key = "building", value = "industrial") %>%
    osmdata_sf()
}
plot(r_sf["Bicycle"], lwd = lwd * 3, reset = FALSE)
plot(st_geometry(buildings$osm_polygons), col = "grey", add = TRUE)
```

The next stage is to aggregate these lines together to create a 'route network'.
This, and many other functions, are described in the [stplanr-paper vignette](https://github.com/ropensci/stplanr/blob/master/vignettes/stplanr-paper.Rmd).

# Motivations

As settlements worldwide have grown and become more complex, the process of planning has had to adapt. Planners today are specialists, in sub-fields such as Emergency, Logistics, Healthcare, Urban and Transport Planning.
The 'art' of planning has become more of a science, with its own array of specialist hardware and software.

Like other types of planning, new technologies are changing and in many ways improving the practice of Transport Planning.
Transport interventions such as new bridges, ports and active travel routes are no longer only decided based on the intuition of public sector or political authorities.
Decisions are now the result of a long socio-technical process involving public consultation, cost-benefit analyses and computer modeling and visualisation.
With the ongoing digital revolution, the importance of this last stage has grown, to the point where transport planning is now a highly technical process, employing dozens of software developers in large planning organizations.
There is now a multi-billion pound global transport planning consultancy industry, to support the decision-making process.
Yet the results of all this labor are unavailable to the vast majority of citizens worldwide.
Transport planning decisions which go against the best available evidence keep getting made.

In this context the aim of **stplanr** is to provide an accessible toolbox for transport planning, with a focus on geographic data.
It is hoped that it will be useful for practitioners and researchers alike, as part of the ongoing transition to open source software taking place in the tech industry.

A further motivation is that the best available [evidence](https://www.nature.com/articles/nclimate2923) suggests the future of civilization depends on our ability to transition away from fossil fuels.
The transport sector is the fastest growing source of emissions by sector, and represents a major roadblock in the path towards a zero-carbon economy. Transport systems are also a major cause of ill health, by enabling sedentary lifestyles and causing numerous road traffic casualties.
Knowledge of these impacts motivated the word 'sustainable' in the package's name: by focusing on active travel and public transport modes, **stplanr** is intended to encourage interventions that reduce dependence on fossil fuels.

# Further resources

**stplanr** is focussed on geographic data.
The reason for this is that almost all transport data, from the spatial distribution of bus stops to the routes that pedestrians take between home and work, contains a spatial element.
Representing this spatial data in a formal class system has many advantages, including sensible defaults for plotting the spatial data on a map and support for a range of geographic operations.

**sf** supports most common geographic data formats used in transport planning (including Shapefiles and GeoJSON files representing points, lines, zones). 
See *stplanr: A package for transport planning* [@lovelace_stplanr_2018] for details.

To get the best out of **stplanr** it helps to have a strong understanding of spatial data in R in general. Chapter 2 of the open source book [*Geocomputation with R*](https://geocompr.robinlovelace.net/) provides an introductory tutorial on the basics of spatial data with R and contains references to more advanced tutorials which may come in handy as your spatial data analysis skills progress.
Further information on geographic data for transport applications can be found in the same book.
See https://geocompr.robinlovelace.net/transport.html.

# Contributing

We welcome your contributions, whether it's filing a bug or feature request in the [issue tracker](https://github.com/ropensci/stplanr/issues), putting in a pull request to improve performance or documentation, or simply letting us know how you're using **stplanr** in your work by citing it or dropping us an email.

# References

---
title: "Origin-destination data with stplanr"
output: rmarkdown::html_vignette
author: "Robin Lovelace and Edward Leigh"
vignette: >
  %\VignetteIndexEntry{Origin-destination data with stplanr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography:
 - references.bib
 - stplanr-citation.bib
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
has_webshot <- "webshot" %in% installed.packages()
```

Note: an updated version of this vignette, called 'od', is available in the `od` package.
View it as follows:

```{r, eval=FALSE}
install.packages("od")
vignette("od")
```


# Introduction: what is OD data?

As the name suggests, origin-destination (OD) data represents movement through geographic space, from an origin (O) to a destination (D).
Sometimes also called '[flow data](https://www.ons.gov.uk/census/2011census/2011censusdata/originanddestinationdata)', OD datasets contain details of trips between two geographic points or, more commonly, zones (which are often represented by a zone centroid).
Most OD datasets refer to start and end locations with 'ID' columns containing character strings such as `zone1`.
These IDs refer to a geographic feature in a separate geographic dataset.
Origin and destination locations are sometimes represented as geographic coordinates.

OD datasets typically contain multiple non geographic attributes.
These usually include, at a minimum, the number of trips that take place from the origin to the destination over a given time period (e.g. a typical work day).
Additional attributes can include breakdown by the mode(s) of transport used for the trips.
Usually only a single mode is captured (trips made by a combination of cycle-train-walk modes are often counted only as 'train' trips).
Additional disaggregations of overall counts may include trip counts at different time periods.

Many OD datasets omit information.
If there is only one time period, then this resides in the metadata for the whole data set.
There is rarely any information about the path taken between the start and end points.
It is typically the job of the analyst to use a routing service (such as [OSRM](https://github.com/riatelab/osrm), [Google Directions API](https://symbolixau.github.io/googleway/articles/googleway-vignette.html#google-directions-api), [CycleStreets.net](https://github.com/Robinlovelace/cyclestreets/) or [OpenRouteService](https://github.com/GIScience/openrouteservice-r/)) or an assignment model (such as those contained in proprietary software such as [SATURN](https://saturnsoftware2.co.uk/) and [Visum](https://www.ptvgroup.com/en/solutions/products/ptv-visum/)) to identify likely routes with reference to shortest path algorithms or generalised cost minimisation algorithms (which account for monetary plus time and quality 'costs').

# The importance of OD data

Despite the rather dull name, OD datasets are a vital part of the modern world: they underpin analysis and models that influence current *and future* transport systems.
Historically, these models, and the OD datasets that drove them, were used to plan for car-dominated cities [@boyce_forecasting_2015].
Now that there is growing evidence of the negative impacts car domination, however, there is a strong argument for transport models being re-purposed.
Origin-destination data can be part of the solution.

From a health perspective transport planning, supported by OD data and analysed primarily using proprietary software and opaque methods, has failed: roads are now the largest cause of death of young people worldwide, killing more than [1 million](https://www.who.int/publications/i/item/9789241565684) people each year [@world_health_organization_global_2018].
Even ignoring problems such as air pollution, obesity and climate change, it is clear that current transport systems are unsustainable.
There are other reasons why transport data analysis and software are important [@lovelace_stplanr_2018].

The purpose of this vignette is to introduce OD data, an important component of many transport planning models, with examples based on data and functions from the stplanr package.
The aim is to enable you to use OD data to inform more sustainable transport plans, for example by identifying 'desire lines' along which policies could cause a modal switch away from cars and towards lower energy modes such as walking, cycling, and public transport.

# An example OD dataset

OD data can be accessed from a range of sources (we will see code that downloads many thousands of OD pairs later in this vignette).
Some 'data carpentry' may be needed before the OD data is ready for analysis.
This vignette does not cover cleaning OD data: we assume you know R and come with 'tidy' data [@wickham_tidy_2014], in which each row represents travel between an origin and a destination (typically zones represented by zone IDs), and each column represents an attribute such as number of trips or vehichle counts by mode or straight line distance.^[
It may be difficult to convert between 'number of trip' and 'number of vehicle' counts for modes in which a single vehicle can contain many people, such as cars (a problem that can be overcome when surveys differentiate between car driver and 'car passenger' as distinct modes), buses and trams if occupancy levels are unknown.
Typically OD data only report single stage trips, but multi-modal trips such as walk-rail-cycle can be represented when such a combination of modes is represented by a new, unique, mode.
]

In simple terms OD data looks like this:

```{r setup, message=FALSE}
library(stplanr)
library(dplyr)
od <- stplanr::od_data_sample %>%
  select(-matches("rail|name|moto|car|tax|home|la_")) %>%
  top_n(n = 14, wt = all)
class(od)
od
```

<!-- The next section discusses this, and other, representations of OD data. -->

<!-- # Representations of OD data -->

Like all data, the object `od`, created in the preceding code chunk, comes from a specific context: the 2011 [UK Census](https://ukdataservice.ac.uk/learning-hub/census/) questions:

- In your main job, what is the address of your workplace?
- How do you usually travel to work (for the longest part, by distance, of your usual journey to work)?
  - Work mainly at or from home
  - Underground, metro, light rail, tram
  - Train
  - ...
  
The object `od` is a data frame containing aggregated answers to these questions (see `?pct::get_od()` for details).
It is *implicitly geographic*: the first two columns refer to geographic entities but do not contain coordinates themselves (OD coordinates are covered below).
Other columns contain attributes associated with each OD pair, typically counting how many people travel by mode of transport.
OD data can be represented in a number of ways, as outlined in the next sections.

# Origin-destination pairs (long form)

The most useful way of representing OD data is the 'long' data frame format described above.
This is increasingly the format used by official statistical agencies, including the UK's Office for National Statistics (ONS), who provide origin destination data as a `.csv` file.
Typically, the first column is the zone code of origin and the second column is the zone code of the destination, as is the case with the object `od`.
Subsequent columns contain attributes such as `all`, meaning trips by all modes, as illustrated below (we will see a matrix representation of this subset of the data in the next section):

```{r}
od[1:3]
```

`geo_code1` refers to the origin, `geo_code2` refers to the destination.

Additional columns can represent addition attributes, such as number of trips by time, mode of travel, type of person, or trip purpose.
The `od` dataset contains column names representing mode of travel (train, bus, bicycle etc), as can be seen with `names(od[-(1:2)])`.
These 'mode' columns contain integers in the example data, but contain characters, dates and other data types, taking advantage of the flexibility of data frames.

# Origin destination matrices

The 'OD matrix' representation of OD data represents each attribute column in the long form as a separate matrix.
Instead of rows representing OD pairs, rows represent all travel from each origin to all destinations (represented as columns).
The **stplanr** function `od_to_odmatrix()` converts between the 'long' to the 'matrix' form on a per column basis, as illustrated below:

```{r}
od_matrix <- od_to_odmatrix(od[1:3])
class(od_matrix)
od_matrix
```

Note that row and column names are now zone codes.
The cell in row 1 and column 2 (`od_matrix[1, 2]`), for example, reports that there are 94 trips from zone `E02002361` to zone `E02002393`.
In the case above, no people travel between the majority of the OD pair combinations, as represented by the `NA`s.
OD matrices are a relatively rudimentary data structure that pre-date R's `data.frame` class.
Typically, they only contained integer counts, providing small and simple datasets that could be used in 20^th^ Century transport modelling software running on limited 20^th^ Century hardware. 

Although 'OD matrix' is still sometimes used informally to refer to any OD datadset, the long OD pair representation is recommended:
OD matrices become unwieldy for large OD datasets, which are likely to be sparse, with many empty cells represented by NAs.
Furthermore, to represent many attributes in matix format, multiple lists of OD matrices or 'OD arrays' must be created.
This is demonstrated in the code chunk below, which represents travel between OD pairs by all modes and by bike:

```{r}
lapply(c("all", "bicycle"), function(x) od_to_odmatrix(od[c("geo_code1", "geo_code2", x)]))
```

The function `odmatrix_to_od()` can converts OD matrices back into the more convenient long form:

```{r}
odmatrix_to_od(od_matrix)
```

# Inter and intra-zonal flows

A common, and sometimes problematic, feature of OD data is 'intra-zonal flows'.
These are trips that start and end in the same zone.
The proportion of travel that is intra-zonal depends largely on the size of the zones used.
It is often useful to separate intra-zonal and inter-zonal flows at the outset, as demonstrated below:

```{r}
(od_inter <- od %>% filter(geo_code1 != geo_code2))
(od_intra <- od %>% filter(geo_code1 == geo_code2))
```

Intra-zonal OD pairs represent short trips (up to the size of the zone within which the trips take place) so are sometimes ignored in OD data analyis.
However, intra-zonal flows can be valuable, for example in measuring the amount of localised transport activity and as a sign of local economies.

# Oneway lines

Another subtly with some ([symetric](https://icaci.org/files/documents/ICC_proceedings/ICC2013/_extendedAbstract/393_proceeding.pdf), where origins and destinations can be the same points) OD data is that oneway flows can hide the extent of bidirectional flows in plots and other types of analysis.
This is illustrated below for a sample of the `od` dataset:

```{r}
(od_min <- od_data_sample[c(1, 2, 9), 1:6])
(od_oneway <- od_oneway(od_min))
```

Note that in the second dataset there are only 2 rows instead of 3.
The function `od_oneway()` aggregates oneway lines to produce bidirectional flows.
By default, it returns the sum of each numeric column for each bidirectional origin-destination pair.

# Desire lines

The previous representations of OD data are all implicitly geographic: their coordinates are not contained in the data, but associated with another object that *is* geographic, typically a zone or a zone centroid.
This is problematic, meaning that multiple objects or files are required to fully represent the same data.
Desire line representations overcome this issue.
They are geographic lines between origin and destination, with the same attributes as in the 'long' representation.

`od2line()` can convert long form OD data to desire lines.
The second argument is a zone or a centroid dataset that contains 'zone IDs' that match the IDs in the first and second columns of the OD data, as illustrated below:

```{r}
z <- zones_sf
class(z)
l <- od2line(flow = od_inter, zones = z)
```

The preceding code chunk created a zones object called `z`, the coordinates of which were used to convert the object `od` into `l`, which are geographic desire lines.
The desire line object is stored in as a geographic simple features object, which has the same number of rows as does the object `od` and one more column:

```{r}
class(l)
nrow(od) - nrow(l)
ncol(l) - ncol(od)
```

The new column is the geometry column, which can be plotted as follows:

```{r}
plot(l$geometry)
```

By default, plotting `l` shows the attributes for each line:

```{r}
plot(l)
```

Because these lines have a coordinate reference system (CRS) inherited from the zones data, they can also be plotted on an interactive map, as follows (result only shown if webshot is installed):

```{r}
library(leaflet)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = l)
```

## Non-matching IDs

Note that in some OD datasets there may be IDs that match no zone.
We can simulate this situation by setting the third origin ID of `od` to `nomatch`, a string that is not in the zones ID:

```{r, error=TRUE}
od$geo_code2[3] <- "nomatch"
od2line(od, z)
```

You should clean your OD data and ensure all ids in the first two columns match the ids in the first column of the zone data before running `od2line()`.


# A larger example: commuter trips in London

The minimal example dataset we've been using so far is fine for demonstrating the key concepts of OD data.
But for more advanced topic, and to get an idea of what is possible with OD data at a city level, it helps to have a larger dataset.

We will use an example dataset representing commuting in London, accessed as follows (note: these code chunks are not evaluated in the vignette because it starts by downloading 2.4 million rows and could take a few minutes to run).
First, we can use the `pct` package to download official data from the UK (note the addition of the % active column):

```{r, eval=FALSE}
library(dplyr)

# get nationwide OD data
od_all <- pct::get_od()
nrow(od_all)
# > 2402201
od_all$Active <- (od_all$bicycle + od_all$foot) /
  od_all$all * 100
centroids_all <- pct::get_centroids_ew() %>% sf::st_transform(4326)
nrow(centroids_all)
# > 7201
london <- pct::pct_regions %>% filter(region_name == "london")
centroids_london <- centroids_all[london, ]
od_london <- od_all %>%
  filter(geo_code1 %in% centroids_london$msoa11cd) %>%
  filter(geo_code2 %in% centroids_london$msoa11cd)
od_london <- od_all[
  od_all$geo_code1 %in% centroids_london$msoa11cd &
    od_all$geo_code2 %in% centroids_london$msoa11cd,
]
```

```{r, eval=FALSE, echo=FALSE}
# aim: create a reproducible OD dataset
od_lnd <- od_london %>%
  select(-matches("rail|name|moto|car|tax|home")) %>%
  filter(geo_code2 == "E02000001") %>%
  top_n(4, wt = all)
z_lnd <- centroids_london %>%
  filter(msoa11cd %in% c(od$geo_code1, od$geo_code2))
```


Now that we have the input OD data (in `od_london`) and zones (population-weighted centroids in `cents_london` in this case), can can convert them to desire lines:

```{r, eval=FALSE}
desire_lines_london <- od2line(od_london, centroids_london)
nrow(desire_lines_london)
# > 352654
```

Even after filering flows to keep only those with origins *and* destinations in London, there are still more than 300k flows. That is a lot to plot.
So we'll further subset them, first so they only contain inter-zonal flows (which are actually lines, intra-zonal flows are lines with length 0, which are essentially points) and second to contain only flows containing above a threshold level of flows:

```{r, eval=FALSE}
min_trips_threshold <- 20
desire_lines_inter <- desire_lines_london %>% filter(geo_code1 != geo_code2)
desire_lines_intra <- desire_lines_london %>% filter(geo_code1 == geo_code2)
desire_lines_top <- desire_lines_inter %>% filter(all >= min_trips_threshold)
nrow(desire_lines_top)
# > 28879
```

If we do any analysis on this dataset, it's important to know how representative it is of all flows.
A crude way to do this is to calculate the proportion of lines and trips that are covered in the dataset:

```{r, eval=FALSE}
nrow(desire_lines_top) / nrow(desire_lines_london)
# > 0.08189046
sum(desire_lines_top$all) / sum(desire_lines_london$all)
# > 0.557343
```

This shows that only 8% of the lines contain more than half (55%) of the total number of trips.

# Plotting origin-destination data

Once you have an OD dataset of a size that can be plotted (20,000 desire lines is quick to plot on most computers) a logical next stage is to plot it, e.g. with `sf`'s `plot()` method:

```{r, eval=FALSE}
plot(desire_lines_top["all"])
```

```{r, echo=FALSE}
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61058906-030a5c80-a3f0-11e9-90b5-d216964e9681.png")
```
You may be disapointed by the result, which is more of a 'hay stack' plot than an intuitive illustration of flows across the city.
To overcome this issue, you can set the aesthetics to emphasize with important flows, e.g. by line width in `sf`'s plotting system:

```{r, eval=FALSE}
lwd <- desire_lines_top$all / mean(desire_lines_top$all) / 10
desire_lines_top$percent_dont_drive <- 100 - desire_lines_top$car_driver / desire_lines_top$all * 100
plot(desire_lines_top["percent_dont_drive"], lwd = lwd, breaks = c(0, 50, 70, 80, 90, 95, 100))
```

```{r, echo=FALSE}
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/62073083-e5ceee00-b237-11e9-9cc7-8bf62d0e9b3f.png")
```

This is better, but is still not ideal: the code was not intuitive to write, and the result is still not publication quality.
Instead, it makes sense to make a dedicated mapping package such **tmap**, as outlined in the [visualisation chapter](https://geocompr.robinlovelace.net/adv-map.html) of the open source book *Geocomputation with R* [@lovelace_geocomputation_2019].
As shown in the transport chapter of that book, OD flows can be visualised with the following code:

```{r, eval=FALSE}
library(tmap)
desire_lines_top <- desire_lines_top %>%
  arrange(Active)
tm_shape(london) + tm_borders() +
  tm_shape(desire_lines_top) +
  tm_lines(
    palette = "plasma", breaks = c(0, 5, 10, 20, 40, 100),
    lwd = "all",
    scale = 9,
    title.lwd = "Number of trips",
    alpha = 0.5,
    col = "Active",
    title = "Active travel (%)",
    legend.lwd.show = FALSE
  ) +
  tm_scale_bar() +
  tm_layout(
    legend.bg.alpha = 0.5,
    legend.bg.color = "white"
  )
```

```{r, echo=FALSE}
# tmap_save(.Last.value, "tmap-london.png")
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61066243-12dc6d80-a3fd-11e9-8805-826a47c553f6.png")
```

The above plot contains much information, providing a visual overview of the transport pattern in the city, telling us that:

- It is a monocentric city, with most flows going to the centre.
- Active transport is geographically dependent, dominating in the central north of the city and with limited uptake on the outskirts of the city.
- Although the city centre dominates, there are many small clusters of flows in the outer region, for example near Heathrow airport, which is located in the far west of the map.

Plotting OD data in this way can tell us much about cities, each of which has a different travel pattern.
You can use the same code to visualise mobility patterns in any city.
See Section [12.4](https://geocompr.robinlovelace.net/transport.html#desire-lines) of *Geocomputation with R* to see results for Bristol, a more polycentric city with a lower average percentage of travel by walking and cycling.

```{r, eval=FALSE, echo=FALSE}
saveRDS(od_all, "od_all.Rds")
piggyback::pb_upload("od_all.Rds")
```

# Summaries by origin and destination

It is possible to group OD data by origin and destination to gain information at the zone level.
The code and resulting plot below, for example, summarises the number of people departing from each zone by mode:

```{r, eval=FALSE}
zones_london <- pct::get_pct_zones("london") %>%
  select("geo_code")
origin_attributes <- desire_lines_top %>%
  sf::st_drop_geometry() %>%
  group_by(geo_code1) %>%
  summarize_if(is.numeric, sum) %>%
  dplyr::rename(geo_code = geo_code1)
# origin_attributes <-
zones_origins <- left_join(zones_london, origin_attributes, by = "geo_code")
plot(zones_origins, border = NA)
```

```{r, echo=FALSE}
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61067619-e7a74d80-a3ff-11e9-8c15-7467717b36ec.png")
```

We can observe a number of features, including that:

- Rail is much more common in the south, reflecting the greater density of the local rail network, with short distances between stops, in the South of the city.
- Cars dominat in the outer fringes, especiall in the West.
- Taxi and motorbike use have intriguing clusters in the West (perhaps around the wealthy Kensington area for taxis).

The pattern is quite different when we calculate the destinations:

```{r, eval=FALSE}
destination_attributes <- desire_lines_top %>%
  sf::st_drop_geometry() %>%
  group_by(geo_code2) %>%
  summarize_if(is.numeric, sum) %>%
  dplyr::rename(geo_code = geo_code2) %>%
  mutate_at(vars(-matches("geo_|all")), funs(. / all)) %>%
  left_join(zones_london, ., by = "geo_code")

plot(destination_attributes, border = NA)
```


```{r, echo=FALSE}
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/61069409-27703400-a404-11e9-9c83-1cd5f2397260.png")
```

# Further reading

Despite the importance of origin-destination datasets for transport research, there are surprisingly few guides dedicated to working with them using open source software.
The following suggestions are based on my own reading  --- if you have any other suggestions of good resources for working with OD data, let me know!

- Section [12.4](https://geocompr.robinlovelace.net/transport.html#desire-lines) of *Geocomputation with R* [@lovelace_geocomputation_2019] puts OD data in the wider context of geographic transport data.
- @martin_origin-destination_2018 describe methods for classifying OD pairs based on demographic data.
- The [kepler.gl](https://kepler.gl/demo/ukcommute) website provides a nifty web application for visualising OD data.
- Documentation for the open source microscopic transport modelling software [SUMO](https://sumo.dlr.de/userdoc/Demand/Importing_O/D_Matrices.html) describes ways of reading-in OD file formats not covered in this vignette.
- An excellent introduction to modelling and visualising OD data in the [introductory vignette](https://github.com/rCarto/flows/blob/master/vignettes/flows.Rmd) of the `flows` R package.

# Summary

In summary, `stplanr` provides many functions for working with OD data.
OD data is an important component in transport planning and modelling that can, if used with creativity and skill, could assist with sustainable transport planning and the global [transition away from fossil fuels](https://www.sei.org/wp-content/uploads/2019/01/realizing-a-just-and-equitable-transition-away-from-fossil-fuels.pdf).
There are many other things can be done with OD data, some of which could be supported by future versions of this package.
To suggest new features of otherwise get in touch, see the `stplanr` issue tracker at [github.com/ropensci/stplanr](https://github.com/ropensci/stplanr/issues).

<!-- # Applications -->

<!-- This is represented in the file `lines_cars.Rds`, representing the top 20,000 desire lines at the MSOA-MSOA level in England and Wales by the number of car km used for travel to work, which can be downloaded, read-in and plotted as follows: -->


```{r, out.width="100%", warning=FALSE, eval=FALSE, echo=FALSE}
u <- "https://github.com/ropensci/stplanr/releases/download/0.2.9/lines_cars.Rds"
f <- file.path(tempdir(), "lines_cars.Rds")
download.file(u, f)
lines_cars <- readRDS(f)
plot(lines_cars["car_km"], lwd = lines_cars$car_km / 1000)
```

<!-- Based on the estimate of the average energy use per km being 2.5 MJ, and that these return trips are made on average 200 times per year, with a circuity of 1.3, we can estimate the total energy use of the 'high energy commutes' as follows: -->

```{r, eval=FALSE, echo=FALSE}
sum(lines_cars$car_km * 2.5 * 200) / 1e9
```

<!-- That represents ~10 petajoules (PJ), only for the top 20,000 most energy intensive commutes. -->
<!-- That may seem like a lot, but represents only a fraction of the UK's total energy use of [~200 Mtoe](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/729451/DUKES_PN.pdf) (8400 PJ). -->


```{r, echo=FALSE, eval=FALSE}
# out-takes and test code
# demonstrate bug/feature in sf
library(sf)
m <- matrix(c(
  0, 0,
  1, 0,
  0, 1,
  0, 0
), ncol = 2)
p <- st_polygon(list(m))

m <- matrix(c(
  0, 0,
  1, 0,
  0, NA,
  0, 0
), ncol = 2)
p <- st_polygon(list(m))
plot(p)

l <- st_linestring(m)
plot(l)
plot(p)
m <- matrix(c(0, 0, 0, NA), ncol = 2)
l <- st_linestring(m)
plot(l)
```

```{r, echo=FALSE, eval=FALSE}
usethis::use_data(od_data_sample)
# aim: get top flows by car use multiplied by distance
# subset flows with more than n people driving:
od_cars <- od_data_all[od_data_all$car_driver >= 50, ]
cents_ew <- pct::get_centroids_ew()
od_cars <- od_cars[
  od_cars$geo_code1 %in% cents_ew$msoa11cd &
    od_cars$geo_code2 %in% cents_ew$msoa11cd,
]
desire_lines_cars <- od2line(od_cars, cents_ew)
plot(desire_lines_cars[1:999, ])
desire_lines_cars$euclidean_distance_m <- as.numeric(sf::st_length(desire_lines_cars)) / 1000
desire_lines_cars$car_km <- desire_lines_cars$car_driver * desire_lines_cars$euclidean_distance_m
lines_cars <- dplyr::top_n(desire_lines_cars, 20000, car_km)
summary(lines_cars$car_driver)
plot(lines_cars["car_km"])
saveRDS(lines_cars, "lines_cars.Rds")
piggyback::pb_upload("lines_cars.Rds")
```

# References
---
title: 'stplanr: A Package for Transport Planning'
author:
  - name: Robin Lovelace
    affiliation: University of Leeds
    address:
    - 34-40 University Road
    - LS2 9JT, UK
    email:  r.lovelace@leeds.ac.uk
  - name: Richard Ellison
    affiliation: University of Sydney
    address:
    - 378 Abercrombie Street
    - Darlington, NSW 2008, Australia
    email:  richard.ellison@sydney.edu.au
abstract: >
  Tools for transport planning should be flexible, scalable and transparent. The **stplanr** package demonstrates and provides a home for such tools, with an emphasis on spatial transport data and non-motorized modes. **stplanr** facilitates common transport planning tasks including: downloading and cleaning transport datasets; creating geographic 'desire lines from origin-destination (OD) data; route assignment, via the `SpatialLinesNetwork` class and interfaces to routing services such as CycleStreets.net; calculation of route segment attributes such as bearing and aggregate flow; and `travel watershed' analysis. This paper demonstrates this functionality using reproducible examples on real transport datasets. More broadly, the experience shows open source software can form the basis of a reproducible transport planning workflow. **stplanr**, alongside other packages and open source projects, could provide a more transparent and democratically accountable alternative to the current approach which is heavily reliant on proprietary and technically and financially inaccessible software.
bibliography:
  - references.bib
  - stplanr-citation.bib
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{stplanr: A Package for Transport Planning}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo=FALSE}
knitr::opts_chunk$set(fig.width = 7, fig.height = 5, eval = FALSE)
```

## Note

This paper has now been peer reviewed and published by the R Journal.
Please see the published version at [journal.r-project.org](https://journal.r-project.org/archive/2018/RJ-2018-053/index.html) and cite it as @lovelace_stplanr_2018.

# Introduction

Transport planning can broadly be defined as the process of designing
and evaluating transport interventions [@willumsen_modelling_2011]
usually with the ultimate aim of improving transport systems from
economic, social and environmental perspectives. This inevitably
involves a degree of subjective judgment and intuition. With the
proliferation of new transport datasets — and the increasing
availability of hardware and software to make sense of them — there is
great potential for the discipline to become more evidence-based and
scientific [@balmer_matsim-t:_2009]. Transport planners have always
undertaken a wide range computational activities
[@boyce_forecasting_2015], but with the digital revolution the demands
have grown beyond the capabilities of a single, monolithic product. The
diversity of tasks , and need for democratic accountability in public
decision making, suggests that future-proof transport planning software
should be:


-   flexible, able to handle a wide range of data formats

-   scalable, able to work at multiple geographic levels from single
    streets to large cities and regions

-   robust and reliable, tested on a range of datasets and able to work
    ’out of the box’ in a range of real-world projects

-   open source and reproducible, ensuring transparency and encouraging
    citizen science

This paper sets out to demonstrate that open source software with a
command-line interface (CLI) can provide a foundation for transport
planning software that meets each of these criteria. R provides a strong
basis for progress in this direction because it already contains
functionality used in common transport planning workflows. , and greatly
improved R’s spatial abilities [@bivand_applied_2013], work that is
being consolidated and extended in the recent package.

Building on these foundations a number of spatial packages have been
developed for applied domains including: disease mapping and modelling,
with packages such as and
[@kim_spatialepi:_2016; @brown_diseasemapping:_2016]; spatial ecology,
with the **adehabitat** family of packages [@calenge_package_2006]; and
visualisation, with packages such as **SpatialEpi**, **diseasemapping** and @RJ-2016-005. However,
there has been little prior work to develop R functionality designed
specifically for transport planning, with the notable exceptions of
[TravelR](https://r-forge.r-project.org/projects/travelr/) (a package on
R-Forge last updated in 2012) and
[gtfsr](https://github.com/ropensci/gtfsr) (a package for handling
General Transit Feed Specification (GTFS) data).

The purpose of **stplanr** is to provide a toolbox rather than a
specific solution for transport planning, with an emphasis on spatial
data and active modes. This emphasis is timely given the recent emphasis
on sustainability [@banister_sustainable_2008] and ‘Big Data’
[@zheng_big_2016] in the wider field of transport planning. A major
motivation was the lack of R packages, and open source software in
general, for transport applications. This may be surprising given the
uqiquity of transport problems;^[
Many people can think of things that could be improved on their local transport networks, especially for walking, cycling and wheel-chairs.
But most lack the evidence to communicate the issues, and potential solutions, to others.
]
R’s proficiency at handling spatial,
temporal and travel survey data that describe transport systems; and the
growing popularity of R in applied domains
[@jalal_overview_2017; @moore_why_2017]. Another motivation is the
growth in open access datasets: the main purpose of early versions of
the package was to process open origin-destination data
[@lovelace_propensity_2017].

R is already used in transport applications, as illustrated by recent
research that applies packages from other domains to transport problems.
For instance, @efthymiou_use_2012 use R to
analyse the data collected from an online survey focused on car-sharing,
bicycle-sharing and electric vehicles. @efthymiou_use_2012
also used R to collect and analyse
transport-related data from Twitter using packages including , and .
These packages were used to download, parse and plot the Twitter data
using a method that can be repeated and the results reproduced or
updated. More general statistical analyses have also been conducted on
transport-related datasets using packages including and
[@diana_studying_2012; @cerin_walking_2013]. Despite the rising use of R
for transport research, there has yet been to be a package for transport
planning.

The design of the R language, with its emphasis on flexibility, data
processing and statistical modelling, suggests it can provide a powerful
environment for transport planning research. There are many quantitative
methods in transport planning, many of which fit into the classic ‘four
stage’ transport model which involves the following steps
[@willumsen_modelling_2011]: (1) trip *generation* to estimate trip
freqency from origins; (2) *distribution* of trips to destinations; (3)
*modal split* of trips between walking, cycling, buses etc.; (4)
*assignment* of trips to the transport route network. To this we would
like to add two more stages for the big data age: (0) data processing
and exploration; and (5) validation. This sequence is not the only way
of transport modelling and some have argued that its dominance has
reduced innovation. However it is certainly a common approach and
provides a useful schema for classifying the kinds of task that
**stplanr** can tackle:

-   Accessing and processing of data on transport infrastructure and
    behaviour (stage 0)

-   Analysis and visualisation of the transport network (0)

-   Analysis of origin-destination (OD) data and the visualisation of
    resulting ‘desire lines’

-   The allocation of desire lines to roads and other guideways via
    routing services

-   The aggregation of routes to estimate total levels of flow on
    segments throughout the transport network

-   Development of models to estimate transport behaviour currently and
    under various scenarios of change

-   The calculation of ‘catchment areas’ affected by transport
    infrastructure

The automation of such tasks can assist researchers and practitioners to
create evidence for decision making. If the data processing and analysis
stages are fast and painless, more time can be dedicated to
visualisation and decision making. This should allow researchers to
focus on problems, rather than on clunky graphical user interfaces
(GUIs), and ad-hoc scripts that could be generalised. Furthermore, if
the process can be made reproducible and accessible (e.g. via online
visualisation packages such as **shiny**), this could help transport planning
move away from reliance on ‘black boxes’ [@waddell_urbansim:_2002] and
empower citizens to challenge decisions made by transport planning
authorities based on the evidence [@hollander_transport_2016]. There are
many advantages of using a scriptable, interactive and open source
language such as R for transport planning. Such an approach enables:
reproducible research; the automation and sharing of code between
researchers; reduced barriers to innovation as anyone can create new
features for the benefit of all planners; easier interaction with non
domain experts (who will lack dedicated software); and integration with
other software systems, as illustrated by the use of to generate
JavaScript for sharing interactive maps for transport planning, as used
in the publicly accessible Propensity to Cycle Tool
[@lovelace_propensity_2017]. Furthermore, R has a strong user community
which can support newcomers (*stplanr* was peer reviewed thanks to the
community surrounding ROpenSci). The advantages of using R specifically
to develop the functionality described in this paper are that it has
excellent geo-statistical capabilities [@pebesma_software_2015],
visualisation packages (e.g. **tmap**, **ggplot2**), support for logit models (which are
useful for modelling modal shift), and support for the many formats that
transport datasets are stored in (e.g. via the **haven**  and **rio** packages).

# Package structure and functionality

The package can be installed and loaded in the usual way (see the package's
[README](https://github.com/ropensci/stplanr) for dependencies and access to development versions):

```{r, eval=FALSE}
install.packages("stplanr")
```


```{r}
library(stplanr)
```

As illustrated by the message emitted when **stplanr** is loaded, it depends on \CRANpkg{sp}. This means that the spatial data classes commonly used in the package will work with generic R functions such as `summary`, `aggregate` and, as illustrated in the figures below, `plot` \citep{bivand_applied_2013}.

## Core functions and classes

The package's core functions are structured around 3 common types of spatial transport data:

- Origin-destination (OD) data, which report the number of people travelling between origin-destination pairs. This type of data is not explicitly spatial (OD datasets are usually represented as data frames) but represents movement over space between points in geographical space. An example is provided in the `flow` dataset.
- Line data, one dimensional linear features on the surface of the Earth. These are typically stored as a `SpatialLinesDataFrame`.
- Route data are special types of lines which have been allocated to the transport network. Routes typically result from the allocation of a straight 'desire line' allocated to the route network with a `route_` function. Route network represent many overlapping routes. All are typically stored as `SpatialLinesDataFrame`.

For ease of use, functions focussed on each data type have been developed with names prefixed with `od_`, `line_` and `route_` respectively. A selection of these is presented in Table 1. Additional 'core functions' could be developed, such as those prefixed with `rn_` (for working with route network data) and `g_` functions for geographic operations such as buffer creation on lat/lon projected data (this function is currently named `buff_geo`). We plan to elicit feedback on such changes before implementing them.

```{r, echo=FALSE, results='asis', message=FALSE}
# stplanr_funs = ls("package:stplanr")
# sel_core = grep(pattern = "od_|^line_|route_", x = stplanr_funs)
# core_funs = stplanr_funs[sel_core]
# args(name = core_funs[1])
fun_table <- read.csv("fun_table.csv", stringsAsFactors = FALSE, check.names = FALSE)
knitr::kable(fun_table, caption = "Selection of functions for working with or generating OD, line and route data types.")
```

With a tip of the hat to the concept of type stability (e.g. as implemented in \CRANpkg{dplyr}), we also plan to make the core functions of **stplanr** more type-stable in future releases.  Core functions, which begin with the prefixes listed above, could follow \CRANpkg{dplyr}'s lead and return only objects with the same class as that of the input. However there are limitations to this approach: it will break existing functionality and mean that output objects have a larger size than necessary (`line_bearing`, for example, does not need to duplicate the spatial data contained in its input). Instead, we plan to continue to name functions around the type of *input* data they take, but are open minded about function input-output data class conventions, especially in the context of the new class system implemented in \CRANpkg{sf}.

A class system has not been developed for each data type (this option is discussed in the final section). The most common data types used in **stplanr** are assumed to be data frames and spatial datasets.

Transport datasets are very diverse. There are therefore many other functions which have more ad-hock names. Rather attempt a systematic description of each of **stplanr**'s functions (which can be gleaned from the online manual) it is more illuminating to see how they work together, as part of a transport planning workflow. As with most workflows, this begins with data access and ends with visualisation.

## Accessing and processing transport data

Gaining access to data is often the first stage in transport research. This is often a long and protracted process which is thankfully becoming easier thanks to the 'open data' movement and packages such as **tigris** for making data access from within R easier \citep{walker_tigris:_2016}.

**stplanr** provides a variety of different functions that facilitate importing common data formats used for transport analysis into R.
Although transport analysis generally requires some transport-specific datasets, it also typically relies heavily on common sources of data including census data.
This being the case, **stplanr** also includes functions that may be useful to those not involved in transport research.
This includes the `read_table_builder` function for importing data from the Australian Bureau of Statistics (ABS) and the UK's Stats19 road traffic casualty dataset. A brief example of the latter is demonstrated below, which begins with downloading the data (warning this downloads ~100 MB of data):

```{r, eval=FALSE}
dl_stats19() # download and extract stats19 road traffic casualty data
```

```
#> [1] "Data saved at: /tmp/RtmpppF3E2/Accidents0514.csv"
#> [2] "Data saved at: /tmp/RtmpppF3E2/Casualties0514.csv"
#> [3] "Data saved at: /tmp/RtmpppF3E2/Vehicles0514.csv"  
```

Once the data has been saved in the default directory, determined by `tempdir()`, it can be read-in and cleaned with the `read_stats19_` functions (note these call `format_stats19_` functions internally to clean the datasets and add correct labels to the variables):

```{r, eval=FALSE}
ac <- read_stats19_ac()
ca <- read_stats19_ca()
ve <- read_stats19_ve()
```

The resulting datasets (representing accident, casualty and vehicle level data, respectively) can be merged and made geographic, as illustrated below:

```{r, eval=FALSE}
library(dplyr)
ca_ac <- inner_join(ca, ac)
ca_cycle <- ca_ac %>%
  filter(Casualty_Severity == "Fatal" & !is.na(Latitude)) %>%
  select(Age = Age_of_Casualty, Mode = Casualty_Type, Longitude, Latitude)
ca_sp <- SpatialPointsDataFrame(coords = ca_cycle[3:4], data = ca_cycle[1:2])
```

Now that this casualty data has been cleaned, subsetted (to only include serious cycle crashes) and converted into a spatial class system, we can analyse them using geographical datasets of the type commonly used by **stplanr**.
The following code, for example, geographically subsets the dataset to include only crashes that occured within the bounding box of a [route network dataset](https://github.com/ropensci/stplanr/blob/master/data/route_network.rda?raw=true) provided by **stplanr** (from version 0.1.7 and beyond) using the function `bb2poly`, which converts a spatial dataset into a box, represented as a rectangular `SpatialPolygonsDataFrame`:

```{r, eval=FALSE}
data("route_network") # devtools::install_github("ropensci/splanr")version 0.1.7
proj4string(ca_sp) <- proj4string(route_network)
bb <- bb2poly(route_network)
proj4string(bb) <- proj4string(route_network)
ca_local <- ca_sp[bb, ]
```

The above code chunk shows the importance of understanding geographical data when working with transport data.
It is only by converting the casualty data into a spatial data class, and adding a coordinate reference system (CRS), that transport planners and researchers can link this important dataset back to the route network.
We can now perform GIS operations on the results.
The next code chunk, for example, finds all the fatalities that took place within 100 m of the route network, using the function `buff_geo`:

```{r, echo=FALSE}
bb <- bb2poly(route_network)
load("reqfiles.RData")
```

```{r, message=FALSE}
rnet_buff_100 <- geo_buffer(route_network, width = 100)
ca_buff <- ca_local[rnet_buff_100, ]
```

These can be visualised using base R graphics, extended by \CRANpkg{sp}, as illustrated in Figure \ref{fig:fats}. This provides a good start for analysis but for publication-quality plots and interactive plots, designed for public engagement, we recommend using dedicated visualisation packages that work with spatial data such as \CRANpkg{tmap}.

```{r fats, fig.cap="Road traffic fatalities in the study area downloaded with with stplanr (crosses). Deaths that happened within 100 m of the route network are represented by circles.", out.width="50%", fig.align="center"}
plot(bb, lty = 4)
plot(rnet_buff_100, col = "grey", add = TRUE)
points(ca_local, pch = 4)
points(ca_buff, cex = 3)
```

## Creating geographic desire lines

Perhaps the most common type of aggregate-level transport information is origin-destination ('OD') data.
This can be presented either as a matrix or (more commonly) a long table
of OD pairs. An example of this type of raw data is provided
below (see `?flow` to see how this dataset was created).  

```{r}
data("flow", package = "stplanr")
head(flow[c(1:3, 12)])
```

Although the flow data displayed above describes movement over
geographical space, it contains no explicitly geographical
information. Instead, the coordinates of the origins and
destinations are linked to a separate geographical dataset
which also must be loaded to analyse the flows. This is a
common problem solved by the function `od2line`.
The geographical data is a set of points representing centroids
of the origin and destinations, saved as a
`SpatialPointsDataFrame`. Geographical data in R is best
represented as such `Spatial*` objects, which use the
`S4` object engine. This explains the close integration of
**stplanr** with R's spatial packages, especially **sp**, which
defines the `S4` spatial object system.

```{r}
data("cents", package = "stplanr")
as.data.frame(cents[1:3, -c(3, 4)])
```

We use `od2line` to combine `flow` and `cents`, to join
the former to the latter. We will visualise the
`l` object created below in the next section. 

```{r, warning=FALSE}
l <- od2line(flow = flow, zones = cents)
```

The data is now in a form that is much easier to analyse. We can plot the
data with the command `plot(l)`, which was not possible before. Because the
`SpatialLinesDataFrame` object also contains data per line, it also helps
with visualisation of the flows, as illustrated in Figure \ref{fig:lines_routes}.

## Allocating flows to the transport network

A common problem faced by transport researchers is network
allocation: converting the 'as the crow flies' lines illustrated in the
figure above into routes. These are the complex, winding
paths that people and
animals make to avoid obstacles such as buildings and to make the journey
faster and more efficient (e.g. by following the route network).

This is difficult (and was until recently near impossible using free software)
because of the size and complexity of transport networks, the complexity
of realistic routing algorithms and need for context-specificity in the routing
engine. Inexperienced cyclists, for example, would take a very different route
than a heavy goods vehicle. **stplanr** tackles this issue by using 3rd party APIs to provide
route-allocation. 

Route allocation is undertaken by \code{route\_}
functions such as \code{route\_cyclestreets}
and \linebreak \code{route\_graphhopper}.
These allocate a single OD pair, represented as a text string to be 'geo-coded', a pair of of coordinates, or two `SpatialPoints` objects, representing origins and destinations.
This is illustrated below with `route_cyclestreet`, which uses the
[CycleStreets.net API](https://www.cyclestreets.net/api/), a routing service "by cyclists for cyclists" that offers a range route strategies (primarily 'fastest', 'quietest' and 'balanced') that are based on a
detailed analysis of cyclist
wayfinding:^[An
API key is needed for this function to work. This can be requested (or purchased for large scale routing) from [cyclestreets.net/api/apply](https://www.cyclestreets.net/api/apply/). See `?route_cyclestreet` for details.
Thanks to Martin Lucas-Smith and Simon Nuttall for making this possible.]

```{r, eval=FALSE}
route_bl <- route_cyclestreets(from = "Bradford", to = "Leeds")
route_c1_c2 <- route_cyclestreets(cents[1, ], cents[2, ])
```

The raw output from routing APIs is usually provided as a JSON or GeoJSON text string. By default, `route_cyclestreet` saves a number of key variables (including length, time, hilliness and busyness variables generated by CycleStreets.net) from the attribute data provided by the API. If the user wants to save the raw output, the `save_raw` argument can be used:

```{r, eval=FALSE}
route_bl_raw <- route_cyclestreets(from = "Bradford", to = "Leeds", save_raw = TRUE)
```

Additional arguments taken by the `route_` functions depend on the routing function in question. By changing the `plan` argument of `route_cyclestreet` to `fastest`, `quietest` or `balanced`, for example, routes favouring speed, quietness or a balance between speed and quietness will be saved, respectively.

To automate the creation of route-allocated lines over many desire lines, the `line2route` function loops over each line, wrapping any `route_` function as an input. The output is a `SpatialLinesDataFrame` with the same number of dimensions as the input dataset (see the right panel in Figure \ref{fig:lines_routes}).

```
routes_fast <- line2route(l = l, route_fun = route_cyclestreet)
```

The result of this 'batch routing' exercise is illustrated in Figure \ref{fig:lines_routes}. The red lines in the left hand panel are very different from the hypothetical straight 'desire lines' often used in transport research, highlighting the importance of this route-allocation functionality. 

```{r lines_routes, out.width='50%', fig.cap='Visualisation of travel desire lines, with width proportional to number of trips between origin and destination (black) and routes allocated to network  (red) in the left-hand panel. The right hand panel shows the route network dataset generated by overline().', fig.show='hold'}
plot(route_network, lwd = 0)
plot(l, lwd = l$All / 10, add = TRUE)
lines(routes_fast, col = "red")
routes_fast$All <- l$All
rnet <- overline(routes_fast, "All", fun = sum)
rnet$flow <- rnet$All / mean(rnet$All) * 3
plot(rnet, lwd = rnet$flow / mean(rnet$flow))
```

To estimate the amount of capacity needed at each segment on the transport network, the `overline` function demonstrated above, is used to divide line geometries into unique segments and aggregate the overlapping values.
The results, illustrated in the right-hand panel of Figure \ref{fig:lines_routes}, can be used to estimate where there is most need to improve the transport network, for example informing the decision of where to build new bicycle paths.

Limitations with the `route_cyclestreet` routing API include its specificity, to one mode (cycling) and a single region (the UK and part of Europe). To overcome these limitations, additional routing APIs were added with the functions `route_graphhopper`, `route_transportapi_public` and `viaroute`. These interface to Graphhopper, TransportAPI and the Open Source Routing Machine (OSRM) routing services, respectively. The great advantage of OSRM is that it allows you to run your own routing services on a local server, greatly increasing the rate of route generation.

A short example of finding the route by car and bike between New York and Oaxaca demonstrates how `route_graphhopper` can collect geographical and other data on routes by various modes, anywhere in the world. The output, shown in Table \ref{tab:xtnyoa}, shows that the function also saves time, distance and (for bike trips)
vertical distance climbed for the trips.

```{r, eval=FALSE, out.width='\\textwidth'}
ny2oaxaca1 <- route_graphhopper("New York", "Oaxaca", vehicle = "bike")
ny2oaxaca2 <- route_graphhopper("New York", "Oaxaca", vehicle = "car")
rbind(ny2oaxaca1@data, ny2oaxaca2@data)
```

```{r, eval=FALSE, echo=FALSE}
nytab <- rbind(ny2oaxaca1@data, ny2oaxaca2@data)
nytab <- cbind(Mode = c("Cycle", "Car"), nytab)
xtnyoa <- xtable(nytab, caption = "Attribute data from the route\\_graphhopper function, from New York to Oaxaca, by cycle and car.", label = "tab:xtnyoa")
print.xtable(xtnyoa, include.rownames = FALSE)
plot(ny2oaxaca1)
plot(ny2oaxaca2, add = TRUE, col = "red")

ny2oaxaca1@data
ny2oaxaca2@data
```

|     time|    dist| change_elev|
|--------:|-------:|-----------:|
| 17522.73| 4885663|    87388.13|
|  2759.89| 4754772|          NA|

<!-- ## Calculating geographic attributes of transport routes -->

## Modelling travel catchment areas

Accessibility to transport services is a particularly important topic when considering public transport or active travel because of the frequent steep reduction in use as distances to access services (or infrastructure) increase.
As a result, the planning for transport services and infrastructure frequently focuses on several measures of accessibility including distance, but also travel times and frequencies and weighted by population.
The functions in **stplanr** are intended to provide a method of estimating these accessibility measures as well as calculating the population that can access specific services (i.e., estimating the catchment area).

Catchment areas in particular are a widely used measure of accessibility that attempts to both quantify the likely target group for a particular service, and visualise the geographic area that is covered by the service.
For instance, passengers are often said to be willing to walk up to 400 metres to a bus stop, or 800 metres to a railway station \citep{el-geneidy_new_2014}.
Although these distances may appear relatively arbitrary and have been found to underestimate the true catchment area of bus stops and railway stations  \citep{el-geneidy_new_2014,daniels_explaining_2013} they nonetheless represent a good, albeit somewhat conservative, starting point from which catchment areas can be determined.

In many cases, catchment areas are calculated on the basis of straight-line (or "as the crow flies") distances.
This is a simplistic, but relatively appealing approach because it requires little additional data and is straight-forward to understand.
**stplanr** provides functionality that calculates catchment areas using straight-line distances with the `calc_catchment` function.
This function takes a `SpatialPolygonsDataFrame` that contains the population (or other) data, typically from a census, and a `Spatial*` layer that contains the geometry of the transport facility.
These two layers are overlayed to calculate statistics for the desired catchments including proportioning polygons to account for the proportion located within the catchment area.

To illustrate how catchment areas can be calculated, **stplanr** contains some sample datasets stored in ESRI Shapefile format (a commonly used format for distributing GIS layers) that can together be used to calculate sample catchment areas.
One of these datasets (`smallsa1`) contains population data for Statistical Area 1 (SA1) zones in Sydney, Australia.
The second contains hypothetical cycleways aligned to streets in Sydney. The code below unzips the datasets and reads in the shapefiles.

```{r loadshapefiles, results='hide',message='hide'}
data_dir <- system.file("extdata", package = "stplanr")
unzip(file.path(data_dir, "smallsa1.zip"))
unzip(file.path(data_dir, "testcycleway.zip"))
sa1income <- as(sf::read_sf("smallsa1.shp"), "Spatial")
testcycleway <- as(sf::read_sf("testcycleway.shp"), "Spatial")
# Remove unzipped files
file.remove(list.files(pattern = "^(smallsa1|testcycleway).*"))
```

Calculating the catchment area is straightforward and in addition to specifying the required datasets, only a vector containing column names to calculate statistics and a distance is required.
Since proportioning the areas assumes projected data, unprojected data are automatically projected to either a common projection (if one is already projected) or a specified projection.
It should be emphasised that the choice of projection is important and has an effect on the results meaning setting a local projection is recommended to achieve the most accurate results.

```{r calccatchment, results='hide'}
catch800m <- calc_catchment(
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  distance = 800,
  projection = "austalbers",
  dissolve = TRUE
)
```

By looking at the data.frame associated with the SpatialPolygonsDataFrame that is returned from the `calc_catchment` function, the total population within the catchment area can be seen to be nearly 40,000 people.
The catchment area can also be plotted as with any other `Spatial*` object using the `plot` function using the code below with the result shown in Figure \ref{fig:catchmentplot}.

```{r catchmentplot, fig.cap='An 800 metre catchment area (red) associated with a cycle path (green) using straight-line distance in Sydney.'}
plot(sa1income, col = "light grey")
plot(catch800m, col = rgb(1, 0, 0, 0.5), add = TRUE)
plot(testcycleway, col = "green", add = TRUE)
```

This simplistic catchment area is useful when the straight-line distance is a reasonable approximation of the route taken to walk (or cycle) to a transport facility.
However, this is often not the case.
The catchment area in Figure \ref{fig:catchmentplot} initially appears reasonable but the red-shaded catchment area includes an area that requires travelling around a bay to access from the (green-coloured) cycleway.
To allow for more realistic catchment areas for most situations, **stplanr** provides the `calc_network_catchment` function that uses the same principle as `calc_catchment` but also takes into account the transport network.

To use `calc_network_catchment`, a transport network needs to be prepared that can be used in conjunction with the previous datasets.
Preparation of the dataset involves using the `SpatialLinesNetwork` function to create a network from a `SpatialLinesDataFrame`.
This function combines a `SpatialLinesDataFrame` with a graph network (using the \CRANpkg{igraph} package) to provide basic routing functionality.
The network is used to calculate the shortest actual paths within the specific catchment distance.
This process involves the following code:

```{r, echo=TRUE, message=FALSE, warning=FALSE, results='hide'}
unzip(file.path(data_dir, "sydroads.zip"))
sydroads <- as(sf::read_sf(".", "roads"), "Spatial")
file.remove(list.files(pattern = "^(roads).*"))
sydnetwork <- SpatialLinesNetwork(sydroads)
```

The network catchment is then calculated using a similar method as with `calc_catchment` but with a few minor changes.
Specifically these are including the `SpatialLinesNetwork`, and using the `maximpedance` parameter to define the distance, with distance being the additional distance from the network.
In contrast to the distance parameter that is based on the straight-line distance in both the `calc_catchment` and `calc_network_catchment` functions, the `maximpedance` parameter is the maximum value in the units of the network's weight attribute.
In practice this is generally distance in metres but can also be travel times, risk or other measures.

```{r, warning=FALSE}
netcatch800m <- calc_network_catchment(
  sln = sydnetwork,
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  maximpedance = 800,
  distance = 100,
  projection = "austalbers"
)
```

Once calculated, the network catchment area can be used just as the straight-line network catchment.
This includes extracting the catchment population of 128,000 and plotting the original catchment area together with the original area with the results shown in Figure \ref{fig:netcatchplot}:

```{r netcatchplot, fig.cap='A 800 metre network catchment are (blue) compared with a catchment area based on Euclidean distance (red) associated with a cycle path (green).'}
plot(sa1income, col = "light grey")
plot(catch800m, col = rgb(1, 0, 0, 0.5), add = TRUE)
plot(netcatch800m, col = rgb(0, 0, 1, 0.5), add = TRUE)
plot(testcycleway, col = "green", add = TRUE)
```

# Modelling and visualisation

## Modelling mode choice

Route-allocated lines allow estimation of *route distance* and 
*cirquity* (route distance divided by Euclidean distance).
These variables can help model the rate of flow between origins and
destination, as illustrated in the left-hand panel of Figure \ref{fig:euclidfastest}.
The code below demonstrates
how objects generated by **stplanr** can be used to undertake such analysis, with the `line_length` function used to find the distance, in meters, of lat/lon data.

```
l$d_euclidean <- line_length(l)
l$d_rf <- routes_fast@data$length
plot(l$d_euclidean, l$d_rf,
  xlab = "Euclidean distance", ylab = "Route distance")
abline(a = 0, b = 1)
abline(a = 0, b = 1.2, col = "green")
abline(a = 0, b = 1.5, col = "red")
```

```{r, echo=FALSE, message=FALSE}
l$d_euclidean <- line_length(l)
l$d_rf <- routes_fast$length
```

The left hand panel of Figure \ref{fig:euclidfastest} shows the expected strong correlation between
Euclidean ($d_E$) and fastest route ($d_{Rf}$) distance. However, some OD pairs
have a proportionally higher route distance than others, as illustrated
by distance from the black line in the above plot: this represents \emph{Circuity ($Q$)}: the ratio of network distance to Euclidean distance \citep{levinson_minimum_2009}:

$$
 Q = \frac{d_{Rf}}{d_E}
$$

An extension to the concept of cirquity is the 'quietness diversion factor' ($QDF$) of a desire line \citep{lovelace_propensity_2016}, the ratio of the route distance of a quiet route option ($d_{Rq}$) to that of the fastest:

$$
 QDF = \frac{d_{Rq}}{d_{Rf}}
$$

Thanks to the 'quietest' route option provided by `route_cyclestreet`, we can estimate average values for both metrics as follows:

```{r, eval=FALSE}
routes_slow <- line2route(l, route_cyclestreet, plan = "quietest")
```

```{r}
l$d_rq <- routes_slow$length # quietest route distance
Q <- mean(l$d_rf / l$d_euclidean, na.rm = TRUE)
QDF <- mean(l$d_rq / l$d_rf, na.rm = TRUE)
Q
QDF
```

The results show that cycle paths are not particularly direct in the study region by international standards \citep{crow_design_2007}.
This is hardly surprisingly given the small size of the sample and the short distances covered: $Q$ tends to decrease at a decaying rate with distance.
What is surprising is that $QDF$ is close to unity, which could imply that the quiet routes are constructed along direct, and therefore sensible routes.
We should caution against such assumptions, however:
It is a small sample of desire lines and, when time is explored, we find that the 'quietness diversion factor with respect to time' ($QDF_t$) is slightly larger:

```{r}
(QDFt <- mean(routes_slow$time / routes_fast$time, na.rm = TRUE))
```

## Models of travel behaviour

There are many ways of estimating flows between origins and destinations,
including spatial interaction models, the four-stage transport model
and gravity models ('distance decay'). **stplanr** aims eventually to facilitate
creation of many types of flow model. 

At present there are no functions for modelling distance decay, but this is something we would like to add in future versions of **stplanr**.
Distance decay is an especially important concept for sustainable transport
planning due to physical limitations on the ability of people to walk and
cycle large distances \citep{iacono_measuring_2010}.

We can explore the relationship between distance and the proportion of trips
made by walking, using the same object `l` generated by **stplanr**.

```{r euclidwalking1, fig.cap='Euclidean distance and walking trips', eval=FALSE}
l$pwalk <- l$On.foot / l$All
plot(l$d_euclidean, l$pwalk,
  cex = l$All / 50,
  xlab = "Euclidean distance (m)", ylab = "Proportion of trips by foot"
)
```

```{r euclidfastest, out.width='100%', fig.cap='Euclidean and fastest route distance of trips in the study area (left) and Euclidean distance vs the proportion of trips made by walking (right).', echo=FALSE}
par(mfrow = c(1, 2))
lgb <- sp::spTransform(l, CRSobj = sp::CRS("+init=epsg:27700"))
l$d_euclidean <- rgeos::gLength(lgb, byid = T)
l$d_rf <- routes_fast@data$length
plot(l$d_euclidean, l$d_rf,
  xlab = "Euclidean distance", ylab = "Route distance"
)
abline(a = 0, b = 1)
abline(a = 0, b = 1.2, col = "green")
abline(a = 0, b = 1.5, col = "red")
l$pwalk <- l$On.foot / l$All
plot(l$d_euclidean, l$pwalk,
  cex = l$All / 50,
  xlab = "Euclidean distance (m)", ylab = "Proportion of trips by foot"
)
```

Based on the right-hand panel in Figure \ref{fig:euclidfastest}, there is a clear negative relationship between distance of trips and the proportion of those trips made by walking.
This is unsurprising: beyond a certain distance (around 1.5km according
the the data presented in the figure above) walking is usually seen as too slow and other modes are considered.
According to the academic literature, this 'distance decay' is non-linear
and there have been a number of functions proposed to fit to distance decay
curves \citep{martinez_new_2013}. From the range of options we test below just two forms.
We will compare the ability of linear and log-square-root functions
to fit the data contained in `l` for walking.

```{r}
lm1 <- lm(pwalk ~ d_euclidean, data = l@data, weights = All)
lm2 <- lm(pwalk ~ d_rf, data = l@data, weights = All)
lm3 <- glm(pwalk ~ d_rf + I(d_rf^0.5),
  data = l@data, weights = All, family = quasipoisson(link = "log")
)
```

The results of these regression models can be seen using `summary()`.
Surprisingly, Euclidean distance was a better predictor of
walking than route distance, but no strong conclusions can be drawn from this finding, with such a small sample of desire lines (n = 42). The results are purely illustrative, of the kind of the possibilities created by using **stplanr** in conjuction with R's modelling capabilities (see Figure \vref{fig:euclidwalking2}).

```{r, echo=FALSE, eval=FALSE}
summary(lm1)
summary(lm2)
summary(lm3)
```

```{r euclidwalking2, fig.cap='Relationship between euclidean distance and walking', out.width="75%", fig.align="center"}
plot(l$d_euclidean, l$pwalk,
  cex = l$All / 50,
  xlab = "Euclidean distance (m)", ylab = "Proportion of trips by foot"
)
l2 <- data.frame(d_euclidean = 1:5000, d_rf = 1:5000)
lm1p <- predict(lm1, l2)
lm2p <- predict(lm2, l2)
lm3p <- predict(lm3, l2)
lines(l2$d_euclidean, lm1p)
lines(l2$d_euclidean, exp(lm2p), col = "green")
lines(l2$d_euclidean, exp(lm3p), col = "red")
```

## Visualisation

Visualisation is an important aspect of any transport study, as it enables researchers to communicate their findings to other researchers, policy-makers and, ultimately, the public. It may therefore come as a surprise that **stplanr** contains no functions for visualisation. Instead, users are encouraged to make use of existing spatial visualisation tools in R, such as **tmap**, **leaflet** and **ggmap** \citep{cheshire_spatial_2015,kahle_ggmap:_2013}.

Furthermore, with the development of online application frameworks such as **shiny**, it is now easier than ever to make the results of transport analysis and modelling projects available to the public. An example is the online interface of the Propensity to Cycle Tool (PCT). The results of the project, generated using **stplanr**, are presented at zone, desire line and Route Network levels \citep{lovelace_propensity_2016}. There is great potential to expand on the principle of publicly accessible transport planning tools via 'web apps', perhaps through new R packages dedicated to visualising transport data.

# Future directions of travel

This paper has demonstrated the great potential for R to be used for transport planning. R's flexibility, powerful GIS capabilities \citep{bivand_applied_2013} and free accessibility makes it well-suited to the needs of transport planners and researchers, especially those wanting to avoid the high costs of market-leading products. Rather than 'reinvent the wheel' (e.g. with a new class system), **stplanr** builds on existing packages and \CRANpkg{sp} classes to work with common transport data formats.

It is useful to see **stplanr**, and R for transport planning in general, as an addition tool in the transport planner's cabinet. It can be understood as one part of a wider movement that is making transport planning a more open and democratic process. Other developments in this movement include the increasing availability of open data \citep{naumova_building_2016} and the rise of open source products for transport modelling, such as [SUMO](https://www.dlr.de/ts/en/desktopdefault.aspx/tabid-9883/16931_read-41000/), [MATSim](https://www.matsim.org/) and MITSIMLAB \citep{saidallah_comparative_2016}.
**stplanr**, with its focus on GIS operations rather than microscopic vehicle-level behaviour, can complement such software and help make better use of new open data sources.

Because transport planning is an inherently spatial activity, **stplanr** occupies an important niche in the transport planning software landscape, with its focus on spatial transport data. There is great potential for development of **stplanr** in many directions.
Desirable developments include the additional of functions for modelling modal split, for examample with functions to create commonly distance decay curves which are commonly found in active travel research \citep{martinez_new_2013} and improving the computational efficiency of existing functions to make the methods more scalable for large databases.
Our priority for **stplanr** however, is to keep the focus on geographic functions for transport planning. There are many opportunities in this direction, including:

- Functions to assess the environment surrounding routes, e.g. via integration with the in-development **osmdata** package.
- Functions to match different GIS routes, perhaps building on the Hausdorf distance algorithm implemented in the \CRANpkg{rgeos} function `gDistance`.
- Additional functions for route-allocation of travel, e.g. via an interface to the OpenTripPlanner API.
- Functions for aggregating very large GPS trace datasets (e.g. into raster cells) for anonymisation and analysis/visualisation purposes.
- The creation of a class system for spatial transport datasets, such as to represent spatial route and a route networks (perhaps with classes named \code{"sr"} and \code{"srn"}). This is not a short-term priority and it would be beneficial to coincide such developments to a migration to \CRANpkg{sf} for spatial classes.

Such spatial data processing capabilities would increase the range of transport planning tasks that **stplanr** can facilitate. For all this planned development activity to be useful, it is vital that new functionality is intuitive. R has a famously steep learning curve. Implementing simple concepts such as consistent naming systems \citep{baath_state_2012} and ensuring 'type stability' can greatly improve the usability of the package. For this reason, much future work in **stplanr** will go into improving documentation and user-friendliness. 

Like much open source software **stplanr** is an open-ended project, a work-in-progress. We have set out clear motivations for developing transport planning capabilities in R and believe that the current version of **stplanr** (0.1.6) provides a major step in that direction compared with what was available a couple of years ago. But there is much more to do. We therefore welcome input on where the package's priorities should lie, how it should evolve in the future and how to ensure it is well-developed and sustained.

<!-- This creates the opportunity for extending existing modelling and data analysis data to better handle typical travel survey formats, e.g. by creating a system to represent travel surveys worldwide by merging trip, individual and household level data. Such non-spatial extensions of R for transport planning could build on the fast grouping and processing functions provided by \CRANpkg{dplyr}. This is an interesting potential direction of travel for **stplanr** or future packages for transport research. -->

<!-- If you would like to contribute or request new features, see -->
<!-- https://github.com/Robinlovelace/stplanr -->

# References
---
title: "Transport routing with stplanr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Transport routing with stplanr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = curl::has_internet()
)
```

```{r setup}
library(stplanr)
```

# Introduction

Routing is the process of identifying routes that enable movement between two geographic locations along the shortest path (based on mode-specific routing profiles) or in some other 'optimal' way, based on route network data.
Most open routing engines rely on OpenStreetMap (OSM) data.

We will use the example of the Isle of Wight to demonstrate routing engines.
To get OSM data for the Isle of Wight you can run the following commands:

```{r, eval=FALSE}
remotes::install_github("itsleeds/geofabrik")
library(geofabrik)
roads_iow = get_geofabrik(name = "Isle of Wight")
f = gf_filename("Isle of Wight")
file.copy(f, "iow.pbf")
options(osrm.server = "https://0.0.0.0:5000/", osrm.profile = "driving")
```


# OSRM

Routing services such as OpenStreetMap Routing Machine (OSRM) require an input network, usually from OSM.

We will use the `osrm` package:

```{r}
library(osrm)
```

In the system terminal run the following commands to make the [OSRM docker image](https://hub.docker.com/r/osrm/osrm-backend/) work for you.

```{r, engine='bash', eval=FALSE}
docker run -t -v "${PWD}:/data" osrm/osrm-backend osrm-extract -p /opt/car.lua /data/iow.pbf
docker run -t -v "${PWD}:/data" osrm/osrm-backend osrm-partition /data/iow.osrm
docker run -t -v "${PWD}:/data" osrm/osrm-backend osrm-customize /data/iow.osrm
docker run -t -i -p 5000:5000 -v "${PWD}:/data" osrm/osrm-backend osrm-routed --algorithm mld /data/iow.osrm
curl "https://127.0.0.1:5000/route/v1/driving/13.388860,52.517037;13.385983,52.496891?steps=true"
```

Now we can do routing in R!

On a single route:

```{r, eval=FALSE}
l = pct::wight_lines_30
p = line2points(l)
r = osrm::osrmRoute(src = p[1, ], dst = p[2, ], returnclass = "sf", overview = "full")
plot(r)
```

```{r, echo=FALSE}
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/86902789-577d1080-c106-11ea-91df-8d0180931562.png")
```

And to find many routes via the `route()` function, resulting in something like the figure below.

```{r, eval=FALSE}
routes_osrm = route(l = l, route_fun = osrmRoute, returnclass = "sf", overview = "full")
rnet_osrm = overline(routes_osrm, attrib = "bicycle")
mapview::mapview(rnet_osrm, lwd = rnet_osrm$bicycle / 10)
```

```{r, eval=FALSE, echo=FALSE}
system.time({
  routes_osrm = route(l = l, route_fun = osrmRoute, returnclass = "sf", overview = "full")
})
30 / 0.9 # around 30 routes per second
saveRDS(routes_osrm, "routes_osrm.Rds")
piggyback::pb_upload("routes_osrm.Rds")
```


```{r, echo=FALSE}
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/86858225-2970df80-c0b8-11ea-8394-07f98f1c8e8a.png")
```

```{r, eval=FALSE}
# tidy up
f = list.files(pattern = "iow")
unlink(x = f, recursive = TRUE)
```

Shut down the docker container.


```{r, eval=FALSE, engine='zsh'}
docker ps
docker stop stupefied_hopper
```


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_cyclestreets.R
\name{route_cyclestreets}
\alias{route_cyclestreets}
\title{Plan a single route with CycleStreets.net}
\usage{
route_cyclestreets(
  from,
  to,
  plan = "fastest",
  silent = TRUE,
  pat = NULL,
  base_url = "https://www.cyclestreets.net",
  reporterrors = TRUE,
  save_raw = "FALSE"
)
}
\arguments{
\item{from}{Text string or coordinates (a numeric vector of
\code{length = 2} representing latitude and longitude) representing a point
on Earth.}

\item{to}{Text string or coordinates (a numeric vector of
\code{length = 2} representing latitude and longitude) representing a point
on Earth. This represents the destination of the trip.}

\item{plan}{Text strong of either "fastest" (default), "quietest" or "balanced"}

\item{silent}{Logical (default is FALSE). TRUE hides request sent.}

\item{pat}{The API key used. By default this is set to NULL and
this is usually aquired automatically through a helper, api_pat().}

\item{base_url}{The base url from which to construct API requests
(with default set to main server)}

\item{reporterrors}{Boolean value (TRUE/FALSE) indicating if cyclestreets (TRUE by default).
should report errors (FALSE by default).}

\item{save_raw}{Boolean value which returns raw list from the json if TRUE (FALSE by default).}
}
\description{
Provides an R interface to the CycleStreets.net cycle planning API,
a route planner made by cyclists for cyclists.
The function returns a SpatialLinesDataFrame object representing the
an estimate of the fastest, quietest or most balance route.
Currently only works for the United Kingdom and part of continental Europe,
though other areas may be requested by contacting CycleStreets.
See \url{https://www.cyclestreets.net/api/}for more information.
}
\details{
This function uses the online routing service
CycleStreets.net to find routes suitable for cyclists
between origins and destinations. Requires an
internet connection, a CycleStreets.net API key
and origins and destinations within the UK (and various areas beyond) to run.

Note that if \code{from} and \code{to} are supplied as
character strings (instead of lon/lat pairs), Google's
geo-coding services are used via \code{geo_code()}.

You need to have an api key for this code to run.
Loading a locally saved copy of the api key text string
before running the function, for example, will ensure it
is available on any computer:

\verb{mytoken <- readLines("~/Dropbox/dotfiles/cyclestreets-api-key-rl") Sys.setenv(CYCLESTREETS = mytoken)}

if you want the API key to be available in future
sessions, set it using the .Renviron file
with \code{usethis::edit_r_environ()}

Read more about the .Renviron here: \code{?.Renviron}
}
\examples{

\dontrun{
from <- c(-1.55, 53.80) # geo_code("leeds")
to <- c(-1.76, 53.80) # geo_code("bradford uk")
json_output <- route_cyclestreets(from = from, to = to, plan = "quietest", save_raw = TRUE)
str(json_output) # what does cyclestreets give you?
rf_lb <- route_cyclestreets(from, to, plan = "fastest")
rf_lb@data
plot(rf_lb)
(rf_lb$length / (1000 * 1.61)) / # distance in miles
  (rf_lb$time / (60 * 60)) # time in hours - average speed here: ~8mph
}

}
\seealso{
line2route
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{odmatrix_to_od}
\alias{odmatrix_to_od}
\title{Convert origin-destination data from wide to long format}
\usage{
odmatrix_to_od(odmatrix)
}
\arguments{
\item{odmatrix}{A matrix with row and columns representing origin and destination zone codes
and cells representing the flow between these zones.}
}
\description{
This function takes a matrix representing travel between origins
(with origin codes in the \code{rownames} of the matrix)
and destinations
(with destination codes in the \code{colnames} of the matrix)
and returns a data frame representing origin-destination pairs.
}
\details{
The function returns a data frame with rows ordered by origin and then destination
zone code values and with names \code{orig}, \code{dest} and \code{flow}.
}
\examples{
odmatrix <- od_to_odmatrix(flow)
odmatrix_to_od(odmatrix)
flow[1:9, 1:3]
odmatrix_to_od(od_to_odmatrix(flow[1:9, 1:3]))
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_funs.R
\name{route_rolling_diff}
\alias{route_rolling_diff}
\title{Return smoothed differences between vector values}
\usage{
route_rolling_diff(x, lag = 1, abs = TRUE)
}
\arguments{
\item{x}{Numeric vector to smooth}

\item{lag}{The window size of the smoothing function. The default, 3, will take
the mean of values before, after and including each value.}

\item{abs}{Should the absolute (always positive) change be returned? True by default}
}
\description{
This function calculates a simple rolling mean in base R. It is useful for
calculating route characteristics such as mean distances of segments and
changes in gradient.
}
\examples{
r1 <- od_data_routes[od_data_routes$route_number == 2, ]
y <- r1$elevations
route_rolling_diff(y, lag = 1)
route_rolling_diff(y, lag = 2)
r1$elevations_diff_1 <- route_rolling_diff(y, lag = 1)
r1$elevations_diff_n <- route_rolling_diff(y, lag = 1, abs = FALSE)
d <- cumsum(r1$distances) - r1$distances / 2
diff_above_mean <- r1$elevations_diff_1 + mean(y)
diff_above_mean_n <- r1$elevations_diff_n + mean(y)
plot(c(0, cumsum(r1$distances)), c(y, y[length(y)]), ylim = c(80, 130))
lines(c(0, cumsum(r1$distances)), c(y, y[length(y)]))
points(d, diff_above_mean)
points(d, diff_above_mean_n, col = "blue")
abline(h = mean(y))
}
\seealso{
Other route_funs: 
\code{\link{route_average_gradient}()},
\code{\link{route_rolling_average}()},
\code{\link{route_rolling_gradient}()},
\code{\link{route_sequential_dist}()},
\code{\link{route_slope_matrix}()},
\code{\link{route_slope_vector}()}
}
\concept{route_funs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{l_poly}
\alias{l_poly}
\title{Line polygon}
\format{
A SpatialPolygon
}
\usage{
data(l_poly)
}
\description{
This dataset represents road width for testing.
}
\examples{
\dontrun{
l <- routes_fast[13, ]
l_poly <- geo_projected(l, rgeos::gBuffer, 8)
plot(l_poly)
plot(routes_fast, add = TRUE)
# allocate road width to relevant line
devtools::use_data(l_poly)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/google-functions.R
\name{dist_google}
\alias{dist_google}
\title{Return travel network distances and time using the Google Maps API}
\usage{
dist_google(
  from,
  to,
  google_api = Sys.getenv("GOOGLEDIST"),
  g_units = "metric",
  mode = c("bicycling", "walking", "driving", "transit"),
  arrival_time = ""
)
}
\arguments{
\item{from}{Two-column matrix or data frame of coordinates representing
latitude and longitude of origins.}

\item{to}{Two-column matrix or data frame of coordinates representing
latitude and longitude of destinations.}

\item{google_api}{String value containing the Google API key to use.}

\item{g_units}{Text string, either metric (default) or imperial.}

\item{mode}{Text string specifying the mode of transport. Can be
bicycling (default), walking, driving or transit}

\item{arrival_time}{Time of arrival in date format.}
}
\description{
Return travel network distances and time using the Google Maps API
}
\details{
Absent authorization, the google API is limited to a maximum of 100
simultaneous queries, and so will, for example, only returns values for up to
10 origins times 10 destinations.
}
\section{Details}{

Estimate travel times accounting for the road network - see \url{https://developers.google.com/maps/documentation/distance-matrix/overview}
Note: Currently returns the json object returned by the Google Maps API and uses the same origins and destinations.
}

\examples{
\dontrun{
# Distances from one origin to one destination
from <- c(-46.3, -23.4)
to <- c(-46.4, -23.4)
dist_google(from = from, to = to, mode = "walking") # not supported on last test
dist_google(from = from, to = to, mode = "driving")
dist_google(from = c(0, 52), to = c(0, 53))
data("cents")
# Distances from between all origins and destinations
dists_cycle <- dist_google(from = cents, to = cents)
dists_drive <- dist_google(cents, cents, mode = "driving")
dists_trans <- dist_google(cents, cents, mode = "transit")
dists_trans_am <- dist_google(cents, cents,
  mode = "transit",
  arrival_time = strptime("2016-05-27 09:00:00",
    format = "\%Y-\%m-\%d \%H:\%M:\%S", tz = "BST"
  )
)
# Find out how much longer (or shorter) cycling takes than walking
summary(dists_cycle$duration / dists_trans$duration)
# Difference between travelling now and for 9am arrival
summary(dists_trans_am$duration / dists_trans$duration)
odf <- points2odf(cents)
odf <- cbind(odf, dists)
head(odf)
flow <- points2flow(cents)
# show the results for duration (thicker line = shorter)
plot(flow, lwd = mean(odf$duration) / odf$duration)
dist_google(c("Hereford"), c("Weobley", "Leominster", "Kington"))
dist_google(c("Hereford"), c("Weobley", "Leominster", "Kington"),
  mode = "transit", arrival_time = strptime("2016-05-27 17:30:00",
    format = "\%Y-\%m-\%d \%H:\%M:\%S", tz = "BST"
  )
)
}
}
\seealso{
Other od: 
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overline.R
\name{lineLabels}
\alias{lineLabels}
\title{Label SpatialLinesDataFrame objects}
\usage{
lineLabels(sl, attrib)
}
\arguments{
\item{sl}{A SpatialLinesDataFrame with overlapping elements}

\item{attrib}{A text string corresponding to a named variable in \code{sl}}
}
\description{
This function adds labels to lines plotted using base graphics. Largely
for illustrative purposes, not designed for publication-quality
graphics.
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\author{
Barry Rowlingson
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overline.R
\name{overline_intersection}
\alias{overline_intersection}
\title{Convert series of overlapping lines into a route network}
\usage{
overline_intersection(sl, attrib, fun = sum)
}
\arguments{
\item{sl}{An \code{sf} \code{LINESTRING} object with overlapping elements}

\item{attrib}{character, column names in sl to be aggregated}

\item{fun}{Named list of functions to summaries the attributes by? \code{sum} is the default.
\code{list(sum = sum, average = mean)} will summarise all \code{attrib}utes by sum and mean.}
}
\description{
This function takes overlapping \code{LINESTRING}s stored in an
\code{sf} object and returns a route network composed of non-overlapping
geometries and aggregated values.
}
\examples{
routes_fast_sf$value <- 1
sl <- routes_fast_sf[4:6, ]
attrib <- c("value", "length")
rnet <- overline_intersection(sl = sl, attrib)
plot(rnet, lwd = rnet$value)
# A larger example
sl <- routes_fast_sf[4:7, ]
rnet <- overline_intersection(sl = sl, attrib = c("value", "length"))
plot(rnet, lwd = rnet$value)
rnet_sf <- overline(routes_fast_sf[4:7, ], attrib = c("value", "length"))
plot(rnet_sf, lwd = rnet_sf$value)

# An even larger example (not shown, takes time to run)
# rnet = overline_intersection(routes_fast_sf, attrib = c("value", "length"))
# rnet_sf <- overline(routes_fast_sf, attrib = c("value", "length"), buff_dist = 10)
# plot(rnet$geometry, lwd = rnet$value * 2, col = "grey")
# plot(rnet_sf$geometry,  lwd = rnet_sf$value, add = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route.R
\name{route_dodgr}
\alias{route_dodgr}
\title{Route on local data using the dodgr package}
\usage{
route_dodgr(from = NULL, to = NULL, l = NULL, net = NULL)
}
\arguments{
\item{from}{An object representing origins
(if lines are provided as the first argument, from is assigned to \code{l})}

\item{to}{An object representing destinations}

\item{l}{Only needed if from and to are empty, in which case this
should be a spatial object representing desire lines}

\item{net}{sf object representing the route network}
}
\description{
Route on local data using the dodgr package
}
\examples{
if (requireNamespace("dodgr")) {
  from <- c(-1.5327, 53.8006) # from <- geo_code("pedallers arms leeds")
  to <- c(-1.5279, 53.8044) # to <- geo_code("gzing")
  # next 4 lines were used to generate `stplanr::osm_net_example`
  # pts <- rbind(from, to)
  # colnames(pts) <- c("X", "Y")
  # net <- dodgr::dodgr_streetnet(pts = pts, expand = 0.1)
  # osm_net_example <- net[c("highway", "name", "lanes", "maxspeed")]
  r <- route_dodgr(from, to, net = osm_net_example)
  plot(osm_net_example$geometry)
  plot(r$geometry, add = TRUE, col = "red", lwd = 5)
}
}
\seealso{
Other routes: 
\code{\link{line2routeRetry}()},
\code{\link{line2route}()},
\code{\link{route_local}()},
\code{\link{route_osrm}()},
\code{\link{route_transportapi_public}()},
\code{\link{route}()}
}
\concept{routes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/node-funs.R
\name{rnet_add_node}
\alias{rnet_add_node}
\title{Add a node to route network}
\usage{
rnet_add_node(rnet, p)
}
\arguments{
\item{rnet}{A route network of the type generated by \code{overline()}}

\item{p}{A point represented by an \code{sf} object the will split the \code{route}}
}
\description{
Add a node to route network
}
\examples{
sample_routes <- routes_fast_sf[2:6, NULL]
sample_routes$value <- rep(1:3, length.out = 5)
rnet <- overline2(sample_routes, attrib = "value")
p <- sf::st_sfc(sf::st_point(c(-1.540, 53.826)), crs = sf::st_crs(rnet))
r_split <- route_split(rnet, p)
plot(rnet$geometry, lwd = rnet$value * 5, col = "grey")
plot(p, cex = 9, add = TRUE)
plot(r_split, col = 1:nrow(r_split), add = TRUE, lwd = r_split$value)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catchmentArea.R
\name{calc_catchment}
\alias{calc_catchment}
\title{Calculate catchment area and associated summary statistics.}
\usage{
calc_catchment(
  polygonlayer,
  targetlayer,
  calccols,
  distance = 500,
  projection = paste0("+proj=aea +lat_1=90 +lat_2=-18.416667 ",
    "+lat_0=0 +lon_0=10 +x_0=0 +y_0=0 +ellps=GRS80",
    " +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"),
  retainAreaProportion = FALSE,
  dissolve = FALSE,
  quadsegs = NULL
)
}
\arguments{
\item{polygonlayer}{A SpatialPolygonsDataFrame containing zones from which
the summary statistics for the catchment variable will be calculated.
Smaller polygons will increase the accuracy of the results.}

\item{targetlayer}{A SpatialPolygonsDataFrame, SpatialLinesDataFrame,
SpatialPointsDataFrame, SpatialPolygons, SpatialLines or SpatialPoints
object containing the specifications of the facility for which the
catchment area is being calculated. If the object contains more than one
facility (e.g., multiple cycle paths) the aggregate catchment area will be
calculated.}

\item{calccols}{A vector of column names containing the variables in the
polygonlayer to be used in the calculation of the summary statistics for
the catchment area. If dissolve = FALSE, all other variables in the
original SpatialPolygonsDataFrame for zones that fall partly or entirely
within the catchment area will be included in the returned
SpatialPolygonsDataFrame but will not be adjusted for the proportion within
the catchment area.}

\item{distance}{Defines the size of the catchment area as the distance
around the targetlayer in the units of the projection
(default = 500 metres)}

\item{projection}{The proj4string used to define the projection to be used
for calculating the catchment areas or a character string 'austalbers' to
use the Australian Albers Equal Area projection. Ignored if the polygonlayer
is projected in which case the targetlayer will be converted to the
projection used by the polygonlayer. In all cases the resulting object will
be reprojected to the original coordinate system and projection of the
polygon layer. Default is an Albers Equal Area projection but for more
reliable results should use a local projection (e.g., Australian Albers
Equal Area project).}

\item{retainAreaProportion}{Boolean value. If TRUE retains a variable in
the resulting SpatialPolygonsDataFrame containing the proportion of the
original area within the catchment area (Default = FALSE).}

\item{dissolve}{Boolean value. If TRUE collapses the underlying zones
within the catchment area into a single region with statistics for the
whole catchment area.}

\item{quadsegs}{Number of line segments to use to approximate a quarter
circle. Parameter passed to buffer functions, default is 5 for sp and
30 for sf.}
}
\description{
Calculate catchment area and associated summary statistics.
}
\section{Details}{

Calculates the catchment area of a facility (e.g., cycle path) using
straight-line distance as well as summary statistics from variables
available in a SpatialPolygonsDataFrame with census tracts or other
zones. Assumes that the frequency of the variable is evenly distributed
throughout the zone. Returns a SpatialPolygonsDataFrame.
}

\examples{
\dontrun{
data_dir <- system.file("extdata", package = "stplanr")
unzip(file.path(data_dir, "smallsa1.zip"))
unzip(file.path(data_dir, "testcycleway.zip"))
sa1income <- as(sf::read_sf("smallsa1.shp"), "Spatial")
testcycleway <- as(sf::read_sf("testcycleway.shp"), "Spatial")
cway_catch <- calc_catchment(
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  distance = 800,
  projection = "austalbers",
  dissolve = TRUE
)
plot(sa1income)
plot(cway_catch, add = TRUE, col = "green")
plot(testcycleway, col = "red", add = TRUE)
sa1income <- sf::read_sf("smallsa1.shp")
testcycleway <- sf::read_sf("testcycleway.shp")
f <- list.files(".", "testcycleway|smallsa1")
file.remove(f)
cway_catch <- calc_catchment(
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  distance = 800,
  projection = "austalbers",
  dissolve = TRUE
)
plot(sa1income$geometry)
plot(testcycleway$geometry, col = "red", add = TRUE)
plot(cway_catch["Total"], add = TRUE)
}
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnet_group.R
\name{rnet_group}
\alias{rnet_group}
\alias{rnet_group.default}
\alias{rnet_group.sfc}
\alias{rnet_group.sf}
\alias{rnet_group.sfNetwork}
\title{Assign segments in a route network to groups}
\usage{
rnet_group(rnet, ...)

\method{rnet_group}{default}(rnet, ...)

\method{rnet_group}{sfc}(
  rnet,
  cluster_fun = igraph::clusters,
  d = NULL,
  as.undirected = TRUE,
  ...
)

\method{rnet_group}{sf}(
  rnet,
  cluster_fun = igraph::clusters,
  d = NULL,
  as.undirected = TRUE,
  ...
)

\method{rnet_group}{sfNetwork}(rnet, cluster_fun = igraph::clusters, ...)
}
\arguments{
\item{rnet}{An sf, sfc, or sfNetwork object representing a route network.}

\item{...}{Arguments passed to other methods.}

\item{cluster_fun}{The clustering function to use. Various clustering functions
are available in the \code{igraph} package. Default: \code{\link[igraph:components]{igraph::clusters()}}.}

\item{d}{Optional distance variable used to classify segments that are
close (within a certain distance specified by \code{d}) to each other but not
necessarily touching}

\item{as.undirected}{Coerce the graph created internally into an undirected
graph with \code{\link[igraph:as.directed]{igraph::as.undirected()}}? TRUE by default, which enables use
of a wider range of clutering functions.}
}
\value{
If the input rnet is an sf/sfc object, it returns an integer vector
reporting the groups of each network element. If the input is an sfNetwork
object, it returns an sfNetwork object with an extra column called
rnet_group representing the groups of each network element. In the latter
case, the connectivity of the spatial object is derived from the sfNetwork
object.
}
\description{
This function assigns linestring features, many of which in an
\code{sf} object can form route networks, into groups.
By default, the function \code{igraph::clusters()} is used to determine
group membership, but any \verb{igraph::cluster*()} function can be used.
See examples and the web page
\href{https://igraph.org/r/doc/communities.html}{igraph.org/r/doc/communities.html}
for more information. From that web page, the following clustering
functions are available:
}
\details{
\verb{cluster_edge_betweenness, cluster_fast_greedy, cluster_label_prop,}
\verb{cluster_leading_eigen, cluster_louvain, cluster_optimal, cluster_spinglass, cluster_walktrap}
}
\examples{
rnet <- rnet_breakup_vertices(stplanr::osm_net_example)
rnet$group <- rnet_group(rnet)
plot(rnet["group"])
# mapview::mapview(rnet["group"])
rnet$group_25m <- rnet_group(rnet, d = 25)
plot(rnet["group_25m"])
rnet$group_walktrap <- rnet_group(rnet, igraph::cluster_walktrap)
plot(rnet["group_walktrap"])
rnet$group_louvain <- rnet_group(rnet, igraph::cluster_louvain)
plot(rnet["group_louvain"])
rnet$group_fast_greedy <- rnet_group(rnet, igraph::cluster_fast_greedy)
plot(rnet["group_fast_greedy"])

# show sfNetwork implementation
sfn <- SpatialLinesNetwork(rnet)
sfn <- rnet_group(sfn)
plot(sfn@sl["rnet_group"])
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route.R
\name{route}
\alias{route}
\title{Plan routes on the transport network}
\usage{
route(
  from = NULL,
  to = NULL,
  l = NULL,
  route_fun = cyclestreets::journey,
  wait = 0,
  n_print = 10,
  list_output = FALSE,
  cl = NULL,
  ...
)
}
\arguments{
\item{from}{An object representing origins
(if lines are provided as the first argument, from is assigned to \code{l})}

\item{to}{An object representing destinations}

\item{l}{Only needed if from and to are empty, in which case this
should be a spatial object representing desire lines}

\item{route_fun}{A routing function to be used for converting the straight lines to routes
\code{\link[=od2line]{od2line()}}}

\item{wait}{How long to wait between routes?
0 seconds by default, can be useful when sending requests to rate limited APIs.}

\item{n_print}{A number specifying how frequently progress updates
should be shown}

\item{list_output}{If FALSE (default) assumes spatial (linestring) object output. Set to TRUE to save output as a list.}

\item{cl}{Cluster}

\item{...}{Arguments passed to the routing function, e.g. \code{\link[=route_cyclestreets]{route_cyclestreets()}}}
}
\description{
Takes origins and destinations, finds the optimal routes between them
and returns the result as a spatial (sf or sp) object.
The definition of optimal depends on the routing function used
}
\examples{
library(sf)
l = od_data_lines[2, ]
\donttest{
if(curl::has_internet()) {
r_walk = route(l = l, route_fun = route_osrm, osrm.profile = "foot")
r_bike = route(l = l, route_fun = route_osrm, osrm.profile = "bike")
plot(r_walk$geometry)
plot(r_bike$geometry, col = "blue", add = TRUE)
# r_bc = route(l = l, route_fun = route_bikecitizens)
# plot(r_bc)
# route(l = l, route_fun = route_bikecitizens, wait = 1)
library(osrm)
r_osrm <- route(
  l = l,
  route_fun = osrmRoute,
  returnclass = "sf"
)
nrow(r_osrm)
plot(r_osrm)
sln <- stplanr::SpatialLinesNetwork(route_network_sf)
# calculate shortest paths
plot(sln)
plot(l$geometry, add = TRUE)
r_local <- stplanr::route(
  l = l,
  route_fun = stplanr::route_local,
  sln = sln
)
plot(r_local["all"], add = TRUE, lwd = 5)
}
}
}
\seealso{
Other routes: 
\code{\link{line2routeRetry}()},
\code{\link{line2route}()},
\code{\link{route_dodgr}()},
\code{\link{route_local}()},
\code{\link{route_osrm}()},
\code{\link{route_transportapi_public}()}
}
\concept{routes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{SpatialLinesNetwork}
\alias{SpatialLinesNetwork}
\title{Create object of class SpatialLinesNetwork or sfNetwork}
\usage{
SpatialLinesNetwork(sl, uselonglat = FALSE, tolerance = 0)
}
\arguments{
\item{sl}{A SpatialLines or SpatialLinesDataFrame containing the lines to
use to create the network.}

\item{uselonglat}{A boolean value indicating if the data should be assumed
to be using WGS84 latitude/longitude coordinates. If \code{FALSE} or not
set, uses the coordinate system specified by the SpatialLines object.}

\item{tolerance}{A numeric value indicating the tolerance (in the units of
the coordinate system) to use as a tolerance with which to match nodes.}
}
\description{
Creates a new SpatialLinesNetwork (for SpatialLines) or sfNetwork (for sf)
object that can be used for routing analysis within R.
}
\section{Details}{

This function is used to create a new SpatialLinesNetwork from an existing
SpatialLines or SpatialLinesDataFrame object. A typical use case is to
represent a transport network for routing and other network analysis
functions. This function and the corresponding SpatialLinesNetwork
class is an implementation of the SpatialLinesNetwork developed by
Edzer Pebesma and presented on \href{https://rpubs.com/edzer/6767}{RPubs}.
The original implementation has been rewritten to better support large
(i.e., detailed city-size) networks and to provide additional methods
useful for conducting transport research following on from the initial
examples provided by \href{https://rpubs.com/janoskaz/10396}{Janoska(2013)}.
}

\examples{
\donttest{
# dont test due to issues with s2 dependency
sln_sf <- SpatialLinesNetwork(route_network_sf)
plot(sln_sf)
shortpath <- sum_network_routes(sln_sf, 1, 50, sumvars = "length")
plot(shortpath$geometry, col = "red", lwd = 4, add = TRUE)
}
}
\references{
Pebesma, E. (2013). Spatial Networks, URL:https://rpubs.com/edzer/6767.

Janoska, Z. (2013). Find shortest path in spatial network,
URL:https://rpubs.com/janoskaz/10396.
}
\seealso{
Other rnet: 
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_cyclestreets.R
\name{api_pat}
\alias{api_pat}
\title{Retrieve personal access token.}
\usage{
api_pat(api_name, force = FALSE)
}
\arguments{
\item{api_name}{Text string of the name of the API you are calling, e.g.
cyclestreets, graphhopper etc.}
}
\description{
Retrieve personal access token.
}
\examples{
\dontrun{
api_pat(api_name = "cyclestreet")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overline.R
\name{overline}
\alias{overline}
\alias{overline2}
\title{Convert series of overlapping lines into a route network}
\usage{
overline(
  sl,
  attrib,
  ncores = 1,
  simplify = TRUE,
  regionalise = 1e+05,
  quiet = ifelse(nrow(sl) < 1000, TRUE, FALSE),
  fun = sum
)

overline2(
  sl,
  attrib,
  ncores = 1,
  simplify = TRUE,
  regionalise = 1e+05,
  quiet = ifelse(nrow(sl) < 1000, TRUE, FALSE),
  fun = sum
)
}
\arguments{
\item{sl}{A spatial object representing routes on a transport network}

\item{attrib}{character, column names in sl to be aggregated}

\item{ncores}{integer, how many cores to use in parallel processing, default = 1}

\item{simplify}{logical, if TRUE group final segments back into lines, default = TRUE}

\item{regionalise}{integer, during simplification regonalisation is used if the number of segments exceeds this value}

\item{quiet}{Should the the function omit messages? \code{NULL} by default,
which means the output will only be shown if \code{sl} has more than 1000 rows.}

\item{fun}{Named list of functions to summaries the attributes by? \code{sum} is the default.
\code{list(sum = sum, average = mean)} will summarise all \code{attrib}utes by sum and mean.}
}
\value{
An \code{sf} object representing a route network
}
\description{
This function takes a series of overlapping lines and converts them into a
single route network.

This function is intended as a replacement for overline() and is significantly faster
especially on large datasets. However, it also uses more memory.
}
\details{
The function can be used to estimate the amount of transport 'flow' at the
route segment level based on input datasets from routing services, for
example linestring geometries created with the \code{route()} function.

The \code{overline()} function breaks each line into many straight
segments and then looks for duplicated segments. Attributes are summed for
all duplicated segments, and if simplify is TRUE the segments with identical
attributes are recombined into linestrings.

The following arguments only apply to the \code{sf} implementation of \code{overline()}:
\itemize{
\item \code{ncores}, the number of cores to use in parallel processing
\item \code{simplify}, should the final segments be converted back into longer lines? The default
setting is \code{TRUE}. \code{simplify = FALSE} results in straight line segments consisting
of only 2 vertices (the start and end point),
resulting in a data frame with many more rows than the simplified results (see examples).
\item \code{regionalise} the threshold number of rows above which
regionalisation is used (see details).
}

For \code{sf} objects Regionalisation breaks the dataset into a 10 x 10 grid and
then performed the simplification across each grid. This significantly
reduces computation time for large datasets, but slightly increases the final
file size. For smaller datasets it increases computation time slightly but
reduces memory usage and so may also be useful.

A known limitation of this method is that overlapping segments of different
lengths are not aggregated. This can occur when lines stop halfway down a
road. Typically these errors are small, but some artefacts may remain within
the resulting data.

For very large datasets nrow(x) > 1000000, memory usage can be significant.
In these cases is is possible to overline subsets of the dataset, rbind the
results together, and then overline again, to produce a final result.

Multicore support is only enabled for the regionalised simplification stage
as it does not help with other stages.
}
\examples{
sl <- routes_fast_sf[2:4, ]
sl$All <- flowlines$All[2:4]
rnet <- overline(sl = sl, attrib = "All")
nrow(sl)
nrow(rnet)
plot(rnet)
rnet_mean <- overline(sl, c("All", "av_incline"), fun = list(mean = mean, sum = sum))
plot(rnet_mean, lwd = rnet_mean$All_sum / mean(rnet_mean$All_sum))
rnet_sf_raw <- overline(sl, attrib = "length", simplify = FALSE)
nrow(rnet_sf_raw)
summary(n_vertices(rnet_sf_raw))
plot(rnet_sf_raw)
rnet_sf_raw$n <- 1:nrow(rnet_sf_raw)
plot(rnet_sf_raw[10:25, ])
# legacy implementation based on sp data
# sl <- routes_fast[2:4, ]
# rnet1 <- overline(sl = sl, attrib = "length")
# rnet2 <- overline(sl = sl, attrib = "length", buff_dist = 1)
# plot(rnet1, lwd = rnet1$length / mean(rnet1$length))
# plot(rnet2, lwd = rnet2$length / mean(rnet2$length))
}
\references{
Morgan M and Lovelace R (2020). Travel flow aggregation: Nationally scalable methods
for interactive and online visualisation of transport behaviour at the road network level.
Environment and Planning B: Urban Analytics and City Science. July 2020.
\doi{10.1177/2399808320942779}.

Rowlingson, B (2015). Overlaying lines and aggregating their values for overlapping
segments. Reproducible question from \url{https://gis.stackexchange.com}. See
\url{https://gis.stackexchange.com/questions/139681/}.
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}

Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\author{
Barry Rowlingson

Malcolm Morgan
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_projected.R
\name{geo_length}
\alias{geo_length}
\title{Calculate line length of line with geographic or projected CRS}
\usage{
geo_length(shp)
}
\arguments{
\item{shp}{A spatial line object}
}
\description{
Takes a line (represented in sf or sp classes)
and returns a numeric value representing distance in meters.
}
\examples{
lib_versions <- sf::sf_extSoftVersion()
lib_versions
if (lib_versions[3] >= "6.3.1") {
  geo_length(routes_fast)
  geo_length(routes_fast_sf)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od_coords2line}
\alias{od_coords2line}
\title{Convert origin-destination coordinates into desire lines}
\usage{
od_coords2line(odc, crs = 4326, remove_duplicates = TRUE)
}
\arguments{
\item{odc}{A data frame or matrix representing the coordinates
of origin-destination data. The first two columns represent the
coordinates of the origin (typically longitude and latitude) points;
the third and fourth columns represent the coordinates of the destination
(in the same CRS). Each row represents travel from origin to destination.}

\item{crs}{A number representing the coordinate reference system
of the result, 4326 by default.}

\item{remove_duplicates}{Should rows with duplicated rows be removed? \code{TRUE} by default.}
}
\description{
Convert origin-destination coordinates into desire lines
}
\examples{
odf <- od_coords(l = flowlines_sf)
odlines <- od_coords2line(odf)
odlines <- od_coords2line(odf, crs = 4326)
plot(odlines)
x_coords <- 1:3
n <- 50
d <- data.frame(lapply(1:4, function(x) sample(x_coords, n, replace = TRUE)))
names(d) <- c("fx", "fy", "tx", "ty")
l <- od_coords2line(d)
plot(l)
nrow(l)
l_with_duplicates <- od_coords2line(d, remove_duplicates = FALSE)
plot(l_with_duplicates)
nrow(l_with_duplicates)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_google.R
\name{route_google}
\alias{route_google}
\title{Find shortest path using Google services}
\usage{
route_google(from, to, mode = "walking", key = Sys.getenv("GOOGLE"), ...)
}
\arguments{
\item{from}{An object representing origins
(if lines are provided as the first argument, from is assigned to \code{l})}

\item{to}{An object representing destinations}

\item{mode}{Mode of transport, walking (default), bicycling, transit, or driving}

\item{key}{Google key. By default it is \code{Sys.getenv("GOOGLE")}. Set it with:
\code{usethis::edit_r_environ()}.}

\item{...}{Arguments passed to the routing function, e.g. \code{\link[=route_cyclestreets]{route_cyclestreets()}}}
}
\description{
Find the shortest path using Google's services.
See the \code{mapsapi} package for details.
}
\examples{
\dontrun{
from <- "university of leeds"
to <- "pedallers arms leeds"
r <- route(from, to, route_fun = cyclestreets::journey)
plot(r)
# r_google <- route(from, to, route_fun = mapsapi::mp_directions) # fails
r_google1 <- route_google(from, to)
plot(r_google1)
r_google <- route(from, to, route_fun = route_google)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od_aggregate_to}
\alias{od_aggregate_to}
\title{Summary statistics of trips arriving at destination zones in OD data}
\usage{
od_aggregate_to(flow, attrib = NULL, FUN = sum, ..., col = 2)
}
\arguments{
\item{flow}{A data frame representing origin-destination data.
The first two columns of this data frame should correspond
to the first column of the data in the zones. Thus in \code{\link[=cents]{cents()}},
the first column is geo_code. This corresponds to the first two columns
of \code{\link[=flow]{flow()}}.}

\item{attrib}{character, column names in sl to be aggregated}

\item{FUN}{A function to summarise OD data by}

\item{...}{Additional arguments passed to \code{FUN}}

\item{col}{The column that the OD dataset is grouped by
(1 by default, the first column usually represents the origin)}
}
\description{
This function takes a data frame of OD data and
returns a data frame reporting summary statistics for each unique zone of destination.
}
\details{
It has some default settings: it assumes the destination ID column is the 2nd
and the default summary statistic is \code{sum()}.
By default, if \code{attrib} is not set, it summarises all numeric columns.
}
\examples{
od_aggregate_to(flow)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnet-clean.R
\name{rnet_breakup_vertices}
\alias{rnet_breakup_vertices}
\title{Break up an sf object with LINESTRING geometry.}
\usage{
rnet_breakup_vertices(rnet, verbose = FALSE)
}
\arguments{
\item{rnet}{An sf or sfc object with LINESTRING geometry representing a route
network.}

\item{verbose}{Boolean. If TRUE, the function prints additional messages.}
}
\value{
An sf or sfc object with LINESTRING geometry created after breaking
up the input object.
}
\description{
This function breaks up a LINESTRING geometry into multiple LINESTRING(s). It
is used mainly for preserving routability of an \code{sfNetwork} object that is
created using Open Street Map data. See details,
\href{https://github.com/ropensci/stplanr/issues/282}{stplanr/issues/282}, and
\href{https://github.com/ropensci/stplanr/issues/416}{stplanr/issues/416}.
}
\details{
A LINESTRING geometry is broken-up when one of the two following conditions
are met:
\enumerate{
\item two or more LINESTRINGS share a POINT which is a boundary point for some
LINESTRING(s), but not all of them (see the rnet_roundabout example);
\item two or more LINESTRINGS share a POINT which is not in the boundary of any
LINESTRING (see the rnet_cycleway_intersection example).
}

The problem with the first example is that, according to algorithm behind
\code{\link[=SpatialLinesNetwork]{SpatialLinesNetwork()}}, two LINESTRINGS are connected if and only if they
share at least one point in their boundaries. The roads and the roundabout
are clearly connected in the "real" world but the corresponding LINESTRING
objects do not share two distinct boundary points. In fact, by Open Street
Map standards, a roundabout is represented as a closed and circular
LINESTRING, and this implies that the roundabout is not connected to the
other roads according to \code{\link[=SpatialLinesNetwork]{SpatialLinesNetwork()}} definition. By the same
reasoning, the roads in the second example are clearly connected in the
"real" world, but they do not share any point in their boundaries. This
function is used to solve this type of problem.
}
\examples{
library(sf)
def_par <- par(no.readonly = TRUE)
par(mar = rep(0, 4))

# Check the geometry of the roundabout example. The dots represent the
# boundary points of the LINESTRINGS. The "isolated" red point in the
# top-left is the boundary point of the roundabout, and it is not shared
# with any other street.
plot(st_geometry(rnet_roundabout), lwd = 2, col = rainbow(nrow(rnet_roundabout)))
boundary_points <- st_geometry(line2points(rnet_roundabout))
points_cols <- rep(rainbow(nrow(rnet_roundabout)), each = 2)
plot(boundary_points, pch = 16, add = TRUE, col = points_cols, cex = 2)

# Clean the roundabout example.
rnet_roundabout_clean <- rnet_breakup_vertices(rnet_roundabout)
plot(st_geometry(rnet_roundabout_clean), lwd = 2, col = rainbow(nrow(rnet_roundabout_clean)))
boundary_points <- st_geometry(line2points(rnet_roundabout_clean))
points_cols <- rep(rainbow(nrow(rnet_roundabout_clean)), each = 2)
plot(boundary_points, pch = 16, add = TRUE, col = points_cols)
# The roundabout is now routable since it was divided into multiple pieces
# (one for each colour), which, according to SpatialLinesNetwork() function,
# are connected.

# Check the geometry of the overpasses example. This example is used to test
# that this function does not create any spurious intersection.
plot(st_geometry(rnet_overpass), lwd = 2, col = rainbow(nrow(rnet_overpass)))
boundary_points <- st_geometry(line2points(rnet_overpass))
points_cols <- rep(rainbow(nrow(rnet_overpass)), each = 2)
plot(boundary_points, pch = 16, add = TRUE, col = points_cols, cex = 2)
# At the moment the network is not routable since one of the underpasses is
# not connected to the other streets.

# Check interactively.
# mapview::mapview(rnet_overpass)

# Clean the network. It should not create any spurious intersection between
# roads located at different heights.
rnet_overpass_clean <- rnet_breakup_vertices(rnet_overpass)
plot(st_geometry(rnet_overpass_clean), lwd = 2, col = rainbow(nrow(rnet_overpass_clean)))
# Check interactively.
# mapview::mapview(rnet_overpass)

# Check the geometry of the cycleway_intersection example. The black dots
# represent the boundary points and we can see that the two roads are not
# connected according to SpatialLinesNetwork() function.
plot(
  rnet_cycleway_intersection$geometry,
  lwd = 2,
  col = rainbow(nrow(rnet_cycleway_intersection)),
  cex = 2
)
plot(st_geometry(line2points(rnet_cycleway_intersection)), pch = 16, add = TRUE)
# Check interactively
# mapview::mapview(rnet_overpass)

# Clean the rnet object and plot the result.
rnet_cycleway_intersection_clean <- rnet_breakup_vertices(rnet_cycleway_intersection)
plot(
  rnet_cycleway_intersection_clean$geometry,
  lwd = 2,
  col = rainbow(nrow(rnet_cycleway_intersection_clean)),
  cex = 2
)
plot(st_geometry(line2points(rnet_cycleway_intersection_clean)), pch = 16, add = TRUE)

par(def_par)
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cents}
\alias{cents}
\alias{cents_sf}
\title{Spatial points representing home locations}
\format{
A spatial dataset with 8 rows and 5 variables
}
\usage{
data(cents)
}
\description{
These points represent population-weighted centroids of Medium Super Output Area (MSOA) zones within a 1 mile radius of of my home when I was writing this package.
}
\details{
\itemize{
\item geo_code the official code of the zone
\item MSOA11NM name zone name
\item percent_fem the percent female
\item avslope average gradient of the zone
}

Cents was generated from the data repository pct-data: https://github.com/npct/pct-data. This data was accessed from within the pct repo: https://github.com/npct/pct, using the following code:
}
\examples{
\dontrun{
cents
plot(cents)
}

}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{route_network}
\alias{route_network}
\alias{route_network_sf}
\title{spatial lines dataset representing a route network}
\format{
A spatial lines dataset 80 rows and 1 column
}
\usage{
data(route_network)
}
\description{
The flow of commuters using different segments of the road network represented in the
\code{\link[=flowlines]{flowlines()}} and \code{\link[=routes_fast]{routes_fast()}} datasets
}
\examples{
\dontrun{
# Generate route network
route_network <- overline(routes_fast, "All", fun = sum)
route_network_sf <- sf::st_as_sf(route_network)
}
}
\seealso{
Other example data: 
\code{\link{destination_zones}},
\code{\link{flow_dests}},
\code{\link{flowlines}},
\code{\link{flow}},
\code{\link{routes_fast}},
\code{\link{routes_slow}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od2odf}
\alias{od2odf}
\title{Extract coordinates from OD data}
\usage{
od2odf(flow, zones)
}
\arguments{
\item{flow}{A data frame representing origin-destination data.
The first two columns of this data frame should correspond
to the first column of the data in the zones. Thus in \code{\link[=cents]{cents()}},
the first column is geo_code. This corresponds to the first two columns
of \code{\link[=flow]{flow()}}.}

\item{zones}{A spatial object representing origins (and destinations
if no separate destinations object is provided) of travel.}
}
\description{
Extract coordinates from OD data
}
\details{
Origin-destination (OD) data is often provided
in the form of 1 line per OD pair, with zone codes of the trip origin in the first
column and the zone codes of the destination in the second column
(see the \href{https://docs.ropensci.org/stplanr/articles/stplanr-od.html}{\code{vignette("stplanr-od")}}) for details.
\code{od2odf()} creates an 'origin-destination data frame', based on a data frame containing
origin and destination cones (\code{flow}) that match the first column in a
a spatial (polygon or point) object (\code{zones}).

The function returns a data frame with coordinates for the origin and destination.
}
\examples{
data(flow)
data(zones)
od2odf(flow[1:2, ], zones)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_code.R
\name{geo_code}
\alias{geo_code}
\title{Convert text strings into points on the map}
\usage{
geo_code(
  address,
  service = "nominatim",
  base_url = "https://maps.google.com/maps/api/geocode/json",
  return_all = FALSE,
  pat = NULL
)
}
\arguments{
\item{address}{Text string representing the address you want to geocode}

\item{service}{Which service to use? Nominatim by default}

\item{base_url}{The base url to query}

\item{return_all}{Should the request return all information returned by Google Maps?
The default is \code{FALSE}: to return only two numbers: the longitude and latitude, in that order}

\item{pat}{The API key used. By default this is set to NULL and
this is usually aquired automatically through a helper, api_pat().}
}
\description{
Generate a lat/long pair from data using Google's geolocation API.
}
\examples{
\dontrun{
geo_code(address = "Hereford")
geo_code("LS7 3HB")
geo_code("hereford", return_all = TRUE)
# needs api key in .Renviron
geo_code("hereford", service = "google", pat = Sys.getenv("GOOGLE"), return_all = TRUE)
}
}
\seealso{
Other nodes: 
\code{\link{nearest_google}()}
}
\concept{nodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_sf_fun.R
\name{as_sf_fun}
\alias{as_sf_fun}
\alias{as_sp_fun}
\title{Convert functions support sf/sp}
\usage{
as_sf_fun(input, FUN, ...)
}
\arguments{
\item{input}{Input object - an sf or sp object}

\item{FUN}{A function that works on sp/sf data}

\item{...}{Arguments passed to \code{FUN}}
}
\description{
Convert functions support sf/sp
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/line_breakup.R
\name{line_breakup}
\alias{line_breakup}
\title{Break up line objects into shorter segments}
\usage{
line_breakup(l, z)
}
\arguments{
\item{l}{An sf object with LINESTRING geometry}

\item{z}{An sf object with \code{POLYGON} geometry or a number representing the
resolution of grid cells used to break up the linestring objects}
}
\value{
An sf object with LINESTRING geometry created after breaking up the
input object.
}
\description{
This function breaks up a LINESTRING geometries into smaller pieces.
}
\examples{
library(sf)
z <- zones_sf$geometry
l <- routes_fast_sf$geometry[2]
l_split <- line_breakup(l, z)
l
l_split
sf::st_length(l)
sum(sf::st_length(l_split))
plot(z)
plot(l, add = TRUE, lwd = 9, col = "grey")
plot(l_split, add = TRUE, col = 1:length(l_split))
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/line_sample.R
\name{n_sample_length}
\alias{n_sample_length}
\title{Sample integer number from given continuous vector of line lengths and probabilities, with total n}
\usage{
n_sample_length(n, l_lengths, weights)
}
\arguments{
\item{n}{Sum of integer values returned}

\item{l_lengths}{Numeric vector of line lengths}

\item{weights}{Relative probabilities of samples on lines}
}
\description{
Sample integer number from given continuous vector of line lengths and probabilities, with total n
}
\examples{
n <- 10
l_lengths <- 1:5
weights <- 9:5
(res <- n_sample_length(n, l_lengths, weights))
sum(res)
n <- 100
l_lengths <- c(12, 22, 15, 14)
weights <- c(38, 10, 44, 34)
(res <- n_sample_length(n, l_lengths, weights))
sum(res)
# more examples:
n_sample_length(5, 1:5, c(0.1, 0.9, 0, 0, 0))
n_sample_length(5, 1:5, c(0.5, 0.3, 0.1, 0, 0))
l <- flowlines[2:6, ]
l_lengths <- line_length(l)
n <- n_sample_length(10, l_lengths, weights = l$All)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quadrants.R
\name{quadrant}
\alias{quadrant}
\title{Split a spatial object into quadrants}
\usage{
quadrant(sp_obj, number_out = FALSE)
}
\arguments{
\item{sp_obj}{Spatial object}

\item{number_out}{Should the output be numbers from 1:4 (FALSE by default)}
}
\description{
Split a spatial object (initially tested on SpatialPolygons) into quadrants.
}
\details{
Returns a character vector of NE, SE, SW, NW corresponding to north-east, south-east
quadrants respectively. If number_out is TRUE, returns numbers from 1:4, respectively.
}
\examples{
data(zones)
sp_obj <- zones
(quads <- quadrant(sp_obj))
plot(sp_obj, col = factor(quads))
points(rgeos::gCentroid(sp_obj), col = "white")
# edge cases (e.g. when using rasters) lead to NAs
sp_obj <- raster::rasterToPolygons(raster::raster(ncol = 3, nrow = 3))
(quads <- quadrant(sp_obj))
plot(sp_obj, col = factor(quads))
}
\seealso{
Other geo: 
\code{\link{bbox_scale}()},
\code{\link{geo_bb_matrix}()},
\code{\link{geo_bb}()},
\code{\link{reproject}()}
}
\concept{geo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{line2points}
\alias{line2points}
\alias{line2pointsn}
\alias{line2vertices}
\title{Convert a spatial (linestring) object to points}
\usage{
line2points(l, ids = rep(1:nrow(l)))

line2pointsn(l)

line2vertices(l)
}
\arguments{
\item{l}{An \code{sf} object or a \code{SpatialLinesDataFrame} from the older \code{sp} package}

\item{ids}{Vector of ids (by default \code{1:nrow(l)})}
}
\description{
The number of points will be double the number of lines with \code{line2points}. A
closely related function, \code{line2pointsn} returns all the points that were
line vertices. The points corresponding with a given line, \code{i}, will be
\code{(2*i):((2*i)+1)}. The last function, \code{line2vertices}, returns all the points
that are vertices but not nodes. If the input \code{l} object is composed by only
1 LINESTRING with 2 POINTS, then it returns an empty \code{sf} object.
}
\examples{
l <- routes_fast_sf[2, ]
lpoints <- line2points(l)
plot(l$geometry)
plot(lpoints, add = TRUE)
# test all vertices:
plot(l$geometry)
lpoints2 <- line2pointsn(l)
plot(lpoints2$geometry, add = TRUE)

# extract only internal vertices
l_internal_vertices <- line2vertices(l)
plot(sf::st_geometry(l), reset = FALSE)
plot(l_internal_vertices, add = TRUE)
# The boundary points are missing
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{find_network_nodes}
\alias{find_network_nodes}
\title{Find graph node ID of closest node to given coordinates}
\usage{
find_network_nodes(sln, x, y = NULL, maxdist = 1000)
}
\arguments{
\item{sln}{SpatialLinesNetwork to search.}

\item{x}{Either the x (longitude) coordinate value, a vector of x values,
a dataframe or matrix with (at least) two columns, the first for coordinate
for x (longitude) values and a second for y (latitude) values, or a named
vector of length two with values of 'lat' and 'lon'. The output of
geo_code() either as a single result or as multiple (using
rbind() ) can also be used.}

\item{y}{Either the y (latitude) coordinate value or a vector of y values.}

\item{maxdist}{The maximum distance within which to match the nodes to
coordinates. If the SpatialLinesNetwork is projected then distance should
be in the same units as the projection. If longlat, then distance is in
metres. Default is 1000.}
}
\value{
An integer value with the ID of the node closest to \verb{(x,y)}
with a value of \code{NA} the closest node is further than \code{maxdist}
from \verb{(x,y)}. If \code{x} is a vector, returns a vector of Node IDs.
}
\description{
Find graph node ID of closest node to given coordinates
}
\section{Details}{

Finds the node ID of the closest point to a single coordinate pair (or a
set of coordinates) from a SpatialLinesNetwork.
}

\examples{
data(routes_fast)
rnet <- overline(routes_fast, attrib = "length")
sln <- SpatialLinesNetwork(rnet)
find_network_nodes(sln, -1.516734, 53.828)
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_funs.R
\name{route_rolling_average}
\alias{route_rolling_average}
\title{Return smoothed averages of vector}
\usage{
route_rolling_average(x, n = 3)
}
\arguments{
\item{x}{Numeric vector to smooth}

\item{n}{The window size of the smoothing function.
The default, 3, will take the mean of values before, after and including
each value.}
}
\description{
This function calculates a simple rolling mean in base R.
It is useful for calculating route characteristics such as mean
distances of segments and changes in gradient.
}
\examples{
y <- od_data_routes$elevations[od_data_routes$route_number == 2]
y
route_rolling_average(y)
route_rolling_average(y, n = 1)
route_rolling_average(y, n = 2)
route_rolling_average(y, n = 3)
}
\seealso{
Other route_funs: 
\code{\link{route_average_gradient}()},
\code{\link{route_rolling_diff}()},
\code{\link{route_rolling_gradient}()},
\code{\link{route_sequential_dist}()},
\code{\link{route_slope_matrix}()},
\code{\link{route_slope_vector}()}
}
\concept{route_funs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\docType{class}
\name{SpatialLinesNetwork-class}
\alias{SpatialLinesNetwork-class}
\title{An S4 class representing a (typically) transport network}
\description{
This class uses a combination of a SpatialLinesDataFrame and an igraph
object to represent transport networks that can be used for routing and
other network analyses.
}
\section{Slots}{

\describe{
\item{\code{sl}}{A SpatialLinesDataFrame with the geometry and other attributes
for each link the in network.}

\item{\code{g}}{The graph network corresponding to \code{sl}.}

\item{\code{nb}}{A list containing vectors of the nodes connected to each node
in the network.}

\item{\code{weightfield}}{A character vector containing the variable (column) name
from the SpatialLinesDataFrame to be used for weighting the network.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overline.R
\name{overline_spatial}
\alias{overline_spatial}
\title{Spatial aggregation of routes represented with sp classes}
\usage{
overline_spatial(sl, attrib, fun = sum, na.zero = FALSE, buff_dist = 0)
}
\arguments{
\item{sl}{SpatialLinesDataFrame with overlapping Lines to split by
number of overlapping features.}

\item{attrib}{character, column names in sl to be aggregated}

\item{fun}{Named list of functions to summaries the attributes by? \code{sum} is the default.
\code{list(sum = sum, average = mean)} will summarise all \code{attrib}utes by sum and mean.}

\item{na.zero}{Sets whether aggregated values with a value of zero are
removed.}

\item{buff_dist}{A number specifying the distance in meters of the buffer to be used to crop lines before running the operation.
If the distance is zero (the default) touching but non-overlapping lines may be aggregated.}
}
\description{
This function, largely superseded by sf implementations, still works
but is not particularly fast.
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_bikecitizens.R
\name{route_bikecitizens}
\alias{route_bikecitizens}
\title{Get a route from the BikeCitizens web service}
\usage{
route_bikecitizens(
  from = NULL,
  to = NULL,
  base_url = "https://map.bikecitizens.net/api/v1/locations/route.json",
  cccode = "gb-leeds",
  routing_profile = "balanced",
  bike_profile = "citybike",
  from_lat = 53.8265,
  from_lon = -1.576195,
  to_lat = 53.80025,
  to_lon = -1.51577
)
}
\arguments{
\item{from}{A numeric vector representing the start point}

\item{to}{A numeric vector representing the end point}

\item{base_url}{The base URL for the routes}

\item{cccode}{The city code for the routes}

\item{routing_profile}{What type of routing to use?}

\item{bike_profile}{What type of bike?}

\item{from_lat}{Latitude of origin}

\item{from_lon}{Longitude of origin}

\item{to_lat}{Latitude of destination}

\item{to_lon}{Longitude of destination}
}
\description{
See \href{https://map.bikecitizens.net/gb-leeds#/!/1/1/53.8265,-1.576195/53.80025,-1.51577}{bikecitizens.net}
for an interactive version of the routing engine used by BikeCitizens.
}
\examples{
\donttest{
if(curl::has_internet()) {
route_bikecitizens()
ldf <- od_coords(stplanr::od_data_lines[2, ])
r <- route_bikecitizens(ldf)
plot(r)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catchmentArea.R
\name{calc_network_catchment}
\alias{calc_network_catchment}
\title{Calculate catchment area and associated summary statistics using network.}
\usage{
calc_network_catchment(
  sln,
  polygonlayer,
  targetlayer,
  calccols,
  maximpedance = 1000,
  distance = 100,
  projection = paste0("+proj=aea +lat_1=90 +lat_2=-18.416667",
    " +lat_0=0 +lon_0=10 +x_0=0 +y_0=0",
    " +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"),
  retainAreaProportion = FALSE,
  dissolve = FALSE
)
}
\arguments{
\item{sln}{The SpatialLinesNetwork to use.}

\item{polygonlayer}{A SpatialPolygonsDataFrame containing zones from which
the summary statistics for the catchment variable will be calculated.
Smaller polygons will increase the accuracy of the results.}

\item{targetlayer}{A SpatialPolygonsDataFrame, SpatialLinesDataFrame or
SpatialPointsDataFrame object containing the specifications of the
facilities and zones for which the catchment areas are being calculated.}

\item{calccols}{A vector of column names containing the variables in the
polygonlayer to be used in the calculation of the summary statistics for
the catchment area. If dissolve = FALSE, all other variables in the
original SpatialPolygonsDataFrame for zones that fall partly or entirely
within the catchment area will be included in the returned
SpatialPolygonsDataFrame but will not be adjusted for the proportion within
the catchment area.}

\item{maximpedance}{The maximum value of the network's weight attribute in
the units of the weight (default = 1000).}

\item{distance}{Defines the additional catchment area around the network
in the units of the projection.
(default = 100 metres)}

\item{projection}{The proj4string used to define the projection to be used
for calculating the catchment areas or a character string 'austalbers' to
use the Australian Albers Equal Area projection. Ignored if the polygonlayer
is projected in which case the targetlayer will be converted to the
projection used by the polygonlayer. In all cases the resulting object will
be reprojected to the original coordinate system and projection of the
polygon layer. Default is an Albers Equal Area projection but for more
reliable results should use a local projection (e.g., Australian Albers
Equal Area project).}

\item{retainAreaProportion}{Boolean value. If TRUE retains a variable in
the resulting SpatialPolygonsDataFrame containing the proportion of the
original area within the catchment area (Default = FALSE).}

\item{dissolve}{Boolean value. If TRUE collapses the underlying zones
within the catchment area into a single region with statistics for the
whole catchment area.}
}
\description{
Calculate catchment area and associated summary statistics using network.
}
\section{Details}{

Calculates the catchment area of a facility (e.g., cycle path) using
network distance (or other weight variable) as well as summary statistics
from variables available in a SpatialPolygonsDataFrame with census tracts
or other zones. Assumes that the frequency of the variable is evenly
distributed throughout the zone. Returns a SpatialPolygonsDataFrame.
}

\examples{
\dontrun{
data_dir <- system.file("extdata", package = "stplanr")
unzip(file.path(data_dir, "smallsa1.zip"), exdir = tempdir())
unzip(file.path(data_dir, "testcycleway.zip"), exdir = tempdir())
unzip(file.path(data_dir, "sydroads.zip"), exdir = tempdir())
sa1income <- readOGR(tempdir(), "smallsa1")
testcycleway <- readOGR(tempdir(), "testcycleway")
sydroads <- readOGR(tempdir(), "roads")
sydnetwork <- SpatialLinesNetwork(sydroads)
calc_network_catchment(
  sln = sydnetwork,
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  maximpedance = 800,
  distance = 200,
  projection = "austalbers",
  dissolve = TRUE
)
}
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope.R
\name{route_sequential_dist}
\alias{route_sequential_dist}
\title{Calculate the sequential distances between sequential coordinate pairs}
\usage{
route_sequential_dist(m, lonlat = TRUE)
}
\arguments{
\item{m}{Matrix containing coordinates and elevations}

\item{lonlat}{Are the coordinates in lon/lat order? \code{TRUE} by default}
}
\description{
Calculate the sequential distances between sequential coordinate pairs
}
\examples{
x <- c(0, 2, 3, 4, 5, 9)
y <- c(0, 0, 0, 0, 0, 1)
m <- cbind(x, y)
route_sequential_dist(m)
}
\seealso{
Other route_funs: 
\code{\link{route_average_gradient}()},
\code{\link{route_rolling_average}()},
\code{\link{route_rolling_diff}()},
\code{\link{route_rolling_gradient}()},
\code{\link{route_slope_matrix}()},
\code{\link{route_slope_vector}()}
}
\concept{route_funs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{points2flow}
\alias{points2flow}
\title{Convert a series of points into geographical flows}
\usage{
points2flow(p)
}
\arguments{
\item{p}{A spatial (point) object}
}
\description{
Takes a series of geographical points and converts them into a spatial (linestring) object
representing the potential flows, or 'spatial interaction', between every combination
of points.
}
\examples{
data(cents)
plot(cents)
flow <- points2flow(cents)
plot(flow, add = TRUE)
flow_sf <- points2flow(cents_sf)
plot(flow_sf)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catchmentArea.R
\name{calc_moving_catchment}
\alias{calc_moving_catchment}
\title{Calculate summary statistics for all features independently.}
\usage{
calc_moving_catchment(
  polygonlayer,
  targetlayer,
  calccols,
  distance = 500,
  projection = "worldalbers",
  retainAreaProportion = FALSE
)
}
\arguments{
\item{polygonlayer}{A SpatialPolygonsDataFrame containing zones from which
the summary statistics for the catchment variable will be calculated.
Smaller polygons will increase the accuracy of the results.}

\item{targetlayer}{A SpatialPolygonsDataFrame, SpatialLinesDataFrame or
SpatialPointsDataFrame object containing the specifications of the
facilities and zones for which the catchment areas are being calculated.}

\item{calccols}{A vector of column names containing the variables in the
polygonlayer to be used in the calculation of the summary statistics for
the catchment areas.}

\item{distance}{Defines the size of the catchment areas as the distance
around the targetlayer in the units of the projection
(default = 500 metres)}

\item{projection}{The proj4string used to define the projection to be used
for calculating the catchment areas or a character string 'austalbers' to
use the Australian Albers Equal Area projection. Ignored if the polygonlayer
is projected in which case the targetlayer will be converted to the
projection used by the polygonlayer. In all cases the resulting object will
be reprojected to the original coordinate system and projection of the
polygon layer. Default is an Albers Equal Area projection but for more
reliable results should use a local projection (e.g., Australian Albers
Equal Area project).}

\item{retainAreaProportion}{Boolean value. If TRUE retains a variable in
the resulting SpatialPolygonsDataFrame containing the proportion of the
original area within the catchment area (Default = FALSE).}
}
\description{
Calculate summary statistics for all features independently.
}
\section{Details}{

Calculates the summary statistics for a catchment area of multiple
facilities or zones using straight-line distance from variables
available in a SpatialPolygonsDataFrame with census tracts or other
zones. Assumes that the frequency of the variable is evenly distributed
throughout the zone. Returns the original source dataframe with additional
columns with summary variables.
}

\examples{
\dontrun{
data_dir <- system.file("extdata", package = "stplanr")
unzip(file.path(data_dir, "smallsa1.zip"))
unzip(file.path(data_dir, "testcycleway.zip"))
sa1income <- readOGR(".", "smallsa1")
testcycleway <- readOGR(".", "testcycleway")
calc_moving_catchment(
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  distance = 800,
  projection = "austalbers"
)
}
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toptail.R
\name{toptail_buff}
\alias{toptail_buff}
\title{Clip the beginning and ends SpatialLines to the edge of SpatialPolygon borders}
\usage{
toptail_buff(l, buff, ...)
}
\arguments{
\item{l}{An sf LINESTRING object}

\item{buff}{An sf POLYGON object to act as the buffer}

\item{...}{Arguments passed to rgeos::gBuffer()}
}
\description{
Takes lines and removes the start and end point, to a distance determined
by the nearest polygon border.
}
\examples{
l <- routes_fast_sf
buff <- zones_sf
r_toptail <- toptail_buff(l, buff)
nrow(l)
nrow(r_toptail)
plot(zones_sf$geometry)
plot(l$geometry, add = TRUE)
plot(r_toptail$geometry, lwd = 5, add = TRUE)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_projected.R
\name{geo_select_aeq}
\alias{geo_select_aeq}
\title{Select a custom projected CRS for the area of interest}
\usage{
geo_select_aeq(shp)
}
\arguments{
\item{shp}{A spatial object with a geographic (WGS84) coordinate system}
}
\description{
This function takes a spatial object with a geographic (WGS84)
CRS and returns a custom projected CRS focussed on the centroid of the object.
This function is especially useful for using units of metres in all directions
for data collected anywhere in the world.
}
\details{
The function is based on this stackexchange answer:
\url{https://gis.stackexchange.com/questions/121489}
}
\examples{
sp::bbox(routes_fast)
new_crs <- geo_select_aeq(routes_fast)
rf_projected <- sp::spTransform(routes_fast, new_crs)
sp::bbox(rf_projected)
line_length <- rgeos::gLength(rf_projected, byid = TRUE)
plot(line_length, rf_projected$length)
shp <- zones_sf
geo_select_aeq(shp)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{od_data_routes}
\alias{od_data_routes}
\title{Example segment-level route data}
\format{
A data frame (tibble) object
}
\description{
See \code{data-raw/generate-data.Rmd} for details on how this was created.
The dataset shows routes between origins and destinations represented in
\code{od_data_lines}
}
\examples{
od_data_routes
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneway.R
\name{od_oneway}
\alias{od_oneway}
\title{Aggregate od pairs they become non-directional}
\usage{
od_oneway(
  x,
  attrib = names(x[-c(1:2)])[vapply(x[-c(1:2)], is.numeric, TRUE)],
  id1 = names(x)[1],
  id2 = names(x)[2],
  stplanr.key = NULL
)
}
\arguments{
\item{x}{A data frame or SpatialLinesDataFrame, representing an OD matrix}

\item{attrib}{A vector of column numbers or names, representing variables to be aggregated.
By default, all numeric variables are selected.
aggregate}

\item{id1}{Optional (it is assumed to be the first column)
text string referring to the name of the variable containing
the unique id of the origin}

\item{id2}{Optional (it is assumed to be the second column)
text string referring to the name of the variable containing
the unique id of the destination}

\item{stplanr.key}{Optional key of unique OD pairs regardless of the order,
e.g., as generated by \code{\link[=od_id_max_min]{od_id_max_min()}} or \code{\link[=od_id_szudzik]{od_id_szudzik()}}}
}
\value{
\code{oneway} outputs a data frame (or \code{sf} data frame) with rows containing
results for the user-selected attribute values that have been aggregated.
}
\description{
For example, sum total travel in both directions.
}
\details{
Flow data often contains movement in two directions: from point A to point B
and then from B to A. This can be problematic for transport planning, because
the magnitude of flow along a route can be masked by flows the other direction.
If only the largest flow in either direction is captured in an analysis, for
example, the true extent of travel will by heavily under-estimated for
OD pairs which have similar amounts of travel in both directions.
Flows in both direction are often represented by overlapping lines with
identical geometries (see \code{\link[=flowlines]{flowlines()}}) which can be confusing
for users and are difficult to plot.
}
\examples{
(od_min <- od_data_sample[c(1, 2, 9), 1:6])
(od_oneway <- od_oneway(od_min))
# (od_oneway_old = onewayid(od_min, attrib = 3:6)) # old implementation
nrow(od_oneway) < nrow(od_min) # result has fewer rows
sum(od_min$all) == sum(od_oneway$all) # but the same total flow
od_oneway(od_min, attrib = "all")
attrib <- which(vapply(flow, is.numeric, TRUE))
flow_oneway <- od_oneway(flow, attrib = attrib)
colSums(flow_oneway[attrib]) == colSums(flow[attrib]) # test if the colSums are equal
# Demonstrate the results from oneway and onewaygeo are identical
flow_oneway_geo <- onewaygeo(flowlines, attrib = attrib)
flow_oneway_sf <- od_oneway(flowlines_sf)
par(mfrow = c(1, 2))
plot(flow_oneway_geo, lwd = flow_oneway_geo$All / mean(flow_oneway_geo$All))
plot(flow_oneway_sf$geometry, lwd = flow_oneway_sf$All / mean(flow_oneway_sf$All))
par(mfrow = c(1, 1))
od_max_min <- od_oneway(od_min, stplanr.key = od_id_character(od_min[[1]], od_min[[2]]))
cor(od_max_min$all, od_oneway$all)
# benchmark performance
# bench::mark(check = FALSE, iterations = 3,
#   onewayid(flowlines_sf, attrib),
#   od_oneway(flowlines_sf)
# )
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od2line}
\alias{od2line}
\alias{od2line2}
\title{Convert origin-destination data to spatial lines}
\usage{
od2line(
  flow,
  zones,
  destinations = NULL,
  zone_code = names(zones)[1],
  origin_code = names(flow)[1],
  dest_code = names(flow)[2],
  zone_code_d = NA,
  silent = FALSE
)

od2line2(flow, zones)
}
\arguments{
\item{flow}{A data frame representing origin-destination data.
The first two columns of this data frame should correspond
to the first column of the data in the zones. Thus in \code{\link[=cents]{cents()}},
the first column is geo_code. This corresponds to the first two columns
of \code{\link[=flow]{flow()}}.}

\item{zones}{A spatial object representing origins (and destinations
if no separate destinations object is provided) of travel.}

\item{destinations}{A spatial object
representing destinations of travel flows.}

\item{zone_code}{Name of the variable in \code{zones} containing the ids of the zone.
By default this is the first column names in the zones.}

\item{origin_code}{Name of the variable in \code{flow} containing the ids of the zone of origin.
By default this is the first column name in the flow input dataset.}

\item{dest_code}{Name of the variable in \code{flow} containing the ids of the zone of destination.
By default this is the second column name in the flow input dataset or the first column name in the
destinations if that is set.}

\item{zone_code_d}{Name of the variable in \code{destinations} containing the ids of the zone.
By default this is the first column names in the destinations.}

\item{silent}{TRUE by default, setting it to TRUE will show you the matching columns}
}
\description{
Origin-destination ('OD') flow data is often provided
in the form of 1 line per flow with zone codes of origin and destination
centroids. This can be tricky to plot and link-up with geographical data.
This function makes the task easier.
}
\details{
Origin-destination (OD) data is often provided
in the form of 1 line per OD pair, with zone codes of the trip origin in the first
column and the zone codes of the destination in the second column
(see the \href{https://docs.ropensci.org/stplanr/articles/stplanr-od.html}{\code{vignette("stplanr-od")}}) for details.
\code{od2line()} creates a spatial (linestring) object representing movement from the origin
to the destination for each OD pair.
It takes data frame containing
origin and destination cones (\code{flow}) that match the first column in a
a spatial (polygon or point) object (\code{zones}).
}
\examples{
od_data <- stplanr::flow[1:20, ]
l <- od2line(flow = od_data, zones = cents_sf)
plot(sf::st_geometry(cents_sf))
plot(l, lwd = l$All / mean(l$All), add = TRUE)
l <- od2line(flow = od_data, zones = cents)
# When destinations are different
head(destinations[1:5])
od_data2 <- flow_dests[1:12, 1:3]
od_data2
flowlines_dests <- od2line(od_data2, cents_sf, destinations = destinations_sf)
flowlines_dests
plot(flowlines_dests)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_projected.R
\name{geo_projected}
\alias{geo_projected}
\alias{gprojected}
\title{Perform GIS functions on a temporary, projected version of a spatial object}
\usage{
geo_projected(shp, fun, crs, silent, ...)
}
\arguments{
\item{shp}{A spatial object with a geographic (WGS84) coordinate system}

\item{fun}{A function to perform on the projected object (e.g. the the rgeos or sf packages)}

\item{crs}{An optional coordinate reference system (if not provided it is set
automatically by \code{\link[=geo_select_aeq]{geo_select_aeq()}})}

\item{silent}{A binary value for printing the CRS details (default: TRUE)}

\item{...}{Arguments to pass to \code{fun}, e.g. \code{byid = TRUE} if the function is \code{rgeos::gLength()}}
}
\description{
This function performs operations on projected data.
}
\examples{
lib_versions <- sf::sf_extSoftVersion()
lib_versions
# fails on some systems (with early versions of PROJ)
if (lib_versions[3] >= "6.3.1") {
  shp <- routes_fast_sf[2:4, ]
  geo_projected(shp, sf::st_buffer, dist = 100)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{sum_network_links}
\alias{sum_network_links}
\title{Summarise links from shortest paths data}
\usage{
sum_network_links(sln, routedata)
}
\arguments{
\item{sln}{The SpatialLinesNetwork or sfNetwork to use.}

\item{routedata}{A dataframe where the first column contains the Node ID(s)
of the start of the routes, the second column indicates the Node ID(s) of
the end of the routes, and any additional columns are summarised by link.
If there are no additional colums, then overlapping routes are counted.}
}
\description{
Summarise links from shortest paths data
}
\section{Details}{

Find the shortest path on the network between specified nodes and returns
a SpatialLinesDataFrame or sf containing the path(s) and summary statistics
of each one.
}

\examples{
sln_sf <- SpatialLinesNetwork(route_network_sf)
plot(sln_sf)
nodes_df <- data.frame(
  start = rep(c(1, 2, 3, 4, 5), each = 4),
  end = rep(c(50, 51, 52, 33), times = 5)
)
weightfield(sln_sf) # field used to determine shortest path
library(sf)
shortpath_sf <- sum_network_links(sln_sf, nodes_df)
plot(shortpath_sf["count"], lwd = shortpath_sf$count, add = TRUE)
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linefuns.R
\name{line_length}
\alias{line_length}
\title{Calculate length of lines in geographic CRS}
\usage{
line_length(l, byid = TRUE)
}
\arguments{
\item{l}{A spatial lines object}

\item{byid}{Logical determining whether the length is returned per object (default is true)}
}
\description{
Calculate length of lines in geographic CRS
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od_dist}
\alias{od_dist}
\title{Quickly calculate Euclidean distances of od pairs}
\usage{
od_dist(flow, zones)
}
\arguments{
\item{flow}{A data frame representing origin-destination data.
The first two columns of this data frame should correspond
to the first column of the data in the zones. Thus in \code{\link[=cents]{cents()}},
the first column is geo_code. This corresponds to the first two columns
of \code{\link[=flow]{flow()}}.}

\item{zones}{A spatial object representing origins (and destinations
if no separate destinations object is provided) of travel.}
}
\description{
It is common to want to know the Euclidean distance between origins and destinations
in OD data. You can calculate this by first converting OD data to SpatialLines data,
e.g. with \code{\link[=od2line]{od2line()}}. However this can be slow and overkill if you just
want to know the distance. This function is a few orders of magnitude faster.
}
\details{
Note: this function assumes that the zones or centroids in \code{cents} have a geographic
(lat/lon) CRS.
}
\examples{
data(flow)
data(cents)
od_dist(flow, cents)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linefuns.R
\name{line_midpoint}
\alias{line_midpoint}
\title{Find the mid-point of lines}
\usage{
line_midpoint(l)
}
\arguments{
\item{l}{A spatial lines object}
}
\description{
This is a wrapper around \code{\link[=SpatialLinesMidPoints]{SpatialLinesMidPoints()}} that allows it to find the midpoint
of lines that are not projected, which have a lat/long CRS.
}
\examples{
data(routes_fast)
line_midpoint(routes_fast[2:5, ])
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\docType{class}
\name{sfNetwork-class}
\alias{sfNetwork-class}
\title{An S4 class representing a (typically) transport network}
\description{
This class uses a combination of a sf layer and an igraph
object to represent transport networks that can be used for routing and
other network analyses.
}
\section{Slots}{

\describe{
\item{\code{sl}}{A sf line layer with the geometry and other attributes
for each link the in network.}

\item{\code{g}}{The graph network corresponding to \code{sl}.}

\item{\code{nb}}{A list containing vectors of the nodes connected to each node
in the network.}

\item{\code{weightfield}}{A character vector containing the variable (column) name
from the SpatialLinesDataFrame to be used for weighting the network.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od_to_odmatrix}
\alias{od_to_odmatrix}
\title{Convert origin-destination data from long to wide format}
\usage{
od_to_odmatrix(flow, attrib = 3, name_orig = 1, name_dest = 2)
}
\arguments{
\item{flow}{A data frame representing flows between origin and destinations}

\item{attrib}{A number or character string representing the column containing the attribute data
of interest from the \code{flow} data frame}

\item{name_orig}{A number or character string representing the zone of origin}

\item{name_dest}{A number or character string representing the zone of destination}
}
\description{
This function takes a data frame representing travel between origins
(with origin codes in \code{name_orig}, typically the 1st column)
and destinations
(with destination codes in \code{name_dest}, typically the second column) and returns a matrix
with cell values (from \code{attrib}, the third column by default) representing travel between
origins and destinations.
}
\examples{
od_to_odmatrix(flow)
od_to_odmatrix(flow[1:9, ])
od_to_odmatrix(flow[1:9, ], attrib = "Bicycle")
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{sum_network_routes}
\alias{sum_network_routes}
\title{Summarise shortest path between nodes on network}
\usage{
sum_network_routes(
  sln,
  start,
  end,
  sumvars = weightfield(sln),
  combinations = FALSE
)
}
\arguments{
\item{sln}{The SpatialLinesNetwork or sfNetwork to use.}

\item{start}{Integer of node indices where route starts.}

\item{end}{Integer of node indices where route ends.}

\item{sumvars}{Character vector of variables for which to calculate
summary statistics. The default value is \code{weightfield(sln)}.}

\item{combinations}{Boolean value indicating if all combinations of start
and ends should be calculated. If TRUE then every start Node ID will be routed
to every end Node ID. This is faster than passing every combination to start
and end. Default is \code{FALSE}.}
}
\description{
Summarise shortest path between nodes on network
}
\section{Details}{

Find the shortest path on the network between specified nodes and returns a
\code{SpatialLinesDataFrame} (or an \code{sf} object with LINESTRING geometry)
containing the path(s) and summary statistics of each one.

The start and end arguments must be integers representing the node index.
To find which node is closest to a geographic point, use \code{find_nearest_node()}.

If the start and end node are identical, the function will return a
degenerate line with just two (identical) points. See
\href{https://github.com/ropensci/stplanr/issues/444}{#444}.
}

\examples{
sln <- SpatialLinesNetwork(route_network)
weightfield(sln) # field used to determine shortest path
shortpath <- sum_network_routes(sln, start = 1, end = 50, sumvars = "length")
plot(shortpath, col = "red", lwd = 4)
plot(sln, add = TRUE)

# with sf objects
sln <- SpatialLinesNetwork(route_network_sf)
weightfield(sln) # field used to determine shortest path
shortpath <- sum_network_routes(sln, start = 1, end = 50, sumvars = "length")
plot(sf::st_geometry(shortpath), col = "red", lwd = 4)
plot(sln, add = TRUE)

# find shortest path between two coordinates
sf::st_bbox(sln@sl)
start_coords <- c(-1.546, 53.826)
end_coords <- c(-1.519, 53.816)
plot(sln)
plot(sf::st_point(start_coords), cex = 3, add = TRUE, col = "red")
plot(sf::st_point(end_coords), cex = 3, add = TRUE, col = "blue")
nodes <- find_network_nodes(sln, rbind(start_coords, end_coords))
shortpath <- sum_network_routes(sln, nodes[1], nodes[2])
plot(sf::st_geometry(shortpath), col = "darkred", lwd = 3, add = TRUE)

# degenerate path
sum_network_routes(sln, start = 1, end = 1)
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo-functions.R
\name{bbox_scale}
\alias{bbox_scale}
\title{Scale a bounding box}
\usage{
bbox_scale(bb, scale_factor)
}
\arguments{
\item{bb}{Bounding box object}

\item{scale_factor}{Numeric vector determining how much the bounding box will grow or shrink.
Two numbers refer to extending the bounding box in x and y dimensions, respectively.
If the value is 1, the output size will be the same as the input.}
}
\description{
Takes a bounding box as an input and outputs a bounding box of a different size, centred at the same point.
}
\examples{
bb <- matrix(c(-1.55, 53.80, -1.50, 53.83), nrow = 2)
bb1 <- bbox_scale(bb, scale_factor = 1.05)
bb2 <- bbox_scale(bb, scale_factor = c(2, 1.05))
bb3 <- bbox_scale(bb, 0.1)
plot(x = bb2[1, ], y = bb2[2, ])
points(bb1[1, ], bb1[2, ])
points(bb3[1, ], bb3[2, ])
points(bb[1, ], bb[2, ], col = "red")
}
\seealso{
Other geo: 
\code{\link{geo_bb_matrix}()},
\code{\link{geo_bb}()},
\code{\link{quadrant}()},
\code{\link{reproject}()}
}
\concept{geo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{line2routeRetry}
\alias{line2routeRetry}
\title{Convert straight spatial (linestring) object from flow data into routes retrying
on connection (or other) intermittent failures}
\usage{
line2routeRetry(lines, pattern = "^Error: ", n_retry = 3, ...)
}
\arguments{
\item{lines}{A spatial (linestring) object}

\item{pattern}{A regex that the error messages must not match to be retried, default
"^Error: " i.e. do not retry errors starting with "Error: "}

\item{n_retry}{Number of times to retry}

\item{...}{Arguments passed to the routing function, e.g. \code{\link[=route_cyclestreets]{route_cyclestreets()}}}
}
\description{
Convert straight spatial (linestring) object from flow data into routes retrying
on connection (or other) intermittent failures
}
\section{Details}{


See \code{\link[=line2route]{line2route()}} for the version that is not retried on errors.
}

\examples{
\dontrun{
data(flowlines)
rf_list <- line2routeRetry(flowlines[1:2, ], pattern = "nonexistanceerror", silent = F)
}
}
\seealso{
Other routes: 
\code{\link{line2route}()},
\code{\link{route_dodgr}()},
\code{\link{route_local}()},
\code{\link{route_osrm}()},
\code{\link{route_transportapi_public}()},
\code{\link{route}()}
}
\concept{routes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo_projected.R
\name{geo_buffer}
\alias{geo_buffer}
\title{Perform a buffer operation on a temporary projected CRS}
\usage{
geo_buffer(shp, dist = NULL, width = NULL, ...)
}
\arguments{
\item{shp}{A spatial object with a geographic CRS (e.g. WGS84)
around which a buffer should be drawn}

\item{dist}{The distance (in metres) of the buffer (when buffering simple features)}

\item{width}{The distance (in metres) of the buffer (when buffering sp objects)}

\item{...}{Arguments passed to the buffer (see \code{?rgeos::gBuffer} or \code{?sf::st_buffer} for details)}
}
\description{
This function solves the problem that buffers will not be circular when used on
non-projected data.
}
\details{
Requires recent version of PROJ (>= 6.3.0).
Buffers on \code{sf} objects with geographic (lon/lat) coordinates can also
be done with the \href{https://r-spatial.github.io/s2/}{\code{s2}} package.
}
\examples{
lib_versions <- sf::sf_extSoftVersion()
lib_versions
if (lib_versions[3] >= "6.3.1") {
  buff_sf <- geo_buffer(routes_fast_sf, dist = 50)
  plot(buff_sf$geometry)
  geo_buffer(routes_fast_sf$geometry, dist = 50)
  # on legacy sp objects (not tested)
  # buff_sp <- geo_buffer(routes_fast, width = 100)
  # class(buff_sp)
  # plot(buff_sp, col = "red")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{routes_slow}
\alias{routes_slow}
\alias{routes_slow_sf}
\title{spatial lines dataset of commuter flows on the travel network}
\format{
A spatial lines dataset 49 rows and 15 columns
}
\usage{
data(routes_slow)
}
\description{
Simulated travel route allocated to the transport network
representing the 'quietest' between \code{\link[=cents]{cents()}}
objects
with \code{\link[=od2line]{od2line()}} (see \code{\link[=flow]{flow()}}).
}
\seealso{
Other example data: 
\code{\link{destination_zones}},
\code{\link{flow_dests}},
\code{\link{flowlines}},
\code{\link{flow}},
\code{\link{route_network}},
\code{\link{routes_fast}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{summary,SpatialLinesNetwork-method}
\alias{summary,SpatialLinesNetwork-method}
\title{Print a summary of a SpatialLinesNetwork}
\usage{
\S4method{summary}{SpatialLinesNetwork}(object, ...)
}
\arguments{
\item{object}{The SpatialLinesNetwork}

\item{...}{Arguments to pass to relevant summary function.}
}
\description{
Print a summary of a SpatialLinesNetwork
}
\examples{
data(routes_fast)
rnet <- overline(routes_fast, attrib = "length")
sln <- SpatialLinesNetwork(rnet)
summary(sln)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/node-funs.R
\name{sln_add_node}
\alias{sln_add_node}
\title{Add node to spatial lines object}
\usage{
sln_add_node(sln, p)
}
\arguments{
\item{sln}{A spatial lines (\code{sfNetwork}) object created by \code{SpatialLinesNetwork}}

\item{p}{A point represented by an \code{sf} object the will split the \code{route}}
}
\description{
Add node to spatial lines object
}
\examples{
sample_routes <- routes_fast_sf[2:6, NULL]
sample_routes$value <- rep(1:3, length.out = 5)
rnet <- overline2(sample_routes, attrib = "value")
sln <- SpatialLinesNetwork(rnet)
p <- sf::st_sfc(sf::st_point(c(-1.540, 53.826)), crs = sf::st_crs(rnet))
sln_nodes <- sln2points(sln)
sln_new <- sln_add_node(sln, p)
route <- route_local(sln_new, p, sln_nodes[9, ])
plot(sln)
plot(sln_nodes, pch = as.character(1:nrow(sln_nodes)), add = TRUE)
plot(route$geometry, lwd = 9, add = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_funs.R
\name{route_rolling_gradient}
\alias{route_rolling_gradient}
\title{Calculate rolling average gradient from elevation data at segment level}
\usage{
route_rolling_gradient(elevations, distances, lag = 1, n = 2, abs = TRUE)
}
\arguments{
\item{elevations}{Elevations, e.g. those provided by the \code{cyclestreets} package}

\item{distances}{Distances, e.g. those provided by the \code{cyclestreets} package}

\item{lag}{The window size of the smoothing function. The default, 3, will take
the mean of values before, after and including each value.}

\item{n}{The window size of the smoothing function.
The default, 3, will take the mean of values before, after and including
each value.}

\item{abs}{Should the absolute (always positive) change be returned? True by default}
}
\description{
Calculate rolling average gradient from elevation data at segment level
}
\examples{
r1 <- od_data_routes[od_data_routes$route_number == 2, ]
y <- r1$elevations
distances <- r1$distances
route_rolling_gradient(y, distances)
route_rolling_gradient(y, distances, abs = FALSE)
route_rolling_gradient(y, distances, n = 3)
route_rolling_gradient(y, distances, n = 4)
r1$elevations_diff_1 <- route_rolling_diff(y, lag = 1)
r1$rolling_gradient <- route_rolling_gradient(y, distances, n = 2)
r1$rolling_gradient3 <- route_rolling_gradient(y, distances, n = 3)
r1$rolling_gradient4 <- route_rolling_gradient(y, distances, n = 4)
d <- cumsum(r1$distances) - r1$distances / 2
diff_above_mean <- r1$elevations_diff_1 + mean(y)
par(mfrow = c(2, 1))
plot(c(0, cumsum(r1$distances)), c(y, y[length(y)]), ylim = c(80, 130))
lines(c(0, cumsum(r1$distances)), c(y, y[length(y)]))
points(d, diff_above_mean)
abline(h = mean(y))
rg <- r1$rolling_gradient
rg[is.na(rg)] <- 0
plot(c(0, d), c(0, rg), ylim = c(0, 0.2))
points(c(0, d), c(0, r1$rolling_gradient3), col = "blue")
points(c(0, d), c(0, r1$rolling_gradient4), col = "grey")
par(mfrow = c(1, 1))
}
\seealso{
Other route_funs: 
\code{\link{route_average_gradient}()},
\code{\link{route_rolling_average}()},
\code{\link{route_rolling_diff}()},
\code{\link{route_sequential_dist}()},
\code{\link{route_slope_matrix}()},
\code{\link{route_slope_vector}()}
}
\concept{route_funs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od_aggregate_from}
\alias{od_aggregate_from}
\title{Summary statistics of trips originating from zones in OD data}
\usage{
od_aggregate_from(flow, attrib = NULL, FUN = sum, ..., col = 1)
}
\arguments{
\item{flow}{A data frame representing origin-destination data.
The first two columns of this data frame should correspond
to the first column of the data in the zones. Thus in \code{\link[=cents]{cents()}},
the first column is geo_code. This corresponds to the first two columns
of \code{\link[=flow]{flow()}}.}

\item{attrib}{character, column names in sl to be aggregated}

\item{FUN}{A function to summarise OD data by}

\item{...}{Additional arguments passed to \code{FUN}}

\item{col}{The column that the OD dataset is grouped by
(1 by default, the first column usually represents the origin)}
}
\description{
This function takes a data frame of OD data and
returns a data frame reporting summary statistics for each unique zone of origin.
}
\details{
It has some default settings: the default summary statistic is \code{sum()} and the
first column in the OD data is assumed to represent the zone of origin.
By default, if \code{attrib} is not set, it summarises all numeric columns.
}
\examples{
od_aggregate_from(flow)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{zones}
\alias{zones}
\alias{zones_sf}
\title{Spatial polygons of home locations for flow analysis.}
\description{
Note: we recommend using the \code{zones_sf} data.
}
\details{
These correspond to the \code{cents_sf} data.

\itemize{
\item geo_code. the official code of the zone
}
}
\examples{
library(sf)
zones_sf
plot(zones_sf)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/line_sample.R
\name{line_sample}
\alias{line_sample}
\title{Sample n points along lines with density proportional to a weight}
\usage{
line_sample(l, n, weights)
}
\arguments{
\item{l}{The SpatialLines object along which to create sample points}

\item{n}{The total number of points to sample}

\item{weights}{The relative probabilities of lines being samples}
}
\description{
Sample n points along lines with density proportional to a weight
}
\examples{
l <- flowlines[2:5, ]
n <- 100
l_lengths <- line_length(l)
weights <- l$All
p <- line_sample(l, 50, weights)
plot(p)
p <- line_sample(l, 50, weights = 1:length(l))
plot(p)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route-transport-api.R
\name{route_transportapi_public}
\alias{route_transportapi_public}
\title{Plan a single route with TransportAPI.com}
\usage{
route_transportapi_public(
  from,
  to,
  silent = FALSE,
  region = "southeast",
  modes = NA,
  not_modes = NA
)
}
\arguments{
\item{from}{Text string or coordinates (a numeric vector of
\code{length = 2} representing latitude and longitude) representing a point
on Earth.}

\item{to}{Text string or coordinates (a numeric vector of
\code{length = 2} representing latitude and longitude) representing a point
on Earth. This represents the destination of the trip.}

\item{silent}{Logical (default is FALSE). TRUE hides request sent.}

\item{region}{String for the active region to use for journey plans.
Possible values are 'southeast' (default) or 'tfl'.}

\item{modes}{Vector of character strings containing modes to use. Default is
to use all modes.}

\item{not_modes}{Vector of character strings containing modes not to use.
Not used if \code{modes} is set.}
}
\description{
Provides an R interface to the TransportAPI.com public transport API.
The function returns a SpatialLinesDataFrame object representing the
public route.
Currently only works for the United Kingdom.
See \url{https://developer.transportapi.com/documentation}for more information.
}
\details{
This function uses the online routing service
TransportAPI.com to find public routes
between origins and destinations. It does not require
any key to access the API.

Note that if \code{from} and \code{to} are supplied as
character strings (instead of lon/lat pairs), Google's
geo-coding services are used via \code{geo_code}.

Note: there is now a dedicated transportAPI package:
https://github.com/ITSLeeds/transportAPI
}
\examples{

\dontrun{
# Plan the 'public' route from Hereford to Leeds
rqh <- route_transportapi_public(from = "Hereford", to = "Leeds")
plot(rq_hfd)
}

# Aim plan public transport routes with transportAPI
}
\seealso{
line2route

Other routes: 
\code{\link{line2routeRetry}()},
\code{\link{line2route}()},
\code{\link{route_dodgr}()},
\code{\link{route_local}()},
\code{\link{route_osrm}()},
\code{\link{route}()}
}
\concept{routes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crs-funs.R
\name{reproject}
\alias{reproject}
\title{Reproject lat/long spatial object so that they are in units of 1m}
\usage{
reproject(shp, crs = geo_select_aeq(shp))
}
\arguments{
\item{shp}{A spatial object with a geographic (WGS84) coordinate system}

\item{crs}{An optional coordinate reference system (if not provided it is set
automatically by \code{\link[=geo_select_aeq]{geo_select_aeq()}}).}
}
\description{
Many GIS functions (e.g. finding the area)
}
\examples{
data(routes_fast)
rf_aeq <- reproject(routes_fast[1:3, ])
rf_osgb <- reproject(routes_fast[1:3, ], 27700)
}
\seealso{
Other geo: 
\code{\link{bbox_scale}()},
\code{\link{geo_bb_matrix}()},
\code{\link{geo_bb}()},
\code{\link{quadrant}()}
}
\concept{geo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{sln2points}
\alias{sln2points}
\title{Generate spatial points representing nodes on a SpatialLinesNetwork
or sfNetwork.}
\usage{
sln2points(sln)
}
\arguments{
\item{sln}{The SpatialLinesNetwork or sfNetwork to use.}
}
\description{
Generate spatial points representing nodes on a SpatialLinesNetwork
or sfNetwork.
}
\examples{
data(routes_fast)
rnet <- overline(routes_fast, attrib = "length")
sln <- SpatialLinesNetwork(rnet)
(sln_nodes <- sln2points(sln))
plot(sln)
plot(sln_nodes, add = TRUE)
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{plot,SpatialLinesNetwork,ANY-method}
\alias{plot,SpatialLinesNetwork,ANY-method}
\title{Plot a SpatialLinesNetwork}
\usage{
\S4method{plot}{SpatialLinesNetwork,ANY}(x, component = "sl", ...)
}
\arguments{
\item{x}{The SpatialLinesNetwork to plot}

\item{component}{The component of the network to plot. Valid values are "sl"
for the geographic (SpatialLines) representation or "graph" for the graph
representation.}

\item{...}{Arguments to pass to relevant plot function.}
}
\description{
Plot a SpatialLinesNetwork
}
\examples{
sln <- SpatialLinesNetwork(route_network)
plot(sln)
plot(sln, component = "graph")
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/node-funs.R
\name{rnet_get_nodes}
\alias{rnet_get_nodes}
\title{Extract nodes from route network}
\usage{
rnet_get_nodes(rnet, p = NULL)
}
\arguments{
\item{rnet}{A route network of the type generated by \code{overline()}}

\item{p}{A point represented by an \code{sf} object the will split the \code{route}}
}
\description{
Extract nodes from route network
}
\examples{
rnet_get_nodes(route_network_sf)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadABS.R
\name{read_table_builder}
\alias{read_table_builder}
\title{Import and format Australian Bureau of Statistics (ABS) TableBuilder files}
\usage{
read_table_builder(dataset, filetype = "csv", sheet = 1, removeTotal = TRUE)
}
\arguments{
\item{dataset}{Either a dataframe containing the original data from
TableBuilder or a character string containing the path of the
unzipped TableBuilder file.}

\item{filetype}{A character string containing the filetype. Valid values
are 'csv', 'legacycsv' and 'xlsx' (default = 'csv'). Required even when
dataset is a dataframe. Use 'legacycsv' for csv files derived from earlier
versions of TableBuilder for which csv outputs were csv versions of the
xlsx files. Current csv output from TableBuilder follow a more standard
csv format.}

\item{sheet}{An integer value containing the index of the sheet in the
xlsx file (default = 1).}

\item{removeTotal}{A boolean value. If TRUE removes the rows and columns
with totals (default = TRUE).}
}
\description{
Import and format Australian Bureau of Statistics (ABS) TableBuilder files
}
\section{Details}{

The Australian Bureau of Statistics (ABS) provides customised tables for
census and other datasets in a format that is difficult to use in R
because it contains rows with additional information.
This function imports the original (unzipped) TableBuilder files in .csv
or .xlsx format before creating an R dataframe with the data.
}

\examples{
data_dir <- system.file("extdata", package = "stplanr")
t1 <- read_table_builder(file.path(data_dir, "SA1Population.csv"))
if (requireNamespace("openxlsx")) {
  t2 <- read_table_builder(file.path(data_dir, "SA1Population.xlsx"),
    filetype = "xlsx", sheet = 1, removeTotal = TRUE
  )
}
f <- file.path(data_dir, "SA1Population.csv")
sa1pop <- read.csv(f, stringsAsFactors = TRUE, header = FALSE)
t3 <- read_table_builder(sa1pop)
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{points2odf}
\alias{points2odf}
\title{Convert a series of points into a dataframe of origins and destinations}
\usage{
points2odf(p)
}
\arguments{
\item{p}{A spatial points object}
}
\description{
Takes a series of geographical points and converts them into a data.frame
representing the potential flows, or 'spatial interaction', between every combination
of points.
}
\examples{
data(cents)
df <- points2odf(cents)
cents_centroids <- rgeos::gCentroid(cents, byid = TRUE)
df2 <- points2odf(cents_centroids)
df3 <- points2odf(cents_sf)
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{?magrittr} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{line2route}
\alias{line2route}
\title{Convert straight OD data (desire lines) into routes}
\usage{
line2route(
  l,
  route_fun = stplanr::route_cyclestreets,
  n_print = 10,
  list_output = FALSE,
  l_id = NA,
  time_delay = 0,
  ...
)
}
\arguments{
\item{l}{A spatial (linestring) object}

\item{route_fun}{A routing function to be used for converting the straight lines to routes
\code{\link[=od2line]{od2line()}}}

\item{n_print}{A number specifying how frequently progress updates
should be shown}

\item{list_output}{If FALSE (default) assumes spatial (linestring) object output. Set to TRUE to save output as a list.}

\item{l_id}{Character string naming the id field from the input lines data,
typically the origin and destination ids pasted together. If absent, the row name of the
straight lines will be used.}

\item{time_delay}{Number or seconds to wait between each query}

\item{...}{Arguments passed to the routing function, e.g. \code{\link[=route_cyclestreets]{route_cyclestreets()}}}
}
\description{
Convert straight OD data (desire lines) into routes
}
\section{Details}{


See \code{\link[=route_cyclestreets]{route_cyclestreets()}} and other route functions for details.

A parallel implementation of this was available until version 0.1.8.
}

\examples{
\dontrun{
# does not run as requires API key
l <- flowlines[2:5, ]
r <- line2route(l)
rq <- line2route(l = l, plan = "quietest", silent = TRUE)
rsc <- line2route(l = l, route_fun = cyclestreets::journey)
plot(r)
plot(r, col = "red", add = TRUE)
plot(rq, col = "green", add = TRUE)
plot(rsc)
plot(l, add = T)
# Plot for a single line to compare 'fastest' and 'quietest' route
n <- 2
plot(l[n, ])
lines(r[n, ], col = "red")
lines(rq[n, ], col = "green")
}
}
\seealso{
Other routes: 
\code{\link{line2routeRetry}()},
\code{\link{route_dodgr}()},
\code{\link{route_local}()},
\code{\link{route_osrm}()},
\code{\link{route_transportapi_public}()},
\code{\link{route}()}
}
\concept{routes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{destination_zones}
\alias{destination_zones}
\alias{destinations}
\alias{destinations_sf}
\title{Example destinations data}
\format{
A spatial dataset with 87 features
}
\usage{
data(destination_zones)
}
\description{
This dataset represents trip destinations on a different geographic
level than the origins stored in the object \code{cents}.
}
\examples{
\dontrun{
# This is how the dataset was constructed - see
# https://cowz.geodata.soton.ac.uk/download/
download.file(
  "https://cowz.geodata.soton.ac.uk/download/files/COWZ_EW_2011_BFC.zip",
  "COWZ_EW_2011_BFC.zip"
)
unzip("COWZ_EW_2011_BFC.zip")
wz <- raster::shapefile("COWZ_EW_2011_BFC.shp")
to_remove <- list.files(pattern = "COWZ", full.names = TRUE, recursive = TRUE)
file.remove(to_remove)
proj4string(wz)
wz <- sp::spTransform(wz, proj4string(zones))
destination_zones <- wz[zones, ]
plot(destination_zones)
devtools::use_data(destination_zones)
head(destination_zones@data)
destinations <- rgeos::gCentroid(destinations, byid = TRUE)
destinations <- sp::SpatialPointsDataFrame(destinations, destination_zones@data)
devtools::use_data(destinations, overwrite = TRUE)
destinations_sf <- sf::st_as_sf(destinations)
devtools::use_data(destinations_sf)
}
}
\seealso{
Other example data: 
\code{\link{flow_dests}},
\code{\link{flowlines}},
\code{\link{flow}},
\code{\link{route_network}},
\code{\link{routes_fast}},
\code{\link{routes_slow}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{weightfield}
\alias{weightfield}
\alias{weightfield<-}
\alias{weightfield,SpatialLinesNetwork-method}
\alias{weightfield,sfNetwork-method}
\alias{weightfield<-,SpatialLinesNetwork,ANY-method}
\alias{weightfield<-,sfNetwork,ANY-method}
\alias{weightfield<-,SpatialLinesNetwork,character-method}
\alias{weightfield<-,sfNetwork,character-method}
\title{Get or set weight field in SpatialLinesNetwork}
\usage{
weightfield(x)

weightfield(x, varname) <- value

weightfield(x, varname) <- value

\S4method{weightfield}{SpatialLinesNetwork}(x)

\S4method{weightfield}{sfNetwork}(x)

\S4method{weightfield}{SpatialLinesNetwork,ANY}(x) <- value

\S4method{weightfield}{sfNetwork,ANY}(x) <- value

\S4method{weightfield}{SpatialLinesNetwork,character}(x, varname) <- value

\S4method{weightfield}{sfNetwork,character}(x, varname) <- value
}
\arguments{
\item{x}{SpatialLinesNetwork to use}

\item{varname}{The name of the variable to set/use.}

\item{value}{Either the name of the variable to use as the weight field or
a dataframe or vector containing the weights to use if \code{varname} is
passed to the replacement function. If the dataframe contains multiple
columns, the column with the same name as \code{varname} is used,
otherwise the first column is used.}
}
\description{
Get or set value of weight field in SpatialLinesNetwork
}
\section{Details}{

These functions manipulate the value of weightfield in a
SpatialLinesNetwork. When changing the value of weightfield, the weights
of the graph network are updated with the values of the corresponding
variables.
}

\examples{
# with sp objects
data(routes_fast)
rnet <- overline(routes_fast, attrib = "length")
sln <- SpatialLinesNetwork(rnet)
weightfield(sln) <- "length"
weightfield(sln, "randomnum") <- sample(1:10, size = nrow(sln@sl), replace = TRUE)
data(routes_fast_sf)
rnet <- overline(routes_fast_sf, attrib = "length")
sln <- SpatialLinesNetwork(rnet)
weightfield(sln) <- "length"
sln@sl$randomnum <- sample(1:10, size = nrow(sln@sl), replace = TRUE)
weightfield(sln) <- "randomnum"
# todo: show the difference that it makes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/node-funs.R
\name{route_split_id}
\alias{route_split_id}
\title{Split route based on the id or coordinates of one of its vertices}
\usage{
route_split_id(r, id = NULL, p = NULL)
}
\arguments{
\item{r}{An \code{sf} object with one feature containing a linestring geometry to be split}

\item{id}{The index of the point on the number to be split}

\item{p}{A point represented by an \code{sf} object the will split the \code{route}}
}
\description{
Split route based on the id or coordinates of one of its vertices
}
\examples{
sample_routes <- routes_fast_sf[2:6, 3]
r <- sample_routes[2, ]
id <- round(n_vertices(r) / 2)
r_split <- route_split_id(r, id = id)
plot(r$geometry, lwd = 9, col = "grey")
plot(r_split, col = c("red", "blue"), add = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{od_data_lines}
\alias{od_data_lines}
\title{Example of desire line representations of origin-destination data from UK Census}
\format{
A data frame (tibble) object
}
\description{
Derived from \code{od_data_sample} showing movement between points represented in \code{cents_sf}
}
\examples{
od_data_lines
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_funs.R
\name{route_average_gradient}
\alias{route_average_gradient}
\title{Return average gradient across a route}
\usage{
route_average_gradient(elevations, distances)
}
\arguments{
\item{elevations}{Elevations, e.g. those provided by the \code{cyclestreets} package}

\item{distances}{Distances, e.g. those provided by the \code{cyclestreets} package}
}
\description{
This function assumes that elevations and distances are in the same units.
}
\examples{
r1 <- od_data_routes[od_data_routes$route_number == 2, ]
elevations <- r1$elevations
distances <- r1$distances
route_average_gradient(elevations, distances) # an average of a 4\% gradient
}
\seealso{
Other route_funs: 
\code{\link{route_rolling_average}()},
\code{\link{route_rolling_diff}()},
\code{\link{route_rolling_gradient}()},
\code{\link{route_sequential_dist}()},
\code{\link{route_slope_matrix}()},
\code{\link{route_slope_vector}()}
}
\concept{route_funs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rnet_roundabout}
\alias{rnet_roundabout}
\title{Example of roundabout data showing problems for SpatialLinesNetwork objects}
\format{
A sf object
}
\description{
See \code{data-raw/rnet_roundabout.R} for details on how this was created.
}
\examples{
rnet_roundabout
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ca_local}
\alias{ca_local}
\title{SpatialPointsDataFrame representing road traffic deaths}
\format{
A SpatialPointsDataFrame with 11 rows and 2 columns
}
\usage{
data(ca_local)
}
\description{
This dataset represents the type of data downloaded and cleaned
using stplanr functions. It represents a very small sample (with most variables stripped)
of open data from the UK's Stats19 dataset.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catchmentArea.R
\name{calc_catchment_sum}
\alias{calc_catchment_sum}
\title{Calculate summary statistics for catchment area.}
\usage{
calc_catchment_sum(
  polygonlayer,
  targetlayer,
  calccols,
  distance = 500,
  projection = paste0("+proj=aea +lat_1=90 +lat_2=-18.416667",
    " +lat_0=0 +lon_0=10 +x_0=0 +y_0=0",
    " +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"),
  retainAreaProportion = FALSE,
  quadsegs = NA
)
}
\arguments{
\item{polygonlayer}{A SpatialPolygonsDataFrame containing zones from which
the summary statistics for the catchment variable will be calculated.
Smaller polygons will increase the accuracy of the results.}

\item{targetlayer}{A SpatialPolygonsDataFrame, SpatialLinesDataFrame,
SpatialPointsDataFrame, SpatialPolygons, SpatialLines or SpatialPoints
object containing the specifications of the facility for which the
catchment area is being calculated. If the object contains more than one
facility (e.g., multiple cycle paths) the aggregate catchment area will be
calculated.}

\item{calccols}{A vector of column names containing the variables in the
polygonlayer to be used in the calculation of the summary statistics for
the catchment area.}

\item{distance}{Defines the size of the catchment area as the distance
around the targetlayer in the units of the projection
(default = 500 metres)}

\item{projection}{The proj4string used to define the projection to be used
for calculating the catchment areas or a character string 'austalbers' to
use the Australian Albers Equal Area projection. Ignored if the polygonlayer
is projected in which case the targetlayer will be converted to the
projection used by the polygonlayer. In all cases the resulting object will
be reprojected to the original coordinate system and projection of the
polygon layer. Default is an Albers Equal Area projection but for more
reliable results should use a local projection (e.g., Australian Albers
Equal Area project).}

\item{retainAreaProportion}{Boolean value. If TRUE retains a variable in
the resulting SpatialPolygonsDataFrame containing the proportion of the
original area within the catchment area (Default = FALSE).}

\item{quadsegs}{Number of line segments to use to approximate a quarter
circle. Parameter passed to buffer functions, default is 5 for sp and
30 for sf.}
}
\description{
Calculate summary statistics for catchment area.
}
\section{Details}{

Calculates the summary statistics for a catchment area of a facility
(e.g., cycle path) using straight-line distance from variables
available in a SpatialPolygonsDataFrame with census tracts or other
zones. Assumes that the frequency of the variable is evenly distributed
throughout the zone. Returns either a single value if calccols is of
length = 1, or a named vector otherwise.
}

\examples{
\dontrun{
data_dir <- system.file("extdata", package = "stplanr")
unzip(file.path(data_dir, "smallsa1.zip"))
unzip(file.path(data_dir, "testcycleway.zip"))
sa1income <- rgdal::readOGR(".", "smallsa1")
testcycleway <- rgdal::readOGR(".", "testcycleway")
calc_catchment_sum(
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  distance = 800,
  projection = "austalbers"
)

calc_catchment_sum(
  polygonlayer = sa1income,
  targetlayer = testcycleway,
  calccols = c("Total"),
  distance = 800,
  projection = "austalbers"
)
}
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/google-functions.R
\name{nearest_google}
\alias{nearest_google}
\title{Generate nearest point on the route network of a point using the Google Maps API}
\usage{
nearest_google(lat, lng, google_api)
}
\arguments{
\item{lat}{Numeric vector containing latitude coordinate for each coordinate
to map. Also accepts dataframe with latitude in the first column and
longitude in the second column.}

\item{lng}{Numeric vector containing longitude coordinate for each
coordinate to map.}

\item{google_api}{String value containing the Google API key to use.}
}
\description{
Generate nearest point on the route network of a point using the Google Maps API
}
\section{Details}{

Retrieve coordinates of the node(s) on the network mapped from coordinates
passed to functions.
}

\examples{
\dontrun{
nearest_google(lat = 50.333, lng = 3.222, google_api = "api_key_here")
}
}
\seealso{
Other nodes: 
\code{\link{geo_code}()}
}
\concept{nodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/node-funs.R
\name{route_split}
\alias{route_split}
\title{Split route in two at point on or near network}
\usage{
route_split(r, p)
}
\arguments{
\item{r}{An \code{sf} object with one feature containing a linestring geometry to be split}

\item{p}{A point represented by an \code{sf} object the will split the \code{route}}
}
\value{
An sf object with 2 feature
}
\description{
Split route in two at point on or near network
}
\examples{
sample_routes <- routes_fast_sf[2:6, NULL]
r <- sample_routes[2, ]
p <- sf::st_sfc(sf::st_point(c(-1.540, 53.826)), crs = sf::st_crs(r))
plot(r$geometry, lwd = 9, col = "grey")
plot(p, add = TRUE)
r_split <- route_split(r, p)
plot(r_split, col = c("red", "blue"), add = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope.R
\name{route_slope_matrix}
\alias{route_slope_matrix}
\title{Calculate the gradient of line segments from a matrix of coordinates}
\usage{
route_slope_matrix(m, e = m[, 3], lonlat = TRUE)
}
\arguments{
\item{m}{Matrix containing coordinates and elevations}

\item{e}{Elevations in same units as x (assumed to be metres)}

\item{lonlat}{Are the coordinates in lon/lat order? \code{TRUE} by default}
}
\description{
Calculate the gradient of line segments from a matrix of coordinates
}
\examples{
x <- c(0, 2, 3, 4, 5, 9)
y <- c(0, 0, 0, 0, 0, 9)
z <- c(1, 2, 2, 4, 3, 1) / 10
m <- cbind(x, y, z)
plot(x, z, ylim = c(-0.5, 0.5), type = "l")
(gx <- route_slope_vector(x, z))
(gxy <- route_slope_matrix(m, lonlat = FALSE))
abline(h = 0, lty = 2)
points(x[-length(x)], gx, col = "red")
points(x[-length(x)], gxy, col = "blue")
title("Distance (in x coordinates) elevation profile",
  sub = "Points show calculated gradients of subsequent lines"
)
}
\seealso{
Other route_funs: 
\code{\link{route_average_gradient}()},
\code{\link{route_rolling_average}()},
\code{\link{route_rolling_diff}()},
\code{\link{route_rolling_gradient}()},
\code{\link{route_sequential_dist}()},
\code{\link{route_slope_vector}()}
}
\concept{route_funs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnet_boundary_points.R
\name{rnet_boundary_points}
\alias{rnet_boundary_points}
\alias{rnet_boundary_df}
\alias{rnet_boundary_unique}
\alias{rnet_boundary_points_lwgeom}
\alias{rnet_duplicated_vertices}
\title{Get points at the beginner and end of linestrings}
\usage{
rnet_boundary_points(rnet)

rnet_boundary_df(rnet)

rnet_boundary_unique(rnet)

rnet_boundary_points_lwgeom(rnet)

rnet_duplicated_vertices(rnet, n = 2)
}
\arguments{
\item{rnet}{An sf or sfc object with LINESTRING geometry representing a route
network.}

\item{n}{The minimum number of time a vertex must be duplicated to be returned}
}
\description{
Get points at the beginner and end of linestrings
}
\examples{
has_sfheaders <- requireNamespace("sfheaders", quietly = TRUE)
if(has_sfheaders) {
rnet <- rnet_roundabout
bp1 <- rnet_boundary_points(rnet)
bp2 <- line2points(rnet) # slower version with lwgeom
bp3 <- rnet_boundary_points_lwgeom(rnet) # slower version with lwgeom
bp4 <- rnet_boundary_unique(rnet)
nrow(bp1)
nrow(bp3)
identical(sort(sf::st_coordinates(bp1)), sort(sf::st_coordinates(bp2)))
identical(sort(sf::st_coordinates(bp3)), sort(sf::st_coordinates(bp4)))
plot(rnet$geometry)
plot(bp3, add = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toptail.R
\name{toptailgs}
\alias{toptailgs}
\title{Clip the first and last n metres of SpatialLines}
\usage{
toptailgs(l, toptail_dist, tail_dist = NULL)
}
\arguments{
\item{l}{A SpatialLines object}

\item{toptail_dist}{The distance (in metres) to top the line by.
Can be either a single value or a vector of the same length as the
SpatialLines object. If tail_dist is missing, is used as the tail distance.}

\item{tail_dist}{The distance (in metres) to tail the line by. Can be
either a single value or a vector of the same length as the SpatialLines
object.}
}
\description{
Takes lines and removes the start and end point, to a distance determined
by the user. Uses the geosphere::distHaversine function and requires
coordinates in WGS84 (lng/lat).
}
\examples{
data("routes_fast")
rf <- routes_fast[2:3, ]
r_toptail <- toptailgs(rf, toptail_dist = 300)
plot(rf, lwd = 3)
plot(r_toptail, col = "red", add = TRUE)
plot(cents, add = TRUE)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/line_via.R
\name{line_via}
\alias{line_via}
\title{Add geometry columns representing a route via intermediary points}
\usage{
line_via(l, p)
}
\arguments{
\item{l}{A spatial lines object}

\item{p}{A spatial points object}
}
\description{
Takes an origin (A) and destination (B), represented by the linestring \code{l},
and generates 3 extra geometries based on points \code{p}:
}
\details{
\enumerate{
\item From A to P1 (P1 being the nearest point to A)
\item From P1 to P2 (P2 being the nearest point to B)
\item From P2 to B
}
}
\examples{
library(sf)
l <- flowlines_sf[2:4, ]
p <- destinations_sf
lv <- line_via(l, p)
lv
# library(mapview)
# mapview(lv) +
#    mapview(lv$leg_orig, col = "red")
plot(lv[3], lwd = 9, reset = FALSE)
plot(lv$leg_orig, col = "red", lwd = 5, add = TRUE)
plot(lv$leg_via, col = "black", add = TRUE)
plot(lv$leg_dest, col = "green", lwd = 5, add = TRUE)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{summary,sfNetwork-method}
\alias{summary,sfNetwork-method}
\title{Print a summary of a sfNetwork}
\usage{
\S4method{summary}{sfNetwork}(object, ...)
}
\arguments{
\item{object}{The sfNetwork}

\item{...}{Arguments to pass to relevant summary function.}
}
\description{
Print a summary of a sfNetwork
}
\examples{
data(routes_fast)
rnet <- overline(routes_fast, attrib = "length")
sln <- SpatialLinesNetwork(rnet)
summary(sln)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rnet_cycleway_intersection}
\alias{rnet_cycleway_intersection}
\title{Example of cycleway intersection data showing problems for SpatialLinesNetwork objects}
\format{
A sf object
}
\description{
See \code{data-raw/rnet_cycleway_intersection} for details on how this was created.
}
\examples{
rnet_cycleway_intersection
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overline.R
\name{islines}
\alias{islines}
\title{Do the intersections between two geometries create lines?}
\usage{
islines(g1, g2)
}
\arguments{
\item{g1}{A spatial object}

\item{g2}{A spatial object}
}
\description{
This is a function required in \code{\link[=overline]{overline()}}. It identifies
whether sets of lines overlap (beyond shared points) or
not.
}
\examples{
\dontrun{
rnet <- overline(routes_fast[c(2, 3, 22), ], attrib = "length")
plot(rnet)
lines(routes_fast[22, ], col = "red") # line without overlaps
islines(routes_fast[2, ], routes_fast[3, ])
islines(routes_fast[2, ], routes_fast[22, ])
# sf implementation
islines(routes_fast_sf[2, ], routes_fast_sf[3, ])
islines(routes_fast_sf[2, ], routes_fast_sf[22, ])
}
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyclestreets.R
\name{nearest_cyclestreets}
\alias{nearest_cyclestreets}
\title{Generate nearest point on the route network of a point using the CycleStreets.net}
\usage{
nearest_cyclestreets(shp = NULL, lat, lng, pat = api_pat("cyclestreet"))
}
\arguments{
\item{shp}{A spatial object}

\item{lat}{Numeric vector containing latitude coordinate for each coordinate
to map. Also accepts dataframe with latitude in the first column and
longitude in the second column.}

\item{lng}{Numeric vector containing longitude coordinate for each
coordinate to map.}

\item{pat}{The API key used. By default this is set to NULL and
this is usually aquired automatically through a helper, api_pat().}
}
\description{
Generate nearest point on the route network of a point using the CycleStreets.net
}
\section{Details}{

Retrieve coordinates of the node(s) on the network mapped from coordinates
passed to functions.

Note: there is now a dedicated cyclestreets package:
https://github.com/Robinlovelace/cyclestreets
}

\examples{
\dontrun{
nearest_cyclestreets(53, 0.02, pat = Sys.getenv("CYCLESTREETS"))
nearest_cyclestreets(cents[1, ], pat = Sys.getenv("CYCLESTREETS"))
nearest_cyclestreets(cents_sf[1, ], pat = Sys.getenv("CYCLESTREETS"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo-functions.R
\name{geo_bb}
\alias{geo_bb}
\alias{bb2poly}
\title{Flexible function to generate bounding boxes}
\usage{
geo_bb(
  shp,
  scale_factor = 1,
  distance = 0,
  output = c("polygon", "points", "bb")
)
}
\arguments{
\item{shp}{Spatial object (from sf or sp packages)}

\item{scale_factor}{Numeric vector determining how much the bounding box will grow or shrink.
Two numbers refer to extending the bounding box in x and y dimensions, respectively.
If the value is 1, the output size will be the same as the input.}

\item{distance}{Distance in metres to extend the bounding box by}

\item{output}{Type of object returned (polygon by default)}
}
\description{
Takes a geographic object or bounding box as an input and outputs a bounding box,
represented as a bounding box, corner points or rectangular polygon.
}
\examples{
# Simple features implementation:
shp <- routes_fast_sf
shp_bb <- geo_bb(shp, distance = 100)
plot(shp_bb, col = "red", reset = FALSE)
plot(geo_bb(routes_fast_sf, scale_factor = 0.8), col = "green", add = TRUE)
plot(geo_bb(routes_fast_sf, output = "points"), add = TRUE)
plot(routes_fast_sf$geometry, add = TRUE)
}
\seealso{
bb_scale

Other geo: 
\code{\link{bbox_scale}()},
\code{\link{geo_bb_matrix}()},
\code{\link{quadrant}()},
\code{\link{reproject}()}
}
\concept{geo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stplanr-package.R
\docType{package}
\name{stplanr-package}
\alias{stplanr-package}
\alias{stplanr}
\title{\strong{stplanr: Sustainable Transport Planning with R}}
\description{
The stplanr package provides functions to access
and analyse data for transportation research, including origin-destination analysis,
route allocation and modelling travel patterns.
}
\section{Interesting functions}{

\itemize{
\item \code{\link[=overline]{overline()}} - Aggregate overlaying route lines and data intelligently
\item \code{\link[=calc_catchment]{calc_catchment()}} - Create a 'catchment area' to show the areas serving a destination
\item \code{\link[=route_cyclestreets]{route_cyclestreets()}} - Finds the fastest routes for cyclists between two places.
}
}

\seealso{
\url{https://github.com/ropensci/stplanr}
}
\author{
Robin Lovelace \email{rob00x@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{od_coords}
\alias{od_coords}
\title{Create matrices representing origin-destination coordinates}
\usage{
od_coords(from = NULL, to = NULL, l = NULL)
}
\arguments{
\item{from}{An object representing origins
(if lines are provided as the first argument, from is assigned to \code{l})}

\item{to}{An object representing destinations}

\item{l}{Only needed if from and to are empty, in which case this
should be a spatial object representing desire lines}
}
\description{
This function takes a wide range of input data types (spatial lines, points or text strings)
and returns a matrix of coordinates representing origin (fx, fy) and destination (tx, ty) points.
}
\examples{
od_coords(from = c(0, 52), to = c(1, 53)) # lon/lat coordinates
od_coords(from = cents[1, ], to = cents[2, ]) # Spatial points
od_coords(cents_sf[1:3, ], cents_sf[2:4, ]) # sf points
# od_coords("Hereford", "Leeds") # geocode locations
od_coords(flowlines[1:3, ])
od_coords(flowlines_sf[1:3, ])
}
\seealso{
Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_dist}()},
\code{\link{od_id}},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{flow}
\alias{flow}
\title{data frame of commuter flows}
\format{
A data frame with 49 rows and 15 columns
}
\usage{
data(flow)
}
\description{
This dataset represents commuter flows (work travel) between origin
and destination zones (see \code{\link[=cents]{cents()}}).
The data is from the UK and is available as open data:
\url{https://wicid.ukdataservice.ac.uk/}.
}
\details{
The variables are as follows:

\itemize{
\item Area.of.residence. id of origin zone
\item Area.of.workplace id of destination zone
\item All. Travel to work flows by all modes
\item \verb{[,4:15]}. Flows for different modes
\item id. unique id of flow
}
Although these variable names are unique to UK data, the data
structure is generalisable and typical of flow data from any source.
The key variables are the origin and destination ids, which link to
the \code{cents} georeferenced spatial objects.
}
\examples{
\dontrun{
# This is how the dataset was constructed - see
# https://github.com/npct/pct - if download to ~/repos
flow <- readRDS("~/repos/pct/pct-data/national/flow.Rds")
data(cents)
o <- flow$Area.of.residence \%in\% cents$geo_code[-1]
d <- flow$Area.of.workplace \%in\% cents$geo_code[-1]
flow <- flow[o & d, ] # subset flows with o and d in study area
library(devtools)
flow$id <- paste(flow$Area.of.residence, flow$Area.of.workplace)
use_data(flow, overwrite = TRUE)

# Convert flows to spatial lines dataset
flowlines <- od2line(flow = flow, zones = cents)
# use_data(flowlines, overwrite = TRUE)

# Convert flows to routes
routes_fast <- line2route(l = flowlines, plan = "fastest")
routes_slow <- line2route(l = flowlines, plan = "quietest")

use_data(routes_fast)
use_data(routes_slow)
routes_fast_sf <- sf::st_as_sf(routes_fast)
routes_slow_sf <- sf::st_as_sf(routes_slow)
}

}
\seealso{
Other example data: 
\code{\link{destination_zones}},
\code{\link{flow_dests}},
\code{\link{flowlines}},
\code{\link{route_network}},
\code{\link{routes_fast}},
\code{\link{routes_slow}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope.R
\name{route_slope_vector}
\alias{route_slope_vector}
\title{Calculate the gradient of line segments from distance and elevation vectors}
\usage{
route_slope_vector(x, e)
}
\arguments{
\item{x}{Vector of locations}

\item{e}{Elevations in same units as x (assumed to be metres)}
}
\description{
Calculate the gradient of line segments from distance and elevation vectors
}
\examples{
x <- c(0, 2, 3, 4, 5, 9)
e <- c(1, 2, 2, 4, 3, 1) / 10
route_slope_vector(x, e)
}
\seealso{
Other route_funs: 
\code{\link{route_average_gradient}()},
\code{\link{route_rolling_average}()},
\code{\link{route_rolling_diff}()},
\code{\link{route_rolling_gradient}()},
\code{\link{route_sequential_dist}()},
\code{\link{route_slope_matrix}()}
}
\concept{route_funs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overline.R
\name{gsection}
\alias{gsection}
\title{Function to split overlapping SpatialLines into segments}
\usage{
gsection(sl, buff_dist = 0)
}
\arguments{
\item{sl}{SpatialLinesDataFrame with overlapping Lines to split by
number of overlapping features.}

\item{buff_dist}{A number specifying the distance in meters of the buffer to be used to crop lines before running the operation.
If the distance is zero (the default) touching but non-overlapping lines may be aggregated.}
}
\description{
Divides SpatialLinesDataFrame objects into separate Lines.
Each new Lines object is the aggregate of a single number
of aggregated lines.
}
\examples{
lib_versions <- sf::sf_extSoftVersion()
lib_versions
# fails on some systems (with early versions of PROJ)
if (lib_versions[3] >= "6.3.1") {
  sl <- routes_fast_sf[2:4, ]
  rsec <- gsection(sl)
  length(rsec) # sections
  plot(rsec, col = seq(length(rsec)))
  rsec <- gsection(sl, buff_dist = 50)
  length(rsec) # 4 features: issue
  plot(rsec, col = seq(length(rsec)))
  # dont test due to issues with sp classes on some set-ups
  # sl <- routes_fast[2:4, ]
  # rsec <- gsection(sl)
  # rsec_buff <- gsection(sl, buff_dist = 1)
  # plot(sl[1], lwd = 9, col = 1:nrow(sl))
  # plot(rsec, col = 5 + (1:length(rsec)), add = TRUE, lwd = 3)
  # plot(rsec_buff, col = 5 + (1:length(rsec_buff)), add = TRUE, lwd = 3)
}
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{plot,sfNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{update_line_geometry}
\alias{update_line_geometry}
\title{Update line geometry}
\usage{
update_line_geometry(l, nl)
}
\arguments{
\item{l}{A SpatialLines object, whose geometry is to be modified}

\item{nl}{A SpatialLines object of the same length as \code{l} to provide the new geometry}
}
\description{
Take two SpatialLines objects and update the geometry of the former with that of the latter,
retaining the data of the former.
}
\examples{
data(flowlines)
l <- flowlines[2:5, ]
nl <- routes_fast
nrow(l)
nrow(nl)
l <- l[!is_linepoint(l), ]
names(l)
names(routes_fast)
l_newgeom <- update_line_geometry(l, nl)
plot(l, lwd = l$All / mean(l$All))
plot(l_newgeom, lwd = l$All / mean(l$All))
names(l_newgeom)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_local.R
\name{route_local}
\alias{route_local}
\title{Plan a route with local data}
\usage{
route_local(sln, from, to, l = NULL, ...)
}
\arguments{
\item{sln}{The SpatialLinesNetwork or sfNetwork to use.}

\item{from}{An object representing origins
(if lines are provided as the first argument, from is assigned to \code{l})}

\item{to}{An object representing destinations}

\item{l}{Only needed if from and to are empty, in which case this
should be a spatial object representing desire lines}

\item{...}{Arguments to pass to \code{sum_network_links}}
}
\description{
This function returns the shortest path between locations
in, or near to, segements on a \code{SpatialLinesNetwork}.
}
\examples{
from <- c(-1.535181, 53.82534)
to <- c(-1.52446, 53.80949)
sln <- SpatialLinesNetwork(route_network_sf)
r <- route_local(sln, from, to)
plot(sln)
plot(r$geometry, add = TRUE, col = "red", lwd = 5)
plot(cents[c(3, 4), ], add = TRUE)
r2 <- route_local(sln = sln, cents_sf[3, ], cents_sf[4, ])
plot(r2$geometry, add = TRUE, col = "blue", lwd = 3)
l <- flowlines_sf[3:5, ]
r3 <- route_local(l = l, sln = sln)
plot(r2$geometry, add = TRUE, col = "blue", lwd = 3)
}
\seealso{
Other routes: 
\code{\link{line2routeRetry}()},
\code{\link{line2route}()},
\code{\link{route_dodgr}()},
\code{\link{route_osrm}()},
\code{\link{route_transportapi_public}()},
\code{\link{route}()}
}
\concept{routes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{line2df}
\alias{line2df}
\title{Convert geographic line objects to a data.frame with from and to coords}
\usage{
line2df(l)
}
\arguments{
\item{l}{A spatial lines object}
}
\description{
This function returns a data frame with fx and fy and tx and ty variables
representing the beginning and end points of spatial line features respectively.
}
\examples{
data(flowlines)
line2df(flowlines[5, ]) # beginning and end of a single straight line
line2df(flowlines) # on multiple lines
line2df(routes_fast[5:6, ]) # beginning and end of routes
line2df(routes_fast_sf[5:6, ]) # beginning and end of routes
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{od_data_sample}
\alias{od_data_sample}
\title{Example of origin-destination data from UK Census}
\format{
A data frame (tibble) object
}
\description{
See \code{data-raw/generate-data.Rmd} for details on how this was created.
}
\examples{
od_data_sample
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{flow_dests}
\alias{flow_dests}
\title{data frame of invented
commuter flows with destinations in a different layer than the origins}
\format{
A data frame with 49 rows and 15 columns
}
\usage{
data(flow_dests)
}
\description{
data frame of invented
commuter flows with destinations in a different layer than the origins
}
\examples{
\dontrun{
# This is how the dataset was constructed
flow_dests <- flow
flow_dests$Area.of.workplace <- sample(x = destinations$WZ11CD, size = nrow(flow))
flow_dests <- dplyr::rename(flow_dests, WZ11CD = Area.of.workplace)
devtools::use_data(flow_dests)
}

}
\seealso{
Other example data: 
\code{\link{destination_zones}},
\code{\link{flowlines}},
\code{\link{flow}},
\code{\link{route_network}},
\code{\link{routes_fast}},
\code{\link{routes_slow}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneway.R
\name{od_id}
\alias{od_id}
\alias{od_id_szudzik}
\alias{od_id_max_min}
\alias{od_id_character}
\title{Combine two ID values to create a single ID number}
\usage{
od_id_szudzik(x, y, ordermatters = FALSE)

od_id_max_min(x, y)

od_id_character(x, y)
}
\arguments{
\item{x}{a vector of numeric, character, or factor values}

\item{y}{a vector of numeric, character, or factor values}

\item{ordermatters}{logical, does the order of values matter to pairing, default = FALSE}
}
\description{
Combine two ID values to create a single ID number
}
\details{
In OD data it is common to have many 'oneway' flows from "A to B" and "B to A".
It can be useful to group these an have a single ID that represents pairs of IDs
with or without directionality, so they contain 'twoway' or bi-directional values.

\verb{od_id*} functions take two vectors of equal length and return a vector of IDs,
which are unique for each combination but the same for twoway flows.
\itemize{
\item the Szudzik pairing function, on two vectors of equal
length. It returns a vector of ID numbers.
}

This function superseeds od_id_order as it is faster on large datasets
}
\examples{
(d <- od_data_sample[2:9, 1:2])
(id <- od_id_character(d[[1]], d[[2]]))
duplicated(id)
od_id_szudzik(d[[1]], d[[2]])
od_id_max_min(d[[1]], d[[2]])
}
\seealso{
od_oneway

Other od: 
\code{\link{dist_google}()},
\code{\link{od2line}()},
\code{\link{od2odf}()},
\code{\link{od_aggregate_from}()},
\code{\link{od_aggregate_to}()},
\code{\link{od_coords2line}()},
\code{\link{od_coords}()},
\code{\link{od_dist}()},
\code{\link{od_oneway}()},
\code{\link{od_to_odmatrix}()},
\code{\link{odmatrix_to_od}()},
\code{\link{points2flow}()},
\code{\link{points2odf}()}
}
\concept{od}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toptail.R
\name{geo_toptail}
\alias{geo_toptail}
\alias{toptail}
\title{Clip the first and last n metres of SpatialLines}
\usage{
geo_toptail(l, toptail_dist, ...)
}
\arguments{
\item{l}{A SpatialLines object}

\item{toptail_dist}{The distance (in metres) to top and tail the line by.
Can either be a single value or a vector of the same length as the
SpatialLines object.}

\item{...}{Arguments passed to rgeos::gBuffer()}
}
\description{
Takes lines and removes the start and end point, to a distance determined
by the user.
}
\details{
Note: \code{\link[=toptailgs]{toptailgs()}} is around 10 times faster, but only works
on data with geographic CRS's due to its reliance on the geosphere
package.
}
\examples{
lib_versions <- sf::sf_extSoftVersion()
lib_versions
# dont test due to issues with sp classes on some set-ups
if (lib_versions[3] >= "6.3.1") {
  # l <- routes_fast[2:4, ] # to run with sp classes
  l <- routes_fast_sf[2:4, ]
  l_top_tail <- geo_toptail(l, 300)
  l_top_tail
  plot(sf::st_geometry(l_top_tail))
  plot(sf::st_geometry(geo_toptail(l, 600)), lwd = 9, add = TRUE)
}
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{flowlines}
\alias{flowlines}
\alias{flowlines_sf}
\title{spatial lines dataset of commuter flows}
\format{
A spatial lines dataset with 49 rows and 15 columns
}
\description{
Flow data after conversion to a spatial format
with \code{\link[=od2line]{od2line()}} (see \code{\link[=flow]{flow()}}).
}
\seealso{
Other example data: 
\code{\link{destination_zones}},
\code{\link{flow_dests}},
\code{\link{flow}},
\code{\link{route_network}},
\code{\link{routes_fast}},
\code{\link{routes_slow}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{sln_clean_graph}
\alias{sln_clean_graph}
\title{Clean spatial network - return an sln with a single connected graph}
\usage{
sln_clean_graph(sln)
}
\arguments{
\item{sln}{A spatial lines (\code{sfNetwork}) object created by \code{SpatialLinesNetwork}}
}
\value{
An sfNetwork object
}
\description{
See https://github.com/ropensci/stplanr/issues/344
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/line_via.R
\name{mats2line}
\alias{mats2line}
\title{Convert 2 matrices to lines}
\usage{
mats2line(mat1, mat2, crs = NA)
}
\arguments{
\item{mat1}{Matrix representing origins}

\item{mat2}{Matrix representing destinations}

\item{crs}{Number representing the coordinate system of the data, e.g. 4326}
}
\description{
Convert 2 matrices to lines
}
\examples{
m1 <- matrix(c(1, 2, 1, 2), ncol = 2)
m2 <- matrix(c(9, 9, 9, 1), ncol = 2)
l <- mats2line(m1, m2)
class(l)
l
lsf <- sf::st_sf(l, crs = 4326)
class(lsf)
plot(lsf)
# mapview::mapview(lsf)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linefuns.R
\name{line_bearing}
\alias{line_bearing}
\title{Find the bearing of straight lines}
\usage{
line_bearing(l, bidirectional = FALSE)
}
\arguments{
\item{l}{A spatial lines object}

\item{bidirectional}{Should the result be returned in a bidirectional format?
Default is FALSE. If TRUE, the same line in the oposite direction would have the same bearing}
}
\description{
This is a simple wrapper around the geosphere function \code{\link[=bearing]{bearing()}} to return the
bearing (in degrees relative to north) of lines.
}
\details{
Returns a boolean vector. TRUE means that the associated line is in fact a point
(has no distance). This can be useful for removing data that will not be plotted.
}
\examples{
lib_versions <- sf::sf_extSoftVersion()
lib_versions
# fails on some systems (with early versions of PROJ)
if (lib_versions[3] >= "6.3.1") {
  bearings_sf_1_9 <- line_bearing(flowlines_sf[1:5, ])
  bearings_sf_1_9 # lines of 0 length have NaN bearing
  bearings_sp_1_9 <- line_bearing(flowlines[1:5, ])
  bearings_sp_1_9
  plot(bearings_sf_1_9, bearings_sp_1_9)
  line_bearing(flowlines_sf[1:5, ], bidirectional = TRUE)
  line_bearing(flowlines[1:5, ], bidirectional = TRUE)
}
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialLinesNetwork.R
\name{plot,sfNetwork,ANY-method}
\alias{plot,sfNetwork,ANY-method}
\title{Plot an sfNetwork}
\usage{
\S4method{plot}{sfNetwork,ANY}(x, component = "sl", ...)
}
\arguments{
\item{x}{The sfNetwork to plot}

\item{component}{The component of the network to plot. Valid values are "sl"
for the geographic (sf) representation or "graph" for the graph
representation.}

\item{...}{Arguments to pass to relevant plot function.}
}
\description{
Plot an sfNetwork
}
\examples{
sln_sf <- SpatialLinesNetwork(route_network_sf)
plot(sln_sf)
}
\seealso{
Other rnet: 
\code{\link{SpatialLinesNetwork}},
\code{\link{calc_catchment_sum}()},
\code{\link{calc_catchment}()},
\code{\link{calc_moving_catchment}()},
\code{\link{calc_network_catchment}()},
\code{\link{find_network_nodes}()},
\code{\link{gsection}()},
\code{\link{islines}()},
\code{\link{lineLabels}()},
\code{\link{overline_spatial}()},
\code{\link{overline}()},
\code{\link{plot,SpatialLinesNetwork,ANY-method}},
\code{\link{rnet_breakup_vertices}()},
\code{\link{rnet_group}()},
\code{\link{sln2points}()},
\code{\link{sum_network_links}()},
\code{\link{sum_network_routes}()}
}
\concept{rnet}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linefuns.R
\name{n_vertices}
\alias{n_vertices}
\title{Retrieve the number of vertices from a SpatialLines or SpatialPolygons object}
\usage{
n_vertices(l)
}
\arguments{
\item{l}{A SpatialLines or SpatalPolygons object}
}
\description{
Returns a vector of the same length as the number of lines,
with the number of vertices per line or polygon.
}
\details{
See \url{https://gis.stackexchange.com/questions/58147/} for more information.
}
\examples{
n_vertices(routes_fast)
n_vertices(routes_fast_sf)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo-functions.R
\name{writeGeoJSON}
\alias{writeGeoJSON}
\title{Write to geojson easily}
\usage{
writeGeoJSON(shp, filename)
}
\arguments{
\item{shp}{Spatial data object}

\item{filename}{File name of the output geojson}
}
\description{
Provides a user-friendly wrapper for \code{sf::st_write()}. Note,
\code{geojson_write} from the geojsonio package
provides the same functionality \url{https://github.com/ropensci/geojsonio}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{osm_net_example}
\alias{osm_net_example}
\title{Example of OpenStreetMap road network}
\format{
An sf object
}
\description{
Example of OpenStreetMap road network
}
\examples{
osm_net_example
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo-functions.R
\name{geo_bb_matrix}
\alias{geo_bb_matrix}
\title{Create matrix representing the spatial bounds of an object}
\usage{
geo_bb_matrix(shp)
}
\arguments{
\item{shp}{Spatial object (from sf or sp packages)}
}
\description{
Converts a range of spatial data formats into a matrix representing the bounding box
}
\examples{
geo_bb_matrix(routes_fast)
geo_bb_matrix(routes_fast_sf)
geo_bb_matrix(cents[1, ])
geo_bb_matrix(c(-2, 54))
geo_bb_matrix(sf::st_coordinates(cents_sf))
}
\seealso{
Other geo: 
\code{\link{bbox_scale}()},
\code{\link{geo_bb}()},
\code{\link{quadrant}()},
\code{\link{reproject}()}
}
\concept{geo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/od-funs.R
\name{points2line}
\alias{points2line}
\title{Convert a series of points, or a matrix of coordinates, into a line}
\usage{
points2line(p)
}
\arguments{
\item{p}{A spatial (points) obect or matrix representing the coordinates of points.}
}
\description{
This is a simple wrapper around \code{\link[=spLines]{spLines()}} that makes the creation of
\code{SpatialLines} objects easy and intuitive
}
\examples{
p <- matrix(1:4, ncol = 2)
library(sp)
l <- points2line(p)
plot(l)
l <- points2line(cents)
plot(l)
p <- line2points(routes_fast)
l <- points2line(p)
plot(l)
l_sf <- points2line(cents_sf)
plot(l_sf)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/node-funs.R
\name{route_nearest_point}
\alias{route_nearest_point}
\title{Find nearest route to a given point}
\usage{
route_nearest_point(r, p, id_out = FALSE)
}
\arguments{
\item{r}{The input route object from which the nearest route is to be found}

\item{p}{The point whose nearest route will be found}

\item{id_out}{Should the index of the matching feature be returned? \code{FALSE} by default}
}
\description{
This function was written as a drop-in replacement for \code{sf::st_nearest_feature()},
which only works with recent versions of GEOS.
}
\examples{
r <- routes_fast_sf[2:6, NULL]
p <- sf::st_sfc(sf::st_point(c(-1.540, 53.826)), crs = sf::st_crs(r))
route_nearest_point(r, p, id_out = TRUE)
r_nearest <- route_nearest_point(r, p)
plot(r$geometry)
plot(p, add = TRUE)
plot(r_nearest, lwd = 5, add = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneway.R
\name{od_id_order}
\alias{od_id_order}
\title{Generate ordered ids of OD pairs so lowest is always first
This function is slow on large datasets, see szudzik_pairing for faster alternative}
\usage{
od_id_order(x, id1 = names(x)[1], id2 = names(x)[2])
}
\arguments{
\item{x}{A data frame or SpatialLinesDataFrame, representing an OD matrix}

\item{id1}{Optional (it is assumed to be the first column)
text string referring to the name of the variable containing
the unique id of the origin}

\item{id2}{Optional (it is assumed to be the second column)
text string referring to the name of the variable containing
the unique id of the destination}
}
\description{
Generate ordered ids of OD pairs so lowest is always first
This function is slow on large datasets, see szudzik_pairing for faster alternative
}
\examples{
x <- data.frame(id1 = c(1, 1, 2, 2, 3), id2 = c(1, 2, 3, 1, 4))
od_id_order(x) # 4th line switches id1 and id2 so stplanr.key is in order
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rnet_overpass}
\alias{rnet_overpass}
\title{Example of overpass data showing problems for SpatialLinesNetwork objects}
\format{
A sf object
}
\description{
See \code{data-raw/rnet_overpass.R} for details on how this was created.
}
\examples{
rnet_overpass
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linefuns.R
\name{angle_diff}
\alias{angle_diff}
\title{Calculate the angular difference between lines and a predefined bearing}
\usage{
angle_diff(l, angle, bidirectional = FALSE, absolute = TRUE)
}
\arguments{
\item{l}{A spatial lines object}

\item{angle}{an angle in degrees relative to North, with 90 being East and -90 being West.
(direction of rotation is ignored).}

\item{bidirectional}{Should the result be returned in a bidirectional format?
Default is FALSE. If TRUE, the same line in the oposite direction would have the same bearing}

\item{absolute}{If TRUE (the default) only positive values can be returned}
}
\description{
This function was designed to find lines that are close to parallel and perpendicular
to some pre-defined route. It can return results that are absolute (contain information
on the direction of turn, i.e. + or - values for clockwise/anticlockwise),
bidirectional (which mean values greater than +/- 90 are impossible).
}
\details{
Building on the convention used in \code{\link[=bearing]{bearing()}} and in many applications,
North is definied as 0, East as 90 and West as -90.
}
\examples{
lib_versions <- sf::sf_extSoftVersion()
lib_versions
# fails on some systems (with early versions of PROJ)
if (lib_versions[3] >= "6.3.1") {
  # Find all routes going North-South
  lines_sf <- od2line(od_data_sample, zones = zones_sf)
  angle_diff(lines_sf[2, ], angle = 0)
  angle_diff(lines_sf[2:3, ], angle = 0)
  a <- angle_diff(flowlines, angle = 0, bidirectional = TRUE, absolute = TRUE)
  plot(flowlines)
  plot(flowlines[a < 15, ], add = TRUE, lwd = 3, col = "red")
  # East-West
  plot(flowlines[a > 75, ], add = TRUE, lwd = 3, col = "green")
}
}
\seealso{
Other lines: 
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overline.R
\name{onewaygeo}
\alias{onewaygeo}
\title{Aggregate flows so they become non-directional (by geometry - the slow way)}
\usage{
onewaygeo(x, attrib)
}
\arguments{
\item{x}{A dataset containing linestring geometries}

\item{attrib}{A text string containing the name of the line's attribute to
aggregate or a numeric vector of the columns to be aggregated}
}
\value{
\code{onewaygeo} outputs a SpatialLinesDataFrame with single lines
and user-selected attribute values that have been aggregated. Only lines
with a distance (i.e. not intra-zone flows) are included
}
\description{
Flow data often contains movement in two directions: from point A to point B
and then from B to A. This can be problematic for transport planning, because
the magnitude of flow along a route can be masked by flows the other direction.
If only the largest flow in either direction is captured in an analysis, for
example, the true extent of travel will by heavily under-estimated for
OD pairs which have similar amounts of travel in both directions.
Flows in both direction are often represented by overlapping lines with
identical geometries (see \code{\link[=flowlines]{flowlines()}}) which can be confusing
for users and are difficult to plot.
}
\details{
This function aggregates directional flows into non-directional flows,
potentially halving the number of lines objects and reducing the number
of overlapping lines to zero.
}
\examples{
plot(flowlines[1:30, ], lwd = flowlines$On.foot[1:30])
singlines <- onewaygeo(flowlines[1:30, ], attrib = which(names(flowlines) == "On.foot"))
plot(singlines, lwd = singlines$On.foot / 2, col = "red", add = TRUE)
\dontrun{
plot(flowlines, lwd = flowlines$All / 10)
singlelines <- onewaygeo(flowlines, attrib = 3:14)
plot(singlelines, lwd = singlelines$All / 20, col = "red", add = TRUE)
sum(singlelines$All) == sum(flowlines$All)
nrow(singlelines)
singlelines_sf <- onewaygeo(flowlines_sf, attrib = 3:14)
sum(singlelines_sf$All) == sum(flowlines_sf$All)
summary(singlelines$All == singlelines_sf$All)
}
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/route_osrm.R
\name{route_osrm}
\alias{route_osrm}
\title{Plan routes on the transport network using the OSRM server}
\usage{
route_osrm(
  from,
  to,
  osrm.server = "https://routing.openstreetmap.de/",
  osrm.profile = "foot"
)
}
\arguments{
\item{from}{An object representing origins
(if lines are provided as the first argument, from is assigned to \code{l})}

\item{to}{An object representing destinations}

\item{osrm.server}{The base URL of the routing server.
getOption("osrm.server") by default.}

\item{osrm.profile}{The routing profile to use, e.g. "car", "bike" or "foot"
(when using the routing.openstreetmap.de test server).
getOption("osrm.profile") by default.}

\item{profile}{Which routing profile to use? One of "foot" (default)
"bike" or "car" for the default open server.}
}
\description{
This function is a simplified and (because it uses GeoJSON not binary polyline format)
slower R interface to OSRM routing services compared with the excellent
\code{\link[osrm:osrmRoute]{osrm::osrmRoute()}} function (which can be used via the \code{\link[=route]{route()}}) function.
}
\examples{
\donttest{
l1 = od_data_lines[49, ]
l1m = od_coords(l1)
from = l1m[, 1:2]
to = l1m[, 3:4]
if(curl::has_internet()) {
r_foot = route_osrm(from, to)
r_bike = route_osrm(from, to, osrm.profile = "bike")
r_car = route_osrm(from, to, osrm.profile = "car")
plot(r_foot$geometry, lwd = 9, col = "grey")
plot(r_bike, col = "blue", add = TRUE)
plot(r_car, col = "red", add = TRUE)
}
}
}
\seealso{
Other routes: 
\code{\link{line2routeRetry}()},
\code{\link{line2route}()},
\code{\link{route_dodgr}()},
\code{\link{route_local}()},
\code{\link{route_transportapi_public}()},
\code{\link{route}()}
}
\concept{routes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{stplanr-deprecated}
\alias{stplanr-deprecated}
\title{Deprecated functions in stplanr}
\description{
These functions are depreciated and will be removed:
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linefuns.R
\name{line_segment}
\alias{line_segment}
\title{Divide SpatialLines dataset into regular segments}
\usage{
line_segment(l, n_segments, segment_length = NA)
}
\arguments{
\item{l}{A spatial lines object}

\item{n_segments}{The number of segments to divide the line into}

\item{segment_length}{The approximate length of segments in the output (overides n_segments if set)}
}
\description{
Divide SpatialLines dataset into regular segments
}
\examples{
data(routes_fast)
l <- routes_fast[2, ]
library(sp)
l_seg2 <- line_segment(l = l, n_segments = 2)
plot(l_seg2, col = l_seg2$group, lwd = 50)
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{is_linepoint}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{routes_fast}
\alias{routes_fast}
\alias{routes_fast_sf}
\title{spatial lines dataset of commuter flows on the travel network}
\format{
A spatial lines dataset with 49 rows and 15 columns
}
\usage{
data(routes_fast)
}
\description{
Simulated travel route allocated to the transport network
representing the 'fastest' between \code{\link[=cents]{cents()}}
objects
with \code{\link[=od2line]{od2line()}} (see \code{\link[=flow]{flow()}}).
}
\seealso{
Other example data: 
\code{\link{destination_zones}},
\code{\link{flow_dests}},
\code{\link{flowlines}},
\code{\link{flow}},
\code{\link{route_network}},
\code{\link{routes_slow}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linefuns.R
\name{is_linepoint}
\alias{is_linepoint}
\title{Identify lines that are points}
\usage{
is_linepoint(l)
}
\arguments{
\item{l}{A spatial lines object}
}
\description{
OD matrices often contain 'intrazonal' flows, where the origin is the same point as the
destination. This function can help identify such intrazonal OD pairs, using 2 criteria:
the total number of vertices (2 or fewer) and whether the origin and destination are the same.
}
\details{
Returns a boolean vector. TRUE means that the associated line is in fact a point
(has no distance). This can be useful for removing data that will not be plotted.
}
\examples{
data(flowlines)
islp <- is_linepoint(flowlines)
nrow(flowlines)
sum(islp)
# Remove invisible 'linepoints'
nrow(flowlines[!islp, ])
}
\seealso{
Other lines: 
\code{\link{angle_diff}()},
\code{\link{geo_toptail}()},
\code{\link{line2df}()},
\code{\link{line2points}()},
\code{\link{line_bearing}()},
\code{\link{line_breakup}()},
\code{\link{line_midpoint}()},
\code{\link{line_sample}()},
\code{\link{line_segment}()},
\code{\link{line_via}()},
\code{\link{mats2line}()},
\code{\link{n_sample_length}()},
\code{\link{n_vertices}()},
\code{\link{onewaygeo}()},
\code{\link{points2line}()},
\code{\link{toptail_buff}()},
\code{\link{toptailgs}()},
\code{\link{update_line_geometry}()}
}
\concept{lines}
