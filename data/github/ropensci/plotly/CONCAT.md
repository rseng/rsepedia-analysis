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

<img src="man/figures/plotly.png" width="200" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/plotly/workflows/R-CMD-check/badge.svg)](https://github.com/plotly/plotly.R/actions)
[![CRAN
Status](https://www.r-pkg.org/badges/version/plotly)](https://cran.r-project.org/package=plotly)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/plotly)](https://www.rpackages.io/package/plotly)
[![monthly](https://cranlogs.r-pkg.org/badges/plotly)](https://www.rpackages.io/package/plotly)
<!-- badges: end -->

An R package for creating interactive web graphics via the open source
JavaScript graphing library
[plotly.js](https://github.com/plotly/plotly.js).

## Installation

Install from CRAN:

``` r
install.packages("plotly")
```

Or install the latest development version (on GitHub) via `{remotes}`:

``` r
remotes::install_github("plotly/plotly")
```

## Getting started

### Web-based ggplot2 graphics

If you use [ggplot2](https://github.com/tidyverse/ggplot2), `ggplotly()`
converts your static plots to an interactive web-based version\!

``` r
library(plotly)
g <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
  xlim(1, 6) + ylim(40, 100)
ggplotly(g)
```

![<https://i.imgur.com/G1rSArP.gifv>](https://i.imgur.com/G1rSArP.gif)

By default, `ggplotly()` tries to replicate the static ggplot2 version
exactly (before any interaction occurs), but sometimes you need greater
control over the interactive behavior. The `ggplotly()` function itself
has some convenient “high-level” arguments, such as `dynamicTicks`,
which tells plotly.js to dynamically recompute axes, when appropriate.
The `style()` function also comes in handy for *modifying* the
underlying trace
attributes (e.g. [hoveron](https://plotly.com/r/reference/#scatter-hoveron)) used to generate the plot:

``` r
gg <- ggplotly(g, dynamicTicks = "y")
style(gg, hoveron = "points", hoverinfo = "x+y+text", hoverlabel = list(bgcolor = "white"))
```

![<https://i.imgur.com/qRvLgea.gifv>](https://imgur.com/qRvLgea.gif)

Moreover, since `ggplotly()` returns a plotly object, you can apply
essentially any function from the R package on that object. Some useful
ones include `layout()` (for [customizing the
layout](https://plotly-r.com/improving-ggplotly.html#modifying-layout)),
`add_traces()` (and its higher-level `add_*()` siblings, for example
`add_polygons()`, for [adding new
traces/data](https://plotly-r.com/improving-ggplotly.html#leveraging-statistical-output)),
`subplot()` (for [combining multiple plotly
objects](https://plotly-r.com/arranging-views.html#arranging-plotly-objects)),
and `plotly_json()` (for inspecting the underlying JSON sent to
plotly.js).

The `ggplotly()` function will also respect some “unofficial”
**ggplot2** aesthetics, namely `text` (for [customizing the
tooltip](https://plotly-r.com/controlling-tooltips.html#tooltip-text-ggplotly)),
`frame` (for [creating
animations](https://plotly-r.com/animating-views.html)),
and `ids` (for ensuring sensible smooth transitions).

### Using plotly without ggplot2

The `plot_ly()` function provides a more direct interface to plotly.js
so you can leverage more specialized chart types (e.g., [parallel
coordinates](https://plotly.com/r/parallel-coordinates-plot/) or
[maps](https://plotly.com/r/maps/)) or even some visualization that the
ggplot2 API won’t ever support (e.g., surface,
[mesh](https://plotly.com/r/3d-mesh/),
[trisurf](https://plotly.com/r/trisurf/), etc).

``` r
plot_ly(z = ~volcano, type = "surface")
```

![<https://plot.ly/~brnvg/1134>](https://plot.ly/~brnvg/1134.png)

## Learn more

To learn more about special features that the plotly R package provides (e.g., [client-side linking](https://plotly-r.com/client-side-linking.html), [**shiny** integration](https://plotly-r.com/linking-views-with-shiny.html), [editing and generating static images](https://plotly-r.com/publish.html), [custom events in JavaScript](https://plotly-r.com/javascript.html), and more), see <https://plotly-r.com>. You may already be familiar with existing plotly documentation (e.g., <https://plotly.com/r/>), which is essentially a language-agnostic how-to guide for learning plotly.js, whereas <https://plotly-r.com> is meant to be more wholistic tutorial written by and for the R user. The package itself ships with a number of demos (list them by running `demo(package = "plotly")`) and shiny/rmarkdown examples (list them by running `plotly_example("shiny")` or `plotly_example("rmd")`). [Carson](https://cpsievert.me) also keeps numerous [slide decks](https://talks.cpsievert.me) with useful examples and concepts.

## Contributing

Please read through our [contributing
guidelines](https://github.com/plotly/plotly.R/blob/master/CONTRIBUTING.md).
Included are directions for opening issues, asking questions,
contributing changes to plotly, and our code of
conduct.

-----

![<https://ropensci.org>](https://www.ropensci.org/public_images/github_footer.png)
# 4.10.0.9000

## New features

* `ggplotly()` now supports the `{ggalluvial}` package. (#2061, thanks @moutikabdessabour)
* `highlight()` now supports `on="plotly_selecting"`, enabling client-side linked brushing via mouse click+drag (no mouse-up event required, as with `on="plotly_selected"`). (#1280)

## Bug fixes

* `ggplotly()` now converts `stat_ecdf()` properly. (#2065)
* `ggplotly()` now correctly handles `geom_tile()` with no `fill` aesthetic. (#2063)
* `ggplotly()` now respects `guide(aes = "none")` (e.g., `guide(fill = "none")`) when constructing legend entries. (#2067)
* Fixed an issue with translating `GGally::ggcorr()` via `ggplotly()`. (#2012)

## Improvements

* `ggplotly()` does not issue warnings with `options(warnPartialMatchArgs = TRUE)` any longer. (#2046, thanks @bersbersbers)

# 4.10.0

## Breaking changes in JavaScript API

* This version of the R package upgrades the version of the underlying plotly.js library from v1.57.1 to v2.5.1. This includes many breaking changes, bug fixes, and improvements to the underlying JavaScript library. Most of the breaking changes are summarized in [this announcement of the 2.0 release](https://community.plotly.com/t/announcing-plotly-js-2-0/53675), but [see here](https://github.com/plotly/plotly.js/blob/HEAD/CHANGELOG.md) for the full changelog.

## Breaking changes in R API

* `ggplotly()` now uses the `layout.legend.title` (instead of `layout.annotations`) plotly.js API to convert guides for discrete scales. (#1961)
* `renderPlotly()` now uses `Plotly.react()` (instead of `Plotly.newPlot()`) to redraw when `layout(transition = )` is specified. This makes it possible/easier to implement a smooth transitions when `renderPlotly()` gets re-executed. (#2001)

## New Features

* Added new functions for static image exporting via the [kaleido python package](https://github.com/plotly/Kaleido), namely `save_image()` and `kaleido()`. See `help(save_image, package = "plotly")` for installation info and example usage. (#1971)

## Improvements

* `ggplotly()` now better positions axis titles for `facet_wrap()`/`facet_grid()`. (#1975)

# 4.9.4.1

## BUG FIXES

* Fixes a bug in `ggplotly()` with `{crosstalk}` and `{ggplot2}` v3.3.4  (#1952).

# 4.9.4

## BUG FIXES

* Duplicate `highlight(selectize=T)` dropdowns are no longer rendered in Shiny (#1936).
* `group_by.plotly()` now properly retains crosstalk information across `{dplyr}` versions (#1920).
* Adds fixes in `ggplotly()` for the upcoming `{ggplot2}` v3.3.4 release (#1952).
* Fixes some issues with `name` and `frames` when both attributes are specified. (#1903 and #1618).

# 4.9.3

## Changes to plotly.js

* This version of the R package upgrades the version of the underlying plotly.js library from v1.52.2 to v1.57.1. This includes many bug fixes and improvements. The [plotly.js release page](https://github.com/plotly/plotly.js/releases) has the full list of changes.

## NEW FEATURES

* `renderPlotly()` now works well with `shiny::bindCache()`, meaning that plotly graphs can now be persistently cached in Shiny apps with `renderPlotly(expr) %>% shiny::bindCache()` (#1879).

* `ggplotly()` now works well with the [**thematic** package](https://rstudio.github.io/thematic/). That is, it can now correctly translate **ggplot2** styling that derives from **thematic**. Note that, in order to use **thematic**'s auto theming in Shiny with `ggplotly()`, you need **shiny** v1.5.0 (or higher) and **htmlwidgets** v1.5.2.9000 (or higher). Relatedly, if these versions are available, one may now also call `getCurrentOutputInfo()` inside `renderPlotly()` to get CSS styles of the output container (#1801 and #1802).

## IMPROVEMENTS

* All HTTP requests are now retried upon failure (#1656, @jameslamb).

* R linebreaks (`\n`) in _factor labels_ are now translated to HTML linebreaks (`<br />`), too. Before, this conversion was only done for colums of type character. ([#1700](https://github.com/plotly/plotly.R/pull/1700), @salim-b).

## BUG FIXES

* When R's `POSIXt` class is serialized to JSON, the time of day is now correctly preserved (in plotly.js expected `'yyyy-mm-dd HH:MM:SS.ssssss'` format). This should fix a whole host of issues where date-times were being rounded. (#1871, @FlukeAndFeather).

* `ggplotly()` now handles discrete axes of a `facet_wrap` and `facet_grid` correctly when there is only one category in panels > 1 (#1577 and #1720).

* `ggplotly()` now correctly accounts for linebreaks in tick label text when computing plot margins (#1791, @trekonom).

* `ggplotly()` now handles `element_blank()` and `factor()` labels in positional scales correctly (#1731 and #1772).

* `ggplotly()` now handles missing `y` aesthetic in `geom_errorbar()` (#1779, @trekonom).

# 4.9.2.1

This is minor patch release with a few minor bug fixes and updates test expectations in anticipation of new R 4.0 defaults.

## BUG FIXES

* Fixes rendering issues non-HTML **rmarkdown** output formats, which was introduced in the 4.9.2 release (#1702).
* Fixes a `ggplotly()` bug in axis tick translation (#1725, #1721).
* `plot_mapbox()` now correctly defaults to a scattermapbox trace (unless z is present, then it defaults to choroplethmapbox) (#1707).
* `ggplotly()` now correctly resolves overlapping axis tick text in `coord_sf()` (#1673).
* A false-positive warning is no longer thrown when attempting to cast `toWebGL()` (#1569).

# 4.9.2

## Changes to plotly.js

* This version of the R package upgrades the version of the underlying plotly.js library from v1.49.4 to v1.52.2. This includes many bug fixes, improvements, as well as 2 new trace types: `treemap` and `image`. The [plotly.js release page](https://github.com/plotly/plotly.js/releases) has the full list of changes.

## IMPROVEMENTS

* The `add_image()` function was added to make it easier to create image traces via `raster` objects.

## BUG FIXES

* `add_sf()`/`geom_sf()` now correctly handle geometry columns that are named something other than `"geometry"` (#1659).
* Specifying an english locale no longer results in error (#1686).

# 4.9.1

## Changes to plotly.js

* This version of the R package upgrades the version of the underlying plotly.js library from v1.46.1 to v1.49.4. The [plotly.js release page](https://github.com/plotly/plotly.js/releases) has the full list of changes.

## IMPROVEMENTS

* `event_data()` gains support for the `plotly_sunburstclick` event (#1648)

## BUG FIXES

* Fixed an issue with correctly capturing the return value of user-expressions to `renderPlotly()` (#1528).
* Fixed a resizing issue where graphs could be incorrectly resized to their initial size in some cases (#1553). 
* `ggplotly()` now positions the x-axis in the last column of a `facet_wrap()` properly (#1501).
* `ggplotly()` now handles `geom_hline()`/`geom_vline()` correctly in conjunction with `coord_flip()` (#1519).
* `event_data()` now correctly relays the `key` attribute for statistical traces (#1610).

# 4.9.0

## Changes to plotly.js

* This version of the R package upgrades the version of the underlying plotly.js library from v1.42.3 to v1.46.1. The [plotly.js release page](https://github.com/plotly/plotly.js/releases) has the full list of changes, but here is summary most pertainent ones for the R package:
  * New trace types: `sunburst`, `waterfall`, `isosurface`.
  * New `hovertemplate` attribute allows for finer-tuned control over tooltip text. See [here](https://plotly-r.com/controlling-tooltips.html#fig:tooltip-format-heatmap) for an example.
  * Providing a string to `title` is now deprecated (but still works). Instead, use `title = list(text = "title")`. This change was made to support a new title placement API (e.g., `title = list(text = "title", xanchor = "left")`). Note that these changes are relevant for `layout.title` as well as `layout.xaxis.title`/`layout.yaxis.title`/etc.

## NEW FEATURES & IMPROVEMENTS

* Several new features and improvements related to accessing plotly.js events in shiny (learn more about them in this RStudio [webinar](https://www.rstudio.com/resources/webinars/accessing-and-responding-to-plotly-events-in-shiny/)):
    * The `event` argument of the `event_data()` function now supports the following events: `plotly_selecting`, `plotly_brushed`, `plotly_brushing`, `plotly_restyle`, `plotly_legendclick`, `plotly_legenddoubleclick`, `plotly_clickannotation`, `plotly_afterplot`, `plotly_doubleclick`, `plotly_deselect`, `plotly_unhover`. For examples, see `plotly_example("shiny", "event_data")`, `plotly_example("shiny", "event_data_legends")`, and  `plotly_example("shiny", "event_data_annotation")`, 
    * New `event_register()` and `event_unregister()` functions for declaring which events to transmit over the wire (i.e., from the browser to the shiny server). Events that are likely to have large overhead are not registered by default, so you'll need to register these: `plotly_selecting`, `plotly_unhover`, `plotly_restyle`, `plotly_legendclick`, and `plotly_legenddoubleclick`.
    * A new `priority` argument. By setting `priority='event'`, the `event` is treated like a true event: any reactive expression using the `event` becomes invalidated (regardless of whether the input values has changed). For an example, see `plotly_example("shiny", "event_priority")`.
    * The `event_data()` function now relays the (official plotly.js) `customdata` attribute in similar fashion to (unofficial) `key` attribute (#1423). Run `plotly_example("shiny", "event_data")` for an example. 
    * `event_data("plotly_selected")` is no longer too eager to clear. That is, it is no longer set to `NULL` when clicking on a plot *after* triggering the "plotly_selected" event (#1121) (#1122).
    
* Several new features and improvements for exporting static graphs with the orca command-line utility:
    * The `orca()` function now supports conversion of much larger figures (#1322) and works without a mapbox api token (#1314).
    * The `orca_serve()` function was added for efficient exporting of many plotly graphs. For examples, see `help(orca_serve)`.
    * The `orca()` function gains new arguments `more_args` and `...` for finer control over the underlying system commands.

* `ggplotly()` now respects horizontal alignment of **ggplot2** titles (e.g., `ggplotly(qplot(1:10) + ggtitle("A title") + theme(plot.title = element_text(hjust = 1)))`).
* **plotly** objects can now be serialized and unserialized in different environments (i.e., you can now use `saveRDS()` to save an object as an rds file and restore it on another machine with `readRDS()`). Note this object is *dynamically* linked to JavaScript libraries, so one should take care to use consistent versions of **plotly** when serializing and unserializing (#1376).
* The `style()` function now supports "partial updates" (i.e. modification of a particular property of an object, rather than the entire object). For example, notice how the first plot retains the original marker shape (a square): `p <- plot_ly(x = 1:10, y = 1:10, symbol = I(15)); subplot(style(p, marker.color = "red"), style(p, marker = list(color = "red")))` (#1342).
* The `method` argument of `plotlyProxyInvoke()` gains support for a `"reconfig"` method. This makes it possible to modify just the configuration of a plot in a **shiny** app. For an example use, see `plotly_example("shiny", "event_data_annotation")`.
* The `plotly_example()` function will now attempt to open the source file(s) used to run the example. Set `edit = FALSE` to prevent the source file(s) from opening.
* An informative warning is now thrown if invalid argument names are supplied to `config()`.

## CHANGES

* If `stroke` is specified, `span` now defaults to `I(1)`. This results in a slightly narrower default span for some trace types (e.g., `box`, `contour`), but it also ensures the `stroke` is always visible when it's relevant (e.g. `plot_ly(x = 1:10, y = 1:10, stroke = I("black"))`), making for a more consistent overall default (#1507).
* The 'collaborate' button no longer appears in the modebar, and is longer supported, so the `config()` function no longer has a `collaborate` argument.
* The `cloud` argument is now deprecated and will be removed in a future version. Use `showSendToCloud` instead.
* `ggplotly()` now translates title information to `layout.title.text` (instead of `layout.title`) and `layout.title.font` (instead of `layout.titlefont`)

## BUG FIXES

* `subplot()` now works much better with annotations, images, and shapes:
  - When `xref`/`yref` references an x/y axis these references are bumped accordingly (#1181).
  - When `xref`/`yref` references paper coordinates, these coordinates are updated accordingly (#1332).
* `subplot()` now repositions shapes with fixed height/width (i.e., `xsizemode`/`ysizemode` of `"pixel"`) correctly (#1494).
* The `colorscale` attribute now correctly handles a wider range of input values (#1432, #1485)
* The colorscale generated via the `color` argument in `plot_ly()` now uses an evenly spaced grid of values instead of quantiles (#1308).
* When using **shinytest** to test a **shiny** that contains **plotly** graph, false positive differences are no longer reported (rstudio/shinytest#174). 
* When the `size` argument maps to `marker.size`, it now converts to an array of appropriate length (#1479).
* The `color` and `stroke` arguments now work as expected for trace types with `fillcolor` but no `fill` attribute (e.g. `box` traces) (#1292).
* Information emitted by in `event_data()` for heatmaps with atomic vectors for `x`/`y`/`z` is now correct (#1141).
* Fixed issue where **dplyr** groups caused a problem in the ordering of data arrays passed to `marker` objects (#1351).
* In some cases, a `ggplotly()` colorbar would cause issues with hover behavior, which is now fixed (#1381).
* An articial marker no longer appears when clearing a crosstalk selection of a plot with a colorbar (#1406).
* Clearing a highlight event via crosstalk no longer deletes all the traces added since initial draw (#1436).
* Recursive attribute validation is now only performed on recursive objects (#1315).
* The `text` attribute is no longer collapsed to a string when `hoveron='fills+points'` (#1448). 
* `layout.[x-y]axis.domain` is no longer supplied a default when `layout.grid` is specified (#1427).
* When uploading charts to a plot.ly account via `api_create()`, layout attributes are no longer incorrectly src-ified, which was causing inconsistencies in local/remote rendering of `ggplotly()` charts (#1197).

# 4.8.0

## NEW FEATURES & IMPROVEMENTS

### plotly.js and `plot_ly()` specific improvements

* Upgraded to plotly.js v1.39.2. A _huge_ amount of features and improvements have been made since v1.29.2 (i.e., the version included in the last CRAN release of the R package - v4.7.1). Highlights include a complete re-write of `scattergl` to make it nearly feature complete with `scatter`, localization of text rendering (i.e., international translations), and six new trace types (`cone`, `scatterpolar`, `scatterpolargl`, `splom`, `table`, & `violin`)! See [here](https://github.com/plotly/plotly.js/releases) for a complete list of plotly.js-specific improvements.
* Support for **sf** (simple feature) data structures was added to `plot_ly()`, `plot_mapbox()`, and `plot_geo()` (via the new `add_sf()` function). See [this blog post](https://blog.cpsievert.me/2018/03/30/visualizing-geo-spatial-data-with-sf-and-plotly) for an overview.
* Better control over the stroke (i.e., outline) appearance of various filled graphical marks via the new "special arguments" (`stroke`, `strokes`, `alpha_stroke`, `span`, and `spans`). For an overview, see the **sf** blog post linked to in the bullet point above and the new package demos (list all demos with `demo(package = "plotly")`).

### `ggplotly()` specific improvements

* `ggplotly()` now supports conversion of **ggplot2**'s `geom_sf()`.
* One may now inform `ggplotly()` about the relevant **shiny** output size via `session$clientData`. This ensures `ggplotly()` sizing is closer to **ggplot2** sizing, even on window resize. For an example, run `plotly_example("shiny", "ggplotly_sizing")`.

### Other improvements relevant for all **plotly** objects

* LaTeX rendering via MathJax is now supported and the new `TeX()` function may be used to flag a character vector as LaTeX (#375). Use the new `mathjax` argument in `config()` to specify either external (`mathjax="cdn"`) or local (`mathjax="local"`) MathJaX. If `"cdn"`, mathjax is loaded externally (meaning an internet connection is needed for TeX rendering). If `"local"`, the PLOTLY_MATHJAX_PATH environment variable must be set to the location (a local file path) of MathJax. IMPORTANT: **plotly** uses SVG-based mathjax rendering which doesn't play nicely with HTML-based rendering (e.g., **rmarkdown** documents and **shiny** apps). To leverage both types of rendering, you must `<iframe>` your plotly graph(s) into the larger document (see [here](https://github.com/plotly/plotly.R/blob/master/inst/examples/rmd/MathJax/index.Rmd) for an **rmarkdown** example  and [here](https://github.com/plotly/plotly.R/blob/master/inst/examples/rmd/MathJax/index.Rmd) for a **shiny** example).
* The selection (i.e., linked-brushing) mode can now switch from 'transient' to 'persistent' by holding the 'shift' key. It's still possible to _force_ persistent selection by setting `persistent = TRUE` in `highlight()`, but `persistent = FALSE` (the default) is now recommended since it allows one to switch between [persistent/transient selection](https://plotly-r.com/client-side-linking.html#fig:txmissing-modes) in the browser, rather than at the command line.
* The `highlight()` function gains a `debounce` argument for throttling the rate at which `on` events may be fired. This is mainly useful for improving user experience when `highlight(on = "plotly_hover")` and mousing over relevant markers at a rapid rate (#1277)
* The new `partial_bundle()` function makes it easy to leverage [partial bundles of plotly.js](https://github.com/plotly/plotly.js#partial-bundles) for reduced file sizes and faster render times.
* The `config()` function gains a `locale` argument for easily changing localization defaults (#1270). This makes it possible localize date axes, and in some cases, modebar buttons (#1270).
* The `plot_geo()` function gains a `offline` argument for rendering `"scattergeo"` traces with or without an internet connection (#356). Leveraging this argument requires the new **plotlyGeoAssets** package.
* Support for async rendering of inside **shiny** apps using the [promises](https://rstudio.github.io/promises/) package (#1209). For an example, run `plotly_example("shiny", "async")`.
* Instead of an error, `ggplotly(NULL, "message")` and `plotly_build(NULL, "message")` now returns `htmltools::div("message")`, making it easier to relay messages in shiny when data isn't yet ready to plot (#1116).
* The `animation_button()` function gains a `label` argument, making it easier to control the label of an animation button generated through the `frame` API (#1205).
* The new `highlight_key()` function provides a wrapper around `crosstalk::SharedData$new()`, making it easier to teach others how to leverage `SharedData` objects with **plotly** and **crosstalk**.

## CHANGES

### `plot_ly()` specific changes

* The `name` attribute is now a "special `plot_ly()` argument" and behaves similar to `split` (it ensures a different trace for every unique value supplied). Although this leads to a breaking change (`name` was previously appended to an automatically generated trace name), it leads to a more flexible and transparent API. Those that wish to have the old behavior back should provide relevant mappings to the `name` attributes (e.g. `plot_ly(mtcars, x = ~wt, y = ~mpg, color = ~factor(vs), name = "a")` should become `plot_ly(mtcars, x = ~wt, y = ~mpg, color = ~factor(vs), name = ~paste(vs, "\na"))`)
* The `color` argument now maps to `fillcolor`, making it much easier to use polygon fills to encode data values (e.g., choropleth maps). For backwards-compatibilty reasons, when `color` maps to `fillcolor`, `alpha` defaults to 0.5 (instead of 1). For an example, `plot_mapbox(mn_res, color = ~INDRESNAME)` or `plot_mapbox(mn_res, split = ~INDRESNAME, color = ~AREA, showlegend = FALSE, stroke = I("black"))`.
* The `color` argument no longer automatically add `"markers"` to the `mode` attribute for scatter/scattergl trace types. Those who wish to have the old behavior back, should add `"markers"` to the `mode` explicity (e.g., change `plot_ly(economics, x = ~pce, y = ~pop, color = ~as.numeric(date), mode = "lines")` to `plot_ly(economics, x = ~pce, y = ~pop, color = ~as.numeric(date), mode = "lines+markers")`).
* The `size` argument now informs a default [error_[x/y].width](https://plotly.com/r/reference/#scatter-error_x-width) (and `span` informs [error_[x/y].thickness](https://plotly.com/r/reference/#scatter-error_x-thickness)). Note you can override the default by specifying directly (e.g. `plot_ly(x = 1:10, y = 1:10, size = I(10), error_x = list(value = 5, width = 0))`).
* `layout.showlegend` now defaults to `TRUE` for a *single* pie trace. This is a more sensible default and matches pure plotly.js behavior.

### Other changes relevant for all **plotly** objects

* All axis objects now default to `automargin = TRUE`. The majority of the time this should make axis labels more readable, but may have un-intended consequences in some rare cases (#1252). 
* The `elementId` field is no longer populated, which fixes the "Ignoring explicitly provided widget ID" warning in shiny applications (#985).

## BUG FIXES

### `ggplotly()` specific fixes

* The default `height`/`width` that `ggplotly()` assumes is now more consistently correct in various context, but it also now requires access to one of the following devices: `Cairo::Cairo()`, `png()`, or `jpg()`. 
* In RStudio, `ggplotly()` was ignoring a specified `height`/`width` (#1190).
* `ggplotly()` now uses fixed heights for facet strips meaning that their height is still correct after a window resize (#1265).

### `plot_ly()` specific fixes

* The `limits` argument of `colorbar()` wasn't being applied to `line.color`/`line.cmin`/`line.cmax` (#1236).
* The `legendgroup` can now properly map data values (#1148).

### Other fixes relevant for all **plotly** objects

* Marker sizes (i.e., `marker.size`) are now _always_ based on the area when `marker.sizemode='area'` (which is the default sizemode when using the `size` argument). Previously, traces with one just one value supplied to `marker.size` were being sized by their diameter (#1133).
* Bug fix for linking views with crosstalk where the source of the selection is an aggregated trace (#1218).
* Resizing behavior, after updating `height`/`width` via **shiny** reactive values, is now correct (#1068).
* Fixed algorithm for coercing the proposed layout to the plot schema (#1156).
* `add_*()` no longer inherits `crosstalk::SharedData` key information when `inherit = FALSE` (#1242).


# 4.7.1

## NEW FEATURES & IMPROVEMENTS

* It is now possible to modify (i.e., update without a full redraw) plotly graphs inside of a shiny app via the new `plotlyProxy()` and `plotlyProxyInvoke()` functions. For examples, see `plotly_example("shiny", "proxy_relayout")` and `plotly_example("shiny", "proxy_mapbox")`. Closes #580. 
* Added a new `plotly_example()` function to make it easier to run shiny/rmarkdown examples included with the package under the `inst/examples` directory.
* The `schema()` function now returns the plot schema (rather just printing it), making it easier to acquire/use values from the official plot schema. See `help(schema)` for an example. Fixes #1038.

## CHANGES

* Upgraded to plotly.js v1.29.2 -- https://github.com/plotly/plotly.js/releases/tag/v1.29.2

## BUG FIXES

* The default sizing in `ggplotly()` is no longer fixed to the device size inside RStudio. Fixes #1033.
* Removed use of `ArrayBuffer.isView()`, which should fix rendering issues on plaforms that don't have a typed array polyfill (e.g., RStudio on Windows). Fixes #1055.
* `event_data("plotly_relayout")` no longer fires `NULL` for any event. Fixes #1039.
* Fixed a bug when using `color` with scattermapbox/scattergeo. Fixes #1038.
* Fixed a highlighting bug when brushing multiple points with `marker.color` as an array. Fixes #1084.


# 4.7.0

## NEW FEATURES & IMPROVEMENTS

* Added support for fixed coordinates (i.e., the aspect ratio component of `coord_equal()`, `coord_fixed()`, `coord_map()`, `coord_quickmap()`).
* Added support for `geom_sf()` and `coord_sf()`.
* The (previously internal) `group2NA()` function is now exported and its performance has been greatly improved thanks to the new **data.table** dependency. Essentially any geom conversion that reduces to a polygon/path should see speed improvements. Similarly, any `plot_ly()` graph that use `group_by()` in conjunction with `add_lines()`, `add_paths()`, `add_segments()`, etc will also see improvements, especially when there is a large number of groups. For details on the speed improvements, see #1022 and #996 (thanks @msummersgill).
* The `api_create()` function gains a new `fileopt` argument, which is inspired from the `fileopt` argument in the (deprecated) `plotly_POST()` function (fixes #976). It currently supports to values: `"new"` and `"overwrite"`. The default, `"overwrite"`, will overwrite existing file(s) with a matching `filename`.
* The `filename` argument in `api_create()` now accepts a character vector of length 2, the first string is used to name the plot, and the second is used to name the grid (i.e., data behind the plot).

## CHANGES

* Upgraded to plotly.js v1.27.1 -- https://github.com/plotly/plotly.js/releases/tag/v1.27.1
* The `traces` argument in the `style()` function now defaults to `NULL` (instead of 1). Meaning that, by default, supplied attributes now modify _every_ trace (instead of the first one).

## Bug fixes 

* Fixes numerous problems with `coord_flip()` (fixes #1012).
* The typed array polyfill is now included *before* the plotly.js bundle, which should fix some rendering issues in some browsers, including RStudio (fixes #1010).
* When creating private plots (via `api_create()`), both the plot and the data behind the plot are private (fixes #976).
* Creating a plot with multiple traces (or frames) via (via `api_create()`) no longer creates multiple grids (fixes #1004).
* The `api_create()` function should now create grid references for all data array attributes (fixes #1014).
* `ggplotly()` no longer opens an (off-screen) graphics device in RStudio for sizing. It will now correctly use the size of the viewer panel when querying the size of the graphics device.
* Margins are no longer always set to `NULL` for pie charts (fixes #1002)
* Fixed bug when highlight multiple 'simple key' traces (fixes #974).

# 4.6.0

## NEW FEATURES & IMPROVEMENTS

* Added a significant amount of support for "multiple linked views". For some relatively basic examples, see the demos (the ones prefixed with "highlight" are most relevant) -- `demo(package = "plotly")`. For a more comprehensive overview, see <https://plotly-r.com/client-side-linking.html>. For some more complex examples, see <https://pedestrians.cpsievert.me/>
* Added the `highlight()` function for configuring selection modes/sequences/options.
* Added support for animation. For some relatively basic examples, see the examples section of `help(animation)`. For a more thorough overview, see <https://plotly-r.com/animating-views.html>
* Added a `frame` argument to `plot_ly()` for creating animations. Also added the `animation_opts()`, `animation_slider()`, and `animation_button()` functions for configuring animation defaults.
* Added a new interface to [v2 of the REST API](https://api.plot.ly/v2). This new interface makes the  `plotly_POST()` and `get_figure()` functions obsolete (use `api_create()` and `api_download_plot()` instead), and thus, are now deprecated, but remain around for backwards-compatibility. For more details, see `help(api)`.
* Added support for conversion of more **ggplot2** geoms via `ggplotly()`: `GeomCol`, `GeomRug`, `GeomCrossbar`, `GeomQuantile`, `GeomSpoke`, `GeomDotplot`, `GeomRasterAnn` (i.e., `annotation_raster()`), and `GeomAnnotationMap` (i.e., `annotation_map()`).
* Added a new function `raster2uri()` which makes it easier to embed raster objects as [images](https://plotly.com/r/reference/#layout-images) via data URIs. For examples, see `help(raster2uri)`.
* `ggplotly()` gains a new argument, `dynamicTicks`, which allows axis ticks to update upon zoom/pan interactions (fixes #485).
* Sensible sizing and positioning defaults are now provided for subplots multiple colorbars.
* R linebreaks are translated to HTML linebreaks (i.e., '\n' translates to '<br />') (fixes #851).
* Added a `plot_dendro()` function for a quick and dirty interactive dendrogram with support for hierarchial selection. For more, see -- <https://plotly-r.com/client-side-linking.html#fig:dendro>
* The `export()` function gains a `selenium` argument for rendering/exporting WebGL plots and exporting to 'svg'/'webp' file formats (via the plotly.js function [Plotly.downloadImage()](https://plotly.com/javascript/plotlyjs-function-reference/#plotlydownloadimage)).
* Better type checking of trace attributes will now automatically reduce a single-valued vector to a constant (when appropriate). This is particularly useful for anchoring multiple traces to a single legend entry via `legendgroup` (see #675, #817, #826).
* The `plotlyOutput()` function gains a `inline` argument which makes it easier to place multiple plots in the same row (in a shiny application).

## CHANGES

* Upgraded to plotly.js v1.26.1 -- https://github.com/plotly/plotly.js/releases/tag/v1.26.1
* `ggplotly()` now applies `format()` to automatically generated hoverinfo. This will allow for finer control over the text displayed (for instance, `options(digits = 4)` can now be used to choose the number of significant digits used). See #834 for an example.
* `HTMLwidgets.renderValue()` should now avoid creating too many active WebGL contexts (thanks @AleksandrIanevski).
* A TypedArray polyfill is now included by default, and the function `remove_typedarray_polyfill()` was added to make it easy to remove it. Fixes #824, #717, #825.
* If no graphics device is already open, `ggplotly()` now tries to open/close a Cairo graphics device, then a bitmap (png/jpeg) device. If neither is available, it errors. This helps to ensure that a *screen* device is never opened by `ggplotly()` (which fixes #829). Furthermore, if `width`/`height` is not specified *and* no graphics device is currently open, a default of 640/480 is used for width/height of the device. 

## BUG FIXES


* Placement of bars (in all cases, even when representing a negative count) should now be correct (applies to `geom_bar()`, `geom_histogram()`, `geom_col()`). Fixes #560, #874, #901, #831.
* Fix for hoverinfo displaying the heights of bars in the translation `geom_bar()` via `ggplotly()`. Fixes #557 and #662.
* `embed_notebook()` now works in *nteract* notebooks (see #768). 
* Axis categories are no longer reordered for matrices (see #863).
* Fix for hoverinfo displaying values after scale transformations (in `ggplotly()`). Fixes #804.
* Font faces for axis titles are now translated (in `ggplotly()`). Fixes #861.

# 4.5.6

## NEW FEATURES

* Added support for the `timezone` argument in __ggplot2__'s `scale_datetime()`. Fixes (#743, thanks @earowang).

## CHANGES

* Now requires  a version of __ggplot2__ higher than 2.1.0 because the new ggproto faceting infrastructure introduced breaking changes.
* A book icon is added to the mode bar, by default, which links to the plotly book. If you want to remove this icon from a plot `p`, do `config(p, modeBarButtonsToRemove = "Collaborate")`
* Specifying height/width in `layout()` is now deprecated. Specify in `ggplotly()` or `plot_ly()`.
* The `ggplotly()` function now preserves all information about the layer mapping. This makes it possible to access input/output data from any layer.

## BUG FIXES

* HTMLwidget.resize() used to ignore specified `layout.width`/`layout.height`.
* When `height`/`width` are specified in `ggplotly()`, relative sizes are now translated correctly. Fixes #489 and #510.
* More careful handling of font when expanding annotation arrays. Fixes #738.
* Ignore data arrays of non-tidy traces. Fixes #737.
* When using `ggplotly()` on a plot with `geom_line` and `group` aesthetic wrong tooltip information was shown. Fixes #774.

# 4.5.5 -- 28 September 2016

## NEW FEATURES

* histogram2d/histogram2dcontour traces now respect the `colors` argument.

## BUG FIX

* Don't traceify by non-existant levels, fixes #735.

# 4.5.4 -- 27 September 2016

## BUG FIX

* Only insert missing values to differentiate groups when it's relevant.

# 4.5.3 -- 27 September 2016

## NEW FEATURES

* The `colorbar()` function gains a new `limits` arguments for controlling the colorscale
limits.

## BUG FIX

* The `z` is now required in `add_heatmap()`. If you want a `z` to be computed, use `add_histogram()`.

# 4.5.2 -- 23 September 2016

## NEW FEATURES

* The new argument, `split`, replaces the old functionality of the now deprecated `group` argument by creating one trace per value.

## BUG FIXES

* Passing plots to `subplot()` without a specified color (once again) match the coloring defaults supplied by plotly.js (see #724).

# 4.5.1 -- 23 September 2016

## NEW FEATURES

* A tibble with a list-column of plotly objects can now be input directly into `subplot()`

## BUG FIXES

* The `colorbar()` function now works on colorbars generated via `z` mapping.

# 4.5.0 -- 22 September 2016

## NEW FEATURES

* Added the `plot_mapbox()` and `plot_geo()` functions, which make it easier to work with the "scattermapbox", "scattergeo", and "choropleth" trace types. See the maps chapter of the plotly book for some examples -- <https://plotly-r.com/maps.html> 
* `subplot()` now accepts, and correctly scales mapbox objects.
* Added the `add_mesh3d()` and `add_pie()` functions as wrappers around the "mesh3d", and "pie" trace types.

## CHANGES

* The `add_scattergeo()` and `add_choropleth()` functions have been deprecated in favor of `plot_geo()`. 
* The `add_area(...)` function changed it's meaning from `add_lines(..., fill = 'tozeroy')` to a wrapper around the area trace <https://plotly.com/r/reference/#area>. This is more consistent with the naming conventions already in place for other `add_()` functions.
* `add_ribbons()` now shows points (instead of fill) on hover.

# 4.4.5 -- 19 September 2016

## NEW FEATURES

* Added a `rangeslider()` function to make it easier to add a range slider to the x-axis.
* Added a `colorbar()` function to make it easier to modify an automatically generated colorbar.

## BUG FIXES

* Bug fix for data arranging (introduced in 4.4.2).
* If the same (discrete) variable is mapped to two different aesthetics, redundant text is no longer generated in the legend entries (e.g., `plot_ly(mpg, x = ~cty, y = ~hwy, symbol = ~factor(cyl), color = ~factor(cyl))`)

# 4.4.4 -- 15 September 2016

## NEW FEATURES

* Added `inherit` argument for all `add_()` functions to avoid inheriting attributes from `plot_ly()`.
* Added the `add_fun()` function to add layers to a plot without modifying the original data associated with a plotly object.
* Added the `add_annotations()` function to make it easier to add annotations.
* Added the `layerData` argument to `ggplotly()` to make it possible to retrieve the data from a particular __ggplot2__ layer.

# 4.4.3 -- 15 September 2016

## CHANGES

* Downgrade to plotly.js v1.16.3 (which is proven to be a stable release and avoids #717) -- https://github.com/plotly/plotly.js/releases/tag/v1.17.2

# 4.4.2 -- 14 September 2016

## BUG FIXES

* Arrange data by trace index _before_ computing groups, fixes #716.

# 4.4.1 -- 14 September 2016

## BUG FIXES

* Restrict to atomic vectors when gathering data for training; otherwise, formulas referencing variables that don't exist in the data, but reference a function can cause problems.

# 4.4.0 -- 13 September 2016

## CHANGES 

* To enhance visibility of small markers, `marker.line.color` is now transparent by default.
* Upgraded to plotly.js v1.17.2 -- https://github.com/plotly/plotly.js/releases/tag/v1.17.2

## BUG FIXES

* It is now possible (again) to set/change attributes of autogenerated `marker.colorbar`.
* The `add_choropleth()` previously wasn't relaying the `z` attribute.
* Factors were being treated as characters in `add_segments()` (resulting in incorrect axis category order).
* No more error in `plot_ly()` when the number of traces is a multiple of ten.

# 4.3.7 -- 11 September 2016

## BUG FIXES

* `event_data()` now works inside shiny modules (#659). For an example, see <https://github.com/plotly/plotly.R/tree/master/inst/examples/shiny/event_data_modules>

# 4.3.6 -- 9 September 2016

## CHANGES

* Upgraded to plotly.js v1.17.0 -- https://github.com/plotly/plotly.js/releases/tag/v1.17.0

## BUG FIXES

* Fix for error handling in `add_bars()`.
* More careful logic for inferring data variables. Fixes #712

# 4.3.5 -- 5 September 2016

## NEW FEATURES

* The internal `as_widget()` function was exported to make it easier to convert a list
(adhering to the plotly spec) to a plotly htmlwidget object. This should only be needed when "manually" editing the underlying list object.
* Warnings about missing attributes now supply information about the relevant trace.

## CHANGES

* vignettes were removed and that documentation will now be hosted at <https://plotly.com/r/>

## BUG FIXES

* Get event_data() working with subplots. Fixes #663

# 4.3.4 -- 31 August 2016

## CHANGES

* Expressions yielding a ggplot2 object can now, once again, be provided to `plotlyOutput()`. In order to make this possible, `ggplotly()` now has a method for plotly objects (the identity function), and `ggplotly()` called on any expression provided to `plotlyOutput()`.

# 4.3.3 -- 29 August 2016

## BUG FIXES

* Bug fix for translation of ggplot2's `geom_text()`.

# 4.3.2 -- 26 August 2016

## NEW FEATURES

* The function `last_plot()` can now be used to retrieve the most recently _printed_ plotly graph. Thanks to this new feature, when `plotly_POST()` is called with no plotly object supplied, the most recently _printed_ plotly graph is sent to the users plotly account.

# 4.3.1 -- 23 August 2016

## CHANGES

* Upgraded to plotly.js v1.16.3 -- https://github.com/plotly/plotly.js/releases/tag/v1.16.3

# 4.3.0 -- 22 August 2016

## NEW FEATURES

* The `colors`/`symbols`/`linetypes` arguments now accept _named_ character vectors.
The names specify the domain (i.e., data values) and the values specify the range
(i.e., visual encodings). This is mainly useful to ensure a particular 
(discrete) data value is mapped to a particular visual attribute (yes, this is similar, in spirit, to ggplot2's `scale_*_manual()`).

## CHANGES

* Symbol and linetype palette defaults are now determined by `scales::shape_pal()` and `scales::linetype_pal()`.
* viridis is the default color scale for ordered factors.
* When mapping a character string to `color`/`symbol`/`linetype`, domain values are 
sorted alphabetically before scales are applied. Also, when mapping a factor to `color`/`symbol`/`linetype`, domain values are sorted according to their factor levels before scales are applied. This leads to more consistent (categorical axis ordering behaves similarly) and predictable (compared to having values sorted in the order in which they appear) behavior.

## BUG FIXES

# 4.2.1 -- 22 August 2016

## BUGFIX

* `alpha` is now applied when `color` isn't specified (fixes #658).

# 4.2.0 -- 11 August 2016

## CHANGES

* `plot_ly()` now orders the categories of a discrete x/y axis according the level ordering (if a factor) or alphabetical order (if a character string). Fixes #578.

# 4.1.1 -- 8 August 2016

## CHANGES

* Upgraded to plotly.js v1.16.1 -- https://github.com/plotly/plotly.js/releases/tag/v1.16.1

# 4.1.0 -- 27 June 2016

## NEW FEATURES

* `ggplotly()` gains a new `originalData` argument which allows one to attach either the original (global) data, or a "scaled"/"trained" version of the data used by __ggplot2__ to draw the graph (for a quick example, `ggplotly(qplot(1:10), originalData = FALSE) %>% plotly_data()`). 
* Hoverinfo is now shown for fill, instead of points, for several geoms (`geom_polygon()`/`geom_hex()`/`geom_rect()`/`geom_map()`). 
* If `stat_identity()` is used, group domain values are preserved and displayed in hoverinfo.
* New functions `hide_guides()`/`hide_legend()` were added (these work similarly to the existing `hide_colorbar()`) to simply the hiding of guides (i.e., legends/colorbars).

## BUG FIXES

* Legend titles (annotations) are no longer generated when no legend is displayed (#635, #607)
* Hoverinfo is no longer displayed if no tooltip variables are present in a layer (#563).
* Facets with 10 or more columns/rows should now render correctly (#640).
* x-axis anchors in `facet_wrap()` should now be correct.

## OTHER CHANGES

* Upgraded to plotly.js v1.15.0 -- https://github.com/plotly/plotly.js/releases/tag/v1.15.0

# 4.0.2 -- 25 June 2016

## BUG FIXES

* Bug fix for formulas evaluating to a logical vector (#650)

# 4.0.1 -- 14 June 2016

## BUG FIXES

* Duplicated values of positional attributes are no longer removed (bug was introduced in v4.0.0).

## OTHER CHANGES

* Upgraded to plotly.js v1.14.2 -- https://github.com/plotly/plotly.js/releases/tag/v1.14.2

# 4.0.0 -- 13 June 2016

## BREAKING CHANGES & IMPROVEMENTS:

* Formulas (instead of plain expressions) are now required when using variable mappings. For example, `plot_ly(mtcars, x = wt, y = mpg, color = vs)` should now be `plot_ly(mtcars, x = ~wt, y = ~mpg, color = ~vs)`. This is a major breaking change, but it is necessary to ensure that evaluation is correct in all contexts (as a result, `evaluate` argument is now deprecated as it is no longer needed). It also has the benefit of being easier to program with (i.e., writing your own custom functions that wrap `plot_ly()`) since it preserves [referential transparency](https://en.wikipedia.org/wiki/Referential_transparency). For more details, see the [lazyeval vignette](https://github.com/hadley/lazyeval/blob/master/vignettes/lazyeval.Rmd)
* The data structure used to represent plotly objects is now an htmlwidget object (instead of a data frame with a special attribute tracking visual mappings). As a result, the `as.widget()` function has deprecated, and [serialization/memory leak problems](https://github.com/rstudio/shiny/issues/1151) are no longer an issue. This change also implies that arbitrary data manipulation functions can no longer be intermingled inside a plot pipeline, but plotly methods for dplyr's data manipulation verbs are now provided (see `?plotly_data` for examples).
* The `group` variable mapping no longer create multiple traces, but instead defines "gaps" within a trace (fixes #418, #381, #577). Groupings should be declared via the new `group_by()` function (see `help(plotly_data)` for examples) instead of the `group` argument (which is now deprecated).
* `plot_ly()` now _initializes_ a plotly object (i.e., won't add a scatter trace by default), meaning that something like `plot_ly(x = 1:10, y = 1:10) %>% add_trace(y = 10:1)` creates one trace, instead of two. That being said, if you manually specify a trace type in `plot_ly()`, it will add a layer with that trace type (e.g. `plot_ly(x = 1:10, y = 1:10, type = "scatter") %>% add_trace(y = 10:1)` draws two scatter traces). If no trace type is provided, a sensible type is inferred from the supplied data, and automatically added (i.e., `plot_ly(x = rnorm(100))` now creates a histogram).
* The `inherit` argument is deprecated. Any arguments/attributes specified in `plot_ly()` will automatically be passed along to additional traces added via `add_trace()` (or any of it's `add_*()` siblings).
* Aesthetic scaling (e.g., `color`, `symbol`, `size`) is applied at the plot-level, instead of the trace level. 
* Size is no longer automatically included in hovertext (closes #549).

## NEW FEATURES & IMPROVEMENTS:

* Added `linetype`/`linetypes` arguments for mapping discrete variables to line types (works very much like the `symbol`/`symbols`).
* Scaling for aesthetics can be avoided via `I()` (closes #428). This is mainly useful for changing default appearance (e.g. `plot_ly(x = 1:10, y = 1:10, color = I("red"))`).
* Symbols and linetypes now recognize `pch` and `lty` values (e.g. `plot_ly(x = 1:25, y = 1:25, symbol = I(0:24))`)
* A new `alpha` argument controls the alpha transparency of `color` (e.g. `plot_ly(x = 1:10, y = 1:10, color = I("red"), alpha = 0.1)`).
* Added a `sizes` argument for controlling the range of marker size scaling.
* New `add_polygons()`/`add_ribbons()`/`add_area()`/`add_segments()`/`add_lines()`/`add_markers()`/`add_paths()`/`add_text()` functions provide a shorthand for common special cases of `add_trace()`.
* New `toWebGL()` function for easy conversion from SVG to WebGL.
* New `export()` function makes it easy to save plots as png/jpeg/pdf (fixes #311).
* Misspecified trace/layout attributes produce a warning.
* New `plotly_data()` function for returning/inspecting data frame(s) associated with a plotly object.
* New `plotly_json()` function for inspecting the data sent to plotly.js (as an R list or JSON).
* `layout()` is now a generic function and uses method dispatch to avoid conflicts with `graphics::layout()` (fixes #464).

## OTHER CHANGES:

* Upgraded to plotly.js v1.14.1 -- https://github.com/plotly/plotly.js/releases/tag/v1.14.1

3.6.5 -- 10 June 2016

IMPROVEMENT:

Multiple rows of facet strips will now be separated by <br> (i.e., line breaks) instead of ,. See #593.

3.6.4 -- 31 May 2016

BUG FIX:

embed_notebook() will no longer use a '.embed' extension in the iframe src attribute. See #613.

3.6.3 -- 24 May 2016

CHANGES:

Provided a better way of reexporting magrittr::`%>%`. See #597.

3.6.2 -- 24 May 2016

CHANGES: 

Removed unnecessary plyr dependency.

3.6.1 -- 23 May 2016

BUG FIX: 

Add a default method for plotly_build. Fixes #592.

3.6.0 -- 16 May 2016

NEW FEATURES & CHANGES:

* Many improvements to the subplot() function:
  * ggplot2 objects are now officially supported (#520).
  * Several new arguments allow one to synchronize x/y axes (#298), height/width (#376), hide/show x/y axis titles.
  * A list of plots can now be passed to the first argument.
  * A new vignette with examples and more explanation can be accessed via `vignette("subplot")`.

* ggplotly() is now a generic function with a method for ggmatrix objects.
* plotly_build() is now a generic function. 

BUG FIX: 

Column facet strips will no longer be drawn when there is only one column.

3.5.7 -- 13 May 2016

CHANGES:

Better defaults for defaultWidth/defaultHeight in the htmlwidget's sizing policy.

BUG FIX:

Pass knitr options to the named argument options. Fixes #582.

3.5.6 -- 12 May 2016

BUG FIX:

Use .embed suffix in iframe src attribute. Fixes #581.

3.5.5 -- 5 May 2016

CHANGES:

ggplotly() will now use plotly's layout.axisid.title (instead of 
layout.annotations) for axis titles on non-faceted plots. 
This will make for a better title placement experience (see #510).

BUG FIX:

Space for interior facet_wrap() strips are now accounted for.

3.5.4 -- 5 May 2016

BUG FIX:

gg2list() now returns an object of class "plotly_built" instead of "plotly"
to ensure a sensible print method is invoked.

3.5.3 -- 3 May 2016

CHANGES:

Upgrade to plotlyjs v1.10.1 -- https://github.com/plotly/plotly.js/releases/tag/v1.10.1

3.5.2 -- 2 May 2016

BUG FIX:

Added missing key properties in ggplotly() converter so selections can be accessible via event_data().

3.5.1 -- 26 Apr 2016

CHANGES:

Upgrade to plotlyjs v1.10.0 -- https://github.com/plotly/plotly.js/releases/tag/v1.10.0

Distinguish between "built" (plotly_built) and "non-built" (plotly_hash) plotly objects. See #562


3.5.0 -- 19 Apr 2016

NEW FEATURES:

The toRGB() function will now respect alpha channels in hex color codes and can recursively apply alpha. 

CHANGES:

The toRGB() function will always output color codes with an alpha channel (e.g. toRGB('black') is now 'rgba(0,0,0,1)' instead of 'rgb(0,0,0)')

3.4.15 -- 18 Apr 2016

BUGFIX:

The alpha in geom_smooth was incorrectly inheriting from other layers. See #551.

3.4.14 -- 15 Apr 2016

CHANGES:

Upgrade to plotlyjs v1.9.0 -- https://github.com/plotly/plotly.js/releases/tag/v1.9.0

3.4.13 -- 6 Apr 2016

BUGFIX:

In some cases, marker color was inheriting from the marker line color when
it shouldn't have. See ##537.

3.4.12 -- 5 Apr 2016

CHANGES:

Upgrade to plotlyjs v1.8.0 -- https://github.com/plotly/plotly.js/releases/tag/v1.8.0

3.4.11 -- 2 Apr 2016

BUGFIX:

Fix bug when altering modebar button defaults

3.4.10 -- 1 Apr 2016

BUGFIX:

Fix a geom_errorbar bug introduced in 3.4.9. See #513.

3.4.9 -- 25 Mar 2016

BUGFIX:

Upgrade to plotlyjs 1.7.0. Fixes #513

3.4.8 -- 23 Mar 2016

BUGFIX:

* Safeguard against null fields in selections. See #530.

3.4.7 -- 19 Mar 2016

BUGFIX:

* Added custom CSS which allows plotly to work nicely in ioslides.

3.4.6 -- 17 Mar 2016

NEW FEATURES:

The 'plotly_relayout' event is now accessible via the event_data() function.

Fixed #514.

3.4.5 -- 17 Mar 2016

BUGFIX:

Fixed #514.

3.4.4 -- 17 Mar 2016

BUGFIX:

Show discrete positional values in tooltip (see #515); better GeomTile conversion; pass plot object into layers2traces.

3.4.3 -- 14 Mar 2016

BUGFIX:

Custom facet labeller functions will now translate correctly. See #507.

3.4.2 -- 14 Mar 2016

BUGFIX:

Automatic resizing will now occur only when layout.autosize is true (the default). See #403.

3.4.1 -- 13 Mar 2016

BUGFIX:

Legend titles are now supported.

3.4.0 -- 12 Mar 2016

NEW FEATURES:

* geom_map() and geom_hex() are now supported.

CHANGES:

* The default value of the fileopt argument was changed from "new" to "overwrite".

BUGFIX:

* Made a number of bugfixes/improvements to hoverinfo & conversion of geom_tile()/geom_point().

3.3.1 -- 10 Mar 2016

CHANGES:

* Changed the mapping argument name to tooltip (which seems like a better name).

BUGFIX:

* Redundant legend entries are no longer shown.

3.2.1 -- 10 Mar 2016

BUGFIX:

* Proper formatting for date tooltips.

3.2.0 -- 10 Mar 2016

CHANGES:

* Legend titles no longer appear in legend entries.
* Tooltips now reflect aesthetic mappings. This makes it easier to decode 
data values from a given visual marking.

NEW FEATURES:

* geom_violin() is now supported.
* ggplotly() gains a mapping argument to control the set of aesthetics to appears in the tooltip as well as their order.

3.1.0 -- 8 Mar 2016

CHANGES:

* The "hidden" sharing option in plotly_POST() was renamed to "secret".
* The default value in the scale argument in plotly_IMAGE() is now 1.

3.0.0 -- 8 Mar 2016

NEW FEATURES:

* ggplotly() is now about 20x faster (it avoids calling ggplot_build() 20+ times). In some cases, it might be even faster since a lot of other redundant computation is avoided.

CHANGES:

* Instead of (trying to) translate both major and minor grid lines, we now translate only major grid lines. This generally produces a result closer to the actual ggplot2 result since ggplot2 doesn't draw ticks on minor grid lines.

BUG FIXES:

* ggplotly() now supports most of scale_*()/theme()/guides(). As a result, this fixes a lot of issues (#482, #481, #479, #476, #473, #460, #456, #454, #453, #447, #443, #434, #422, #421, #399, #379, #378, #357, #318, #316, #242, #232, #211, #203, #185, #184, #161). In order to support all of scale_x_*() an scale_y_*(), we always use linear axis types, and supply ticktext/tickvals to plotly.js. This has some unfortunate consequences on hoverformatting, which may be addressed in future releases of plotly.js -- https://github.com/plotly/plotly.js/issues/320

2.5.0 -- 1 Mar 2016

NEW FEATURES

* New event_data() function provides easy access to plotly events in shiny.
For an example, see https://github.com/ropensci/plotly/tree/master/inst/examples/plotlyEvents

* plot_ly() and ggplotly() gain a source argument to differentiate between 
plotly events in shiny apps with multiple plots. ggplotly() also gains width 
and height arguments.

CHANGES

The arguments filename, fileopt, world_readable in ggplotly() were removed as
they should be provided to plotly_POST() instead. 

2.4.4 -- 13 Feb 2016

as.widget() now returns htmlwidget objects untouched. See #449.

2.4.3 -- 11 Feb 2016

Ensure that we always return HTTPS links. Fixes #455

2.4.2 -- 9 Feb 2016

Fix for on-premise domain configuration. 

2.4.1 -- 2 Feb 2016

Attach base_url in as.widget() so it works in multiple contexts

2.4.0 -- 1 Feb 2016

* Pass plot configuration using ... to avoid conflicts in defaults/documentation
* Upgrade to plotly.js 1.5.1

2.3.4 -- 1 Feb 2016

Added a plotly_api_domain environment variable for configuring the API domain. Fixes #441

2.3.3 -- 27 Jan 2016

Bump axis number for each trace matching a panel number. fixes #318

2.3.2 -- 25 Jan 2016

More accurate list of data_array properties. Fixes #415

2.3.1 -- 25 Jan 2016

More accurate conversion of path width. Fixes #373.

2.3.0 -- 19 Jan 2016

Add sharing argument and deprecate world_readable. Fixes #332

2.2.4 -- 18 Jan 2016

Fix for error in embed_notebook(). See #409.

2.2.3 -- 18 Jan 2016

Fix for geom_vline(). See #402.

2.2.2 -- 18 Jan 2016

Fix bar orientation when we detect geom_bar() + coord_flip() in ggplotly(). Fixes #390.

2.2.1 -- 18 Jan 2016

Search for axis title in scene object. fixes #393.

2.2.0 -- 13 Jan 2016

The default for layout.hovermode is now 'closest' for non-line scatter traces

2.1.3 -- 12 Jan 2016

Fix size and alpha translation for geom_point. Fixes #386

2.1.2 -- 11 Jan 2016

Upgraded to plotlyjs 1.4.1. For a list of changes, see https://github.com/plotly/plotly.js/releases/tag/v1.4.1

2.1.1 -- 11 Jan 2016

Upgraded to plotlyjs 1.4. For a list of changes, see https://github.com/plotly/plotly.js/releases/tag/v1.4.0

2.1.0 -- 29 Dec 2015

plot_ly() now defaults to inherit=FALSE and plotly_build() is now idempotent. Fixes #280 and #277. See #368 for details.

2.0.19 -- 23 Dec 2015

Added as.widget() function for conveniency in converting plotly object to htmlwidget objects. See #294.

2.0.18 -- 22 Dec 2015

Fix #365

2.0.17 -- 22 Dec 2015

Fix #358

2.0.16 -- 18 Dec 2015

Require ggplot2 2.0.0 or higher. For details, see #269.

2.0.15 -- 13 Dec 2015

Fix #346

2.0.14 -- 13 Dec 2015

Fix #212

2.0.13 -- 12 Dec 2015

Fix #286

2.0.12 -- 11 Dec 2015

Fix #221

2.0.11 -- 11 Dec 2015

Fix #250

2.0.10 -- 10 Dec 2015

Fix #225

2.0.9 -- 10 Dec 2015

Fix #333

2.0.8 -- 10 Dec 2015

Fix a bug with geom_segment (see #321 & #228) 

2.0.7 -- 10 Dec 2015

Fix #233

2.0.6 -- 2 Dec 2015

Upgrade to plotlyjs 1.1.1. Fixes #319.

2.0.5 -- 1 Dec 2015

Fix for legend names. See #236.

2.0.4 -- 28 Nov 2015

Fix #313.

2.0.3 -- 18 Nov 2015

Fixed bug causing knitr options to be ignored. Also added VignetteBuilder to DESCRIPTION to vignette is available.

2.0.2 -- 17 Nov 2015

Using plotly_build() on a ggplot object should always return a plotly object

2.0.1 -- 17 Nov 2015

Better printing of server figures. Documentation and other fixes for initial CRAN release!

2.0.0 -- 2 Nov 2015

Added a dependency on htmlwidgets and 'offline' plots are now the default. If you want to create a figure on a plotly server, you need to use `plotly_POST()`. Also added a `config()` function to control the default appearance of the interactive plot

1.0.10 -- 3 Nov 2015

Fixed #292.

1.0.9 -- 28 Sep 2015

Fixed filename, fileopt arguments in plot_ly. Specifying the same filename will now overwrite the plot if it exists.

1.0.8 -- 14 Sep 2015

Added the plotly_IMAGES() function which interfaces to the images endpoint https://api.plot.ly/v2/#images

Details -> https://github.com/ropensci/plotly/pull/279

1.0.7 -- 26 Aug 2015

See https://github.com/ropensci/plotly/pull/275

1.0.6 -- 25 Aug 2015

Fix a bug with subplot domain calculations (see https://github.com/ropensci/plotly/pull/274)

1.0.5 -- 20 Aug 2015

Fix issue converting plotly offline markdown documents to HTML when using `markdown::markdownToHTML`

1.0.4 -- 14 Aug 2015

Bug fix for subplot. See #265

1.0.3 -- 7 Aug 2015

Improved legend positioning. See #241

1.0.2 -- 2 Aug 2015

* last_plot() will now look for the last plotly object; if not found, it will try to find the last ggplot object.
* Officially added the filename, fileopt, and world_readable arguments to plot_ly() and ggplotly().
* If plotly offline is not available, the shiny.launch.browser option is changed to open a web brower. See #245.
* Various namespace/documentation improvements for R CMD check.

1.0.1 -- 2 Aug 2015

Removed the stream() function as it wasn't ready to be included.

1.0.0 -- 31 July 2015

A major reworking of package internals which includes a few backwards incompatible changes.

Major changes include:

(1) New high-level grammar for expressing Plotly graphs from R (see the `plot_ly()`, `add_trace()`, `layout()`, and `style()` functions).
(2) New print methods which make it easier to create, modify, and embed Plotly graphs.
(3) Added a `subplot()` function for putting several graphs on a single page.
(4) Added the `renderPlotly()` and `plotlyOutput()` functions for embedding plotly graphs in shiny applications.
(5) Added `offline()` function for creating standalone HTML pages via Plotly Offline (see http://purchasing.plot.ly/)

For more details, see the new vignettes with `browseVignettes(package = "plotly")` and/or the pull request -> https://github.com/ropensci/plotly/pull/226

0.6.3 -- 2 June 2015

Add new tests inspired by the R Cookbook distributions #214

0.6.2 -- 19 May 2015

In geom_bar(stat = "identity"), sum y values if multiple for a given x.

0.6.1 -- 5 May 2015

Add test-cookbook-lines.R and fix bugs that showed up in those tests.

0.6 -- 4 May 2015

Let gg2list() return a figure object (backwards incompatible change).

0.5.29 -- 16 April 2015

geom_density() as filled area chart #202

0.5.28 -- 15 April 2015

Let ggplot handle histogram binning. Fix #198

0.5.27 -- 19 Mar 2015

Reimplement geom_ribbon as a basic polygon. Fix #191. Fix #192.

0.5.26 -- 18 Mar 2015

Implemented geom_rect #178

0.5.25 -- 10 March 2015

Implemented geom_smooth() #183

0.5.24 -- 10 March 2015

Implemented facet_wrap(scales="free") #167

0.5.23 -- 10 March 2015.

geom_ribbon() now respects alpha transparency

0.5.22 -- 2 March 2015.

Fixes for ylim() #171.

0.5.21 -- 23 February 2015.

Fixes for error bars and tick marks.

0.5.20 -- 9 February 2015.

Add alpha transparency to fill conversion.
Let geom_area support colour and fill aesthetics.

0.5.19 -- 23 January 2015.

Support class conversion such as as.Date() within ggplot code.

0.5.18 -- 22 January 2015.

Return proper filepath when filename contains directories.

0.5.17 -- 30 December 2014.

Support date-time binning in histograms.

0.5.16 -- 29 December 2014.

Support colour aesthetic in geom_text().

0.5.15 -- 19 December 2014.

Use proper RCurlOptions in get_figure() method.

0.5.14 -- 1 December 2014.

Make layers geom_line + geom_point only one trace in Plotly.

0.5.13 -- 27 November 2014.

Rename translation file and server endpoint parameter to be hip.

0.5.12 -- 12 November 2014.

Improve legend title position.

0.5.11 -- 11 November 2014.

Show legend title.

0.5.10 -- 7 November 2014.

Improve showlegend and fix legend’s `x` position.

0.5.9 -- 3 November 2014.

Default colours for geom_polygon().

0.5.8 -- 30 October 2014.

Support hline over a factor x range.
Default colours for geom_boxplot().

0.5.7 -- 29 October 2014.

Default colours for geom_area() and geom_ribbon().

0.5.6 -- 28 October 2014.

Convert line size faithfully.

0.5.5 -- 24 October 2014.

Support category histograms (with factors).

0.5.4 -- 22 October 2014.

Support conversion of geom_vline().

0.5.3 -- 21 October 2014.

Support conversion of geom_bar() with position_dodge().

0.5.2 -- 18 October 2014.

Support aesthetic shape in geom_path() and, hence, geom_line() (conversion).

0.5.1 -- 15 October 2014.

Do not show zero lines by default (as in ggplot2 plots).

0.5.0 -- 15 October 2014.

From now on, version numbers are meaningful again...
Many changes meanwhile, especially support for more geoms.

0.4 -- 7 April 2014.

Re-write geom to trace conversion code.

0.3.8 -- 21 March 2014.

ggplotly takes the last_plot() by default.

Support for ggplotly layout elements title, tickcolor, gridcolor,
showlegend, plot_bgcolor, paper_bgcolor, tickangle, axis titles, plot
border colors.

0.3.7 -- 14 March 2014.

For ggplotly:

- if on the command line, open a web browser (as before).

- if in knitr/Rmd in a chunk with plotly=TRUE, embed the plot.

0.3.6 -- 10 March 2014.

Merge ggplotly code.

0.3.5
# Contributing to plotly.js

## Opening issues

See the [opening issues template](https://github.com/ropensci/plotly/blob/master/.github/ISSUE_TEMPLATE.md)

## Development guidelines

If you'd like to contribute changes to plotly, we use [the GitHub flow](https://guides.github.com/introduction/flow/index.html) for proposing, submitting, reviewing, and accepting changes. If you aren't familiar with git and/or GitHub, we recommend studying these excellent free resources by [Hadley Wickham](http://r-pkgs.had.co.nz/git.html) and [Jenny Bryan](http://happygitwithr.com). We also prefer R coding style that adheres to <http://style.tidyverse.org/>.

If your pull request fixes a bug and/or implements a new feature, please [write a test via testthat](http://r-pkgs.had.co.nz/tests.html) to demonstrate it's working. Tests should generally check the return value of `plotly_build()`, but can also use **vdiffr**'s `expect_doppelganger()` to add a visual test! By default, `devtools::tests()` won't run the visual tests, but you *can* run visual tests on your machine via `Sys.setenv("VDIFFR" = "true"); devtools::test()`. That being said, false positives are likely to occur simply due to superfluous differences in your system environment, so we recommend running visual tests via docker.

### Running visual tests via docker

To ensure a consistent and reproducible environment, visual tests should be run against the [cpsievert/plotly-orca](https://hub.docker.com/r/cpsievert/plotly-orca/) docker image. If you add a new visual test and/or expect any differences, run a container like so:

```shell
git clone https://github.com/ropensci/plotly.git
cd plotly
docker run -v $(pwd):/home/plotly --privileged -p 3838:3838 cpsievert/plotly-orca
```

This will launch a shiny app for inspecting and validating any visual differences. To see the shiny app, open your browser to <http://0.0.0.0/3838>. If there are differences that look 'good', you should validate them via the shiny app. This will automatically copy over the new "baseline" figures over to your host machine (so that you can git add/commit/push the new baselines). If, for some reason, you want to just run the visual tests to see if they'll pass, do:

```shell
docker run -e VMODE="ci" -v $(pwd):/home/plotly --privileged cpsievert/plotly-orca
```

## Code of Conduct

We want to encourage a warm, welcoming, and safe environment for contributing to this project. See the [code of conduct](https://github.com/ropensci/plotly/blob/master/CONDUCT.md) for more information.
Please briefly describe your problem and what output you expect. If you have a question, please don't use this form, but instead ask on the community forum <http://community.plot.ly/c/api/r> or stackoverflow <http://stackoverflow.com>.

Please include a minimal reprex. The goal of a reprex is to make it as easy as possible for me to recreate your problem so that I can fix it. If you've never heard of a reprex before, start by reading <https://github.com/tidyverse/reprex#readme>, and follow the advice further down the page. Do NOT include session info unless it's explicitly asked for, or you've used `reprex::reprex(..., si = TRUE)` to hide it away.  

Delete these instructions once you have read them.

---

Brief description of the problem

```r
# insert reprex here
```
---
title: "Printing plotly objects in a knitr/rmarkdown doc"
output: html_document
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  message = FALSE,
  fig.width = 10,
  fig.height = 4,
  comment = "#>",
  collapse = TRUE
)
```

Printing anything created via the R package **plotly** should "just work" in a knitr/rmarkdown document -- including representations of things on your plotly account. However, the default print method that the package provides may not work for your purposes, so this document is designed to help you go beyond those defaults. This is especially useful for objects representing "remote files".

## Remote files

For example, if you create a plot or grid (via `api_create()`), then print the result, you get an HTML iframe pointing to that object on your account.

```{r}
library(plotly)
d <- api_create(mtcars)
d
```



As it turns out, the `api_create()` function returns a bunch of metadata about that file on your account. A nice way to inspect that information is to leverage the `jsonedit()` function from the **listviewer** package.

```{r}
listviewer::jsonedit(d)
```

Storing this information is a good idea since now you can [modify the "remote file"](https://api.plot.ly/v2/files#partial_update) at a later point. Let's rename the file using the "low-level" `api()` interface.

```{r}
nm <- paste(sample(LETTERS, 20), collapse = "-")
d2 <- api(
  file.path("files", d$fid), "PATCH", list(filename = nm)
)
identical(d2$filename, nm)
```

The `api_create()` function also understands how to "upload" ggplot2/plotly objects to the web platform. Printing in this case will again produce an HTML iframe pointing to the plot as it appears on the platform. 


```{r}
p <- api_create(qplot(1:10, 1:10))
p
```

The metadata returned for a plot is structured very much like a grid (i.e., what we saw previously), but now we can leverage some attributes unique to a plot, such as image urls.

```{r}
library(htmltools)
tags$img(src = p$image_urls$default)
```

## Downloading files

You can also download a plot or a grid already hosted on plotly's web platform (assuming they're public, or you have the proper crendentials). When you download a plot, it is converted to an htmlwidget, meaning that when you print it, the plot will render entirely locally (i.e., you don't need internet access for it to render).

```{r}
p <- api_download_plot(200, "cpsievert")
layout(
  p, title ="An htmlwidget version of <a href='https://plot.ly/~cpsievert/200'>this</a> plot"
)
```

There is no guarantee a plotly grid will _always_ map back to an R data frame, so `api_download_grid()` returns the abstract list representation of the data, but will _try to_ convert it to a data frame when printing.

```{r}
g <- api_download_grid(14681, "cpsievert")
g
```

```{r}
# note how the actual data is inside the 'preview' element
listviewer::jsonedit(g)
```



---
title: "Render plotly MathJax inside rmarkdown"
output: html_document
---

Some HTML-based MathJax:

$$ \alpha+\beta $$

You _could_ print this **plotly** graph with SVG-based rendering, but it would break the HTML-based rendering of **rmarkdown**!

```{r message = FALSE}
library(plotly)

p <- plotly_empty() %>%
  add_trace(x = 1, y = 1, text = TeX("\\alpha"), mode = "text", size = I(1000)) %>%
  config(mathjax = "cdn")
```

Instead, use something like the **widgetframe** package to create a responsive iframe which ensure the SVG-based rendering that plotly requires is done independently of **rmarkdown**'s HTML-based rendering.

```{r}
widgetframe::frameableWidget(p)
```

Or, do it the old-fashioned way: save your plotly graph to an HTML file via `htmlwidgets::saveWidget()` then use an HTML `<iframe>`

```{r}
htmlwidgets::saveWidget(p, "my-plotly-plot.html")
```


<iframe src="my-plotly-plot.html" width="100%" height="400" id="igraph" scrolling="no" seamless="seamless" frameBorder="0"> </iframe>
---
title: "Using plotly with onRender"
author: "Carson Sievert"
date: "`r Sys.Date()`"
output: html_document
---

```{r, message = FALSE, warning = FALSE}
library(plotly)
library(htmlwidgets)

set.seed(1056)

nPatients <- 50
nVisits <- 10

df <- data.frame(
  fev1_perc = rnorm(n = nPatients * nVisits, mean = 100, sd = 10),
  uin = rep(seq(nPatients), each = nVisits),
  visit = rep(seq(nVisits), nPatients)
)
c1 <- list(color = toRGB("steelblue", 0.5))
c2 <- list(color = toRGB("orange", 0.5))

# The color mapping is used only to generate multiple traces 
# (which makes it easier to highlight lines on the JS side)
df %>%
  plot_ly(
    x = ~visit, y = ~fev1_perc, split = ~factor(uin), marker = c1, line = c2
  ) %>% 
  layout(hovermode = "closest", showlegend = FALSE) %>%
  onRender('
    function(el, x) { 
      var graphDiv = document.getElementById(el.id);
      // reduce the opacity of every trace except for the hover one
      el.on("plotly_hover", function(e) { 
        var traces = [];
        for (var i = 0; i < x.data.length; i++) {
          if (i !== e.points[0].curveNumber) traces.push(i);
        }
        Plotly.restyle(graphDiv, "opacity", 0.2, traces);
      })
     el.on("plotly_unhover", function(e) { 
       var traces = [];
       for (var i = 0; i < x.data.length; i++) traces.push(i);
       Plotly.restyle(graphDiv, "opacity", 1, traces);
     })
    } 
  ')
```
---
title: "Flex Dashboard"
output: 
  flexdashboard::flex_dashboard:
    orientation: rows
---


```{r setup, include=FALSE}
library(plotly)
library(maps)
knitr::opts_chunk$set(message = FALSE)
```

Rows {data-height:600}
------------------------------------------------------------------------------

### Chart A

```{r}
# This example modifies code from Hadley Wickham -- https://gist.github.com/hadley/233134
# It also uses data from Nathan Yau's flowingdata site -- http://flowingdata.com/
unemp <- read.csv("http://datasets.flowingdata.com/unemployment09.csv")
names(unemp) <- c("id", "state_fips", "county_fips", "name", "year", 
                  "?", "?", "?", "rate")
unemp$county <- tolower(gsub(" County, [A-Z]{2}", "", unemp$name))
unemp$state <- gsub("^.*([A-Z]{2}).*$", "\\1", unemp$name)
county_df <- map_data("county")
names(county_df) <- c("long", "lat", "group", "order", "state_name", "county")
county_df$state <- state.abb[match(county_df$state_name, tolower(state.name))]
county_df$state_name <- NULL
state_df <- map_data("state")
choropleth <- merge(county_df, unemp, by = c("state", "county"))
choropleth <- choropleth[order(choropleth$order), ]
choropleth$rate_d <- cut(choropleth$rate, breaks = c(seq(0, 10, by = 2), 35))

# provide a custom tooltip to plotly with the county name and actual rate
choropleth$text <- with(choropleth, paste0("County: ", name, "<br>Rate: ", rate))
p <- ggplot(choropleth, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = rate_d, text = text), 
               colour = alpha("white", 1/2), size = 0.2) + 
  geom_polygon(data = state_df, colour = "white", fill = NA) +
  scale_fill_brewer(palette = "PuRd") + theme_void()
# just show the text aesthetic in the tooltip
ggplotly(p, tooltip = "text")
```

### Chart B

```{r}
crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
crimesm <- tidyr::gather(crimes, variable, value, -state)
states_map <- map_data("state")
g <- ggplot(crimesm, aes(map_id = state)) +
  geom_map(aes(fill = value), map = states_map) +
  expand_limits(x = states_map$long, y = states_map$lat) +
  facet_wrap( ~ variable) + theme_void()
ggplotly(g)
```

Rows {data-height:400}
------------------------------------------------------------------------------


### Chart C

```{r}
m <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
  stat_density_2d() + xlim(0.5, 6) + ylim(40, 110)
ggplotly(m)
```

### Chart D

```{r}
m <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
  xlim(0.5, 6) + ylim(40, 110)
ggplotly(m)
```


### Chart E

```{r}
m <- ggplot(faithful, aes(x = eruptions, y = waiting)) + geom_hex() 
ggplotly(m)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layers2traces.R
\name{geom2trace}
\alias{geom2trace}
\title{Convert a "basic" geoms to a plotly.js trace.}
\usage{
geom2trace(data, params, p)
}
\arguments{
\item{data}{the data returned by \code{plotly::to_basic}.}

\item{params}{parameters for the geom, statistic, and 'constant' aesthetics}

\item{p}{a ggplot2 object (the conversion may depend on scales, for instance).}
}
\description{
This function makes it possible to convert ggplot2 geoms that
are not included with ggplot2 itself. Users shouldn't need to use
this function. It exists purely to allow other package authors to write
their own conversion method(s).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{plotly_POST}
\alias{plotly_POST}
\title{Create/Modify plotly graphs}
\usage{
plotly_POST(
  x = last_plot(),
  filename = NULL,
  fileopt = "overwrite",
  sharing = c("public", "private", "secret"),
  ...
)
}
\arguments{
\item{x}{either a ggplot object, a plotly object, or a list.}

\item{filename}{character string describing the name of the plot in your
plotly account. Use / to specify directories. If a directory path does not
exist it will be created. If this argument is not specified and the title
of the plot exists, that will be used for the filename.}

\item{fileopt}{character string describing whether to create a "new" plotly,
"overwrite" an existing plotly, "append" data to existing plotly,
or "extend" it.}

\item{sharing}{If 'public', anyone can view this graph. It will appear in
your profile and can appear in search engines. You do not need to be
logged in to Plotly to view this chart.
If 'private', only you can view this plot. It will not appear in the
Plotly feed, your profile, or search engines. You must be logged in to
Plotly to view this graph. You can privately share this graph with other
Plotly users in your online Plotly account and they will need to be logged
in to view this plot.
If 'secret', anyone with this secret link can view this chart. It will
not appear in the Plotly feed, your profile, or search engines.
If it is embedded inside a webpage or an IPython notebook, anybody who is
viewing that page will be able to view the graph.
You do not need to be logged in to view this plot.}

\item{...}{not used}
}
\description{
Deprecated: see \code{\link[=api_create]{api_create()}}.
}
\seealso{
\code{\link[=plot_ly]{plot_ly()}}, \code{\link[=signup]{signup()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highlight.R
\name{attrs_selected}
\alias{attrs_selected}
\title{Specify attributes of selection traces}
\usage{
attrs_selected(opacity = 1, ...)
}
\arguments{
\item{opacity}{a number between 0 and 1 specifying the overall opacity of
the selected trace}

\item{...}{other trace attributes attached to the selection trace.}
}
\description{
By default the name of the selection trace derives from the selected values.
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplotly.R
\name{ggplotly}
\alias{ggplotly}
\title{Convert ggplot2 to plotly}
\usage{
ggplotly(
  p = ggplot2::last_plot(),
  width = NULL,
  height = NULL,
  tooltip = "all",
  dynamicTicks = FALSE,
  layerData = 1,
  originalData = TRUE,
  source = "A",
  ...
)
}
\arguments{
\item{p}{a ggplot object.}

\item{width}{Width of the plot in pixels (optional, defaults to automatic sizing).}

\item{height}{Height of the plot in pixels (optional, defaults to automatic sizing).}

\item{tooltip}{a character vector specifying which aesthetic mappings to show
in the tooltip. The default, "all", means show all the aesthetic mappings
(including the unofficial "text" aesthetic). The order of variables here will
also control the order they appear. For example, use
\code{tooltip = c("y", "x", "colour")} if you want y first, x second, and
colour last.}

\item{dynamicTicks}{should plotly.js dynamically generate axis tick labels?
Dynamic ticks are useful for updating ticks in response to zoom/pan
interactions; however, they can not always reproduce labels as they
would appear in the static ggplot2 image.}

\item{layerData}{data from which layer should be returned?}

\item{originalData}{should the "original" or "scaled" data be returned?}

\item{source}{a character string of length 1. Match the value of this string
with the source argument in \code{\link[=event_data]{event_data()}} to retrieve the
event data corresponding to a specific plot (shiny apps can have multiple plots).}

\item{...}{arguments passed onto methods.}
}
\description{
This function converts a \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}} object to a
plotly object.
}
\details{
Conversion of relative sizes depends on the size of the current
graphics device (if no device is open, width/height of a new (off-screen)
device defaults to 640/480). In other words, \code{height} and
\code{width} must be specified at runtime to ensure sizing is correct.
For examples on how to specify the output container's \code{height}/\code{width} in a
shiny app, see \code{plotly_example("shiny", "ggplotly_sizing")}.
}
\examples{
\dontrun{
# simple example
ggpenguins <- qplot(bill_length_mm , body_mass_g, 
data = palmerpenguins::penguins, color = species)
ggplotly(ggpenguins)

data(canada.cities, package = "maps")
viz <- ggplot(canada.cities, aes(long, lat)) +
  borders(regions = "canada") +
  coord_equal() +
  geom_point(aes(text = name, size = pop), colour = "red", alpha = 1/2)
ggplotly(viz, tooltip = c("text", "size"))

# linked scatterplot brushing
d <- highlight_key(mtcars)
qplot(data = d, x = mpg, y = wt) \%>\%
  subplot(qplot(data = d, x = mpg, y = vs)) \%>\% 
  layout(title = "Click and drag to select points") \%>\%
  highlight("plotly_selected")


# more brushing (i.e. highlighting) examples
demo("crosstalk-highlight-ggplotly", package = "plotly")

# client-side linked brushing in a scatterplot matrix
highlight_key(palmerpenguins::penguins) \%>\%
  GGally::ggpairs(aes(colour = Species), columns = 1:4) \%>\%
  ggplotly(tooltip = c("x", "y", "colour")) \%>\%
  highlight("plotly_selected")
}

}
\references{
\url{https://plotly.com/ggplot2/}
}
\seealso{
\code{\link[=plot_ly]{plot_ly()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{rangeslider}
\alias{rangeslider}
\title{Add a range slider to the x-axis}
\usage{
rangeslider(p, start = NULL, end = NULL, ...)
}
\arguments{
\item{p}{plotly object.}

\item{start}{a start date/value.}

\item{end}{an end date/value.}

\item{...}{these arguments are documented here
\url{https://plotly.com/r/reference/#layout-xaxis-rangeslider}}
}
\description{
Add a range slider to the x-axis
}
\examples{

plot_ly(x = time(USAccDeaths), y = USAccDeaths) \%>\% 
  add_lines() \%>\%
  rangeslider()
  
d <- tibble::tibble(
  time = seq(as.Date("2016-01-01"), as.Date("2016-08-31"), by = "days"),
  y = rnorm(seq_along(time))
 )
 
plot_ly(d, x = ~time, y = ~y) \%>\%
  add_lines() \%>\%
  rangeslider(d$time[5], d$time[50])
  

}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.api_grid}
\alias{print.api_grid}
\title{Print a plotly grid object}
\usage{
\method{print}{api_grid}(x, ...)
}
\arguments{
\item{x}{a plotly grid object}

\item{...}{additional arguments (currently ignored)}
}
\description{
Print a plotly grid object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kaleido.R
\name{print.kaleidoScope}
\alias{print.kaleidoScope}
\title{Print method for kaleido}
\usage{
\method{print}{kaleidoScope}(x, ...)
}
\arguments{
\item{x}{a \code{\link[=kaleido]{kaleido()}} object.}

\item{...}{currently unused.}
}
\description{
S3 method for \code{\link[=kaleido]{kaleido()}}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{plotly_empty}
\alias{plotly_empty}
\title{Create a complete empty plotly graph.}
\usage{
plotly_empty(...)
}
\arguments{
\item{...}{arguments passed onto \code{\link[=plot_ly]{plot_ly()}}}
}
\description{
Useful when used with \code{\link[=subplot]{subplot()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layers2traces.R
\name{to_basic}
\alias{to_basic}
\title{Convert a geom to a "basic" geom.}
\usage{
to_basic(data, prestats_data, layout, params, p, ...)
}
\arguments{
\item{data}{the data returned by \code{ggplot2::ggplot_build()}.}

\item{prestats_data}{the data before statistics are computed.}

\item{layout}{the panel layout.}

\item{params}{parameters for the geom, statistic, and 'constant' aesthetics}

\item{p}{a ggplot2 object (the conversion may depend on scales, for instance).}

\item{...}{currently ignored}
}
\description{
This function makes it possible to convert ggplot2 geoms that
are not included with ggplot2 itself. Users shouldn't need to use
this function. It exists purely to allow other package authors to write
their own conversion method(s).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{event_register}
\alias{event_register}
\title{Register a shiny input value}
\usage{
event_register(p, event = NULL)
}
\arguments{
\item{p}{a plotly object.}

\item{event}{The type of plotly event. All supported events are listed in the
function signature above (i.e., the usage section).}
}
\description{
Register a shiny input value
}
\seealso{
\link{event_data}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_exports.R
\name{api_create}
\alias{api_create}
\alias{api_create.plotly}
\alias{api_create.ggplot}
\alias{api_create.data.frame}
\alias{api_download_plot}
\alias{api_download_grid}
\alias{api}
\title{Tools for working with plotly's REST API (v2)}
\usage{
api_create(
  x = last_plot(),
  filename = NULL,
  fileopt = c("overwrite", "new"),
  sharing = c("public", "private", "secret"),
  ...
)

\method{api_create}{plotly}(
  x = last_plot(),
  filename = NULL,
  fileopt = "overwrite",
  sharing = "public",
  ...
)

\method{api_create}{ggplot}(
  x = last_plot(),
  filename = NULL,
  fileopt = "overwrite",
  sharing = "public",
  ...
)

\method{api_create}{data.frame}(x, filename = NULL, fileopt = "overwrite", sharing = "public", ...)

api_download_plot(id, username)

api_download_grid(id, username)

api(endpoint = "/", verb = "GET", body = NULL, ...)
}
\arguments{
\item{x}{An R object to hosted on plotly's web platform.
Can be a plotly/ggplot2 object or a \link{data.frame}.}

\item{filename}{character vector naming file(s). If \code{x} is a plot,
can be a vector of length 2 naming both the plot AND the underlying grid.}

\item{fileopt}{character string describing whether to "overwrite" existing
files or ensure "new" file(s) are always created.}

\item{sharing}{If 'public', anyone can view this graph. It will appear in
your profile and can appear in search engines. You do not need to be
logged in to Plotly to view this chart.
If 'private', only you can view this plot. It will not appear in the
Plotly feed, your profile, or search engines. You must be logged in to
Plotly to view this graph. You can privately share this graph with other
Plotly users in your online Plotly account and they will need to be logged
in to view this plot.
If 'secret', anyone with this secret link can view this chart. It will
not appear in the Plotly feed, your profile, or search engines.
If it is embedded inside a webpage or an IPython notebook, anybody who is
viewing that page will be able to view the graph.
You do not need to be logged in to view this plot.}

\item{...}{For \code{api()}, these arguments are passed onto
\code{\link[httr:RETRY]{httr::RETRY()}}. For \code{api_create()}, these arguments are
included in the body of the HTTP request.}

\item{id}{a filename id.}

\item{username}{a plotly username.}

\item{endpoint}{the endpoint (i.e., location) for the request.
To see a list of all available endpoints, call \code{api()}.
Any relevant query parameters should be included here (see examples).}

\item{verb}{name of the HTTP verb to use (as in, \code{\link[httr:RETRY]{httr::RETRY()}}).}

\item{body}{body of the HTTP request(as in, \code{\link[httr:RETRY]{httr::RETRY()}}).
If this value is not already converted to JSON
(via \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}), it uses the internal \code{to_JSON()}
to ensure values are "automatically unboxed" (i.e., vec.}
}
\description{
Convenience functions for working with version 2 of plotly's REST API.
Upload R objects to a plotly account via \code{api_create()} and download
plotly objects via \code{api_download_plot()}/\code{api_download_grid()}.
For anything else, use \code{api()}.
}
\examples{

\dontrun{

# ------------------------------------------------------------
# api_create() makes it easy to upload ggplot2/plotly objects
# and/or data frames to your plotly account
# ------------------------------------------------------------

# A data frame creates a plotly "grid". Printing one will take you 
# to the it's web address so you can start creating!
(m <- api_create(mtcars))

# A plotly/ggplot2 object create a plotly "plot".
p <- plot_ly(mtcars, x = ~factor(vs))
(r <- api_create(p))

# api_create() returns metadata about the remote "file". Here is
# one way you could use that metadata to download a plot for local use:
fileID <- strsplit(r$file$fid, ":")[[1]]
layout(
  api_download_plot(fileID[2], fileID[1]),
  title = sprintf("Local version of <a href='\%s'>this</a> plot", r$file$web_url)
)

------------------------------------------------------------
# The api() function provides a low-level interface for performing 
# any action at any endpoint! It always returns a list.
# ------------------------------------------------------------

# list all the endpoints
api()

# search the entire platform!
# see https://api.plot.ly/v2/search
api("search?q=overdose")
api("search?q=plottype:pie trump fake")

# these examples will require a user account
usr <- Sys.getenv("plotly_username", NA)
if (!is.na(usr)) {
  # your account info https://api.plot.ly/v2/#users
  api(sprintf("users/\%s", usr))
  # your folders/files https://api.plot.ly/v2/folders#user
  api(sprintf("folders/home?user=\%s", usr))
}

# Retrieve a specific file https://api.plot.ly/v2/files#retrieve
api("files/cpsievert:14681")

# change the filename https://api.plot.ly/v2/files#update
# (note: this won't work unless you have proper credentials to the relevant account)
api("files/cpsievert:14681", "PATCH", list(filename = "toy file")) 

# Copy a file https://api.plot.ly/v2/files#lookup
api("files/cpsievert:14681/copy", "POST")

# Create a folder https://api.plot.ly/v2/folders#create
api("folders", "POST", list(path = "/starts/at/root/and/ends/here"))

}

}
\references{
\url{https://api.plot.ly/v2}
}
\seealso{
\code{\link[=signup]{signup()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{hide_legend}
\alias{hide_legend}
\title{Hide legend}
\usage{
hide_legend(p)
}
\arguments{
\item{p}{a plotly object.}
}
\description{
Hide legend
}
\examples{

p <- plot_ly(mtcars, x = ~wt, y = ~cyl, color = ~factor(cyl))
hide_legend(p)
}
\seealso{
\code{\link[=hide_colorbar]{hide_colorbar()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.api_grid_local}
\alias{print.api_grid_local}
\title{Print a plotly grid object}
\usage{
\method{print}{api_grid_local}(x, ...)
}
\arguments{
\item{x}{a plotly grid object}

\item{...}{additional arguments (currently ignored)}
}
\description{
Print a plotly grid object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dev.R
\name{plotly_json}
\alias{plotly_json}
\title{Inspect JSON sent to plotly.js}
\usage{
plotly_json(p = last_plot(), jsonedit = interactive(), pretty = TRUE, ...)
}
\arguments{
\item{p}{a plotly or ggplot object.}

\item{jsonedit}{use \link[listviewer:jsonedit]{listviewer::jsonedit} to view the JSON?}

\item{pretty}{adds indentation whitespace to JSON output. Can be TRUE/FALSE
or a number specifying the number of spaces to indent. See \link[jsonlite:prettify]{jsonlite::prettify}.}

\item{...}{other options passed onto \link[listviewer:jsonedit]{listviewer::jsonedit}}
}
\description{
This function is useful for obtaining/viewing/debugging JSON
sent to plotly.js.
}
\examples{
  
plotly_json(plot_ly())
plotly_json(plot_ly(), FALSE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly.R
\name{plot_dendro}
\alias{plot_dendro}
\title{Plot an interactive dendrogram}
\usage{
plot_dendro(d, set = "A", xmin = -50, height = 500, width = 500, ...)
}
\arguments{
\item{d}{a dendrogram object}

\item{set}{defines a crosstalk group}

\item{xmin}{minimum of the range of the x-scale}

\item{height}{height}

\item{width}{width}

\item{...}{arguments supplied to \code{\link[=subplot]{subplot()}}}
}
\description{
This function takes advantage of nested key selections to implement an
interactive dendrogram. Selecting a node selects all the labels (i.e. leafs)
under that node.
}
\examples{

\dontrun{
hc <- hclust(dist(USArrests), "ave")
dend1 <- as.dendrogram(hc)
plot_dendro(dend1, height = 600) \%>\% 
  hide_legend() \%>\% 
  highlight(persistent = TRUE, dynamic = TRUE)
}

}
\seealso{
\code{\link[=plot_ly]{plot_ly()}}, \code{\link[=plot_mapbox]{plot_mapbox()}}, \code{\link[=ggplotly]{ggplotly()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highlight.R
\name{highlight}
\alias{highlight}
\title{Query graphical elements in multiple linked views}
\usage{
highlight(
  p,
  on = "plotly_click",
  off,
  persistent = getOption("persistent", FALSE),
  dynamic = FALSE,
  color = NULL,
  selectize = FALSE,
  defaultValues = NULL,
  opacityDim = getOption("opacityDim", 0.2),
  selected = attrs_selected(),
  debounce = 0,
  ...
)
}
\arguments{
\item{p}{a plotly visualization.}

\item{on}{turn on a selection on which event(s)? To disable on events
altogether, use \code{NULL}. Currently the following are supported:
\itemize{
\item \code{'plotly_click'}
\item \code{'plotly_hover'}
\item \code{'plotly_selected'}: triggered through rectangular
(layout.dragmode = 'select') or lasso (layout.dragmode = 'lasso') brush.
}}

\item{off}{turn off a selection on which event(s)? To disable off
events altogether, use \code{NULL}. Currently the following are supported:
\itemize{
\item \code{'plotly_doubleclick'}: triggered on a double mouse click while
(layout.dragmode = 'zoom') or (layout.dragmode = 'pan')
\item \code{'plotly_deselect'}: triggered on a double mouse click while
(layout.dragmode = 'select') or (layout.dragmode = 'lasso')
\item \code{'plotly_relayout'}: triggered whenever axes are rescaled
(i.e., clicking the home button in the modebar) or whenever the height/width
of the plot changes.
}}

\item{persistent}{should selections persist (i.e., accumulate)? We often
refer to the default (\code{FALSE}) as a 'transient' selection mode;
which is recommended, because one may switch from 'transient' to
'persistent' selection by holding the shift key.}

\item{dynamic}{should a widget for changing selection colors be included?}

\item{color}{character string of color(s) to use for
highlighting selections. See \code{\link[=toRGB]{toRGB()}} for valid color
specifications. If \code{NULL} (the default), the color of selected marks
are not altered.}

\item{selectize}{provide a selectize.js widget for selecting keys? Note that
the label used for this widget derives from the groupName of the SharedData object.}

\item{defaultValues}{a vector of values for setting a "default selection".
These values should match the key attribute.}

\item{opacityDim}{a number between 0 and 1 used to reduce the
opacity of non-selected traces (by multiplying with the existing opacity).}

\item{selected}{attributes of the selection, see \code{\link[=attrs_selected]{attrs_selected()}}.}

\item{debounce}{amount of time to wait before firing an event (in milliseconds).
The default of 0 means do not debounce at all.
Debouncing is mainly useful when \code{on = "plotly_hover"} to avoid firing too many events
when users clickly move the mouse over relevant graphical marks.}

\item{...}{currently not supported.}
}
\description{
This function sets a variety of options for brushing (i.e., highlighting)
multiple plots. These options are primarily designed for linking
multiple plotly graphs, and may not behave as expected when linking
plotly to another htmlwidget package via crosstalk. In some cases,
other htmlwidgets will respect these options, such as persistent selection
in leaflet (see \code{demo("highlight-leaflet", package = "plotly")}).
}
\examples{

# These examples are designed to show you how to highlight/brush a *single*
# view. For examples of multiple linked views, see `demo(package = "plotly")` 

d <- highlight_key(txhousing, ~city)
p <- ggplot(d, aes(date, median, group = city)) + geom_line()
gg <- ggplotly(p, tooltip = "city") 
highlight(gg, dynamic = TRUE)

# supply custom colors to the brush 
cols <- toRGB(RColorBrewer::brewer.pal(3, "Dark2"), 0.5)
highlight(gg, on = "plotly_hover", color = cols, dynamic = TRUE)

# Use attrs_selected() for complete control over the selection appearance
# note any relevant colors you specify here should override the color argument
s <- attrs_selected(
  showlegend = TRUE,
  mode = "lines+markers",
  marker = list(symbol = "x")
)

highlight(layout(gg, showlegend = TRUE), selected = s)

}
\references{
\url{https://plotly-r.com/client-side-linking.html}
}
\seealso{
\code{\link[=attrs_selected]{attrs_selected()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{hide_guides}
\alias{hide_guides}
\title{Hide guides (legends and colorbars)}
\usage{
hide_guides(p)
}
\arguments{
\item{p}{a plotly object.}
}
\description{
Hide guides (legends and colorbars)
}
\seealso{
\code{\link[=hide_legend]{hide_legend()}}, \code{\link[=hide_colorbar]{hide_colorbar()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{config}
\alias{config}
\title{Set the default configuration for plotly}
\usage{
config(
  p,
  ...,
  cloud = FALSE,
  showSendToCloud = cloud,
  locale = NULL,
  mathjax = NULL
)
}
\arguments{
\item{p}{a plotly object}

\item{...}{these arguments are documented at
\url{https://github.com/plotly/plotly.js/blob/master/src/plot_api/plot_config.js}}

\item{cloud}{deprecated. Use \code{showSendToCloud} instead.}

\item{showSendToCloud}{include the send data to cloud button?}

\item{locale}{locale to use. See \href{https://github.com/plotly/plotly.js/tree/master/dist#to-include-localization}{here} for more info.}

\item{mathjax}{add \href{https://github.com/plotly/plotly.js/tree/master/dist#to-support-mathjax}{MathJax rendering support}.
If \code{"cdn"}, mathjax is loaded externally (meaning an internet connection is needed for
TeX rendering). If \code{"local"}, the PLOTLY_MATHJAX_PATH environment variable must be
set to the location (a local file path) of MathJax. IMPORTANT: \strong{plotly} uses SVG-based
mathjax rendering which doesn't play nicely with HTML-based rendering
(e.g., \strong{rmarkdown} documents and \strong{shiny} apps). To leverage both types of rendering,
you must \verb{<iframe>} your plotly graph(s) into the larger document
(see \href{https://github.com/plotly/plotly.R/blob/master/inst/examples/rmd/MathJax/index.Rmd}{here}
for an \strong{rmarkdown} example and
\href{https://github.com/plotly/plotly.R/blob/master/inst/examples/rmd/MathJax/index.Rmd}{here} for a \strong{shiny} example).}
}
\description{
Set the default configuration for plotly
}
\examples{

# remove the plotly logo and collaborate button from modebar
config(plot_ly(), displaylogo = FALSE, collaborate = FALSE)

# enable mathjax
# see more examples at https://plotly.com/r/LaTeX/
plot_ly(x = c(1, 2, 3, 4), y = c(1, 4, 9, 16)) \%>\%
  layout(title = TeX("\\\\text{Some mathjax: }\\\\alpha+\\\\beta x")) \%>\%
  config(mathjax = "cdn")

# change the language used to render date axes and on-graph text 
# (e.g., modebar buttons)
today <- Sys.Date()
x <- seq.Date(today, today + 360, by = "day")
p <- plot_ly(x = x, y = rnorm(length(x))) \%>\%
  add_lines()

# japanese
config(p, locale = "ja")
# german
config(p, locale = "de")
# spanish
config(p, locale = "es")
# chinese
config(p, locale = "zh-CN")

}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplotly.R
\name{bbox}
\alias{bbox}
\title{Estimate bounding box of a rotated string}
\usage{
bbox(txt = "foo", angle = 0, size = 12)
}
\arguments{
\item{txt}{a character string of length 1}

\item{angle}{sets the angle of the tick labels with respect to the
horizontal (e.g., \code{tickangle} of -90 draws the tick labels vertically)}

\item{size}{vertical size of a character}
}
\description{
Estimate bounding box of a rotated string
}
\references{
https://www.dropbox.com/s/nc6968prgw8ne4w/bbox.pdf?dl=0
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add_data}
\alias{add_data}
\title{Add data to a plotly visualization}
\usage{
add_data(p, data = NULL)
}
\arguments{
\item{p}{a plotly visualization}

\item{data}{a data frame.}
}
\description{
Add data to a plotly visualization
}
\examples{

plot_ly() \%>\% add_data(economics) \%>\% add_trace(x = ~date, y = ~pce)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partial_bundles.R
\name{partial_bundle}
\alias{partial_bundle}
\title{Use a partial bundle of plotly.js}
\usage{
partial_bundle(p, type = "auto", local = TRUE, minified = TRUE)
}
\arguments{
\item{p}{a plotly object.}

\item{type}{name of the (partial) bundle. The default, \code{'auto'}, attempts to
find the smallest single bundle that can render \code{p}. If no single partial bundle
can render \code{p}, then the full bundle is used.}

\item{local}{whether or not to download the partial bundle so that it can be
viewed later without an internet connection.}

\item{minified}{whether or not to use a minified js file (non-minified file can be useful for debugging plotly.js)}
}
\description{
Leveraging plotly.js' partial bundles can lead to smaller file sizes
and faster rendering. The full list of available bundles, and the
trace types that they support, are available
\href{https://github.com/plotly/plotly.js/blob/master/dist/README.md#partial-bundles}{here}
}
\details{
WARNING: use this function with caution when rendering multiple
plotly graphs on a single website. That's because, if multiple plotly.js
bundles are used, the most recent bundle will override the other bundles.
See the examples section for an example.
}
\examples{

# ----------------------------------------------------------------------
# This function is always safe to use when rendering a single 
# plotly graph. In this case, we get a 3x file reduction.
# ----------------------------------------------------------------------

\dontrun{
library(plotly)
p <- plot_ly(x = 1:10, y = 1:10) \%>\% add_markers()
save_widget <- function(p, f) {
  owd <- setwd(dirname(f))
  on.exit(setwd(owd))
  htmlwidgets::saveWidget(p, f)
  mb <- round(file.info(f)$size / 1e6, 3)
  message("File is: ", mb," MB")
}
f1 <- tempfile(fileext = ".html")
f2 <- tempfile(fileext = ".html")
save_widget(p, f1)
save_widget(partial_bundle(p), f2)

# ----------------------------------------------------------------------
# But, since plotly.js bundles override one another, 
# be careful when putting multiple graphs in a larger document!
# Note how the surface (part of the gl3d bundle) renders, but the 
# heatmap (part of the cartesian bundle) doesn't...
# ----------------------------------------------------------------------

library(htmltools)
p1 <- plot_ly(z = ~volcano) \%>\% 
  add_heatmap() \%>\%
  partial_bundle()
p2 <- plot_ly(z = ~volcano) \%>\% 
  add_surface() \%>\%
  partial_bundle()
browsable(tagList(p1, p2))
}

}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{plotly-shiny}
\alias{plotly-shiny}
\alias{plotlyOutput}
\alias{renderPlotly}
\title{Shiny bindings for plotly}
\usage{
plotlyOutput(
  outputId,
  width = "100\%",
  height = "400px",
  inline = FALSE,
  reportTheme = TRUE
)

renderPlotly(expr, env = parent.frame(), quoted = FALSE)
}
\arguments{
\item{outputId}{output variable to read from}

\item{width, height}{Must be a valid CSS unit (like \code{"100\%"},
\code{"400px"}, \code{"auto"}) or a number, which will be coerced to a
string and have \code{"px"} appended. Note that, for height, using "auto"
or "100\%" generally will not work as expected, because of how
height is computed with HTML/CSS.}

\item{inline}{use an inline (\code{span()}) or block container
(\code{div()}) for the output}

\item{reportTheme}{whether or not to report CSS styles (if a sufficient
version of shiny and htmlwidgets is available).}

\item{expr}{An expression that generates a plotly}

\item{env}{The environment in which to evaluate \code{expr}.}

\item{quoted}{Is \code{expr} a quoted expression (with \code{quote()})? This
is useful if you want to save an expression in a variable.}
}
\description{
Output and render functions for using plotly within Shiny
applications and interactive Rmd documents.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplotly.R
\name{gg2list}
\alias{gg2list}
\title{Convert a ggplot to a list.}
\usage{
gg2list(
  p,
  width = NULL,
  height = NULL,
  tooltip = "all",
  dynamicTicks = FALSE,
  layerData = 1,
  originalData = TRUE,
  source = "A",
  ...
)
}
\arguments{
\item{p}{ggplot2 plot.}

\item{width}{Width of the plot in pixels (optional, defaults to automatic sizing).}

\item{height}{Height of the plot in pixels (optional, defaults to automatic sizing).}

\item{tooltip}{a character vector specifying which aesthetic tooltips to show in the
tooltip. The default, "all", means show all the aesthetic tooltips
(including the unofficial "text" aesthetic).}

\item{dynamicTicks}{accepts the following values: \code{FALSE}, \code{TRUE}, \code{"x"}, or \code{"y"}.
Dynamic ticks are useful for updating ticks in response to zoom/pan/filter
interactions; however, there is no guarantee they reproduce axis tick text
as they would appear in the static ggplot2 image.}

\item{layerData}{data from which layer should be returned?}

\item{originalData}{should the "original" or "scaled" data be returned?}

\item{source}{a character string of length 1. Match the value of this string
with the source argument in \code{\link[=event_data]{event_data()}} to retrieve the
event data corresponding to a specific plot (shiny apps can have multiple plots).}

\item{...}{currently not used}
}
\value{
a 'built' plotly object (list with names "data" and "layout").
}
\description{
Convert a ggplot to a list.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/last_plot.R
\name{last_plot}
\alias{last_plot}
\title{Retrieve the last plot to be modified or created.}
\usage{
last_plot()
}
\description{
Retrieve the last plot to be modified or created.
}
\seealso{
\code{\link[ggplot2:last_plot]{ggplot2::last_plot()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/style.R
\name{style}
\alias{style}
\title{Modify trace(s)}
\usage{
style(p, ..., traces = NULL)
}
\arguments{
\item{p}{A plotly visualization.}

\item{...}{Visual properties.}

\item{traces}{numeric vector. Which traces should be modified? By default,
attributes place in \code{...} will be applied to every trace.}
}
\description{
Modify trace(s) of an existing plotly visualization. Useful when used in
conjunction with \code{\link[=get_figure]{get_figure()}}.
}
\examples{

# style() is especially useful in conjunction with ggplotly()
# It allows you to leverage the underlying plotly.js library to change 
# the return result of ggplotly()
(p <- ggplotly(qplot(data = mtcars, wt, mpg, geom = c("point", "smooth"))))

# removes hoverinfo for the line/ribbon traces (use `plotly_json()` to verify!)
style(p, hoverinfo = "none", traces = c(2, 3))

# another example with plot_ly() instead of ggplotly()
marker <- list(
  color = "red",
  line = list(
    width = 20, 
    color = "black"
 )
)
(p <- plot_ly(x = 1:10, y = 1:10, marker = marker))

# note how the entire (marker) object is replaced if a list is provided
style(p, marker = list(line = list(color = "blue")))

# similar to plotly.js, you can update a particular attribute like so 
# https://github.com/plotly/plotly.js/issues/1866#issuecomment-314115744
style(p, marker.line.color = "blue") 
# this clobbers the previously supplied marker.line.color
style(p, marker.line = list(width = 2.5), marker.size = 10)

}
\seealso{
\code{\link[=api_download_plot]{api_download_plot()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toRGB.R
\name{toRGB}
\alias{toRGB}
\title{Convert R colours to RGBA hexadecimal colour values}
\usage{
toRGB(x, alpha = 1)
}
\arguments{
\item{x}{see the \code{col} argument in \code{col2rgb} for valid specifications}

\item{alpha}{alpha channel on 0-1 scale}
}
\value{
hexadecimal colour value (if is.na(x), return "transparent" for compatibility with Plotly)
}
\description{
Convert R colours to RGBA hexadecimal colour values
}
\examples{

toRGB("steelblue") 
# [1] "rgba(70,130,180,1)"

m <- list(
  color = toRGB("red"),
  line = list(
    color = toRGB("black"), 
    width = 19
  )
)

plot_ly(x = 1, y = 1, marker = m)

}
\seealso{
\code{\link[=showRGB]{showRGB()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.api_plot}
\alias{print.api_plot}
\title{Print a plot on plotly's platform}
\usage{
\method{print}{api_plot}(x, ...)
}
\arguments{
\item{x}{a plotly figure object}

\item{...}{additional arguments (currently ignored)}
}
\description{
Print a plot on plotly's platform
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{toWebGL}
\alias{toWebGL}
\title{Convert trace types to WebGL}
\usage{
toWebGL(p)
}
\arguments{
\item{p}{a plotly or ggplot object.}
}
\description{
Convert trace types to WebGL
}
\examples{

# currently no bargl trace type
toWebGL(ggplot() + geom_bar(aes(1:10)))
toWebGL(qplot(1:10, 1:10))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export.R
\name{export}
\alias{export}
\title{Export a plotly graph to a static file}
\usage{
export(p = last_plot(), file = "plotly.png", selenium = NULL, ...)
}
\arguments{
\item{p}{a plotly or ggplot object.}

\item{file}{a filename. The file type is inferred from the file extension.
Valid extensions include 'jpeg' | 'png' | 'webp' | 'svg' | 'pdf'}

\item{selenium}{used only when \code{p} is a WebGL plot or the output
format is 'webp' or 'svg'. Should be an object of class "rsClientServer"
returned by \code{RSelenium::rsDriver}.}

\item{...}{if \code{p} is non-WebGL and the output file format is
jpeg/png/pdf arguments are passed along to \code{webshot::webshot()}.
Otherwise, they are ignored.}
}
\description{
This function is in the process of being deprecated (use \link{orca} instead).
}
\details{
For SVG plots, a screenshot is taken via \code{webshot::webshot()}.
Since \code{phantomjs} (and hence \code{webshot}) does not support WebGL,
the RSelenium package is used for exporting WebGL plots.
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{hobbs}
\alias{hobbs}
\title{Hobbs data}
\format{
A data frame with three variables: \code{r}, \code{t},
\code{nms}.
}
\usage{
hobbs
}
\description{
Description TBD.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/embed.R
\name{embed_notebook}
\alias{embed_notebook}
\title{Embed a plot as an iframe into a Jupyter Notebook}
\usage{
embed_notebook(x, width = NULL, height = NULL, file = NULL)
}
\arguments{
\item{x}{a plotly object}

\item{width}{attribute of the iframe. If \code{NULL}, the width in
\code{plot_ly} is used. If that is also \code{NULL}, '100\%' is the default.}

\item{height}{attribute of the iframe. If \code{NULL}, the height in
\code{plot_ly} is used. If that is also \code{NULL}, '400px' is the default.}

\item{file}{deprecated.}
}
\description{
Embed a plot as an iframe into a Jupyter Notebook
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly_example.R
\name{plotly_example}
\alias{plotly_example}
\title{Run a plotly example(s)}
\usage{
plotly_example(type = c("demo", "shiny", "rmd"), name, edit = TRUE, ...)
}
\arguments{
\item{type}{the type of example}

\item{name}{the name of the example (valid names depend on \code{type}).}

\item{edit}{whether to open the relevant source files using \link{file.edit}. Only relevant if \code{type} is \code{"shiny"} or \code{"rmd"}.}

\item{...}{arguments passed onto the suitable method.}
}
\description{
Provides a unified interface for running demos, shiny apps, and Rmd documents
which are bundled with the package.
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{raster2uri}
\alias{raster2uri}
\title{Encode a raster object as a data URI}
\usage{
raster2uri(r, ...)
}
\arguments{
\item{r}{an object coercable to a raster object via \code{\link[=as.raster]{as.raster()}}}

\item{...}{arguments passed onto \code{\link[=as.raster]{as.raster()}}.}
}
\description{
Encode a raster object as a data URI, which is suitable for
use with \code{layout()} \href{https://plotly.com/r/reference/#layout-images}{images}.
This is especially convenient for embedding raster images on a plot in
a self-contained fashion (i.e., so they don't depend on external URL links).
}
\examples{

# a red gradient (from ?as.raster)
r <- as.raster(matrix(hcl(0, 80, seq(50, 80, 10)), nrow = 4, ncol = 5))
plot(r)

# embed the raster as an image
plot_ly(x = 1, y = 1) \%>\% 
  layout(
    images = list(list(
     source = raster2uri(r),
     xref = "paper", 
     yref = "paper", 
     x = 0, y = 0, 
     sizex = 0.5, sizey = 0.5, 
     xanchor = "left", yanchor = "bottom"
  ))
 ) 
}
\references{
\url{https://plotly-r.com/embedding-images.html}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{layout}
\alias{layout}
\title{Modify the layout of a plotly visualization}
\usage{
layout(p, ..., data = NULL)
}
\arguments{
\item{p}{A plotly object.}

\item{...}{Arguments to the layout object. For documentation,
see \url{https://plotly.com/r/reference/#Layout_and_layout_style_objects}}

\item{data}{A data frame to associate with this layout (optional). If not
provided, arguments are evaluated using the data frame in \code{\link[=plot_ly]{plot_ly()}}.}
}
\description{
Modify the layout of a plotly visualization
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/animate.R
\name{animation_opts}
\alias{animation_opts}
\alias{animation}
\alias{animation_slider}
\alias{animation_button}
\title{Animation configuration options}
\usage{
animation_opts(
  p,
  frame = 500,
  transition = frame,
  easing = "linear",
  redraw = TRUE,
  mode = "immediate"
)

animation_slider(p, hide = FALSE, ...)

animation_button(p, ..., label)
}
\arguments{
\item{p}{a plotly object.}

\item{frame}{The amount of time between frames (in milliseconds).
Note that this amount should include the \code{transition}.}

\item{transition}{The duration of the smooth transition between
frames (in milliseconds).}

\item{easing}{The type of transition easing. See the list of options here
\url{https://github.com/plotly/plotly.js/blob/master/src/plots/animation_attributes.js}}

\item{redraw}{Trigger a redraw of the plot at completion of the transition?
A redraw may significantly impact performance, but may be necessary to
update graphical elements that can't be transitioned.}

\item{mode}{Describes how a new animate call interacts with currently-running
animations. If \code{immediate}, current animations are interrupted and
the new animation is started. If \code{next}, the current frame is allowed
to complete, after which the new animation is started. If \code{afterall}
all existing frames are animated to completion before the new animation
is started.}

\item{hide}{remove the animation slider?}

\item{...}{for \code{animation_slider}, attributes are passed to a special
layout.sliders object tied to the animation frames.
The definition of these attributes may be found here
\url{https://github.com/plotly/plotly.js/blob/master/src/components/sliders/attributes.js}
For \code{animation_button}, arguments are passed to a special
layout.updatemenus button object tied to the animation
\url{https://github.com/plotly/plotly.js/blob/master/src/components/updatemenus/attributes.js}}

\item{label}{a character string used for the animation button's label}
}
\description{
Animations can be created by either using the \code{frame} argument in
\code{\link[=plot_ly]{plot_ly()}} or the (unofficial) \code{frame} ggplot2 aesthetic in
\code{\link[=ggplotly]{ggplotly()}}. By default, animations populate a play button
and slider component for controlling the state of the animation
(to pause an animation, click on a relevant location on the slider bar).
Both the play button and slider component transition between frames according
rules specified by \code{\link[=animation_opts]{animation_opts()}}.
}
\examples{

df <- data.frame(
  x = c(1, 2, 2, 1, 1, 2),
  y = c(1, 2, 2, 1, 1, 2),
  z = c(1, 1, 2, 2, 3, 3)
)
plot_ly(df) \%>\%
  add_markers(x = 1.5, y = 1.5) \%>\%
  add_markers(x = ~x, y = ~y, frame = ~z)

# it's a good idea to remove smooth transitions when there is
# no relationship between objects in each view
plot_ly(mtcars, x = ~wt, y = ~mpg, frame = ~cyl) \%>\%
  animation_opts(transition = 0)

# works the same way with ggplotly
if (interactive()) {
  p <- ggplot(txhousing, aes(month, median)) +
    geom_line(aes(group = year), alpha = 0.3) +
    geom_smooth() +
    geom_line(aes(frame = year, ids = month), color = "red") +
    facet_wrap(~ city)
 
  ggplotly(p, width = 1200, height = 900) \%>\%
    animation_opts(1000)
}

  
#' # for more, see https://plotly.com/r/animating-views.html

}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signup.R
\name{signup}
\alias{signup}
\title{Create a new plotly account.}
\usage{
signup(username, email, save = TRUE)
}
\arguments{
\item{username}{Desired username.}

\item{email}{Desired email.}

\item{save}{If request is successful, should the username & API key be
automatically stored as an environment variable in a .Rprofile?}
}
\value{
\itemize{
\item api_key key to use with the api
\item tmp_pw temporary password to access your plotly account
}
}
\description{
A sign up interface to plotly through the R Console.
}
\examples{
\dontrun{
# You need a plotly username and API key to communicate with the plotly API.

# If you don't already have an API key, you can obtain one with a valid
# username and email via signup().
s <- signup('anna.lyst', 'anna.lyst@plot.ly')

# If you already have a username and API key, please create the following
# environment variables:
Sys.setenv("plotly_username" = "me")
Sys.setenv("plotly_api_key" = "mykey")
# You can also change the default domain if you have a plotly server.
Sys.setenv("plotly_domain" = "http://mydomain.com")

# If you want to automatically load these environment variables when you
# start R, you can put them inside your ~/.Rprofile 
# (see help(.Rprofile) for more details)

}
}
\references{
https://plotly.com/rest/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subplots.R
\name{subplot}
\alias{subplot}
\title{View multiple plots in a single view}
\usage{
subplot(
  ...,
  nrows = 1,
  widths = NULL,
  heights = NULL,
  margin = 0.02,
  shareX = FALSE,
  shareY = FALSE,
  titleX = shareX,
  titleY = shareY,
  which_layout = "merge"
)
}
\arguments{
\item{...}{One of the following
\itemize{
\item any number of plotly/ggplot2 objects.
\item a list of plotly/ggplot2 objects.
\item a tibble with one list-column of plotly/ggplot2 objects.
}}

\item{nrows}{number of rows for laying out plots in a grid-like structure.
Only used if no domain is already specified.}

\item{widths}{relative width of each column on a 0-1 scale. By default all
columns have an equal relative width.}

\item{heights}{relative height of each row on a 0-1 scale. By default all
rows have an equal relative height.}

\item{margin}{either a single value or four values (all between 0 and 1).
If four values are provided, the first is used as the left margin, the second
is used as the right margin, the third is used as the top margin, and the
fourth is used as the bottom margin.
If a single value is provided, it will be used as all four margins.}

\item{shareX}{should the x-axis be shared amongst the subplots?}

\item{shareY}{should the y-axis be shared amongst the subplots?}

\item{titleX}{should x-axis titles be retained?}

\item{titleY}{should y-axis titles be retained?}

\item{which_layout}{adopt the layout of which plot? If the default value of
"merge" is used, layout options found later in the sequence of plots will
override options found earlier in the sequence. This argument also accepts a
numeric vector specifying which plots to consider when merging.}
}
\value{
A plotly object
}
\description{
View multiple plots in a single view
}
\examples{

# pass any number of plotly objects to subplot()
p1 <- plot_ly(economics, x = ~date, y = ~uempmed)
p2 <- plot_ly(economics, x = ~date, y = ~unemploy)
subplot(p1, p2, p1, p2, nrows = 2, margin = 0.05)

#'  # anchor multiple traces on the same legend entry
 p1 <- add_lines(p1, color = I("black"), name = "1st", legendgroup = "1st")
 p2 <- add_lines(p2, color = I("red"), name = "2nd", legendgroup = "2nd")
 
 subplot(
   p1, style(p1, showlegend = FALSE),
   p2, style(p2, showlegend = FALSE),
   nrows = 2, margin = 0.05
 )

# or pass a list
economics_long \%>\%
  split(.$variable) \%>\%
  lapply(function(d) plot_ly(d, x = ~date, y = ~value)) \%>\%
  subplot(nrows = NROW(.), shareX = TRUE)
  
# or pass a tibble with a list-column of plotly objects
economics_long \%>\%
  group_by(variable) \%>\%
  do(p = plot_ly(., x = ~date, y = ~value)) \%>\%
  subplot(nrows = NROW(.), shareX = TRUE)
  
# learn more at https://plotly.com/r/subplot.html

}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{event_data}
\alias{event_data}
\title{Access plotly user input event data in shiny}
\usage{
event_data(
  event = c("plotly_hover", "plotly_unhover", "plotly_click", "plotly_doubleclick",
    "plotly_selected", "plotly_selecting", "plotly_brushed", "plotly_brushing",
    "plotly_deselect", "plotly_relayout", "plotly_restyle", "plotly_legendclick",
    "plotly_legenddoubleclick", "plotly_clickannotation", "plotly_afterplot",
    "plotly_sunburstclick"),
  source = "A",
  session = shiny::getDefaultReactiveDomain(),
  priority = c("input", "event")
)
}
\arguments{
\item{event}{The type of plotly event. All supported events are listed in the
function signature above (i.e., the usage section).}

\item{source}{a character string of length 1. Match the value of this string
with the \code{source} argument in \code{\link[=plot_ly]{plot_ly()}} (or \code{\link[=ggplotly]{ggplotly()}}) to respond to
events emitted from that specific plot.}

\item{session}{a shiny session object (the default should almost always be used).}

\item{priority}{the priority of the corresponding shiny input value.
If equal to \code{"event"}, then \code{\link[=event_data]{event_data()}} always triggers re-execution,
instead of re-executing only when the relevant shiny input value changes
(the default).}
}
\description{
This function must be called within a reactive shiny context.
}
\examples{
\dontrun{
plotly_example("shiny", "event_data")
}
}
\references{
\itemize{
\item \url{https://plotly-r.com/linking-views-with-shiny.html#shiny-plotly-inputs}
\item \url{https://plotly.com/javascript/plotlyjs-function-reference/}
}
}
\seealso{
\link{event_register}, \link{event_unregister}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly_build.R
\name{plotly_build}
\alias{plotly_build}
\title{'Build' (i.e., evaluate) a plotly object}
\usage{
plotly_build(p, registerFrames = TRUE)
}
\arguments{
\item{p}{a ggplot object, or a plotly object, or a list.}

\item{registerFrames}{should a frame trace attribute be interpreted as frames in an animation?}
}
\description{
This generic function creates the list object sent to plotly.js
for rendering. Using this function can be useful for overriding defaults
provided by \code{ggplotly}/\code{plot_ly} or for debugging rendering
errors.
}
\examples{

p <- plot_ly(economics, x = ~date, y = ~pce)
# the unevaluated plotly object
str(p)
# the evaluated data
str(plotly_build(p)$x$data)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group2NA.R
\name{group2NA}
\alias{group2NA}
\title{Separate groups with missing values}
\usage{
group2NA(
  data,
  groupNames = "group",
  nested = NULL,
  ordered = NULL,
  retrace.first = inherits(data, "GeomPolygon")
)
}
\arguments{
\item{data}{a data frame.}

\item{groupNames}{character vector of grouping variable(s)}

\item{nested}{other variables that group should be nested
(i.e., ordered) within.}

\item{ordered}{a variable to arrange by (within nested & groupNames). This
is useful primarily for ordering by x}

\item{retrace.first}{should the first row of each group be appended to the
last row? This is useful for enclosing polygons with lines.}
}
\value{
a data.frame with rows ordered by: \code{nested},
then \code{groupNames}, then \code{ordered}. As long as \code{groupNames}
contains valid variable names, new rows will also be inserted to separate
the groups.
}
\description{
This function is used internally by plotly, but may also be useful to some
power users. The details section explains when and why this function is useful.
}
\details{
If a group of scatter traces share the same non-positional characteristics
(i.e., color, fill, etc), it is more efficient to draw them as a single trace
with missing values that separate the groups (instead of multiple traces),
In this case, one should also take care to make sure
\href{https://plotly.com/r/reference/#scatter-connectgaps}{connectgaps}
is set to \code{FALSE}.
}
\examples{

# note the insertion of new rows with missing values 
group2NA(mtcars, "vs", "cyl")

# need to group lines by city somehow!
plot_ly(txhousing, x = ~date, y = ~median) \%>\% add_lines()

# instead of using group_by(), you could use group2NA()
tx <- group2NA(txhousing, "city")
plot_ly(tx, x = ~date, y = ~median) \%>\% add_lines()

# add_lines() will ensure paths are sorted by x, but this is equivalent
tx <- group2NA(txhousing, "city", ordered = "date")
plot_ly(tx, x = ~date, y = ~median) \%>\% add_paths()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly.R
\name{plot_geo}
\alias{plot_geo}
\title{Initiate a plotly-geo object}
\usage{
plot_geo(data = data.frame(), ..., offline = FALSE)
}
\arguments{
\item{data}{A data frame (optional).}

\item{...}{arguments passed along to \code{\link[=plot_ly]{plot_ly()}}.}

\item{offline}{whether or not to include geo assets so that the map
can be viewed with or without an internet connection. The plotlyGeoAssets
package is required for this functionality.}
}
\description{
Use this function instead of \code{\link[=plot_ly]{plot_ly()}} to initialize
a plotly-geo object. This enforces the entire plot so use
the scattergeo trace type, and enables higher level geometries
like \code{\link[=add_polygons]{add_polygons()}} to work
}
\examples{

map_data("world", "canada") \%>\%
  group_by(group) \%>\%
  plot_geo(x = ~long, y = ~lat) \%>\%
  add_markers(size = I(1))

}
\seealso{
\code{\link[=plot_ly]{plot_ly()}}, \code{\link[=plot_mapbox]{plot_mapbox()}}, \code{\link[=ggplotly]{ggplotly()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add_fun}
\alias{add_fun}
\title{Apply function to plot, without modifying data}
\usage{
add_fun(p, fun, ...)
}
\arguments{
\item{p}{a plotly object.}

\item{fun}{a function. Should take a plotly object as input and return a
modified plotly object.}

\item{...}{arguments passed to \code{fun}.}
}
\description{
Useful when you need two or more layers that apply a summary statistic
to the original data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly_data.R
\name{highlight_key}
\alias{highlight_key}
\title{Highlight/query data based on primary key}
\usage{
highlight_key(x, ...)
}
\arguments{
\item{x}{a plotly visualization or a \code{data.frame}.}

\item{...}{arguments passed to \code{crosstalk::SharedData$new()}}
}
\value{
An object of class \link[crosstalk:SharedData]{crosstalk::SharedData}
}
\description{
This function simply creates an object of class \link[crosstalk:SharedData]{crosstalk::SharedData}.
The reason it exists is to make it easier to teach others how to leverage
its functionality in plotly. It also makes it more discoverable if one
is already aware of \link{highlight}.
}
\seealso{
\link{highlight}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add_trace}
\alias{add_trace}
\alias{add_markers}
\alias{add_text}
\alias{add_paths}
\alias{add_lines}
\alias{add_segments}
\alias{add_polygons}
\alias{add_sf}
\alias{add_table}
\alias{add_ribbons}
\alias{add_image}
\alias{add_area}
\alias{add_pie}
\alias{add_bars}
\alias{add_histogram}
\alias{add_histogram2d}
\alias{add_histogram2dcontour}
\alias{add_heatmap}
\alias{add_contour}
\alias{add_boxplot}
\alias{add_surface}
\alias{add_mesh}
\alias{add_scattergeo}
\alias{add_choropleth}
\title{Add trace(s) to a plotly visualization}
\usage{
add_trace(p, ..., data = NULL, inherit = TRUE)

add_markers(p, x = NULL, y = NULL, z = NULL, ..., data = NULL, inherit = TRUE)

add_text(
  p,
  x = NULL,
  y = NULL,
  z = NULL,
  text = NULL,
  ...,
  data = NULL,
  inherit = TRUE
)

add_paths(p, x = NULL, y = NULL, z = NULL, ..., data = NULL, inherit = TRUE)

add_lines(p, x = NULL, y = NULL, z = NULL, ..., data = NULL, inherit = TRUE)

add_segments(
  p,
  x = NULL,
  y = NULL,
  xend = NULL,
  yend = NULL,
  ...,
  data = NULL,
  inherit = TRUE
)

add_polygons(p, x = NULL, y = NULL, ..., data = NULL, inherit = TRUE)

add_sf(p, ..., x = ~x, y = ~y, data = NULL, inherit = TRUE)

add_table(p, ..., rownames = TRUE, data = NULL, inherit = TRUE)

add_ribbons(
  p,
  x = NULL,
  ymin = NULL,
  ymax = NULL,
  ...,
  data = NULL,
  inherit = TRUE
)

add_image(p, z = NULL, colormodel = NULL, ..., data = NULL, inherit = TRUE)

add_area(p, r = NULL, theta = NULL, t = NULL, ..., data = NULL, inherit = TRUE)

add_pie(p, values = NULL, labels = NULL, ..., data = NULL, inherit = TRUE)

add_bars(p, x = NULL, y = NULL, ..., data = NULL, inherit = TRUE)

add_histogram(p, x = NULL, y = NULL, ..., data = NULL, inherit = TRUE)

add_histogram2d(
  p,
  x = NULL,
  y = NULL,
  z = NULL,
  ...,
  data = NULL,
  inherit = TRUE
)

add_histogram2dcontour(
  p,
  x = NULL,
  y = NULL,
  z = NULL,
  ...,
  data = NULL,
  inherit = TRUE
)

add_heatmap(p, x = NULL, y = NULL, z = NULL, ..., data = NULL, inherit = TRUE)

add_contour(p, z = NULL, ..., data = NULL, inherit = TRUE)

add_boxplot(p, x = NULL, y = NULL, ..., data = NULL, inherit = TRUE)

add_surface(p, z = NULL, ..., data = NULL, inherit = TRUE)

add_mesh(p, x = NULL, y = NULL, z = NULL, ..., data = NULL, inherit = TRUE)

add_scattergeo(p, ...)

add_choropleth(p, z = NULL, ..., data = NULL, inherit = TRUE)
}
\arguments{
\item{p}{a plotly object}

\item{...}{Arguments (i.e., attributes) passed along to the trace \code{type}.
See \code{\link[=schema]{schema()}} for a list of acceptable attributes for a given trace \code{type}
(by going to \code{traces} -> \code{type} -> \code{attributes}). Note that attributes
provided at this level may override other arguments
(e.g. \code{plot_ly(x = 1:10, y = 1:10, color = I("red"), marker = list(color = "blue"))}).}

\item{data}{A data frame (optional) or \link[crosstalk:SharedData]{crosstalk::SharedData} object.}

\item{inherit}{inherit attributes from \code{\link[=plot_ly]{plot_ly()}}?}

\item{x}{the x variable.}

\item{y}{the y variable.}

\item{z}{a numeric matrix (unless \code{\link[=add_image]{add_image()}}, which wants a raster object, see \code{\link[=as.raster]{as.raster()}}).}

\item{text}{textual labels.}

\item{xend}{"final" x position (in this context, x represents "start")}

\item{yend}{"final" y position (in this context, y represents "start")}

\item{rownames}{whether or not to display the rownames of \code{data}.}

\item{ymin}{a variable used to define the lower boundary of a polygon.}

\item{ymax}{a variable used to define the upper boundary of a polygon.}

\item{colormodel}{Sets the colormodel for image traces if \code{z} is not a raster object.
If \code{z} is a raster object (see \code{\link[=as.raster]{as.raster()}}), the \code{'rgba'} colormodel is always used.}

\item{r}{Sets the radial coordinates.}

\item{theta}{Sets the angular coordinates.}

\item{t}{Deprecated. Use \code{theta} instead.}

\item{values}{the value to associated with each slice of the pie.}

\item{labels}{the labels (categories) corresponding to \code{values}.}
}
\description{
Add trace(s) to a plotly visualization
}
\examples{

# the `plot_ly()` function initiates an object, and if no trace type
# is specified, it sets a sensible default
p <- plot_ly(economics, x = ~date, y = ~uempmed)
p

# some `add_*()` functions are a specific case of a trace type
# for example, `add_markers()` is a scatter trace with mode of markers
add_markers(p)

# scatter trace with mode of text
add_text(p, text = "\%")

# scatter trace with mode of lines 
add_paths(p)

# like `add_paths()`, but ensures points are connected according to `x`
add_lines(p)

# if you prefer to work with plotly.js more directly, can always
# use `add_trace()` and specify the type yourself
add_trace(p, type = "scatter", mode = "markers+lines")

# mappings provided to `plot_ly()` are "global", but can be overwritten
plot_ly(economics, x = ~date, y = ~uempmed, color = I("red"), showlegend = FALSE) \%>\% 
  add_lines() \%>\%
  add_markers(color = ~pop)

# a number of `add_*()` functions are special cases of the scatter trace
plot_ly(economics, x = ~date) \%>\% 
  add_ribbons(ymin = ~pce - 1e3, ymax = ~pce + 1e3)

# use `group_by()` (or `group2NA()`) to apply visual mapping
# once per group (e.g. one line per group)
txhousing \%>\% 
  group_by(city) \%>\% 
  plot_ly(x = ~date, y = ~median) \%>\%
  add_lines(color = I("black"))

\dontrun{
# use `add_sf()` or `add_polygons()` to create geo-spatial maps
# http://blog.cpsievert.me/2018/03/30/visualizing-geo-spatial-data-with-sf-and-plotly/
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  plot_ly() \%>\% add_sf(data = nc)
}

# univariate summary statistics
plot_ly(mtcars, x = ~factor(vs), y = ~mpg) \%>\% 
  add_boxplot()
plot_ly(mtcars, x = ~factor(vs), y = ~mpg) \%>\% 
  add_trace(type = "violin")
  
# `add_histogram()` does binning for you...
mtcars \%>\%
  plot_ly(x = ~factor(vs)) \%>\%
  add_histogram()
  
# ...but you can 'pre-compute' bar heights in R
mtcars \%>\%
  dplyr::count(vs) \%>\%
  plot_ly(x = ~vs, y = ~n) \%>\%
  add_bars()

# the 2d analogy of add_histogram() is add_histogram2d()/add_histogram2dcontour()
library(MASS)
(p <- plot_ly(geyser, x = ~waiting, y = ~duration))
add_histogram2d(p)
add_histogram2dcontour(p)

# the 2d analogy of add_bars() is add_heatmap()/add_contour()
# (i.e., bin counts must be pre-specified)
den <- kde2d(geyser$waiting, geyser$duration)
p <- plot_ly(x = den$x, y = den$y, z = den$z)
add_heatmap(p)
add_contour(p)

# `add_table()` makes it easy to map a data frame to the table trace type
plot_ly(economics) \%>\% 
  add_table()

# pie charts!
ds <- data.frame(labels = c("A", "B", "C"), values = c(10, 40, 60))
plot_ly(ds, labels = ~labels, values = ~values) \%>\%
  add_pie() \%>\%
  layout(title = "Basic Pie Chart using Plotly")
  
data(wind)
plot_ly(wind, r = ~r, t = ~t) \%>\% 
  add_area(color = ~nms) \%>\%
  layout(radialaxis = list(ticksuffix = "\%"), orientation = 270)

# ------------------------------------------------------------
# 3D chart types
# ------------------------------------------------------------
plot_ly(z = ~volcano) \%>\% 
  add_surface()
plot_ly(x = c(0, 0, 1), y = c(0, 1, 0), z = c(0, 0, 0)) \%>\% 
  add_mesh()
}

}
\references{
\url{https://plotly-r.com/overview.html}

\url{https://plotly.com/r/}

\url{https://plotly.com/r/reference/}
}
\seealso{
\code{\link[=plot_ly]{plot_ly()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{offline}
\alias{offline}
\title{Plotly Offline}
\usage{
offline(p, height, width, out_dir, open_browser)
}
\arguments{
\item{p}{a plotly object}

\item{height}{A valid CSS unit. (like "100\\%", "600px", "auto") or a number,
which will be coerced to a string and have "px" appended.}

\item{width}{A valid CSS unit. (like "100\\%", "600px", "auto") or a number,
which will be coerced to a string and have "px" appended.}

\item{out_dir}{a directory to place the visualization.
If \code{NULL}, a temporary directory is used when the offline object is printed.}

\item{open_browser}{open the visualization after creating it?}
}
\value{
a plotly object of class "offline"
}
\description{
Deprecated in version 2.0 (offline plots are now the default)
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly_IMAGE.R
\name{plotly_IMAGE}
\alias{plotly_IMAGE}
\title{Create a static image}
\usage{
plotly_IMAGE(
  x,
  width = 1000,
  height = 500,
  format = "png",
  scale = 1,
  out_file,
  ...
)
}
\arguments{
\item{x}{either a plotly object or a list.}

\item{width}{Image width in pixels}

\item{height}{Image height in pixels}

\item{format}{The desired image format 'png', 'jpeg', 'svg', 'pdf', 'eps', or 'webp'}

\item{scale}{Both png and jpeg formats will be scaled beyond the specified width and height by this number.}

\item{out_file}{A filename for writing the image to a file.}

\item{...}{arguments passed onto \code{httr::RETRY}}
}
\description{
The images endpoint turns a plot (which may be given in multiple forms)
into an image of the desired format.
}
\examples{
\dontrun{
p <- plot_ly(x = 1:10)
Png <- plotly_IMAGE(p, out_file = "plotly-test-image.png")
Jpeg <- plotly_IMAGE(p, format = "jpeg", out_file = "plotly-test-image.jpeg")
Svg <- plotly_IMAGE(p, format = "svg",  out_file = "plotly-test-image.svg")
Pdf <- plotly_IMAGE(p, format = "pdf",  out_file = "plotly-test-image.pdf")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mathjax.R
\name{TeX}
\alias{TeX}
\title{Render TeX in a plotly graph using MathJax}
\usage{
TeX(x)
}
\arguments{
\item{x}{a character vector}
}
\description{
This function makes it slightly easier to render TeX in a plotly graph --
it ensures that MathJax is included with the final result and also
ensures the provided string is surrounded with \code{$} (this is what plotly.js
uses to declare a string as TeX).
}
\examples{

plot_ly(x = c(1, 2, 3, 4), y = c(1, 4, 9, 16)) \%>\%
  layout(title = TeX("\\\\text{Some mathjax: }\\\\alpha+\\\\beta x")) \%>\%
  config(mathjax = "cdn")
}
\seealso{
\link{config}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{knit_print.api_grid}
\alias{knit_print.api_grid}
\title{Embed a plotly grid as an iframe in a knitr doc}
\usage{
knit_print.api_grid(x, options, ...)
}
\arguments{
\item{x}{a plotly figure object}

\item{options}{knitr options.}

\item{...}{placeholder.}
}
\description{
Embed a plotly grid as an iframe in a knitr doc
}
\references{
https://github.com/yihui/knitr/blob/master/vignettes/knit_print.Rmd
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mic}
\alias{mic}
\title{Mic data}
\format{
A data frame with three variables: \code{r}, \code{t},
\code{nms}.
}
\usage{
mic
}
\description{
Description TBD.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{knit_print.api_plot}
\alias{knit_print.api_plot}
\title{Embed a plotly figure as an iframe in a knitr doc}
\usage{
knit_print.api_plot(x, options, ...)
}
\arguments{
\item{x}{a plotly figure object}

\item{options}{knitr options.}

\item{...}{placeholder.}
}
\description{
Embed a plotly figure as an iframe in a knitr doc
}
\references{
https://github.com/yihui/knitr/blob/master/vignettes/knit_print.Rmd
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R, R/pipe.R
\docType{import}
\name{mutate}
\alias{mutate}
\alias{mutate_}
\alias{transmute}
\alias{transmute_}
\alias{select}
\alias{select_}
\alias{rename}
\alias{rename_}
\alias{group_by}
\alias{group_by_}
\alias{groups}
\alias{ungroup}
\alias{summarise}
\alias{summarise_}
\alias{do}
\alias{do_}
\alias{arrange}
\alias{arrange_}
\alias{distinct}
\alias{distinct_}
\alias{slice}
\alias{slice_}
\alias{filter}
\alias{filter_}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{dplyr}{\code{\link[dplyr]{arrange}}, \code{\link[dplyr:se-deprecated]{arrange_}}, \code{\link[dplyr]{distinct}}, \code{\link[dplyr:se-deprecated]{distinct_}}, \code{\link[dplyr]{do}}, \code{\link[dplyr:se-deprecated]{do_}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr:se-deprecated]{filter_}}, \code{\link[dplyr]{group_by}}, \code{\link[dplyr:se-deprecated]{group_by_}}, \code{\link[dplyr:group_data]{groups}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr:se-deprecated]{mutate_}}, \code{\link[dplyr]{rename}}, \code{\link[dplyr:se-deprecated]{rename_}}, \code{\link[dplyr]{select}}, \code{\link[dplyr:se-deprecated]{select_}}, \code{\link[dplyr]{slice}}, \code{\link[dplyr:se-deprecated]{slice_}}, \code{\link[dplyr]{summarise}}, \code{\link[dplyr:se-deprecated]{summarise_}}, \code{\link[dplyr:mutate]{transmute}}, \code{\link[dplyr:se-deprecated]{transmute_}}, \code{\link[dplyr:group_by]{ungroup}}}

  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{knit_print.api_grid_local}
\alias{knit_print.api_grid_local}
\title{Embed a plotly grid as an iframe in a knitr doc}
\usage{
knit_print.api_grid_local(x, options, ...)
}
\arguments{
\item{x}{a plotly figure object}

\item{options}{knitr options.}

\item{...}{placeholder.}
}
\description{
Embed a plotly grid as an iframe in a knitr doc
}
\references{
https://github.com/yihui/knitr/blob/master/vignettes/knit_print.Rmd
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly.R
\name{remove_typedarray_polyfill}
\alias{remove_typedarray_polyfill}
\title{Remove TypedArray polyfill}
\usage{
remove_typedarray_polyfill(p)
}
\arguments{
\item{p}{a plotly object}
}
\description{
By default, plotly.js' TypedArray polyfill is included as a dependency, so
printing "just works" in any context. Many users won't need this polyfill,
so this function may be used to remove it and thus reduce the size of the page.
}
\details{
The polyfill seems to be only relevant for those rendering plots
via phantomjs and RStudio on some Windows platforms.
}
\examples{

\dontrun{
p1 <- plot_ly()
p2 <- remove_typedarray_polyfill(p1)
t1 <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(p1, t1)
file.info(t1)$size
htmlwidgets::saveWidget(p2, t1)
file.info(t1)$size
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.api}
\alias{print.api}
\title{Print method for a 'generic' API response}
\usage{
\method{print}{api}(x, ...)
}
\arguments{
\item{x}{a list.}

\item{...}{additional arguments (currently ignored)}
}
\description{
Print method for a 'generic' API response
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{wind}
\alias{wind}
\title{Wind data}
\format{
A data frame with three variables: \code{r}, \code{t},
\code{nms}.
}
\usage{
wind
}
\description{
Description TBD.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kaleido.R
\name{save_image}
\alias{save_image}
\alias{kaleido}
\title{Save plot as a static image}
\usage{
save_image(p, file, ..., width = NULL, height = NULL, scale = NULL)

kaleido(...)
}
\arguments{
\item{p}{a plot object.}

\item{file}{a file path with a suitable file extension (png, jpg, jpeg,
webp, svg, or pdf).}

\item{...}{not currently used.}

\item{width, height}{The width/height of the exported image in layout
pixels. If \code{scale} is 1, this will also be the width/height of the exported
image in physical pixels.}

\item{scale}{The scale factor to use when exporting
the figure. A scale factor larger than 1.0 will increase the image
resolution with respect to the figure's layout pixel dimensions. Whereas as
scale factor of less than 1.0 will decrease the image resolution.}
}
\value{
For \code{save_image()}, the generated \code{file}. For \code{kaleido()}, an environment that contains:
\itemize{
\item \code{transform()}: a function to convert plots objects into static images. This function has the same signature (i.e., arguments) as \code{save_image()}
\item \code{shutdown()}: a function for shutting down any currently running subprocesses
that were launched via \code{transform()}
\item \code{scope}: a reference to the underlying \code{kaleido.scopes.plotly.PlotlyScope}
python object. Modify this object to customize the underlying Chromium
subprocess and/or configure other details such as URL to plotly.js, MathJax, etc.
}
}
\description{
Static image exporting via \href{https://github.com/plotly/Kaleido/}{the kaleido python package}. \code{kaleido()} imports
kaleido into a \pkg{reticulate}d Python session and returns a \verb{$transform()}
method for converting R plots into static images. \code{save_image()} provides a convenience wrapper around \code{kaleido()$transform()}.
}
\section{Installation}{


\code{kaleido()} requires \href{https://github.com/plotly/Kaleido/}{the kaleido python package} to be usable via the \pkg{reticulate} package. Here is a recommended way to do the installation:\preformatted{install.packages('reticulate')
reticulate::install_miniconda()
reticulate::conda_install('r-reticulate', 'python-kaleido')
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')
}
}

\examples{

\dontrun{
  # Save a single image
  p <- plot_ly(x = 1:10)
  tmp <- tempfile(fileext = ".png")
  save_image(p, tmp)
  file.show(tmp)

  # Efficiently save multiple images
  scope <- kaleido()
  for (i in 1:5) {
    scope$transform(p, tmp)
  }
  # Remove and garbage collect to remove 
  # R/Python objects and shutdown subprocesses
  rm(scope); gc()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toRGB.R
\name{showRGB}
\alias{showRGB}
\title{View colors already formatted by toRGB()}
\usage{
showRGB(x, ...)
}
\arguments{
\item{x}{character string specifying color(s).}

\item{...}{arguments passed along to \code{scales::show_col}.}
}
\description{
Useful for viewing colors after they've been converted to plotly.js'
color format -- "rgba(255, 255, 255, 1)"
}
\examples{

showRGB(toRGB(colors()), labels = FALSE)
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{colorbar}
\alias{colorbar}
\title{Modify the colorbar}
\usage{
colorbar(p, ..., limits = NULL, which = 1)
}
\arguments{
\item{p}{a plotly object}

\item{...}{arguments are documented here
\url{https://plotly.com/r/reference/#scatter-marker-colorbar}.}

\item{limits}{numeric vector of length 2. Set the extent of the colorbar scale.}

\item{which}{colorbar to modify? Should only be relevant for subplots with
multiple colorbars.}
}
\description{
Modify the colorbar
}
\examples{

p <- plot_ly(mtcars, x = ~wt, y = ~mpg, color = ~cyl)

# pass any colorbar attribute -- 
# https://plotly.com/r/reference/#scatter-marker-colorbar
colorbar(p, len = 0.5)

# Expand the limits of the colorbar
colorbar(p, limits = c(0, 20))
# values outside the colorbar limits are considered "missing"
colorbar(p, limits = c(5, 6))

# also works on colorbars generated via a z value
corr <- cor(diamonds[vapply(diamonds, is.numeric, logical(1))])
plot_ly(x = rownames(corr), y = colnames(corr), z = corr) \%>\%
 add_heatmap() \%>\%
 colorbar(limits = c(-1, 1))
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{as.widget}
\alias{as.widget}
\title{Convert a plotly object to an htmlwidget object}
\usage{
as.widget(x, ...)
}
\arguments{
\item{x}{a plotly object.}

\item{...}{other options passed onto \code{htmlwidgets::createWidget}}
}
\description{
This function was deprecated in  4.0.0,
as plotly objects are now htmlwidget objects,
so there is no need to convert them.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly.R
\name{as_widget}
\alias{as_widget}
\title{Convert a list to a plotly htmlwidget object}
\usage{
as_widget(x, ...)
}
\arguments{
\item{x}{a plotly object.}

\item{...}{other options passed onto \code{htmlwidgets::createWidget}}
}
\description{
Convert a list to a plotly htmlwidget object
}
\examples{

trace <- list(x = 1, y = 1)
obj <- list(data = list(trace), layout = list(title = "my plot"))
as_widget(obj)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{event_unregister}
\alias{event_unregister}
\title{Un-register a shiny input value}
\usage{
event_unregister(p, event = NULL)
}
\arguments{
\item{p}{a plotly object.}

\item{event}{The type of plotly event. All supported events are listed in the
function signature above (i.e., the usage section).}
}
\description{
Un-register a shiny input value
}
\seealso{
\link{event_data}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly.R
\name{plot_ly}
\alias{plot_ly}
\title{Initiate a plotly visualization}
\usage{
plot_ly(
  data = data.frame(),
  ...,
  type = NULL,
  name,
  color,
  colors = NULL,
  alpha = NULL,
  stroke,
  strokes = NULL,
  alpha_stroke = 1,
  size,
  sizes = c(10, 100),
  span,
  spans = c(1, 20),
  symbol,
  symbols = NULL,
  linetype,
  linetypes = NULL,
  split,
  frame,
  width = NULL,
  height = NULL,
  source = "A"
)
}
\arguments{
\item{data}{A data frame (optional) or \link[crosstalk:SharedData]{crosstalk::SharedData} object.}

\item{...}{Arguments (i.e., attributes) passed along to the trace \code{type}.
See \code{\link[=schema]{schema()}} for a list of acceptable attributes for a given trace \code{type}
(by going to \code{traces} -> \code{type} -> \code{attributes}). Note that attributes
provided at this level may override other arguments
(e.g. \code{plot_ly(x = 1:10, y = 1:10, color = I("red"), marker = list(color = "blue"))}).}

\item{type}{A character string specifying the trace type (e.g. \code{"scatter"}, \code{"bar"}, \code{"box"}, etc).
If specified, it \emph{always} creates a trace, otherwise}

\item{name}{Values mapped to the trace's name attribute. Since a trace can
only have one name, this argument acts very much like \code{split} in that it
creates one trace for every unique value.}

\item{color}{Values mapped to relevant 'fill-color' attribute(s)
(e.g. \href{https://plotly.com/r/reference/#scatter-fillcolor}{fillcolor},
\href{https://plotly.com/r/reference/#scatter-marker-color}{marker.color},
\href{https://plotly.com/r/reference/#scatter-textfont-color}{textfont.color}, etc.).
The mapping from data values to color codes may be controlled using
\code{colors} and \code{alpha}, or avoided altogether via \code{\link[=I]{I()}} (e.g., \code{color = I("red")}).
Any color understood by \code{\link[grDevices:col2rgb]{grDevices::col2rgb()}} may be used in this way.}

\item{colors}{Either a colorbrewer2.org palette name (e.g. "YlOrRd" or "Blues"),
or a vector of colors to interpolate in hexadecimal "#RRGGBB" format,
or a color interpolation function like \code{colorRamp()}.}

\item{alpha}{A number between 0 and 1 specifying the alpha channel applied to \code{color}.
Defaults to 0.5 when mapping to \href{https://plotly.com/r/reference/#scatter-fillcolor}{fillcolor} and 1 otherwise.}

\item{stroke}{Similar to \code{color}, but values are mapped to relevant 'stroke-color' attribute(s)
(e.g., \href{https://plotly.com/r/reference/#scatter-marker-line-color}{marker.line.color}
and \href{https://plotly.com/r/reference/#scatter-line-color}{line.color}
for filled polygons). If not specified, \code{stroke} inherits from \code{color}.}

\item{strokes}{Similar to \code{colors}, but controls the \code{stroke} mapping.}

\item{alpha_stroke}{Similar to \code{alpha}, but applied to \code{stroke}.}

\item{size}{(Numeric) values mapped to relevant 'fill-size' attribute(s)
(e.g., \href{https://plotly.com/r/reference/#scatter-marker-size}{marker.size},
\href{https://plotly.com/r/reference/#scatter-textfont-size}{textfont.size},
and \href{https://plotly.com/r/reference/#scatter-error_x-width}{error_x.width}).
The mapping from data values to symbols may be controlled using
\code{sizes}, or avoided altogether via \code{\link[=I]{I()}} (e.g., \code{size = I(30)}).}

\item{sizes}{A numeric vector of length 2 used to scale \code{size} to pixels.}

\item{span}{(Numeric) values mapped to relevant 'stroke-size' attribute(s)
(e.g.,
\href{https://plotly.com/r/reference/#scatter-marker-line-width}{marker.line.width},
\href{https://plotly.com/r/reference/#scatter-line-width}{line.width} for filled polygons,
and \href{https://plotly.com/r/reference/#scatter-error_x-thickness}{error_x.thickness})
The mapping from data values to symbols may be controlled using
\code{spans}, or avoided altogether via \code{\link[=I]{I()}} (e.g., \code{span = I(30)}).}

\item{spans}{A numeric vector of length 2 used to scale \code{span} to pixels.}

\item{symbol}{(Discrete) values mapped to \href{https://plotly.com/r/reference/#scatter-marker-symbol}{marker.symbol}.
The mapping from data values to symbols may be controlled using
\code{symbols}, or avoided altogether via \code{\link[=I]{I()}} (e.g., \code{symbol = I("pentagon")}).
Any \link{pch} value or \href{https://plotly.com/r/reference/#scatter-marker-symbol}{symbol name} may be used in this way.}

\item{symbols}{A character vector of \link{pch} values or \href{https://plotly.com/r/reference/#scatter-marker-symbol}{symbol names}.}

\item{linetype}{(Discrete) values mapped to \href{https://plotly.com/r/reference/#scatter-line-dash}{line.dash}.
The mapping from data values to symbols may be controlled using
\code{linetypes}, or avoided altogether via \code{\link[=I]{I()}} (e.g., \code{linetype = I("dash")}).
Any \code{lty} (see \link{par}) value or \href{https://plotly.com/r/reference/#scatter-line-dash}{dash name} may be used in this way.}

\item{linetypes}{A character vector of \code{lty} values or \href{https://plotly.com/r/reference/#scatter-line-dash}{dash names}}

\item{split}{(Discrete) values used to create multiple traces (one trace per value).}

\item{frame}{(Discrete) values used to create animation frames.}

\item{width}{Width in pixels (optional, defaults to automatic sizing).}

\item{height}{Height in pixels (optional, defaults to automatic sizing).}

\item{source}{a character string of length 1. Match the value of this string
with the source argument in \code{\link[=event_data]{event_data()}} to retrieve the
event data corresponding to a specific plot (shiny apps can have multiple plots).}
}
\description{
This function maps R objects to \href{https://plotly.com/javascript/}{plotly.js},
an (MIT licensed) web-based interactive charting library. It provides
abstractions for doing common things (e.g. mapping data values to
fill colors (via \code{color}) or creating \link{animation}s (via \code{frame})) and sets
some different defaults to make the interface feel more 'R-like'
(i.e., closer to \code{\link[=plot]{plot()}} and \code{\link[ggplot2:qplot]{ggplot2::qplot()}}).
}
\details{
Unless \code{type} is specified, this function just initiates a plotly
object with 'global' attributes that are passed onto downstream uses of
\code{\link[=add_trace]{add_trace()}} (or similar). A \link{formula} must always be used when
referencing column name(s) in \code{data} (e.g. \code{plot_ly(mtcars, x = ~wt)}).
Formulas are optional when supplying values directly, but they do
help inform default axis/scale titles
(e.g., \code{plot_ly(x = mtcars$wt)} vs \code{plot_ly(x = ~mtcars$wt)})
}
\examples{
\dontrun{

# plot_ly() tries to create a sensible plot based on the information you 
# give it. If you don't provide a trace type, plot_ly() will infer one.
plot_ly(economics, x = ~pop)
plot_ly(economics, x = ~date, y = ~pop)
# plot_ly() doesn't require data frame(s), which allows one to take 
# advantage of trace type(s) designed specifically for numeric matrices
plot_ly(z = ~volcano)
plot_ly(z = ~volcano, type = "surface")

# plotly has a functional interface: every plotly function takes a plotly
# object as it's first input argument and returns a modified plotly object
add_lines(plot_ly(economics, x = ~date, y = ~unemploy/pop))

# To make code more readable, plotly imports the pipe operator from magrittr
economics \%>\% plot_ly(x = ~date, y = ~unemploy/pop) \%>\% add_lines()

# Attributes defined via plot_ly() set 'global' attributes that 
# are carried onto subsequent traces, but those may be over-written
plot_ly(economics, x = ~date, color = I("black")) \%>\%
 add_lines(y = ~uempmed) \%>\%
 add_lines(y = ~psavert, color = I("red"))

# Attributes are documented in the figure reference -> https://plotly.com/r/reference
# You might notice plot_ly() has named arguments that aren't in this figure
# reference. These arguments make it easier to map abstract data values to
# visual attributes.
p <- plot_ly(palmerpenguins::penguins, x = ~bill_length_mm, y = ~body_mass_g)
add_markers(p, color = ~bill_depth_mm, size = ~bill_depth_mm)
add_markers(p, color = ~species)
add_markers(p, color = ~species, colors = "Set1")
add_markers(p, symbol = ~species)
add_paths(p, linetype = ~species)

}

}
\references{
\url{https://plotly-r.com/overview.html}
}
\seealso{
\itemize{
\item For initializing a plotly-geo object: \code{\link[=plot_geo]{plot_geo()}}
\item For initializing a plotly-mapbox object: \code{\link[=plot_mapbox]{plot_mapbox()}}
\item For translating a ggplot2 object to a plotly object: \code{\link[=ggplotly]{ggplotly()}}
\item For modifying any plotly object: \code{\link[=layout]{layout()}}, \code{\link[=add_trace]{add_trace()}}, \code{\link[=style]{style()}}
\item For linked brushing: \code{\link[=highlight]{highlight()}}
\item For arranging multiple plots: \code{\link[=subplot]{subplot()}}, \code{\link[crosstalk:bscols]{crosstalk::bscols()}}
\item For inspecting plotly objects: \code{\link[=plotly_json]{plotly_json()}}
\item For quick, accurate, and searchable plotly.js reference: \code{\link[=schema]{schema()}}
}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{hide_colorbar}
\alias{hide_colorbar}
\title{Hide color bar(s)}
\usage{
hide_colorbar(p)
}
\arguments{
\item{p}{a plotly object.}
}
\description{
Hide color bar(s)
}
\examples{

p <- plot_ly(mtcars, x = ~wt, y = ~cyl, color = ~cyl)
hide_colorbar(p)
  
}
\seealso{
\code{\link[=hide_legend]{hide_legend()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{get_figure}
\alias{get_figure}
\title{Request a figure object}
\usage{
get_figure(username, id)
}
\arguments{
\item{username}{corresponding username for the figure.}

\item{id}{of the Plotly figure.}
}
\description{
Deprecated: see \code{\link[=api_download_plot]{api_download_plot()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dev.R
\name{schema}
\alias{schema}
\title{Acquire (and optionally display) plotly's plot schema}
\usage{
schema(jsonedit = interactive(), ...)
}
\arguments{
\item{jsonedit}{use \code{listviewer::jsonedit} to view the JSON?}

\item{...}{other options passed onto \code{listviewer::jsonedit}}
}
\description{
The schema contains valid attributes names, their value type,
default values (if any), and min/max values (if applicable).
}
\examples{
s <- schema()

# retrieve acceptable `layout.mapbox.style` values
if (!is.na(Sys.getenv('MAPBOX_TOKEN', NA))) {
  styles <- s$layout$layoutAttributes$mapbox$style$values
  subplot(
    plot_mapbox() \%>\% layout(mapbox = list(style = styles[3])),
    plot_mapbox() \%>\% layout(mapbox = list(style = styles[5]))
  )
}



}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proxy.R
\name{plotlyProxy}
\alias{plotlyProxy}
\alias{plotlyProxyInvoke}
\title{Modify a plotly object inside a shiny app}
\usage{
plotlyProxy(
  outputId,
  session = shiny::getDefaultReactiveDomain(),
  deferUntilFlush = TRUE
)

plotlyProxyInvoke(p, method, ...)
}
\arguments{
\item{outputId}{single-element character vector indicating the output ID
map to modify (if invoked from a Shiny module, the namespace will be added
automatically)}

\item{session}{the Shiny session object to which the map belongs; usually the
default value will suffice.}

\item{deferUntilFlush}{indicates whether actions performed against this
instance should be carried out right away, or whether they should be held
until after the next time all of the outputs are updated.}

\item{p}{a plotly proxy object (created with \code{plotlyProxy})}

\item{method}{a plotlyjs method to invoke. For a list of options,
visit \url{https://plotly.com/javascript/plotlyjs-function-reference/}}

\item{...}{unnamed arguments passed onto the plotly.js method}
}
\description{
Modify a plotly object inside a shiny app
}
\examples{


if (require("shiny") && interactive()) {
  plotly_example("shiny", "proxy_relayout")
  plotly_example("shiny", "proxy_mapbox")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add_annotations}
\alias{add_annotations}
\title{Add an annotation(s) to a plot}
\usage{
add_annotations(p, text = NULL, ..., data = NULL, inherit = TRUE)
}
\arguments{
\item{p}{a plotly object}

\item{text}{annotation text (required).}

\item{...}{these arguments are documented at
\url{https://github.com/plotly/plotly.js/blob/master/src/components/annotations/attributes.js}}

\item{data}{a data frame.}

\item{inherit}{inherit attributes from \code{\link[=plot_ly]{plot_ly()}}?}
}
\description{
Add an annotation(s) to a plot
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orca.R
\name{orca}
\alias{orca}
\alias{orca_serve}
\title{Static image exporting via orca}
\usage{
orca(
  p,
  file = "plot.png",
  format = tools::file_ext(file),
  scale = NULL,
  width = NULL,
  height = NULL,
  mathjax = FALSE,
  parallel_limit = NULL,
  verbose = FALSE,
  debug = FALSE,
  safe = FALSE,
  more_args = NULL,
  ...
)

orca_serve(
  port = 5151,
  mathjax = FALSE,
  safe = FALSE,
  request_limit = NULL,
  keep_alive = TRUE,
  window_max_number = NULL,
  quiet = FALSE,
  debug = FALSE,
  more_args = NULL,
  ...
)
}
\arguments{
\item{p}{a plotly object.}

\item{file}{output filename.}

\item{format}{the output format (png, jpeg, webp, svg, pdf, eps).}

\item{scale}{Sets the image scale. Applies to all output images.}

\item{width}{Sets the image width. If not set, defaults to \code{layout.width} value.
Applies to all output images.}

\item{height}{Sets the image height. If not set, defaults to \code{layout.height} value.
Applies to all output images.}

\item{mathjax}{whether or not to include MathJax (required to render \link{TeX}).
If \code{TRUE}, the PLOTLY_MATHJAX_PATH environment variable must be set and point
to the location of MathJax (this variable is also used to render \link{TeX} in
interactive graphs, see \link{config}).}

\item{parallel_limit}{Sets the limit of parallel tasks run.}

\item{verbose}{Turn on verbose logging on stdout.}

\item{debug}{Starts app in debug mode and turn on verbose logs on stdout.}

\item{safe}{Turns on safe mode: where figures likely to make browser window
hang during image generating are skipped.}

\item{more_args}{additional arguments to pass along to system command. This is useful
for specifying display and/or electron options, such as \code{--enable-webgl} or \code{--disable-gpu}.}

\item{...}{for \code{orca()}, additional arguments passed along to \code{processx::run}. For
\code{orca_serve()}, additional arguments passed along to \code{processx::process}.}

\item{port}{Sets the server's port number.}

\item{request_limit}{Sets a request limit that makes orca exit when reached.}

\item{keep_alive}{Turn on keep alive mode where orca will (try to) relaunch server if process unexpectedly exits.}

\item{window_max_number}{Sets maximum number of browser windows the server can keep open at a given time.}

\item{quiet}{Suppress all logging info.}
}
\description{
Superseded by \code{\link[=kaleido]{kaleido()}}.
}
\section{Methods}{


The \code{orca_serve()} function returns an object with two methods:

\describe{
\item{\code{export(p, file = "plot.png", format = tools::file_ext(file), scale = NULL, width = NULL, height = NULL)}}{
Export a static image of a plotly graph. Arguments found here are the same as those found in \code{orca()}
}
\item{\code{close()}}{Close down the orca server and kill the underlying node process.}
}
}

\section{Fields}{


The \code{orca_serve()} function returns an object with two fields:

\describe{
\item{\code{port}}{The port number that the server is listening to.}
\item{\code{process}}{An R6 class for controlling and querying the underlying node process.}
}
}

\examples{

\dontrun{
# NOTE: in a headless environment, you may need to set `more_args="--enable-webgl"`
# to export webgl correctly
p <- plot_ly(z = ~volcano) \%>\% add_surface()
orca(p, "surface-plot.svg")

#' # launch the server
server <- orca_serve()

# export as many graphs as you'd like
server$export(qplot(1:10), "test1.pdf")
server$export(plot_ly(x = 1:10, y = 1:10), "test2.pdf")

# the underlying process is exposed as a field, so you
# have full control over the external process
server$process$is_alive()

# convenience method for closing down the server
server$close()

# remove the exported files from disk
unlink("test1.pdf")
unlink("test2.pdf")
}

}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{res_mn}
\alias{res_mn}
\title{Minnesotan Indian Reservation Lands}
\format{
An sf data frame with 13 features and 5 fields
}
\usage{
res_mn
}
\description{
Minnesotan Indian Reservation Lands
}
\references{
\url{https://www.dot.state.mn.us/maps/gdma/gis-data.html}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly_data.R
\name{plotly_data}
\alias{plotly_data}
\alias{groups.plotly}
\alias{ungroup.plotly}
\alias{group_by.plotly}
\alias{mutate.plotly}
\alias{do.plotly}
\alias{summarise.plotly}
\alias{arrange.plotly}
\alias{select.plotly}
\alias{filter.plotly}
\alias{distinct.plotly}
\alias{slice.plotly}
\alias{rename.plotly}
\alias{transmute.plotly}
\alias{group_by_.plotly}
\alias{mutate_.plotly}
\alias{do_.plotly}
\alias{summarise_.plotly}
\alias{arrange_.plotly}
\alias{select_.plotly}
\alias{filter_.plotly}
\alias{distinct_.plotly}
\alias{slice_.plotly}
\alias{rename_.plotly}
\alias{transmute_.plotly}
\title{Obtain data associated with a plotly graph}
\usage{
plotly_data(p, id = p$x$cur_data)

\method{groups}{plotly}(x)

\method{ungroup}{plotly}(x, ...)

\method{group_by}{plotly}(.data, ...)

\method{mutate}{plotly}(.data, ...)

\method{do}{plotly}(.data, ...)

\method{summarise}{plotly}(.data, ...)

\method{arrange}{plotly}(.data, ...)

\method{select}{plotly}(.data, ...)

\method{filter}{plotly}(.data, ...)

\method{distinct}{plotly}(.data, ...)

\method{slice}{plotly}(.data, ...)

\method{rename}{plotly}(.data, ...)

\method{transmute}{plotly}(.data, ...)

\method{group_by_}{plotly}(.data, ...)

\method{mutate_}{plotly}(.data, ...)

\method{do_}{plotly}(.data, ...)

\method{summarise_}{plotly}(.data, ...)

\method{arrange_}{plotly}(.data, ...)

\method{select_}{plotly}(.data, ...)

\method{filter_}{plotly}(.data, ...)

\method{distinct_}{plotly}(.data, ...)

\method{slice_}{plotly}(.data, ...)

\method{rename_}{plotly}(.data, ...)

\method{transmute_}{plotly}(.data, ...)
}
\arguments{
\item{p}{a plotly visualization.}

\item{id}{a character string or number referencing an "attribute layer".}

\item{x}{a plotly visualization.}

\item{...}{arguments passed onto the relevant method.}

\item{.data}{a plotly visualization.}
}
\description{
\code{plotly_data()} returns data associated with
a plotly visualization (if there are multiple data frames, by default,
it returns the most recent one).
}
\examples{

# use group_by() to define groups of visual markings
p <- txhousing \%>\%
  group_by(city) \%>\%
  plot_ly(x = ~date, y = ~sales)
p

# plotly objects preserve data groupings 
groups(p)
plotly_data(p)

# dplyr verbs operate on plotly objects as if they were data frames
p <- economics \%>\%
  plot_ly(x = ~date, y = ~unemploy / pop) \%>\%
  add_lines() \%>\%
  mutate(rate = unemploy / pop) \%>\% 
  filter(rate == max(rate))
plotly_data(p)
add_markers(p)
layout(p, annotations = list(x = ~date, y = ~rate, text = "peak"))

# use group_by() + do() + subplot() for trellis displays 
d <- group_by(mpg, drv)
plots <- do(d, p = plot_ly(., x = ~cty, name = ~drv))
subplot(plots[["p"]], nrows = 3, shareX = TRUE)

# arrange displays by their mean
means <- summarise(d, mn = mean(cty, na.rm = TRUE))
means \%>\%
  dplyr::left_join(plots) \%>\%
  arrange(mn) \%>\%
  subplot(nrows = NROW(.), shareX = TRUE)
  
# more dplyr verbs applied to plotly objects
p <- mtcars \%>\%
  plot_ly(x = ~wt, y = ~mpg, name = "scatter trace") \%>\%
  add_markers()
p \%>\% slice(1) \%>\% plotly_data()
p \%>\% slice(1) \%>\% add_markers(name = "first observation")
p \%>\% filter(cyl == 4) \%>\% plotly_data()
p \%>\% filter(cyl == 4) \%>\% add_markers(name = "four cylinders")


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly.R
\name{plot_mapbox}
\alias{plot_mapbox}
\title{Initiate a plotly-mapbox object}
\usage{
plot_mapbox(data = data.frame(), ...)
}
\arguments{
\item{data}{A data frame (optional).}

\item{...}{arguments passed along to \code{\link[=plot_ly]{plot_ly()}}. They should be
valid scattermapbox attributes - \url{https://plotly.com/r/reference/#scattermapbox}.
Note that x/y can also be used in place of lat/lon.}
}
\description{
Use this function instead of \code{\link[=plot_ly]{plot_ly()}} to initialize
a plotly-mapbox object. This enforces the entire plot so use
the scattermapbox trace type, and enables higher level geometries
like \code{\link[=add_polygons]{add_polygons()}} to work
}
\examples{
\dontrun{

plot_mapbox(res_mn)
plot_mapbox(res_mn, color = ~INDRESNAME)

map_data("world", "canada") \%>\%
  group_by(group) \%>\%
  plot_mapbox(x = ~long, y = ~lat) \%>\%
  add_polygons() \%>\%
  layout(
    mapbox = list(
      center = list(lat = ~median(lat), lon = ~median(long))
    )
  )
}

}
\seealso{
\code{\link[=plot_ly]{plot_ly()}}, \code{\link[=plot_geo]{plot_geo()}}, \code{\link[=ggplotly]{ggplotly()}}
}
\author{
Carson Sievert
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{plotly}
\alias{plotly}
\title{Main interface to plotly}
\usage{
plotly(username, key)
}
\arguments{
\item{username}{plotly username}

\item{key}{plotly API key}
}
\description{
This function is now deprecated. It used to provide a way to store plotly
account credentials, but that can now be done with environment variables.
For more details and examples, see \url{https://plotly.com/r/getting-started/}.

If you're here looking for an intro/overview of the package, see the
\href{https://github.com/plotly/plotly.R/#getting-started}{readme}
}
\seealso{
\code{\link[=ggplotly]{ggplotly()}}, \code{\link[=plot_ly]{plot_ly()}}, \code{\link[=signup]{signup()}}
}
\keyword{internal}
