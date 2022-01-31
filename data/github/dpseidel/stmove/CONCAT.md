
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis Build
Status](https://travis-ci.com/dpseidel/stmove.svg?token=ZVrezsGfh5uSAe6FpgAU&branch=master)](https://travis-ci.com/dpseidel/stmove)
[![Codecov test
coverage](https://codecov.io/gh/dpseidel/stmove/branch/master/graph/badge.svg?token=A1gUYaWSSY)](https://codecov.io/gh/dpseidel/stmove)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
status](https://www.r-pkg.org/badges/version/stmove)](https://cran.r-project.org/package=stmove)

# stmove

The goal of stmove is to make more accessible and transparent standard
spatio-temporal approaches to interpreting movement data before getting
into more challenging aspects of deconstructing movement trajectories.
Analyses in this package expect “clean, regular, movement data time
series” which consists of a sequence of points `(x, y, time)` where
`time = 0, 1, 2, 3,..., T` and all missing points have been interpolated
and filled in.

For a detailed review of stmove’s motivation and functionality, please
see [our preprint available on
Biorxiv](http://biorxiv.org/cgi/content/short/758987v1).

## Installation

stmove is in active development and not yet available on CRAN. To
download the current version of the package use the following code:

``` r
# install.packages("remotes")
remotes::install_github("dpseidel/stmove")
```

If you encounter bugs or have features you would like to see
incorporated in future version of stmove, please open an issue.

## Usage

The primary function of this package is `build_report` which, given a
clean regularized trajectory, will deliver a .Rmd and .pdf report
including the results of the following computations.

### Basic path distributions:

1.  Generate step size time series \(S={(t,s(t)) | over [0,T]}\) and
    plot step-size histogram
2.  Generate turning angle time series \(A={(t,a(t)) | over [0,T]}\) and
    plot turning angle distribution

### Basic path statistics:

1.  If data is sub hourly, plot running means, standard deviations, and
    auto correlations for s(t) and a(t), plus cross-correlation
    s(t)-a(t) (basic stats collection: BSC) using a 3 hour (or 6 hour if
    data is only hourly) “sliding time window” (STW)
2.  Generate and plot 12-hourly BSC using a 12 hour “jumping time
    window” (JTW)
3.  Generate and plot 14/15 day BSC following a new-full-new moon
    sequence using a JTW
4.  Generate and plot seasonal BSC for how ever many seasons are
    available using a JTW

### Basic path visualizations.

1.  Plot the trajectory over space
2.  Generate a wavelet plot that allows us to visualize possible
    periodic components in s(t) and a(t), auto and cross correlation
    coefficients.

### Basic space constructions:

1.  Construct 25%, 50% and 95% home range isopleths (optionally by
    season) using two methods:
    1.  `k-LoCoH` hull sets
    2.  autocorrelated utilization distribution analysis implemented
        with `ctmm::akde`

Once these are done, one can then pursue various kinds of analysis that
address questions of interest (e.g., GLM models of location and
landscape, HMM modeling, step section analysis), but this package
strives to set a standard for what is minimally needed before embarking
on such analyses.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis Build Status](https://travis-ci.com/dpseidel/stmove.svg?token=ZVrezsGfh5uSAe6FpgAU&branch=master)](https://travis-ci.com/dpseidel/stmove)
[![Codecov test coverage](https://codecov.io/gh/dpseidel/stmove/branch/master/graph/badge.svg?token=A1gUYaWSSY)](https://codecov.io/gh/dpseidel/stmove)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN status](https://www.r-pkg.org/badges/version/stmove)](https://cran.r-project.org/package=stmove)


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# stmove

The goal of stmove is to make more accessible and transparent standard spatio-temporal
approaches to interpreting movement data before getting into more challenging 
aspects of deconstructing movement trajectories.  Analyses in this package expect
“clean, regular, movement data time series” which consists of a sequence of
points `(x, y, time)` where `time = 0, 1, 2, 3,..., T` and all missing points 
have been interpolated and filled in.


For a detailed review of stmove's motivation and functionality, please see [our preprint available on Biorxiv](http://biorxiv.org/cgi/content/short/758987v1).


## Installation

stmove is in active development and not yet available on CRAN. To download the current version of the package use the following code:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("dpseidel/stmove")
```

If you encounter bugs or have features you would like to see incorporated in future version of stmove, please open an issue. 

## Usage

The primary function of this package is `build_report` which, given a clean
regularized trajectory, will deliver a .Rmd and .pdf report including the 
results of the following computations. 

### Basic path distributions:
1. Generate step size time series $S={(t,s(t)) | over [0,T]}$ and plot step-size histogram
2. Generate turning angle time series $A={(t,a(t)) | over [0,T]}$ and plot turning angle distribution

### Basic path statistics:
1. If data is sub hourly, plot running means, standard deviations, and auto 
correlations for s(t) and a(t), plus cross-correlation s(t)-a(t) (basic stats 
collection: BSC) using a 3 hour (or 6 hour if data is only hourly) "sliding time window” (STW)
2. Generate and plot 12-hourly BSC  using a 12 hour “jumping time window” (JTW)
3. Generate and plot 14/15 day BSC following a new-full-new moon sequence using a JTW
4. Generate and plot seasonal BSC for how ever many seasons are available using a JTW

### Basic path visualizations.
1.  Plot the trajectory over space
2.  Generate a wavelet plot that allows us to visualize possible periodic 
components in s(t) and a(t), auto and cross correlation coefficients.

### Basic space constructions:
1.  Construct 25%, 50%  and 95% home range isopleths (optionally by season) using two methods: 
    1. `k-LoCoH` hull sets
    2. autocorrelated utilization distribution analysis implemented with `ctmm::akde`

Once these are done, one can then pursue various kinds of analysis that address 
questions of interest (e.g., GLM models of location and landscape, HMM modeling,
step section analysis), but this package strives to set a standard for what is 
minimally needed before embarking on such analyses.
---
title: "Population Trajectory Analyses"
subtitle: "created with pkg stmove" 
date: "`r format(Sys.time(), '%d %B, %Y')`"
output: pdf_document
params:
  df: NA
  proj4: NA
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
```

```{r pop, messages = F,  fig.height = 4}
dist <- dist_map(params$df, params$proj4)

plot_timeline(params$df)
```
---
title: "Movement Report for `r params$df$id[1]`"
subtitle: "created with pkg stmove" 
date: "`r format(Sys.time(), '%d %B, %Y')`"
output: pdf_document
geometry: margin=2cm
params:
  df: NA
  stats: NA
  seas: NA
  construct: NA
  proj4: NA
  wavelet: NA
---

```{r setup, include = F}
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

library(stmove)
library(patchwork)

construct <- (!is.null(params$construct))
klocoh <- "klocoh" %in% tolower(params$construct)
akde <- "akde" %in% tolower(params$construct)
rolling <- ("rolling" %in% params$stats)
diurnal <- ("diurnal" %in% params$stats)
lunar <- ("lunar" %in% params$stats)
seasonal <- ("seasonal" %in% params$stats)
wavelet <- (!is.null(params$wavelet))
```

# Visualize the trajectory

```{r}
ggplot(params$df, aes(x, y)) + geom_point() + labs(x = "Easting", y = "Northing") +
  theme_minimal()
```

# Step Length & Turning Angle Distributions:

```{r, fig.height = 4}
ss <- ss_dist(x = params$df)
```

```{r, fig.height = 4}
ta <- ta_dist(x = params$df)
```

`r if (any(c(rolling, diurnal, lunar, seasonal))) '# Basic Statistical Summary'`

`r if (rolling) '## Fine Scale Sliding Window'`

```{r eval = rolling, messages = F, fig.height = 4}
roll <<- rolling_stats(params$df)

ggplot(roll, aes(date, y = mean_dist)) +
  geom_path(na.rm = T) +
  facet_wrap(~ parse_factor(paste(
    lubridate::month(date, label = TRUE),
    lubridate::year(date)
  ), ordered = T),
  scales = "free_x"
  ) +
  labs(
    title = "Mean Step Size on a Rolling interval",
    y = "Mean Step Size (m)"
  ) + theme_minimal() +
  scale_x_datetime("Day of the month", date_labels = "%d")


ggplot(roll, aes(date, y = mean_ang)) +
  geom_path(na.rm = T) +
  facet_wrap(~ parse_factor(paste(
    lubridate::month(date, label = TRUE),
    lubridate::year(date)
  ), ordered = T),
  scales = "free_x"
  ) +
  labs(
    title = "Mean Turning Angle on a Rolling interval",
    y = "Mean Angle (radians)"
  ) + theme_minimal() +
  scale_x_datetime("Day of the month", date_labels = "%d")


suppressMessages(ggplot(roll, aes(date)) +
  geom_smooth(aes(y = acf_ang, color = "Relative Angle"), na.rm = T) +
  geom_smooth(aes(y = acf_dist, color = "Step Size"), na.rm = T) +
  geom_smooth(aes(y = ccf, color = "Cross Correlation"), na.rm = T) +
  ylab("Correlation") +
  ggtitle("Lag 1 Auto and Cross Correlations in Turning Angle & Step Size",
    subtitle = "Rolling Window"
  )) + theme_minimal()
```

`r if (diurnal) '## Diurnal Cycle'`

```{r eval = diurnal, message = FALSE}
diurnal_stats <<- interval_stats(params$df)

# We need to make some decisions about how to plot these things...
plots <- purrr::map(
  names(diurnal_stats)[c(3, 6, 5, 8, 9)],
  function(x) {
    y <- sym(x)
    p <- ggplot(
      diurnal_stats,
      aes(interval_start, y = !!y, color = TOD)
    ) + geom_point(na.rm = T, size = .3) +
      geom_smooth(na.rm = T) +
      xlab(NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 25, vjust = .5))

    if (x == "acf_ang") {
      return(p + scale_color_discrete("Time of Day") +
        theme(legend.position = "bottom"))
    } else {
      return(p + guides(color = "none"))
    }
  }
)

suppressMessages({
  plots[[1]] + plots[[2]]
} - {
  plots[[3]] + plots[[4]] + plots[[5]]
} + plot_layout(nrow = 2) +
  plot_annotation(title = "Diurnal Interval Statistics"))
```

`r if (lunar) '## Lunar Cycle'`

```{r eval = lunar}
lunar_stats <<- interval_stats(params$df, type = "lunar")

# We need to make some decisions about how to plot these things...
plots <- purrr::map(
  names(lunar_stats)[c(3, 6, 5, 8, 9)],
  function(x) {
    y <- sym(x)
    p <- ggplot(
      lunar_stats,
      aes(x = as.numeric(interval_start), y = !!y, color = phase)
    ) +
      geom_step(aes(group = 1)) +
      # scale_x_continuous(breaks = 1:nlevels(lunar_stats$interval_start),
      #                   labels = levels(lunar_stats$interval_start)) +
      xlab(NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 25, vjust = .5))

    if (x == "acf_ang") {
      return(p + scale_color_discrete("Lunar Phase") +
        theme(legend.position = "bottom"))
    } else {
      return(p + guides(color = "none"))
    }
  }
)

{
  plots[[1]] + plots[[2]]
} - {
  plots[[3]] + plots[[4]] + plots[[5]]
} + plot_layout(nrow = 2) +
  plot_annotation(title = "Lunar Interval Statistics")
```

`r if (seasonal) '## Seasonal Cycle'`

```{r eval = seasonal}
seasonal_stats <<- interval_stats(params$df, type = "seasonal", seas = params$seas)

p1 <- seasonal_stats %>%
  # tidyr::gather("stat", "estimate", -1) %>%
  ggplot(., aes(x = as.numeric(interval_start), y = mean_dist)) +
  geom_col() +
  geom_errorbar(aes(
    ymin = mean_dist,
    ymax = mean_dist + sd_dist
  )) + labs(x = "", y = "Average Step Size (+ 1 SD)")
p2 <- seasonal_stats %>%
  ggplot(., aes(x = as.numeric(interval_start), y = mean_ang)) +
  geom_col() +
  labs(x = "Season", y = "Mean Turning Angle")

p3 <- seasonal_stats %>%
  select(-mean_dist, -sd_dist, -mean_ang, -sd_ang) %>%
  tidyr::gather("stat", "estimate", -1) %>%
  ggplot(., aes(
    x = as.factor(as.numeric(interval_start)), y = estimate,
    group = stat, color = stat
  )) +
  geom_point() +
  geom_line() +
  lims(y = c(-1, 1)) +
  labs(x = "")

p1 + p2 + p3 +
  plot_annotation(
    title = "Interval Statistics Across Seasons",
    theme = theme(plot.title = element_text(hjust = .5))
  ) &
  theme_minimal()
```

`r if (wavelet) '# Visualizing Periodicity'`

```{r eval = wavelet}
wave <<- wavelet(params$df, stats = params$wavelet, useRaster = TRUE)
```

`r if (construct) '# Space Use Constructions'`

`r if (klocoh) '## A k-LoCoh hull set'`

```{r, eval = klocoh}
locoh <<- construct(params$df, type = "klocoh", proj4 = params$proj4)
```

`r if (akde) '## An autocorrelated kernel density estimation'`
```{r, eval = akde}
AKDE <<- construct(params$df, type = "akde", proj4 = params$proj4)
```

\newpage
# Suggested Next Steps

This report provides a first-pass look at your trajectories and patterns you see here
may be further investigated in any number of ways. Below we suggest some common 
next questions and steps and associated packages that could help you begin to address them. 
Also see suggestions from Seidel et al. 2019, Appendix C of Seidel et al. 2018,
or tutorials from the 2018 Movement Ecology in R workshop hosted [here](https://www.danaseidel.com/MovEco-R-Workshop/).

- Habitat Selection Analysis
    - the `sf` package offers tools for ready manipulation of spatial vector data. 
    - the `raster` and `velox` packages offer tools for raster data manipulation. 
    - the `lme4` package is useful for fitting generalized linear models with fixed or random effects
- Hidden Markov Models can be fit or Behavioral Change Point Analysis done with tools from:
    - `momentuHMM`
    - `moveHmm`
    - `bcpa`
- Dynamic Interaction between individuals may be explored use
    - `wildlifeDI`
    
Also consider the numerous analyses made available by the extentive "adehabitat" packages:
`adehabitatLT`, `adehabitatHR`, `adehabitatMA`, `adehabitatHS` and their indepth vignettes. 
---
title: "Getting started with stmove"
author: "Eric Dougherty"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

The `stmove` package builds upon several other analytics and visualization packages to offer the user a straight-forward user interface and a more digestible report of a set of high-level analyses frequently undertaken to describe animal movement data. To take full advantage of the `build_report` function, the input data needs to be a regularized trajectory, meaning that geospatial fixes need to occur at regular intervals throughout the period of interest. Often, movement data collected in the field will not satisfy this requirement. This means that some pre-processing may be necessary, and the `stmove` package offers a few utility functions to aid in the preparation of the movement tracks.

This tutorial will cover these utility functions, taking a raw trajectory from its original form to a regularized path. Using the result of this pre-processing, we will cover some of the functionality offered by the `build_report` function and walk through the final report that it creates. The tutorial is conducted in RStudio to demonstrate the user interface elements of the package, but it can be run without this feature in native R.

```{r, warning=FALSE, message=FALSE}
library(stmove)
library(sf)
library(dplyr)
library(ggplot2)
```

## Example Path

For this tutorial, we will use an elephant trajectory from Southern Africa. The individual (AG195) was tracked between July 2009 and January 2012, resulting in nearly 130,000 relocation points. For purpose of demonstration, we will focus on the year-long period between January 1, 2010 and December 31, 2010. The original data collection procedure involved alternative fix intervals of 1 minute and 19 minutes. This means that every other point will be separated by approximately 20 minutes. We will eliminate the intevening short fix intervals for the high-level analyses here by making use of the `regularize` function. For convenience, the stmove package bundles the raw 2010 trajectory for AG195 as `AG195_raw`.

```{r}
head(AG195_raw)
```

The spatial position of the animal in this case has been collected in the form of latlong coordinates, but it will be much more useful in a projection like UTM where the units are more easily interpretable. Since we know that this data was collected in UTM Zone 33S, we can easily change its projection during preprocessing.

The easiest way to perform such reprojections is with the `sf` package, which works well with the tidyverse tools, like those used above for manipulating the data frame. We'll want to use the `st_as_sf` function, passing the first two columns as our `coords` argument and 4326 (the EPSG value for the latlong projection in WGS84) as the `crs` argument. 

```{r}
AG195_sf <- AG195_raw %>%
  st_as_sf(coords = 1:2, crs = 4236) %>%
  st_transform(crs = 32733)
```

If we have a general idea of what the track should look like in space, we can plot the result of that function, which created a new geometry column, in order to verify that the points are correctly projected as latlong coordinates.

```{r}
ggplot() +
  geom_sf(data = AG195_sf) +
  coord_sf(datum = st_crs(32733))
```

That appears to be correct, and we can print out the head of our new dataframe to see that the geometry column is in UTM units rather than decimal degrees:

```{r}
head(AG195_sf)
```

Perfect! The next thing we'll notice about the dataset here is the fact that the positional fixes are not exactly regular in the normal sense. Even the much-reduced view offered by the `head` function demonstrates that fixes were actually captured at two intervals, 1 minute and 19 minutes. This means that every other point is separated by 20 minutes, and that will be the regular interval we aim for as we move forward with this analysis. There may be some very interesting information in those 1 minute intervals, but for the sake of the broad-scale analyses we will be conducting with `stmove`, it is very reasonable to eliminate these intervening points.

One way to go about regularizing this path is the `regularize` utility function provided by `stmove`. In order to use this function, we'll need to revert from the simple features geometry column that we've just created to x and y coordinates. The data frame that we'll need to feed into the `regularize` function consists of those two correctly projected coordinate columns, a datetime column (renamed "date" to appease the inner workings of the function), and the id column (not particualrly relevant here, but it would be important to maintain if we were preprocessing/analyzing multiple individuals simultaneously). We'll also eliminate the geometry column by setting it to NULL. This final data frame will have 51832 rows and 4 columns.

```{r}
AG195_reg <- AG195_sf %>%
  mutate(
    x = st_coordinates(AG195_sf)[, 1],
    y = st_coordinates(AG195_sf)[, 2]
  ) %>%
  mutate(date = datetime) %>%
  select(x, y, date, id)
st_geometry(AG195_reg) <- NULL
AG195_reg <- stmove::regularize(AG195_reg, dt = 20, units = "min")
```

The results is a much more informative, temporally-regularized path. In addition to the coordinates, timestamps, and id columns we fed into the function, we have several additional columns that were created when the trajectory was converted to an ltraj object in the `regularize` function. One last step we'll need to undertake is to determine if that regularization process (or the data collection process itself) resulted in any gaps in the data. The regularization function will fill in NA values for the x and y coordinates at timestamps that we expect to have a reading, but where none is available. Using a table, we can ascertain the frequency of such gaps in our path:

```{r}
table(is.na(AG195_reg$x))
```

It appears we have only 31 missing points! That is astoundingly few for a trajectory that has 26209 points! We'll still need to fill in those gaps to make sure our downstream analyses work. One way to do that is with the `kalman` utility function provided by `stmove`.

```{r}
AG195_final <- stmove::kalman(AG195_reg)
table(is.na(AG195$x))
```

This function may take some time to run, but its doing some extraordinary things under the hood. Using all of the known points and the distribution of step lengths between them, the algorithm interpolates any missing points that we have. This process should be used with some discernment, however. In this case, the proportion of missing points to known points is so low that we should have no problem filling in those gaps, but other trajectories may be riddled with missing points once they have been temporally regularized. In such cases, the interpolation process may lose some of its integrity, as the model has less information with which to predict those unknown values. We would not recommend relying upon the Kalman smoothing approach implemented in `stmove` if the proportion of missing points is around 5%.

We can take a look at our final data frame (which should look a lot like it used to, but perhaps slightly sparser). We should be wary of any peculiar outliers that may have appeared, as these may be examples of the Kalman smoother misbehaving. Note that the plot below may be somewhat distorted compared to the plot we made earlier. The reason is that we are no longer plotting projected geographic points here, we are simply plotting x and y values. Even so, we should be able to pick out any outliers, which is the main reason for building this intermediate plot following the Kalman smoothing process.

```{r}
ggplot() +
  geom_point(aes(x = x, y = y), data = AG195_final)
```

These steps have provided a suitable path for running the `build_report` function offered by `stmove`. The result of this is a PDF file with several interesting broad-scale analyses of the path in question. We simply designate the data frame on which the analyses will be conducted, the name or filepath of the output file, and the projection (in our case "+init=epsg:32733"), and the report will be created with rolling and diurnal statistics, a klocoh home range, and a wavelet diagram. Each of these provides potentially interesting insights that we will discuss below.

The build process might take a little while, and we would expect that time to scale with the length of the trajectory, but even with the intensive wavelet and home range calculations, the report is done within a matter of minutes! Here, we have set the `wavelet` argument to TRUE, but the default is NULL, so if you want to see that, you'll need to add that argument on top of the required `df`, `path`, and `proj4` arguments. The default arguments for `stats` are both "rolling"" and "diurnal"", but we will add "lunar" to the set, but we won't add "seasonal" because we would also need to input a vector of the start dates of each season. For our purposes here, the other statistics will be sufficient. The default argument for `construct` is "klocoh" (for the k local convex hull method), which will also be perfect for our analyses here. One could alternatively (or additionally) select "akde" for the adaptive kernel density estimation approach to constructing the home range.

```{r, eval=FALSE}
stmove::build_report(AG195_final,
  path = "~/Desktop/",
  proj4 = "+init=epsg:32733",
  stats = c("rolling", "diurnal", "lunar"),
  construct = c("klocoh"),
  wavelet = "dist"
)
```

 

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stmove.R
\docType{package}
\name{stmove}
\alias{stmove}
\alias{stmove-package}
\title{The \code{stmove} package}
\description{
See the README on \href{https://github.com/dpseidel/stmove#readme}{GitHub}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{create_telemetry}
\alias{create_telemetry}
\title{Build telemetry objects using projected data}
\usage{
create_telemetry(df, proj4)
}
\arguments{
\item{df}{a dataframe containing columns x, y, date representing relocations in space and time.}

\item{proj4}{a character string indicating the proj.4 definition of the
coordinate reference system defining the relocations}
}
\description{
A helper function to create telemetry objects required by \code{ctmm}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_report.r
\name{build_report}
\alias{build_report}
\title{Build a report with basic spatio-temporal movement computations}
\usage{
build_report(df, path = here::here(), proj4, stats = c("rolling",
  "diurnal"), construct = c("klocoh"), seas = NULL, wavelet = NULL)
}
\arguments{
\item{df}{a dataframe containing columns "x", "y", and "date"}

\item{path}{absolute file path for output .pdf, defaults to project working directory.
see also \link[here]{here}.}

\item{proj4}{a character string indicating the proj.4 definition of the
coordinate reference system defining the relocations}

\item{stats}{a character vector of stats to calculate, options include: "rolling",
"diurnal", "lunar", and "seasonal"}

\item{construct}{a character vector indicating which spacetime construction methods to use.
options include "klocoh" and "akde"}

\item{seas}{a character vector including the start date of each season interval in the format
\%Y-\%m-\%d, e.g. "2015-01-01".  Required if \code{type == "seasonal"}.}

\item{wavelet}{a character vector indicating which variables to calculate wavelet analysis on,
options in "dist", "rel.angle", "acf_dist", "acf_ang", and "ccf".}
}
\description{
Build a report with basic spatio-temporal movement computations
}
\examples{
\donttest{
build_report(AG195, proj4 = "+proj=utm +zone=33 +south +datum=WGS84 +units=m +no_defs")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AG268}
\alias{AG268}
\title{AG268: elephant trajectory}
\format{A data frame with 53940 rows and 10 variables:
\describe{
\item{x}{x coordinate}
\item{y}{y coordinate}
\item{date}{timestamp of relocation}
\item{id}{animal id}
\item{real}{logical indicating empirical or interpolated position}
}}
\usage{
AG268
}
\description{
A dataframe containing the clean regularized relocation points for
African elephant (\emph{Loxodonta africana}) AG195. All points are from
2010 and were collected by Getz et al. x and y coordinates are projected in
WGS 84 / UTM zone 33S (espg:32733). The fix rate for this individual is 1 fix
every 15 minutes. Columns are as follows:
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructions.r
\name{construct}
\alias{construct}
\title{Home Range Construction by Two Methods}
\usage{
construct(df, type = c("klocoh", "akde"), proj4)
}
\arguments{
\item{df}{a dataframe containing columns x, y, date representing relocations in space and time.}

\item{type}{a character string indicating the type/s of constructions to build, "klocoh" or "akde"}

\item{proj4}{a character string indicating the proj.4 definition of the
coordinate reference system defining the relocations}
}
\description{
Home Range Construction by Two Methods
}
\examples{
\donttest{
klocoh <- construct(AG195, type = "klocoh", proj = "+proj=utm +zone=33 +south")
akde <- construct(AG195, type = "akde", proj = "+proj=utm +zone=33 +south")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.r
\name{interval_stats}
\alias{interval_stats}
\title{Interval Summary Statistics}
\usage{
interval_stats(df, type = "diurnal", seas = NULL)
}
\arguments{
\item{df}{a dataframe containing columns "x", "y", and "date"}

\item{type}{a character string indicating the intervals to calculate, one of
"diurnal", "lunar", or "seasonal".}

\item{seas}{a character vector including the start date of each season interval in the format
\%Y-\%m-\%d, e.g. "2015-01-01".  Required if \code{type == "seasonal"}.}
}
\description{
Interval Summary Statistics
}
\details{
Note: if your trajectory includes dates before the first date in your \code{seas} vector, they
will automatically be considered a separate season.
}
\examples{
# diurnal is default
interval_stats(AG195)
interval_stats(AG195, type = "lunar")
# for seasonal, include y-m-d formatted `seas` vector
interval_stats(AG195, type = "seasonal", seas = c("2010-03-01", "2010-06-01", "2010-9-01"))
}
\seealso{
\link[lunar]{lunar.phase}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{kalman}
\alias{kalman}
\title{Kalman Smoothing with Structural Time Series}
\usage{
kalman(df, warn = TRUE)
}
\arguments{
\item{df}{a dataframe containing columns x, y, date representing relocations in space and time.}

\item{warn}{logical, should warnings be issued if interpolation is above 5\%?}
}
\value{
a dataframe with additional binary column \code{real} flagging those points whose positions were interpolated.
}
\description{
Interpolate missing values in a trajectory.
}
\details{
The replacement points are generated using a structural time series model fitted by maximum likelihood.
}
\seealso{
\link[imputeTS]{na.kalman}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{regularize}
\alias{regularize}
\title{Regularize a near-regular trajectory}
\usage{
regularize(df, dt, units = "min", ref = NULL)
}
\arguments{
\item{df}{a dataframe containing columns x, y, date representing relocations in space and time.}

\item{dt}{the expected time lag between relocations}

\item{units}{a character string indicating the time units for dt and tol}

\item{ref}{a datetime from which to start the fixes. Default \code{ref=NULL} will
calculate ref value by rounding the first timestamp in \code{df}}
}
\description{
Regularize a near-regular trajectory
}
\seealso{
\link[adehabitatLT]{setNA} \link[adehabitatLT]{sett0} \link[adehabitatLT]{subsample}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AG195}
\alias{AG195}
\title{AG195: elephant trajectory}
\format{A data frame with 53940 rows and 10 variables:
\describe{
\item{x}{x coordinate}
\item{y}{y coordinate}
\item{date}{timestamp of relocation}
\item{id}{animal id}
\item{real}{logical indicating empirical or interpolated position}
}}
\usage{
AG195
}
\description{
A dataframe containing the clean regularized relocation points for
African elephant (\emph{Loxodonta africana}) AG195. All points are from
2010 and were collected by Getz et al. x and y coordinates are projected in
WGS 84 / UTM zone 33S (espg:32733). The fix rate for this individual is 1 fix
every 20 minutes. Columns are as follows:
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.r
\name{ss_dist}
\alias{ss_dist}
\alias{ta_dist}
\title{Basic Distributions}
\usage{
ss_dist(x, plot = T)

ta_dist(x, plot = T)
}
\arguments{
\item{x}{a dataframe with columns: x, y, date, and id (optional)}

\item{plot}{a logical indicating whether or not to return a histogram of the distribution}
}
\value{
a numeric vector and (optionally) a plot of the turning angle or step length distribution
}
\description{
The first step to any movement analysis is to extract the step length and turning angle
distributions of trajectory. These functions will calculate step length and relative
turning angles for a given trajectory and plot associated histograms.
}
\examples{
\donttest{
ss <- ss_dist(AG195)
ta <- ta_dist(AG195)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population.R
\name{dist_map}
\alias{dist_map}
\title{Spatial Temporal Distribution Plot}
\usage{
dist_map(df, proj4, labels = TRUE)
}
\arguments{
\item{df}{a dataframe containing columns "x", "y", "date", and "id"}

\item{proj4}{a character string indicating the proj.4 definition of the
coordinate reference system defining the relocations}

\item{labels}{logical, indicating whether or not to label points with IDs, implemented using
\link[ggrepel]{geom_text_repel}}
}
\description{
Spatial Temporal Distribution Plot
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.r
\name{rolling_stats}
\alias{rolling_stats}
\title{Rolling Summary Statistics}
\usage{
rolling_stats(df, n_roll = NULL)
}
\arguments{
\item{df}{a dataframe containing columns "x", "y", and "date}

\item{n_roll}{numeric, number of fixes (e.g. width of the window) over which to
calculate rolling statistics. If NULL, the default is to roll over 3 hours for fix-rates
of 1 hour or less, and over 6 hours for fix rates occuring at intervals greater than 1 hour.
If your fix rate is every 4 hours or greater and n_roll is not set, the function
will warn you to set \code{n_roll}.}
}
\description{
Rolling Summary Statistics
}
\details{
some detailed discussion of how auto and cross correlation are calculated and
na handling
}
\examples{
roll <- rolling_stats(AG195)
}
\seealso{
\link[RcppRoll]{roll_mean} \link[RcppRoll]{roll_sdr} \link[TTR]{runCor}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population.R
\name{plot_timeline}
\alias{plot_timeline}
\title{Population Timeline Plot}
\usage{
plot_timeline(df)
}
\arguments{
\item{df}{a dataframe containing columns "x", "y", "date", and "id"}
}
\description{
Population Timeline Plot
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wavelet.R
\name{wavelet}
\alias{wavelet}
\title{Wavelet analyses}
\usage{
wavelet(df, stats = c("dist", "rel.angle", "acf_dist", "acf_ang", "ccf"),
  plot = T, ...)
}
\arguments{
\item{df}{a dataframe with columns: x, y, date}

\item{stats}{a character vector indicating which variables to calculate wavelet analysis on,
options in "dist", "rel.angle", "acf_dist", "acf_ang", and "ccf". By default, analysis is run on all 5.}

\item{plot}{a logical indicating whether or not to return a histogram of the distribution}

\item{...}{Arguments passed to \link[dplR]{wavelet.plot}. Only relevant when \code{plot = T}. Useful
for changing plot defaults e.g. axis labels.}
}
\description{
Build morelet wavelet diagrams for step size, turning angle, auto and cross correlation
coefficients through time.
}
\examples{
\donttest{
dist_wave <- wavelet(AG195, stats = "dist", UseRaster = T)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_report_addin.r
\name{report_addin}
\alias{report_addin}
\title{Build a Report (RStudio Add-in)}
\usage{
report_addin()
}
\description{
\code{report_addin()} opens an \href{https://shiny.rstudio.com/articles/gadgets.html}{RStudio gadget} and
\href{http://rstudio.github.io/rstudioaddins/}{addin} that allows you to
choose which analyses to run and build a report.
Appears as "Build report" in the RStudio Addins menu.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AG195_raw}
\alias{AG195_raw}
\title{AG195_raw: raw elephant trajectory}
\format{A data frame with 51964 rows and 4 variables:
\describe{
\item{longitude}{x coordinate}
\item{latitude}{y coordinate}
\item{datetime}{timestamp of relocation}
\item{id}{animal id}
}}
\usage{
AG195_raw
}
\description{
A dataframe containing the raw irregular relocation points for
African elephant (\emph{Loxodonta africana}) AG195. All points are from
2010 and were collected by Getz et al. Coordinates reflect the WGS84 geographic
coordinate system. The fix rate for this individual alternates
between every 1 and 19 minutes.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{df_check}
\alias{df_check}
\title{Check dataframes conform to stmove style}
\usage{
df_check(df)
}
\arguments{
\item{df}{a dataframe}
}
\description{
An (internal) helper function to check dataframes and provide helpful errors
}
\keyword{internal}
