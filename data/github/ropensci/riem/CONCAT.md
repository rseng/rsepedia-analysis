riem
====
<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/riem)](http://cran.r-project.org/package=riem)
[![codecov](https://app.codecov.io/gh/ropensci/riem/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/riem)
[![](https://badges.ropensci.org/39_status.svg)](https://github.com/ropensci/software-review/issues/39)
  [![R-CMD-check](https://github.com/ropensci/riem/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/riem/actions)
  <!-- badges: end -->


This package allows to get weather data from ASOS stations (airports)
via the great website of the [Iowa Environment
Mesonet](https://mesonet.agron.iastate.edu/request/download.phtml?network=IN__ASOS).

Installation
============

Install the package with:

``` r
install.packages("riem")
```

Or install the development version using
[devtools](https://github.com/r-lib/devtools) with:

``` r
library("devtools")
install_github("ropensci/riem")
```

# Documentation

Please refer to the [`pkgdown` website](https://docs.ropensci.org/riem/) to read docs, in particular the [reference](https://docs.ropensci.org/riem/reference/index.html) and the [vignettes](https://docs.ropensci.org/riem/articles/index.html).

# Use cases in the wild

Submit your use cases by opening [a new issue](https://github.com/ropensci/riem/issues/new)!

* [@maelle](https://github.com/maelle/)'s [Where to live in the US](https://masalmon.eu/2017/11/16/wheretoliveus/)

* [@Ryo-N7](https://github.com/Ryo-N7)'s [Where to live in Japan: XKCD-themed climate plots and maps!](https://ryo-n7.github.io/2017-11-22-japan-xkcd-weather-index/)

* [@sharlagelfand](https://github.com/sharlagelfand)'s [Oh, the cold places I've lived](https://sharlagelfand.netlify.app/posts/oh-the-cold-places-ive-lived/)

* In a course: [@byuidatascience](https://github.com/byuidatascience)'s [M335 Task 15: How is the weather?](https://byuistats.github.io/M335/backgrounds.html#task_15:_how_is_the_weather)

Meta
----

-   Please [report any issues or
    bugs](https://github.com/ropensci/riem/issues).
-   License: GPL
-   Get citation information for `ropenaq` in R doing
    `citation(package = 'riem')`
-   Please note that this project is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree
    to abide by its terms.

[![ropensci\_footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
# riem (development version)

# riem 0.2.0

* Switches to newer IEM metadata web services (#35, @akrherz)

# riem 0.1.1

* Eliminates a few dependencies (dplyr, lazyeval, readr) to make installation easier.

* Now the default end date for `riem_measures` is the current date as given by `Sys.Date()`.

# riem 0.1.0

* Added a `NEWS.md` file to track changes to the package.



## Test environments
* local Ubuntu install
* win-builder
* GitHub Actions

## R CMD check results

R CMD check results
0 errors | 0 warnings | 0 notes

## Release summary

Switches to newer IEM metadata web services.

## Reverse dependencies

There are no reverse dependencies.
---
title: "riem"
---

The riem package allows to get weather data from ASOS stations (airports) via the awesome website of the [Iowa Environment Mesonet](https://mesonet.agron.iastate.edu/request/download.phtml?network=IN__ASOS).


# Installation

Install the package with:

```{r eval = FALSE}
install.packages("riem")
```

Or install the development version using [devtools](https://github.com/hadley/devtools) with:

```{r, eval = FALSE}
library("devtools")
install_github("ropenscilabs/riem")
```

# Get available networks

```{r, warning = FALSE, message = FALSE}
library("riem")
riem_networks() 
```

# Get available stations for one network

```{r}
riem_stations(network = "IN__ASOS") 
```


# Get measures for one station

Possible variables are (copied from [here](https://mesonet.agron.iastate.edu/request/download.phtml), see also the [ASOS user guide](http://www.nws.noaa.gov/asos/pdfs/aum-toc.pdf))

* station: three or four character site identifier

* valid: timestamp of the observation (UTC)

* tmpf: Air Temperature in Fahrenheit, typically @ 2 meters

* dwpf: Dew Point Temperature in Fahrenheit, typically @ 2 meters

* relh: Relative Humidity in \%

* drct: Wind Direction in degrees from north

* sknt: Wind Speed in knots

* p01i: One hour precipitation for the period from the observation time to the time of the previous hourly precipitation reset. This varies slightly by site. Values are in inches. This value may or may not contain frozen precipitation melted by some device on the sensor or estimated by some other means. Unfortunately, we do not know of an authoritative database denoting which station has which sensor.

* alti: Pressure altimeter in inches

* mslp: Sea Level Pressure in millibar

* vsby: Visibility in miles

* gust: Wind Gust in knots

* skyc1: Sky Level 1 Coverage

* skyc2: Sky Level 2 Coverage

* skyc3: Sky Level 3 Coverage

* skyc4: Sky Level 4 Coverage

* skyl1: Sky Level 1 Altitude in feet

* skyl2: Sky Level 2 Altitude in feet

* skyl3: Sky Level 3 Altitude in feet

* skyl4: Sky Level 4 Altitude in feet

* presentwx: Present Weather Codes (space seperated), see e.g. [this manual](http://www.ofcm.gov/fmh-1/pdf/H-CH8.pdf) for further explanations.

* feel: Apparent Temperature (Wind Chill or Heat Index) in degF

* ice_accretion_1hr: Ice Accretion over 1 Hour in inch

* ice_accretion_3hr: Ice Accretion over 3 Hour in inch

* ice_accretion_6hr: Ice Accretion over 6 Hour in inch

* relh: Relative Humidity in % 

* metar: unprocessed reported observation in METAR format

* peak_wind_gust: Wind gust in knots from the METAR PK WND remark, this value may be different than the value found in the gust field. The gust field is derived from the standard METAR wind report.

* peak_wind_drct: The wind direction in degrees North denoted in the METAR PK WND remark.

* peak_wind_time: The timestamp of the PK WND value in the same timezone as the valid field and controlled by the tz parameter.


```{r}
measures <- riem_measures(station = "VOHY", date_start = "2000-01-01", date_end = "2016-04-22") 
head(measures)
```

For conversion of wind speed or temperature into other units, see [this package](https://github.com/geanders/weathermetrics/).

---
title: "Forecasting temperature at your vacation destination"
author: "M. Salmon"
---

For Christmas I'll travel to Marseille. What temperatures should I expect there? I could of course open a weather app, but in this vignette I want to give an example using the `riem` and  `forecast` packages. 

# Find airport for Marseille

The name of the network for France is "FR__ASOS". I already know there's only one airport near the city.

```{r, cache = FALSE, message = FALSE, warning = FALSE}
library("riem")
library("dplyr")
france_airports <- riem_stations(network = "FR__ASOS")
marseilles_airport <- filter(france_airports, grepl("MARSEILLE", name) | grepl("Marseille", name))
marseilles_airport

```

# Get time series of temperature for Marseille airport

We'll transform it to daily average, and convert Fahrenheit to Celsius thanks to the `weathermetrics` package. We impute the missing values and remove outliers via the use of `forecast::tsclean`.


```{r, cache = FALSE, fig.width=8, fig.height=6, message = FALSE, warning = FALSE}
marseille <- riem_measures(station = marseilles_airport$id,
                           date_start = "2010 01 01")

marseille <- group_by(marseille, day = as.Date(valid))
marseille <- summarize(marseille, temperature = mean(tmpf))
marseille <- mutate(marseille, temperature = weathermetrics::fahrenheit.to.celsius(temperature))
library("ggplot2")


library("forecast")
marseille_ts = ts(as.vector(tsclean(marseille$temperature)), freq=365.25, start=c(2010, 1))
autoplot(marseille_ts) +
  ylab("Daily average temperature in Marseille airport (ºC)") +
  xlab("Time (days)")

```

# Forecast for Marseille

For this we use the `forecast` package. We use the `stlm` because our time series obviously present yearly seasonality.

```{r, cache = FALSE, fig.width=8, fig.height=6}

fit <- stlm(marseille_ts)
pred <- forecast(fit, h = 7)
# plot
theme_set(theme_gray(base_size = 14))
autoplot(pred) +
  ylab("Daily average temperature in Marseille airport (ºC)") +
  xlab("Time (days)") +
  ggtitle("How cold will I be during the holidays?",
          subtitle = "Data accessed via the rOpenSci riem package and forecasted with forecast")

```

Mmh I don't see anything, but `autoplot.forecast` has an `include` parameters, so I'll only plot the last 31 values.


```{r, cache = FALSE, fig.width=8, fig.height=6}

autoplot(pred, include = 31) +
  ylab("Daily average temperature in Marseille airport (ºC)") +
  xlab("Time (days)") +
  ggtitle("How cold will I be during the holidays?",
          subtitle = "Data accessed via the rOpenSci riem package and forecasted with forecast")
```

# What if I went somewhere else?

Ok, what if I had travelled to, say, Hyderabad in India?

```{r, echo = FALSE, cache = FALSE, fig.width=8, fig.height=6}

hyderabad <- riem_measures(station = "VOHY",
                           date_start = "2010 01 01")

hyderabad <- group_by(hyderabad, day = as.Date(valid))
hyderabad <- summarize(hyderabad, temperature = mean(tmpf))
hyderabad <- mutate(hyderabad, temperature = weathermetrics::fahrenheit.to.celsius(temperature))

hyderabad_ts = ts(as.vector(tsclean(hyderabad$temperature)), freq=365.25, start=c(2010, 1))
fit <- stlm(hyderabad_ts)
pred <- forecast(fit, h = 7)
autoplot(pred, include = 31) +
  ylab("Daily average temperature in Hyderabad airport (ºC)") +
  xlab("Time (days)") +
  ggtitle("Hyderabad forecasted temperatures in December 2016",
          subtitle = "Data accessed via the rOpenSci riem package and forecasted with forecast")

```

Without surprise, we forecast I'd have enjoyed warmer weather. 

I wouldn't advise you to really use such code to forecast temperature, but I'd recommend you to use `riem` for getting weather airport data quite easily and to dig more deeply into `forecast` functionalities if you're interested in time series forecasting. And stay warm!
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures.R
\name{riem_measures}
\alias{riem_measures}
\title{Function for getting weather data from one station}
\usage{
riem_measures(
  station = "VOHY",
  date_start = "2014-01-01",
  date_end = as.character(Sys.Date())
)
}
\arguments{
\item{station}{station ID, see riem_stations()}

\item{date_start}{date of start of the desired data, e.g. "2000-01-01"}

\item{date_end}{date of end of the desired data, e.g. "2016-04-22"}
}
\value{
a data.frame (tibble tibble) with measures, the number of columns can vary from station to station,
but possible variables are
\itemize{
\item station: three or four character site identifier
\item valid: timestamp of the observation (UTC)
\item tmpf: Air Temperature in Fahrenheit, typically @ 2 meters
\item dwpf: Dew Point Temperature in Fahrenheit, typically @ 2 meters
\item relh: Relative Humidity in %
\item drct: Wind Direction in degrees from north
\item sknt: Wind Speed in knots
\item p01i: One hour precipitation for the period from the observation time to the time of the previous hourly precipitation reset. This varies slightly by site. Values are in inches. This value may or may not contain frozen precipitation melted by some device on the sensor or estimated by some other means. Unfortunately, we do not know of an authoritative database denoting which station has which sensor.
\item alti: Pressure altimeter in inches
\item mslp: Sea Level Pressure in millibar
\item vsby: Visibility in miles
\item gust: Wind Gust in knots
\item skyc1: Sky Level 1 Coverage
\item skyc2: Sky Level 2 Coverage
\item skyc3: Sky Level 3 Coverage
\item skyc4: Sky Level 4 Coverage
\item skyl1: Sky Level 1 Altitude in feet
\item skyl2: Sky Level 2 Altitude in feet
\item skyl3: Sky Level 3 Altitude in feet
\item skyl4: Sky Level 4 Altitude in feet
\item presentwx: Present Weather Codes (space seperated),
 see e.g. Chapter 8 of [this manual](https://www.ofcm.gov/publications/fmh/FMH1/FMH1.pdf) for further explanations.
\item feel: Apparent Temperature (Wind Chill or Heat Index) in degF
\item ice_accretion_1hr: Ice Accretion over 1 Hour in inch
\item ice_accretion_3hr: Ice Accretion over 3 Hour in inch
\item ice_accretion_6hr: Ice Accretion over 6 Hour in inch
\item relh: Relative Humidity in %
\item metar: unprocessed reported observation in METAR format
\item peak_wind_gust: Wind gust in knots from the METAR PK WND remark, this value may be different than the value found in the gust field. The gust field is derived from the standard METAR wind report.
\item peak_wind_drct: The wind direction in degrees North denoted in the METAR PK WND remark.
\item peak_wind_time: The timestamp of the PK WND value in the same timezone as the valid field and controlled by the tz parameter.
}
}
\description{
Function for getting weather data from one station
}
\details{
The data is queried through \url{https://mesonet.agron.iastate.edu/request/download.phtml}.
}
\examples{
\dontrun{
riem_measures(station = "VOHY", date_start = "2000-01-01", date_end = "2016-04-22")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stations.R
\name{riem_stations}
\alias{riem_stations}
\title{Function for getting stations of an ASOS network}
\usage{
riem_stations(network = NULL)
}
\arguments{
\item{network}{A single network code, see riem_networks() for finding the code corresponding to a name.}
}
\value{
a data.frame (tibble tibble) with the id, name, longitude (lon) and latitude (lat) of each station in the network.
}
\description{
Function for getting stations of an ASOS network
}
\details{
You can see a map of stations in a network at \url{https://mesonet.agron.iastate.edu/request/download.phtml}.
}
\examples{
\dontrun{
riem_stations(network = "IN__ASOS")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networks.R
\name{riem_networks}
\alias{riem_networks}
\title{Function for getting ASOS and AWOS networks}
\usage{
riem_networks()
}
\value{
a data.frame (tibble tibble) with the names and codes of available networks.
}
\description{
Function for getting ASOS and AWOS networks
}
\examples{
\dontrun{
riem_networks()
}
}
