pendulum
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/pendulum/workflows/R-check/badge.svg)](https://github.com/ropensci/pendulum/actions?query=workflow%3AR-check)

Time classes for R, w/ time mocking capability via [timefuzz][]

**Why?** - Some programming languages support [monkey patching][monkey]. However, R does not officially allow this (there's a way to do it, but it's a hack). Therefore, the only way that I can think of to "mock time" is to use a different set of time classes/objects. That is, we can't just mock time using `Sys.Date()` directly as there's no way to monkey patch the code behind that function.  Some may say it's good that R doesn't officially support monkey patching; there's definitely upsides; but a downside is that mocking things like time becomes difficult.

- `clock()` is designed for setting a specific time
- `now()` is designed for getting the current time AND is mockable
- `sys_time()` and `sys_date()` are drop in replacements for `Sys.time()` and `Sys.Date()` AND are mockable

See the docs to get started: https://docs.ropensci.org/pendulum/

## Installation


```r
remotes::install_github("ropensci/pendulum")
```


```r
library(pendulum)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pendulum/issues).
* License: MIT
* Get citation information for `pendulum` in R doing `citation(package = 'pendulum')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[timefuzz]: https://github.com/ropensci/timefuzz
[monkey]: https://en.wikipedia.org/wiki/Monkey_patch
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

* Submit an issue on the [Issues page](https://github.com/ropensci/pendulum/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/pendulum.git`
* Make sure to track progress upstream (i.e., on our version of `pendulum` at `ropensci/pendulum`) by doing `git remote add upstream https://github.com/ropensci/pendulum.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/pendulum`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
pendulum
=====

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/pendulum/workflows/R-check/badge.svg)](https://github.com/ropensci/pendulum/actions?query=workflow%3AR-check)

Time classes for R, w/ time mocking capability via [timefuzz][]

**Why?** - Some programming languages support [monkey patching][monkey]. However, R does not officially allow this (there's a way to do it, but it's a hack). Therefore, the only way that I can think of to "mock time" is to use a different set of time classes/objects. That is, we can't just mock time using `Sys.Date()` directly as there's no way to monkey patch the code behind that function.  Some may say it's good that R doesn't officially support monkey patching; there's definitely upsides; but a downside is that mocking things like time becomes difficult.

- `clock()` is designed for setting a specific time
- `now()` is designed for getting the current time AND is mockable
- `sys_time()` and `sys_date()` are drop in replacements for `Sys.time()` and `Sys.Date()` AND are mockable

See the docs to get started: https://docs.ropensci.org/pendulum/

## Installation

```{r eval=FALSE}
remotes::install_github("ropensci/pendulum")
```

```{r}
library(pendulum)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pendulum/issues).
* License: MIT
* Get citation information for `pendulum` in R doing `citation(package = 'pendulum')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[timefuzz]: https://github.com/ropensci/timefuzz
[monkey]: https://en.wikipedia.org/wiki/Monkey_patch
---
title: "pendulum"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{pendulum}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Time classes for R, w/ time mocking capability via [timefuzz][]

```{r setup}
library(pendulum)
```

## clock

```{r}
clock(2009)$time
clock(2009, 3, 13)$time
clock(2009, 3, 13, 1, 4, 53)$time
```

```{r}
x <- clock(2009, 3, 13)
x$utc
x$date
```

## now

```{r}
now()$now()
now("UTC")$now()
```

## sys_date

Mockable [Sys.Date()] drop in replacement. Uses `now()` internally so can be mocked as shown below.

```{r}
Sys.Date()
sys_date()
```

## sys_time

Mockable [Sys.time()] drop in replacement. Uses `now()` internally so can be mocked as shown below.

```{r}
Sys.time()
sys_time()
```

## use in a function

```{r}
todays_date <- function() sys_time()
todays_date()
```

now let's mock time

```{r}
library(timefuzz)
x <- time_fuzz$new()
## set to today + 5 days
x$freeze(Sys.Date() + 5)
todays_date()
# if you run it again, you get the same EXACT time
todays_date()
```

"unfreeze" and you're back to real time

```{r}
x$unfreeze()
todays_date()
Sys.sleep(1)
todays_date()
Sys.sleep(1)
todays_date()
```

[timefuzz]: https://github.com/ropensci/timefuzz
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sys_time.R
\name{sys_time}
\alias{sys_time}
\title{Mockable \code{\link[=Sys.time]{Sys.time()}} drop in replacement}
\usage{
sys_time()
}
\value{
\code{POSIXct}, matching what \code{\link[=Sys.time]{Sys.time()}} returns
}
\description{
Mockable \code{\link[=Sys.time]{Sys.time()}} drop in replacement
}
\examples{
Sys.time()
sys_time()

clock_mock()
library(timefuzz)
x <- time_fuzz$new()
x$freeze(Sys.time() + 500)
sys_time()
clock_mock(FALSE)
sys_time()
}
\seealso{
\code{\link[=sys_date]{sys_date()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pendulum-package.R
\docType{package}
\name{pendulum-package}
\alias{pendulum-package}
\title{pendulum}
\description{
Time Classes
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Now.R
\name{now}
\alias{now}
\title{now}
\usage{
now(tzone = "")
}
\arguments{
\item{tzone}{(character) a valid time zone, defaults to your local timezone}
}
\description{
now
}
\examples{
x <- now()
x
x$time
x$date
x$utc
x$unix_epoch
x$min
x$sec

# mocking
Sys.Date()
clock_mock()
library(timefuzz)
mock <- time_fuzz$new()
## set to today + 10 days
mock$freeze(Sys.Date() + 10)
z <- now()
z
z$time
z$date
z$utc
z$unix_epoch
z$min
z$sec
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_mock.R
\name{clock_mock}
\alias{clock_mock}
\title{Mock time}
\usage{
clock_mock(on = TRUE)
}
\arguments{
\item{on}{(logical) turn mocking on with \code{TRUE} or turn off with \code{FALSE}.
By default is \code{FALSE}}
}
\description{
Mock time
}
\details{
\code{timefuzz} package required for mocking behavior
}
\examples{
\dontrun{

if (interactive()) {
  # load timefuzz
  library(timefuzz)
  library(clock)

  # turn on mocking
  clock_mock()

  # do stuff
  x <- time_fuzz$new()
  ## set to today + 435 days
  x$freeze(Sys.Date() + 435)

  # get time again
  now()

  # turn off mocking
  clock_mock(FALSE)

  # get time again
  now()
}

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock.R
\name{clock}
\alias{clock}
\title{Clock class}
\usage{
clock(...)
}
\arguments{
\item{...}{date parts, must be supplied in order (i.e., if you want to give
seconds, you must give all others first): year, month, day, hour, minutes,
seconds. uses \code{\link[=ISOdate]{ISOdate()}} internally to convert to a \code{POSIXct}. optional;
if nothing supplied, the current time is used.}
}
\value{
a \link{POSIXct} object
}
\description{
Clock class
}
\details{
\strong{Methods}
These must be called with parens
\describe{
\item{\code{now(tzone = "")}}{
Get the date time right now
param \code{tzone}: a valid time zone
}
\item{\code{date()}}{
Same as $now() but as a \code{Date} class
}
\item{\code{utc()}}{
Same as $now() but with time zone set to UTC
}
}

\strong{Active Methods}
These are called without parens
\describe{
\item{\code{time}}{
Get the current time or the time you set on initialize
}
\item{\code{year}}{
Get the year component of the time
}
\item{\code{month}}{
Get the month component of the time
}
\item{\code{day}}{
Get the day component of the time
}
\item{\code{hour}}{
Get the hour component of the time
}
\item{\code{min}}{
Get the minute component of the time
}
\item{\code{sec}}{
Get the seconds component of the time
}
\item{\code{utc_offset}}{
Get the utc offset
}
\item{\code{unix_epoch}}{
Get the time in seconds since unix epoch
}
}
}
\examples{
# if no time set when initialized
## uses current time when you call functions with parens
x <- clock(2009)
x
## and uses current time when you call active functions
x$utc
x$date
x$year
x$month
 
# if time set when initialized
## ignores that time when calling functions with parens
z <- clock(2009, 3, 13)
z
z$utc
## and uses the user set time when you call active functions
z$year
z$month
z$utc_offset
z$unix_epoch

# 
clock(2009)$time
clock(2009, 3, 13)$time
clock(2009, 3, 13, 1, 4, 53)$time

now()
now("UTC")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sys_date.R
\name{sys_date}
\alias{sys_date}
\title{Mockable \code{\link[=Sys.Date]{Sys.Date()}} drop in replacement}
\usage{
sys_date()
}
\value{
\code{Date}, matching what \code{\link[=Sys.Date]{Sys.Date()}} returns
}
\description{
Mockable \code{\link[=Sys.Date]{Sys.Date()}} drop in replacement
}
\examples{
Sys.Date()
sys_date()

clock_mock()
library(timefuzz)
x <- time_fuzz$new()
x$freeze(Sys.Date() + 5)
sys_date()
clock_mock(FALSE)
sys_date()
}
\seealso{
\code{\link[=sys_time]{sys_time()}}
}
