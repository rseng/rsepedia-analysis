timefuzz
========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/timefuzz/workflows/R-check/badge.svg)](https://github.com/ropensci/timefuzz/actions?query=workflow%3AR-check)

Time travel to test time dependent code - a port of Ruby's [timecop](https://github.com/travisjeffery/timecop)

**BEWARE: VERY ALPHA**

Package API:

 - `time_fuzz`
 - `TimeStackItem`
 - `Time`

Features supported:

- freeze: freeze time to a specific point

Hope to support soon:

- travel: travel back to a specific point in time, but allow time to continue moving forward from there
- scale: scale time by a given scaling factor that will cause time to move at an accelerated pace

## Installation


```r
remotes::install_github("ropensci/timefuzz")
```


```r
library("timefuzz")
```

## freeze

### freeze with a code block


```r
library(timefuzz)
library(pendulum)
library(testthat)
```

`book_due()` is a toy function that tells us if a book is due


```r
book_due <- function(due_date = Sys.Date() + 10) {
  sys_date() > due_date
}
```

Given the due date of 2021-01-29 the book is not due


```r
expect_false(book_due()) # FALSE
```

Create a `time_fuzz` object


```r
x <- time_fuzz$new()
x
#> <time_fuzz> 
#>   date:
```

Call `freez()`, passing the date you want to freeze time to, and then a code block to run
in the context of that frozen time. Here we'll freeze time to today + 450 days


```r
x$freeze(Sys.Date() + 450, {
  expect_true(book_due())
})
#> Error: book_due() is not TRUE
#> 
#> `actual`:   FALSE
#> `expected`: TRUE
```

`book_due()` results in `TRUE` now, whereas it was `FALSE` above in real time

### freeze without a code block


```r
x <- time_fuzz$new()
## set to today + 450 days
x$freeze(Sys.Date() + 450)
```

We're in the freezed date. So any time based actions using the [pendulum][] package 
are now using your frozen time context.


```r
sys_time()
#> [1] "2021-01-19 10:29:39 PST"
```

call `$unfreeze()` to unfreeze


```r
x$unfreeze()
```

now we're back in current time


```r
sys_time()
#> [1] "2021-01-19 10:29:39 PST"
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/timefuzz/issues).
* License: MIT
* Get citation information for `timefuzz` in R doing `citation(package = 'timefuzz')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
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

* Submit an issue on the [Issues page](https://github.com/ropensci/timefuzz/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/brranching.git`
* Make sure to track progress upstream (i.e., on our version of `brranching` at `ropensci/timefuzz`) by doing `git remote add upstream https://github.com/ropensci/timefuzz.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/timefuzz`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
timefuzz
========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/timefuzz/workflows/R-check/badge.svg)](https://github.com/ropensci/timefuzz/actions?query=workflow%3AR-check)

Time travel to test time dependent code - a port of Ruby's [timecop](https://github.com/travisjeffery/timecop)

**BEWARE: VERY ALPHA**

Package API:

```{r echo=FALSE, comment=NA, results="asis"}
cat(paste(" -", paste(sprintf("`%s`", getNamespaceExports("timefuzz")), collapse = "\n - ")))
```

Features supported:

- freeze: freeze time to a specific point

Hope to support soon:

- travel: travel back to a specific point in time, but allow time to continue moving forward from there
- scale: scale time by a given scaling factor that will cause time to move at an accelerated pace

## Installation

```{r eval=FALSE}
remotes::install_github("ropensci/timefuzz")
```

```{r}
library("timefuzz")
```

## freeze

### freeze with a code block

```{r}
library(timefuzz)
library(pendulum)
library(testthat)
```

`book_due()` is a toy function that tells us if a book is due

```{r}
book_due <- function(due_date = Sys.Date() + 10) {
  sys_date() > due_date
}
```

Given the due date of `r Sys.Date() + 10` the book is not due

```{r}
expect_false(book_due()) # FALSE
```

Create a `time_fuzz` object

```{r}
x <- time_fuzz$new()
x
```

Call `freez()`, passing the date you want to freeze time to, and then a code block to run
in the context of that frozen time. Here we'll freeze time to today + 450 days

```{r}
x$freeze(Sys.Date() + 450, {
  expect_true(book_due())
})
```

`book_due()` results in `TRUE` now, whereas it was `FALSE` above in real time

### freeze without a code block

```{r}
x <- time_fuzz$new()
## set to today + 450 days
x$freeze(Sys.Date() + 450)
```

We're in the freezed date. So any time based actions using the [pendulum][] package 
are now using your frozen time context.

```{r}
sys_time()
```

call `$unfreeze()` to unfreeze

```{r}
x$unfreeze()
```

now we're back in current time

```{r}
sys_time()
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/timefuzz/issues).
* License: MIT
* Get citation information for `timefuzz` in R doing `citation(package = 'timefuzz')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time.R
\name{Time}
\alias{Time}
\title{Time}
\description{
Time
}
\examples{
\dontrun{
x <- Time$new()
x$now_without_mock_time()
x$now_with_mock_time()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time_stack_item.R
\name{TimeStackItem}
\alias{TimeStackItem}
\title{TimeStackItem}
\description{
TimeStackItem
}
\examples{
\dontrun{
x <- TimeStackItem$new(mock_type = "freeze", date = "2019-02-18")
x
x$time
x$time$time
x$time$year
x$scaling_factor_
x$travel_offset_
x$year()
x$min()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time_fuzz.R
\name{time_fuzz}
\alias{time_fuzz}
\title{time_fuzz}
\value{
an \link{time_fuzz} object
}
\description{
time_fuzz
}
\details{
\strong{Methods}
\describe{
\item{\code{freeze(...)}}{
Freeze time
}
}
}
\examples{
\dontrun{
x <- time_fuzz$new()
x
x$date

x <- time_fuzz$new()
x
x$freeze
x$freeze("2019-01-29", {
  5 + 5
})

x$scale
x$scale({
  5 + 5
})

library(timefuzz)
library(clock)
library(testthat)
# pendulum::clock_mock()
pendulum:::clock_opts$mock
(cl <- clock())
book_due <- function(due_date = "2020-08-25") {
  as.POSIXct(pendulum::clock()$date()) > as.POSIXct(due_date)
}
expect_false(book_due()) # FALSE
x <- time_fuzz$new()
x
x$freeze(Sys.Date() + 60, {
  # pendulum::clock_mock()
  # cat(pendulum:::clock_opts$mock, sep = "\n")
  # cat(as.character(pendulum::clock()$date()), sep = "\n")
  # cat(pendulum:::clock_opts$mock, sep = "\n")
  expect_true(book_due())
})

# no block passed
clock_mock()
clock()$now()
Sys.time()
x <- time_fuzz$new()
## set to today + 435 days
x$freeze(Sys.Date() + 435)
pendulum:::clock_opts$mock
clock()$now()
x$unfreeze()
pendulum:::clock_opts$mock
clock()$now()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timefuzz-package.R
\docType{package}
\name{timefuzz-package}
\alias{timefuzz-package}
\alias{timefuzz}
\title{timefuzz}
\description{
Time Travel to Test Time Dependent Code
}
\author{
Scott Chamberlain
}
\keyword{package}
