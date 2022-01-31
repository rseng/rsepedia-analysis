conditionz
==========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/conditionz/workflows/R-check/badge.svg)](https://github.com/ropensci/conditionz/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/conditionz/coverage.svg?branch=master)](https://codecov.io/github/ropensci/conditionz?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/conditionz)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/conditionz)](https://cran.r-project.org/package=conditionz)

control how many times conditions are thrown

docs: https://ropensci.github.io/conditionz/

Package API:

 - `handle_messages`
 - `handle_conditions`
 - `ConditionKeeper`
 - `handle_warnings`
 - `capture_message`
 - `capture_warning`

Use cases for `conditionz` functions:

- `ConditionKeeper` is what you want to use if you want to keep track of conditions inside a
function being applied many times, either in a for loop or lapply style fashion.
- `handle_conditions`/`handle_messages`/`handle_warnings` is what you want to use if the multiple
conditions are happening within a single function or code block
- `capture_message`/`capture_warning` are meant for capturing messages/warnings into a useable
list

## Installation

The CRAN version:


```r
install.packages("conditionz")
```

Or the development version:


```r
remotes::install_github("ropensci/conditionz")
```


```r
library("conditionz")
```

## ConditionKeeper

`ConditionKeeper` is the internal R6 class that handles keeping track of
conditions and lets us determine if conditions have been encountered,
how many times, etc.


```r
x <- ConditionKeeper$new(times = 4)
x
#> ConditionKeeper
#>  id: 98e020c2-e805-4da0-a369-d1939dc7c998
#>  times: 4
#>  messages: 0
x$get_id()
#> [1] "98e020c2-e805-4da0-a369-d1939dc7c998"
x$add("one")
x$add("two")
x
#> ConditionKeeper
#>  id: 98e020c2-e805-4da0-a369-d1939dc7c998
#>  times: 4
#>  messages: 2
#>   one  two
x$thrown_already("one")
#> [1] TRUE
x$thrown_already("bears")
#> [1] FALSE
x$not_thrown_yet("bears")
#> [1] TRUE

x$add("two")
x$add("two")
x$add("two")
x$thrown_times("two")
#> [1] 4
x$thrown_enough("two")
#> [1] TRUE
x$thrown_enough("one")
#> [1] FALSE
```

## basic usage

A simple function that throws messages


```r
squared <- function(x) {
  stopifnot(is.numeric(x))
  y <- x^2
  if (y > 20) message("woops, > than 20! check your numbers")
  return(y)
}
foo <- function(x) {
  vapply(x, function(z) squared(z), numeric(1))
}
bar <- function(x, times = 1) {
  y <- ConditionKeeper$new(times = times)
  on.exit(y$purge())
  vapply(x, function(z) y$handle_conditions(squared(z)), numeric(1))
}
```

Running the function normally throws many messages


```r
foo(1:10)
#> woops, > than 20! check your numbers
#> woops, > than 20! check your numbers
#> woops, > than 20! check your numbers
#> woops, > than 20! check your numbers
#> woops, > than 20! check your numbers
#> woops, > than 20! check your numbers
#>  [1]   1   4   9  16  25  36  49  64  81 100
```

Using in `ConditionKeeper` allows you to control how many messages
are thrown


```r
bar(x = 1:10)
#> woops, > than 20! check your numbers
#>  [1]   1   4   9  16  25  36  49  64  81 100
```


```r
bar(1:10, times = 3)
#> woops, > than 20! check your numbers
#> 
#> woops, > than 20! check your numbers
#> 
#> woops, > than 20! check your numbers
#>  [1]   1   4   9  16  25  36  49  64  81 100
```

## benchmark

definitely need to work on performance


```r
library(microbenchmark)
microbenchmark::microbenchmark(
  normal = suppressMessages(foo(1:10)),
  with_conditionz = suppressMessages(bar(1:10)),
  times = 100
)
#> Unit: microseconds
#>             expr      min        lq     mean    median       uq      max neval
#>           normal  902.391  939.6685 1071.314  957.6225 1076.281 3251.181   100
#>  with_conditionz 1719.505 1746.8675 1927.367 1788.8385 1993.597 4090.570   100
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/conditionz/issues).
* License: MIT
* Get citation information for `conditionz` in R doing `citation(package = 'conditionz')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
conditionz 0.1.0
================

### NEW FEATURES

* Released to CRAN
## Test environments

* local OS X install, R 3.5.2 patched
* ubuntu 14.04 (on travis-ci), R 3.5.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

New submission

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.

---

This is a new release. I have read and agree to the the 
CRAN policies at https://cran.r-project.org/web/packages/policies.html

Thanks!
Scott Chamberlain
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

* Submit an issue on the [Issues page](https://github.com/ropensci/conditionz/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/conditionz.git`
* Make sure to track progress upstream (i.e., on our version of `conditionz` at `ropensci/conditionz`) by doing `git remote add upstream https://github.com/ropensci/conditionz.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/conditionz`
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
conditionz
==========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/conditionz/workflows/R-check/badge.svg)](https://github.com/ropensci/conditionz/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/conditionz/coverage.svg?branch=master)](https://codecov.io/github/ropensci/conditionz?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/conditionz)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/conditionz)](https://cran.r-project.org/package=conditionz)

control how many times conditions are thrown

docs: https://ropensci.github.io/conditionz/

Package API:

```{r echo=FALSE, comment=NA, results="asis"}
cat(paste(" -", paste(sprintf("`%s`", getNamespaceExports("conditionz")), collapse = "\n - ")))
```

Use cases for `conditionz` functions:

- `ConditionKeeper` is what you want to use if you want to keep track of conditions inside a
function being applied many times, either in a for loop or lapply style fashion.
- `handle_conditions`/`handle_messages`/`handle_warnings` is what you want to use if the multiple
conditions are happening within a single function or code block
- `capture_message`/`capture_warning` are meant for capturing messages/warnings into a useable
list

## Installation

The CRAN version:

```{r eval=FALSE}
install.packages("conditionz")
```

Or the development version:

```{r eval=FALSE}
remotes::install_github("ropensci/conditionz")
```

```{r}
library("conditionz")
```

## ConditionKeeper

`ConditionKeeper` is the internal R6 class that handles keeping track of
conditions and lets us determine if conditions have been encountered,
how many times, etc.

```{r}
x <- ConditionKeeper$new(times = 4)
x
x$get_id()
x$add("one")
x$add("two")
x
x$thrown_already("one")
x$thrown_already("bears")
x$not_thrown_yet("bears")

x$add("two")
x$add("two")
x$add("two")
x$thrown_times("two")
x$thrown_enough("two")
x$thrown_enough("one")
```

## basic usage

A simple function that throws messages

```{r}
squared <- function(x) {
  stopifnot(is.numeric(x))
  y <- x^2
  if (y > 20) message("woops, > than 20! check your numbers")
  return(y)
}
foo <- function(x) {
  vapply(x, function(z) squared(z), numeric(1))
}
bar <- function(x, times = 1) {
  y <- ConditionKeeper$new(times = times)
  on.exit(y$purge())
  vapply(x, function(z) y$handle_conditions(squared(z)), numeric(1))
}
```

Running the function normally throws many messages

```{r}
foo(1:10)
```

Using in `ConditionKeeper` allows you to control how many messages
are thrown

```{r}
bar(x = 1:10)
```

```{r}
bar(1:10, times = 3)
```

## benchmark

definitely need to work on performance

```{r}
library(microbenchmark)
microbenchmark::microbenchmark(
  normal = suppressMessages(foo(1:10)),
  with_conditionz = suppressMessages(bar(1:10)),
  times = 100
)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/conditionz/issues).
* License: MIT
* Get citation information for `conditionz` in R doing `citation(package = 'conditionz')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConditionKeeper.R
\name{ConditionKeeper}
\alias{ConditionKeeper}
\title{ConditionKeeper}
\description{
R6 class with methods for handling conditions
}
\examples{
x <- ConditionKeeper$new(times = 4)
x
x$get_id()
x$add("one")
x$add("two")
x
x$thrown_already("one")
x$thrown_already("bears")
x$not_thrown_yet("bears")

x$add("two")
x$add("two")
x$add("two")
x$thrown_times("two")
x$thrown_enough("two")
x$thrown_enough("one")

foo <- function(x) {
  message("you gave: ", x)
  return(x)
}
foo('a')
x$handle_conditions(foo('a'))

x <- ConditionKeeper$new(times = 4, condition = "warning")
x
x$add("one")
x$add("two")
x
}
\seealso{
\code{\link[=handle_conditions]{handle_conditions()}}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{bucket}}{list holding conditions}

\item{\code{times}}{number of times}

\item{\code{condition}}{(character) type of condition, message or warning}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ConditionKeeper$new()}}
\item \href{#method-print}{\code{ConditionKeeper$print()}}
\item \href{#method-add}{\code{ConditionKeeper$add()}}
\item \href{#method-remove}{\code{ConditionKeeper$remove()}}
\item \href{#method-purge}{\code{ConditionKeeper$purge()}}
\item \href{#method-thrown_already}{\code{ConditionKeeper$thrown_already()}}
\item \href{#method-not_thrown_yet}{\code{ConditionKeeper$not_thrown_yet()}}
\item \href{#method-thrown_times}{\code{ConditionKeeper$thrown_times()}}
\item \href{#method-thrown_enough}{\code{ConditionKeeper$thrown_enough()}}
\item \href{#method-get_id}{\code{ConditionKeeper$get_id()}}
\item \href{#method-handle_conditions}{\code{ConditionKeeper$handle_conditions()}}
\item \href{#method-clone}{\code{ConditionKeeper$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{ConditionKeeper} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$new(times = 1, condition = "message")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{times}}{(integer) number of times to throw condition. required.
default: 1}

\item{\code{condition}}{(character) which condition, one of "message" (default)
or "warning"}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{ConditionKeeper} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{ConditionKeeper} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add"></a>}}
\if{latex}{\out{\hypertarget{method-add}{}}}
\subsection{Method \code{add()}}{
add a condition to internal storage
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$add(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a condition}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
self
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remove"></a>}}
\if{latex}{\out{\hypertarget{method-remove}{}}}
\subsection{Method \code{remove()}}{
remove the first condition from internal storage;
returns that condition so you know what you removed
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$remove()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the condition removed
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-purge"></a>}}
\if{latex}{\out{\hypertarget{method-purge}{}}}
\subsection{Method \code{purge()}}{
removes all conditions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$purge()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
NULL
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-thrown_already"></a>}}
\if{latex}{\out{\hypertarget{method-thrown_already}{}}}
\subsection{Method \code{thrown_already()}}{
has the condition been thrown already?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$thrown_already(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a condition}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-not_thrown_yet"></a>}}
\if{latex}{\out{\hypertarget{method-not_thrown_yet}{}}}
\subsection{Method \code{not_thrown_yet()}}{
has the condition NOT been thrown yet?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$not_thrown_yet(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a condition}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-thrown_times"></a>}}
\if{latex}{\out{\hypertarget{method-thrown_times}{}}}
\subsection{Method \code{thrown_times()}}{
number of times the condition has been thrown
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$thrown_times(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a condition}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
numeric
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-thrown_enough"></a>}}
\if{latex}{\out{\hypertarget{method-thrown_enough}{}}}
\subsection{Method \code{thrown_enough()}}{
has the condition been thrown enough? "enough" being:
thrown number of times equal to what you specified in the \code{times}
parameter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$thrown_enough(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{a condition}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_id"></a>}}
\if{latex}{\out{\hypertarget{method-get_id}{}}}
\subsection{Method \code{get_id()}}{
get the internal ID for the ConditionKeeper object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$get_id()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a UUID (character)
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-handle_conditions"></a>}}
\if{latex}{\out{\hypertarget{method-handle_conditions}{}}}
\subsection{Method \code{handle_conditions()}}{
pass a code block or function and handle conditions
within it
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$handle_conditions(expr)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{expr}}{an expression}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
the result of calling the expression
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ConditionKeeper$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/capture_condition.R
\name{capture_message}
\alias{capture_message}
\alias{capture_warning}
\alias{capture_condition}
\title{capture condition}
\usage{
capture_message(expr)

capture_warning(expr)
}
\description{
capture condition
}
\examples{
foom <- function(x) {
  message("its too bad")
  return(x)
}
capture_message(foom(4))

foow <- function(x) {
  warning("its too bad")
  return(x)
}
capture_warning(foow(4))
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_conditions.R
\name{handle_conditions}
\alias{handle_conditions}
\alias{handle_messages}
\alias{handle_warnings}
\title{Handle conditions}
\usage{
handle_conditions(expr, condition = "message", times = 1)

handle_messages(expr, times = 1)

handle_warnings(expr, times = 1)
}
\arguments{
\item{expr}{an expression}

\item{condition}{(character) one of "message" or "warning"}

\item{times}{(integer) max. times a condition should be thrown.
default: 1}
}
\value{
whatever the \code{expr} returns
}
\description{
Handle conditions
}
\details{
Uses \link{ConditionKeeper} internally
}
\examples{
foo <- function(x) {
  message("you gave: ", x)
  return(x)
}

foo('a')
capture_message(foo('a'))
handle_conditions(foo('a'))
suppressMessages(handle_conditions(foo('a')))
handle_conditions(foo('a'), "message")

bar <- function(x) {
  for (i in x) message("you gave: ", i)
  return(x)
}
bar(1:5)
handle_conditions(bar(1:5))

handle_messages(foo('a'))

hello <- function(x) {
  warning("you gave: ", x)
  return(x)
}
handle_warnings(hello('a'))

# code block
handle_warnings({
  as.numeric(letters[1:3])
  as.numeric(letters[4:6])
  as.numeric(letters[7:9])
})
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditionz-package.R
\docType{package}
\name{conditionz-package}
\alias{conditionz-package}
\alias{conditionz}
\title{conditionz}
\description{
condition control
}
\author{
Scott Chamberlain \href{mailto:myrmecocystus@gmail.com}{myrmecocystus@gmail.com}
}
\keyword{package}
