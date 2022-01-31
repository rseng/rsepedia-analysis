hoardr
======



[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/hoardr)](https://cranchecks.info/pkgs/hoardr)
[![R-check](https://github.com/ropensci/hoardr/workflows/R-check/badge.svg)](https://github.com/ropensci/hoardr/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/hoardr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/hoardr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/hoardr)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/hoardr)](https://cran.r-project.org/package=hoardr)


`hoard` - manage cached files

Exposes a single `R6` object so that when the package is imported in another
package for managing cached files, you don't need to pollute the NAMESPACE
with a bunch of functions. (you can always just `hoardr::fxn`, but
with a single object there are other benefits as well [maintaining state, e.g.]).

## install

stable


```r
install.packages("hoardr")
```

dev version


```r
remotes::install_github("ropensci/hoardr")
```


```r
library(hoardr)
```

## usage

initialize client


```r
(x <- hoardr::hoard())
#> <hoard> 
#>   path: 
#>   cache path:
```

set cache path


```r
x$cache_path_set("foobar", type = 'tempdir')
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpWAAkp3/R/foobar"
```

make the directory if doesn't exist


```r
x$mkdir()
```

put a file in the cache


```r
cat("hello world", file = file.path(x$cache_path_get(), "foo.txt"))
```

list the files


```r
x$list()
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpWAAkp3/R/foobar/foo.txt"
```

details


```r
x$details()
#> <cached files>
#>   directory: /var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmpWAAkp3/R/foobar
#> 
#>   file: /foo.txt
#>   size: 0 mb
```

delete by file name


```r
x$delete("foo.txt")
x$list()
#> character(0)
```

## todo

see [issue 1](https://github.com/ropensci/hoardr/issues/1)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/hoardr/issues).
* License: MIT
* Get citation information for `hoardr` in R doing `citation(package = 'hoardr')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/hoardr/blob/master/CODE_OF_CONDUCT.md
hoardr 0.5.2
============

### BUG FIXES

* Important fix: `HoardClient`, called by `hoardr()` function, was storing the cache path in an environment inside the R6 class. If multiple instances of `HoardClient` exist in the same R session, the cache path for any one then affects all others. Fixed by storing as a private variable int he R6 class instead of in an environment (#14)


hoardr 0.5.0
============

### NEW FEATURES

* Gains new method on the `HoardClient` object to check if one or more files exist, returning a data.frame (#10)
* `cache_path_set()` method on `HoardClient` gains new parameter `full_path` to make the base cache path directly with a full path rather than using the three other parameters (`path`, `type`, and `prefix`) (#12)


hoardr 0.2.0
============

### CHANGES

* Compliance with CRAN policies about writing to users disk (#6)

### MINOR IMPROVEMENTS

* Improved documentation (#7)

### BUG FIXES

* Change `key()` and `keys()` to use `file=TRUE` (#8)
* Fix R6 import warning (#5)


hoardr 0.1.0
============

### NEW FEATURES

* released to CRAN
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
## Test environments

* local OS X install, R 3.5.1 patched
* ubuntu 14.04 (on travis-ci), R 3.5.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

## Reverse dependencies

* I have run R CMD check on the 11 downstream dependencies
(<https://github.com/ropensci/hoardr/blob/master/revdep/README.md>).
No problems were found related to this package.

---

This submission includes a fix so that many objects in the same R session don't share variables anymore.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/hoardr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/hoardr.git`
* Make sure to track progress upstream (i.e., on our version of `hoardr` at `ropensci/hoardr`) by doing `git remote add upstream https://github.com/ropensci/hoardr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/hoardr`

### Check out our [discussion forum](https://discuss.ropensci.org)
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.5.1 Patched (2018-11-18 r75627) |
|os       |macOS Mojave 10.14.1                        |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2018-12-01                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|hoardr  |0.5.0 |0.5.2 |*  |

# Revdeps

## All (11)

|package                                  |version |error |warning |note |
|:----------------------------------------|:-------|:-----|:-------|:----|
|bomrang                                  |0.4.0   |      |        |     |
|crminer                                  |0.2.0   |      |        |     |
|[finch](problems.md#finch)               |0.2.0   |      |        |1    |
|fulltext                                 |1.1.0   |      |        |     |
|[getCRUCLdata](problems.md#getcrucldata) |0.2.5   |1     |        |     |
|rcoreoa                                  |0.3.0   |      |        |     |
|rdpla                                    |0.2.0   |      |        |     |
|rerddap                                  |0.4.2   |      |        |     |
|rnoaa                                    |0.7.0   |      |        |     |
|taxizedb                                 |0.1.4   |      |        |     |
|traits                                   |0.3.0   |      |        |     |

# finch

Version: 0.2.0

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘rappdirs’
      All declared Imports should be used.
    ```

# getCRUCLdata

Version: 0.2.5

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      1/1 mismatches
      [1] 12.9 - 12.9 == 5.72e-07
      
      ── 2. Failure: Test that create_stack creates tmn if requested (@test-create_CRU
      raster::maxValue(CRU_stack_list[[1]][[1]]) not equal to 4.3.
      1/1 mismatches
      [1] 4.3 - 4.3 == 1.91e-07
      
      ══ testthat results  ═══════════════════════════════════════════════════════════
      OK: 637 SKIPPED: 23 FAILED: 2
      1. Failure: Test that create_stack creates tmx if requested (@test-create_CRU_stack.R#868) 
      2. Failure: Test that create_stack creates tmn if requested (@test-create_CRU_stack.R#1233) 
      
      Error: testthat unit tests failed
      Execution halted
    ```

hoardr
======

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/hoardr)](https://cranchecks.info/pkgs/hoardr)
[![R-check](https://github.com/ropensci/hoardr/workflows/R-check/badge.svg)](https://github.com/ropensci/hoardr/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/hoardr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/hoardr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/hoardr)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/hoardr)](https://cran.r-project.org/package=hoardr)


`hoard` - manage cached files

Exposes a single `R6` object so that when the package is imported in another
package for managing cached files, you don't need to pollute the NAMESPACE
with a bunch of functions. (you can always just `hoardr::fxn`, but
with a single object there are other benefits as well [maintaining state, e.g.]).

## install

stable

```{r eval=FALSE}
install.packages("hoardr")
```

dev version

```{r eval=FALSE}
remotes::install_github("ropensci/hoardr")
```

```{r}
library(hoardr)
```

## usage

initialize client

```{r}
(x <- hoardr::hoard())
```

set cache path

```{r}
x$cache_path_set("foobar", type = 'tempdir')
```

make the directory if doesn't exist

```{r}
x$mkdir()
```

put a file in the cache

```{r}
cat("hello world", file = file.path(x$cache_path_get(), "foo.txt"))
```

list the files

```{r}
x$list()
```

details

```{r}
x$details()
```

delete by file name

```{r}
x$delete("foo.txt")
x$list()
```

## todo

see [issue 1](https://github.com/ropensci/hoardr/issues/1)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/hoardr/issues).
* License: MIT
* Get citation information for `hoardr` in R doing `citation(package = 'hoardr')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/hoardr/blob/master/CODE_OF_CONDUCT.md
---
title: "Introduction to the hoardr package"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
vignette: >
  %\VignetteIndexEntry{Introduction to the hoardr package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

hoardr introduction
===================

`hoardr` is a package for managing cached files.

The benefit of using `hoardr` vs. raw `rapddirs` is that `hoardr` exposes
an easy to use `R6` class that has variables and functions within it,
so you don't have to import function foo or bar, etc. Just a single
object.

You can easily wrap `hoardr` with user facing functions in your own
package to manage cached files.

If you find any bugs or have any feature requests get in touch at
<https://github.com/ropensci/hoardr>.

## Install

Stable from CRAN

```{r eval=FALSE}
install.packages("hoardr")
```

Dev version

```{r eval=FALSE}
devtools::install_github(c("ropensci/hoardr"))
```

```{r}
library("hoardr")
```

## initialize client

```{r}
(x <- hoardr::hoard())
```

## set cache path

```{r}
x$cache_path_set("foobar", type = 'tempdir')
```

## make the directory if doesn't exist

```{r}
x$mkdir()
```

## put a file in the cache

```{r}
cat("hello world", file = file.path(x$cache_path_get(), "foo.txt"))
```

## list the files

```{r}
x$list()
```

## details

```{r}
x$details()
```

## delete by file name

```{r}
x$delete("foo.txt")
x$list()
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hoard-package.R
\docType{package}
\name{hoardr-package}
\alias{hoardr-package}
\alias{hoardr}
\title{hoardr}
\description{
Manage Cached Files
}
\details{
\code{hoardr} is a tiny package with just a single export \code{\link[=hoard]{hoard()}}.
The package goal is to make it easy to setup, use, and manage cached
files for another R package. In addition, you can export functions in
your own package using \code{hoardr} for users to manage their cached
files.
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hoard_client.R
\name{hoard}
\alias{hoard}
\title{hoardr class}
\arguments{
\item{path}{(character) a path to cache files in. required}

\item{type}{(character) type of cache. One of "user_cache_dir" (default),
"user_log_dir", "user_data_dir", "user_config_dir", "site_data_dir",
"site_config_dir". Can also pass in any function that gives a path to a
directory, e.g., \code{tempdir()}. required.}
}
\description{
hoardr class
}
\details{
For the purposes of caching, you'll likely want to stick with
\code{user_cache_dir}, but you can change the type of cache with the \code{type}
parameter.

\code{hoard} is just a tiny wrapper around \code{HoardClient$new()}, which isn't
itself exported, but you can use it if you want via \code{:::}

\strong{Methods}
\describe{
\item{\code{cache_path_get()}}{
Get the cache path
\strong{return}: (character) path to the cache directory
}
\item{\code{cache_path_set(path = NULL, type = "user_cache_dir", prefix = "R", full_path = NULL)}}{
Set the cache path. By default, we set cache path to
\code{file.path(user_cache_dir, prefix, path)}. Note that this does not
actually make the directory, but just sets the path to it.
\itemize{
\item path (character) the path to be appended to the cache path set
by \code{type}
\item type (character) the type of cache, see \link{rappdirs}
\item prefix (character) prefix to the \code{path} value. Default: "R"
\item full_path (character) instead of using \code{path}, \code{type}, and \code{prefix}
just set the full path with this parameter
}
\strong{return}: (character) path to the cache directory just set
}
\item{\code{list()}}{
List files in the directory (full file paths)
\strong{return}: (character) vector of file paths for files in the cache
}
\item{\code{mkdir()}}{
Make the directory if doesn't exist already
\strong{return}: \code{TRUE}, invisibly
}
\item{\code{delete(files, force = TRUE)}}{
Delete files by name
\itemize{
\item files (character) vector/list of file paths
\item force (logical) force deletion? Default: \code{TRUE}
}
\strong{return}: nothing
}
\item{\code{delete_all(force = TRUE)}}{
Delete all files
\itemize{
\item force (logical) force deletion? Default: \code{FALSE}
}
\strong{return}: nothing
}
\item{\code{details(files = NULL)}}{
Get file details
\itemize{
\item files (character) vector/list of file paths
}
\strong{return}: objects of class \code{cache_info}, each with brief summary
info including file path and file size
}
\item{\code{keys(algo = "md5")}}{
Get a hash for all files. Note that these keys may not be unique
if the files are identical, leading to identical hashes
\strong{return}: (character) hashes for the files
}
\item{\code{key(x, algo = "md5")}}{
Get a hash for a single file. Note that these keys may not be unique
if the files are identical, leading to identical hashes
\itemize{
\item x (character) path to a file
\item algo (character) the algorithm to be used, passed on to
\code{\link[digest:digest]{digest::digest()}}, choices: md5 (default), sha1, crc32, sha256,
sha512, xxhash32, xxhash64 and murmur32.
}
\strong{return}: (character) hash for the file
}
\item{\code{files()}}{
Get all files as HoardFile objects
\strong{return}: (character) paths to the files
}
\item{\code{compress()}}{
Compress files into a zip file - leaving only the zip file
\strong{return}: (character) path to the cache directory
}
\item{\code{uncompress()}}{
Uncompress all files and remove zip file
\strong{return}: (character) path to the cache directory
}
\item{\code{exists(files)}}{
Check if files exist
\itemize{
\item files: (character) one or more files, paths are optional
}
\strong{return}: (data.frame) with two columns:
\itemize{
\item files: (character) file path
\item exists: (boolean) does it exist or not
}
}
}
}
\examples{
(x <- hoard())
x$cache_path_set(path = "foobar", type = 'tempdir')
x
x$path
x$cache_path_get()

# Or you can set the full path directly with `full_path`
mydir <- file.path(tempdir(), "foobar")
x$cache_path_set(full_path = mydir)
x
x$path
x$cache_path_get()

# make the directory if doesn't exist already
x$mkdir()

# list files in dir
x$list()
cat(1:10000L, file = file.path(x$cache_path_get(), "foo.txt"))
x$list()

# add more files
cat(letters, file = file.path(x$cache_path_get(), "foo2.txt"))
cat(LETTERS, file = file.path(x$cache_path_get(), "foo3.txt"))

# see if files exist
x$exists("foo.txt") # exists
x$exists(c("foo.txt", "foo3.txt")) # both exist
x$exists(c("foo.txt", "foo3.txt", "stuff.txt")) # one doesn't exist

# cache details
x$details()

# delete files by name - we prepend the base path for you
x$delete("foo.txt")
x$list()
x$details()

# delete all files
cat("one\ntwo\nthree", file = file.path(x$cache_path_get(), "foo.txt"))
cat("asdfasdf asd fasdf", file = file.path(x$cache_path_get(), "bar.txt"))
x$delete_all()
x$list()

# make/get a key for a file
cat(1:10000L, file = file.path(x$cache_path_get(), "foo.txt"))
x$keys()
x$key(x$list()[1])

# as files
Map(function(z) z$exists(), x$files())

# compress and uncompress
x$compress()
x$uncompress()

# reset cache path
x$cache_path_set(path = "stuffthings", type = "tempdir")
x
x$cache_path_get()
x$list()

# cleanup
unlink(x$cache_path_get())
}
