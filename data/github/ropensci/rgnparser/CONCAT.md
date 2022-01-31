rgnparser
=========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rgnparser)](https://cranchecks.info/pkgs/rgnparser)
[![R-check](https://github.com/ropensci/rgnparser/workflows/R-check/badge.svg)](https://github.com/ropensci/rgnparser/actions/)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rgnparser)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rgnparser)](https://cran.r-project.org/package=rgnparser)

**rgnparser**: Parse Scientific Names

Docs: https://ropensci.github.io/rgnparser/

## Installation


```r
install.packages("rgnparser")
# OR
remotes::install_github("ropensci/rgnparser")
```


```r
library("rgnparser")
```

## Install gnparser

The command line tool written in Go, gnparser, is required to use this package.

If you want to install gnparser on your own, instructions can be found at the
gnparser repo (https://github.com/gnames/gnparser)

There is a helper function in **rgnparser** for downloading and installing
gnparser on major operating systems (macOS, Windows, Linux):


```r
rgnparser::install_gnparser()
```

It installs the latest gnparser version by default, but you can specify which 
version to install. You can also install gnparser outside of R yourself
(see above).


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rgnparser/issues).
* License: MIT
* Get citation information for `rgnparser` in R doing `citation(package = 'rgnparser')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rgnparser 0.2.0
===============

### NEW FEATURES

* A new gnparser version (v1) is out. In addition, gnparser has moved from Gitlab to Github; which also required changes to `install_gnparser()` because we download Go binaries from the gnparser source repository (#7)
* As part of new gnparser version, arguments have changed in `gn_parse()` and `gn_parse_tidy()`: `format` has been removed. `batch_size` and `ignore_tags` were added to both functions, while `details` was added to `gn_parse()` only. See docs for details.  (#11)
* gnparser v1 or greater is now required (#10)

### DEFUNCT

* `gn_debug()` is now defunct. the gnparser command for this function was removed in gnparser v1 (#9)

### BUG FIXES

* `gn_version()` was broken with the new gnparser version, fixed now (#8)
* xxx (#xx)


rgnparser 0.1.0
===============

### NEW FEATURES

* First submission to CRAN.
## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (release and devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

-----

This submission makes appropriate changes for a new major version of the gnparser Go command line tool.

Thanks!
Scott Chamberlain
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

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

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the `pegax` project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
rgnparser
=========

```{r, echo=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rgnparser)](https://cranchecks.info/pkgs/rgnparser)
[![R-check](https://github.com/ropensci/rgnparser/workflows/R-check/badge.svg)](https://github.com/ropensci/rgnparser/actions/)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rgnparser)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rgnparser)](https://cran.r-project.org/package=rgnparser)

**rgnparser**: Parse Scientific Names

Docs: https://ropensci.github.io/rgnparser/

## Installation

```{r eval=FALSE}
install.packages("rgnparser")
# OR
remotes::install_github("ropensci/rgnparser")
```

```{r}
library("rgnparser")
```

## Install gnparser

The command line tool written in Go, gnparser, is required to use this package.

If you want to install gnparser on your own, instructions can be found at the
gnparser repo (https://github.com/gnames/gnparser)

There is a helper function in **rgnparser** for downloading and installing
gnparser on major operating systems (macOS, Windows, Linux):

```{r eval=FALSE}
rgnparser::install_gnparser()
```

It installs the latest gnparser version by default, but you can specify which 
version to install. You can also install gnparser outside of R yourself
(see above).


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rgnparser/issues).
* License: MIT
* Get citation information for `rgnparser` in R doing `citation(package = 'rgnparser')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "rgnparser"
author: "Scott Chamberlain"
date: "2021-01-21"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to rgnparser}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



rgnparser: Parse Scientific Names

An R interface to `gnparser` at https://github.com/gnames/gnparser

## Installation


```r
install.packages("rgnparser")
# OR
remotes::install_github("ropensci/rgnparser")
```


```r
library(rgnparser)
```

## Install gnparser

The command line tool written in Go, gnparser, is required to use this package.

If you want to install gnparser on your own, instructions can be found at the [gnparser repo][gnparser].

There is a helper function in **rgnparser** for downloading and installing gnparser on major
operating systems (macOS, Windows, Linux):


```r
rgnparser::install_gnparser()
```

It installs the latest gnparser version by default, but you can specify which 
version to install. You can also install gnparser outside of R yourself (see above).

## gn_version

gnparser version


```r
gn_version()
#> $version
#> [1] "v1.0.0"
#> 
#> $build
#> [1] "2021-01-19_14:45:28UTC"
```


## gn_parse_tidy

output a data.frame with more minimal information


```r
x <- c("Quadrella steyermarkii (Standl.) Iltis &amp; Cornejo",
  "Parus major Linnaeus, 1788", "Helianthus annuus var. texanus")
gn_parse_tidy(x)
#> # A tibble: 3 x 9
#>   id    verbatim cardinality canonicalstem canonicalsimple canonicalfull
#>   <chr> <chr>          <dbl> <chr>         <chr>           <chr>        
#> 1 3e33… Quadrel…           2 Quadrella st… Quadrella stey… Quadrella st…
#> 2 e4e1… Parus m…           2 Parus maior   Parus major     Parus major  
#> 3 e571… Heliant…           3 Helianthus a… Helianthus ann… Helianthus a…
#> # … with 3 more variables: authorship <chr>, year <dbl>, quality <dbl>
```

It's pretty fast, thanks to gnparser of course


```r
n <- 10000L
# get random scientific names from taxize
spp <- taxize::names_list(rank = "species", size = n)
timed <- system.time(gn_parse_tidy(spp))
timed
#>    user  system elapsed 
#>   1.225   0.113   0.555
```

Just 0.555 sec for 10000 names

## gn_parse

output a list of lists with more detailed information


```r
x <- c("Quadrella steyermarkii (Standl.) Iltis &amp; Cornejo",
  "Parus major Linnaeus, 1788", "Helianthus annuus var. texanus")
gn_parse(x)
#> [[1]]
#> [[1]]$parsed
#> [1] TRUE
#> 
#> [[1]]$quality
#> [1] 3
#> 
#> [[1]]$qualityWarnings
#>   quality                           warning
#> 1       3 HTML tags or entities in the name
#> 
#> [[1]]$verbatim
#> [1] "Quadrella steyermarkii (Standl.) Iltis & Cornejo"
#> 
#> [[1]]$normalized
...
```

[gnparser]: https://github.com/gnames/gnparser
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gn_parse_tidy.R
\name{gn_parse_tidy}
\alias{gn_parse_tidy}
\title{gn_parse_tidy}
\usage{
gn_parse_tidy(x, threads = 4, batch_size = NULL, ignore_tags = FALSE)
}
\arguments{
\item{x}{(character) vector of scientific names. required}

\item{threads}{(integer/numeric) number of threads to run. CPU's
threads number is the default. default: \code{4}}

\item{batch_size}{(integer/numeric) maximum number of names in a
batch send for processing. default: \code{NULL}}

\item{ignore_tags}{(logical) ignore HTML entities and tags when
parsing. default: \code{FALSE}}
}
\value{
a data.frame
}
\description{
extract names using gnparser into a tidy tibble
}
\details{
This function focuses on a data.frame result that's easy
to munge downstream - note that this function does not do additional
details as does \code{\link[=gn_parse]{gn_parse()}}.
}
\examples{
trys <- function(x) try(x, silent=TRUE)
if (interactive()) {
x <- c("Quadrella steyermarkii (Standl.) Iltis & Cornejo",
  "Parus major Linnaeus, 1788", "Helianthus annuus var. texanus")
trys(gn_parse_tidy(x))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install_gnparser.R
\name{install_gnparser}
\alias{install_gnparser}
\title{Install gnparser}
\usage{
install_gnparser(version = "latest", force = FALSE)
}
\arguments{
\item{version}{The gnparser version number, e.g., \verb{1.0.0}; the default
\code{latest} means the latest version (fetched from GitLab releases).
Alternatively, this argument can take a file path of the zip archive or
tarball of gnparser that has already been downloaded from GitLab,
in which case it will not be downloaded again. The minimum version
is \code{v1.0.0} because gnparser v1 introduced breaking changes - and
we don't support older versions of gnparser here.}

\item{force}{Whether to install gnparser even if it has already been
installed. This may be useful when upgrading gnparser.}
}
\description{
Downloads the appropriate gnparser executable for your platform and
tries to copy it to a system directory so \pkg{rgnparser} can run the
\code{gnparser} command.
}
\details{
This function tries to install gnparser to \code{Sys.getenv('APPDATA')} on
Windows, \verb{~/Library/Application Support} on macOS, and \verb{~/bin/} on
other platforms (such as Linux). If these directories are not writable, the
package directory \code{gnparser} of \pkg{rgnparser} will be used. If it
still fails, you have to install gnparser by yourself and make sure it can
be found via the environment variable \code{PATH}.

This is just a helper function and may fail to choose the correct gnparser
executable for your operating system, especially if you are not on Windows
or Mac or a major Linux distribution. When in doubt, read the gnparser
documentation and install it yourself:
https://github.com/gnames/gnparser#installation
}
\note{
modified from \code{blogdown::install_hugo}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gn_parse.R
\name{gn_parse}
\alias{gn_parse}
\title{gn_parse}
\usage{
gn_parse(
  x,
  threads = NULL,
  batch_size = NULL,
  ignore_tags = FALSE,
  details = FALSE
)
}
\arguments{
\item{x}{(character) vector of scientific names. required}

\item{threads}{(integer/numeric) number of threads to run. CPU's
threads number is the default. default: \code{4}}

\item{batch_size}{(integer/numeric) maximum number of names in a
batch send for processing. default: \code{NULL}}

\item{ignore_tags}{(logical) ignore HTML entities and tags when
parsing. default: \code{FALSE}}

\item{details}{(logical) Return more details for a parsed name}
}
\value{
a list
}
\description{
extract names using gnparser
}
\examples{
trys <- function(x) try(x, silent=TRUE)
if (interactive()) {
x <- c("Quadrella steyermarkii (Standl.) Iltis & Cornejo",
  "Parus major Linnaeus, 1788", "Helianthus annuus var. texanus")
trys(gn_parse(x[1]))
trys(gn_parse(x[2]))
trys(gn_parse(x[3]))
trys(gn_parse(x))
# details
w <- trys(gn_parse(x, details = TRUE))
w[[1]]$details # details for one name
lapply(w, "[[", "details") # details for all names
z <- trys(gn_parse(x, details = FALSE)) # compared to regular
z
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gn_version.R
\name{gn_version}
\alias{gn_version}
\title{gn_version}
\usage{
gn_version()
}
\value{
named list, with \code{version} and \code{build}
}
\description{
get gnparser version information
}
\examples{
trys <- function(x) try(x, silent=TRUE)
if (interactive()) {
trys(gn_version())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgnparser-package.R
\docType{package}
\name{rgnparser-package}
\alias{rgnparser-package}
\alias{rgnparser}
\title{rgnparser}
\description{
Parse scientific names using gnparser
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gn_debug.R
\name{gn_debug}
\alias{gn_debug}
\title{gn_debug}
\usage{
gn_debug(...)
}
\arguments{
\item{...}{ignored}
}
\description{
DEFUNCT
}
