rerddap
=====



[![cran checks](https://cranchecks.info/badges/worst/rerddap)](https://cranchecks.info/pkgs/rerddap)
[![R-check](https://github.com/ropensci/rerddap/workflows/R-check/badge.svg)](https://github.com/ropensci/rerddap/actions)
[![codecov.io](https://codecov.io/github/ropensci/rerddap/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rerddap?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rerddap)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rerddap)](https://cran.r-project.org/package=rerddap)

`rerddap` is a general purpose R client for working with ERDDAP servers.

Package Docs: <https://docs.ropensci.org/rerddap/>

## Installation

From CRAN


```r
install.packages("rerddap")
```

Or development version from GitHub


```r
remotes::install_github("ropensci/rerddap")
```


```r
library("rerddap")
```

Some users may experience an installation error, stating to install 1 or more 
packages, e.g., you may need `DBI`, in which case do, for example, 
`install.packages("DBI")` before installing `rerddap`.

## Background

ERDDAP is a server built on top of OPenDAP, which serves some NOAA data. You can get gridded data (griddap (<https://upwell.pfeg.noaa.gov/erddap/griddap/documentation.html>)), which lets you query from gridded datasets, or table data (tabledap (<https://upwell.pfeg.noaa.gov/erddap/tabledap/documentation.html>)) which lets you query from tabular datasets. In terms of how we interface with them, there are similarities, but some differences too. We try to make a similar interface to both data types in `rerddap`.

## NetCDF

`rerddap` supports NetCDF format, and is the default when using the `griddap()` function. NetCDF is a binary file format, and will have a much smaller footprint on your disk than csv. The binary file format means it's harder to inspect, but the `ncdf4` package makes it easy to pull data out and write data back into a NetCDF file. Note the the file extension for NetCDF files is `.nc`. Whether you choose NetCDF or csv for small files won't make much of a difference, but will with large files.

## Caching

Data files downloaded are cached in a single directory on your machine determined by the `hoardr` package. When you use `griddap()` or `tabledap()` functions, we construct a MD5 hash from the base URL, and any query parameters - this way each query is separately cached. Once we have the hash, we look in the cache directory for a matching hash. If there's a match we use that file on disk - if no match, we make a http request for the data to the ERDDAP server you specify.

## ERDDAP servers

You can get a data.frame of ERDDAP servers using the function `servers()`. Most I think serve some kind of NOAA data, but there are a few that aren't NOAA data.  If you know of more ERDDAP servers, send a pull request, or let us know.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rerddap/issues).
* License: MIT
* Get citation information for `rerddap` in R doing `citation(package = 'rerddap')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rerddap 0.8.0
=============

* Added global search function
* fixed bug when dataset has a decreasing coordinate that
is not latitude or longitude


rerddap 0.7.6
=============

### MINOR IMPROVEMENTS

* fixed a bug in dealing with trailing slashes in URLs

rerddap 0.7.4
=============

### MINOR IMPROVEMENTS

* fix a broken test

rerddap 0.7.0
=============

### MINOR IMPROVEMENTS

* vignettes only on package documentation site now  (#87)
* `server()` (to fetch known ERDDAP server URLs) now uses the list maintained by `irishmarineinstitute/awesome-erddap` on GitHub (#86)
* better error handling for `griddap()`: if no dimension arguments passed, we error saying so (and no http requests made); in addition, if a dataset is passed to `griddap()`, to which the output of `info()` was also passed, then we can check if the dataset has griddap data or not, and fail saying so if not (#91)
* `griddap()` and `tabledap()`: if `info()` output passed to these two funcitons, we will now use the url within that info output, and use a message telling the user we are doing so; now you don't have to set the url if you pass info output  (#92)


rerddap 0.6.5
=============

### BUG FIXES

* fix a `convert_units` test that was failing because remote service had changed the response


rerddap 0.6.4
=============

### BUG FIXES

* fix to internal fxn `err_handle()` for handling http errors - ERDDAP servers changed to some weird JSON-ish type format (#85)


rerddap 0.6.0
=============

### MINOR IMPROVEMENTS

* change all `tibble::as_data_frame`/`tibble::data_frame` to `tibble::as_tibble` (#79)
* `info()` gains new element in its output list, `base_url`, the base url for the ERDDAP server under consideration (#80)
* improved docs for `griddap()` with respect to what's returned from the function  (#81)
* fix some test fixtures to use preserve exact bytes so that cran checks on debian clang devel don't fail (#83)
* add .github files: contributing, issue template, pull request template

### BUG FIXES

* fix for lat/lon parsing within `griddap()` to account for cases when min and max are reversed from the order they should be in (#78)
* fix to `griddap()` to parse additioanl dimensions returned; previously we were only returning time, lat, and lon, plus one more (#82) thanks @afredstonhermann


rerddap 0.5.0
=============

### MINOR IMPROVEMENTS

* added new `Caching` section to package level manual file (`?rerddap`) about caching  (#52)
* use markdown docs in package (#75)
* replace `httr` with `crul` (#54)
* cache most tests with HTTP requests using `vcr` (#76)
* add test for `read` parameter in `griddap()` (#47)
* use default url via `eurl()`; used as default in main functions; set default url with env vars, see `?eurl`  (#41)
* improve handling and reporting back to user of ERDDAP server errors (#70) (#73)
* change to `griddap()`: when nc format gridded datasets have latitude and longitude we "melt" them into a data.frame for easy downstream consumption. When nc format gridded datasets do not have latitude and longitude components, we do not read in the data, throw a warning saying so. You can readin the nc file yourself with the file path (#74)
* for for `griddap()` to support cases in wihch lat/lon runs north to south and south to north (#68)

### BUG FIXES

* `memory()` usage in `griddap()` wasn't working. fixed now (#77)


rerddap 0.4.2
=============

### NEW FEATURES

* Now using `hoardr` to manage caching paths and such (#60). Also
now asking users where they want to cache files, either in a 
`rappdirs` user cache dir or a temp directory. Now on tests and examples
we use temp dirs.
* Related to above, new functions `cache_info()` to get cache path and 
number of cached files, and `cache_setup()` to set cache path.
* Related to above, `cache_details()`, `cache_list()`, and `cache_delete()`
lose their `cache_path` parameter - now cache path is set package wide and 
we use the same cache path, so no need to set in the fxn call.

### MINOR IMPROVEMENTS

* Fixes to a number of `griddap()` and `tabledap()` examples to use 
datasets that still exist (previous examples used datasets that are no
gone)


rerddap 0.4.0
=============

### NEW FEATURES

* New vignette added that goes in to much more depth than 
the original vignette (#51) thx to @rmendels
* `info()` function gains new attribute `url` with the 
base url for the ERDDAP server used (#42)
* Replaced usage of internal compact data.frame code to 
use `tibble` package (#45)

### MINOR IMPROVEMENTS

* Added another ERDDAP server to `servers()` function (#49)
* Changed base URLs for default ERDDAP server from `http` 
to `https`  (#50)
* Added note to docs for `griddap()` and `tabledap()` for how
to best deal with 500 server errors (#48)
* Replaced all `dplyr::rbind_all` uses with `dplyr::bind_rows` (#46)


rerddap 0.3.4
=============

### MINOR IMPROVEMENTS

* Removed use of `ncdf` package, which has been taken off CRAN.
Using `ncdf4` now for all NetCDF file manipulation. (#35)
* Failing better now with custom error catching (#31)
* Added many internal checks for parameter inputs, warning or
stopping as necessary - ERDDAP servers silently drop with no
informative messages (#32)

### BUG FIXES

* Using now `file.info()$size` instead of `file.size()` to be
backwards compatible with R versions < 3.2


rerddap 0.3.0
=============

### NEW FEATURES

* Cache functions accept the outputs of `griddap()` and `tabledap()`
so that the user can easily see cache details or delete the file from
the cache without having to manually get the file name. (#30)

### MINOR IMPROVEMENTS

* All package dependencies now use `importFrom` so we only import
functions we need instead of their global namespaces.

### BUG FIXES

* Fixed bug in parsing data from netcdf files, affected the
`griddap()` function (#28)


rerddap 0.2.0
=============

### NEW FEATURES

* Added a suite of functions to manage local cached files (#17)

### MINOR IMPROVEMENTS

* Added new ERDDAP server to list of servers in the `servers()` function (#21)

### BUG FIXES

* Fixed a few cases across a number of functions in which an empty list
passed to `query` parmaeter in `httr::GET()` caused an error (#23)
* Fixed retrieval of path to file written to disk by `httr::write_disk()` (#24)
* `last` is a value accepted by ERDDAP servers, but internal functions
weren't checking correctly, fixed now. (#25)
* `as.info()` wasn't passing on the `url` parameter to the `info()` function.
fixed now. (#26)


rerddap 0.1.0
=============

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.1.1
* rhub ubuntu, debian, solaris, Apple Silicon (M1), macOS 11.6 Big Sur
* win-builder (devel and release)

## R CMD check results

OK from all checks

## Reverse dependencies

* No probelms with reverse dependencies.

---

This version fixes a bug dealing with coordinates that
are not lat-lon and are in decreasing order.  Adds a new
search function.

Thanks! 
Roy Mendelssohn
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

*  CORRECT: you edit a roxygen comment in a `.R` file below `R/`.
*  INCORRECT: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes. The `README` should contain badges for any continuous integration services used by the package.
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the `rerddadp` project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html) for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) DO NOT SHARE SCREENSHOTS OF CODE; SHARE ACTUAL CODE IN TEXT FORMAT -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.4 Patched (2021-02-17 r80031) |
|os       |macOS Big Sur 10.16                         |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-03-05                                  |

# Dependencies

|package   |old   |new      |Δ  |
|:---------|:-----|:--------|:--|
|rerddap   |0.7.0 |0.7.1.92 |*  |
|cli       |NA    |2.3.1    |*  |
|crul      |NA    |1.1.0    |*  |
|dplyr     |NA    |1.0.4    |*  |
|lifecycle |NA    |1.0.0    |*  |
|pillar    |NA    |1.5.1    |*  |
|tibble    |NA    |3.1.0    |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*