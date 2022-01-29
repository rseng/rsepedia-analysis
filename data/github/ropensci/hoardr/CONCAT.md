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

