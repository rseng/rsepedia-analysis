cchecks
=======



[![R-CMD-check](https://github.com/ropensci/cchecks/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cchecks/actions/)
[![codecov.io](https://codecov.io/github/ropensci/cchecks/coverage.svg?branch=master)](https://codecov.io/github/ropensci/cchecks?branch=master)

R client for the CRAN checks API at <https://cranchecks.info>

[CRAN checks API docs][docs]

authentication is only needed for the CRAN checks API for the functions that start with `cchn`

See https://docs.ropensci.org/cchecks for full documentation on `cchecks`

## Install


```r
remotes::install_github("ropensci/cchecks")
```


```r
library("cchecks")
```

## heartbeat

- `cch_heartbeat()`

## packages

- current day package check data: `cch_pkgs()` or `cch_pkgs("packagename")`
- historical package check data (30 days back): `cch_pkgs_history()`
- historical package check data for all packages by day: `cch_history()`
- search historical data: `cch_pkgs_search()`

There's an important shortcoming of historical data. The links in the historical
data in the `checks` field are not date specific. If you go to a link in historical
data, for example for April 2nd, 2020, links in that set of data link to whatever
the current check data is for that package. The `check_details` field is
date specific though; the text is scraped from the package checks page each day
and stored, so you can count on that to be date specific. There are sometimes
links to further checks, often of compiled packages on various types of checks
that CRAN runs; we do not have those check results - we could get them but 
have not take the time to sort that out.


## maintainers

- all maintainers: `cch_maintainers()`
- maintainers by email: `cch_maintainers("maelle.salmon_at_yahoo.se")`

## notifications

Functions for working with notifications are all prefixed with `cchn`. 

`cchn` functions are designed to be used from within an R package directory. The functions
look for the package name and maintainer email address. Functions copy heavily from 
https://github.com/r-hub/rhub

The functions

- `cchn_register()`: registration
- `cchn_pkg_rule_list()`/`cchn_rule_list()`: list your own rules
- `cchn_pkg_rule_get()`/`cchn_rule_get()`: get a rule by id
- `cchn_pkg_rule_add()`/`cchn_rule_add()`: create a rule
- `cchn_pkg_rule_delete()`/`cchn_rule_delete()`: delete a rule by id (get id from `cchn_pkg_rule_list`/`cchn_rule_list`)

Functions prefixed with `cchn_pkg_` operate within a package 
directory. That is, your current working directory is an R
package, and is the package for which you want to handle CRAN checks 
notifications. These functions make sure that you are inside of 
an R package, and use the email address and package name based
on the directory you're in.

Functions prefixed with just `cchn_` do not operate within a package.
These functions do not guess package name at all, but require the user
to supply a package name (for those functions that require a package name);
and instead of guessing an email address from your package, we guess email
from the cached cchecks email file (see `?cchn_register`).

The first thing to do is to register an email address. In an R session in a working directory for 
one of your packages that is on CRAN, run `cchn_register()`. This function:

- registers your email address
- a validation email is sent to you right away with your token
- paste the token from the amil into the R session
- the email and token are then saved in a file on your machine
- all `cchn_rule*` functions use this cached token
- you don't need to pass the token in any `cchn` function calls

If you run `cchn_register()` in the same package directory (with the same email address), 
you'll be issued a new token, which will be updated in your cached token file.

It's entirely possible to have more than one email address you use across different R packages. 
If you run `cchn_register()` in a different package directory (with a different email address
from a previous run of `cchn_register()`), you'll be issued a different token associated 
with that new email address.

See `?cchn_rules` for details on how the rules work and many examples of adding rules.

Note that you can only manage your own rules. You can not list, get, or delete rules 
of other users.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/cchecks/issues).
* License: MIT
* Get citation information for `cchecks` in R doing `citation(package = 'cchecks')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[docs]: https://cranchecks.info/docs
[coc]: https://github.com/ropensci/cchecks/blob/master/CODE_OF_CONDUCT.md
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

* Submit an issue on the [Issues page](https://github.com/ropensci/cchecks/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/cchecks.git`
* Make sure to track progress upstream (i.e., on our version of `cchecks` at `ropensci/cchecks`) by doing `git remote add upstream https://github.com/ropensci/cchecks.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/cchecks`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
cchecks
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![R-CMD-check](https://github.com/ropensci/cchecks/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cchecks/actions/)
[![codecov.io](https://codecov.io/github/ropensci/cchecks/coverage.svg?branch=master)](https://codecov.io/github/ropensci/cchecks?branch=master)

R client for the CRAN checks API at <https://cranchecks.info>

[CRAN checks API docs][docs]

authentication is only needed for the CRAN checks API for the functions that start with `cchn`

See https://docs.ropensci.org/cchecks for full documentation on `cchecks`

## Install

```{r eval=FALSE}
remotes::install_github("ropensci/cchecks")
```

```{r}
library("cchecks")
```

## heartbeat

- `cch_heartbeat()`

## packages

- current day package check data: `cch_pkgs()` or `cch_pkgs("packagename")`
- historical package check data (30 days back): `cch_pkgs_history()`
- historical package check data for all packages by day: `cch_history()`
- search historical data: `cch_pkgs_search()`

There's an important shortcoming of historical data. The links in the historical
data in the `checks` field are not date specific. If you go to a link in historical
data, for example for April 2nd, 2020, links in that set of data link to whatever
the current check data is for that package. The `check_details` field is
date specific though; the text is scraped from the package checks page each day
and stored, so you can count on that to be date specific. There are sometimes
links to further checks, often of compiled packages on various types of checks
that CRAN runs; we do not have those check results - we could get them but 
have not take the time to sort that out.


## maintainers

- all maintainers: `cch_maintainers()`
- maintainers by email: `cch_maintainers("maelle.salmon_at_yahoo.se")`

## notifications

Functions for working with notifications are all prefixed with `cchn`. 

`cchn` functions are designed to be used from within an R package directory. The functions
look for the package name and maintainer email address. Functions copy heavily from 
https://github.com/r-hub/rhub

The functions

- `cchn_register()`: registration
- `cchn_pkg_rule_list()`/`cchn_rule_list()`: list your own rules
- `cchn_pkg_rule_get()`/`cchn_rule_get()`: get a rule by id
- `cchn_pkg_rule_add()`/`cchn_rule_add()`: create a rule
- `cchn_pkg_rule_delete()`/`cchn_rule_delete()`: delete a rule by id (get id from `cchn_pkg_rule_list`/`cchn_rule_list`)

Functions prefixed with `cchn_pkg_` operate within a package 
directory. That is, your current working directory is an R
package, and is the package for which you want to handle CRAN checks 
notifications. These functions make sure that you are inside of 
an R package, and use the email address and package name based
on the directory you're in.

Functions prefixed with just `cchn_` do not operate within a package.
These functions do not guess package name at all, but require the user
to supply a package name (for those functions that require a package name);
and instead of guessing an email address from your package, we guess email
from the cached cchecks email file (see `?cchn_register`).

The first thing to do is to register an email address. In an R session in a working directory for 
one of your packages that is on CRAN, run `cchn_register()`. This function:

- registers your email address
- a validation email is sent to you right away with your token
- paste the token from the amil into the R session
- the email and token are then saved in a file on your machine
- all `cchn_rule*` functions use this cached token
- you don't need to pass the token in any `cchn` function calls

If you run `cchn_register()` in the same package directory (with the same email address), 
you'll be issued a new token, which will be updated in your cached token file.

It's entirely possible to have more than one email address you use across different R packages. 
If you run `cchn_register()` in a different package directory (with a different email address
from a previous run of `cchn_register()`), you'll be issued a different token associated 
with that new email address.

See `?cchn_rules` for details on how the rules work and many examples of adding rules.

Note that you can only manage your own rules. You can not list, get, or delete rules 
of other users.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/cchecks/issues).
* License: MIT
* Get citation information for `cchecks` in R doing `citation(package = 'cchecks')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[docs]: https://cranchecks.info/docs
[coc]: https://github.com/ropensci/cchecks/blob/master/CODE_OF_CONDUCT.md
---
title: Client for the cranchecks.info API
author: Scott Chamberlain
date: "2020-05-19"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Client for the cranchecks.info API}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



`cchecks` is a client for the <https://cranchecks.info> API

## Install


```r
remotes::install_github("ropensci/cchecks")
```


```r
library("cchecks")
```

## heartbeat


```r
cch_heartbeat()
#> $routes
#>  [1] "/"                               "/docs"                          
#>  [3] "/heartbeat/?"                    "/pkgs"                          
#>  [5] "/pkgs/:name"                     "/maintainers"                   
#>  [7] "/maintainers/:email"             "/badges/:type/:package"         
#>  [9] "/badges/flavor/:flavor/:package" "/pkgs/:name/history"            
#> [11] "/history/:date"                  "/search"                        
#> [13] "/notifications/rules"            "/notifications/rules/:id"
```

## packages

all


```r
cch_pkgs(limit = 1)
#> $found
#> [1] 16617
#> 
#> $count
#> [1] 1
#> 
#> $offset
#> [1] 0
#> 
#> $error
#> NULL
#> 
#> $data
#>   package                                                               url
#> 1 localIV https://cloud.r-project.org/web/checks/check_results_localIV.html
#>               date_updated summary.any summary.ok summary.note summary.warn
#> 1 2020-05-20T00:04:00.609Z       FALSE         12            0            0
#>   summary.error summary.fail
#> 1             0            0
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       checks
#> 1 r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-debian-gcc, r-devel-linux-x86_64-fedora-clang, r-devel-linux-x86_64-fedora-gcc, r-devel-windows-ix86+x86_64, r-patched-linux-x86_64, r-patched-solaris-x86, r-release-linux-x86_64, r-release-osx-x86_64, r-release-windows-ix86+x86_64, r-oldrel-osx-x86_64, r-oldrel-windows-ix86+x86_64, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 0.3.0, 2.62, 2.12, 0, 0, 16, 2.72, 0, 1.9, 0, 30, 0, 5, 42.55, 33.16, 0, 0, 53, 42.57, 0, 42.1, 0, 73, 0, 54, 45.17, 35.28, 54.87, 55, 69, 45.29, 82.6, 44, 0, 103, 0, 59, OK, OK, OK, OK, OK, OK, OK, OK, OK, OK, OK, OK, https://www.R-project.org/nosvn/R.check/r-devel-linux-x86_64-debian-clang/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-devel-linux-x86_64-debian-gcc/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-devel-linux-x86_64-fedora-clang/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-devel-linux-x86_64-fedora-gcc/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-devel-windows-ix86+x86_64/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-patched-linux-x86_64/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-patched-solaris-x86/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-release-linux-x86_64/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-release-osx-x86_64/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-release-windows-ix86+x86_64/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-oldrel-osx-x86_64/localIV-00check.html, https://www.R-project.org/nosvn/R.check/r-oldrel-windows-ix86+x86_64/localIV-00check.html
#>   check_details
#> 1            NA
```

by name


```r
x <- cch_pkgs(c("geojsonio", "leaflet", "MASS"))
lapply(x, "[[", c("data", "summary"))
#> [[1]]
#> [[1]]$any
#> [1] TRUE
#> 
#> [[1]]$ok
#> [1] 11
#> 
#> [[1]]$note
#> [1] 1
#> 
#> [[1]]$warn
#> [1] 0
#> 
#> [[1]]$error
#> [1] 0
#> 
#> [[1]]$fail
#> [1] 0
#> 
#> 
#> [[2]]
#> [[2]]$any
#> [1] TRUE
#> 
#> [[2]]$ok
#> [1] 11
#> 
#> [[2]]$note
#> [1] 1
#> 
#> [[2]]$warn
#> [1] 0
#> 
#> [[2]]$error
#> [1] 0
#> 
#> [[2]]$fail
#> [1] 0
#> 
#> 
#> [[3]]
#> [[3]]$any
#> [1] FALSE
#> 
#> [[3]]$ok
#> [1] 12
#> 
#> [[3]]$note
#> [1] 0
#> 
#> [[3]]$warn
#> [1] 0
#> 
#> [[3]]$error
#> [1] 0
#> 
#> [[3]]$fail
#> [1] 0
```

historical data


```r
x <- cch_pkgs_history(x = "geojsonio")
x$data$history
#> # A tibble: 30 x 4
#>    date_updated summary$any   $ok $note $warn $error $fail checks
#>    <chr>        <lgl>       <int> <int> <int>  <int> <int> <list>
#>  1 2020-05-19T… TRUE           11     1     0      0     0 <df[,…
#>  2 2020-05-18T… TRUE           11     1     0      0     0 <df[,…
#>  3 2020-05-17T… TRUE           11     1     0      0     0 <df[,…
#>  4 2020-05-16T… TRUE           10     1     1      0     0 <df[,…
#>  5 2020-05-15T… TRUE           11     1     0      0     0 <df[,…
#>  6 2020-05-14T… TRUE           11     1     0      0     0 <df[,…
#>  7 2020-05-13T… TRUE           11     1     0      0     0 <df[,…
#>  8 2020-05-12T… TRUE           11     1     0      0     0 <df[,…
#>  9 2020-05-11T… TRUE           11     1     0      0     0 <df[,…
#> 10 2020-05-10T… TRUE           11     1     0      0     0 <df[,…
#> # … with 20 more rows, and 6 more variables: check_details$version <chr>,
#> #   $check <chr>, $result <chr>, $output <chr>, $flavors <list>,
#> #   $additional_issues <list>
```

search historical data


```r
x <- cch_pkgs_search(q = "memory")
x$data
#> # A tibble: 30 x 5
#>    package date_updated summary$any   $ok $note $warn $error $fail checks
#>    <chr>   <chr>        <lgl>       <int> <int> <int>  <int> <int> <list>
#>  1 tidyft  2020-04-20T… TRUE            9     0     2      0     0 <df[,…
#>  2 tidyft  2020-04-21T… TRUE           10     0     2      0     0 <df[,…
#>  3 tidyft  2020-04-22T… TRUE           10     0     2      0     0 <df[,…
#>  4 tidyft  2020-04-23T… TRUE           10     0     2      0     0 <df[,…
#>  5 tidyft  2020-04-24T… TRUE           10     0     2      0     0 <df[,…
#>  6 tidyft  2020-04-25T… TRUE           10     0     2      0     0 <df[,…
#>  7 tidyft  2020-04-26T… TRUE           10     0     2      0     0 <df[,…
#>  8 tidyft  2020-04-27T… TRUE           10     0     2      0     0 <df[,…
#>  9 tidyft  2020-04-28T… TRUE           10     0     2      0     0 <df[,…
#> 10 tidyft  2020-04-29T… TRUE           10     0     2      0     0 <df[,…
#> # … with 20 more rows, and 6 more variables: check_details$version <chr>,
#> #   $check <chr>, $result <chr>, $output <chr>, $flavors <list>,
#> #   $additional_issues <list>
```

## maintainers

all


```r
cch_maintainers(limit = 1)
#> $found
#> [1] 9450
#> 
#> $count
#> [1] 1
#> 
#> $offset
#> [1] 0
#> 
#> $error
#> NULL
#> 
#> $data
#>                    email             name
#> 1 f.briatte_at_gmail.com François Briatte
#>                                                                                url
#> 1 https://cloud.r-project.org/web/checks/check_results_f.briatte_at_gmail.com.html
#>               date_updated                       table
#> 1 2020-05-20T00:02:21.528Z ggnetwork, TRUE, 7, 5, 0, 0
#>                                                                                             packages
#> 1 ggnetwork, https://cloud.r-project.org/web/checks/check_results_ggnetwork.html, NOTE, OK, 5, 7, NA
```

by name


```r
cch_maintainers(c("maelle.salmon_at_yahoo.se", "13268259225_at_163.com"))
#> [[1]]
#> [[1]]$error
#> NULL
#> 
#> [[1]]$data
#> [[1]]$data$email
#> [1] "maelle.salmon_at_yahoo.se"
#> 
#> [[1]]$data$name
#> [1] "Maëlle Salmon"
#> 
#> [[1]]$data$url
#> [1] "https://cloud.r-project.org/web/checks/check_results_maelle.salmon_at_yahoo.se.html"
#> 
#> [[1]]$data$date_updated
#> [1] "2020-05-20T00:02:21.579Z"
#> 
#> [[1]]$data$table
#>       package   any ok note warn error
#> 1 monkeylearn  TRUE  7    5    0     0
#> 2    opencage  TRUE  6    6    0     0
#> 3        riem FALSE 12    0    0     0
#> 4     ropenaq FALSE 12    0    0     0
#> 5 rtimicropem  TRUE  7    5    0     0
#> 
#> [[1]]$data$packages
#>       package
#> 1 monkeylearn
#> 2    opencage
#> 3        riem
#> 4     ropenaq
#> 5 rtimicropem
#>                                                                     url
#> 1 https://cloud.r-project.org/web/checks/check_results_monkeylearn.html
#> 2    https://cloud.r-project.org/web/checks/check_results_opencage.html
#> 3        https://cloud.r-project.org/web/checks/check_results_riem.html
#> 4     https://cloud.r-project.org/web/checks/check_results_ropenaq.html
#> 5 https://cloud.r-project.org/web/checks/check_results_rtimicropem.html
#>     check_result version
#> 1 NOTE, OK, 5, 7      NA
#> 2 NOTE, OK, 6, 6      NA
#> 3         OK, 12      NA
#> 4         OK, 12      NA
#> 5 NOTE, OK, 5, 7      NA
#> 
#> 
#> 
#> [[2]]
#> [[2]]$error
#> NULL
#> 
#> [[2]]$data
#> [[2]]$data$email
#> [1] "13268259225_at_163.com"
#> 
#> [[2]]$data$name
#> [1] "Xiao-Ping You"
#> 
#> [[2]]$data$url
#> [1] "https://cloud.r-project.org/web/checks/check_results_13268259225_at_163.com.html"
#> 
#> [[2]]$data$date_updated
#> [1] "2020-05-20T00:02:22.034Z"
#> 
#> [[2]]$data$table
#>   package  any ok note warn error
#> 1    XHWE TRUE  0   12    0     0
#> 
#> [[2]]$data$packages
#>   package                                                            url
#> 1    XHWE https://cloud.r-project.org/web/checks/check_results_XHWE.html
#>   check_result version
#> 1         NULL      NA
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cch_maintainers.R
\name{cch_maintainers}
\alias{cch_maintainers}
\title{Get maintainer based checks}
\usage{
cch_maintainers(x = NULL, limit = 10, offset = 0, ...)
}
\arguments{
\item{x}{email slug name, optional, if you pass in more than one
we'll do async}

\item{limit}{number of records to return. Default: 10}

\item{offset}{record number to start at. Default: 0}

\item{...}{Curl options passed to \code{\link[crul:HttpClient]{crul::HttpClient()}} or
\code{\link[crul:Async]{crul::Async()}}}
}
\value{
list of info about a package(s)
}
\description{
Get maintainer based checks
}
\examples{
\dontrun{
x <- cch_maintainers()
x$data
x$data$summary
x$data$checks
x$data$date_updated
x$data$package
x$data$url

cch_maintainers("00gerhard_at_gmail.com")
cch_maintainers(c("123saga_at_gmail.com", "13268259225_at_163.com",
	 "csardi.gabor_at_gmail.com"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cch_pkgs.R
\name{cch_pkgs}
\alias{cch_pkgs}
\title{Get package checks data}
\usage{
cch_pkgs(x = NULL, limit = 10, offset = 0, ...)
}
\arguments{
\item{x}{package name, optional, if you pass in more than one
we'll do async}

\item{limit}{number of records to return. Default: 10}

\item{offset}{record number to start at. Default: 0}

\item{...}{Curl options passed to \code{\link[crul:HttpClient]{crul::HttpClient()}} or
\code{\link[crul:Async]{crul::Async()}}}
}
\value{
list of info about a package(s)
}
\description{
Get package checks data
}
\note{
this function only gets the current days checks; see
\code{\link[=cch_pkgs_history]{cch_pkgs_history()}} for historical data
}
\examples{
\dontrun{
x <- cch_pkgs()
x$data
x$data$summary
x$data$checks
x$data$date_updated
x$data$package
x$data$url

cch_pkgs("geojsonio")
cch_pkgs(c("geojsonio", "leaflet", "MASS"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cch_heartbeat.R
\name{cch_heartbeat}
\alias{cch_heartbeat}
\title{Get heartbeat}
\usage{
cch_heartbeat(...)
}
\arguments{
\item{...}{Curl options passed to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
list of routes
}
\description{
Get heartbeat
}
\examples{
\dontrun{
cch_heartbeat()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cch_history.R
\name{cch_history}
\alias{cch_history}
\title{Get historical check data for all packages by date}
\usage{
cch_history(date, ...)
}
\arguments{
\item{date}{(character) a date of the form \code{YYYY-MM-DD}. required}

\item{...}{Curl options passed to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a tibble with columns:
\itemize{
\item package: character vector of package names
\item summary: character vector of JSON hash's of check summary data
\item checks: character vector of JSON hash's of checks performed
\item check_details: character vector of check details. if no check
details the string will be "null"; if details given, then
a JSON hash of details
\item date_updated: character vector of dates, the date the check was
performed on
}
}
\description{
Get historical check data for all packages by date
}
\details{
This function gets historical data for all packages for a single
day; see \code{\link[=cch_pkgs_history]{cch_pkgs_history()}} for last 30 days history for particular
packages

You have to do a bit of data wrangling to get this data into a
easily sortable/filterable/etc. form
}
\examples{
\dontrun{
x <- cch_history(date = "2020-04-01")
str(x)
lapply(x$summary[1:3], jsonlite::fromJSON)
}
}
\seealso{
\code{\link[=cch_pkgs_history]{cch_pkgs_history()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cchecks-package.R
\docType{package}
\name{cchecks-package}
\alias{cchecks-package}
\alias{cchecks}
\title{cchecks}
\description{
Client for the cranchecks.info API
}
\section{Docs}{

See \url{https://cranchecks.info/docs}
}

\section{Historical data}{

There's an important shortcoming of historical data. The links in the
historical data in the \code{checks} field are not date specific. If you go to a
link in historical data, for example for April 2nd, 2020, links in that set
of data link to whatever the current check data is for that package. The
\code{check_details} field is date specific though; the text is scraped from the
package checks page each day and stored, so you can count on that to be date
specific. There are sometimes links to further checks, often of compiled
packages on various types of checks that CRAN runs; we do not have those
check results - we could get them but have not take the time to sort that
out.
}

\author{
Scott Chamberlain \email{sckott@protonmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cchn_rules.R, R/rule-add.R, R/rule-delete.R,
%   R/rule-get.R, R/rule-list.R
\name{cchn_rules}
\alias{cchn_rules}
\alias{cchn_rules_add}
\alias{cchn_pkg_rule_add}
\alias{cchn_rule_add}
\alias{cchn_pkg_rule_delete}
\alias{cchn_rule_delete}
\alias{cchn_pkg_rule_get}
\alias{cchn_rule_get}
\alias{cchn_pkg_rule_list}
\alias{cchn_rule_list}
\title{Notifications: add, list, get, delete notification rules}
\usage{
cchn_rules_add(rules, email, quiet = FALSE, ...)

cchn_pkg_rule_add(
  status = NULL,
  platform = NULL,
  time = NULL,
  regex = NULL,
  package = NULL,
  email = NULL,
  path = ".",
  quiet = FALSE,
  ...
)

cchn_rule_add(
  package,
  status = NULL,
  platform = NULL,
  time = NULL,
  regex = NULL,
  email = NULL,
  quiet = FALSE,
  ...
)

cchn_pkg_rule_delete(id, email = NULL, path = ".", ...)

cchn_rule_delete(id, email = NULL, quiet = FALSE, ...)

cchn_pkg_rule_get(id, email = NULL, path = ".", ...)

cchn_rule_get(id, email = NULL, quiet = FALSE, ...)

cchn_pkg_rule_list(email = NULL, path = ".", ...)

cchn_rule_list(email = NULL, quiet = FALSE, ...)
}
\arguments{
\item{rules}{(list) a list of rules. each element in the list must
be a named list, in which each must have a "package", and a set of
rules (see below)}

\item{email}{(character) email address to use for interaction with
the CRAN checks API. we use the email address in the maintainers slot
of the DESCRIPTION file of your working directory. you can supply an
email address instead}

\item{quiet}{(logical) suppress messages? default: \code{FALSE}}

\item{...}{Curl options passed to \link[crul:verb-GET]{crul::verb-GET}, \link[crul:verb-POST]{crul::verb-POST}, or
\link[crul:verb-DELETE]{crul::verb-DELETE}}

\item{status}{(character) a check status, one of: error, warn, note, fail}

\item{platform}{(character) a platform, a string to match against the
platform strings used by cran checks. e.g., "osx" would match any osx
platform check results, whereas you could limit the rule to just a single
specific platform by using the target platforms exact string
"r-oldrel-osx-x86_64". Leave as \code{NULL} (the default) to match all
platforms.}

\item{time}{(integer) number of days}

\item{regex}{(character) a regex string}

\item{package}{(character) a package name. if \code{NULL}, we attempt to
get the package name from the working directory, and fail out if there's
not a valid package structure/package name}

\item{path}{(character) path to a directory containing an R package}

\item{id}{(integer) a rule id. note that you can not get or delete
rules that are not yours. required}
}
\value{
\itemize{
\item \code{cchn_pkg_rule_add()}/\code{cchn_rule_add()}/\code{cchn_rules_add()}: message about
the rule added, and a note about using \code{cchn_rule_list()} to list your rules
\item \code{cchn_pkg_rule_get()}/\code{cchn_rule_get()}: list with elements \code{error} and
\code{data} (a list of the parts of the rule)
\item \code{cchn_pkg_rule_list()}/\code{cchn_rule_list()}: list with elements \code{error} and
\code{data} (a data.frame of all the rules associated with your email)
\item \code{cchn_pkg_rule_delete()}/\code{cchn_rule_delete()}: if deletion works, a
message saying "ok"
}
}
\description{
Notifications: add, list, get, delete notification rules
}
\details{
Functions prefixed with \code{cchn_pkg_} operate within a package
directory. That is, your current working directory is an R
package, and is the package for which you want to handle CRAN checks
notifications. These functions make sure that you are inside of
an R package, and use the email address and package name based
on the directory you're in.

Functions prefixed with just \code{cchn_} do not operate within a package.
These functions do not guess package name at all, but require the user
to supply a package name (for those functions that require a package name);
and instead of guessing an email address from your package, we guess email
from the cached cchecks email file.
\itemize{
\item \code{cchn_pkg_rule_add()}/\code{cchn_rule_add()}: add a rule, one rule per
function call
\item \code{cchn_rules_add()}: add many rules at once; no option for package context
\item \code{cchn_pkg_rule_get()}/\code{cchn_rule_get()}: get a rule by rule id (see
\code{cchn_pkg_rule_list()}/\code{cchn_rule_list()} to get ids; can only get rules for
the authenticated user)
\item \code{cchn_pkg_rule_list()}/\code{cchn_rule_list()}: list rules for the
authenticated user - \code{cchn_pkg_rule_list()} lists rules only for the package
in question, while \code{cchn_rule_list()} lists all rules for the user (email)
\item \code{cchn_pkg_rule_delete()}/\code{cchn_rule_delete()}: delete a rule by rule id
(only those for the authenticated user)
}
}
\section{example rules}{


Note that the first parameter \code{package} is left out for brevity
\itemize{
\item ERROR for at least 1 day across all platforms
\itemize{
\item \code{cchn_rule_add(status = 'error')}
}
\item ERROR for 3 days in a row across 2 or more platforms
\itemize{
\item \code{cchn_rule_add(status = 'error', time = 3, platform = 2)}
}
\item ERROR for 2 days in a row on all osx platforms
\itemize{
\item \code{cchn_rule_add(status = 'error', time = 2, platform = "osx")}
}
\item ERROR for 2 days in a row on all release R versions
\itemize{
\item \code{cchn_rule_add(status = 'error', time = 2, platform = "release")}
}
\item WARN for 4 days in a row on any platform except Solaris
\itemize{
\item \code{cchn_rule_add(status = 'warn', time = 4, platform = "-solaris")}
}
\item WARN for 2 days in a row across 9 or more platforms
\itemize{
\item \code{cchn_rule_add(status = 'warn', time = 2, platform = 10)}
}
\item NOTE on any osx platform
\itemize{
\item \code{cchn_rule_add(status = 'note', platform = "osx")}
}
\item NOTE on any platform
\itemize{
\item \code{cchn_rule_add(status = 'note')}
}
\item error details contain regex 'install'
\itemize{
\item \code{cchn_rule_add(regex = "install")}
}
}
}

\examples{
\dontrun{
## Workflow 1: within a package directory
# (x <- cchn_pkg_rule_list())
# if (length(x$data$id) > 0) {
#  cchn_pkg_rule_get(x$data$id[1])
#  cchn_pkg_rule_delete(x$data$id[1])
#  cchn_pkg_rule_get(id = x$data$id[1])
# }

## Workflow 2: not in a package directory
# (x <- cchn_rule_list())
# if (length(x$data$id) > 0) {
#  cchn_rule_get(x$data$id[1])
#  cchn_rule_delete(x$data$id[1])
#  cchn_rule_get(id = x$data$id[1])
# }

## cchn_pkg_rule_add: add a rule - guesses the package name
## you can specify the package name instead
# cchn_pkg_rule_add(status = "note", platform = 3,
#  email = "some-email")
## cchn_rule_add: add a rule - not in package context, must 
## specify the package name
# cchn_rule_add(package = "foobar", status = "note", platform = 3,
#  email = "some-email")

## cchn_pkg_rule_add: should guess package name and email
# cchn_pkg_rule_add(status = "note", platform = 3)

## cchn_rule_add: package name must be supplied. takes first email
## from cached emails.csv file, see `?cchn_register` for more
# cchn_rule_add(package = "foobar", status = "warn", platform = 2)

## cchn_rules_add: add many rules at once
## no package context here, email and package names must be given
# pkg <- "charlatan"
# rules <- list(
#   list(package = pkg, status = "warn"),
#   list(package = pkg, status = "error", time = 4)
# )
# cchn_rules_add(rules, "your-email", verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cch_pkgs_history.R
\name{cch_pkgs_history}
\alias{cch_pkgs_history}
\title{Get historical check data for packages}
\usage{
cch_pkgs_history(x, limit = 30, offset = 0, ...)
}
\arguments{
\item{x}{package name, required, if you pass in more than one
we'll do async}

\item{limit}{number of records to return. Default: 10}

\item{offset}{record number to start at. Default: 0}

\item{...}{Curl options passed to \code{\link[crul:HttpClient]{crul::HttpClient()}} or
\code{\link[crul:Async]{crul::Async()}}}
}
\value{
list of info about a package(s)
}
\description{
Get historical check data for packages
}
\details{
this function gets historical data; for current day check data only
see \code{\link[=cch_pkgs]{cch_pkgs()}}

data is only available for 30 days prior to today's date, see the
\code{\link[=cch_history]{cch_history()}} function for older data
}
\examples{
\dontrun{
x <- cch_pkgs_history(x = "geojsonio")
x
x$data
x$data$package
x$data$history
x$data$history$summary
x$data$history$summary$any
x$data$history$check_details

# many packages
res <- cch_pkgs_history(c("geojsonio", "leaflet", "MASS"))
res

# pagination
cch_pkgs_history(x = "geojsonio", limit = 3)
cch_pkgs_history(x = "geojsonio", limit = 3, offset = 4)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cchn_register.R
\name{cchn_register}
\alias{cchn_register}
\title{Notifications: register your email address and get a token}
\usage{
cchn_register(email = NULL, token = NULL, ...)
}
\arguments{
\item{email}{(character) email address to use for interaction with
the CRAN checks API. If no email address is given, we go through a few
steps: check for the cached file mentioned below for any emails; use
\code{\link[whoami:email_address]{whoami::email_address()}}; if the location you are running this in is
a package, we look for the maintainer's email address in the package.}

\item{token}{(character) your CRAN checks API token. you shouldn't need
to pass a token here. if you used \code{\link[=cchn_register]{cchn_register()}} your token should
be cached}

\item{...}{Curl options passed to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
\code{NULL} - nothing returned
}
\description{
Notifications: register your email address and get a token
}
\details{
We cache a file with email addresses and tokens at the path
\code{file.path(rappdirs::user_data_dir("cranchecks", "cchecks"), "emails.csv")}
You can run that in R to get the path for the file on your machine.

To get a new token for an email address that was previously registered,
go to the file above and delete the line with the email address and token
for the email address in question; remember to save the change.
Then when you run \code{cchn_register()} again for that email you can get
a new token.

To add an email address that was validated before (probably on
another machine), to the configuration file, call this function
with the ‘email’ and ‘token’ arguments.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cch_pkgs_search.R
\name{cch_pkgs_search}
\alias{cch_pkgs_search}
\title{Search historical package check data}
\usage{
cch_pkgs_search(
  q,
  package = NULL,
  one_each = FALSE,
  fields = NULL,
  limit = 30,
  offset = 0,
  ...
)
}
\arguments{
\item{q}{(character) full text query string}

\item{package}{(character) a package name. limit results to a single
package, e.g, \code{package="taxize"}}

\item{one_each}{(logical) if \code{TRUE}, return a single result for each package;
useful if you want to find out what packages match a particular query, and
don't care which day that match happened. default: \code{FALSE}}

\item{fields}{(character) vector of fields to return, e.g.,
\code{fields=c("package", "check_details")}}

\item{limit}{number of records to return. Default: 10}

\item{offset}{record number to start at. Default: 0}

\item{...}{Curl options passed to \code{\link[crul:HttpClient]{crul::HttpClient()}} or
\code{\link[crul:Async]{crul::Async()}}}
}
\value{
list of info about a package(s)
}
\description{
Search historical package check data
}
\examples{
\dontrun{
x <- cch_pkgs_search(q = "memory")
x
x$data
x$data$package
x$data
x$data$summary
x$data$summary$any
x$data$check_details

# restrict returned fields
res <- cch_pkgs_search("memory", fields = c("package", "check_details"))
res$data$check_details$output[1]
grepl('memory', res$data$check_details$output[1])

# one each, one record per package
res <- cch_pkgs_search("memory", one_each = TRUE,
  fields = c("package", "date_updated"))
res

# pagination
cch_pkgs_search("memory", limit = 3)
cch_pkgs_search("memory", limit = 3, offset = 4)
}
}
