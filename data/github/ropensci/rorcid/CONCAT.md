rorcid
======



<!-- badges: start -->
[![R build status](https://github.com/ropensci/rorcid/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rorcid/actions)
[![cran checks](https://cranchecks.info/badges/worst/rorcid)](https://cranchecks.info/pkgs/rorcid)
[![codecov.io](https://codecov.io/github/ropensci/rorcid/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rorcid?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rorcid?color=2ED968)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rorcid)](https://cran.r-project.org/package=rorcid)
<!-- badges: end -->

`rorcid` is an R programmatic interface to the Orcid public API. `rorcid` is not a product developed or distributed by ORCID®.

rorcid docs: <https://docs.ropensci.org/rorcid/>

Orcid API docs:

* Public API docs: <https://members.orcid.org/api>
* Swagger docs: https://pub.orcid.org/v3.0/#/Development_Public_API_v3.0

The package now works with the `v3.0` ORCID API. It's too complicated to allow users to work with different versions of the API, so it's hard-coded to `v3.0`.

## Authentication

There are three ways to authenticate with `rorcid`:

- Interactively login with OAuth. This doesn't require any input on
your part. We use a client id and client secret key to ping ORCID.org;
at which point you log in with your username/password; then we get back
a token (same as the above option). We don't know your username or
password, only the token that we get back. We cache that token locally
in a hidden file in whatever working directory you're in. If you delete
that file, or run the code from a new working directory, then we
re-authorize.
- Use a `client_id` and `client_secret` to do 2-legged OAuth. 
ORCID docs at https://members.orcid.org/api/oauth/2legged-oauth and
https://members.orcid.org/api/post-oauthtoken-reading-public-data
This requires you to register a "client application". See 
https://orcid.org/content/register-client-application-2 for instructions
- Use a token as a result of either of the two above approaches. The token
is a alphanumeric UUID, e.g. `dc0a6b6b-b4d4-4276-bc89-78c1e9ede56e`. You
can get this token by running `orcid_auth()`, then storing that key
(the uuid alone, not the "Bearer " part) either as en environment
variable in your `.Renviron` file in your home directory (with the name
`ORCID_TOKEN`), or as an R option in your `.Rprofile` file (with the name
`orcid_token`). See [Startup] for more information.
Either an environment variable or R option work. If we don't find
either we do the next option.

We recommend the 3rd option if possible, specifically, storing the token
as an environment variable permanently.

If authentication fails, you can still use `rorcid`. ORCID does not require
authentication at this point, but may in the future - this prepares you
for when that happens.

See <https://info.orcid.org/documentation/integration-guide/getting-started-with-your-orcid-integration/#easy-faq-2569> for more about ORCID 
OAuth Scopes.

## Computing environments without browsers

One pitfall is when you are using `rorcid` on a server, and you're ssh'ed
in, so that there's no way to open a browser to do the OAuth browser
flow. Similarly for any other situation in which a browser can not be
opened. In this case, run `orcid_auth()` on another machine in which you do
have the ability to open a browser, then collect the info that's ouptput
from `orcid_auth()` and store it as an environment variable (see above).

## Installation

Stable version


```r
install.packages("rorcid")
```

Development version


```r
remotes::install_github("ropensci/rorcid")
```


```r
library('rorcid')
```

## Docs

Get started with rorcid at <https://docs.ropensci.org/rorcid/>

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rorcid/issues)
* License: MIT
* Get citation information for `rorcid` in R doing `citation(package = 'rorcid')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rorcid 0.7.0
============

### MINOR IMPROVEMENTS

* `orcid_search()` changes: now returns number of records found as an attribute; get it by doing `attr(x, "found")` if `x` is the result of a call to `orcid_search()`. In addition, added an example of retrieving number of records found, as well as added more documentation to the function page on pagination. (#86)
* Auth gains a 3rd method: 2 legged OAuth. Updated `?orcid_auth` docs on the 3 different auth methods, as well as the README and auth vignette. (#87)
* `orcid_citations()` produced wrongly formatted BibTeX strings under particular circumstances; fixed. (#92) 


rorcid 0.6.4
============

### MINOR IMPROVEMENTS

* removed `recursive` parameter in `orcid()` function as it wasn't used internally (#65)
* fix doc title for `orcid_external_identifiers()` (#81)
* in `orcid_search()` removed parameters `current_prim_inst` and `patent_number` as those have been removed from ORCID; added additional fields to query on in the `orcid()` function: `peer-review-type`, `peer-review-role`, `peer-review-group-id`, `biography`, and `external-id-type-and-value` (#82)

### BUG FIXES 

* remove unused param definition for `...` in `works()` - via R-devel checks  (#84)


rorcid 0.6.0
============

### NEW FEATURES

* `orcid_citations()` gains new internal function `extract_bibtex` to attempt to better parse bibtex (#74) thanks @RLumSK
* `orcid_search()` gains new parameter `affiliation_org` to search by affiliation organization name, and parses `affiliation_org` output (#75)
* `orcid_search()` gains new parameters `ringgold_org_id` and `grid_org_id`, as well as parsing those outputs if present (#76)

### MINOR IMPROVEMENTS

* clarify names of variables to store the ORCID API token (#71) thanks @fmichonneau
* dont run tests if ORCID token not found (#79)

### BUG FIXES 

* fix to `orcid_search()`: parsing issues were fixed by giving default NA's when no results found for certain result sections (#72)
* tests now using cached requests/responses via vcr (#73)
* `orcid_search()` fix: `orcid-identifier.path` was not returned when no results found; now we check for that in the results, and return an empty tibble if not found (#75)
* in `orcid_search()`, keywords have to be passed as multiple instances of `keyword` rather than as `keywords`; fixed now; user facing still just uses `keywords` as the input param (#77)

rorcid 0.5.2
============

### DEPLOYMENT

* Switched from using TravisCI and Appveyor to using only GitHub Actions for
  _R_.

rorcid 0.5.0
============

`rorcid` now works with the v3 ORCID API (#63) (#68) (#70)

### NEW FEATURES

* with ORCID v3 API, new functions added: `orcid_distinctions()`, `orcid_invited_positions()`, `orcid_memberships()`, `orcid_qualifications()`, `orcid_research_resources()`, and `orcid_services()`
* gains new fxn `orcid_citations()` for getting citations for an ORCID ID in user specified formats - leverages `rcrossref` and `handlr` packages (#51) (#69)
* new function added `orcid_search()`, a wrapper around `orcid()` function as an easier interface than `orcid()` - see ropensci/codemetar issue 83 for discussion (#54)
* requires R 3.5 or greater

### MINOR IMPROVEMENTS

* dataset added to the package `issn_title`, a named vector, with values as journal names and names as their ISSN values (sourced from Crossref). see `?orcid_peer_reviews` examples for an example of using the dataset to gather journal titles from jorunal ISSN's (#52)
* added documentation on ORCID authentication to README (#60) thanks @maelle
* use `fauxpas::find_error_class` method instead of internal hack (#61)
* Added more examples to the vignette (#56) thanks @bomeara
* fix many typos (#59) thanks @maelle 
* add section to `?orcid_auth` documentation about "Computing evironments without browsers" - you can't do OAuth flow in a non-interactive session (#55) thanks @pkraker for the find
* changes for `orcid_works()`: `put_code` parameter now accepts up to 50 put codes; significant changes internally to make it easier to combine results into a data.frame  (#44) thanks @gorkang

### BUG FIXES 

* `httpuv` package added to Suggests and used inside only the `orcid_auth()` function when doing the OAuth flow because out of band (OOB) OAuth doesn't work without httpuv (#67) thanks @ciakovx for finding that
* fix to `identifiers()` function - was failing on results that gave zero length lists (#40) thanks @agbarnett


rorcid 0.4.0
============

Most changes in this version are to update the package to work with the new ORCID API (`v2.1`). (#37) (#40)

### NEW FEATURES

* `rorcid` now support OAuth authentication. We still recommend to not use OAuth, but to get a token and store that as an environment variable. See `?orcid_auth` for help (#26)
* To work with the new ORCID API, we've introduced new functions: `orcid_activities()`, `orcid_address()`, `orcid_auth()`, `orcid_bio()`, `orcid_educations()`, `orcid_email()`, `orcid_employments()`, `orcid_external_identifiers()`, `orcid_fundings()`, `orcid_keywords()`, `orcid_other_names()`, `orcid_peer_reviews()`, `orcid_person()`, `orcid_ping()`, `orcid_researcher_urls()`, `orcid_works()`
* we've introduced a new packcage import `data.table` for binding lists together into a data.frame

### MINOR IMPROVEMENTS

* Fixes to `identifiers()` for new API. Includes better failure behavior on classes it doesn't support (#34) (#39)
* No changes here other than adding some examples, but ORCID API now allows ot search on some works metadata (e.g., titles) (#33)
* Examples added for searching on specific fields (#36)
* `orcid_id()` changed internally; now wraps the new function `orcid_person()` (#41)
* using Markdown docs now (#35)
* replaced `httr` with `crul` for HTTP requests. we have retained `httr` only to do OAuth (#32)
* `orcid_id()` loses its `profile` parameter due to the ORCID API change. it does pass on parameters to `orcid_person()`, so see that man file
* `works()` now returns a tibble/data.frame instead of a list of items
* A much updated package level man file with lots of docs


rorcid 0.3.0
============

### NEW FEATURES

* Added a vignette (#20)
* `orcid_id()` function gains output for employment and funding (#24) (#29)

### MINOR IMPROVEMENTS

* change all `is()` calls to `inherits()` (#30)
* using `tibble` package now for compact data.frame outputs instead
of internal code. an associated change in the output of both `orcid()`
and `orcid_doi()` is that we now return a tibble (data.frame) instead of
a data.frame as a slot in a list. we add how many results are returned from 
your search as an attribute on the data.frame. Access it like 
`attr(out, "found")` (#25)
* base ORCID API URL changed from `http` to `https` scheme
* genereal improvements to documentation throughout package

### DEPRECATED AND DEFUNCT

* `summary.or_id()` is now defunct. see `?rorcid-defunct`


rorcid 0.2.2
============

### NEW FEATURES

* Require `httr >= v1.1.0` (#23)

### MINOR IMPROVEMENTS

* Updated `dplyr` tidy data.frame internal code (#21)
* Changes to internal use of `httr::content()` to parse to text, then read JSON
manually using `jsonlite` & to always set `encoding` explicitly in the same calls (#22)

### BUG FIXES

* Fix to `as.orcid()` and presumably other function calls by requiring
`httr >= v1.1.0` because older versions cause a problem when parsing
responses (#23) thanks @ericwatt


rorcid 0.2.0
============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.3 Patched
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 reverse dependency. (Summary at <https://github.com/ropensci/rorcid/blob/master/revdep/README.md>). No problems were found.

---

This version makes minor improvements, including improved documentation on authentication.

Thanks!
Scott Chamberlain
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rorcid/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rorcid.git`
* Make sure to track progress upstream (i.e., on our version of `rorcid` at `ropensci/rorcid`) by doing `git remote add upstream https://github.com/ropensci/rorcid.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rorcid`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Email

Do not email the maintainer - open an issue instead.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>

# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.3 Patched (2021-01-08 r79819) |
|os       |macOS Catalina 10.15.7                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-01-20                                  |

# Dependencies

|package |old   |new      |Δ  |
|:-------|:-----|:--------|:--|
|rorcid  |0.6.4 |0.6.4.99 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*rorcid
======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

<!-- badges: start -->
[![R build status](https://github.com/ropensci/rorcid/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rorcid/actions)
[![cran checks](https://cranchecks.info/badges/worst/rorcid)](https://cranchecks.info/pkgs/rorcid)
[![codecov.io](https://codecov.io/github/ropensci/rorcid/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rorcid?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rorcid?color=2ED968)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rorcid)](https://cran.r-project.org/package=rorcid)
<!-- badges: end -->

`rorcid` is an R programmatic interface to the Orcid public API. `rorcid` is not a product developed or distributed by ORCID®.

rorcid docs: <https://docs.ropensci.org/rorcid/>

Orcid API docs:

* Public API docs: <https://members.orcid.org/api>
* Swagger docs: https://pub.orcid.org/v3.0/#/Development_Public_API_v3.0

The package now works with the `v3.0` ORCID API. It's too complicated to allow users to work with different versions of the API, so it's hard-coded to `v3.0`.

## Authentication

There are three ways to authenticate with `rorcid`:

- Interactively login with OAuth. This doesn't require any input on
your part. We use a client id and client secret key to ping ORCID.org;
at which point you log in with your username/password; then we get back
a token (same as the above option). We don't know your username or
password, only the token that we get back. We cache that token locally
in a hidden file in whatever working directory you're in. If you delete
that file, or run the code from a new working directory, then we
re-authorize.
- Use a `client_id` and `client_secret` to do 2-legged OAuth. 
ORCID docs at https://members.orcid.org/api/oauth/2legged-oauth and
https://members.orcid.org/api/post-oauthtoken-reading-public-data
This requires you to register a "client application". See 
https://orcid.org/content/register-client-application-2 for instructions
- Use a token as a result of either of the two above approaches. The token
is a alphanumeric UUID, e.g. `dc0a6b6b-b4d4-4276-bc89-78c1e9ede56e`. You
can get this token by running `orcid_auth()`, then storing that key
(the uuid alone, not the "Bearer " part) either as en environment
variable in your `.Renviron` file in your home directory (with the name
`ORCID_TOKEN`), or as an R option in your `.Rprofile` file (with the name
`orcid_token`). See [Startup] for more information.
Either an environment variable or R option work. If we don't find
either we do the next option.

We recommend the 3rd option if possible, specifically, storing the token
as an environment variable permanently.

If authentication fails, you can still use `rorcid`. ORCID does not require
authentication at this point, but may in the future - this prepares you
for when that happens.

See <https://info.orcid.org/documentation/integration-guide/getting-started-with-your-orcid-integration/#easy-faq-2569> for more about ORCID 
OAuth Scopes.

## Computing environments without browsers

One pitfall is when you are using `rorcid` on a server, and you're ssh'ed
in, so that there's no way to open a browser to do the OAuth browser
flow. Similarly for any other situation in which a browser can not be
opened. In this case, run `orcid_auth()` on another machine in which you do
have the ability to open a browser, then collect the info that's ouptput
from `orcid_auth()` and store it as an environment variable (see above).

## Installation

Stable version

```{r eval=FALSE}
install.packages("rorcid")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/rorcid")
```

```{r}
library('rorcid')
```

## Docs

Get started with rorcid at <https://docs.ropensci.org/rorcid/>

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rorcid/issues)
* License: MIT
* Get citation information for `rorcid` in R doing `citation(package = 'rorcid')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: rorcid authentication
author: Scott Chamberlain
date: "2021-01-20"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{rorcid authentication}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

## Authentication

There are three ways to authenticate with `rorcid`:

- Interactively login with OAuth. This doesn't require any input on
your part. We use a client id and client secret key to ping ORCID.org;
at which point you log in with your username/password; then we get back
a token (same as the above option). We don't know your username or
password, only the token that we get back. We cache that token locally
in a hidden file in whatever working directory you're in. If you delete
that file, or run the code from a new working directory, then we
re-authorize.
- Use a `client_id` and `client_secret` to do 2-legged OAuth. 
ORCID docs at https://members.orcid.org/api/oauth/2legged-oauth and
https://members.orcid.org/api/post-oauthtoken-reading-public-data
This requires you to register a "client application". See 
https://orcid.org/content/register-client-application-2 for instructions
- Use a token as a result of either of the two above approaches. The token
is a alphanumeric UUID, e.g. `dc0a6b6b-b4d4-4276-bc89-78c1e9ede56e`. You
can get this token by running `orcid_auth()`, then storing that key
(the uuid alone, not the "Bearer " part) either as en environment
variable in your `.Renviron` file in your home directory (with the name
`ORCID_TOKEN`), or as an R option in your `.Rprofile` file (with the name
`orcid_token`). See [Startup] for more information.
Either an environment variable or R option work. If we don't find
either we do the next option.

We recommend the 3rd option if possible, specifically, storing the token
as an environment variable permanently.

If authentication fails, you can still use `rorcid`. ORCID does not require
authentication at this point, but may in the future - this prepares you
for when that happens.

See <https://info.orcid.org/documentation/integration-guide/getting-started-with-your-orcid-integration/#easy-faq-2569> for more about ORCID 
OAuth Scopes.

## Computing environments without browsers

One pitfall is when you are using `rorcid` on a server, and you're ssh'ed
in, so that there's no way to open a browser to do the OAuth browser
flow. Similarly for any other situation in which a browser can not be
opened. In this case, run `orcid_auth()` on another machine in which you do
have the ability to open a browser, then collect the info that's ouptput
from `orcid_auth()` and store it as an environment variable (see above).
---
title: rorcid introduction
author: Scott Chamberlain
date: "2021-01-20"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{rorcid introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



`rorcid` is an R programmatic interface to the Orcid public API. `rorcid` is not a product developed or distributed by ORCID®.

Orcid API docs: http://members.orcid.org/api

The package now works with the `v3.0` ORCID API. It's too complicated to allow users to work with different versions of the API, so it's hard-coded to `v3.0`.


## Load rorcid


```r
library("rorcid")
```

## as.orcid

There's a function `as.orcid()` in this package to help coerce an Orcid ID to an object that holds details for that Orcid ID, prints a nice summary, and you can browse easily to that profile. E.g.


```r
as.orcid(x = "0000-0002-1642-628X")
#> <ORCID> 0000-0002-1642-628X
#>   Name: Boettiger, Carl
#>   URL (first): http://www.carlboettiger.info
#>   Country: US
#>   Keywords: Ecology, Evolution, Regime Shifts, Stochastic Dynamics
```

Or you can pass in many IDs


```r
as.orcid(c("0000-0003-1620-1408", "0000-0002-9341-7985"))
#> [[1]]
#> <ORCID> 0000-0003-1620-1408
#>   Name: Johnson, Thomas
#>   URL (first): 
#>   Country: US
#>   Keywords: 
#> 
#> [[2]]
#> <ORCID> 0000-0002-9341-7985
#>   Name: Binfield, Peter
#>   URL (first): 
#>   Country: US
#>   Keywords:
```

The `browse()` function lets you browser to a profile easily with a single function call


```r
browse(as.orcid("0000-0002-1642-628X"))
```

![profile](http://f.cl.ly/items/3d1o0k1X3R1U110C0u1u/Screen%20Shot%202015-02-10%20at%207.21.57%20PM.png)

## Get works

The `works()` function helps get works data from an orcid data object. The output of `works()` uses a print method to just print citations for each work.


```r
(out <- works(orcid_id("0000-0002-0233-1757")))
#> # A tibble: 6 x 32
#>   `put-code` url   type  `journal-title` visibility path  `display-index`
#> *      <int> <lgl> <chr> <lgl>           <chr>      <chr> <chr>          
#> 1    5296064 NA    jour… NA              public     /000… 0              
#> 2    5296065 NA    jour… NA              public     /000… 0              
#> 3    5296066 NA    jour… NA              public     /000… 0              
#> 4    9012984 NA    jour… NA              public     /000… 0              
#> 5    9012985 NA    jour… NA              public     /000… 0              
#> 6    9012986 NA    jour… NA              public     /000… 0              
#> # … with 25 more variables: `created-date.value` <dbl>,
#> #   `last-modified-date.value` <dbl>, `source.source-client-id` <lgl>,
#> #   `source.assertion-origin-orcid` <lgl>,
#> #   `source.assertion-origin-client-id` <lgl>,
#> #   `source.assertion-origin-name` <lgl>, `source.source-orcid.uri` <chr>,
#> #   `source.source-orcid.path` <chr>, `source.source-orcid.host` <chr>,
#> #   `source.source-name.value` <chr>, title.subtitle <lgl>,
#> #   `title.translated-title` <lgl>, title.title.value <chr>,
#> #   `external-ids.external-id` <list>, `publication-date.year.value` <chr>,
#> #   `publication-date.month.value` <chr>, `publication-date.day.value` <chr>,
#> #   `publication-date.day` <lgl>, `external-ids` <lgl>,
#> #   `publication-date` <lgl>, `source.source-orcid` <lgl>,
#> #   `source.source-client-id.uri` <chr>, `source.source-client-id.path` <chr>,
#> #   `source.source-client-id.host` <chr>, url.value <chr>
```

## Search Orcid

Get a list of names and Orcid IDs matching a name query


```r
orcid(query = "carl boettiger")
#> # A tibble: 1,000 x 3
#>    `orcid-identifier.uri`            `orcid-identifier.pa… `orcid-identifier.ho…
#>  * <chr>                             <chr>                 <chr>                
#>  1 https://orcid.org/0000-0002-1642… 0000-0002-1642-628X   orcid.org            
#>  2 https://orcid.org/0000-0002-5951… 0000-0002-5951-4503   orcid.org            
#>  3 https://orcid.org/0000-0002-3554… 0000-0002-3554-5196   orcid.org            
#>  4 https://orcid.org/0000-0003-1853… 0000-0003-1853-1574   orcid.org            
#>  5 https://orcid.org/0000-0002-7462… 0000-0002-7462-1956   orcid.org            
#>  6 https://orcid.org/0000-0002-7790… 0000-0002-7790-5102   orcid.org            
#>  7 https://orcid.org/0000-0003-1021… 0000-0003-1021-5374   orcid.org            
#>  8 https://orcid.org/0000-0001-8899… 0000-0001-8899-7850   orcid.org            
#>  9 https://orcid.org/0000-0001-7084… 0000-0001-7084-5402   orcid.org            
#> 10 https://orcid.org/0000-0003-0241… 0000-0003-0241-3557   orcid.org            
#> # … with 990 more rows
```

You can string together many search terms


```r
orcid(query = "johnson cardiology houston")
#> # A tibble: 1,000 x 3
#>    `orcid-identifier.uri`            `orcid-identifier.pa… `orcid-identifier.ho…
#>  * <chr>                             <chr>                 <chr>                
#>  1 https://orcid.org/0000-0002-0897… 0000-0002-0897-2301   orcid.org            
#>  2 https://orcid.org/0000-0002-5281… 0000-0002-5281-4466   orcid.org            
#>  3 https://orcid.org/0000-0001-8188… 0000-0001-8188-0078   orcid.org            
#>  4 https://orcid.org/0000-0002-2862… 0000-0002-2862-4160   orcid.org            
#>  5 https://orcid.org/0000-0002-0682… 0000-0002-0682-9982   orcid.org            
#>  6 https://orcid.org/0000-0002-4334… 0000-0002-4334-4001   orcid.org            
#>  7 https://orcid.org/0000-0002-6401… 0000-0002-6401-4597   orcid.org            
#>  8 https://orcid.org/0000-0003-4792… 0000-0003-4792-0149   orcid.org            
#>  9 https://orcid.org/0000-0001-9667… 0000-0001-9667-1615   orcid.org            
#> 10 https://orcid.org/0000-0003-0945… 0000-0003-0945-6138   orcid.org            
#> # … with 990 more rows
```

And you can use start and rows arguments to do pagination


```r
orcid("johnson cardiology houston", start = 2, rows = 3)
#> # A tibble: 3 x 3
#>   `orcid-identifier.uri`            `orcid-identifier.pat… `orcid-identifier.ho…
#> * <chr>                             <chr>                  <chr>                
#> 1 https://orcid.org/0000-0001-8188… 0000-0001-8188-0078    orcid.org            
#> 2 https://orcid.org/0000-0002-2862… 0000-0002-2862-4160    orcid.org            
#> 3 https://orcid.org/0000-0002-0682… 0000-0002-0682-9982    orcid.org
```

## Search by Orcid ID


```r
out <- orcid_id(orcid = "0000-0002-9341-7985")
out$`0000-0002-9341-7985`$addresses
#> $`last-modified-date`
#> $`last-modified-date`$value
#> [1] 1.465227e+12
#> 
#> 
#> $address
#>   visibility                               path put-code display-index
#> 1     public /0000-0002-9341-7985/address/89515    89515             0
#>   created-date.value last-modified-date.value source.source-client-id
#> 1       1.453659e+12             1.465227e+12                      NA
#>   source.assertion-origin-orcid source.assertion-origin-client-id
#> 1                            NA                                NA
#>   source.assertion-origin-name               source.source-orcid.uri
#> 1                           NA https://orcid.org/0000-0002-9341-7985
#>   source.source-orcid.path source.source-orcid.host source.source-name.value
#> 1      0000-0002-9341-7985                orcid.org           Peter Binfield
#>   country.value
#> 1            US
#> 
#> $path
#> [1] "/0000-0002-9341-7985/address"
```


## Search by DOIs

There is a helper function `check_dois()` that uses a regex checker to see if your DOIs are likely good or likely bad:

All good DOIs


```r
dois <- c("10.1371/journal.pone.0025995","10.1371/journal.pone.0053712",
       "10.1371/journal.pone.0054608","10.1371/journal.pone.0055937")
check_dois(dois)
#> $good
#> [1] "10.1371/journal.pone.0025995" "10.1371/journal.pone.0053712"
#> [3] "10.1371/journal.pone.0054608" "10.1371/journal.pone.0055937"
#> 
#> $bad
#> NULL
```

Some good, some bad


```r
dois <- c("10.1016/j.medpal.2008.12.005","10.1080/00933104.2000.10505926","10.1037/a0024480",
        "10.1002/anie.196603172","2344","asdf","232","asdf","23dd")
check_dois(dois)
#> $good
#> [1] "10.1016/j.medpal.2008.12.005"   "10.1080/00933104.2000.10505926"
#> [3] "10.1037/a0024480"               "10.1002/anie.196603172"        
#> 
#> $bad
#> [1] "2344" "asdf" "232"  "asdf" "23dd"
```

Basic search


```r
orcid_doi(dois = "10.1087/20120404")
#> [[1]]
#> # A tibble: 13 x 3
#>    `orcid-identifier.uri`            `orcid-identifier.pa… `orcid-identifier.ho…
#>  * <chr>                             <chr>                 <chr>                
#>  1 https://orcid.org/0000-0002-3181… 0000-0002-3181-3456   orcid.org            
#>  2 https://orcid.org/0000-0002-1290… 0000-0002-1290-9735   orcid.org            
#>  3 https://orcid.org/0000-0001-7343… 0000-0001-7343-9784   orcid.org            
#>  4 https://orcid.org/0000-0003-3738… 0000-0003-3738-1487   orcid.org            
#>  5 https://orcid.org/0000-0001-5727… 0000-0001-5727-2427   orcid.org            
#>  6 https://orcid.org/0000-0003-1603… 0000-0003-1603-8743   orcid.org            
#>  7 https://orcid.org/0000-0002-2123… 0000-0002-2123-6317   orcid.org            
#>  8 https://orcid.org/0000-0002-0860… 0000-0002-0860-8606   orcid.org            
#>  9 https://orcid.org/0000-0003-0493… 0000-0003-0493-0128   orcid.org            
#> 10 https://orcid.org/0000-0003-3188… 0000-0003-3188-6273   orcid.org            
#> 11 https://orcid.org/0000-0002-5993… 0000-0002-5993-8592   orcid.org            
#> 12 https://orcid.org/0000-0003-1419… 0000-0003-1419-2405   orcid.org            
#> 13 https://orcid.org/0000-0001-5109… 0000-0001-5109-3700   orcid.org            
#> 
#> attr(,"class")
#> [1] "orcid_doi"
```

This DOI is not a real one, but a partial DOI, then we can fuzzy search


```r
orcid_doi(dois = "10.1087/2", fuzzy = TRUE, rows = 5)
#> [[1]]
#> # A tibble: 5 x 3
#>   `orcid-identifier.uri`            `orcid-identifier.pat… `orcid-identifier.ho…
#> * <chr>                             <chr>                  <chr>                
#> 1 https://orcid.org/0000-0001-8701… 0000-0001-8701-615X    orcid.org            
#> 2 https://orcid.org/0000-0001-6081… 0000-0001-6081-0708    orcid.org            
#> 3 https://orcid.org/0000-0001-5919… 0000-0001-5919-8670    orcid.org            
#> 4 https://orcid.org/0000-0002-0634… 0000-0002-0634-3376    orcid.org            
#> 5 https://orcid.org/0000-0001-8959… 0000-0001-8959-3380    orcid.org            
#> 
#> attr(,"class")
#> [1] "orcid_doi"
```

Function is vectorized, search for many DOIs


```r
dois <- c("10.1371/journal.pone.0025995","10.1371/journal.pone.0053712",
       "10.1371/journal.pone.0054608","10.1371/journal.pone.0055937")
res <- orcid_doi(dois = dois, fuzzy = TRUE)
res[[1]]
#> # A tibble: 1,000 x 3
#>    `orcid-identifier.uri`            `orcid-identifier.pa… `orcid-identifier.ho…
#>  * <chr>                             <chr>                 <chr>                
#>  1 https://orcid.org/0000-0002-0051… 0000-0002-0051-973X   orcid.org            
#>  2 https://orcid.org/0000-0001-7767… 0000-0001-7767-3500   orcid.org            
#>  3 https://orcid.org/0000-0001-6577… 0000-0001-6577-5526   orcid.org            
#>  4 https://orcid.org/0000-0001-9637… 0000-0001-9637-4127   orcid.org            
#>  5 https://orcid.org/0000-0002-3461… 0000-0002-3461-9990   orcid.org            
#>  6 https://orcid.org/0000-0002-0153… 0000-0002-0153-9708   orcid.org            
#>  7 https://orcid.org/0000-0003-1403… 0000-0003-1403-7093   orcid.org            
#>  8 https://orcid.org/0000-0002-5591… 0000-0002-5591-5721   orcid.org            
#>  9 https://orcid.org/0000-0003-0551… 0000-0003-0551-2550   orcid.org            
#> 10 https://orcid.org/0000-0003-4410… 0000-0003-4410-6539   orcid.org            
#> # … with 990 more rows
```

## Get formatted citations for an ORCID ID

One workflow is to get publications associated with an ORCID profile. The following will extract all the works with a DOI, then use the `rcrossref` package to get nicely formatted references for each, and then export them to a bibtex file


```r
my_dois <- rorcid::identifiers(rorcid::works("0000-0002-0337-5997"))
pubs <- rcrossref::cr_cn(dois = my_dois, format = "bibtex")
invisible(lapply(pubs, write, "pubs.bib", append=TRUE))
```


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_person.R
\name{orcid_person}
\alias{orcid_person}
\title{Get personal data for a person}
\usage{
orcid_person(orcid, details = FALSE, ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{details}{(logical). also get details. Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get personal data for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_person(orcid = "0000-0002-9341-7985")
res$`0000-0002-9341-7985`
names(res$`0000-0002-9341-7985`)
res$`0000-0002-9341-7985`$`last-modified`
res$`0000-0002-9341-7985`$`keywords`
res$`0000-0002-9341-7985`$`biography`
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_keywords.R
\name{orcid_keywords}
\alias{orcid_keywords}
\title{Get education information for a person}
\usage{
orcid_keywords(orcid, put_code = NULL, format = "application/json", ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get education information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all data
res <- orcid_keywords(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
names(res$`0000-0002-1642-628X`)
res$`0000-0002-1642-628X`$`keyword`

# individual ones
orcid_keywords("0000-0002-1642-628X", 31202)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_other_names.R
\name{orcid_other_names}
\alias{orcid_other_names}
\title{Get education information for a person}
\usage{
orcid_other_names(orcid, put_code = NULL, format = "application/json", ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get education information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all data
res <- orcid_other_names(orcid = "0000-0001-7893-4389")
res$`0000-0001-7893-4389`
names(res$`0000-0001-7893-4389`)
res$`0000-0001-7893-4389`$`other-name`

# individual ones
orcid_other_names("0000-0001-7893-4389", 239534)

# formats
orcid_other_names("0000-0001-7893-4389", format = "application/xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/browse.R
\name{browse}
\alias{browse}
\title{Navigate to an ORCID profile in your default browser}
\usage{
browse(orcid)
}
\arguments{
\item{orcid}{An \code{or_cid} class object}
}
\description{
Navigate to an ORCID profile in your default browser
}
\examples{
\dontrun{
browse(as.orcid("0000-0002-1642-628X"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_fundings.R
\name{orcid_fundings}
\alias{orcid_fundings}
\title{Get funding information for a person}
\usage{
orcid_fundings(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get funding summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get funding information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all funding data
res <- orcid_fundings(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
names(res$`0000-0002-1642-628X`)
res$`0000-0002-1642-628X`$`group`

# individual funding records
orcid_fundings(orcid = "0000-0002-1642-628X", 385627)

# funding summary information
orcid_fundings(orcid = "0000-0002-1642-628X", 385627, summary = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.orcid.R
\name{as.orcid}
\alias{as.orcid}
\title{Convert an ORCID or something like an ORCID object}
\usage{
as.orcid(x, ...)
}
\arguments{
\item{x}{An ORCID id, passed to \code{print}}

\item{...}{Further args passed on to \code{\link[=orcid_id]{orcid_id()}}}
}
\value{
an S3 object of class \code{or_cid}, which pretty prints
for brevity
}
\description{
Convert an ORCID or something like an ORCID object
}
\examples{
\dontrun{
as.orcid(x="0000-0002-1642-628X")
out <- orcid("text:English", rows = 20)
as.orcid(out$`orcid-identifier.path`[1])

# Passon further args to orcid_id()
as.orcid("0000-0002-1642-628X", verbose = TRUE)

# Browse to a profile
# browse(as.orcid("0000-0002-1642-628X"))

# many ORCIDs as a character vector
ids <- c("0000-0002-1642-628X", "0000-0002-9341-7985")
as.orcid(ids)

# many in a list via orcid_id()
(x <- lapply(ids, orcid_id))
as.orcid(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identifiers.R
\name{identifiers}
\alias{identifiers}
\alias{identifiers.works}
\alias{identifiers.list}
\alias{identifiers.orcid_id}
\alias{identifiers.orcid}
\alias{identifiers.orcid_doi}
\title{Get identifiers}
\usage{
identifiers(x, type = "doi", ...)

\method{identifiers}{works}(x, type = "doi", ...)

\method{identifiers}{list}(x, type = "doi", ...)

\method{identifiers}{orcid_id}(x, type = "doi", ...)

\method{identifiers}{orcid}(x, type = "doi", ...)

\method{identifiers}{orcid_doi}(x, type = "doi", ...)
}
\arguments{
\item{x}{An object of class works, orcid, orcid_id, orcid_doi, or a list
that contains any number of the previous objects.}

\item{type}{(character) One of doi (default), pmid, pmc, eid, other_id,
orcid, scopus, researcherid. The orcid's here are for works, not
individuals. This parameter is ignored for classes \code{orcid} and \code{orcid_doi}
both of which would go down a rabbit hole of getting works for all
ORCIDs which could take a while.}

\item{...}{Ignored.}
}
\value{
(character) vector of identifiers, or NULL if none found
}
\description{
This function aims to pluck out just identifiers into a vector
for easy use downstream (e.g., use DOIs to fetch article metadata). You
can still manually fetch additional data from outputs of functions in
this package.
}
\examples{
\dontrun{
# Result of call to works()
x <- works(orcid_id("0000-0001-8607-8025"))
# doi by default
identifiers(x)
# orcids
identifiers(x, "orcid")
# pmid
identifiers(x, "pmid")
# pmc 
identifiers(x, "pmc") 
# other_id
identifiers(x, "other_id")

# Result of call to orcid_id()
x <- orcid_id(orcid = "0000-0002-9341-7985")
identifiers(x, "doi")
identifiers(x, "eid")

# Result of call to orcid()
x <- orcid(query="carl+boettiger")
identifiers(x)

# Result of call to orcid_doi()
x <- orcid_doi(dois="10.1087/20120404", fuzzy=TRUE)
identifiers(x)
}
}
\references{
list of identifiers
https://pub.qa.orcid.org/v2.0/identifiers?locale=en
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_ping.R
\name{orcid_ping}
\alias{orcid_ping}
\title{Check if ORCID API is up and running}
\usage{
orcid_ping(...)
}
\arguments{
\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
a text string
}
\description{
Check if ORCID API is up and running
}
\examples{
\dontrun{
orcid_ping()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/works.R
\name{works}
\alias{works}
\title{Get works data}
\usage{
works(x)
}
\arguments{
\item{x}{Anything that can be coerced via \code{\link[=as.orcid]{as.orcid()}}, see
\code{\link[=as.orcid]{as.orcid()}} for help}
}
\value{
A tibble (data.frame)
}
\description{
Get works data
}
\details{
This function gets works using the function \link{orcid_works}
and packages up the data in a data.frame for easier processing
}
\examples{
\dontrun{
out <- works(orcid_id("0000-0002-9341-7985"))
out
out$type
out$path

works( orcid_id("0000-0002-1642-628X") )
works( orcid_id("0000-0003-1444-9135") )
works( orcid_id("0000-0003-1419-2405") )

out <- orcid(query="keyword:ecology")
works(orcid_id(out$`orcid-identifier.path`[7]))
works(orcid_id(out$`orcid-identifier.path`[8]))
works(orcid_id(out$`orcid-identifier.path`[9]))
works(orcid_id(out$`orcid-identifier.path`[10]))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_address.R
\name{orcid_address}
\alias{orcid_address}
\title{Get address information for a person}
\usage{
orcid_address(orcid, put_code = NULL, format = "application/json", ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get address information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all addresses
res <- orcid_address(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
names(res$`0000-0002-1642-628X`)
res$`0000-0002-1642-628X`$`address`

# individual address
orcid_address(orcid = "0000-0002-1642-628X", 288064)

# format
orcid_address(orcid = "0000-0002-1642-628X", 288064, "application/xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_dois.R
\name{check_dois}
\alias{check_dois}
\title{Verify DOI's are likely good}
\usage{
check_dois(x)
}
\arguments{
\item{x}{One or more DOIs}
}
\value{
A list of length two, one slot for good DOIs, one for bad
}
\description{
Verify DOI's are likely good
}
\examples{
\dontrun{
check_dois("10.1087/20120404")

dois=c("10.1371/journal.pone.0025995","10.1371/journal.pone.0053712",
       "10.1371/journal.pone.0054608","10.1371/journal.pone.0055937")
check_dois(dois)

dois=c("10.1016/j.medpal.2008.12.005","10.1080/00933104.2000.10505926",
       "10.1037/a0024480", "10.1002/anie.196603172","2344","asdf","232",
       "asdf","23dd")
check_dois(dois)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rorcid-package.R
\name{rorcid-defunct}
\alias{rorcid-defunct}
\title{Defunct functions in rorcid}
\description{
\itemize{
\item \code{\link[=summary.or_cid]{summary.or_cid()}}: Function is gone. Deemed not really
that useful, and hard to maintain given other changes in the package.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rorcid-package.R
\docType{data}
\name{issn_title}
\alias{issn_title}
\title{Lookup vector for journal titles by ISSN}
\description{
named vector of journal titles. the values are journal titles and
the names are ISSN's.
}
\details{
length: 57,968

data collected on 2018-06-13 from Crossref
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_researcher_urls.R
\name{orcid_researcher_urls}
\alias{orcid_researcher_urls}
\title{Get researcher urls for a person}
\usage{
orcid_researcher_urls(orcid, put_code = NULL, format = "application/json", ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get researcher urls for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all data
res <- orcid_researcher_urls(orcid = "0000-0003-1444-9135")
res$`0000-0003-1444-9135`
names(res$`0000-0003-1444-9135`)
res$`0000-0003-1444-9135`$`researcher-url`

# individual ones
orcid_researcher_urls("0000-0003-1444-9135", 304093)
orcid_researcher_urls("0000-0003-1444-9135", c(332241, 304093))

# formats
orcid_researcher_urls("0000-0003-1444-9135", 304093, 
  format = "application/xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_citations.R
\name{orcid_citations}
\alias{orcid_citations}
\title{Get citations}
\usage{
orcid_citations(
  orcid,
  put_code = NULL,
  cr_format = "bibtex",
  cr_style = "apa",
  cr_locale = "en-US",
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s) of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{cr_format}{Used in Crossref queries only. Name of the format. One of
"rdf-xml", "turtle",
"citeproc-json", "citeproc-json-ish", "text", "ris", "bibtex" (default),
"crossref-xml", "datacite-xml","bibentry", or "crossref-tdm". The format
"citeproc-json-ish" is a format that is not quite proper citeproc-json.
passed to \code{rcrossref::cr_cn}. The special "citeproc2bibtex" value asks
for citeproc-json from Crossref, then converts it into bibtex format
using \link[handlr:HandlrClient]{handlr::HandlrClient}}

\item{cr_style}{Used in Crossref queries only. A CSL style (for text
format only). See 'get_styles()' for options. Default: apa.
passed to \code{rcrossref::cr_cn}}

\item{cr_locale}{Used in Crossref queries only. Language locale.
See \link{Sys.getlocale}, passed to \code{rcrossref::cr_cn}}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID

data.frame, with the columns:
\itemize{
\item put: ORCID PUT code, identifying the work identifier in ORCID's records
\item id: the external identifier
\item id_type: the type of external identifier
\item format: the citation format retrieved
\item citation: the citation as JSON
}
}
\description{
Get citations
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.

This function is focused on getting citations only.
You can get all citations for an ORCID, or for certain works
using a PUT code, or for many PUT codes.

We attempt to get citations via Crossref using \pkg{rcrossref}
whenever possible as they are the most flexible and don't have as
many mistakes in the text. If there is no DOI, we fetch the
citation from ORCID.

Right now we get JSON citations back. We'd like to support bibtex
format. DOI.org supports this but not ORCID.
}
\examples{
\dontrun{
(res <- orcid_citations(orcid = "0000-0002-9341-7985"))
(res2 <- orcid_citations(orcid = "0000-0002-1642-628X"))
(res2 <- orcid_citations(orcid = c("0000-0002-9341-7985", "0000-0002-1642-628X")))

# get individual works
## a single put code
(a <- orcid_citations(orcid = "0000-0002-9341-7985", put_code = 5011717))
## many put codes
(b <- orcid_citations(orcid = "0000-0002-9341-7985",
   put_code = c(5011717, 15536016)))

# request other formats, Crossref only
orcid_citations(orcid = "0000-0002-9341-7985", cr_format = "turtle")

# parse citation data if you wish
# for parsing bibtex can use bibtex package or others
(res <- orcid_citations(orcid = "0000-0002-9341-7985"))
lapply(res[res$format == "csl-json", "citation"][[1]], jsonlite::fromJSON)

# lots of citations
orcid_citations(orcid = "0000-0001-8642-6325")

# example with no external identifier, returns NA's
orcid_citations(orcid = "0000-0001-8642-6325", 26222265)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_id.r
\name{orcid_id}
\alias{orcid_id}
\title{Get data for particular ORCID's}
\usage{
orcid_id(orcid, ...)
}
\arguments{
\item{orcid}{(character) A single Orcid identifier, of the
form XXXX-XXXX-XXXX-XXXX}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A named list of results - from a call to \code{\link[=orcid_person]{orcid_person()}}
}
\description{
Get data for particular ORCID's
}
\examples{
\dontrun{
res <- orcid_id(orcid = "0000-0002-9341-7985")
res$`0000-0002-9341-7985`
res$`0000-0002-9341-7985`$`name`
res$`0000-0002-9341-7985`$`other-names`
res$`0000-0002-9341-7985`$`biography`
res$`0000-0002-9341-7985`$`researcher-urls`
res$`0000-0002-9341-7985`$`emails`
res$`0000-0002-9341-7985`$`addresses`
res$`0000-0002-9341-7985`$`keywords`
res$`0000-0002-9341-7985`$`external-identifiers`
res$`0000-0002-9341-7985`$`emails`

ids <- c("0000-0003-1620-1408", "0000-0002-9341-7985")
res <- lapply(ids, orcid_id)
vapply(res, function(x) x[[1]]$name$`family-name`$value, "")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_works.R
\name{orcid_works}
\alias{orcid_works}
\title{Get works for a person}
\usage{
orcid_works(orcid, put_code = NULL, format = "application/json", ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get works for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# get all works
res <- orcid_works(orcid = "0000-0002-9341-7985")
res$`0000-0002-9341-7985`
res$`0000-0002-9341-7985`$works
res$`0000-0002-9341-7985`$works$type
str(res$`0000-0002-9341-7985`)

# get individual works
a <- orcid_works(orcid = "0000-0002-9341-7985", put_code = 5011717)
a$`0000-0002-9341-7985`
a$`0000-0002-9341-7985`$works
b <- orcid_works(orcid = "0000-0002-9341-7985", 
   put_code = c(5011717, 15536016))
b$`0000-0002-9341-7985`

# change formats
orcid_works("0000-0002-9341-7985", 5011717, "application/json")
orcid_works("0000-0002-9341-7985", 5011717, "application/xml")
orcid_works("0000-0002-9341-7985", 5011717, 
  "application/vnd.orcid+xml; qs=5")
orcid_works("0000-0002-9341-7985", 5011717, 
  "application/vnd.citationstyles.csl+json")

# get citations
id <- "0000-0001-7678-8656"
x <- orcid_works(id)
wks <- orcid_works(id, put_code = x[[1]]$works$`put-code`)
wks[[1]]$works$`work.citation.citation-value`

## or send many put codes at once, will be split into chunks of 50 each
id <- "0000-0001-6758-5101"
z <- orcid_works(id)
pcodes <- z[[1]]$works$`put-code`
length(pcodes)
res <- orcid_works(orcid = id, put_code = pcodes)
head(res$`0000-0001-6758-5101`$works)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_external_identifiers.R
\name{orcid_external_identifiers}
\alias{orcid_external_identifiers}
\title{Get external identifiers for a person}
\usage{
orcid_external_identifiers(
  orcid,
  put_code = NULL,
  format = "application/json",
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get external identifiers for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all data
res <- orcid_external_identifiers(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
names(res$`0000-0002-1642-628X`)
res$`0000-0002-1642-628X`$`external-identifier`

# individual records
orcid_external_identifiers(orcid = "0000-0002-1642-628X", 141736)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_distinctions.R
\name{orcid_distinctions}
\alias{orcid_distinctions}
\title{Get distinction for a person}
\usage{
orcid_distinctions(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get education summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get distinction for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_distinctions(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
res$`0000-0002-1642-628X`$`created-date`
res$`0000-0002-1642-628X`$`affiliation-group`
res$`0000-0002-1642-628X`$path
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_qualifications.R
\name{orcid_qualifications}
\alias{orcid_qualifications}
\title{Get qualifications for a person}
\usage{
orcid_qualifications(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get peer review summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get qualifications for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_services.R
\name{orcid_services}
\alias{orcid_services}
\title{Get services}
\usage{
orcid_services(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get education summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get services
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_services(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
res$`0000-0002-1642-628X`$`last-modified-date`
res$`0000-0002-1642-628X`$`affiliation-group`
res$`0000-0002-1642-628X`$path
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_employments.R
\name{orcid_employments}
\alias{orcid_employments}
\title{Get employment information for a person}
\usage{
orcid_employments(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get employment summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get employment information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all employment data
res <- orcid_employments(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
names(res$`0000-0002-1642-628X`)
res$`0000-0002-1642-628X`$`employment-summary`

# individual employment records
orcid_employments(orcid = "0000-0002-1642-628X", 1115445)
orcid_employments(orcid = "0000-0002-1642-628X", 148496)

# employment summary information
orcid_employments(orcid = "0000-0002-1642-628X", 1115445, summary = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_activities.R
\name{orcid_activities}
\alias{orcid_activities}
\title{Get activities for a person}
\usage{
orcid_activities(orcid, ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get activities for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_activities(orcid = "0000-0002-9341-7985")
res$`0000-0002-9341-7985`
names(res$`0000-0002-9341-7985`)
res$`0000-0002-9341-7985`$`last-modified`
res$`0000-0002-9341-7985`$`educations`
res$`0000-0002-9341-7985`$`fundings`
res$`0000-0002-9341-7985`$`peer-reviews`
res$`0000-0002-9341-7985`$`works`
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_memberships.R
\name{orcid_memberships}
\alias{orcid_memberships}
\title{Get memberships for a person}
\usage{
orcid_memberships(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get education summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get memberships for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_memberships(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
res$`0000-0002-1642-628X`$`created-date`
res$`0000-0002-1642-628X`$`affiliation-group`
res$`0000-0002-1642-628X`$memberships
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{orcid_auth}
\alias{orcid_auth}
\alias{rorcid-auth}
\title{ORCID authorization}
\usage{
orcid_auth(
  scope = "/authenticate",
  reauth = FALSE,
  redirect_uri = getOption("rorcid.redirect_uri"),
  client_id = NULL,
  client_secret = NULL
)
}
\arguments{
\item{scope}{(character) one or more scopes. default: \code{"/authenticate"}.
see "ORCID OAuth Scopes" section below for other scope options}

\item{reauth}{(logical) Force re-authorization? default: \code{FALSE}}

\item{redirect_uri}{(character) a redirect URI. optional. set by passing
to this parameter or using the R option \code{rorcid.redirect_uri}}

\item{client_id}{(character) a client id. optional}

\item{client_secret}{(character) a client secret. optional}
}
\value{
a character string with the access token prefixed with "Bearer "
}
\description{
ORCID authorization
}
\details{
There are three ways to authorise with \pkg{rorcid}:
\itemize{
\item Interactively login with OAuth. This doesn't require any input on
your part. We use a client id and client secret key to ping ORCID.org;
at which point you log in with your username/password; then we get back
a token (same as the above option). We don't know your username or
password, only the token that we get back. We cache that token locally
in a hidden file in whatever working directory you're in. If you delete
that file, or run the code from a new working directory, then we
re-authorize.
\item Use a \code{client_id} and \code{client_secret} to do 2-legged OAuth.
ORCID docs at https://members.orcid.org/api/oauth/2legged-oauth and
https://members.orcid.org/api/post-oauthtoken-reading-public-data
This requires you to register a "client application". See
https://orcid.org/content/register-client-application-2 for instructions
\item Use a token as a result of either of the two above approaches. The token
is a alphanumeric UUID, e.g. \verb{dc0a6b6b-b4d4-4276-bc89-78c1e9ede56e}. You
can get this token by running \code{orcid_auth()}, then storing that key
(the uuid alone, not the "Bearer " part) either as en environment
variable in your \code{.Renviron} file in your home directory (with the name
\code{ORCID_TOKEN}), or as an R option in your \code{.Rprofile} file (with the name
\code{orcid_token}). See \link{Startup} for more information.
Either an environment variable or R option work. If we don't find
either we do the next option.
}

We recommend the 3rd option if possible, specifically, storing the token
as an environment variable permanently.

If authentication fails, you can still use \pkg{rorcid}. ORCID does not require
authentication at this point, but may in the future - this prepares you
for when that happens :)
}
\note{
This function is used within \pkg{rorcid} to get/do authentication.
}
\section{ORCID OAuth Scopes}{

https://info.orcid.org/faq/what-is-an-oauth-scope-and-which-scopes-does-orcid-support/
}

\section{Computing environments without browsers}{

One pitfall is when you are using rorcid on a server, and you're ssh'ed
in, so that there's no way to open a browser to do the OAuth browser
flow. Similarly for any other situation in which a browser can not be
opened. In this case, run \code{orcid_auth()} on another machine in which you do
have the ability to open a browser, then collect the info that's ouptput
from \code{orcid_auth()} and store it as an environment variable (see above).
}

\examples{
\dontrun{
x <- orcid_auth()
orcid_auth(reauth = TRUE)
# orcid_auth(scope = "/read-public", reauth = TRUE)

# supply client_id AND client_secret to avoid 3 legged, interactive OAuth
# orcid_auth(client_id = "---", client_secret = "---")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_email.R
\name{orcid_email}
\alias{orcid_email}
\title{Get education information for a person}
\usage{
orcid_email(orcid, ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get education information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_email(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
names(res$`0000-0002-1642-628X`)
res$`0000-0002-1642-628X`$`email`
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_research_resources.R
\name{orcid_research_resources}
\alias{orcid_research_resources}
\title{Get research resources for a person}
\usage{
orcid_research_resources(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get education summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get research resources for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_research_resources(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
res$`0000-0002-1642-628X`$`last-modified-date`
res$`0000-0002-1642-628X`$`group`
res$`0000-0002-1642-628X`$path
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated-defunct.R
\name{summary.or_cid}
\alias{summary.or_cid}
\title{This function is defunct.}
\usage{
\method{summary}{or_cid}(object, ...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_search.R
\name{orcid_search}
\alias{orcid_search}
\title{Orcid search - more user friendly than \code{\link[=orcid]{orcid()}}}
\usage{
orcid_search(
  given_name = NULL,
  family_name = NULL,
  past_inst = NULL,
  current_inst = NULL,
  affiliation_org = NULL,
  ringgold_org_id = NULL,
  grid_org_id = NULL,
  credit_name = NULL,
  other_name = NULL,
  email = NULL,
  digital_object_ids = NULL,
  work_title = NULL,
  grant_number = NULL,
  keywords = NULL,
  text = NULL,
  rows = 10,
  start = NULL,
  ...
)
}
\arguments{
\item{given_name}{(character) given name}

\item{family_name}{(character) family name}

\item{past_inst}{(character) past institution}

\item{current_inst}{(character) current institution}

\item{affiliation_org}{(character) affiliation organization name}

\item{ringgold_org_id}{(character) ringgold organization id}

\item{grid_org_id}{(character) grid organization id}

\item{credit_name}{(character) credit name}

\item{other_name}{(character) other name}

\item{email}{(character) email}

\item{digital_object_ids}{(character) digital object ids}

\item{work_title}{(character) work title}

\item{grant_number}{(character) grant number}

\item{keywords}{(character) keywords to search. character vector, one
or more keywords}

\item{text}{(character) text to search}

\item{rows}{(integer) number of records to return}

\item{start}{(integer) record number to start at}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
a \code{data.frame} with three columns:
\itemize{
\item first: given name
\item last: family name
\item orcid: ORCID identifier
}

If no results are found, an empty (0 rows) data.frame
is returned
}
\description{
Orcid search - more user friendly than \code{\link[=orcid]{orcid()}}
}
\details{
The goal of this function is to make a human friendly
way to search ORCID.

Thus, internally we map the parameters given to this function
to the actual parameters that ORCID wants that are not
so human friendly.

We don't include all possible fields you could search against
here - for that use \code{\link[=orcid]{orcid()}}

Importantly, we return the first 10 results, following the default
setting for the \code{rows} parameter in \code{\link[=orcid]{orcid()}}. You can set the rows
parameter in this function to a max of 200. The maximum is an
upper bound set by the ORCID API. You can get the number of results
found programatically by fetching the \code{found} attribute on the ouput
of this function, e.g., \code{attr(x, "found")}.
}
\note{
\code{current_prim_inst} and \code{patent_number} parameters have been
removed as ORCID has removed them
}
\section{How parameters are combined}{

We combine multiple parameters with \code{AND}, such that
e.g., \code{given_name="Jane"} and \code{family_name="Doe"} gets passed
to ORCID as \verb{given-names:Jane AND family-name:Doe}
}

\examples{
\dontrun{
orcid_search(given_name = "carl", family_name = "boettiger")
orcid_search(given_name = "carl")
orcid_search(given_name = "carl", rows = 2)
orcid_search(keywords = c("birds", "turtles"))
orcid_search(affiliation_org = '("Boston University" OR BU)')
orcid_search(ringgold_org_id = '1438')
orcid_search(grid_org_id = 'grid.5509.9')
orcid_search(email = '*@orcid.org')
orcid_search(given_name = "carl", verbose = TRUE)
# get number of results found
x <- orcid_search(ringgold_org_id = '1438')
attr(x, "found")
}
}
\references{
\url{https://members.orcid.org/api/tutorial/search-orcid-registry}
}
\seealso{
\code{\link[=orcid]{orcid()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid.r
\name{orcid}
\alias{orcid}
\title{Search for ORCID ID's.}
\usage{
orcid(
  query = NULL,
  start = NULL,
  rows = NULL,
  defType = NULL,
  q.alt = NULL,
  qf = NULL,
  mm = NULL,
  qs = NULL,
  pf = NULL,
  ps = NULL,
  pf2 = NULL,
  ps2 = NULL,
  pf3 = NULL,
  ps3 = NULL,
  tie = NULL,
  bq = NULL,
  bf = NULL,
  boost = NULL,
  uf = NULL,
  lowercaseOperators = NULL,
  fuzzy = FALSE,
  recursive = FALSE,
  ...
)
}
\arguments{
\item{query}{Search terms. You can do quite complicated queries using the
SOLR syntax. See examples below. For all possible fields to query, do
\code{data(fields)}}

\item{start}{Result number to start on. Keep in mind that pages start at 0.
Default: 0}

\item{rows}{Numer of results to return. Default: 10. Max: 200}

\item{defType}{Query syntax. One of edismax or X. See Details for more.}

\item{q.alt}{If specified, this query will be used (and parsed by default
using standard query parsing syntax) when the main query string is not
specified or blank. This comes in handy when you need something like a
match-all-docs query (don't forget &rows=0 for that one!) in order to get
collection-wise faceting counts.}

\item{qf}{(Query Fields) List of fields and the "boosts" to associate with
each of them when building DisjunctionMaxQueries from the user's query}

\item{mm}{(Minimum 'Should' Match).}

\item{qs}{(Query Phrase Slop) Amount of slop on phrase queries explicitly
included in the user's query string (in qf fields; affects matching).}

\item{pf}{(Phrase Fields) Once the list of matching documents has been
identified using the "fq" and "qf" params, the "pf" param can be used to
"boost" the score of documents in cases where all of the terms in the "q"
param appear in close proximity}

\item{ps}{(Phrase Slop) Default amount of slop on phrase queries built with
"pf", "pf2" and/or "pf3" fields (affects boosting).}

\item{pf2}{(Phrase bigram fields) As with 'pf' but chops the input into
bi-grams, e.g. "the brown fox jumped" is queried as "the brown" "brown fox"
"fox jumped"}

\item{ps2}{(Phrase bigram slop) As with 'ps' but sets default slop factor for
'pf2'. If not specified, 'ps' will be used.}

\item{pf3}{(Phrase trigram fields) As with 'pf' but chops the input into
tri-grams, e.g. "the brown fox jumped" is queried as "the brown fox"
"brown fox jumped"}

\item{ps3}{(Phrase trigram slop) As with 'ps' but sets default slop factor
for 'pf3'. If not specified, 'ps' will be used.}

\item{tie}{(Tie breaker) Float value to use as tiebreaker in
DisjunctionMaxQueries (should be something much less than 1)}

\item{bq}{(Boost Query) A raw query string (in the SolrQuerySyntax) that will
be included with the user's query to influence the score. See references}

\item{bf}{(Boost Function, additive) Functions (with optional boosts) that
will be included in the user's query to influence the score. Any function
supported natively by Solr can be used, along with a boost value, e.g.:
recip(rord(myfield),1,2,3)^1.5}

\item{boost}{(Boost Function, multiplicative) As for 'bf' but multiplies the
boost into the  score}

\item{uf}{(User Fields) Specifies which schema fields the end user shall be
allowed to query for explicitly. This parameter supports wildcards.}

\item{lowercaseOperators}{This param controls whether to try to interpret
lowercase words as boolean operators such as "and", "not" and "or".
Set &lowercaseOperators=true to allow this. Default is "false".}

\item{fuzzy}{Use fuzzy matching on input DOIs. Defaults to FALSE. If FALSE,
we stick "digital-object-ids" before the DOI so that the search sent to
ORCID is for that exact DOI. If TRUE, we use some regex to find the DOI.}

\item{recursive}{DEFUNCT}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
a data.frame (tibble). You can access number of results found like
\code{attr(result, "found")}. Note that with ORCID API v2 and greater,
results here are only the identifiers. To get other metadata/data
you can take the identifiers and use other functions in this package.
}
\description{
Search for ORCID ID's.
}
\details{
All query syntaxes available in SOLR 3.6 (
\url{https://lucene.apache.org/solr/guide/8_7/the-standard-query-parser.html})
are supported, including Lucene with Solr extensions (default), DisMax,
and Extended Dismax.

You can use any of the following within the query statement:
given-names, family-name, credit-name, other-names, email, grant-number,
patent-number, keyword, worktitle, digital-objectids, current-institution,
affiliation-name, current-primary-institution, text, past-institution,
peer-review-type, peer-review-role, peer-review-group-id, biography,
external-id-type-and-value

For more complicated queries the ORCID API supports using ExtendedDisMax.
See the documentation on the web here:
\url{https://lucene.apache.org/solr/guide/8_7/the-extended-dismax-query-parser.html}

Note that when constructing queries, you don't need to use syntax like
\code{+}, etc., \code{crul}, the http client we use internally, will do that
for you. For example, instead of writing \code{johnson+cardiology}, just
write \verb{johnson cardiology}, and instead of writing
\code{johnson+AND+cardiology}, write \verb{johnson AND cardiology}. Though,
you still need to use \code{AND}, \code{OR}, etc. to join term/queries
together.
}
\examples{
\dontrun{
# Get a list of names and Orcid IDs matching a name query
orcid(query="carl+boettiger")
orcid(query="given-names:carl AND family-name:boettiger")

# by email
orcid(query="email:cboettig@berkeley.edu")

# You can string together many search terms
orcid(query="johnson cardiology houston")

# peer review group id
orcid("peer-review-group-id:1996-3068")

# And use boolean operators
orcid("johnson AND(caltech OR 'California Institute of Technology')")

# And you can use start and rows arguments to do pagination
orcid("johnson cardiology houston", start = 2, rows = 3)

# Use search terms, here family name
orcid("family-name:Sanchez", start = 4, rows = 6)

# Use search terms, here...
orcid(query="Raymond", start=0, rows=10, defType="edismax")

# Search using keywords
orcid(query="keyword:ecology")

# Search by DOI
orcid(query="10.1087/20120404")

# Note the difference between the first wrt the second and third
## See also orcid_doi() function for searching by DOIs
orcid("10.1087/20120404")
orcid('"10.1087/20120404"')
## doi
orcid('digital-object-ids:"10.1087/20120404"')
## doi prefix
orcid('digital-object-ids:"10.1087/*"')

# search by work titles
orcid('work-titles:Modern developments in holography and its materials')
orcid('pmc:PMC3901677')

## Using more complicated SOLR queries

# Use the qf parameter to "boost" query fields so they are ranked higher
# 	See how it is different than the second query without using "qf"
orcid(defType = "edismax", query = "Raymond",
   qf = "given-names^1.0 family-name^2.0", start = 0, rows = 10)
orcid(query = "Raymond", start = 0, rows = 10)

# Use other SOLR parameters as well, here mm. Using the "mm" param, 1 and
# 2 word queries require that all of the optional clauses match, but for
# queries with three or more clauses one missing clause is allowed...
# See for more: http://bit.ly/1uyMLDQ
orcid(defType = "edismax",
      query="keyword:ecology OR evolution OR conservation",
      mm = 2, rows = 20)
}
}
\references{
\url{https://members.orcid.org/api/tutorial/search-orcid-registry}
\url{https://lucene.apache.org/solr/guide/8_7/the-extended-dismax-query-parser.html}
}
\seealso{
\code{\link[=orcid_doi]{orcid_doi()}} \code{\link[=orcid_id]{orcid_id()}} \code{\link[=orcid_search]{orcid_search()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rorcid-package.R
\docType{package}
\name{rorcid-package}
\alias{rorcid-package}
\title{A programmatic R interface the Orcid.org API}
\description{
A R interface to the Orcid public API. \pkg{rorcid} is not a product
developed or distributed by ORCID.

ORCID website: https://orcid.org/

Orcid API docs: http://members.orcid.org/api

Some key \pkg{rorcid} function:
\itemize{
\item \code{\link[=as.orcid]{as.orcid()}} - coerce various inputs to ORCID class
\item \code{\link[=browse]{browse()}} - browse to a profile in your default browser
\item \code{\link[=check_dois]{check_dois()}} - check that strings are likely to be DOIs
\item \code{\link[=identifiers]{identifiers()}} - grab identifiers out of various objects
\item \code{\link[=orcid]{orcid()}} and \code{\link[=orcid_search]{orcid_search()}} - Search for ORCID id's
\item \code{\link[=orcid_doi]{orcid_doi()}} - Search by DOI
\item \code{\link[=orcid_id]{orcid_id()}} - Search by ORCID id, and get either bio,
profile, or works
\item \code{\link[=works]{works()}} - Parse out works from various objects
}
}
\section{API routes not implemented}{


Not quite sure what these do so haven't messed with them.
\itemize{
\item \verb{/\{orcid\}/notification-permission/\{id\}}
\item \verb{/client/\{client_id\}}
\item \verb{/group-id-record}
\item \verb{/group-id-record/\{putCode\}}
}
}

\section{Rate Limits}{

Definitions:
\itemize{
\item Request a second - Number of request that can be made a second.
Value: 8 per second (24 with API v2rc+) - Haven't been able to find up
to date values for API v3 (so assume they are the same I guess)
\item Burst - Number of request we will allow to be queued before rejecting.
The request in the queue are slowed down at the request a second rate.
Value: 40 (same with API v2rc+) - Haven't been able to find up to date
values for API v3 (so assume they are the same I guess)
}

If you exceed the burst, you'll get a 503 responses. Developers should do
their best to avoid approaching those limits.
}

\seealso{
\link{rorcid-auth} for Authentication information
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_doi.r
\name{orcid_doi}
\alias{orcid_doi}
\title{Search for ORCID ID's using DOIs}
\usage{
orcid_doi(dois = NULL, start = NULL, rows = NULL, fuzzy = FALSE, ...)
}
\arguments{
\item{dois}{(character) Digital object identifier (DOI), a vector fo DOIs.}

\item{start}{(integer) Result number to start on. Keep in mind that pages
start at 0.}

\item{rows}{(integer) Numer of results to return.}

\item{fuzzy}{(logical) Use fuzzy matching on input DOIs. Defaults to FALSE.
If FALSE, we stick "digital-object-ids" before the DOI so that the search
sent to ORCID is for that exact DOI. If TRUE, we use some regex to find
the DOI.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Search for ORCID ID's using DOIs
}
\examples{
\dontrun{
orcid_doi(dois="10.1087/20120404", fuzzy=TRUE)

# fuzzy is FALSE by default
orcid_doi(dois="10.1087/20120404", fuzzy=FALSE)

# This DOI is not a real one, but a partial DOI, then we can fuzzy search
# get more than defualt 10 records (or rows)
orcid_doi(dois="10.1087/2", fuzzy=TRUE, rows=20) 

# If you don't input proper DOIs, the function will get mad
dois <- c("10.1371/journal.pone.0025995","10.1371/journal.pone.0053712",
       "10.1371/journal.pone.0054608","10.1371/journal.pone.0055937")
orcid_doi(dois=dois)

# dois <- c("10.1016/j.medpal.2008.12.005","10.1080/00933104.2000.10505926",
#          "10.1037/a0024480", "10.1002/anie.196603172","2344","asdf","232",
#          "asdf","23dd")
# orcid_doi(dois=dois)

orcid_doi(dois="10.1087/20120404", fuzzy=FALSE) 
orcid_doi(dois="10.1371/journal.pone.0025995", fuzzy=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_peer_reviews.R
\name{orcid_peer_reviews}
\alias{orcid_peer_reviews}
\title{Get peer review information for a person}
\usage{
orcid_peer_reviews(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get peer review summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get peer review information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all peer review data
res <- orcid_peer_reviews(orcid = "0000-0001-7678-8656")
res$`0000-0001-7678-8656`
names(res$`0000-0001-7678-8656`)
res$`0000-0001-7678-8656`$`group`

# get individual works
orcid_peer_reviews("0000-0003-1444-9135", 75565)

# summary
orcid_peer_reviews("0000-0003-1444-9135", 75565, summary = TRUE)

# get Journal titles via ISSN's provided in results, using the 
# provided issn_title dataset
x <- orcid_peer_reviews("0000-0001-7678-8656", put_code = "220419")
issn <- strsplit(x[[1]]$`review-group-id`, ":")[[1]][[2]]
issn_title[[issn]]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rorcid-package.R
\docType{data}
\name{fields}
\alias{fields}
\title{Lookup table for search fields}
\description{
Lookup table for search fields
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_bio.R
\name{orcid_bio}
\alias{orcid_bio}
\title{Get biography data for a person}
\usage{
orcid_bio(orcid, format = "application/json", ...)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get biography data for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_bio(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
res$`0000-0002-1642-628X`$`created-date`
res$`0000-0002-1642-628X`$`last-modified-date`
res$`0000-0002-1642-628X`$`content`
res$`0000-0002-1642-628X`$`visibility`
res$`0000-0002-1642-628X`$path
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_invited_positions.R
\name{orcid_invited_positions}
\alias{orcid_invited_positions}
\title{Get invited positions for a person}
\usage{
orcid_invited_positions(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get education summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get invited positions for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
res <- orcid_invited_positions(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
res$`0000-0002-1642-628X`$`created-date`
res$`0000-0002-1642-628X`$`affiliation-group`
res$`0000-0002-1642-628X`$path
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orcid_educations.R
\name{orcid_educations}
\alias{orcid_educations}
\title{Get education information for a person}
\usage{
orcid_educations(
  orcid,
  put_code = NULL,
  format = "application/json",
  summary = FALSE,
  ...
)
}
\arguments{
\item{orcid}{(character) Orcid identifier(s), of the form
XXXX-XXXX-XXXX-XXXX. required.}

\item{put_code}{(character/integer) one or more put codes. up to
50. optional}

\item{format}{(character) Name of the content-type format. One of
"application/vnd.orcid+xml; qs=5", "application/orcid+xml; qs=3",
"application/xml", "application/vnd.orcid+json; qs=4",
"application/orcid+json; qs=2", "application/json"
"application/vnd.citationstyles.csl+json". optional}

\item{summary}{(logical) get education summary for a put code.
Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
A list of results for each Orcid ID passed in, with each element
named by the Orcid ID
}
\description{
Get education information for a person
}
\details{
This function is vectorized, so you can pass in many ORCID's, and
there's an element returned for each ORCID you put in.
}
\examples{
\dontrun{
# all education data
res <- orcid_educations(orcid = "0000-0002-1642-628X")
res$`0000-0002-1642-628X`
names(res$`0000-0002-1642-628X`)
res$`0000-0002-1642-628X`$`education-summary`

# individual education records
orcid_educations(orcid = "0000-0002-1642-628X", 148494)

# education summary information
orcid_educations(orcid = "0000-0002-1642-628X", 148494, summary = TRUE)
}
}
