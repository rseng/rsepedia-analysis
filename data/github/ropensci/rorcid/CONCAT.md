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

*Wow, no problems at all. :)**Wow, no problems at all. :)*