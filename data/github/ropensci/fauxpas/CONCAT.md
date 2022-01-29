fauxpas
=======



[![R-check](https://github.com/ropensci/fauxpas/workflows/R-check/badge.svg)](https://github.com/ropensci/fauxpas/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/fauxpas)](https://cranchecks.info/pkgs/fauxpas)
[![cran version](http://www.r-pkg.org/badges/version/fauxpas)](https://cran.r-project.org/package=fauxpas)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/fauxpas)](https://github.com/metacran/cranlogs.app)

`fauxpas` does http errors

* HTTP error classes more in line with Ruby/Python/Etc.
* An error class for each HTTP status in case a user wants to
be specific to an HTTP status code, and general purpose handlers
for any error
* Work with any of the major R http clients: `crul`, `curl`, `httr`, (maybe
`RCurl` later)
* Provide flexiblity for what to do on an HTTP error, including
custom functions and message templates

Info Links:

* HTTP on wikipedia: <https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol>
* HTTP2 on wikipedia: <https://en.wikipedia.org/wiki/HTTP/2>
* HTTP/1.1: <https://www.w3.org/Protocols/rfc2616/rfc2616-sec10.html>
* HTTP/2: <https://tools.ietf.org/html/rfc7540>

## Install

CRAN version


```r
install.packages("fauxpas")
```

Dev version


```r
remotes::install_github("ropensci/fauxpas")
```


```r
library("fauxpas")
```

## use with crul


```r
library("crul")
cli <- HttpClient$new("https://httpbin.org/status/414")
res <- cli$get()
http(res)
#> Error: Request-URI Too Long (HTTP 414).
http414(res)
#> Error: Request-URI Too Long (HTTP 414).
```


```r
x <- HTTPRequestURITooLong$new()
x$do_verbose(res)
#> Error: Request-URI Too Long (HTTP 414).
#> - The server is refusing to service the request because the Request-URI is
#>    longer than the server is willing to interpret. This rare condition is only likely
#>    to occur when a client has improperly converted a POST request to a GET request
#>    with long query information, when the client has descended into a URI black hole
#>    of redirection (e.g., a redirected URI prefix that points to a suffix of itself),
#>    or when the server is under attack by a client attempting to exploit security
#>    holes present in some servers using fixed-length buffers for reading or
#>    manipulating the Request-URI.
```

## use with curl


```r
library("curl")
h <- curl::new_handle()
curl::handle_setopt(h)
resp <- curl::curl_fetch_memory("https://httpbin.org/status/404", h)
http(resp)
#> Error: Not Found (HTTP 404).
http404(resp)
#> Error: Not Found (HTTP 404).
```


```r
x <- HTTPNotFound$new()
x$do_verbose(resp)
#> Error:  Not Found (HTTP 404).
#>  - The server has not found anything matching the Request-URI. No indication is
#> given of whether the condition is temporary or permanent. The 410 (Gone) status
#> code SHOULD be used if the server knows, through some internally configurable
#> mechanism, that an old resource is permanently unavailable and has no forwarding
#> address. #> This status code is commonly used when the server does not wish to
#> reveal exactly why the request has been refused, or when no other response is
#> applicable.
```

## use with httr


```r
library("httr")
res <- GET("https://httpbin.org/status/405")
http405(res)
#> Error: Method Not Allowed (HTTP 405).
```


```r
x <- HTTPMethodNotAllowed$new()
x$do_verbose(res)
#> Error: Method Not Allowed (HTTP 405).
#>  - The method specified in the Request-Line is not allowed for the resource
#> identified by the Request-URI. The response MUST include an Allow header
#> containing a list of valid methods for the requested resource.
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/fauxpas/issues)
* License: MIT
* Get citation information for `fauxpas` in R doing `citation(package = 'fauxpas')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
fauxpas 0.5.0
=============

### NEW FEATURES

* changed default for parameter `behavior` in function `http()` to `auto` instead of `stop` (#12)
* `http()` gains new parameter `muffle` (logical) to optionally muffle any warning/message from http status codes in series  (#12)


fauxpas 0.2.0
=============

### NEW FEATURES

* New function `find_error_class()` for a user to get a `HTTP*` class object with just the status code (#10)
* Added a vignette (#4)
* Added an additional object type that fauxpas responds to: the `Response` from the `webmockr` package and the `VcrResponse` from the `vcr` package (#9)

### MINOR IMPROVEMENTS

* `behavior_type` variable in the `Error` class object is now a private variable, and can be set with a new method `$set_behavior()` (#7)
* For `http*()` S3 list methods added additional checks for response objects from `curl` package that are simple lists, so are challenging to verify (#5)
* Added tests for all `http*` methods via a generator (#11)

### BUG FIXES

* Fixed `message_template` behavior, see above related comment about new message template param (#8)


fauxpas 0.1.0
=============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local macOS install, R 3.6.3 patched
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 6 downstream dependencies
(<https://github.com/ropensci/fauxpas/blob/master/revdep/README.md>).
No problems were found.

---

This version includes a couple of new features.

Thanks! Scott Chamberlain
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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/fauxpas/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/fauxpas.git`
* Make sure to track progress upstream (i.e., on our version of `fauxpas` at `ropensci/fauxpas`) by doing `git remote add upstream https://github.com/ropensci/fauxpas.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/fauxpas`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.6.3 Patched (2020-03-11 r78147) |
|os       |macOS Catalina 10.15.4                      |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-04-13                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|fauxpas |0.2.0 |0.5.0 |*  |

# Revdeps

*Wow, no problems at all. :)*# Check times

|   |package |version | check_time|
|:--|:-------|:-------|----------:|
|1  |crul    |0.5.2   |       35.4|
|3  |rorcid  |0.4.0   |       22.8|
|2  |nneo    |0.1.0   |       12.2|


*Wow, no problems at all. :)*