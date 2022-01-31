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


*Wow, no problems at all. :)*fauxpas
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

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

```{r eval=FALSE}
install.packages("fauxpas")
```

Dev version

```{r eval=FALSE}
remotes::install_github("ropensci/fauxpas")
```

```{r}
library("fauxpas")
```

## use with crul

```{r eval=FALSE}
library("crul")
cli <- HttpClient$new("https://httpbin.org/status/414")
res <- cli$get()
http(res)
#> Error: Request-URI Too Long (HTTP 414).
http414(res)
#> Error: Request-URI Too Long (HTTP 414).
```

```{r eval=FALSE}
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

```{r eval=FALSE}
library("curl")
h <- curl::new_handle()
curl::handle_setopt(h)
resp <- curl::curl_fetch_memory("https://httpbin.org/status/404", h)
http(resp)
#> Error: Not Found (HTTP 404).
http404(resp)
#> Error: Not Found (HTTP 404).
```

```{r eval=FALSE}
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

```{r eval=FALSE}
library("httr")
res <- GET("https://httpbin.org/status/405")
http405(res)
#> Error: Method Not Allowed (HTTP 405).
```

```{r eval=FALSE}
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
---
title: "Introduction to the fauxpas package"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to the fauxpas package}
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

`fauxpas` does http errors!

* HTTP error classes more in line with Ruby/Python/Etc. (psssttt: that is a good thing)
* An error class for each HTTP status in case a user wants to
be specific to an HTTP status code, and general purpose handlers
for any error
* Work with any of the major R http clients: `crul`, `curl`, `httr` (maybe
`RCurl` later)
* Provide flexiblity for what to do on an HTTP error, including
custom functions and message templates

To see example of `fauxpas` usage in other packages, check out the code within any of:

- [nneo](https://github.com/ropenscilabs/nneo)
- [rorcid](https://github.com/ropensci/rorcid)

There are __a lot__ of functions in this package. Don't be scared away by that. 
There are a few major sets of functions that are easy to understand and use:

- There's a function `http()` that is a general purpose function that handles any 
HTTP status codes. There's also a set of functions that follow the pattern `http*()`
where `*` is a HTTP status code (e.g. `http201()`). You can use the general 
purpose function `http()` or the one specific to the status code you are interested
in.
- There's a lower level `R6` class `Error` that is wrapped by `http()` that like 
`http()` is a general purpose function that handles any HTTP status codes. 
There's also a set of `R6` classes that follow the pattern `HTTP*()`
where `*` is a HTTP name (e.g. `HTTPGatewayTimeout`). You can use the general 
purpose class `Error()` or the one specific to the status you are interested
in.

`fauxpas` uses [httpcode][] under the hood - which holds all the info about HTTP status 
codes. 

## Install

CRAN version

```{r eval=FALSE}
install.packages("fauxpas")
```

Dev version

```{r eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/fauxpas")
```

```{r}
library("fauxpas")
```

## Find classes for a HTTP status code

When using `fauxpas` in another package I personally like to use the `HTTP*` classes.

Since these are named with their HTTP reason (e.g. `HTTPGatewayTimeout`), it's not super
straightforward to find them with a status code, which is what one always has as a response
from a server. A good way to find these is with a new function `find_error_class()`. Just 
pass the function a HTTP status code and it will find the matching class. 

```{r}
(x <- find_error_class(418))
```

Which returns the matching error class object, with which you can initialize a new object:

```{r}
x$new()
```

If you pass it a status code it doesn't know about it errors

```{r eval=FALSE}
find_error_class(999)
#> Error in find_error_class(999) : no method found for 999
```


## Use fauxpas with HTTP clients

### crul

```{r eval=FALSE}
library("crul")
cli <- HttpClient$new("https://httpbin.org/status/414")
res <- cli$get()
http(res)
#> Error: Request-URI Too Long (HTTP 414).
http414(res)
#> Error: Request-URI Too Long (HTTP 414).
```

```{r eval=FALSE}
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

### curl

```{r eval=FALSE}
library("curl")
h <- curl::new_handle()
curl::handle_setopt(h)
resp <- curl::curl_fetch_memory("https://httpbin.org/status/404", h)
http(resp)
#> Error: Not Found (HTTP 404).
http404(resp)
#> Error: Not Found (HTTP 404).
```

```{r eval=FALSE}
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

### httr

```{r eval=FALSE}
library("httr")
res <- GET("https://httpbin.org/status/405")
http405(res)
#> Error: Method Not Allowed (HTTP 405).
```

```{r eval=FALSE}
x <- HTTPMethodNotAllowed$new()
x$do_verbose(res)
#> Error: Method Not Allowed (HTTP 405).
#>  - The method specified in the Request-Line is not allowed for the resource
#> identified by the Request-URI. The response MUST include an Allow header
#> containing a list of valid methods for the requested resource.
```

[httpcode]: https://github.com/sckott/httpcode
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors-children.R
\name{Error-Classes}
\alias{Error-Classes}
\alias{HTTPContinue}
\alias{HTTPSwitchProtocol}
\alias{HTTPProcessing}
\alias{HTTPOK}
\alias{HTTPCreated}
\alias{HTTPAccepted}
\alias{HTTPNonAuthoritativeInformation}
\alias{HTTPNoContent}
\alias{HTTPResetContent}
\alias{HTTPPartialContent}
\alias{HTTPMultiStatus}
\alias{HTTPAlreadyReported}
\alias{HTTPImUsed}
\alias{HTTPMultipleChoices}
\alias{HTTPMovedPermanently}
\alias{HTTPFound}
\alias{HTTPSeeOther}
\alias{HTTPNotModified}
\alias{HTTPUseProxy}
\alias{HTTPSwitchProxy}
\alias{HTTPTemporaryRedirect}
\alias{HTTPPermanentRedirect}
\alias{HTTPBadRequest}
\alias{HTTPUnauthorized}
\alias{HTTPPaymentRequired}
\alias{HTTPForbidden}
\alias{HTTPNotFound}
\alias{HTTPMethodNotAllowed}
\alias{HTTPNotAcceptable}
\alias{HTTPProxyAuthenticationRequired}
\alias{HTTPRequestTimeout}
\alias{HTTPConflict}
\alias{HTTPGone}
\alias{HTTPLengthRequired}
\alias{HTTPPreconditionFailed}
\alias{HTTPRequestEntityTooLarge}
\alias{HTTPRequestURITooLong}
\alias{HTTPUnsupportedMediaType}
\alias{HTTPRequestRangeNotSatisfiable}
\alias{HTTPExpectationFailed}
\alias{HTTPTeaPot}
\alias{HTTPAuthenticationTimeout}
\alias{HTTPMethodFailure}
\alias{HTTPMisdirectedRequest}
\alias{HTTPUnprocessableEntity}
\alias{HTTPLocked}
\alias{HTTPFailedDependency}
\alias{HTTPUnorderedCollection}
\alias{HTTPUpgradeRequired}
\alias{HTTPPreconditionRequired}
\alias{HTTPTooManyRequests}
\alias{HTTPRequestHeaderFieldsTooLarge}
\alias{HTTPLoginTimeout}
\alias{HTTPNoResponse}
\alias{HTTPRetryWith}
\alias{HTTPBlockedByWindowsParentalControls}
\alias{HTTPUnavailableForLegalReasons}
\alias{HTTPRequestHeaderTooLarge}
\alias{HTTPCertError}
\alias{HTTPNoCert}
\alias{HTTPHTTPToHTTPS}
\alias{HTTPTokenExpiredInvalid}
\alias{HTTPClientClosedRequest}
\alias{HTTPInternalServerError}
\alias{HTTPNotImplemented}
\alias{HTTPBadGateway}
\alias{HTTPServiceUnavailable}
\alias{HTTPGatewayTimeout}
\alias{HTTPHTTPVersionNotSupported}
\alias{HTTPVariantAlsoNegotiates}
\alias{HTTPInsufficientStorage}
\alias{HTTPLoopDetected}
\alias{HTTPBandwidthLimitExceeded}
\alias{HTTPNotExtended}
\alias{HTTPNetworkAuthenticationRequired}
\alias{HTTPNetworkReadTimeoutError}
\alias{HTTPNetworkConnectTimeoutError}
\alias{HTTPWebServerReturnedUnknownError}
\alias{HTTPWebServerIsDown}
\alias{HTTPConnectionTimedOut}
\alias{HTTPOriginIsUnreachable}
\alias{HTTPATimeoutOccurred}
\alias{HTTPSSLHandshakeFailed}
\alias{HTTPInvalidSSLCertificate}
\alias{HTTPRailgunError}
\title{Individual error classes}
\description{
These error classes are for each HTTP error, and inherit from
the \code{\link{Error}} class in this package.
}
\details{
In addition to what's available in \code{\link{Error}},
these classes have a single variable \code{mssg} that is the very
verbose complete message describing the HTTP condition in detail.
You can include that message in your condition by using \code{do_verbose}
(see below)

\strong{Methods}

In addition to the methods documented in \code{\link{Error}}, these methods
also have:
\itemize{
  \item \code{do_verbose(response, template)} {

  Execute condition, whether it be message, warning, or error.

  \itemize{
   \item response: is any response from \pkg{crul}, \pkg{curl}, or \pkg{httr}
  Execute condition, whether it be message, warning, error, or your
  own custom function. This method uses \code{message_template_verbose},
  and uses it's default value.
   \item template: a template to use for the verbose message, see \code{\link{Error}}
  for details
  }
  }
}
}
\examples{
if (requireNamespace("crul")) {

 library("crul")
 res <- HttpClient$new("https://httpbin.org/status/414")$get()
 x <- HTTPRequestURITooLong$new()
 x
 \dontrun{
 x$do(res)
 x$do_verbose(res)
 }

 # behavior
 x <- HTTPRequestURITooLong$new(behavior = "warning")
 \dontrun{
 x$do(res)
 x$do_verbose(res)
 }

 x <- HTTPRequestURITooLong$new(behavior = "message")
 \dontrun{
 x$do(res)
 x$do_verbose(res)
 }

 # with message template
 (x <- HTTPRequestURITooLong$new(
   message_template = "{{reason}} ............ {{status}}",
   message_template_verbose = "{{reason}} .>.>.>.>.>.> {{status}}\n {{message}}"
 ))
 \dontrun{
 x$do(res)
 x$do_verbose(res)
 }
}

}
\seealso{
\code{\link[fauxpas]{Error}}, \code{\link[fauxpas]{http}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error_http_generated.R, R/s3-classes.R
\name{http100}
\alias{http100}
\alias{http101}
\alias{http102}
\alias{http200}
\alias{http201}
\alias{http202}
\alias{http203}
\alias{http204}
\alias{http205}
\alias{http206}
\alias{http207}
\alias{http208}
\alias{http226}
\alias{http300}
\alias{http301}
\alias{http302}
\alias{http303}
\alias{http304}
\alias{http305}
\alias{http306}
\alias{http307}
\alias{http308}
\alias{http400}
\alias{http401}
\alias{http402}
\alias{http403}
\alias{http404}
\alias{http405}
\alias{http406}
\alias{http407}
\alias{http408}
\alias{http409}
\alias{http410}
\alias{http411}
\alias{http412}
\alias{http413}
\alias{http414}
\alias{http415}
\alias{http416}
\alias{http417}
\alias{http418}
\alias{http419}
\alias{http420}
\alias{http422}
\alias{http423}
\alias{http424}
\alias{http425}
\alias{http426}
\alias{http428}
\alias{http429}
\alias{http431}
\alias{http440}
\alias{http444}
\alias{http449}
\alias{http450}
\alias{http451}
\alias{http494}
\alias{http495}
\alias{http496}
\alias{http497}
\alias{http498}
\alias{http499}
\alias{http500}
\alias{http501}
\alias{http502}
\alias{http503}
\alias{http504}
\alias{http505}
\alias{http506}
\alias{http507}
\alias{http508}
\alias{http509}
\alias{http510}
\alias{http511}
\alias{http598}
\alias{http599}
\alias{http}
\title{higher level error wrappers}
\usage{
http100(response, behavior = "auto", message_template, muffle = FALSE)

http101(response, behavior = "auto", message_template, muffle = FALSE)

http102(response, behavior = "auto", message_template, muffle = FALSE)

http200(response, behavior = "auto", message_template, muffle = FALSE)

http201(response, behavior = "auto", message_template, muffle = FALSE)

http202(response, behavior = "auto", message_template, muffle = FALSE)

http203(response, behavior = "auto", message_template, muffle = FALSE)

http204(response, behavior = "auto", message_template, muffle = FALSE)

http205(response, behavior = "auto", message_template, muffle = FALSE)

http206(response, behavior = "auto", message_template, muffle = FALSE)

http207(response, behavior = "auto", message_template, muffle = FALSE)

http208(response, behavior = "auto", message_template, muffle = FALSE)

http226(response, behavior = "auto", message_template, muffle = FALSE)

http300(response, behavior = "auto", message_template, muffle = FALSE)

http301(response, behavior = "auto", message_template, muffle = FALSE)

http302(response, behavior = "auto", message_template, muffle = FALSE)

http303(response, behavior = "auto", message_template, muffle = FALSE)

http304(response, behavior = "auto", message_template, muffle = FALSE)

http305(response, behavior = "auto", message_template, muffle = FALSE)

http306(response, behavior = "auto", message_template, muffle = FALSE)

http307(response, behavior = "auto", message_template, muffle = FALSE)

http308(response, behavior = "auto", message_template, muffle = FALSE)

http400(response, behavior = "auto", message_template, muffle = FALSE)

http401(response, behavior = "auto", message_template, muffle = FALSE)

http402(response, behavior = "auto", message_template, muffle = FALSE)

http403(response, behavior = "auto", message_template, muffle = FALSE)

http404(response, behavior = "auto", message_template, muffle = FALSE)

http405(response, behavior = "auto", message_template, muffle = FALSE)

http406(response, behavior = "auto", message_template, muffle = FALSE)

http407(response, behavior = "auto", message_template, muffle = FALSE)

http408(response, behavior = "auto", message_template, muffle = FALSE)

http409(response, behavior = "auto", message_template, muffle = FALSE)

http410(response, behavior = "auto", message_template, muffle = FALSE)

http411(response, behavior = "auto", message_template, muffle = FALSE)

http412(response, behavior = "auto", message_template, muffle = FALSE)

http413(response, behavior = "auto", message_template, muffle = FALSE)

http414(response, behavior = "auto", message_template, muffle = FALSE)

http415(response, behavior = "auto", message_template, muffle = FALSE)

http416(response, behavior = "auto", message_template, muffle = FALSE)

http417(response, behavior = "auto", message_template, muffle = FALSE)

http418(response, behavior = "auto", message_template, muffle = FALSE)

http419(response, behavior = "auto", message_template, muffle = FALSE)

http420(response, behavior = "auto", message_template, muffle = FALSE)

http422(response, behavior = "auto", message_template, muffle = FALSE)

http423(response, behavior = "auto", message_template, muffle = FALSE)

http424(response, behavior = "auto", message_template, muffle = FALSE)

http425(response, behavior = "auto", message_template, muffle = FALSE)

http426(response, behavior = "auto", message_template, muffle = FALSE)

http428(response, behavior = "auto", message_template, muffle = FALSE)

http429(response, behavior = "auto", message_template, muffle = FALSE)

http431(response, behavior = "auto", message_template, muffle = FALSE)

http440(response, behavior = "auto", message_template, muffle = FALSE)

http444(response, behavior = "auto", message_template, muffle = FALSE)

http449(response, behavior = "auto", message_template, muffle = FALSE)

http450(response, behavior = "auto", message_template, muffle = FALSE)

http451(response, behavior = "auto", message_template, muffle = FALSE)

http494(response, behavior = "auto", message_template, muffle = FALSE)

http495(response, behavior = "auto", message_template, muffle = FALSE)

http496(response, behavior = "auto", message_template, muffle = FALSE)

http497(response, behavior = "auto", message_template, muffle = FALSE)

http498(response, behavior = "auto", message_template, muffle = FALSE)

http499(response, behavior = "auto", message_template, muffle = FALSE)

http500(response, behavior = "auto", message_template, muffle = FALSE)

http501(response, behavior = "auto", message_template, muffle = FALSE)

http502(response, behavior = "auto", message_template, muffle = FALSE)

http503(response, behavior = "auto", message_template, muffle = FALSE)

http504(response, behavior = "auto", message_template, muffle = FALSE)

http505(response, behavior = "auto", message_template, muffle = FALSE)

http506(response, behavior = "auto", message_template, muffle = FALSE)

http507(response, behavior = "auto", message_template, muffle = FALSE)

http508(response, behavior = "auto", message_template, muffle = FALSE)

http509(response, behavior = "auto", message_template, muffle = FALSE)

http510(response, behavior = "auto", message_template, muffle = FALSE)

http511(response, behavior = "auto", message_template, muffle = FALSE)

http598(response, behavior = "auto", message_template, muffle = FALSE)

http599(response, behavior = "auto", message_template, muffle = FALSE)

http(response, behavior = "auto", message_template, muffle = FALSE)
}
\arguments{
\item{response}{The result of a call via \pkg{crul}, \pkg{curl}, or
\pkg{httr}}

\item{behavior}{Behavior of the error. default: auto. See Details}

\item{message_template}{A message template. optional. use whisker
templating. names to use are: reason and status. use in template
like \code{\{\{reason\}\}} and  \code{\{\{status\}\}}. Note that
\code{\{\{message\}\}} that is used in \code{message_template_verbose}
will be ignored here.}

\item{muffle}{(logical) whether to not respond when status codes
in 1xx-3xx series. Default: \code{FALSE}}
}
\description{
higher level error wrappers
}
\note{
These \code{http*} methods only use \code{$do} and not
\code{$do_verbose}.
}
\section{behavior parameter options}{

\itemize{
 \item stop - use \code{stop}
 \item warning - use \code{warning}
 \item message - use \code{message}
 \item auto - toggle between \code{stop} and \code{message} depending
 on the HTTP status code series. Defaults will be:
 \itemize{
  \item 1xx: \code{message}
  \item 2xx: \code{message}
  \item 3xx: \code{message}
  \item 4xx: \code{stop}
  \item 5xx: \code{stop}
 }
}

Of course, you can always override the defaults.
}

\section{using package \pkg{curl}}{

curl reponses are simple lists, so we have little to go on to make sure
it's a response from the \pkg{curl} package. We check for list names
internally but of course you could pass in a list with the right named
elements, while the values are complete nonsense, in which case
we'll probably fail badly. There's not much we can do.
}

\examples{
if (requireNamespace("crul")) {
 library("crul")
 res <- HttpClient$new("https://httpbin.org/status/418")$get()
 \dontrun{http(res)}
 http(res, behavior = "warning")
 http(res, behavior = "message")

 res <- HttpClient$new("https://httpbin.org/status/414")$get()
 \dontrun{http414(res)}
 http(res, behavior = "warning")
 http(res, behavior = "message")

 res <- HttpClient$new("https://httpbin.org/asdfafadsf")$get()
 \dontrun{http404(res)}
 http(res, behavior = "warning")
 http(res, behavior = "message")
}

if (requireNamespace("curl")) {
 library("curl")
 h <- curl::new_handle()
 curl::handle_setopt(h)
 res <- curl::curl_fetch_memory("https://httpbin.org/status/418", h)
 \dontrun{http(res)}
 http(res, behavior = "warning")
 http(res, behavior = "message")
}

if (requireNamespace("httr")) {
 library("httr")
 res <- GET("https://httpbin.org/status/418")
 \dontrun{http(res)}
 http(res, behavior = "warning")
 http(res, behavior = "message")
}

# muffle responses
if (requireNamespace("crul")) {
 library("crul")
 res201 <- HttpClient$new("https://httpbin.org")$get("status/201")
 res404 <- HttpClient$new("https://httpbin.org")$get("status/404")
 # status codes < 300 CAN be muffled - i.e., return the http response object
 http(res201, muffle = TRUE)
 # status codes > 300 CAN NOT be muffled - i.e., return the http response object
 \dontrun{
 http(res404, muffle = TRUE)
 }
}
}
\seealso{
\code{\link{Error}}, \code{\link[fauxpas]{Error-Classes}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_error_class.R
\name{find_error_class}
\alias{find_error_class}
\title{Find error classes}
\usage{
find_error_class(status_code)
}
\arguments{
\item{status_code}{(numeric,integer) A status code}
}
\value{
an object of class \code{R6ClassGenerator}. call \code{$new()}
to initialize a new instance
}
\description{
Find error classes
}
\examples{
find_error_class(414)
find_error_class(418)
find_error_class(505)

# initialize the class
find_error_class(418)$new()

# not found
\dontrun{find_error_class(999)}
}
\seealso{
\code{\link[fauxpas]{Error}}, \code{\link[fauxpas]{Error-Classes}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{Error}
\alias{Error}
\title{Error class}
\arguments{
\item{behavior}{Behavior of the error. default: auto. See Details}

\item{message_template}{A message template. optional. use whisker
templating. names to use are: reason and status. use in template
like \code{\{\{reason\}\}} and  \code{\{\{status\}\}}. Note that
\code{\{\{message\}\}} that is used in \code{message_template_verbose}
will be ignored here.}

\item{call.}{(logical) indicating if the call should become part
of the error message. Default: \code{FALSE}}

\item{message_template_verbose}{A verbose message template. optional.
use whisker templating. names to use are: reason, status, message.
use in template like \code{\{\{reason\}\}}, \code{\{\{status\}\}}, and
\code{\{\{message\}\}}. Note that this is ignored here, but is
used in the \code{HTTP*} methods (e.g. \code{HTTPBadRequest})}

\item{muffle}{(logical) whether to not respond when status codes
in 1xx-3xx series. Default: \code{FALSE}}
}
\description{
Error class
}
\details{
\strong{Methods}
\itemize{
  \item \code{do(response, mssg)} {

  Execute condition, whether it be message, warning, or error.

  \itemize{
   \item response: is any response from \pkg{crul}, \pkg{curl}, or \pkg{httr}
  Execute condition, whether it be message, warning, error, or your
  own custom function. This method uses \code{message_template_verbose},
  and uses it's default value.
   \item mssg: character string message to include in call. ignored if
  template does not have a \code{message} entry
   }
  }

  \item \code{set_behavior(behavior)}

  Set behavior, same as setting behavior on initializing with \code{$new()}
}
}
\section{behavior parameter options}{

\itemize{
 \item stop - use \code{stop}
 \item warning - use \code{warning}
 \item message - use \code{message}
 \item auto - toggle between \code{stop} and \code{message} depending
 on the HTTP status code series. Defaults will be:
 \itemize{
  \item 1xx: \code{message}
  \item 2xx: \code{message}
  \item 3xx: \code{message}
  \item 4xx: \code{stop}
  \item 5xx: \code{stop}
 }
}

Of course, you can always override the defaults.
}

\examples{
Error$new()
# reset behavior
(z <- Error$new())
z$set_behavior("warning")
z

if (requireNamespace("crul")) {
 library("crul")
 res <- HttpClient$new("https://httpbin.org/status/418")$get()

 # stop
 (x <- Error$new(behavior = "stop"))
 \dontrun{x$do(res)}

 # warn
 (x <- Error$new(behavior = "warning"))
 x$do(res)

 # do vs. do_verbose
 x <- HTTPRequestURITooLong$new(behavior = "stop")
 res <- HttpClient$new("https://httpbin.org/status/414")$get()
 \dontrun{
 http414(res)
 ## with template
 http414(res, message_template = "{{status}}\n  --> {{reason}}")
 x$do(res)
 x$do_verbose(res)
 }

 # service unavailable
 x <- HTTPServiceUnavailable$new(behavior = "stop")
 res <- HttpClient$new("https://httpbin.org/status/503")$get()
 \dontrun{
 x$do(res)
 x$do_verbose(res)
 }

 # message template
 y <- Error$new(message_template = "{{reason}} ............ {{status}}")
 res <- HttpClient$new("https://httpbin.org/status/418")$get()
 \dontrun{
 y$do(res)
 }

 yy <- Error$new(message_template = "{{status}}\n  --> {{reason}}")
 yy$message_template
 \dontrun{
 yy$do(res)
 }

 ## with verbose message
 library(crul)
 res <- HttpClient$new("https://httpbin.org/status/401")$get()
 yy <- HTTPUnauthorized$new()
 zz <- HTTPUnauthorized$new(
   message_template = "HTTP({{status}}): {{reason}}\n  {{message}}"
 )
 yy$message_template; zz$message_template
 \dontrun{
 yy$do(res)
 zz$do(res)
 yy$do_verbose(res)
 zz$do_verbose(res)
 }

 yy <- Error$new(
   message_template = "HTTP({{status}}): {{reason}}\n  {{message}}"
 )
 yy$message_template
 \dontrun{yy$do(res)}

 # muffle responses
 (x <- Error$new(muffle = TRUE))
 res <- crul::HttpClient$new("https://httpbin.org/status/226")$get()
 z <- x$do(res)
 z
}
}
\seealso{
\code{\link[fauxpas]{http}}, \code{\link[fauxpas]{Error-Classes}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fauxpas-package.R
\docType{package}
\name{fauxpas-package}
\alias{fauxpas-package}
\alias{fauxpas}
\title{fauxpas}
\description{
HTTP Error Helpers
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
