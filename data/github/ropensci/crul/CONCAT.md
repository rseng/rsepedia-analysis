# crul <img src="man/figures/logo.png" align="right" alt="" width="120">



[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/crul/workflows/R-check/badge.svg)](https://github.com/ropensci/crul/actions/)
[![codecov](https://codecov.io/gh/ropensci/crul/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/crul)
[![cran checks](https://cranchecks.info/badges/worst/crul)](https://cranchecks.info/pkgs/crul)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/crul)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/crul)](https://cran.r-project.org/package=crul)

An HTTP client, taking inspiration from Ruby's [faraday](https://rubygems.org/gems/faraday) and Python's `requests`

Package documentation: <https://docs.ropensci.org/crul/>

Some Features:

* `HttpClient` - Main interface to making HTTP requests. Synchronous requests only.
* `HttpResponse` - HTTP response object, used for all responses across the
different clients.
* `Paginator` - Auto-paginate through requests - supports a subset of all possible
pagination scenarios - will fill out more scenarios soon
* `Async` - Asynchronous HTTP requests - a simple interface for many URLS -
whose interface is similar to `HttpClient` - all URLs are treated the same.
* `AsyncVaried` - Asynchronous HTTP requests - accepts any number of `HttpRequest`
objects - with a different interface than `HttpClient`/`Async` due to the nature
of handling requests with different HTTP methods, options, etc.
* set curl options globally: `set_auth()`, `set_headers()`, and more
* Writing to disk and streaming: available with both synchronous requests
as well as async requests
* Hooks on requests and responses are available in the `HttpClient` method only, 
and allow you to trigger functions to run on requests or responses, or both.
See `?hooks` for the details and examples
* Mocking: `crul` integrates with [webmockr](https://github.com/ropensci/webmockr) to mock
HTTP requests. Checkout the [http testing book][book]
* Test caching: `crul` also integrates with [vcr](https://github.com/ropensci/vcr) to cache http requests/responses. Checkout the [http testing book][book]

## Installation

CRAN version


```r
install.packages("crul")
```

Latest binaries from rOpenSci


```r
install.packages("crul", repos = "https://dev.ropensci.org")
```

Dev version from GitHub


```r
remotes::install_github("ropensci/crul")
```


```r
library("crul")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/crul/issues).
* License: MIT
* Get citation information for `crul` in R doing `citation(package = 'crul')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
* Where does the package name come from? It was a play on "curl", the popular command line client.
* Where does the sticker design come from? The sticker idea [arose from a tweet](https://github.com/ropensci/crul/issues/42) - crul is close (ish) to Krull, a 1980's movie with a "mystical five-pointed weapon". The association with Krull was not known before naming the package.

[book]: https://books.ropensci.org/http-testing/
crul 1.2
========

### DOCUMENTATION

* fix example in `AsyncQueue` docs (#146) thanks @johnbaums !
* update `HttpClient` docs to state that it's an R6 class, and give some details on what an R6 class is and links to more info (#155)

### NEW FEATURES

* `AsyncQueue` gains methods: `parse`, `status_code`, `status`, `content`, and `times` (#156)
* `$responses()` method now returns an S3 class with an associated print method to prevent printing a lot of results to the screen; print method pritns a summary of results, and at most 10 results, just status code and url (#157)


### MINOR IMPROVEMENTS

* parsing response headers gains a check for whether encoding is valid, and if not tries to set Latin1 encoding, and if that doesn't work, fails out with message (#163) (#164) thanks @FlukeAndFeather


crul 1.1
========

### NEW FEATURES

* `Paginator` gains support for query parameter combination `page`/`per_page`  to automatically paginate (#145)

### MINOR IMPROVEMENTS

* fix typo (#149) thanks @dpprdan
* Change to how numbers are handled in query parameters. We unfortunately hadn't tested this package with large numbers, which were being converted to scientific notation with a certain number of digits before a decimal. Fixed handling of query parameters to avoid this problem. Fix for `Paginator` as well as for `HttpClient` (#151) (#152) (#153) thanks @ateucher

### BUG FIXES

* sometimes weird response headers are returned in an HTTP response that can not be easily parsed; `crul` would raise an error when this header parsing happens, but now we raise a warning instead (#150)


crul 1.0
========

### ok related changes

* `ok()` can now accept more than 1 status code so that you can check if the status of a url is within a set of status codes rather than equal to 1 status code (#124)
* `ok()` gains a parameter `verb` to use either head or get requests. in addition added more documentation (#125) to the function on how to get the "right answer" for whether a url is ok/up  (#123) (#127)
* `ok()` gains parameter `ua_random`, which if `TRUE`, will use a random user agent string pulled from a vector of 50 user agent strings generated from `charlatan::UserAgentProvider` (#138)

### NEW FEATURES

* gains new async class `AsyncQueue` for doing async requests with rate limits (#139)
* gains new functions `curl_verbose()` and `set_verbose()`. `curl_verbose()` can be set by passing to the initialize step (e.g., `HttpClient$new(url, verbose=curl_verbose())`), and gets more compact verbose curl output, while also getting request body information (and response body optionally). `set_verbose()` is sets `curl_verbose()` globally (#141)
* gains new vignette "How to choose a client" for choosing which crul class to use (e.g., `HttpClient` vs. `Async`)  (#133) (#143)

### MINOR IMPROVEMENTS

* package sticker done, shown in README (#42)
* improve function/class reference page in docs site (#131)
* improvements to the best practices vignette (#132)
* removed unused private variable in the `AsyncVaried` class (#140)
* fix inaccuracy in documentation for the RETRY method (#130)
* `HttpRequest` now adds the query (if present) to the printed url in the print method for the class (it was absent before now) (#128)
* use new roxygen2 support for R6 classes (#126)
* removed `delete-requests` and `post-requests` manual files - mostly redundant with other documentation


crul 0.9.0
==========

### NEW FEATURES

* `HttpResponse` response object gains new methods for checking response content types, includes: `raise_for_ct`, `raise_for_ct_html`, `raise_for_ct_json`, `raise_for_ct_xml`. these behave similarly to `raise_for_status`, and can behave as a warning or raise an error through stop (#119) (#120)

### MINOR IMPROVEMENTS

* fix to prep_body internal function to handle various body inputs; now avoids warning about `as.character.form_file` when both httr and crul are loaded (#112)
* finish off "Failing with fauxpas" section of the "API package best practices" vignette (#121)

### BUG FIXES

* the `head()` verb on `HttpClient` was no capturing `auth` when set on initialization (#122)


crul 0.8.4
==========

### MINOR IMPROVEMENTS

* `jsonlite` package moved to Imports (#112)
* the `parse()` method in the `HttpResponse` object now checks whether the response raw bytes can be converted to character, and if not just returns raw bytes (#115) (#116)
* give vignettes titles (#113) (#114)

### BUG FIXES

* no longer setting `cainfo` curl option, fixes problem arising from change in recent libcurl version (#117)


crul 0.8.0
==========

### NEW FEATURES

* you can now pass on parameters through the `parse()` method of an `HttpResponse` class to the internally called function `iconv()` to finely control the usage of `iconv` for cases in which normal encoding conversion doesn't work  (#110)

### MINOR IMPROVEMENTS

* use `httpcode` package instead of `fauxpas` package within `ok()` function (#108) (#109) thanks @maelle !
* fix links to http testing book - ropensci -> ropenscilabs (#111)


crul 0.7.4
==========

### NEW FEATURES

* event hooks added to `HttpClient`. both request and response hooks supported. not supported in async methods for now (#76) (#107)

### MINOR IMPROVEMENTS

* improve `$parse()` behavior (in the `HttpResponse` object) when using disk or stream. `$parse()` was throwing a warning when using disk and an error when using stream. and improves behavior when doing async requests (#104)
* `Paginator` gains optional progress bar through the new `progress` parameter. In addition, the `cat()` calls inside the method were removed, so as not to insert newlines with each page and to not print "OK" when done (#106) thanks @boshek

### BUG FIXES

* passing on opts/headers now works with `Async` (#101) (#103)
* streaming was broken in `AsyncVaried` with `curl` of a certain version, works now (#102) (#103)


crul 0.7.0
==========

### NEW FEATURES

* `HttpClient` gains a `retry` method: retries any request verb until successful (HTTP response status < 400) or a condition for giving up is met. (#89) (#95) thanks @hlapp
* `HttpClient`, `HttpRequest`, and `Async` classes gain `verb` method for doing HTTP requests specifying any of the supported HTTP verbs  (#97)
* `HttpClient` and `Paginator` gain a `url_fetch` method: get the URL that would be sent in an HTTP request without sending the HTTP request. Useful for getting the URL before executing an HTTP request if you need to check something about the URL first. (#92)
* new vignette for "API package best practices" (#65)
* Package gains manual files for each HTTP verb to facilitate linking to package documentation for information on each HTTP verb (#98)
* Intermediate headers (e.g., those in redirect chains) are now given back in a new slot in the `HttpResponse` class as `$response_headers_all` as an unnamed list, with each element a named list of headers; the last list in the set is the final response headers that match those given in the `$response_headers` slot (#60) (#99)

### BUG FIXES

* some dangling file connections were left open - now fixed (#93) (#95)
* fix `url_parse`: lacked check that input was a string, and that it was length 1 - this PR fixed that (#100) thanks @aaronwolen

### DEFUNCT

* `HttpStubbedResponse` was removed from the package - it may have been used at some point, but is not used in the package anymore (#88)

crul 0.6.0
==========

### NEW FEATURES

* `Async` and `AsyncVaried` now support simple auth, see `?auth`  (#70)
* gains new function `ok()` to ping a URL to see if it's up or not, returns a single boolean (#71) (#73)
* `HttpClient` and `HttpRequest` gain new parameter `progress` that accepts a function to use to construct a progress bar. For now accepts `httr::progress()` but will accept other options in the future (#20) (#81)
* gains a new vignette for curl options (#7)
* can now set curl options globally using new functions `set_auth()`, `set_headers()`, `set_opts()`, `set_proxy()`, and `crul_settings()` (#48) (#85)

### MINOR IMPROVEMENTS

* explicitly import `httpcode::http_code` (#80)
* fix vignette names to make them more clear and add numbers to order them (#64)
* change print function for `Async` and `AsyncVaried` to print max of 10 and tell user how many total and remaining not shown (#72)
* added support to `proxy()` for socks, e.g. to use with TOR (#79)
* now when `Async` and `AsyncVaried` requests fail, they don't error but instead we capture the error and pass it back in the result. this way any failure requests don't stop progress of the entire async request suite (#74) (#84)


crul 0.5.2
==========

### MINOR IMPROVEMENTS

* Fixed handling of user agent: you can pass a UA string 
as a curl option or a header. Previously, we were wrongly overwriting
the user input UA if given as a curl option - but were not doing 
so if given as a header. This is fixed now.  (#63) thx to @maelle and @dpprdan

### BUG FIXES

* Fix to `Paginator` - it wasn't handling pagination correctly. 
In addition, fixed to hopefully handle all scenarios now. added
more tests (#62)
* Fixed handling of query parameters. We were using `urltools::url_encode`
to encode strings, but it wasn't encoding correctly in some locales. Using `curl::curl_escape` fixes the problem. Encoding is done on query values and names  (#67) (#68)

crul 0.5.0
==========

### NEW FEATURES

* Gains a new R6 class `Paginator` to help users automatically paginate through multiple requests. It only supports query parameter based paginating for now. We'll add support later for other types including cursors (e.g., used in Solr servers), and for link headers (e.g., used in the GitHub API). Please get in touch if you find any problems with `Paginator`. (#56)
* Async classes `Async` and `Asyncvaried` gain ability to write to disk and stream data (to disk or elsewhere, e.g. R console or to an R object) (#46) thanks @artemklevtsov for the push to do this

### MINOR IMPROVEMENTS

* Improved documentation for `auth` to indicate that `user` and `pwd` are indeed required - and to further indicate that one can pass in `NULL` to those parameters (similar to an empty string `""` in `httr::authenticate`) when one e.g. may want to use `gssnegotiate` method (#43)
* Fixed query builder so that one can now protect query parameters by wrapping them in `I()` (#55)

### BUG FIXES

* Fixed bug in `head` requests with `HttpClient` when passing `query` parameter - it was failing previously. Added `query` parameter back. (#52)


crul 0.4.0
==========

### NEW FEATURES

* file uploads now work, see new function `upload()` and examples (#25)

### MINOR IMPROVEMENTS

* fixes to reused curl handles - within a connection object only,
not across connection objects (#45)
* `crul` now drops any options passed in to `opts` or to `...` that 
are not in set of allowed curl options, see `curl::curl_options()` (#49)
* cookies should now be persisted across requests within 
a connection object, see new doc `?cookies` for how to set cookies (#44)
* gather cainfo and use in curl options when applicable (#51)
* remove `disk` and `stream` from `head` method in `HttpClient` 
and `HttpRequest` as no body returned in a HEAD request


crul 0.3.8
==========

### BUG FIXES

* Fixed `AsyncVaried` to return async responses in the order that
they were passed in. This also fixes this exact same behavior in 
`Async` because `Async` uses `AsyncVaried` internally. (#41)
thanks @dirkschumacher for reporting



crul 0.3.6
==========

* Note: This version gains support for integration with 
`webmockr`, which is now on CRAN.

### NEW FEATURES

* New function `auth()` to do simple authentication (#33)
* New function `HttpStubbedResponse` for making a stubbed 
response object for the `webmockr` integration (#4)
* New function `mock()` to turn on mocking - it's off by 
default. If `webmockr` is not installed but user attempts 
to use mocking we error with message to install 
`webmockr` (#4)

### MINOR IMPROVEMENTS

* Use `gzip-deflate` by deafult for each request 
to make sure gzip compression is used if the server 
can do it (#34)
* Change `useragent` to `User-Agent` as default user 
agent header (#35)
* Now we make sure that user supplied headers override the 
default headers if they are of the same name (#36)



crul 0.3.4
==========

### NEW FEATURES

* New utility functions `url_build` and `url_parse` (#31)

### MINOR IMPROVEMENTS

* Now using markdown for documentation (#32)
* Better documentation for `AsyncVaried` (#30)
* New vignette on how to use `crul` in realistic
scenarios rather than brief examples to demonstrate
individual features (#29)
* Better documentation for `HttpRequest` (#28)
* Included more tests

### BUG FIXES

* Fixed put/patch/delete as weren't passing body
correctly in `HttpClient` (#26)
* DRY out code for preparing requests - simplify to
use helper functions (#27)


crul 0.3.0
==========

### NEW FEATURES

* Added support for asynchronous HTTP requests, including two new
R6 classes: `Async` and `AsyncVaried`. The former being a simpler
interface treating all URLs with same options/HTTP method, and the latter
allowing any type of request through the new R6 class `HttpRequest` (#8) (#24)
* New R6 class `HttpRequest` to support `AsyncVaried` - this method
only defines a request, but does not execute it. (#8)

### MINOR IMPROVEMENTS

* Added support for proxies (#22)

### BUG FIXES

* Fixed parsing of headers from FTP servers (#21)





crul 0.2.0
==========

### MINOR IMPROVEMENTS

* Created new manual files for various tasks to document
usage better (#19)
* URL encode paths - should fix any bugs where spaces between words
caused errors previously (#17)
* URL encode query parameters - should fix any bugs where spaces between words
caused errors previously (#11)
* request headers now passed correctly to response object (#13)
* response headers now parsed to a list for easier access (#14)
* Now supporting multiple query parameters of the same name, wasn't
possible in last version (#15)




crul 0.1.6
==========

### NEW FEATURES

* Improved options for using curl options. Can manually add
to list of curl options or pass in via `...`. And we
check that user doesn't pass in prohibited options
(`curl` package takes care of checking that options
are valid) (#5)
* Incorporated `fauxpas` package for dealing with HTTP
conditions. It's a Suggest, so only used if installed (#6)
* Added support for streaming via `curl::curl_fetch_stream`.
`stream` param defaults to `NULL` (thus ignored), or pass in a
function to use streaming. Only one of memory, streaming or
disk allowed. (#9)
* Added support for streaming via `curl::curl_fetch_disk`.
`disk` param defaults to `NULL` (thus ignored), or pass in a
path to write to disk instead of use memory. Only one of memory,
streaming or disk allowed. (#12)

### MINOR IMPROVEMENTS

* Added missing `raise_for_status()` method on the
`HttpResponse` class (#10)

### BUG FIXES

* Was importing `httpcode` but wasn't using it in the package.
Now using the package in `HttpResponse`






crul 0.1.0
==========

### NEW FEATURES

* Released to CRAN.
## Test environments

* local macOS install, R 4.1.2
* ubuntu 18.04 (on github actions), R 4.1.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 70 downstream dependencies. No problems were found related to this package.

---

This version includes some minor improvements.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/crul/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/crul.git`
* Make sure to track progress upstream (i.e., on our version of `crul` at `ropensci/crul`) by doing `git remote add upstream https://github.com/ropensci/crul.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/crul`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-06-30 r78761) |
|os       |macOS Catalina 10.15.6                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-07-30                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|crul    |0.9.0 |1.0.0 |*  |

# Revdeps

## Failed to check (3)

|package                              |version |error  |warning |note |
|:------------------------------------|:-------|:------|:-------|:----|
|[microdemic](failures.md#microdemic) |0.5.0   |__+1__ |        |     |
|mindicador                           |0.1.5   |1      |        |     |
|mstrio                               |?       |       |        |     |

# microdemic

<details>

* Version: 0.5.0
* Source code: https://github.com/cran/microdemic
* URL: https://github.com/ropensci/microdemic (devel), https://docs.ropensci.org/microdemic (website)
* BugReports: https://github.com/ropensci/microdemic/issues
* Date/Publication: 2020-01-27 22:00:02 UTC
* Number of recursive dependencies: 50

Run `revdep_details(,"microdemic")` for more info

</details>

## Newly broken

*   R CMD check timed out
    

# Check times

|   |package          |version | check_time|
|:--|:----------------|:-------|----------:|
|29 |ropenaq          |0.2.5   |      191.7|
|28 |rnoaa            |0.7.0   |      119.2|
|7  |finch            |0.2.0   |       97.3|
|31 |rplos            |0.8.0   |       94.8|
|8  |fulltext         |1.0.1   |       88.8|
|38 |spocc            |0.7.0   |       80.6|
|20 |rcrossref        |0.8.0   |       77.1|
|34 |rvertnet         |0.6.2   |       76.8|
|35 |seaaroundus      |1.2.0   |       76.2|
|2  |ccafs            |0.1.0   |       64.4|
|37 |solrium          |1.0.0   |       63.7|
|33 |rtimes           |0.5.0   |       63.3|
|39 |traits           |0.3.0   |         63|
|24 |rgbif            |0.9.9   |       51.9|
|36 |sofa             |0.3.0   |       46.2|
|21 |rdatacite        |0.3.0   |       45.9|
|18 |rbison           |0.5.4   |       43.9|
|3  |crminer          |0.1.4   |       43.6|
|41 |wikitaxa         |0.2.0   |       42.9|
|30 |rorcid           |0.4.0   |       42.6|
|32 |rredlist         |0.4.0   |         37|
|6  |fauxpas          |0.1.0   |       35.5|
|23 |rdryad           |0.3.0   |       33.7|
|22 |rdpla            |0.2.0   |       33.2|
|40 |webmockr         |0.1.0   |       33.2|
|13 |openadds         |0.2.0   |       31.7|
|26 |ritis            |0.7.0   |       31.5|
|17 |rbhl             |0.8.0   |       28.9|
|14 |pangaear         |0.6.0   |       28.8|
|10 |microdemic       |0.2.0   |       26.2|
|1  |bold             |0.5.0   |         26|
|9  |jaod             |0.1.0   |       25.6|
|27 |rjsonapi         |0.1.0   |       22.9|
|19 |rcoreoa          |0.1.0   |       22.6|
|11 |natserv          |0.1.4   |       21.5|
|42 |worrms           |0.2.0   |       21.1|
|4  |discgolf         |0.2.0   |       18.3|
|5  |duckduckr        |1.0.0   |       17.7|
|25 |rif              |0.2.0   |       17.7|
|12 |nneo             |0.1.0   |       17.6|
|15 |pleiades         |0.2.0   |       16.8|
|16 |postlightmercury |1.2     |       13.8|


Hi,

This is an automated email to let you know about the release of {{{ my_package }}}, which I'll submit to CRAN on {{{ date }}}.

To check for potential problems, I ran `R CMD check` on your package {{{your_package}}} (v{{{your_version}}}).

I found: {{{your_summary}}}.

{{#you_have_problems}}
{{{your_results}}}

If I got an ERROR because I couldn't install your package (or one of it's dependencies), my apologies. You'll have to run the checks yourself (unfortunately I don't have the time to diagnose installation failures as I have to run checks on hundreds of packages).

Otherwise, please carefully look at the results, and let me know if I've introduced a bug in {{{ my_package }}}. If I have, I'd really appreciate a minimal reproducible example that uses only {{{ my_package }}} functions. That way I can find and fix the bug as quickly as possible.

If it doesn't look like a bug in {{{ my_package }}}, please prepare an update for CRAN. Ideally you'll tweak your package so it works with both the released and development versions of crul. Otherwise, be prepared to submit your package to CRAN soon after I let you know that I've submitted.

To get the development version of {{{ my_package }}} so you can run the checks yourself, you can run:

    # install.packages("devtools")
    devtools::install_github("{{my_github}}")

To see what's changed visit <https://github.com/{{{my_github}}}/blob/main/NEWS.md>.

{{/you_have_problems}}
{{^you_have_problems}}
It looks like everything is ok, so you don't need to take any action, but you might want to read the NEWS, <https://github.com/{{{my_github}}}/blob/main/NEWS.md>, to see what's changed.
{{/you_have_problems}}


If you have any questions about this email, please feel free to respond directly.

Regards,

{{{ me }}}
# microdemic

<details>

* Version: 0.5.0
* Source code: https://github.com/cran/microdemic
* URL: https://github.com/ropensci/microdemic (devel), https://docs.ropensci.org/microdemic (website)
* BugReports: https://github.com/ropensci/microdemic/issues
* Date/Publication: 2020-01-27 22:00:02 UTC
* Number of recursive dependencies: 50

Run `revdep_details(,"microdemic")` for more info

</details>

## Newly broken

*   R CMD check timed out
    

# mindicador

<details>

* Version: 0.1.5
* Source code: https://github.com/cran/mindicador
* URL: https://github.com/pachamaltese/mindicador
* BugReports: https://github.com/pachamaltese/mindicador/issues
* Date/Publication: 2020-06-02 09:30:02 UTC
* Number of recursive dependencies: 84

Run `revdep_details(,"mindicador")` for more info

</details>

## In both

*   R CMD check timed out
    

# mstrio

<details>

* Version: 
* Source code: ???
* URL: https://docs.ropensci.org/crul (website) https://github.com/ropensci/crul (devel) https://books.ropensci.org/http-testing/ (user manual)
* BugReports: https://github.com/ropensci/crul/issues
* Number of recursive dependencies: 0

Run `revdep_details(,"")` for more info

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
