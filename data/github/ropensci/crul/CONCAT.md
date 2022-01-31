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
# crul <img src="man/figures/logo.png" align="right" alt="" width="120">

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE
)
```

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

```{r eval=FALSE}
install.packages("crul")
```

Latest binaries from rOpenSci

```{r eval=FALSE}
install.packages("crul", repos = "https://dev.ropensci.org")
```

Dev version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/crul")
```

```{r}
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
---
title: 3. async with crul
author: Scott Chamberlain
date: "2020-07-09"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{3. async with crul}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



Asynchronous requests with `crul`.

There are two interfaces to asynchronous requests in `crul`:

1. Simple async: any number of URLs, all treated with the same curl options,
headers, etc., and only one HTTP method type at a time.
2. Varied request async: build any type of request and execute all asynchronously.

The first option takes less thinking, less work, and is good solution when you
just want to hit a bunch of URLs asynchronously.

The second option is ideal when you want to set curl options/headers on each
request and/or want to do different types of HTTP methods on each request.

One thing to think about before using async is whether the data provider is
okay with it. It's possible that a data provider's service may be brought down
if you do too many async requests.


```r
library("crul")
```

## simple async

Build request object with 1 or more URLs




```r
(cc <- Async$new(
  urls = c(
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?a=5&b=6',
    'https://httpbin.org/ip'
  )
))
#> <crul async connection> 
#>   curl options: 
#>   proxies: 
#>   auth: 
#>   headers: 
#>   urls: (n: 3)
#>    https://httpbin.org/get?a=5
#>    https://httpbin.org/get?a=5&b=6
#>    https://httpbin.org/ip
```

Make request with any HTTP method


```r
(res <- cc$get())
#> [[1]]
#> <crul response> 
#>   url: https://httpbin.org/get?a=5
#>   request_headers: 
#>   response_headers: 
#>     status: HTTP/2 200 
#>     date: Thu, 09 Jul 2020 20:14:45 GMT
#>     content-type: application/json
#>     content-length: 400
#>     server: gunicorn/19.9.0
#>     access-control-allow-origin: *
#>     access-control-allow-credentials: true
#>   params: 
#>     a: 5
#>   status: 200
#> 
#> [[2]]
#> <crul response> 
#>   url: https://httpbin.org/get?a=5&b=6
#>   request_headers: 
#>   response_headers: 
#>     status: HTTP/2 200 
#>     date: Thu, 09 Jul 2020 20:14:45 GMT
#>     content-type: application/json
#>     content-length: 419
#>     server: gunicorn/19.9.0
#>     access-control-allow-origin: *
#>     access-control-allow-credentials: true
#>   params: 
#>     a: 5
#>     b: 6
#>   status: 200
#> 
#> [[3]]
#> <crul response> 
#>   url: https://httpbin.org/ip
#>   request_headers: 
#>   response_headers: 
#>     status: HTTP/2 200 
#>     date: Thu, 09 Jul 2020 20:14:45 GMT
#>     content-type: application/json
#>     content-length: 31
#>     server: gunicorn/19.9.0
#>     access-control-allow-origin: *
#>     access-control-allow-credentials: true
#>   status: 200
```

You get back a list matching length of the number of input URLs

Access object variables and methods just as with `HttpClient` results, here just one at a time.


```r
res[[1]]$url
#> [1] "https://httpbin.org/get?a=5"
res[[1]]$success()
#> [1] TRUE
res[[1]]$parse("UTF-8")
#> [1] "{\n  \"args\": {\n    \"a\": \"5\"\n  }, \n  \"headers\": {\n    \"Accept\": \"application/json, text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\": \"R (4.0.2 x86_64-apple-darwin17.0 x86_64 darwin17.0)\", \n    \"X-Amzn-Trace-Id\": \"Root=1-5f077ab5-686fd087166386a10fe80f25\"\n  }, \n  \"origin\": \"24.21.229.59\", \n  \"url\": \"https://httpbin.org/get?a=5\"\n}\n"
```

Or apply access/method calls across many results, e.g., parse all results


```r
lapply(res, function(z) z$parse("UTF-8"))
#> [[1]]
#> [1] "{\n  \"args\": {\n    \"a\": \"5\"\n  }, \n  \"headers\": {\n    \"Accept\": \"application/json, text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\": \"R (4.0.2 x86_64-apple-darwin17.0 x86_64 darwin17.0)\", \n    \"X-Amzn-Trace-Id\": \"Root=1-5f077ab5-686fd087166386a10fe80f25\"\n  }, \n  \"origin\": \"24.21.229.59\", \n  \"url\": \"https://httpbin.org/get?a=5\"\n}\n"
#> 
#> [[2]]
#> [1] "{\n  \"args\": {\n    \"a\": \"5\", \n    \"b\": \"6\"\n  }, \n  \"headers\": {\n    \"Accept\": \"application/json, text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\": \"R (4.0.2 x86_64-apple-darwin17.0 x86_64 darwin17.0)\", \n    \"X-Amzn-Trace-Id\": \"Root=1-5f077ab5-0c5caed1b0cac8d405458e4e\"\n  }, \n  \"origin\": \"24.21.229.59\", \n  \"url\": \"https://httpbin.org/get?a=5&b=6\"\n}\n"
#> 
#> [[3]]
#> [1] "{\n  \"origin\": \"24.21.229.59\"\n}\n"
```

## varied request async


```r
req1 <- HttpRequest$new(
  url = "https://httpbin.org/get?a=5",
  opts = list(
    verbose = TRUE
  )
)
req1$get()
#> <crul http request> get
#>   url: https://httpbin.org/get?a=5
#>   curl options: 
#>     verbose: TRUE
#>   proxies: 
#>   auth: 
#>   headers: 
#>   progress: FALSE

req2 <- HttpRequest$new(
  url = "https://httpbin.org/post?a=5&b=6"
)
req2$post(body = list(a = 5))
#> <crul http request> post
#>   url: https://httpbin.org/post?a=5&b=6
#>   curl options: 
#>   proxies: 
#>   auth: 
#>   headers: 
#>   progress: FALSE

(res <- AsyncVaried$new(req1, req2))
#> <crul async varied connection>
#>   requests: (n: 2)
#>    get: https://httpbin.org/get?a=5 
#>    post: https://httpbin.org/post?a=5&b=6
```

Make requests asynchronously


```r
res$request()
```

Parse all results


```r
res$parse()
#> [1] "{\n  \"args\": {\n    \"a\": \"5\"\n  }, \n  \"headers\": {\n    \"Accept\": \"application/json, text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\": \"R (4.0.2 x86_64-apple-darwin17.0 x86_64 darwin17.0)\", \n    \"X-Amzn-Trace-Id\": \"Root=1-5f077ab6-74ff987c9d00dd500c8ad920\"\n  }, \n  \"origin\": \"24.21.229.59\", \n  \"url\": \"https://httpbin.org/get?a=5\"\n}\n"                                                                                                                                                                                                                                                       
#> [2] "{\n  \"args\": {\n    \"a\": \"5\", \n    \"b\": \"6\"\n  }, \n  \"data\": \"\", \n  \"files\": {}, \n  \"form\": {\n    \"a\": \"5\"\n  }, \n  \"headers\": {\n    \"Accept\": \"application/json, text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\", \n    \"Content-Length\": \"137\", \n    \"Content-Type\": \"multipart/form-data; boundary=------------------------cb67be38c1a8a084\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\": \"libcurl/7.64.1 r-curl/4.3 crul/0.9.2.93\", \n    \"X-Amzn-Trace-Id\": \"Root=1-5f077ab6-d206c2e86f429778e6cb9810\"\n  }, \n  \"json\": null, \n  \"origin\": \"24.21.229.59\", \n  \"url\": \"https://httpbin.org/post?a=5&b=6\"\n}\n"
```


```r
lapply(res$parse(), jsonlite::prettify)
#> [[1]]
#> {
#>     "args": {
#>         "a": "5"
#>     },
#>     "headers": {
#>         "Accept": "application/json, text/xml, application/xml, */*",
#>         "Accept-Encoding": "gzip, deflate",
#>         "Host": "httpbin.org",
#>         "User-Agent": "R (4.0.2 x86_64-apple-darwin17.0 x86_64 darwin17.0)",
#>         "X-Amzn-Trace-Id": "Root=1-5f077ab6-74ff987c9d00dd500c8ad920"
#>     },
#>     "origin": "24.21.229.59",
#>     "url": "https://httpbin.org/get?a=5"
#> }
#>  
#> 
#> [[2]]
#> {
#>     "args": {
#>         "a": "5",
#>         "b": "6"
#>     },
#>     "data": "",
#>     "files": {
#> 
#>     },
#>     "form": {
#>         "a": "5"
#>     },
#>     "headers": {
#>         "Accept": "application/json, text/xml, application/xml, */*",
#>         "Accept-Encoding": "gzip, deflate",
#>         "Content-Length": "137",
#>         "Content-Type": "multipart/form-data; boundary=------------------------cb67be38c1a8a084",
#>         "Host": "httpbin.org",
#>         "User-Agent": "libcurl/7.64.1 r-curl/4.3 crul/0.9.2.93",
#>         "X-Amzn-Trace-Id": "Root=1-5f077ab6-d206c2e86f429778e6cb9810"
#>     },
#>     "json": null,
#>     "origin": "24.21.229.59",
#>     "url": "https://httpbin.org/post?a=5&b=6"
#> }
#> 
```

Status codes


```r
res$status_code()
#> [1] 200 200
```
---
title: 4. curl options
author: Scott Chamberlain
date: "2021-02-05"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{4. curl options}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



> adapted in part from the blog post [Curling - exploring web request options][post]


Most times you request data from the web, you should have no problem. However, you eventually will run into problems. In addition, there are advanced things you can do modifying requests to web resources that fall in the _advanced stuff_ category.

Requests to web resources are served over the `http` protocol via [curl][]. `curl` _is a command line tool and library for transferring data with URL syntax, supporting (lots of protocols)_ . `curl` has many options that you may not know about.

I'll go over some of the common and less commonly used curl options, and try to explain why you may want to use some of them.

## Discover curl options

You can go to the source, that is the curl book at <https://everything.curl.dev/>. In R: `curl::curl_options()` for finding curl options. which gives information for each curl option, including the libcurl variable name (e.g., `CURLOPT_CERTINFO`) and the type of variable (e.g., logical).

## Other ways to use curl besides R

Perhaps the canonical way to use curl is on the command line. You can get curl for your operating system at <https://curl.se/download.html>, though hopefully you already have curl. Once you have curl, you can have lots of fun. For example, get the contents of the Google landing page:

```sh
curl https://www.google.com
```

* If you like that you may also like [httpie][], a Python command line tool that is a little more convenient than curl (e.g., JSON output is automatically parsed and colorized).
* Alot of data from the web is in JSON format. A great command line tool to pair with `curl` is [jq][].

> Note: if you are on windows you may require extra setup if you want to play with curl on the command line. OSX and linux have it by default. On Windows 8, installing the latest version from here https://curl.se/download.html#Win64 worked for me.


## general info

With `crul` you have to set curl options per each object, so not globally across all HTTP requests. We may allow the global curl option setting in the future.


## using curl options in other packages

We recommend using `...` to allow users to pass in curl options. For example, lets say you have a function in a package

```r
foo <- function() {
  z <- crul::HttpClient$new(url = yoururl)
  z$get()
}
```

To make it easy for users to pass in curl options use an `...`

```r
foo <- function(...) {
  z <- crul::HttpClient$new(url = yoururl, opts = list(...))
  z$get()
}
```

Then we can pass in any combination of acceptable curl options:

```r
foo(verbose = TRUE)
#> verbose curl output
```

You can instead make users pass in a list, e.g.:

```r
foo <- function(opts = list()) {
  z <- crul::HttpClient$new(url = yoururl, opts = opts)
  z$get()
}
```

Then a user has to pass curl options like: 

```r
foo(opts = list(verbose = TRUE))
```


## timeout

> Set a timeout for a request. If request exceeds timeout, request stops.

relevant commands:

* `timeout_ms=<integer>`

```r
HttpClient$new("https://www.google.com/search", 
  opts = list(timeout_ms = 1))$get()
#> Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) :
#>  Timeout was reached: Operation timed out after 35 milliseconds with 0 bytes received
```

_Why use this?_ You sometimes are working with a web resource that is somewhat unreliable. For example, if you want to run a script on a server that may take many hours, and the web resource could be down at some point during that time, you could set the timeout and error catch the response so that the script doesn't hang on a server that's not responding. Another example could be if you call a web resource in an R package. In your test suite, you may want to test that a web resource is responding quickly, so you could set a timeout, and not test if that fails.

## verbose

> Print detailed info on a curl call

relevant commands:

* `verbose=<boolean>`

Just do a `HEAD` request so we don't have to deal with big output

```r
HttpClient$new("https://httpbin.org", 
  opts = list(verbose = TRUE))$head()
#> > HEAD / HTTP/1.1
#> Host: httpbin.org
#> User-Agent: libcurl/7.54.0 r-curl/3.2 crul/0.5.4.9521
#> Accept: */*
#> Accept-Encoding: gzip, deflate
#> 
#> < HTTP/1.1 200 OK
#> < Connection: keep-alive
#> < Server: gunicorn/19.8.1
#> < Date: Fri, 06 Jul 2018 17:56:50 GMT
#> < Content-Type: text/html; charset=utf-8
#> < Content-Length: 8344
#> < Access-Control-Allow-Origin: *
#> < Access-Control-Allow-Credentials: true
#> < Via: 1.1 vegur
```

_Why use this?_ As you can see verbose output gives you lots of information that may be useful for debugging a request. You typically don't need verbose output unless you want to inspect a request.

## headers

> Add headers to modify requests, including authentication, setting content-type, accept type, etc.

relevant commands:

* `HttpClient$new(headers = list(...))`


```r
x <- HttpClient$new("https://httpbin.org", 
  headers = list(
    Accept = "application/json", 
    foo = "bar"
  ), 
  opts = list(verbose = TRUE)
)
x$head()
#> > HEAD / HTTP/1.1
#> Host: httpbin.org
#> User-Agent: libcurl/7.54.0 r-curl/3.2 crul/0.5.4.9521
#> Accept-Encoding: gzip, deflate
#> Accept: application/json
#> foo: bar
#> 
#> < HTTP/1.1 200 OK
#> < Connection: keep-alive
#> < Server: gunicorn/19.8.1
#> < Date: Fri, 06 Jul 2018 17:59:15 GMT
#> < Content-Type: text/html; charset=utf-8
#> < Content-Length: 8344
#> < Access-Control-Allow-Origin: *
#> < Access-Control-Allow-Credentials: true
#> < Via: 1.1 vegur
```

_Why use this?_ For some web resources, using headers is mandatory, and `httr` makes including them quite easy. Headers are nice too because e.g., passing authentication in the header instead of the URL string means your private data is not as exposed to prying eyes.


## authenticate

> Set authentication details for a resource

relevant commands:

* `auth()`

`auth()` for basic username/password authentication

```r
auth(user = "foo", pwd = "bar")
#> $userpwd
#> [1] "foo:bar"
#> 
#> $httpauth
#> [1] 1
#> 
#> attr(,"class")
#> [1] "auth"
#> attr(,"type")
#> [1] "basic"
```

To use an API key, this depends on the data provider. They may request it one or either of the header


```r
HttpClient$new("https://httpbin.org/get", headers = list(Authorization = "Bearer 234kqhrlj2342"))
```

or as a query parameter (which is passed in the URL string)

```r
HttpClient$new("https://httpbin.org/get", query = list(api_key = "<your key>"))
```

Another authentication option is OAuth. OAuth is not supported in `crul` yet. You can always do OAuth with 
`httr` and then take your token and pass it in as a header/etc. with `crul`.

## cookies

> Set or get cookies.

relevant commands:

* `auth()`

Set cookies (just showing response headers)

```r
x <- HttpClient$new(url = "https://www.google.com", opts = list(verbose = TRUE))
res <- x$get()
#> < HTTP/1.1 200 OK
#> < Date: Fri, 06 Jul 2018 23:25:49 GMT
#> < Expires: -1
#> < Cache-Control: private, max-age=0
#> < Content-Type: text/html; charset=ISO-8859-1
#> < P3P: CP="This is not a P3P policy! See g.co/p3phelp for more info."
#> < Content-Encoding: gzip
#> < Server: gws
#> < X-XSS-Protection: 1; mode=block
#> < X-Frame-Options: SAMEORIGIN
#> * Added cookie 1P_JAR="2018-07-06-23" for domain google.com, path /, expire 1533511549
#> < Set-Cookie: 1P_JAR=2018-07-06-23; expires=Sun, 05-Aug-2018 23:25:49 GMT; path=/; domain=.google.com
#> * Added cookie NID="134=yt47WC-2mhTgQpkSCMz_ySTig3MCJD5Bx_lNj_aVLAwKu8SPMX-CCowKfU8zSv2cJ2OjiX2LTrYnhWMGvIDieCC419v0VHvlm4Hl9xln9-r4MZwcnqwTZQPT72VNE0uA" for domain google.com, path /, expire 1546730749
#> < Set-Cookie: NID=134=yt47WC-2mhTgQpkSCMz_ySTig3MCJD5Bx_lNj_aVLAwKu8SPMX-CCowKfU8zSv2cJ2OjiX2LTrYnhWMGvIDieCC419v0VHvlm4Hl9xln9-r4MZwcnqwTZQPT72VNE0uA; expires=Sat, 05-Jan-2019 23:25:49 GMT; path=/; domain=.google.com; HttpOnly
#> < Alt-Svc: quic=":443"; ma=2592000; v="43,42,41,39,35"
#> < Transfer-Encoding: chunked
```

If there are cookies in a response, you can access them with `curl::handle_cookies` like:

```r
curl::handle_cookies(res$handle)
#>                  domain flag path secure          expiration   name
#> 1           .google.com TRUE    /  FALSE 2018-08-05 16:25:16 1P_JAR
#> 2 #HttpOnly_.google.com TRUE    /  FALSE 2019-01-05 15:25:16    NID
#>   value
#> 1 2018-07-06-23
#> 2 134=4E_Zo-cY8hRLNSj47jRJQML0CPQ8Ip__ ...
```

## progress

> Print curl progress

relevant commands:

* `HttpClient$new(progress = fxn)`

```r
x <- HttpClient$new("https://httpbin.org/get", progress = httr::progress())
#> |==================================| 100%
```

_Why use this?_ As you could imagine, this is increasingly useful as a request for a web resource takes longer and longer. For very long requests, this will help you know approximately when a request will finish.


## proxies

> When behind a proxy, give authentication details for your proxy.

relevant commands:

* `HttpClient$new(proxies = proxy("http://97.77.104.22:3128", "foo", "bar"))`


```r
prox <- proxy("125.39.66.66", port = 80, username = "username", password = "password")
HttpClient$new("http://www.google.com/search", proxies = prox)
```

_Why use this?_ Most of us likely don't need to worry about this. However, if you are in a work place, or maybe in certain geographic locations, you may have to use a proxy. I haven't personally used a proxy in R, so any feedback on this is great.


## user agent

> Some resources require a user-agent string.

relevant commands:

* `HttpClient$new(headers = list(`User-Agent` = "foobar"))`
OR 
* `HttpClient$new(opts = list(useragent = "foobar"))`

both result in the same thing


_Why use this?_ This is set by default in a http request, as you can see in the first example above for user agent. Some web APIs require that you set a specific user agent. For example, the [GitHub API](https://developer.github.com/v3/#user-agent-required) requires that you include a user agent string in the header of each request that is your username or the name of your application so they can contact you if there is a problem.


[curl]: https://curl.se/
[jq]: https://stedolan.github.io/jq/
[httpie]: https://github.com/httpie/httpie
[gbif]: https://www.gbif.org/
[post]: https://ropensci.org/blog/2014/12/18/curl-options/
---
title: 6. Choosing a HTTP request class
author: Scott Chamberlain
date: "2020-07-13"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{6. Choosing a HTTP request class}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



There are a number of different crul classes that do HTTP requests. The following table compares features across the classes.

|Class|HTTP verbs|Asynchronous?|Packges using|
|---|---|---|---|
|HttpClient|all|no|[bcdata][],[chirps][],[duckduckr][],[gfonts][],[mindicador][],[nsapi][],[tradestatistics][],[viafr][]|
|Paginator|all except retry|no|---|
|ok|head,get|no|[dams][]|
|Async|all except retry|yes|[fulltext][],[rdryad][],[rnoaa][]|
|AsyncVaried|all except retry|yes|[ropenaq][],[rcrossref][],[taxize][],[rcitoid][],[mstrio][]|
|AsyncQueue|all except retry|yes|---|

**HttpClient** is the main class for doing synchronous HTTP requests. It supports all HTTP verbs including retry. It was the first class in this package. See also `crul::ok()`, which builds on this class.

**Paginator** is a class for automating pagination. It requires an instance of `HttpClient` as it's first parameter. It does not handle asynchronous requests at this time, but may in the future. `Paginator` may be the right class to use when you don't know the total number of results. Beware however, that if there are A LOT of results (and a lot depends on your internet speed and the server response time) the requests may take a long time to finish - just plan wisely to fit your needs.

**ok** is a convienence function light wrapper around `HttpClient`. It's use case is for determining if a URL is "up", or "okay". You don't have to instantiate an R6 class as you do with the other classes discussed here, but you can pass an `HttpClient` class to it if you like.

With **Async** you can make HTTP requests in parallel. `Async` does not at this time support retry. It's targeted at the use case where you don't mind having request settings the same for all requests - you just pass in any number of URLs and all requests get the same headers, auth, curl options applied, if any.

**AsyncVaried** allows you to customize each request using `HttpRequest` (See below); that is, every HTTP request run asynchronously can have its own curl options, headers, etc. 

**AsyncQueue** is the newest class, inheriting from `AsyncVaried`, introduced in curl v1. The motivation behind this class is: sometimes you may want to do HTTP requests in parallel, but there's rate limiting or some other reason to want to not simply send off all requests immediately. `AsyncQueue` sets up a queue, splitting up requests into buckets, and executes requests based on a `sleep` or `req_per_min` (requests per minute) setting.

## HttpRequest

`HttpRequest` is related here, but not in the table above because it doesn't do actual HTTP requests, but is used to build HTTP requests to pass in to `AsyncVaried`. The simplified class `Async` relative to `AsyncVaried` uses `HttpRequest` internally to build requests. 

## More Async

See the [async with crul](async.html) vignette for more details on asynchronous requests.

[bcdata]: https://bcgov.github.io/bcdata/
[chirps]: https://docs.ropensci.org/chirps/
[duckduckr]: https://github.com/dirkschumacher/duckduckr
[gfonts]: https://dreamrs.github.io/gfonts/
[mindicador]: https://pacha.dev/mindicador/
[nsapi]: https://rmhogervorst.nl/nsapi/
[tradestatistics]: https://docs.ropensci.org/tradestatistics/
[viafr]: https://github.com/stefanieschneider/viafr
[fulltext]:https://docs.ropensci.org/fulltext/
[rdryad]: https://docs.ropensci.org/rdryad/
[rnoaa]: https://docs.ropensci.org/rnoaa/
[ropenaq]: https://docs.ropensci.org/ropenaq/
[rcrossref]: https://docs.ropensci.org/rcrossref/
[taxize]: https://docs.ropensci.org/taxize/
[rcitoid]: https://docs.ropensci.org/rcitoid/
[mstrio]: https://cran.r-project.org/package=mstrio
[dams]: https://jsta.github.io/dams/
---
title: 2. crul workflows
author: Scott Chamberlain
date: "2020-07-09"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{2. crul workflows}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



The following aims to help you decide how to use `crul` in different 
scenarios.

`crul` is aimed a bit more at developers than at the casual 
user doing HTTP requests. That is, `crul` is probably a better fit 
for an R package developer, mainly because it heavily uses `R6` - 
an interface that's very unlike the interface in `httr` but very
similar to interacting with classes in Ruby/Python.

```r
library("crul")
```

## A simple HTTP request function

Most likely you'll want to do a `GET` request - so let's start with that - 
though the details are not much different for other HTTP verbs.

And in most cases you'll likely not want to do asynchronous requests - though
see below if you do.

You'll probably want to write a small function, like so (annotated for 
clarity)



```r
make_request <- function(url) {
  # create a HttpClient object, defining the url
  cli <- crul::HttpClient$new(url = url)
  # do a GET request
  res <- cli$get()
  # check to see if request failed or succeeded
  # - if succeeds this will return nothing and proceeds to next step
  res$raise_for_status()
  # parse response to plain text (JSON in this case) - most likely you'll 
  # want UTF-8 encoding
  txt <- res$parse("UTF-8")
  # parse the JSON to an R list
  jsonlite::fromJSON(txt)
}
```

Use the function



```r
make_request("https://httpbin.org/get")
#> $args
#> named list()
#> 
#> $headers
#> $headers$Accept
#> [1] "application/json, text/xml, application/xml, */*"
#> 
#> $headers$`Accept-Encoding`
#> [1] "gzip, deflate"
#> 
#> $headers$Host
#> [1] "httpbin.org"
#> 
#> $headers$`User-Agent`
#> [1] "libcurl/7.64.1 r-curl/4.3 crul/0.9.2.93"
#> 
#> $headers$`X-Amzn-Trace-Id`
#> [1] "Root=1-5f077ab0-7142dc982d077268af3d2a28"
#> 
#> 
#> $origin
#> [1] "24.21.229.59"
#> 
#> $url
#> [1] "https://httpbin.org/get"
```

Now you can use the `make_request` function in your script or package.

## More customized function

Once you get more familiar (or if you're already familiar with HTTP) you may
want to have more control, toggle more switches.

In the next function, we'll allow for users to pass in curl options, use 
a custom HTTP status checker, and xxx.



```r
make_request2 <- function(url, ...) {
  # create a HttpClient object, defining the url
  cli <- crul::HttpClient$new(url = url)
  # do a GET request, allow curl options to be passed in
  res <- cli$get(...)
  # check to see if request failed or succeeded
  # - a custom approach this time combining status code, 
  #   explanation of the code, and message from the server
  if (res$status_code > 201) {
    mssg <- jsonlite::fromJSON(res$parse("UTF-8"))$message$message
    x <- res$status_http()
    stop(
      sprintf("HTTP (%s) - %s\n  %s", x$status_code, x$explanation, mssg),
      call. = FALSE
    )
  }
  # parse response
  txt <- res$parse("UTF-8")
  # parse the JSON to an R list
  jsonlite::fromJSON(txt)
}
```

Use the function



```r
make_request2("https://api.crossref.org/works?rows=0")
#> $status
#> [1] "ok"
#> 
#> $`message-type`
#> [1] "work-list"
#> 
#> $`message-version`
#> [1] "1.0.0"
#> 
#> $message
#> $message$facets
#> named list()
#> 
#> $message$`total-results`
#> [1] 115235329
#> 
#> $message$items
#> list()
#> 
#> $message$`items-per-page`
#> [1] 0
#> 
#> $message$query
#> $message$query$`start-index`
#> [1] 0
#> 
#> $message$query$`search-terms`
#> NULL
```

No different from the first function (besides the URL). However, now we can 
pass in curl options:


```r
make_request2("https://api.crossref.org/works?rows=0", verbose = TRUE)
make_request2("https://api.crossref.org/works?rows=0", timeout_ms = 1)
```

We can also pass named parameters supported in the `get` method, including
`query`, `disk`, and `stream`.



```r
make_request2("https://api.crossref.org/works", query = list(rows = 0))
#> $status
#> [1] "ok"
#> 
#> $`message-type`
#> [1] "work-list"
#> 
#> $`message-version`
#> [1] "1.0.0"
#> 
#> $message
#> $message$facets
#> named list()
#> 
#> $message$`total-results`
#> [1] 115235329
#> 
#> $message$items
#> list()
#> 
#> $message$`items-per-page`
#> [1] 0
#> 
#> $message$query
#> $message$query$`start-index`
#> [1] 0
#> 
#> $message$query$`search-terms`
#> NULL
```

In addition, the failure behavior is different, and customized to the 
specific web resource we are working with


```r
make_request2("https://api.crossref.org/works?rows=asdf")
#> Error: HTTP (400) - Bad request syntax or unsupported method
#>   Integer specified as asdf but must be a positive integer less than or equal to 1000.
```

## Asynchronous requests

You may want to use asynchronous HTTP requests when any one HTTP request 
takes "too long". This is of course all relative. You may be dealing with a 
server that responds very slowly, or other circumstances. 

See the __async with crul__ vignette for more details on asynchronous requests.
---
title: 5. API package best practices
author: Scott Chamberlain
date: "2021-08-16"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{5. API package best practices}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---


The `crul` package documentation mostly documents how to work
with any particular function or class, but does not detail 
how you would use the package in a more realistic context. This
vignette outlines what we think of as best practices for using
`crul` in scripts or an R package.

## Importing crul

In most cases you'll only need to import one thing from `crul`: 
`HttpClient`. Add crul to `Imports` in your `DESCRIPTION` file,
and add an entry like `@importFrom crul HttpClient` somewhere in 
your package documentation, for example:

```r
#' Some function
#' 
#' @export
#' @importFrom crul HttpClient
#' ...
```

## HTTP request function

If you have more than one function that needs to make an HTTP request
it's probably useful to have a function for doing HTTP requests. The
following is an example of a function.


```r
xGET <- function(url, path, args = list(), ...) {
  cli <- crul::HttpClient$new(url, opts = list(...))
  res <- cli$get(path = path, query = args)
  res$raise_for_status()
  res$raise_for_ct_json()
  res$parse("UTF-8")
}
```

There's some features to note in the above function:

* `url`: this really depends on your setup. In some cases the base URL 
doesn't change, so you can remove the `url` parameter and define the 
url in the `crul::HttpClient$new()` call. 
* `path`: this likely will hold anything after the base path
* `args`: named list of query arguments. the default of `list()` 
means you can then use the function and not have to pass `args` 
in cases where no query args are needed.
* `...`: it's called an ellipsis. see example and discussion below.

You can use the function like:


```r
x <- xGET("https://httpbin.org", "get", args = list(foo = "bar"))
# parse the JSON to a list
jsonlite::fromJSON(x)
# more parsing
```

Because we used an ellipsis, anyone can pass in curl options like:


```r
xGET("https://xxx.org", args = list(foo = "bar"), verbose = TRUE)
```

Curl options are important for digging into the details of HTTP 
requests, and go a long way towards users being able to sort out 
their own problems, and help you diagnose problems as well.

Alternatively, you can just do the HTTP request in your `xGET` function
and return the response object - and line by line, or with 
another function, parse results as needed.

## Failing with fauxpas

[fauxpas][] is in Suggests in this package. If you don't have it
installed, no worries, but if you do have it installed, we use
fauxpas.

There is not much difference with the default `raise_for_status()`
between using fauxpas and not using it. 

However, you can construct your own replacement with fauxpas that 
gives you more flexibility in how you deal with HTTP status codes.

First, make an HTTP request:


```r
con <- HttpClient$new("https://httpbin.org/status/404")
res <- con$get()
```

Then use `fauxpas::find_error_class` to get the correct R6 error 
class for the status code, in this case `404`


```r
x <- fauxpas::find_error_class(res$status_code)$new()
#> <HTTPNotFound>
#>  behavior: stop
#>  message_template: {{reason}} (HTTP {{status}})
#>  message_template_verbose: {{reason}} (HTTP {{status}}).\n - {{message}}
```

We can then do one of two things: use `$do()` or `$do_verbose()`. `$do()` 
is simpler and gives you thhe same thing `$raise_for_status()` gives, but
allows you to change behavior (stop vs. warning vs. message), and how the 
message is formatted. By default we get:


```r
x$do(res)
#> Error: Not Found (HTTP 404)
```

We can change the template using `whisker` templating


```r
x$do(res, template = "{{status}}\n  --> {{reason}}")
#> Error: 404
#>  --> Not Found
```

`$do_verbose()` gives you a lot more detail about the status code, possibly more
than you want:


```r
x$do_verbose(res)
#> Error: Not Found (HTTP 404).
#>  - The server has not found anything matching the Request-URI. No indication
#>  is given of whether the condition is temporary or permanent. The 410 (Gone)
#>  status code SHOULD be used if the server knows, through some internally configurable
#>  mechanism, that an old resource is permanently unavailable and has no forwarding
#>  address. This status code is commonly used when the server does not wish to
#>  reveal exactly why the request has been refused, or when no other response
#> is applicable.
```

You can change behavior to either `warning` or `message`:


```r
x$behavior <- "warning"
x$do(res)
#> Warning message:
#> Not Found (HTTP 404)
x$behavior <- "message"
x$do(res)
#> Not Found (HTTP 404)
```

## Retrying requests

In some failure scenarios it may make sense to retry the same request. 
For example, if a 429 "Too many requests" http status is returned, you
can retry the request after a certain amount of time (that time should
be supplied by the server). We suggest using RETRY if you are in these
scenarios. See [`HttpClient$retry()`](https://docs.ropensci.org/crul/reference/HttpClient.html#method-retry)
for more information.


## Mocking with webmockr

[webmockr][] is a package for stubbing and setting expectations on 
HTTP requests. It has support for working with two HTTP request 
packages: [crul][] and [httr][]. 

There are a variety of use cases for `webmockr`. 

* Use it in an interactive R session where you're working 
on a project and want to mock HTTP requests and set certain responses. 
* You can be on a plane and still allow requests to work without an 
internet connection by setting a response to give back. 
* Test hard to test scenarios in your code or package. `webmockr` 
allows you to give back exact responses just as you describe and 
even fail with certain HTTP conditions. Getting certain failures 
to happen with a remote server can sometimes be difficult. 
* Package test suite: you can use `webmockr` in a test suite, 
although the next section covers `vcr` which builds on top of 
`webmockr` and is ideal for tests.

See the book [HTTP mocking and testing in R][book] for more.


## Testing with vcr

[vcr][] records and replays HTTP requests. Its main use case is for 
caching HTTP requests in test suites in R packages. It has support 
for working with two HTTP request packages: [crul][] and [httr][]. 

To use `vcr` for testing the setup is pretty easy. 

1. Add `vcr` to Suggests in your DESCRIPTION file
2. Make a file in your `tests/testthat/` directory called 
`helper-yourpackage.R` (or skip if as similar file already exists). 
In that file use the following lines to setup your path for storing 
cassettes (change path to whatever you want):

```r
library("vcr")
invisible(vcr::vcr_configure())
```

3. In your tests, for whichever tests you want to use `vcr`, wrap 
the tests in a `vcr::use_cassette()` call like:

```r
library(testthat)
test_that("my test", {
  vcr::use_cassette("rl_citation", {
    aa <- rl_citation()

    expect_is(aa, "character")
    expect_match(aa, "IUCN")
    expect_match(aa, "www.iucnredlist.org")
  })
})
```

That's it! Just run your tests are you normally would and any HTTP 
requests done by `crul` or `httr` will be cached on the first test run
then the cached responses used every time thereafter. 

See the book [HTTP mocking and testing in R][book] for more.

## What else?

Let us know if there's anything else you'd like to see in this 
document and/or if there's anything that can be explained better.

Last, the httr package has a similar article on best practices, see 
<https://httr.r-lib.org/articles/api-packages.html>


[crul]: https://github.com/ropensci/crul
[webmockr]: https://github.com/ropensci/webmockr
[vcr]: https://github.com/ropensci/vcr
[httr]: https://github.com/r-lib/httr
[fauxpas]: https://github.com/ropensci/fauxpas
[book]: https://books.ropensci.org/http-testing/
---
title: 1. crul introduction
author: Scott Chamberlain
date: "2020-07-10"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{1. crul introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



`crul` is an HTTP client for R.

## Install

Stable CRAN version


```r
install.packages("crul")
```

Dev version


```r
remotes::install_github("ropensci/crul")
```


```r
library("crul")
```

## HttpClient - the main interface

`HttpClient` is where to start


```r
(x <- HttpClient$new(
  url = "https://httpbin.org",
  opts = list(
    timeout = 1
  ),
  headers = list(
    a = "hello world"
  )
))
#> <crul connection> 
#>   url: https://httpbin.org
#>   curl options: 
#>     timeout: 1
#>   proxies: 
#>   auth: 
#>   headers: 
#>     a: hello world
#>   progress: FALSE
#>   hooks:
```

Makes a R6 class, that has all the bits and bobs you'd expect for doing HTTP
requests. When it prints, it gives any defaults you've set. As you update
the object you can see what's been set


```r
x$opts
#> $timeout
#> [1] 1
```


```r
x$headers
#> $a
#> [1] "hello world"
```

## Do some HTTP requests

The client object created above has http methods that you can call,
and pass paths to, as well as query parameters, body values, and any other
curl options.

Here, we'll do a __GET__ request on the route `/get` on our base url
`https://httpbin.org` (the full url is then `https://httpbin.org/get`)


```r
res <- x$get("get")
```

The response from a http request is another R6 class `HttpResponse`, which
has slots for the outputs of the request, and some functions to deal with
the response:

Status code


```r
res$status_code
#> [1] 200
```

The content


```r
res$content
#>   [1] 7b 0a 20 20 22 61 72 67 73 22 3a 20 7b 7d 2c 20 0a 20 20 22 68 65 61 64 65
#>  [26] 72 73 22 3a 20 7b 0a 20 20 20 20 22 41 22 3a 20 22 68 65 6c 6c 6f 20 77 6f
#>  [51] 72 6c 64 22 2c 20 0a 20 20 20 20 22 41 63 63 65 70 74 22 3a 20 22 61 70 70
#>  [76] 6c 69 63 61 74 69 6f 6e 2f 6a 73 6f 6e 2c 20 74 65 78 74 2f 78 6d 6c 2c 20
#> [101] 61 70 70 6c 69 63 61 74 69 6f 6e 2f 78 6d 6c 2c 20 2a 2f 2a 22 2c 20 0a 20
#> [126] 20 20 20 22 41 63 63 65 70 74 2d 45 6e 63 6f 64 69 6e 67 22 3a 20 22 67 7a
#> [151] 69 70 2c 20 64 65 66 6c 61 74 65 22 2c 20 0a 20 20 20 20 22 48 6f 73 74 22
#> [176] 3a 20 22 68 74 74 70 62 69 6e 2e 6f 72 67 22 2c 20 0a 20 20 20 20 22 55 73
#> [201] 65 72 2d 41 67 65 6e 74 22 3a 20 22 6c 69 62 63 75 72 6c 2f 37 2e 36 34 2e
#> [226] 31 20 72 2d 63 75 72 6c 2f 34 2e 33 20 63 72 75 6c 2f 30 2e 39 2e 34 2e 39
#> [251] 31 22 2c 20 0a 20 20 20 20 22 58 2d 41 6d 7a 6e 2d 54 72 61 63 65 2d 49 64
#> [276] 22 3a 20 22 52 6f 6f 74 3d 31 2d 35 66 30 38 64 36 63 65 2d 61 61 32 30 39
#> [301] 37 30 64 63 62 63 31 33 64 30 61 37 65 38 66 32 35 65 36 22 0a 20 20 7d 2c
#> [326] 20 0a 20 20 22 6f 72 69 67 69 6e 22 3a 20 22 32 34 2e 32 31 2e 32 32 39 2e
#> [351] 35 39 22 2c 20 0a 20 20 22 75 72 6c 22 3a 20 22 68 74 74 70 73 3a 2f 2f 68
#> [376] 74 74 70 62 69 6e 2e 6f 72 67 2f 67 65 74 22 0a 7d 0a
```

HTTP method


```r
res$method
#> [1] "get"
```

Request headers


```r
res$request_headers
#> $`User-Agent`
#> [1] "libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91"
#> 
#> $`Accept-Encoding`
#> [1] "gzip, deflate"
#> 
#> $Accept
#> [1] "application/json, text/xml, application/xml, */*"
#> 
#> $a
#> [1] "hello world"
```

Response headers


```r
res$response_headers
#> $status
#> [1] "HTTP/2 200 "
#> 
#> $date
#> [1] "Fri, 10 Jul 2020 20:59:58 GMT"
#> 
#> $`content-type`
#> [1] "application/json"
#> 
#> $`content-length`
#> [1] "393"
#> 
#> $server
#> [1] "gunicorn/19.9.0"
#> 
#> $`access-control-allow-origin`
#> [1] "*"
#> 
#> $`access-control-allow-credentials`
#> [1] "true"
```

All response headers, including intermediate headers, if any


```r
res$response_headers_all
#> [[1]]
#> [[1]]$status
#> [1] "HTTP/2 200 "
#> 
#> [[1]]$date
#> [1] "Fri, 10 Jul 2020 20:59:58 GMT"
#> 
#> [[1]]$`content-type`
#> [1] "application/json"
#> 
#> [[1]]$`content-length`
#> [1] "393"
#> 
#> [[1]]$server
#> [1] "gunicorn/19.9.0"
#> 
#> [[1]]$`access-control-allow-origin`
#> [1] "*"
#> 
#> [[1]]$`access-control-allow-credentials`
#> [1] "true"
```

And you can parse the content with a provided function:


```r
res$parse()
#> [1] "{\n  \"args\": {}, \n  \"headers\": {\n    \"A\": \"hello world\", \n    \"Accept\": \"application/json, text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\": \"libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91\", \n    \"X-Amzn-Trace-Id\": \"Root=1-5f08d6ce-aa20970dcbc13d0a7e8f25e6\"\n  }, \n  \"origin\": \"24.21.229.59\", \n  \"url\": \"https://httpbin.org/get\"\n}\n"
jsonlite::fromJSON(res$parse())
#> $args
#> named list()
#> 
#> $headers
#> $headers$A
#> [1] "hello world"
#> 
#> $headers$Accept
#> [1] "application/json, text/xml, application/xml, */*"
#> 
#> $headers$`Accept-Encoding`
#> [1] "gzip, deflate"
#> 
#> $headers$Host
#> [1] "httpbin.org"
#> 
#> $headers$`User-Agent`
#> [1] "libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91"
#> 
#> $headers$`X-Amzn-Trace-Id`
#> [1] "Root=1-5f08d6ce-aa20970dcbc13d0a7e8f25e6"
#> 
#> 
#> $origin
#> [1] "24.21.229.59"
#> 
#> $url
#> [1] "https://httpbin.org/get"
```

With the `HttpClient` object, which holds any configuration stuff
we set, we can make other HTTP verb requests. For example, a `HEAD`
request:


```r
x$post(
  path = "post", 
  body = list(hello = "world")
)
```


## write to disk


```r
x <- HttpClient$new(url = "https://httpbin.org")
f <- tempfile()
res <- x$get(disk = f)
# when using write to disk, content is a path
res$content 
#> [1] "/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T//RtmppSPqIf/file20b5651a4a3a"
```

Read lines


```r
readLines(res$content, n = 10)
#>  [1] "<!DOCTYPE html>"                                                                                                               
#>  [2] "<html lang=\"en\">"                                                                                                            
#>  [3] ""                                                                                                                              
#>  [4] "<head>"                                                                                                                        
#>  [5] "    <meta charset=\"UTF-8\">"                                                                                                  
#>  [6] "    <title>httpbin.org</title>"                                                                                                
#>  [7] "    <link href=\"https://fonts.googleapis.com/css?family=Open+Sans:400,700|Source+Code+Pro:300,600|Titillium+Web:400,600,700\""
#>  [8] "        rel=\"stylesheet\">"                                                                                                   
#>  [9] "    <link rel=\"stylesheet\" type=\"text/css\" href=\"/flasgger_static/swagger-ui.css\">"                                      
#> [10] "    <link rel=\"icon\" type=\"image/png\" href=\"/static/favicon.ico\" sizes=\"64x64 32x32 16x16\" />"
```

## stream data


```r
(x <- HttpClient$new(url = "https://httpbin.org"))
#> <crul connection> 
#>   url: https://httpbin.org
#>   curl options: 
#>   proxies: 
#>   auth: 
#>   headers: 
#>   progress: FALSE
#>   hooks:
res <- x$get('stream/5', stream = function(x) cat(rawToChar(x)))
#> {"url": "https://httpbin.org/stream/5", "args": {}, "headers": {"Host": "httpbin.org", "X-Amzn-Trace-Id": "Root=1-5f08d6cf-871364ce53d065c8d32eefd4", "User-Agent": "libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91", "Accept-Encoding": "gzip, deflate", "Accept": "application/json, text/xml, application/xml, */*"}, "origin": "24.21.229.59", "id": 0}
#> {"url": "https://httpbin.org/stream/5", "args": {}, "headers": {"Host": "httpbin.org", "X-Amzn-Trace-Id": "Root=1-5f08d6cf-871364ce53d065c8d32eefd4", "User-Agent": "libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91", "Accept-Encoding": "gzip, deflate", "Accept": "application/json, text/xml, application/xml, */*"}, "origin": "24.21.229.59", "id": 1}
#> {"url": "https://httpbin.org/stream/5", "args": {}, "headers": {"Host": "httpbin.org", "X-Amzn-Trace-Id": "Root=1-5f08d6cf-871364ce53d065c8d32eefd4", "User-Agent": "libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91", "Accept-Encoding": "gzip, deflate", "Accept": "application/json, text/xml, application/xml, */*"}, "origin": "24.21.229.59", "id": 2}
#> {"url": "https://httpbin.org/stream/5", "args": {}, "headers": {"Host": "httpbin.org", "X-Amzn-Trace-Id": "Root=1-5f08d6cf-871364ce53d065c8d32eefd4", "User-Agent": "libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91", "Accept-Encoding": "gzip, deflate", "Accept": "application/json, text/xml, application/xml, */*"}, "origin": "24.21.229.59", "id": 3}
#> {"url": "https://httpbin.org/stream/5", "args": {}, "headers": {"Host": "httpbin.org", "X-Amzn-Trace-Id": "Root=1-5f08d6cf-871364ce53d065c8d32eefd4", "User-Agent": "libcurl/7.64.1 r-curl/4.3 crul/0.9.4.91", "Accept-Encoding": "gzip, deflate", "Accept": "application/json, text/xml, application/xml, */*"}, "origin": "24.21.229.59", "id": 4}
# when streaming, content is NULL
res$content 
#> NULL
```

## Learn more 

Learn more with the other vignettes:

- [crul workflows](how-to-use-crul.html)
- [async with crul](async.html)
- [curl options](curl-options.html)
- [API package best practices](best-practices-api-packages.html)
- [Choosing a HTTP request class](choosing-a-client.html)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/response.R
\name{HttpResponse}
\alias{HttpResponse}
\title{Base HTTP response object}
\description{
Class with methods for handling HTTP responses
}
\details{
\strong{Additional Methods}
\describe{
\item{\code{raise_for_ct(type, charset = NULL, behavior = "stop")}}{
Check response content-type; stop or warn if not matched. Parameters:
\itemize{
\item type: (character) a mime type to match against; see
\link[mime:mimemap]{mime::mimemap} for allowed values
\item charset: (character) if a charset string given, we check that
it matches the charset in the content type header. default: NULL
\item behavior: (character) one of stop (default) or warning
}
}
\item{\code{raise_for_ct_html(charset = NULL, behavior = "stop")}}{
Check that the response content-type is \code{text/html}; stop or warn if
not matched. Parameters: see \code{raise_for_ct()}
}
\item{\code{raise_for_ct_json(charset = NULL, behavior = "stop")}}{
Check that the response content-type is \code{application/json}; stop or
warn if not matched. Parameters: see \code{raise_for_ct()}
}
\item{\code{raise_for_ct_xml(charset = NULL, behavior = "stop")}}{
Check that the response content-type is \code{application/xml}; stop or warn if
not matched. Parameters: see \code{raise_for_ct()}
}
}
}
\section{R6 classes}{

This is an R6 class from the package \pkg{R6}. Find out more
about R6 at \url{https://r6.r-lib.org/}. After creating an instance of an R6
class (e.g., \code{x <- HttpClient$new(url = "https://httpbin.org")}) you can
access values and methods on the object \code{x}.
}

\examples{
\dontrun{
x <- HttpResponse$new(method = "get", url = "https://httpbin.org")
x$url
x$method

x <- HttpClient$new(url = 'https://httpbin.org')
(res <- x$get('get'))
res$request_headers
res$response_headers
res$parse()
res$status_code
res$status_http()
res$status_http()$status_code
res$status_http()$message
res$status_http()$explanation
res$success()

x <- HttpClient$new(url = 'https://httpbin.org/status/404')
(res <- x$get())
# res$raise_for_status()

x <- HttpClient$new(url = 'https://httpbin.org/status/414')
(res <- x$get())
# res$raise_for_status()
}
}
\seealso{
\link{content-types}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{method}}{(character) one or more URLs}

\item{\code{url}}{(character) one or more URLs}

\item{\code{opts}}{(character) one or more URLs}

\item{\code{handle}}{(character) one or more URLs}

\item{\code{status_code}}{(character) one or more URLs}

\item{\code{request_headers}}{(character) one or more URLs}

\item{\code{response_headers}}{(character) one or more URLs}

\item{\code{response_headers_all}}{(character) one or more URLs}

\item{\code{modified}}{(character) one or more URLs}

\item{\code{times}}{(character) one or more URLs}

\item{\code{content}}{(character) one or more URLs}

\item{\code{request}}{(character) one or more URLs}

\item{\code{raise_for_ct}}{for ct method (general)}

\item{\code{raise_for_ct_html}}{for ct method (html)}

\item{\code{raise_for_ct_json}}{for ct method (json)}

\item{\code{raise_for_ct_xml}}{for ct method (xml)}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{HttpResponse$print()}}
\item \href{#method-new}{\code{HttpResponse$new()}}
\item \href{#method-parse}{\code{HttpResponse$parse()}}
\item \href{#method-success}{\code{HttpResponse$success()}}
\item \href{#method-status_http}{\code{HttpResponse$status_http()}}
\item \href{#method-raise_for_status}{\code{HttpResponse$raise_for_status()}}
\item \href{#method-clone}{\code{HttpResponse$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for HttpResponse objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpResponse$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new HttpResponse object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpResponse$new(
  method,
  url,
  opts,
  handle,
  status_code,
  request_headers,
  response_headers,
  response_headers_all,
  modified,
  times,
  content,
  request
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{(character) HTTP method}

\item{\code{url}}{(character) A url, required}

\item{\code{opts}}{(list) curl options}

\item{\code{handle}}{A handle}

\item{\code{status_code}}{(integer) status code}

\item{\code{request_headers}}{(list) request headers, named list}

\item{\code{response_headers}}{(list) response headers, named list}

\item{\code{response_headers_all}}{(list) all response headers, including
intermediate redirect headers, unnamed list of named lists}

\item{\code{modified}}{(character) modified date}

\item{\code{times}}{(vector) named vector}

\item{\code{content}}{(raw) raw binary content response}

\item{\code{request}}{request object, with all details}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse"></a>}}
\if{latex}{\out{\hypertarget{method-parse}{}}}
\subsection{Method \code{parse()}}{
Parse the raw response content to text
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpResponse$parse(encoding = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{encoding}}{(character) A character string describing the
current encoding. If left as \code{NULL}, we attempt to guess the
encoding. Passed to \code{from} parameter in \code{iconv}}

\item{\code{...}}{additional parameters passed on to \code{iconv} (options: sub,
mark, toRaw). See \code{?iconv} for help}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-success"></a>}}
\if{latex}{\out{\hypertarget{method-success}{}}}
\subsection{Method \code{success()}}{
Was status code less than or equal to 201
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpResponse$success()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status_http"></a>}}
\if{latex}{\out{\hypertarget{method-status_http}{}}}
\subsection{Method \code{status_http()}}{
Get HTTP status code, message, and explanation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpResponse$status_http(verbose = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{verbose}}{(logical) whether to get verbose http status description,
default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
object of class "http_code", a list with slots for status_code,
message, and explanation
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-raise_for_status"></a>}}
\if{latex}{\out{\hypertarget{method-raise_for_status}{}}}
\subsection{Method \code{raise_for_status()}}{
Check HTTP status and stop with appropriate
HTTP error code and message if >= 300. otherwise use \pkg{httpcode}.
If you have \code{fauxpas} installed we use that.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpResponse$raise_for_status()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
stop or warn with message
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpResponse$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/progress.R
\name{progress}
\alias{progress}
\title{progress bars}
\description{
progress bars
}
\details{
pass \code{httr::progress()} to \code{progress} param in \link{HttpClient},
which pulls out relevant info to pass down to \pkg{curl}

if file sizes known you get progress bar; if file sizes not
known you get bytes downloaded

See the README for examples
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curl_options.R
\name{curl_verbose}
\alias{curl_verbose}
\title{curl verbose method}
\usage{
curl_verbose(data_out = TRUE, data_in = FALSE, info = FALSE, ssl = FALSE)
}
\arguments{
\item{data_out}{Show data sent to the server}

\item{data_in}{Show data recieved from the server}

\item{info}{Show informational text from curl. This is mainly useful for
debugging https and auth problems, so is disabled by default}

\item{ssl}{Show even data sent/recieved over SSL connections?}
}
\description{
curl verbose method
}
\details{
line prefixes:
\itemize{
\item \code{*} informative curl messages
\item \verb{=>} headers sent (out)
\item \code{>} data sent (out)
\item \verb{*>} ssl data sent (out)
\item \code{<=} headers received (in)
\item \code{<} data received (in)
\item \verb{<*} ssl data received (in)
}
}
\note{
adapted from \code{httr::verbose}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/content-types.R
\name{content-types}
\alias{content-types}
\title{Working with content types}
\description{
The \link{HttpResponse} class holds all the responses elements for an HTTP
request. This document details how to work specifically with the
content-type of the response headers
}
\section{Content types}{

The "Content-Type" header in HTTP responses gives the media type of the
response. The media type is both the data format and how the data is
intended to be processed by a recipient. (modified from rfc7231)
}

\section{Behavior of the parameters HttpResponse raise_for_ct* methods}{

\itemize{
\item type: (only applicable for the \code{raise_for_ct()} method): instead of
using one of the three other content type methods for html, json, or xml,
you can specify a mime type to check, any of those in \link[mime:mimemap]{mime::mimemap}
\item charset: if you don't give a value to this parameter, we only
check that the content type is what you expect; that is, the charset,
if given, is ignored.
\item behavior: by default when you call this method, and the content type
does not match what the method expects, then we run \code{stop()} with a
message. Instead of stopping, you can choose \code{behavior="warning"}
and we'll throw a warning instead, allowing any downstream processing
to proceed.
}
}

\examples{
\dontrun{
(x <- HttpClient$new(url = "https://httpbin.org"))
(res <- x$get())

## see the content type
res$response_headers

## check that the content type is text/html
res$raise_for_ct_html()

## it's def. not json
# res$raise_for_ct_json()

## give custom content type
res$raise_for_ct("text/html")
# res$raise_for_ct("application/json")
# res$raise_for_ct("foo/bar")

## check charset in addition to the media type
res$raise_for_ct_html(charset = "utf-8")
# res$raise_for_ct_html(charset = "utf-16")

# warn instead of stop
res$raise_for_ct_json(behavior = "warning")
}
}
\references{
spec for content types:
\url{https://tools.ietf.org/html/rfc7231#section-3.1.1.5}

spec for media types:
\url{https://tools.ietf.org/html/rfc7231#section-3.1.1.1}
}
\seealso{
\link{HttpResponse}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_url.R
\name{url_build}
\alias{url_build}
\alias{url_parse}
\title{Build and parse URLs}
\usage{
url_build(url, path = NULL, query = NULL)

url_parse(url)
}
\arguments{
\item{url}{(character) a url, length 1}

\item{path}{(character) a path, length 1}

\item{query}{(list) a named list of query parameters}
}
\value{
\code{url_build} returns a character string URL; \code{url_parse}
returns a list with URL components
}
\description{
Build and parse URLs
}
\examples{
url_build("https://httpbin.org")
url_build("https://httpbin.org", "get")
url_build("https://httpbin.org", "post")
url_build("https://httpbin.org", "get", list(foo = "bar"))

url_parse("httpbin.org")
url_parse("http://httpbin.org")
url_parse(url = "https://httpbin.org")
url_parse("https://httpbin.org/get")
url_parse("https://httpbin.org/get?foo=bar")
url_parse("https://httpbin.org/get?foo=bar&stuff=things")
url_parse("https://httpbin.org/get?foo=bar&stuff=things[]")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crul-package.r
\docType{package}
\name{crul-package}
\alias{crul-package}
\alias{crul}
\title{crul}
\description{
\strong{HTTP R client}
}
\section{Package API}{

\itemize{
\item \code{\link[=HttpClient]{HttpClient()}} - create a connection client, set all
your http options, make http requests
\item \code{\link[=HttpResponse]{HttpResponse()}} - mostly for internal use, handles
http responses
\item \code{\link[=Paginator]{Paginator()}} - auto-paginate through requests
\item \code{\link[=Async]{Async()}} - asynchronous requests
\item \code{\link[=AsyncVaried]{AsyncVaried()}} - varied asynchronous requests
\item \code{\link[=HttpRequest]{HttpRequest()}} - generate an HTTP request, mostly for
use in building requests to be used in \code{Async} or \code{AsyncVaried}
\item \code{\link[=mock]{mock()}} - Turn on/off mocking, via \code{webmockr}
\item \code{\link[=auth]{auth()}} - Simple authentication helper
\item \code{\link[=proxy]{proxy()}} - Proxy helper
\item \code{\link[=upload]{upload()}} - File upload helper
\item set curl options globally: \code{\link[=set_auth]{set_auth()}}, \code{\link[=set_headers]{set_headers()}},
\code{\link[=set_opts]{set_opts()}}, \code{\link[=set_proxy]{set_proxy()}}, and \code{\link[=crul_settings]{crul_settings()}}
}
}

\section{HTTP verbs (or HTTP request methods)}{


See \link{verb-GET}, \link{verb-POST}, \link{verb-PUT}, \link{verb-PATCH}, \link{verb-DELETE},
\link{verb-HEAD} for details.
\itemize{
\item \link{HttpClient} is the main interface for making HTTP requests,
and includes methods for each HTTP verb
\item \link{HttpRequest} allows you to prepare a HTTP payload for use with
\link{AsyncVaried}, which provides asynchronous requests for varied
HTTP methods
\item \link{Async} provides asynchronous requests for a single HTTP method
at a time
\item the \code{verb()} method can be used on all the above to request
a specific HTTP verb
}
}

\section{Checking HTTP responses}{


\code{\link[=HttpResponse]{HttpResponse()}} has helpers for checking and raising warnings/errors.
\itemize{
\item \link{content-types} details the various options for checking content
types and throwing a warning or error if the response content
type doesn't match what you expect. Mis-matched content-types are
typically a good sign of a bad response. There's methods built
in for json, xml and html, with the ability to set any
custom content type
\item \code{raise_for_status()} is a method on \code{\link[=HttpResponse]{HttpResponse()}} that checks
the HTTP status code, and errors with the appropriate message for
the HTTP status code, optionally using the package \code{fauxpas}
if it's installed.
}
}

\section{HTTP conditions}{

We use \code{fauxpas} if you have it installed for handling HTTP
conditions but if it's not installed we use \pkg{httpcode}
}

\section{Mocking}{

Mocking HTTP requests is supported via the \pkg{webmockr}
package. See \link{mock} for guidance, and
\url{https://books.ropensci.org/http-testing/}
}

\section{Caching}{

Caching HTTP requests is supported via the \pkg{vcr}
package. See \url{https://books.ropensci.org/http-testing/}
}

\section{Links}{


Source code: \url{https://github.com/ropensci/crul}

Bug reports/feature requests: \url{https://github.com/ropensci/crul/issues}
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verbs.R
\name{verb-DELETE}
\alias{verb-DELETE}
\title{HTTP verb info: DELETE}
\description{
The DELETE method deletes the specified resource.
}
\section{The DELETE method}{

The DELETE method requests that the origin server remove the
association between the target resource and its current
functionality.  In effect, this method is similar to the rm command
in UNIX: it expresses a deletion operation on the URI mapping of the
origin server rather than an expectation that the previously
associated information be deleted.

See \url{https://tools.ietf.org/html/rfc7231#section-4.3.5} for further
details.
}

\examples{
\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")
x$delete(path = 'delete')

## a list
(res1 <- x$delete('delete', body = list(hello = "world"), verbose = TRUE))
jsonlite::fromJSON(res1$parse("UTF-8"))

## a string
(res2 <- x$delete('delete', body = "hello world", verbose = TRUE))
jsonlite::fromJSON(res2$parse("UTF-8"))

## empty body request
x$delete('delete', verbose = TRUE)
}

}
\references{
\url{https://tools.ietf.org/html/rfc7231#section-4.3.5}
}
\seealso{
\link{crul-package}

Other verbs: 
\code{\link{verb-GET}},
\code{\link{verb-HEAD}},
\code{\link{verb-PATCH}},
\code{\link{verb-POST}},
\code{\link{verb-PUT}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set.R
\name{crul-options}
\alias{crul-options}
\alias{set_opts}
\alias{set_verbose}
\alias{set_proxy}
\alias{set_auth}
\alias{set_headers}
\alias{crul_settings}
\title{Set curl options, proxy, and basic auth}
\usage{
set_opts(...)

set_verbose()

set_proxy(x)

set_auth(x)

set_headers(...)

crul_settings(reset = FALSE)
}
\arguments{
\item{...}{For \code{set_opts()} any curl option in the set
\code{\link[curl:curl_options]{curl::curl_options()}}. For \code{set_headers()} a named list of headers}

\item{x}{For \code{set_proxy()} a \code{proxy} object made with \code{\link[=proxy]{proxy()}}. For
\code{set_auth()} an \code{auth} object made with \code{\link[=auth]{auth()}}}

\item{reset}{(logical) reset all settings (aka, delete them).
Default: \code{FALSE}}
}
\description{
Set curl options, proxy, and basic auth
}
\details{
\itemize{
\item \code{set_opts()}: set curl options; supports any options in
\code{\link[curl:curl_options]{curl::curl_options()}}
\item \code{set_verbose()}: set custom curl verbose; sets \code{verbose=TRUE}
and \code{debugfunction} to the callback result from \code{\link[=curl_verbose]{curl_verbose()}}
\item \code{set_proxy()}: set proxy settings, accepts \code{\link[=proxy]{proxy()}}
\item \code{set_auth()}: set authorization, accepts \code{\link[=auth]{auth()}}
\item \code{set_headers()}: set request headers, a named list
\item \code{crul_settings()}: list all settigns set via these functions
}
}
\note{
the \code{mock} option will be seen in output of \code{crul_settings()}
but is set via the function \code{\link[=mock]{mock()}}
}
\examples{
if (interactive()) {
# get settings
crul_settings()

# curl options
set_opts(timeout_ms = 1000)
crul_settings()
set_opts(timeout_ms = 4000)
crul_settings()
set_opts(verbose = TRUE)
crul_settings()
\dontrun{
HttpClient$new('https://httpbin.org')$get('get')
}
# set_verbose - sets: `verbose=TRUE`, and `debugfunction` to 
# result of call to `curl_verbose()`, see `?curl_verbose`
set_verbose()
crul_settings()

# basic authentication
set_auth(auth(user = "foo", pwd = "bar", auth = "basic"))
crul_settings()

# proxies
set_proxy(proxy("http://97.77.104.22:3128"))
crul_settings()

# headers
crul_settings(TRUE) # reset first
set_headers(foo = "bar")
crul_settings()
set_headers(`User-Agent` = "hello world")
crul_settings()
\dontrun{
set_opts(verbose = TRUE)
HttpClient$new('https://httpbin.org')$get('get')
}

# reset
crul_settings(TRUE)
crul_settings()

# works with async functions
## Async
set_opts(verbose = TRUE)
cc <- Async$new(urls = c(
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?foo=bar'))
(res <- cc$get())

## AsyncVaried
set_opts(verbose = TRUE)
set_headers(stuff = "things")
reqlist <- list(
  HttpRequest$new(url = "https://httpbin.org/get")$get(),
  HttpRequest$new(url = "https://httpbin.org/post")$post())
out <- AsyncVaried$new(.list = reqlist)
out$request()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mocking.R
\name{mock}
\alias{mock}
\title{Mocking HTTP requests}
\usage{
mock(on = TRUE)
}
\arguments{
\item{on}{(logical) turn mocking on with \code{TRUE} or turn off with \code{FALSE}.
By default is \code{FALSE}}
}
\description{
Mocking HTTP requests
}
\details{
\code{webmockr} package required for mocking behavior
}
\examples{
\dontrun{

if (interactive()) {
  # load webmockr
  library(webmockr)
  library(crul)

  URL <- "https://httpbin.org"

  # turn on mocking
  crul::mock()

  # stub a request
  stub_request("get", file.path(URL, "get"))
  webmockr:::webmockr_stub_registry

  # create an HTTP client
  (x <- HttpClient$new(url = URL))

  # make a request - matches stub - no real request made
  x$get('get')

  # allow net connect
  webmockr::webmockr_allow_net_connect()
  x$get('get', query = list(foo = "bar"))
  webmockr::webmockr_disable_net_connect()
  x$get('get', query = list(foo = "bar"))
}

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/async-queue.R
\name{AsyncQueue}
\alias{AsyncQueue}
\title{AsyncQueue}
\description{
An AsyncQueue client
}
\section{R6 classes}{

This is an R6 class from the package \pkg{R6}. Find out more
about R6 at \url{https://r6.r-lib.org/}. After creating an instance of an R6
class (e.g., \code{x <- HttpClient$new(url = "https://httpbin.org")}) you can
access values and methods on the object \code{x}.
}

\examples{
\dontrun{
# Using sleep
reqlist <- list(
  HttpRequest$new(url = "https://httpbin.org/get")$get(),
  HttpRequest$new(url = "https://httpbin.org/post")$post(),
  HttpRequest$new(url = "https://httpbin.org/put")$put(),
  HttpRequest$new(url = "https://httpbin.org/delete")$delete(),
  HttpRequest$new(url = "https://httpbin.org/get?g=5")$get(),
  HttpRequest$new(
    url = "https://httpbin.org/post")$post(body = list(y = 9)),
  HttpRequest$new(
    url = "https://httpbin.org/get")$get(query = list(hello = "world")),
  HttpRequest$new(url = "https://ropensci.org")$get(),
  HttpRequest$new(url = "https://ropensci.org/about")$get(),
  HttpRequest$new(url = "https://ropensci.org/packages")$get(),
  HttpRequest$new(url = "https://ropensci.org/community")$get(),
  HttpRequest$new(url = "https://ropensci.org/blog")$get(),
  HttpRequest$new(url = "https://ropensci.org/careers")$get()
)
out <- AsyncQueue$new(.list = reqlist, bucket_size = 5, sleep = 3)
out
out$bucket_size # bucket size
out$requests() # list requests
out$request() # make requests
out$responses() # list responses

# Using requests per minute
if (interactive()) {
x="https://raw.githubusercontent.com/ropensci/roregistry/gh-pages/registry.json"
z <- HttpClient$new(x)$get()
urls <- jsonlite::fromJSON(z$parse("UTF-8"))$packages$url
repos = Filter(length, regmatches(urls, gregexpr("ropensci/[A-Za-z]+", urls)))
repos = unlist(repos)
auth <- list(Authorization = paste("token", Sys.getenv('GITHUB_PAT')))
reqs <- lapply(repos[1:50], function(w) {
  HttpRequest$new(paste0("https://api.github.com/repos/", w), headers = auth)$get()
})

out <- AsyncQueue$new(.list = reqs, req_per_min = 30)
out
out$bucket_size
out$requests()
out$request()
out$responses()
}}
}
\seealso{
Other async: 
\code{\link{AsyncVaried}},
\code{\link{Async}},
\code{\link{HttpRequest}}
}
\concept{async}
\section{Super class}{
\code{\link[crul:AsyncVaried]{crul::AsyncVaried}} -> \code{AsyncQueue}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{bucket_size}}{(integer) number of requests to send at once}

\item{\code{sleep}}{(integer) number of seconds to sleep between each bucket}

\item{\code{req_per_min}}{(integer) requests per minute}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{AsyncQueue$print()}}
\item \href{#method-new}{\code{AsyncQueue$new()}}
\item \href{#method-request}{\code{AsyncQueue$request()}}
\item \href{#method-responses}{\code{AsyncQueue$responses()}}
\item \href{#method-parse}{\code{AsyncQueue$parse()}}
\item \href{#method-status_code}{\code{AsyncQueue$status_code()}}
\item \href{#method-status}{\code{AsyncQueue$status()}}
\item \href{#method-content}{\code{AsyncQueue$content()}}
\item \href{#method-times}{\code{AsyncQueue$times()}}
\item \href{#method-clone}{\code{AsyncQueue$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="crul" data-topic="AsyncVaried" data-id="requests">}\href{../../crul/html/AsyncVaried.html#method-requests}{\code{crul::AsyncVaried$requests()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for AsyncQueue objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{AsyncQueue} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$new(
  ...,
  .list = list(),
  bucket_size = 5,
  sleep = NULL,
  req_per_min = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{..., .list}}{Any number of objects of class \code{\link[=HttpRequest]{HttpRequest()}},
must supply inputs to one of these parameters, but not both}

\item{\code{bucket_size}}{(integer) number of requests to send at once.
default: 5. See Details.}

\item{\code{sleep}}{(integer) seconds to sleep between buckets.
default: NULL (not set)}

\item{\code{req_per_min}}{(integer) maximum number of requests per minute.
if \code{NULL} (default), its ignored}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Must set either \code{sleep} or \code{req_per_min}. If you set
\code{req_per_min} we calculate a new \code{bucket_size} when \verb{$new()} is
called
}

\subsection{Returns}{
A new \code{AsyncQueue} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-request"></a>}}
\if{latex}{\out{\hypertarget{method-request}{}}}
\subsection{Method \code{request()}}{
Execute asynchronous requests
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$request()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing, responses stored inside object, though will print
messages if you choose verbose output
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-responses"></a>}}
\if{latex}{\out{\hypertarget{method-responses}{}}}
\subsection{Method \code{responses()}}{
List responses
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$responses()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list of \code{HttpResponse} objects, empty list before
requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse"></a>}}
\if{latex}{\out{\hypertarget{method-parse}{}}}
\subsection{Method \code{parse()}}{
parse content
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$parse(encoding = "UTF-8")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{encoding}}{(character) the encoding to use in parsing.
default:"UTF-8"}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character vector, empty character vector before
requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status_code"></a>}}
\if{latex}{\out{\hypertarget{method-status_code}{}}}
\subsection{Method \code{status_code()}}{
Get HTTP status codes for each response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$status_code()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
numeric vector, empty numeric vector before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status"></a>}}
\if{latex}{\out{\hypertarget{method-status}{}}}
\subsection{Method \code{status()}}{
List HTTP status objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$status()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list of \code{http_code} objects, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-content"></a>}}
\if{latex}{\out{\hypertarget{method-content}{}}}
\subsection{Method \code{content()}}{
Get raw content for each response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$content()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
raw list, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-times"></a>}}
\if{latex}{\out{\hypertarget{method-times}{}}}
\subsection{Method \code{times()}}{
curl request times
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$times()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list of named numeric vectors, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncQueue$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verbs.R
\name{verb-PATCH}
\alias{verb-PATCH}
\title{HTTP verb info: PATCH}
\description{
The PATCH method is used to apply partial modifications to a resource.
}
\section{The PATCH method}{

The PATCH method requests that a set of changes described in the
request entity be applied to the resource identified by the Request-
URI.  The set of changes is represented in a format called a "patch
document" identified by a media type.  If the Request-URI does not
point to an existing resource, the server MAY create a new resource,
depending on the patch document type (whether it can logically modify
a null resource) and permissions, etc.

See \url{https://tools.ietf.org/html/rfc5789#section-2} for further
details.
}

\examples{
\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")
x$patch(path = 'patch', body = list(hello = "mars"))
}

}
\references{
\url{https://tools.ietf.org/html/rfc5789}
}
\seealso{
\link{crul-package}

Other verbs: 
\code{\link{verb-DELETE}},
\code{\link{verb-GET}},
\code{\link{verb-HEAD}},
\code{\link{verb-POST}},
\code{\link{verb-PUT}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http-headers.R
\name{http-headers}
\alias{http-headers}
\title{Working with HTTP headers}
\description{
Working with HTTP headers
}
\examples{
\dontrun{
(x <- HttpClient$new(url = "https://httpbin.org"))

# set headers
(res <- HttpClient$new(
  url = "https://httpbin.org",
  opts = list(
    verbose = TRUE
  ),
  headers = list(
    a = "stuff",
    b = "things"
  )
))
res$headers
# reassign header value
res$headers$a <- "that"
# define new header
res$headers$c <- "what"
# request
res$get('get')

## setting content-type via headers
(res <- HttpClient$new(
  url = "https://httpbin.org",
  opts = list(
    verbose = TRUE
  ),
  headers = list(`Content-Type` = "application/json")
))
res$get('get')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verbs.R
\name{verb-POST}
\alias{verb-POST}
\title{HTTP verb info: POST}
\description{
The POST method is used to submit an entity to the specified resource,
often causing a change in state or side effects on the server.
}
\section{The POST method}{

If one or more resources has been created on the origin server as a
result of successfully processing a POST request, the origin server
SHOULD send a 201 (Created) response containing a Location header
field that provides an identifier for the primary resource created
(Section 7.1.2 \url{https://tools.ietf.org/html/rfc7231#section-7.1.2})
and a representation that describes the status of the
request while referring to the new resource(s).

See \url{https://tools.ietf.org/html/rfc7231#section-4.3.3} for further
details.
}

\examples{
\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")

# a named list
x$post(path='post', body = list(hello = "world"))

# a string
x$post(path='post', body = "hello world")

# an empty body request
x$post(path='post')

# encode="form"
res <- x$post(path="post",
  encode = "form",
  body = list(
    custname = 'Jane',
    custtel = '444-4444',
    size = 'small',
    topping = 'bacon',
    comments = 'make it snappy'
  )
)
jsonlite::fromJSON(res$parse("UTF-8"))

# encode="json"
res <- x$post("post",
  encode = "json",
  body = list(
    genus = 'Gagea',
    species = 'pratensis'
  )
)
jsonlite::fromJSON(res$parse())
}

}
\references{
\url{https://tools.ietf.org/html/rfc7231#section-4.3.3}
}
\seealso{
\link{crul-package}

Other verbs: 
\code{\link{verb-DELETE}},
\code{\link{verb-GET}},
\code{\link{verb-HEAD}},
\code{\link{verb-PATCH}},
\code{\link{verb-PUT}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paginator.R
\name{Paginator}
\alias{Paginator}
\title{Paginator client}
\value{
a list, with objects of class \code{\link[=HttpResponse]{HttpResponse()}}.
Responses are returned in the order they are passed in.
}
\description{
A client to help you paginate
}
\details{
See \code{\link[=HttpClient]{HttpClient()}} for information on parameters
}
\section{R6 classes}{

This is an R6 class from the package \pkg{R6}. Find out more
about R6 at \url{https://r6.r-lib.org/}. After creating an instance of an R6
class (e.g., \code{x <- HttpClient$new(url = "https://httpbin.org")}) you can
access values and methods on the object \code{x}.
}

\section{Methods to paginate}{


Supported now:
\itemize{
\item \code{limit_offset}: the most common way (in my experience), so is the default.
This method involves setting how many records and what record to start at
for each request. We send these query parameters for you.
\item \code{page_perpage}: set the page to fetch and (optionally) how many records
to get per page
}

Supported later, hopefully:
\itemize{
\item \code{link_headers}: link headers are URLS for the next/previous/last
request given in the response header from the server. This is relatively
uncommon, though is recommended by JSONAPI and is implemented by a
well known API (GitHub).
\item \code{cursor}: this works by a single string given back in each response, to
be passed in the subsequent response, and so on until no more records
remain. This is common in Solr
}
}

\examples{
\dontrun{
if (interactive()) {
# limit/offset approach
con <- HttpClient$new(url = "https://api.crossref.org")
cc <- Paginator$new(client = con, limit_param = "rows",
   offset_param = "offset", limit = 50, chunk = 10)
cc
cc$get('works')
cc
cc$responses()
cc$status()
cc$status_code()
cc$times()
# cc$content()
cc$parse()
lapply(cc$parse(), jsonlite::fromJSON)

# page/per page approach (with no per_page param allowed)
conn <- HttpClient$new(url = "https://discuss.ropensci.org")
cc <- Paginator$new(client = conn, by = "page_perpage",
 page_param = "page", per_page_param = "per_page", limit = 90, chunk = 30)
cc
cc$get('c/usecases/l/latest.json')
cc$responses()
lapply(cc$parse(), jsonlite::fromJSON)

# page/per_page
conn <- HttpClient$new('https://api.inaturalist.org')
cc <- Paginator$new(conn, by = "page_perpage", page_param = "page",
 per_page_param = "per_page", limit = 90, chunk = 30)
cc
cc$get('v1/observations', query = list(taxon_name="Helianthus"))
cc$responses()
res <- lapply(cc$parse(), jsonlite::fromJSON)
res[[1]]$total_results
vapply(res, "[[", 1L, "page")
vapply(res, "[[", 1L, "per_page")
vapply(res, function(w) NROW(w$results), 1L)
## another
ccc <- Paginator$new(conn, by = "page_perpage", page_param = "page",
 per_page_param = "per_page", limit = 500, chunk = 30, progress = TRUE)
ccc
ccc$get('v1/observations', query = list(taxon_name="Helianthus"))
res2 <- lapply(ccc$parse(), jsonlite::fromJSON)
vapply(res2, function(w) NROW(w$results), 1L)

# progress bar
(con <- HttpClient$new(url = "https://api.crossref.org"))
cc <- Paginator$new(client = con, limit_param = "rows",
   offset_param = "offset", limit = 50, chunk = 10,
   progress = TRUE)
cc
cc$get('works')
}}

## ------------------------------------------------
## Method `Paginator$url_fetch`
## ------------------------------------------------

\dontrun{
cli <- HttpClient$new(url = "https://api.crossref.org")
cc <- Paginator$new(client = cli, limit_param = "rows",
   offset_param = "offset", limit = 50, chunk = 10)
cc$url_fetch('works')
cc$url_fetch('works', query = list(query = "NSF"))
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{http_req}}{an object of class \code{HttpClient}}

\item{\code{by}}{(character) how to paginate. Only 'limit_offset' supported
for now. In the future will support 'link_headers' and 'cursor'.
See Details.}

\item{\code{chunk}}{(numeric/integer) the number by which to chunk
requests, e.g., 10 would be be each request gets 10 records.
number is passed through \code{\link[=format]{format()}} to prevent larger numbers
from being scientifically formatted}

\item{\code{limit_param}}{(character) the name of the limit parameter.
Default: limit}

\item{\code{offset_param}}{(character) the name of the offset parameter.
Default: offset}

\item{\code{limit}}{(numeric/integer) the maximum records wanted.
number is passed through \code{\link[=format]{format()}} to prevent larger numbers
from being scientifically formatted}

\item{\code{page_param}}{(character) the name of the page parameter.
Default: NULL}

\item{\code{per_page_param}}{(character) the name of the per page parameter.
Default: NULL}

\item{\code{progress}}{(logical) print a progress bar, using \link[utils:txtProgressBar]{utils::txtProgressBar}.
Default: \code{FALSE}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{Paginator$print()}}
\item \href{#method-new}{\code{Paginator$new()}}
\item \href{#method-get}{\code{Paginator$get()}}
\item \href{#method-post}{\code{Paginator$post()}}
\item \href{#method-put}{\code{Paginator$put()}}
\item \href{#method-patch}{\code{Paginator$patch()}}
\item \href{#method-delete}{\code{Paginator$delete()}}
\item \href{#method-head}{\code{Paginator$head()}}
\item \href{#method-responses}{\code{Paginator$responses()}}
\item \href{#method-status_code}{\code{Paginator$status_code()}}
\item \href{#method-status}{\code{Paginator$status()}}
\item \href{#method-parse}{\code{Paginator$parse()}}
\item \href{#method-content}{\code{Paginator$content()}}
\item \href{#method-times}{\code{Paginator$times()}}
\item \href{#method-url_fetch}{\code{Paginator$url_fetch()}}
\item \href{#method-clone}{\code{Paginator$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for \code{Paginator} objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Paginator} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$new(
  client,
  by = "limit_offset",
  limit_param = NULL,
  offset_param = NULL,
  limit = NULL,
  chunk = NULL,
  page_param = NULL,
  per_page_param = NULL,
  progress = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{client}}{an object of class \code{HttpClient}, from a call to \link{HttpClient}}

\item{\code{by}}{(character) how to paginate. Only 'limit_offset' supported for
now. In the future will support 'link_headers' and 'cursor'. See Details.}

\item{\code{limit_param}}{(character) the name of the limit parameter.
Default: limit}

\item{\code{offset_param}}{(character) the name of the offset parameter.
Default: offset}

\item{\code{limit}}{(numeric/integer) the maximum records wanted}

\item{\code{chunk}}{(numeric/integer) the number by which to chunk requests,
e.g., 10 would be be each request gets 10 records}

\item{\code{page_param}}{(character) the name of the page parameter.}

\item{\code{per_page_param}}{(character) the name of the per page parameter.}

\item{\code{progress}}{(logical) print a progress bar, using \link[utils:txtProgressBar]{utils::txtProgressBar}.
Default: \code{FALSE}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Paginator} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get"></a>}}
\if{latex}{\out{\hypertarget{method-get}{}}}
\subsection{Method \code{get()}}{
make a paginated GET request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$get(path = NULL, query = list(), ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-post"></a>}}
\if{latex}{\out{\hypertarget{method-post}{}}}
\subsection{Method \code{post()}}{
make a paginated POST request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$post(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-put"></a>}}
\if{latex}{\out{\hypertarget{method-put}{}}}
\subsection{Method \code{put()}}{
make a paginated PUT request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$put(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-patch"></a>}}
\if{latex}{\out{\hypertarget{method-patch}{}}}
\subsection{Method \code{patch()}}{
make a paginated PATCH request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$patch(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-delete"></a>}}
\if{latex}{\out{\hypertarget{method-delete}{}}}
\subsection{Method \code{delete()}}{
make a paginated DELETE request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$delete(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-head"></a>}}
\if{latex}{\out{\hypertarget{method-head}{}}}
\subsection{Method \code{head()}}{
make a paginated HEAD request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$head(path = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
not sure if this makes any sense or not yet
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-responses"></a>}}
\if{latex}{\out{\hypertarget{method-responses}{}}}
\subsection{Method \code{responses()}}{
list responses
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$responses()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list of \code{HttpResponse} objects, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status_code"></a>}}
\if{latex}{\out{\hypertarget{method-status_code}{}}}
\subsection{Method \code{status_code()}}{
Get HTTP status codes for each response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$status_code()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
numeric vector, empty numeric vector before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status"></a>}}
\if{latex}{\out{\hypertarget{method-status}{}}}
\subsection{Method \code{status()}}{
List HTTP status objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$status()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list of \code{http_code} objects, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse"></a>}}
\if{latex}{\out{\hypertarget{method-parse}{}}}
\subsection{Method \code{parse()}}{
parse content
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$parse(encoding = "UTF-8")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{encoding}}{(character) the encoding to use in parsing.
default:"UTF-8"}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character vector, empty character vector before
requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-content"></a>}}
\if{latex}{\out{\hypertarget{method-content}{}}}
\subsection{Method \code{content()}}{
Get raw content for each response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$content()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
raw list, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-times"></a>}}
\if{latex}{\out{\hypertarget{method-times}{}}}
\subsection{Method \code{times()}}{
curl request times
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$times()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list of named numeric vectors, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-url_fetch"></a>}}
\if{latex}{\out{\hypertarget{method-url_fetch}{}}}
\subsection{Method \code{url_fetch()}}{
get the URL that would be sent (i.e., before executing
the request) the only things that change the URL are path and query
parameters; body and any curl options don't change the URL
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$url_fetch(path = NULL, query = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
URLs (character)
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
cli <- HttpClient$new(url = "https://api.crossref.org")
cc <- Paginator$new(client = cli, limit_param = "rows",
   offset_param = "offset", limit = 50, chunk = 10)
cc$url_fetch('works')
cc$url_fetch('works', query = list(query = "NSF"))
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Paginator$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle.R
\name{handle}
\alias{handle}
\title{Make a handle}
\usage{
handle(url, ...)
}
\arguments{
\item{url}{(character) A url. required.}

\item{...}{options passed on to \code{\link[curl:handle]{curl::new_handle()}}}
}
\description{
Make a handle
}
\examples{
handle("https://httpbin.org")

# handles - pass in your own handle
\dontrun{
h <- handle("https://httpbin.org")
(res <- HttpClient$new(handle = h))
out <- res$get("get")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curl-options.R
\name{curl-options}
\alias{curl-options}
\alias{user-agent}
\alias{verbose}
\alias{timeout}
\title{curl options}
\description{
With the \code{opts} parameter you can pass in various
curl options, including user agent string, whether to get verbose
curl output or not, setting a timeout for requests, and more. See
\code{\link[curl:curl_options]{curl::curl_options()}} for all the options you can use. Note that
you need to give curl options exactly as given in
\code{\link[curl:curl_options]{curl::curl_options()}}.
}
\examples{
\dontrun{
url <- "https://httpbin.org"

# set curl options on client initialization
(res <- HttpClient$new(url = url, opts = list(verbose = TRUE)))
res$opts
res$get('get')

# or set curl options when performing HTTP operation
(res <- HttpClient$new(url = url))
res$get('get', verbose = TRUE)
res$get('get', stuff = "things")

# set a timeout
(res <- HttpClient$new(url = url, opts = list(timeout_ms = 1)))
# res$get('get')

# set user agent either as a header or an option
HttpClient$new(url = url,
  headers = list(`User-Agent` = "hello world"),
  opts = list(verbose = TRUE)
)$get('get')

HttpClient$new(url = url,
  opts = list(verbose = TRUE, useragent = "hello world")
)$get('get')

# You can also set custom debug function via the verbose 
# parameter when calling `$new()`
res <- HttpClient$new(url, verbose=curl_verbose())
res
res$get("get")
res <- HttpClient$new(url, verbose=curl_verbose(data_in=TRUE))
res$get("get")
res <- HttpClient$new(url, verbose=curl_verbose(info=TRUE))
res$get("get")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/upload.R
\name{upload}
\alias{upload}
\title{upload file}
\usage{
upload(path, type = NULL)
}
\arguments{
\item{path}{(character) a single path, file must exist}

\item{type}{(character) a file type, guessed by \link[mime:guess_type]{mime::guess_type} if
not given}
}
\description{
upload file
}
\examples{
\dontrun{
# image
path <- file.path(Sys.getenv("R_DOC_DIR"), "html/logo.jpg")
(x <- HttpClient$new(url = "https://eu.httpbin.org"))
res <- x$post(path = "post", body = list(y = upload(path)))
res$content

# text file, in a list
file <- upload(system.file("CITATION"))
res <- x$post(path = "post", body = list(y = file))
jsonlite::fromJSON(res$parse("UTF-8"))

# text file, as data
res <- x$post(path = "post", body = file)
jsonlite::fromJSON(res$parse("UTF-8"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hooks.R
\name{hooks}
\alias{hooks}
\title{Event Hooks}
\description{
Trigger functions to run on requests and/or responses.
See Details for more.
}
\details{
Functions passed to \code{request} are run \strong{before} the
request occurs. The meaning of triggering a function on the
request is that you can do things to the request object.

Functions passed to \code{response} are run \strong{once} the
request is done, and the response object is created.
The meaning of triggering a function on the
response is to do things on the response object.

The above for request and response applies the same
whether you make real HTTP requests or mock with
\code{webmockr}.
}
\note{
Only supported on \link{HttpClient} for now
}
\examples{
\dontrun{
# hooks on the request
fun_req <- function(request) {
  cat(paste0("Requesting: ", request$url$url), sep = "\n")
}
(x <- HttpClient$new(url = "https://httpbin.org",
  hooks = list(request = fun_req)))
x$hooks
x$hooks$request
r1 <- x$get('get')

captured_req <- list()
fun_req2 <- function(request) {
  cat("Capturing Request", sep = "\n")
  captured_req <<- request
}
(x <- HttpClient$new(url = "https://httpbin.org",
  hooks = list(request = fun_req2)))
x$hooks
x$hooks$request
r1 <- x$get('get')
captured_req



# hooks on the response
fun_resp <- function(response) {
  cat(paste0("status_code: ", response$status_code), sep = "\n")
}
(x <- HttpClient$new(url = "https://httpbin.org",
  hooks = list(response = fun_resp)))
x$url
x$hooks
r1 <- x$get('get')

# both
(x <- HttpClient$new(url = "https://httpbin.org",
  hooks = list(request = fun_req, response = fun_resp)))
x$get("get")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writing-options.R
\name{writing-options}
\alias{writing-options}
\title{Writing data options}
\description{
Writing data options
}
\examples{
\dontrun{
# write to disk
(x <- HttpClient$new(url = "https://httpbin.org"))
f <- tempfile()
res <- x$get("get", disk = f)
res$content # when using write to disk, content is a path
readLines(res$content)
close(file(f))

# streaming response
(x <- HttpClient$new(url = "https://httpbin.org"))
res <- x$get('stream/50', stream = function(x) cat(rawToChar(x)))
res$content # when streaming, content is NULL


## Async
(cc <- Async$new(
  urls = c(
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?foo=bar',
    'https://httpbin.org/get?b=4',
    'https://httpbin.org/get?stuff=things',
    'https://httpbin.org/get?b=4&g=7&u=9&z=1'
  )
))
files <- replicate(5, tempfile())
(res <- cc$get(disk = files, verbose = TRUE))
lapply(files, readLines)

## Async varied
### disk
f <- tempfile()
g <- tempfile()
req1 <- HttpRequest$new(url = "https://httpbin.org/get")$get(disk = f)
req2 <- HttpRequest$new(url = "https://httpbin.org/post")$post(disk = g)
req3 <- HttpRequest$new(url = "https://httpbin.org/get")$get()
(out <- AsyncVaried$new(req1, req2, req3))
out$request()
out$content()
readLines(f)
readLines(g)
out$parse()
close(file(f))
close(file(g))

### stream - to console
fun <- function(x) print(x)
req1 <- HttpRequest$new(url = "https://httpbin.org/get"
)$get(query = list(foo = "bar"), stream = fun)
req2 <- HttpRequest$new(url = "https://httpbin.org/get"
)$get(query = list(hello = "world"), stream = fun)
(out <- AsyncVaried$new(req1, req2))
out$request()
out$content()

### stream - to an R object
lst <- list()
fun <- function(x) lst <<- append(lst, list(x))
req1 <- HttpRequest$new(url = "https://httpbin.org/get"
)$get(query = list(foo = "bar"), stream = fun)
req2 <- HttpRequest$new(url = "https://httpbin.org/get"
)$get(query = list(hello = "world"), stream = fun)
(out <- AsyncVaried$new(req1, req2))
out$request()
lst
cat(vapply(lst, function(z) rawToChar(z$content), ""), sep = "\n")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client.R
\name{HttpClient}
\alias{HttpClient}
\title{HTTP client}
\value{
an \link{HttpResponse} object
}
\description{
Create and execute HTTP requests
}
\note{
A little quirk about \code{crul} is that because user agent string can
be passed as either a header or a curl option (both lead to a \code{User-Agent}
header being passed in the HTTP request), we return the user agent
string in the \code{request_headers} list of the response even if you
pass in a \code{useragent} string as a curl option. Note that whether you pass
in as a header like \code{User-Agent} or as a curl option like \code{useragent},
it is returned as \code{request_headers$User-Agent} so at least accessing
it in the request headers is consistent.
}
\section{R6 classes}{

This is an R6 class from the package \pkg{R6}. Find out more
about R6 at \url{https://r6.r-lib.org/}. After creating an instance of an R6
class (e.g., \code{x <- HttpClient$new(url = "https://httpbin.org")}) you can
access values and methods on the object \code{x}.
}

\section{handles}{

curl handles are re-used on the level of the connection object, that is,
each \code{HttpClient} object is separate from one another so as to better
separate connections.

If you don't pass in a curl handle to the \code{handle} parameter,
it gets created when a HTTP verb is called. Thus, if you try to get \code{handle}
after creating a \code{HttpClient} object only passing \code{url} parameter, \code{handle}
will be \code{NULL}. If you pass a curl handle to the \code{handle} parameter, then
you can get the handle from the \code{HttpClient} object. The response from a
http verb request does have the handle in the \code{handle} slot.
}

\examples{
\dontrun{
# set your own handle
(h <- handle("https://httpbin.org"))
(x <- HttpClient$new(handle = h))
x$handle
x$url
(out <- x$get("get"))
x$handle
x$url
class(out)
out$handle
out$request_headers
out$response_headers
out$response_headers_all

# if you just pass a url, we create a handle for you
#  this is how most people will use HttpClient
(x <- HttpClient$new(url = "https://httpbin.org"))
x$url
x$handle # is empty, it gets created when a HTTP verb is called
(r1 <- x$get('get'))
x$url
x$handle
r1$url
r1$handle
r1$content
r1$response_headers
r1$parse()

(res_get2 <- x$get('get', query = list(hello = "world")))
res_get2$parse()
library("jsonlite")
jsonlite::fromJSON(res_get2$parse())

# post request
(res_post <- x$post('post', body = list(hello = "world")))

## empty body request
x$post('post')

# put request
(res_put <- x$put('put'))

# delete request
(res_delete <- x$delete('delete'))

# patch request
(res_patch <- x$patch('patch'))

# head request
(res_head <- x$head())

# query params are URL encoded for you, so DO NOT do it yourself
## if you url encode yourself, it gets double encoded, and that's bad
(x <- HttpClient$new(url = "https://httpbin.org"))
res <- x$get("get", query = list(a = 'hello world'))

# access intermediate headers in response_headers_all
x <- HttpClient$new("https://doi.org/10.1007/978-3-642-40455-9_52-1")
bb <- x$get()
bb$response_headers_all
}

## ------------------------------------------------
## Method `HttpClient$verb`
## ------------------------------------------------

\dontrun{
(x <- HttpClient$new(url = "https://httpbin.org"))
x$verb('get')
x$verb('GET')
x$verb('GET', query = list(foo = "bar"))
x$verb('retry', 'GET', path = "status/400")
}

## ------------------------------------------------
## Method `HttpClient$retry`
## ------------------------------------------------

\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")

# retry, by default at most 3 times
(res_get <- x$retry("GET", path = "status/400"))

# retry, but not for 404 NOT FOUND
(res_get <- x$retry("GET", path = "status/404", terminate_on = c(404)))

# retry, but only for exceeding rate limit (note that e.g. Github uses 403)
(res_get <- x$retry("GET", path = "status/429", retry_only_on = c(403, 429)))
}

## ------------------------------------------------
## Method `HttpClient$url_fetch`
## ------------------------------------------------

x <- HttpClient$new(url = "https://httpbin.org")
x$url_fetch()
x$url_fetch('get')
x$url_fetch('post')
x$url_fetch('get', query = list(foo = "bar"))
}
\seealso{
\link{http-headers}, \link{writing-options}, \link{cookies}, \link{hooks}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{url}}{(character) a url}

\item{\code{opts}}{(list) named list of curl options}

\item{\code{proxies}}{a \code{\link[=proxy]{proxy()}} object}

\item{\code{auth}}{an \code{\link[=auth]{auth()}} object}

\item{\code{headers}}{(list) named list of headers, see \link{http-headers}}

\item{\code{handle}}{a \code{\link[=handle]{handle()}}}

\item{\code{progress}}{only supports \code{httr::progress()}, see \link{progress}}

\item{\code{hooks}}{a named list, see \link{hooks}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{HttpClient$print()}}
\item \href{#method-new}{\code{HttpClient$new()}}
\item \href{#method-get}{\code{HttpClient$get()}}
\item \href{#method-post}{\code{HttpClient$post()}}
\item \href{#method-put}{\code{HttpClient$put()}}
\item \href{#method-patch}{\code{HttpClient$patch()}}
\item \href{#method-delete}{\code{HttpClient$delete()}}
\item \href{#method-head}{\code{HttpClient$head()}}
\item \href{#method-verb}{\code{HttpClient$verb()}}
\item \href{#method-retry}{\code{HttpClient$retry()}}
\item \href{#method-handle_pop}{\code{HttpClient$handle_pop()}}
\item \href{#method-url_fetch}{\code{HttpClient$url_fetch()}}
\item \href{#method-clone}{\code{HttpClient$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for \code{HttpClient} objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new HttpClient object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$new(
  url,
  opts,
  proxies,
  auth,
  headers,
  handle,
  progress,
  hooks,
  verbose
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{url}}{(character) A url. One of \code{url} or \code{handle} required.}

\item{\code{opts}}{any curl options}

\item{\code{proxies}}{a \code{\link[=proxy]{proxy()}} object}

\item{\code{auth}}{an \code{\link[=auth]{auth()}} object}

\item{\code{headers}}{named list of headers, see \link{http-headers}}

\item{\code{handle}}{a \code{\link[=handle]{handle()}}}

\item{\code{progress}}{only supports \code{httr::progress()}, see \link{progress}}

\item{\code{hooks}}{a named list, see \link{hooks}}

\item{\code{verbose}}{a special handler for verbose curl output,
accepts a function only. default is \code{NULL}. if used, \code{verbose}
and \code{debugfunction} curl options are ignored if passed to \code{opts}
on \verb{$new()} and ignored if \code{...} passed to a http method call}

\item{\code{urls}}{(character) one or more URLs}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{HttpClient} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get"></a>}}
\if{latex}{\out{\hypertarget{method-get}{}}}
\subsection{Method \code{get()}}{
Make a GET request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$get(path = NULL, query = list(), disk = NULL, stream = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-post"></a>}}
\if{latex}{\out{\hypertarget{method-post}{}}}
\subsection{Method \code{post()}}{
Make a POST request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$post(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-put"></a>}}
\if{latex}{\out{\hypertarget{method-put}{}}}
\subsection{Method \code{put()}}{
Make a PUT request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$put(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-patch"></a>}}
\if{latex}{\out{\hypertarget{method-patch}{}}}
\subsection{Method \code{patch()}}{
Make a PATCH request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$patch(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-delete"></a>}}
\if{latex}{\out{\hypertarget{method-delete}{}}}
\subsection{Method \code{delete()}}{
Make a DELETE request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$delete(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-head"></a>}}
\if{latex}{\out{\hypertarget{method-head}{}}}
\subsection{Method \code{head()}}{
Make a HEAD request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$head(path = NULL, query = list(), ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-verb"></a>}}
\if{latex}{\out{\hypertarget{method-verb}{}}}
\subsection{Method \code{verb()}}{
Use an arbitrary HTTP verb supported on this class
Supported verbs: "get", "post", "put", "patch", "delete", "head". Also
supports retry
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$verb(verb, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{verb}}{an HTTP verb supported on this class: "get",
"post", "put", "patch", "delete", "head". Also supports retry.}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
(x <- HttpClient$new(url = "https://httpbin.org"))
x$verb('get')
x$verb('GET')
x$verb('GET', query = list(foo = "bar"))
x$verb('retry', 'GET', path = "status/400")
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-retry"></a>}}
\if{latex}{\out{\hypertarget{method-retry}{}}}
\subsection{Method \code{retry()}}{
Retry a request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$retry(
  verb,
  ...,
  pause_base = 1,
  pause_cap = 60,
  pause_min = 1,
  times = 3,
  terminate_on = NULL,
  retry_only_on = NULL,
  onwait = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{verb}}{an HTTP verb supported on this class: "get",
"post", "put", "patch", "delete", "head". Also supports retry.}

\item{\code{...}}{For \code{retry}, the options to be passed on to the method
implementing the requested verb, including curl options. Otherwise,
curl options, only those in the acceptable set from \code{\link[curl:curl_options]{curl::curl_options()}}
except the following: httpget, httppost, post, postfields, postfieldsize,
and customrequest}

\item{\code{pause_base, pause_cap, pause_min}}{basis, maximum, and minimum for
calculating wait time for retry. Wait time is calculated according to the
exponential backoff with full jitter algorithm. Specifically, wait time is
chosen randomly between \code{pause_min} and the lesser of \code{pause_base * 2} and
\code{pause_cap}, with \code{pause_base} doubling on each subsequent retry attempt.
Use \code{pause_cap = Inf} to not terminate retrying due to cap of wait time
reached.}

\item{\code{times}}{the maximum number of times to retry. Set to \code{Inf} to
not stop retrying due to exhausting the number of attempts.}

\item{\code{terminate_on, retry_only_on}}{a vector of HTTP status codes. For
\code{terminate_on}, the status codes for which to terminate retrying, and for
\code{retry_only_on}, the status codes for which to retry the request.}

\item{\code{onwait}}{a callback function if the request will be retried and
a wait time is being applied. The function will be passed two parameters,
the response object from the failed request, and the wait time in seconds.
Note that the time spent in the function effectively adds to the wait time,
so it should be kept simple.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
Retries the request given by \code{verb} until successful
(HTTP response status < 400), or a condition for giving up is met.
Automatically recognizes \code{Retry-After} and \code{X-RateLimit-Reset} headers
in the response for rate-limited remote APIs.
}

\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")

# retry, by default at most 3 times
(res_get <- x$retry("GET", path = "status/400"))

# retry, but not for 404 NOT FOUND
(res_get <- x$retry("GET", path = "status/404", terminate_on = c(404)))

# retry, but only for exceeding rate limit (note that e.g. Github uses 403)
(res_get <- x$retry("GET", path = "status/429", retry_only_on = c(403, 429)))
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-handle_pop"></a>}}
\if{latex}{\out{\hypertarget{method-handle_pop}{}}}
\subsection{Method \code{handle_pop()}}{
reset your curl handle
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$handle_pop()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-url_fetch"></a>}}
\if{latex}{\out{\hypertarget{method-url_fetch}{}}}
\subsection{Method \code{url_fetch()}}{
get the URL that would be sent (i.e., before executing
the request) the only things that change the URL are path and query
parameters; body and any curl options don't change the URL
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$url_fetch(path = NULL, query = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list. any numeric values are
passed through \code{\link[=format]{format()}} to prevent larger numbers from being
scientifically formatted}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
URL (character)
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{x <- HttpClient$new(url = "https://httpbin.org")
x$url_fetch()
x$url_fetch('get')
x$url_fetch('post')
x$url_fetch('get', query = list(foo = "bar"))
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpClient$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verbs.R
\name{verb-GET}
\alias{verb-GET}
\title{HTTP verb info: GET}
\description{
The GET method requests a representation of the specified resource.
Requests using GET should only retrieve data.
}
\section{The GET method}{

The GET method requests transfer of a current selected
representation for the target resource.  GET is the primary
mechanism of information retrieval and the focus of almost all
performance optimizations. Hence, when people speak of retrieving
some identifiable information via HTTP, they are generally referring
to making a GET request.

It is tempting to think of resource identifiers as remote file system
pathnames and of representations as being a copy of the contents of
such files.  In fact, that is how many resources are implemented (see
Section 9.1 (\url{https://tools.ietf.org/html/rfc7231#section-9.1})
for related security considerations).  However, there are
no such limitations in practice.  The HTTP interface for a resource
is just as likely to be implemented as a tree of content objects, a
programmatic view on various database records, or a gateway to other
information systems.  Even when the URI mapping mechanism is tied to
a file system, an origin server might be configured to execute the
files with the request as input and send the output as the
representation rather than transfer the files directly.  Regardless,
only the origin server needs to know how each of its resource identifiers
corresponds to an implementation and how each implementation manages
to select and send a current representation of the target resource
in a response to GET.

A client can alter the semantics of GET to be a "range request",
requesting transfer of only some part(s) of the selected
representation, by sending a Range header field in the request
(RFC7233: \url{https://tools.ietf.org/html/rfc7233}).

A payload within a GET request message has no defined semantics;
sending a payload body on a GET request might cause some existing
implementations to reject the request.

The response to a GET request is cacheable; a cache MAY use it to
satisfy subsequent GET and HEAD requests unless otherwise indicated
by the Cache-Control header field (Section 5.2 of RFC7234:
\url{https://tools.ietf.org/html/rfc7234#section-5.2}).
}

\examples{
\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")
x$get(path = 'get')
}

}
\references{
\url{https://tools.ietf.org/html/rfc7231#section-4.3.1}
}
\seealso{
\link{crul-package}

Other verbs: 
\code{\link{verb-DELETE}},
\code{\link{verb-HEAD}},
\code{\link{verb-PATCH}},
\code{\link{verb-POST}},
\code{\link{verb-PUT}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ok.R
\name{ok}
\alias{ok}
\title{check if a url is okay}
\usage{
ok(x, status = 200L, info = TRUE, verb = "head", ua_random = FALSE, ...)
}
\arguments{
\item{x}{either a URL as a character string, or an object of
class \link{HttpClient}}

\item{status}{(integer) one or more HTTP status codes, must be integers.
default: \code{200L}, since this is the most common signal
that a URL is okay, but there may be cases in which your URL
is okay if it's a \code{201L}, or some other status code.}

\item{info}{(logical) in the case of an error, do you want a
\code{message()} about it? Default: \code{TRUE}}

\item{verb}{(character) use "head" (default) or "get" HTTP verb
for the request. note that "get" will take longer as it returns a
body. however, "verb=get" may be your only option if a url
blocks head requests}

\item{ua_random}{(logical) use a random user agent string?
default: \code{TRUE}. if you set \code{useragent} curl option it will override
this setting. The random user agent string is pulled from a vector of
50 user agent strings generated from \code{charlatan::UserAgentProvider}
(by executing \code{replicate(30, UserAgentProvider$new()$user_agent())})}

\item{...}{args passed on to \link{HttpClient}}
}
\value{
a single boolean, if \code{TRUE} the URL is up and okay,
if \code{FALSE} it is down; but, see Details
}
\description{
check if a url is okay
}
\details{
We internally verify that status is an integer and
in the known set of HTTP status codes, and that info is a boolean

You may have to fiddle with the parameters to \code{ok()} as well as
curl options to get the "right answer". If you think you are getting
incorrectly getting \code{FALSE}, the first thing to do is to pass in
\code{verbose=TRUE} to \code{ok()}. That will give you verbose curl output and will
help determine what the issue may be. Here's some different scenarios:
\itemize{
\item the site blocks head requests: some sites do this, try \code{verb="get"}
\item it will be hard to determine a site that requires this, but it's
worth trying a random useragent string, e.g., \code{ok(useragent = "foobar")}
\item some sites are up and reachable but you could get a 403 Unauthorized
error, there's nothing you can do in this case other than having access
\item its possible to get a weird HTTP status code, e.g., LinkedIn gives
a 999 code, they're trying to prevent any programmatic access
}

A \code{FALSE} result may be incorrect depending on the use case. For example,
if you want to know if curl based scraping will work without fiddling with
curl options, then the \code{FALSE} is probably correct, but if you want to
fiddle with curl options, then first step would be to send \code{verbose=TRUE}
to see whats going on with any redirects and headers. You can set headers,
user agent strings, etc. to get closer to the request you want to know
about. Note that a user agent string is always passed by default, but it
may not be the one you want.
}
\examples{
\dontrun{
# 200
ok("https://www.google.com") 
# 200
ok("https://httpbin.org/status/200")
# more than one status
ok("https://www.google.com", status = c(200L, 202L))
# 404
ok("https://httpbin.org/status/404")
# doesn't exist
ok("https://stuff.bar")
# doesn't exist
ok("stuff")

# use get verb instead of head
ok("http://animalnexus.ca")
ok("http://animalnexus.ca", verb = "get")

# some urls will require a different useragent string
# they probably regex the useragent string
ok("https://doi.org/10.1093/chemse/bjq042")
ok("https://doi.org/10.1093/chemse/bjq042", verb = "get", useragent = "foobar")

# with random user agent's
## here, use a request hook to print out just the user agent string so 
## we can see what user agent string is being sent off
fun_ua <- function(request) {
  message(paste0("User-agent: ", request$options$useragent), sep = "\n")
}
z <- crul::HttpClient$new("https://doi.org/10.1093/chemse/bjq042", 
 hooks = list(request = fun_ua))
z
replicate(5, ok(z, ua_random=TRUE), simplify=FALSE)
## if you set useragent option it will override ua_random=TRUE
ok("https://doi.org/10.1093/chemse/bjq042", useragent="foobar", ua_random=TRUE)

# with HttpClient
z <- crul::HttpClient$new("https://httpbin.org/status/404", 
 opts = list(verbose = TRUE))
ok(z)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proxies.R
\name{proxies}
\alias{proxies}
\alias{proxy}
\title{proxy options}
\usage{
proxy(url, user = NULL, pwd = NULL, auth = "basic")
}
\arguments{
\item{url}{(character) URL, with scheme (http/https), domain and
port (must be numeric). required.}

\item{user}{(character) username, optional}

\item{pwd}{(character) password, optional}

\item{auth}{(character) authentication type, one of basic (default),
digest, digest_ie, gssnegotiate, ntlm, any or \code{NULL}. optional}
}
\description{
proxy options
}
\details{
See https://www.hidemyass.com/proxy for a list of proxies you
can use
}
\examples{
proxy("http://97.77.104.22:3128")
proxy("97.77.104.22:3128")
proxy("http://97.77.104.22:3128", "foo", "bar")
proxy("http://97.77.104.22:3128", "foo", "bar", auth = "digest")
proxy("http://97.77.104.22:3128", "foo", "bar", auth = "ntlm")

# socks
proxy("socks5://localhost:9050/", auth = NULL)

\dontrun{
# with proxy (look at request/outgoing headers)
# (res <- HttpClient$new(
#   url = "http://www.google.com",
#   proxies = proxy("http://97.77.104.22:3128")
# ))
# res$proxies
# res$get(verbose = TRUE)

# vs. without proxy (look at request/outgoing headers)
# (res2 <- HttpClient$new(url = "http://www.google.com"))
# res2$get(verbose = TRUE)


# Use authentication
# (res <- HttpClient$new(
#   url = "http://google.com",
#   proxies = proxy("http://97.77.104.22:3128", user = "foo", pwd = "bar")
# ))

# another example
# (res <- HttpClient$new(
#   url = "http://ip.tyk.nu/",
#   proxies = proxy("http://200.29.191.149:3128")
# ))
# res$get()$parse("UTF-8")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/async.R
\name{Async}
\alias{Async}
\title{Simple async client}
\value{
a list, with objects of class \code{\link[=HttpResponse]{HttpResponse()}}.
Responses are returned in the order they are passed in. We print the
first 10.
}
\description{
An async client to work with many URLs, but all with the same HTTP method
}
\details{
See \code{\link[=HttpClient]{HttpClient()}} for information on parameters.
}
\section{Failure behavior}{

HTTP requests mostly fail in ways that you are probably familiar with,
including when there's a 400 response (the URL not found), and when the
server made a mistake (a 500 series HTTP status code).

But requests can fail sometimes where there is no HTTP status code, and
no agreed upon way to handle it other than to just fail immediately.

When a request fails when using synchronous requests (see \link{HttpClient})
you get an error message that stops your code progression
immediately saying for example:
\itemize{
\item "Could not resolve host: https://foo.com"
\item "Failed to connect to foo.com"
\item "Resolving timed out after 10 milliseconds"
}

However, for async requests we don't want to fail immediately because
that would stop the subsequent requests from occurring. Thus, when
we find that a request fails for one of the reasons above we
give back a \link{HttpResponse} object just like any other response, and:
\itemize{
\item capture the error message and put it in the \code{content} slot of the
response object (thus calls to \code{content} and \code{parse()} work correctly)
\item give back a \code{0} HTTP status code. we handle this specially when testing
whether the request was successful or not with e.g., the \code{success()}
method
}
}

\section{R6 classes}{

This is an R6 class from the package \pkg{R6}. Find out more
about R6 at \url{https://r6.r-lib.org/}. After creating an instance of an R6
class (e.g., \code{x <- HttpClient$new(url = "https://httpbin.org")}) you can
access values and methods on the object \code{x}.
}

\examples{
\dontrun{
cc <- Async$new(
  urls = c(
    'https://httpbin.org/',
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?foo=bar'
  )
)
cc
(res <- cc$get())
res[[1]]
res[[1]]$url
res[[1]]$success()
res[[1]]$status_http()
res[[1]]$response_headers
res[[1]]$method
res[[1]]$content
res[[1]]$parse("UTF-8")
lapply(res, function(z) z$parse("UTF-8"))

# curl options/headers with async
urls = c(
 'https://httpbin.org/',
 'https://httpbin.org/get?a=5',
 'https://httpbin.org/get?foo=bar'
)
cc <- Async$new(urls = urls, 
  opts = list(verbose = TRUE),
  headers = list(foo = "bar")
)
cc
(res <- cc$get())

# using auth with async
dd <- Async$new(
  urls = rep('https://httpbin.org/basic-auth/user/passwd', 3),
  auth = auth(user = "foo", pwd = "passwd"),
  opts = list(verbose = TRUE)
)
dd
res <- dd$get()
res
vapply(res, function(z) z$status_code, double(1))
vapply(res, function(z) z$success(), logical(1))
lapply(res, function(z) z$parse("UTF-8"))

# failure behavior
## e.g. when a URL doesn't exist, a timeout, etc.
urls <- c("http://stuffthings.gvb", "https://foo.com", 
  "https://httpbin.org/get")
conn <- Async$new(urls = urls)
res <- conn$get()
res[[1]]$parse("UTF-8") # a failure
res[[2]]$parse("UTF-8") # a failure
res[[3]]$parse("UTF-8") # a success
}

## ------------------------------------------------
## Method `Async$get`
## ------------------------------------------------

\dontrun{
(cc <- Async$new(urls = c(
    'https://httpbin.org/',
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?foo=bar'
  )))
(res <- cc$get())
}

## ------------------------------------------------
## Method `Async$verb`
## ------------------------------------------------

\dontrun{
cc <- Async$new(
  urls = c(
    'https://httpbin.org/',
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?foo=bar'
  )
)
(res <- cc$verb('get'))
lapply(res, function(z) z$parse("UTF-8"))
}
}
\seealso{
Other async: 
\code{\link{AsyncQueue}},
\code{\link{AsyncVaried}},
\code{\link{HttpRequest}}
}
\concept{async}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{urls}}{(character) one or more URLs}

\item{\code{opts}}{any curl options}

\item{\code{proxies}}{named list of headers}

\item{\code{auth}}{an object of class \code{auth}}

\item{\code{headers}}{named list of headers}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{Async$print()}}
\item \href{#method-new}{\code{Async$new()}}
\item \href{#method-get}{\code{Async$get()}}
\item \href{#method-post}{\code{Async$post()}}
\item \href{#method-put}{\code{Async$put()}}
\item \href{#method-patch}{\code{Async$patch()}}
\item \href{#method-delete}{\code{Async$delete()}}
\item \href{#method-head}{\code{Async$head()}}
\item \href{#method-verb}{\code{Async$verb()}}
\item \href{#method-clone}{\code{Async$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for Async objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Async object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$new(urls, opts, proxies, auth, headers)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{urls}}{(character) one or more URLs}

\item{\code{opts}}{any curl options}

\item{\code{proxies}}{a \code{\link[=proxy]{proxy()}} object}

\item{\code{auth}}{an \code{\link[=auth]{auth()}} object}

\item{\code{headers}}{named list of headers}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Async} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get"></a>}}
\if{latex}{\out{\hypertarget{method-get}{}}}
\subsection{Method \code{get()}}{
execute the \code{GET} http verb for the \code{urls}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$get(path = NULL, query = list(), disk = NULL, stream = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) URL path, appended to the base URL}

\item{\code{query}}{(list) query terms, as a named list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
\code{NULL} (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
(cc <- Async$new(urls = c(
    'https://httpbin.org/',
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?foo=bar'
  )))
(res <- cc$get())
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-post"></a>}}
\if{latex}{\out{\hypertarget{method-post}{}}}
\subsection{Method \code{post()}}{
execute the \code{POST} http verb for the \code{urls}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$post(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  disk = NULL,
  stream = NULL,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) URL path, appended to the base URL}

\item{\code{query}}{(list) query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
\code{NULL} (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-put"></a>}}
\if{latex}{\out{\hypertarget{method-put}{}}}
\subsection{Method \code{put()}}{
execute the \code{PUT} http verb for the \code{urls}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$put(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  disk = NULL,
  stream = NULL,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) URL path, appended to the base URL}

\item{\code{query}}{(list) query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
\code{NULL} (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-patch"></a>}}
\if{latex}{\out{\hypertarget{method-patch}{}}}
\subsection{Method \code{patch()}}{
execute the \code{PATCH} http verb for the \code{urls}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$patch(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  disk = NULL,
  stream = NULL,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) URL path, appended to the base URL}

\item{\code{query}}{(list) query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
\code{NULL} (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-delete"></a>}}
\if{latex}{\out{\hypertarget{method-delete}{}}}
\subsection{Method \code{delete()}}{
execute the \code{DELETE} http verb for the \code{urls}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$delete(
  path = NULL,
  query = list(),
  body = NULL,
  encode = "multipart",
  disk = NULL,
  stream = NULL,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) URL path, appended to the base URL}

\item{\code{query}}{(list) query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
\code{NULL} (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-head"></a>}}
\if{latex}{\out{\hypertarget{method-head}{}}}
\subsection{Method \code{head()}}{
execute the \code{HEAD} http verb for the \code{urls}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$head(path = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) URL path, appended to the base URL}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-verb"></a>}}
\if{latex}{\out{\hypertarget{method-verb}{}}}
\subsection{Method \code{verb()}}{
execute any supported HTTP verb
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$verb(verb, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{verb}}{(character) a supported HTTP verb: get, post, put, patch, delete,
head.}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
cc <- Async$new(
  urls = c(
    'https://httpbin.org/',
    'https://httpbin.org/get?a=5',
    'https://httpbin.org/get?foo=bar'
  )
)
(res <- cc$verb('get'))
lapply(res, function(z) z$parse("UTF-8"))
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Async$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/httprequest.R
\name{HttpRequest}
\alias{HttpRequest}
\title{HTTP request object}
\description{
Create HTTP requests
}
\details{
This R6 class doesn't do actual HTTP requests as does
\code{\link[=HttpClient]{HttpClient()}} - it is for building requests to use for async HTTP
requests in \code{\link[=AsyncVaried]{AsyncVaried()}}

Note that you can access HTTP verbs after creating an \code{HttpRequest}
object, just as you can with \code{HttpClient}. See examples for usage.

Also note that when you call HTTP verbs on a \code{HttpRequest} object you
don't need to assign the new object to a variable as the new details
you've added are added to the object itself.

See \code{\link[=HttpClient]{HttpClient()}} for information on parameters.
}
\section{R6 classes}{

This is an R6 class from the package \pkg{R6}. Find out more
about R6 at \url{https://r6.r-lib.org/}. After creating an instance of an R6
class (e.g., \code{x <- HttpClient$new(url = "https://httpbin.org")}) you can
access values and methods on the object \code{x}.
}

\examples{
\dontrun{
x <- HttpRequest$new(url = "https://httpbin.org/get")
## note here how the HTTP method is shown on the first line to the right
x$get()

## assign to a new object to keep the output
z <- x$get()
### get the HTTP method
z$method()

(x <- HttpRequest$new(url = "https://httpbin.org/get")$get())
x$url
x$payload

(x <- HttpRequest$new(url = "https://httpbin.org/post"))
x$post(body = list(foo = "bar"))

HttpRequest$new(
  url = "https://httpbin.org/get",
  headers = list(
    `Content-Type` = "application/json"
  )
)
}

## ------------------------------------------------
## Method `HttpRequest$verb`
## ------------------------------------------------

z <- HttpRequest$new(url = "https://httpbin.org/get")
res <- z$verb('get', query = list(hello = "world"))
res$payload
}
\seealso{
\link{http-headers}, \link{writing-options}

Other async: 
\code{\link{AsyncQueue}},
\code{\link{AsyncVaried}},
\code{\link{Async}}
}
\concept{async}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{url}}{(character) a url}

\item{\code{opts}}{(list) named list of curl options}

\item{\code{proxies}}{a \code{\link[=proxy]{proxy()}} object}

\item{\code{auth}}{an \code{\link[=auth]{auth()}} object}

\item{\code{headers}}{(list) named list of headers, see \link{http-headers}}

\item{\code{handle}}{a \code{\link[=handle]{handle()}}}

\item{\code{progress}}{only supports \code{httr::progress()}, see \link{progress}}

\item{\code{payload}}{resulting payload after request}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{HttpRequest$print()}}
\item \href{#method-new}{\code{HttpRequest$new()}}
\item \href{#method-get}{\code{HttpRequest$get()}}
\item \href{#method-post}{\code{HttpRequest$post()}}
\item \href{#method-put}{\code{HttpRequest$put()}}
\item \href{#method-patch}{\code{HttpRequest$patch()}}
\item \href{#method-delete}{\code{HttpRequest$delete()}}
\item \href{#method-head}{\code{HttpRequest$head()}}
\item \href{#method-verb}{\code{HttpRequest$verb()}}
\item \href{#method-method}{\code{HttpRequest$method()}}
\item \href{#method-clone}{\code{HttpRequest$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for \code{HttpRequest} objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{HttpRequest} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$new(url, opts, proxies, auth, headers, handle, progress)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{url}}{(character) A url. One of \code{url} or \code{handle} required.}

\item{\code{opts}}{any curl options}

\item{\code{proxies}}{a \code{\link[=proxy]{proxy()}} object}

\item{\code{auth}}{an \code{\link[=auth]{auth()}} object}

\item{\code{headers}}{named list of headers, see \link{http-headers}}

\item{\code{handle}}{a \code{\link[=handle]{handle()}}}

\item{\code{progress}}{only supports \code{httr::progress()}, see \link{progress}}

\item{\code{urls}}{(character) one or more URLs}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{HttpRequest} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get"></a>}}
\if{latex}{\out{\hypertarget{method-get}{}}}
\subsection{Method \code{get()}}{
Define a GET request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$get(path = NULL, query = list(), disk = NULL, stream = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-post"></a>}}
\if{latex}{\out{\hypertarget{method-post}{}}}
\subsection{Method \code{post()}}{
Define a POST request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$post(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-put"></a>}}
\if{latex}{\out{\hypertarget{method-put}{}}}
\subsection{Method \code{put()}}{
Define a PUT request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$put(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-patch"></a>}}
\if{latex}{\out{\hypertarget{method-patch}{}}}
\subsection{Method \code{patch()}}{
Define a PATCH request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$patch(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-delete"></a>}}
\if{latex}{\out{\hypertarget{method-delete}{}}}
\subsection{Method \code{delete()}}{
Define a DELETE request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$delete(
  path = NULL,
  query = list(),
  body = NULL,
  disk = NULL,
  stream = NULL,
  encode = "multipart",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{query}}{query terms, as a named list}

\item{\code{body}}{body as an R list}

\item{\code{disk}}{a path to write to. if NULL (default), memory used.
See \code{\link[curl:curl_fetch]{curl::curl_fetch_disk()}} for help.}

\item{\code{stream}}{an R function to determine how to stream data. if
NULL (default), memory used. See \code{\link[curl:curl_fetch]{curl::curl_fetch_stream()}}
for help}

\item{\code{encode}}{one of form, multipart, json, or raw}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-head"></a>}}
\if{latex}{\out{\hypertarget{method-head}{}}}
\subsection{Method \code{head()}}{
Define a HEAD request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$head(path = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{URL path, appended to the base URL}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-verb"></a>}}
\if{latex}{\out{\hypertarget{method-verb}{}}}
\subsection{Method \code{verb()}}{
Use an arbitrary HTTP verb supported on this class
Supported verbs: get, post, put, patch, delete, head
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$verb(verb, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{verb}}{an HTTP verb supported on this class: get,
post, put, patch, delete, head. Also supports retry.}

\item{\code{...}}{curl options, only those in the acceptable set from
\code{\link[curl:curl_options]{curl::curl_options()}} except the following: httpget, httppost, post,
postfields, postfieldsize, and customrequest}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{z <- HttpRequest$new(url = "https://httpbin.org/get")
res <- z$verb('get', query = list(hello = "world"))
res$payload
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-method"></a>}}
\if{latex}{\out{\hypertarget{method-method}{}}}
\subsection{Method \code{method()}}{
Get the HTTP method (if defined)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$method()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(character) the HTTP method
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpRequest$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asyncvaried.R
\name{AsyncVaried}
\alias{AsyncVaried}
\title{Async client for different request types}
\value{
An object of class \code{AsyncVaried} with variables and methods.
\link{HttpResponse} objects are returned in the order they are passed in.
We print the first 10.
}
\description{
An async client to do many requests, each with different URLs, curl options,
etc.
}
\section{Failure behavior}{

HTTP requests mostly fail in ways that you are probably familiar with,
including when there's a 400 response (the URL not found), and when the
server made a mistake (a 500 series HTTP status code).

But requests can fail sometimes where there is no HTTP status code, and
no agreed upon way to handle it other than to just fail immediately.

When a request fails when using synchronous requests (see \link{HttpClient})
you get an error message that stops your code progression
immediately saying for example:
\itemize{
\item "Could not resolve host: https://foo.com"
\item "Failed to connect to foo.com"
\item "Resolving timed out after 10 milliseconds"
}

However, for async requests we don't want to fail immediately because
that would stop the subsequent requests from occurring. Thus, when
we find that a request fails for one of the reasons above we
give back a \link{HttpResponse} object just like any other response, and:
\itemize{
\item capture the error message and put it in the \code{content} slot of the
response object (thus calls to \code{content} and \code{parse()} work correctly)
\item give back a \code{0} HTTP status code. we handle this specially when testing
whether the request was successful or not with e.g., the \code{success()}
method
}
}

\section{R6 classes}{

This is an R6 class from the package \pkg{R6}. Find out more
about R6 at \url{https://r6.r-lib.org/}. After creating an instance of an R6
class (e.g., \code{x <- HttpClient$new(url = "https://httpbin.org")}) you can
access values and methods on the object \code{x}.
}

\examples{
\dontrun{
# pass in requests via ...
req1 <- HttpRequest$new(
  url = "https://httpbin.org/get",
  opts = list(verbose = TRUE),
  headers = list(foo = "bar")
)$get()
req2 <- HttpRequest$new(url = "https://httpbin.org/post")$post()

# Create an AsyncVaried object
out <- AsyncVaried$new(req1, req2)

# before you make requests, the methods return empty objects
out$status()
out$status_code()
out$content()
out$times()
out$parse()
out$responses()

# make requests
out$request()

# access various parts
## http status objects
out$status()
## status codes
out$status_code()
## content (raw data)
out$content()
## times
out$times()
## parsed content
out$parse()
## response objects
out$responses()

# use $verb() method to select http verb
method <- "post"
req1 <- HttpRequest$new(
  url = "https://httpbin.org/post",
  opts = list(verbose = TRUE),
  headers = list(foo = "bar")
)$verb(method)
req2 <- HttpRequest$new(url = "https://httpbin.org/post")$verb(method)
out <- AsyncVaried$new(req1, req2)
out
out$request()
out$responses()

# pass in requests in a list via .list param
reqlist <- list(
  HttpRequest$new(url = "https://httpbin.org/get")$get(),
  HttpRequest$new(url = "https://httpbin.org/post")$post(),
  HttpRequest$new(url = "https://httpbin.org/put")$put(),
  HttpRequest$new(url = "https://httpbin.org/delete")$delete(),
  HttpRequest$new(url = "https://httpbin.org/get?g=5")$get(),
  HttpRequest$new(
    url = "https://httpbin.org/post")$post(body = list(y = 9)),
  HttpRequest$new(
    url = "https://httpbin.org/get")$get(query = list(hello = "world"))
)

out <- AsyncVaried$new(.list = reqlist)
out$request()
out$status()
out$status_code()
out$content()
out$times()
out$parse()

# using auth with async
url <- "https://httpbin.org/basic-auth/user/passwd"
auth <- auth(user = "user", pwd = "passwd")
reqlist <- list(
  HttpRequest$new(url = url, auth = auth)$get(),
  HttpRequest$new(url = url, auth = auth)$get(query = list(a=5)),
  HttpRequest$new(url = url, auth = auth)$get(query = list(b=3))
)
out <- AsyncVaried$new(.list = reqlist)
out$request()
out$status()
out$parse()

# failure behavior
## e.g. when a URL doesn't exist, a timeout, etc.
reqlist <- list(
  HttpRequest$new(url = "http://stuffthings.gvb")$get(),
  HttpRequest$new(url = "https://httpbin.org")$head(),
  HttpRequest$new(url = "https://httpbin.org", 
   opts = list(timeout_ms = 10))$head()
)
(tmp <- AsyncVaried$new(.list = reqlist))
tmp$request()
tmp$responses()
tmp$parse("UTF-8")

# access intemediate redirect headers
dois <- c("10.7202/1045307ar", "10.1242/jeb.088898", "10.1121/1.3383963")
reqlist <- list(
  HttpRequest$new(url = paste0("https://doi.org/", dois[1]))$get(),
  HttpRequest$new(url = paste0("https://doi.org/", dois[2]))$get(),
  HttpRequest$new(url = paste0("https://doi.org/", dois[3]))$get()
)
tmp <- AsyncVaried$new(.list = reqlist)
tmp$request()
tmp
lapply(tmp$responses(), "[[", "response_headers_all")
}
}
\seealso{
Other async: 
\code{\link{AsyncQueue}},
\code{\link{Async}},
\code{\link{HttpRequest}}
}
\concept{async}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{AsyncVaried$print()}}
\item \href{#method-new}{\code{AsyncVaried$new()}}
\item \href{#method-request}{\code{AsyncVaried$request()}}
\item \href{#method-responses}{\code{AsyncVaried$responses()}}
\item \href{#method-requests}{\code{AsyncVaried$requests()}}
\item \href{#method-parse}{\code{AsyncVaried$parse()}}
\item \href{#method-status_code}{\code{AsyncVaried$status_code()}}
\item \href{#method-status}{\code{AsyncVaried$status()}}
\item \href{#method-content}{\code{AsyncVaried$content()}}
\item \href{#method-times}{\code{AsyncVaried$times()}}
\item \href{#method-clone}{\code{AsyncVaried$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for AsyncVaried objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new AsyncVaried object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$new(..., .list = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{..., .list}}{Any number of objects of class \code{\link[=HttpRequest]{HttpRequest()}},
must supply inputs to one of these parameters, but not both}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{AsyncVaried} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-request"></a>}}
\if{latex}{\out{\hypertarget{method-request}{}}}
\subsection{Method \code{request()}}{
Execute asynchronous requests
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$request()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing, responses stored inside object, though will print
messages if you choose verbose output
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-responses"></a>}}
\if{latex}{\out{\hypertarget{method-responses}{}}}
\subsection{Method \code{responses()}}{
List responses
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$responses()}\if{html}{\out{</div>}}
}

\subsection{Details}{
An S3 print method is used to summarise results. \link{unclass}
the output to see the list, or index to results, e.g., \verb{[1]}, \verb{[1:3]}
}

\subsection{Returns}{
a list of \code{HttpResponse} objects, empty list before
requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-requests"></a>}}
\if{latex}{\out{\hypertarget{method-requests}{}}}
\subsection{Method \code{requests()}}{
List requests
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$requests()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list of \code{HttpRequest} objects, empty list before
requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse"></a>}}
\if{latex}{\out{\hypertarget{method-parse}{}}}
\subsection{Method \code{parse()}}{
parse content
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$parse(encoding = "UTF-8")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{encoding}}{(character) the encoding to use in parsing.
default:"UTF-8"}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character vector, empty character vector before
requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status_code"></a>}}
\if{latex}{\out{\hypertarget{method-status_code}{}}}
\subsection{Method \code{status_code()}}{
Get HTTP status codes for each response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$status_code()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
numeric vector, empty numeric vector before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-status"></a>}}
\if{latex}{\out{\hypertarget{method-status}{}}}
\subsection{Method \code{status()}}{
List HTTP status objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$status()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list of \code{http_code} objects, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-content"></a>}}
\if{latex}{\out{\hypertarget{method-content}{}}}
\subsection{Method \code{content()}}{
Get raw content for each response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$content()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
raw list, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-times"></a>}}
\if{latex}{\out{\hypertarget{method-times}{}}}
\subsection{Method \code{times()}}{
curl request times
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$times()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list of named numeric vectors, empty list before requests made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AsyncVaried$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cookies.R
\name{cookies}
\alias{cookies}
\title{Working with cookies}
\description{
Working with cookies
}
\examples{
\dontrun{
x <- HttpClient$new(
  url = "https://httpbin.org",
  opts = list(
    cookie = "c=1;f=5",
    verbose = TRUE
  )
)
x

# set cookies
(res <- x$get("cookies"))
jsonlite::fromJSON(res$parse("UTF-8"))

(x <- HttpClient$new(url = "https://httpbin.org"))
res <- x$get("cookies/set", query = list(foo = 123, bar = "ftw"))
jsonlite::fromJSON(res$parse("UTF-8"))
curl::handle_cookies(handle = res$handle)

# reuse handle
res2 <- x$get("get", query = list(hello = "world"))
jsonlite::fromJSON(res2$parse("UTF-8"))
curl::handle_cookies(handle = res2$handle)

# DOAJ
x <- HttpClient$new(url = "https://doaj.org")
res <- x$get("api/v1/journals/f3f2e7f23d444370ae5f5199f85bc100",
  verbose = TRUE)
res$response_headers$`set-cookie`
curl::handle_cookies(handle = res$handle)
res2 <- x$get("api/v1/journals/9abfb36b06404e8a8566e1a44180bbdc",
  verbose = TRUE)

## reset handle
x$handle_pop()
## cookies no longer sent, as handle reset
res2 <- x$get("api/v1/journals/9abfb36b06404e8a8566e1a44180bbdc",
  verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{auth}
\alias{auth}
\title{Authentication}
\usage{
auth(user, pwd, auth = "basic")
}
\arguments{
\item{user}{(character) username, required. see Details.}

\item{pwd}{(character) password, required. see Details.}

\item{auth}{(character) authentication type, one of basic (default),
digest, digest_ie, gssnegotiate, ntlm, or any. required}
}
\description{
Authentication
}
\details{
Only supporting simple auth for now, OAuth later maybe.

For \code{user} and \code{pwd} you are required to pass in some value.
The value can be \code{NULL} to - which is equivalent to passing in an
empty string like \code{""} in \code{httr::authenticate}. You may want to pass
in \code{NULL} for both \code{user} and \code{pwd} for example if you are using
\code{gssnegotiate} auth type. See example below.
}
\examples{
auth(user = "foo", pwd = "bar", auth = "basic")
auth(user = "foo", pwd = "bar", auth = "digest")
auth(user = "foo", pwd = "bar", auth = "ntlm")
auth(user = "foo", pwd = "bar", auth = "any")

# gssnegotiate auth
auth(NULL, NULL, "gssnegotiate")

\dontrun{
# with HttpClient
(res <- HttpClient$new(
  url = "https://httpbin.org/basic-auth/user/passwd",
  auth = auth(user = "user", pwd = "passwd")
))
res$auth
x <- res$get()
jsonlite::fromJSON(x$parse("UTF-8"))

# with HttpRequest
(res <- HttpRequest$new(
  url = "https://httpbin.org/basic-auth/user/passwd",
  auth = auth(user = "user", pwd = "passwd")
))
res$auth
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verbs.R
\name{verb-HEAD}
\alias{verb-HEAD}
\title{HTTP verb info: HEAD}
\description{
The HEAD method asks for a response identical to that of a GET request,
but without the response body.
}
\section{The HEAD method}{

The HEAD method is identical to GET except that the
server MUST NOT send a message body in the response (i.e., the
response terminates at the end of the header section).  The server
SHOULD send the same header fields in response to a HEAD request as
it would have sent if the request had been a GET, except that the
payload header fields MAY be omitted.  This method can
be used for obtaining metadata about the selected representation
without transferring the representation data and is often used for
testing hypertext links for validity, accessibility, and recent
modification.

See \url{https://tools.ietf.org/html/rfc7231#section-4.3.2} for further
details.
}

\examples{
\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")
x$head()
}

}
\references{
\url{https://tools.ietf.org/html/rfc7231#section-4.3.2}
}
\seealso{
\link{crul-package}

Other verbs: 
\code{\link{verb-DELETE}},
\code{\link{verb-GET}},
\code{\link{verb-PATCH}},
\code{\link{verb-POST}},
\code{\link{verb-PUT}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verbs.R
\name{verb-PUT}
\alias{verb-PUT}
\title{HTTP verb info: PUT}
\description{
The PUT method replaces all current representations of the target
resource with the request payload.
}
\section{The PUT method}{

The PUT method requests that the state of the target resource be
created or replaced with the state defined by the representation
enclosed in the request message payload.  A successful PUT of a given
representation would suggest that a subsequent GET on that same
target resource will result in an equivalent representation being
sent in a 200 (OK) response. However, there is no guarantee that
such a state change will be observable, since the target resource
might be acted upon by other user agents in parallel, or might be
subject to dynamic processing by the origin server, before any
subsequent GET is received.  A successful response only implies that
the user agent's intent was achieved at the time of its processing by
the origin server.

If the target resource does not have a current representation and the
PUT successfully creates one, then the origin server MUST inform the
user agent by sending a 201 (Created) response.  If the target
resource does have a current representation and that representation
is successfully modified in accordance with the state of the enclosed
representation, then the origin server MUST send either a 200 (OK) or
a 204 (No Content) response to indicate successful completion of the
request.

See \url{https://tools.ietf.org/html/rfc7231#section-4.3.4} for further
details.
}

\examples{
\dontrun{
x <- HttpClient$new(url = "https://httpbin.org")
x$put(path = 'put', body = list(foo = "bar"))
}

}
\references{
\url{https://tools.ietf.org/html/rfc7231#section-4.3.4}
}
\seealso{
\link{crul-package}

Other verbs: 
\code{\link{verb-DELETE}},
\code{\link{verb-GET}},
\code{\link{verb-HEAD}},
\code{\link{verb-PATCH}},
\code{\link{verb-POST}}
}
\concept{verbs}
