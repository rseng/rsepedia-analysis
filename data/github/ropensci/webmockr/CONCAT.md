---
title: "none"
output: html_fragment
---

The very very short version is: `webmockr` helps you stub HTTP requests so you 
don't have to repeat yourself.

**More details**

You tell `webmockr` what HTTP request you want to match against and if it sees a 
request matching your criteria it doesn't actually do the HTTP request. Instead,
it gives back the same object you would have gotten back with a real request, but 
only with the bits it knows about. For example, we can't give back the actual 
data you'd get from a real HTTP request as the request wasn't performed.

In addition, if you set an expectation of what `webmockr` should return, we 
return that. For example, if you expect a request to return a 418 error 
(I'm a Teapot), then that's what you'll get.

**What you can match against**

* HTTP method (required)

Plus any single or combination of the following:

* URI
    * Right now, we can match directly against URI's, and with regex URI patterns. 
  Eventually, we will support RFC 6570 URI templates. 
    * We normalize URI paths so that URL encoded things match 
  URL un-encoded things (e.g. `hello world` to `hello%20world`)
* Query parameters
    * We normalize query parameter values so that URL encoded things match 
  URL un-encoded things (e.g. `message = hello world` to 
  `message = hello%20world`)
* Request headers
    * We normalize headers and treat all forms of same headers as equal. For 
  example, the following two sets of headers are equal:
        * `list(H1 = "value1", content_length = 123, X_CuStOm_hEAder = "foo")`
        * `list(h1 = "value1", "Content-Length" = 123, "x-cuSTOM-HeAder" = "foo")`
* Request body

**Real HTTP requests**

There's a few scenarios to think about when using `webmockr`:

After doing

```r
library(webmockr)
```

`webmockr` is loaded but not turned on. At this point `webmockr` doesn't 
change anythning.

Once you turn on `webmockr` like 

```r
webmockr::enable()
```

`webmockr` will now by default not allow real HTTP requests from the http 
libraries that adapters are loaded for (right now only `crul`).

You can optionally allow real requests via `webmockr_allow_net_connect()`, and
disallow real requests via `webmockr_disable_net_connect()`. You can check 
whether you are allowing real requests with `webmockr_net_connect_allowed()`.

Certain kinds of real HTTP requests allowed: We don't suppoprt this yet, 
but you can allow localhost HTTP requests with the `allow_localhost` parameter
in the `webmockr_configure()` function. 

**Storing actual HTTP responses**

`webmockr` doesn't do that. Check out [vcr][]

[vcr]: https://github.com/ropensci/vcr
webmockr
========



<!-- NOTE: run `make readme` to generate the README.md -->

[![cran checks](https://cranchecks.info/badges/worst/webmockr)](https://cranchecks.info/pkgs/webmockr)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/webmockr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/webmockr/actions/)
[![codecov](https://codecov.io/gh/ropensci/webmockr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/webmockr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/webmockr)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/webmockr)](https://cran.r-project.org/package=webmockr)


R library for stubbing and setting expectations on HTTP requests.

Port of the Ruby gem [webmock](https://github.com/bblimke/webmock)

<details> <summary><strong>How it works in detail</strong></summary> <p>





<p>The very very short version is: <code>webmockr</code> helps you stub HTTP requests so you don’t have to repeat yourself.</p>
<p><strong>More details</strong></p>
<p>You tell <code>webmockr</code> what HTTP request you want to match against and if it sees a request matching your criteria it doesn’t actually do the HTTP request. Instead, it gives back the same object you would have gotten back with a real request, but only with the bits it knows about. For example, we can’t give back the actual data you’d get from a real HTTP request as the request wasn’t performed.</p>
<p>In addition, if you set an expectation of what <code>webmockr</code> should return, we return that. For example, if you expect a request to return a 418 error (I’m a Teapot), then that’s what you’ll get.</p>
<p><strong>What you can match against</strong></p>
<ul>
<li>HTTP method (required)</li>
</ul>
<p>Plus any single or combination of the following:</p>
<ul>
<li>URI
<ul>
<li>Right now, we can match directly against URI’s, and with regex URI patterns. Eventually, we will support RFC 6570 URI templates.</li>
<li>We normalize URI paths so that URL encoded things match URL un-encoded things (e.g. <code>hello world</code> to <code>hello%20world</code>)</li>
</ul></li>
<li>Query parameters
<ul>
<li>We normalize query parameter values so that URL encoded things match URL un-encoded things (e.g. <code>message = hello world</code> to <code>message = hello%20world</code>)</li>
</ul></li>
<li>Request headers
<ul>
<li>We normalize headers and treat all forms of same headers as equal. For example, the following two sets of headers are equal:
<ul>
<li><code>list(H1 = &quot;value1&quot;, content_length = 123, X_CuStOm_hEAder = &quot;foo&quot;)</code></li>
<li><code>list(h1 = &quot;value1&quot;, &quot;Content-Length&quot; = 123, &quot;x-cuSTOM-HeAder&quot; = &quot;foo&quot;)</code></li>
</ul></li>
</ul></li>
<li>Request body</li>
</ul>
<p><strong>Real HTTP requests</strong></p>
<p>There’s a few scenarios to think about when using <code>webmockr</code>:</p>
<p>After doing</p>
<pre class="r"><code>library(webmockr)</code></pre>
<p><code>webmockr</code> is loaded but not turned on. At this point <code>webmockr</code> doesn’t change anythning.</p>
<p>Once you turn on <code>webmockr</code> like</p>
<pre class="r"><code>webmockr::enable()</code></pre>
<p><code>webmockr</code> will now by default not allow real HTTP requests from the http libraries that adapters are loaded for (right now only <code>crul</code>).</p>
<p>You can optionally allow real requests via <code>webmockr_allow_net_connect()</code>, and disallow real requests via <code>webmockr_disable_net_connect()</code>. You can check whether you are allowing real requests with <code>webmockr_net_connect_allowed()</code>.</p>
<p>Certain kinds of real HTTP requests allowed: We don’t suppoprt this yet, but you can allow localhost HTTP requests with the <code>allow_localhost</code> parameter in the <code>webmockr_configure()</code> function.</p>
<p><strong>Storing actual HTTP responses</strong></p>
<p><code>webmockr</code> doesn’t do that. Check out <a href="https://github.com/ropensci/vcr">vcr</a></p>

</p></details>

## Features

* Stubbing HTTP requests at low http client lib level
* Setting and verifying expectations on HTTP requests
* Matching requests based on method, URI, headers and body
* Support for `testthat` via [vcr][]
* Can be used for testing or outside of a testing context

## Supported HTTP libraries

* [crul](https://github.com/ropensci/crul)
* [httr](https://github.com/r-lib/httr)

## Install

from cran


```r
install.packages("webmockr")
```

Dev version


```r
remotes::install_github("ropensci/webmockr")
```


```r
library(webmockr)
```

## Enable webmockr


```r
webmockr::enable()
#> CrulAdapter enabled!
#> HttrAdapter enabled!
```

## Inside a test framework


```r
library(crul)
library(testthat)

# make a stub
stub_request("get", "https://httpbin.org/get") %>%
   to_return(body = "success!", status = 200)
#> <webmockr stub> 
#>   method: get
#>   uri: https://httpbin.org/get
#>   with: 
#>     query: 
#>     body: 
#>     request_headers: 
#>   to_return: 
#>   - status: 200
#>     body: success!
#>     response_headers: 
#>     should_timeout: FALSE
#>     should_raise: FALSE

# check that it's in the stub registry
stub_registry()
#> <webmockr stub registry> 
#>  Registered Stubs
#>    GET: https://httpbin.org/get   | to_return:   with body "success!"  with status 200

# make the request
z <- crul::HttpClient$new(url = "https://httpbin.org")$get("get")

# run tests (nothing returned means it passed)
expect_is(z, "HttpResponse")
expect_equal(z$status_code, 200)
expect_equal(z$parse("UTF-8"), "success!")
```




## Outside a test framework


```r
library(crul)
```

### Stubbed request based on uri only and with the default response


```r
stub_request("get", "https://httpbin.org/get")
#> <webmockr stub> 
#>   method: get
#>   uri: https://httpbin.org/get
#>   with: 
#>     query: 
#>     body: 
#>     request_headers: 
#>   to_return:
```


```r
x <- HttpClient$new(url = "https://httpbin.org")
x$get('get')
#> <crul response> 
#>   url: https://httpbin.org/get
#>   request_headers: 
#>     User-Agent: libcurl/7.64.1 r-curl/4.3 crul/1.1.0
#>     Accept-Encoding: gzip, deflate
#>     Accept: application/json, text/xml, application/xml, */*
#>   response_headers: 
#>   status: 200
```

set return objects


```r
stub_request("get", "https://httpbin.org/get") %>%
  wi_th(
    query = list(hello = "world")) %>%
    to_return(status = 418)
#> <webmockr stub> 
#>   method: get
#>   uri: https://httpbin.org/get
#>   with: 
#>     query: hello=world
#>     body: 
#>     request_headers: 
#>   to_return: 
#>   - status: 418
#>     body: 
#>     response_headers: 
#>     should_timeout: FALSE
#>     should_raise: FALSE
```


```r
x$get('get', query = list(hello = "world"))
#> <crul response> 
#>   url: https://httpbin.org/get
#>   request_headers: 
#>     User-Agent: libcurl/7.64.1 r-curl/4.3 crul/1.1.0
#>     Accept-Encoding: gzip, deflate
#>     Accept: application/json, text/xml, application/xml, */*
#>   response_headers: 
#>   status: 418
```

### Stubbing requests based on method, uri and query params


```r
stub_request("get", "https://httpbin.org/get") %>%
  wi_th(query = list(hello = "world"), 
        headers = list('User-Agent' = 'libcurl/7.51.0 r-curl/2.6 crul/0.3.6', 
                       'Accept-Encoding' = "gzip, deflate"))
#> <webmockr stub> 
#>   method: get
#>   uri: https://httpbin.org/get
#>   with: 
#>     query: hello=world
#>     body: 
#>     request_headers: User-Agent=libcurl/7.51.0 r-cur..., Accept-Encoding=gzip, deflate
#>   to_return:
```


```r
stub_registry()
#> <webmockr stub registry> 
#>  Registered Stubs
#>    GET: https://httpbin.org/get 
#>    GET: https://httpbin.org/get?hello=world   | to_return:    with status 418 
#>    GET: https://httpbin.org/get?hello=world   with headers {"User-Agent":"libcurl/7.51.0 r-curl/2.6 crul/0.3.6","Accept-Encoding":"gzip, deflate"}
```


```r
x <- HttpClient$new(url = "https://httpbin.org")
x$get('get', query = list(hello = "world"))
#> <crul response> 
#>   url: https://httpbin.org/get
#>   request_headers: 
#>     User-Agent: libcurl/7.64.1 r-curl/4.3 crul/1.1.0
#>     Accept-Encoding: gzip, deflate
#>     Accept: application/json, text/xml, application/xml, */*
#>   response_headers: 
#>   status: 418
```

### Stubbing requests and set expectation of a timeout


```r
stub_request("post", "https://httpbin.org/post") %>% to_timeout()
#> <webmockr stub> 
#>   method: post
#>   uri: https://httpbin.org/post
#>   with: 
#>     query: 
#>     body: 
#>     request_headers: 
#>   to_return: 
#>   - status: 
#>     body: 
#>     response_headers: 
#>     should_timeout: TRUE
#>     should_raise: FALSE
x <- HttpClient$new(url = "https://httpbin.org")
x$post('post')
#> Error: Request Timeout (HTTP 408).
#>  - The client did not produce a request within the time that the server was prepared to wait. The client MAY repeat the request without modifications at any later time.
```

### Stubbing requests and set HTTP error expectation


```r
library(fauxpas)
stub_request("get", "https://httpbin.org/get?a=b") %>% to_raise(HTTPBadRequest)
#> <webmockr stub> 
#>   method: get
#>   uri: https://httpbin.org/get?a=b
#>   with: 
#>     query: 
#>     body: 
#>     request_headers: 
#>   to_return: 
#>   - status: 
#>     body: 
#>     response_headers: 
#>     should_timeout: FALSE
#>     should_raise: HTTPBadRequest
x <- HttpClient$new(url = "https://httpbin.org")
x$get('get', query = list(a = "b"))
#> Error: Bad Request (HTTP 400).
#>  - The request could not be understood by the server due to malformed syntax. The client SHOULD NOT repeat the request without modifications.
```

## httr integration


```r
library(webmockr)
library(httr)
#> 
#> Attaching package: 'httr'
#> The following object is masked from 'package:crul':
#> 
#>     handle

# turn on httr mocking
httr_mock()
```


```r
# no stub found
GET("https://httpbin.org/get")
#> Error: Real HTTP connections are disabled.
#> Unregistered request:
#>   GET https://httpbin.org/get   with headers {Accept: application/json, text/xml, application/xml, */*}
#> 
#> You can stub this request with the following snippet:
#> 
#>    stub_request('get', uri = 'https://httpbin.org/get') %>%
#>      wi_th(
#>        headers = list('Accept' = 'application/json, text/xml, application/xml, */*')
#>      )
#> ============================================================
```

make a stub


```r
stub_request('get', uri = 'https://httpbin.org/get') %>%
  wi_th(
    headers = list('Accept' = 'application/json, text/xml, application/xml, */*')
  ) %>%
  to_return(status = 418, body = "I'm a teapot!!!", headers = list(im_a = "teapot"))
#> <webmockr stub> 
#>   method: get
#>   uri: https://httpbin.org/get
#>   with: 
#>     query: 
#>     body: 
#>     request_headers: Accept=application/json, te...
#>   to_return: 
#>   - status: 418
#>     body: I'm a teapot!!!
#>     response_headers: im_a=teapot
#>     should_timeout: FALSE
#>     should_raise: FALSE
```

now returns mocked response



```r
(res <- GET("https://httpbin.org/get"))
res$status_code
#> [1] 418
res$headers
#> $im_a
#> [1] "teapot"
```

## Writing to disk

Write to a file before mocked request




```r
## make a temp file
f <- tempfile(fileext = ".json")
## write something to the file
cat("{\"hello\":\"world\"}\n", file = f)
readLines(f)
#> [1] "{\"hello\":\"world\"}"
## make the stub
invisible(stub_request("get", "https://httpbin.org/get") %>% 
  to_return(body = file(f)))
## make a request
out <- HttpClient$new("https://httpbin.org/get")$get(disk = f)
readLines(file(f))
#> [1] "{\"hello\":\"world\"}"
```

OR - you can use `mock_file()` to have `webmockr` handle file and contents


```r
g <- tempfile(fileext = ".json")
## make the stub
invisible(stub_request("get", "https://httpbin.org/get") %>% 
  to_return(body = mock_file(g, "{\"hello\":\"mars\"}\n")))
## make a request
out <- crul::HttpClient$new("https://httpbin.org/get")$get(disk = g)
readLines(out$content)
#> [1] "{\"hello\":\"world\"}"
```

Writing to disk is supported in both `crul` and `httr`

## Many requests in a row

e.g., many redirects, then a final successful request


```r
webmockr::enable()
library(crul)
library(fauxpas)

z <- stub_request("get", "https://httpbin.org/get")
to_return(z, status = 200, body = "foobar", headers = list(a = 5))
to_return(z, status = 200, body = "bears", headers = list(b = 6))
to_raise(z, HTTPBadRequest)
z

con <- crul::HttpClient$new(url = "https://httpbin.org")
# the first to_return()
first <- con$get("get")
first
first$parse("UTF-8")
# the second to_return()
second <- con$get("get")
second
second$parse("UTF-8")
# the third to_return() - fails as specified
third <- con$get("get")
```

Note that subsequent requests past the number of responses given with `to_return()`/etc.
simply gives the last response you specified. Although if you set a `to_timeout` or 
`to_raise` this feature won't happen since you fail out.


## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Aaron Wolen](https://github.com/aaronwolen)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/webmockr/issues).
* License: MIT
* Get citation information for `webmockr` in R doing `citation(package = 'webmockr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[vcr]: https://github.com/ropensci/vcr
webmockr 0.8.0
==============

### NEW FEATURES

* `enable()` and the `enable()` method on the `Adapter` R6 class gain new parameter `quiet` to toggle whether messages are printed or not  (#112)

### MINOR IMPROVEMENTS

* to re-create http response objects for both httr and crul we were using the url from the request object; now we use the url from the response object, BUT if there is no url in the response object we fall back to using the url from the request object (#110) (#113)
* improve docs: add further explanation to manual files for both `to_raise()` and `to_return()` to explain the differenc between them and when you may want to use them (#100)


webmockr 0.7.4
==============

### MINOR IMPROVEMENTS

* to support vcr being able to recreate httr objects fully (see github issue ropensci/vcr#132) we needed to handle additional parts of httr request objects: fields and output - with this change vcr should return objects much closer to what real httr requests return (#109)

### BUG FIXES

* bug fix + improvement: fixes for simple authentication - `wi_th()` now supports `basic_auth` to mock basic authentication either with `crul::auth()` or `httr::authenticate()` (#108)


webmockr 0.7.0
==============

### NEW FEATURES

* Gains ability to define more than 1 returned HTTP response, and the order in which the HTTP responses are returned. The idea is from the Ruby webmock library, but the implementation is different because the Ruby and R languages are very different. You can give more than one `to_return()` one creating a stub, or if you want to return the same response each time, you can use the new `times` parameter within `to_return()`. As a related use case (#31) you can mock http retry's using this new feature (#10) (#32) (#101)
* Gains new function `webmockr_reset()` to be able to reset stub registry and request registry in one function call (#97) (#101)
* Gains support for mocking simple authentication. `wi_th()` now accepts `basic_auth` in addition to query, body, and headers. Note that authentication type is ignored (#103)

### MINOR IMPROVEMENTS

* change to how URI's are matched in `stub_request()`: we weren't allowing matching URI's without schemes; you can now do that. In addition, webmockr can match URI's without the "http" scheme, but does not match if the scheme is "https". See `UriPattern` for more (#102)
* another change to how URI's are matched: now query params compared separately to the URI; note that regex not allowed in query params (#104) - And now query parameters are compared with the same code both when regex uri is used and when it is not (#107)
* URI matching for stubs is now done only on the URI's themselves; that is, query parameters are removed before comparison, so only the base url with http scheme, plus paths, are compared (#107)
* wasn't sure `write_disk_path` behavior was correct when using httr, seems to be working, added tests for it (#79)
* values for query parameters given to `wi_th()` are now all coerced to character class to make sure that all comparisons of stubs and requests are done with the same class (character) (#107)

### BUG FIXES

* fix for `uri_regex` usage in `stub_request()`: no longer curl escape the `uri_regex` given, only escape a non-regex uri (#106)


webmockr 0.6.2
==============

* change to `CrulAdapter`: do not use `normalizePath` on the `write_disk_path` path so that relative paths are not changed to full paths - added tests for this (#95) (#96)


webmockr 0.6.0
==============

### NEW FEATURES

* new `Adapter` class to consolidate common code for the `HttrAdapter` and `CrulAdapter` classes, which inherit from `Adapter`; not a user facing change (#87)
* pkgdown documentation site gains grouping of functions to help the user navigate the package: see https://docs.ropensci.org/webmockr/reference/ (#93)

### MINOR IMPROVEMENTS

* now correctly fails with informative message when `write_disk_path` is `NULL` when the user is trying to write to disk while using webmockr (#78)
* improve README construction; use html child for the details section (#81)
* fix matching stub matching for bodies when bodies are JSON encoded (#82) 
* when vcr was loaded real HTTP requests were being performed twice when they should have only been performed once (#91) (#92)

### BUG FIXES

* fix for `set_body()` method in the `Response` class - handle cases where user writing to disk and not, and handle raw bytes correctly (#80)
* fix to `to_s()` method in `StubbedRequest` class - was formatting query parameters incorrectly (#83)
* fix to `BodyPattern` class to handle upload objects in a list; related issue fixed where `wi_th()` parameter `body` was not handling upload objects (#84) (#85)
* httr requests were failing when vcr loaded, but with no cassette inserted; fixed `handle_request()` to skip vcr-related code unless a cassette is inserted (#86) (#88)


webmockr 0.5.0
==============

### NEW FEATURES

* `webmockr` now supports mocking writing to disk. TLDR: see `?mocking-disk-writing` to get started - That is, both of the major high level http clients in R, crul and httr, support writing directly to disk (rather than the user manually getting the http response and writing it to disk). supporting this required quite a bit of work, both in code and in thinking about how to support the various scenarios in which users can find themselves when dealing with writing to disk - Please get in touch if you have problems with this (#57) (#76)
* gains `request_registry_clear()` method to easily clear all requests in the request registry (#75)

### MINOR IMPROVEMENTS

* better docs for R6 classes with R6 support in new roxygen2 version on cran (#77)
* httr simple auth was being ignored - its now supported (simple auth with crul already worked) (#74)

### BUG FIXES

* fix to handle raw responses that can not be converted to character, such as images; needed due to issue https://github.com/ropensci/vcr/issues/112 (#72) (#73)


webmockr 0.4.0
==============

### MINOR IMPROVEMENTS

* fix link to http testing book, change ropensci to ropenscilabs (#67)
* fixes to request matching: single match types working now (e.g., just match on query, or just on headers); in addition, header matching now works; added examples of single match types (#68) (#69)

### BUG FIXES

* fix stub specification within crul and httr adapters; typo in setting headers (#70)


webmockr 0.3.4
==============

### DEFUNCT

* underscore methods `to_return_()` and `wi_th_()` are defunct (#60) (#64)

### NEW FEATURES

* `to_return()` gains parameter `.list` (#60) (#64)

### MINOR IMPROVEMENTS

* typo fixes (#62) thanks @Bisaloo !
* improved the print method for stubs, found in `StubbedRequest`, to have better behavior for very long strings such as in headers and bodies (#63)

### BUG FIXES

* fix date in mocked `httr` response object to match the date format that `httr` uses in real HTTP requests (#58) (#61) via <https://github.com/ropensci/vcr/issues/91>
* fix response headers in mocked `httr` response objects. `httr` makes the list of headers insensitive to case, so we now use that function from the package (#59) (#61)
* `to_return()` and `wi_th()` drop use of the `lazyeval` package and fall back to using the simple `list(...)` - fixes problem where creating stubs was failing within `test_that()` blocks due to some weird lazy eval conflicts (i think) (#60) (#64) thanks @karawoo !


webmockr 0.3.0
==============

### MINOR IMPROVEMENTS

* returned mocked response headers were retaining case that the user gave - whereas they should be all lowercased to match the output in `crul` and `httr`. now fixed. (#49) thanks @hlapp
* returned mocked response headers were not all of character class, but depended on what class was given by the user on creating the stub. this is now fixed, returning all character class values for response headers (#48) thanks @hlapp
* skip tests that require `vcr` if `vcr` is not available (#53)
* internal change to crul adapter to produce the same http response as a new version of crul returns - adds a `response_headers_all` slot  (#51) (#54)


webmockr 0.2.9
==============

### MINOR IMPROVEMENTS

* make `request_registry()` and `stub_registry()` print methods more similar to avoid confusion for users (#35)
* update docs for `enable`/`disable` to indicate that `crul` and `httr` supported (#46) (related to #45)
* wrap httr adapter examples in `requireNamespace` so only run when httr available
* clean up `.onLoad` call, removing commented out code, and add note about creating adapter objects does not load crul and httr packages

### BUG FIXES

* fix to `enable()` and `disable()` methods. even though `httr` is in Suggests, we were loading all adapters (crul, httr) with `stop` when the package was not found. We now give a message and skip when a package not installed. In addition, we `enable()` and `disable()` gain an `adapter` parameter to indicate which package you want to enable or disable. If `adapter` not given we attempt all adapters. Note that this bug shouldn't have affected `vcr` users as `httr` is in Imports in that package, so you'd have to have `httr` installed   (#45) thanks to @maelle for uncovering the problem


webmockr 0.2.8
==============

### NEW FEATURES

* Added support for integration with package `httr`; see `HttrAdapter` for the details; `webmockr` now integrates with two HTTP R packages: `crul` and `httr` (#43) (#44)
* Along with `httr` integration is a new method `httr_mock()` to turn on mocking for `httr`; and two methods `build_httr_response` and `build_httr_request` meant for internal use


webmockr 0.2.6
==============

### NEW FEATURES

* Added support for integration with package `vcr` (now on CRAN) for doing HTTP request caching


webmockr 0.2.4
==============

### NEW FEATURES

* New function `enabled()` to ask if `webmockr` is enabled, gives a
boolean
* `wi_th()` gains new parameter `.list` as an escape hatch to avoid
NSE. examples added in the `wi_th` man file to clarify its use

### MINOR IMPROVEMENTS

* matching by request body was not supported, it now is; added examples
of matching on request body, see `?stub_request`  (#36)
* make sure that the adapter for `crul` handles all types of matches (#29)
* removed all internal usage of pipes in the package. still exporting
pipe for users (#30)
* fixed internals to give vcr error when vcr loaded - for future release
with vcr support (#34)
* require newest `crul` version

### BUG FIXES

* Error messages with the suggest stub were not giving bodies. They 
now give bodies if needed along with method, uri, headers, query (#37)
* Fixed `Response` class that was not dealing with capitalization 
correctly


webmockr 0.2.0
==============

### NEW FEATURES

* New function `to_raise()` to say that a matched response should return a certain exception, currently `to_raise` accepts error classes from the `fauxpas` package (#9)
* New function `to_timeout()` to say that a matched response should return a timeout. This is a special case of `to_raise` to easily do a timeout expectation (#11)
* New function `request_registry()` to list requests in the request registry (#23)
* package `crul` moved to Imports from Suggests as it's the only http client supported for now. will move back to Suggests once we support at least one other http client
* `webmockr_configure()` changes: `turn_on` has been removed; `allow_net_connect` and `allow_localhost` were ignored before, but are now used and are now set to `FALSE` by default; fixed usage of `allow` which now accepts character vector of URLs instead of a boolean; the following correctly marked as being ignored for now until fixed `net_http_connect_on_start`, `show_stubbing_instructions`, `query_values_notation`, `show_body_diff` (#19) (#21)
* `webmockr_disable_net_connect()` now accepts an `allow` parameter to disable all other connections except those URLs given in `allow`
* `webmockr_net_connect_allowed()` now accepts a `uri` parameter to test if a URI/URL is allowed

### MINOR IMPROVEMENTS

* Fixed printed stub statement when printed to the console - we weren't including headers accurately (#18)
* Added examples to the `stub_registry()` and `stub_registry_clea()` manual files (#24)
* internal methods `build_crul_request` and `build_crul_response` moved outside of the `CrulAdapter` class so that they can be accesed like `webmockr::` in other packages
* `enable()` and `disable()` now return booleans invisibly
* General improvements to documentation throughout
* Added linting of user inputs to the `to_return()` method, and docs details on what to input to the method
* Added linting of user inputs to the `wi_th()` method, and docs details on what to input to the method

### BUG FIXES

* Fixed option `allow_localhost`, which wasn't actually workin before (#25)

### DEPRECATED AND DEFUNCT

* `webmockr_enable()` and `webmockr_disable` are now defunct. Use `webmockr::enable()` and `webmockr::disable()` instead



webmockr 0.1.0
==============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local macOS install, R 4.0.4 Patched
* ubuntu (on GitHub Actions), R 4.0.4
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

I have checked the 13 reverse dependencies, and no problems were found. See (<https://github.com/ropensci/webmockr/tree/master/revdep>).

---

In this version, a function gains a parameter, a fix is made and documentation is improved.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/webmockr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/webmockr.git`
* Make sure to track progress upstream (i.e., on our version of `webmockr` at `ropensci/webmockr`) by doing `git remote add upstream https://github.com/ropensci/webmockr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/webmockr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

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
|date     |2021-03-12                                  |

# Dependencies

|package  |old   |new      |Δ  |
|:--------|:-----|:--------|:--|
|webmockr |0.7.4 |0.7.5.95 |*  |
|crul     |NA    |1.1.0    |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*