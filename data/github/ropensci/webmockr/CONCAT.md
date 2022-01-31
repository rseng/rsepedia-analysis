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

*Wow, no problems at all. :)**Wow, no problems at all. :)*webmockr
========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE
)
```

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

```{r details, child="details.html"}
```

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

```{r eval=FALSE}
install.packages("webmockr")
```

Dev version

```{r eval=FALSE}
remotes::install_github("ropensci/webmockr")
```

```{r}
library(webmockr)
```

## Enable webmockr

```{r}
webmockr::enable()
```

## Inside a test framework

```{r}
library(crul)
library(testthat)

# make a stub
stub_request("get", "https://httpbin.org/get") %>%
   to_return(body = "success!", status = 200)

# check that it's in the stub registry
stub_registry()

# make the request
z <- crul::HttpClient$new(url = "https://httpbin.org")$get("get")

# run tests (nothing returned means it passed)
expect_is(z, "HttpResponse")
expect_equal(z$status_code, 200)
expect_equal(z$parse("UTF-8"), "success!")
```

```{r echo=FALSE}
stub_registry_clear()
```


## Outside a test framework

```{r}
library(crul)
```

### Stubbed request based on uri only and with the default response

```{r}
stub_request("get", "https://httpbin.org/get")
```

```{r}
x <- HttpClient$new(url = "https://httpbin.org")
x$get('get')
```

set return objects

```{r}
stub_request("get", "https://httpbin.org/get") %>%
  wi_th(
    query = list(hello = "world")) %>%
    to_return(status = 418)
```

```{r}
x$get('get', query = list(hello = "world"))
```

### Stubbing requests based on method, uri and query params

```{r}
stub_request("get", "https://httpbin.org/get") %>%
  wi_th(query = list(hello = "world"), 
        headers = list('User-Agent' = 'libcurl/7.51.0 r-curl/2.6 crul/0.3.6', 
                       'Accept-Encoding' = "gzip, deflate"))
```

```{r}
stub_registry()
```

```{r}
x <- HttpClient$new(url = "https://httpbin.org")
x$get('get', query = list(hello = "world"))
```

### Stubbing requests and set expectation of a timeout

```{r error=TRUE}
stub_request("post", "https://httpbin.org/post") %>% to_timeout()
x <- HttpClient$new(url = "https://httpbin.org")
x$post('post')
```

### Stubbing requests and set HTTP error expectation

```{r error=TRUE}
library(fauxpas)
stub_request("get", "https://httpbin.org/get?a=b") %>% to_raise(HTTPBadRequest)
x <- HttpClient$new(url = "https://httpbin.org")
x$get('get', query = list(a = "b"))
```

## httr integration

```{r}
library(webmockr)
library(httr)

# turn on httr mocking
httr_mock()
```

```{r eval=FALSE}
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

```{r}
stub_request('get', uri = 'https://httpbin.org/get') %>%
  wi_th(
    headers = list('Accept' = 'application/json, text/xml, application/xml, */*')
  ) %>%
  to_return(status = 418, body = "I'm a teapot!!!", headers = list(im_a = "teapot"))
```

now returns mocked response


```{r eval=FALSE} 
(res <- GET("https://httpbin.org/get"))
res$status_code
#> [1] 418
res$headers
#> $im_a
#> [1] "teapot"
```

## Writing to disk

Write to a file before mocked request

```{r echo=FALSE}
stub_registry_clear()
request_registry_clear()
```

```{r}
## make a temp file
f <- tempfile(fileext = ".json")
## write something to the file
cat("{\"hello\":\"world\"}\n", file = f)
readLines(f)
## make the stub
invisible(stub_request("get", "https://httpbin.org/get") %>% 
  to_return(body = file(f)))
## make a request
out <- HttpClient$new("https://httpbin.org/get")$get(disk = f)
readLines(file(f))
```

OR - you can use `mock_file()` to have `webmockr` handle file and contents

```{r}
g <- tempfile(fileext = ".json")
## make the stub
invisible(stub_request("get", "https://httpbin.org/get") %>% 
  to_return(body = mock_file(g, "{\"hello\":\"mars\"}\n")))
## make a request
out <- crul::HttpClient$new("https://httpbin.org/get")$get(disk = g)
readLines(out$content)
```

Writing to disk is supported in both `crul` and `httr`

## Many requests in a row

e.g., many redirects, then a final successful request

```{r eval=FALSE}
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock_file.R
\name{mock_file}
\alias{mock_file}
\title{Mock file}
\usage{
mock_file(path, payload)
}
\arguments{
\item{path}{(character) a file path. required}

\item{payload}{(character) string to be written to the file given
at \code{path} parameter. required}
}
\value{
a list with S3 class \code{mock_file}
}
\description{
Mock file
}
\examples{
mock_file(path = tempfile(), payload = "{\"foo\": \"bar\"}")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{webmockr_enable}
\alias{webmockr_enable}
\title{This function is defunct.}
\usage{
webmockr_enable(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stub_request.R
\name{stub_request}
\alias{stub_request}
\title{Stub an http request}
\usage{
stub_request(method = "get", uri = NULL, uri_regex = NULL)
}
\arguments{
\item{method}{(character) HTTP method, one of "get", "post", "put", "patch",
"head", "delete", "options" - or the special "any" (for any method)}

\item{uri}{(character) The request uri. Can be a full or partial uri.
\pkg{webmockr} can match uri's without the "http" scheme, but does
not match if the scheme is "https". required, unless \code{uri_regex} given.
See \link{UriPattern} for more. See the "uri vs. uri_regex" section}

\item{uri_regex}{(character) A URI represented as regex. required, if \code{uri}
not given. See examples and the "uri vs. uri_regex" section}
}
\value{
an object of class \code{StubbedRequest}, with print method describing
the stub.
}
\description{
Stub an http request
}
\details{
Internally, this calls \link{StubbedRequest} which handles the logic

See \code{\link[=stub_registry]{stub_registry()}} for listing stubs, \code{\link[=stub_registry_clear]{stub_registry_clear()}}
for removing all stubs and \code{\link[=remove_request_stub]{remove_request_stub()}} for removing specific
stubs

If multiple stubs match the same request, we use the first stub. So if you
want to use a stub that was created after an earlier one that matches,
remove the earlier one(s).

Note on \code{wi_th()}: If you pass \code{query} values are coerced to character
class in the recorded stub. You can pass numeric, integer, etc., but
all will be coerced to character.

See \code{\link[=wi_th]{wi_th()}} for details on request body/query/headers and
\code{\link[=to_return]{to_return()}} for details on how response status/body/headers
are handled
}
\note{
Trailing slashes are dropped from stub URIs before matching
}
\section{uri vs. uri_regex}{

When you use \code{uri}, we compare the URIs without query params AND
also the query params themselves without the URIs.

When you use \code{uri_regex} we don't compare URIs and query params;
we just use your regex string defined in \code{uri_regex} as the pattern
for a call to \link{grepl}
}

\section{Mocking writing to disk}{

See \link{mocking-disk-writing}
}

\examples{
\dontrun{
# basic stubbing
stub_request("get", "https://httpbin.org/get")
stub_request("post", "https://httpbin.org/post")

# any method, use "any"
stub_request("any", "https://httpbin.org/get")

# list stubs
stub_registry()

# request headers
stub_request("get", "https://httpbin.org/get") \%>\%
   wi_th(headers = list('User-Agent' = 'R'))

# request body
stub_request("post", "https://httpbin.org/post") \%>\%
   wi_th(body = list(foo = 'bar'))
stub_registry()
library(crul)
x <- crul::HttpClient$new(url = "https://httpbin.org")
crul::mock()
x$post('post', body = list(foo = 'bar'))

# add expectation with to_return
stub_request("get", "https://httpbin.org/get") \%>\%
  wi_th(
    query = list(hello = "world"),
    headers = list('User-Agent' = 'R')) \%>\%
  to_return(status = 200, body = "stuff", headers = list(a = 5))

# list stubs again
stub_registry()

# regex
stub_request("get", uri_regex = ".+ample\\\\..")

# set stub an expectation to timeout
stub_request("get", "https://httpbin.org/get") \%>\% to_timeout()
x <- crul::HttpClient$new(url = "https://httpbin.org")
res <- x$get('get')

# raise exception
library(fauxpas)
stub_request("get", "https://httpbin.org/get") \%>\% to_raise(HTTPAccepted)
stub_request("get", "https://httpbin.org/get") \%>\% to_raise(HTTPAccepted, HTTPGone)

x <- crul::HttpClient$new(url = "https://httpbin.org")
stub_request("get", "https://httpbin.org/get") \%>\% to_raise(HTTPBadGateway)
crul::mock()
x$get('get')

# pass a list to .list
z <- stub_request("get", "https://httpbin.org/get")
wi_th(z, .list = list(query = list(foo = "bar")))

# just body
stub_request("any", uri_regex = ".+") \%>\%
   wi_th(body = list(foo = 'bar'))
## with crul
library(crul)
x <- crul::HttpClient$new(url = "https://httpbin.org")
crul::mock()
x$post('post', body = list(foo = 'bar'))
x$put('put', body = list(foo = 'bar'))
## with httr
library(httr)
httr_mock()
POST('https://example.com', body = list(foo = 'bar'))
PUT('https://google.com', body = list(foo = 'bar'))


# just headers
headers <- list(
  'Accept-Encoding' = 'gzip, deflate', 
  'Accept' = 'application/json, text/xml, application/xml, */*')
stub_request("any", uri_regex = ".+") \%>\% wi_th(headers = headers)
library(crul)
x <- crul::HttpClient$new(url = "https://httpbin.org", headers = headers)
crul::mock()
x$post('post')
x$put('put', body = list(foo = 'bar'))
x$get('put', query = list(stuff = 3423234L))

# many responses
## the first response matches the first to_return call, and so on
stub_request("get", "https://httpbin.org/get") \%>\% 
  to_return(status = 200, body = "foobar", headers = list(a = 5)) \%>\% 
  to_return(status = 200, body = "bears", headers = list(b = 6))
con <- crul::HttpClient$new(url = "https://httpbin.org")
con$get("get")$parse("UTF-8")
con$get("get")$parse("UTF-8")

## OR, use times with to_return() to repeat the same response many times
library(fauxpas)
stub_request("get", "https://httpbin.org/get") \%>\% 
  to_return(status = 200, body = "apple-pie", times = 2) \%>\% 
  to_raise(HTTPUnauthorized)
con <- crul::HttpClient$new(url = "https://httpbin.org")
con$get("get")$parse("UTF-8")
con$get("get")$parse("UTF-8")
con$get("get")$parse("UTF-8")

# clear all stubs
stub_registry()
stub_registry_clear()
}
}
\seealso{
\code{\link[=wi_th]{wi_th()}}, \code{\link[=to_return]{to_return()}}, \code{\link[=to_timeout]{to_timeout()}}, \code{\link[=to_raise]{to_raise()}},
\code{\link[=mock_file]{mock_file()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapter-httr.R
\name{build_httr_request}
\alias{build_httr_request}
\title{Build a httr request}
\usage{
build_httr_request(x)
}
\arguments{
\item{x}{an unexecuted httr request object}
}
\value{
a httr request
}
\description{
Build a httr request
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RequestPattern.R
\name{RequestPattern}
\alias{RequestPattern}
\title{RequestPattern class}
\description{
class handling all request matchers
}
\examples{
\dontrun{
(x <- RequestPattern$new(method = "get", uri = "httpbin.org/get"))
x$body_pattern
x$headers_pattern
x$method_pattern
x$uri_pattern
x$to_s()

# make a request signature
rs <- RequestSignature$new(method = "get", uri = "http://httpbin.org/get")

# check if it matches
x$matches(rs)

# regex uri
(x <- RequestPattern$new(method = "get", uri_regex = ".+ossref.org"))
x$uri_pattern
x$uri_pattern$to_s()
x$to_s()

# uri with query parameters
(x <- RequestPattern$new(
    method = "get", uri = "https://httpbin.org/get",
    query = list(foo = "bar")
))
x$to_s()
## query params included in url, not separately
(x <- RequestPattern$new(
  method = "get", uri = "https://httpbin.org/get?stuff=things"
))
x$to_s()
x$query_params

# just headers (via setting method=any & uri_regex=.+)
headers <- list(
  'User-Agent' = 'Apple',
  'Accept-Encoding' = 'gzip, deflate', 
  'Accept' = 'application/json, text/xml, application/xml, */*')
x <- RequestPattern$new(
   method = "any",
   uri_regex = ".+",
   headers = headers)
x$to_s()
rs <- RequestSignature$new(method = "any", uri = "http://foo.bar", 
  options = list(headers = headers))
rs
x$matches(rs)

# body
x <- RequestPattern$new(method = "post", uri = "httpbin.org/post",
  body = list(y = crul::upload(system.file("CITATION"))))
x$to_s()
rs <- RequestSignature$new(method = "post", uri = "http://httpbin.org/post",
  options = list(
     body = list(y = crul::upload(system.file("CITATION")))))
rs
x$matches(rs)
}
}
\seealso{
pattern classes for HTTP method \link{MethodPattern}, headers
\link{HeadersPattern}, body \link{BodyPattern}, and URI/URL \link{UriPattern}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{method_pattern}}{xxx}

\item{\code{uri_pattern}}{xxx}

\item{\code{body_pattern}}{xxx}

\item{\code{headers_pattern}}{xxx}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{RequestPattern$new()}}
\item \href{#method-matches}{\code{RequestPattern$matches()}}
\item \href{#method-to_s}{\code{RequestPattern$to_s()}}
\item \href{#method-clone}{\code{RequestPattern$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{RequestPattern} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestPattern$new(
  method,
  uri = NULL,
  uri_regex = NULL,
  query = NULL,
  body = NULL,
  headers = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{the HTTP method (any, head, options, get, post, put,
patch, trace, or delete). "any" matches any HTTP method. required.}

\item{\code{uri}}{(character) request URI. required or uri_regex}

\item{\code{uri_regex}}{(character) request URI as regex. required or uri}

\item{\code{query}}{(list) query parameters, optional}

\item{\code{body}}{(list) body request, optional}

\item{\code{headers}}{(list) headers, optional}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{RequestPattern} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-matches"></a>}}
\if{latex}{\out{\hypertarget{method-matches}{}}}
\subsection{Method \code{matches()}}{
does a request signature match the selected matchers?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestPattern$matches(request_signature)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request_signature}}{a \link{RequestSignature} object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
Print pattern for easy human consumption
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestPattern$to_s()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestPattern$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/zzz.R
\name{webmockr_crul_fetch}
\alias{webmockr_crul_fetch}
\title{execute a curl request}
\usage{
webmockr_crul_fetch(x)
}
\arguments{
\item{x}{an object}
}
\value{
a curl response
}
\description{
execute a curl request
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapter-httr.R
\name{httr_mock}
\alias{httr_mock}
\title{Turn on httr mocking
Sets a callback that routes httr request through webmockr}
\usage{
httr_mock(on = TRUE)
}
\arguments{
\item{on}{(logical) set to \code{TRUE} to turn on, and \code{FALSE}
to turn off. default: \code{TRUE}}
}
\value{
Silently returns \code{TRUE} when enabled and \code{FALSE} when disabled.
}
\description{
Turn on httr mocking
Sets a callback that routes httr request through webmockr
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{webmockr_disable}
\alias{webmockr_disable}
\title{This function is defunct.}
\usage{
webmockr_disable(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RequestPattern.R
\name{MethodPattern}
\alias{MethodPattern}
\title{MethodPattern}
\description{
method matcher
}
\details{
Matches regardless of case. e.g., POST will match to post
}
\examples{
(x <- MethodPattern$new(pattern = "post"))
x$pattern
x$matches(method = "post")
x$matches(method = "POST")

# all matches() calls should be TRUE
(x <- MethodPattern$new(pattern = "any"))
x$pattern
x$matches(method = "post")
x$matches(method = "GET")
x$matches(method = "HEAD")
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{pattern}}{(character) an http method}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{MethodPattern$new()}}
\item \href{#method-matches}{\code{MethodPattern$matches()}}
\item \href{#method-to_s}{\code{MethodPattern$to_s()}}
\item \href{#method-clone}{\code{MethodPattern$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{MethodPattern} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MethodPattern$new(pattern)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{pattern}}{(character) a HTTP method, lowercase}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{MethodPattern} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-matches"></a>}}
\if{latex}{\out{\hypertarget{method-matches}{}}}
\subsection{Method \code{matches()}}{
test if the pattern matches a given http method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MethodPattern$matches(method)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{(character) a HTTP method, lowercase}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
Print pattern for easy human consumption
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MethodPattern$to_s()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{MethodPattern$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/to_timeout.R
\name{to_timeout}
\alias{to_timeout}
\title{Set timeout as an expected return on a match}
\usage{
to_timeout(.data)
}
\arguments{
\item{.data}{input. Anything that can be coerced to a \code{StubbedRequest} class
object}
}
\value{
an object of class \code{StubbedRequest}, with print method describing
the stub
}
\description{
Set timeout as an expected return on a match
}
\note{
see examples in \code{\link[=stub_request]{stub_request()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stub_registry_clear.R
\name{stub_registry_clear}
\alias{stub_registry_clear}
\title{stub_registry_clear}
\usage{
stub_registry_clear()
}
\value{
an empty list invisibly
}
\description{
Clear all stubs in the stub registry
}
\examples{
(x <- stub_request("get", "https://httpbin.org/get"))
stub_registry()
stub_registry_clear()
stub_registry()
}
\seealso{
Other stub-registry: 
\code{\link{StubRegistry}},
\code{\link{remove_request_stub}()},
\code{\link{stub_registry}()}
}
\concept{stub-registry}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RequestRegistry.R
\name{RequestRegistry}
\alias{RequestRegistry}
\title{RequestRegistry}
\description{
keeps track of HTTP requests
}
\examples{
x <- RequestRegistry$new()
z1 <- RequestSignature$new("get", "http://scottchamberlain.info")
z2 <- RequestSignature$new("post", "https://httpbin.org/post")
x$register_request(request = z1)
x$register_request(request = z1)
x$register_request(request = z2)
# print method to list requests
x

# more complex requests
w <- RequestSignature$new(
  method = "get",
  uri = "https:/httpbin.org/get",
  options = list(headers = list(`User-Agent` = "foobar", stuff = "things"))
)
w$to_s()
x$register_request(request = w)
x


# hashes, and number of times each requested
x$request_signatures$hash

# times_executed method
pat <- RequestPattern$new(
  method = "get",
  uri = "https:/httpbin.org/get",
  headers = list(`User-Agent` = "foobar", stuff = "things")
)
pat$to_s()
x$times_executed(pat)
z <- RequestPattern$new(method = "get", uri = "http://scottchamberlain.info")
x$times_executed(z)
w <- RequestPattern$new(method = "post", uri = "https://httpbin.org/post")
x$times_executed(w)

## pattern with no matches - returns 0 (zero)
pat <- RequestPattern$new(
  method = "get",
  uri = "http://recology.info/"
)
pat$to_s()
x$times_executed(pat)

# reset the request registry
x$reset()
}
\seealso{
\code{\link[=stub_registry]{stub_registry()}} and \link{StubRegistry}

Other request-registry: 
\code{\link{HashCounter}},
\code{\link{request_registry}()}
}
\concept{request-registry}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{request_signatures}}{a HashCounter object}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{RequestRegistry$print()}}
\item \href{#method-reset}{\code{RequestRegistry$reset()}}
\item \href{#method-register_request}{\code{RequestRegistry$register_request()}}
\item \href{#method-times_executed}{\code{RequestRegistry$times_executed()}}
\item \href{#method-clone}{\code{RequestRegistry$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{RequestRegistry} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestRegistry$print(x, ...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-reset"></a>}}
\if{latex}{\out{\hypertarget{method-reset}{}}}
\subsection{Method \code{reset()}}{
Reset the registry to no registered requests
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestRegistry$reset()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned; ressets registry to no requests
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-register_request"></a>}}
\if{latex}{\out{\hypertarget{method-register_request}{}}}
\subsection{Method \code{register_request()}}{
Register a request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestRegistry$register_request(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{a character string of the request, serialized from
a \code{RequestSignature$new(...)$to_s()}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; registers the request
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-times_executed"></a>}}
\if{latex}{\out{\hypertarget{method-times_executed}{}}}
\subsection{Method \code{times_executed()}}{
How many times has a request been made
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestRegistry$times_executed(request_pattern)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request_pattern}}{an object of class \code{RequestPattern}}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
if no match is found for the request pattern, 0 is returned
}

\subsection{Returns}{
integer, the number of times the request has been made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestRegistry$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/webmockr-opts.R
\name{webmockr_configure}
\alias{webmockr_configure}
\alias{webmockr_configure_reset}
\alias{webmockr_configuration}
\alias{webmockr_allow_net_connect}
\alias{webmockr_disable_net_connect}
\alias{webmockr_net_connect_allowed}
\title{webmockr configuration}
\usage{
webmockr_configure(
  allow_net_connect = FALSE,
  allow_localhost = FALSE,
  allow = NULL,
  net_http_connect_on_start = FALSE,
  show_stubbing_instructions = FALSE,
  query_values_notation = FALSE,
  show_body_diff = FALSE
)

webmockr_configure_reset()

webmockr_configuration()

webmockr_allow_net_connect()

webmockr_disable_net_connect(allow = NULL)

webmockr_net_connect_allowed(uri = NULL)
}
\arguments{
\item{allow_net_connect}{(logical) Default: \code{FALSE}}

\item{allow_localhost}{(logical) Default: \code{FALSE}}

\item{allow}{(character) one or more URI/URL to allow (and by extension
all others are not allowed)}

\item{net_http_connect_on_start}{(logical) Default: \code{FALSE}. ignored for
now}

\item{show_stubbing_instructions}{(logical) Default: \code{FALSE}. ignored for
now}

\item{query_values_notation}{(logical) Default: \code{FALSE}. ignored for
now}

\item{show_body_diff}{(logical) Default: \code{FALSE}. ignored for
now}

\item{uri}{(character) a URI/URL as a character string - to determine
whether or not it is allowed}
}
\description{
webmockr configuration
}
\section{webmockr_allow_net_connect}{

If there are stubs found for a request, even if net connections are
allowed (by running \code{webmockr_allow_net_connect()}) the stubbed
response will be returned. If no stub is found, and net connections
are allowed, then a real HTTP request can be made.
}

\examples{
\dontrun{
webmockr_configure()
webmockr_configure(
 allow_localhost = TRUE
)
webmockr_configuration()
webmockr_configure_reset()

webmockr_allow_net_connect()
webmockr_net_connect_allowed()

# disable net connect for any URIs
webmockr_disable_net_connect()
### gives NULL with no URI passed
webmockr_net_connect_allowed()
# disable net connect EXCEPT FOR given URIs
webmockr_disable_net_connect(allow = "google.com")
### is a specific URI allowed?
webmockr_net_connect_allowed("google.com")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_request_stub.R
\name{remove_request_stub}
\alias{remove_request_stub}
\title{Remove a request stub}
\usage{
remove_request_stub(stub)
}
\arguments{
\item{stub}{a request stub, of class \code{StubbedRequest}}
}
\value{
logical, \code{TRUE} if removed, \code{FALSE} if not removed
}
\description{
Remove a request stub
}
\examples{
(x <- stub_request("get", "https://httpbin.org/get"))
stub_registry()
remove_request_stub(x)
stub_registry()
}
\seealso{
Other stub-registry: 
\code{\link{StubRegistry}},
\code{\link{stub_registry_clear}()},
\code{\link{stub_registry}()}
}
\concept{stub-registry}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wi_th.R
\name{wi_th}
\alias{wi_th}
\title{Set additional parts of a stubbed request}
\usage{
wi_th(.data, ..., .list = list())
}
\arguments{
\item{.data}{input. Anything that can be coerced to a \code{StubbedRequest} class
object}

\item{...}{Comma separated list of named variables. accepts the following:
\code{query}, \code{body}, \code{headers}, \code{basic_auth}. See Details.}

\item{.list}{named list, has to be one of \code{query}, \code{body},
\code{headers} and/or \code{basic_auth}. An alternative to passing in via \code{...}.
Don't pass the same thing to both, e.g. don't pass 'query' to \code{...}, and
also 'query' to this parameter}
}
\value{
an object of class \code{StubbedRequest}, with print method describing
the stub
}
\description{
Set query params, request body, request headers and/or basic_auth
}
\details{
\code{with} is a function in the \code{base} package, so we went with
\code{wi_th}

Values for query, body, headers, and basic_auth:
\itemize{
\item query: (list) a named list. values are coerced to character
class in the recorded stub. You can pass numeric, integer, etc., but
all will be coerced to character.
\item body: various, including character string, list, raw, numeric,
upload (\code{crul::upload} or \code{httr::upload_file}, they both create the
same object in the end)
\item headers: (list) a named list
\item basic_auth: (character) a length two vector, username and password.
authentication type (basic/digest/ntlm/etc.) is ignored. that is,
mocking authenciation right now does not take into account the
authentication type. We don't do any checking of the username/password
except to detect edge cases where for example, the username/password
were probably not set by the user on purpose (e.g., a URL is
picked up by an environment variable)
}

Note that there is no regex matching on query, body, or headers. They
are tested for matches in the following ways:
\itemize{
\item query: compare stubs and requests with \code{identical()}. this compares
named lists, so both list names and values are compared
\item body: varies depending on the body format (list vs. character, etc.)
\item headers: compare stub and request values with \code{==}. list names are
compared with \code{\%in\%}. \code{basic_auth} is included in headers (with the name
Authorization)
}
}
\note{
see more examples in \code{\link[=stub_request]{stub_request()}}
}
\examples{
# first, make a stub object
req <- stub_request("post", "https://httpbin.org/post")

# add body
# list
wi_th(req, body = list(foo = "bar"))
# string
wi_th(req, body = '{"foo": "bar"}')
# raw
wi_th(req, body = charToRaw('{"foo": "bar"}'))
# numeric
wi_th(req, body = 5)
# an upload
wi_th(req, body = crul::upload(system.file("CITATION")))
# wi_th(req, body = httr::upload_file(system.file("CITATION")))

# add query - has to be a named list
wi_th(req, query = list(foo = "bar"))

# add headers - has to be a named list
wi_th(req, headers = list(foo = "bar"))
wi_th(req, headers = list(`User-Agent` = "webmockr/v1", hello="world"))

# .list - pass in a named list instead
wi_th(req, .list = list(body = list(foo = "bar")))

# basic authentication
wi_th(req, basic_auth = c("user", "pass"))
wi_th(req, basic_auth = c("user", "pass"), headers = list(foo = "bar"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RequestPattern.R
\name{BodyPattern}
\alias{BodyPattern}
\title{BodyPattern}
\description{
body matcher
}
\examples{
# make a request signature
bb <- RequestSignature$new(
  method = "get",
  uri = "https:/httpbin.org/get",
  options = list(
    body = list(foo = "bar", a = 5)
  )
)

# make body pattern object
## FALSE
z <- BodyPattern$new(pattern = list(foo = "bar"))
z$pattern
z$matches(bb$body)
## TRUE
z <- BodyPattern$new(pattern = list(foo = "bar", a = 5))
z$pattern
z$matches(bb$body)

# uploads in bodies
## upload NOT in a list
bb <- RequestSignature$new(
  method = "post", uri = "https:/httpbin.org/post",
  options = list(body = crul::upload(system.file("CITATION"))))
bb$body
z <- BodyPattern$new(pattern = 
  crul::upload(system.file("CITATION")))
z$pattern
z$matches(bb$body)

## upload in a list
bb <- RequestSignature$new(
  method = "post", uri = "https:/httpbin.org/post",
  options = list(body = list(y = crul::upload(system.file("CITATION")))))
bb$body
z <- BodyPattern$new(pattern =
  list(y = crul::upload(system.file("CITATION"))))
z$pattern
z$matches(bb$body)
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{pattern}}{a list}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{BodyPattern$new()}}
\item \href{#method-matches}{\code{BodyPattern$matches()}}
\item \href{#method-to_s}{\code{BodyPattern$to_s()}}
\item \href{#method-clone}{\code{BodyPattern$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{BodyPattern} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BodyPattern$new(pattern)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{pattern}}{(list) a body object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{BodyPattern} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-matches"></a>}}
\if{latex}{\out{\hypertarget{method-matches}{}}}
\subsection{Method \code{matches()}}{
Match a list of headers against that stored
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BodyPattern$matches(body, content_type = "")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{body}}{(list) the body}

\item{\code{content_type}}{(character) content type}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
Print pattern for easy human consumption
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BodyPattern$to_s()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{BodyPattern$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/webmockr.R
\docType{package}
\name{webmockr-package}
\alias{webmockr-package}
\alias{webmockr}
\title{webmockr}
\description{
Stubbing and setting expectations on HTTP requests
}
\section{Features}{

\itemize{
\item Stubbing HTTP requests at low http client lib level
\item Setting and verifying expectations on HTTP requests
\item Matching requests based on method, URI, headers and body
\item Supports multiple HTTP libraries, including \pkg{crul} and
\pkg{httr}
\item Integration with HTTP test caching library \pkg{vcr}
}
}

\examples{
library(webmockr)
stub_request("get", "https://httpbin.org/get")
stub_request("post", "https://httpbin.org/post")
stub_registry()
}
\author{
Scott Chamberlain \email{myrmecocystus+r@gmail.com}

Aaron Wolen
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RequestRegistry.R
\name{HashCounter}
\alias{HashCounter}
\title{HashCounter}
\description{
hash with counter, to store requests, and count each time
it is used
}
\examples{
x <- HashCounter$new()
x$hash
z <- RequestSignature$new(method = "get", uri = "https:/httpbin.org/get")
x$put(z)
x$hash
x$get(z)
x$put(z)
x$get(z)
}
\seealso{
Other request-registry: 
\code{\link{RequestRegistry}},
\code{\link{request_registry}()}
}
\concept{request-registry}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{hash}}{(list) a list for internal use only, with elements
\code{key}, \code{sig}, and \code{count}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-put}{\code{HashCounter$put()}}
\item \href{#method-get}{\code{HashCounter$get()}}
\item \href{#method-clone}{\code{HashCounter$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-put"></a>}}
\if{latex}{\out{\hypertarget{method-put}{}}}
\subsection{Method \code{put()}}{
Register a request by it's key
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HashCounter$put(req_sig)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{req_sig}}{an object of class \code{RequestSignature}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; registers request and iterates
internal counter
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get"></a>}}
\if{latex}{\out{\hypertarget{method-get}{}}}
\subsection{Method \code{get()}}{
Get a request by key
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HashCounter$get(req_sig)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{req_sig}}{an object of class \code{RequestSignature}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
(integer) the count of how many times the request has been made
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HashCounter$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/RequestPattern.R
\name{HeadersPattern}
\alias{HeadersPattern}
\title{HeadersPattern}
\description{
headers matcher
}
\details{
\code{webmockr} normalises headers and treats all forms of same headers as equal:
i.e the following two sets of headers are equal:
\code{list(Header1 = "value1", content_length = 123, X_CuStOm_hEAder = "foo")}
and
\code{list(header1 = "value1", "Content-Length" = 123, "x-cuSTOM-HeAder" = "foo")}
}
\examples{
(x <- HeadersPattern$new(pattern = list(a = 5)))
x$pattern
x$matches(list(a = 5))

# different cases
(x <- HeadersPattern$new(pattern = list(Header1 = "value1")))
x$pattern
x$matches(list(header1 = "value1"))
x$matches(list(header1 = "value2"))

# different symbols
(x <- HeadersPattern$new(pattern = list(`Hello_World` = "yep")))
x$pattern
x$matches(list(`hello-world` = "yep"))
x$matches(list(`hello-worlds` = "yep"))

headers <- list(
  'User-Agent' = 'Apple',
  'Accept-Encoding' = 'gzip, deflate', 
  'Accept' = 'application/json, text/xml, application/xml, */*')
(x <- HeadersPattern$new(pattern = headers))
x$to_s()
x$pattern
x$matches(headers)
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{pattern}}{a list}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{HeadersPattern$new()}}
\item \href{#method-matches}{\code{HeadersPattern$matches()}}
\item \href{#method-empty_headers}{\code{HeadersPattern$empty_headers()}}
\item \href{#method-to_s}{\code{HeadersPattern$to_s()}}
\item \href{#method-clone}{\code{HeadersPattern$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{HeadersPattern} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HeadersPattern$new(pattern)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{pattern}}{(list) a pattern, as a named list, must be named,
e.g,. \code{list(a = 5, b = 6)}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{HeadersPattern} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-matches"></a>}}
\if{latex}{\out{\hypertarget{method-matches}{}}}
\subsection{Method \code{matches()}}{
Match a list of headers against that stored
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HeadersPattern$matches(headers)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{headers}}{(list) named list of headers, e.g,. \code{list(a = 5, b = 6)}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-empty_headers"></a>}}
\if{latex}{\out{\hypertarget{method-empty_headers}{}}}
\subsection{Method \code{empty_headers()}}{
Are headers empty? tests if null or length==0
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HeadersPattern$empty_headers(headers)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{headers}}{named list of headers}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
Print pattern for easy human consumption
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HeadersPattern$to_s()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HeadersPattern$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/HttpLibAdapterRegistry.R
\name{HttpLibAdapaterRegistry}
\alias{HttpLibAdapaterRegistry}
\title{HttpLibAdapaterRegistry}
\description{
http lib adapter registry
}
\examples{
x <- HttpLibAdapaterRegistry$new()
x$register(CrulAdapter$new())
x
x$adapters
x$adapters[[1]]$name
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{adapters}}{list}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{HttpLibAdapaterRegistry$print()}}
\item \href{#method-register}{\code{HttpLibAdapaterRegistry$register()}}
\item \href{#method-clone}{\code{HttpLibAdapaterRegistry$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{HttpLibAdapaterRegistry} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpLibAdapaterRegistry$print(x, ...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-register"></a>}}
\if{latex}{\out{\hypertarget{method-register}{}}}
\subsection{Method \code{register()}}{
Register an http library adapter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpLibAdapaterRegistry$register(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{an http lib adapter, e.g., \link{CrulAdapter}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing, registers the library adapter
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttpLibAdapaterRegistry$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/defunct.R
\name{to_return_}
\alias{to_return_}
\title{This function is defunct.}
\usage{
to_return_(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pluck_body.R
\name{pluck_body}
\alias{pluck_body}
\title{Extract the body from an HTTP request}
\usage{
pluck_body(x)
}
\arguments{
\item{x}{an unexecuted crul \emph{or} httr request object}
}
\value{
one of the following:
\itemize{
\item \code{NULL} if the request is not associated with a body
\item \code{NULL} if an upload is used not in a list
\item list containing the multipart-encoded body
\item character vector with the JSON- or raw-encoded body, or upload form file
}
}
\description{
Returns an appropriate representation of the data contained within a request
body based on its encoding.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapter-crul.R
\name{build_crul_request}
\alias{build_crul_request}
\title{Build a crul request}
\usage{
build_crul_request(x)
}
\arguments{
\item{x}{an unexecuted crul request object}
}
\value{
a crul request
}
\description{
Build a crul request
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RequestSignature.R
\name{RequestSignature}
\alias{RequestSignature}
\title{RequestSignature}
\description{
General purpose request signature builder
}
\examples{
# make request signature
x <- RequestSignature$new(method = "get", uri = "https:/httpbin.org/get")
# method
x$method
# uri
x$uri
# request signature to string
x$to_s()

# headers
w <- RequestSignature$new(
  method = "get",
  uri = "https:/httpbin.org/get",
  options = list(headers = list(`User-Agent` = "foobar", stuff = "things"))
)
w
w$headers
w$to_s()

# headers and body
bb <- RequestSignature$new(
  method = "get",
  uri = "https:/httpbin.org/get",
  options = list(
    headers = list(`User-Agent` = "foobar", stuff = "things"),
    body = list(a = "tables")
  )
)
bb
bb$headers
bb$body
bb$to_s()

# with disk path
f <- tempfile()
bb <- RequestSignature$new(
  method = "get",
  uri = "https:/httpbin.org/get",
  options = list(disk = f)
)
bb
bb$disk
bb$to_s()
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{method}}{(character) an http method}

\item{\code{uri}}{(character) a uri}

\item{\code{body}}{(various) request body}

\item{\code{headers}}{(list) named list of headers}

\item{\code{proxies}}{(list) proxies as a named list}

\item{\code{auth}}{(list) authentication details, as a named list}

\item{\code{url}}{internal use}

\item{\code{disk}}{(character) if writing to disk, the path}

\item{\code{fields}}{(various) request body details}

\item{\code{output}}{(various) request output details, disk, memory, etc}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{RequestSignature$new()}}
\item \href{#method-print}{\code{RequestSignature$print()}}
\item \href{#method-to_s}{\code{RequestSignature$to_s()}}
\item \href{#method-clone}{\code{RequestSignature$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{RequestSignature} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestSignature$new(method, uri, options = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{the HTTP method (any, head, options, get, post, put,
patch, trace, or delete). "any" matches any HTTP method. required.}

\item{\code{uri}}{(character) request URI. required.}

\item{\code{options}}{(list) options. optional. See Details.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{RequestSignature} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{RequestSignature} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestSignature$print()}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
Request signature to a string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestSignature$to_s()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a character string representation of the request signature
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestSignature$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/adapter-crul.R, R/adapter-httr.R, R/adapter.R
\name{CrulAdapter}
\alias{CrulAdapter}
\alias{HttrAdapter}
\alias{Adapter}
\title{Adapters for Modifying HTTP Requests}
\description{
\code{Adapter} is the base parent class used to implement
\pkg{webmockr} support for different HTTP clients. It should not be used
directly. Instead, use one of the client-specific adapters that webmockr
currently provides:
\itemize{
\item \code{CrulAdapter} for \pkg{crul}
\item \code{HttrAdapter} for \pkg{httr}
}
}
\details{
Note that the documented fields and methods are the same across all
client-specific adapters.
}
\examples{
\dontrun{
if (requireNamespace("httr", quietly = TRUE)) {
# library(httr)

# normal httr request, works fine
# real <- GET("https://httpbin.org/get")
# real

# with webmockr
# library(webmockr)
## turn on httr mocking
# httr_mock()
## now this request isn't allowed
# GET("https://httpbin.org/get")
## stub the request
# stub_request('get', uri = 'https://httpbin.org/get') \%>\%
#   wi_th(
#     headers = list('Accept' = 'application/json, text/xml, application/xml, */*')
#   ) \%>\%
#   to_return(status = 418, body = "I'm a teapot!", headers = list(a = 5))
## now the request succeeds and returns a mocked response
# (res <- GET("https://httpbin.org/get"))
# res$status_code
# rawToChar(res$content)

# allow real requests while webmockr is loaded
# webmockr_allow_net_connect()
# webmockr_net_connect_allowed()
# GET("https://httpbin.org/get?animal=chicken")
# webmockr_disable_net_connect()
# webmockr_net_connect_allowed()
# GET("https://httpbin.org/get?animal=chicken")

# httr_mock(FALSE)
}
}
}
\section{Super class}{
\code{\link[webmockr:Adapter]{webmockr::Adapter}} -> \code{CrulAdapter}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{client}}{HTTP client package name}

\item{\code{name}}{adapter name}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clone}{\code{CrulAdapter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="disable">}\href{../../webmockr/html/Adapter.html#method-disable}{\code{webmockr::Adapter$disable()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="enable">}\href{../../webmockr/html/Adapter.html#method-enable}{\code{webmockr::Adapter$enable()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="handle_request">}\href{../../webmockr/html/Adapter.html#method-handle_request}{\code{webmockr::Adapter$handle_request()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="initialize">}\href{../../webmockr/html/Adapter.html#method-initialize}{\code{webmockr::Adapter$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="remove_stubs">}\href{../../webmockr/html/Adapter.html#method-remove_stubs}{\code{webmockr::Adapter$remove_stubs()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CrulAdapter$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
\section{Super class}{
\code{\link[webmockr:Adapter]{webmockr::Adapter}} -> \code{HttrAdapter}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{client}}{HTTP client package name}

\item{\code{name}}{adapter name}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clone}{\code{HttrAdapter$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="disable">}\href{../../webmockr/html/Adapter.html#method-disable}{\code{webmockr::Adapter$disable()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="enable">}\href{../../webmockr/html/Adapter.html#method-enable}{\code{webmockr::Adapter$enable()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="handle_request">}\href{../../webmockr/html/Adapter.html#method-handle_request}{\code{webmockr::Adapter$handle_request()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="initialize">}\href{../../webmockr/html/Adapter.html#method-initialize}{\code{webmockr::Adapter$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="webmockr" data-topic="Adapter" data-id="remove_stubs">}\href{../../webmockr/html/Adapter.html#method-remove_stubs}{\code{webmockr::Adapter$remove_stubs()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HttrAdapter$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{client}}{HTTP client package name}

\item{\code{name}}{adapter name}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Adapter$new()}}
\item \href{#method-enable}{\code{Adapter$enable()}}
\item \href{#method-disable}{\code{Adapter$disable()}}
\item \href{#method-handle_request}{\code{Adapter$handle_request()}}
\item \href{#method-remove_stubs}{\code{Adapter$remove_stubs()}}
\item \href{#method-clone}{\code{Adapter$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Adapter object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Adapter$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-enable"></a>}}
\if{latex}{\out{\hypertarget{method-enable}{}}}
\subsection{Method \code{enable()}}{
Enable the adapter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Adapter$enable(quiet = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{quiet}}{(logical) suppress messages? default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{TRUE}, invisibly
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-disable"></a>}}
\if{latex}{\out{\hypertarget{method-disable}{}}}
\subsection{Method \code{disable()}}{
Disable the adapter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Adapter$disable(quiet = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{quiet}}{(logical) suppress messages? default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{FALSE}, invisibly
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-handle_request"></a>}}
\if{latex}{\out{\hypertarget{method-handle_request}{}}}
\subsection{Method \code{handle_request()}}{
All logic for handling a request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Adapter$handle_request(req)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{req}}{a request}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
various outcomes
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remove_stubs"></a>}}
\if{latex}{\out{\hypertarget{method-remove_stubs}{}}}
\subsection{Method \code{remove_stubs()}}{
Remove all stubs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Adapter$remove_stubs()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned; removes all request stubs
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Adapter$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/webmockr_reset.R
\name{webmockr_reset}
\alias{webmockr_reset}
\title{webmockr_reset}
\usage{
webmockr_reset()
}
\value{
nothing
}
\description{
Clear all stubs and the request counter
}
\details{
this function runs \code{\link[=stub_registry_clear]{stub_registry_clear()}} and
\code{\link[=request_registry_clear]{request_registry_clear()}} - so you can run those two yourself
to achieve the same thing
}
\examples{
# webmockr_reset()
}
\seealso{
\code{\link[=stub_registry_clear]{stub_registry_clear()}} \code{\link[=request_registry_clear]{request_registry_clear()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{webmockr-defunct}
\alias{webmockr-defunct}
\title{Defunct functions in \pkg{webmockr}}
\description{
\itemize{
\item \code{\link[=webmockr_enable]{webmockr_enable()}}: Function removed, see \code{\link[=enable]{enable()}}
\item \code{\link[=webmockr_disable]{webmockr_disable()}}: Function removed, see \code{\link[=disable]{disable()}}
\item \link{to_return_}: Only \code{\link[=to_return]{to_return()}} is available now
\item \link{wi_th_}: Only \code{\link[=wi_th]{wi_th()}} is available now
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
Pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapter-httr.R
\name{build_httr_response}
\alias{build_httr_response}
\title{Build a httr response}
\usage{
build_httr_response(req, resp)
}
\arguments{
\item{req}{a request}

\item{resp}{a response}
}
\value{
a httr response
}
\description{
Build a httr response
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_raise.R
\name{to_raise}
\alias{to_raise}
\title{Set raise error condition}
\usage{
to_raise(.data, ...)
}
\arguments{
\item{.data}{input. Anything that can be coerced to a \code{StubbedRequest}
class object}

\item{...}{One or more HTTP exceptions from the \pkg{fauxpas} package. Run
\code{grep("HTTP*", getNamespaceExports("fauxpas"), value = TRUE)} for a list of
possible exceptions}
}
\value{
an object of class \code{StubbedRequest}, with print method describing
the stub
}
\description{
Set raise error condition
}
\details{
The behavior in the future will be:

When multiple exceptions are passed, the first is used on the first
mock, the second on the second mock, and so on. Subsequent mocks use the
last exception

But for now, only the first exception is used until we get that fixed
}
\note{
see examples in \code{\link[=stub_request]{stub_request()}}
}
\section{Raise vs. Return}{

\code{to_raise()} always raises a stop condition, while \code{to_return(status=xyz)} only
sets the status code on the returned HTTP response object. So if you want to
raise a stop condition then \code{to_raise()} is what you want. But if you
don't want to raise a stop condition use \code{to_return()}. Use cases for each
vary. For example, in a unit test you may have a test expecting a 503 error;
in this case \code{to_raise()} makes sense. In another case, if a unit test
expects to test some aspect of an HTTP response object that httr or crul
typically returns, then you'll want \code{to_return()}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_registry.R
\name{request_registry}
\alias{request_registry}
\alias{request_registry_clear}
\title{List or clear requests in the request registry}
\usage{
request_registry()

request_registry_clear()
}
\value{
an object of class \code{RequestRegistry}, print method gives the
requests in the registry and the number of times each one has been
performed
}
\description{
List or clear requests in the request registry
}
\details{
\code{request_registry()} lists the requests that have been made
that webmockr knows about; \code{request_registry_clear()} resets the
request registry (removes all recorded requests)
}
\examples{
webmockr::enable()
stub_request("get", "https://httpbin.org/get") \%>\%
  to_return(body = "success!", status = 200)

# nothing in the request registry
request_registry()

# make the request
z <- crul::HttpClient$new(url = "https://httpbin.org")$get("get")

# check the request registry - the request was made 1 time
request_registry()

# do the request again
z <- crul::HttpClient$new(url = "https://httpbin.org")$get("get")

# check the request registry - now it's been made 2 times, yay!
request_registry()

# clear the request registry
request_registry_clear()
webmockr::disable()
}
\seealso{
Other request-registry: 
\code{\link{HashCounter}},
\code{\link{RequestRegistry}}
}
\concept{request-registry}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapter-crul.R
\name{build_crul_response}
\alias{build_crul_response}
\title{Build a crul response}
\usage{
build_crul_response(req, resp)
}
\arguments{
\item{req}{a request}

\item{resp}{a response}
}
\value{
a crul response
}
\description{
Build a crul response
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StubRegistry.R
\name{StubRegistry}
\alias{StubRegistry}
\title{StubRegistry}
\description{
stub registry to keep track of \link{StubbedRequest} stubs
}
\examples{
\dontrun{
# Make a stub
stub1 <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
stub1$with(headers = list('User-Agent' = 'R'))
stub1$to_return(status = 200, body = "foobar", headers = list())
stub1

# Make another stub
stub2 <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
stub2

# Put both stubs in the stub registry
reg <- StubRegistry$new()
reg$register_stub(stub = stub1)
reg$register_stub(stub = stub2)
reg
reg$request_stubs
}
}
\seealso{
Other stub-registry: 
\code{\link{remove_request_stub}()},
\code{\link{stub_registry_clear}()},
\code{\link{stub_registry}()}
}
\concept{stub-registry}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{request_stubs}}{(list) list of request stubs}

\item{\code{global_stubs}}{(list) list of global stubs}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{StubRegistry$print()}}
\item \href{#method-register_stub}{\code{StubRegistry$register_stub()}}
\item \href{#method-find_stubbed_request}{\code{StubRegistry$find_stubbed_request()}}
\item \href{#method-request_stub_for}{\code{StubRegistry$request_stub_for()}}
\item \href{#method-remove_request_stub}{\code{StubRegistry$remove_request_stub()}}
\item \href{#method-remove_all_request_stubs}{\code{StubRegistry$remove_all_request_stubs()}}
\item \href{#method-is_registered}{\code{StubRegistry$is_registered()}}
\item \href{#method-clone}{\code{StubRegistry$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{StubRegistry} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$print(x, ...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-register_stub"></a>}}
\if{latex}{\out{\hypertarget{method-register_stub}{}}}
\subsection{Method \code{register_stub()}}{
Register a stub
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$register_stub(stub)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{stub}}{an object of type \link{StubbedRequest}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; registers the stub
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-find_stubbed_request"></a>}}
\if{latex}{\out{\hypertarget{method-find_stubbed_request}{}}}
\subsection{Method \code{find_stubbed_request()}}{
Find a stubbed request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$find_stubbed_request(req)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{req}}{an object of class \link{RequestSignature}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
an object of type \link{StubbedRequest}, if matched
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-request_stub_for"></a>}}
\if{latex}{\out{\hypertarget{method-request_stub_for}{}}}
\subsection{Method \code{request_stub_for()}}{
Find a stubbed request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$request_stub_for(request_signature)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request_signature}}{an object of class \link{RequestSignature}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
logical, 1 or more
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remove_request_stub"></a>}}
\if{latex}{\out{\hypertarget{method-remove_request_stub}{}}}
\subsection{Method \code{remove_request_stub()}}{
Remove a stubbed request by matching request signature
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$remove_request_stub(stub)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{stub}}{an object of type \link{StubbedRequest}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; removes the stub from the registry
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remove_all_request_stubs"></a>}}
\if{latex}{\out{\hypertarget{method-remove_all_request_stubs}{}}}
\subsection{Method \code{remove_all_request_stubs()}}{
Remove all request stubs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$remove_all_request_stubs()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned; removes all request stubs
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_registered"></a>}}
\if{latex}{\out{\hypertarget{method-is_registered}{}}}
\subsection{Method \code{is_registered()}}{
Find a stubbed request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$is_registered(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{an object of class \link{RequestSignature}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; registers the stub
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubRegistry$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/stub_registry.R
\name{stub_registry}
\alias{stub_registry}
\title{List stubs in the stub registry}
\usage{
stub_registry()
}
\value{
an object of class \code{StubRegistry}, print method gives the
stubs in the registry
}
\description{
List stubs in the stub registry
}
\examples{
# make a stub
stub_request("get", "https://httpbin.org/get") \%>\%
  to_return(body = "success!", status = 200)

# check the stub registry, there should be one in there
stub_registry()

# make another stub
stub_request("get", "https://httpbin.org/get") \%>\%
  to_return(body = "woopsy", status = 404)

# check the stub registry, now there are two there
stub_registry()

# to clear the stub registry
stub_registry_clear()
}
\seealso{
Other stub-registry: 
\code{\link{StubRegistry}},
\code{\link{remove_request_stub}()},
\code{\link{stub_registry_clear}()}
}
\concept{stub-registry}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mocking-disk-writing.R
\name{mocking-disk-writing}
\alias{mocking-disk-writing}
\title{Mocking writing to disk}
\description{
Mocking writing to disk
}
\examples{
\dontrun{
# enable mocking
enable()

# Write to a file before mocked request

# crul
library(crul)
## make a temp file
f <- tempfile(fileext = ".json")
## write something to the file
cat("{\"hello\":\"world\"}\n", file = f)
readLines(f)
## make the stub
stub_request("get", "https://httpbin.org/get") \%>\% 
  to_return(body = file(f))
## make a request
(out <- HttpClient$new("https://httpbin.org/get")$get(disk = f))
out$content
readLines(out$content)

# httr
library(httr)
## make a temp file
f <- tempfile(fileext = ".json")
## write something to the file
cat("{\"hello\":\"world\"}\n", file = f)
readLines(f)
## make the stub
stub_request("get", "https://httpbin.org/get") \%>\% 
  to_return(body = file(f), 
   headers = list('content-type' = "application/json"))
## make a request
## with httr, you must set overwrite=TRUE or you'll get an errror
out <- GET("https://httpbin.org/get", write_disk(f, overwrite=TRUE))
out
out$content
content(out, "text", encoding = "UTF-8")


# Use mock_file to have webmockr handle file and contents

# crul
library(crul)
f <- tempfile(fileext = ".json")
## make the stub
stub_request("get", "https://httpbin.org/get") \%>\% 
  to_return(body = mock_file(f, "{\"hello\":\"mars\"}\n"))
## make a request
(out <- crul::HttpClient$new("https://httpbin.org/get")$get(disk = f))
out$content
readLines(out$content)

# httr
library(httr)
## make a temp file
f <- tempfile(fileext = ".json")
## make the stub
stub_request("get", "https://httpbin.org/get") \%>\% 
  to_return(
    body = mock_file(path = f, payload = "{\"foo\": \"bar\"}"),
    headers = list('content-type' = "application/json")
  )
## make a request
out <- GET("https://httpbin.org/get", write_disk(f))
out
## view stubbed file content
out$content
readLines(out$content)
content(out, "text", encoding = "UTF-8")

# disable mocking
disable()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{wi_th_}
\alias{wi_th_}
\title{This function is defunct.}
\usage{
wi_th_(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StubbedRequest.R
\name{StubbedRequest}
\alias{StubbedRequest}
\title{StubbedRequest}
\description{
stubbed request class underlying \code{\link[=stub_request]{stub_request()}}
}
\examples{
\dontrun{
x <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
x$method
x$uri
x$with(headers = list('User-Agent' = 'R', apple = "good"))
x$to_return(status = 200, body = "foobar", headers = list(a = 5))
x
x$to_s()

# many to_return's
x <- StubbedRequest$new(method = "get", uri = "httpbin.org")
x$to_return(status = 200, body = "foobar", headers = list(a = 5))
x$to_return(status = 200, body = "bears", headers = list(b = 6))
x
x$to_s()

# raw body
x <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
x$to_return(status = 200, body = raw(0), headers = list(a = 5))
x$to_s()
x

x <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
x$to_return(status = 200, body = charToRaw("foo bar"),
  headers = list(a = 5))
x$to_s()
x

# basic auth
x <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
x$with(basic_auth = c("foo", "bar"))
x$to_s()
x

# file path
x <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
f <- tempfile()
x$to_return(status = 200, body = file(f), headers = list(a = 5))
x
x$to_s()
unlink(f)

# to_file(): file path and payload to go into the file
#   payload written to file during mocked response creation
x <- StubbedRequest$new(method = "get", uri = "api.crossref.org")
f <- tempfile()
x$to_return(status = 200, body = mock_file(f, "{\"foo\": \"bar\"}"),
  headers = list(a = 5))
x
x$to_s()
unlink(f)

# uri_regex
(x <- StubbedRequest$new(method = "get", uri_regex = ".+ossref.org"))
x$method
x$uri_regex
x$to_s()

# to timeout
(x <- StubbedRequest$new(method = "get", uri_regex = ".+ossref.org"))
x$to_s()
x$to_timeout()
x$to_s()
x

# to raise
library(fauxpas)
(x <- StubbedRequest$new(method = "get", uri_regex = ".+ossref.org"))
x$to_s()
x$to_raise(HTTPBadGateway)
x$to_s()
x
}
}
\seealso{
\code{\link[=stub_request]{stub_request()}}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{method}}{(xx) xx}

\item{\code{uri}}{(xx) xx}

\item{\code{uri_regex}}{(xx) xx}

\item{\code{uri_parts}}{(xx) xx}

\item{\code{host}}{(xx) xx}

\item{\code{query}}{(xx) xx}

\item{\code{body}}{(xx) xx}

\item{\code{basic_auth}}{(xx) xx}

\item{\code{request_headers}}{(xx) xx}

\item{\code{response_headers}}{(xx) xx}

\item{\code{responses_sequences}}{(xx) xx}

\item{\code{status_code}}{(xx) xx}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{StubbedRequest$new()}}
\item \href{#method-print}{\code{StubbedRequest$print()}}
\item \href{#method-with}{\code{StubbedRequest$with()}}
\item \href{#method-to_return}{\code{StubbedRequest$to_return()}}
\item \href{#method-to_timeout}{\code{StubbedRequest$to_timeout()}}
\item \href{#method-to_raise}{\code{StubbedRequest$to_raise()}}
\item \href{#method-to_s}{\code{StubbedRequest$to_s()}}
\item \href{#method-clone}{\code{StubbedRequest$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{StubbedRequest} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$new(method, uri = NULL, uri_regex = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{the HTTP method (any, head, get, post, put,
patch, or delete). "any" matches any HTTP method. required.}

\item{\code{uri}}{(character) request URI. either this or \code{uri_regex}
required. \pkg{webmockr} can match uri's without the "http" scheme,
but does not match if the scheme is "https". required, unless
\code{uri_regex} given. See \link{UriPattern} for more.}

\item{\code{uri_regex}}{(character) request URI as regex. either this or \code{uri}
required}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{StubbedRequest} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{StubbedRequest} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$print(x, ...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-with"></a>}}
\if{latex}{\out{\hypertarget{method-with}{}}}
\subsection{Method \code{with()}}{
Set expectations for what's given in HTTP request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$with(
  query = NULL,
  body = NULL,
  headers = NULL,
  basic_auth = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{query}}{(list) request query params, as a named list. optional}

\item{\code{body}}{(list) request body, as a named list. optional}

\item{\code{headers}}{(list) request headers as a named list. optional.}

\item{\code{basic_auth}}{(character) basic authentication. optional.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets only
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_return"></a>}}
\if{latex}{\out{\hypertarget{method-to_return}{}}}
\subsection{Method \code{to_return()}}{
Set expectations for what's returned in HTTP response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$to_return(status, body, headers)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{status}}{(numeric) an HTTP status code}

\item{\code{body}}{(list) response body, one of: \code{character}, \code{json},
\code{list}, \code{raw}, \code{numeric}, \code{NULL}, \code{FALSE}, or a file connection
(other connetion types not supported)}

\item{\code{headers}}{(list) named list, response headers. optional.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets whats to be returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_timeout"></a>}}
\if{latex}{\out{\hypertarget{method-to_timeout}{}}}
\subsection{Method \code{to_timeout()}}{
Response should time out
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$to_timeout()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_raise"></a>}}
\if{latex}{\out{\hypertarget{method-to_raise}{}}}
\subsection{Method \code{to_raise()}}{
Response should raise an exception \code{x}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$to_raise(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{(character) an exception message}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
Response as a character string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$to_s()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(character) the response as a string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{StubbedRequest$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/Response.R
\name{Response}
\alias{Response}
\title{Response}
\description{
custom webmockr http response class
}
\examples{
\dontrun{
(x <- Response$new())

x$set_url("https://httpbin.org/get")
x

x$set_request_headers(list('Content-Type' = "application/json"))
x
x$request_headers

x$set_response_headers(list('Host' = "httpbin.org"))
x
x$response_headers

x$set_status(404)
x
x$get_status()

x$set_body("hello world")
x
x$get_body()
# raw body
x$set_body(charToRaw("hello world"))
x
x$get_body()

x$set_exception("exception")
x
x$get_exception()
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{url}}{(character) a url}

\item{\code{body}}{(various) list, character, etc}

\item{\code{content}}{(various) response content/body}

\item{\code{request_headers}}{(list) a named list}

\item{\code{response_headers}}{(list) a named list}

\item{\code{options}}{(character) list}

\item{\code{status_code}}{(integer) an http status code}

\item{\code{exception}}{(character) an exception message}

\item{\code{should_timeout}}{(logical) should the response timeout?}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Response$new()}}
\item \href{#method-print}{\code{Response$print()}}
\item \href{#method-set_url}{\code{Response$set_url()}}
\item \href{#method-get_url}{\code{Response$get_url()}}
\item \href{#method-set_request_headers}{\code{Response$set_request_headers()}}
\item \href{#method-get_request_headers}{\code{Response$get_request_headers()}}
\item \href{#method-set_response_headers}{\code{Response$set_response_headers()}}
\item \href{#method-get_respone_headers}{\code{Response$get_respone_headers()}}
\item \href{#method-set_body}{\code{Response$set_body()}}
\item \href{#method-get_body}{\code{Response$get_body()}}
\item \href{#method-set_status}{\code{Response$set_status()}}
\item \href{#method-get_status}{\code{Response$get_status()}}
\item \href{#method-set_exception}{\code{Response$set_exception()}}
\item \href{#method-get_exception}{\code{Response$get_exception()}}
\item \href{#method-clone}{\code{Response$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Response} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$new(options = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{options}}{(list) a list of options}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Response} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for the \code{Response} class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$print(x, ...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-set_url"></a>}}
\if{latex}{\out{\hypertarget{method-set_url}{}}}
\subsection{Method \code{set_url()}}{
set the url for the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$set_url(url)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{url}}{(character) a url}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets url
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_url"></a>}}
\if{latex}{\out{\hypertarget{method-get_url}{}}}
\subsection{Method \code{get_url()}}{
get the url for the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$get_url()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(character) a url
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_request_headers"></a>}}
\if{latex}{\out{\hypertarget{method-set_request_headers}{}}}
\subsection{Method \code{set_request_headers()}}{
set the request headers for the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$set_request_headers(headers, capitalize = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{headers}}{(list) named list}

\item{\code{capitalize}}{(logical) whether to capitalize first letters of
each header; default: \code{TRUE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets request headers on the response
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_request_headers"></a>}}
\if{latex}{\out{\hypertarget{method-get_request_headers}{}}}
\subsection{Method \code{get_request_headers()}}{
get the request headers for the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$get_request_headers()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(list) request headers, a named list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_response_headers"></a>}}
\if{latex}{\out{\hypertarget{method-set_response_headers}{}}}
\subsection{Method \code{set_response_headers()}}{
set the response headers for the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$set_response_headers(headers, capitalize = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{headers}}{(list) named list}

\item{\code{capitalize}}{(logical) whether to capitalize first letters of
each header; default: \code{TRUE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets response headers on the response
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_respone_headers"></a>}}
\if{latex}{\out{\hypertarget{method-get_respone_headers}{}}}
\subsection{Method \code{get_respone_headers()}}{
get the response headers for the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$get_respone_headers()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(list) response headers, a named list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_body"></a>}}
\if{latex}{\out{\hypertarget{method-set_body}{}}}
\subsection{Method \code{set_body()}}{
set the body of the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$set_body(body, disk = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{body}}{(various types)}

\item{\code{disk}}{(logical) whether its on disk; default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets body on the response
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_body"></a>}}
\if{latex}{\out{\hypertarget{method-get_body}{}}}
\subsection{Method \code{get_body()}}{
get the body of the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$get_body()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
various
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_status"></a>}}
\if{latex}{\out{\hypertarget{method-set_status}{}}}
\subsection{Method \code{set_status()}}{
set the http status of the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$set_status(status)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{status}}{(integer) the http status}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets the http status of the response
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_status"></a>}}
\if{latex}{\out{\hypertarget{method-get_status}{}}}
\subsection{Method \code{get_status()}}{
get the http status of the response
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$get_status()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(integer) the http status
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_exception"></a>}}
\if{latex}{\out{\hypertarget{method-set_exception}{}}}
\subsection{Method \code{set_exception()}}{
set an exception
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$set_exception(exception)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{exception}}{(character) an exception string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned; sets an exception
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_exception"></a>}}
\if{latex}{\out{\hypertarget{method-get_exception}{}}}
\subsection{Method \code{get_exception()}}{
get the exception, if set
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$get_exception()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(character) an exception
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Response$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/RequestPattern.R
\name{UriPattern}
\alias{UriPattern}
\title{UriPattern}
\description{
uri matcher
}
\examples{
# trailing slash
(z <- UriPattern$new(pattern = "http://foobar.com"))
z$matches("http://foobar.com") # TRUE
z$matches("http://foobar.com/") # TRUE

# without scheme
## matches http by default: does not match https by default
(z <- UriPattern$new(pattern = "foobar.com"))
z$matches("http://foobar.com") # TRUE
z$matches("http://foobar.com/") # TRUE
z$matches("https://foobar.com") # FALSE
z$matches("https://foobar.com/") # FALSE
## to match https, you'll have to give the complete url
(z <- UriPattern$new(pattern = "https://foobar.com"))
z$matches("https://foobar.com/") # TRUE
z$matches("http://foobar.com/") # FALSE

# default ports
(z <- UriPattern$new(pattern = "http://foobar.com"))
z$matches("http://foobar.com:80") # TRUE
z$matches("http://foobar.com:80/") # TRUE
z$matches("http://foobar.com:443") # TRUE
z$matches("http://foobar.com:443/") # TRUE

# user info - FIXME, not sure we support this yet
(z <- UriPattern$new(pattern = "http://foobar.com"))
z$matches("http://user:pass@foobar.com")

# regex
(z <- UriPattern$new(regex_pattern = ".+ample\\\\.."))
z$matches("http://sample.org") # TRUE
z$matches("http://example.com") # TRUE
z$matches("http://tramples.net") # FALSE

# add query parameters
(z <- UriPattern$new(pattern = "http://foobar.com"))
z$add_query_params(list(pizza = "cheese", cheese = "cheddar"))
z
z$pattern
z$matches("http://foobar.com?pizza=cheese&cheese=cheddar") # TRUE
z$matches("http://foobar.com?pizza=cheese&cheese=swiss") # FALSE

# query parameters in the uri
(z <- UriPattern$new(pattern = "https://httpbin.org/get?stuff=things"))
z$add_query_params() # have to run this method to gather query params
z$matches("https://httpbin.org/get?stuff=things") # TRUE
z$matches("https://httpbin.org/get?stuff2=things") # FALSE

# regex add query parameters
(z <- UriPattern$new(regex_pattern = "https://foobar.com/.+/order"))
z$add_query_params(list(pizza = "cheese"))
z
z$pattern
z$matches("https://foobar.com/pizzas/order?pizza=cheese") # TRUE
z$matches("https://foobar.com/pizzas?pizza=cheese") # FALSE

# query parameters in the regex uri
(z <- UriPattern$new(regex_pattern = "https://x.com/.+/order\\\\?fruit=apple"))
z$add_query_params() # have to run this method to gather query params
z$matches("https://x.com/a/order?fruit=apple") # TRUE
z$matches("https://x.com/a?fruit=apple") # FALSE

# any pattern
(z <- UriPattern$new(regex_pattern = "stuff\\\\.com.+"))
z$regex
z$pattern
z$matches("http://stuff.com") # FALSE
z$matches("https://stuff.com/stff") # TRUE
z$matches("https://stuff.com/apple?bears=brown&bats=grey") # TRUE
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{pattern}}{(character) pattern holder}

\item{\code{regex}}{a logical}

\item{\code{query_params}}{a list, or \code{NULL} if empty}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{UriPattern$new()}}
\item \href{#method-matches}{\code{UriPattern$matches()}}
\item \href{#method-pattern_matches}{\code{UriPattern$pattern_matches()}}
\item \href{#method-query_params_matches}{\code{UriPattern$query_params_matches()}}
\item \href{#method-extract_query}{\code{UriPattern$extract_query()}}
\item \href{#method-add_query_params}{\code{UriPattern$add_query_params()}}
\item \href{#method-to_s}{\code{UriPattern$to_s()}}
\item \href{#method-clone}{\code{UriPattern$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{UriPattern} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$new(pattern = NULL, regex_pattern = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{pattern}}{(character) a uri, as a character string. if scheme
is missing, it is added (we assume http)}

\item{\code{regex_pattern}}{(character) a uri as a regex character string,
see \link[base:regex]{base::regex}. if scheme is missing, it is added (we assume
http)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{UriPattern} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-matches"></a>}}
\if{latex}{\out{\hypertarget{method-matches}{}}}
\subsection{Method \code{matches()}}{
Match a uri against a pattern
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$matches(uri)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{uri}}{(character) a uri}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pattern_matches"></a>}}
\if{latex}{\out{\hypertarget{method-pattern_matches}{}}}
\subsection{Method \code{pattern_matches()}}{
Match a URI
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$pattern_matches(uri)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{uri}}{(character) a uri}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-query_params_matches"></a>}}
\if{latex}{\out{\hypertarget{method-query_params_matches}{}}}
\subsection{Method \code{query_params_matches()}}{
Match query parameters of a URI
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$query_params_matches(uri)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{uri}}{(character) a uri}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a boolean
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-extract_query"></a>}}
\if{latex}{\out{\hypertarget{method-extract_query}{}}}
\subsection{Method \code{extract_query()}}{
Extract query parameters as a named list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$extract_query(uri)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{uri}}{(character) a uri}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
named list, or \code{NULL} if no query parameters
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add_query_params"></a>}}
\if{latex}{\out{\hypertarget{method-add_query_params}{}}}
\subsection{Method \code{add_query_params()}}{
Add query parameters to the URI
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$add_query_params(query_params)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{query_params}}{(list|character) list or character}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned, updates uri pattern
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
Print pattern for easy human consumption
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$to_s()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UriPattern$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/flipswitch.R
\name{enable}
\alias{enable}
\alias{enabled}
\alias{disable}
\title{Enable or disable webmockr}
\usage{
enable(adapter = NULL, options = list(), quiet = FALSE)

enabled(adapter = "crul")

disable(adapter = NULL, options = list(), quiet = FALSE)
}
\arguments{
\item{adapter}{(character) the adapter name, 'crul' or 'httr'.
one or the other. if none given, we attempt to enable both
adapters}

\item{options}{list of options - ignored for now.}

\item{quiet}{(logical) suppress messages? default: \code{FALSE}}
}
\value{
\code{enable()} and \code{disable()} invisibly returns booleans for
each adapter, as a result of running enable or disable, respectively,
on each \link{HttpLibAdapaterRegistry} object. \code{enabled} returns a
single boolean
}
\description{
Enable or disable webmockr
}
\details{
\code{enable()} enables \pkg{webmockr} for all adapters.
\code{disable()} disables \pkg{webmockr} for all adapters.  \code{enabled()}
answers whether \pkg{webmockr} is enabled for a given adapter
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_return.R
\name{to_return}
\alias{to_return}
\title{Expectation for what's returned from a stubbed request}
\usage{
to_return(.data, ..., .list = list(), times = 1)
}
\arguments{
\item{.data}{input. Anything that can be coerced to a \code{StubbedRequest} class
object}

\item{...}{Comma separated list of named variables. accepts the following:
\code{status}, \code{body}, \code{headers}. See Details for more.}

\item{.list}{named list, has to be one of 'status', 'body',
and/or 'headers'. An alternative to passing in via \code{...}. Don't pass the
same thing to both, e.g. don't pass 'status' to \code{...}, and also 'status' to
this parameter}

\item{times}{(integer) number of times the given response should be
returned; default: 1. value must be greater than or equal to 1. Very large
values probably don't make sense, but there's no maximum value. See
Details.}
}
\value{
an object of class \code{StubbedRequest}, with print method describing
the stub
}
\description{
Set response status code, response body, and/or response headers
}
\details{
Values for status, body, and headers:
\itemize{
\item status: (numeric/integer) three digit status code
\item body: various: \code{character}, \code{json}, \code{list}, \code{raw}, \code{numeric},
\code{NULL}, \code{FALSE}, a file connection (other connetion types
not supported), or a \code{mock_file} function call (see \code{\link[=mock_file]{mock_file()}})
\item headers: (list) a named list, must be named
}

response headers are returned with all lowercase names and the values
are all of type character. if numeric/integer values are given
(e.g., \code{to_return(headers = list(a = 10))}), we'll coerce any
numeric/integer values to character.
}
\note{
see more examples in \code{\link[=stub_request]{stub_request()}}
}
\section{multiple \code{to_return()}}{

You can add more than one \code{to_return()} to a webmockr stub (including
\code{\link[=to_raise]{to_raise()}}, \code{\link[=to_timeout]{to_timeout()}}). Each one is a HTTP response returned.
That is, you'll match to an HTTP request based on \code{stub_request()} and
\code{wi_th()}; the first time the request is made, the first response
is returned; the second time the reqeust is made, the second response
is returned; and so on.

Be aware that webmockr has to track number of requests
(see \code{\link[=request_registry]{request_registry()}}), and so if you use multiple \code{to_return()}
or the \code{times} parameter, you must clear the request registry
in order to go back to mocking responses from the start again.
\code{\link[=webmockr_reset]{webmockr_reset()}} clears the stub registry and  the request registry,
after which you can use multiple responses again (after creating
your stub(s) again of course)
}

\section{Raise vs. Return}{

\code{to_raise()} always raises a stop condition, while \code{to_return(status=xyz)} only
sets the status code on the returned HTTP response object. So if you want to
raise a stop condition then \code{to_raise()} is what you want. But if you
don't want to raise a stop condition use \code{to_return()}. Use cases for each
vary. For example, in a unit test you may have a test expecting a 503 error;
in this case \code{to_raise()} makes sense. In another case, if a unit test
expects to test some aspect of an HTTP response object that httr or crul
typically returns, then you'll want \code{to_return()}.
}

\examples{
# first, make a stub object
foo <- function() {
  stub_request("post", "https://httpbin.org/post")
}

# add status, body and/or headers
foo() \%>\% to_return(status = 200)
foo() \%>\% to_return(body = "stuff")
foo() \%>\% to_return(body = list(a = list(b = "world")))
foo() \%>\% to_return(headers = list(a = 5))
foo() \%>\% 
  to_return(status = 200, body = "stuff", headers = list(a = 5))

# .list - pass in a named list instead
foo() \%>\% to_return(.list = list(body = list(foo = "bar")))

# multiple responses using chained `to_return()`
foo() \%>\% to_return(body = "stuff") \%>\% to_return(body = "things")

# many of the same response using the times parameter
foo() \%>\% to_return(body = "stuff", times = 3)
}
