
## A ‘robots.txt’ Parser and ‘Webbot’/‘Spider’/‘Crawler’ Permissions Checker

[![ropensci\_footer](https://raw.githubusercontent.com/ropensci/robotstxt/master/logo/github_footer.png)](https://ropensci.org)

**Status**

*lines of R code:* 1007, *lines of test code:* 1758

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![](https://badges.ropensci.org/25_status.svg)](https://github.com/ropensci/software-review/issues/25)
<a href="https://travis-ci.org/ropensci/robotstxt"><img src="https://api.travis-ci.org/ropensci/robotstxt.svg?branch=master"><a/>
<a href="https://cran.r-project.org/package=robotstxt"><img src="http://www.r-pkg.org/badges/version/robotstxt"></a>
[![cran
checks](https://cranchecks.info/badges/summary/reshape)](https://cran.r-project.org/web/checks/check_results_reshape.html)
<a href="https://codecov.io/gh/ropensci/robotstxt"><img src="https://codecov.io/gh/ropensci/robotstxt/branch/master/graph/badge.svg" alt="Codecov" /></a>
<img src="http://cranlogs.r-pkg.org/badges/grand-total/robotstxt">
<img src="http://cranlogs.r-pkg.org/badges/robotstxt">

**Development version**

0.7.13 - 2020-08-19 / 20:39:24

**Description**

Provides functions to download and parse ‘robots.txt’ files. Ultimately
the package makes it easy to check if bots (spiders, crawler, scrapers,
…) are allowed to access specific resources on a domain.

**License**

MIT + file LICENSE <br>Peter Meissner \[aut, cre\], Kun Ren \[aut, cph\]
(Author and copyright holder of list\_merge.R.), Oliver Keys \[ctb\]
(original release code review), Rich Fitz John \[ctb\] (original release
code review)

**Citation**

``` r
citation("robotstxt")
```

**BibTex for citing**

``` r
toBibtex(citation("robotstxt"))
```

**Contribution - AKA The-Think-Twice-Be-Nice-Rule**

Please note that this project is released with a Contributor Code of
Conduct. By participating in this project you agree to abide by its
terms:

> As contributors and maintainers of this project, we pledge to respect
> all people who contribute through reporting issues, posting feature
> requests, updating documentation, submitting pull requests or patches,
> and other activities.
> 
> We are committed to making participation in this project a
> harassment-free experience for everyone, regardless of level of
> experience, gender, gender identity and expression, sexual
> orientation, disability, personal appearance, body size, race,
> ethnicity, age, or religion.
> 
> Examples of unacceptable behavior by participants include the use of
> sexual language or imagery, derogatory comments or personal attacks,
> trolling, public or private harassment, insults, or other
> unprofessional conduct.
> 
> Project maintainers have the right and responsibility to remove, edit,
> or reject comments, commits, code, wiki edits, issues, and other
> contributions that are not aligned to this Code of Conduct. Project
> maintainers who do not follow the Code of Conduct may be removed from
> the project team.
> 
> Instances of abusive, harassing, or otherwise unacceptable behavior
> may be reported by opening an issue or contacting one or more of the
> project maintainers.
> 
> This Code of Conduct is adapted from the Contributor Covenant
> (<https://www.contributor-covenant.org/>), version 1.0.0, available at
> <https://www.contributor-covenant.org/version/1/0/0/code-of-conduct/>

## Installation

**Installation and start - stable version**

``` r
install.packages("robotstxt")
library(robotstxt)
```

**Installation and start - development version**

``` r
devtools::install_github("ropensci/robotstxt")
library(robotstxt)
```

## Usage

**Robotstxt class documentation**

``` r
?robotstxt
```

Simple path access right checking (the functional way) …

``` r
library(robotstxt)
options(robotstxt_warn = FALSE)


paths_allowed(
  paths  = c("/api/rest_v1/?doc", "/w/"), 
  domain = "wikipedia.org", 
  bot    = "*"
)
##  wikipedia.org
## [1]  TRUE FALSE

paths_allowed(
  paths = c(
    "https://wikipedia.org/api/rest_v1/?doc", 
    "https://wikipedia.org/w/"
  )
)
##  wikipedia.org                       wikipedia.org
## [1]  TRUE FALSE
```

… or (the object oriented way) …

``` r
library(robotstxt)
options(robotstxt_warn = FALSE)

rtxt <- 
  robotstxt(domain = "wikipedia.org")

rtxt$check(
  paths = c("/api/rest_v1/?doc", "/w/"), 
  bot   = "*"
)
## [1]  TRUE FALSE
```

### Retrieval

Retrieving the robots.txt file for a domain:

``` r
# retrieval
rt <- 
  get_robotstxt("https://petermeissner.de")

# printing
rt
## [robots.txt]
## --------------------------------------
## 
## # just do it - punk
```

### Interpretation

Checking whether or not one is supposadly allowed to access some
resource from a web server is - unfortunately - not just a matter of
downloading and parsing a simple robots.txt file.

First there is no official specification for robots.txt files so every
robots.txt file written and every robots.txt file read and used is an
interpretation. Most of the time we all have a common understanding on
how things are supposed to work but things get more complicated at the
edges.

Some interpretation problems:

  - finding no robots.txt file at the server (e.g. HTTP status code 404)
    implies that everything is allowed
  - subdomains should have there own robots.txt file if not it is
    assumed that everything is allowed
  - redirects involving protocol changes - e.g. upgrading from http to
    https - are followed and considered no domain or subdomain change -
    so whatever is found at the end of the redirect is considered to be
    the robots.txt file for the original domain
  - redirects from subdomain www to the doamin is considered no domain
    change - so whatever is found at the end of the redirect is
    considered to be the robots.txt file for the subdomain originally
    requested

### Event Handling

Because the interpretation of robots.txt rules not just depends on the
rules specified within the file, the package implements an event handler
system that allows to interpret and re-interpret events into rules.

Under the hood the `rt_request_handler()` function is called within
`get_robotstxt()`. This function takes an {httr} request-response object
and a set of event handlers. Processing the request and the handlers it
checks for various events and states around getting the file and reading
in its content. If an event/state happened the event handlers are passed
on to the `request_handler_handler()` along for problem resolution and
collecting robots.txt file transformations:

  - rule priorities decide if rules are applied given the current state
    priority
  - if rules specify signals those are emitted (e.g. error, message,
    warning)
  - often rules imply overwriting the raw content with a suitable
    interpretation given the circumstances the file was (or was not)
    retrieved

Event handler rules can either consist of 4 items or can be functions -
the former being the usual case and that used throughout the package
itself. Functions like `paths_allowed()` do have parameters that allow
passing along handler rules or handler functions.

Handler rules are lists with the following items:

  - `over_write_file_with`: if the rule is triggered and has higher
    priority than those rules applied beforehand (i.e. the new priority
    has an higher value than the old priority) than the robots.txt file
    retrieved will be overwritten by this character vector
  - `signal`: might be `"message"`, `"warning"`, or `"error"` and will
    use the signal function to signal the event/state just handled.
    Signaling a warning or a message might be suppressed by setting the
    function paramter `warn = FALSE`.
  - `cache` should the package be allowed to cache the results of the
    retrieval or not
  - `priority` the priority of the rule specified as numeric value,
    rules with higher priority will be allowed to overwrite robots.txt
    file content changed by rules with lower priority

The package knows the following rules with the following defaults:

  - `on_server_error` :
  - given a server error - the server is unable to serve a file - we
    assume that something is terrible wrong and forbid all paths for the
    time being but do not cache the result so that we might get an
    updated file later on

<!-- end list -->

``` r
on_server_error_default
## $over_write_file_with
## [1] "User-agent: *\nDisallow: /"
## 
## $signal
## [1] "error"
## 
## $cache
## [1] FALSE
## 
## $priority
## [1] 20
```

  - `on_client_error` :
  - client errors encompass all HTTP status 4xx status codes except 404
    which is handled directly
  - despite the fact that there are a lot of codes that might indicate
    that the client has to take action (authentication, billing, … see:
    <https://de.wikipedia.org/wiki/HTTP-Statuscode>) in the case of
    retrieving robots.txt with simple GET request things should just
    work and any client error is treated as if there is no file
    available and thus scraping is generally allowed

<!-- end list -->

``` r
on_client_error_default
## $over_write_file_with
## [1] "User-agent: *\nAllow: /"
## 
## $signal
## [1] "warning"
## 
## $cache
## [1] TRUE
## 
## $priority
## [1] 19
```

  - `on_not_found` :
  - HTTP status code 404 has its own handler but is treated the same
    ways other client errors: if there is no file available and thus
    scraping is generally allowed

<!-- end list -->

``` r
on_not_found_default
## $over_write_file_with
## [1] "User-agent: *\nAllow: /"
## 
## $signal
## [1] "warning"
## 
## $cache
## [1] TRUE
## 
## $priority
## [1] 1
```

  - `on_redirect` :
  - redirects are ok - often redirects redirect from HTTP schema to
    HTTPS - robotstxt will use whatever content it has been redirected
    to

<!-- end list -->

``` r
on_redirect_default
## $cache
## [1] TRUE
## 
## $priority
## [1] 3
```

  - `on_domain_change` :
  - domain changes are handled as if the robots.txt file did not exist
    and thus scraping is generally allowed

<!-- end list -->

``` r
on_domain_change_default
## $signal
## [1] "warning"
## 
## $cache
## [1] TRUE
## 
## $priority
## [1] 4
```

  - `on_file_type_mismatch` :
  - if {robotstxt} gets content with content type other than text it
    probably is not a robotstxt file, this situation is handled as if no
    file was provided and thus scraping is generally allowed

<!-- end list -->

``` r
on_file_type_mismatch_default
## $over_write_file_with
## [1] "User-agent: *\nAllow: /"
## 
## $signal
## [1] "warning"
## 
## $cache
## [1] TRUE
## 
## $priority
## [1] 6
```

  - `on_suspect_content` :
  - if {robotstxt} cannot parse it probably is not a robotstxt file,
    this situation is handled as if no file was provided and thus
    scraping is generally allowed

<!-- end list -->

``` r
on_suspect_content_default
## $over_write_file_with
## [1] "User-agent: *\nAllow: /"
## 
## $signal
## [1] "warning"
## 
## $cache
## [1] TRUE
## 
## $priority
## [1] 7
```

### Design Map for Event/State Handling

**from version 0.7.x onwards**

While previous releases were concerned with implementing parsing and
permission checking and improving performance the 0.7.x release will be
about robots.txt retrieval foremost. While retrieval was implemented
there are corner cases in the retrieval stage that very well influence
the interpretation of permissions granted.

**Features and Problems handled:**

  - now handles corner cases of retrieving robots.txt files
  - e.g. if no robots.txt file is available this basically means “you
    can scrape it all”
  - but there are further corner cases (what if there is a server error,
    what if redirection takes place, what is redirection takes place to
    different domains, what if a file is returned but it is not
    parsable, or is of format HTML or JSON, …)

**Design Decisions**

1.  the whole HTTP request-response-chain is checked for certain
    event/state types
      - server error
      - client error
      - file not found (404)
      - redirection
      - redirection to another domain
2.  the content returned by the HTTP is checked against
      - mime type / file type specification mismatch
      - suspicious content (file content does seem to be JSON, HTML, or
        XML instead of robots.txt)
3.  state/event handler define how these states and events are handled
4.  a handler handler executes the rules defined in individual handlers
5.  handler can be overwritten
6.  handler defaults are defined that they should always do the right
    thing
7.  handler can …
      - overwrite the content of a robots.txt file (e.g. allow/disallow
        all)
      - modify how problems should be signaled: error, warning, message,
        none
      - if robots.txt file retrieval should be cached or not
8.  problems (no matter how they were handled) are attached to the
    robots.txt’s as attributes, allowing for …
      - transparency
      - reacting post-mortem to the problems that occured
9.  all handler (even the actual execution of the HTTP-request) can be
    overwritten at runtime to inject user defined behaviour beforehand

### Warnings

By default all functions retrieving robots.txt files will warn if there
are

  - any HTTP events happening while retrieving the file (e.g. redirects)
    or
  - the content of the file does not seem to be a valid robots.txt file.

The warnings in the following example can be turned of in three ways:

(example)

``` r
library(robotstxt)

paths_allowed("petermeissner.de")
##  petermeissner.de
## [1] TRUE
```

(solution 1)

``` r
library(robotstxt)

suppressWarnings({
  paths_allowed("petermeissner.de")
})
##  petermeissner.de
## [1] TRUE
```

(solution 2)

``` r
library(robotstxt)

paths_allowed("petermeissner.de", warn = FALSE)
##  petermeissner.de
## [1] TRUE
```

(solution 3)

``` r
library(robotstxt)

options(robotstxt_warn = FALSE)

paths_allowed("petermeissner.de")
##  petermeissner.de
## [1] TRUE
```

### Inspection and Debugging

The robots.txt files retrieved are basically mere character vectors:

``` r
rt <- 
  get_robotstxt("petermeissner.de")

as.character(rt)
## [1] "# just do it - punk\n"

cat(rt)
## # just do it - punk
```

The last HTTP request is stored in an object

``` r
rt_last_http$request
## Response [https://petermeissner.de/robots.txt]
##   Date: 2020-09-03 19:05
##   Status: 200
##   Content-Type: text/plain
##   Size: 20 B
## # just do it - punk
```

But they also have some additional information stored as attributes.

``` r
names(attributes(rt))
## [1] "problems" "cached"   "request"  "class"
```

Events that might change the interpretation of the rules found in the
robots.txt file:

``` r
attr(rt, "problems")
## $on_redirect
## $on_redirect[[1]]
## $on_redirect[[1]]$status
## [1] 301
## 
## $on_redirect[[1]]$location
## [1] "https://petermeissner.de/robots.txt"
## 
## 
## $on_redirect[[2]]
## $on_redirect[[2]]$status
## [1] 200
## 
## $on_redirect[[2]]$location
## NULL
```

The {httr} request-response object that allwos to dig into what exactly
was going on in the client-server exchange.

``` r
attr(rt, "request")
## Response [https://petermeissner.de/robots.txt]
##   Date: 2020-09-03 19:05
##   Status: 200
##   Content-Type: text/plain
##   Size: 20 B
## # just do it - punk
```

… or lets us retrieve the original content given back by the server:

``` r
httr::content(
  x        = attr(rt, "request"), 
  as       = "text",
  encoding = "UTF-8"
)
## [1] "# just do it - punk\n"
```

… or have a look at the actual HTTP request issued and all response
headers given back by the server:

``` r
# extract request-response object
rt_req <- 
  attr(rt, "request")

# HTTP request
rt_req$request
## <request>
## GET http://petermeissner.de/robots.txt
## Output: write_memory
## Options:
## * useragent: libcurl/7.64.1 r-curl/4.3 httr/1.4.1
## * ssl_verifypeer: 1
## * httpget: TRUE
## Headers:
## * Accept: application/json, text/xml, application/xml, */*
## * user-agent: R version 3.6.3 (2020-02-29)

# response headers
rt_req$all_headers
## [[1]]
## [[1]]$status
## [1] 301
## 
## [[1]]$version
## [1] "HTTP/1.1"
## 
## [[1]]$headers
## $server
## [1] "nginx/1.10.3 (Ubuntu)"
## 
## $date
## [1] "Thu, 03 Sep 2020 19:05:45 GMT"
## 
## $`content-type`
## [1] "text/html"
## 
## $`content-length`
## [1] "194"
## 
## $connection
## [1] "keep-alive"
## 
## $location
## [1] "https://petermeissner.de/robots.txt"
## 
## attr(,"class")
## [1] "insensitive" "list"       
## 
## 
## [[2]]
## [[2]]$status
## [1] 200
## 
## [[2]]$version
## [1] "HTTP/1.1"
## 
## [[2]]$headers
## $server
## [1] "nginx/1.10.3 (Ubuntu)"
## 
## $date
## [1] "Thu, 03 Sep 2020 19:05:45 GMT"
## 
## $`content-type`
## [1] "text/plain"
## 
## $`content-length`
## [1] "20"
## 
## $`last-modified`
## [1] "Thu, 03 Sep 2020 15:33:01 GMT"
## 
## $connection
## [1] "keep-alive"
## 
## $etag
## [1] "\"5f510cad-14\""
## 
## $`accept-ranges`
## [1] "bytes"
## 
## attr(,"class")
## [1] "insensitive" "list"
```

### Transformation

For convenience the package also includes a `as.list()` method for
robots.txt files.

``` r
as.list(rt)
## $content
## [1] "# just do it - punk\n"
## 
## $robotstxt
## [1] "# just do it - punk\n"
## 
## $problems
## $problems$on_redirect
## $problems$on_redirect[[1]]
## $problems$on_redirect[[1]]$status
## [1] 301
## 
## $problems$on_redirect[[1]]$location
## [1] "https://petermeissner.de/robots.txt"
## 
## 
## $problems$on_redirect[[2]]
## $problems$on_redirect[[2]]$status
## [1] 200
## 
## $problems$on_redirect[[2]]$location
## NULL
## 
## 
## 
## 
## $request
## Response [https://petermeissner.de/robots.txt]
##   Date: 2020-09-03 19:05
##   Status: 200
##   Content-Type: text/plain
##   Size: 20 B
## # just do it - punk
```

### Caching

The retrieval of robots.txt files is cached on a per R-session basis.
Restarting an R-session will invalidate the cache. Also using the the
function parameter `froce = TRUE` will force the package to re-retrieve
the robots.txt file.

``` r
paths_allowed("petermeissner.de/I_want_to_scrape_this_now", force = TRUE, verbose = TRUE)
##  petermeissner.de                      rt_robotstxt_http_getter: force http get
## [1] TRUE
paths_allowed("petermeissner.de/I_want_to_scrape_this_now",verbose = TRUE)
##  petermeissner.de                      rt_robotstxt_http_getter: cached http get
## [1] TRUE
```

## More information

  - <https://www.robotstxt.org/norobots-rfc.txt>
  - [Have a look at the vignette at
    https://cran.r-project.org/package=robotstxt/vignettes/using\_robotstxt.html](https://cran.r-project.org/package=robotstxt/vignettes/using_robotstxt.html)
  - [Google on
    robots.txt](https://developers.google.com/search/reference/robots_txt?hl=en)
  - <https://wiki.selfhtml.org/wiki/Grundlagen/Robots.txt>
  - <https://support.google.com/webmasters/answer/6062608?hl=en>
  - <https://www.robotstxt.org/robotstxt.html>
NEWS robotstxt
==========================================================================


0.7.13 | 2020-09-03
--------------------------------------------------------------------------

- CRAN compliance - prevent URL forwarding (HTTP 301): add www to URLs


0.7.12 | 2020-09-03
--------------------------------------------------------------------------

- CRAN compliance - prevent URL forwarding (HTTP 301): add trailing slashes to URLs



0.7.11 | 2020-09-02
--------------------------------------------------------------------------

- CRAN compliance - LICENCE file wording; prevent URL forwarding (HTTP 301)




0.7.10 | 2020-08-19
--------------------------------------------------------------------------

- fix problem in parse_robotstxt() - comment in last line of robots.txt file would lead to errornous parsing - reported by @gittaca, https://github.com/ropensci/robotstxt/pull/59 and https://github.com/ropensci/robotstxt/issues/60





0.7.9 | 2020-08-02
--------------------------------------------------------------------------

- fix problem is_valid_robotstxt() - robots.txt validity check was to lax - reported by @gittaca, https://github.com/ropensci/robotstxt/issues/58





0.7.8 | 2020-07-22
--------------------------------------------------------------------------

- fix problem with domain name extraction - reported by @gittaca, https://github.com/ropensci/robotstxt/issues/57
- fix problem with vArYING CasE in robots.txt field names - reported by @steffilazerte, https://github.com/ropensci/robotstxt/issues/55





0.7.7 | 2020-06-17
--------------------------------------------------------------------------

- fix problem in rt_request_handler - reported by @MHWauben https://github.com/dmi3kno/polite/issues/28 - patch by @dmi3kno





0.7.6 | 2020-06-13
--------------------------------------------------------------------------

- make info whether or not results were cached available - requested by @dmi3kno, https://github.com/ropensci/robotstxt/issues/53




0.7.5 | 2020-06-07
--------------------------------------------------------------------------

- **fix** passing through more parameters from robotstxt() to get_robotstxt() - reported and implemented by @dmi3kno





0.7.3 | 2020-05-29
--------------------------------------------------------------------------

- **minor** : improve printing of robots.txt
- add request data as attribute to robots.txt
- add `as.list()` method for robots.txt
- adding several paragrpahs to the README file
- **major** : finishing handlers - quality check, documentation
- **fix** : Partial matching warnings #51 - reported by @mine-cetinkaya-rundel





0.7.2 | 2020-05-04
--------------------------------------------------------------------------

- **minor** : changes in dependencies were introducing errors when no scheme/protocoll was provided in URL -- fixed https://github.com/ropensci/robotstxt/issues/50





0.7.1 | 2018-01-09
--------------------------------------------------------------------------

- **minor** : modifying robots.txt parser to be more robust against different formatting of robots.txt files -- fixed https://github.com/ropensci/robotstxt/issues/48





0.7.0 | 2018-11-27
--------------------------------------------------------------------------

- **major** : introducing http handler to allow for better interpretation of robots.txt files in case of certain events: redirects, server error, client error, suspicous content, ...



0.6.4 | 2018-09-14
--------------------------------------------------------------------------

- **minor** : pass through of parameter for content encoding 



0.6.3 | 2018-09-14
--------------------------------------------------------------------------

- **minor** : introduced parameter encoding to `get_robotstxt()` that defaults to "UTF-8" which does the content function anyways - but now it will not complain about it
- **minor** : added comment to help files specifying use of trailing slash in paths pointing to folders in `paths_allowed` and `robotstxt`.




0.6.2 | 2018-07-18
--------------------------------------------------------------------------

- **minor** : changed from `future::future_lapply()` to `future.apply::future_lapply()` to make package compatible with versions of future after 1.8.1




0.6.1 | 2018-05-30
--------------------------------------------------------------------------

- **minor** : package was moved to other repo location and project status badge was added



0.6.0 | 2018-02-10
--------------------------------------------------------------------------

- **change/fix** check function paths_allowed() would not return correct result in some edge cases, indicating that spiderbar/rep-cpp check method is more reliable and shall be the default and only  method: [see 1](https://github.com/ropensci/robotstxt/issues/22), [see 2](https://github.com/hrbrmstr/spiderbar/issues/2), [see 3](https://github.com/seomoz/rep-cpp/issues/33)




0.5.2 | 2017-11-12
--------------------------------------------------------------------------

- **fix** : rt_get_rtxt() would break on Windows due trying to readLines() from folder




0.5.1 | 2017-11-11
--------------------------------------------------------------------------

- **change** : spiderbar is now non-default second (experimental) check method
- **fix** : there were warnings in case of multiple domain guessing



0.5.0 | 2017-10-07
--------------------------------------------------------------------------

- **feature** : spiderbar's can_fetch() was added, now one can choose which check method to use for checking access rights 
- **feature** : use futures (from package future) to speed up retrieval and parsing
- **feature** : now there is a `get_robotstxts()` function wich is a 'vectorized' version of `get_robotstxt()`
- **feature** : `paths_allowed()` now allows checking via either robotstxt parsed robots.txt files or via functionality provided by the spiderbar package (the latter should be faster by approximatly factor 10)
- **feature** : various functions now have a ssl_verifypeer option (analog to CURL option https://curl.haxx.se/libcurl/c/CURLOPT_SSL_VERIFYPEER.html) which might help with robots.txt file retrieval in some cases
- **change** : user_agent for robots.txt file retrieval will now default to: `sessionInfo()$R.version$version.string` 
- **change** : robotstxt now assumes it knows how to parse --> if it cannot parse it assumes that it got no valid robots.txt file meaning that there are no restrictions
- **fix** : valid_robotstxt would not accept some actual valid robotstxt files



0.4.1 | 2017-08-20
--------------------------------------------------------------------------

- **restructure** : put each function in separate file
- **fix** : parsing would go bonkers for robots.txt of cdc.gov (e.g. combining all robots with all permissions) due to errornous handling of carriage return character (reported by @hrbrmstr - thanks)



0.4.0 | 2017-07-14
--------------------------------------------------------------------------

- **user_agent** parameter **added** to robotstxt() and paths_allowed to allow for user defined HTTP user-agent send when retrieving robots.txt file from domain



0.3.4 | 2017-07-08
--------------------------------------------------------------------------

- **fix** : non robots.txt files (e.g. html files returned by server instead of the requested robots.txt / facebook.com) would be handled as if it were non existent / empty files (reported by @simonmunzert - thanks)
- **fix** : UTF-8 encoded robots.txt with BOM (byte order mark) would break parsing although files were otherwise valid robots.txt files




0.3.3 | 2016-12-10
--------------------------------------------------------------------------

- updating NEWS file and switching to NEWS.md





0.3.2 | 2016-04-28 
--------------------------------------------------------------------------

- CRAN publication





0.3.1 | 2016-04-27 
--------------------------------------------------------------------------

- get_robotstxt() tests for HTTP errors and handles them, warnings might be suppressed while un-plausible HTTP status codes will lead to stoping the function https://github.com/ropenscilabs/robotstxt#5

- dropping R6 dependency and use list implementation instead https://github.com/ropenscilabs/robotstxt#6

- use caching for get_robotstxt() https://github.com/ropenscilabs/robotstxt#7 / https://github.com/ropenscilabs/robotstxt/commit/90ad735b8c2663367db6a9d5dedbad8df2bc0d23

- make explicit, less error prone usage of httr::content(rtxt) https://github.com/ropenscilabs/robotstxt#

- replace usage of missing for parameter check with explicit NULL as default value for parameter https://github.com/ropenscilabs/robotstxt#9

- partial match useragent / useragents https://github.com/ropenscilabs/robotstxt#10

- explicit declaration encoding: encoding="UTF-8" in httr::content() https://github.com/ropenscilabs/robotstxt#11





version 0.1.2 // 2016-02-08 ...
--------------------------------------------------------------------------

- first feature complete version on CRAN





This is a re-submission.

## Your comments

You still have not followed the moved content in

   Found the following (possibly) invalid URLs:
     URL: https://contributor-covenant.org/ (moved to
https://www.contributor-covenant.org/)
       From: README.md
       Status: 200
       Message: OK
     URL: https://contributor-covenant.org/version/1/0/0/ (moved to
https://www.contributor-covenant.org/version/1/0/0/)
       From: README.md
       Status: 200
       Message: OK

Please fix and resubmit.




## My actions

I changed the URLs.

Sorry, that this is such a mess - I am really trying to get them fixed and
minimize workload on behalf of CRAN. 

(HTTP 301 forwards are however also a very common in the ever changing internet 
and in contrast to 4xx and 5xx HTTP status codes no errors as such - 
"301 moved permanetly" can mean anything really including any kind of 
redirect for any reason. 
It would be great if those test would be 
part of the normal R CMD check?. ...So I (and all the others) 
can fetch them locally or via the CI pipelines set up and ready to catch problems 
beforhand). 





## Test environments

- Ubuntu precise (16.04) (on travis-ci: old, current, devel; https://travis-ci.org/github/ropensci/robotstxt ) --> OK

- Win10 lokal with R 3.6.3 --> OK
- win-builder   - devel    --> OK
- win-builder   - release  --> OK


## R CMD check results

0 errors | 0 warnings | 0 notes



## Reverse Dependency Checks

Seems ok.



---
output: github_document
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##",
  fig.path = "README-"
)
```


```{r, include=FALSE}
options("width"=110)
tmp <- packageDescription( basename(getwd()) )
```


```{r, results='asis', echo=FALSE}
cat("##", tmp$Title)
```



```{r, include=FALSE}
filelist.R   <- list.files("R", recursive = TRUE, pattern="\\.R$", ignore.case = TRUE, full.names = TRUE)
filelist.tests   <- list.files("tests", recursive = TRUE, pattern="\\.R$", ignore.case = TRUE, full.names = TRUE)
filelist.cpp <- list.files("src", recursive = TRUE, pattern="\\.cpp$", ignore.case = TRUE, full.names = TRUE)
lines.R      <- unlist(lapply(filelist.R, readLines, warn = FALSE))
lines.tests  <- unlist(lapply(filelist.tests, readLines, warn = FALSE))
lines.cpp    <- unlist(lapply(filelist.cpp, readLines, warn = FALSE))
length.R     <- length(grep("(^\\s*$)|(^\\s*#)|(^\\s*//)", lines.R,  value = TRUE, invert = TRUE))
length.tests <- length(grep("(^\\s*$)|(^\\s*#)|(^\\s*//)", lines.tests,  value = TRUE, invert = TRUE))
length.cpp   <- length(grep("(^\\s*$)|(^\\s*#)|(^\\s*//)", lines.cpp,  value = TRUE, invert = TRUE))
```




[![ropensci\_footer](https://raw.githubusercontent.com/ropensci/robotstxt/master/logo/github_footer.png)](https://ropensci.org)



**Status**


*lines of R code:* `r length.R`, *lines of test code:* `r length.tests`

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![](https://badges.ropensci.org/25_status.svg)](https://github.com/ropensci/software-review/issues/25)
<a href="https://travis-ci.org/ropensci/robotstxt"><img src="https://api.travis-ci.org/ropensci/robotstxt.svg?branch=master"><a/>
<a href="https://cran.r-project.org/package=robotstxt"><img src="http://www.r-pkg.org/badges/version/robotstxt"></a>
[![cran checks](https://cranchecks.info/badges/summary/reshape)](https://cran.r-project.org/web/checks/check_results_reshape.html)
<a href="https://codecov.io/gh/ropensci/robotstxt"><img src="https://codecov.io/gh/ropensci/robotstxt/branch/master/graph/badge.svg" alt="Codecov" /></a>
<img src="http://cranlogs.r-pkg.org/badges/grand-total/robotstxt">
<img src="http://cranlogs.r-pkg.org/badges/robotstxt">





**Development version**

```{r, include=FALSE}
source_files <- 
  grep(
    "/R/|/src/|/tests/",
    list.files(recursive = TRUE, full.names = TRUE), 
    value = TRUE
  )
last_change <- 
  as.character(
    format(max(file.info(source_files)$mtime), tz="UTC")
  )
```


```{r, results='asis', echo=FALSE}
cat(tmp$Version)
cat(" - ")
cat(stringr::str_replace(last_change, " ", " / "))
```

**Description**

```{r, results='asis', echo=FALSE}
cat(tmp$Description)
```


**License**

```{r, results='asis', echo=FALSE}
cat(tmp$License, "<br>")
cat(tmp$Author)
```




**Citation**


```{r, results='asis',  eval=FALSE}
citation("robotstxt")
```

**BibTex for citing**

```{r, eval=FALSE}
toBibtex(citation("robotstxt"))
```




**Contribution - AKA The-Think-Twice-Be-Nice-Rule**

Please note that this project is released with a Contributor Code of Conduct. By participating in this project you agree to abide by its terms:

> As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.
> 
> We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.
> 
> Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.
> 
> Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.
> 
> Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.
> 
> This Code of Conduct is adapted from the Contributor Covenant 
(https://www.contributor-covenant.org/), version 1.0.0, available at 
https://www.contributor-covenant.org/version/1/0/0/code-of-conduct/



## Installation

**Installation and start - stable version**

```{r, eval=FALSE}
install.packages("robotstxt")
library(robotstxt)
```


**Installation and start - development version**

```{r, eval=FALSE}
devtools::install_github("ropensci/robotstxt")
library(robotstxt)
```



## Usage

**Robotstxt class documentation**

```{r, eval=FALSE}
?robotstxt
```


Simple path access right checking (the functional way) ... 

```{r}
library(robotstxt)
options(robotstxt_warn = FALSE)


paths_allowed(
  paths  = c("/api/rest_v1/?doc", "/w/"), 
  domain = "wikipedia.org", 
  bot    = "*"
)

paths_allowed(
  paths = c(
    "https://wikipedia.org/api/rest_v1/?doc", 
    "https://wikipedia.org/w/"
  )
)
```

... or (the object oriented way) ...

```{r}
library(robotstxt)
options(robotstxt_warn = FALSE)

rtxt <- 
  robotstxt(domain = "wikipedia.org")

rtxt$check(
  paths = c("/api/rest_v1/?doc", "/w/"), 
  bot   = "*"
)
```


### Retrieval

Retrieving the robots.txt file for a domain: 

```{r}
# retrieval
rt <- 
  get_robotstxt("https://petermeissner.de")

# printing
rt
```



###  Interpretation

Checking whether or not one is supposadly allowed to access some resource from a 
web server is - unfortunately - not just a matter of downloading and parsing a 
simple robots.txt file. 

First there is no official specification for robots.txt files so every robots.txt 
file written and every robots.txt file read and used is an interpretation. Most of
the time we all have a common understanding on how things are supposed to work 
but things get more complicated at the edges. 

Some interpretation problems:

- finding no robots.txt file at the server (e.g. HTTP status code 404) implies that everything is allowed
- subdomains should have there own robots.txt file if not it is assumed that everything is allowed
- redirects involving protocol changes - e.g. upgrading from http to https - are followed and considered no domain or subdomain change - so whatever is found at the end of the redirect is considered to be the robots.txt file for the original domain
- redirects from subdomain www to the doamin is considered no domain change  - so whatever is found at the end of the redirect is considered to be the robots.txt file for the subdomain originally requested



### Event Handling

Because the interpretation of robots.txt rules not just depends on the rules specified within the file,
the package implements an event handler system that allows to interpret and re-interpret events into rules. 

Under the hood the `rt_request_handler()` function is called within `get_robotstxt()`. 
This function takes an {httr} request-response object and a set of event handlers. 
Processing the request and the handlers it checks for various events and states 
around getting the file and reading in its content. If an  event/state happened
the event handlers are passed on to the `request_handler_handler()` along for 
problem resolution and collecting robots.txt file transformations:

- rule priorities decide if rules are applied given the current state priority
- if rules specify signals those are emitted (e.g. error, message, warning)
- often rules imply overwriting the raw content with a suitable interpretation given the circumstances the file was (or was not) retrieved


Event handler rules can either consist of 4 items or can be functions - the former being the usual case and that used throughout the package itself. 
Functions like `paths_allowed()` do have parameters that allow passing 
along handler rules or handler functions.

Handler rules are lists with the following items: 

- `over_write_file_with`: if the rule is triggered and has higher priority than those rules applied beforehand (i.e. the new priority has an higher value than the old priority) than the robots.txt file retrieved will be overwritten by this character vector
- `signal`: might be `"message"`, `"warning"`, or `"error"` and will use the signal function to signal the event/state just handled. Signaling a warning or a message might be suppressed by setting the function paramter `warn = FALSE`.
- `cache` should the package be allowed to cache the results of the retrieval or not
- `priority` the priority of the rule specified as numeric value, rules with higher priority will be allowed to overwrite robots.txt file content changed by rules with lower priority



The package knows the following rules with the following defaults: 

- `on_server_error       ` : 
- given a server error - the server is unable to serve a file - we assume that something is terrible wrong and forbid all paths for the time being but do not cache the result so that we might get an updated file later on

```{r}
on_server_error_default
```


- `on_client_error       ` : 
- client errors encompass all HTTP status 4xx status codes except 404 which is handled directly
- despite the fact that there are a lot of codes that might indicate that the client has to take action (authentication, billing, ... see: https://de.wikipedia.org/wiki/HTTP-Statuscode) in the case of retrieving robots.txt with simple GET request things should just work and any client error is treated as if there is no file available and thus scraping is generally allowed

```{r}
on_client_error_default
```

- `on_not_found          ` : 
- HTTP status code 404 has its own handler but is treated the same ways other client errors: if there is no file available and thus scraping is generally allowed


```{r}
on_not_found_default
```

- `on_redirect           ` : 
- redirects are ok - often redirects redirect from HTTP schema to HTTPS - robotstxt will use whatever content it has been redirected to

```{r}
on_redirect_default
```

- `on_domain_change      ` : 
- domain changes are handled as if the robots.txt file did not exist and thus scraping is generally allowed

```{r}
on_domain_change_default
```

- `on_file_type_mismatch ` : 
- if {robotstxt} gets content with content type other than text it probably is not a robotstxt file, this situation is handled as if no file was provided and thus scraping is generally allowed

```{r}
on_file_type_mismatch_default
```

- `on_suspect_content    ` :
- if {robotstxt} cannot parse it probably is not a robotstxt file, this situation is handled as if no file was provided and thus scraping is generally allowed


```{r}
on_suspect_content_default
```



### Design Map for Event/State Handling

**from version 0.7.x onwards**

While previous releases were concerned with implementing parsing and permission checking and improving performance the 0.7.x release will be about robots.txt retrieval foremost. While retrieval was implemented there are corner cases in the retrieval stage that very well influence the interpretation of permissions granted.


**Features and Problems handled:**

- now handles corner cases of retrieving robots.txt files
- e.g. if no robots.txt file is available this basically means "you can scrape it all" 
- but there are further corner cases (what if there is a server error, what if redirection takes place, what is redirection takes place to different domains, what if a file is returned but it is not parsable, or is of format HTML or JSON, ...)


**Design Decisions**

1. the whole HTTP request-response-chain is checked for certain event/state types
    - server error
    - client error
    - file not found (404)
    - redirection
    - redirection to another domain
2. the content returned by the HTTP is checked against 
    - mime type / file type specification mismatch 
    - suspicious content (file content does seem to be JSON, HTML, or XML instead of  robots.txt)
3. state/event handler define how these states and events are handled
4. a handler handler executes the rules defined in individual handlers
5. handler can be overwritten
6. handler defaults are defined that they should always do the right thing
7. handler can ...
    - overwrite the content of a robots.txt file (e.g. allow/disallow all)
    - modify how problems should be signaled: error, warning, message, none
    - if robots.txt file retrieval should be cached or not
8. problems (no matter how they were handled) are attached to the robots.txt's as attributes, allowing for ... 
    - transparency
    - reacting post-mortem to the problems that occured
9. all handler (even the actual execution of the HTTP-request) can be overwritten at runtime to inject user defined behaviour beforehand









### Warnings

By default all functions retrieving robots.txt files will warn if there are 

- any HTTP events happening while retrieving the file (e.g. redirects) or
- the content of the file does not seem to be a valid robots.txt file.

The warnings in the following example can be turned of in three ways: 


```{r, include = FALSE}
options("robotstxt_warn" = TRUE)
```

(example)
```{r}
library(robotstxt)

paths_allowed("petermeissner.de")
```



(solution 1)
```{r}
library(robotstxt)

suppressWarnings({
  paths_allowed("petermeissner.de")
})
```


(solution 2)
```{r}
library(robotstxt)

paths_allowed("petermeissner.de", warn = FALSE)
```



(solution 3)
```{r}
library(robotstxt)

options(robotstxt_warn = FALSE)

paths_allowed("petermeissner.de")
```



### Inspection and Debugging

The robots.txt files retrieved are basically mere character vectors:

```{r}
rt <- 
  get_robotstxt("petermeissner.de")

as.character(rt)

cat(rt)
```

The last HTTP request is stored in an object 

```{r}
rt_last_http$request
```



But they also have some additional information stored as attributes.


```{r}
names(attributes(rt))
```

Events that might change the interpretation of the rules found in the robots.txt file:

```{r}
attr(rt, "problems")
```

The {httr} request-response object that allwos to dig into what exactly was going on in the client-server exchange. 

```{r}
attr(rt, "request")
```

... or lets us retrieve the original content given back by the server:

```{r}
httr::content(
  x        = attr(rt, "request"), 
  as       = "text",
  encoding = "UTF-8"
)
```


... or have a look at the actual HTTP request issued and all response headers given back by the server:

```{r}
# extract request-response object
rt_req <- 
  attr(rt, "request")

# HTTP request
rt_req$request

# response headers
rt_req$all_headers
```






### Transformation


For convenience the package also includes a `as.list()` method for robots.txt files. 

```{r}
as.list(rt)
```



### Caching 

The retrieval of robots.txt files is cached on a per R-session basis. 
Restarting an R-session will invalidate the cache. Also using the the 
function parameter `froce = TRUE` will force the package to re-retrieve the 
robots.txt file. 

```{r}
paths_allowed("petermeissner.de/I_want_to_scrape_this_now", force = TRUE, verbose = TRUE)
paths_allowed("petermeissner.de/I_want_to_scrape_this_now",verbose = TRUE)
```



## More information

- https://www.robotstxt.org/norobots-rfc.txt
- [Have a look at the vignette at https://cran.r-project.org/package=robotstxt/vignettes/using_robotstxt.html ](https://cran.r-project.org/package=robotstxt/vignettes/using_robotstxt.html)
- [Google on robots.txt](https://developers.google.com/search/reference/robots_txt?hl=en)
- https://wiki.selfhtml.org/wiki/Grundlagen/Robots.txt
- https://support.google.com/webmasters/answer/6062608?hl=en
- https://www.robotstxt.org/robotstxt.html



---
title: "Using Robotstxt"
author: "Peter Meissner"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    toc: true
    css: style.css
vignette: >
  %\VignetteIndexEntry{using_robotstxt}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Description 

The package provides a simple 'robotstxt' class and accompanying methods to parse and check 'robots.txt' files. Data fields are provided as data frames and vectors. Permissions can be checked by providing path character vectors and optional bot names.
    
    
# Robots.txt files

Robots.txt files are a way to kindly ask webbots, spiders, crawlers, wanderers and the like to access or not access certain parts of a webpage. The de facto  'standard' never made it beyond a informal ["Network Working Group INTERNET DRAFT"](http://www.robotstxt.org/norobots-rfc.txt). Nonetheless, the use of  robots.txt files is widespread (e.g. https://en.wikipedia.org/robots.txt, https://www.google.com/robots.txt) and bots from Google, Yahoo and the like will adhere to the rules defined in robots.txt files - although, their interpretation of those rules might differ (e.g. [rules for googlebot ](https://developers.google.com/search/reference/robots_txt)). 

As the name of the files already suggests robots.txt files are plain text and always found at the root of a domain. The syntax of the files in essence follows a `fieldname: value` scheme with optional preceding `user-agent: ...` lines to indicate the scope of the following rule block. Blocks are separated by blank lines and the omission of a user-agent field (which directly corresponds to the HTTP user-agent field) is seen as referring to all bots. `#` serves to comment lines and parts of lines. Everything after `#` until the end of line is regarded a comment. Possible field names are: user-agent, disallow, allow, crawl-delay, sitemap, and host. 


Let us have an example file to get an idea how a robots.txt file might look like. The file below starts with a comment line followed by a line disallowing access to any content -- everything that is contained in root ("`/`") -- for all bots. The next block concerns GoodBot and NiceBot. Those two get the previous permissions lifted by being disallowed nothing. The third block is for PrettyBot. PrettyBot likes shiny stuff and therefor gets a special permission for everything contained in the "`/shinystuff/`" folder while all other restrictions still hold. In the last block all bots are asked to pause at least 5 seconds between two visits. 


```robots.txt
# this is a comment 
# a made up example of an robots.txt file

Disallow: /

User-agent: GoodBot # another comment
User-agent: NiceBot
Disallow: 

User-agent: PrettyBot
Allow: /shinystuff/

Crawl-Delay: 5
```

For more information have a look at: http://www.robotstxt.org/norobots-rfc.txt, where the robots.txt file 'standard' is described formally. Valuable introductions can be found at http://www.robotstxt.org/robotstxt.html as well as at https://en.wikipedia.org/wiki/Robots_exclusion_standard - of cause. 

# Fast food usage for the uninterested

```{r, message=FALSE}
library(robotstxt)
paths_allowed("http://google.com/")
paths_allowed("http://google.com/search")
```



# Example Usage 

First, let us load the package. In addition we load the dplyr package to be able to use the magrittr pipe operator `%>%` and some easy to read and remember data manipulation functions.

```{r, message=FALSE}
library(robotstxt)
library(dplyr)
```

## object oriented style

The first step is to create an instance of the robotstxt class provided by the package. The instance has to be initiated via providing either domain or the actual text of the robots.txt file. If only the domain is provided, the robots.txt file will be downloaded automatically. Have a look at `?robotstxt` for descriptions of all data fields and methods as well as their parameters. 


```{r, include=FALSE}
rtxt <- 
  robotstxt(
    domain = "wikipedia.org", 
    text   = robotstxt:::rt_get_rtxt("robots_wikipedia.txt")
  )
```

```{r, eval=FALSE}
rtxt <- robotstxt(domain="wikipedia.org")
```

`rtxt` is of class `robotstxt`.

```{r}
class(rtxt)
```

Printing the object lets us glance at all data fields and methods in `rtxt` - we have access to the text as well as all common fields. Non-standard fields are collected in `other`.

```{r}
rtxt
```

Checking permissions works via `rtxt`'s `check` method by providing one or more paths. If no bot name is provided `"*"` - meaning any bot - is assumed. 


```{r}
# checking for access permissions
rtxt$check(paths = c("/","api/"), bot = "*")
rtxt$check(paths = c("/","api/"), bot = "Orthogaffe")
rtxt$check(paths = c("/","api/"), bot = "Mediapartners-Google*  ")
```



## functional style

While working with the robotstxt class is recommended the checking can be done with functions only as well. In the following we (1) download the robots.txt file; (2) parse it and (3) check permissions.

```{r, include=FALSE}
r_text <- robotstxt:::rt_get_rtxt("robots_new_york_times.txt")
```

```{r, eval=FALSE}
r_text        <- get_robotstxt("nytimes.com") 
```

```{r}
r_parsed <- parse_robotstxt(r_text)
r_parsed
```

```{r}
paths_allowed(
  paths  = c("images/","/search"), 
  domain = c("wikipedia.org", "google.com"),
  bot    = "Orthogaffe"
)
``` 






% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{named_list}
\alias{named_list}
\title{make automatically named list}
\usage{
named_list(...)
}
\arguments{
\item{...}{things to be put in list}
}
\description{
make automatically named list
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_robotstxt.R
\name{get_robotstxt}
\alias{get_robotstxt}
\title{downloading robots.txt file}
\usage{
get_robotstxt(
  domain,
  warn = getOption("robotstxt_warn", TRUE),
  force = FALSE,
  user_agent = utils::sessionInfo()$R.version$version.string,
  ssl_verifypeer = c(1, 0),
  encoding = "UTF-8",
  verbose = FALSE,
  rt_request_handler = robotstxt::rt_request_handler,
  rt_robotstxt_http_getter = robotstxt::get_robotstxt_http_get,
  on_server_error = on_server_error_default,
  on_client_error = on_client_error_default,
  on_not_found = on_not_found_default,
  on_redirect = on_redirect_default,
  on_domain_change = on_domain_change_default,
  on_file_type_mismatch = on_file_type_mismatch_default,
  on_suspect_content = on_suspect_content_default
)
}
\arguments{
\item{domain}{domain from which to download robots.txt file}

\item{warn}{warn about being unable to download domain/robots.txt because of}

\item{force}{if TRUE instead of using possible cached results the function
will re-download the robotstxt file HTTP response status 404. If this
happens,}

\item{user_agent}{HTTP user-agent string to be used to retrieve robots.txt
file from domain}

\item{ssl_verifypeer}{analog to CURL option
\url{https://curl.haxx.se/libcurl/c/CURLOPT_SSL_VERIFYPEER.html} -- and
might help with robots.txt file retrieval in some cases}

\item{encoding}{Encoding of the robots.txt file.}

\item{verbose}{make function print out more information}

\item{rt_request_handler}{handler function that handles request according to
the event handlers specified}

\item{rt_robotstxt_http_getter}{function that executes HTTP request}

\item{on_server_error}{request state handler for any 5xx status}

\item{on_client_error}{request state handler for any 4xx HTTP status that is
not 404}

\item{on_not_found}{request state handler for HTTP status 404}

\item{on_redirect}{request state handler for any 3xx HTTP status}

\item{on_domain_change}{request state handler for any 3xx HTTP status where
domain did change as well}

\item{on_file_type_mismatch}{request state handler for content type other
than 'text/plain'}

\item{on_suspect_content}{request state handler for content that seems to be
something else than a robots.txt file (usually a JSON, XML or HTML)}
}
\description{
downloading robots.txt file
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_robotstxt.R
\name{parse_robotstxt}
\alias{parse_robotstxt}
\title{function parsing robots.txt}
\usage{
parse_robotstxt(txt)
}
\arguments{
\item{txt}{content of the robots.txt file}
}
\value{
a named list with useragents, comments, permissions, sitemap
}
\description{
function parsing robots.txt
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paths_allowed.R
\name{paths_allowed}
\alias{paths_allowed}
\title{check if a bot has permissions to access page(s)}
\usage{
paths_allowed(
  paths = "/",
  domain = "auto",
  bot = "*",
  user_agent = utils::sessionInfo()$R.version$version.string,
  check_method = c("spiderbar"),
  warn = getOption("robotstxt_warn", TRUE),
  force = FALSE,
  ssl_verifypeer = c(1, 0),
  use_futures = TRUE,
  robotstxt_list = NULL,
  verbose = FALSE,
  rt_request_handler = robotstxt::rt_request_handler,
  rt_robotstxt_http_getter = robotstxt::get_robotstxt_http_get,
  on_server_error = on_server_error_default,
  on_client_error = on_client_error_default,
  on_not_found = on_not_found_default,
  on_redirect = on_redirect_default,
  on_domain_change = on_domain_change_default,
  on_file_type_mismatch = on_file_type_mismatch_default,
  on_suspect_content = on_suspect_content_default
)
}
\arguments{
\item{paths}{paths for which to check bot's permission, defaults to "/". Please, note that path to a folder should end with a trailing slash ("/").}

\item{domain}{Domain for which paths should be checked. Defaults to "auto".
If set to "auto" function will try to guess the domain by parsing the paths
argument. Note however, that these are educated guesses which might utterly
fail. To be on the safe side, provide appropriate domains manually.}

\item{bot}{name of the bot, defaults to "*"}

\item{user_agent}{HTTP user-agent string to be used to retrieve robots.txt
file from domain}

\item{check_method}{at the moment only kept for backward compatibility reasons - do not use parameter anymore --> will let the function simply use the default}

\item{warn}{suppress warnings}

\item{force}{if TRUE instead of using possible cached results the function
will re-download the robotstxt file HTTP response status 404. If this
happens,}

\item{ssl_verifypeer}{analog to CURL option
\url{https://curl.haxx.se/libcurl/c/CURLOPT_SSL_VERIFYPEER.html} -- and
might help with robots.txt file retrieval in some cases}

\item{use_futures}{Should future::future_lapply be used for possible
parallel/async retrieval or not. Note: check out help
pages and vignettes of package future on how to set up
plans for future execution because the robotstxt package
does not do it on its own.}

\item{robotstxt_list}{either NULL -- the default -- or a list of character
vectors with one vector per path to check}

\item{verbose}{make function print out more information}

\item{rt_request_handler}{handler function that handles request according to
the event handlers specified}

\item{rt_robotstxt_http_getter}{function that executes HTTP request}

\item{on_server_error}{request state handler for any 5xx status}

\item{on_client_error}{request state handler for any 4xx HTTP status that is
not 404}

\item{on_not_found}{request state handler for HTTP status 404}

\item{on_redirect}{request state handler for any 3xx HTTP status}

\item{on_domain_change}{request state handler for any 3xx HTTP status where
domain did change as well}

\item{on_file_type_mismatch}{request state handler for content type other
than 'text/plain'}

\item{on_suspect_content}{request state handler for content that seems to be
something else than a robots.txt file (usually a JSON, XML or HTML)}
}
\description{
check if a bot has permissions to access page(s)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paths_allowed_worker_spiderbar.R
\name{paths_allowed_worker_spiderbar}
\alias{paths_allowed_worker_spiderbar}
\title{paths_allowed_worker spiderbar flavor}
\usage{
paths_allowed_worker_spiderbar(domain, bot, paths, robotstxt_list)
}
\arguments{
\item{domain}{Domain for which paths should be checked. Defaults to "auto".
If set to "auto" function will try to guess the domain by parsing the paths
argument. Note however, that these are educated guesses which might utterly
fail. To be on the safe side, provide appropriate domains manually.}

\item{bot}{name of the bot, defaults to "*"}

\item{paths}{paths for which to check bot's permission, defaults to "/". Please, note that path to a folder should end with a trailing slash ("/").}

\item{robotstxt_list}{either NULL -- the default -- or a list of character
vectors with one vector per path to check}
}
\description{
paths_allowed_worker spiderbar flavor
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robotstxt.R
\name{robotstxt}
\alias{robotstxt}
\title{Generate a representations of a robots.txt file}
\usage{
robotstxt(
  domain = NULL,
  text = NULL,
  user_agent = NULL,
  warn = getOption("robotstxt_warn", TRUE),
  force = FALSE,
  ssl_verifypeer = c(1, 0),
  encoding = "UTF-8",
  verbose = FALSE,
  on_server_error = on_server_error_default,
  on_client_error = on_client_error_default,
  on_not_found = on_not_found_default,
  on_redirect = on_redirect_default,
  on_domain_change = on_domain_change_default,
  on_file_type_mismatch = on_file_type_mismatch_default,
  on_suspect_content = on_suspect_content_default
)
}
\arguments{
\item{domain}{Domain for which to generate a representation. If text equals to NULL,
the function will download the file from server - the default.}

\item{text}{If automatic download of the robots.txt is not preferred, the text can be
supplied directly.}

\item{user_agent}{HTTP user-agent string to be used to retrieve robots.txt
file from domain}

\item{warn}{warn about being unable to download domain/robots.txt because of}

\item{force}{if TRUE instead of using possible cached results the function
will re-download the robotstxt file HTTP response status 404. If this
happens,}

\item{ssl_verifypeer}{analog to CURL option
\url{https://curl.haxx.se/libcurl/c/CURLOPT_SSL_VERIFYPEER.html} -- and
might help with robots.txt file retrieval in some cases}

\item{encoding}{Encoding of the robots.txt file.}

\item{verbose}{make function print out more information}

\item{on_server_error}{request state handler for any 5xx status}

\item{on_client_error}{request state handler for any 4xx HTTP status that is
not 404}

\item{on_not_found}{request state handler for HTTP status 404}

\item{on_redirect}{request state handler for any 3xx HTTP status}

\item{on_domain_change}{request state handler for any 3xx HTTP status where
domain did change as well}

\item{on_file_type_mismatch}{request state handler for content type other
than 'text/plain'}

\item{on_suspect_content}{request state handler for content that seems to be
something else than a robots.txt file (usually a JSON, XML or HTML)}
}
\value{
Object (list) of class robotstxt with parsed data from a
  robots.txt (domain, text, bots, permissions, host, sitemap, other) and one
  function to (check()) to check resource permissions.
}
\description{
The function generates a list that entails data resulting from parsing a robots.txt file
as well as a function called check that enables to ask the representation if bot (or
particular bots) are allowed to access a resource on the domain.
}
\section{Fields}{

\describe{
\item{\code{domain}}{character vector holding domain name for which the robots.txt
file is valid; will be set to NA if not supplied on initialization}

\item{\code{text}}{character vector of text of robots.txt file; either supplied on
initialization or automatically downloaded from domain supplied on
initialization}

\item{\code{bots}}{character vector of bot names mentioned in robots.txt}

\item{\code{permissions}}{data.frame of bot permissions found in robots.txt file}

\item{\code{host}}{data.frame of host fields found in robots.txt file}

\item{\code{sitemap}}{data.frame of sitemap fields found in robots.txt file}

\item{\code{other}}{data.frame of other - none of the above - fields found in
robots.txt file}

\item{\code{check()}}{Method to check for bot permissions. Defaults to the
domains root and no bot in particular. check() has two arguments:
paths and bot. The first is for supplying the paths for which to check
permissions and the latter to put in the name of the bot.
Please, note that path to a folder should end with a trailing slash ("/").}
}}

\examples{
\dontrun{
rt <- robotstxt(domain="google.com")
rt$bots
rt$permissions
rt$check( paths = c("/", "forbidden"), bot="*")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt_request_handler.R,
%   R/rt_request_handler_defaults.R
\docType{data}
\name{rt_request_handler}
\alias{rt_request_handler}
\alias{on_server_error_default}
\alias{on_client_error_default}
\alias{on_not_found_default}
\alias{on_redirect_default}
\alias{on_domain_change_default}
\alias{on_sub_domain_change_default}
\alias{on_file_type_mismatch_default}
\alias{on_suspect_content_default}
\title{rt_request_handler}
\format{
An object of class \code{list} of length 4.

An object of class \code{list} of length 4.

An object of class \code{list} of length 4.

An object of class \code{list} of length 2.

An object of class \code{list} of length 3.

An object of class \code{list} of length 2.

An object of class \code{list} of length 4.

An object of class \code{list} of length 4.
}
\usage{
rt_request_handler(
  request,
  on_server_error = on_server_error_default,
  on_client_error = on_client_error_default,
  on_not_found = on_not_found_default,
  on_redirect = on_redirect_default,
  on_domain_change = on_domain_change_default,
  on_sub_domain_change = on_sub_domain_change_default,
  on_file_type_mismatch = on_file_type_mismatch_default,
  on_suspect_content = on_suspect_content_default,
  warn = TRUE,
  encoding = "UTF-8"
)

on_server_error_default

on_client_error_default

on_not_found_default

on_redirect_default

on_domain_change_default

on_sub_domain_change_default

on_file_type_mismatch_default

on_suspect_content_default
}
\arguments{
\item{request}{result of an HTTP request (e.g. httr::GET())}

\item{on_server_error}{request state handler for any 5xx status}

\item{on_client_error}{request state handler for any 4xx HTTP status that is
not 404}

\item{on_not_found}{request state handler for HTTP status 404}

\item{on_redirect}{request state handler for any 3xx HTTP status}

\item{on_domain_change}{request state handler for any 3xx HTTP status where
domain did change as well}

\item{on_sub_domain_change}{request state handler for any 3xx HTTP status where
domain did change but only to www-sub_domain}

\item{on_file_type_mismatch}{request state handler for content type other
than 'text/plain'}

\item{on_suspect_content}{request state handler for content that seems to be
something else than a robots.txt file (usually a JSON, XML or HTML)}

\item{warn}{suppress warnings}

\item{encoding}{The text encoding to assume if no encoding is provided in the
headers of the response}
}
\value{
a list with three items following the following schema: \cr \code{
  list( rtxt = "", problems = list( "redirect" = list( status_code = 301 ),
  "domain" = list(from_url = "...", to_url = "...") ) ) }
}
\description{
A helper function for get_robotstxt() that will extract the robots.txt file
from the HTTP request result object. furthermore it will inform
get_robotstxt() if the request should be cached and which problems occured.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_domain.R
\name{remove_domain}
\alias{remove_domain}
\title{function to remove domain from path}
\usage{
remove_domain(x)
}
\arguments{
\item{x}{path aka URL from which to first infer domain and then remove it}
}
\description{
function to remove domain from path
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_robotstxt.R
\name{print.robotstxt}
\alias{print.robotstxt}
\title{printing robotstxt}
\usage{
\method{print}{robotstxt}(x, ...)
}
\arguments{
\item{x}{robotstxt instance to be printed}

\item{...}{goes down the sink}
}
\description{
printing robotstxt
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{rt_list_rtxt}
\alias{rt_list_rtxt}
\title{list robots.txt files saved along with the package}
\usage{
rt_list_rtxt()
}
\description{
list robots.txt files saved along with the package:
these functions ar very handy for testing (not used otherwise)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http_was_redirected.R
\name{http_was_redirected}
\alias{http_was_redirected}
\title{http_was_redirected}
\usage{
http_was_redirected(response)
}
\arguments{
\item{response}{an httr response object, e.g. from a call to httr::GET()}
}
\value{
logical of length 1 indicating whether or not any redirect happened
  during the HTTP request
}
\description{
http_was_redirected
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_robotstxt_http_get.R
\docType{data}
\name{rt_last_http}
\alias{rt_last_http}
\alias{get_robotstxt_http_get}
\title{storage for http request response objects}
\format{
An object of class \code{environment} of length 1.
}
\usage{
rt_last_http

get_robotstxt_http_get(
  domain,
  user_agent = utils::sessionInfo()$R.version$version.string,
  ssl_verifypeer = 1
)
}
\arguments{
\item{domain}{the domain to get tobots.txt. file for}

\item{user_agent}{the user agent to use for HTTP request header}

\item{ssl_verifypeer}{analog to CURL option
\url{https://curl.haxx.se/libcurl/c/CURLOPT_SSL_VERIFYPEER.html}
-- and might help with robots.txt file retrieval in some cases}
}
\description{
storage for http request response objects

get_robotstxt() worker function to execute HTTP request
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt_get_comments.R
\name{rt_get_comments}
\alias{rt_get_comments}
\title{extracting comments from robots.txt}
\usage{
rt_get_comments(txt)
}
\arguments{
\item{txt}{content of the robots.txt file}
}
\description{
extracting comments from robots.txt
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http_domain_changed.R
\name{http_domain_changed}
\alias{http_domain_changed}
\title{http_domain_changed}
\usage{
http_domain_changed(response)
}
\arguments{
\item{response}{an httr response object, e.g. from a call to httr::GET()}
}
\value{
logical of length 1 indicating whether or not any domain change
    happened during the HTTP request
}
\description{
http_domain_changed
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{re-export magrittr pipe operator}
\description{
re-export magrittr pipe operator
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_robotstxts.R
\name{get_robotstxts}
\alias{get_robotstxts}
\title{function to get multiple robotstxt files}
\usage{
get_robotstxts(
  domain,
  warn = TRUE,
  force = FALSE,
  user_agent = utils::sessionInfo()$R.version$version.string,
  ssl_verifypeer = c(1, 0),
  use_futures = FALSE,
  verbose = FALSE,
  rt_request_handler = robotstxt::rt_request_handler,
  rt_robotstxt_http_getter = robotstxt::get_robotstxt_http_get,
  on_server_error = on_server_error_default,
  on_client_error = on_client_error_default,
  on_not_found = on_not_found_default,
  on_redirect = on_redirect_default,
  on_domain_change = on_domain_change_default,
  on_file_type_mismatch = on_file_type_mismatch_default,
  on_suspect_content = on_suspect_content_default
)
}
\arguments{
\item{domain}{domain from which to download robots.txt file}

\item{warn}{warn about being unable to download domain/robots.txt because of}

\item{force}{if TRUE instead of using possible cached results the function
will re-download the robotstxt file HTTP response status 404. If this
happens,}

\item{user_agent}{HTTP user-agent string to be used to retrieve robots.txt
file from domain}

\item{ssl_verifypeer}{analog to CURL option
\url{https://curl.haxx.se/libcurl/c/CURLOPT_SSL_VERIFYPEER.html}
-- and might help with robots.txt file retrieval in some cases}

\item{use_futures}{Should future::future_lapply be used for possible
parallel/async retrieval or not. Note: check out help
pages and vignettes of package future on how to set up
plans for future execution because the robotstxt package
does not do it on its own.}

\item{verbose}{make function print out more information}

\item{rt_request_handler}{handler function that handles request according to
the event handlers specified}

\item{rt_robotstxt_http_getter}{function that executes HTTP request}

\item{on_server_error}{request state handler for any 5xx status}

\item{on_client_error}{request state handler for any 4xx HTTP status that is
not 404}

\item{on_not_found}{request state handler for HTTP status 404}

\item{on_redirect}{request state handler for any 3xx HTTP status}

\item{on_domain_change}{request state handler for any 3xx HTTP status where
domain did change as well}

\item{on_file_type_mismatch}{request state handler for content type other
than 'text/plain'}

\item{on_suspect_content}{request state handler for content that seems to be
something else than a robots.txt file (usually a JSON, XML or HTML)}
}
\description{
function to get multiple robotstxt files
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{rt_get_rtxt}
\alias{rt_get_rtxt}
\title{load robots.txt files saved along with the package}
\usage{
rt_get_rtxt(name = sample(rt_list_rtxt(), 1))
}
\arguments{
\item{name}{name of the robots.txt files, defaults to a random drawn file ;-)}
}
\description{
load robots.txt files saved along with the package:
these functions are very handy for testing (not used otherwise)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_suspect_robotstxt.R
\name{is_suspect_robotstxt}
\alias{is_suspect_robotstxt}
\title{is_suspect_robotstxt}
\usage{
is_suspect_robotstxt(text)
}
\arguments{
\item{text}{content of a robots.txt file provides as character vector}
}
\description{
function that checks if file is valid / parsable robots.txt file
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http_subdomain_changed.R
\name{http_subdomain_changed}
\alias{http_subdomain_changed}
\title{http_subdomain_changed}
\usage{
http_subdomain_changed(response)
}
\arguments{
\item{response}{an httr response object, e.g. from a call to httr::GET()}
}
\value{
logical of length 1 indicating whether or not any domain change
    happened during the HTTP request
}
\description{
http_subdomain_changed
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt_cache.R
\docType{data}
\name{rt_cache}
\alias{rt_cache}
\title{get_robotstxt() cache}
\format{
An object of class \code{environment} of length 0.
}
\usage{
rt_cache
}
\description{
get_robotstxt() cache
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_valid_robotstxt.R
\name{is_valid_robotstxt}
\alias{is_valid_robotstxt}
\title{function that checks if file is valid / parsable robots.txt file}
\usage{
is_valid_robotstxt(text, check_strickt_ascii = FALSE)
}
\arguments{
\item{text}{content of a robots.txt file provides as character vector}

\item{check_strickt_ascii}{whether or not to check if content does adhere to the specification of RFC to use plain text aka ASCII}
}
\description{
function that checks if file is valid / parsable robots.txt file
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sanitize_path.R
\name{sanitize_path}
\alias{sanitize_path}
\title{making paths uniform}
\usage{
sanitize_path(path)
}
\arguments{
\item{path}{path to be sanitized}
}
\value{
sanitized path
}
\description{
making paths uniform
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fix_url.R
\name{fix_url}
\alias{fix_url}
\title{fix_url}
\usage{
fix_url(url)
}
\arguments{
\item{url}{a character string containing a single URL}
}
\description{
fix_url
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_robotstxt_text.R
\name{print.robotstxt_text}
\alias{print.robotstxt_text}
\title{printing robotstxt_text}
\usage{
\method{print}{robotstxt_text}(x, ...)
}
\arguments{
\item{x}{character vector aka robotstxt$text to be printed}

\item{...}{goes down the sink}
}
\description{
printing robotstxt_text
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt_get_fields.R
\name{rt_get_fields}
\alias{rt_get_fields}
\title{extracting permissions from robots.txt}
\usage{
rt_get_fields(txt, regex = "", invert = FALSE)
}
\arguments{
\item{txt}{content of the robots.txt file}

\item{regex}{regular expression specify field}

\item{invert}{invert selection made via regex?}
}
\description{
extracting permissions from robots.txt
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_merge.R
\name{list_merge}
\alias{list_merge}
\title{Merge a number of named lists in sequential order}
\usage{
list_merge(...)
}
\arguments{
\item{...}{named lists}
}
\description{
Merge a number of named lists in sequential order
}
\details{
List merging is usually useful in the merging of program
settings or configuraion with multiple versions across time,
or multiple administrative levels. For example, a program
settings may have an initial version in which most keys are
defined and specified. In later versions, partial modifications
are recorded. In this case, list merging can be useful to merge
all versions of settings in release order of these versions. The
result is an fully updated settings with all later modifications
applied.
}
\author{
Kun Ren <mail@renkun.me>


The function merges a number of lists in sequential order
by \code{modifyList}, that is, the later list always
modifies the former list and form a merged list, and the
resulted list is again being merged with the next list.
The process is repeated until all lists in \code{...} or
\code{list} are exausted.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/null_to_default.R
\name{null_to_defeault}
\alias{null_to_defeault}
\title{null_to_defeault}
\usage{
null_to_defeault(x, d)
}
\arguments{
\item{x}{value to check and return}

\item{d}{value to return in case x is NULL}
}
\description{
null_to_defeault
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt_get_useragent.R
\name{rt_get_useragent}
\alias{rt_get_useragent}
\title{extracting HTTP useragents from robots.txt}
\usage{
rt_get_useragent(txt)
}
\arguments{
\item{txt}{content of the robots.txt file}
}
\description{
extracting HTTP useragents from robots.txt
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/guess_domain.R
\name{guess_domain}
\alias{guess_domain}
\title{function guessing domain from path}
\usage{
guess_domain(x)
}
\arguments{
\item{x}{path aka URL from which to infer domain}
}
\description{
function guessing domain from path
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_url.R
\name{parse_url}
\alias{parse_url}
\title{parse_url}
\usage{
parse_url(url)
}
\arguments{
\item{url}{url to parse into its components}
}
\value{
data.frame with columns protocol, domain, path
}
\description{
parse_url
}
\examples{

\dontrun{
url <-
c(
  "google.com",
  "google.com/",
  "www.google.com",
  "http://google.com",
  "https://google.com",
  "sub.domain.whatever.de"
  "s-u-b.dom-ain.what-ever.de"
)

parse_url(url)
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt_get_fields_worker.R
\name{rt_get_fields_worker}
\alias{rt_get_fields_worker}
\title{extracting robotstxt fields}
\usage{
rt_get_fields_worker(txt, type = "all", regex = NULL, invert = FALSE)
}
\arguments{
\item{txt}{content of the robots.txt file}

\item{type}{name or names of the fields to be returned, defaults to all
fields}

\item{regex}{subsetting field names via regular expressions}

\item{invert}{field selection}
}
\description{
extracting robotstxt fields
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_list.R
\name{as.list.robotstxt_text}
\alias{as.list.robotstxt_text}
\title{Method as.list() for class robotstxt_text}
\usage{
\method{as.list}{robotstxt_text}(x, ...)
}
\arguments{
\item{x}{class robotstxt_text object to be transformed into list}

\item{...}{further arguments (inherited from \code{base::as.list()})}
}
\description{
Method as.list() for class robotstxt_text
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_handler_handler.R
\name{request_handler_handler}
\alias{request_handler_handler}
\title{request_handler_handler}
\usage{
request_handler_handler(request, handler, res, info = TRUE, warn = TRUE)
}
\arguments{
\item{request}{the request object returned by call to httr::GET()}

\item{handler}{the handler either a character string entailing various options or a function producing a specific list, see return.}

\item{res}{a list a list with elements '[handler names], ...', 'rtxt', and 'cache'}

\item{info}{info to add to problems list}

\item{warn}{if FALSE warnings and messages are suppressed}
}
\value{
a list with elements '[handler name]', 'rtxt', and 'cache'
}
\description{
Helper function to handle robotstxt handlers.
}
