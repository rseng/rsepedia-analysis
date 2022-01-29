
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



