vcr
===



<!-- README.md is generated from README.Rmd. Please edit that file -->

[![cran checks](https://cranchecks.info/badges/worst/vcr)](https://cranchecks.info/pkgs/vcr)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/vcr/workflows/R-check/badge.svg)](https://github.com/ropensci/vcr/actions/)
[![codecov](https://codecov.io/gh/ropensci/vcr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/vcr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/vcr)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/vcr)](https://cran.r-project.org/package=vcr)


Easier HTTP testing! Record HTTP requests and responses on disk and replay them for the unit tests of your R package, to make them independent from any connection, faster, and more complete. An R port of the Ruby gem [vcr](https://github.com/vcr/vcr)

## Elevator pitch


* **Setup vcr for your package with `vcr::use_vcr()`**
* Tweak the configuration to protect your secrets
* **Sprinkle your tests with `vcr::use_cassette()` to save HTTP interactions to disk in "cassettes" files**
* If you want to test for package behavior when the API returns e.g. a 404 or 503 code, edit the cassettes, or use [webmockr](https://docs.ropensci.org/webmockr/)

Now your tests can work without any internet connection!

[Demo of adding vcr testing to an R package](https://github.com/maelle/exemplighratia/pull/2/files)

## Installation

CRAN version:


```r
install.packages("vcr")
```

Development version:


```r
remotes::install_github("ropensci/vcr")
```


```r
library("vcr")
library("crul")
```


## Docs

Check out the [HTTP testing book](https://books.ropensci.org/http-testing) and the [vcr vignettes](https://docs.ropensci.org/vcr/articles/).

## Supported HTTP libraries

* [crul](https://docs.ropensci.org/crul)
* [httr](https://httr.r-lib.org/)

## Getting Started


The docs assume you are using testthat for your unit tests.

### `use_vcr`

You can then set up your package to use `vcr` with:

```r
vcr::use_vcr()
```

This will:

* put `vcr` into the `DESCRIPTION`
* check that `testthat` is setup
* setup `testthat` if not
* set the recorded cassettes to be saved in and sourced from `tests/fixtures`
* setup a config file for `vcr`
* add an example test file for `vcr`
* make a `.gitattributes` file with settings for `vcr` 
* make a `./tests/testthat/setup-vcr.R` file

What you will see in the R console:

```
◉ Using package: vcr.example  
◉ assuming fixtures at: tests/fixtures  
✓ Adding vcr to Suggests field in DESCRIPTION  
✓ Creating directory: ./tests/testthat  
◉ Looking for testthat.R file or similar  
✓ tests/testthat.R: added  
✓ Adding vcr config to tests/testthat/setup-vcr.example.R  
✓ Adding example test file tests/testthat/test-vcr_example.R  
✓ .gitattributes: added  
◉ Learn more about `vcr`: https://books.ropensci.org/http-testing
```

### Protecting secrets

Secrets often turn up in API work. A common example is an API key. 
`vcr` saves responses from APIs as YAML files, and this will include your secrets unless you indicate to `vcr` what they are and how to protect them.
The `vcr_configure` function has the `filter_sensitive_data` argument function for just this situation. 
The `filter_sensitive_data` argument takes a named list where the _name_ of the list is the string that will be used in the recorded cassettes _instead of_ the secret, which is the list _item_. 
`vcr` will manage the replacement of that for you, so all you need to do is to edit your `setup-vcr.R` file like this:

```r
library("vcr") # *Required* as vcr is set up on loading
invisible(vcr::vcr_configure(
  dir = "../fixtures"
))
vcr::check_cassette_names()
```

Use the `filter_sensitive_data` argument in the `vcr_configure` function to show `vcr` how to keep your secret. The best way to store secret information is to have it in a `.Renviron` file. Assuming that that is already in place, supply a named list to the `filter_sensitive_data` argument.

```r
library("vcr")
invisible(vcr::vcr_configure(
  filter_sensitive_data = list("<<<my_api_key>>>" = Sys.getenv('APIKEY')),  # add this
  dir = "../fixtures"
))
vcr::check_cassette_names()
```

Notice we wrote `Sys.getenv('APIKEY')` and not the API key directly, otherwise you'd have written your API key to a file that might end up in a public repo.

The will get your secret information from the environment, and make sure that whenever `vcr` records a new cassette, it will replace the secret information with `<<<my_api_key>>>`. You can find out more about this in the [HTTP testing book](https://books.ropensci.org/http-testing/) chapter on security.

The addition of the line above will instruct `vcr` to replace any string in cassettes it records that are equivalent to your string which is stored as the `APIKEY` environmental variable with the masking string `<<<my_api_key>>>`. In practice, you might get a `YAML` that looks a little like this:

```yaml
http_interactions:
- request:
    method: post
    ...
    headers:
      Accept: application/json, text/xml, application/xml, */*
      Content-Type: application/json
      api-key: <<<my_api_key>>>
    ...
```
Here, my `APIKEY` environmental variable would have been stored as the `api-key` value, but `vcr` has realised this and recorded the string `<<<my_api_key>>>` instead.

Once the cassette is recorded, `vcr` no longer needs the API key as no real requests will be made.
Furthermore, as by default requests matching does not include the API key, things will work.

**Now, how to ensure tests work in the absence of a real API key?**

E.g. to have tests pass on continuous integration for external pull requests to your code repository.

* vcr does not need an actual API key for requests once the cassettes are created, as no real requests will be made.
* you still need to fool your _package_ into believing there is an API key as it will construct requests with it. So add the following lines to a testthat setup file (e.g. `tests/testthat/setup-vcr.R`)

```r
if (!nzchar(Sys.getenv("APIKEY"))) {
  Sys.setenv("APIKEY" = "foobar")
}
```

#### Using an `.Renviron`

A simple way to manage local environmental variables is to use an [`.Renviron` file](https://rstats.wtf/r-startup.html#renviron). 
Your `.Renviron` file might look like this:

```sh
APIKEY="mytotallysecretkey"
```

You can have this set at a project or user level, and `usethis` has the [`usethis::edit_r_environ()`](https://usethis.r-lib.org/reference/edit.html) function to help edit the file.

## Usage





### In tests

In your tests, for whichever tests you want to use `vcr`, wrap them in a `vcr::use_cassette()` call like:

```r
library(testthat)
vcr::use_cassette("rl_citation", {
  test_that("my test", {
    aa <- rl_citation()

    expect_is(aa, "character")
    expect_match(aa, "IUCN")
    expect_match(aa, "www.iucnredlist.org")
  })
})
```

OR put the `vcr::use_cassette()` block on the inside, but put `testthat` expectations outside of
the `vcr::use_cassette()` block:

```r
library(testthat)
test_that("my test", {
  vcr::use_cassette("rl_citation", {
    aa <- rl_citation()
  })

  expect_is(aa, "character")
  expect_match(aa, "IUCN")
  expect_match(aa, "www.iucnredlist.org")
})
```

Don't wrap the `use_cassette()` block inside your  `test_that()` block with `testthat` expectations inside the `use_cassette()` block, as you'll only get the line number that the `use_cassette()` block starts on on failures.

The first time you run the tests, a "cassette" i.e. a file with recorded HTTP interactions, is created at `tests/fixtures/rl_citation.yml`.
The times after that, the cassette will be used.
If you change your code and more HTTP interactions are needed in the code wrapped by `vcr::use_cassette("rl_citation"`, delete `tests/fixtures/rl_citation.yml` and run the tests again for re-recording the cassette.

### Outside of tests

If you want to get a feel for how vcr works, although you don't need too.




```r
library(vcr)
library(crul)

cli <- crul::HttpClient$new(url = "https://eu.httpbin.org")
system.time(
  use_cassette(name = "helloworld", {
    cli$get("get")
  })
)
```

The request gets recorded, and all subsequent requests of the same form used the cached HTTP response, and so are much faster


```r
system.time(
  use_cassette(name = "helloworld", {
    cli$get("get")
  })
)
```



Importantly, your unit test deals with the same inputs and the same outputs - but behind the scenes you use a cached HTTP response - thus, your tests run faster.

The cached response looks something like (condensed for brevity):

```yaml
http_interactions:
- request:
    method: get
    uri: https://eu.httpbin.org/get
    body:
      encoding: ''
      string: ''
    headers:
      User-Agent: libcurl/7.54.0 r-curl/3.2 crul/0.5.2
  response:
    status:
      status_code: '200'
      message: OK
      explanation: Request fulfilled, document follows
    headers:
      status: HTTP/1.1 200 OK
      connection: keep-alive
    body:
      encoding: UTF-8
      string: "{\n  \"args\": {}, \n  \"headers\": {\n    \"Accept\": \"application/json,
        text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\",
        \n    \"Connection\": \"close\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\":
        \"libcurl/7.54.0 r-curl/3.2 crul/0.5.2\"\n  }, \n  \"origin\": \"111.222.333.444\",
        \n  \"url\": \"https://eu.httpbin.org/get\"\n}\n"
  recorded_at: 2018-04-03 22:55:02 GMT
  recorded_with: vcr/0.1.0, webmockr/0.2.4, crul/0.5.2
```

All components of both the request and response are preserved, so that the HTTP client (in this case `crul`) can reconstruct its own response just as it would if it wasn't using `vcr`.


### Less basic usage

For tweaking things to your needs, make sure to read the docs about [configuration](https://docs.ropensci.org/vcr/articles/configuration.html) (e.g., where are the fixtures saved? can they be re-recorded automatically regulary?) and [request matching](https://docs.ropensci.org/vcr/articles/request_matching.html) (how does vcr match a request to a recorded interaction?)


## Terminology


* _vcr_: the name comes from the idea that we want to record something and play it back later, like a vcr
* _cassette_: A _thing_ to record HTTP interactions to. Right now the only option is the file system (writing to files), but in the future could be other things, e.g. a key-value store like Redis
* _fixture_: A fixture is something used to consistently test a piece of software. In this case, a cassette (just defined above) is a fixture - used in unit tests. If you use our setup function `vcr_setup()` the default directory created to hold cassettes is called `fixtures/` as a signal as to what the folder contains.
* Persisters: how to save requests - currently only option is the file system
* _serialize_: translating data into a format that can be stored; here, translate HTTP request and response data into a representation on disk to read back later
* Serializers: how to serialize the HTTP response - currently only option is YAML; other options in the future could include e.g. JSON
* _insert cassette_: create a cassette (all HTTP interactions will be recorded to this cassette)
* _eject cassette_: eject the cassette (no longer recording to that cassette)
* _replay_: refers to using a cached result of an http request that was recorded earlier

## Workflows



### vcr for tests

* See [usage section](#usage)

* When running tests or checks of your whole package, note that some users have found different results with
`devtools::check()` vs. `devtools::test()`. It's not clear why this would make a difference. Do let us know
if you run into this problem.

### vcr in your R project

You can use `vcr` in any R project as well.

* Load `vcr` in your project
* Similar to the above example, use `use_cassette` to run code that does HTTP requests.
* The first time a real request is done, and after that the cached response will be used.

### How it works in lots of detail

See the [vignette about internals](https://docs.ropensci.org/vcr/articles/internals.html)

### Just want to mock and not store on disk?

You're looking for [webmockr][], that vcr itself uses. 
`webmockr` only matches requests based on criteria you choose, but does not cache HTTP interactions to disk as `vcr` does.

## Configuration


See also the [configuration vignette](https://docs.ropensci.org/vcr/articles/configuration.html).

We set the following defaults:

* dir = `"."`
 * record = `"once"`
 * match_requests_on = `"c("method", "uri")"`
 * allow_unused_http_interactions = `TRUE`
 * serialize_with = `"yaml"`
 * json_pretty = `FALSE`
 * persist_with = `"FileSystem"`
 * ignore_hosts = `NULL`
 * ignore_localhost = `FALSE`
 * ignore_request = `NULL`
 * uri_parser = `"crul::url_parse"`
 * preserve_exact_body_bytes = `FALSE`
 * turned_off = `FALSE`
 * re_record_interval = `NULL`
 * clean_outdated_http_interactions = `FALSE`
 * allow_http_connections_when_no_cassette = `FALSE`
 * cassettes = `list()`
 * linked_context = `NULL`
 * log = `FALSE`
 * log_opts = `list(file = "vcr.log", log_prefix = "Cassette", date = TRUE)`
 * filter_sensitive_data = `NULL`
 * filter_request_headers = `NULL`
 * filter_response_headers = `NULL`
 * write_disk_path = `NULL`
 * verbose_errors = `FALSE`


You can get the defaults programmatically with

```r
vcr_config_defaults()
```

You can change all the above defaults with `vcr_configure()`:

```r
vcr_configure()
```

Calling `vcr_configuration()` gives you some of the more important configuration parameters in a nice tidy print out


```r
vcr_configuration()
#> <vcr configuration>
#>   Cassette Dir: .
#>   Record: once
#>   Serialize with: yaml
#>   URI Parser: crul::url_parse
#>   Match Requests on: method, uri
#>   Preserve Bytes?: FALSE
#>   Logging?: FALSE
#>   ignored hosts: 
#>   ignore localhost?: FALSE
#>   Write disk path:
```


<!-- `use_cassette()` is an easier approach. An alternative is to use
`insert_cassett()` + `eject_cassette()`.

`use_cassette()` does both insert and eject operations for you, but
you can instead do them manually by using the above functions. You do have
to eject the cassette after using insert. -->

For more details refer to the [configuration vignette](https://docs.ropensci.org/vcr/articles/configuration.html)

## Matching/Matchers


`vcr` looks for similarity in your HTTP requests to cached requests. You
can set what is examined about the request with one or more of the
following options:

* `body`
* `headers`
* `host`
* `method`
* `path`
* `query`
* `uri`

By default, we use `method` (HTTP method, e.g., `GET`) and `uri` (test for exact match against URI, e.g., `http://foo.com`).

You can set your own options by tweaking the `match_requests_on` parameter:




```r
use_cassette(name = "one", {
    cli$post("post", body = list(a = 5))
  },
  match_requests_on = c('method', 'headers', 'body')
)
```

For more details refer to the [request matching vignette](https://docs.ropensci.org/vcr/articles/request_matching.html).

## vcr in other languages

The canonical `vcr` (in Ruby) lists ports in other languages at <https://github.com/vcr/vcr>

## Note about missing features


There's a number of features in this package that are not yet supported, but for which their parameters are found in the package. 

We've tried to make sure the parameters that are ignored are marked as such. Keep an eye out for package updates for changes in these parameters, and/or let us know you want it and we can move it up in the priority list.

## Example packages using vcr

* [rgbif][]
* [rredlist][]
* [bold][]
* [wikitaxa][]
* [worrms][]
* [microdemic][]
* [zbank][]
* [rplos][]
* [ritis][]

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Aaron Wolen](https://github.com/aaronwolen)
* [Maëlle Salmon](https://github.com/maelle)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/vcr/issues)
* License: MIT
* Get citation information for `vcr` in R doing `citation(package = 'vcr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[webmockr]: https://docs.ropensci.org/webmockr
[crul]: https://docs.ropensci.org/crul
[rgbif]: https://github.com/ropensci/rgbif
[rredlist]: https://github.com/ropensci/rredlist
[bold]: https://github.com/ropensci/bold
[wikitaxa]: https://github.com/ropensci/wikitaxa
[worrms]: https://github.com/ropensci/worrms
[microdemic]: https://github.com/ropensci/microdemic
[zbank]: https://github.com/ropenscilabs/zbank
[rplos]: https://github.com/ropensci/rplos
[ritis]: https://github.com/ropensci/ritis
vcr 1.0.2
=========

### BUG FIXES

* fix to `vcr_test_path()` to find root package path correctly (#235) (#236)

vcr 1.0.0
=========

### NEW FEATURES

* `check_cassette_names()` gains `allowed_duplicates` parameter to allow duplicate cassette names; we typically advise users not to use duplicate cassette names, but there are cases where you may want to share cassettes across tests (#227)
* `vcr_configure()` gains `filter_query_parameters` parameter for filtering out query parameters so they don't show up in the recorded request on disk (#212)
* `use_vcr()`: now sets a mimimum vcr version, which is usually the latest (stable) version on CRAN. You can of course easily remove or change the version requirement yourself after running it (#214)
* `vcr_configure()` gains `warn_on_empty_cassette` parameter: Should a warning be thrown when an empty cassette is detected? Empty cassettes are cleaned up (deleted) either way (#224) thanks @llrs and @dpprdan
* `vcr_configure()` gains `quiet` parameter: suppress any messages from both vcr and webmockr (#226) (#25)
* `vcr_configure()` gains new option `filter_sensitive_data_regex`; now `filter_sensitive_data` is for fixed string matching, while `filter_sensitive_data_regex` is for regex based matching (#222) thanks @tomsing1 for reporting
* gains package import `rprojroot`

### MINOR IMPROVEMENTS

* `filter_sensitive_data` option now strips leading and trailing single and double quotes from strings before being used IN CASE a user accidentally quotes a secret - logic being that even though a secret may have a single or double quote in it, its very unlikely that it would have both a leading and trailing quote (single or double) (#221)

### Documentation

* new vignette explaining the design of the vcr package (also can be found in the HTTP Testing book) (#232) (#233)
* no user facing change - but vignettes moved into man/rmdhunks so that they can be pulled into the HTTP Testing book easily (#209) (#216)
* fix in configuration vignette to clarify a `filter_request_headers` example  (#215) thanks @maelle
* docs update (#33) (#217)

### BUG FIXES

* `filter_request_headers` was unfortunately adding a request header to the request written to disk when the header did not exist; now fixed (#213)
* bug in internal function `is_base64()`; `strsplit()` needed `useBytes=TRUE` (#219)
* `filter_sensitive_data` was not working when strings contained regex characters; fixed, and see also above new config variable for regex specific filtering  (#222) thanks @tomsing1 for reporting
* `vcr_test_path()` should now correctly set paths (#225) (#228) (#229) (#230)


vcr 0.6.0
=========

### NEW FEATURES

* We have a new vcr contributor! @maelle (#198)
* Gains a new serializer: JSON. You can use this serializer by setting globally `vcr_configure(serialize_with="json")` or per cassette `use_cassette(..., serialize_with="json")`. The JSON serializer uses `jsonlite` under the hood. Note that we by default do not write JSON to disk preserving newlines; that is the JSON is all on one line. You can use pretty printing by setting `json_pretty` in `vcr_configure()`. As part of this change, factored out new R6 class `Serializer` from which both JSON and YAML serializers inherit (#32)
* Gains two new configuration options for managing secrets: `filter_request_headers` and `filter_response_headers`. These are implemented differently than `filter_sensitive_data`. The two new filters do simple value replacement or complete removal of request or response headers, whereas `filter_sensitive_data` uses regex to replace strings anywhere in the stored request/response. See the "Configure vcr" vignette for details (#182)
* request matching: `host` and `path` now work (#177) (see also #70)
* In previous versions of vcr the `insert_cassette()`/`eject_cassette()` workflow did not work because the webmockr triggers required only worked when using `use_cassette()`. This has been fixed now so you can use `use_cassette()`, passing a code block to it, or run `insert_cassette()` then run any code, then when finished run `eject_cassette()`.  (#24) thanks @Robsteranium for the nudge, may not have fixed this without it
* improve debugging experience: new vignette "Debugging your tests that use vcr", including new function `vcr_test_path()` - which is now used in `use_vcr()` so that the correct path to tests is used when running tests both interactively and non-interactively (#192) (#193)
* Dependencies: dropped `lazyeval` from Imports; `withr` added to Suggests; minimum `webmockr` version now `v0.7.4`
* In README, point to rOpenSci code of conduct rather than file in repo
* Gains function `skip_if_vcr_off()` to use in tests to skip a test if vcr is turned off (#191) (#195)

### MINOR IMPROVEMENTS

* slight factor out of some code in YAML serializer to use elsewhere (#203) (#204)
* serializers: drop `$deserialize_string()` method as was not used - rename `$deserialize_path()` method to just `$deserialize()` (#189)
* serializers: with the new JSON serializer, documentation added to `?vcr_configure` and `?use_cassette` stating that you can have different cassettes with the same name as long as they use different serializers (and then have different file extensions). if you want to change serializers but do not want to keep the old cassette with the old serializer make sure to clean up the old file (#188)
* now using GitHub Actions - remove Travis-CI and Appveyor (#175)
* fixes for tests not being idempotent (#174) thanks @alex-gable
* clean up `UnhandledHTTPRequestError` - remove unused variable `cassette` in `$initialize()` method (always use `current_cassette()` to get the cassette being used) (#163) tip from @aaronwolen
* A change in latest webmockr release (`v0.7.4`) allowed for changes here to return an httr `response` object that more closely matches what httr returns in a real HTTP request. Before this the major problem was, assuming `x` is a httr `response` object, `x$request` was a `RequestSignature` object (from `webmockr`), whereas the class in a real httr response object is `request`  (#132)
* Re-factor of `Cassette` class greatly simplifying webmockr HTTP request stubbing (#98) (#173) big thanks to @alex-gable !
* `HTTPInteractionList` improvement: in checking for request matches against those on disk we were checking all requesets in a cassette - faster to check and stop when a match found. Using new factored out function to do this checking that stops when first match found. Many more tests added to check this behavior (#69)
* base64 encoded output in cassettes when using YAML serializer are now wrapped to approximately 80 character width (triggered when `preserve_exact_body_bytes=TRUE`) - this makes cassettes longer however. Implementing this brought in use of `cpp11` (first use of C++ in vcr). This makes base64 encoded response body recording consistent with how vcr's in other programming languages do it  (#41)
* `decode_compressed_response` option removed from `Cassette` class - wasn't being used and won't be used (#30)
* add additional examples to `VcrResponse` docs showing what the `update_content_length_header()` does (#29)
* `use_vcr()` changes: 1) now creates a test helper file called `setup-pkgname.R` instead of `helper-pkgname.R`; 2) now by default sets directory for fixtures using `dir = vcr_test_path("fixtures")` instead of `dir = "../fixtures"`. See other news item about `vcr_test_path`

### DOCUMENTATION

* better description of vcr at top of README (#198)
* delete unused docs folder in repository (docs built elsewhere) (#210)
* tell users that explicitly loading vcr is required in your test setup (#185) (#186) thanks @KevCaz
* added explanation of where and how `webmockr` is integrated in `Cassette` class - see section "Points of webmockr integration" in `?Cassette` (#176) (see also #173)
* improved getting started and protecting secrets sections in the introduction vignette (#170) (#172) thanks @DaveParr
* add to introduction vignette a section titled "how to ensure tests work in the absence of a real API key" (#137) (#194)



vcr 0.5.4
=========

### NEW FEATURES

* Error messages when tests using vcr fail are now simpler, primarily to reduce the space error messages take up. The user can toggle whether they get the new simplified error messages or the older format more verbose messages using the `verbose_errors` setting in the `vcr_configure()` function. In addition, `vcr_last_error()` gives the last full error, but that doesn't help in non-interactive mode; if in non-interactive mode, which most users will be in when running the entire test suite for a package, you can set an environment variable (`VCR_VERBOSE_ERRORS`) to toggle this setting (e.g.,
`Sys.setenv(VCR_VERBOSE_ERRORS=TRUE); devtools::test()`) (#121) (#154)

### MINOR IMPROVEMENTS

* changed `write_disk_path` handling internally to not run it through `normalizePath` before recording it to the cassette; passing the path through `normalizePath` was leading to the full path recorded in the cassette, which means in a package testing context that a test that uses a file on disk will (likely) only work on the machine the cassette was first created on. with relative paths in a package context, a test that has a file written on disk should now work in different testing contexts (locally, and various continuous integration platforms) (#135) (#166)
* added a bit of documentation about large files created when using vcr, and how to ignore them if needed within `.Rinstignore` and/or `.Rbuildignore` (#164)

vcr 0.5.0
=========

### NEW FEATURES

* new function `check_cassette_names` to use in your `helper-pkgname.R` file in your test suite; it checks for duplicated cassette names only. Any use of `insert_cassette()` (thereby, any use of `use_cassette()`) uses a revamped version of an internal fxn that checks for an improved list of potential problems in cassette names (#116) (#159)
* `use_vcr()` adds gitignore cassette diffs via the addition of a `gitattributes` file (#109)
* `vcr_configure()` overhaul: function no longer has each setting as a parameter; rather, it has an ellipsis (`...`), and internally we check parameters passed in. The documentation (`?vcr_configure`) lists the details for each available parameter. Importantly, each call to `vcr_configure()` now only changes the vcr settings for parameters passed in to the function; to reset all vcr settings, run `vcr_configure_reset()`  (#136) (#141)
* `insert_cassette()` and `use_cassette()` now inherit any vcr settings set by `vcr_configure()`; this wasn't happening consistently before. Most default parameter values in `insert_cassette/use_cassette` set to `NULL`, in which case they inherit from whatever values are set by `vcr_configure()`, but can be overriden (#151) (#153)

### MINOR IMPROVEMENTS

* define _serialize_, _cassette_, and _fixture_ in the README (#138) (#139)
* fix `filter_sensitive_data` parameter description in `vcr_configure` docs  (#129)
* move higher up in README a brief description of what this package does (#140)
* import `utils::getParseData` so its in namespace (#142)
* better cleanup of some stray test files left on disk (#148)
* `use_vcr()` no longer uses `context()` in example test file (#144)
* improved documentation of functions and environment variables for turning vcr on and off and when to use each of them - documentation mostly in the HTTP Testing book at https://books.ropensci.org/http-testing/lightswitch.html (#131)
* fix a `use_cassette` test (#133)
* Add assertions to `vcr_configure()` when parameters are set by the user to fail early (#156)

### BUG FIXES

* fix for handling of http requests that request image data AND do not write that data to disk; in addition, fix usage of `preserve_exact_body_bytes` when image data is in the response body (#128) thanks @Rekyt
* vcr now should handle request bodies correctly on POST requests (#143)
* Request matching was failing for empty bodies when "body" was one of the matchers (#157) (#161)
* fix to `sensitive_remove()` internal function used when the user sets `filter_sensitive_data` in `vcr_configure()`; when an env var is missing in the `filter_sensitive_data` list, `sensitive_remove()` was causing C stack errors in some cases (#160) thanks @zachary-foster
* fix for recording JSON-encoded bodies; vcr wasn't handling HTTP requests when the user set the body to be encoded as JSON (e.g., `encode="json"` with crul or httr) (#130)


vcr 0.4.0
=========

### NEW FEATURES

* vcr now can handle requests from both `crul` and `httr` that write to disk; `crul` supports this with the `disk` parameter and `httr` through the `write_disk()` function; see the section on mocking writing to disk in the http testing book https://books.ropensci.org/http-testing/vcr-usage.html#vcr-disk; also see `?mocking-disk-writing` within `webmockr` for mocking writing to disk without using `vcr`, and the section in the http testing book https://books.ropensci.org/http-testing/webmockr-stubs.html#webmockr-disk  (#81) (#125)
* vcr gains ability to completely turn off vcr for your test suite even if you're using `vcr::use_cassette`/`vcr::insert_cassette`; this is helpful if you want to run tests both with and without vcr; workflows are supported both for setting env vars on the command line as well as working interactively within R; see `?lightswitch` for details (#37)
* ignoring requests now works, with some caveats: it only works for now with `crul` (not `httr`), and works for ignoring specifc hosts, and localhosts, but not for custom callbacks. See the vcr configuration vignette https://docs.ropensci.org/vcr/articles/configuration.html#ignoring-some-requests for discussion and examples (#127)

### MINOR IMPROVEMENTS

* documentation for R6 classes should be much better now; roxygen2 now officially supports R6 classes (#123)
* added minimal cassette name checking; no spaces allowed and no file extensions allowed; more checks may be added later (#106)

### BUG FIXES

* fix handling of http response bodies that are images; we were converting raw class bodies into character, which was causing images to error, which can't be converted to character; we now check if a body can be converted to character or not and if not, leave it as is (#112) (#119) thanks @Rekyt for the report
* simple auth with package `httr` wasn't working (`htrr::authenticate()`); we were not capturing use of `authenticate`; it's been solved now (#113)
* we were not properly capturing request bodies with package `httr` requests; that's been fixed (#122)
* httr adapter was failing on second run, reading a cached response. fixed now (#124)
* `response_summary()` fixed; this function prints a summary of the http response body; sometimes this function would fail with multibyte string error because the `gsub` call would change the encoding, then would fail on the `substring` call; we now set `useBytes = TRUE` in the `gsub` call to avoid this problem (#126)


vcr 0.3.0
=========

### NEW FEATURES

* new internal method `up_to_date_interactions` in `cassette_class` now allows filtering cassettes by user specified date (#96) (#104)
* re-recording now works - see new `use_casette()` parameters `re_record_interval` and `clean_outdated_http_interactions`; you can now set a re-record interval (in seconds) so that you can for example always re-record cassettes if you don't want cassettes to be more than X days old; depends on new internal method `up_to_date_interactions` (#104) (#105)

### MINOR IMPROVEMENTS

* fix link to HTTP Testing Book: ropensci -> ropenscilabs (#100)
* add new section to HTTP Testing Book on "vcr enabled testing" with sub-sections on check vs. test, your package on CRAN, and your package on continuous integration sites (#102)

### BUG FIXES

* fix request body matching - partly through fixes to `webmockr` package (requires v0.4 or greater); more generally, makes single type request matching (e.g., just HTTP method, or just URL) possible, it was not working before, but is now working; added examples of doing single type matching (#70) (#76) (#108)
* fixed type in `cassette_class` where typo lead to not setting headers correctly in the `webmockr::wi_th()` call (#107)


vcr 0.2.6
=========

## NEW FEATURES

* gains function `use_vcr()` to setup `vcr` for your package. This requires 3 pkgs all in Suggests; so are not required if you don't need to use `use_vcr()` (#52) (#95) thanks @maelle for the feedback!
* `vcr` actually supports all four recording modes: `none`, `once`, `new_episodes`, and `all`. `once` is what's used by default. See `?recording` for description of the recording modes. For now [the test file test-ause_cassette_record_modes.R](https://github.com/ropensci/vcr/blob/master/tests/testthat/test-ause_cassette_record_modes.R) gives some examples and what to expect for each record mode; in the future the http testing book will have much more information in the _Record modes_ chapter <https://books.ropensci.org/http-testing/record-modes.html> ([commit](https://github.com/ropensci/vcr/commit/04aa5f784b18308d8f62d1b6b0be2f3e140f2a5a))

### MINOR IMPROVEMENTS

* lots of tidying for better/consistent style
* fix for a partial argument call in `as.list()`: `all` to `all.names` ([commit](https://github.com/ropensci/vcr/commit/b20a2d5ffd0f65175dee4d84aa9573f3652df1d2))

### BUG FIXES

* error thrown with `httr` due to wrong date format. the problem was in the `webmockr` package. see [ropensci/webmockr#58](https://github.com/ropensci/webmockr/issues/58) (#91) thanks @Bisaloo
* fix for `use_cassette()` when using `httr`: we weren't collecting `status_code` and storing it with the cassette (#92) thanks @Bisaloo
* fixes for `use_cassette()` for `httr`: was working fine with a single httr request, but not with 2 or more (#93) (#94) thanks @Rekyt
* in error blocks with `use_cassette()` the URL is presented from the request, and if there's a secret (API key) in the URL as a query parameter (or in any other place in the URL) then that secret is shown to the world (including if the error block happens on CI on the public web). This is fixed now; we use directives from your `filter_sensitive_data` call in `vcr_configure()` to mask secrets in error messages (#89) (#90)


vcr 0.2.2
=========

### MINOR IMPROVEMENTS

* typo fixes (#85) thanks @Rekyt
* added to docs: at least one person has reported different results using `vcr` with `devtools::check` vs. `devtools::test` (#83)
* changed suggested usage of `vcr` in test suites from `use_cassette` block wrapped in `test_that` to the other way around; leads to `testthat` pointing to the actual test line that failed rather than pointing to the start of the `use_cassette` block (#86)

### BUG FIXES

* Fix for `%||%` internal function. Was incorrectly doing logical comparison; when headers list was passed one or more of the tests in the if statement had length > 1. Dev R is testing for this (#87)


vcr 0.2.0
=========

### NEW FEATURES

* gains support for the `httr` package. `vcr` now supports `crul` and `httr`. Some of the integration for `httr` is via `webmockr`, while some of the tooling resides here in `vcr`  (#73) (#79)

### BUG FIXES

* fix handling of response bodies when not raw type (#77) (#78)


vcr 0.1.0
=========

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS, R 4.1.0
* ubuntu 16.04 (on GitHub Actions), R 4.1.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 47 reverse dependencies. Summary at https://github.com/ropensci/vcr/actions?query=workflow%3Arevdep Problems were found in checks, but were not related to vcr.

--------

This version fixes an issue introduced in v1 released a few weeks ago on CRAN that caused problems in 6 packages that depend on vcr.

Thanks very much,
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

* Submit an issue on the [Issues page](https://github.com/ropensci/vcr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g., `git clone https://github.com/<yourgithubusername>/vcr.git`
* Make sure to track progress upstream (i.e., on our version of `vcr` at `ropensci/vcr`) by doing `git remote add upstream https://github.com/ropensci/vcr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/vcr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
## revdepcheck results

We checked 52 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.1.0 (2021-05-18) |
|os       |macOS Big Sur 10.16          |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |US/Pacific                   |
|date     |2021-05-31                   |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|vcr     |1.0.0 |1.0.2 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*