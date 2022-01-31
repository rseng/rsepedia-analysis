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

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output: github_document
---

vcr
===

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![cran checks](https://cranchecks.info/badges/worst/vcr)](https://cranchecks.info/pkgs/vcr)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/vcr/workflows/R-check/badge.svg)](https://github.com/ropensci/vcr/actions/)
[![codecov](https://codecov.io/gh/ropensci/vcr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/vcr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/vcr)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/vcr)](https://cran.r-project.org/package=vcr)


Easier HTTP testing! Record HTTP requests and responses on disk and replay them for the unit tests of your R package, to make them independent from any connection, faster, and more complete. An R port of the Ruby gem [vcr](https://github.com/vcr/vcr)

## Elevator pitch

```{r child='man/rmdhunks/elevator-pitch.Rmd'} 
```

## Installation

CRAN version:

```{r eval=FALSE}
install.packages("vcr")
```

Development version:

```{r eval=FALSE}
remotes::install_github("ropensci/vcr")
```

```{r}
library("vcr")
library("crul")
```


## Docs

Check out the [HTTP testing book](https://books.ropensci.org/http-testing) and the [vcr vignettes](https://docs.ropensci.org/vcr/articles/).

## Supported HTTP libraries

* [crul](https://docs.ropensci.org/crul/)
* [httr](https://httr.r-lib.org/)

## Getting Started

```{r child='man/rmdhunks/setup.Rmd'} 
```

## Usage

```{r child='man/rmdhunks/basic-usage.Rmd'} 
```


## Terminology

```{r child='man/rmdhunks/glossary.Rmd'} 
```

## Workflows

```{r child='man/rmdhunks/workflows.Rmd'} 
```

### How it works in lots of detail

See the [vignette about internals](https://docs.ropensci.org/vcr/articles/internals.html)

### Just want to mock and not store on disk?

You're looking for [webmockr][], that vcr itself uses. 
`webmockr` only matches requests based on criteria you choose, but does not cache HTTP interactions to disk as `vcr` does.

## Configuration

```{r child='man/rmdhunks/configuration.Rmd'} 
```

## Matching/Matchers

```{r child='man/rmdhunks/matching.Rmd'} 
```

## vcr in other languages

The canonical `vcr` (in Ruby) lists ports in other languages at <https://github.com/vcr/vcr>

## Note about missing features

```{r child='man/rmdhunks/missing-features.Rmd'} 
```

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

[webmockr]: https://docs.ropensci.org/webmockr/
[crul]: https://docs.ropensci.org/crul/
[rgbif]: https://github.com/ropensci/rgbif
[rredlist]: https://github.com/ropensci/rredlist
[bold]: https://github.com/ropensci/bold
[wikitaxa]: https://github.com/ropensci/wikitaxa
[worrms]: https://github.com/ropensci/worrms
[microdemic]: https://github.com/ropensci/microdemic
[zbank]: https://github.com/ropenscilabs/zbank
[rplos]: https://github.com/ropensci/rplos
[ritis]: https://github.com/ropensci/ritis
---
title: "Introduction to vcr"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{1. vcr introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
	comment = "#>",
	collapse = TRUE,
	warning = FALSE,
	message = FALSE,
  eval = FALSE
)
```

vcr introduction
================

`vcr` is an R port of the Ruby gem [VCR](https://github.com/vcr/vcr) (i.e., a translation, there's no Ruby here :))

`vcr` helps you stub and record HTTP requests so you don't have to repeat HTTP requests.

The main use case is for unit tests, but you can use it outside of the unit test use case.

`vcr` works with the `crul` and `httr` HTTP request packages.

Check out the [HTTP testing book](https://books.ropensci.org/http-testing/) for a lot more documentation on `vcr`, `webmockr`, and `crul`, and other packages.

## Elevator pitch

```{r child='../man/rmdhunks/elevator-pitch.Rmd', eval=TRUE} 
```

## Installation

CRAN

```{r eval=FALSE}
install.packages("vcr")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/vcr")
```

```{r}
library("vcr")
```

## Getting Started

```{r child='../man/rmdhunks/setup.Rmd', eval=TRUE} 
```

## Basic usage

```{r child='../man/rmdhunks/basic-usage.Rmd', eval=TRUE} 
```


## Terminology

```{r child='../man/rmdhunks/glossary.Rmd', eval=TRUE} 
```

## Workflows

```{r child='../man/rmdhunks/workflows.Rmd', eval=TRUE} 
```

## Configuration

```{r child='../man/rmdhunks/configuration.Rmd', eval=TRUE} 
```

## Matching/Matchers

```{r child='../man/rmdhunks/matching.Rmd', eval=TRUE} 
```

## Note about missing features

```{r child='../man/rmdhunks/missing-features.Rmd', eval=TRUE} 
```


---
title: "Configure vcr request matching"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{3. request matching}
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

Request matching
================

```{r child='../man/rmdhunks/request-matching-vignette.Rmd', eval=TRUE} 
```

## More documentation

Check out the [http testing book](https://books.ropensci.org/http-testing/) for a lot more documentation on `vcr`, `webmockr`, and `crul`

---
title: "Debugging your tests that use vcr"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{4. vcr tests debugging}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(vcr)
```

```{r child='../man/rmdhunks/debugging-vignette.Rmd', eval=TRUE} 
```
---
title: "Design of vcr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Design of vcr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r child='../man/rmdhunks/vcr-design.Rmd', eval=TRUE} 
```
---
title: "How vcr works, in a lot of details"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{internals}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(vcr)
```

```{r child='../man/rmdhunks/internals-vignette.Rmd', eval=TRUE} 
```



---
title: "Configure vcr"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{2. vcr configuration}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
	comment = "#>",
	collapse = TRUE,
	warning = FALSE,
	message = FALSE,
	eval = FALSE
)
```

vcr configuration
=================

`vcr` configuration

```{r}
library("vcr")
```

```{r child='../man/rmdhunks/configuration-vignette.Rmd', eval=TRUE} 
```

## More documentation

Check out the [http testing book](https://books.ropensci.org/http-testing/) for a lot more documentation on `vcr`, `webmockr`, and `crul`
---
title: "Record modes"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{vcr record modes}
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

Request matching
================

```{r child='../man/rmdhunks/record-modes.Rmd', eval=TRUE} 
```

## More documentation

Check out the [http testing book](https://books.ropensci.org/http-testing/) for a lot more documentation on `vcr`, `webmockr`, and `crul`

---
title: "Mocking writing to disk"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Writing to disk}
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

Request matching
================

```{r setup}
library("vcr")
```

```{r child='../man/rmdhunks/write-to-disk.Rmd', eval=TRUE} 
```

## More documentation

Check out the [http testing book](https://books.ropensci.org/http-testing/) for a lot more documentation on `vcr`, `webmockr`, and `crul`

---
title: "Why and how to edit your vcr cassettes?"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{5. cassette manual editing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(vcr)
```


```{r child='../man/rmdhunks/cassette-editing-vignette.Rmd', eval=TRUE} 
```

## Conclusion

In this vignette we saw why and how to edit your vcr cassettes.
We also presented approaches that use webmockr instead of vcr for mocking API responses.
We mentioned editing cassettes by hand but you could also write a script using the `yaml` 
or `jsonlite` package to edit your YAML/JSON cassettes programmatically.
---
title: "Turning vcr on and off"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Turning vcr on and off}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(vcr)
```

```{r child='../man/rmdhunks/lightswitch.Rmd', eval=TRUE} 
```

```{r, echo=FALSE}

vcr::turn_on()
```
## Why edit cassettes?

By design vcr is very good at recording HTTP interactions that actually took place.
Now sometimes when testing/demo-ing your package you will want to use _fake_ HTTP interactions.
For instance:

* What happens if the web API returns a 503 code? Is there an informative error?
* What happens if it returns a 503 and then a 200 code? Does the retry work?
* What if the API returns too much data for even simple queries and you want to make your cassettes smaller?

In all these cases, you can edit your cassettes as long as you are aware of the risks!

## Risks related to cassette editing

* If you use a vcr cassette where you replace a 200 code with a 503 code, and vcr is turned off, 
the test will fail because the API will probably not return an error. Use `vcr::skip_if_vcr_off()`.
* If you edit cassettes by hand you can't re-record them easily, you'd need to re-record them then re-apply your edits.

Therefore you'll need to develop a good workflow.

## Example 1: test using an edited cassette with a 503

First, write your test e.g.

```r
vcr::use_cassette("api-error", {
  testhat("Errors are handled well", {
    vcr::skip_if_vcr_off()
    expect_error(call_my_api()), "error message")
  })
})

```

Then run your tests a first time.

1. It will fail
2. It will have created a cassette under `tests/fixtures/api-error.yml` that looks 
something like

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

You can edit to (new status code)

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
      status_code: '503'
```

And run your test again, it should pass!
Note the use of `vcr::skip_if_vcr_off()`: if vcr is turned off, there is a real 
API request and most probably this request won't get a 503 as a status code.

### The same thing with webmockr

The advantage of the approach involving editing cassettes is that you only learn
 one thing, which is vcr.
Now, by using the webmockr directly in your tests, you can also test for the 
behavior of your package in case of errors.
Below we assume `api_url()` returns the URL `call_my_api()` calls.

```r
testhat("Errors are handled well", {
  webmockr::enable()
  stub <- webmockr::stub_request("get", api_url())
  webmockr::to_return(stub, status = 503)
  expect_error(call_my_api()), "error message")
  webmockr::disable()

})

```

A big pro of this approach is that it works even when vcr is turned off.
A con is that it's quite different from the vcr syntax.

## Example 2: test using an edited cassette with a 503 then a 200

Here we assume your package contains some sort of [retry](https://blog.r-hub.io/2020/04/07/retry-wheel/).

First, write your test e.g.

```r
vcr::use_cassette("api-error", {
  testhat("Errors are handled well", {
    vcr::skip_if_vcr_off()
    expect_message(thing <- call_my_api()), "retry message")
    expect_s4_class(thing, "data.frame")
  })
})

```

Then run your tests a first time.

1. It will fail
2. It will have created a cassette under `tests/fixtures/api-error.yml` that looks 
something like

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

You can duplicate the HTTP interaction, and make the first one return a 503 status code.
vcr will first use the first interaction, then the second one, when making the same request.

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
      status_code: '503'
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

And run your test again, it should pass!
Note the use of `vcr::skip_if_vcr_off()`: if vcr is turned off, there is a real 
API request and most probably this request won't get a 503 as a status code.

### The same thing with webmockr

The advantage of the approach involving editing cassettes is that you only learn
 one thing, which is vcr.
Now, by using the webmockr directly in your tests, you can also test for the 
behavior of your package in case of errors.
Below we assume `api_url()` returns the URL `call_my_api()` calls.

```r
testhat("Errors are handled well", {
  webmockr::enable()
  stub <- webmockr::stub_request("get", api_url())
  stub %>%
  to_return(status = 503)  %>%
  to_return(status = 200, body = "{\n  \"args\": {}, \n  \"headers\": {\n    \"Accept\": \"application/json,
        text/xml, application/xml, */*\", \n    \"Accept-Encoding\": \"gzip, deflate\",
        \n    \"Connection\": \"close\", \n    \"Host\": \"httpbin.org\", \n    \"User-Agent\":
        \"libcurl/7.54.0 r-curl/3.2 crul/0.5.2\"\n  }, \n  \"origin\": \"111.222.333.444\",
        \n  \"url\": \"https://eu.httpbin.org/get\"\n}\n", headers = list(b = 6))
  expect_message(thing <- call_my_api()), "retry message")
    expect_s4_class(thing, "data.frame")
  webmockr::disable()

})

```

The pro of this approach is the elegance of the stubbing, with the two different responses.
Each webmockr function like `to_return()` even has an argument `times` indicating the 
number of times the given response should be returned.

The con is that on top of being different from vcr, in this case where we also needed 
a good response in the end (the one with a 200 code, and an actual body), writing the
mock is much more cumbersome than just recording a vcr cassette.
To match previously recorded requests, `vcr` has to try to match new 
HTTP requests to a previously recorded one. By default, we match on 
HTTP method (e.g., `GET`) and URI (e.g., `http://foo.com`), following 
Ruby's VCR gem.

You can customize how we match requests with one or more of the 
following options, some of which are on by default, some of
which can be used together, and some alone.

* `method`: Use the **method** request matcher to match requests on the HTTP method
(i.e. GET, POST, PUT, DELETE, etc). You will generally want to use
this matcher. The **method** matcher is used (along with the **uri** matcher)
by default if you do not specify how requests should match.
* `uri`: Use the **uri** request matcher to match requests on the request URI. The
**uri** matcher is used (along with the **method** matcher) by default
if you do not specify how requests should match.
* `host`: Use the **host** request matcher to match requests on the request host.
You can use this (alone, or in combination with **path**) as an
alternative to **uri** so that non-deterministic portions of the URI
are not considered as part of the request matching.
* `path`: Use the **path** request matcher to match requests on the path portion
of the request URI. You can use this (alone, or in combination with **host**)
as an alternative to **uri** so that non-deterministic portions of the URI
* `query`: Use the **query** request matcher to match requests on the query string
portion of the request URI. You can use this (alone, or in combination with
others) as an alternative to **uri** so that non-deterministic portions of
the URI are not considered as part of the request matching.
* `body`: Use the **body** request matcher to match requests on the request body.
* `headers`: Use the **headers** request matcher to match requests on the request headers.

You can set your own options by tweaking the `match_requests_on` parameter in 
`use_cassette()`:

```{r}
library(vcr)
```


```{r echo=FALSE, eval=FALSE}
unlink(file.path(cassette_path(), "foo_bar.yml"))
```

```{r eval = FALSE}
use_cassette(name = "foo_bar", {
    cli$post("post", body = list(a = 5))
  }, 
  match_requests_on = c('method', 'headers', 'body')
)
```

```{r echo=FALSE, eval=FALSE}
unlink(file.path(cassette_path(), "foo_bar.yml"))
```

## Matching

### headers

```{r echo=FALSE, eval=FALSE}
unlink(file.path(cassette_path(), "nothing_new.yml"))
```

```{r eval = FALSE}
library(crul)
library(vcr)
cli <- crul::HttpClient$new("https://httpbin.org/get", 
  headers = list(foo = "bar"))
use_cassette(name = "nothing_new", {
    one <- cli$get()
  }, 
  match_requests_on = 'headers'
)
cli$headers$foo <- "stuff"
use_cassette(name = "nothing_new", {
    two <- cli$get()
  }, 
  match_requests_on = 'headers'
)
one$request_headers
two$request_headers
```

```{r echo=FALSE, eval=FALSE}
unlink(file.path(cassette_path(), "nothing_new.yml"))
```

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

```{r echo=FALSE}
unlink(file.path(cassette_path(), "one.yml"))
```

```{r eval = FALSE}
use_cassette(name = "one", {
    cli$post("post", body = list(a = 5))
  },
  match_requests_on = c('method', 'headers', 'body')
)
```

For more details refer to the [request matching vignette](https://docs.ropensci.org/vcr/articles/request_matching.html).

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
* _vcr_: the name comes from the idea that we want to record something and play it back later, like a vcr
* _cassette_: A _thing_ to record HTTP interactions to. Right now the only option is the file system (writing to files), but in the future could be other things, e.g. a key-value store like Redis
* _fixture_: A fixture is something used to consistently test a piece of software. In this case, a cassette (just defined above) is a fixture - used in unit tests. If you use our setup function `vcr_setup()` the default directory created to hold cassettes is called `fixtures/` as a signal as to what the folder contains.
* Persisters: how to save requests - currently only option is the file system
* _serialize_: translating data into a format that can be stored; here, translate HTTP request and response data into a representation on disk to read back later
* Serializers: how to serialize the HTTP response - currently only option is YAML; other options in the future could include e.g. JSON
* _insert cassette_: create a cassette (all HTTP interactions will be recorded to this cassette)
* _eject cassette_: eject the cassette (no longer recording to that cassette)
* _replay_: refers to using a cached result of an http request that was recorded earlier

```{r, echo=FALSE, message=FALSE}
library("vcr")
```

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

```{r echo=FALSE, results='hide', eval=identical(Sys.getenv("IN_PKGDOWN"), "true")}
suppressPackageStartupMessages(require(vcr, quietly = TRUE))
unlink(file.path(cassette_path(), "helloworld.yml"))
vcr_configure(dir = tempdir())
```

```{r, eval=identical(Sys.getenv("IN_PKGDOWN"), "true")}
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

```{r, eval=identical(Sys.getenv("IN_PKGDOWN"), "true")}
system.time(
  use_cassette(name = "helloworld", {
    cli$get("get")
  })
)
```

```{r echo=FALSE, eval=identical(Sys.getenv("IN_PKGDOWN"), "true")}
unlink(file.path(cassette_path(), "helloworld.yml"))
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
See also the [configuration vignette](https://docs.ropensci.org/vcr/articles/configuration.html).

We set the following defaults:

```{r, echo = FALSE, results='asis', collapse=TRUE}
defaults <- rev(vcr_config_defaults())
defaults[unlist(lapply(defaults, is.character))] <- paste0('"', defaults[unlist(lapply(defaults, is.character))], '"')
cat(sprintf("* %s = `%s`\n", names(defaults), defaults))
```


You can get the defaults programmatically with

```r
vcr_config_defaults()
```

You can change all the above defaults with `vcr_configure()`:

```r
vcr_configure()
```

Calling `vcr_configuration()` gives you some of the more important configuration parameters in a nice tidy print out

```{r}
vcr_configuration()
```


<!-- `use_cassette()` is an easier approach. An alternative is to use
`insert_cassett()` + `eject_cassette()`.

`use_cassette()` does both insert and eject operations for you, but
you can instead do them manually by using the above functions. You do have
to eject the cassette after using insert. -->

For more details refer to the [configuration vignette](https://docs.ropensci.org/vcr/articles/configuration.html)
**The steps**

1. Use either `vcr::use_cassette()` or `vcr::insert_cassette()`
  a. If you use `vcr::insert_cassette()`, make sure to run `vcr::eject_cassette()` when you're done to stop recording
2. When you first run a request with `vcr` there's no cached data to use, so we allow HTTP requests until your request is done.
3. Before we run the real HTTP request, we "stub" the request with `webmockr` so that future requests will match the stub.
This stub is an R6 class with details of the interaction (request + response), but is not on disk.
4. After the stub is made, we run the real HTTP request.
5. We then disallow HTTP requests so that if the request is done again we use the cached response
6. The last thing we do is write the HTTP interaction to disk in a mostly human readable form.

When you run that request again using `vcr::use_cassette()` or `vcr::insert_cassette()`:

* We use `webmockr` to match the request to cached requests, and since we stubbed the request the first time we used the cached response.

Of course if you do a different request, even slightly (but depending on which matching format you decided to use), then
the request will have no matching stub and no cached response, and then a real HTTP request is done - we then cache it, then subsequent requests will pull from that cached response.

`webmockr` has adapters for each R client (crul and httr) - so that we actually intercept HTTP requests when `webmockr` is loaded and the user turns it on. So, `webmockr` doesn't actually require an internet or localhost connection at all, but can do its thing just fine by matching on whatever the user requests to match on. In fact, `webmockr` doesn't allow real HTTP requests by default, but can be toggled off of course.

The main use case we are going for in `vcr` is to deal with real HTTP requests and responses, so we allow real HTTP requests when we need to, and turn it off when we don't.

This gives us a very flexible and powerful framework where we can support `webmockr` and `vcr` integration for any number of R clients for HTTP requests and support many different formats serialized to disk.
There's a number of features in this package that are not yet supported, but for which their parameters are found in the package. 

We've tried to make sure the parameters that are ignored are marked as such. Keep an eye out for package updates for changes in these parameters, and/or let us know you want it and we can move it up in the priority list.
* **Setup vcr for your package with `vcr::use_vcr()`**
* Tweak the configuration to protect your secrets
* **Sprinkle your tests with `vcr::use_cassette()` to save HTTP interactions to disk in "cassettes" files**
* If you want to test for package behavior when the API returns e.g. a 404 or 503 code, edit the cassettes, or use [webmockr](https://docs.ropensci.org/webmockr/)

Now your tests can work without any internet connection!

[Demo of adding vcr testing to an R package](https://github.com/maelle/exemplighratia/pull/2/files), [corresponding narrative](https://books.ropensci.org/http-testing/vcr.html).
This section explains `vcr`'s internal design and architecture.

### Where vcr comes from and why R6

`vcr` was "ported" from the Ruby gem (aka, library) of the same name[^1]. 
Because it was ported from Ruby, an object-oriented programming language
I thought it would be easier to use an object system in R that most
closely resemble that used in Ruby (at least in my opinion). This 
thinking lead to choosing [R6][]. The exported functions users interact
with are not R6 classes, but are rather normal R functions. However,
most of the internal code in the package uses R6. Thus, familiarity
with R6 is important for people that may want to contribute to `vcr`,
but not required at all for `vcr` users.

### Principles

#### An easy to use interface hides complexity

As described above, `vcr` uses R6 internally, but users interact with
normal R functions. Internal functions that are quite complicated are
largely R6 while exported, simpler functions users interact with are
normal R functions.

#### Class/function names are inherited from Ruby vcr

Since R `vcr` was ported from Ruby, we kept most of the names of 
functions/classes and variables. So if you're wondering about why 
a function, class, or variable has a particular name, its derivation 
can not be found out in this package, for the most part that is.

#### Hooks into HTTP clients

Perhaps the most fundamental thing about that this package work is
how it knows what HTTP requests are being made. This stumped me for 
quite a long time. When looking at Ruby vcr, at first I thought it
must be "listening" for HTTP requests somehow. Then I found out about
[monkey patching][mp] (see also [this blog post as a good 
summary][mpblog]); that's how it's achieved in Ruby. That is, the Ruby
vcr package literally overrides certain methods in Ruby HTTP clients, 
hijacking internals of the HTTP clients. 

However, monkey patching is not allowed in R. Thus, in R we have to
somehow have "hooks" into HTTP clients in R. Fortunately, Scott is the
maintainer of one of the HTTP clients, `crul`, so was able to quickly
create a hook. Very fortunately, there was already a hook mechanism
in the `httr` package. 

The actual hooks are not in `vcr`, but in `webmockr`. `vcr` depends on
`webmockr` for hooking into HTTP clients `httr` and `crul`. 

### Internal classes

An overview of some of the more important aspects of vcr.

#### Configuration

An internal object (`vcr_c`) is created when `vcr` is loaded with
the default vcr configuration options inside of an R6 class `VCRConfig` -
see <https://github.com/ropensci/vcr/blob/master/R/onLoad.R>. This
class is keeps track of default and user specified configuration options.
You can access `vcr_c` using triple namespace `:::`, though it is not
intended for general use. Whenever you make calls to `vcr_configure()`
or other configuration functions, `vcr_c` is affected.

#### Cassette class

`Cassette` is an R6 class that handles internals/state for each cassette.
Each time you run `use_cassette()` this class is used. The class has quite
a few methods in it, so there's a lot going on in the class. Ideally
the class would be separated into subclasses to handle similar sets
of logic, but there's not an easy way to do that with R6. 

Of note in `Cassette` is that when called, within the `initialize()` 
call `webmockr` is used to create webmockr stubs. 

#### How HTTP requests are handled

Within `webmockr`, there are calls to the vcr class `RequestHandler`, which
has child classes `RequestHandlerCrul` and `RequestHandlerHttr` for 
`crul` and `httr`, respectively. These classes determine what to do with
each HTTP request. The options for each HTTP request include:

- **Ignored** You can ignore HTTP requests under certain rules using the 
configuration options `ignore_hosts` and `ignore_localhost`
- **Stubbed by vcr** This is an HTTP request for which a match is found
in the cassette defined in the `use_cassette()`/`insert_cassette()` call.
In this case the matching request/response from the cassette is returned
with no real HTTP request allowed.
- **Recordable** This is an HTTP request for which no match is found
in the cassette defined in the `use_cassette()`/`insert_cassette()` call.
In this case a real HTTP request is allowed, and the request/response is
recorded to the cassette.
- **Unhandled** This is a group of cases, all of which cause an error
to be thrown with a message trying to help the user figure out how to
fix the problem.

If you use [vcr logging][logging] you'll see these categories in your logs.

#### Serializers

Serializers handle in what format cassettes are written to files on disk.
The current options are YAML (default) and JSON. YAML was implemented first
in `vcr` because that's the default option in Ruby vcr.

An R6 class `Serializer` is the parent class for all serializer types;
`YAML` and `JSON` are both R6 classes that inherit from `Serializer`. Both
`YAML` and `JSON` define just two methods: `serialize()` and `deserialize()`
for converting R structures to yaml or json, and converting yaml or json back
to R structures, respectively.


### Environments

#### Logging

An internal environment (`vcr_log_env`) is used when you use logging.
At this point it only keeps track of one variable - `file` - to be able
to refer to what file is used for logging across many classes/functions
that need to write to the log file.

#### A bit of housekeeping

Another internal environment (`vcr__env`) is used to keep track of a
few items, including the current cassette in use, and the last vcr error.

#### Lightswitch

Another internal environment (`light_switch`) is used to keep track of users
turning on and off `vcr`. See `?lightswitch`.



[^1]: The first version of Ruby's vcr was released in February 2010
https://rubygems.org/gems/vcr/versions/0.1.0. Ruby vcr source code:
https://github.com/vcr/vcr/




[R6]: https://adv-r.hadley.nz/r6.html
[mp]: https://en.wikipedia.org/wiki/Monkey_patch
[mpblog]: https://culttt.com/2015/06/17/what-is-monkey-patching-in-ruby/
[logging]: https://docs.ropensci.org/vcr/articles/debugging.html?q=logging#logging-1
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
Record modes dictate under what circumstances http requests/responses are
recorded to cassettes (disk). Set the recording mode with the parameter 
`record` in the `use_cassette()` and `insert_cassette()` functions.

## once

The `once` record mode will:

- Replay previously recorded interactions.
- Record new interactions if there is no cassette file.
- Cause an error to be raised for new requests if there is a cassette file.


It is similar to the `new_episodes` record mode, but will prevent new,
unexpected requests from being made (i.e. because the request URI changed
or whatever).

`once` is the default record mode, used when you do not set one.

## none

The `none` record mode will:

- Replay previously recorded interactions.
- Cause an error to be raised for any new requests.


This is useful when your code makes potentially dangerous
HTTP requests.  The `none` record mode guarantees that no
new HTTP requests will be made.

## new_episodes

The `new_episodes` record mode will:

- Record new interactions.
- Replay previously recorded interactions.


It is similar to the `once` record mode, but will **always** record new
interactions, even if you have an existing recorded one that is similar
(but not identical, based on the `match_request_on` option).

## all

The `all` record mode will:

- Record new interactions.
- Never replay previously recorded interactions.


This can be temporarily used to force vcr to re-record
a cassette (i.e. to ensure the responses are not out of date)
or can be used when you simply want to log all HTTP requests.

```{r}
library("vcr")
```

You can also get the default configuration variables via `vcr_config_defaults()`

```{r}
vcr_config_defaults()
```

These defaults are set when you load `vcr` - you can override any of them as described below.

## Set configuration variables

Use `vcr_configure()` to set configuration variables.

For example, set a single variable:

```{r}
vcr_configure(
  dir = "foobar/vcr_cassettes"
)
```

Or many at once:

```{r}
vcr_configure(
  dir = "foobar/vcr_cassettes",
  record = "all"
)
```

## Re-set to defaults

```{r}
vcr_configure_reset()
```

## Details on some of the config options

### dir

Directory where cassettes are stored

```{r}
vcr_configure(dir = "new/path")
```

### record

The record mode

One of: 'all', 'none', 'new_episodes', 'once'. See `?recording` for info on the options

```{r}
vcr_configure(record = "new_episodes")
```

### match_requests_on 

Customize how `vcr` matches requests

```{r}
vcr_configure(match_requests_on = c('query', 'headers'))
```

### allow_unused_http_interactions 

Allow HTTP connections when no cassette

Default is `TRUE`, and thus does not error when http interactions are unused. You 
can set to `FALSE` in which case vcr errors when a cassette is ejected and 
not all http interactions have been used.

```{r}
vcr_configure(allow_unused_http_interactions = FALSE)
```

### serialize_with

Which serializer to use. Right now only option is "yaml"

```{r}
vcr_configure(serialize_with = "yaml")
```

### persist_with

Which persister to use. Right now only option is "FileSystem"

```{r}
vcr_configure(persist_with = "FileSystem")
```

### ignoring some requests

**ignore_hosts**

Specify particular hosts to ignore. By ignore, we mean that 
real HTTP requests to the ignored host will be allowed to occur, while
all others will not.

```{r}
vcr_configure(ignore_hosts = "google.com")
```

**ignore_localhost**

Ignore all localhost requests

```{r}
vcr_configure(ignore_localhost = TRUE)
```

**ignore_request**

THIS DOESN'T WORK YET

**How to ignore requests**

For ignoring requests, you can for example, have real http requests go through (ignored by `vcr`) while other requests are handled by `vcr`. For example, let's say you want requests to `google.com` to be ignored:

```{r eval=FALSE}
vcr_configure(ignore_hosts = "google.com")
use_cassette("foo_bar", {
  crul::HttpClient$new("https://httpbin.org/get")$get()
  crul::HttpClient$new("https://google.com")$get()
})
```

The request to httpbin.org will be handled by `vcr`, a cassette created for the request/response to that url, while the google.com request will be ignored and not cached at all.

Note: ignoring requests only works for the `crul` package for now; it should work for `httr` in a later `vcr` version.

### uri_parse

Which uri parser to use

By default we use `crul::url_parse`, but you can use a different one. Remember 
to pass in the function quoted, and namespaced.

```{r}
vcr_configure(uri_parser = "urltools::url_parse")
```

### preserve_exact_body_bytes

Some HTTP servers are not well-behaved and respond with invalid data. Set 
`preserve_exact_body_bytes` to `TRUE` to base64 encode the result body in 
order to preserve the bytes exactly as-is. `vcr` does not do this by
default, since base64-encoding the string removes the human readability 
of the cassette.

```{r}
vcr_configure(preserve_exact_body_bytes = TRUE)
```

### filter_sensitive_data 

A named list of values to replace. Sometimes your package or script is 
working with sensitive tokens/keys, which you do not want to accidentally 
share with the world.

Before recording (writing to a cassette) we do the replacement and then when
reading from the cassette we do the reverse replacement to get back
to the real data.

```r
vcr_configure(
  filter_sensitive_data = list("<some_api_key>" = Sys.getenv('MY_API_KEY'))
)
```

Before recording to disk, the env var `MY_API_KEY` is retrieved from your machine, 
and we find instances of it, and replace with `<some_api_key>`. When replaying 
to create the HTTP response object we put the real value of the env var 
back in place.

To target specific request or response headers see `filter_request_headers`
and `filter_response_headers`.

### filter_request_headers

Expects a character vector or a named list. If a character vector, or any
unnamed element in a list, the request header is removed before being 
written to the cassette.

If a named list is passed, the name is the header and the value is the
value with which to replace the real value.

A request header you set to remove or replace is only removed/replaced
from the cassette, and any requests using a cassette, but will still be in
your crul or httr response objects on a real request that creates the
cassette.

Examples:

```r
vcr_configure(
  filter_request_headers = "Authorization"
)
vcr_configure(
  filter_request_headers = c("Authorization", "User-Agent")
)
vcr_configure(
  filter_request_headers = list(Authorization = "<<<not-my-bearer-token>>>")
)
```

### filter_response_headers

Expects a character vector or a named list. If a character vector, or any
unnamed element in a list, the response header is removed before being
written to the cassette.

If a named list is passed, the name is the header and the value is the
value with which to replace the real value.

A response header you set to remove or replace is only removed/replaced
from the cassette, and any requests using a cassette, but will still be in
your crul or httr response objects on a real request that creates the
cassette.

Examples:

```r
vcr_configure(
  filter_response_headers = "server"
)
vcr_configure(
  filter_response_headers = c("server", "date")
)
vcr_configure(
  filter_response_headers = list(server = "fake-server")
)
```

### filter_query_parameters

Expects a character vector or a named list. If a character vector, or any
unnamed element in a list, the query parameter is removed (both parameter
name and value) before being written to the cassette.

If a named list is passed, the name is the query parameter name and the
value is the value with which to replace the real value.

A response header you set to remove or replace is only removed/replaced
from the cassette, and any requests using a cassette, but will still be in
your crul or httr response objects on a real request that creates the
cassette.

Beware of your `match_requests_on` option when using this filter. If you
filter out a query parameter it's probably a bad idea to match on `query`
given that there is no way for vcr to restore the exact http request
from your cassette after one or more query parameters is removed or changed.
One way you could filter a query parameter and still match on query or 
at least on the complete uri is to use replacement behavior (a named list),
but instead of `list(a="b")` use two values `list(a=c("b","c"))`, where
"c" is the string to be stored in the cassette. You could of course replace
those values with values from environment variables so that you obscure 
the real values if your code is public.

Examples:

```r
# completely drop parameter "user"
vcr_configure(
  filter_query_parameters = "user"
)
# completely drop parameters "user" and "api_key"
vcr_configure(
  filter_query_parameters = c("user", "api_key")
)
# replace the value of parameter "api_key" with "fake-api-key"
# NOTE: in this case there's no way to put back any value on
# subsequent requests, so we have to match by dropping this
# parameter value before comparing URIs
vcr_configure(
  filter_query_parameters = list(api_key = "fake-api-key")
)
# replace the value found at Sys.getenv("MY_API_KEY") of parameter
# "api_key" with the value "foo". When using a cassette on subsequent
# requests, we can replace "foo" with the value at Sys.getenv("MY_API_KEY")
# before doing the URI comparison
vcr_configure(
  filter_query_parameters = list(api_key = c(Sys.getenv("MY_API_KEY"), "foo"))
)
```
Cassette names:

- Should be meaningful so that it is obvious to you what test/function
they relate to. Meaningful names are important so that you can quickly
determine to what test file or test block a cassette belongs. Note that
vcr cannot check that your cassette names are meaningful.
- Should not be duplicated. Duplicated cassette names would lead to
a test using the wrong cassette.
- Should not have spaces. Spaces can lead to problems in using file paths.
- Should not include a file extension. vcr handles file extensions for
the user.
- Should not have illegal characters that can lead to problems in using
file paths: `/`, `?`, `<`, `>`, `\\`, `:`, `*`, `|`, and `\"`
- Should not have control characters, e.g., `\n`
- Should not have just dots, e.g., `.` or `..`
- Should not have Windows reserved words, e.g., `com1`
- Should not have trailing dots
- Should not be longer than 255 characters

`vcr::check_cassette_names()` is meant to be run during your tests, from a
`setup-pkgname.R` file inside the `tests/testthat` directory. It only
checks that cassette names are not duplicated. Note that if you do 
need to have duplicated cassette names you can do so by using the
`allowed_duplicates` parameter in `check_cassette_names()`. A helper function
`check_cassette_names()` runs inside [insert_cassette()] that checks
that cassettes do not have: spaces, file extensions, unaccepted
characters (slashes)
If you have http requests for which you write the response to disk, then
use `vcr_configure()` to set the `write_disk_path` option. See more about 
the write_disk_path configuration option in `vignette("configuration", package = "vcr")`.

Here, we create a temporary directory, then set the fixtures

```{r}
tmpdir <- tempdir()
vcr_configure(
  dir = file.path(tmpdir, "fixtures"),
  write_disk_path = file.path(tmpdir, "files")
)
```

Then pass a file path (that doesn't exist yet) to crul's `disk` parameter.
`vcr` will take care of handling writing the response to that file in
addition to the cassette.

```{r}
library(crul)
## make a temp file
f <- tempfile(fileext = ".json")
## make a request
cas <- use_cassette("test_write_to_disk", {
  out <- HttpClient$new("https://httpbin.org/get")$get(disk = f)
})
file.exists(out$content)
out$parse()
```

This also works with `httr`. The only difference is that you write to disk
with a function `httr::write_disk(path)` rather than a parameter.

Note that when you write to disk when using `vcr`, the cassette is slightly
changed. Instead of holding the http response body itself, the cassette
has the file path with the response body.

```yaml
http_interactions:
- request:
    method: get
    uri: https://httpbin.org/get
  response:
    headers:
      status: HTTP/1.1 200 OK
      access-control-allow-credentials: 'true'
    body:
      encoding: UTF-8
      file: yes
      string: /private/var/folders/fc/n7g_vrvn0sx_st0p8lxb3ts40000gn/T/Rtmp5W4olr/files/file177e2e5d97ec.json
```

And the file has the response body that otherwise would have been in the `string`
yaml field above:

```json
{
  "args": {}, 
  "headers": {
    "Accept": "application/json, text/xml, application/xml, */*", 
    "Accept-Encoding": "gzip, deflate", 
    "Host": "httpbin.org", 
    "User-Agent": "libcurl/7.54.0 r-curl/4.3 crul/0.9.0"
  }, 
  "origin": "24.21.229.59, 24.21.229.59", 
  "url": "https://httpbin.org/get"
}
```

```{r echo=FALSE}
invisible(vcr_configure_reset())
```
Sometimes your tests using a vcr cassette will fail and you will want to debug them.

## An HTTP request has been made that vcr does not know how to handle

If you get an error starting with "An HTTP request has been made that vcr does not know how to handle:" when running your tests, 
it means that the code in your test makes an HTTP request 
for which there is no matching information in the cassette you are using. 
You might have added a request, or changed one slightly.

The easy fix is: delete the cassette and re-run the test to re-record the cassette. 
Run the test a second time to ensure all is well. 
If not, escalate to the next paragraph.

Maybe you didn't actually want to change the request you are making.
Make sure the requests do not contain something random, or something related to e.g. what time it is now, in the URI (`http://foo.com?time=13`).
To make sure things are not varying, you might want to use mocking (of e.g. a function returning the current time), setting a random seed, using [`withr`](https://withr.r-lib.org/) (for e.g. setting an option to a certain value in your test).

### Actual debugging

Ideally you will want to run the code of the tests as if it were run inside tests, 
in particular, using the same vcr cassette.

### Prepare your debugging environment

You will first need to load either the vcr helper `tests/testthat/helper-vcr.R` (e.g. via `devtools::load_all()`) 
or source the vcr setup file `tests/testthat/setup-vcr.R` i.e. the file with these lines (and maybe others)

```r
library("vcr")
invisible(vcr::vcr_configure(
  dir = vcr::vcr_test_path("fixtures"),
  filter_sensitive_data = list("<<github_api_token>>" = Sys.getenv('GITHUB_PAT'))
))
vcr::check_cassette_names()

```

If instead of `vcr::vcr_test_path("fixtures")` you see `"../fixtures"`, 
replace `"../fixtures"` with `vcr::vcr_test_path("fixtures")`,
as `vcr::vcr_test_path()` is a function that is meant to help exactly what you will want: 
have the path to `tests/fixtures/` work from tests and from the root (which is where you will be running the code to debug it).

So that is one step (loading the vcr helper or sourcing the vcr setup file), 
or maybe two (if you also had to replace `"../fixtures"` with `vcr::vcr_test_path("fixtures")`).

### Debugging itself

Now look at the test whose code you are trying to debug e.g.

```r
foo <- function() crul::ok('https://httpbin.org/get')

test_that("foo works", {
  vcr::use_cassette("testing", {
    x <- foo()
  })
  expect_true(x)
})
```

If you want to run the code as if you were in the test,

```r
foo <- function() crul::ok('https://httpbin.org/get')
vcr::insert_cassette("testing") # it will be created if needed
x <- foo()
x
# further interactive debugging and fixes
vcr::eject_cassette("testing")
```

### Logging

You can use vcr's built in logging to help in your debugging process. To configure logging,
use the `vcr_configure()` function, and set `log=TRUE` and set options for logging on the
`log_opts` parameter as a named list. See `?vcr_configure` for details.

Here, we are setting our log file to be a temporary file that will be cleaned up at the end
of the R session. Here, the file extension is `.log`, but the file extension does not matter.

```r
vcr::vcr_configure(
  dir = vcr::vcr_test_path("fixtures"),
  log = TRUE,
  log_opts = list(file = file.path(tempdir(), "vcr.log"))
)
```

With `log=TRUE` you can continue with debugging. Open the log file you set in a text
editor or other location; or examine in your shell/terminal.

As an example, after running the block above

```r
foo <- function() crul::ok('https://httpbin.org/get')

test_that("foo works", {
  vcr::use_cassette("testing", {
    x <- foo()
  })
  expect_true(x)
})
```

If we open the log file we'll see the logs for each step vcr takes in handling an HTTP request.
The logs have information on what cassette was used, what exact time it was recorded, what 
matchers were in use, the cassette options, and how a request is handled.

```
[Cassette: 'testing'] - 2020-11-24 16:05:17 - Init. HTTPInteractionList w/ request matchers [method, uri] & 0 interaction(s): {  }
[Cassette: 'testing'] - 2020-11-24 16:05:17 - Initialized with options: {name: testing, record: once, serialize_with: yaml, persist_with: FileSystem, match_requests_on: c("method", "uri"), update_content_length_header: FALSE, allow_playback_repeats: FALSE, preserve_exact_body_bytes: FALSE}
[Cassette: 'testing'] - 2020-11-24 16:05:17 - Handling request: head https://httpbin.org/get (disabled: FALSE)
[Cassette: 'testing'] - 2020-11-24 16:05:17 - Identified request type: (recordable) for head https://httpbin.org/get
[Cassette: 'testing'] - 2020-11-24 16:05:17 -    Recorded HTTP interaction: head https://httpbin.org/get => 200 
```

Logging isn't meant to be turned on all the time - rather only for debugging/informational purposes.

### Return to normal development

Make sure you ejected the cassette you were using!

Unless your vcr helper/setup file tweaked more things than you would like, 
you do not even need to re-start R, but you could, just to be on the safe side.

Sometimes you may need to turn off `vcr`, either for individual
function calls, individual test blocks, whole test files, or
for the entire package. The following attempts to break down
all the options.

`vcr` has the following four exported functions:

- `turned_off()` - Turns vcr off for the duration of a code block
- `turn_off()` - Turns vcr off completely, so that it no longer handles every
HTTP request
- `turn_on()` - turns vcr on; the opposite of `turn_off()`
- `turned_on()` - Asks if vcr is turned on, returns a boolean

Instead of using the above four functions, you could use environment
variables to achieve the same thing. This way you could enable/disable
`vcr` in non-interactive environments such as continuous integration,
Docker containers, or running R non-interactively from the command line.
The full set of environment variables `vcr` uses, all of which accept
only `TRUE` or `FALSE`:

- `VCR_TURN_OFF`: turn off vcr altogether; set to `TRUE` to skip any vcr
usage; default: `FALSE`
- `VCR_TURNED_OFF`: set the `turned_off` internal package setting; this
does not turn off vcr completely as does `VCR_TURN_OFF` does, but rather
is looked at together with `VCR_IGNORE_CASSETTES`
- `VCR_IGNORE_CASSETTES`: set the `ignore_cassettes` internal package
setting; this is looked at together with `VCR_TURNED_OFF`

## turned_off {#turned-off}

`turned_off()` lets you temporarily make a real HTTP request without completely turning
`vcr` off, unloading it, etc.

What happens internally is we turn off `vcr`, run your code block, then on exit
turn `vcr` back on - such that `vcr` is only turned off for the duration of your
code block. Even if your code block errors, `vcr` will be turned back on
due to use of `on.exit(turn_on())`

```r
library(vcr)
library(crul)
turned_off({
  con <- HttpClient$new(url = "https://httpbin.org/get")
  con$get()
})
```

```r
#> <crul response>
#>   url: https://httpbin.org/get
#>   request_headers:
#>     User-Agent: libcurl/7.54.0 r-curl/4.3 crul/0.9.0
#>     Accept-Encoding: gzip, deflate
#>     Accept: application/json, text/xml, application/xml, */*
#>   response_headers:
#>     status: HTTP/1.1 200 OK
#>     date: Fri, 14 Feb 2020 19:44:46 GMT
#>     content-type: application/json
#>     content-length: 365
#>     connection: keep-alive
#>     server: gunicorn/19.9.0
#>     access-control-allow-origin: *
#>     access-control-allow-credentials: true
#>   status: 200
```

## turn_off/turn_on {#turn-off-on}

`turn_off()` is different from `turned_off()` in that `turn_off()` is not aimed
at a single call block, but rather it turns `vcr` off for the entire package.
`turn_off()` does check first before turning `vcr` off that there is not currently
a cassette in use. `turn_off()` is meant to make R ignore `vcr::insert_cassette()`
and `vcr::use_cassette()` blocks in your test suite - letting the code in the block
run as if they were not wrapped in `vcr` code - so that all you have to do to run
your tests with cached requests/responses AND with real HTTP requests is toggle
a single R function or environment variable.

```r
library(vcr)
vcr_configure(dir = tempdir())
# real HTTP request works - vcr is not engaged here
crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
# wrap HTTP request in use_cassette() - vcr is engaged here
use_cassette("foo_bar", {
  crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
})
# turn off & ignore cassettes - use_cassette is ignored, real HTTP request made
turn_off(ignore_cassettes = TRUE)
use_cassette("foo_bar", {
  crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
})
# if you turn off and don't ignore cassettes, error thrown
turn_off(ignore_cassettes = FALSE)
use_cassette("foo_bar", {
  res2=crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
})
# vcr back on - now use_cassette behaves as before
turn_on()
use_cassette("foo_bar3", {
  res2=crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
})
```

## turned_on {#turned-on}

`turned_on()` does what it says on the tin - it tells you if `vcr` is turned on
or not. 

```{r}
library(vcr)
turn_on()
turned_on()
turn_off()
turned_on()
```

## Environment variables {#lightswitch-env-vars}

The `VCR_TURN_OFF` environment variable can be used within R or on the command line
to turn off `vcr`. For example, you can run tests for a package that uses `vcr`, but
ignore any `use_cassette`/`insert_cassette` usage, by running this on the command
line in the root of your package:

```
VCR_TURN_OFF=true Rscript -e "devtools::test()"
```

Or, similarly within R:

```r
Sys.setenv(VCR_TURN_OFF = TRUE)
devtools::test()
```

The `VCR_TURNED_OFF` and `VCR_IGNORE_CASSETTES` environment variables can be used
in combination to achieve the same thing as `VCR_TURN_OFF`:

```
VCR_TURNED_OFF=true VCR_IGNORE_CASSETTES=true Rscript -e "devtools::test()"
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_handler-httr.R
\name{RequestHandlerHttr}
\alias{RequestHandlerHttr}
\title{RequestHandlerHttr}
\description{
Methods for the httr package, building on \link{RequestHandler}
}
\examples{
\dontrun{
vcr_configure(
 dir = tempdir(),
 record = "once"
)

# GET request
library(httr)
load("~/httr_req.rda")
req
x <- RequestHandlerHttr$new(req)
# x$handle()

# POST request
library(httr)
webmockr::httr_mock()
mydir <- file.path(tempdir(), "testing_httr")
invisible(vcr_configure(dir = mydir))
use_cassette(name = "testing2", {
  res <- POST("https://httpbin.org/post", body = list(foo = "bar"))
}, match_requests_on = c("method", "uri", "body"))

load("~/httr_req_post.rda")
insert_cassette("testing3")
httr_req_post
x <- RequestHandlerHttr$new(httr_req_post)
x
# x$handle()
self=x

}
}
\section{Super class}{
\code{\link[vcr:RequestHandler]{vcr::RequestHandler}} -> \code{RequestHandlerHttr}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{RequestHandlerHttr$new()}}
\item \href{#method-clone}{\code{RequestHandlerHttr$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="vcr" data-topic="RequestHandler" data-id="handle">}\href{../../vcr/html/RequestHandler.html#method-handle}{\code{vcr::RequestHandler$handle()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{RequestHandlerHttr} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestHandlerHttr$new(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{The request from an object of class \code{HttpInteraction}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{RequestHandlerHttr} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestHandlerHttr$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/as-casette.R
\name{as.cassette}
\alias{as.cassette}
\alias{as.cassettepath}
\title{Coerce names, etc. to cassettes}
\usage{
as.cassette(x, ...)

as.cassettepath(x)
}
\arguments{
\item{x}{Input, a cassette name (character), or something that
can be coerced to a cassette}

\item{...}{further arguments passed on to \code{\link[=cassettes]{cassettes()}} or
[read_cassette_meta()}
}
\value{
a cassette of class \code{Cassette}
}
\description{
Coerce names, etc. to cassettes

Coerce to a cassette path
}
\examples{
\dontrun{
vcr_configure(dir = tempfile())
insert_cassette("foobar")
cassettes(on_disk = FALSE)
cassettes(on_disk = TRUE)
as.cassette("foobar", on_disk = FALSE)
eject_cassette() # eject the current cassette

# cleanup
unlink(file.path(tempfile(), "foobar.yml"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logger.R
\name{vcr_logging}
\alias{vcr_logging}
\alias{vcr_log_file}
\alias{vcr_log_info}
\title{vcr log file setup}
\usage{
vcr_log_file(file, overwrite = TRUE)

vcr_log_info(message, include_date = TRUE)
}
\arguments{
\item{file}{(character) a file path, required}

\item{overwrite}{(logical) whether or not to overwrite the file at
'file' if it already exists. Default: \code{TRUE}}

\item{message}{(character) a message to log}

\item{include_date}{(logical) include date and time in each log entry.
Default: \code{FALSE}}
}
\description{
vcr log file setup
}
\examples{
# user workflow
vcr_configuration()
logfile <- file.path(tempdir(), "vcr.log")
vcr_configure(dir = tempdir(), log = TRUE, log_opts = list(file = logfile))

readLines(logfile) # empty

# log messages
vcr_log_info("hello world!")
readLines(logfile)
vcr_log_info("foo bar")
readLines(logfile)
## many messages
vcr_log_info(c("brown cow", "blue horse"))
readLines(logfile)
vcr_log_info(c("brown cow", "blue horse", "green goat"))
readLines(logfile)

# standalone workflow
# set a file to log to
vcr_log_file((f <- tempfile()))
readLines(f) # empty

# log messages
vcr_log_info("hello world!")
readLines(logfile)
vcr_log_info("foo bar")
readLines(logfile)

# cleanup
unlink(f)
unlink(logfile)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/real_connections.R
\name{real_http_connections_allowed}
\alias{real_http_connections_allowed}
\title{Are real http connections allowed?}
\usage{
real_http_connections_allowed()
}
\value{
boolean, \code{TRUE} if real HTTP requests allowed; \code{FALSE} if not
}
\description{
Are real http connections allowed?
}
\examples{
real_http_connections_allowed()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serializers.R
\name{Serializers}
\alias{Serializers}
\alias{serializer_fetch}
\title{Cassette serializers}
\usage{
serializer_fetch(x = "yaml", name)
}
\description{
Keeps track of the cassette serializers in a hash-like object
}
\details{
\strong{Private Methods}
\describe{
\item{\code{serialize_get()}}{
Gets a named serializer. This is also run on \code{Serializers$new()}
}
}
}
\examples{
\dontrun{
(aa <- Serializers$new())
aa$name
aa$serializers
yaml_serializer <- aa$serializers$new()
yaml_serializer

x <- Serializers$new(name = "json")
x$serializers$new()
json_serializer <- x$serializers$new()
json_serializer
}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{serializers}}{(list) list of serializer names}

\item{\code{name}}{(character) Name of a serializer. "yaml" (default) or "json"}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Serializers$new()}}
\item \href{#method-clone}{\code{Serializers$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Serializers object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Serializers$new(serializers = list(), name = "yaml")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{serializers}}{(list) list of serializer names}

\item{\code{name}}{(character) Name of a serializer. "yaml" (default) or "json"}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Serializers} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Serializers$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/skip-vcr-off.R
\name{skip_if_vcr_off}
\alias{skip_if_vcr_off}
\title{Skip tests if vcr is off}
\usage{
skip_if_vcr_off()
}
\value{
Nothing, skip test.
}
\description{
Custom testthat skipper to skip tests if vcr is turned off via the
environment variable \code{VCR_TURN_OFF}.
}
\details{
This might be useful if your test will fail with real requests:
when the cassette was e.g. edited (a real request produced a 200 status code
but you made it a 502 status code for testing the behavior of your code
when the API errors)
or if the tests are very specific (e.g. testing a date was correctly parsed,
but making a real request would produce a different date).
}
\seealso{
\code{\link[=turn_off]{turn_off()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/response_class.R
\name{VcrResponse}
\alias{VcrResponse}
\title{The response of an HTTPInteraction}
\description{
Custom vcr http response object
}
\examples{
\dontrun{
vcr_configure(dir = tempdir())

# basic example of VcrResponse use
url <- "https://google.com"
(cli <- crul::HttpClient$new(url = url))
(res <- cli$get("get", query = list(q = "stuff")))
(x <- VcrResponse$new(res$status_http(), res$response_headers,
   res$parse("UTF-8"), res$response_headers$status))
x$body
x$status
x$headers
x$http_version
x$to_hash()
x$from_hash(x$to_hash())

# update content length header
## example 1
### content-length header present, but no change
url <- "https://fishbase.ropensci.org"
cli <- crul::HttpClient$new(url = url, headers = list(`Accept-Encoding` = '*'))
res <- cli$get("species/34")
x <- VcrResponse$new(res$status_http(), res$response_headers,
   res$parse("UTF-8"), res$response_headers$status)
x$headers$`content-length`
x$update_content_length_header()
x$headers$`content-length`

## example 2
### no content-length header b/c a transfer-encoding header is included
### and no content-length header allowed if transfer-encoding header
### used (via rfc7230)
url <- "https://google.com"
cli <- crul::HttpClient$new(url = url)
res <- cli$get()
x <- VcrResponse$new(res$status_http(), res$response_headers,
   rawToChar(res$content), res$response_headers$status)
x$headers$`content-length` # = NULL
x$update_content_length_header() # no change, b/c header doesn't exist
x$headers$`content-length` # = NULL

## example 3
### content-length header present, and does change
body <- " Hello World "
x <- VcrResponse$new(200, list('content-length'=nchar(body)),
  body, "HTTP/2")
x$headers$`content-length` # = 13
x$body <- gsub("^\\\\s|\\\\s$", "", x$body)
x$headers$`content-length` # = 13
x$update_content_length_header()
x$headers$`content-length` # = 11

# check if body is compressed
url <- "https://fishbase.ropensci.org"
(cli <- crul::HttpClient$new(url = url))
(res <- cli$get("species/3"))
res$response_headers
(x <- VcrResponse$new(res$status_http(), res$response_headers,
   res$parse("UTF-8"), res$response_headers$status))
x$content_encoding()
x$is_compressed()

# with disk
url <- "https://google.com"
(cli <- crul::HttpClient$new(url = url))
f <- tempfile()
(res <- cli$get("get", query = list(q = "stuff"), disk = f))
(x <- VcrResponse$new(res$status_http(), res$response_headers,
   f, res$response_headers$status, disk = TRUE))
}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{status}}{the status of the response}

\item{\code{headers}}{the response headers}

\item{\code{body}}{the response body}

\item{\code{http_version}}{the HTTP version}

\item{\code{opts}}{a list}

\item{\code{adapter_metadata}}{Additional metadata used by a specific VCR adapter}

\item{\code{hash}}{a list}

\item{\code{disk}}{a boolean}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{VcrResponse$new()}}
\item \href{#method-to_hash}{\code{VcrResponse$to_hash()}}
\item \href{#method-from_hash}{\code{VcrResponse$from_hash()}}
\item \href{#method-update_content_length_header}{\code{VcrResponse$update_content_length_header()}}
\item \href{#method-get_header}{\code{VcrResponse$get_header()}}
\item \href{#method-edit_header}{\code{VcrResponse$edit_header()}}
\item \href{#method-delete_header}{\code{VcrResponse$delete_header()}}
\item \href{#method-content_encoding}{\code{VcrResponse$content_encoding()}}
\item \href{#method-is_compressed}{\code{VcrResponse$is_compressed()}}
\item \href{#method-clone}{\code{VcrResponse$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new VcrResponse object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$new(
  status,
  headers,
  body,
  http_version,
  opts,
  adapter_metadata = NULL,
  disk
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{status}}{the status of the response}

\item{\code{headers}}{the response headers}

\item{\code{body}}{the response body}

\item{\code{http_version}}{the HTTP version}

\item{\code{opts}}{a list}

\item{\code{adapter_metadata}}{Additional metadata used by a specific VCR adapter}

\item{\code{disk}}{boolean, is body a file on disk}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{VcrResponse} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_hash"></a>}}
\if{latex}{\out{\hypertarget{method-to_hash}{}}}
\subsection{Method \code{to_hash()}}{
Create a hash
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$to_hash()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-from_hash"></a>}}
\if{latex}{\out{\hypertarget{method-from_hash}{}}}
\subsection{Method \code{from_hash()}}{
Get a hash back to an R list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$from_hash(hash)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{hash}}{a list}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
an \code{VcrResponse} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-update_content_length_header"></a>}}
\if{latex}{\out{\hypertarget{method-update_content_length_header}{}}}
\subsection{Method \code{update_content_length_header()}}{
Updates the Content-Length response header so that
it is accurate for the response body
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$update_content_length_header()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
no return; modifies the content length header
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_header"></a>}}
\if{latex}{\out{\hypertarget{method-get_header}{}}}
\subsection{Method \code{get_header()}}{
Get a header by name
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$get_header(key)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{key}}{(character) header name to get}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
the header value (if it exists)
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-edit_header"></a>}}
\if{latex}{\out{\hypertarget{method-edit_header}{}}}
\subsection{Method \code{edit_header()}}{
Edit a header
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$edit_header(key, value = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{key}}{(character) header name to edit}

\item{\code{value}}{(character) new value to assign}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; modifies the header in place
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-delete_header"></a>}}
\if{latex}{\out{\hypertarget{method-delete_header}{}}}
\subsection{Method \code{delete_header()}}{
Delete a header
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$delete_header(key)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{key}}{(character) header name to delete}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; the header is deleted if it exists
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-content_encoding"></a>}}
\if{latex}{\out{\hypertarget{method-content_encoding}{}}}
\subsection{Method \code{content_encoding()}}{
Get the content-encoding header value
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$content_encoding()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(character) the content-encoding value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_compressed"></a>}}
\if{latex}{\out{\hypertarget{method-is_compressed}{}}}
\subsection{Method \code{is_compressed()}}{
Checks if the encoding is one of "gzip" or "deflate"
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$is_compressed()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{VcrResponse$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/request_class.R
\name{Request}
\alias{Request}
\title{The request of an HTTPInteraction}
\description{
object that handled all aspects of a request
}
\examples{
url <- "https://eu.httpbin.org/post"
body <- list(foo = "bar")
headers <- list(
  `User-Agent` = "libcurl/7.54.0 r-curl/3.2 crul/0.5.2",
  `Accept-Encoding` = "gzip, deflate",
  Accept = "application/json, text/xml, application/xml, */*"
)

(x <- Request$new("POST", url, body, headers))
x$body
x$method
x$uri
x$host
x$path
x$headers
h <- x$to_hash()
x$from_hash(h)
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{method}}{(character) http method}

\item{\code{uri}}{(character) a uri}

\item{\code{scheme}}{(character) scheme (http or https)}

\item{\code{host}}{(character) host (e.g., stuff.org)}

\item{\code{path}}{(character) path (e.g., foo/bar)}

\item{\code{query}}{(character) query params, named list}

\item{\code{body}}{(character) named list}

\item{\code{headers}}{(character) named list}

\item{\code{skip_port_stripping}}{(logical) whether to strip thhe port}

\item{\code{hash}}{(character) a named list - internal use}

\item{\code{opts}}{(character) options - internal use}

\item{\code{disk}}{(logical) xx}

\item{\code{fields}}{(various) request body details}

\item{\code{output}}{(various) request output details, disk, memory, etc}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Request$new()}}
\item \href{#method-to_hash}{\code{Request$to_hash()}}
\item \href{#method-from_hash}{\code{Request$from_hash()}}
\item \href{#method-clone}{\code{Request$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Request} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Request$new(method, uri, body, headers, opts, disk, fields, output)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{method}}{(character) the HTTP method (i.e. head, options, get,
post, put, patch or delete)}

\item{\code{uri}}{(character) request URI}

\item{\code{body}}{(character) request body}

\item{\code{headers}}{(named list) request headers}

\item{\code{opts}}{(named list) options internal use}

\item{\code{disk}}{(boolean), is body a file on disk}

\item{\code{fields}}{(various) post fields}

\item{\code{output}}{(various) output details}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Request} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_hash"></a>}}
\if{latex}{\out{\hypertarget{method-to_hash}{}}}
\subsection{Method \code{to_hash()}}{
Convert the request to a list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Request$to_hash()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-from_hash"></a>}}
\if{latex}{\out{\hypertarget{method-from_hash}{}}}
\subsection{Method \code{from_hash()}}{
Convert the request to a list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Request$from_hash(hash)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{hash}}{a list}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new \code{Request} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Request$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/lightswitch.R
\name{lightswitch}
\alias{lightswitch}
\alias{turned_off}
\alias{turn_on}
\alias{turned_on}
\alias{turn_off}
\title{Turn vcr on and off, check on/off status, and turn off for a given http call}
\usage{
turned_off(..., ignore_cassettes = FALSE)

turn_on()

turned_on()

turn_off(ignore_cassettes = FALSE)
}
\arguments{
\item{...}{Any block of code to run, presumably an http request}

\item{ignore_cassettes}{(logical) Controls what happens when a cassette is
inserted while vcr is turned off. If \code{TRUE} is passed, the cassette
insertion will be ignored; otherwise an error will be raised.
Default: \code{FALSE}}
}
\description{
Turn vcr on and off, check on/off status, and turn off for a given http call
}
\details{
Sometimes you may need to turn off \code{vcr}, either for individual function
calls, individual test blocks, whole test files, or for the entire
package. The following attempts to break down all the options.

\code{vcr} has the following four exported functions:
\itemize{
\item \code{turned_off()} - Turns vcr off for the duration of a code block
\item \code{turn_off()} - Turns vcr off completely, so that it no longer
handles every HTTP request
\item \code{turn_on()} - turns vcr on; the opposite of \code{turn_off()}
\item \code{turned_on()} - Asks if vcr is turned on, returns a boolean
}

Instead of using the above four functions, you could use environment
variables to achieve the same thing. This way you could enable/disable
\code{vcr} in non-interactive environments such as continuous integration,
Docker containers, or running R non-interactively from the command line.
The full set of environment variables \code{vcr} uses, all of which accept
only \code{TRUE} or \code{FALSE}:
\itemize{
\item \code{VCR_TURN_OFF}: turn off vcr altogether; set to \code{TRUE} to skip any
vcr usage; default: \code{FALSE}
\item \code{VCR_TURNED_OFF}: set the \code{turned_off} internal package setting;
this does not turn off vcr completely as does \code{VCR_TURN_OFF} does,
but rather is looked at together with \code{VCR_IGNORE_CASSETTES}
\item \code{VCR_IGNORE_CASSETTES}: set the \code{ignore_cassettes} internal package
setting; this is looked at together with \code{VCR_TURNED_OFF}
}
\subsection{turned_off}{

\code{turned_off()} lets you temporarily make a real HTTP request without
completely turning \code{vcr} off, unloading it, etc.

What happens internally is we turn off \code{vcr}, run your code block, then
on exit turn \code{vcr} back on - such that \code{vcr} is only turned off for the
duration of your code block. Even if your code block errors, \code{vcr} will
be turned back on due to use of \code{on.exit(turn_on())}\if{html}{\out{<div class="r">}}\preformatted{library(vcr)
library(crul)
turned_off(\{
  con <- HttpClient$new(url = "https://httpbin.org/get")
  con$get()
\})
}\if{html}{\out{</div>}}\if{html}{\out{<div class="r">}}\preformatted{#> <crul response>
#>   url: https://httpbin.org/get
#>   request_headers:
#>     User-Agent: libcurl/7.54.0 r-curl/4.3 crul/0.9.0
#>     Accept-Encoding: gzip, deflate
#>     Accept: application/json, text/xml, application/xml, */*
#>   response_headers:
#>     status: HTTP/1.1 200 OK
#>     date: Fri, 14 Feb 2020 19:44:46 GMT
#>     content-type: application/json
#>     content-length: 365
#>     connection: keep-alive
#>     server: gunicorn/19.9.0
#>     access-control-allow-origin: *
#>     access-control-allow-credentials: true
#>   status: 200
}\if{html}{\out{</div>}}
}

\subsection{turn_off/turn_on}{

\code{turn_off()} is different from \code{turned_off()} in that \code{turn_off()} is
not aimed at a single call block, but rather it turns \code{vcr} off for the
entire package. \code{turn_off()} does check first before turning \code{vcr} off
that there is not currently a cassette in use. \code{turn_off()} is meant to
make R ignore \code{vcr::insert_cassette()} and \code{vcr::use_cassette()} blocks
in your test suite - letting the code in the block run as if they were
not wrapped in \code{vcr} code - so that all you have to do to run your tests
with cached requests/responses AND with real HTTP requests is toggle a
single R function or environment variable.\if{html}{\out{<div class="r">}}\preformatted{library(vcr)
vcr_configure(dir = tempdir())
# real HTTP request works - vcr is not engaged here
crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
# wrap HTTP request in use_cassette() - vcr is engaged here
use_cassette("foo_bar", \{
  crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
\})
# turn off & ignore cassettes - use_cassette is ignored, real HTTP request made
turn_off(ignore_cassettes = TRUE)
use_cassette("foo_bar", \{
  crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
\})
# if you turn off and don't ignore cassettes, error thrown
turn_off(ignore_cassettes = FALSE)
use_cassette("foo_bar", \{
  res2=crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
\})
# vcr back on - now use_cassette behaves as before
turn_on()
use_cassette("foo_bar3", \{
  res2=crul::HttpClient$new(url = "https://eu.httpbin.org/get")$get()
\})
}\if{html}{\out{</div>}}
}

\subsection{turned_on}{

\code{turned_on()} does what it says on the tin - it tells you if \code{vcr} is
turned on or not.\if{html}{\out{<div class="r">}}\preformatted{library(vcr)
turn_on()
turned_on()
}\if{html}{\out{</div>}}\preformatted{## [1] TRUE
}\if{html}{\out{<div class="r">}}\preformatted{turn_off()
}\if{html}{\out{</div>}}\preformatted{## vcr turned off; see ?turn_on to turn vcr back on
}\if{html}{\out{<div class="r">}}\preformatted{turned_on()
}\if{html}{\out{</div>}}\preformatted{## [1] FALSE
}
}

\subsection{Environment variables}{

The \code{VCR_TURN_OFF} environment variable can be used within R or on the
command line to turn off \code{vcr}. For example, you can run tests for a
package that uses \code{vcr}, but ignore any \code{use_cassette}/\code{insert_cassette}
usage, by running this on the command line in the root of your package:\preformatted{VCR_TURN_OFF=true Rscript -e "devtools::test()"
}

Or, similarly within R:\if{html}{\out{<div class="r">}}\preformatted{Sys.setenv(VCR_TURN_OFF = TRUE)
devtools::test()
}\if{html}{\out{</div>}}

The \code{VCR_TURNED_OFF} and \code{VCR_IGNORE_CASSETTES} environment variables
can be used in combination to achieve the same thing as \code{VCR_TURN_OFF}:\preformatted{VCR_TURNED_OFF=true VCR_IGNORE_CASSETTES=true Rscript -e "devtools::test()"
}
}
}
\examples{
\dontrun{
vcr_configure(dir = tempdir())

turn_on()
turned_on()
turn_off()

# turn off for duration of a block
library(crul)
turned_off({
 res <- HttpClient$new(url = "https://eu.httpbin.org/get")$get()
})
res

# turn completely off
turn_off()
library(webmockr)
crul::mock()
# HttpClient$new(url = "https://eu.httpbin.org/get")$get(verbose = TRUE)
turn_on()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/persisters-file.R
\name{FileSystem}
\alias{FileSystem}
\title{File system persister}
\description{
The only built-in cassette persister. Persists cassettes
to the file system.
}
\details{
\strong{Private Methods}
\describe{
\item{\code{storage_location()}}{
Get storage location
}
\item{\code{absolute_path_to_file()}}{
Get absolute path to the \code{storage_location}
}
}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{file_name}}{(character) the file name, not whole path}

\item{\code{write_fxn}}{(character) fxn to use for writing to disk}

\item{\code{content}}{(character) content to record to a cassette}

\item{\code{path}}{(character) storage directory for cassettes}

\item{\code{write2disk}}{(character) write to disk or make a new FileSystem}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{FileSystem$new()}}
\item \href{#method-get_cassette}{\code{FileSystem$get_cassette()}}
\item \href{#method-is_empty}{\code{FileSystem$is_empty()}}
\item \href{#method-set_cassette}{\code{FileSystem$set_cassette()}}
\item \href{#method-clone}{\code{FileSystem$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{FileSystem} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileSystem$new(
  file_name = NULL,
  write_fxn = NULL,
  content = NULL,
  path = NULL,
  write2disk = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file_name}}{(character) the file name, not whole path}

\item{\code{write_fxn}}{(character) fxn to use for writing to disk}

\item{\code{content}}{(character) content to record to a cassette}

\item{\code{path}}{(character) storage directory for cassettes}

\item{\code{write2disk}}{(logical) write to disk or just make a new FileSystem
object. Default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{FileSystem} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_cassette"></a>}}
\if{latex}{\out{\hypertarget{method-get_cassette}{}}}
\subsection{Method \code{get_cassette()}}{
Gets the cassette for the given storage key (file name)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileSystem$get_cassette(file_name = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file_name}}{(character) the file name, not whole path}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
named list, from \code{yaml::yaml.load_file}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_empty"></a>}}
\if{latex}{\out{\hypertarget{method-is_empty}{}}}
\subsection{Method \code{is_empty()}}{
Checks if a cassette is empty or not
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileSystem$is_empty()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_cassette"></a>}}
\if{latex}{\out{\hypertarget{method-set_cassette}{}}}
\subsection{Method \code{set_cassette()}}{
Sets the cassette for the given storage key (file name)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileSystem$set_cassette(file_name = NULL, content)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file_name}}{(character) the file name, not whole path}

\item{\code{content}}{(character) content to record to a cassette}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; writes to disk
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FileSystem$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/request_handler-crul.R
\name{RequestHandlerCrul}
\alias{RequestHandlerCrul}
\title{RequestHandlerCrul}
\description{
Methods for the crul package, building on \link{RequestHandler}
}
\examples{
\dontrun{
vcr_configure(
 dir = tempdir(),
 record = "once"
)

data(crul_request)
crul_request$url$handle <- curl::new_handle()
crul_request
x <- RequestHandlerCrul$new(crul_request)
# x$handle()

# body matching
library(vcr)
library(crul)
vcr_configure(dir = tempdir(), log = TRUE,
 log_opts = list(file = file.path(tempdir(), "vcr.log")))
cli <- HttpClient$new(url = "https://httpbin.org")

## testing, same uri and method, changed body in 2nd block
use_cassette(name = "apple7", {
  resp <- cli$post("post", body = list(foo = "bar"))
}, match_requests_on = c("method", "uri", "body"))
## should error, b/c record="once"
if (interactive()) {
  use_cassette(name = "apple7", {
    resp <- cli$post("post", body = list(foo = "bar"))
    resp2 <- cli$post("post", body = list(hello = "world"))
  }, match_requests_on = c("method", "uri", "body"))
}
cas <- insert_cassette(name = "apple7", 
  match_requests_on = c("method", "uri", "body"))
resp2 <- cli$post("post", body = list(foo = "bar"))
eject_cassette("apple7")

## testing, same body, changed method in 2nd block
if (interactive()) {
use_cassette(name = "apple8", {
  x <- cli$post("post", body = list(hello = "world"))
}, match_requests_on = c("method", "body"))
use_cassette(name = "apple8", {
  x <- cli$get("post", body = list(hello = "world"))
}, match_requests_on = c("method", "body"))
}

## testing, same body, changed uri in 2nd block
# use_cassette(name = "apple9", {
#   x <- cli$post("post", body = list(hello = "world"))
#   w <- cli$post("get", body = list(hello = "world"))
# }, match_requests_on = c("method", "body"))
# use_cassette(name = "apple9", {
#   NOTHING HERE
# }, match_requests_on = c("method", "body"))
# unlink(file.path(vcr_configuration()$dir, "apple9.yml"))
}
}
\section{Super class}{
\code{\link[vcr:RequestHandler]{vcr::RequestHandler}} -> \code{RequestHandlerCrul}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-clone}{\code{RequestHandlerCrul$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="vcr" data-topic="RequestHandler" data-id="handle">}\href{../../vcr/html/RequestHandler.html#method-handle}{\code{vcr::RequestHandler$handle()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="vcr" data-topic="RequestHandler" data-id="initialize">}\href{../../vcr/html/RequestHandler.html#method-initialize}{\code{vcr::RequestHandler$initialize()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestHandlerCrul$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/eject_cassette.R
\name{eject_cassette}
\alias{eject_cassette}
\title{Eject a cassette}
\usage{
eject_cassette(
  cassette = NULL,
  options = list(),
  skip_no_unused_interactions_assertion = NULL
)
}
\arguments{
\item{cassette}{(character) a single cassette names to eject; see \code{name}
parameter definition in \code{\link[=insert_cassette]{insert_cassette()}} for cassette name rules}

\item{options}{(list) a list of options to apply to the eject process}

\item{skip_no_unused_interactions_assertion}{(logical) If \code{TRUE}, this will
skip the "no unused HTTP interactions" assertion enabled by the
\code{allow_unused_http_interactions = FALSE} cassette option. This is intended
for use when your test has had an error, but your test framework has
already handled it - IGNORED FOR NOW}
}
\value{
The ejected cassette if there was one
}
\description{
Eject a cassette
}
\examples{
vcr_configure(dir = tempdir())
insert_cassette("hello")
(x <- current_cassette())

# by default does current cassette
x <- eject_cassette()
x
# can also select by cassette name
# eject_cassette(cassette = "hello")
}
\seealso{
\code{\link[=use_cassette]{use_cassette()}}, \code{\link[=insert_cassette]{insert_cassette()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_matcher_registry.R
\name{RequestMatcherRegistry}
\alias{RequestMatcherRegistry}
\title{RequestMatcherRegistry}
\description{
handles request matchers
}
\note{
r1=from new request; r2=from recorded interaction
}
\examples{
\dontrun{
(x <- RequestMatcherRegistry$new())
x$default_matchers
x$registry
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{registry}}{initialze registry list with a request, or leave empty}

\item{\code{default_matchers}}{request matchers to use. default: method, uri}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{RequestMatcherRegistry$new()}}
\item \href{#method-register}{\code{RequestMatcherRegistry$register()}}
\item \href{#method-register_built_ins}{\code{RequestMatcherRegistry$register_built_ins()}}
\item \href{#method-try_to_register_body_as_json}{\code{RequestMatcherRegistry$try_to_register_body_as_json()}}
\item \href{#method-clone}{\code{RequestMatcherRegistry$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new RequestMatcherRegistry object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestMatcherRegistry$new(
  registry = list(),
  default_matchers = list("method", "uri")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{registry}}{initialze registry list with a request, or leave empty}

\item{\code{default_matchers}}{request matchers to use. default: method, uri}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{RequestMatcherRegistry} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-register"></a>}}
\if{latex}{\out{\hypertarget{method-register}{}}}
\subsection{Method \code{register()}}{
Register a custom matcher
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestMatcherRegistry$register(name, func)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{matcher name}

\item{\code{func}}{function that describes a matcher, should return
a single boolean}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; registers the matcher
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-register_built_ins"></a>}}
\if{latex}{\out{\hypertarget{method-register_built_ins}{}}}
\subsection{Method \code{register_built_ins()}}{
Register all built in matchers
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestMatcherRegistry$register_built_ins()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
no return; registers all built in matchers
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-try_to_register_body_as_json"></a>}}
\if{latex}{\out{\hypertarget{method-try_to_register_body_as_json}{}}}
\subsection{Method \code{try_to_register_body_as_json()}}{
Try to register body as JSON
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestMatcherRegistry$try_to_register_body_as_json(r1, r2)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{r1, r2}}{\link{Request} class objects}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; registers the matcher
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestMatcherRegistry$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/recording.R
\name{recording}
\alias{recording}
\title{vcr recording options}
\description{
vcr recording options
}
\details{
Record modes dictate under what circumstances http requests/responses
are recorded to cassettes (disk). Set the recording mode with the
parameter \code{record} in the \code{use_cassette()} and \code{insert_cassette()}
functions.
\subsection{once}{

The \code{once} record mode will:
\itemize{
\item Replay previously recorded interactions.
\item Record new interactions if there is no cassette file.
\item Cause an error to be raised for new requests if there is a cassette
file.
}

It is similar to the \code{new_episodes} record mode, but will prevent new,
unexpected requests from being made (i.e. because the request URI
changed or whatever).

\code{once} is the default record mode, used when you do not set one.
}

\subsection{none}{

The \code{none} record mode will:
\itemize{
\item Replay previously recorded interactions.
\item Cause an error to be raised for any new requests.
}

This is useful when your code makes potentially dangerous HTTP requests.
The \code{none} record mode guarantees that no new HTTP requests will be
made.
}

\subsection{new_episodes}{

The \code{new_episodes} record mode will:
\itemize{
\item Record new interactions.
\item Replay previously recorded interactions.
}

It is similar to the \code{once} record mode, but will \strong{always} record new
interactions, even if you have an existing recorded one that is similar
(but not identical, based on the \code{match_request_on} option).
}

\subsection{all}{

The \code{all} record mode will:
\itemize{
\item Record new interactions.
\item Never replay previously recorded interactions.
}

This can be temporarily used to force vcr to re-record a cassette
(i.e. to ensure the responses are not out of date) or can be used when
you simply want to log all HTTP requests.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http_interaction_list.R
\name{HTTPInteractionList}
\alias{HTTPInteractionList}
\title{HTTPInteractionList class}
\description{
keeps track of all \link{HTTPInteraction} objects
}
\details{
\strong{Private Methods}
\describe{
\item{\code{has_unused_interactions()}}{
Are there any unused interactions? returns boolean
}
\item{\code{matching_interaction_index_for()}}{
asdfadf
}
\item{\code{matching_used_interaction_for(request)}}{
asdfadfs
}
\item{\code{interaction_matches_request(request, interaction)}}{
Check if a request matches an interaction (logical)
}
\item{\code{from_hash()}}{
Get a hash back.
}
\item{\code{request_summary(z)}}{
Get a request summary (character)
}
\item{\code{response_summary(z)}}{
Get a response summary (character)
}
}
}
\examples{
\dontrun{
vcr_configure(
 dir = tempdir(),
 record = "once"
)

# make interactions
## make the request
### turn off mocking
crul::mock(FALSE)
url <- "https://eu.httpbin.org/post"
cli <- crul::HttpClient$new(url = url)
res <- cli$post(body = list(a = 5))

## request
(request <- Request$new("POST", url, list(a = 5), res$headers))
## response
(response <- VcrResponse$new(
   res$status_http(),
   res$response_headers,
   res$parse("UTF-8"),
   res$response_headers$status))
## make an interaction
(inter <- HTTPInteraction$new(request = request, response = response))

# make an interactionlist
(x <- HTTPInteractionList$new(
   interactions = list(inter),
   request_matchers = vcr_configuration()$match_requests_on
))
x$interactions
x$request_matchers
x$parent_list
x$parent_list$response_for()
x$parent_list$has_interaction_matching()
x$parent_list$has_used_interaction_matching()
x$parent_list$remaining_unused_interaction_count()
x$used_interactions
x$allow_playback_repeats
x$interactions
x$response_for(request)
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{interactions}}{(list) list of interaction class objects}

\item{\code{request_matchers}}{(character) vector of request matchers}

\item{\code{allow_playback_repeats}}{whether to allow playback repeats}

\item{\code{parent_list}}{A list for empty objects, see \code{NullList}}

\item{\code{used_interactions}}{(list) Interactions that have been used}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{HTTPInteractionList$new()}}
\item \href{#method-response_for}{\code{HTTPInteractionList$response_for()}}
\item \href{#method-has_interaction_matching}{\code{HTTPInteractionList$has_interaction_matching()}}
\item \href{#method-has_used_interaction_matching}{\code{HTTPInteractionList$has_used_interaction_matching()}}
\item \href{#method-remaining_unused_interaction_count}{\code{HTTPInteractionList$remaining_unused_interaction_count()}}
\item \href{#method-assert_no_unused_interactions}{\code{HTTPInteractionList$assert_no_unused_interactions()}}
\item \href{#method-clone}{\code{HTTPInteractionList$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{HTTPInteractionList} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteractionList$new(
  interactions,
  request_matchers,
  allow_playback_repeats = FALSE,
  parent_list = NullList$new(),
  used_interactions = list()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{interactions}}{(list) list of interaction class objects}

\item{\code{request_matchers}}{(character) vector of request matchers}

\item{\code{allow_playback_repeats}}{whether to allow playback repeats or not}

\item{\code{parent_list}}{A list for empty objects, see \code{NullList}}

\item{\code{used_interactions}}{(list) Interactions that have been used. That is,
interactions that are on disk in the current cassette, and a
request has been made that matches that interaction}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{HTTPInteractionList} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-response_for"></a>}}
\if{latex}{\out{\hypertarget{method-response_for}{}}}
\subsection{Method \code{response_for()}}{
Check if there's a matching interaction, returns a
response object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteractionList$response_for(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{The request from an object of class \code{HTTPInteraction}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-has_interaction_matching"></a>}}
\if{latex}{\out{\hypertarget{method-has_interaction_matching}{}}}
\subsection{Method \code{has_interaction_matching()}}{
Check if has a matching interaction
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteractionList$has_interaction_matching(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{The request from an object of class \code{HTTPInteraction}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-has_used_interaction_matching"></a>}}
\if{latex}{\out{\hypertarget{method-has_used_interaction_matching}{}}}
\subsection{Method \code{has_used_interaction_matching()}}{
check if has used interactions matching a given request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteractionList$has_used_interaction_matching(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{The request from an object of class \code{HTTPInteraction}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-remaining_unused_interaction_count"></a>}}
\if{latex}{\out{\hypertarget{method-remaining_unused_interaction_count}{}}}
\subsection{Method \code{remaining_unused_interaction_count()}}{
Number of unused interactions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteractionList$remaining_unused_interaction_count()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
integer
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-assert_no_unused_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-assert_no_unused_interactions}{}}}
\subsection{Method \code{assert_no_unused_interactions()}}{
Checks if there are no unused interactions left.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteractionList$assert_no_unused_interactions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
various
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteractionList$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/request_response.R
\name{request_response}
\alias{request_response}
\alias{request_summary}
\alias{response_summary}
\title{request and response summary methods}
\usage{
request_summary(request, request_matchers = "")

response_summary(response)
}
\arguments{
\item{request}{a \link{Request} object}

\item{request_matchers}{(character) a vector of matchers.
Default: \code{""}}

\item{response}{a \link{VcrResponse} object}
}
\value{
character string, of either request or response
}
\description{
request and response summary methods
}
\details{
By default, method and uri are included
in the request summary - if body and/or headers are
specified in \code{request_matchers}, then they are also
included

HTTP status code and response body are included in the
response summary. The response body is truncated to a
max of 80 characters

In \code{response_summary()} we use \link{gsub} with \code{useBytes=TRUE} to avoid
problems sometimes seen with multibyte strings - this shouldn't affect
your data/etc. as this is only for printing a summary of the response
}
\examples{
# request
url <- "https://httpbin.org"
body <- list(foo = "bar")
headers <- list(
  `User-Agent` = "r-curl/3.2",
  `Accept-Encoding` = "gzip, deflate",
  Accept = "application/json"
)

(x <- Request$new("POST", url, body, headers))
request_summary(request = x)
request_summary(request = x, c('method', 'uri'))
request_summary(request = x, c('method', 'uri', 'body'))
request_summary(request = x, c('method', 'uri', 'headers'))
request_summary(request = x, c('method', 'uri', 'body', 'headers'))

# response
status <- list(status_code = 200, message = "OK",
  explanation = "Request fulfilled, document follows")
headers <- list(
  status = "HTTP/1.1 200 OK",
  connection = "keep-alive",
  date = "Tue, 24 Apr 2018 04:46:56 GMT"
)
response_body <- 
"{\"args\": {\"q\": \"stuff\"}, \"headers\": {\"Accept\": \"text/html\"}}\n"
(x <- VcrResponse$new(status, headers,
   response_body, "HTTP/1.1 200 OK"))
response_summary(x)

## with binary body
# path <- "tests/testthat/png_eg.rda"
# load(path)
# (x <- VcrResponse$new(status, headers, png_eg, "HTTP/1.1 200 OK"))
# response_summary(x)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cassette_class.R
\name{Cassette}
\alias{Cassette}
\title{Cassette handler}
\value{
an object of class \code{Cassette}
}
\description{
Main R6 class that is called from the main user facing
function \code{\link[=use_cassette]{use_cassette()}}
}
\section{Points of webmockr integration}{

\itemize{
\item \code{initialize()}: webmockr is used in the \code{initialize()} method to
create webmockr stubs. stubs are created on call to \code{Cassette$new()}
within \code{insert_cassette()}, but then on exiting \code{use_cassette()},
or calling \code{eject()} on \code{Cassette} class from \code{insert_cassette()},
stubs are cleaned up.
\item \code{eject()} method: \code{\link[webmockr:enable]{webmockr::disable()}} is called before exiting
eject to disable webmock so that webmockr does not affect any HTTP
requests that happen afterwards
\item \code{call_block()} method: call_block is used in the \code{\link[=use_cassette]{use_cassette()}}
function to evaluate whatever code is passed to it; within call_block
\code{\link[webmockr:webmockr_configure]{webmockr::webmockr_allow_net_connect()}} is run before we evaluate
the code block to allow real HTTP requests, then
\code{\link[webmockr:webmockr_configure]{webmockr::webmockr_disable_net_connect()}} is called after evalulating
the code block to disallow real HTTP requests
\item \code{make_http_interaction()} method: \code{\link[webmockr:pluck_body]{webmockr::pluck_body()}} utility
function is used to pull the request body out of the HTTP request
\item \code{serialize_to_crul()} method: method: \link[webmockr:RequestSignature]{webmockr::RequestSignature} and
\link[webmockr:Response]{webmockr::Response} are used to build a request and response,
respectively, then passed to \code{\link[webmockr:build_crul_response]{webmockr::build_crul_response()}}
to make a complete \code{crul} HTTP response object
}
}

\examples{
library(vcr)
vcr_configure(dir = tempdir())

res <- Cassette$new(name = "bob")
res$file()
res$originally_recorded_at()
res$recording()
res$serializable_hash()
res$eject()
res$should_remove_matching_existing_interactions()
res$storage_key()
res$match_requests_on

# record all requests
res <- Cassette$new("foobar", record = "all")
res$eject()

# cleanup
unlink(file.path(tempdir(), c("bob.yml", "foobar.yml")))

library(vcr)
vcr_configure(dir = tempdir())
res <- Cassette$new(name = "jane")
library(crul)
HttpClient$new("https://httpbin.org")$get("get")
}
\seealso{
\code{\link[=vcr_configure]{vcr_configure()}}, \code{\link[=use_cassette]{use_cassette()}}, \code{\link[=insert_cassette]{insert_cassette()}}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{name}}{(character) cassette name}

\item{\code{record}}{(character) record mode}

\item{\code{manfile}}{(character) cassette file path}

\item{\code{recorded_at}}{(character) date/time recorded at}

\item{\code{serialize_with}}{(character) serializer to use (yaml|json)}

\item{\code{serializer}}{(character) serializer to use (yaml|json)}

\item{\code{persist_with}}{(character) persister to use (FileSystem only)}

\item{\code{persister}}{(character) persister to use (FileSystem only)}

\item{\code{match_requests_on}}{(character) matchers to use
default: method & uri}

\item{\code{re_record_interval}}{(numeric) the re-record interval}

\item{\code{tag}}{ignored, not used right now}

\item{\code{tags}}{ignored, not used right now}

\item{\code{root_dir}}{root dir, gathered from \code{\link[=vcr_configuration]{vcr_configuration()}}}

\item{\code{update_content_length_header}}{(logical) Whether to overwrite the
\code{Content-Length} header}

\item{\code{allow_playback_repeats}}{(logical) Whether to allow a single HTTP
interaction to be played back multiple times}

\item{\code{allow_unused_http_interactions}}{(logical) ignored, not used right now}

\item{\code{exclusive}}{(logical) ignored, not used right now}

\item{\code{preserve_exact_body_bytes}}{(logical) Whether to base64 encode the
bytes of the requests and responses}

\item{\code{args}}{(list) internal use}

\item{\code{http_interactions_}}{(list) internal use}

\item{\code{new_recorded_interactions}}{(list) internal use}

\item{\code{clean_outdated_http_interactions}}{(logical) Should outdated interactions
be recorded back to file}

\item{\code{to_return}}{(logical) internal use}

\item{\code{cassette_opts}}{(list) various cassette options}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Cassette$new()}}
\item \href{#method-print}{\code{Cassette$print()}}
\item \href{#method-call_block}{\code{Cassette$call_block()}}
\item \href{#method-eject}{\code{Cassette$eject()}}
\item \href{#method-file}{\code{Cassette$file()}}
\item \href{#method-recording}{\code{Cassette$recording()}}
\item \href{#method-is_empty}{\code{Cassette$is_empty()}}
\item \href{#method-originally_recorded_at}{\code{Cassette$originally_recorded_at()}}
\item \href{#method-serializable_hash}{\code{Cassette$serializable_hash()}}
\item \href{#method-interactions_to_record}{\code{Cassette$interactions_to_record()}}
\item \href{#method-merged_interactions}{\code{Cassette$merged_interactions()}}
\item \href{#method-up_to_date_interactions}{\code{Cassette$up_to_date_interactions()}}
\item \href{#method-should_re_record}{\code{Cassette$should_re_record()}}
\item \href{#method-should_stub_requests}{\code{Cassette$should_stub_requests()}}
\item \href{#method-should_remove_matching_existing_interactions}{\code{Cassette$should_remove_matching_existing_interactions()}}
\item \href{#method-storage_key}{\code{Cassette$storage_key()}}
\item \href{#method-raw_cassette_bytes}{\code{Cassette$raw_cassette_bytes()}}
\item \href{#method-make_dir}{\code{Cassette$make_dir()}}
\item \href{#method-deserialized_hash}{\code{Cassette$deserialized_hash()}}
\item \href{#method-previously_recorded_interactions}{\code{Cassette$previously_recorded_interactions()}}
\item \href{#method-write_recorded_interactions_to_disk}{\code{Cassette$write_recorded_interactions_to_disk()}}
\item \href{#method-record_http_interaction}{\code{Cassette$record_http_interaction()}}
\item \href{#method-any_new_recorded_interactions}{\code{Cassette$any_new_recorded_interactions()}}
\item \href{#method-make_args}{\code{Cassette$make_args()}}
\item \href{#method-write_metadata}{\code{Cassette$write_metadata()}}
\item \href{#method-http_interactions}{\code{Cassette$http_interactions()}}
\item \href{#method-make_http_interaction}{\code{Cassette$make_http_interaction()}}
\item \href{#method-serialize_to_crul}{\code{Cassette$serialize_to_crul()}}
\item \href{#method-clone}{\code{Cassette$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Cassette} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$new(
  name,
  record,
  serialize_with,
  persist_with,
  match_requests_on,
  re_record_interval,
  tag,
  tags,
  update_content_length_header,
  allow_playback_repeats,
  allow_unused_http_interactions,
  exclusive,
  preserve_exact_body_bytes,
  clean_outdated_http_interactions
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{The name of the cassette. vcr will sanitize this to ensure it
is a valid file name.}

\item{\code{record}}{The record mode. Default: "once". In the future we'll support
"once", "all", "none", "new_episodes". See \link{recording} for more information}

\item{\code{serialize_with}}{(character) Which serializer to use.
Valid values are "yaml" (default), the only one supported for now.}

\item{\code{persist_with}}{(character) Which cassette persister to
use. Default: "file_system". You can also register and use a
custom persister.}

\item{\code{match_requests_on}}{List of request matchers
to use to determine what recorded HTTP interaction to replay. Defaults to
\verb{["method", "uri"]}. The built-in matchers are "method", "uri",
"headers" and "body" ("host" and "path" not supported yet, but should
be in a future version)}

\item{\code{re_record_interval}}{(numeric) When given, the cassette will be
re-recorded at the given interval, in seconds.}

\item{\code{tag, tags}}{tags ignored, not used right now}

\item{\code{update_content_length_header}}{(logical) Whether or
not to overwrite the \code{Content-Length} header of the responses to
match the length of the response body. Default: \code{FALSE}}

\item{\code{allow_playback_repeats}}{(logical) Whether or not to
allow a single HTTP interaction to be played back multiple times.
Default: \code{FALSE}.}

\item{\code{allow_unused_http_interactions}}{(logical) ignored, not used right now}

\item{\code{exclusive}}{(logical) ignored, not used right now}

\item{\code{preserve_exact_body_bytes}}{(logical) Whether or not
to base64 encode the bytes of the requests and responses for
this cassette when serializing it. See also \code{preserve_exact_body_bytes}
in \code{\link[=vcr_configure]{vcr_configure()}}. Default: \code{FALSE}}

\item{\code{clean_outdated_http_interactions}}{(logical) Should outdated interactions
be recorded back to file. Default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Cassette} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for \code{Cassette} objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$print(x, ...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-call_block"></a>}}
\if{latex}{\out{\hypertarget{method-call_block}{}}}
\subsection{Method \code{call_block()}}{
run code
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$call_block(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{pass in things to be evaluated}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
various
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eject"></a>}}
\if{latex}{\out{\hypertarget{method-eject}{}}}
\subsection{Method \code{eject()}}{
ejects the current cassette
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$eject()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
self
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-file"></a>}}
\if{latex}{\out{\hypertarget{method-file}{}}}
\subsection{Method \code{file()}}{
get the file path for the cassette
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$file()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-recording"></a>}}
\if{latex}{\out{\hypertarget{method-recording}{}}}
\subsection{Method \code{recording()}}{
is the cassette in recording mode?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$recording()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_empty"></a>}}
\if{latex}{\out{\hypertarget{method-is_empty}{}}}
\subsection{Method \code{is_empty()}}{
is the cassette on disk empty
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$is_empty()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-originally_recorded_at"></a>}}
\if{latex}{\out{\hypertarget{method-originally_recorded_at}{}}}
\subsection{Method \code{originally_recorded_at()}}{
timestamp the cassette was originally recorded at
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$originally_recorded_at()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
POSIXct date
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serializable_hash"></a>}}
\if{latex}{\out{\hypertarget{method-serializable_hash}{}}}
\subsection{Method \code{serializable_hash()}}{
Get a list of the http interactions to record + recorded_with
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$serializable_hash()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-interactions_to_record"></a>}}
\if{latex}{\out{\hypertarget{method-interactions_to_record}{}}}
\subsection{Method \code{interactions_to_record()}}{
Get the list of http interactions to record
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$interactions_to_record()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-merged_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-merged_interactions}{}}}
\subsection{Method \code{merged_interactions()}}{
Get interactions to record
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$merged_interactions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-up_to_date_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-up_to_date_interactions}{}}}
\subsection{Method \code{up_to_date_interactions()}}{
Cleans out any old interactions based on the
re_record_interval and clean_outdated_http_interactions settings
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$up_to_date_interactions(interactions)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{interactions}}{list of http interactions, of class \link{HTTPInteraction}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
list of interactions to record
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-should_re_record"></a>}}
\if{latex}{\out{\hypertarget{method-should_re_record}{}}}
\subsection{Method \code{should_re_record()}}{
Should re-record interactions?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$should_re_record()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-should_stub_requests"></a>}}
\if{latex}{\out{\hypertarget{method-should_stub_requests}{}}}
\subsection{Method \code{should_stub_requests()}}{
Is record mode NOT "all"?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$should_stub_requests()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-should_remove_matching_existing_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-should_remove_matching_existing_interactions}{}}}
\subsection{Method \code{should_remove_matching_existing_interactions()}}{
Is record mode "all"?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$should_remove_matching_existing_interactions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-storage_key"></a>}}
\if{latex}{\out{\hypertarget{method-storage_key}{}}}
\subsection{Method \code{storage_key()}}{
Get the serializer path
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$storage_key()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-raw_cassette_bytes"></a>}}
\if{latex}{\out{\hypertarget{method-raw_cassette_bytes}{}}}
\subsection{Method \code{raw_cassette_bytes()}}{
Get character string of entire cassette; bytes is a misnomer
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$raw_cassette_bytes()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-make_dir"></a>}}
\if{latex}{\out{\hypertarget{method-make_dir}{}}}
\subsection{Method \code{make_dir()}}{
Create the directory that holds the cassettes, if not present
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$make_dir()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
no return; creates a directory recursively, if missing
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialized_hash"></a>}}
\if{latex}{\out{\hypertarget{method-deserialized_hash}{}}}
\subsection{Method \code{deserialized_hash()}}{
get http interactions from the cassette via the serializer
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$deserialized_hash()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-previously_recorded_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-previously_recorded_interactions}{}}}
\subsection{Method \code{previously_recorded_interactions()}}{
get all previously recorded interactions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$previously_recorded_interactions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-write_recorded_interactions_to_disk"></a>}}
\if{latex}{\out{\hypertarget{method-write_recorded_interactions_to_disk}{}}}
\subsection{Method \code{write_recorded_interactions_to_disk()}}{
write recorded interactions to disk
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$write_recorded_interactions_to_disk()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-record_http_interaction"></a>}}
\if{latex}{\out{\hypertarget{method-record_http_interaction}{}}}
\subsection{Method \code{record_http_interaction()}}{
record an http interaction (doesn't write to disk)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$record_http_interaction(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{an crul or httr response object, with the request at \verb{$request}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-any_new_recorded_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-any_new_recorded_interactions}{}}}
\subsection{Method \code{any_new_recorded_interactions()}}{
Are there any new recorded interactions?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$any_new_recorded_interactions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-make_args"></a>}}
\if{latex}{\out{\hypertarget{method-make_args}{}}}
\subsection{Method \code{make_args()}}{
make list of all options
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$make_args()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-write_metadata"></a>}}
\if{latex}{\out{\hypertarget{method-write_metadata}{}}}
\subsection{Method \code{write_metadata()}}{
write metadata to the cassette
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$write_metadata()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-http_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-http_interactions}{}}}
\subsection{Method \code{http_interactions()}}{
make \link{HTTPInteractionList} object, assign to http_interactions_ var
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$http_interactions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
nothing returned
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-make_http_interaction"></a>}}
\if{latex}{\out{\hypertarget{method-make_http_interaction}{}}}
\subsection{Method \code{make_http_interaction()}}{
Make an \code{HTTPInteraction} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$make_http_interaction(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{an crul or httr response object, with the request at \verb{$request}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
an object of class \link{HTTPInteraction}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_to_crul"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_to_crul}{}}}
\subsection{Method \code{serialize_to_crul()}}{
Make a crul response object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$serialize_to_crul()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a crul response
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cassette$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/serializers-json.R
\name{JSON}
\alias{JSON}
\title{The JSON serializer}
\description{
class with methods for serializing via \pkg{jsonlite}
}
\keyword{internal}
\section{Super class}{
\code{\link[vcr:Serializer]{vcr::Serializer}} -> \code{JSON}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{JSON$new()}}
\item \href{#method-serialize}{\code{JSON$serialize()}}
\item \href{#method-deserialize}{\code{JSON$deserialize()}}
\item \href{#method-clone}{\code{JSON$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{JSON} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JSON$new(path = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) full path to the yaml file}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{JSON} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize"></a>}}
\if{latex}{\out{\hypertarget{method-serialize}{}}}
\subsection{Method \code{serialize()}}{
Serializes the given hash using internal fxn write_json
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JSON$serialize(x, path, bytes)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{(list) the object to serialize}

\item{\code{path}}{(character) the file path}

\item{\code{bytes}}{(logical) whether to preserve exact body bytes or not}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
(character) the json string to write to disk
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize}{}}}
\subsection{Method \code{deserialize()}}{
Deserializes the content at the file path using
jsonlite::fromJSON
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JSON$deserialize()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(list) the deserialized object, an R list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{JSON$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/serializer.R
\name{Serializer}
\alias{Serializer}
\title{Serializer class - base class for JSON/YAML serializers}
\description{
Serializer class - base class for JSON/YAML serializers

Serializer class - base class for JSON/YAML serializers
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{file_extension}}{(character) A file extension}

\item{\code{path}}{(character) full path to the yaml file}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Serializer$new()}}
\item \href{#method-serialize}{\code{Serializer$serialize()}}
\item \href{#method-deserialize}{\code{Serializer$deserialize()}}
\item \href{#method-clone}{\code{Serializer$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new YAML object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Serializer$new(file_extension = NULL, path = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file_extension}}{(character) A file extension}

\item{\code{path}}{(character) path to the cassette, excluding the cassette
directory and the file extension}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{YAML} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize"></a>}}
\if{latex}{\out{\hypertarget{method-serialize}{}}}
\subsection{Method \code{serialize()}}{
Serializes a hash - REPLACED BY YAML/JSON METHODS
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Serializer$serialize(x, path, bytes)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{(list) the object to serialize}

\item{\code{path}}{(character) the file path}

\item{\code{bytes}}{(logical) whether to preserve exact body bytes or not}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
(character) the YAML or JSON string to write to disk
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize}{}}}
\subsection{Method \code{deserialize()}}{
Serializes a file - REPLACED BY YAML/JSON METHODS
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Serializer$deserialize()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Serializer$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/vcr-package.R
\docType{package}
\name{vcr}
\alias{vcr}
\alias{vcr-package}
\title{vcr: Record HTTP Calls to Disk}
\description{
\pkg{vcr} records test suite 'HTTP' requests and replay them during
future runs.
}
\details{
Check out the \href{https://books.ropensci.org/http-testing/}{http testing book}
for a lot more documentation on \code{vcr}, \code{webmockr}, and \code{crul}
}
\section{Backstory}{

A Ruby gem of the same name (\code{VCR}, \url{https://github.com/vcr/vcr}) was
created many years ago and is the original. Ports in many languages
have been done. Check out that GitHub repo for all the details on
how the canonical version works.
}

\section{Main functions}{

The \link{use_cassette} function is most likely what you'll want to use. It
sets the cassette you want to record to, inserts the cassette, and then
ejects the cassette, recording the interactions to the cassette.

Instead, you can use \link{insert_cassette}, but then you have to make sure
to use \link{eject_cassette}.
}

\section{vcr configuration}{

\link{vcr_configure} is the function to use to set R session wide settings.
See it's manual file for help.
}

\section{Record modes}{

See \link{recording} for help on record modes.
}

\section{Request matching}{

See \link{request-matching} for help on the many request matching options.
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http_interactions.R
\name{http_interactions}
\alias{http_interactions}
\title{Get the http interactions of the current cassette}
\usage{
http_interactions()
}
\value{
object of class \code{HTTPInteractionList} if there is a current
cassette in use, or \code{NullList} if no cassette in use
}
\description{
Get the http interactions of the current cassette
}
\examples{
\dontrun{
vcr_configure(dir = tempdir())
insert_cassette("foo_bar")
webmockr::webmockr_allow_net_connect()
library(crul)
cli <- crul::HttpClient$new("https://eu.httpbin.org/get")
one <- cli$get(query = list(a = 5))
z <- http_interactions()
z
z$interactions
z$used_interactions
# on eject, request written to the cassette
eject_cassette("foo_bar")

# insert cassette again
insert_cassette("foo_bar")
# now interactions will be present 
z <- http_interactions()
z$interactions
z$used_interactions
invisible(cli$get(query = list(a = 5)))
z$used_interactions

# cleanup
unlink(file.path(tempdir(), "foo_bar.yml"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_ignorer.R
\name{RequestIgnorer}
\alias{RequestIgnorer}
\title{Request ignorer}
\description{
request ignorer methods
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{LOCALHOST_ALIASES}}{A constant with values: 'localhost', '127.0.0.1',
and '0.0.0.0'}

\item{\code{ignored_hosts}}{vector of ignored host URI's}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{RequestIgnorer$new()}}
\item \href{#method-ignore_request}{\code{RequestIgnorer$ignore_request()}}
\item \href{#method-ignore_localhost}{\code{RequestIgnorer$ignore_localhost()}}
\item \href{#method-ignore_localhost_value}{\code{RequestIgnorer$ignore_localhost_value()}}
\item \href{#method-ignore_hosts}{\code{RequestIgnorer$ignore_hosts()}}
\item \href{#method-should_be_ignored}{\code{RequestIgnorer$should_be_ignored()}}
\item \href{#method-clone}{\code{RequestIgnorer$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{RequestIgnorer} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestIgnorer$new()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A new \code{RequestIgnorer} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ignore_request"></a>}}
\if{latex}{\out{\hypertarget{method-ignore_request}{}}}
\subsection{Method \code{ignore_request()}}{
Will ignore any request for which the given function
returns \code{TRUE}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestIgnorer$ignore_request()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
no return; defines request ignorer hook
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ignore_localhost"></a>}}
\if{latex}{\out{\hypertarget{method-ignore_localhost}{}}}
\subsection{Method \code{ignore_localhost()}}{
ignore all localhost values (localhost, 127.0.0.1, 0.0.0.0)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestIgnorer$ignore_localhost()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
no return; sets to ignore all localhost aliases
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ignore_localhost_value"></a>}}
\if{latex}{\out{\hypertarget{method-ignore_localhost_value}{}}}
\subsection{Method \code{ignore_localhost_value()}}{
ignore a specific named localhost
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestIgnorer$ignore_localhost_value(value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{(character) A localhost value to ignore, e.g., 'localhost'}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; defines request ignorer hook
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ignore_hosts"></a>}}
\if{latex}{\out{\hypertarget{method-ignore_hosts}{}}}
\subsection{Method \code{ignore_hosts()}}{
ignore any named host
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestIgnorer$ignore_hosts(hosts)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{hosts}}{(character) vector of hosts to ignore}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; adds host to ignore
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-should_be_ignored"></a>}}
\if{latex}{\out{\hypertarget{method-should_be_ignored}{}}}
\subsection{Method \code{should_be_ignored()}}{
method to determine whether to ignore a request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestIgnorer$should_be_ignored(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{request to ignore}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; defines request ignorer hook
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestIgnorer$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/vcr_setup.R
\name{use_vcr}
\alias{use_vcr}
\title{Setup vcr for a package}
\usage{
use_vcr(dir = ".", verbose = TRUE)
}
\arguments{
\item{dir}{(character) path to package root. default's to
current directory}

\item{verbose}{(logical) print progress messages. default: \code{TRUE}}
}
\value{
only messages about progress, returns invisible()
}
\description{
Setup vcr for a package
}
\details{
Sets a mimimum vcr version, which is usually the latest
(stable) version on CRAN. You can of course easily remove or change
the version requirement yourself after running this function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcr_test_path.R
\name{vcr_test_path}
\alias{vcr_test_path}
\title{Locate file in tests directory}
\usage{
vcr_test_path(...)
}
\arguments{
\item{...}{Character vectors giving path component. each character string
gets added on to the path, e.g., \code{vcr_test_path("a", "b")} becomes
\code{tests/a/b} relative to the root of the package.}
}
\value{
A character vector giving the path
}
\description{
This function, similar to \code{testthat::test_path()}, is designed to work both
interactively and during tests, locating files in the \verb{tests/} directory.
}
\note{
Beware if you have more than one \code{tests} directories in your package
root. This may not work as intended if that is the case.
}
\examples{
if (interactive()) {
vcr_test_path("fixtures")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cassettes.R
\name{cassettes}
\alias{cassettes}
\alias{current_cassette}
\alias{cassette_path}
\title{List cassettes, get current cassette, etc.}
\usage{
cassettes(on_disk = TRUE, verb = FALSE)

current_cassette()

cassette_path()
}
\arguments{
\item{on_disk}{(logical) Check for cassettes on disk + cassettes in session
(\code{TRUE}), or check for only cassettes in session (\code{FALSE}). Default: \code{TRUE}}

\item{verb}{(logical) verbose messages}
}
\description{
List cassettes, get current cassette, etc.
}
\details{
\itemize{
\item \code{cassettes()}: returns cassettes found in your R session, you can toggle
whether we pull from those on disk or not
\item \code{current_cassette()}: returns an empty list when no cassettes are in use,
while it returns the current cassette (a \code{Cassette} object) when one is
in use
\item \code{cassette_path()}: just gives you the current directory path where
cassettes will be stored
}
}
\examples{
vcr_configure(dir = tempdir())

# list all cassettes
cassettes()
cassettes(on_disk = FALSE)

# list the currently active cassette
insert_cassette("stuffthings")
current_cassette()
eject_cassette()

cassettes()
cassettes(on_disk = FALSE)

# list the path to cassettes
cassette_path()
vcr_configure(dir = file.path(tempdir(), "foo"))
cassette_path()

vcr_configure_reset()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/persisters.R
\name{Persisters}
\alias{Persisters}
\alias{persister_fetch}
\title{Cassette persisters}
\usage{
persister_fetch(x = "FileSystem", file_name)
}
\description{
Keeps track of the cassette persisters in a hash-like object
}
\details{
There's only one option: \code{FileSystem}
\strong{Private Methods}
\describe{
\item{\code{persister_get()}}{
Gets and sets a named persister
}
}
}
\examples{
(aa <- Persisters$new())
aa$name
aa$persisters
yaml_serializer <- aa$persisters$new()
yaml_serializer
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{persisters}}{(list) internal use, holds persister object}

\item{\code{name}}{(character)}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Persisters$new()}}
\item \href{#method-clone}{\code{Persisters$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Persisters} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Persisters$new(persisters = list(), name = "FileSystem")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{persisters}}{(list) a list}

\item{\code{name}}{(character) Persister name, only option right now
is "FileSystem"}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Persisters} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Persisters$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/configuration.R
\name{vcr_configure}
\alias{vcr_configure}
\alias{vcr_configure_reset}
\alias{vcr_configuration}
\alias{vcr_config_defaults}
\title{Global Configuration Options}
\usage{
vcr_configure(...)

vcr_configure_reset()

vcr_configuration()

vcr_config_defaults()
}
\arguments{
\item{...}{configuration settings used to override defaults. See below for a
complete list of valid arguments.}
}
\description{
Configurable options that define vcr's default behavior.
}
\section{Configurable settings}{

\subsection{vcr options}{
\subsection{File locations}{
\itemize{
\item \code{dir} Cassette directory
\item \code{write_disk_path} (character) path to write files to
for any requests that write responses to disk. by default this parameter
is \code{NULL}. For testing a package, you'll probably want this path to
be in your \verb{tests/} directory, perhaps next to your cassettes
directory, e.g., where your cassettes are in \code{tests/fixtures}, your
files from requests that write to disk are in \code{tests/files}.
If you want to ignore these files in your installed package,
add them to \code{.Rinstignore}. If you want these files ignored on build
then add them to \code{.Rbuildignore} (though if you do, tests that depend
on these files probably will not work because they won't be found; so
you'll likely have to skip the associated tests as well).
}
}

\subsection{Contexts}{
\itemize{
\item \code{turned_off} (logical) VCR is turned on by default. Default:
\code{FALSE}
\item \code{allow_unused_http_interactions} (logical) Default: \code{TRUE}
\item \code{allow_http_connections_when_no_cassette} (logical) Determines how vcr
treats HTTP requests that are made when no vcr cassette is in use. When
\code{TRUE}, requests made when there is no vcr cassette in use will be allowed.
When \code{FALSE} (default), an \link{UnhandledHTTPRequestError} error will be raised
for any HTTP request made when there is no cassette in use
}
}

\subsection{Filtering}{
\itemize{
\item \code{ignore_hosts} (character) Vector of hosts to ignore. e.g., localhost, or
google.com. These hosts are ignored and real HTTP requests allowed to go
through
\item \code{ignore_localhost} (logical) Default: \code{FALSE}
\item \code{ignore_request} List of requests to ignore. NOT USED RIGHT NOW, sorry
\item \code{filter_sensitive_data} named list of values to replace. Format is:\preformatted{list(thing_to_replace_it_with = thing_to_replace)
}

We replace all instances of \code{thing_to_replace} with
\code{thing_to_replace_it_with}. Uses \code{\link[=gsub]{gsub()}} internally, with \code{fixed=TRUE};
so does exact matches. Before recording (writing to a cassette) we do
the replacement and then when reading from the cassette we do the reverse
replacement to get back to the real data. Before record replacement happens
in internal function \code{write_interactions()}, while before playback
replacement happens in internal function \code{YAML$deserialize()}
\item \code{filter_sensitive_data_regex} named list of values to replace. Follows
\code{filter_sensitive_data} format, except uses \code{fixed=FALSE} in the \code{\link[=gsub]{gsub()}}
function call; this means that the value in \code{thing_to_replace} is a regex
pattern.
\item \code{filter_request_headers} (character/list) \strong{request} headers to filter.
A character vector of request headers to remove - the headers will not be
recorded to disk. Alternatively, a named list similar to
\code{filter_sensitive_data} instructing vcr with what value to replace the
real value of the request header.
\item \code{filter_response_headers} (named list) \strong{response} headers to filter.
A character vector of response headers to remove - the headers will not be
recorded to disk. Alternatively, a named list similar to
\code{filter_sensitive_data} instructing vcr with what value to replace the
real value of the response header.
\item \code{filter_query_parameters} (named list) query parameters to filter.
A character vector of query parameters to remove - the query parameters
will not be recorded to disk. Alternatively, a named list similar to
\code{filter_sensitive_data} instructing vcr with what value to replace the
real value of the query parameter.
}
}

}

\subsection{Errors}{
\itemize{
\item \code{verbose_errors} Do you want more verbose errors or less verbose
errors when cassette recording/usage fails? Default is \code{FALSE}, that is,
less verbose errors. If \code{TRUE}, error messages will include more details
about what went wrong and suggest possible solutions. For testing
in an interactive R session, if \code{verbose_errors=FALSE}, you can run
\code{vcr_last_error()} to get the full error. If in non-interactive mode,
which most users will be in when running the entire test suite for a
package, you can set an environment variable (\code{VCR_VERBOSE_ERRORS})
to toggle this setting (e.g.,
\verb{Sys.setenv(VCR_VERBOSE_ERRORS=TRUE); devtools::test()})
}
\subsection{Internals}{
\itemize{
\item \code{cassettes} (list) don't use
\item \code{linked_context} (logical) linked context
\item \code{uri_parser} the uri parser, default: \code{crul::url_parse()}
}
}

\subsection{Logging}{
\itemize{
\item \code{log} (logical) should we log important vcr things? Default: \code{FALSE}
\item \code{log_opts} (list) Additional logging options:
\itemize{
\item 'file' either \code{"console"} or a file path to log to
\item 'log_prefix' default: "Cassette". We insert the cassette name after
that prefix, then the rest of the message.
\item More to come...
}
}
}

}

\subsection{Cassette Options}{

These settings can be configured globally, using \code{vcr_configure()}, or
locally, using either \code{use_cassette()} or \code{insert_cassette()}. Global
settings are applied to \emph{all} cassettes but are overridden by settings
defined locally for individual cassettes.
\itemize{
\item \code{record} (character) One of 'all', 'none', 'new_episodes', or 'once'.
See \link{recording}
\item \code{match_requests_on} vector of matchers. Default: (\code{method}, \code{uri})
See \link{request-matching} for details.
\item \code{serialize_with}: (character) "yaml" or "json". Note that you can have
multiple cassettes with the same name as long as they use different
serializers; so if you only want one cassette for a given cassette name,
make sure to not switch serializers, or clean up files you no longer need.
\item \code{json_pretty}: (logical) want JSON to be newline separated to be easier
to read? Or remove newlines to save disk space? default: FALSE
\item \code{persist_with} (character) only option is "FileSystem"
\item \code{preserve_exact_body_bytes} (logical) preserve exact body bytes for
\item \code{re_record_interval} (numeric) When given, the cassette will be
re-recorded at the given interval, in seconds.
\item \code{clean_outdated_http_interactions} (logical) Should outdated interactions
be recorded back to file. Default: \code{FALSE}
\item \code{quiet} (logical) Suppress any messages from both vcr and webmockr.
Default: \code{TRUE}
\item \code{warn_on_empty_cassette} (logical) Should a warning be thrown when an
empty cassette is detected? Empty cassettes are cleaned up (deleted) either
way. This option only determines whether a warning is thrown or not.
Default: \code{FALSE}
}
}
}

\examples{
vcr_configure(dir = tempdir())
vcr_configure(dir = tempdir(), record = "all")
vcr_configuration()
vcr_config_defaults()
vcr_configure(dir = tempdir(), ignore_hosts = "google.com")
vcr_configure(dir = tempdir(), ignore_localhost = TRUE)


# logging
vcr_configure(dir = tempdir(), log = TRUE,
  log_opts = list(file = file.path(tempdir(), "vcr.log")))
vcr_configure(dir = tempdir(), log = TRUE, log_opts = list(file = "console"))
vcr_configure(dir = tempdir(), log = TRUE,
 log_opts = list(
   file = file.path(tempdir(), "vcr.log"),
   log_prefix = "foobar"
))
vcr_configure(dir = tempdir(), log = FALSE)

# filter sensitive data
vcr_configure(dir = tempdir(),
  filter_sensitive_data = list(foo = "<bar>")
)
vcr_configure(dir = tempdir(),
  filter_sensitive_data = list(foo = "<bar>", hello = "<world>")
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_matching.R
\name{request-matching}
\alias{request-matching}
\title{vcr request matching}
\description{
There are a number of options, some of which are on by default, some of
which can be used together, and some alone.
}
\details{
To match previously recorded requests, \code{vcr} has to try to match new
HTTP requests to a previously recorded one. By default, we match on HTTP
method (e.g., \code{GET}) and URI (e.g., \verb{http://foo.com}), following Ruby’s
VCR gem.

You can customize how we match requests with one or more of the
following options, some of which are on by default, some of which can be
used together, and some alone.
\itemize{
\item \code{method}: Use the \strong{method} request matcher to match requests on
the HTTP method (i.e. GET, POST, PUT, DELETE, etc). You will
generally want to use this matcher. The \strong{method} matcher is used
(along with the \strong{uri} matcher) by default if you do not specify
how requests should match.
\item \code{uri}: Use the \strong{uri} request matcher to match requests on the
request URI. The \strong{uri} matcher is used (along with the \strong{method}
matcher) by default if you do not specify how requests should match.
\item \code{host}: Use the \strong{host} request matcher to match requests on the
request host. You can use this (alone, or in combination with
\strong{path}) as an alternative to \strong{uri} so that non-deterministic
portions of the URI are not considered as part of the request
matching.
\item \code{path}: Use the \strong{path} request matcher to match requests on the
path portion of the request URI. You can use this (alone, or in
combination with \strong{host}) as an alternative to \strong{uri} so that
non-deterministic portions of the URI
\item \code{query}: Use the \strong{query} request matcher to match requests on the
query string portion of the request URI. You can use this (alone, or
in combination with others) as an alternative to \strong{uri} so that
non-deterministic portions of the URI are not considered as part of
the request matching.
\item \code{body}: Use the \strong{body} request matcher to match requests on the
request body.
\item \code{headers}: Use the \strong{headers} request matcher to match requests on
the request headers.
}

You can set your own options by tweaking the \code{match_requests_on}
parameter in \code{use_cassette()}:\if{html}{\out{<div class="r">}}\preformatted{library(vcr)
}\if{html}{\out{</div>}}\if{html}{\out{<div class="r">}}\preformatted{use_cassette(name = "foo_bar", \{
    cli$post("post", body = list(a = 5))
  \}, 
  match_requests_on = c('method', 'headers', 'body')
)
}\if{html}{\out{</div>}}
\subsection{Matching}{
\subsection{headers}{\if{html}{\out{<div class="r">}}\preformatted{library(crul)
library(vcr)
cli <- crul::HttpClient$new("https://httpbin.org/get", 
  headers = list(foo = "bar"))
use_cassette(name = "nothing_new", \{
    one <- cli$get()
  \}, 
  match_requests_on = 'headers'
)
cli$headers$foo <- "stuff"
use_cassette(name = "nothing_new", \{
    two <- cli$get()
  \}, 
  match_requests_on = 'headers'
)
one$request_headers
two$request_headers
}\if{html}{\out{</div>}}
}

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request_handler.R
\name{RequestHandler}
\alias{RequestHandler}
\title{RequestHandler}
\description{
Base handler for http requests, deciding whether a
request is stubbed, to be ignored, recordable, or unhandled
}
\details{
\strong{Private Methods}
\describe{
\item{\code{request_type(request)}}{
Get the request type
}
\item{\code{externally_stubbed()}}{
just returns FALSE
}
\item{\code{should_ignore()}}{
should we ignore the request, depends on request ignorer
infrastructure that's not working yet
}
\item{\code{has_response_stub()}}{
Check if there is a matching response stub in the
http interaction list
}
\item{\code{get_stubbed_response()}}{
Check for a response and get it
}
\item{\code{request_summary(request)}}{
get a request summary
}
\item{\code{on_externally_stubbed_request(request)}}{
on externally stubbed request do nothing
}
\item{\code{on_ignored_request(request)}}{
on ignored request, do something
}
\item{\code{on_recordable_request(request)}}{
on recordable request, record the request
}
\item{\code{on_unhandled_request(request)}}{
on unhandled request, run UnhandledHTTPRequestError
}
}
}
\examples{
\dontrun{
# record mode: once
vcr_configure(
 dir = tempdir(),
 record = "once"
)

data(crul_request)
crul_request$url$handle <- curl::new_handle()
crul_request
x <- RequestHandler$new(crul_request)
# x$handle()

# record mode: none
vcr_configure(
 dir = tempdir(),
 record = "none"
)
data(crul_request)
crul_request$url$handle <- curl::new_handle()
crul_request
insert_cassette("testing_record_mode_none", record = "none")
#file.path(vcr_c$dir, "testing_record_mode_none.yml")
x <- RequestHandlerCrul$new(crul_request)
# x$handle()
crul_request$url$url <- "https://api.crossref.org/works/10.1039/c8sm90002g/"
crul_request$url$handle <- curl::new_handle()
z <- RequestHandlerCrul$new(crul_request)
# z$handle()
eject_cassette("testing_record_mode_none")
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{request_original}}{original, before any modification}

\item{\code{request}}{the request, after any modification}

\item{\code{vcr_response}}{holds \link{VcrResponse} object}

\item{\code{stubbed_response}}{the stubbed response}

\item{\code{cassette}}{the cassette holder}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{RequestHandler$new()}}
\item \href{#method-handle}{\code{RequestHandler$handle()}}
\item \href{#method-clone}{\code{RequestHandler$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{RequestHandler} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestHandler$new(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{The request from an object of class \code{HttpInteraction}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{RequestHandler} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-handle"></a>}}
\if{latex}{\out{\hypertarget{method-handle}{}}}
\subsection{Method \code{handle()}}{
Handle the request (\code{request} given in \verb{$initialize()})
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestHandler$handle()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
handles a request, outcomes vary
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RequestHandler$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/hooks.R
\name{Hooks}
\alias{Hooks}
\title{Hooks class}
\description{
Helps define new hooks, hold hooks, and accessors to get and
use hooks.
}
\details{
\strong{Private Methods}
\describe{
\item{\code{make_hook(x, plac, fun)}}{
Make a hook.
- x (character) Hook name
- plac Placement, one of "start" or "end"
- fun a function/callback
}
}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{hooks}}{intenal use}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-invoke_hook}{\code{Hooks$invoke_hook()}}
\item \href{#method-clear_hooks}{\code{Hooks$clear_hooks()}}
\item \href{#method-define_hook}{\code{Hooks$define_hook()}}
\item \href{#method-clone}{\code{Hooks$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-invoke_hook"></a>}}
\if{latex}{\out{\hypertarget{method-invoke_hook}{}}}
\subsection{Method \code{invoke_hook()}}{
invoke a hook
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Hooks$invoke_hook(hook_type, args)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{hook_type}}{(character) Hook name}

\item{\code{args}}{(named list) Args passed when invoking a hook}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
executes hook
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clear_hooks"></a>}}
\if{latex}{\out{\hypertarget{method-clear_hooks}{}}}
\subsection{Method \code{clear_hooks()}}{
clear all hooks
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Hooks$clear_hooks()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
no return
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-define_hook"></a>}}
\if{latex}{\out{\hypertarget{method-define_hook}{}}}
\subsection{Method \code{define_hook()}}{
define a hook
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Hooks$define_hook(hook_type, fun, prepend = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{hook_type}}{(character) Hook name}

\item{\code{fun}}{A function}

\item{\code{prepend}}{(logical) Whether to prepend or add to the end
of the string. Default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
no return; defines hook internally
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Hooks$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/insert_cassette.R
\name{insert_cassette}
\alias{insert_cassette}
\title{Insert a cassette to record HTTP requests}
\usage{
insert_cassette(
  name,
  record = NULL,
  match_requests_on = NULL,
  update_content_length_header = FALSE,
  allow_playback_repeats = FALSE,
  serialize_with = NULL,
  persist_with = NULL,
  preserve_exact_body_bytes = NULL,
  re_record_interval = NULL,
  clean_outdated_http_interactions = NULL
)
}
\arguments{
\item{name}{The name of the cassette. vcr will check this to ensure it
is a valid file name. Not allowed: spaces, file extensions, control
characters (e.g., \verb{\\n}), illegal characters ('/', '?', '<', '>', '\\', ':',
'*', '|', and '\"'), dots alone (e.g., '.', '..'), Windows reserved
words (e.g., 'com1'), trailing dots (can cause problems on Windows),
names longer than 255 characters. See section "Cassette names"}

\item{record}{The record mode (default: \code{"once"}). See \link{recording} for a
complete list of the different recording modes.}

\item{match_requests_on}{List of request matchers
to use to determine what recorded HTTP interaction to replay. Defaults to
\verb{["method", "uri"]}. The built-in matchers are "method", "uri", "host",
"path", "headers", "body" and "query"}

\item{update_content_length_header}{(logical) Whether or
not to overwrite the \code{Content-Length} header of the responses to
match the length of the response body. Default: \code{FALSE}}

\item{allow_playback_repeats}{(logical) Whether or not to
allow a single HTTP interaction to be played back multiple times.
Default: \code{FALSE}.}

\item{serialize_with}{(character) Which serializer to use.
Valid values are "yaml" (default) and "json". Note that you can have
multiple cassettes with the same name as long as they use different
serializers; so if you only want one cassette for a given cassette name,
make sure to not switch serializers, or clean up files you no longer need.}

\item{persist_with}{(character) Which cassette persister to
use. Default: "file_system". You can also register and use a
custom persister.}

\item{preserve_exact_body_bytes}{(logical) Whether or not
to base64 encode the bytes of the requests and responses for
this cassette when serializing it. See also \code{preserve_exact_body_bytes}
in \code{\link[=vcr_configure]{vcr_configure()}}. Default: \code{FALSE}}

\item{re_record_interval}{(integer) How frequently (in seconds) the
cassette should be re-recorded. default: \code{NULL} (not re-recorded)}

\item{clean_outdated_http_interactions}{(logical) Should outdated
interactions be recorded back to file? default: \code{FALSE}}
}
\value{
an object of class \code{Cassette}
}
\description{
Insert a cassette to record HTTP requests
}
\details{
Cassette names:
\itemize{
\item Should be meaningful so that it is obvious to you what test/function
they relate to. Meaningful names are important so that you can
quickly determine to what test file or test block a cassette
belongs. Note that vcr cannot check that your cassette names are
meaningful.
\item Should not be duplicated. Duplicated cassette names would lead to a
test using the wrong cassette.
\item Should not have spaces. Spaces can lead to problems in using file
paths.
\item Should not include a file extension. vcr handles file extensions for
the user.
\item Should not have illegal characters that can lead to problems in
using file paths: \code{/}, \verb{?}, \code{<}, \code{>}, \verb{\\}, \code{:}, \code{*}, \code{|}, and \verb{\"}
\item Should not have control characters, e.g., \verb{\n}
\item Should not have just dots, e.g., \code{.} or \code{..}
\item Should not have Windows reserved words, e.g., \code{com1}
\item Should not have trailing dots
\item Should not be longer than 255 characters
}

\code{vcr::check_cassette_names()} is meant to be run during your tests, from
a \code{setup-pkgname.R} file inside the \code{tests/testthat} directory. It only
checks that cassette names are not duplicated. Note that if you do need
to have duplicated cassette names you can do so by using the
\code{allowed_duplicates} parameter in \code{check_cassette_names()}. A helper
function \code{check_cassette_names()} runs inside
\code{\link[=insert_cassette]{insert_cassette()}} that checks that cassettes
do not have: spaces, file extensions, unaccepted characters (slashes)
}
\section{Cassette options}{


Default values for arguments controlling cassette behavior are
inherited from vcr's global configuration. See \code{\link[=vcr_configure]{vcr_configure()}} for a
complete list of options and their default settings. You can override these
options for a specific cassette by changing an argument's value to something
other than \code{NULL} when calling either \code{insert_cassette()} or
\code{use_cassette()}.
}

\examples{
\dontrun{
library(vcr)
library(crul)
vcr_configure(dir = tempdir())
webmockr::webmockr_allow_net_connect()

(x <- insert_cassette(name = "leo5"))
current_cassette()
x$new_recorded_interactions
x$previously_recorded_interactions()
cli <- crul::HttpClient$new(url = "https://httpbin.org")
cli$get("get")
x$new_recorded_interactions # 1 interaction
x$previously_recorded_interactions() # empty
webmockr::stub_registry() # not empty
# very important when using inject_cassette: eject when done
x$eject() # same as eject_cassette("leo5")
x$new_recorded_interactions # same, 1 interaction
x$previously_recorded_interactions() # now not empty
## stub_registry now empty, eject() calls webmockr::disable(), which
## calls the disable method for each of crul and httr adadapters,
## which calls webmockr's remove_stubs() method for each adapter
webmockr::stub_registry()

# cleanup
unlink(file.path(tempdir(), "leo5.yml"))
}
}
\seealso{
\code{\link[=use_cassette]{use_cassette()}}, \code{\link[=eject_cassette]{eject_cassette()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_cassette.R
\name{use_cassette}
\alias{use_cassette}
\title{Use a cassette to record HTTP requests}
\usage{
use_cassette(
  name,
  ...,
  record = NULL,
  match_requests_on = NULL,
  update_content_length_header = FALSE,
  allow_playback_repeats = FALSE,
  serialize_with = NULL,
  persist_with = NULL,
  preserve_exact_body_bytes = NULL,
  re_record_interval = NULL,
  clean_outdated_http_interactions = NULL
)
}
\arguments{
\item{name}{The name of the cassette. vcr will check this to ensure it
is a valid file name. Not allowed: spaces, file extensions, control
characters (e.g., \verb{\\n}), illegal characters ('/', '?', '<', '>', '\\', ':',
'*', '|', and '\"'), dots alone (e.g., '.', '..'), Windows reserved
words (e.g., 'com1'), trailing dots (can cause problems on Windows),
names longer than 255 characters. See section "Cassette names"}

\item{...}{a block of code containing one or more requests (required). Use
curly braces to encapsulate multi-line code blocks. If you can't pass a code
block use \code{\link[=insert_cassette]{insert_cassette()}} instead.}

\item{record}{The record mode (default: \code{"once"}). See \link{recording} for a
complete list of the different recording modes.}

\item{match_requests_on}{List of request matchers
to use to determine what recorded HTTP interaction to replay. Defaults to
\verb{["method", "uri"]}. The built-in matchers are "method", "uri", "host",
"path", "headers", "body" and "query"}

\item{update_content_length_header}{(logical) Whether or
not to overwrite the \code{Content-Length} header of the responses to
match the length of the response body. Default: \code{FALSE}}

\item{allow_playback_repeats}{(logical) Whether or not to
allow a single HTTP interaction to be played back multiple times.
Default: \code{FALSE}.}

\item{serialize_with}{(character) Which serializer to use.
Valid values are "yaml" (default) and "json". Note that you can have
multiple cassettes with the same name as long as they use different
serializers; so if you only want one cassette for a given cassette name,
make sure to not switch serializers, or clean up files you no longer need.}

\item{persist_with}{(character) Which cassette persister to
use. Default: "file_system". You can also register and use a
custom persister.}

\item{preserve_exact_body_bytes}{(logical) Whether or not
to base64 encode the bytes of the requests and responses for
this cassette when serializing it. See also \code{preserve_exact_body_bytes}
in \code{\link[=vcr_configure]{vcr_configure()}}. Default: \code{FALSE}}

\item{re_record_interval}{(integer) How frequently (in seconds) the
cassette should be re-recorded. default: \code{NULL} (not re-recorded)}

\item{clean_outdated_http_interactions}{(logical) Should outdated
interactions be recorded back to file? default: \code{FALSE}}
}
\value{
an object of class \code{Cassette}
}
\description{
Use a cassette to record HTTP requests
}
\details{
A run down of the family of top level \pkg{vcr} functions
\itemize{
\item \code{use_cassette} Initializes a cassette. Returns the inserted
cassette.
\item \code{insert_cassette} Internally used within \code{use_cassette}
\item \code{eject_cassette} ejects the current cassette. The cassette
will no longer be used. In addition, any newly recorded HTTP interactions
will be written to disk.
}
}
\section{Cassette options}{


Default values for arguments controlling cassette behavior are
inherited from vcr's global configuration. See \code{\link[=vcr_configure]{vcr_configure()}} for a
complete list of options and their default settings. You can override these
options for a specific cassette by changing an argument's value to something
other than \code{NULL} when calling either \code{insert_cassette()} or
\code{use_cassette()}.
}

\section{Behavior}{

This function handles a few different scenarios:
\itemize{
\item when everything runs smoothly, and we return a \code{Cassette} class object
so you can inspect the cassette, and the cassette is ejected
\item when there is an invalid parameter input on cassette creation,
we fail with a useful message, we don't return a cassette, and the
cassette is ejected
\item when there is an error in calling your passed in code block,
we return with a useful message, and since we use \code{on.exit()}
the cassette is still ejected even though there was an error,
but you don't get an object back
\item whenever an empty cassette (a yml/json file) is found, we delete it
before returning from the \code{use_cassette()} function call. we achieve
this via use of \code{on.exit()} so an empty cassette is deleted even
if there was an error in the code block you passed in
}
}

\section{Cassettes on disk}{

Note that \emph{"eject"} only means that the R session cassette is no longer
in use. If any interactions were recorded to disk, then there is a file
on disk with those interactions.
}

\section{Using with tests (specifically \pkg{testthat})}{

There's a few ways to get correct line numbers for failed tests and
one way to not get correct line numbers:

\emph{Correct}: Either wrap your \code{test_that()} block inside your \code{use_cassette()}
block, OR if you put your \code{use_cassette()} block inside your \code{test_that()}
block put your \code{testthat} expectations outside of the \code{use_cassette()}
block.

\emph{Incorrect}: By wrapping the \code{use_cassette()} block inside your
\code{test_that()} block with your \pkg{testthat} expectations inside the
\code{use_cassette()} block, you'll only get the line number that the
\code{use_cassette()} block starts on.
}

\examples{
\dontrun{
library(vcr)
library(crul)
vcr_configure(dir = tempdir())

use_cassette(name = "apple7", {
  cli <- HttpClient$new(url = "https://httpbin.org")
  resp <- cli$get("get")
})
readLines(file.path(tempdir(), "apple7.yml"))

# preserve exact body bytes - records in base64 encoding
use_cassette("things4", {
  cli <- crul::HttpClient$new(url = "https://httpbin.org")
  bbb <- cli$get("get")
}, preserve_exact_body_bytes = TRUE)
## see the body string value in the output here
readLines(file.path(tempdir(), "things4.yml"))

# cleanup
unlink(file.path(tempdir(), c("things4.yml", "apple7.yml")))


# with httr
library(vcr)
library(httr)
vcr_configure(dir = tempdir(), log = TRUE, log_opts = list(file = file.path(tempdir(), "vcr.log")))

use_cassette(name = "stuff350", {
  res <- GET("https://httpbin.org/get")
})
readLines(file.path(tempdir(), "stuff350.yml"))

use_cassette(name = "catfact456", {
  res <- GET("https://catfact.ninja/fact")
})

# record mode: none
library(crul)
vcr_configure(dir = tempdir())

## make a connection first
conn <- crul::HttpClient$new("https://eu.httpbin.org")
## this errors because 'none' disallows any new requests
# use_cassette("none_eg", (res2 <- conn$get("get")), record = "none")
## first use record mode 'once' to record to a cassette
one <- use_cassette("none_eg", (res <- conn$get("get")), record = "once")
one; res
## then use record mode 'none' to see it's behavior
two <- use_cassette("none_eg", (res2 <- conn$get("get")), record = "none")
two; res2
}
}
\seealso{
\code{\link[=insert_cassette]{insert_cassette()}}, \code{\link[=eject_cassette]{eject_cassette()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcr-package.R
\docType{data}
\name{crul_request}
\alias{crul_request}
\title{An HTTP request as prepared by the \pkg{crul} package}
\format{
A list
}
\description{
The object is a list, and is the object that is passed on to
\pkg{webmockr} and \pkg{vcr} instead of routing through
\pkg{crul} as normal. Used in examples/tests.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http_interaction.R
\name{HTTPInteraction}
\alias{HTTPInteraction}
\title{HTTPInteraction class}
\description{
object holds request and response objects
}
\details{
\strong{Methods}
\describe{
\item{\code{to_hash()}}{
Create a hash from the HTTPInteraction object
}
\item{\code{from_hash(hash)}}{
Create a HTTPInteraction object from a hash
}
}
}
\examples{
\dontrun{
# make the request
library(vcr)
url <- "https://eu.httpbin.org/post"
body <- list(foo = "bar")
cli <- crul::HttpClient$new(url = url)
res <- cli$post(body = body)

# build a Request object
(request <- Request$new("POST", uri = url,
  body = body, headers = res$response_headers))
# build a VcrResponse object
(response <- VcrResponse$new(
   res$status_http(),
   res$response_headers,
   res$parse("UTF-8"),
   res$response_headers$status))

# make HTTPInteraction object
(x <- HTTPInteraction$new(request = request, response = response))
x$recorded_at
x$to_hash()

# make an HTTPInteraction from a hash with the object already made
x$from_hash(x$to_hash())

# Make an HTTPInteraction from a hash alone
my_hash <- x$to_hash()
HTTPInteraction$new()$from_hash(my_hash)
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{request}}{A \code{Request} class object}

\item{\code{response}}{A \code{VcrResponse} class object}

\item{\code{recorded_at}}{(character) Time http interaction recorded at}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{HTTPInteraction$new()}}
\item \href{#method-to_hash}{\code{HTTPInteraction$to_hash()}}
\item \href{#method-from_hash}{\code{HTTPInteraction$from_hash()}}
\item \href{#method-clone}{\code{HTTPInteraction$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{HTTPInteraction} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteraction$new(request, response, recorded_at)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{A \code{Request} class object}

\item{\code{response}}{A \code{VcrResponse} class object}

\item{\code{recorded_at}}{(character) Time http interaction recorded at}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{HTTPInteraction} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_hash"></a>}}
\if{latex}{\out{\hypertarget{method-to_hash}{}}}
\subsection{Method \code{to_hash()}}{
Create a hash from the HTTPInteraction object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteraction$to_hash()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a named list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-from_hash"></a>}}
\if{latex}{\out{\hypertarget{method-from_hash}{}}}
\subsection{Method \code{from_hash()}}{
Create a HTTPInteraction object from a hash
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteraction$from_hash(hash)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{hash}}{a named list}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new \code{HttpInteraction} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HTTPInteraction$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/errors.R
\name{UnhandledHTTPRequestError}
\alias{UnhandledHTTPRequestError}
\alias{vcr_last_error}
\title{UnhandledHTTPRequestError}
\usage{
vcr_last_error()
}
\description{
Handle http request errors
}
\details{
How this error class is used:
If \code{record="once"} we trigger this.

Users can use vcr in the context of both \code{\link[=use_cassette]{use_cassette()}}
and \code{\link[=insert_cassette]{insert_cassette()}}

For the former, all requests go through the call_block
But for the latter, requests go through webmockr.

Where is one place where we can put UnhandledHTTPRequestError
that will handle both use_cassette and insert_cassette?
}
\section{Error situations where this is invoked}{

\itemize{
\item \code{record=once} AND there's a new request that doesn't match
the one in the cassette on disk
\itemize{
\item in webmockr: if no stub found and there are recorded
interactions on the cassette, and record = once, then
error with UnhandledHTTPRequestError
\itemize{
\item but if record != once, then allow it, unless record == none
}
}
\item others?
}
}

\examples{
vcr_configure(dir = tempdir())
cassettes()
insert_cassette("turtle")
request <- Request$new("post", 'https://eu.httpbin.org/post?a=5',
  "", list(foo = "bar"))

err <- UnhandledHTTPRequestError$new(request)
err$request_description()
err$current_matchers()
err$match_request_on_headers()
err$match_request_on_body()
err$formatted_headers()
cat(err$formatted_headers(), "\n")
cat(err$cassettes_description(), "\n")
cat(err$cassettes_list(), "\n")
err$formatted_suggestions()
cat(err$format_bullet_point('foo bar', 1), "\n")
err$suggestion_for("use_new_episodes")
err$suggestions()
err$no_cassette_suggestions()
err$record_mode_suggestion()
err$has_used_interaction_matching()
err$match_requests_on_suggestion()

# err$construct_message()

# cleanup
eject_cassette("turtle")
unlink(tempdir())
\dontrun{
# vcr_last_error()
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{request}}{a \link{Request} object}

\item{\code{cassette}}{a cassette name}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{UnhandledHTTPRequestError$new()}}
\item \href{#method-run}{\code{UnhandledHTTPRequestError$run()}}
\item \href{#method-construct_message}{\code{UnhandledHTTPRequestError$construct_message()}}
\item \href{#method-request_description}{\code{UnhandledHTTPRequestError$request_description()}}
\item \href{#method-current_matchers}{\code{UnhandledHTTPRequestError$current_matchers()}}
\item \href{#method-match_request_on_headers}{\code{UnhandledHTTPRequestError$match_request_on_headers()}}
\item \href{#method-match_request_on_body}{\code{UnhandledHTTPRequestError$match_request_on_body()}}
\item \href{#method-formatted_headers}{\code{UnhandledHTTPRequestError$formatted_headers()}}
\item \href{#method-cassettes_description}{\code{UnhandledHTTPRequestError$cassettes_description()}}
\item \href{#method-cassettes_list}{\code{UnhandledHTTPRequestError$cassettes_list()}}
\item \href{#method-get_help}{\code{UnhandledHTTPRequestError$get_help()}}
\item \href{#method-formatted_suggestions}{\code{UnhandledHTTPRequestError$formatted_suggestions()}}
\item \href{#method-format_bullet_point}{\code{UnhandledHTTPRequestError$format_bullet_point()}}
\item \href{#method-format_foot_note}{\code{UnhandledHTTPRequestError$format_foot_note()}}
\item \href{#method-suggestion_for}{\code{UnhandledHTTPRequestError$suggestion_for()}}
\item \href{#method-suggestions}{\code{UnhandledHTTPRequestError$suggestions()}}
\item \href{#method-no_cassette_suggestions}{\code{UnhandledHTTPRequestError$no_cassette_suggestions()}}
\item \href{#method-record_mode_suggestion}{\code{UnhandledHTTPRequestError$record_mode_suggestion()}}
\item \href{#method-has_used_interaction_matching}{\code{UnhandledHTTPRequestError$has_used_interaction_matching()}}
\item \href{#method-match_requests_on_suggestion}{\code{UnhandledHTTPRequestError$match_requests_on_suggestion()}}
\item \href{#method-clone}{\code{UnhandledHTTPRequestError$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{UnhandledHTTPRequestError} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$new(request)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{request}}{(Request) a \link{Request} object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{UnhandledHTTPRequestError} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-run"></a>}}
\if{latex}{\out{\hypertarget{method-run}{}}}
\subsection{Method \code{run()}}{
Run unhandled request handling
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$run()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
various
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-construct_message"></a>}}
\if{latex}{\out{\hypertarget{method-construct_message}{}}}
\subsection{Method \code{construct_message()}}{
Construct and execute stop message for why request failed
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$construct_message()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a stop message
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-request_description"></a>}}
\if{latex}{\out{\hypertarget{method-request_description}{}}}
\subsection{Method \code{request_description()}}{
construct request description
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$request_description()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-current_matchers"></a>}}
\if{latex}{\out{\hypertarget{method-current_matchers}{}}}
\subsection{Method \code{current_matchers()}}{
get current request matchers
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$current_matchers()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-match_request_on_headers"></a>}}
\if{latex}{\out{\hypertarget{method-match_request_on_headers}{}}}
\subsection{Method \code{match_request_on_headers()}}{
are headers included in current matchers?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$match_request_on_headers()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-match_request_on_body"></a>}}
\if{latex}{\out{\hypertarget{method-match_request_on_body}{}}}
\subsection{Method \code{match_request_on_body()}}{
is body includled in current matchers?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$match_request_on_body()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-formatted_headers"></a>}}
\if{latex}{\out{\hypertarget{method-formatted_headers}{}}}
\subsection{Method \code{formatted_headers()}}{
get request headers
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$formatted_headers()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cassettes_description"></a>}}
\if{latex}{\out{\hypertarget{method-cassettes_description}{}}}
\subsection{Method \code{cassettes_description()}}{
construct description of current or lack thereof cassettes
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$cassettes_description()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cassettes_list"></a>}}
\if{latex}{\out{\hypertarget{method-cassettes_list}{}}}
\subsection{Method \code{cassettes_list()}}{
cassette details
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$cassettes_list()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_help"></a>}}
\if{latex}{\out{\hypertarget{method-get_help}{}}}
\subsection{Method \code{get_help()}}{
get help message for non-verbose error
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$get_help()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-formatted_suggestions"></a>}}
\if{latex}{\out{\hypertarget{method-formatted_suggestions}{}}}
\subsection{Method \code{formatted_suggestions()}}{
make suggestions for what to do
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$formatted_suggestions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-format_bullet_point"></a>}}
\if{latex}{\out{\hypertarget{method-format_bullet_point}{}}}
\subsection{Method \code{format_bullet_point()}}{
add bullet point to beginning of a line
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$format_bullet_point(lines, index)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lines}}{(character) vector of strings}

\item{\code{index}}{(integer) a number}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-format_foot_note"></a>}}
\if{latex}{\out{\hypertarget{method-format_foot_note}{}}}
\subsection{Method \code{format_foot_note()}}{
make a foot note
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$format_foot_note(url, index)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{url}}{(character) a url}

\item{\code{index}}{(integer) a number}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-suggestion_for"></a>}}
\if{latex}{\out{\hypertarget{method-suggestion_for}{}}}
\subsection{Method \code{suggestion_for()}}{
get a suggestion by key
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$suggestion_for(key)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{key}}{(character) a character string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-suggestions"></a>}}
\if{latex}{\out{\hypertarget{method-suggestions}{}}}
\subsection{Method \code{suggestions()}}{
get all suggestions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$suggestions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-no_cassette_suggestions"></a>}}
\if{latex}{\out{\hypertarget{method-no_cassette_suggestions}{}}}
\subsection{Method \code{no_cassette_suggestions()}}{
get all no cassette suggestions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$no_cassette_suggestions()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-record_mode_suggestion"></a>}}
\if{latex}{\out{\hypertarget{method-record_mode_suggestion}{}}}
\subsection{Method \code{record_mode_suggestion()}}{
get the appropriate record mode suggestion
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$record_mode_suggestion()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-has_used_interaction_matching"></a>}}
\if{latex}{\out{\hypertarget{method-has_used_interaction_matching}{}}}
\subsection{Method \code{has_used_interaction_matching()}}{
are there any used interactions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$has_used_interaction_matching()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
logical
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-match_requests_on_suggestion"></a>}}
\if{latex}{\out{\hypertarget{method-match_requests_on_suggestion}{}}}
\subsection{Method \code{match_requests_on_suggestion()}}{
match requests on suggestion
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$match_requests_on_suggestion()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UnhandledHTTPRequestError$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/check_cassette_names.R
\name{check_cassette_names}
\alias{check_cassette_names}
\title{Check cassette names}
\usage{
check_cassette_names(
  pattern = "test-",
  behavior = "stop",
  allowed_duplicates = NULL
)
}
\arguments{
\item{pattern}{(character) regex pattern for file paths to check.
this is done inside of \verb{tests/testthat/}. default: "test-"}

\item{behavior}{(character) "stop" (default) or "warning". if "warning",
we use \code{immediate.=TRUE} so the warning happens at the top of your
tests rather than you seeing it after tests have run (as would happen
by default)}

\item{allowed_duplicates}{(character) cassette names that can be duplicated}
}
\description{
Check cassette names
}
\details{
Cassette names:
\itemize{
\item Should be meaningful so that it is obvious to you what test/function
they relate to. Meaningful names are important so that you can
quickly determine to what test file or test block a cassette
belongs. Note that vcr cannot check that your cassette names are
meaningful.
\item Should not be duplicated. Duplicated cassette names would lead to a
test using the wrong cassette.
\item Should not have spaces. Spaces can lead to problems in using file
paths.
\item Should not include a file extension. vcr handles file extensions for
the user.
\item Should not have illegal characters that can lead to problems in
using file paths: \code{/}, \verb{?}, \code{<}, \code{>}, \verb{\\}, \code{:}, \code{*}, \code{|}, and \verb{\"}
\item Should not have control characters, e.g., \verb{\n}
\item Should not have just dots, e.g., \code{.} or \code{..}
\item Should not have Windows reserved words, e.g., \code{com1}
\item Should not have trailing dots
\item Should not be longer than 255 characters
}

\code{vcr::check_cassette_names()} is meant to be run during your tests, from
a \code{setup-pkgname.R} file inside the \code{tests/testthat} directory. It only
checks that cassette names are not duplicated. Note that if you do need
to have duplicated cassette names you can do so by using the
\code{allowed_duplicates} parameter in \code{check_cassette_names()}. A helper
function \code{check_cassette_names()} runs inside
\code{\link[=insert_cassette]{insert_cassette()}} that checks that cassettes
do not have: spaces, file extensions, unaccepted characters (slashes)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serializers-yaml.R
\name{YAML}
\alias{YAML}
\title{The YAML serializer}
\description{
class with methods for serializing via the \pkg{yaml} package
}
\keyword{internal}
\section{Super class}{
\code{\link[vcr:Serializer]{vcr::Serializer}} -> \code{YAML}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{YAML$new()}}
\item \href{#method-serialize}{\code{YAML$serialize()}}
\item \href{#method-deserialize}{\code{YAML$deserialize()}}
\item \href{#method-clone}{\code{YAML$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new YAML object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{YAML$new(path = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{(character) path to the cassette, excluding the cassette
directory and the file extension}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{YAML} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize"></a>}}
\if{latex}{\out{\hypertarget{method-serialize}{}}}
\subsection{Method \code{serialize()}}{
Serializes the given hash using internal fxn write_yaml
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{YAML$serialize(x, path, bytes)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{(list) the object to serialize}

\item{\code{path}}{(character) the file path}

\item{\code{bytes}}{(logical) whether to preserve exact body bytes or not}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
(character) the YAML string to write to disk
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize}{}}}
\subsection{Method \code{deserialize()}}{
Deserializes the content at the path using
yaml::yaml.load_file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{YAML$deserialize()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
(list) the deserialized object, an R list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{YAML$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
