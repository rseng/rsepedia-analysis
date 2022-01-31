rcrossref: R interface to CrossRef APIs
=======================================



[![cran checks](https://cranchecks.info/badges/worst/rcrossref)](https://cranchecks.info/pkgs/rcrossref)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/rcrossref/workflows/R-check/badge.svg)](https://github.com/ropensci/rcrossref/actions/)
[![codecov](https://codecov.io/gh/ropensci/rcrossref/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rcrossref)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcrossref)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcrossref)](https://cran.r-project.org/package=rcrossref)

## CrossRef documentation

* Crossref API: https://github.com/CrossRef/rest-api-doc#readme
* Crossref's API issue tracker: https://gitlab.com/crossref/issues
* Crossref metadata search API: https://www.crossref.org/labs/crossref-metadata-search/
* Crossref DOI Content Negotiation: https://citation.crosscite.org/docs.html
* Crossref Text and Data Mining (TDM) Services: https://www.crossref.org/education/retrieve-metadata/rest-api/text-and-data-mining/

## Installation

Stable version from CRAN


```r
install.packages("rcrossref")
```

Or development version from GitHub


```r
remotes::install_github("ropensci/rcrossref")
```

Load `rcrossref`


```r
library('rcrossref')
```

## Register for the Polite Pool

If you are intending to access Crossref regularly you will want to send your email address with your queries. This has the advantage that queries are placed in the polite pool of servers. Including your email address is good practice as described in the Crossref documentation under Good manners (https://github.com/CrossRef/rest-api-doc#good-manners--more-reliable-service). The second advantage is that Crossref can contact you if there is a problem with a query.

Details on how to register your email in a call can be found at `?rcrossref-package`. To pass your email address to Crossref, simply store it as an environment variable in .Renviron like this:

Open file: `file.edit("~/.Renviron")`

Add email address to be shared with Crossref `crossref_email= "name@example.com"`

Save the file and restart your R session

To stop sharing your email when using rcrossref simply delete it from your .Renviron file. 

## Documentation

See https://docs.ropensci.org/rcrossref/ to get started

## Meta

* Please report any issues or bugs: https://github.com/ropensci/rcrossref/issues
* License: MIT
* Get citation information for `rcrossref` in R doing `citation(package = 'rcrossref')`
* Please note that this package is released with a Contributor Code of Conduct (https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

---

This package is part of a richer suite called fulltext (https://github.com/ropensci/fulltext), along with several other packages, that provides the ability to search for and retrieve full text of open access scholarly articles.

---
rcrossref 1.1.0
===============

### MINOR IMPROVEMENTS

* move `bibtex` package to Suggests as it has been orphaned, use it conditionally  (#209)
* change all uses of `tibble::tbl_df` to `tibble::as_tibble` (#206)
* additional fields now parsed in the `cr_works()` output: short-container-title, references-count, is-referenced-by-count, language, content-domain, and update-to (#208)


rcrossref 1.0.0
===============

### MINOR IMPROVEMENTS

* docs page vignette has names (#197)
* `query.title` field query is no longer supported by Crossref, removed from package (#198)
* fix some links in readme (#199)
* improve documentation for how to do deep paging (new section "Deep paging"), and how to use the cursor specifically (#190)

### BUG FIXES

* fix error when progress used and when limit param not used  (#190)
* encoding fix for the Crossref RStudio Addin (#194) (#201)


rcrossref 0.9.2
===============

## NEW FEATURES

* `cr_funders()`, `cr_journals()`, `cr_licenses()`, `cr_members()`, `cr_prefixes()`, `cr_types()`, `cr_works()` gain ability to show a progress bar when using deep pagination when using `works=TRUE` (#186) (#188)

### MINOR IMPROVEMENTS

* fix a test fixture that had non-ascii characters (#187)
* `cr_works()` returned a tibble in the `$data` slot in each case except for when a single DOI was passed to the `doi` parameter. now fixed  (#184) thanks @martinjhnhadley


rcrossref 0.9.0
===============

### NEW FEATURES

* big Crossref Addin update (#171) all work by @haozhu233
* Async HTTP introduced for `cr_works()`, `cr_works_()`, and `cr_citation_count()`. See the new parameter `async` (logical) in those functions. For `cr_citation_count()`, it now accepts more than 1 DOI, and the output has changed from a numeric value to a data.frame (columns: `doi` and `count`). With `async=TRUE` for `cr_works()` you get a list of data.frame's; while for `cr_works_()` you get a list of JSON's (#121) (#160) (#182)

### MINOR IMPROVEMENTS

* package tests now using `vcr` for HTTP request/response caching (#178) (#179)
* in `works` data, now returning `published-print` and `published-online` fields (#181)
* add new filters for `/works` routes (when `works = TRUE`): `isbn`, `reference_visibility`, `has_content_domain`, and `has_domain_restriction` (#176) (#177)
* add new return element in `/works` routes (when `works = TRUE`): `$reference` gives the references cited in the article (these are not the articles citing the target article, sorry) (#176)
* improve documentation on the "polite pool". If you supply an email address Crossref will put you in the polite pool (#173) (#175) thanks @poldham
* added example to `cr_abstract()` of handling many DOIs while allowing for failures without stopping progress (#174) thanks @zackbatist

### BUG FIXES

* fix internal parsing of funder data in `cr_works()` (#180) thanks @nicholasmfraser for the bug report
* fix for `cr_citation_count()`: when given a bad/invalid/malformed DOI, throw warning and give back an `NA` (#164) thanks for the report @chreman
* Fix for instances in the package where we incorrectly were doing logical comparison's; detected via _R_CHECK_LENGTH_1_LOGIC2_


rcrossref 0.8.4
===============

### NEW FEATURES

* The RStudio Addin now includes ability to search by article metadata, e.g., title (#148) (#152) 
* you can now set a different base url for `cr_cn()`  (#168)

### MINOR IMPROVEMENTS

* better behavior for `cr_cn()` when bibentry not valid, that is not parseable. before the package we use to parse `bibtex` would stop on invalid bibentry data, but now we get around the invalid bits and give back the bibentry (#147)
* updated and improved documentation for filters (#161)
* improved description of the `dois` parameter (#162) thanks @ms609
* Much improved error behavior for `cr_*` functions. We now give errors like `404 (client error): /works/blblbl - Resource not found.`, which includes HTTP status code, major class of error (success, redirect, client, server), the route requested, and the error message (#163) 

### BUG FIXES

* fixed bug in `cr_works()` in which field queries should have been possible for title and affiliation, but were not. Fixed now.  (#149)
* `cr_journals()` was not correctly parsing data when more than 1 ISSN given and `works` set to `TRUE`, fixed now (#156)
* `cr_journals()` was not correctly parsing data when no ISSN found and `works` set to `TRUE`, fixed now (#150)
* `cr_journals()` was not correctly handling queries with multiple ISSN's and `works` set to `FALSE`, fixed now (#151)
* fixed bug in `cr_works()`, and any other `cr_*` function that set `works=TRUE`. the `license` slot can have more than 1 result, and we were only giving the first back. fixed the parsing on this to give back all license results (#170)


rcrossref 0.8.0
===============

### NEW FEATURES

* Include documentation and internals to support passing a user 
supplied email to Crossref to get put into higher rate limit
pool (#145) thanks @njahn82
* Functions that operate on the `/works` route gain a `select`
parameter to select certain fields to return (#146)

### MINOR IMPROVEMENTS

* Updated docs for additional options for the `sort` parameter (#142)
* Updated docs for additional options for doing field queries (the 
`flq` parameter) (#143)
* Filters: New filters added to the package, and some removed 
(e.g, `publisher-name`). You can see the filters with the 
functions `filter_details`/`filter_names`. Beware, some filters
error sometimes with the Crossref API - they may not work, but they may, let me know at <https://github.com/ropensci/rcrossref/issues> or let Crossref know at <https://github.com/CrossRef/rest-api-doc/issues> (#136) (#139) (#141)


rcrossref 0.7.0
===============

### NEW FEATURES

* All text mining functionality moved into a new package: `crminer`
<https://github.com/ropensci/crminer> . Functions that did
text mining stuff now defunct, see `?rcrossref-defunct` (#122)
* All Crossref API requests now using `https` instead of `http` (#133)

### MINOR IMPROVEMENTS

* replace `xml2::xml_find_one` with `xml2::xml_find_first` (#128)
thanks @njahn82
* replaced `httr` with `crul` for HTTP requests (#132)
* Now using markdown for documentation (#135)
* Added more documentation to `cr_journals` and `cr_works`
about what the returned data fields `backfile_dois` and `current_dois`
really mean (#105) thanks @SteveViss

### BUG FIXES

* Fix to `cr_prefixes` to not fail when no results found (#130)
thanks @globbestael
* Fixed `cr_works` to allow queries like `facet = license:*` to be
passed to `facet` parameter (was always allowed by Crossref, but we
neglected to allow it - previously only allowed a boolean) (#129)
* Fixed `cr_funders` and `cr_journals` to give back facet data along
with other data (#134)
* Fix to `cr_*` functions to check for a missing content-type
headers and instead of failing, we continue anyway and try to parse
data as sometimes Crossref doesn't give back a content type header
at all (#127)


rcrossref 0.6.0
===============

### NEW FEATURES

* Added support for field queries (see
<https://github.com/CrossRef/rest-api-doc/blob/master/rest_api.md#field-queries>
for information on them). Means that you can query on specific fields rather
than using `query` parameter which queries across all fields (#111)
* Now using `rappdirs` for local storage and caching for `cr_ft_text`  (#106)

### MINOR IMPROVEMENTS

* Added to man files where appropriate new 10K max value for the
`offset` parameter (#126)
* Added to pkg level man file new rate limit headers included,
and how users can get to those, via `config=verbose()` call (#124)
* Better failure modes on input parameters, still more work to do
surely (#101)
* sleeping now between tests to avoid making crossref rate
limit gate keepers mad (#125)
* `cr_search` and `cr_search_free` are now defunct. They were marked
deprecated in previous version, and warned of defunct, and now
they are defunct. Similar functionality can be done with e.g., `cr_works()`
(#102)
* `crosscite` is now defunct. The functionality of this function can be
achieved with `cr_cn()` (#82)
* `cr_fundref` is now defunct. Crossref changed their name `fundref`
to `funders`, so we've changed our function, see `cr_funders()` (#83)
* parameter `sample` maximum value is now 100, was previously 1000.
documentation updated. (#118)
* New filters `has-clinical-trial-number` and `has-abstract` added to
the package, see `?filters` for help (#120)



rcrossref 0.5.8
===============

### NEW FEATURES

* Added an RStudio Addin for searching for citations. See `?rcrossref` for
more. Addin authored by Hao Zhu @haozhu233 (#114)
* New function `cr_abstract()` that tries to get an abstract via XML provided by
Crossref - NOTE: an abstract is rarely available though (#116)

### BUG FIXES

* Fixed bug in `cr_cn()` where DOIs with no minting agency found were
failing because we were previously stopping when no agency found.
Now, we just assume Crossref and move on from there. (#117)
thanks @dfalster !
* Fix to `cr_r()` when number requested > 100. Actual fix is in
`cr_works()`. Max for sample used to be 1000, asked this on the
Crossref API forum,
see <https://github.com/CrossRef/rest-api-doc/issues/146> (#115)
* Fix to `cr_journals()` in internal parsing, was failing in cases
where `ISSN` array was of length zero

rcrossref 0.5.4
===============

### MINOR IMPROVEMENTS

* Improved documentation for `cr_citation_count()` to remove PLOS
reference as the function isn't only for PLOS works (#108)
* Changed use of `dplyr::rbind_all()` to `dplyr::bind_rows()` (#113)

rcrossref 0.5.2
===============

### NEW FEATURES

* User-agent string is passed with every request to Crossref with
names and versions of this package, and its HTTP dependency packages: `httr`
and `curl` (which `httr` depends on). Will potentially be useful to
Crossref to know how many requests come from this R client (#100)

### DEPRECATED

* `cr_search()` and `cr_search_free()` use old Crossref web services, so
are now marked deprecated, and will throw a deprecation message, but can
still be used. They will both be defunct in `v0.6` of this package (#99)

### MINOR IMPROVEMENTS

* `XML` replaced with `xml2` (#98)
* `httr::content()` calls: all parse to text then parse content
manually. in addition, encoding explicitly set to `UTF-8` on
`httr::content()` calls (#98)

### BUG FIXES

* Bug fix to `cr_journals()` - fix to parse correctly on some failed requests
(#97) thanks @nkorf
* Bug fix to `cr_fundref()/cr_funders()` - parsing wasn't working correctly in
all cases

rcrossref 0.5.0
===============

Skipped `v0.4` to `v0.5` because of many changes - as described below.

### NEW FEATURES

* Support added for 'deep paging' with the newer Crossref search API. Two new params added to each function: `cursor`, which accepts a cursor alphanumeric string or the special `*`, which indicates that you want to initiate deep paging; `cursor_max`, which is not in the Crossref API, but just used here in this package to indicate where to stop - otherwise, you'd get all results, even if there was 70 million, for example. A new internal `R6` class used to make cursor requests easy (#77)
* New function `id_converter()` to get a PMID from a DOI and vice versa (#49)
* Added a Code of Conduct.
* New function `cr_types()`, along with its low level equivalent `cr_types_()` for when you just want a list or json back (#92)
* New suite of functions composing a low-level API for the Crossref search API. These functions only do data request, and return json or a list, and do not attempt to parse to a data.frame. The new functions: `cr_funders_()`, `cr_journals_()`, `cr_licenses_()`, `cr_members_()`, `cr_prefixes_()`, `cr_types_()`, `cr_works_()`. These functions are a bit faster, and aren't subject to parsing errors in the event of a change in the Crossref API. (#93)
* Added new `filter_names()` and `filter_details()` functions to get information on what filters
are available, the expected values, and what they mean.

### MINOR IMPROVEMENTS

* Added documentation for new filter types, and added them to list of filters for `filter_names()` and `filter_details()` (#73)
* `cr_funders()` alias added to `cr_fundref()` (#74)
* Made note in docs that funders without IDs don't show up on the `/funders` route,s in `cr_funders()` (#79)
* Made note in docs that `sample` parameter ignored unless `works=TRUE` (#81)
* Added note to docs that only what is returned in the API is what is searched when you search the Crossref API - that is, abstracts and full text aren't searched (#91)
* `cr_cn()` now checks that the user supplied content-type is supported for the DOI minting agency associated with the DOI (#88) (thanks @njahn82)
* Removed `.progress` parameter use internally where it wasn't applicable.
* `sample` parameter dropped from `cr_licenses()`.
* `cr_works()` parsing changed. We now don't attempt to flatten nested arrays, but instead give them back as data.frame's nested within the main data.frame. For example, `author` often has many entries, so we return that as a single column, but indexing to that column gives back a data.frame with a row for each author, and N number of columns. Hopefully this doesn't break too much code downstream :)
* Additional text added to the package level man file (`?rcrossref`) to explain: what you're actually searching when you search; deprecated and defunct functions; and explanation of high vs. low level API.

### BUG FIXES

* Fix to `cr_members()` to warn on error instead of stop during parsing (#68)
* Fix to internal parser for `cr_works()` to output links data, for full text links (#70)
* Minor fix in `cr_cn()` example that didn't work (#80)
* Fixed parsing of `affiliation` data inside `author` object in Crossref search API returned data (#84)
* Fixed parsing of funder `award` slot in Crossref search API returned data (#90)

### DEPRECATED

* `crosscite()` deprecated, will be removed in a future version of this package (#78)
* `cr_fundref()` now has a deprecated message, and will be removed in the next version (#74)

rcrossref 0.3.4
===============

### NEW FEATURES

* Added new function `crosscite()` to work with the
Citeproc service (http://crosscite.org/citeproc/) (#60)

### BUG FIXES

* Fixed problems related to `httr` `v1` (#65)
* Import all non-base R functions (#64)
* The agency route was down used by the `cr_agency()` function,
back up and fixed now (#63)

rcrossref 0.3.0
===============

### NEW FEATURES

* New function `extract_pdf()` to extract text from pdfs
* New function `cr_ft_links()` to get links for full text content of an article (#10)
* New function `cr_ft_text()` to get links for full text content of an article. In addition,
`cr_ft_pdf()`, `cr_ft_plain()`, and `cr_ft_xml()` are convenience functions that will get
the format pdf, plain text, or xml, respectively. You can of course specify format in the
`cr_ft_text()` function with the `type` parameter (#10) (#42)

### MINOR IMPROVEMENTS

* Filled out more tests (#45)

### BUG FIXES

* No longer assign queried doi to the `data.frame` in `cr_works()`, which caused failure if
a non-Crossref DOI included (#52)
* `pmid2doi()` and `doi2pmid()` functions removed temporarily as the web service is down
temporarily, but will be online again soon from Crossref (#48)

rcrossref 0.2.1
===============

### MINOR IMPROVEMENTS

* Fixes for man file examples. (#35)
* `cr_citation()` is deprecated (stil useable, but will be removed in a
future version of the package). use `cr_cn()` instead. (#34)

rcrossref 0.2.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.2 Patched
* ubuntu 16.04 LTS (on Travis-CI), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 8 downstream dependencies. There were no problems. See the summary at <https://github.com/ropensci/rcrossref/tree/master/revdep>

-------

This submission moves the bibtex package to Suggests as it has been orphaned, and a few minor improvements.

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

## Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/rcrossref/issues) - be sure to include R session information and a reproducible example.

## Code contributions

### Broad overview of contributing workflow

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rcrossref.git`
* Make sure to track progress upstream (i.e., on our version of `rcrossref` at `ropensci/rcrossref`) by doing `git remote add upstream https://github.com/ropensci/rcrossref.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs (see Tests below)
* Push up to your account
* Submit a pull request to home base at `ropensci/rcrossref`

### Tests

To add tests, go to the folder `tests/testthat/`. Tests are generally organized as individual files for each exported function from the package (that is, listed as an export in the `NAMESPACE` file). If you are adding a new exported function, add a new test file. If you are changing an existing function, work in the tests file for that function, unless it doesn't have tests, in which case make a new test file.

The book R packages book provides [a chapter on testing in general](http://r-pkgs.had.co.nz/tests.html). Do consult that first if you aren't familiar with testing in R.

The easiest set up to run tests is from within an R session:

```r
library(devtools)
library(testthat)
# loads the package
load_all()
```

To test an individual test file

```r
test_file("tests/testthat/test-foobar.R")
```

To run all tests

```r
devtools::test()
```

If you are running tests that have `skip_on_cran()` in them, set `Sys.setenv(NOT_CRAN = "true")` prior to running tests.


### Making changes

In addition to changing the code, do make sure to udpate the documentation if applicable. The R packages book book has a [chapter on documentation](http://r-pkgs.had.co.nz/man.html) you should read if you aren't familiar.

After code and documentation has been changed, update documentation by running either `devtools::document()` or `roxygen2::roxygenise()`.

Make sure if you change what packages or even functions within packages are imported, most likely add the package to Imports in the DESCRIPTION file and list what functions are imported in the `rcrossref-package.R` file.

Be conservative about adding new dependencies.


### Style

* Make sure code, documentation, and comments are no more than 80 characters in width.
* Use `<-` instead of `=` for assignment
* Always use `snake_case` (and all lowercase) instead of `camelCase`

## Check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.2 Patched (2020-09-01 r79114) |
|os       |macOS Catalina 10.15.6                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-09-30                                  |

# Dependencies

|package   |old     |new   |Δ  |
|:---------|:-------|:-----|:--|
|rcrossref |1.0.0   |1.1.0 |*  |
|bibtex    |0.4.2.3 |NA    |*  |

# Revdeps

*Wow, no problems at all. :)*Hi,

This is an automated email to let you know about the release of {{{ my_package }}}, which I'll submit to CRAN on ({{{ date }}}).

To check for potential problems, I ran `R CMD check` on your package {{{your_package}}} ({{{your_version}}}).

I found: {{{your_summary}}}.

{{#you_have_problems}}
{{{your_results}}}

If I got an ERROR because I couldn't install your package (or one of it's dependencies), my apologies. You'll have to run the checks yourself (unfortunately I don't have the time to diagnose installation failures as I have to run checks on hundreds of packages).

Otherwise, please carefully look at the results, and let me know if I've introduced a bug in {{{ my_package }}}. If I have, I'd really appreciate a minimal reproducible example that uses only {{{ my_package }}} functions. That way I can find and fix the bug as quickly as possible.

If it doesn't look like a bug in {{{ my_package }}}, please prepare an update for CRAN. Ideally you'll tweak your package so it works with both the released and development versions of dplyr. Otherwise, be prepared to submit your package to CRAN soon after I let you know that I've submitted.

To get the development version of {{{ my_package }}} so you can run the checks yourself, you can run:

    # install.packages("devtools")
    devtools::install_github("{{my_github}}")

To see what's changed visit <https://github.com/{{{my_github}}}/blob/master/NEWS.md>.

{{/you_have_problems}}
{{^you_have_problems}}
It looks like everything is ok, so you don't need to take any action, but you might want to read the NEWS, <https://github.com/{{{my_github}}}/blob/master/NEWS.md>, to see what's changed.
{{/you_have_problems}}

If you have any questions about this email, please feel free to respond directly.

Regards,

{{{ me }}}
*Wow, no problems at all. :)*rcrossref: R interface to CrossRef APIs
=======================================

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![cran checks](https://cranchecks.info/badges/worst/rcrossref)](https://cranchecks.info/pkgs/rcrossref)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/rcrossref/workflows/R-check/badge.svg)](https://github.com/ropensci/rcrossref/actions/)
[![codecov](https://codecov.io/gh/ropensci/rcrossref/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rcrossref)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcrossref)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcrossref)](https://cran.r-project.org/package=rcrossref)

## CrossRef documentation

* Crossref API: https://github.com/CrossRef/rest-api-doc#readme
* Crossref's API issue tracker: https://gitlab.com/crossref/issues
* Crossref metadata search API: https://www.crossref.org/labs/crossref-metadata-search/
* Crossref DOI Content Negotiation: https://citation.crosscite.org/docs.html
* Crossref Text and Data Mining (TDM) Services: https://www.crossref.org/education/retrieve-metadata/rest-api/text-and-data-mining/

## Installation

Stable version from CRAN

```{r eval=FALSE}
install.packages("rcrossref")
```

Or development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/rcrossref")
```

Load `rcrossref`

```{r}
library('rcrossref')
```

## Register for the Polite Pool

If you are intending to access Crossref regularly you will want to send your email address with your queries. This has the advantage that queries are placed in the polite pool of servers. Including your email address is good practice as described in the Crossref documentation under Good manners (https://github.com/CrossRef/rest-api-doc#good-manners--more-reliable-service). The second advantage is that Crossref can contact you if there is a problem with a query.

Details on how to register your email in a call can be found at `?rcrossref-package`. To pass your email address to Crossref, simply store it as an environment variable in .Renviron like this:

Open file: `file.edit("~/.Renviron")`

Add email address to be shared with Crossref `crossref_email= "name@example.com"`

Save the file and restart your R session

To stop sharing your email when using rcrossref simply delete it from your .Renviron file. 

## Documentation

See https://docs.ropensci.org/rcrossref/ to get started

## Meta

* Please report any issues or bugs: https://github.com/ropensci/rcrossref/issues
* License: MIT
* Get citation information for `rcrossref` in R doing `citation(package = 'rcrossref')`
* Please note that this package is released with a Contributor Code of Conduct (https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

---

This package is part of a richer suite called fulltext (https://github.com/ropensci/fulltext), along with several other packages, that provides the ability to search for and retrieve full text of open access scholarly articles.

---
---
title: rcrossref introduction 
author: Scott Chamberlain
date: "2021-08-19"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{rcrossref introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



## Installation

Install stable version from CRAN


```r
install.packages("rcrossref")
```

Or development version from GitHub


```r
remotes::install_github("ropensci/rcrossref")
```


```r
library("rcrossref")
```

## Citation search

CrossRef's DOI Content Negotiation service, where you can citations back in various formats, including `apa`


```r
cr_cn(dois = "10.1371/journal.pone.0112608", format = "text", style = "apa")
#> [1] "Wang, Q., & Taylor, J. E. (2014). Quantifying Human Mobility Perturbation and Resilience in Hurricane Sandy. PLoS ONE, 9(11), e112608. doi:10.1371/journal.pone.0112608"
```

There are a lot more styles. We include a dataset as a character vector within the package, accessible via the `get_styles()` function, e.g.,


```r
get_styles()[1:5]
#> [1] "academy-of-management-review"        
#> [2] "accident-analysis-and-prevention"    
#> [3] "aci-materials-journal"               
#> [4] "acm-sig-proceedings-long-author-list"
#> [5] "acm-sig-proceedings"
```

`bibtex`


```r
cat(cr_cn(dois = "10.1126/science.169.3946.635", format = "bibtex"))
#> @article{Frank_1970,
#> 	doi = {10.1126/science.169.3946.635},
#> 	url = {https://doi.org/10.1126%2Fscience.169.3946.635},
#> 	year = 1970,
#> 	month = {aug},
#> 	publisher = {American Association for the Advancement of Science ({AAAS})},
#> 	volume = {169},
#> 	number = {3946},
#> 	pages = {635--641},
#> 	author = {H. S. Frank},
#> 	title = {The Structure of Ordinary Water: New data and interpretations are yielding new insights into this fascinating substance},
#> 	journal = {Science}
#> }
```

`bibentry`


```r
cr_cn(dois = "10.6084/m9.figshare.97218", format = "bibentry")
#> $doi
#> [1] "10.6084/M9.FIGSHARE.97218"
#> 
#> $url
#> [1] "https://figshare.com/articles/thesis/Regime_shifts_in_ecology_and_evolution_(PhD_Dissertation)/97218"
#> 
#> $author
#> [1] "Boettiger, Carl"
#> 
#> $keywords
#> [1] "Evolutionary Biology, FOS: Biological sciences, FOS: Biological sciences, Ecology"
#> 
#> $title
#> [1] "Regime shifts in ecology and evolution (PhD Dissertation)"
#> 
#> $publisher
#> [1] "figshare"
#> 
#> $year
#> [1] "2012"
#> 
#> $copyright
#> [1] "Creative Commons Attribution 4.0 International"
#> 
#> $key
#> [1] "https://doi.org/10.6084/m9.figshare.97218"
#> 
#> $entry
#> [1] "article"
```

## Citation count

Citation count, using OpenURL


```r
cr_citation_count(doi="10.1371/journal.pone.0042793")
#>                            doi count
#> 1 10.1371/journal.pone.0042793    43
```

## Search Crossref metadata API

The following functions all use the CrossRef API.

### Look up funder information


```r
cr_funders(query="NSF")
#> $meta
#>   total_results search_terms start_index items_per_page
#> 1            39          NSF           0             20
#> 
#> $data
#> # A tibble: 20 × 6
#>    id           location       name         alt.names       uri     tokens      
#>    <chr>        <chr>          <chr>        <chr>           <chr>   <chr>       
#>  1 100000001    United States  National Sc… "USA NSF, NSF,… http:/… national, s…
#>  2 100015388    United States  Kansas NSF … "KNE, NSF EPSC… http:/… kansas, nsf…
#>  3 100016323    United States  Arkansas NS… "Arkansas EPSC… http:/… arkansas, n…
#>  4 100003187    United States  National Sl… "NSF"           http:/… national, s…
#>  5 501100000930 Australia      National St… "NSF"           http:/… national, s…
#>  6 501100004190 Norway         Norsk Sykep… "NSF, Norwegia… http:/… norsk, syke…
#>  7 501100020414 United Kingdom Neuroscienc… "The Neuroscie… http:/… neuroscienc…
#>  8 100000154    United States  Division of… "IOS, NSF Divi… http:/… division, o…
#>  9 100017338    China          Key Program… ""              http:/… key, progra…
#> 10 100000084    United States  Directorate… "NSF Directora… http:/… directorate…
#> 11 100016620    United States  Nick Simons… "NSF, The Nick… http:/… nick, simon…
#> 12 100006445    United States  Center for … "CHM, NSF, Uni… http:/… center, for…
#> 13 100017325    United States  Engineering… "ERC, The NSF … http:/… engineering…
#> 14 100008367    Denmark        Statens Nat… "Danish Nation… http:/… statens, na…
#> 15 100000179    United States  Office of t… "NSF Office of… http:/… office, of,…
#> 16 501100019492 China          National Na… "NSFC-General … http:/… national, n…
#> 17 501100011002 China          National Na… "NSFC-Yunnan J… http:/… national, n…
#> 18 501100008982 Sri Lanka      National Sc… "National Scie… http:/… national, s…
#> 19 501100001809 China          National Na… "NNSF of China… http:/… national, n…
#> 20 501100014220 China          National Na… "NSFC-Henan Jo… http:/… national, n…
#> 
#> $facets
#> NULL
```

### Check the DOI minting agency


```r
cr_agency(dois = '10.13039/100000001')
#> $DOI
#> [1] "10.13039/100000001"
#> 
#> $agency
#> $agency$id
#> [1] "crossref"
#> 
#> $agency$label
#> [1] "Crossref"
```

### Search works (i.e., articles, books, etc.)


```r
cr_works(filter=c(has_orcid=TRUE, from_pub_date='2004-04-04'), limit=1)
#> $meta
#>   total_results search_terms start_index items_per_page
#> 1       6556799           NA           0              1
#> 
#> $data
#> # A tibble: 1 × 26
#>   container.title created  deposited  published.online doi   indexed issn  issue
#>   <chr>           <chr>    <chr>      <chr>            <chr> <chr>   <chr> <chr>
#> 1 Belügyi Szemle  2021-04… 2021-04-26 2021-04-26       10.3… 2021-0… 2677… 4    
#> # … with 18 more variables: issued <chr>, member <chr>, page <chr>,
#> #   prefix <chr>, publisher <chr>, score <chr>, source <chr>,
#> #   reference.count <chr>, references.count <chr>,
#> #   is.referenced.by.count <chr>, title <chr>, type <chr>, url <chr>,
#> #   volume <chr>, abstract <chr>, short.container.title <chr>, author <list>,
#> #   link <list>
#> 
#> $facets
#> NULL
```

### Search journals


```r
cr_journals(issn=c('1803-2427','2326-4225'))
#> $data
#> # A tibble: 2 × 53
#>   title    publisher   issn   last_status_che… deposits_abstra… deposits_orcids…
#>   <chr>    <chr>       <chr>  <date>           <lgl>            <lgl>           
#> 1 Journal… De Gruyter… 1803-… 2021-08-18       TRUE             TRUE            
#> 2 Journal… American S… 2326-… 2021-08-18       FALSE            FALSE           
#> # … with 47 more variables: deposits <lgl>,
#> #   deposits_affiliations_backfile <lgl>,
#> #   deposits_update_policies_backfile <lgl>,
#> #   deposits_similarity_checking_backfile <lgl>,
#> #   deposits_award_numbers_current <lgl>,
#> #   deposits_resource_links_current <lgl>, deposits_articles <lgl>,
#> #   deposits_affiliations_current <lgl>, deposits_funders_current <lgl>, …
#> 
#> $facets
#> NULL
```

### Search license information


```r
cr_licenses(query = 'elsevier')
#> $meta
#>   total_results search_terms start_index items_per_page
#> 1            52           NA          NA             NA
#> 
#> $data
#> # A tibble: 20 × 2
#>    URL                                                                    work.count
#>    <chr>                                                                       <int>
#>  1 http://aspb.org/publications/aspb-journals/open-articles                        1
#>  2 http://creativecommons.org/licenses/by-nc-nd/3.0                                3
#>  3 http://creativecommons.org/licenses/by-nc-nd/3.0/                              11
#>  4 http://creativecommons.org/licenses/by-nc-nd/4.0/                              21
#>  5 http://creativecommons.org/licenses/by-nc/4.0/                                  6
#>  6 http://creativecommons.org/licenses/by-sa/4.0                                   1
#>  7 http://creativecommons.org/licenses/by/2.0                                      2
#>  8 http://creativecommons.org/licenses/by/3.0                                      2
#>  9 http://creativecommons.org/licenses/by/3.0/                                     2
#> 10 http://creativecommons.org/licenses/by/3.0/igo/                                 1
#> 11 http://creativecommons.org/licenses/by/4.0                                     11
#> 12 http://creativecommons.org/licenses/by/4.0/                                    22
#> 13 http://doi.wiley.com/10.1002/tdm_license_1                                    136
#> 14 http://doi.wiley.com/10.1002/tdm_license_1.1                                 2255
#> 15 http://iopscience.iop.org/info/page/text-and-data-mining                        2
#> 16 http://iopscience.iop.org/page/copyright                                        2
#> 17 http://journals.iucr.org/services/copyrightpolicy.html                         11
#> 18 http://journals.iucr.org/services/copyrightpolicy.html#TDM                     11
#> 19 http://journals.sagepub.com/page/policies/text-and-data-mining-license        365
#> 20 http://onlinelibrary.wiley.com/termsAndConditions                              63
```

### Search based on DOI prefixes


```r
cr_prefixes(prefixes=c('10.1016','10.1371','10.1023','10.4176','10.1093'))
#> $meta
#> NULL
#> 
#> $data
#>                               member                                    name
#> 1   http://id.crossref.org/member/78                             Elsevier BV
#> 2  http://id.crossref.org/member/340        Public Library of Science (PLoS)
#> 3  http://id.crossref.org/member/297 Springer Science and Business Media LLC
#> 4 http://id.crossref.org/member/1989                    Co-Action Publishing
#> 5  http://id.crossref.org/member/286           Oxford University Press (OUP)
#>                                  prefix
#> 1 http://id.crossref.org/prefix/10.1016
#> 2 http://id.crossref.org/prefix/10.1371
#> 3 http://id.crossref.org/prefix/10.1023
#> 4 http://id.crossref.org/prefix/10.4176
#> 5 http://id.crossref.org/prefix/10.1093
#> 
#> $facets
#> list()
```

### Search CrossRef members


```r
cr_members(query='ecology', limit = 5)
#> $meta
#>   total_results search_terms start_index items_per_page
#> 1            24      ecology           0              5
#> 
#> $data
#> # A tibble: 5 × 56
#>      id primary_name    location     last_status_che… current.dois backfile.dois
#>   <int> <chr>           <chr>        <date>           <chr>        <chr>        
#> 1  4302 Immediate Scie… Toronto, ON… 2021-08-18       0            6            
#> 2  6933 Knowledge Ecol… Washington,… 2021-08-18       0            1            
#> 3  1950 Journal of Vec… United Stat… 2021-08-18       0            0            
#> 4  2899 Association fo… Eugene, OR,… 2021-08-18       0            0            
#> 5  7745 Institute of A… Makhachkala… 2021-08-18       172          678          
#> # … with 50 more variables: total.dois <chr>, prefixes <chr>,
#> #   coverge.affiliations.current <chr>,
#> #   coverge.similarity.checking.current <chr>, coverge.funders.backfile <chr>,
#> #   coverge.licenses.backfile <chr>, coverge.funders.current <chr>,
#> #   coverge.affiliations.backfile <chr>, coverge.resource.links.backfile <chr>,
#> #   coverge.orcids.backfile <chr>, coverge.update.policies.current <chr>,
#> #   coverge.open.references.backfile <chr>, coverge.orcids.current <chr>, …
#> 
#> $facets
#> NULL
```

### Get N random DOIs

`cr_r()` uses the function `cr_works()` internally.


```r
cr_r()
#>  [1] "10.1111/biot.1990.4.issue-2"           
#>  [2] "10.7717/peerj.5714/table-5"            
#>  [3] "10.1061/(asce)0733-947x(2008)134:1(34)"
#>  [4] "10.1016/j.clml.2015.04.108"            
#>  [5] "10.1371/journal.pntd.0009063.s016"     
#>  [6] "10.31616/asj.2020.0425"                
#>  [7] "10.15446/revfacmed.v65n3.61313"        
#>  [8] "10.1586/ecp.09.36"                     
#>  [9] "10.1007/978-3-8348-9174-7_1"           
#> [10] "10.1080/09511920601160171"
```

You can pass in the number of DOIs you want back (default is 10)


```r
cr_r(2)
#> [1] "10.5194/gmd-2021-94-supplement" "10.1108/prt.2009.12938bad.009"
```

## Get full text

Publishers can optionally provide links in the metadata they provide to Crossref for full text of the work, but that data is often missing. Find out more about it at https://support.crossref.org/hc/en-us/articles/215750183-Crossref-Text-and-Data-Mining-Services

Get some DOIs for articles that provide full text, and that have `CC-BY 3.0` licenses (i.e., more likely to actually be open)


```r
out <-
  cr_works(filter = list(has_full_text = TRUE,
    license_url = "http://creativecommons.org/licenses/by/3.0/"))
(dois <- out$data$doi)
#>  [1] "10.1155/jamsa/2006/42542"        "10.1016/s0370-2693(02)01651-9"  
#>  [3] "10.1016/s0370-2693(02)01624-6"   "10.1155/ijmms/2006/89545"       
#>  [5] "10.1016/s0370-2693(01)01058-9"   "10.1016/s0370-2693(01)01257-6"  
#>  [7] "10.1016/s0370-2693(01)01287-4"   "10.1016/s0370-2693(01)01385-5"  
#>  [9] "10.1002/cfg.80"                  "10.1002/cfg.118"                
#> [11] "10.1002/cfg.166"                 "10.1002/cfg.59"                 
#> [13] "10.1002/cfg.108"                 "10.1088/1742-6596/1147/1/012044"
#> [15] "10.1088/1742-6596/1147/1/012077" "10.1088/1742-6596/578/1/012006" 
#> [17] "10.1088/1755-1315/222/1/012020"  "10.1002/ecs2.2575"              
#> [19] "10.1088/1755-1315/214/1/012004"  "10.1088/1755-1315/214/1/012021"
```

From the output of `cr_works` we can get full text links if we know where to look:


```r
do.call("rbind", out$data$link)
#> # A tibble: 55 × 4
#>    URL                           content.type  content.version intended.applica…
#>    <chr>                         <chr>         <chr>           <chr>            
#>  1 http://downloads.hindawi.com… application/… vor             text-mining      
#>  2 http://downloads.hindawi.com… unspecified   vor             similarity-check…
#>  3 https://api.elsevier.com/con… text/xml      vor             text-mining      
#>  4 https://api.elsevier.com/con… text/plain    vor             text-mining      
#>  5 https://api.elsevier.com/con… text/xml      vor             text-mining      
#>  6 https://api.elsevier.com/con… text/plain    vor             text-mining      
#>  7 http://downloads.hindawi.com… application/… vor             text-mining      
#>  8 http://downloads.hindawi.com… unspecified   vor             similarity-check…
#>  9 https://api.elsevier.com/con… text/xml      vor             text-mining      
#> 10 https://api.elsevier.com/con… text/plain    vor             text-mining      
#> # … with 45 more rows
```

From there, you can grab your full text, but because most links require
authentication, enter another package: `crminer`.

You'll need package `crminer` for the rest of the work.

Onc we have DOIs, get URLs to full text content


```r
if (!requireNamespace("crminer")) {
  install.packages("crminer")
}
```


```r
library(crminer)
#> Error in library(crminer): there is no package called 'crminer'
(links <- crm_links("10.1155/2014/128505"))
#> Error in crm_links("10.1155/2014/128505"): could not find function "crm_links"
```

Then use those URLs to get full text


```r
crm_pdf(links)
#> <document>/Users/sckott/Library/Caches/R/crminer/128505.pdf
#>   Pages: 1
#>   No. characters: 1565
#>   Created: 2014-09-15
```

See also fulltext (https://github.com/ropensci/fulltext) for getting scholarly text 
for text mining.
---
title: Crossref filters
author: Scott Chamberlain
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Crossref filters}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Crossref Filter Names

These are associated with the Fundref API

Since these Crossref options use dashes in the names in the API, we provide a crosswalk to these so that you can use underscores wherever they use dashes. Also, filters that have dots in them are undscores instead, so `license.url` becomes `license_url`. Additionally, using colons is not very R like, so if Crossref requires `has-funder:true`, just do `has-funder=TRUE` in R. Pass these options in to the `cr_*()` functions in the filter parameter as a list, e.g., `filter=list(has_funder=TRUE, has_full_text=TRUE)`

Filters allow you to narrow queries. All filter results are lists.  The following filters are supported:

| filter     | possible values | description|
|:-----------|:----------------|:-----------|
| `has_funder` | `{logical}` | metadata which includes one or more funder entry |
| `funder` | `{funder_id}` | metadata which include the `{funder_id}` in FundRef data |
| `prefix` | `{owner_prefix}` | metadata belonging to a DOI owner prefix `{owner_prefix}` (e.g. `10.1016` ) |
| `member` | `{member_id}` | metadata belonging to a CrossRef member |
| `from_index_date` | `{date}` | metadata indexed since (inclusive) `{date}` |
| `until_index_date` | `{date}` | metadata indexed before (inclusive) `{date}` |
| `from_deposit_date` | `{date}` | metadata last (re)deposited since (inclusive) `{date}` |
| `until_deposit_date` | `{date}` | metadata last (re)deposited before (inclusive) `{date}` |
| `from_update_date` | `{date}` | Metadata updated since (inclusive) `{date}`. Currently the same as `from_deposit_date`. |
| `until_update_date` | `{date}` | Metadata updated before (inclusive) `{date}`. Currently the same as `until_deposit_date`. |
| `from_first_deposit_date` | `{date}` | metadata first deposited since (inclusive) `{date}` [^*] |
| `until_first_deposit_date` | `{date}` | metadata first deposited before (inclusive) `{date}` [^*] |
| `from_pub_date` | `{date}` | metadata where published date is since (inclusive) `{date}` |
| `until_pub_date` | `{date}` | metadata where published date is before (inclusive)  `{date}` |
| `has_license` | `{logical}` | metadata that includes any `<license_ref>` elements. |
| `license_url` | `{url}` | metadata where `<license_ref>` value equals `{url}` |
| `license_version` | `{string}` | metadata where the `<license_ref>`'s `applies_to` attribute  is `{string}`|
| `license_delay` | `{integer}` | metadata where difference between publication date and the `<license_ref>`'s `start_date` attribute is <= `{integer}` (in days)|
| `has_full_text` | `{logical}` | metadata that includes any full text `<resource>` elements. |
| `full_text_version` | `{string}`  | metadata where `<resource>` element's `content_version` attribute is `{string}`. |
| `full_text_type` | `{mime_type}`  | metadata where `<resource>` element's `content_type` attribute is `{mime_type}` (e.g. `application/pdf`). |
| `public_references` | | metadata where publishers allow references to be distributed publicly. [^*] |
| `has_references` | `{logical}` | metadata for works that have a list of references |
| `has_archive` | `{logical}` | metadata which include name of archive partner |
| `archive` | `{string}` | metadata which where value of archive partner is `{string}` |
| `has_orcid` | `{logical}` | metadata which includes one or more ORCIDs |
| `orcid` | `{orcid}` | metadata where `<orcid>` element's value = `{orcid}` |
| `issn` | `{issn}` | metadata where record has an ISSN = `{issn}`. Format is `xxxx-xxxx`. |
| `type` | `{type}` | metadata records whose type = `{type}`. Type must be an ID value from the list of types returned by the `/types` resource |
| `directory` | `{directory}` | metadata records whose article or serial are mentioned in the given `{directory}`. Currently the only supported value is `doaj`. |
| `doi` | `{doi}` | metadata describing the DOI `{doi}` |
| `updates` | `{doi}` | metadata for records that represent editorial updates to the DOI `{doi}` |
| `is_update` | `{logical}` | metadata for records that represent editorial updates |
| `has_update_policy` | `{logical}` | metadata for records that include a link to an editorial update policy |
| `container_title` | | metadata for records with a publication title exactly with an exact match |
| `publisher_name` | | metadata for records with an exact matching publisher name |
| `category_name` | |metadata for records with an exact matching category label |
| `type_name` | | metadata for records with an exactly matching type label |
| `award_number` |{`award_number`} |metadata for records with a matching award number. Optionally combine with `award_funder` |
| `award_funder` | {`funder doi or id`} | metadata for records with an award with matching funder. Optionally combine with `award_number` |
| `from_created_date` | {`date`} | metadata first deposited since (inclusive) `{date}` |
| `until_created_date` | {`date`} | metadata first deposited before (inclusive) `{date}` |
| `affiliation` | | metadata for records with at least one contributor with the given affiliation |
| `has_affiliation`| | metadata for records that have any affiliation information |
| `article_number` | | metadata for records with a given article number |
| `alternative_id` | | metadata for records with the given alternative ID, which may be a publisher_specific ID, or any other identifier a publisher may have provided |
| `assertion_group` | | metadata for records with an assertion in a particular group |
| `assertion` | | metadata for records with a particular named assertion |

[^*]: Not implemented yet.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_cn.r
\name{cr_cn}
\alias{cr_cn}
\title{Get citations in various formats from CrossRef.}
\usage{
cr_cn(
  dois,
  format = "bibtex",
  style = "apa",
  locale = "en-US",
  raw = FALSE,
  .progress = "none",
  url = NULL,
  ...
)
}
\arguments{
\item{dois}{Search by a single DOI or many DOIs.}

\item{format}{Name of the format. One of "rdf-xml", "turtle",
"citeproc-json", "citeproc-json-ish", "text", "ris", "bibtex" (default),
"crossref-xml", "datacite-xml","bibentry", or "crossref-tdm". The format
"citeproc-json-ish" is a format that is not quite proper citeproc-json.
Note that the package bibtex is required when \code{format="bibtex"}; the
package is in Suggests so is not required when installing rcrossref}

\item{style}{a CSL style (for text format only). See \code{\link[=get_styles]{get_styles()}}
for options. Default: 'apa'. If there's a style that CrossRef doesn't support
you'll get a  \verb{(500) Internal Server Error}}

\item{locale}{Language locale. See \code{?Sys.getlocale}}

\item{raw}{(logical) Return raw text in the format given by \code{format}
parameter. Default: \code{FALSE}}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{url}{(character) Base URL for the content negotiation request.
Default: "https://doi.org"}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}
}
\description{
Get citations in various formats from CrossRef.
}
\details{
See http://citation.crosscite.org/docs.html for more info
on the Crossref Content Negotiation API service.

DataCite DOIs: Some values of the \code{format} parameter won't work with
DataCite DOIs, i.e. "citeproc-json", "crossref-xml", "crossref-tdm",
"onix-xml".

MEDRA DOIs only work with "rdf-xml", "turtle", "citeproc-json-ish", "ris",
"bibtex", "bibentry", "onix-xml".

See examples below.

See \code{\link[=cr_agency]{cr_agency()}}

Note that the format type \code{citeproc-json} uses the CrossRef API at
\code{api.crossref.org}, while all others are content negotiated via
\verb{http://data.crossref.org}, \verb{http://data.datacite.org} or
\verb{http://data.medra.org}. DOI agency is checked first (see
\code{\link[=cr_agency]{cr_agency()}}).
}
\examples{
\dontrun{
cr_cn(dois="10.1126/science.169.3946.635")
cr_cn(dois="10.1126/science.169.3946.635", "citeproc-json")
cr_cn(dois="10.1126/science.169.3946.635", "citeproc-json-ish")
cr_cn("10.1126/science.169.3946.635", "rdf-xml")
cr_cn("10.1126/science.169.3946.635", "crossref-xml")
cr_cn("10.1126/science.169.3946.635", "text")

# return an R bibentry type
cr_cn("10.1126/science.169.3946.635", "bibentry")
cr_cn("10.6084/m9.figshare.97218", "bibentry")

# return an apa style citation
cr_cn("10.1126/science.169.3946.635", "text", "apa")
cr_cn("10.1126/science.169.3946.635", "text", "harvard3")
cr_cn("10.1126/science.169.3946.635", "text", "elsevier-harvard")
cr_cn("10.1126/science.169.3946.635", "text", "ecoscience")
cr_cn("10.1126/science.169.3946.635", "text", "heredity")
cr_cn("10.1126/science.169.3946.635", "text", "oikos")

# example with many DOIs
dois <- cr_r(2)
cr_cn(dois, "text", "apa")

# Cycle through random styles - print style on each try
stys <- get_styles()
foo <- function(x){
 cat(sprintf("<Style>:\%s\n", x), sep = "\n\n")
 cat(cr_cn("10.1126/science.169.3946.635", "text", style=x))
}
foo(sample(stys, 1))

# Using DataCite DOIs
## some formats don't work
# cr_cn("10.5284/1011335", "crossref-xml")
# cr_cn("10.5284/1011335", "crossref-tdm")
## But most do work
cr_cn("10.5284/1011335", "text")
cr_cn("10.5284/1011335", "datacite-xml")
cr_cn("10.5284/1011335", "rdf-xml")
cr_cn("10.5284/1011335", "turtle")
cr_cn("10.5284/1011335", "citeproc-json-ish")
cr_cn("10.5284/1011335", "ris")
cr_cn("10.5284/1011335", "bibtex")
cr_cn("10.5284/1011335", "bibentry")

# Using Medra DOIs
cr_cn("10.3233/ISU-150780", "onix-xml")

# Get raw output
cr_cn(dois = "10.1002/app.27716", format = "citeproc-json", raw = TRUE)

# sometimes messy DOIs even work
## in this case, a DOI minting agency can't be found
## but we proceed anyway, just assuming it's "crossref"
cr_cn("10.1890/0012-9615(1999)069[0569:EDILSA]2.0.CO;2")

# Use a different base url
cr_cn("10.1126/science.169.3946.635", "text", url = "https://data.datacite.org")
cr_cn("10.1126/science.169.3946.635", "text", url = "http://dx.doi.org")
cr_cn("10.1126/science.169.3946.635", "text", "heredity", url = "http://dx.doi.org")
cr_cn("10.5284/1011335", url = "https://citation.crosscite.org/format", 
   style = "oikos")
cr_cn("10.5284/1011335", url = "https://citation.crosscite.org/format", 
   style = "plant-cell-and-environment")
cr_cn("10.5284/1011335", url = "https://data.datacite.org", 
   style = "plant-cell-and-environment")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_defunct.R
\name{cr_fundref}
\alias{cr_fundref}
\title{fundref}
\usage{
cr_fundref(...)
}
\description{
fundref
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_types.R
\name{cr_types}
\alias{cr_types}
\alias{cr_types_}
\title{Search CrossRef types}
\usage{
cr_types(
  types = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  flq = NULL,
  select = NULL,
  ...
)

cr_types_(
  types = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  parse = FALSE,
  flq = NULL,
  select = NULL,
  ...
)
}
\arguments{
\item{types}{(character) Type identifier, e.g., journal}

\item{query}{Query terms}

\item{filter}{Filter options. See examples for usage examples
and \code{\link{filters}} for what filters are available.
\code{filter} is available for use with \code{cr_works} and
other \code{crossref} family functions with \code{works=TRUE}}

\item{offset}{Number of record to start at. Minimum: 1. For
\code{\link{cr_works}}, and any function setting \code{works = TRUE},
the maximum offset value is 10000. For larger requests use \code{cursor}.}

\item{limit}{Number of results to return in the query. Not relavant when
searching with specific dois. Default: 20. Max: 1000}

\item{sample}{(integer) Number of random results to return. when you use
the sample parameter, the rows and offset parameters are ignored.
Ignored unless \code{works=TRUE}. Max: 100}

\item{sort}{Field to sort on. Acceptable set of fields to sort on:
\itemize{
\item \code{score} OR \code{relevance} - Sort by relevance score
\item \code{updated} - Sort by date of most recent change to metadata.
Currently the same as deposited.
\item \code{deposited} - Sort by time of most recent deposit
\item \code{indexed} - Sort by time of most recent index
\item \code{published} - Sort by publication date
\item \code{published-print} - Sort by print publication date
\item \code{published-online} - Sort by online publication date
\item \code{issued} - Sort by issued date (earliest known publication date)
\item \code{is-referenced-by-count} - Sort by number of times this DOI is
referenced by other Crossref DOIs
\item \code{references-count} - Sort by number of references included in
the references section of the document identified by this DOI
}}

\item{order}{(character) Sort order, one of 'asc' or 'desc'}

\item{facet}{(logical) Include facet results. Boolean or string with
field to facet on. Valid fields are *, affiliation, funder-name,
funder-doi, orcid, container-title, assertion, archive, update-type,
issn, published, source, type-name, publisher-name, license,
category-name, assertion-group. Default: \code{FALSE}}

\item{works}{(logical) If \code{TRUE}, works returned as well, if not then not.}

\item{cursor}{(character) Cursor character string to do deep paging.
Default is None. Pass in '*' to start deep paging. Any combination of
query, filters and facets may be used with deep paging cursors.
While the \code{limit} parameter may be specified along with cursor,
offset and sample cannot be used with the cursor. See
https://github.com/CrossRef/rest-api-doc#deep-paging-with-cursors}

\item{cursor_max}{(integer) Max records to retrieve. Only used when
cursor param used. Because deep paging can result in continuous requests
until all are retrieved, use this parameter to set a maximum number of
records. Of course, if there are less records found than this value,
you will get only those found. When cursor pagination is being used
the \code{limit} parameter sets the chunk size per request.}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{flq}{field queries. One or more field queries. Acceptable set of
field query parameters are:
\itemize{
\item \code{query.container-title}    - Query container-title aka.
publication name
\item \code{query.author} - Query author first and given names
\item \code{query.editor} - Query editor first and given names
\item \code{query.chair}    - Query chair first and given names
\item \code{query.translator} - Query translator first and given names
\item \code{query.contributor} - Query author, editor, chair and
translator first and given names
\item \code{query.bibliographic} - Query bibliographic information, useful
for citation lookup. Includes titles, authors, ISSNs and publication years
\item \code{query.affiliation} - Query contributor affiliations
}

Note: \code{query.title} has been removed - use \code{query.bibliographic}
as a replacement}

\item{select}{(character) One or more field to return (only those fields
are returned)}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}

\item{parse}{(logical) Whether to output json \code{FALSE} or parse to
list \code{TRUE}. Default: \code{FALSE}}
}
\description{
Search CrossRef types
}
\details{
BEWARE: The API will only work for CrossRef DOIs.
}
\note{
See the "Rate limiting" seciton in \link{rcrossref} to get
into the "fast lane"
}
\section{Deep paging (using the cursor)}{

When using the cursor, a character string called \code{next-cursor} is
returned from Crossref that we use to do the next request, and so on. We
use a while loop to get number of results up to the value of
\code{cursor_max}. Since we are doing each request for you, you may not
need the \code{next-cursor} string, but if you do want it, you can get
to it by indexing into the result like \code{x$meta$next_cursor}

Note that you can pass in curl options when using cursor, via \code{"..."}
}

\examples{
\dontrun{
cr_types()
cr_types("monograph")
cr_types("monograph", works=TRUE)
cr_types(c('monograph', 'book-set', 'book', 'book-track'))
cr_types(c('monograph', 'book-set'), works=TRUE)

## get facets back
cr_types("journal-article", works=TRUE, facet=TRUE)$facets
cr_types("monograph", works=TRUE, facet="license:*", limit = 0)
cr_types(c('monograph', 'book-set'), works=TRUE, facet=TRUE)

# Use the cursor for deep paging
cr_types("journal-article", works = TRUE, cursor = "*",
   cursor_max = 500, limit = 100)
cr_types(c('monograph', 'book-set'), works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
## with optional progress bar
cr_types("journal-article", works = TRUE, cursor = "*",
   cursor_max = 500, limit = 100, .progress = TRUE)

# query doesn't work unless using works=TRUE
### you get results, but query param is silently dropped
cr_types(query = "ecology")

# print progress - only works when passing more than one type
cr_types(c('monograph', 'book-set'), works=TRUE, .progress='text')

# Low level function - does no parsing to data.frame, get json or a list
cr_types_('monograph')
cr_types_('monograph', parse = TRUE)
cr_types_("journal-article", works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)

# field queries
## query.container-title
cr_types("journal-article", works = TRUE,
  flq = c(`query.container-title` = 'Ecology'))

# select only certain fields to return
res <- cr_types("journal-article", works = TRUE, select = c('DOI', 'title'))
names(res$data)
}
}
\references{
https://github.com/CrossRef/rest-api-doc
}
\seealso{
Other crossref: 
\code{\link{cr_funders}()},
\code{\link{cr_journals}()},
\code{\link{cr_licenses}()},
\code{\link{cr_members}()},
\code{\link{cr_prefixes}()},
\code{\link{cr_works}()}
}
\concept{crossref}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_ft_text.R
\name{cr_ft_xml}
\alias{cr_ft_xml}
\title{Get full text xml from a DOI}
\usage{
cr_ft_xml(...)
}
\description{
Get full text xml from a DOI
}
\note{
see \code{crminer::crm_xml}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_cn.r
\name{GET_agency_id}
\alias{GET_agency_id}
\title{Get doi agency id to identify resource location}
\usage{
GET_agency_id(x, ...)
}
\arguments{
\item{x}{doi}
}
\description{
Get doi agency id to identify resource location
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filters}
\alias{filters}
\alias{filter_names}
\alias{filter_details}
\title{Get filter details and names.}
\usage{
filter_names()

filter_details(x = NULL)
}
\arguments{
\item{x}{(character) Optional filter name. If not given, all filters
returned.}
}
\description{
Get filter details and names.
}
\details{
Note that all filter names in this package have periods
and dashes replaced with underscores as they both cause problems
in an R console.
}
\examples{
filter_names()
filter_details()
filter_details()$has_authenticated_orcid
filter_details()$has_authenticated_orcid$possible_values
filter_details()$has_authenticated_orcid$description
filter_details("issn")
filter_details("iss")
filter_details(c("issn", "alternative_id"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_styles.R
\name{get_styles}
\alias{get_styles}
\title{Get list of styles from github.com/citation-style-language/styles}
\usage{
get_styles(...)
}
\arguments{
\item{...}{Named parameters passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Get list of styles from github.com/citation-style-language/styles
}
\examples{
\dontrun{
x <- get_styles()
x[1:10]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_citation.R
\name{cr_citation}
\alias{cr_citation}
\title{This function is defunct}
\usage{
cr_citation(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_defunct.R
\name{rcrossref-defunct}
\alias{rcrossref-defunct}
\title{Defunct functions in rcrossref}
\description{
These functions are gone, no longer available.
}
\details{
\itemize{
\item \code{\link[=cr_citation]{cr_citation()}}: Crossref is trying to sunset their
OpenURL API, which  this function uses. So this function is now removed.
See the function \code{\link[=cr_cn]{cr_cn()}}, which does the same things, but with
more functionality, using the new Crossref API.
\item \code{\link[=pmid2doi]{pmid2doi()}} and \code{\link[=doi2pmid]{doi2pmid()}}: The API behind
these functions is down for good, see \code{\link[=id_converter]{id_converter()}} for
similar functionality.
\item \code{\link[=cr_search]{cr_search()}}: The functionality of this function can
be achieved with the new Crossref API. See functions \code{\link[=cr_works]{cr_works()}}
et al.
\item \code{\link[=cr_search_free]{cr_search_free()}}: The functionality of this function can
be achieved with the new Crossref API. See functions \code{\link[=cr_works]{cr_works()}}
et al.
\item \code{\link[=crosscite]{crosscite()}}: The functionality of this function can be
achieved with \code{\link[=cr_cn]{cr_cn()}}
\item \code{\link[=cr_fundref]{cr_fundref()}}: Crossref changed their name "fundref"
to "funders", so we've changed our function, see \code{\link[=cr_funders]{cr_funders()}}
\item \code{\link[=cr_ft_text]{cr_ft_text()}}: This function and other text mining
functions are incorporated in a new package \code{crminer}.
\item \code{\link[=cr_ft_links]{cr_ft_links()}}: This function and other text mining
functions are incorporated in a new package \code{crminer}.
\item \code{\link[=cr_ft_pdf]{cr_ft_pdf()}}: This function and other text mining
functions are  incorporated in a new package \code{crminer}.
\item \code{\link[=cr_ft_plain]{cr_ft_plain()}}: This function and other text mining
functions are incorporated in a new package \code{crminer}.
\item \code{\link[=cr_ft_text]{cr_ft_text()}}: This function and other text mining
functions are incorporated in a new package \code{crminer}.
\item \code{\link[=cr_ft_xml]{cr_ft_xml()}}: This function and other text mining
functions are incorporated in a new package \code{crminer}.
\item \code{\link[=as.tdmurl]{as.tdmurl()}}: This function and other text mining
functions are incorporated in a new package \code{crminer}.
\item \code{\link[=extract_xpdf]{extract_xpdf()}}: This function and other text mining
functions are incorporated in a new package \code{crminer}.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_defunct.R
\name{rcrossref-deprecated}
\alias{rcrossref-deprecated}
\title{Deprecated functions in rcrossref}
\description{
None at the moment
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_pdf.R
\name{extract_xpdf}
\alias{extract_xpdf}
\title{Extract text from a single pdf document.}
\usage{
extract_xpdf(...)
}
\description{
Extract text from a single pdf document.
}
\note{
see \code{crminer::crm_extract}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/id_converter.R
\name{id_converter}
\alias{id_converter}
\title{Get a PMID from a DOI, and vice versa.}
\usage{
id_converter(x, type = NULL, ...)
}
\arguments{
\item{x}{(character) One or more of: doi, pmid, pmcid, or manuscript id,
see examples. required.}

\item{type}{(character) one of doi, pmid, pmcid, or manuscript id}

\item{...}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Get a PMID from a DOI, and vice versa.
}
\examples{
\dontrun{
# get a pmid/pmcid from a doi
id_converter("10.1038/ng.590")

# pmid to doi/pmcid
id_converter("20495566", "pmid")
id_converter("20495566")
# id_converter("20495566", "doi") #error

# pmcid to doi/pmid
id_converter("PMC2883744", "pmcid")
id_converter("PMC2883744")

# manuscript id
id_converter("NIHMS311352")

# more than 1 ID
id_converter(c("PMC3531190","PMC3245039"))

# error, wrong type passed for id given
# id_converter("PMC2883744", "doi")
# error, 200 ids or less
# ids <- cr_r(100)
# id_converter(c(ids, ids, ids))
}
}
\references{
https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_search.r
\name{cr_search}
\alias{cr_search}
\title{Search the CrossRef Metadata API.}
\usage{
cr_search(...)
}
\description{
Search the CrossRef Metadata API.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_ft_links.R
\name{cr_ft_links}
\alias{cr_ft_links}
\title{Get full text links from a DOI}
\usage{
cr_ft_links(...)
}
\description{
Get full text links from a DOI
}
\note{
see \code{crminer::crm_links}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_ft_text.R
\name{cr_ft_plain}
\alias{cr_ft_plain}
\title{Get full text plain from a DOI}
\usage{
cr_ft_plain(...)
}
\description{
Get full text plain from a DOI
}
\note{
see \code{crminer::crm_plain}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_members.r
\name{cr_members}
\alias{cr_members}
\alias{cr_members_}
\title{Search CrossRef members}
\usage{
cr_members(
  member_ids = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  flq = NULL,
  select = NULL,
  ...
)

cr_members_(
  member_ids = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  parse = FALSE,
  flq = NULL,
  select = NULL,
  ...
)
}
\arguments{
\item{member_ids}{One or more member ids. See examples. Alternatively,
you can query for them using the query parameter.}

\item{query}{Query terms}

\item{filter}{Filter options. See examples for usage examples
and \code{\link{filters}} for what filters are available.
\code{filter} is available for use with \code{cr_works} and
other \code{crossref} family functions with \code{works=TRUE}}

\item{offset}{Number of record to start at. Minimum: 1. For
\code{\link{cr_works}}, and any function setting \code{works = TRUE},
the maximum offset value is 10000. For larger requests use \code{cursor}.}

\item{limit}{Number of results to return in the query. Not relavant when
searching with specific dois. Default: 20. Max: 1000}

\item{sample}{(integer) Number of random results to return. when you use
the sample parameter, the rows and offset parameters are ignored.
Ignored unless \code{works=TRUE}. Max: 100}

\item{sort}{Field to sort on. Acceptable set of fields to sort on:
\itemize{
\item \code{score} OR \code{relevance} - Sort by relevance score
\item \code{updated} - Sort by date of most recent change to metadata.
Currently the same as deposited.
\item \code{deposited} - Sort by time of most recent deposit
\item \code{indexed} - Sort by time of most recent index
\item \code{published} - Sort by publication date
\item \code{published-print} - Sort by print publication date
\item \code{published-online} - Sort by online publication date
\item \code{issued} - Sort by issued date (earliest known publication date)
\item \code{is-referenced-by-count} - Sort by number of times this DOI is
referenced by other Crossref DOIs
\item \code{references-count} - Sort by number of references included in
the references section of the document identified by this DOI
}}

\item{order}{(character) Sort order, one of 'asc' or 'desc'}

\item{facet}{(logical) Include facet results. Boolean or string with
field to facet on. Valid fields are *, affiliation, funder-name,
funder-doi, orcid, container-title, assertion, archive, update-type,
issn, published, source, type-name, publisher-name, license,
category-name, assertion-group. Default: \code{FALSE}}

\item{works}{(logical) If \code{TRUE}, works returned as well, if not then not.}

\item{cursor}{(character) Cursor character string to do deep paging.
Default is None. Pass in '*' to start deep paging. Any combination of
query, filters and facets may be used with deep paging cursors.
While the \code{limit} parameter may be specified along with cursor,
offset and sample cannot be used with the cursor. See
https://github.com/CrossRef/rest-api-doc#deep-paging-with-cursors}

\item{cursor_max}{(integer) Max records to retrieve. Only used when
cursor param used. Because deep paging can result in continuous requests
until all are retrieved, use this parameter to set a maximum number of
records. Of course, if there are less records found than this value,
you will get only those found. When cursor pagination is being used
the \code{limit} parameter sets the chunk size per request.}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{flq}{field queries. One or more field queries. Acceptable set of
field query parameters are:
\itemize{
\item \code{query.container-title}    - Query container-title aka.
publication name
\item \code{query.author} - Query author first and given names
\item \code{query.editor} - Query editor first and given names
\item \code{query.chair}    - Query chair first and given names
\item \code{query.translator} - Query translator first and given names
\item \code{query.contributor} - Query author, editor, chair and
translator first and given names
\item \code{query.bibliographic} - Query bibliographic information, useful
for citation lookup. Includes titles, authors, ISSNs and publication years
\item \code{query.affiliation} - Query contributor affiliations
}

Note: \code{query.title} has been removed - use \code{query.bibliographic}
as a replacement}

\item{select}{(character) One or more field to return (only those fields
are returned)}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}

\item{parse}{(logical) Whether to output json \code{FALSE} or parse to
list \code{TRUE}. Default: \code{FALSE}}
}
\description{
Search CrossRef members
}
\details{
BEWARE: The API will only work for CrossRef DOIs.
}
\note{
See the "Rate limiting" seciton in \link{rcrossref} to get
into the "fast lane"
}
\section{Deep paging (using the cursor)}{

When using the cursor, a character string called \code{next-cursor} is
returned from Crossref that we use to do the next request, and so on. We
use a while loop to get number of results up to the value of
\code{cursor_max}. Since we are doing each request for you, you may not
need the \code{next-cursor} string, but if you do want it, you can get
to it by indexing into the result like \code{x$meta$next_cursor}

Note that you can pass in curl options when using cursor, via \code{"..."}
}

\examples{
\dontrun{
cr_members(member_ids=98)
cr_members(member_ids=340)

cr_members(member_ids=98, works=TRUE)
cr_members(member_ids=c(10,98,45,1,9))
cr_members(member_ids=c(10,98,45,1,9), works=TRUE)

cr_members(query='hindawi')
cr_members(query='ecology')

# facets
cr_members(member_ids=98, works=TRUE, facet=TRUE, limit = 0)
cr_members(member_ids=98, works=TRUE, facet="license:*", limit = 0)

# curl options
cr_members(member_ids=98, verbose = TRUE)

# Use the cursor for deep paging
cr_members(member_ids=98, works = TRUE, cursor = "*",
   cursor_max = 500, limit = 100)
cr_members(member_ids=c(10, 98, 45), works = TRUE, cursor = "*",
   cursor_max = 200, limit = 100)
## with optional progress bar
cr_members(member_ids=98, works = TRUE, cursor = "*",
   cursor_max = 500, limit = 100, .progress = TRUE)

# data not found
# cr_members(query="adfdf")
# cr_members(member_ids=c(323234343434,3434343434), works=TRUE, facet=TRUE)
# cr_members(member_ids=c(323234343434,3434343434,98), works=TRUE,facet=TRUE)

# Low level function - does no parsing to data.frame, get json or a list
cr_members_(query = 'hindawi')
cr_members_(member_ids = 98)
cr_members_(query = 'hindawi', parse=TRUE)
cr_members_(member_ids = 98, works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
cr_members_(member_ids = 98, works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100, parse=TRUE)

# field queries
## query.container-title
cr_members(98, works = TRUE, flq = c(`query.container-title` = 'Ecology'))


# select only certain fields to return
res <- cr_members(98, works = TRUE, select = c('DOI', 'title'))
names(res$data)
}
}
\references{
https://github.com/CrossRef/rest-api-doc
}
\seealso{
Other crossref: 
\code{\link{cr_funders}()},
\code{\link{cr_journals}()},
\code{\link{cr_licenses}()},
\code{\link{cr_prefixes}()},
\code{\link{cr_types}()},
\code{\link{cr_works}()}
}
\concept{crossref}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_ft_text.R
\name{cr_ft_pdf}
\alias{cr_ft_pdf}
\title{Get full text pdf from a DOI}
\usage{
cr_ft_pdf(...)
}
\description{
Get full text pdf from a DOI
}
\note{
see \code{crminer::crm_pdf}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_ft_text.R
\name{cr_ft_text}
\alias{cr_ft_text}
\title{Get full text from a DOI}
\usage{
cr_ft_text(...)
}
\description{
Get full text from a DOI
}
\note{
see \code{crminer::crm_text}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_fundref.r
\name{cr_funders}
\alias{cr_funders}
\alias{cr_funders_}
\title{Search the CrossRef Fundref API}
\usage{
cr_funders(
  dois = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  flq = NULL,
  select = NULL,
  ...
)

cr_funders_(
  dois = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  parse = FALSE,
  flq = NULL,
  select = NULL,
  ...
)
}
\arguments{
\item{dois}{Search by a single DOI or many DOIs.}

\item{query}{Query terms}

\item{filter}{Filter options. See examples for usage examples
and \code{\link{filters}} for what filters are available.
\code{filter} is available for use with \code{cr_works} and
other \code{crossref} family functions with \code{works=TRUE}}

\item{offset}{Number of record to start at. Minimum: 1. For
\code{\link{cr_works}}, and any function setting \code{works = TRUE},
the maximum offset value is 10000. For larger requests use \code{cursor}.}

\item{limit}{Number of results to return in the query. Not relavant when
searching with specific dois. Default: 20. Max: 1000}

\item{sample}{(integer) Number of random results to return. when you use
the sample parameter, the rows and offset parameters are ignored.
Ignored unless \code{works=TRUE}. Max: 100}

\item{sort}{Field to sort on. Acceptable set of fields to sort on:
\itemize{
\item \code{score} OR \code{relevance} - Sort by relevance score
\item \code{updated} - Sort by date of most recent change to metadata.
Currently the same as deposited.
\item \code{deposited} - Sort by time of most recent deposit
\item \code{indexed} - Sort by time of most recent index
\item \code{published} - Sort by publication date
\item \code{published-print} - Sort by print publication date
\item \code{published-online} - Sort by online publication date
\item \code{issued} - Sort by issued date (earliest known publication date)
\item \code{is-referenced-by-count} - Sort by number of times this DOI is
referenced by other Crossref DOIs
\item \code{references-count} - Sort by number of references included in
the references section of the document identified by this DOI
}}

\item{order}{(character) Sort order, one of 'asc' or 'desc'}

\item{facet}{(logical) Include facet results. Boolean or string with
field to facet on. Valid fields are *, affiliation, funder-name,
funder-doi, orcid, container-title, assertion, archive, update-type,
issn, published, source, type-name, publisher-name, license,
category-name, assertion-group. Default: \code{FALSE}}

\item{works}{(logical) If \code{TRUE}, works returned as well, if not then not.}

\item{cursor}{(character) Cursor character string to do deep paging.
Default is None. Pass in '*' to start deep paging. Any combination of
query, filters and facets may be used with deep paging cursors.
While the \code{limit} parameter may be specified along with cursor,
offset and sample cannot be used with the cursor. See
https://github.com/CrossRef/rest-api-doc#deep-paging-with-cursors}

\item{cursor_max}{(integer) Max records to retrieve. Only used when
cursor param used. Because deep paging can result in continuous requests
until all are retrieved, use this parameter to set a maximum number of
records. Of course, if there are less records found than this value,
you will get only those found. When cursor pagination is being used
the \code{limit} parameter sets the chunk size per request.}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{flq}{field queries. One or more field queries. Acceptable set of
field query parameters are:
\itemize{
\item \code{query.container-title}    - Query container-title aka.
publication name
\item \code{query.author} - Query author first and given names
\item \code{query.editor} - Query editor first and given names
\item \code{query.chair}    - Query chair first and given names
\item \code{query.translator} - Query translator first and given names
\item \code{query.contributor} - Query author, editor, chair and
translator first and given names
\item \code{query.bibliographic} - Query bibliographic information, useful
for citation lookup. Includes titles, authors, ISSNs and publication years
\item \code{query.affiliation} - Query contributor affiliations
}

Note: \code{query.title} has been removed - use \code{query.bibliographic}
as a replacement}

\item{select}{(character) One or more field to return (only those fields
are returned)}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}

\item{parse}{(logical) Whether to output json \code{FALSE} or parse to
list \code{TRUE}. Default: \code{FALSE}}
}
\description{
Search the CrossRef Fundref API
}
\details{
BEWARE: The API will only work for CrossRef DOIs.

This function name changing to \code{cr_funders} in the next version -
both work for now
}
\note{
See the "Rate limiting" seciton in \link{rcrossref} to get
into the "fast lane"
}
\section{Deep paging (using the cursor)}{

When using the cursor, a character string called \code{next-cursor} is
returned from Crossref that we use to do the next request, and so on. We
use a while loop to get number of results up to the value of
\code{cursor_max}. Since we are doing each request for you, you may not
need the \code{next-cursor} string, but if you do want it, you can get
to it by indexing into the result like \code{x$meta$next_cursor}

Note that you can pass in curl options when using cursor, via \code{"..."}
}

\section{NOTE}{

Funders without IDs don't show up on the /funders route, and in this
function. Some funders don't have assigned IDs in Crossref's system,
so won't show up in searches.
}

\examples{
\dontrun{
cr_funders(query="NSF", limit=1)
cr_funders(query="NSF")
cr_funders(dois='10.13039/100000001')
out <- cr_funders(dois=c('10.13039/100000001','10.13039/100000015'))
out['10.13039/100000001']
out[['10.13039/100000001']]

cr_funders(dois='10.13039/100000001')
cr_funders(dois='10.13039/100000001', works=TRUE, limit=5)

cr_funders(dois=c('10.13039/100000001','10.13039/100000015'))
cr_funders(dois=c('10.13039/100000001','10.13039/100000015'), works=TRUE)

## get facets back
cr_funders("10.13039/100000001", works=TRUE, facet=TRUE, limit = 0)
cr_funders("10.13039/100000001", works=TRUE, facet="license:*", limit = 0)
cr_funders('100000001', works = TRUE, cursor = "*", cursor_max = 500,
   limit = 100, facet=TRUE)

# Curl options
cr_funders(dois='10.13039/100000001', verbose = TRUE)

# If not found, and > 1 DOI given, those not found dropped
cr_funders(dois=c("adfadfaf","asfasf"))
cr_funders(dois=c("adfadfaf","asfasf"), works=TRUE)
cr_funders(dois=c("10.13039/100000001","asfasf"))
cr_funders(dois=c("10.13039/100000001","asfasf"), works=TRUE)

# Use the cursor for deep paging
cr_funders('100000001', works = TRUE, cursor = "*", cursor_max = 500,
   limit = 100)
cr_funders(c('100000001', '100000002'), works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
## with optional progress bar
cr_funders('100000001', works = TRUE, cursor = "*", cursor_max = 500,
   limit = 100, .progress = TRUE)

# Low level function - does no parsing to data.frame, get json or a list
cr_funders_(query = 'nsf')
cr_funders_('10.13039/100000001')
cr_funders_(query = 'science', parse=TRUE)
cr_funders_('10.13039/100000001', works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
cr_funders_('10.13039/100000001', works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100, parse = TRUE)

# field queries
## query.container-title
cr_funders('10.13039/100000001', works = TRUE,
  flq = c(`query.container-title` = 'Ecology'))

# select only certain fields to return
res <- cr_funders('10.13039/100000001', works = TRUE, 
  select = c('DOI', 'title'))
names(res$data)
}
}
\references{
https://github.com/CrossRef/rest-api-doc
}
\seealso{
Other crossref: 
\code{\link{cr_journals}()},
\code{\link{cr_licenses}()},
\code{\link{cr_members}()},
\code{\link{cr_prefixes}()},
\code{\link{cr_types}()},
\code{\link{cr_works}()}
}
\concept{crossref}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_search_free.r
\name{cr_search_free}
\alias{cr_search_free}
\title{Search the CrossRef Metadata for DOIs using free form references.}
\usage{
cr_search_free(...)
}
\description{
Search the CrossRef Metadata for DOIs using free form references.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_citation_count.r
\name{cr_citation_count}
\alias{cr_citation_count}
\title{Get a citation count via CrossRef OpenURL}
\usage{
cr_citation_count(
  doi,
  url = "http://www.crossref.org/openurl/",
  key = "cboettig@ropensci.org",
  async = FALSE,
  ...
)
}
\arguments{
\item{doi}{(character) One or more digital object identifiers. If
\code{async=FALSE} we do synchronous HTTP requests in an \code{lapply} call, but
if \code{async=TRUE}, we do asynchronous HTTP requests.}

\item{url}{(character) the url for the function (should be left to default)}

\item{key}{your Crossref OpenURL email address, either enter, or loads
from \code{.Rprofile}. We use a default, so you don't need to pass this.}

\item{async}{(logical) use async HTTP requests. Default: \code{FALSE}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
a data.frame, with columns \code{doi} and \code{count}. The count column
has numeric values that are the citation count for that DOI, or \code{NA} if
not found or no count available
}
\description{
Get a citation count via CrossRef OpenURL
}
\details{
See https://www.crossref.org/labs/openurl/ for more info on this
Crossref API service.

This number is also known as \strong{cited-by}

Note that this number may be out of sync/may not match that that the
publisher is showing (if they show it) for the same DOI/article.

We've contacted Crossref about this, and they have confirmed this.
Unfortunately, we can not do anything about this.

I would imagine it's best to use this data instead of from the publishers,
and this data you can get programatically :)
}
\section{failure behavior}{

When a DOI does not exist, we may not get a proper HTTP status code
to throw a proper stop status, so we grep on the text itself, and throw
a stop if only one DOI passed and not using async, or warning if more
than one DOI passed or if using async.
}

\examples{
\dontrun{
cr_citation_count(doi="10.1371/journal.pone.0042793")
cr_citation_count(doi="10.1016/j.fbr.2012.01.001")
## many
dois <- c("10.1016/j.fbr.2012.01.001", "10.1371/journal.pone.0042793")
cr_citation_count(doi = dois)
# DOI not found
cr_citation_count(doi="10.1016/j.fbr.2012")

# asyc
dois <- c("10.1016/j.fbr.2012.01.001", "10.1371/journal.pone.0042793", 
 "10.1016/j.fbr.2012", "10.1109/tsp.2006.874779", "10.1007/bf02231542", 
 "10.1007/s00277-016-2782-z", "10.1002/9781118339893.wbeccp020", 
 "10.1177/011542659200700105", "10.1002/chin.197444438", 
 "10.1002/9781118619599.ch4", "10.1007/s00466-012-0724-8", 
 "10.1017/s0376892900029477", "10.1167/16.12.824")
res <- cr_citation_count(doi = dois, async = TRUE)
## verbose curl
res <- cr_citation_count(doi = dois, async = TRUE, verbose = TRUE)
res
## time comparison
system.time(cr_citation_count(doi = dois, async = TRUE))
system.time(cr_citation_count(doi = dois, async = FALSE))

# from a set of random DOIs
cr_citation_count(cr_r(50), async = TRUE)
}
}
\seealso{
\code{\link[=cr_search]{cr_search()}}, \code{\link[=cr_r]{cr_r()}}
}
\author{
Carl Boettiger \email{cboettig@gmail.com},
Scott Chamberlain
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_prefixes.r
\name{cr_prefixes}
\alias{cr_prefixes}
\alias{cr_prefixes_}
\title{Search CrossRef prefixes}
\usage{
cr_prefixes(
  prefixes,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  flq = NULL,
  select = NULL,
  ...
)

cr_prefixes_(
  prefixes,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  parse = FALSE,
  flq = NULL,
  select = NULL,
  ...
)
}
\arguments{
\item{prefixes}{(character) Publisher prefixes, one or more in a vector or
list. Required.}

\item{query}{Query terms}

\item{filter}{Filter options. See examples for usage examples
and \code{\link{filters}} for what filters are available.
\code{filter} is available for use with \code{cr_works} and
other \code{crossref} family functions with \code{works=TRUE}}

\item{offset}{Number of record to start at. Minimum: 1. For
\code{\link{cr_works}}, and any function setting \code{works = TRUE},
the maximum offset value is 10000. For larger requests use \code{cursor}.}

\item{limit}{Number of results to return in the query. Not relavant when
searching with specific dois. Default: 20. Max: 1000}

\item{sample}{(integer) Number of random results to return. when you use
the sample parameter, the rows and offset parameters are ignored.
Ignored unless \code{works=TRUE}. Max: 100}

\item{sort}{Field to sort on. Acceptable set of fields to sort on:
\itemize{
\item \code{score} OR \code{relevance} - Sort by relevance score
\item \code{updated} - Sort by date of most recent change to metadata.
Currently the same as deposited.
\item \code{deposited} - Sort by time of most recent deposit
\item \code{indexed} - Sort by time of most recent index
\item \code{published} - Sort by publication date
\item \code{published-print} - Sort by print publication date
\item \code{published-online} - Sort by online publication date
\item \code{issued} - Sort by issued date (earliest known publication date)
\item \code{is-referenced-by-count} - Sort by number of times this DOI is
referenced by other Crossref DOIs
\item \code{references-count} - Sort by number of references included in
the references section of the document identified by this DOI
}}

\item{order}{(character) Sort order, one of 'asc' or 'desc'}

\item{facet}{(logical) Include facet results. Boolean or string with
field to facet on. Valid fields are *, affiliation, funder-name,
funder-doi, orcid, container-title, assertion, archive, update-type,
issn, published, source, type-name, publisher-name, license,
category-name, assertion-group. Default: \code{FALSE}}

\item{works}{(logical) If \code{TRUE}, works returned as well, if not then not.}

\item{cursor}{(character) Cursor character string to do deep paging.
Default is None. Pass in '*' to start deep paging. Any combination of
query, filters and facets may be used with deep paging cursors.
While the \code{limit} parameter may be specified along with cursor,
offset and sample cannot be used with the cursor. See
https://github.com/CrossRef/rest-api-doc#deep-paging-with-cursors}

\item{cursor_max}{(integer) Max records to retrieve. Only used when
cursor param used. Because deep paging can result in continuous requests
until all are retrieved, use this parameter to set a maximum number of
records. Of course, if there are less records found than this value,
you will get only those found. When cursor pagination is being used
the \code{limit} parameter sets the chunk size per request.}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{flq}{field queries. One or more field queries. Acceptable set of
field query parameters are:
\itemize{
\item \code{query.container-title}    - Query container-title aka.
publication name
\item \code{query.author} - Query author first and given names
\item \code{query.editor} - Query editor first and given names
\item \code{query.chair}    - Query chair first and given names
\item \code{query.translator} - Query translator first and given names
\item \code{query.contributor} - Query author, editor, chair and
translator first and given names
\item \code{query.bibliographic} - Query bibliographic information, useful
for citation lookup. Includes titles, authors, ISSNs and publication years
\item \code{query.affiliation} - Query contributor affiliations
}

Note: \code{query.title} has been removed - use \code{query.bibliographic}
as a replacement}

\item{select}{(character) One or more field to return (only those fields
are returned)}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}

\item{parse}{(logical) Whether to output json \code{FALSE} or parse to
list \code{TRUE}. Default: \code{FALSE}}
}
\description{
Search CrossRef prefixes
}
\details{
BEWARE: The API will only work for CrossRef DOIs.

Note that any one publisher can have more than one DOI. If you want to
search on all DOIs for a publisher, pass in all DOIs, or see
\code{\link[=cr_members]{cr_members()}}, and pass in the \code{member_ids} parameter.

Notes from CrossRef (quoting them):

The prefix of a CrossRef DOI does NOT indicate who currently owns the DOI.
It only reflects who originally registered the DOI. CrossRef metadata has
an \code{owner_prefix} element that records the current owner of the
CrossRef DOI in question.

CrossRef also has member IDs for depositing organisations. A single member
may control multiple owner prefixes, which in turn may control a number of
DOIs. When looking at works published by a certain organisation, member
IDs and the member routes should be used.
}
\note{
See the "Rate limiting" seciton in \link{rcrossref} to get
into the "fast lane"
}
\section{Deep paging (using the cursor)}{

When using the cursor, a character string called \code{next-cursor} is
returned from Crossref that we use to do the next request, and so on. We
use a while loop to get number of results up to the value of
\code{cursor_max}. Since we are doing each request for you, you may not
need the \code{next-cursor} string, but if you do want it, you can get
to it by indexing into the result like \code{x$meta$next_cursor}

Note that you can pass in curl options when using cursor, via \code{"..."}
}

\examples{
\dontrun{
cr_prefixes(prefixes="10.1016")
cr_prefixes(prefixes="10.1016", works=TRUE)
cr_prefixes(prefixes=c('10.1016','10.1371','10.1023','10.4176','10.1093'))
cr_prefixes(prefixes=c('10.1016','10.1371'), works=TRUE)
cr_prefixes(prefixes="10.1016", works=TRUE, filter=c(has_full_text=TRUE), 
  limit=5)
cr_prefixes(prefixes="10.1016", works=TRUE, query='ecology', limit=4)
cr_prefixes(prefixes="10.1016", works=TRUE, query='ecology', limit=4)

# facets - only avail. when works=TRUE
cr_prefixes(prefixes="10.1016", works=TRUE, facet=TRUE)
cr_prefixes(prefixes="10.1016", works=TRUE, facet="license:*", limit=0)
cr_prefixes(prefixes=c('10.1016','10.1371'), works=TRUE, facet=TRUE,
  limit=0)

# Use the cursor for deep paging
cr_prefixes("10.1016", works = TRUE, cursor = "*", cursor_max = 500,
   limit = 100)
cr_prefixes(c('10.1016', '10.1371'), works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
## with optional progress bar
cr_prefixes("10.1016", works = TRUE, cursor = "*", cursor_max = 500,
   limit = 100, .progress = TRUE)

# Low level function - does no parsing to data.frame, get json or a list
cr_prefixes_("10.1016")
cr_prefixes_(c('10.1016', '10.1371'))
cr_prefixes_("10.1016", works = TRUE, query = 'ecology', limit = 10)
cr_prefixes_("10.1016", works = TRUE, query = 'ecology', parse=TRUE,
   limit = 10)
cr_prefixes_("10.1016", works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
cr_prefixes_("10.1016", works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100, parse = TRUE)

# field queries
## query.container-title
cr_prefixes("10.1016", works = TRUE,
  flq = c(`query.container-title` = 'Ecology'))

# select only certain fields to return
res <- cr_prefixes("10.1016", works = TRUE, select = c('DOI', 'title'))
names(res$data)
}
}
\references{
https://github.com/CrossRef/rest-api-doc
}
\seealso{
Other crossref: 
\code{\link{cr_funders}()},
\code{\link{cr_journals}()},
\code{\link{cr_licenses}()},
\code{\link{cr_members}()},
\code{\link{cr_types}()},
\code{\link{cr_works}()}
}
\concept{crossref}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_works.R
\name{cr_works}
\alias{cr_works}
\alias{cr_works_}
\title{Search CrossRef works (articles)}
\usage{
cr_works(
  dois = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  flq = NULL,
  select = NULL,
  async = FALSE,
  ...
)

cr_works_(
  dois = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  parse = FALSE,
  flq = NULL,
  select = NULL,
  async = FALSE,
  ...
)
}
\arguments{
\item{dois}{Search by a single DOI or many DOIs.  Note that using this
parameter at the same time as the \code{query}, \code{limit}, \code{select} or \code{flq}
parameter will result in an error.}

\item{query}{Query terms}

\item{filter}{Filter options. See examples for usage examples
and \code{\link{filters}} for what filters are available.
\code{filter} is available for use with \code{cr_works} and
other \code{crossref} family functions with \code{works=TRUE}}

\item{offset}{Number of record to start at. Minimum: 1. For
\code{\link{cr_works}}, and any function setting \code{works = TRUE},
the maximum offset value is 10000. For larger requests use \code{cursor}.}

\item{limit}{Number of results to return in the query. Not relavant when
searching with specific dois. Default: 20. Max: 1000}

\item{sample}{(integer) Number of random results to return. when you use
the sample parameter, the rows and offset parameters are ignored.
Ignored unless \code{works=TRUE}. Max: 100}

\item{sort}{Field to sort on. Acceptable set of fields to sort on:
\itemize{
\item \code{score} OR \code{relevance} - Sort by relevance score
\item \code{updated} - Sort by date of most recent change to metadata.
Currently the same as deposited.
\item \code{deposited} - Sort by time of most recent deposit
\item \code{indexed} - Sort by time of most recent index
\item \code{published} - Sort by publication date
\item \code{published-print} - Sort by print publication date
\item \code{published-online} - Sort by online publication date
\item \code{issued} - Sort by issued date (earliest known publication date)
\item \code{is-referenced-by-count} - Sort by number of times this DOI is
referenced by other Crossref DOIs
\item \code{references-count} - Sort by number of references included in
the references section of the document identified by this DOI
}}

\item{order}{(character) Sort order, one of 'asc' or 'desc'}

\item{facet}{(logical) Include facet results. Boolean or string with
field to facet on. Valid fields are *, affiliation, funder-name,
funder-doi, orcid, container-title, assertion, archive, update-type,
issn, published, source, type-name, publisher-name, license,
category-name, assertion-group. Default: \code{FALSE}}

\item{cursor}{(character) Cursor character string to do deep paging.
Default is None. Pass in '*' to start deep paging. Any combination of
query, filters and facets may be used with deep paging cursors.
While the \code{limit} parameter may be specified along with cursor,
offset and sample cannot be used with the cursor. See
https://github.com/CrossRef/rest-api-doc#deep-paging-with-cursors}

\item{cursor_max}{(integer) Max records to retrieve. Only used when
cursor param used. Because deep paging can result in continuous requests
until all are retrieved, use this parameter to set a maximum number of
records. Of course, if there are less records found than this value,
you will get only those found. When cursor pagination is being used
the \code{limit} parameter sets the chunk size per request.}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{flq}{field queries. One or more field queries. Acceptable set of
field query parameters are:
\itemize{
\item \code{query.container-title}    - Query container-title aka.
publication name
\item \code{query.author} - Query author first and given names
\item \code{query.editor} - Query editor first and given names
\item \code{query.chair}    - Query chair first and given names
\item \code{query.translator} - Query translator first and given names
\item \code{query.contributor} - Query author, editor, chair and
translator first and given names
\item \code{query.bibliographic} - Query bibliographic information, useful
for citation lookup. Includes titles, authors, ISSNs and publication years
\item \code{query.affiliation} - Query contributor affiliations
}

Note: \code{query.title} has been removed - use \code{query.bibliographic}
as a replacement}

\item{select}{(character) One or more field to return (only those fields
are returned)}

\item{async}{(logical) use async HTTP requests. Default: \code{FALSE}}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}

\item{parse}{(logical) Whether to output json \code{FALSE} or parse to
list \code{TRUE}. Default: \code{FALSE}}
}
\description{
Search CrossRef works (articles)
}
\note{
See the "Rate limiting" seciton in \link{rcrossref} to get
into the "fast lane"
}
\section{Deep paging (using the cursor)}{

When using the cursor, a character string called \code{next-cursor} is
returned from Crossref that we use to do the next request, and so on. We
use a while loop to get number of results up to the value of
\code{cursor_max}. Since we are doing each request for you, you may not
need the \code{next-cursor} string, but if you do want it, you can get
to it by indexing into the result like \code{x$meta$next_cursor}

Note that you can pass in curl options when using cursor, via \code{"..."}
}

\section{Beware}{

The API will only work for CrossRef DOIs.
}

\section{Functions}{

\itemize{
\item \code{cr_works()} - Does data request and parses to data.frame for
easy downstream consumption
\item \code{cr_works_()} - Does data request, and gives back json (default)
or lists, with no attempt to parse to data.frame's
}
}

\section{Explanation of some data fields}{

\itemize{
\item score: a term frequency, inverse document frequency score that
comes from the Crossref Solr backend, based on bibliographic metadata
fields title, publication title, authors, ISSN, publisher, and
date of publication.
}
}

\examples{
\dontrun{
# Works funded by the NSF
cr_works(query="NSF")

# Works that include renear but not ontologies
cr_works(query="renear+-ontologies")

# Filter
cr_works(query="global state", filter=c(has_orcid=TRUE), limit=3)
# Filter by multiple fields
cr_works(filter=c(has_orcid=TRUE, from_pub_date='2004-04-04'))
# Only full text articles
cr_works(filter=c(has_full_text = TRUE))
# has affilitation data
cr_works(filter=c(has_affiliation = TRUE))
# has abstract
cr_works(filter=c(has_abstract = TRUE))
# has clinical trial number
cr_works(filter=c(has_clinical_trial_number = TRUE))

# Querying dois
cr_works(dois='10.1063/1.3593378')
cr_works('10.1371/journal.pone.0033693')
cr_works(dois='10.1007/12080.1874-1746')
cr_works(dois=c('10.1007/12080.1874-1746','10.1007/10452.1573-5125',
   '10.1111/(issn)1442-9993'))

# progress bar
cr_works(dois=c('10.1007/12080.1874-1746','10.1007/10452.1573-5125'),
   .progress="text")

# Include facetting in results
cr_works(query="NSF", facet=TRUE)
## Get facets only, by setting limit=0
cr_works(query="NSF", facet=TRUE, limit=0)
## you can also set facet to a query
cr_works(facet = "license:*", limit=0)

# Sort results
cr_works(query="ecology", sort='relevance', order="asc")
res <- cr_works(query="ecology", sort='score', order="asc")
res$data$score
cr_works(query="ecology", sort='published')
x=cr_works(query="ecology", sort='published-print')
x=cr_works(query="ecology", sort='published-online')

# Get a random number of results
cr_works(sample=1)
cr_works(sample=10)

# You can pass in dot separated fields to filter on specific fields
cr_works(filter=c(award.number='CBET-0756451',
   award.funder='10.13039/100000001'))

# Use the cursor for deep paging
cr_works(query="NSF", cursor = "*", cursor_max = 300, limit = 100)
cr_works(query="NSF", cursor = "*", cursor_max = 300, limit = 100,
   facet = TRUE)
## with optional progress bar
x <- cr_works(query="NSF", cursor = "*", cursor_max = 1200, limit = 200, 
  .progress = TRUE)

# Low level function - does no parsing to data.frame, get json or a list
cr_works_(query = "NSF")
cr_works_(query = "NSF", parse=TRUE)
cr_works_(query="NSF", cursor = "*", cursor_max = 300, limit = 100)
cr_works_(query="NSF", cursor = "*", cursor_max = 300, limit = 100,
   parse=TRUE)

# field queries
## query.author
res <- cr_works(query = "ecology", flq = c(query.author = 'Boettiger'))

## query.container-title
res <- cr_works(query = "ecology",
  flq = c(`query.container-title` = 'Ecology'))

## query.author and query.bibliographic
res <- cr_works(query = "ecology",
  flq = c(query.author = 'Smith', query.bibliographic = 'cell'))

# select only certain fields to return
res <- cr_works(query = "NSF", select = c('DOI', 'title'))
names(res$data)

# asyc
queries <- c("ecology", "science", "cellular", "birds", "European",
  "bears", "beets", "laughter", "hapiness", "funding")
res <- cr_works(query = queries, async = TRUE)
res_json <- cr_works_(query = queries, async = TRUE)
unname(vapply(res_json, class, ""))
jsonlite::fromJSON(res_json[[1]])

queries <- c("ecology", "science", "cellular")
res <- cr_works(query = queries, async = TRUE, verbose = TRUE)
res

# time
queries <- c("ecology", "science", "cellular", "birds", "European",
  "bears", "beets", "laughter", "hapiness", "funding")
system.time(cr_works(query = queries, async = TRUE))
system.time(lapply(queries, function(z) cr_works(query = z)))
}
}
\references{
https://github.com/CrossRef/rest-api-doc
}
\seealso{
Other crossref: 
\code{\link{cr_funders}()},
\code{\link{cr_journals}()},
\code{\link{cr_licenses}()},
\code{\link{cr_members}()},
\code{\link{cr_prefixes}()},
\code{\link{cr_types}()}
}
\concept{crossref}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_abstract.R
\name{cr_abstract}
\alias{cr_abstract}
\title{Get abstract}
\usage{
cr_abstract(doi, ...)
}
\arguments{
\item{doi}{(character) a DOI, required.}

\item{...}{Named parameters passed on to \code{\link[crul]{HttpClient}}}
}
\description{
Get abstract
}
\examples{
\dontrun{
# abstract found
cr_abstract('10.1109/TASC.2010.2088091')
cr_abstract("10.1175//2572.1")
cr_abstract("10.1182/blood.v16.1.1039.1039")

# doi not found
# cr_abstract(doi = '10.5284/1011335')

# abstract not found, throws error
# cr_abstract(doi = '10.1126/science.169.3946.635')

# a random DOI
# cr_abstract(cr_r(1))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdmurl.R
\name{as.tdmurl}
\alias{as.tdmurl}
\title{Coerce a url to a tdmurl with a specific type}
\usage{
as.tdmurl(...)
}
\description{
Coerce a url to a tdmurl with a specific type
}
\note{
see \code{crminer::as_tdmurl}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_r.r
\name{cr_r}
\alias{cr_r}
\title{Get a random set of DOI's through CrossRef.}
\usage{
cr_r(sample = 10, ...)
}
\arguments{
\item{sample}{The number of returned random DOIs. Maximum: 100. Default: 20.}

\item{...}{Further args passed on to \code{\link[=cr_works]{cr_works()}}}
}
\value{
A character vector of DOIs
}
\description{
Get a random set of DOI's through CrossRef.
}
\examples{
\dontrun{
# Default search gets 10 random DOIs
cr_r()

# Get 30 DOIs
cr_r(30)
}
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_journals.r
\name{cr_journals}
\alias{cr_journals}
\alias{cr_journals_}
\title{Search CrossRef journals}
\usage{
cr_journals(
  issn = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  flq = NULL,
  select = NULL,
  ...
)

cr_journals_(
  issn = NULL,
  query = NULL,
  filter = NULL,
  offset = NULL,
  limit = NULL,
  sample = NULL,
  sort = NULL,
  order = NULL,
  facet = FALSE,
  works = FALSE,
  cursor = NULL,
  cursor_max = 5000,
  .progress = "none",
  parse = FALSE,
  flq = NULL,
  select = NULL,
  ...
)
}
\arguments{
\item{issn}{(character) One or more ISSN's. Format: XXXX-XXXX.}

\item{query}{Query terms}

\item{filter}{Filter options. See examples for usage examples
and \code{\link{filters}} for what filters are available.
\code{filter} is available for use with \code{cr_works} and
other \code{crossref} family functions with \code{works=TRUE}}

\item{offset}{Number of record to start at. Minimum: 1. For
\code{\link{cr_works}}, and any function setting \code{works = TRUE},
the maximum offset value is 10000. For larger requests use \code{cursor}.}

\item{limit}{Number of results to return in the query. Not relavant when
searching with specific dois. Default: 20. Max: 1000}

\item{sample}{(integer) Number of random results to return. when you use
the sample parameter, the rows and offset parameters are ignored.
Ignored unless \code{works=TRUE}. Max: 100}

\item{sort}{Field to sort on. Acceptable set of fields to sort on:
\itemize{
\item \code{score} OR \code{relevance} - Sort by relevance score
\item \code{updated} - Sort by date of most recent change to metadata.
Currently the same as deposited.
\item \code{deposited} - Sort by time of most recent deposit
\item \code{indexed} - Sort by time of most recent index
\item \code{published} - Sort by publication date
\item \code{published-print} - Sort by print publication date
\item \code{published-online} - Sort by online publication date
\item \code{issued} - Sort by issued date (earliest known publication date)
\item \code{is-referenced-by-count} - Sort by number of times this DOI is
referenced by other Crossref DOIs
\item \code{references-count} - Sort by number of references included in
the references section of the document identified by this DOI
}}

\item{order}{(character) Sort order, one of 'asc' or 'desc'}

\item{facet}{(logical) Include facet results. Boolean or string with
field to facet on. Valid fields are *, affiliation, funder-name,
funder-doi, orcid, container-title, assertion, archive, update-type,
issn, published, source, type-name, publisher-name, license,
category-name, assertion-group. Default: \code{FALSE}}

\item{works}{(logical) If \code{TRUE}, works returned as well, if not then not.}

\item{cursor}{(character) Cursor character string to do deep paging.
Default is None. Pass in '*' to start deep paging. Any combination of
query, filters and facets may be used with deep paging cursors.
While the \code{limit} parameter may be specified along with cursor,
offset and sample cannot be used with the cursor. See
https://github.com/CrossRef/rest-api-doc#deep-paging-with-cursors}

\item{cursor_max}{(integer) Max records to retrieve. Only used when
cursor param used. Because deep paging can result in continuous requests
until all are retrieved, use this parameter to set a maximum number of
records. Of course, if there are less records found than this value,
you will get only those found. When cursor pagination is being used
the \code{limit} parameter sets the chunk size per request.}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{flq}{field queries. One or more field queries. Acceptable set of
field query parameters are:
\itemize{
\item \code{query.container-title}    - Query container-title aka.
publication name
\item \code{query.author} - Query author first and given names
\item \code{query.editor} - Query editor first and given names
\item \code{query.chair}    - Query chair first and given names
\item \code{query.translator} - Query translator first and given names
\item \code{query.contributor} - Query author, editor, chair and
translator first and given names
\item \code{query.bibliographic} - Query bibliographic information, useful
for citation lookup. Includes titles, authors, ISSNs and publication years
\item \code{query.affiliation} - Query contributor affiliations
}

Note: \code{query.title} has been removed - use \code{query.bibliographic}
as a replacement}

\item{select}{(character) One or more field to return (only those fields
are returned)}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}

\item{parse}{(logical) Whether to output json \code{FALSE} or parse to
list \code{TRUE}. Default: \code{FALSE}}
}
\description{
Search CrossRef journals
}
\details{
BEWARE: The API will only work for CrossRef DOIs.

Note that some parameters are ignored unless \code{works=TRUE}: sample, sort,
order, filter
}
\note{
See the "Rate limiting" seciton in \link{rcrossref} to get
into the "fast lane"
}
\section{Deep paging (using the cursor)}{

When using the cursor, a character string called \code{next-cursor} is
returned from Crossref that we use to do the next request, and so on. We
use a while loop to get number of results up to the value of
\code{cursor_max}. Since we are doing each request for you, you may not
need the \code{next-cursor} string, but if you do want it, you can get
to it by indexing into the result like \code{x$meta$next_cursor}

Note that you can pass in curl options when using cursor, via \code{"..."}
}

\section{Explanation of some data fields}{

\itemize{
\item backfile_dois: Back file records have a publication date older than
two years ago.
\item current_dois: Current records are anything published in the last
two years.
}
}

\examples{
\dontrun{
cr_journals(issn="2167-8359")
cr_journals()
cr_journals(issn="2167-8359", works=TRUE)
cr_journals(issn=c('1803-2427','2326-4225'))
cr_journals(query="ecology")
cr_journals(issn="2167-8359", query='ecology', works=TRUE,
   sort='score', order="asc")
cr_journals(issn="2167-8359", query='ecology', works=TRUE, sort='score',
   order="desc")
cr_journals(issn="2167-8359", works=TRUE,
   filter=c(from_pub_date='2014-03-03'))
cr_journals(query="peerj")
cr_journals(issn='1803-2427', works=TRUE)
cr_journals(issn='1803-2427', works=TRUE, sample=1)
cr_journals(limit=2)

## get facets back
cr_journals('1803-2427', works=TRUE, facet=TRUE)
cr_journals('1803-2427', works=TRUE, facet="published:*", limit = 0)
cr_journals(issn=c('1803-2427','2326-4225'), works=TRUE,
  facet="published:*", limit = 10)

# Use the cursor for deep paging
cr_journals(issn='1932-6203', works = TRUE, cursor = "*", cursor_max = 500,
   limit = 100)
cr_journals(c('1932-6203', '0028-0836'), works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
## with optional progress bar
cr_journals(issn='1932-6203', works = TRUE, cursor = "*", cursor_max = 90,
   limit = 30, .progress = TRUE)

# fails, if you want works, you must give an ISSN
# cr_journals(query = "ecology", filter=c(has_full_text = TRUE),
#    works = TRUE)

# Low level function - does no parsing to data.frame, get json or a list
cr_journals_(query = 'ecology')
cr_journals_("2167-8359")
cr_journals_(query = 'ecology', parse=TRUE)
cr_journals_("2167-8359", works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100)
cr_journals_("2167-8359", works = TRUE, cursor = "*",
   cursor_max = 300, limit = 100, parse = TRUE)

# field queries
## query.author
cr_journals("2167-8359", works = TRUE, flq = c(`query.author` = 'Jane'))

# select only certain fields to return
res <- cr_journals('2167-8359', works = TRUE, 
  select = c('DOI', 'title'))
names(res$data)
}
}
\references{
https://github.com/CrossRef/rest-api-doc
}
\seealso{
Other crossref: 
\code{\link{cr_funders}()},
\code{\link{cr_licenses}()},
\code{\link{cr_members}()},
\code{\link{cr_prefixes}()},
\code{\link{cr_types}()},
\code{\link{cr_works}()}
}
\concept{crossref}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcrossref-package.R
\docType{package}
\name{rcrossref-package}
\alias{rcrossref-package}
\alias{rcrossref}
\title{rcrossref}
\description{
R Client for Various CrossRef APIs.
}
\section{Crossref APIs}{

rcrossref interacts with the main Crossref metadata search API at
https://github.com/CrossRef/rest-api-doc,
the old metadata search API at http://search.labs.crossref.org/, their
DOI Content Negotiation service at
http://citation.crosscite.org/docs.html, and
the \emph{Text and Data Mining} project http://tdmsupport.crossref.org/

Crossref's API issue tracker lives at https://gitlab.com/crossref/issues
it's a good place to go ask them about things related to their API
that go beyond the R interface here.
}

\section{Defunct}{

See \link{rcrossref-deprecated} and \link{rcrossref-defunct}
for details.
}

\section{What am I actually searching?}{

When you use the \verb{cr_*()} functions in this package, you are using
the Crossref search API described at
https://github.com/CrossRef/rest-api-doc
When you search with query terms, on Crossref servers they are not searching
full text, or even abstracts of articles, but only what is available in the
data that is returned to you. That is, they search article titles, authors,
etc. For some discussion on this, see
https://gitlab.com/crossref/issues/issues/101
}

\section{Rate limiting}{

From time to time Crossref needs to impose rate limits to ensure that
the free API is usable by all. Any rate limits that are in effect will
be advertised in the \code{X-Rate-Limit-Limit} and
\code{X-Rate-Limit-Interval} HTTP headers.

This boils down to: they allow X number of requests per some time period.
The numbers can change so we can't give a rate limit that will always
be in effect. If you're curious pass in \code{verbose = TRUE} to
your function call, and you'll get headers that will display these rate
limits.

\strong{Be nice and share your email with Crossref}

The Crossref team encourage requests with appropriate contact information
and will forward you to a dedicated API cluster for improved performance when
you share your email address with them.
https://github.com/CrossRef/rest-api-doc#good-manners--more-reliable-service

To pass your email address to Crossref via this client, simply store it
as environment variable in \code{.Renviron} like this:
\enumerate{
\item Open file:
\code{file.edit("~/.Renviron")}
\item Add email address to be shared with Crossref
\code{crossref_email = name@example.com}
\item Save the file and restart your R session
}

Don't wanna share your email any longer? Simply delete it from
\verb{~/.Renviron}
}

\section{Text mining}{

All Crossref specific text mining functions are now deprecated, and
moved to a new package \code{crminer}.

Another package \pkg{fulltext} is designed solely to do general purpose
text mining involving Crossref and other sources of scholarly metadata
and full text.
}

\section{High and Low Level APIs}{

For the Crossref search API (the functions \code{\link[=cr_funders]{cr_funders()}},
\code{\link[=cr_journals]{cr_journals()}}, \code{\link[=cr_licenses]{cr_licenses()}},
\code{\link[=cr_members]{cr_members()}}, \code{\link[=cr_prefixes]{cr_prefixes()}}, \code{\link[=cr_types]{cr_types()}},
\code{\link[=cr_works]{cr_works()}}), there is a high level API and a low level. The
high level is accessible through those functions just listed (e.g.,
\code{\link[=cr_works]{cr_works()}}), whereas the low level is accessible via the same
fxn name with an underscore (e.g., \code{\link[=cr_works_]{cr_works_()}}). The high level
API does data requests, and parses to data.frame's. Since the high level
API functions have been around a while, we didn't want to break their
behavior, so the low level API functions are separate, and only do the data
request, giving back json or a list, with no attempt to parse any further. The
low level API functions will be faster because there's much less parsing, and
therefore less prone to potential errors due to changes in the Crossref API
that could cause parsing errors. Note that cursor feature works with both
high and low level.
}

\section{RStudio Addin}{

On installation of \pkg{rcrossref} you get an RStudio Addin. To use the Addin,
go to the top toolbar > Tools > Addins > Add Crossref Citations. You'll get a
window pop up that you can put in DOIs for. If the DOI is found, the bibtex
citations will be added to a file called \code{crossref.bib}. New citations
will be appended to that file. Addin authored by Hao Zhu
https://github.com/haozhu233
}

\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cr_agency.r
\name{cr_agency}
\alias{cr_agency}
\title{Check the DOI minting agency on one or more dois}
\usage{
cr_agency(dois = NULL, .progress = "none", ...)
}
\arguments{
\item{dois}{(character) One or more article or organization dois.}

\item{.progress}{Show a \code{plyr}-style progress bar? Options are "none", "text",
"tk", "win", and "time".  See \code{\link[plyr]{create_progress_bar}} for details
of each. Only used when passing in multiple ids (e.g., multiple DOIs, DOI prefixes,
etc.), or when using the \code{cursor} param. When using the \code{cursor} param,
this argument only accept a boolean, either \code{TRUE} or \code{FALSE}; any
non-boolean is coerced to \code{FALSE}.}

\item{...}{Named parameters passed on to \code{\link[crul]{verb-GET}}}
}
\description{
Check the DOI minting agency on one or more dois
}
\examples{
\dontrun{
cr_agency(dois = '10.13039/100000001')
cr_agency(
  dois = c('10.13039/100000001','10.13039/100000015','10.5284/1011335'))
}
}
\references{
https://github.com/CrossRef/rest-api-doc
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmid2doi.R
\name{pmid2doi}
\alias{pmid2doi}
\alias{doi2pmid}
\title{Get a PMID from a DOI, and vice versa.}
\usage{
pmid2doi(...)

doi2pmid(...)
}
\description{
Get a PMID from a DOI, and vice versa.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crosscite.R
\name{crosscite}
\alias{crosscite}
\title{Crosscite - citation formatter}
\usage{
crosscite(...)
}
\description{
Crosscite - citation formatter
}
\keyword{internal}
