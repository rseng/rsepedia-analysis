

fulltext
========

[![cran checks](https://cranchecks.info/badges/flavor/release/fulltext)](https://cranchecks.info/pkgs/fulltext)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/fulltext/workflows/R-check/badge.svg)](https://github.com/ropensci/fulltext/actions/)
[![codecov](https://codecov.io/gh/ropensci/fulltext/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/fulltext)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/fulltext)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/fulltext)](https://cran.r-project.org/package=fulltext)

__Get full text research articles__

Checkout the [package docs][docs] and the [fulltext manual][ftbook] to get started.

-----

rOpenSci has a number of R packages to get either full text, metadata, or both from various publishers. The goal of `fulltext` is to integrate these packages to create a single interface to many data sources.

`fulltext` makes it easy to do text-mining by supporting the following steps:

* Search for articles - `ft_search`
* Fetch articles - `ft_get`
* Get links for full text articles (xml, pdf) - `ft_links`
* Extract text from articles / convert formats - `ft_extract`
* Collect all texts into a data.frame - `ft_table`

Previously supported use cases, extracted out to other packages:

* Collect bits of articles that you actually need - moved to package `pubchunks`
* Supplementary data from papers has been moved to the `suppdata`


It's easy to go from the outputs of `ft_get` to text-mining packages such as 
`tm` and `quanteda`

Data sources in `fulltext` include:

* Crossref - via the `rcrossref` package
* Public Library of Science (PLOS) - via the `rplos` package
* Biomed Central
* arXiv - via the `aRxiv` package
* bioRxiv - via the `biorxivr` package
* PMC/Pubmed via Entrez - via the `rentrez` package
* Scopus - internal tooling
* Semantic Scholar - internal tooling
* Many more are supported via the above sources (e.g., _Royal Society Open Science_ is
available via Pubmed)
* We __will__ add more, as publishers open up, and as we have time...See the issues

Authentication: A number of publishers require authentication via API key, and some even more
draconian authentication processes involving checking IP addresses. We are working on supporting
all the various authentication things for different publishers, but of course all the OA content
is already easily available. See the **Authentication** section in `?fulltext-package` after 
loading the package.

We'd love your feedback. Let us know what you think in the issue tracker


## Installation

Stable version from CRAN


```r
install.packages("fulltext")
```

Development version from GitHub


```r
remotes::install_github("ropensci/fulltext")
```

Load library


```r
library('fulltext')
```

## Interoperability with other packages downstream

Note: this example not included in vignettes as that would require the two below packages in Suggests here. To see many examples and documentation see the [package docs][docs] and the [fulltext manual][ftbook].


```r
cache_options_set(path = (td <- 'foobar'))
res <- ft_get(c('10.7554/eLife.03032', '10.7554/eLife.32763'), type = "pdf")
library(readtext)
x <- readtext::readtext(file.path(cache_options_get()$path, "*.pdf"))
```


```r
library(quanteda)
quanteda::corpus(x)
```

## Contributors

* Scott Chamberlain
* Will Pearse
* Katrin Leinweber

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/fulltext/issues).
* License: MIT
* Get citation information for `fulltext`: `citation(package = 'fulltext')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.


[docs]: https://docs.ropensci.org/fulltext/
[ftbook]: https://books.ropensci.org/fulltext/
fulltext 2.0
==============

### MINOR IMPROVEMENTS

* `ft_abstract()`: change default `from` value from "plos" to "crossref" (#235)

### BUG FIXES

* fix PLOS internals for `ft_abstract()` (#234)
* fix how Miscrosoft Academic query was structured before being sent off - affecting `ft_search()` (#236)


fulltext 1.7.0
==============

### NEW FEATURES

* The Cross Text and Data Mining program has ended. It ended officially at the end of 2020, but it took Wiley until early February to implement a TDM service. This impacts `ft_get` only, and impacts retrieval of only Wiley and Elsevier works (which is a lot since they're very large publishers). You now have to get a TDM key from both Wiley and Elsevier. See the "Authentication" section in `?fulltext` for details. See also `?"ft_get-warnings"`: there's a new warning if we detect that you have a `CROSSREF_TDM` env var   (#224)
* Added new `ft_get` plugin for Trans Tech Publications (#232)
* Import dependency `crminer` dropped. A few `crminer` functions brought over into this package and renamed. This is due to Crossref TDM program ending (#233)

### MINOR IMPROVEMENTS

* new `ft_search()` docs section "Pagination" describing what parameter names to use for each data source to do pagination

### BUG FIXES

* `ft_search()` fix: some data source specific options (e.g., `plosopts`) were being overridden by the global equivalent options (e.g., in `ft_search(limit = 5, from = c("plos", "entrez"), plosopts = list(limit=10))` the limit sent to PLOS was 5 instead of 10). Fixed now so that data source specfiic options override the matching global option (#229)
* Internal API Fatcat changed its parameter name to limit what fields are returned. Fixed now (#231)


fulltext 1.6.0
==============

### NEW FEATURES

* HTTP requests previously made to the ftdoi.org API are now done locally with data/files cached locally; speeds up `ft_get` usage, especially if you're on a slow internet connection (#225) (#226)

### MINOR IMPROVEMENTS

* all vignettes only available on docs site now: https://docs.ropensci.org/fulltext/ (#227)

### BUG FIXES

* `ft_links()` opts parameters dropped (e.g., `plosopts`); no options can be passed down to each data source (#223)


fulltext 1.5.0
==============

### NEW FEATURES

* `ft_abstract()` gains new data source option: Semantic Scholar. Uses async http requests. Beware of rate limits. See 'Rate limits' in `?fulltext-package` (#218)

### MINOR IMPROVEMENTS

* using vcr for some tests (#212)
* improved documentation on proxy settings for Scopus data source added to the 'Authentication' section of the package level manual file, see `?fulltext-package` (#220) thanks @behrica
* for `ft_get()` added new plugin helper for publisher Company of Biologists; and a fix in the Company of Biologists metadata returned from the ftdoi.org API (#216)
* cleaned out old code for brief data.frame printing, use tibble instead

### BUG FIXES

* `from` parameter in `ft_get()` for specifying what data source to use should have been accepting more than 1 value, but was only accepting 1 value. Now `from` accepts 1 or more values (#222)


fulltext 1.4.0
==============

### NEW FEATURES

* `ft_get()` gains a new publisher plugin: Cambridge. its not available in the `from` parameter, but if you pass in a Cambridge DOI, a plugin now exists to attempt to get the pdf (xml not provided) (#195)

### MINOR IMPROVEMENTS

* Wiley now provides xml (in addition to pdf) for at least some articles. So `ft_get()` for Wiley now supports `type="xml"` in addition to `type="pdf"` (#209)
* reduce duplicated code in `ft_get` plugins for checking `type` parameter (#210)
* changed Wiley plugin for `ft_get()`: using different URLs for attempting to fetch articles, those with base url `api.wiley.com` instead of `onlinelibrary.wiley.com` (#191)

### BUG FIXES

* `ft_get()` wasn't allowing `type="pdf"` for PLOS data source; now allowed in addition to xml; also plos requests for articles now using http instead of http (#205) thanks @clbti for the report
* fix to internal fxn `fat_cat_search_one()` used inside of `ft_get()`: fixed data.frame subsetting error due to sometimes missing columns (#206)
* fix to `ft_links()`: was erroring when publisher not supported yet; now gives back no data (#207)
* `ft_browse()` was failing when `what="macrodocs"` selected; macrodocs.org is now dead; `what` parameter now defunct  (#208)
* fix `ft_table()` - was failing on reading malformed files (e.g., an xml file with pdf content inside)  (#211)


fulltext 1.3.0
==============

### NEW FEATURES

* `ft_get()` gains new data source: ScienceDirect (#196) thanks @knh11545
* performance improvement to `ft_get()`: when using internal fxn `get_unknown`(when publisher is not initially knowsn) we run another internal fxn `get_publisher` which is slow cause can only throw one DOI per request to Crossref API; We now use the FatCat (https://fatcat.wiki/) API to do more efficient publisher lookup (#192) (#201)

### MINOR IMPROVEMENTS

* curl options can be passed into most functions now (#199)
* change to `ft_get()` when using Elsevier: we now by default detect if the downloaded file is abstact only and treat that as you dont have access - you can toggle this within `elsevieropts` - see `?ft_get` for help (#200)
* fix link for fulltext manual, ropensci->ropenscilabs (#204)

### BUG FIXES

* fix to `ft_abstract()`: when using plos, we were passing on identifiers to the wrong Solr parameter passed to PLOS Solr API (#202) thanks @mrweiler !
* fix to `ft_links()` with data from Entrez was failing when NULLs were encountered (#203)


fulltext 1.2.0
==============

### NEW FEATURES

* `ft_get()` gains a `progress` parameter, `TRUE` printing a progress bar and `FALSE` not. By default we do not print a progress bar to be consistent with the behavior of previous versions (#140) (#190)
* `cache_options_set()` gains new parameter `full_path` to set the entire path to the cache, use like `cache_options_set(full_path = yourpath)` (#185) thanks @bomeara for the feature request

### DEFUNCT

* `ft_chunks()` and `ft_tabularize()` were deprecated in the previous version, and are now defunct. See the new package https://github.com/ropensci/pubchunks for the same (and improved) functionality (#146) (#181)
* `ft_get_si()` is defunct. It's been pulled out into a new package suppdata  (#186) (#188)

### BUG FIXES

* Fix to `eupmc_search()` internal function. At some piont Europe PMC changed the way they do paging, what parameters are used, etc. Fixed now. See examples for how to do paging; uses a cursor model instead of a rows/offset model  (#184) thanks @jshleap for the report
* two hopefully fixes for Wiley URLs: a) try a 2nd url pattern if the first fails in general for Wiley links; and b) if highwire links found, replace with a url pattern that should work (#189) thanks @bomeara for the report

### MINOR IMPROVEMENTS

* remaining `httr` code removed, now using `crul` for all HTTP requests (#187)
* filled out `Scopus` and `Crossref TDM` parts in the Authentication section of the package level manual file `?fulltext-package`; and add more details on authentication in the `?ft_get` manual file  (#182) (#183)


fulltext 1.1.0
==============

### NEW FEATURES

* gains new function `cache_file_info()` to get information on possibly bad files in your cache - which you can use to remove files as you see fit (#142) (#174) thx @lucymerobinson for the push
* gains new function `as.ft_data()` to create the same output as `ft_get()` returns, but instead pulls all files from your cache (#142) (#172) thanks @lucymerobinson
* `ft_get()` gains new attribute of a data.frame in the `errors` slot with information on each article and what error we collected or `NA_character_` if none; should help with sorting out problems across all requests (#176)
* scopus option in `ft_search()` gains support for facets; see `?scopus_search` (#170) thanks @lucymerobinson

### BUG FIXES

* fixed bug in `ft_search()` for microsoft academic plugin (#154)
* fixed bug in `ft_search()` for scopus plugin - we weren't looping over requests correctly (#161)
* fix bug in `ft_links()` - when result of `ft_search()` passed to `ft_links` with bad urls or with more than 1 then `ft_links` was failing; fix by filtering on `intended-application` field from Crossref via fix in dependency package `crminer` (#173)
* additional check added in `ft_get()` to check for an invalid file that gave a 200 status code (so passes the status code check) (#175)
* bring back support for `start` parameter in `ft_search()` for Scopus to offset what record a query starts at (#180)
* fix to `ft_get()` for entrez data source to loop internally when more than 50 ids requested to avoid 414 http error (URI too long) (#167) thanks @t088

### MINOR IMPROVEMENTS

* change base url for EuropePMC to `https`
* use `https` for all `doi.org` requests (#155) thanks @katrinleinweber 
* better explanation of the may `opts` parameters in `ft_search()` (#161)
* after many attempts, finally (I think) sorted out Microbiology Society for the `ft_get()` function (#163) thanks to @low-decarie
* mention `suppdata` package in the README (#164)
* clarify language in `ft_collect()` that we are saving files to disk locally only (#165)
* fix typos (#166) (#168) thanks @maelle @AugustT
* removed `whisker` package dependency (#156)
* scopus gains support for more parameters; see `?scopus_search`  (#152)
* updated docs with information on using NCBI Entrez API keys (#159)
* plugin for publisher Instituto de Investigaciones Filologicas to `ft_get()` added (#117) thanks @andreifoldes 
* clarified in `ft_search()` docs that for some sources we loop internally to get whatever number of records the user wants, while others we can not (#162)

### DEPRECATED

Continuing to focus the scope of this package the functions `ft_chunks()` and `ft_tabularize()` are now deprecated, and will be removed (defunct) in a future version of this package. See the new package pubchunks for the same and better functionality. (#181)


fulltext 1.0.1
==============

### BUG FIXES

* Fix bug in internal function `get_ext()` which parses either xml, pdf, or plain text from files on disk - it was failing on Linux maxchines due to a faulty regex (#151)

### MINOR IMPROVEMENTS

* Updated formats vignette, three publishers that used to provide XML no longer do (#150)


fulltext 1.0
============

Check out the fulltext manual (https://books.ropensci.org/fulltext/) for detailed documentation.

`fulltext` has undergone a re-organization, which includes a bump in the major version to `v1` to reinforce the large changes the package has undergone. Changes include:

- Function name standardization with the `ft_` prefix. e,g, `chunks` is now `ft_chunks`
- `ft_get` has undergone major re-organization - biggest of which may be that all full text XML/plain text/PDF goes to disk to simplify the user interface.
- `storr` is now imported to manage mapping between real DOIs and file paths that include normalized DOIs - and aids in the function `ft_table()` for creating a data.frame of text results
- Note that with the `ft_get()` overhaul, the only option is to write to disk. Before we attempted to provide many different options for saving XML and PDF data, but it was too complicated. This has implications for using the output of `ft_get()` - the output is only the paths to the files - use `ft_collect()`  to collect the text if you want to use `ft_chunks()` or other `fulltext` functions downstream.

### NEW FEATURES

* `chunks` changed to `ft_chunks` (#139)
* `collect` changed to `ft_collect` (#139)
* `tabularize` changed to `ft_tabularize` (#139)
* `get_text` changed to `ft_text` (#139)
* `ft_get()` gains new parameter `try_unknown` that attempts to try to find full text for a given DOI/ID even if we don't have code plugins specifically for that publisher. This includes trying to get a full text link from Crossref and the ftdoi.org API (#137)
* Gains function `ft_table` that outputs a data.frame of all downloaded articles with DOIs, filenames, and the text of each article, similar to the `readtext` package (#134)
* Gains function `ft_abstract` for fetching abstracts, including support for getting abstracts from Scopus, Microsoft Academic, Crossref, and PLOS (#98) (#115)
* Microsoft Academic added as another data source both in `ft_abstract` and in `ft_search` via the `microdemic` package (#99) (#115)
* `ft_get()` gains an S3 method for `ft_links` - that is, you can pass the output of `ft_links()` to `ft_get()` (#103)
* `ft_get()` gains many new plugins, including for: Informa, Scientific Societies, Europe PMC, Elsevier, Wiley, xxx (#121) (#112) (#52) (#96) (#120) (#xxx)
* Gains new functions to list available plugins for each of `ft_get()`/`ft_links()`/`ft_search()`: `ft_get_ls()`/`ft_links_ls()`/`ft_search_ls()` (#122)
* `ft_chunks()` gains support for Elsevier XML (#116) (#118)
* Scopus added as a new datasource in both `ft_search()` and `ft_abstract` (#95)
* gains new object `ftxt_cache` for managing/listing details of cached files from `ft_get()`
* `ft_search()` gains Scopus and Microsoft Academic options
* `ft_serialize()` loses `file`, `rcache`, and `redis` options, retaining options for converting between list, JSON, and XML formats.
* The package level manual file at `?fulltext-package` gains new sections on authentication and rate limiting
* gains new manual file for how to interpret warnings when using `ft_get()`; see `?ft_get-warnings`

### DEFUNCT

With the re-focusing of the package these functions seemed out of scope, so have been removed:

* `pdfx`, `pdfx_html`, and `pdfx_targz` are now defunct (#145)
* `ft_extract_corpus` is now defunct

The following functions have changed names. Usually I'd mark functions as deprecated in a version, then defunct in the next version, but since we're moving to `v1` here, it made sense to rip the bandade off and make the old function names defunct.

* `chunks` is now defunct - function name changed to `ft_chunks`
* `tabularize` is now defunct - function name changed to `ft_tabularize`
* `collect` is now defunct - function name changed to `ft_collect`
* `get_text` is now defunct - function name changed to `ft_text`

Other defunct:

* `cache_clear` was never working anyway, and is now removed.

Along with the above changes and others the packages `R.cache`, `rredis`, `digest`, and `tm` have been removed from Imports

### MINOR IMPROVEMENTS

* Replaced `httr` with `crul` mostly throughtout the package (#104)
* We are now using DOIs in file names written to disk from `ft_get()`. DOIs are normalized before using to create file paths (#138)
* Output of `ft_get()` should now be correctly named lists after the publisher and the DOI/ID (#126)
* Switched for eLife DOIs to attempt to get full text links from Crossref API instead of constructing by hand (#127)
* Now using `hoardr` package for managing cached files (#124)
* sending user agent string for this R pkg to the ftdoi.org API when calling it now (#141)
* Documented option to include your email address with Crossref API queries to get in their "fast lane" (#123)
* Documented how to get rate limit headers/information for Elsevier/Scopus requests (#109)
* For text extraction from PDFs - only using `pdftools` now - done away with other options (#82)
* `biorxiv_search` is now exported but the man file is hidden (using `roxygen2` directive `@internal`), so you can still get to the manual file doing `?biorxiv_search`

### BUG FIXES

* New internal PLOS API method instead of a function from `rplos` 
because we need to write to disk instead of return parsed XML 
output (#148)
* `ft_get()` now appropriately using cached version of file if found (#130)
* The `type` parameter in `ft_get()` was ignored previously, now definitely
used. (#128)
* Fix `biorxiv_search` (#97) (#113)


fulltext 0.1.8
===============

### MINOR IMPROVEMENTS

* require newest `rcrossref` and `rplos` versions that use `dplyr::bind_rows()` instead
of `dplyr::rbind_all()` to avoid errors/warnings (#89) (#90)

fulltext 0.1.6
===============

### MINOR IMPROVEMENTS

* More documentation added for the `from` parameter for `ft_get_si()` to clarify its use, and fails better when used inappropriately (#68) (#77)
* `ft_get_si()` now gives file type information as attributes so that downstream uses can access that information instead of having to guess file types (#69)

### BUG FIXES

* Fixes to `ft_get_si()` to work with changes in the publisher Wiley's URL changes (#71) (#73)

fulltext 0.1.4
===============

### NEW FEATURES

* New function `ft_get_si()` to grab supplementary files for any
article (#61) (#62) - thanks @willpearse
* New function `ft_links()` to grab links for the full text version 
of an article from entire output of `ft_search()`, their individual 
components (i.e., data sources), or from character vector of DOIs (#36)

### MINOR IMPROVEMENTS

* Lowercased all column headers from `ft_search()` (#63)
* DRYed out plugins for `ft_get()` to reduce code duplication (#48)

### BUG FIXES

* Fixed bug in `ft_search()` where limit integer was too big (#57)
* Fix to `ft_get()` to create a directory if it doesn't exist 
already (#56)
* Fixed bug in `ft_search()` that caused problems in some scenarios (#55)
* Better error message when pdf passed to `pdfx()` function when 
pdf is not useable, that is, e.g., a scanned pdf (#53)

fulltext 0.1.0
===============

* Released to CRAN.
## Test environments

* local macOS install, R 4.1.0
* ubuntu 16.04 (on github actions), R 4.1.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 3 downstream dependencies. No errors were found. Summary at <https://github.com/ropensci/fulltext/blob/master/revdep/README.md>

--------

This version fixes some bugs.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/fulltext/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/fulltext.git`
* Make sure to track progress upstream (i.e., on our version of `fulltext` at `ropensci/fulltext`) by doing `git remote add upstream https://github.com/ropensci/fulltext.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/fulltext`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.3 Patched (2021-02-02 r79929) |
|os       |macOS Catalina 10.15.7                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-02-12                                  |

# Dependencies

|package   |old   |new        |Δ  |
|:---------|:-----|:----------|:--|
|fulltext  |1.6.0 |1.7.0      |*  |
|htmltools |NA    |0.5.1.9000 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

fulltext
========

[![cran checks](https://cranchecks.info/badges/flavor/release/fulltext)](https://cranchecks.info/pkgs/fulltext)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/fulltext/workflows/R-check/badge.svg)](https://github.com/ropensci/fulltext/actions/)
[![codecov](https://codecov.io/gh/ropensci/fulltext/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/fulltext)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/fulltext)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/fulltext)](https://cran.r-project.org/package=fulltext)

__Get full text research articles__

Checkout the [package docs][docs] and the [fulltext manual][ftbook] to get started.

-----

rOpenSci has a number of R packages to get either full text, metadata, or both from various publishers. The goal of `fulltext` is to integrate these packages to create a single interface to many data sources.

`fulltext` makes it easy to do text-mining by supporting the following steps:

* Search for articles - `ft_search`
* Fetch articles - `ft_get`
* Get links for full text articles (xml, pdf) - `ft_links`
* Extract text from articles / convert formats - `ft_extract`
* Collect all texts into a data.frame - `ft_table`

Previously supported use cases, extracted out to other packages:

* Collect bits of articles that you actually need - moved to package `pubchunks`
* Supplementary data from papers has been moved to the `suppdata`


It's easy to go from the outputs of `ft_get` to text-mining packages such as 
`tm` and `quanteda`

Data sources in `fulltext` include:

* Crossref - via the `rcrossref` package
* Public Library of Science (PLOS) - via the `rplos` package
* Biomed Central
* arXiv - via the `aRxiv` package
* bioRxiv - via the `biorxivr` package
* PMC/Pubmed via Entrez - via the `rentrez` package
* Scopus - internal tooling
* Semantic Scholar - internal tooling
* Many more are supported via the above sources (e.g., _Royal Society Open Science_ is
available via Pubmed)
* We __will__ add more, as publishers open up, and as we have time...See the issues

Authentication: A number of publishers require authentication via API key, and some even more
draconian authentication processes involving checking IP addresses. We are working on supporting
all the various authentication things for different publishers, but of course all the OA content
is already easily available. See the **Authentication** section in `?fulltext-package` after 
loading the package.

We'd love your feedback. Let us know what you think in the issue tracker


## Installation

Stable version from CRAN

```{r eval=FALSE}
install.packages("fulltext")
```

Development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/fulltext")
```

Load library

```{r}
library('fulltext')
```

## Interoperability with other packages downstream

Note: this example not included in vignettes as that would require the two below packages in Suggests here. To see many examples and documentation see the [package docs][docs] and the [fulltext manual][ftbook].

```{r eval=FALSE}
cache_options_set(path = (td <- 'foobar'))
res <- ft_get(c('10.7554/eLife.03032', '10.7554/eLife.32763'), type = "pdf")
library(readtext)
x <- readtext::readtext(file.path(cache_options_get()$path, "*.pdf"))
```

```{r eval=FALSE}
library(quanteda)
quanteda::corpus(x)
```

## Contributors

* Scott Chamberlain
* Will Pearse
* Katrin Leinweber

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/fulltext/issues).
* License: MIT
* Get citation information for `fulltext`: `citation(package = 'fulltext')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.


[docs]: https://docs.ropensci.org/fulltext/
[ftbook]: https://books.ropensci.org/fulltext/
---
title: Getting full text
author: Scott Chamberlain
date: "2020-07-01"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Getting full text}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



The main interface to fetching full text is through `ft_get()`.


```r
library("fulltext")
```

## Structure of the returned object from ft_get()

Simple call, pass in a DOI and say where you want to get data from (by default, it's _plos_)


```r
res <- ft_get('10.1371/journal.pone.0086169', from = 'plos')
```

The article text and metadata is stored in the output object.

The `res` object is a list, with slots for each of the data sources, b/c you can request 
data from more than 1 data source.


```r
names(res)
#> [1] "plos"          "entrez"        "elife"         "pensoft"      
#> [5] "arxiv"         "biorxiv"       "elsevier"      "sciencedirect"
#> [9] "wiley"
```

Let's dig into the `plos` source object, which is another list, including metadata the 
text data itself (in the `data` slot).


```r
res$plos
#> $found
#> [1] 1
#> 
#> $dois
#> [1] "10.1371/journal.pone.0086169"
#> 
#> $data
#> $data$backend
#> [1] "ext"
#> 
#> $data$cache_path
#> [1] "/Users/sckott/Library/Caches/R/fulltext"
#> 
#> $data$path
#> $data$path$`10.1371/journal.pone.0086169`
#> $data$path$`10.1371/journal.pone.0086169`$path
#> [1] "/Users/sckott/Library/Caches/R/fulltext/10_1371_journal_pone_0086169.xml"
#> 
#> $data$path$`10.1371/journal.pone.0086169`$id
#> [1] "10.1371/journal.pone.0086169"
#> 
#> $data$path$`10.1371/journal.pone.0086169`$type
#> [1] "xml"
#> 
#> $data$path$`10.1371/journal.pone.0086169`$error
#> NULL
#> 
#> 
#> 
#> $data$data
#> NULL
#> 
#> 
#> $opts
#> $opts$doi
#> [1] "10.1371/journal.pone.0086169"
#> 
#> $opts$type
#> [1] "xml"
#> 
#> $opts$progress
#> [1] FALSE
#> 
#> 
#> $errors
#>                             id error
#> 1 10.1371/journal.pone.0086169  <NA>
```

Indexing to the `data` slot takes us to another list with metadata and the article


```r
res$plos$data
#> $backend
#> [1] "ext"
#> 
#> $cache_path
#> [1] "/Users/sckott/Library/Caches/R/fulltext"
#> 
#> $path
#> $path$`10.1371/journal.pone.0086169`
#> $path$`10.1371/journal.pone.0086169`$path
#> [1] "/Users/sckott/Library/Caches/R/fulltext/10_1371_journal_pone_0086169.xml"
#> 
#> $path$`10.1371/journal.pone.0086169`$id
#> [1] "10.1371/journal.pone.0086169"
#> 
#> $path$`10.1371/journal.pone.0086169`$type
#> [1] "xml"
#> 
#> $path$`10.1371/journal.pone.0086169`$error
#> NULL
#> 
#> 
#> 
#> $data
#> NULL
```

Going down one more index gets us the data object. There is no actual text as we have to 
collect it from the file on disk. See `ft_collect()` to get the text.


```r
res$plos$data
#> $backend
#> [1] "ext"
#> 
#> $cache_path
#> [1] "/Users/sckott/Library/Caches/R/fulltext"
#> 
#> $path
#> $path$`10.1371/journal.pone.0086169`
#> $path$`10.1371/journal.pone.0086169`$path
#> [1] "/Users/sckott/Library/Caches/R/fulltext/10_1371_journal_pone_0086169.xml"
#> 
#> $path$`10.1371/journal.pone.0086169`$id
#> [1] "10.1371/journal.pone.0086169"
#> 
#> $path$`10.1371/journal.pone.0086169`$type
#> [1] "xml"
#> 
#> $path$`10.1371/journal.pone.0086169`$error
#> NULL
#> 
#> 
#> 
#> $data
#> NULL
```

## Fetching many articles

You can get a bunch of DOIs first, e.g., from PLOS using the `rplos` package


```r
library("rplos")
(dois <- searchplos(q = "*:*", fl = 'id',
   fq = list('doc_type:full', "article_type:\"research article\""), limit = 5)$data$id)
#> [1] "10.1371/journal.pone.0020843" "10.1371/journal.pone.0022257"
#> [3] "10.1371/journal.pone.0023139" "10.1371/journal.pone.0023138"
#> [5] "10.1371/journal.pone.0023119"
ft_get(dois, from = 'plos')
#> <fulltext text>
#> [Docs] 5 
#> [Source] ext - /Users/sckott/Library/Caches/R/fulltext 
#> [IDs] 10.1371/journal.pone.0020843 10.1371/journal.pone.0022257
#>      10.1371/journal.pone.0023139 10.1371/journal.pone.0023138
#>      10.1371/journal.pone.0023119 ...
```

## Different data sources

### Articles from eLife

One article


```r
ft_get('10.7554/eLife.04300', from = 'elife')
#> <fulltext text>
#> [Docs] 1 
#> [Source] ext - /Users/sckott/Library/Caches/R/fulltext 
#> [IDs] 10.7554/eLife.04300 ...
```

Many articles


```r
ft_get(c('10.7554/eLife.04300','10.7554/eLife.03032'), from = 'elife')
#> <fulltext text>
#> [Docs] 2 
#> [Source] ext - /Users/sckott/Library/Caches/R/fulltext 
#> [IDs] 10.7554/eLife.04300 10.7554/eLife.03032 ...
```

### Articles from Frontiers in Pharmacology (publisher: Frontiers)


```r
doi <- '10.3389/fphar.2014.00109'
ft_get(doi, from = "entrez")
#> <fulltext text>
#> [Docs] 1 
#> [Source] ext - /Users/sckott/Library/Caches/R/fulltext 
#> [IDs] 4050532 ...
```

## Search using ft_search()

For example, search entrez, get some DOIs, then fetch some articles


```r
(res <- ft_search(query = 'ecology', from = 'entrez'))
#> Query:
#>   [ecology] 
#> Found:
#>   [PLoS: 0; BMC: 0; Crossref: 0; Entrez: 203062; arxiv: 0; biorxiv: 0; Europe PMC: 0; Scopus: 0; Microsoft: 0] 
#> Returned:
#>   [PLoS: 0; BMC: 0; Crossref: 0; Entrez: 10; arxiv: 0; biorxiv: 0; Europe PMC: 0; Scopus: 0; Microsoft: 0]
res$entrez$data$doi
#>  [1] "10.1111/tesg.12439"               "10.1002/tie.22161"               
#>  [3] "10.1182/bloodadvances.2019001144" "10.1186/s13063-020-04513-w"      
#>  [5] "10.1007/s10109-020-00330-6"       "10.3390/ma13112563"              
#>  [7] "10.3390/molecules25112618"        "10.3390/molecules25112634"       
#>  [9] "10.3390/molecules25112471"        "10.1186/s40168-020-00861-6"
```

Get articles


```r
ft_get(res$entrez$data$doi[1:3], from = 'entrez')
#> <fulltext text>
#> [Docs] 3 
#> [Source] ext - /Users/sckott/Library/Caches/R/fulltext 
#> [IDs] 7323205 7323122 7322952 ...
```

## Collect full text from file on disk

When using `ft_get()` you write the files to disk, and you have to pull text out of them as a 
separate step.


```r
(res <- ft_get('10.1371/journal.pone.0086169', from = 'plos'))
#> <fulltext text>
#> [Docs] 1 
#> [Source] ext - /Users/sckott/Library/Caches/R/fulltext 
#> [IDs] 10.1371/journal.pone.0086169 ...
```

One way to do that is with `ft_collect()`. Before running `ft_collect()` the `data` slot is `NULL`.


```r
res$plos$data$data
#> NULL
```

Run `ft_collect()`


```r
res <- res %>% ft_collect
```

After running `ft_collect()` the `data` slot has the text. If there's more than one article they are named
by the identifier


```r
res$plos$data$data
#> $`10.1371/journal.pone.0086169`
#> {xml_document}
#> <article article-type="research-article" dtd-version="3.0" lang="en" xmlns:mml="http://www.w3.org/1998/Math/MathML" xmlns:xlink="http://www.w3.org/1999/xlink">
#> [1] <front>\n  <journal-meta>\n    <journal-id journal-id-type="nlm-ta">PLoS  ...
#> [2] <body>\n  <sec id="s1">\n    <title>Introduction</title>\n    <p>Since th ...
#> [3] <back>\n  <ack>\n    <p>We thank Joan Silk, Julienne Rutherford, and two  ...
```

---
title: Article formats
author: Scott Chamberlain
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Article formats}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

## Information on article formats

There are various formats in which articles are provided by publishers, including pdf, plain text, xml, json, and more. The following is non-exhaustive table of formats provided by publisher or journal. Examples are included for each publisher if they support the format - click on the `Y` to get an example file.

> Note: many of these journals are also in PMC, where some formats are provided: PDF, ePub.

Publisher | pdf | xml | epub | Notes |
----------| ----| --- | ---- | ----- |
AIRCC | [Y][airccp] | N |  N |... |
arXiv | [Y][arxivp] | N |  N |... |
Bentham | [Y][bep] | N |  N |... |
BiomedCentral | Y | N | N |... |
bioRxiv | [Y][biorxivp] | N | N |... |
CogentOA | Y | N | N |... |
Copernicus | [Y][copp] | [Y][copx] ^[1] | N |... |
De Gruyter | Y | N | N |... |
Dovepress | [Y][dovep] | N | N |... |
eLife | [Y][ep] | [Y][ex] | N |... |
FrontiersIn | Y | Y ^[1] | Y | ReadCube in browser |
Hindawi | [Y][hp] | [Y][hx] | [Y][he] | Don't show XML link on page  |
Karger | Y | N | N |... |
MDPI | [Y][mdpip] | Y | N |... |
Nature | [Y][naturep] | N | N |... |
PeerJ | [Y][peerjp] | [Y][peerjx] | N |... |
Pensoft | [Y][pep] | [Y][pex] | N |... |
PLoS | [Y][plp] | [Y][plx] ^[1] | N |... |
Sage | Y | N | N |... |
Scielo | [Y][scielop] | [Y][scielox] | N | May only be some journals |
SERSC | Y | N | N |... |
Springer Open | Y | N | N |... |

^1: NLM-DTD XML schema - http://dtd.nlm.nih.gov/
^2: Wiley and Elsevier do have a few open access journals each,
which provide PDFs, but no XML. Elsevier has text mining web services
but they are so painful to use that we will not support it here. Do
put pressure on these two enormous publishers to give XML of articles, at
least for their open access journals.

[plp]: http://www.plosone.org/article/fetchObject.action?uri=info%3Adoi%2F10.1371%2Fjournal.pone.0107510&representation=PDF
[peerjp]: https://peerj.com/articles/1142.pdf
[pep]: http://zookeys.pensoft.net/lib/ajax_srv/article_elements_srv.php?action=download_pdf&item_id=4351
[ep]: https://cdn.elifesciences.org/articles/31770/elife-31770-v1.pdf
[hp]: http://downloads.hindawi.com/journals/crid/2014/246965.pdf
[kp]: https://www.karger.com/Article/Pdf/370302
[copp]: http://www.biogeosciences.net/11/7331/2014/bg-11-7331-2014.pdf
[bep]: http://benthamopen.com/contents/pdf/TONEUJ/TONEUJ-9-21.pdf
[arxivp]: http://arxiv.org/pdf/1507.08559v1.pdf?
[biorxivp]: http://biorxiv.org/content/biorxiv/early/2015/07/26/023259.full.pdf
[degrutp]: https://www.degruyter.com/downloadpdf/j/biolet.2014.51.issue-2/biolet-2015-0008/biolet-2015-0008.pdf
[mdpip]: http://www.mdpi.com/1999-4915/7/8/2817/pdf
[airccp]: http://airccse.org/journal/cnc/7115cnc04.pdf
[naturep]: http://www.nature.com/articles/srep12550.pdf
[dovep]: http://www.dovepress.com/getfile.php?fileID=24696
[scielop]: http://www.scielo.br/pdf/cbab/v14n1/04.pdf


[plx]: http://www.plosone.org/article/fetchObjectAttachment.action?uri=info%3Adoi%2F10.1371%2Fjournal.pone.0107510&representation=XML
[peerjx]: https://peerj.com/articles/1142.xml
[pex]: http://zookeys.pensoft.net/lib/ajax_srv/article_elements_srv.php?action=download_xml&item_id=4351
[ex]: https://cdn.elifesciences.org/articles/31770/elife-31770-v1.xml
[hx]: http://downloads.hindawi.com/journals/tswj/2014/649260.xml
[copx]: http://www.biogeosciences.net/11/7331/2014/bg-11-7331-2014.xml
[scielox]: http://www.scielo.br/scieloOrg/php/articleXML.php?pid=S1984-70332014000100004&lang=en

[he]: http://downloads.hindawi.com/journals/crid/2014/246965.epub
---
title: fulltext introduction
author: Scott Chamberlain
date: "2020-07-01"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{fulltext introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



`fulltext` is a package to facilitate text mining. It focuses on open access journals. This package makes it easier to search for articles, download those articles in full text if available, convert pdf format to plain text, and extract text chunks for vizualization/analysis. We are planning to add bits for analysis in future versions. The steps in bullet form:

* Search
* Retrieve
* Convert
* Text
* Extract

## Load fulltext


```r
library("fulltext")
```

## Search for articles

Search for the term _ecology_ in PLOS journals.


```r
(res1 <- ft_search(query = 'ecology', from = 'plos'))
#> Query:
#>   [ecology] 
#> Found:
#>   [PLoS: 52809; BMC: 0; Crossref: 0; Entrez: 0; arxiv: 0; biorxiv: 0; Europe PMC: 0; Scopus: 0; Microsoft: 0] 
#> Returned:
#>   [PLoS: 10; BMC: 0; Crossref: 0; Entrez: 0; arxiv: 0; biorxiv: 0; Europe PMC: 0; Scopus: 0; Microsoft: 0]
```

Each publisher/search-engine has a slot with metadata and data


```r
res1$plos
#> Query: [ecology] 
#> Records found, returned: [52809, 10] 
#> License: [CC-BY] 
#> # A tibble: 10 x 1
#>    id                          
#>    <chr>                       
#>  1 10.1371/journal.pone.0001248
#>  2 10.1371/journal.pone.0059813
#>  3 10.1371/journal.pone.0080763
#>  4 10.1371/journal.pone.0220747
#>  5 10.1371/journal.pone.0155019
#>  6 10.1371/journal.pone.0175014
#>  7 10.1371/journal.pone.0150648
#>  8 10.1371/journal.pone.0208370
#>  9 10.1371/journal.pcbi.1003594
#> 10 10.1371/journal.pone.0102437
```

## Get full text links

`ft_links()` - get links for articles (xml and pdf).


```r
w <- ft_search(query = 'biology', from = 'entrez', limit = 5)
ft_links(w)
#> <fulltext links>
#> [Found] 5 
#> [IDs] ID_32559462 ID_32460012 ID_31748751 ID_29060747 ID_27680492 ...
```

Or pass in DOIs directly


```r
ft_links(w$entrez$data$doi, from = "entrez")
#> <fulltext links>
#> [Found] 5 
#> [IDs] ID_32559462 ID_32460012 ID_31748751 ID_29060747 ID_27680492 ...
```

## Get full text

Using the results from `ft_search()` we can grab full text of some articles


```r
search_res <- ft_search(query = 'ecology', from = 'plos')
(out <- ft_get(search_res))
#> <fulltext text>
#> [Docs] 10 
#> [Source] ext - /Users/sckott/Library/Caches/R/fulltext 
#> [IDs] 10.1371/journal.pone.0001248 10.1371/journal.pone.0059813
#>      10.1371/journal.pone.0080763 10.1371/journal.pone.0220747
#>      10.1371/journal.pone.0155019 10.1371/journal.pone.0175014
#>      10.1371/journal.pone.0150648 10.1371/journal.pone.0208370
#>      10.1371/journal.pcbi.1003594 10.1371/journal.pone.0102437 ...
```

Dig in to the PLOS data


```r
out$plos$data$path[[1]]
#> $path
#> [1] "/Users/sckott/Library/Caches/R/fulltext/10_1371_journal_pone_0001248.xml"
#> 
#> $id
#> [1] "10.1371/journal.pone.0001248"
#> 
#> $type
#> [1] "xml"
#> 
#> $error
#> NULL
```

## Extract text from pdfs

Ideally for text mining you have access to XML or other text based formats. However, 
sometimes you only have access to PDFs. In this case you want to extract text 
from PDFs. `fulltext` can help with that. 

You can extract from any pdf from a file path, like:


```r
path <- system.file("examples", "example1.pdf", package = "fulltext")
ft_extract(path)
#> <document>/Library/Frameworks/R.framework/Versions/4.0/Resources/library/fulltext/examples/example1.pdf
#>   Title: Suffering and mental health among older people living in nursing homes---a mixed-methods study
#>   Producer: pdfTeX-1.40.10
#>   Creation date: 2015-07-17
```

Let's search for articles from arXiv, a preprint service. Here, get pdf from 
an article with ID `cond-mat/9309029`:


```r
res <- ft_get('cond-mat/9309029', from = "arxiv")
res2 <- ft_extract(res)
# the first page
res2$arxiv$data$data$`cond-mat/9309029`[1]
#> [1] "                                                                                                                          FERMILAB-PUB-93/15-T\n                                                                                                                                       March 1993\n                                                                                                                             Revised: January 1994\n                                           The Thermodynamics and Economics of Waste\n                                           Dallas C. Kennedy, Research Associate,\narXiv:cond-mat/9309029v8 26 Jan 1994\n                                           Fermi National Accelerator Laboratory,\n                                           P.O. Box 500 MS106, Batavia, Illinois 60510 USA∗\n                                                                                         Abstract\n                                       The increasingly relevant problem of natural resource use and waste production, disposal, and reuse is\n                                       examined from several viewpoints: economic, technical, and thermodynamic. Alternative economies are\n                                       studied, with emphasis on recycling of waste to close the natural resource cycle. The physical nature of\n                                       human economies and constraints on recycling and energy efficiency are stated in terms of entropy and\n                                       complexity.\n                                           What is a cynic? A man who knows the price of everything, and the value of nothing.\n                                                                                                                                      Oscar Wilde [1]\n                                            Our planet is finite in size and, except for a few energy and matter flows, its biosphere forms a closed\n                                       system. The envelope that makes life possible extends a small distance below the surface of the Earth and\n                                       less than a hundred miles into the atmosphere. While almost closed, the biosphere is not static, but is\n                                       constantly changing, moving flows of energy, air, water, soil, and life around in a shifting, never-repeating\n                                       pattern. The combination of general physical laws and the specific properties of the Earth places important\n                                       constraints on the activities of life, some embodied in the metabolism and forms of living creatures, others\n                                       imprinted into their genes by selective effects. All life needs sources of energy for sustenance and imposes\n                                       a burden of waste on its environment. Since this waste is usually harmful to the creatures emitting it, the\n                                       environment must, if these creatures are to continue living, break the waste down into less toxic forms and\n                                       possibly reuse it.\n                                            The growing dominance of humankind over the planet, both by technological power and by numbers,\n                                       imposes certain costs on the biosphere, sometimes a result of conscious attempts at controlling Nature, but\n                                       more frequently by unwitting influence. Moreover, the burden of carrying the activities of human sustenance,\n                                       unlike that of other animals, cannot be understood by considering the physical and biological activities of each\n                                       person in isolation. Because of their unique position at the top of the food chain, their tool-making abilities,\n                                       and the co-operative character of human activities (the division of labor or specialization), the physical and\n                                       biological aspects must be considered together in the context of the peculiarities of economic life [2]. The\n                                       economic aspects take on an independent importance because of the absence of any automatic, given means\n                                       of human subsistence and the extension of individual self-sufficiency by surplus production and trading [3].\n                                       This point of view is necessary for comprehending the ecological significance of all economies more elaborate\n                                       than the simplest subsistence or household economies, up to and including the most sophisticated systems\n                                       of technology and trade. On the other hand, the formulation of economic theory has generally taken place,\n                                         * Present address: Department of Physics, University of Florida, Gainesville, Florida 32611 USA\n                                                                                              1\n"
```

## Extract text chunks

Requires the pubchunks library.

We have a few functions to help you pull out certain parts of an article. 
For example, perhaps you want to get just the authors from your articles, 
or just the abstracts. 

Here, we'll search for some PLOS articles, then get their full text, then
extract various parts of each article with `pub_chunks()`.


```r
if (requireNamespace("pubchunks")) {
library(pubchunks)

# Search
res <- ft_search(query = "ecology", from = "plos")
(x <- ft_get(res))

# Extract DOIs
x %>% ft_collect() %>% pub_chunks("doi")

# Extract DOIs and categories
x %>% ft_collect() %>% pub_chunks(c("doi", "title"))


# `pub_tabularize` attempts to help you put the data that comes out of 
# `pub_chunks()` in to a `data.frame`, that we all know and love. 
x %>% ft_collect() %>% pub_chunks(c("doi", "history")) %>% pub_tabularize()

}
#> $plos
#> $plos$`10.1371/journal.pone.0001248`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0001248       2007-07-02       2007-11-06       plos
#> 
#> $plos$`10.1371/journal.pone.0059813`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0059813       2012-09-16       2013-02-19       plos
#> 
#> $plos$`10.1371/journal.pone.0080763`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0080763       2013-08-15       2013-10-16       plos
#> 
#> $plos$`10.1371/journal.pone.0220747`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0220747       2019-04-17       2019-07-21       plos
#> 
#> $plos$`10.1371/journal.pone.0155019`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0155019       2015-09-22       2016-04-22       plos
#> 
#> $plos$`10.1371/journal.pone.0175014`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0175014       2016-09-23       2017-03-20       plos
#> 
#> $plos$`10.1371/journal.pone.0150648`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0150648       2015-09-19       2016-02-16       plos
#> 
#> $plos$`10.1371/journal.pone.0208370`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0208370       2018-11-13       2019-01-15       plos
#> 
#> $plos$`10.1371/journal.pcbi.1003594`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pcbi.1003594       2014-01-09       2014-03-14       plos
#> 
#> $plos$`10.1371/journal.pone.0102437`
#>                            doi history.received history.accepted .publisher
#> 1 10.1371/journal.pone.0102437       2013-11-27       2014-06-19       plos
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ma_utils.R
\name{microsoft-internals}
\alias{microsoft-internals}
\alias{microsoft_search}
\alias{microsoft_links}
\title{Microsoft Academic search}
\usage{
microsoft_search(
  query,
  count = 10,
  offset = 0,
  orderby = NULL,
  atts = c("Id", "DN", "VFN", "DOI", "D"),
  key = NULL,
  ...
)

microsoft_links(
  query,
  count = 10,
  offset = 0,
  orderby = NULL,
  atts = c("Id", "AA.AuN", "J.JN", "Ti", "Y", "E", "CC"),
  key = NULL,
  ...
)
}
\arguments{
\item{query}{(character) query terms}

\item{count}{(integer) number of records to return. default: 10}

\item{offset}{(integer) record to start at. default: 0}

\item{orderby}{(character) field to sort results by}

\item{atts}{(character) character vector of fields to return}

\item{key}{(character) microsoft academic API key, see
\code{Authentication} section in \link{fulltext-package}}
}
\description{
Wraps \code{microdemic::ma_evaluate}
}
\examples{
\dontrun{
microsoft_search2(query = "Y='19'...",
  key = Sys.getenv("MICROSOFT_ACADEMIC_KEY"))
}
}
\references{
https://academic.microsoft.com/
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_browse.R
\name{ft_browse}
\alias{ft_browse}
\title{Browse an article in your default browser}
\usage{
ft_browse(x, browse = TRUE)
}
\arguments{
\item{x}{An object of class \code{ft_data} - the output from a call to
\code{\link[=ft_get]{ft_get()}}}

\item{browse}{(logical) Whether to browse (default) or not. If \code{FALSE},
return the url.}
}
\description{
Browse an article in your default browser
}
\examples{
\dontrun{
x <- ft_get('10.7554/eLife.04300', from='elife')
ft_browse(x)
ft_browse(x, browse=FALSE)

ft_browse( ft_get('10.3389/fphar.2014.00109', from="entrez") )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_table.R
\name{ft_table}
\alias{ft_table}
\title{Collect metadata and text into a data.frame}
\usage{
ft_table(path = NULL, type = NULL, encoding = NULL, xml_extract_text = TRUE)
}
\arguments{
\item{path}{a directory path, must exist}

\item{type}{(character) type of files to get. Default is \code{NULL} which gets all types.
Can be one of pdf, xml, or plain (file extensions: pdf, xml, and txt, respectively)}

\item{encoding}{(character) encoding, if \code{NULL} we get it from \code{getOption("encoding")}}

\item{xml_extract_text}{(logical) for XML, should we extract the text (\code{TRUE}) or
return a string as XML (\code{FALSE}). Default: \code{TRUE}}
}
\description{
Facilitates downstream processing with text mining packages
by providing metadata and full text in a tidy data.frame format
}
\details{
You can alternatively use \code{readtext::readtext()} or similar functions
to achieve a similar outcome.
}
\examples{
\dontrun{
if (interactive()) {
## from a directory path
x <- ft_table()
x

## only xml
ft_table(type = "xml")

## only pdf
ft_table(type = "pdf")

## don't pull text out of xml, just give back the xml please
x <- ft_table(xml_extract_text = FALSE)
x
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_clean.R
\name{cache_file_info}
\alias{cache_file_info}
\title{Get information on possibly bad files in your cache}
\usage{
cache_file_info()
}
\value{
list, with three elements:
\itemize{
\item xml_not_valid: xml files that could not be read in with
\code{xml2::read_xml()}
\item xml_abstract_only: xml files that only have abstracts.
you can of choose to retain these if you like
\item pdf_not_valid: pdf files that could not be read in with
\code{pdftools::pdf_info()}
}
}
\description{
Get information on possibly bad files in your cache
}
\details{
This function only identifies possibly bad files.
You have to remove/delete them yourself. See example for
how to do so. You can also open up your cache folder and
delete them that way as well.
}
\examples{
\dontrun{
# identify likely bad files
res <- cache_file_info()

# you can remove them yourself, e.g.,
# invisible(lapply(res$xml_abstract_only, unlink))
}
}
\seealso{
Other caching-functions: 
\code{\link{cache}},
\code{\link{ftxt_cache}}
}
\concept{caching-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{get_text}
\alias{get_text}
\title{This function is defunct.}
\usage{
get_text(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_search.R
\name{ft_search}
\alias{ft_search}
\alias{ft_search_ls}
\title{Search for full text}
\usage{
ft_search(
  query,
  from = "plos",
  limit = 10,
  start = 0,
  plosopts = list(),
  bmcopts = list(),
  crossrefopts = list(),
  entrezopts = list(),
  arxivopts = list(),
  biorxivopts = list(),
  euroopts = list(),
  scopusopts = list(),
  maopts = list(),
  ...
)

ft_search_ls()
}
\arguments{
\item{query}{(character) Query terms}

\item{from}{(character) Source to query, one or more of \code{"plos"}, \code{"bmc"},
\code{"crossref"}, \code{"entrez"}, \code{"arxiv"}, \code{"biorxiv"}, \code{"europmc"}, \code{"scopus"},
or \code{"ma"}}

\item{limit}{(integer) Number of records to return. default: 10. See also
Pagination section below.}

\item{start}{(integer) Record number to start at. Only used for
'scopus' right now. default: 0. Note that with some data sources we loop
internally to get all the results you want with the \code{limit} parameter, so
\code{start} in those cases will be ignored. See \strong{Looping} section below.}

\item{plosopts}{(list) PLOS options, a named list. See \code{\link[rplos:searchplos]{rplos::searchplos()}}}

\item{bmcopts}{(list) BMC options, a named list. See \code{\link[=bmc_search]{bmc_search()}}}

\item{crossrefopts}{(list) Crossref options, a named list.
See \code{\link[rcrossref:cr_works]{rcrossref::cr_works()}}}

\item{entrezopts}{(list) Entrez options, a named list.
See \code{\link[rentrez:entrez_search]{rentrez::entrez_search()}}}

\item{arxivopts}{(list) arxiv options, a named list.
See \code{\link[aRxiv:arxiv_search]{aRxiv::arxiv_search()}}}

\item{biorxivopts}{(list) biorxiv options, a named list.
See \code{\link[=biorxiv_search]{biorxiv_search()}}}

\item{euroopts}{(list) Euro PMC options, a named list. See \code{\link[=eupmc_search]{eupmc_search()}}}

\item{scopusopts}{(list) Scopus options, a named list.
See \code{\link[=scopus_search]{scopus_search()}}}

\item{maopts}{(list) Microsoft Academic options, a named list. See
\code{\link[=microsoft_search]{microsoft_search()}}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}, see
examples below. curl options are ignored for: arxiv (however, you
can wrap your call to arxiv in \code{httr::with_config})}
}
\value{
An object of class \code{ft}, and objects of class \code{ft_ind}
within each source. You can access each data source with \code{$}
}
\description{
\code{ft_search} is a one stop shop for searching for articles
across many publishers and repositories. We currently support search for
PLOS via the  \pkg{rplos} package, Crossref via the \pkg{rcrossref}
package, Entrez via the \pkg{rentrez} package, arXiv via the \pkg{aRxiv}
package, and BMC, Biorxiv, EuropePMC, and Scopus via internal helper
functions in this package.

Many publishers' content is searchable via Crossref and Entrez - of course
this doesn't mean that we can get full text for those articles. In the
output objects of this function, we attempt to help by indicating what
license is used for articles.
}
\details{
Each of \code{plosopts}, \code{scopusopts}, etc. expect
a named list.

See \strong{Rate Limits} and \strong{Authentication} in
\link{fulltext-package} for rate limiting and authentication information,
respectively

See https://dev.elsevier.com/sc_search_tips.html for help/tips
on searching with Scopus
}
\note{
for all \verb{*opts} parameters, see the function linked to in
the parameter definition for what you can pass to it.
}
\section{Looping}{

Note that we necessarily have to treat different sources/publishers
differently internally in this function. Some we can search and get
back as many results as desired automatically, while with others you'd
have to manually iterate through to get all your results.
Notes on different sources:
\itemize{
\item PLOS: \code{\link[rplos:searchplos]{rplos::searchplos()}} used and includes internal looping of
requests
\item BMC: using internal function \code{bmc_search} that does not
loop, so you have to iterate through requests manually
\item Crossref: \code{\link[rcrossref:cr_works]{rcrossref::cr_works()}} used, but does not include
internal looping of requests, but the max limit for one request
is relatively high at 1000
\item Entrez: \code{\link[rentrez:entrez_search]{rentrez::entrez_search()}} used, but does not include
internal looping of requests
\item arXiv: \code{\link[aRxiv:arxiv_search]{aRxiv::arxiv_search()}} used and includes internal looping
of requests
\item BiorXiv: using internal function \code{biorxiv_search} that does not
loop, so you have to iterate through requests manually
\item Europe BMC: using internal function \code{eupmc_search} that does not
loop, so you have to iterate through requests manually
\item Scopus: using internal function \code{scopus_search_loop} that does
include internal looping of requests
\item Microsoft AR: using internal function \code{microsoft_search} that does not
loop, so you have to iterate through requests manually
}
}

\section{Pagination}{

For each data source you can pass named parameters to a list matching
that data source name, e.g., \code{plosopts} for PLOS. If you pass pagination
parameters per data source they will override the global pagination
parameters set in \code{ft_search()}. They are for each data source
(limit parameter name/offset parameter name):
\itemize{
\item PLOS: \code{limit}/\code{start}
\item Crossref: \code{limit}/\code{offset}
\item BMC: \code{limit}/\code{offset}
\item Entrez: \code{retmax}/\code{retstart}
\item Europe PMC: \code{per_page}/(see \code{\link[=eupmc_search]{eupmc_search()}})
\item arXiv: \code{limit}/\code{start}
\item BioRxiv: \code{limit}/none
\item Scopus: \code{count}/\code{start}
\item Microsoft Academic: \code{count}/\code{offset}
}
}

\examples{
# List publishers included
ft_search_ls()

\dontrun{
# Plos
(res1 <- ft_search(query='ecology', from='plos'))
res1$plos
ft_search(query='climate change', from='plos', limit=500, 
  plosopts=list(
   fl=c('id','author','eissn','journal','counter_total_all',
   'alm_twitterCount')))

# Crossref
(res2 <- ft_search(query='ecology', from='crossref'))
res2$crossref

# BioRxiv
(res <- ft_search(query='owls', from='biorxiv'))
res$biorxiv

# Entrez
(res <- ft_search(query='ecology', from='entrez'))
res$entrez

# arXiv
(res <- ft_search(query='ecology', from='arxiv'))
res$arxiv

# BMC - can be very slow
(res <- ft_search(query='ecology', from='bmc'))
res$bmc

# Europe PMC
(res <- ft_search(query='ecology', from='europmc'))
res$europmc
## get the next batch of results, using the cursorMark result
ft_search(query='ecology', from='europmc', 
  euroopts = list(cursorMark = res$europmc$cursorMark))

# Scopus
(res <- ft_search(query = 'ecology', from = 'scopus', limit = 100,
   scopusopts = list(key = Sys.getenv('ELSEVIER_SCOPUS_KEY'))))
res$scopus
## pagination
(res <- ft_search(query = 'ecology', from = 'scopus', 
   scopusopts = list(key = Sys.getenv('ELSEVIER_SCOPUS_KEY')), limit = 5))
## lots of results
(res <- ft_search(query = "ecology community elk cow", from = 'scopus', 
   limit = 100,
   scopusopts = list(key = Sys.getenv('ELSEVIER_SCOPUS_KEY'))))
res$scopus
## facets
## FIXME: apparently I don't have access to facets anymore
# (res <- ft_search(query = 'ecology', from = 'scopus', 
#   scopusopts = list(
#     key = Sys.getenv('ELSEVIER_SCOPUS_KEY'), 
#     facets = "subjarea(count=5)"
#   ), limit = 5))
# res$scopus

# PLOS, Crossref, and arxiv
(res <- ft_search(query='ecology', from=c('plos','crossref','arxiv')))
res$plos
res$arxiv
res$crossref

# Microsoft academic search
key <- Sys.getenv("MICROSOFT_ACADEMIC_KEY")
(res <- ft_search("Y='19'...", from = "microsoft", maopts = list(key = key)))
res$ma$data$DOI

# curl options
ft_search(query='ecology', from='plos', verbose = TRUE)
ma_key <- Sys.getenv("MICROSOFT_ACADEMIC_KEY")
ft_search("Y='19'...", from='microsoft', maopts = list(key = ma_key),
  verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.ft_data.R
\name{as.ft_data}
\alias{as.ft_data}
\title{Coerce directory of papers to ft_data object}
\usage{
as.ft_data(path = NULL)
}
\arguments{
\item{path}{cache path. if not given, we use the default
cache path. Default: \code{NULL}}
}
\value{
an object of class \code{ft_data}
}
\description{
create the same object that \code{\link[=ft_get]{ft_get()}} outputs
from your cached files - without having to run
\code{\link[=ft_get]{ft_get()}} again
}
\details{
We use an internal store of identifiers to keep track of files.
These identifiers are in the output of \code{\link[=ft_get]{ft_get()}} and you can see them
in that output.  If a file does not have a matching entry in our index
of files (e.g., if you drop a file into the cache location as in the
example below), then we assign it an index based on the file path; we'd
ideally use an article DOI or similar but we can not safely retrieve it
with just a file path.
}
\examples{
# put a file in the cache in case there aren't any
dir <- file.path(tempdir(), "testing")
dir.create(dir)
file <- system.file("examples", "elife.xml", package = "fulltext")
writeLines(readLines(file), tempfile(tmpdir = dir, fileext = ".xml"))

# call as.ft_data
x <- as.ft_data(path = dir)

# output lives underneath a special list index "cached" 
#   representing already present files
x$cached

\dontrun{
# collect chunks
if (requireNamespace("pubchunks")) {
  library(pubchunks)
  res <- ft_collect(x)
  pub_chunks(res, c("doi", "title")) \%>\% pub_tabularize()
}
}
}
\seealso{
\code{\link[=ft_get]{ft_get()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_extract.R
\name{ft_extract}
\alias{ft_extract}
\title{Extract text from a single pdf document}
\usage{
ft_extract(x)
}
\arguments{
\item{x}{Path to a pdf file, or an object of class \code{ft_data}, the
output from \code{\link[=ft_get]{ft_get()}}}
}
\value{
An object of class \code{pdft_char} in the case of character input,
or of class \code{ft_data} in the case of \code{ft_data} input
}
\description{
\code{ft_extract} attemps to make it easy to extract text from
PDFs, using \pkg{pdftools}. Inputs can be either paths to PDF
files, or the output of \code{\link[=ft_get]{ft_get()}} (class \code{ft_data}).
}
\examples{
\dontrun{
path <- system.file("examples", "example1.pdf", package = "fulltext")
(res <- ft_extract(path))

# use on output of ft_get() to extract pdf to text
## arxiv
res <- ft_get('cond-mat/9309029', from = "arxiv")
res2 <- ft_extract(res)
res$arxiv$data
res2$arxiv$data

## biorxiv
res <- ft_get('10.1101/012476')
res2 <- ft_extract(res)
res$biorxiv$data
res2$biorxiv$data
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_links.R
\name{ft_links}
\alias{ft_links}
\alias{ft_links_ls}
\title{Get full text links}
\usage{
ft_links(x, from = NULL, ...)

ft_links_ls()
}
\arguments{
\item{x}{One of \code{ft}, \code{ft_ind}, or a character string of DOIs.}

\item{from}{Source to query. Ignored when \code{ft_ind} class passed.}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient} (plos, bmc,
crossref) or \code{httr::GET()} (entrez), see examples below}
}
\value{
An object of class ft_links, with either a list or data.frame for
each DOI, with links for XML and PDF links (typically).
}
\description{
Get full text links
}
\details{
Inputs can be an object of class \code{ft}, \code{ft_ind}, or a
character string of DOIs. You can specify a specific source for four sources
(PLOS, BMC, Crossref, and Entrez), but any other publishers we guess the
publisher form the input DOI(s), then attempt to generate full text links
based on the publisher (if found). Of course, guessing the publisher makes
things slower as it requires an HTTP request.

Strategy varies by publisher. For some we can construct XML and PDF links
only from the DOI. For others, we need to make an HTTP request to the
publisher to get additional information - this of course makes things slower.

See \strong{Rate Limits} and \strong{Authentication} in
\link{fulltext-package} for rate limiting and authentication information,
respectively
}
\examples{
# List publishers included
ft_links_ls()

\dontrun{
# Entrez
(res1 <- ft_search(query='ecology', from='entrez'))
res1$entrez$data$doi
## directly from ft_search output
(out <- ft_links(res1))
out$entrez
out$entrez$data[[1]]
## directly individual elements of ft_search output
(out <- ft_links(res1$entrez))
out$entrez
## from character vector of DOIs
x <- c("10.1371/journal.pone.0086169", "10.1016/j.ympev.2010.07.013")
(out2 <- ft_links(x, from = "entrez"))
out2$entrez

# Crossref
(res2 <- ft_search(query='ecology', from='crossref'))
res2$crossref$data$doi
## directly from ft_search output
(out <- ft_links(res2))
out$crossref
out$crossref$data[[1]]
## directly individual elements of ft_search output
(out <- ft_links(res2$crossref))
out$crossref
## from character vector of DOIs
x <- c("10.1016/s1754-5048(14)00139-1", 
       "10.1016/b978-0-12-378260-1.50017-8")
(out2 <- ft_links(x, from = "crossref"))
out2$crossref

# PLOS
(res3 <- ft_search(query='ecology', from='plos', plosopts=list(
   fl=c('id','author','eissn','journal','counter_total_all',
        'alm_twitterCount'))))
res3$plos$data$id
## directly from ft_search output
(out <- ft_links(res3))
out$plos
out$plos$data[[1]]
## directly individual elements of ft_search output
(out <- ft_links(res3$plos))
out$plos
## from character vector of DOIs
x <- c("10.1371/journal.pone.0017342", "10.1371/journal.pone.0091497")
out3 <- ft_links(x)
out3$plos

# BMC
(res <- ft_search(query='ecology', from='bmc'))
res$bmc
## directly from ft_search output
(out <- ft_links(res))
out$bmc
out$bmc$data[[1]]
## directly individual elements of ft_search output
(out <- ft_links(res$bmc))
out$bmc

# Character input
out4 <- ft_links('10.1371/journal.pone.0086169')
out4$plos

# other publishers
## elife
res <- ft_links(c('10.7554/eLife.03032', '10.7554/eLife.02747'))
res$elife

## peerj
ft_links('10.7717/peerj.228')
ft_links(c('10.7717/peerj.228', '10.7717/peerj.1200'))

## wiley
res <- ft_links('10.1006/asle.2001.0035', from = "crossref")
res$crossref$data[[1]]$url

## informa
res <- ft_links('10.1174/02134749660569378', from = "crossref")
res$crossref$data[[1]]$url

## frontiersin
(res <- ft_links('10.3389/fphar.2014.00109'))
res$frontiersin

## copernicus
(res <- ft_links('10.5194/angeo-31-2157-2013'))
res$copernicus

## cogent
(res <- ft_links('10.1080/23311916.2014.938430'))
res$informa

## bmc
(res <- ft_links('10.1186/2049-2618-2-7'))
res$springer
(res <- ft_links('10.1186/2049-2618-2-7', from = "bmc"))

## Many publishers, elife and peerj
res <- ft_links(c('10.7554/eLife.03032', '10.7717/peerj.228'))
res$elife
res$peerj


# curl options
ft_links("10.2458/v17i1.21696", from = "crossref", verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chunks.R
\name{ft_tabularize}
\alias{ft_tabularize}
\title{Extract chunks of data from articles}
\usage{
ft_tabularize(...)
}
\description{
Extract chunks of data from articles
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmc_utils.R
\name{bmc_search}
\alias{bmc_search}
\title{Search for gene sequences available for a species from NCBI.}
\usage{
bmc_search(query, limit = 10, offset = 1, key = NULL, ...)
}
\arguments{
\item{query}{Search terms.}

\item{limit}{Number of records to return. Default 10.}

\item{offset}{Record number to start at. Default: 1}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list of length 2
}
\description{
Search for gene sequences available for a species from NCBI.
}
\examples{
\dontrun{
bmc_search(query='ecology')
bmc_search('fire', limit=3)
bmc_search('fire', limit=2, page=1)
bmc_search('fire', limit=2, page=2)

# Search, then get full text
out <- bmc_search('ecology')
(urls <- vapply(out$records$url, "[[", "", 'value'))
browseURL(urls[1])

# curl debugging help
bmc_search('ecology', verbose = TRUE)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_ftdmurl.R
\name{as_ftdmurl}
\alias{as_ftdmurl}
\alias{as_ftdmurl.ftdmurl}
\alias{as_ftdmurl.character}
\title{Coerce a url to a tdmurl with a specific type}
\usage{
as_ftdmurl(url, type, doi = NULL, member = NULL, intended_application = NULL)

\method{as_ftdmurl}{ftdmurl}(url, type, doi = NULL, member = NULL, intended_application = NULL)

\method{as_ftdmurl}{character}(url, type, doi = NULL, member = NULL, intended_application = NULL)
}
\arguments{
\item{url}{(character) A URL.}

\item{type}{(character) One of 'xml' (default), 'html', 'plain', 'pdf',
'unspecified', or 'all'}

\item{doi}{(character) A DOI, optional, default: \code{NULL}}

\item{member}{(character) Crossref member id. optional}

\item{intended_application}{(character) intended application string,
optional}
}
\description{
A tmd url is just a URL with some attributes to make it easier
to handle within other functions in this package.
}
\examples{
as_ftdmurl("http://downloads.hindawi.com/journals/bmri/2014/201717.xml",
   "xml")
as_ftdmurl("http://downloads.hindawi.com/journals/bmri/2014/201717.pdf",
   "pdf")
out <-
 as_ftdmurl("http://downloads.hindawi.com/journals/bmri/2014/201717.pdf",
   "pdf", "10.1155/2014/201717")
attributes(out)
identical(attr(out, "type"), "pdf")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftdoi_other.R
\name{prefix_local}
\alias{prefix_local}
\title{prefix local}
\usage{
prefix_local(doi)
}
\arguments{
\item{doi}{(characte) a doi}
}
\value{
a named list with: prefix, member, name
}
\description{
prefix local
}
\examples{
\dontrun{
prefix_local('10.3390/ani4010082')
}
}
\seealso{
Other ftdoi: 
\code{\link{ftd_doi}()},
\code{\link{ftd_fetch_patterns}()},
\code{\link{ftd_members}()},
\code{\link{ftd_prefixes}()},
\code{\link{ftdoi_cache}}
}
\concept{ftdoi}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{ft_browse_sections}
\alias{ft_browse_sections}
\title{This function is defunct.}
\usage{
ft_browse_sections(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftd_doi.R
\name{ftd_doi}
\alias{ftd_doi}
\title{DOI}
\usage{
ftd_doi(doi, ...)
}
\arguments{
\item{doi}{(character) one or more DOIs. required}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
data.frame of rows equal to number of DOIs supplied, with columns:
\itemize{
\item doi: the doi
\item url: url for the article
\item content_type: content type of the article format
\item issn: ISSN for the journal containing the DOI
\item member_name: Crossref member name
\item member_url: Crossref member url
}
}
\description{
DOI
}
\examples{
\dontrun{
# pensoft
ftd_doi(doi = '10.3897/zookeys.594.8768')
ftd_doi(doi = '10.3897/mycokeys.54.34571')
ftd_doi(doi = '10.3897/phytokeys.99.26489')
ftd_doi(doi = '10.3897/subtbiol.13.6719')
# plos
ftd_doi(doi = '10.1371/journal.pgen.1006546')
ftd_doi(c('10.1371/journal.pgen.1006546', '10.1371/journal.pbio.1001809'))
# mdpi
ftd_doi('10.3390/ani4010082')
ftd_doi(doi = c('10.3390/ani4010082', "10.3390/ijms19040965",
  "10.3390/rs9010083"))
# frontiers
ftd_doi(doi = '10.3389/fmed.2015.00081')
# informa
ftd_doi(doi = '10.4324/9780203538333')
# thieme
ftd_doi(doi = '10.1055/s-0042-103414')
# peerj
ftd_doi(doi = '10.7717/peerj.991')
ftd_doi(doi = '10.7717/peerj-cs.39')
# American Phyiscal Society
ftd_doi(doi = '10.1103/physreve.68.067402')
# Royal Society of Chemistry
ftd_doi(doi = '10.1039/c4ra04415k')
# Karger
ftd_doi(doi = '10.1159/000360225')
ftd_doi(doi = c("10.1159/000094345","10.1159/000086754"))
# Trans Tech Publications
ftd_doi(doi = '10.4028/www.scientific.net/msf.702-703.774')
# Emerald
ftd_doi(doi = '10.1108/00251740210413370')
# mixed publishers
ftd_doi(doi = c("10.1371/journal.pgen.1006546","10.1159/000086754"))
# Pleiades
ftd_doi(doi = '10.1134/s1063784215120075')
# Instituto de Investigaciones Filologicas
ftd_doi(doi = '10.1016/s0185-3082(14)70398-0')
ftd_doi(c('10.1016/s0185-2574(13)71376-5', '10.19130/iifl.nt.1997.15.0.650'))
# Sage
ftd_doi(doi = '10.1177/0267659117690248')
ftd_doi('10.1177/002193470003000403')
# SPIE
ftd_doi(c("10.1117/12.59493", "10.1117/12.460027",
  "10.1117/1.jei.27.3.033002"))
# PNAS
ftd_doi(c("10.1073/pnas.93.19.10405", "10.1073/pnas.88.4.1182",
  "10.1073/pnas.87.24.9794"))
# Springer
ftd_doi("10.1007/s10107-017-1136-5")
ftd_doi(c("10.1007/s10107-017-1136-5", "10.1007/978-94-017-8625-6",
  "10.1016/s0952-8733(00)00008-8"))
# American Society of Clinical Oncology
ftd_doi(c("10.1200/JCO.20.01121", "10.1200/JCO.19.02959",
  "10.1200/JCO.20.01002"))
# AIP: American Institute of Physics
ftd_doi(c("10.1063/1.5046187", "10.1063/1.4973652", "10.1063/1.5080806"))
# ACS
ftd_doi(c("10.1021/am508843z", "10.1021/acs.analchem.8b05115",
  "10.1021/acs.jchemed.5b00997"))
# The Royal Society
ftd_doi(c("10.1098/rspa.2007.1849", "10.1098/rstb.1970.0037",
  "10.1098/rsif.2006.0142"))
# Company of Biologists
ftd_doi("10.1242/jeb.00137")
ftd_doi(c("10.1242/dev.00905", "10.1242/dev.00915"))
ftd_doi("10.1242/bio.042192")
# Hindawi
ftd_doi("10.1155/2017/4724852")
ftd_doi("10.1155/2020/6914878")
# IOP
ftd_doi("10.1088/2043-6262/7/2/025018")
# AAAS
# z <- rcrossref::cr_members(221, works=TRUE)
# dois <- z$data$doi
# ftd_doi(dois[12:20])
# ftd_doi(dois[2])
ftd_doi("10.1126/science.276.5312.548")
# Oxford
# z <- rcrossref::cr_members(286, works=TRUE)
# dois <- z$data$doi
# ftd_doi(dois[1:5])
ftd_doi("10.1016/s0895-7061(01)02279-8")
# CDC
# z <- rcrossref::cr_members(1822, works=TRUE)
# dois <- z$data$doi
# ftd_doi(dois[1:5])
ftd_doi("10.3201/eid1209.051606")
# Elsevier
## a cc-by3 paper
ftd_doi(doi="10.1016/j.jsamd.2019.02.002")
ftd_doi(c("10.1016/j.nuclphysbps.2015.09.127", "10.1016/j.nuclphysb.2011.09.011", 
 "10.1016/j.eurpolymj.2018.07.009", "10.1016/j.jsamd.2019.02.002",
 "10.1016/j.physletb.2015.11.072"))
# American Society for Microbiology
ftd_doi(doi="10.1128/jcm.39.12.4344-4348.2001")
ftd_doi(c("10.1128/jcm.42.6.2623-2628.2004",
  "10.1128/jcm.42.9.4147-4153.2004",
  "10.1128/jcm.40.10.3826-3830.2002", 
  "10.1128/jcm.41.3.943-947.2003"))
## some DOIs we just can't easily make URLs for, returns NA
ftd_doi(c("10.1128/mcb.11.10.4966", "10.1128/cmr.7.1.14"))
# Walter de Gruyter
ftd_doi(doi="10.1515/geo-2020-0173")
ftd_doi(doi="10.1515/ci.2013.35.2.19b")
ftd_doi(c("10.1515/geo-2020-0173", "10.1515/ci.2013.35.2.bm", 
  "10.2478/jvetres-2020-0058", "10.2478/acmy-2020-0008"))
# Biorxiv
ftd_doi(doi='10.1101/012476')
}
}
\seealso{
Other ftdoi: 
\code{\link{ftd_fetch_patterns}()},
\code{\link{ftd_members}()},
\code{\link{ftd_prefixes}()},
\code{\link{ftdoi_cache}},
\code{\link{prefix_local}()}
}
\concept{ftdoi}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_providers.R
\name{ft_providers}
\alias{ft_providers}
\title{Search for information on journals or publishers.}
\usage{
ft_providers(journal = NULL, publisher = NULL, limit = 10, ...)
}
\arguments{
\item{journal}{Query terms}

\item{publisher}{Source to query}

\item{limit}{Number of records to return.}

\item{...}{Further args passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
An object of class ft_p
}
\description{
Search for information on journals or publishers.
}
\examples{
\dontrun{
ft_providers(journal="Stem Cells International")
ft_providers(publisher="hindawi")
ft_providers(publisher="journal")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{collect}
\alias{collect}
\title{This function is defunct.}
\usage{
collect(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_get.R
\name{ft_get}
\alias{ft_get}
\alias{ft_get_ls}
\title{Download full text articles}
\usage{
ft_get(
  x,
  from = NULL,
  type = "xml",
  try_unknown = TRUE,
  bmcopts = list(),
  entrezopts = list(),
  elifeopts = list(),
  elsevieropts = list(),
  sciencedirectopts = list(),
  wileyopts = list(),
  crossrefopts = list(),
  progress = FALSE,
  ...
)

ft_get_ls()
}
\arguments{
\item{x}{Either identifiers for papers, either DOIs (or other ids) as a
list of character strings, or a character vector, OR an object of class \code{ft},
as returned from \code{\link[=ft_search]{ft_search()}}}

\item{from}{Source to query. Optional.}

\item{type}{(character) one of xml (default), pdf, or plain (Elsevier
and ScienceDirect only). We choose to go with xml as the default as it has
structure that a machine can reason about, but you are of course free to try
to get xml, pdf, or plain (in the case of Elsevier and ScienceDirect).}

\item{try_unknown}{(logical) if publisher plugin not already known, we try
to fetch full text link either using \pkg{ftdoi} package or from Crossref.
If not found at \pkg{ftdoi} or at Crossref we skip with a warning.
If found with \pkg{ftdoi} or Crossref we attempt to download. Only
applicable in \code{character} and \code{list} S3 methods. Default: \code{TRUE}}

\item{bmcopts}{BMC options. parameter DEPRECATED}

\item{entrezopts}{Entrez options, a named list. See
\code{\link[rentrez:entrez_search]{rentrez::entrez_search()}} and \code{\link[=entrez_fetch]{entrez_fetch()}}}

\item{elifeopts}{eLife options, a named list.}

\item{elsevieropts}{Elsevier options, a named list. Use \code{retain_non_ft=TRUE}
to retain files that do not actually have full text but likely only have an
abstract. By default we set \code{retain_non_ft=FALSE} so that if we detect
that you only got an abstract back, we delete it and report an error
that you likely don't have access.}

\item{sciencedirectopts}{Elsevier ScienceDirect options, a named list.}

\item{wileyopts}{Wiley options, a named list.}

\item{crossrefopts}{Crossref options, a named list.}

\item{progress}{(logical) whether to show progress bar or not.
default: \code{FALSE}. if \code{TRUE}, we use \code{utils::txtProgressBar()} and
\code{utils::setTxtProgressBar()}
to create the progress bar; and each progress bar connection is closed
on function exit. A progress bar is run for each data source.
Works for all S3 methods except \code{ft_get.links}. When articles are not
already downloaded you see the progress bar. If articles are already
downloaded/cached, normally we throw messages saying so, but if a
progress bar is requested, then the messages are suppressed to
not interrupt the progress bar.}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}, see examples below}
}
\value{
An object of class \code{ft_data} (of type \code{S3}) with slots for
each of the publishers. The returned object is split up by publishers
because the full text format is the same within publisher - which should
facilitate text mining downstream as different steps may be needed for
each publisher's content.

Note that we have a print method for \code{ft_data} so you see something
like this:\preformatted{<fulltext text>
[Docs] 4
[Source] ext - /Users/foobar/Library/Caches/R/fulltext
[IDs] 10.2307/1592482 10.2307/1119209 10.1037/11755-024 ...
}

Within each publisher there is a list with the elements:
\itemize{
\item \code{found}: number of full text articles found
\item \code{dois}: the DOIs given and searched for
\item \code{data}
\itemize{
\item \code{backend}: the backend. right now only \code{ext} for "by file extension",
we may add other backends in the future, thus we retain this
\item \code{cache_path}: the base directory path for file caching
\item \code{path}: if file retrieved the full path to the file. if file not
retrived this is \code{NULL}
\item \code{data}: if text extracted (see \code{\link[=ft_collect]{ft_collect()}}) the text will be here,
but until then this is \code{NULL}
}
\item \code{opts}: the options given like article type, dois
\item \code{errors}: data.frame of errors, with two columns for article id and error
}
}
\description{
\code{ft_get} is a one stop shop to fetch full text of articles,
either XML or PDFs. We have specific support for PLOS via the
\pkg{rplos} package, Entrez via the \pkg{rentrez} package, and arXiv via the
\pkg{aRxiv} package. For other publishers, we have helpers to \code{ft_get} to
sort out links for full text based on user input. Articles are saved on
disk. See \code{Details} for help on how to use this function.
}
\details{
There are various ways to use \code{ft_get}:
\itemize{
\item Pass in only DOIs - leave \code{from} parameter \code{NULL}. This route will
first query Crossref API for the publisher of the DOI, then we'll use
the appropriate method to fetch full text from the publisher. If a publisher
is not found for the DOI, then we'll throw back a message telling you a
publisher was not found.
\item Pass in DOIs (or other pub IDs) and use the \code{from} parameter. This route
means we don't have to make an extra API call to Crossref (thus, this route
is faster) to determine the publisher for each DOI. We go straight to
getting full text based on the publisher.
\item Use \code{\link[=ft_search]{ft_search()}} to search for articles. Then pass that output to
this function, which will use info in that object. This behaves the same
as the previous option in that each DOI has publisher info so we know how to
get full text for each DOI.
}

Note that some publishers are available via Entrez, but often not recent
articles, where "recent" may be a few months to a year or so. In that case,
make sure to specify the publisher, or else you'll get back no data.
}
\section{Important Access Notes}{

See \strong{Rate Limits} and \strong{Authentication} in
\link{fulltext-package} for rate limiting and authentication information,
respectively.

In particular, take note that when fetching full text from Wiley, the only
way that's done is through the Crossref Text and Data Mining service. See
the Authenticaiton section of \link{fulltext-package} for all the details.

When fetching articles from Elsevier, the only way that used to be done
was through the Crossref TDM flow. However, Crossref TDM is going away.
See \strong{Authentication} in \link{fulltext-package} for details.
}

\section{Notes on the \code{type} parameter}{

Type is sometimes ignored, sometimes used. For certain data sources,
they only accept one type. By data source/publisher:
\itemize{
\item PLOS: pdf and xml
\item Entrez: only xml
\item eLife: pdf and xml
\item Pensoft: pdf and xml
\item arXiv: only pdf
\item BiorXiv: only pdf
\item Elsevier: xml and plain
\item Elsevier ScienceDirect: xml and plain
\item Wiley: pdf and xml
\item Peerj: pdf and xml
\item Informa: only pdf
\item FrontiersIn: pdf and xml
\item Copernicus: pdf and xml
\item Scientific Societies: only pdf
\item Cambridge: only pdf
\item Crossref: depends on the publisher
\item other data sources/publishers: there are too many to cover here - will
try to make a helper in the future for what is covered by different
publishers
}
}

\section{How data is stored}{

\code{ft_get} used to have many options for "backends". We have simplified this
to one option. That one option is that all full text files are written
to disk on your machine. You can choose where these files are stored.

In addition, files are named by their IDs (usually DOIs), and the file
extension for the full text type (pdf or xml usually). This makes inspecting
the files easy.
}

\section{Data formats}{

xml full text is stored in \code{.xml} files. pdf is stored in \code{.pdf} files.
And plain text is stored in \code{.txt} files.
}

\section{Reusing cached articles}{

All files are written to disk and we check for a file matching the given
DOI/ID on each request - if found we use it and throw message saying so.
}

\section{Caching}{

Previously, you could set caching options in each \code{ft_get} function call.
We've simplified this to only setting caching options through the function
\code{\link[=cache_options_set]{cache_options_set()}} - and you can get your cache options using
\code{\link[=cache_options_get]{cache_options_get()}}. See those docs for help on caching.
}

\section{Notes on specific publishers}{

\itemize{
\item arXiv: The IDs passed are not actually DOIs, though they look
similar. Thus, there's no way to not pass in the \code{from} parameter
as we can't determine unambiguously that the IDs passed in are from
arXiv.org.
\item bmc: Is a hot mess since the Springer acquisition. It's been
removed as an officially supported plugin, some DOIs from them may
still work when passed in here, who knows, it's a mess.
}
}

\section{Warnings}{

You will see warnings thrown in the R shell or in the resulting object.
See \link{ft_get-warnings} for more information on what warnings mean.
}

\examples{
# List publishers included
ft_get_ls()

\dontrun{
# If you just have DOIs and don't know the publisher
## PLOS
ft_get('10.1371/journal.pone.0086169')

# Collect all errors from across papers
#   similarly can combine from different publishers as well
res <- ft_get(c('10.7554/eLife.03032', '10.7554/eLife.aaaa'), from = "elife")
res$elife$errors

## PeerJ
ft_get('10.7717/peerj.228')
ft_get('10.7717/peerj.228', type = "pdf")

## eLife
### xml
ft_get('10.7554/eLife.03032')
res <- ft_get(c('10.7554/eLife.03032', '10.7554/eLife.32763'),
  from = "elife")
res$elife
respdf <- ft_get(c('10.7554/eLife.03032', '10.7554/eLife.32763'),
  from = "elife", type = "pdf")
respdf$elife

elife_xml <- ft_get('10.7554/eLife.03032', from = "elife")
library(magrittr)
elife_xml \%<>\% ft_collect()
elife_xml$elife
### pdf
elife_pdf <- ft_get(c('10.7554/eLife.03032', '10.7554/eLife.32763'),
  from = "elife", type = "pdf")
elife_pdf$elife
elife_pdf \%<>\% ft_collect()
elife_pdf \%>\% ft_extract()

## some BMC DOIs will work, but some may not, who knows
ft_get(c('10.1186/2049-2618-2-7', '10.1186/2193-1801-3-7'), from = "entrez")

## FrontiersIn
res <- ft_get(c('10.3389/fphar.2014.00109', '10.3389/feart.2015.00009'))
res
res$frontiersin

## Hindawi - via Entrez
res <- ft_get(c('10.1155/2014/292109','10.1155/2014/162024',
'10.1155/2014/249309'))
res
res$hindawi
res$hindawi$data$path
res \%>\% ft_collect() \%>\% .$hindawi

## F1000Research - via Entrez
x <- ft_get('10.12688/f1000research.6522.1')
## Two different publishers via Entrez - retains publisher names
res <- ft_get(c('10.1155/2014/292109', '10.12688/f1000research.6522.1'))
res$hindawi
res$f1000research

## Thieme -
### coverage is hit and miss, it's not great
ft_get('10.1055/s-0032-1316462')

## Pensoft
ft_get('10.3897/mycokeys.22.12528')

## Copernicus
out <- ft_get(c('10.5194/angeo-31-2157-2013', '10.5194/bg-12-4577-2015'))
out$copernicus

## arXiv - only pdf, you have to pass in the from parameter
res <- ft_get(x='cond-mat/9309029', from = "arxiv")
res$arxiv
res \%>\% ft_extract  \%>\% .$arxiv

## bioRxiv - only pdf
res <- ft_get(x='10.1101/012476')
res$biorxiv

## AAAS - only pdf
res <- ft_get(x='10.1126/science.276.5312.548')
res$aaas

# The Royal Society
res <- ft_get("10.1098/rspa.2007.1849")
ft_get(c("10.1098/rspa.2007.1849", "10.1098/rstb.1970.0037",
  "10.1098/rsif.2006.0142"))

## Karger Publisher
(x <- ft_get('10.1159/000369331'))
x$karger

## MDPI Publisher
(x <- ft_get('10.3390/nu3010063'))
x$mdpi
ft_get('10.3390/nu7085279')
ft_get(c('10.3390/nu3010063', '10.3390/nu7085279'))

# Scientific Societies
## this is a paywall article, you may not have access or you may
x <- ft_get("10.1094/PHYTO-04-17-0144-R")
x$scientificsocieties

# Informa
x <- ft_get("10.1080/03088839.2014.926032")
ft_get("10.1080/03088839.2013.863435")

## CogentOA - part of Inform/Taylor Francis now
ft_get('10.1080/23311916.2014.938430')

library(rplos)
(dois <- searchplos(q="*:*", fl='id',
   fq=list('doc_type:full',"article_type:\"research article\""),
   limit=5)$data$id)
ft_get(dois)
ft_get(c('10.7717/peerj.228','10.7717/peerj.234'))

# elife
ft_get('10.7554/eLife.04300', from='elife')
ft_get(c('10.7554/eLife.04300', '10.7554/eLife.03032'), from='elife')
## search for elife papers via Entrez
dois <- ft_search("elife[journal]", from = "entrez")
ft_get(dois)

# Frontiers in Pharmacology (publisher: Frontiers)
doi <- '10.3389/fphar.2014.00109'
ft_get(doi, from="entrez")

# Hindawi Journals
ft_get(c('10.1155/2014/292109','10.1155/2014/162024','10.1155/2014/249309'),
  from='entrez')

# Frontiers Publisher - Frontiers in Aging Nueroscience
res <- ft_get("10.3389/fnagi.2014.00130", from='entrez')
res$entrez

# Search entrez, get some DOIs
(res <- ft_search(query='ecology', from='entrez'))
res$entrez$data$doi
ft_get(res$entrez$data$doi[1], from='entrez')
ft_get(res$entrez$data$doi[1:3], from='entrez')

# Search entrez, and pass to ft_get()
(res <- ft_search(query='ecology', from='entrez'))
ft_get(res)

# elsevier, ugh
## set the environment variable Sys.setenv(ELSEVIER_TDM_KEY = "your key")
### an open access article
ft_get(x = "10.1016/j.trac.2016.01.027", from = "elsevier")
### non open access article
#### If you don't have access, by default you get abstract only, and we
##### treat it as an error as we assume you want full text
ft_get(x = "10.1016/j.trac.2016.05.027", from = "elsevier")
#### If you want to retain whatever Elsevier gives you
##### set "retain_non_ft = TRUE"
ft_get(x = "10.1016/j.trac.2016.05.027", from = "elsevier",
  elsevieropts = list(retain_non_ft = TRUE))

# sciencedirect
## set the environment variable Sys.setenv(ELSEVIER_TDM_KEY = "your key")
ft_get(x = "10.1016/S0140-6736(13)62329-6", from = "sciencedirect")

# wiley, ugh
## set the environment variable Sys.setenv(WILEY_TDM_KEY = "your key")
ft_get(x = "10.1006/asle.2001.0035", from = "wiley", type = "pdf")
## xml
ft_get(x = "10.1111/evo.13812", from = "wiley")

## highwire fiasco paper
ft_get(x = "10.3732/ajb.1300053", from = "wiley")
ft_get(x = "10.3732/ajb.1300053", from = "wiley", type = "pdf")

# IEEE, ugh
ft_get('10.1109/TCSVT.2012.2221191', type = "pdf")

# AIP Publishing
ft_get('10.1063/1.4967823', try_unknown = TRUE)

# PNAS
ft_get('10.1073/pnas.1708584115', try_unknown = TRUE)

# American Society for Microbiology
ft_get('10.1128/cvi.00178-17')

# American Society of Clinical Oncology
ft_get('10.1200/JCO.18.00454')

# American Institute of Physics
ft_get('10.1063/1.4895527')

# American Chemical Society
ft_get(c('10.1021/la903074z', '10.1021/jp048806z'))

# Royal Society of Chemistry
ft_get('10.1039/c8cc06410e')


# From ft_links output
## Crossref
(res2 <- ft_search(query = 'ecology', from = 'crossref', limit = 3,
  crossrefopts = list(filter = list(has_full_text=TRUE, member=98))))
(out <- ft_links(res2))
(ress <- ft_get(x = out, type = "pdf"))
ress$crossref

(x <- ft_links("10.1111/2041-210X.12656", "crossref"))
(y <- ft_get(x))

## Cambridge
x <- ft_get("10.1017/s0922156598230305")
x$cambridge
z <- ft_get("10.1017/jmo.2019.20")
z$cambridge
m <- ft_get("10.1017/S0266467419000270")
m$cambridge

## No publisher plugin provided yet
ft_get('10.1037/10740-005')
### no link available for this DOI
res <- ft_get('10.1037/10740-005', try_unknown = TRUE)
res[[1]]

# Get a progress bar - off by default
library(rplos)
(dois <- searchplos(q="*:*", fl='id',
   fq=list('doc_type:full',"article_type:\"research article\""),
   limit=5)$data$id)
## when articles not already downloaded you see the progress bar
b <- ft_get(dois, progress = TRUE)
## if articles already downloaded/cached, normally we through messages
## saying so
b <- ft_get(dois, progress = FALSE)
## but if a progress bar is requested, then the messages are suppressed
b <- ft_get(dois, progress = TRUE)

# curl options
ft_get("10.1371/journal.pcbi.1002487", verbose = TRUE)
ft_get('10.3897/mycokeys.22.12528', from = "pensoft", verbose = TRUE)
}
}
\seealso{
\code{\link[=as.ft_data]{as.ft_data()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{chunks}
\alias{chunks}
\title{This function is defunct.}
\usage{
chunks(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chunks.R
\name{ft_chunks}
\alias{ft_chunks}
\title{Extract chunks of data from articles}
\usage{
ft_chunks(...)
}
\description{
Extract chunks of data from articles
}
\keyword{internal}
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
% Please edit documentation in R/scopus_utils.R
\name{scopus_search}
\alias{scopus_search}
\alias{scopus_search_loop}
\title{Scopus search}
\usage{
scopus_search(
  query = NULL,
  count = 25,
  start = 0,
  type = "search",
  search_type = "scopus",
  facets = NULL,
  view = NULL,
  date = NULL,
  sort = NULL,
  content = NULL,
  subj = NULL,
  key = NULL,
  ...
)

scopus_search_loop(
  query = NULL,
  count = 25,
  start = 0,
  type = "search",
  search_type = "scopus",
  facets = NULL,
  view = NULL,
  date = NULL,
  sort = NULL,
  content = NULL,
  subj = NULL,
  key = NULL,
  ...
)
}
\arguments{
\item{query}{(character) query terms, as a single character vector}

\item{count}{(integer/numeric) results to return: default: 25}

\item{start}{(integer/numeric) offset value, default: 0}

\item{type}{(character) type of search, default: search}

\item{search_type}{(character) search type, default: scopus}

\item{facets}{(list) facets, see
https://dev.elsevier.com/tecdoc_api_facets.html for how to construct
facet queries}

\item{view}{the fields to return, see
https://dev.elsevier.com/guides/ScopusSearchViews.htm}

\item{date}{Represents the date range associated with the search,
with the lowest granularity being year. e.g. 2002-2007}

\item{sort}{Represents the sort field name and order. A plus in front of
the sort field name indicates ascending order, a minus indicates
descending order. If sort order is not specified (i.e. no + or -) then
the order defaults to ascending (ASC). Up to three fields can be
specified, each delimited by a comma. The precedence is determined by
their order (i.e. first is primary, second is secondary, and
third is tertiary). . e.g., "overDate,-title"}

\item{content}{filter specific categories of content that should be
searched/returned. one of: core, dummy, all (default)}

\item{subj}{the subject area code associated with the content category
desired. Note that these subject code mapping vary based upon the
environment in which the request is executed. See Details for choices.}

\item{key}{(character) api key. get a key at
https://dev.elsevier.com/index.html}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Scopus search
}
\details{
Rate limits for search are 20,000 per every 7 days. You likely
won't make that many requests in 7 days, but if you do e.g., make 20K in
5 days, then you have to wait 2 days for the clock to reset, than you'll
be able to make 20K again.

See https://dev.elsevier.com/api_key_settings.html for rate
limit information.

See https://dev.elsevier.com/sc_search_tips.html for help/tips
on searching
}
\section{subj choices include}{

\itemize{
\item AGRI: Agricultural and Biological Sciences
\item ARTS: Arts and Humanities
\item BIOC: Biochemistry, Genetics and Molecular Biology
\item BUSI: Business, Management and Accounting
\item CENG: Chemical Engineering
\item CHEM: Chemistry
\item COMP: Computer Science
\item DECI: Decision Sciences
\item DENT: Dentistry
\item EART: Earth and Planetary Sciences
\item ECON: Economics, Econometrics and Finance
\item ENER: Energy
\item ENGI: Engineering
\item ENVI: Environmental Science
\item HEAL: Health Professions
\item IMMU: Immunology and Microbiology
\item MATE: Materials Science
\item MATH: Mathematics
\item MEDI: Medicine
\item NEUR: Neuroscience
\item NURS: Nursing
\item PHAR: Pharmacology, Toxicology and Pharmaceutics
\item PHYS: Physics and Astronomy
\item PSYC: Psychology
\item SOCI: Social Sciences
\item VETE: Veterinary
\item MULT: Multidisciplinary
}
}

\examples{
\dontrun{
res <- scopus_search(query = "ecology")
res

#scopus_search(query = x, type = "abstract")

# looping through
res <- scopus_search_loop(query = "ecology community elk cow")

# using facets
## scopus_search
# res <- scopus_search(query = "ecology", facets = "subjarea(count=5)")
# res
# res$`search-results`$link
# res$`search-results`$entry
# res$`search-results`$facet

## more examples
# x <- scopus_search(query = "ecology", facets = "language(count=4)",
#   count = 1)
# x$`search-results`$facet
# x <- scopus_search(query = "ecology",
#   facets = "pubyear(count=3);doctype();language(count=4)")
# x$`search-results`$facet

## scopus_search_loop
# res <- scopus_search_loop(query = "ecology", facets = "subjarea(count=5)",
#   count = 200)
# res$found
# head(res$results)
# NROW(res$results)
# res$facets

# sort
x <- scopus_search(query = "ecology", sort = "-title")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fulltext-package.R
\docType{package}
\name{fulltext-package}
\alias{fulltext-package}
\alias{fulltext}
\title{Fulltext search and retrieval of scholarly texts.}
\description{
fulltext is a single interface to many sources of scholarly
texts. In practice, this means only ones that are legally useable. We will
support sources that require authentication on a case by case basis - that
is, if more than just a few people will use it, and it's not too
burdensome to include, then we can include that source.
}
\section{Manual}{

See https://books.ropensci.org/fulltext/ for a longer form
manual for using \pkg{fulltext}.
}

\section{What's included}{

We currently include support for search and full text retrieval for a variety
of publishers. See \code{\link[=ft_search]{ft_search()}} for what we include for search,
and  \code{\link[=ft_get]{ft_get()}} for what we include for full text retrieval.
}

\section{Use cases}{

The following are tasks/use cases supported:
\itemize{
\item search - \code{\link[=ft_search]{ft_search()}}
\item get texts - \code{\link[=ft_get]{ft_get()}}
\item get full text links - \code{\link[=ft_links]{ft_links()}}
\item get abstracts - \code{\link[=ft_abstract]{ft_abstract()}}
\item extract text from pdfs - \code{\link[=ft_extract]{ft_extract()}}
\item serialize to different data formats - \code{\link[=ft_serialize]{ft_serialize()}}
\item extract certain article sections (e.g., authors) - moved to \pkg{pubchunks}
}
}

\section{DOI delays}{

Beware that DOIs are not searchable via Crossref/Entrez immediately. The
delay may be as much as a few days, though should be less than a day. This
delay should become shorter as services improve. The point of this is that
you man not find a match for a relatively new DOI (e.g., for an article
published the same day). We've tried to account for this for some publishers.
For example, for Crossref we search Crossref for a match for a DOI, and if
none is found we attempt to retrieve the full text from the publisher
directly.
}

\section{Rate limits}{

\strong{Scopus}: 20,000 per 7 days. See
https://dev.elsevier.com/api_key_settings.html for rate
limit information. To see what your personal rate limit details are,
request verbose HTTP request output - this will vary on the function
you are using - see the docs for the function. See the response
headers \code{X-RateLimit-Limit}, \code{X-RateLimit-Remaining}, and
\code{X-RateLimit-Reset} (your limit, those requests remaining, and UTC
date/time it will reset)

\strong{Microsoft}: 10,000 per month, and 1 per second. There are no rate
limit headers, sorry :(

\strong{PLOS}: There are no known rate limits for PLOS, though if you do
hit something let us know.

\strong{Crossref}: From time to time Crossref needs to impose rate limits
to ensure that the free API is usable by all. Any rate limits that are in
effect will be advertised in the \code{X-Rate-Limit-Limit} and
\code{X-Rate-Limit-Interval} HTTP headers. This boils down to: they allow X
number of requests per some time period. The numbers can change so we
can't give a rate limit that will always be in effect. If you're curious
pass in \code{verbose = TRUE} to your function call, and you'll get headers
that will display these rate limits. See also \strong{Authentication}.

\strong{Semantic Scholar}: Not documented in their docs, and no response
headers given. At time of this writing (2020-07-01) the rate limit is:
100 requests per 5-minutes per IP address. or 20 requests per min. Note
that this rate limit may change.
}

\section{Authentication}{


\strong{BMC}: BMC is integrated into Springer Publishers now,
and that API requires an API key.  Get your key by signing up at
https://dev.springer.com/, then you'll get a key. Pass the key to a
named parameter \code{key} to \code{bmcopts}. Or, save your key in your \code{.Renviron}
file as \code{SPRINGER_KEY}, and we'll read it in for you, and you don't
have to pass in anything.

\strong{Scopus}: Scopus requires two things: an API key and your institution must
have access. For the API key, go to https://dev.elsevier.com/index.html,
register for an account, then when you're in your account, create an API key.
Pass in as variable \code{key} to \code{scopusopts}, or store your key under the name
\code{ELSEVIER_SCOPUS_KEY} as an environment variable in \code{.Renviron}, and
we'll read it in for you. See \link{Startup} for help. For the institution access
go to a browser and see if you have access to the journal(s) you want.
If you don't have access in a browser you probably won't have access via
this package. If you aren't physically at your institution you will likely
need to be on a VPN or similar and eventually require correct proxy settings,
so that your IP address is in the range
that the two publishers are accepting for that institution.
It might be, that the API access seems to work even while
in the wrong IP range or have wrong proxy settings,
but you are not able to see the abstracts, they will be empty.
By using the currect curl options into the calls to ft_search or ft_abstracts even
the most complex proxy including authentication should work. As an example:

\preformatted{
opts <- list(key="your-scopus-key")
ft_abstract(x = dois, from = "scopus", scopusopts = opts,
  proxy="proxy-ip-address",
  proxyport=your-proxy-port,
  proxyuserpwd="username:password", # often the same as your windows login
  proxyauth=8) # ntlm - authentication
}

\strong{Elsevier/ScienceDirect}: Elsevier ScienceDirect requires two things: an
API key and your institution must have access. For the API key,
go to https://dev.elsevier.com/index.html,
register for an account, then when you're in your account, create an API key
that is allowed to access the TDM API (must accept their TDM policy).
Pass in as variable \code{key} to \code{elsevieropts}/\code{sciencedirectopts}, or store
your key under the name \code{ELSEVIER_TDM_KEY} as an environment variable in
\code{.Renviron}, and we'll read it in for you. See \link{Startup} for help. For
the institution access
go to a browser and see if you have access to the journal(s) you want.
If you don't have access in a browser you probably won't have access via
this package. If you aren't physically at your institution you will likely
need to be on a VPN or similar so that your IP address is in the range
that the publisher is accepting for that institution.

\strong{Wiley}: Replacing Crossref TDM service as of February 2021, Wiley
now requires you get a Wiley TDM key. Get one at
https://onlinelibrary.wiley.com/library-info/resources/text-and-datamining
Pass in as variable \code{key} to \code{wileyopts}, or preferably store
your key under the name \code{WILEY_TDM_KEY} as an environment variable in
\code{.Renviron}, and we'll read it in for you. See \link{Startup} for help. Some
notes about Wiley's TDM service:
\itemize{
\item They always respond initially with a redirect to a server dedicated to
the serving of binary resources - fulltext takes care of this
\item Wiley uses rate-limiting: no more than 3 requests per second. you
may get 429 errors if making too requests too rapidly
}

\strong{Microsoft}: Get a key by creating an Azure account
then request a key for \strong{Academic Knowledge API} within
\strong{Cognitive Services}. Store it as an environment variable in your
\code{.Renviron} file - see \link{Startup} for help. Pass your
API key into \code{maopts} as a named element in a list like
\code{list(key = Sys.getenv('MICROSOFT_ACADEMIC_KEY'))}

\strong{Crossref}: Crossref encourages requests with contact information
(an email address) and will forward you to a dedicated API cluster
for improved performance when you share your email address with them. This
is called the "Polite Pool".
https://github.com/CrossRef/rest-api-doc#good-manners--more-reliable-service
To pass your email address to Crossref via this client, store it
as an environment variable in \code{.Renviron} like
\code{crossref_email=name@example.com}, or \code{CROSSREF_EMAIL=name@example.com}.
Save the file and restart your R session. To stop sharing your email when
using rcrossref simply delete it from your \code{.Renviron} file OR to temporarily
not use your email unset it for the session
like \code{Sys.unsetenv('crossref_email')}. To be sure your in the polite pool
use curl verbose by e.g., \code{ft_cr_links(doi = "10.5555/515151", verbose = TRUE)}

\strong{Crossref TDM}: TDM = "Text and Data Mining". This used to apply to just
two publishers - Wiley and Elsevier - This service officially shut down at
the end of 2020. For Elsevier, see the
"Elsevier/ScienceDirect" section above. For Wiley, see the "Wiley" section
above.

\strong{Entrez}: NCBI limits users to making only 3 requests per second. But, users
who register for an API key are able to make up to ten requests per second.
Getting a key is simple; register for a "my ncbi" account then click on a
button in the account settings page. Once you have an API key, you can pass it
as the argument \code{api_key} to \code{entrezopts} in both \code{\link[=ft_get]{ft_get()}} and \code{\link[=ft_search]{ft_search()}}.
However, we advise you use environment variables instead as they are more secure.
To do that you can set an environment variable for the current R session like
\code{Sys.setenv(ENTREZ_KEY="yourkey")} OR better yet set it in your \code{.Renviron}
or equivalent file with an entry like \code{ENTREZ_KEY=yourkey} so that it is
used across R sessions.

No authentication needed for \strong{PLOS}, \strong{eLife}, \strong{arxiv}, \strong{biorxiv},
\strong{Euro PMC}

Open an issue if you run into trouble with authentication.
}

\section{Feedback}{

Let us know what you think at https://github.com/ropensci/fulltext/issues
}

\author{
Scott Chamberlain

Will Pearse

Helge Knüttel
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{fulltext-defunct}
\alias{fulltext-defunct}
\title{Defunct functions in fulltext}
\description{
\itemize{
\item \link{ft_extract_corpus} Function removed. As part of focusing scope of the
package we're trying to limit dependencies, so downstream use of \code{tm} can
still easily be done.
\item \link{pdfx}: Function removed. We're trying to focus the scope of the
package - and this function is more out of scope now.
\item \link{chunks}: Function name changed to \code{\link[=ft_chunks]{ft_chunks()}}
\item \link{tabularize}: Function name changed to \code{\link[=ft_tabularize]{ft_tabularize()}}
\item \link{collect}: Function name changed to \code{\link[=ft_collect]{ft_collect()}}
\item \link{get_text}: Function name changed to \code{\link[=ft_text]{ft_text()}}
\item \code{cache_clear} was never working anyway, and is now removed
\item \link{ft_browse_sections}: no sign that function used, and allows
to remove a dependency
\item \link{ft_get_si}: moved to package \code{suppdata}
\item \link{ft_chunks}: moved to package \code{pubchunks}
\item \link{ft_tabularize}: moved to package \code{pubchunks}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftdoi_cache.R
\name{ftdoi_cache}
\alias{ftdoi_cache}
\title{Caching}
\description{
Manage cached \code{ftdoi} files with \pkg{hoardr}
}
\details{
The dafault cache directory is
\code{paste0(rappdirs::user_cache_dir(), "/R/ftdoi")}, but you can set
your own path using \code{cache_path_set()}

\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}
\section{Useful user functions}{

\itemize{
\item \code{ftdoi_cache$cache_path_get()} get cache path
\item \code{ftdoi_cache$cache_path_set()} set cache path. You can set the
entire path directly via the \code{full_path} arg like
\code{ftdoi_cache$cache_path_set(full_path = "your/path")}
\item \code{ftdoi_cache$list()} returns a character vector of full
path file names
\item \code{ftdoi_cache$files()} returns file objects with metadata
\item \code{ftdoi_cache$details()} returns files with details
\item \code{ftdoi_cache$delete()} delete specific files
\item \code{ftdoi_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
ftdoi_cache

# list files in cache
ftdoi_cache$list()

# delete certain database files
# ftdoi_cache$delete("file path")
# ftdoi_cache$list()

# delete all files in cache
# ftdoi_cache$delete_all()
# ftdoi_cache$list()
}
}
\seealso{
Other ftdoi: 
\code{\link{ftd_doi}()},
\code{\link{ftd_fetch_patterns}()},
\code{\link{ftd_members}()},
\code{\link{ftd_prefixes}()},
\code{\link{prefix_local}()}
}
\concept{ftdoi}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{ft_extract_corpus}
\alias{ft_extract_corpus}
\title{This function is defunct.}
\usage{
ft_extract_corpus(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bio_utils.R
\name{biorxiv_search}
\alias{biorxiv_search}
\title{Biorxiv search}
\usage{
biorxiv_search(query, limit = 10, date_from = NULL, date_to = NULL, ...)
}
\arguments{
\item{query}{query terms}

\item{limit}{records to return. default: 10}

\item{date_from, date_to}{date begin and end, of form YYYY-MM-DD}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Biorxiv search
}
\details{
We search Biorxiv first, get DOIs, then search Crossref -
one consequence of this is that you may get back less than the number of
results you requested even if Biorxiv found equal to or more than
the amount you requested - BECAUSE we take the DOIs from the results and
go out to Crossref to get what we think is better metadata than what
Biorxiv has.
}
\examples{
\dontrun{
biorxiv_search(query = "owls")
biorxiv_search(query = "owls", date_from = "2016-01-01", 
  date_to = "2016-12-30", limit = 10)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftdoi_patterns.R
\name{ftd_fetch_patterns}
\alias{ftd_fetch_patterns}
\title{Download patterns files}
\usage{
ftd_fetch_patterns()
}
\value{
character vector of file paths
}
\description{
Does various checks to see if patterns files alrady downloaded,
out of date, if some/all are deleted and in need of an update
}
\seealso{
Other ftdoi: 
\code{\link{ftd_doi}()},
\code{\link{ftd_members}()},
\code{\link{ftd_prefixes}()},
\code{\link{ftdoi_cache}},
\code{\link{prefix_local}()}
}
\concept{ftdoi}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{tabularize}
\alias{tabularize}
\title{This function is defunct.}
\usage{
tabularize(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftdoi_other.R
\name{ftd_members}
\alias{ftd_members}
\title{Members}
\usage{
ftd_members(id = NULL)
}
\arguments{
\item{id}{(character) a Crossref member ID. Default is \code{NULL} which
gets all members}
}
\description{
Members
}
\examples{
\dontrun{
ftd_members()
ftd_members(221)
ftd_members(1965)

# not found
# ftd_members(999)
}
}
\seealso{
Other ftdoi: 
\code{\link{ftd_doi}()},
\code{\link{ftd_fetch_patterns}()},
\code{\link{ftd_prefixes}()},
\code{\link{ftdoi_cache}},
\code{\link{prefix_local}()}
}
\concept{ftdoi}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftxt_cache.R
\name{ftxt_cache}
\alias{ftxt_cache}
\title{Inspect and manage cached files}
\description{
Inspect and manage cached files
}
\section{Useful user functions for managing cached files}{

\itemize{
\item \code{ftxt_cache$list()} returns a character vector of full
path file names
\item \code{ftxt_cache$files()} returns file objects with metadata
\item \code{ftxt_cache$details()} returns files with details
\item \code{ftxt_cache$delete()} delete specific files
\item \code{ftxt_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
ftxt_cache

# list files in cache
ftxt_cache$list()

# list details of files in cache
ftxt_cache$details()

# delete certain database files
# ftxt_cache$delete("file path")
# ftxt_cache$list()

# delete all files in cache
# ftxt_cache$delete_all()
# ftxt_cache$list()
}
}
\seealso{
\link{cache}, \code{\link[=cache_file_info]{cache_file_info()}}

Other caching-functions: 
\code{\link{cache_file_info}()},
\code{\link{cache}}
}
\concept{caching-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_cr_links.R
\name{ft_cr_links}
\alias{ft_cr_links}
\title{Get Crossref full text links from a DOI}
\usage{
ft_cr_links(doi, type = "all", ...)
}
\arguments{
\item{doi}{(character) A Digital Object Identifier (DOI). required.}

\item{type}{(character) One of 'xml', 'html', 'plain', 'pdf',
'unspecified', or 'all' (default). required.}

\item{...}{Named parameters passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
\code{NULL} if no full text links given; a list of tdmurl objects if
links found. a tdmurl object is an S3 class wrapped around a simple list,
with attributes for:
\itemize{
\item type: type, matchin type passed to the function
\item doi: DOI
\item member: Crossref member ID
\item intended_application: intended application, e.g., text-mining
}
}
\description{
Get Crossref full text links from a DOI
}
\details{
Note that this function is not vectorized.

Some links returned will not in fact lead you to full text
content as you would understandbly think and expect. That is, if you
use the \code{filter} parameter with e.g., \code{\link[rcrossref:cr_works]{rcrossref::cr_works()}}
and filter to only full text content, some links may actually give back
only metadata for an article. Elsevier is perhaps the worst offender,
for one because they have a lot of entries in Crossref TDM, but most
of the links that are apparently full text are not in fact full text,
but only metadata. You can get full text if you are part of a subscribing
institution to that specific Elsever content, but otherwise, you're SOL.

Note that there are still some bugs in the data returned form CrossRef.
For example, for the publisher eLife, they return a single URL with
content-type application/pdf, but the URL is not for a PDF, but for both
XML and PDF, and content-type can be set with that URL as either XML or
PDF to get that type.

In another example, all Elsevier URLs at time of writing are have
\code{http} scheme, while those don't actually work, so we have a
custom fix in this function for that publisher. Anyway, expect changes...
}
\section{Register for the Polite Pool}{

See of 'Authentication' setion of the \link{fulltext-package} manual page
}

\examples{
\dontrun{
dois <- c("10.1245/s10434-016-5211-6",
"10.17159/2413-3108/2016/v0i55a49", "10.17159/2413-3108/2015/v0i53a455",
"10.17159/2413-3108/2006/v0i18a982", "10.1007/s10665-016-9845-y", 
"10.1016/j.ad.2015.06.020", "10.1016/j.medipa.2014.03.002")

# pdf link
ft_cr_links(doi = "10.5555/515151", "pdf")

# xml and plain text links
ft_cr_links(dois[1], "pdf")
ft_cr_links(dois[6], "xml")
ft_cr_links(dois[7], "plain")
ft_cr_links(dois[1]) # all is the default

# pdf link
ft_cr_links(doi = "10.5555/515151", "pdf")
ft_cr_links(doi = "10.3897/phytokeys.52.5250", "pdf")

# many calls, use e.g., lapply
lapply(dois[1:3], ft_cr_links)

# elsevier
## DOI that is open acccess
ft_cr_links('10.1016/j.physletb.2010.10.049')
## DOI that is not open acccess
ft_cr_links('10.1006/jeth.1993.1066')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftdoi_other.R
\name{ftd_prefixes}
\alias{ftd_prefixes}
\title{Prefixes}
\usage{
ftd_prefixes(id = NULL)
}
\arguments{
\item{id}{(character) a DOI prefix. Default is \code{NULL}, which
gets all}
}
\value{
named list of details of the publisher for the DOI prefix
}
\description{
Prefixes
}
\examples{
\dontrun{
ftd_prefixes()
ftd_prefixes(id = '10.1080')

# doesn't work
# ftd_prefixes(id = '10.9999')
}
}
\seealso{
Other ftdoi: 
\code{\link{ftd_doi}()},
\code{\link{ftd_fetch_patterns}()},
\code{\link{ftd_members}()},
\code{\link{ftdoi_cache}},
\code{\link{prefix_local}()}
}
\concept{ftdoi}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_serialize.R
\name{ft_serialize}
\alias{ft_serialize}
\alias{ft_get_keys}
\title{Serialize raw text to other formats, including to disk}
\usage{
ft_serialize(x, to = "xml", from = NULL, ...)

ft_get_keys(x)
}
\arguments{
\item{x}{Input object, output from a call to \code{ft_get}. Required.}

\item{to}{(character) Format to serialize to. One of list,
xml, or json. Required. Output to xml returns object of
class XMLInternalDocument.}

\item{from}{(character) Format \code{x} is currently in. Function attempts
to use metadata provided, or guess from data itself. Optional.
CURRENTLY IGNORED.}

\item{...}{Further args passed on to \code{xml2::read_xml()} or
\code{jsonlite::toJSON()}}
}
\value{
An object of class \code{ft_parsed}
}
\description{
\code{ft_serialize} helps you convert to various data formats. If
your data is in unparsed XML (i.e., character class), you can convert to
parsed XML. If in XML, you can convert to (ugly-ish) JSON, or a list.
}
\examples{
\dontrun{
res <- ft_get('10.7717/peerj.228')

# if articles in xml format, parse the XML
(out <- ft_serialize(ft_collect(res), to='xml'))
out$peerj$data$data[[1]] # the xml

# From XML to JSON
(out <- ft_serialize(ft_collect(res), to='json'))
out$peerj$data$data$`10.7717/peerj.228` # the json
jsonlite::fromJSON(out$peerj$data$data$`10.7717/peerj.228`)

# To a list
out <- ft_serialize(ft_collect(res), to='list')
out$peerj$data$data
out$peerj$data$data[[1]]$body$sec$title
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collect.R
\name{ft_collect}
\alias{ft_collect}
\alias{ft_text}
\alias{ft_text.default}
\alias{ft_text.ft_data}
\title{Collect article text from local files}
\usage{
ft_collect(x, ...)

ft_text(x, ...)

\method{ft_text}{default}(x, ...)

\method{ft_text}{ft_data}(x, ...)
}
\arguments{
\item{x}{Input. An object of class \code{ft_data}}

\item{...}{Further args, ignored.}
}
\value{
an object of class \code{ft_data}, but the \code{data} slot should have
character string of text from the XML/plain text/PDF file
}
\description{
\code{ft_collect} grabs full text data from file paths in your
\code{ft_data} object (result of call to \code{ft_get()}). \code{ft_text} is a
convenience function to grab the nested text data and bring it up in
the list for easier access
}
\details{
The result of this call is actual text you can read
}
\examples{
\dontrun{
# Get some data
x <- ft_get('10.1371/journal.pone.0086169')

# note that the data is not in the object, gives NULL
x$plos$data$data

# Collect data from the .xml file
y <- x \%>\% ft_collect()

# note how the data is now in the object
y$plos$data$data

# Let's get the actual 
## ft_collect() alone, replaces file pointers with parsed text, 
##  maintaining object structure
x \%>\% ft_collect() 
## pulls the text out of the object
x \%>\% ft_collect() \%>\% ft_text()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{pdfx}
\alias{pdfx}
\title{This function is defunct.}
\usage{
pdfx(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ft_abstract.R
\name{ft_abstract}
\alias{ft_abstract}
\alias{ft_abstract_ls}
\title{Get abstracts}
\usage{
ft_abstract(
  x,
  from = "crossref",
  plosopts = list(),
  scopusopts = list(),
  maopts = list(),
  crossrefopts = list(),
  ...
)

ft_abstract_ls()
}
\arguments{
\item{x}{(character) DOIs as a character vector. See Details.}

\item{from}{Source to query. One or more of crossref (default), plos,
scopus, microsoft, or semanticscholar}

\item{plosopts}{PLOS options, a named list.}

\item{scopusopts}{Scopus options, a named list.}

\item{maopts}{Microsoft Academic options, a named list.}

\item{crossrefopts}{Crossref options, a named list.}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}, see
examples below}
}
\value{
An object of class \code{ft_abstract}
}
\description{
Get abstracts
}
\details{
See \strong{Rate Limits} and \strong{Authentication} in
\link{fulltext-package} for rate limiting and authentication information,
respectively. In particular take int account Semantic Scholar rate limits
because we do asynchronous requests to Semantic Scholar, which means
you can get data fast, but you'll hit your rate limit fast too.

There's no options to pass on when \code{from="semanticscholar"}, other than
curl options via \code{...}

When \code{from="semanticscholar"}, ids passed to \code{x} can be various types: DOI,
S2 paper id (Semantic Scholar id), arXiv id, MAG id, ACL id, PubMed id,
or Corpus id. If you use DOIs or S2 paper ids you can pass them to \code{x}
as is. However, if you use other id types you need to prefix each id
with the name of the type of id, options are: "arXiv", "MAG", "ACL",
"PMID", "CorpusID"
}
\examples{
# List publishers included
ft_abstract_ls()

\dontrun{
# PLOS
## search
(res <- ft_search(query = 'biology', from = 'plos', limit = 25,
   plosopts = list(fq = list('doc_type:full', '-article_type:correction',
                  '-article_type:viewpoints'))))
## get abstracts
dois <- res$plos$data$id
(out <- ft_abstract(x = dois, from = "plos"))
out$plos

# Semantic Scholar
(out <- ft_abstract(x = dois, from = "semanticscholar"))
out$semanticscholar
## using arxiv ids
arxiv_ids <- c("0710.3491", "0804.0713", "0810.4821", "1003.0315")
(out <- ft_abstract(x = paste0("arXiv:", arxiv_ids), from = "semanticscholar"))
out$semanticscholar

# Scopus
opts <- list(key = Sys.getenv('ELSEVIER_SCOPUS_KEY'))

## search
(res <- ft_search(query = 'biology', from = 'scopus', scopusopts = opts,
  limit = 25))
## get abstract
dois <- na.omit(res$scopus$data$`prism:doi`)
out <- ft_abstract(x = dois[1:3], from = "scopus", scopusopts = opts)
out
out$scopus

# use scopus Ids
(res <- ft_search(query = 'biology', from = 'scopus', scopusopts = opts,
  limit = 50))
ids <- fulltext:::strextract(res$scopus$data$`dc:identifier`, "[0-9]+")
(out <- ft_abstract(x = ids[1:4], from = 'scopus',
  scopusopts = list(
    key = Sys.getenv('ELSEVIER_SCOPUS_KEY'),
    id_type = "scopus_id"
  )
))

# Microsoft
key <- Sys.getenv("MICROSOFT_ACADEMIC_KEY")
(res <- ft_search("Y=[2010, 2012)", from = "microsoft",
  maopts = list(key = key)))
ids <- res$ma$data$Id
(out <- ft_abstract(x = ids, from = "microsoft",
  maopts = list(
    key = Sys.getenv('MICROSOFT_ACADEMIC_KEY')
  )
))
out$ma
cat(unlist(lapply(out$ma, "[[", "abstract")), sep = "\n\n")

# Crossref
(res <- ft_search("ecology", from = "crossref",
  crossrefopts = list(filter = c(has_abstract = TRUE))))
ids <- res$crossref$data$doi
(out <- ft_abstract(x = ids, from = "crossref"))
out$crossref

# curl options
ft_abstract("10.2458/v17i1.21696", from = "crossref", verbose = TRUE)
ft_abstract("10.1371/journal.pcbi.1002487", from = "plos", verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fulltext-warnings.R
\name{ft_get-warnings}
\alias{ft_get-warnings}
\title{fulltext warnings details}
\description{
What can you do about the various warnings?
}
\details{
This document is in relation to the function \code{\link[=ft_get]{ft_get()}}
}
\section{No plugin}{


For the warning "no plugin for Crossref ...", this is
what happened internally:

This happens when we don't have a hard coded
plugin for that specific publisher within this packge
(use \code{ft_get_ls()} to see what hard coded publisher plugins
we have), but we do have generic functions for Crossref and
ftdoi that are also tried and may get a result. You
are welcome to open up an issue at
https://github.com/ropensci/fulltext/issues to discuss
publisher specific plugins.
}

\section{Access or an error}{


For the warning "you may not have access to ... or an error occurred"
we've likely tried to get the full text but either an error
occurred (which can be a lot of things), or you don't have access
to the full text.

If you think the problem may be that you don't have access,
check whether you are on an IP address that has access to the
full text, and if you're not, get on one that does - most
likely by being on campus/etc. or through a VPN.
}

\section{Part of an article}{


For the warning "... was not found or may be a DOI for a part of an article"
this happens for certain publishers (e.g., PLOS) that issue DOIs for
parts of articles (e.g., abstract, body, supplements) - in which case it
doesn't make sense to get full text of, for example, supplements.
}

\section{No Crossref link}{


For the warning "no link found from Crossref", this happens when
we've gone through the route of searching for a full text URL from
the Crossref API, and there wasn't one, so we stop searching and
give that warning.
}

\section{CROSSREF_TDM env var detected}{


The Crossref Text and Data Mining service ended at the end of 2020.
We only check for the \code{CROSSREF_TDM} environment variable when you request
articles from Elsevier and Wiley as those two publishers were the only
major publishers that used the Crossref TDM service. If we detect the key
we'll throw a warning - you should probably remove your \code{CROSSREF_TDM}
key as the service is no longer operational.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/europe_pmc_utils.R
\name{eupmc}
\alias{eupmc}
\alias{eupmc_search}
\alias{eupmc_fields}
\alias{eupmc_xml}
\title{Europe PMC utilities}
\usage{
eupmc_search(
  query,
  resulttype = "lite",
  synonym = FALSE,
  per_page = 25,
  cursorMark = "*",
  ...
)

eupmc_fields(...)

eupmc_xml(id, ...)
}
\arguments{
\item{query}{(character) Search terms. Required. See Details}

\item{resulttype}{(character) The result type can either be idlist,
core or lite. This parameter determines the fields returned by XML and
JSON formats, but it has no effect on the DC format. See Details.}

\item{synonym}{(boolean) Synonym searches are not made by default
(Default: \code{FALSE}), however queries can be expanded using MeSH
terminology and the UniProt synonyms list.  For example aspirin, a synonym
would be acetylsalicylic acid;this could be included in
the search by setting the parameter value to \code{TRUE}}

\item{per_page}{(integer) Number of records to return. Max: 1000.
Default: 25}

\item{cursorMark}{(character) cursor string, default: \code{*}}

\item{...}{curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{id}{A single Europe PMC article identifier, begins with "PMC",
followed by numbers, e.g.,  "PMC3257301"}

\item{sort}{(character) The default sort order is relevance. Specify the
sort field and sort order. This parameter provides "asc" or "desc" order
for every single-valued field: P_PDATE_D, AUTH_FIRST, CITED etc. For
example sorting by CITED in ascending order: CITED asc}
}
\value{
\code{eupmc_search} returns a list with results. \code{eupmc_fields} returns
a data.frame. \code{eupmc_xml} returns an object of class \code{xml_document}
}
\description{
Europe PMC utilities
}
\section{\code{query} parameter options}{

\itemize{
\item a keyword or combination of keywords (e.g. HPV virus).
\item a phrase with enclosing speech marks (e.g. "human malaria").
\item a fielded search (e.g. auth:stoehr). Available search fields are listed
in the Appendix 1 of the Reference Guide or can be retrieved using the
fields module of the API.
\item a specific publication (e.g. ext_id:781840 src:med) Specify ext_id as the
article identifier, and src as the source database. List of the data sources
can be found on the help pages or in section 3 of the Reference Guide.
}
}

\section{\code{resulttype} parameter options}{

\itemize{
\item idlist - returns a list of IDs and sources for the given search terms
\item lite - returns key metadata for the given search terms; this is the
default value if the parameter is unspecified.
\item core - returns full metadata for a given publication ID; including
abstract, full text links, and MeSH terms.
}
}

\examples{
\dontrun{
# search
eupmc_search(query = 'ecology')
eupmc_search(query = 'human malaria')
eupmc_search(query = '"human malaria"')
eupmc_search(query = 'auth:stoehr')
eupmc_search(query = 'journal:pnas')
eupmc_search(query = 'journal:pnas')
eupmc_search(query = 'ext_id:781840 src:med')
eupmc_search(query = 'ext_id:IND43783977 src:agr')

# list indexed search fields
x <- eupmc_fields()
NROW(x)
head(x)

# get full text XML
eupmc_xml('PMC3257301')
}
}
\references{
https://europepmc.org/RestfulWebService
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{ft_get_si}
\alias{ft_get_si}
\title{This function is defunct.}
\usage{
ft_get_si(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{cache}
\alias{cache}
\alias{cache_options_set}
\alias{cache_options_get}
\title{Set or get cache options}
\usage{
cache_options_set(
  path = "fulltext",
  backend = "ext",
  overwrite = FALSE,
  full_path = NULL
)

cache_options_get()
}
\arguments{
\item{path}{(character) End of directory path. Default: "fulltext".
See Details.}

\item{backend}{(character) Only "ext" supported for now.}

\item{overwrite}{(logical) overwrite cached file or not. Default: \code{FALSE}}

\item{full_path}{(character) Full cache path. If given \code{path} is ignored.}
}
\value{
both functions return the cache options in a named list
}
\description{
Set or get cache options
}
\section{Managing cached files}{

The dafault cache directory is \code{paste0(rappdirs::user_cache_dir(), "/R/fulltext")},
but you can set your own path using \code{cache_path_set()}

You can alternatively set the entire cache path with the \code{full_path}
parameter.

You can only pass \code{path} or \code{full_path} - but not both.

\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}

\examples{
\dontrun{
cache_options_get()
cache_options_set(path = "foobar")
cache_options_get()

# set full path
path <- tempdir()
cache_options_set(full_path = path)
cache_options_get()
}
}
\seealso{
\link{ftxt_cache}, \code{\link[=cache_file_info]{cache_file_info()}}

Other caching-functions: 
\code{\link{cache_file_info}()},
\code{\link{ftxt_cache}}
}
\concept{caching-functions}
