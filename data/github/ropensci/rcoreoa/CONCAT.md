rcoreoa
=======



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/rcoreoa/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rcoreoa/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/rcoreoa/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rcoreoa?branch=master)
[![cran checks](https://cranchecks.info/badges/worst/rcoreoa)](https://cranchecks.info/pkgs/rcoreoa)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcoreoa)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcoreoa)](https://cran.r-project.org/package=rcoreoa)

CORE API R client

CORE API docs: https://core.ac.uk/docs/

rcoreoa docs: https://docs.ropensci.org/rcoreoa/

Get an API key at https://core.ac.uk/api-keys/register. You'll need one, 
so do this now if you haven't yet. Once you have the key, you can pass it 
into the `key` parameter, or as a much better option store your key as an 
environment variable with the name `CORE_KEY` or an R option as `core_key`. 
See `?Startup` for how to work with env vars and R options

## About CORE

CORE's tagline is: "Aggregating the world's open access research papers"

CORE offers seamless access to millions of open access research papers, enrich
the collected data for text-mining and provide unique services to the research
community.

For more infos on CORE, see https://core.ac.uk/about

## Install


```r
install.packages("rcoreoa")
```

Development version


```r
remotes::install_github("ropensci/rcoreoa")
```


```r
library("rcoreoa")
```

## Get started

Get started with an introduction to the package: https://docs.ropensci.org/rcoreoa/articles/rcoreoa

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Aristotelis Charalampous](https://github.com/aristotelisxs)
* [Simon Goring](https://github.com/SimonGoring)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rcoreoa/issues).
* License: MIT
* Get citation information for `rcoreoa` in R doing `citation(package = 'rcoreoa')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rcoreoa 0.4.0
=============

### NEW FEATURES

* `core_articles_dedup()` function added for article deduplication (#23)
* `core_advanced_search()` overhauled because the CORE advanced search fields and syntax have completely changed. The interface changed from passing in a data.frame because we needed more flexibility to let users define how to combine multiple queries  (#20)

### MINOR IMPROVEMENTS

* using vcr for test suite http request caching (#21)
* add max value of `limit` (number of records per request) parameter to documentation (#18)

rcoreoa 0.3.0
=============

### NEW FEATURES

* gains new methods `core_articles_search()` and `core_articles_search_()` for searching for articles (#15)
* gains new manual file `?core_cache` for an overview of caching
* package `hoardr` now used for managing file caching (#13)
* function `core_advanced_search()` gains support for the `language` filter (#16) thanks @chreman

### MINOR IMPROVEMENTS

* support for passing in more than one journal/article/etc. identifier for `core_articles_pdf()`/`core_articles_history()` - already supported in other functions (#10)
* `overwrite` parameter was being ignored in `core_articles_pdf()`, now is passed on internally (#12)
* `parse` parameter dropped in `core_articles_pdf()`, only used internally (#11)
* Improved docs on how to get and use API keys (#14)
* Improve documentation about what a 404 error response means - that the thing reqeusted does not exist (#9)


rcoreoa 0.1.0
=============

### NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.2 Patched
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version adds a function, changes a function interface due to changes in the remove service, and improves documentation.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rcoreoa/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rcoreoa.git`
* Make sure to track progress upstream (i.e., on our version of `rcoreoa` at `ropensci/rcoreoa`) by doing `git remote add upstream https://github.com/ropensci/rcoreoa.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rcoreoa`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) 

Make sure not to share any secrets. -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
rcoreoa
=======

```{r echo=FALSE}
library("knitr")
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/rcoreoa/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rcoreoa/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/ropensci/rcoreoa/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rcoreoa?branch=master)
[![cran checks](https://cranchecks.info/badges/worst/rcoreoa)](https://cranchecks.info/pkgs/rcoreoa)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rcoreoa)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rcoreoa)](https://cran.r-project.org/package=rcoreoa)

CORE API R client

CORE API docs: https://core.ac.uk/docs/

rcoreoa docs: https://docs.ropensci.org/rcoreoa/

Get an API key at https://core.ac.uk/api-keys/register. You'll need one, 
so do this now if you haven't yet. Once you have the key, you can pass it 
into the `key` parameter, or as a much better option store your key as an 
environment variable with the name `CORE_KEY` or an R option as `core_key`. 
See `?Startup` for how to work with env vars and R options

## About CORE

CORE's tagline is: "Aggregating the world's open access research papers"

CORE offers seamless access to millions of open access research papers, enrich
the collected data for text-mining and provide unique services to the research
community.

For more infos on CORE, see https://core.ac.uk/about

## Install

```{r eval=FALSE}
install.packages("rcoreoa")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/rcoreoa")
```

```{r}
library("rcoreoa")
```

## Get started

Get started with an introduction to the package: https://docs.ropensci.org/rcoreoa/articles/rcoreoa

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Aristotelis Charalampous](https://github.com/aristotelisxs)
* [Simon Goring](https://github.com/SimonGoring)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rcoreoa/issues).
* License: MIT
* Get citation information for `rcoreoa` in R doing `citation(package = 'rcoreoa')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Introduction to the rcoreoa package"
author: "Scott Chamberlain, Simon Goring"
date: "2020-07-06"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to the rcoreoa package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`rcoreoa` - Client for the CORE API (https://core.ac.uk/docs/).
CORE (https://core.ac.uk) aggregates open access research
outputs from repositories and journals worldwide and make them
available to the public.

## Installation

The package can be installed directly from the CRAN repository:


```r
install.packages("rcoreoa")
```

Or you can install the development version:


```r
remotes::install_github("ropensci/rcoreoa")
```

Once the package is installed, simply run: 


```r
library("rcoreoa")
```

## Obtaining an API key

The Core API requires an API key, and, as such, requires you to register for the key on the Core Website (https://core.ac.uk/api-keys/register).  Once you register with your email address you will be sent an API key that looks a little like this:

`thisISyourAPIkeyITlooksLIKEaLONGstringOFletters`

## Using the API Key

Best practice is to set the API key as an environment variable for your system, and then call it in R using `Sys.getenv()`.  If you set the parameter in `.Renviron` it is permanently available to your R sessions.  See `?Startup`.  Be aware that if you are using version control you do not want to commit the `.Renviron` file in your local directory.  Either edit your global `.Renviron` file, or make sure that `.Renviron` is added to your `.gitignore` file.

Within the `.Renviron` file you will add:

```
CORE_KEY=thisISyourAPIkeyITlooksLIKEaLONGstringOFletters
```

The key may also be included in a file such as a `.bash_profile` file, or elsewhere.  Users may decide which works best for them.  Once you've added the API key, restart your R session and test to make sure the key has been added using the command:


```r
Sys.getenv("CORE_KEY")
```

If you get this to work, you're doing great and we can move on to the next section.  If this is still not working for you, check to make sure you have saved the `.Renviron` file, that it is in the same directory as your current project's working directory, and that the name you have given the variable in the `.Renviron` file is the same as the name you are calling in `Sys.getenv()`.

## An Introduction to the Functions

The `rcoreoa` package accesses CORE's API to facilitate open text searching.  The API allows an individual to search for articles based on text string searches using `core_search()`.  Given a set of article IDs from the `core_search()`, users can then find more bibliographic information on the article (`core_articles()`) and article publishing history (`core_articles_history()`), on the journals in which the article was published (`core_journals()`).

All of the functions return structured R objects, but can return JSON character strings by appending an unerscore (`_`) to the function name.  We will illustrate the difference:


```r
core_journals(id = '2167-8359')
#> $status
#> [1] "OK"
#> 
#> $data
#> $data$title
#> [1] "PeerJ"
#> 
#> $data$identifiers
#> [1] "oai:doaj.org/journal:576e4d34b8bf461bb586f1e90d80d7cc"
#> [2] "issn:2167-8359"                                       
#> [3] "url:https://doaj.org/toc/2167-8359"                   
#> 
#> $data$subjects
#> [1] "biomedical" "health"     "genetics"   "ecology"    "biology"   
#> [6] "Medicine"   "R"         
#> 
#> $data$language
#> [1] "EN"
#> 
#> $data$rights
#> [1] "CC BY"
#> 
#> $data$publisher
#> [1] "PeerJ Inc."
```

And with the underscore:


```r
core_journals_(id = '2167-8359')
#> [1] "{\"status\":\"OK\",\"data\":{\"title\":\"PeerJ\",\"identifiers\":[\"oai:doaj.org\\/journal:576e4d34b8bf461bb586f1e90d80d7cc\",\"issn:2167-8359\",\"url:https:\\/\\/doaj.org\\/toc\\/2167-8359\"],\"subjects\":[\"biomedical\",\"health\",\"genetics\",\"ecology\",\"biology\",\"Medicine\",\"R\"],\"language\":\"EN\",\"rights\":\"CC BY\",\"publisher\":\"PeerJ Inc.\"}}"
```

Through this Vignette we will illustrate some of the tools available as part of the package within a workflow that seeks to perform some basic bibliometric analysis.

## Pagination

Note that you are limited to a maximum of 100 results for the search functions;
use combination of `page` and `limit` parameters to paginate through results. 
For example:


```r
x1 <- core_search(query = 'ecology', limit = 100, page = 1)
x2 <- core_search(query = 'ecology', limit = 100, page = 2)
head(x1$data[,1:3])
#>                _index   _type       _id
#> 1 articles_2019_06_05 article 102353278
#> 2 articles_2019_06_05 article  20955435
#> 3 articles_2019_06_05 article 101526846
#> 4 articles_2019_06_05 article 103034034
#> 5 articles_2019_06_05 article 101912806
#> 6 articles_2019_06_05 article 103746797
head(x2$data[,1:3])
#>                _index   _type      _id
#> 1 articles_2019_06_05 article 79029462
#> 2 articles_2019_06_05 article 79046038
#> 3 articles_2019_06_05 article 79019557
#> 4 articles_2019_06_05 article 24465080
#> 5 articles_2019_06_05 article 79035703
#> 6 articles_2019_06_05 article 79024454
```

## high- vs. low-level interfaces

Each function has a higher level interface that does HTTP request for data and parses
the JSON using `jsonlite`. This is meant for those who want everything done for them,
but there's a time penalty for as the parsing adds extra time. If you just want raw JSON
unparsed text, you can use the low level interface.

The low level version of each function has `_` at the end (e.g., `core_search_`), while the
high level version doesn't have the `_` (e.g., `core_search`).

The high level version of each function uses the low level method, and the low level method
does all the logic and HTTP requesting, whereas the high level simply parses the output.

## Search


```r
res <- core_search(query = 'ecology', limit = 12)
tibble::as_tibble(res$data)
#> # A tibble: 12 x 5
#>    `_index` `_type` `_id` `_score` `_source`$id $authors $citations
#>    <chr>    <chr>   <chr>    <dbl> <chr>        <list>   <list>    
#>  1 article… article 1023…     18.5 102353278    <chr [2… <list [0]>
#>  2 article… article 2095…     18.5 20955435     <chr [1… <list [0]>
#>  3 article… article 1015…     17.5 101526846    <chr [1… <list [0]>
#>  4 article… article 1030…     17.5 103034034    <chr [1… <list [0]>
#>  5 article… article 1019…     17.5 101912806    <chr [1… <list [0]>
#>  6 article… article 1037…     17.5 103746797    <chr [1… <list [0]>
#>  7 article… article 1038…     17.5 103888859    <chr [1… <list [0]>
#>  8 article… article 1545…     17.5 154520564    <chr [1… <list [0]>
#>  9 article… article 1545…     17.5 154520563    <chr [1… <list [0]>
#> 10 article… article 1004…     17.5 100453958    <chr [1… <list [0]>
#> 11 article… article 1011…     17.5 101147175    <chr [1… <list [0]>
#> 12 article… article 1017…     17.5 101760738    <chr [1… <list [0]>
#> # … with 50 more variables: $contributors <list>, $datePublished <chr>,
#> #   $deleted <chr>, $description <chr>, $fullText <lgl>,
#> #   $fullTextIdentifier <chr>, $identifiers <list>, $journals <lgl>,
#> #   $language <lgl>, $duplicateId <lgl>, $publisher <chr>, $rawRecordXml <chr>,
#> #   $relations <list>, $repositories <list>,
#> #   $repositoryDocument$pdfStatus <int>, $$textStatus <int>,
#> #   $$metadataAdded <dbl>, $$metadataUpdated <dbl>, $$timestamp <dbl>,
#> #   $$depositedDate <dbl>, $$indexed <int>, $$deletedStatus <chr>,
#> #   $$pdfSize <int>, $$tdmOnly <lgl>, $$pdfOrigin <chr>, $similarities <lgl>,
#> #   $subjects <list>, $title <chr>, $topics <list>, $types <list>,
#> #   $urls <list>, $year <int>, $doi <lgl>, $oai <chr>, $downloadUrl <chr>,
#> #   $pdfHashValue <lgl>, $documentType <lgl>, $documentTypeConfidence <lgl>,
#> #   $citationCount <lgl>, $estimatedCitationCount <lgl>, $acceptedDate <lgl>,
#> #   $depositedDate <dbl>, $publishedDate <lgl>, $issn <lgl>,
#> #   $crossrefDocument <lgl>, $magDocument <lgl>, $attachmentCount <int>,
#> #   $repositoryPublicReleaseDate <lgl>, $extendedMetadataAttributes <lgl>,
#> #   $orcidAuthors <lgl>
```

## Advanced Search


```r
res <- core_advanced_search(core_query(identifiers='"oai:aura.abdn.ac.uk:2164/3837"',
  identifiers='"oai:aura.abdn.ac.uk:2164/3843"', op="OR"))
tibble::as_tibble(res$data[[1]])
#> # A tibble: 2 x 5
#>   `_index` `_type` `_id` `_score` `_source`$id $authors $citations $contributors
#>   <chr>    <chr>   <chr>    <dbl> <chr>        <list>   <list>     <list>       
#> 1 article… article 2549…     15.7 25497597     <chr [5… <list [0]> <chr [4]>    
#> 2 article… article 2549…     15.7 25497603     <chr [8… <list [0]> <chr [3]>    
#> # … with 50 more variables: $datePublished <chr>, $deleted <chr>,
#> #   $description <chr>, $fullText <chr>, $fullTextIdentifier <lgl>,
#> #   $identifiers <list>, $journals <list>, $language <lgl>, $duplicateId <lgl>,
#> #   $publisher <lgl>, $rawRecordXml <chr>, $relations <list>,
#> #   $repositories <list>, $repositoryDocument$pdfStatus <int>,
#> #   $$textStatus <int>, $$metadataAdded <dbl>, $$metadataUpdated <dbl>,
#> #   $$timestamp <dbl>, $$depositedDate <dbl>, $$indexed <int>,
#> #   $$deletedStatus <chr>, $$pdfSize <int>, $$tdmOnly <lgl>, $$pdfOrigin <lgl>,
#> #   $similarities <lgl>, $subjects <list>, $title <chr>, $topics <list>,
#> #   $types <list>, $urls <list>, $year <int>, $doi <chr>, $oai <chr>,
#> #   $downloadUrl <chr>, $pdfHashValue <chr>, $documentType <chr>,
#> #   $documentTypeConfidence <int>, $citationCount <lgl>,
#> #   $estimatedCitationCount <lgl>, $acceptedDate <lgl>, $depositedDate <dbl>,
#> #   $publishedDate <lgl>, $issn <lgl>, $attachmentCount <int>,
#> #   $repositoryPublicReleaseDate <lgl>,
#> #   $extendedMetadataAttributes$attachmentCount <int>,
#> #   $$publicReleaseDate <lgl>, $crossrefDocument <lgl>, $magDocument <lgl>,
#> #   $orcidAuthors <lgl>
```

# Articles


```r
core_articles(id = 21132995)
#> $status
#> [1] "OK"
#> 
#> $data
#> $data$id
#> [1] "21132995"
#> 
#> $data$authors
#> list()
#> 
#> $data$contributors
#> [1] "The Pennsylvania State University CiteSeerX Archives"
#> 
#> $data$datePublished
#> [1] "2010-02-17"
#> 
#> $data$description
#> [1] "Abstract. This paper discusses the potential contribution of an eco-theology to the management of marine resources. The claim of the Christian gospel is that God has a plan for everything in the universe and we are to live to bring it about. The Hebrew/Christian world view is significantly different from naturalism and humanism, the prevailing Greek world views in natural resource management. Three biblical paradigms are examined with insights into key elements in the management of fisheries: dominion; regulation and valuation and caring. In dominion we see the strength of mankind&apos;s rule over other species, including fish, misused. Fisheries management generally fails to reign in this driving force, rewarding greed while calling for restraint. Regulation and its impact on mindsets and behaviour, is a theme evident in the Old and New Testaments- entering the promised land, keeping the law and caring for others, both humans and fish. The biblical view also emphasises life, death and resurrection as the process seen in nature. New fishery management paradigms may only develop after old ways and attitudes have died. Significant attitudinal change is essential to improve fisheries management and to achieve new management arrangements. Improved fishery stewardship may require &quot;a new fisher&quot;, relationally mature and societally accountable, to achieve the goals of sustainable fishery management through a variety of policy paradigms. The Christian world view promotes such attitudinal change to improve stewardship through reconciling issues in our relationship with God, neighbours and nature. It is worthy of further investigation"
#> 
#> $data$identifiers
#> [1] "oai:CiteSeerX.psu:10.1.1.152.9831" NA                                 
#> 
#> $data$relations
#> list()
#> 
#> $data$repositories
#>    id openDoarId      name uri urlHomepage urlOaipmh uriJournals physicalName
#> 1 145          0 CiteSeerX  NA          NA        NA          NA       noname
#>   source software metadataFormat description journal roarId pdfStatus nrUpdates
#> 1     NA       NA             NA          NA      NA      0        NA         0
#>   disabled lastUpdateTime repositoryLocation
#> 1    FALSE             NA                 NA
#> 
#> $data$repositoryDocument
#> $data$repositoryDocument$pdfStatus
#> [1] 0
#> 
#> $data$repositoryDocument$metadataAdded
#> [1] 1.413993e+12
#> 
#> $data$repositoryDocument$metadataUpdated
#> [1] 1.529412e+12
#> 
#> $data$repositoryDocument$depositedDate
#> [1] 1.266365e+12
#> 
#> 
#> $data$subjects
#> [1] "text"
#> 
#> $data$topics
#> [1] "Ecology"              "Theology"             "Eco-theology"        
#> [4] "Fisheries Management"
#> 
#> $data$types
#> list()
#> 
#> $data$year
#> [1] 2010
#> 
#> $data$oai
#> [1] "oai:CiteSeerX.psu:10.1.1.152.9831"
#> 
#> $data$downloadUrl
#> [1] ""
```

# Article history


```r
core_articles_history(id = '21132995')
#> $status
#> [1] "OK"
#> 
#> $data
#>              datetime
#> 1 2016-08-03 00:13:41
#> 2 2014-10-22 16:42:14
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               metadata
#> 1                                                                                                               <record><header><identifier>\n    \n      \n        oai:CiteSeerX.psu:10.1.1.152.9831</identifier><datestamp>\n        2010-02-17</datestamp>\n      </header><metadata><oai_dc:dc xmlns:oai_dc="http://www.openarchives.org/OAI/2.0/oai_dc/" xmlns:dc="http://purl.org/dc/elements/1.1/" xsi:schemaLocation="http://www.openarchives.org/OAI/2.0/oai_dc/ http://www.openarchives.org/OAI/2.0/oai_dc.xsd" ><dc:title>\n      \n        \n          </dc:title><dc:subject>\n          Ecology</dc:subject><dc:subject>\n          Theology</dc:subject><dc:subject>\n          Eco-theology</dc:subject><dc:subject>\n          Fisheries Management</dc:subject><dc:description>\n          Abstract. This paper discusses the potential contribution of an eco-theology to the management of marine resources. The claim of the Christian gospel is that God has a plan for everything in the universe and we are to live to bring it about. The Hebrew/Christian world view is significantly different from naturalism and humanism, the prevailing Greek world views in natural resource management. Three biblical paradigms are examined with insights into key elements in the management of fisheries: dominion; regulation and valuation and caring. In dominion we see the strength of mankind&apos;s rule over other species, including fish, misused. Fisheries management generally fails to reign in this driving force, rewarding greed while calling for restraint. Regulation and its impact on mindsets and behaviour, is a theme evident in the Old and New Testaments- entering the promised land, keeping the law and caring for others, both humans and fish. The biblical view also emphasises life, death and resurrection as the process seen in nature. New fishery management paradigms may only develop after old ways and attitudes have died. Significant attitudinal change is essential to improve fisheries management and to achieve new management arrangements. Improved fishery stewardship may require &quot;a new fisher&quot;, relationally mature and societally accountable, to achieve the goals of sustainable fishery management through a variety of policy paradigms. The Christian world view promotes such attitudinal change to improve stewardship through reconciling issues in our relationship with God, neighbours and nature. It is worthy of further investigation.</dc:description><dc:contributor>\n          The Pennsylvania State University CiteSeerX Archives</dc:contributor><dc:publisher>\n          </dc:publisher><dc:date>\n          2010-02-17</dc:date><dc:date>\n          2010-02-17</dc:date><dc:format>\n          application/pdf</dc:format><dc:type>\n          text</dc:type><dc:identifier>\n          http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.152.9831</dc:identifier><dc:source>\n          http://www.orst.edu/Dept/IIFET/2000/papers/mcilgorm.pdf</dc:source><dc:language>\n          en</dc:language><dc:rights>\n          Metadata may be used without restrictions as long as the oai identifier remains attached to it.</dc:rights>\n        </oai_dc:dc>\n      </metadata>\n    </record>
#> 2 <record><header><identifier>\n        \n      \n    \n    \n      \n        oai:CiteSeerX.psu:10.1.1.152.9831</identifier><datestamp>\n        2010-02-17</datestamp>\n      </header><metadata><oai_dc:dc xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:oai_dc="http://www.openarchives.org/OAI/2.0/oai_dc/" xsi:schemaLocation="http://www.openarchives.org/OAI/2.0/oai_dc/ http://www.openarchives.org/OAI/2.0/oai_dc.xsd" ><dc:title>\n      \n      \n        \n          </dc:title><dc:subject>\n      \n      \n        \n          \n          Ecology</dc:subject><dc:subject>\n          Theology</dc:subject><dc:subject>\n          Eco-theology</dc:subject><dc:subject>\n          Fisheries Management</dc:subject><dc:description>\n          Abstract. This paper discusses the potential contribution of an eco-theology to the management of marine resources. The claim of the Christian gospel is that God has a plan for everything in the universe and we are to live to bring it about. The Hebrew/Christian world view is significantly different from naturalism and humanism, the prevailing Greek world views in natural resource management. Three biblical paradigms are examined with insights into key elements in the management of fisheries: dominion; regulation and valuation and caring. In dominion we see the strength of mankind&apos;s rule over other species, including fish, misused. Fisheries management generally fails to reign in this driving force, rewarding greed while calling for restraint. Regulation and its impact on mindsets and behaviour, is a theme evident in the Old and New Testaments- entering the promised land, keeping the law and caring for others, both humans and fish. The biblical view also emphasises life, death and resurrection as the process seen in nature. New fishery management paradigms may only develop after old ways and attitudes have died. Significant attitudinal change is essential to improve fisheries management and to achieve new management arrangements. Improved fishery stewardship may require &quot;a new fisher&quot;, relationally mature and societally accountable, to achieve the goals of sustainable fishery management through a variety of policy paradigms. The Christian world view promotes such attitudinal change to improve stewardship through reconciling issues in our relationship with God, neighbours and nature. It is worthy of further investigation.</dc:description><dc:contributor>\n          The Pennsylvania State University CiteSeerX Archives</dc:contributor><dc:publisher>\n          </dc:publisher><dc:date>\n          \n          2010-02-17</dc:date><dc:date>\n          2010-02-17</dc:date><dc:format>\n          application/pdf</dc:format><dc:type>\n          text</dc:type><dc:identifier>\n          http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.152.9831</dc:identifier><dc:source>\n          http://www.orst.edu/Dept/IIFET/2000/papers/mcilgorm.pdf</dc:source><dc:language>\n          en</dc:language><dc:rights>\n          Metadata may be used without restrictions as long as the oai identifier remains attached to it.</dc:rights>\n        </oai_dc:dc>\n        \n      </metadata>\n        \n      \n    </record>
```

# Journals


```r
core_journals(id = '2220-721X')
```

# Get PDFs

The `_` for these methods means that you get a file path back to the PDF, while the
high level version without the `_` parses the pdf to text for you.


```r
core_articles_pdf_(11549557)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcoreoa-package.R
\docType{package}
\name{rcoreoa-package}
\alias{rcoreoa-package}
\alias{rcoreoa}
\title{rcoreoa - CORE R client}
\description{
CORE is a web service for metadata on scholarly journal articles. Find
CORE at \url{https://core.ac.uk} and their API docs at
\url{https://core.ac.uk/docs/}.
}
\section{Package API}{

Each API endpoint has two functions that interface with it - a higher
level interface and a lower level interface. The lower level functions
have an underscore (\verb{_}) at the end of the function name, while their
corresponding higher level companions do not. The higher level functions
parse to list/data.frame's (as tidy as possible). Lower level
functions give back JSON (character class) thus are slightly faster not
spending time on parsing to R structures.
\itemize{
\item \code{\link[=core_articles]{core_articles()}} / \code{\link[=core_articles_]{core_articles_()}} - get article metadata
\item \code{\link[=core_articles_history]{core_articles_history()}} / \code{\link[=core_articles_history_]{core_articles_history_()}} - get
article history metadata
\item \code{\link[=core_articles_pdf]{core_articles_pdf()}} / \code{\link[=core_articles_pdf_]{core_articles_pdf_()}} - download
article PDF, and optionally extract text
\item \code{\link[=core_journals]{core_journals()}} / \code{\link[=core_journals_]{core_journals_()}} - get journal metadata
\item \code{\link[=core_repos]{core_repos()}} / \code{\link[=core_repos_]{core_repos_()}} - get repository metadata
\item \code{\link[=core_repos_search]{core_repos_search()}} / \code{\link[=core_repos_search_]{core_repos_search_()}} - search for
repositories
\item \code{\link[=core_search]{core_search()}} / \code{\link[=core_search_]{core_search_()}} - search articles
\item \code{\link[=core_advanced_search]{core_advanced_search()}} -  advanced search of articles
}
}

\section{Authentication}{

You'll need a CORE API token/key to use this package. Get one at
\url{https://core.ac.uk/api-keys/register}
}

\section{Pagination}{

Note that you are limited to a maximum of 100 results for the search
functions; use combination of \code{page} and \code{limit} parameters to
paginate through results. For example:\preformatted{x1 <- core_search(query = 'ecology', limit = 100, page = 1)
x2 <- core_search(query = 'ecology', limit = 100, page = 2)
}
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}

Aristotelis Charalampous
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_journals.R, R/core_repos.R
\name{core_journals}
\alias{core_journals}
\alias{core_journals_}
\alias{core_repos_}
\title{Get journal via its ISSN}
\usage{
core_journals(id, key = NULL, method = "GET", parse = TRUE, ...)

core_journals_(id, key = NULL, method = "GET", ...)

core_repos_(id, key = NULL, method = "GET", ...)
}
\arguments{
\item{id}{(integer) One or more journal ISSNs. Required}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{method}{(character) one of 'GET' (default) or 'POST'}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Get journal via its ISSN
}
\details{
\code{core_journals} does the HTTP request and parses, while
\code{core_journals_} just does the HTTP request, gives back JSON as a
character string

These functions take one article ID at a time. Use lapply/loops/etc for
many ids
}
\examples{
\dontrun{
core_journals(id = '2167-8359')

ids <- c("2167-8359", "2050-084X")
res <- lapply(ids, core_journals)
vapply(res, "[[", "", c("data", "title"))

# just http request, get text back
core_journals_('2167-8359')

# post request, ideal for lots of ISSNs
if (requireNamespace("rcrossref", quietly = TRUE)) {
 library(rcrossref)
 res <- lapply(c("bmc", "peerj", "elife", "plos", "frontiers"), function(z)
    cr_journals(query = z))
 ids <- na.omit(unlist(lapply(res, function(b) b$data$issn)))
 out <- core_journals(ids, method = "POST")
 head(out)
}

}
}
\references{
\url{https://core.ac.uk/docs/#!/journals/getJournalByIssnBatch}
\url{https://core.ac.uk/docs/#!/journals/getJournalByIssn}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_articles.R
\name{core_articles}
\alias{core_articles}
\alias{core_articles_}
\title{Get articles}
\usage{
core_articles(
  id,
  metadata = TRUE,
  fulltext = FALSE,
  citations = FALSE,
  similar = FALSE,
  duplicate = FALSE,
  urls = FALSE,
  extractedUrls = FALSE,
  faithfulMetadata = FALSE,
  key = NULL,
  method = "GET",
  parse = TRUE,
  ...
)

core_articles_(
  id,
  metadata = TRUE,
  fulltext = FALSE,
  citations = FALSE,
  similar = FALSE,
  duplicate = FALSE,
  urls = FALSE,
  extractedUrls = FALSE,
  faithfulMetadata = FALSE,
  key = NULL,
  method = "GET",
  ...
)
}
\arguments{
\item{id}{(integer) CORE ID of the article that needs to be fetched.
Required}

\item{metadata}{Whether to retrieve the full article metadata or only the
ID. Default: \code{TRUE}}

\item{fulltext}{Whether to retrieve full text of the article. Default:
\code{FALSE}}

\item{citations}{Whether to retrieve citations found in the article.
Default: \code{FALSE}}

\item{similar}{Whether to retrieve a list of similar articles.
Default: \code{FALSE}
Because the similar articles are calculated on demand, setting this
parameter to true might slightly slow down the response time    query}

\item{duplicate}{Whether to retrieve a list of CORE IDs of different
versions of the article. Default: \code{FALSE}}

\item{urls}{Whether to retrieve a list of URLs from which the article can
be downloaded. This can include links to PDFs as well as HTML pages.
Default: \code{FALSE}}

\item{extractedUrls}{Whether to retrieve a list of URLs which were extracted
from the article fulltext. Default: \code{FALSE}. This parameter is not
available in CORE API v2.0 beta}

\item{faithfulMetadata}{Whether to retrieve the raw XML metadata of the
article. Default: \code{FALSE}}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{method}{(character) one of 'GET' (default) or 'POST'}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Get articles
}
\details{
\code{core_articles} does the HTTP request and parses, while
\code{core_articles_} just does the HTTP request, gives back JSON as a
character string

These functions take one article ID at a time. Use lapply/loops/etc for
many ids
}
\examples{
\dontrun{
core_articles(id = 21132995)
core_articles(id = 21132995, fulltext = TRUE)
core_articles(id = 21132995, citations = TRUE)

# when passing >1 article ID
ids <- c(20955435, 21132995, 21813171, 22815670, 23828884,
   23465055, 23831838, 23923390, 22559733)
## you can use method="GET" with lapply or similar
res <- lapply(ids, core_articles)
vapply(res, "[[", "", c("data", "datePublished"))

## or use method="POST" passing all at once
res <- core_articles(ids, method = "POST")
head(res$data)
res$data$authors

# just http request, get text back
core_articles_(id = '21132995')
## POST, can pass many at once
core_articles_(id = ids, method = "POST")
}
}
\references{
\url{https://core.ac.uk/docs/#!/articles/getArticleByCoreIdBatch}
\url{https://core.ac.uk/docs/#!/articles/getArticleByCoreId}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_articles_search.R
\name{core_articles_search}
\alias{core_articles_search}
\alias{core_articles_search_}
\title{Search CORE articles}
\usage{
core_articles_search(
  query,
  metadata = TRUE,
  fulltext = FALSE,
  citations = FALSE,
  similar = FALSE,
  duplicate = FALSE,
  urls = FALSE,
  faithfulMetadata = FALSE,
  page = 1,
  limit = 10,
  key = NULL,
  parse = TRUE,
  ...
)

core_articles_search_(
  query,
  metadata = TRUE,
  fulltext = FALSE,
  citations = FALSE,
  similar = FALSE,
  duplicate = FALSE,
  urls = FALSE,
  faithfulMetadata = FALSE,
  page = 1,
  limit = 10,
  key = NULL,
  ...
)
}
\arguments{
\item{query}{(character) query string, required}

\item{metadata}{(logical) Whether to retrieve the full article metadata or
only the ID. Default: \code{TRUE}}

\item{fulltext}{(logical) Whether to retrieve full text of the article.
Default: \code{FALSE}}

\item{citations}{(logical) Whether to retrieve citations found in the
article. Default: \code{FALSE}}

\item{similar}{(logical) Whether to retrieve a list of similar articles.
Default: \code{FALSE}. Because the similar articles are calculated on demand,
setting this parameter to true might slightly slow down the response time}

\item{duplicate}{(logical) Whether to retrieve a list of CORE IDs of
different versions of the article. Default: \code{FALSE}}

\item{urls}{(logical) Whether to retrieve a list of URLs from which the
article can be downloaded. This can include links to PDFs as well as
HTML pages. Default: \code{FALSE}}

\item{faithfulMetadata}{(logical) Returns the records raw XML metadata
from the original repository. Default: \code{FALSE}}

\item{page}{(character) page number (default: 1), optional}

\item{limit}{(character) records to return (default: 10, minimum: 10,
maximum: 100), optional}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Search CORE articles
}
\details{
\code{core_articles_search} does the HTTP request and parses, while
\code{core_articles_search_} just does the HTTP request, gives back JSON as a character
string
}
\examples{
\dontrun{
core_articles_search(query = 'ecology')
core_articles_search(query = 'ecology', parse = FALSE)
core_articles_search(query = 'ecology', limit = 12)
out = core_articles_search(query = 'ecology', fulltext = TRUE)

core_articles_search_(query = 'ecology')
jsonlite::fromJSON(core_articles_search_(query = 'ecology'))

# post request
query <- c('data mining', 'semantic web')
res <- core_articles_search(query)
head(res$data)
res$data[[2]]$doi
}
}
\references{
\url{https://core.ac.uk/docs/#!/all/search}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_repos_search.R, R/core_search.R
\name{core_repos_search_}
\alias{core_repos_search_}
\alias{core_search}
\alias{core_search_}
\title{Search CORE}
\usage{
core_repos_search_(query, page = 1, limit = 10, key = NULL, ...)

core_search(query, page = 1, limit = 10, key = NULL, parse = TRUE, ...)

core_search_(query, page = 1, limit = 10, key = NULL, ...)
}
\arguments{
\item{query}{(character) query string, required}

\item{page}{(character) page number (default: 1), optional}

\item{limit}{(character) records to return (default: 10, minimum: 10,
maximum: 100), optional}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}
}
\description{
Search CORE
}
\details{
\code{core_search} does the HTTP request and parses, while
\code{core_search_} just does the HTTP request, gives back JSON as a character
string
}
\examples{
\dontrun{
core_search(query = 'ecology')
core_search(query = 'ecology', parse = FALSE)
core_search(query = 'ecology', limit = 12)

core_search_(query = 'ecology')
library("jsonlite")
jsonlite::fromJSON(core_search_(query = 'ecology'))

# post request
query <- c('data mining', 'machine learning', 'semantic web')
res <- core_search(query)
res
res$status
res$totalHits
res$data
head(res$data[[1]])
}
}
\references{
https://core.ac.uk/docs/#!/all/search
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_advanced_search.R
\name{core_advanced_search}
\alias{core_advanced_search}
\alias{core_query}
\title{Advanced Search CORE}
\usage{
core_advanced_search(
  ...,
  page = 1,
  limit = 10,
  key = NULL,
  parse = TRUE,
  .list = list()
)

core_query(..., op = "AND")
}
\arguments{
\item{...}{for \code{core_query()}, query fields, see Details. for
\code{core_advanced_search()} any number of queries, results of calling
\code{core_query()}. Required. See Details.}

\item{page}{(character) page number (default: 1), optional}

\item{limit}{(character) records to return (default: 10, minimum: 10,
maximum: 100), optional}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{.list}{alternative to passing \code{core_query()} calls to \code{...};
just create a list of them and pass to this parameter; easier for
programming with}

\item{op}{(character) operator to combine multiple fields. options:
\code{AND}, \code{OR}}
}
\value{
data.frame with the following columns:
\itemize{
\item \code{status}: string, which will be 'OK' or 'Not found' or
'Too many queries' or 'Missing parameter' or 'Invalid parameter' or
'Parameter out of bounds'
\item \code{totalHits}: integer, Total number of items matching the search criteria
\item \code{data}: list, a list of relevant resource
}
}
\description{
Advanced Search CORE
}
\details{
\code{query} should be one or more calls to \code{core_query()},
(at least one is required):
\itemize{
\item \code{title}
\item \code{description}
\item \code{fullText}
\item \code{authors}
\item \code{publisher}: string, to be used as an absolute match against the
publisher name metadata field
\item \code{repositories.id}: repository id
\item \code{repositories.name}: repository name
\item \code{doi}: string, to be used as an absolute match against the repository
name metadata field (all other fields will be ignored if included)
\item \code{oai}
\item \code{identifiers}
\item \code{language.name}: string, to filter against languages specified in
https://en.wikipedia.org/wiki/ISO_639-1
\item \code{year}: year to filter to
}

\code{core_advanced_search} does the HTTP request and parses, while
\code{core_advanced_search_} just does the HTTP request, gives back JSON as a
character string
}
\examples{
\dontrun{
## compose queries
core_query(title="psychology", year=2014)
core_query(title="psychology", year=2014, op="OR")
core_query(identifiers='"oai:aura.abdn.ac.uk:2164/3837"',
  identifiers='"oai:aura.abdn.ac.uk:2164/3843"', op="OR")

## do actual searches
core_advanced_search(
  core_query(identifiers='"oai:aura.abdn.ac.uk:2164/3837"',
    identifiers='"oai:aura.abdn.ac.uk:2164/3843"', op="OR")
)

res <- core_advanced_search(
  core_query(title="psychology"),
  core_query(doi='"10.1186/1471-2458-6-309"'),
  core_query(language.name="german")
)
res
}
}
\references{
https://core.ac.uk/docs/#!/all/searchBatch
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caching.R
\name{core_cache}
\alias{core_cache}
\title{Caching}
\description{
Manage cached \code{rcoreoa} files with \pkg{hoardr}
}
\details{
The dafault cache directory is
\code{paste0(rappdirs::user_cache_dir(), "/R/rcoreoa")}, but you can set
your own path using \code{cache_path_set()}

\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}
\section{Useful user functions}{

\itemize{
\item \code{core_cache$cache_path_get()} get cache path
\item \code{core_cache$cache_path_set()} set cache path
\item \code{core_cache$list()} returns a character vector of full
path file names
\item \code{core_cache$files()} returns file objects with metadata
\item \code{core_cache$details()} returns files with details
\item \code{core_cache$delete()} delete specific files
\item \code{core_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
core_cache

# list files in cache
core_cache$list()

# delete certain database files
# core_cache$delete("file path")
# core_cache$list()

# delete all files in cache
# core_cache$delete_all()
# core_cache$list()

# set a different cache path from the default
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_articles_pdf.R
\name{core_articles_pdf}
\alias{core_articles_pdf}
\alias{core_articles_pdf_}
\title{Download article pdf}
\usage{
core_articles_pdf(id, key = NULL, overwrite = FALSE, ...)

core_articles_pdf_(id, key = NULL, overwrite = FALSE, ...)
}
\arguments{
\item{id}{(integer) CORE ID of the article that needs to be fetched.
One or more. Required}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{overwrite}{(logical) overwrite file or not if already
on disk. Default: \code{FALSE}}

\item{...}{Curl options passed to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\value{
\code{core_articles_pdf_} returns a file path on success.
When many IDs passed to \code{core_articles_pdf} it returns a list (equal to
length of IDs) where each element is a character vector of length equal
to number of pages in the PDF; but on failure throws warning and returns
NULL. When single ID apssed to \code{core_articles_pdf} it returns a character
vector of length equal to number of pages in the PDF, but on failure
stops with message
}
\description{
Download article pdf
}
\details{
\code{core_articles_pdf} does the HTTP request and parses
PDF to text, while \code{core_articles_pdf_} just does the HTTP request
and gives back the path to the file

If you get a message like \verb{Error: Not Found (HTTP 404)}, that means
a PDF was not found. That is, it does not exist. That is, there is no
PDF associated with the article ID you searched for. This is the
correct behavior, and nothing is wrong with this function or this
package. We could do another web request to check if the id you
pass in has a PDF or not first, but that's another request, slowing
this function down.
}
\examples{
\dontrun{
# just http request, get file path back
core_articles_pdf_(11549557)

# get paper and parse to text
core_articles_pdf(11549557)

ids <- c(11549557, 385071)
res <- core_articles_pdf(ids)
cat(res[[1]][1])
cat(res[[2]][1])
}
}
\references{
\url{https://core.ac.uk/docs/#!/articles/getArticlePdfByCoreId}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_articles_history.R
\name{core_articles_history}
\alias{core_articles_history}
\alias{core_articles_history_}
\title{Get article history}
\usage{
core_articles_history(id, page = 1, limit = 10, key = NULL, parse = TRUE, ...)

core_articles_history_(id, page = 1, limit = 10, key = NULL, ...)
}
\arguments{
\item{id}{(integer) CORE ID of the article that needs to be fetched.
One or more. Required}

\item{page}{(character) page number (default: 1), optional}

\item{limit}{(character) records to return (default: 10, minimum: 10,
maximum: 100), optional}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
\code{core_articles_history_} returns a JSON string on success.
\code{core_articles_history} returns a list (equal to \code{id} length) where each
element is a list of length two with elements for data and status of
the request; on failure the data slot is NULL.
}
\description{
Get article history
}
\details{
\code{core_articles_history} does the HTTP request and parses,
while \code{core_articles_history_} just does the HTTP request, gives back JSON
as a character string
}
\examples{
\dontrun{
core_articles_history(id = 21132995)

ids <- c(20955435, 21132995, 21813171, 22815670, 14045109, 23828884,
   23465055, 23831838, 23923390, 22559733)
res <- core_articles_history(ids)
vapply(res, function(z) z$data$datetime[1], "")

# just http request, get text back
core_articles_history_(21132995)
}
}
\references{
\url{https://core.ac.uk/docs/#!/articles/getArticleHistoryByCoreId}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_repos_search.R
\name{core_repos_search}
\alias{core_repos_search}
\title{Search CORE repositories}
\usage{
core_repos_search(query, page = 1, limit = 10, key = NULL, parse = TRUE, ...)
}
\arguments{
\item{query}{(character) query string, required}

\item{page}{(character) page number (default: 1), optional}

\item{limit}{(character) records to return (default: 10, minimum: 10,
maximum: 100), optional}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Search CORE repositories
}
\details{
\code{core_repos_search} does the HTTP request and parses, while
\code{core_repos_search_} just does the HTTP request, gives back JSON as
a character string

A POST method is allowed on this route, but it's not supported yet.
}
\examples{
\dontrun{
core_repos_search(query = 'mathematics')
core_repos_search(query = 'physics', parse = FALSE)
core_repos_search(query = 'pubmed')

core_repos_search_(query = 'pubmed')
library("jsonlite")
jsonlite::fromJSON(core_repos_search_(query = 'pubmed'))
}
}
\references{
\url{https://core.ac.uk/docs/#!/repositories/search}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_articles_dedup.R
\name{core_articles_dedup}
\alias{core_articles_dedup}
\alias{core_articles_dedup_}
\title{Article deduplication}
\usage{
core_articles_dedup(
  doi = NULL,
  title = NULL,
  year = NULL,
  description = NULL,
  fulltext = NULL,
  identifier = NULL,
  repositoryId = NULL,
  key = NULL,
  parse = TRUE,
  ...
)

core_articles_dedup_(
  doi = NULL,
  title = NULL,
  year = NULL,
  description = NULL,
  fulltext = NULL,
  identifier = NULL,
  repositoryId = NULL,
  key = NULL,
  ...
)
}
\arguments{
\item{doi}{(character) the DOI for which the duplicates will be identified.
optional}

\item{title}{(character) title to match when looking for duplicate articles.
Either \code{year} or \code{description} should also be supplied if this parameter is
supplied. optional}

\item{year}{(character) year the article was published. Only used in
combination with the value for \code{title} parameter. optional}

\item{description}{(character) abstract for an article based on which its
duplicates will be found. This should be more than 500 characters. Value
for the \code{title} parameter should also be supplied if this is supplied.
optional}

\item{fulltext}{(character) Full text for an article based on which its
duplicates will be found. optional}

\item{identifier}{(character) CORE ID of the article for which the
duplicates will be identified. optional}

\item{repositoryId}{(character) Limit the duplicates search to
particular repository id. optional}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Article deduplication
}
\examples{
\dontrun{
core_articles_dedup(title = "Managing exploratory innovation", year = 2010)
core_articles_dedup_(title = "Managing exploratory innovation", year = 2010)

ab = 'Neonicotinoid seed dressings have caused concern world-wide. We use
large field experiments to assess the effects of neonicotinoid-treated crops
on three bee species across three countries (Hungary, Germany, and the
United Kingdom). Winter-sown oilseed rape was grown commercially with
either seed coatings containing neonicotinoids (clothianidin or
thiamethoxam) or no seed treatment (control). For honey bees, we found
both negative (Hungary and United Kingdom) and positive (Germany)
effects during crop flowering. In Hungary, negative effects on honey
bees (associated with clothianidin) persisted over winter and resulted
in smaller colonies in the following spring (24\% declines). In wild
bees (Bombus terrestris and Osmia bicornis), reproduction was
negatively correlated with neonicotinoid residues. These findings
point to neonicotinoids causing a reduced capacity of bee species
to establish new populations in the year following exposure.'
core_articles_dedup(
  title = "Country-specific effects of neonicotinoid pesticides on honey bees and wild bees",
  description = ab, verbose = TRUE)
}
}
\references{
\url{https://core.ac.uk/docs/#!/articles/nearDuplicateArticles}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_repos.R
\name{core_repos}
\alias{core_repos}
\title{Get repositories via their repository IDs}
\usage{
core_repos(id, key = NULL, method = "GET", parse = TRUE, ...)
}
\arguments{
\item{id}{(integer) One or more repository IDs. Required}

\item{key}{A CORE API key. Get one at
\url{https://core.ac.uk/api-keys/register}. Once you have the key,
you can pass it into this parameter, or as a much better option,
store your key as an environment variable with the name
\code{CORE_KEY} or an R option as \code{core_key}. See \code{?Startup}
for how to work with env vars and R options}

\item{method}{(character) one of 'GET' (default) or 'POST'}

\item{parse}{(logical) Whether to parse to list (\code{FALSE}) or
data.frame (\code{TRUE}). Default: \code{TRUE}}

\item{...}{Curl options passed to \code{\link[crul]{HttpClient}}}
}
\description{
Get repositories via their repository IDs
}
\details{
\code{core_repos} does the HTTP request and parses, while
\code{core_repos_} just does the HTTP request, gives back JSON as a
character string

These functions take one article ID at a time. Use lapply/loops/etc for
many ids
}
\examples{
\dontrun{
core_repos(id = 507)
core_repos(id = 444)

ids <- c(507, 444, 70)
res <- lapply(ids, core_repos)
vapply(res, "[[", "", c("data", "name"))

# just http request, get json as character vector back
core_repos_(507)
}
}
\references{
\url{https://core.ac.uk/docs/#!/repositories/getRepositoryById}
\url{https://core.ac.uk/docs/#!/repositories/getRepositoryByIdBatch}
}
