rplos
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rplos)](https://cranchecks.info/pkgs/rplos)
[![R-check](https://github.com/ropensci/rplos/workflows/R-check/badge.svg)](https://github.com/ropensci/rplos/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/rplos/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rplos?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rplos)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rplos)](https://cran.r-project.org/package=rplos)

## Install

You can get this package at CRAN [here](https://cran.r-project.org/package=rplos), or install it within R by doing


```r
install.packages("rplos")
```

Or install the development version from GitHub


```r
remotes::install_github("ropensci/rplos")
```


```r
library("rplos")
```

## What is this?

`rplos` is a package for accessing full text articles from the Public Library of Science journals using their API.

## Information

You used to need a key to use `rplos` - you no longer do as of 2015-01-13 (or `v0.4.5.999`).

rplos vignetttes: <https://docs.ropensci.org/rplos/>

PLOS API documentation: <http://api.plos.org/>

PLOS Solr schema is at <https://gist.github.com/openAccess/9e76aa7fa6135be419968b1372c86957> but is 1.5 years old so may not be up to date.

Crossref API documentation can be found at <https://github.com/CrossRef/rest-api-doc>. See also [rcrossref](https://github.com/ropensci/rcrossref) ([on CRAN](https://cran.r-project.org/package=rcrossref)) with a much fuller implementation of R functions for all Crossref endpoints.

## Throttling

Beware, PLOS recently has started throttling requests. That is,
they will give error messages like "(503) Service Unavailable -
The server cannot process the request due to a high load", which
means you've done too many requests in a certain time period. Here's
[what they say](http://api.plos.org/solr/faq/#solr_api_recommended_usage) on the matter:

> Please limit your API requests to 7200 requests a day, 300 per hour, 10 per minute and allow 5 seconds for your search to return results. If you exceed this threshold, we will lock out your IP address. If you're a high-volume user of the PLOS Search API and need more API requests a day, please contact us at api@plos.org to discuss your options. We currently limit API users to no more than five concurrent connections from a single IP address.

## Quick start

### Search

Search for the term ecology, and return id (DOI) and publication date, limiting to 5 items


```r
searchplos('ecology', 'id,publication_date', limit = 5)
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1    55873     0
#> 
#> $data
#> # A tibble: 5 x 2
#>   id                           publication_date    
#>   <chr>                        <chr>               
#> 1 10.1371/journal.pone.0001248 2007-11-28T00:00:00Z
#> 2 10.1371/journal.pone.0059813 2013-04-24T00:00:00Z
#> 3 10.1371/journal.pone.0080763 2013-12-10T00:00:00Z
#> 4 10.1371/journal.pone.0246749 2021-02-08T00:00:00Z
#> 5 10.1371/journal.pone.0220747 2019-08-01T00:00:00Z
```

Get DOIs for full article in PLoS One


```r
searchplos(q="*:*", fl='id', fq=list('journal_key:PLoSONE',
   'doc_type:full'), limit=5)
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1   246927     0
#> 
#> $data
#> # A tibble: 5 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pone.0002399
#> 2 10.1371/journal.pone.0002401
#> 3 10.1371/journal.pone.0002403
#> 4 10.1371/journal.pone.0002405
#> 5 10.1371/journal.pone.0002407
```

Query to get some PLOS article-level metrics, notice difference between two outputs


```r
out <- searchplos(q="*:*", fl=c('id','counter_total_all','alm_twitterCount'), fq='doc_type:full')
out_sorted <- searchplos(q="*:*", fl=c('id','counter_total_all','alm_twitterCount'),
   fq='doc_type:full', sort='counter_total_all desc')
head(out$data)
#> # A tibble: 6 x 3
#>   id                           alm_twitterCount counter_total_all
#>   <chr>                                   <int>             <int>
#> 1 10.1371/journal.pcbi.0020071                0             17862
#> 2 10.1371/journal.pbio.1000152                0              7180
#> 3 10.1371/journal.pbio.1000153                5             27702
#> 4 10.1371/journal.pbio.1000159                0             14256
#> 5 10.1371/journal.pbio.1000165                8             30629
#> 6 10.1371/journal.pbio.1000166                0              6884
head(out_sorted$data)
#> # A tibble: 6 x 3
#>   id                           alm_twitterCount counter_total_all
#>   <chr>                                   <int>             <int>
#> 1 10.1371/journal.pmed.0020124             3890           3226572
#> 2 10.1371/journal.pone.0133079              312           2481307
#> 3 10.1371/journal.pcbi.1003149              216           1766230
#> 4 10.1371/journal.pmed.1000376                9           1213573
#> 5 10.1371/journal.pone.0141854             3804            996693
#> 6 10.1371/journal.pcbi.0030102               65            986338
```

A list of articles about social networks that are popular on a social network


```r
searchplos(q="*:*",fl=c('id','alm_twitterCount'),
   fq=list('doc_type:full','subject:"Social networks"','alm_twitterCount:[100 TO 10000]'),
   sort='counter_total_month desc')
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1       64     0
#> 
#> $data
#> # A tibble: 10 x 2
#>    id                           alm_twitterCount
#>    <chr>                                   <int>
#>  1 10.1371/journal.pmed.1000316             1199
#>  2 10.1371/journal.pone.0069841              896
#>  3 10.1371/journal.pone.0073791             1884
#>  4 10.1371/journal.pcbi.1005399              604
#>  5 10.1371/journal.pcbi.1007513              122
#>  6 10.1371/journal.pone.0236517              123
#>  7 10.1371/journal.pbio.1002195              778
#>  8 10.1371/journal.pone.0118093              133
#>  9 10.1371/journal.pone.0155923              200
#> 10 10.1371/journal.pone.0162678              172
```

Show all articles that have these two words less then about 15 words apart


```r
searchplos(q='everything:"sports alcohol"~15', fl='title', fq='doc_type:full', limit=3)
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1      163     0
#> 
#> $data
#> # A tibble: 3 x 1
#>   title                                                                         
#>   <chr>                                                                         
#> 1 Alcohol Advertising in Sport and Non-Sport TV in Australia, during Children’s…
#> 2 Alcohol intoxication at Swedish football matches: A study using biological sa…
#> 3 Symptoms of Insomnia and Sleep Duration and Their Association with Incident S…
```

Narrow results to 7 words apart, changing the ~15 to ~7


```r
searchplos(q='everything:"sports alcohol"~7', fl='title', fq='doc_type:full', limit=3)
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1       91     0
#> 
#> $data
#> # A tibble: 3 x 1
#>   title                                                                         
#>   <chr>                                                                         
#> 1 Alcohol Advertising in Sport and Non-Sport TV in Australia, during Children’s…
#> 2 Alcohol intoxication at Swedish football matches: A study using biological sa…
#> 3 Symptoms of Insomnia and Sleep Duration and Their Association with Incident S…
```

Remove DOIs for annotations (i.e., corrections) and Viewpoints articles


```r
searchplos(q='*:*', fl=c('id','article_type'),
   fq=list('-article_type:correction','-article_type:viewpoints'), limit=5)
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1  2434439     0
#> 
#> $data
#> # A tibble: 5 x 2
#>   id                                                  article_type    
#>   <chr>                                               <chr>           
#> 1 10.1371/journal.pbio.1000146/title                  Unsolved Mystery
#> 2 10.1371/journal.pbio.1000146/abstract               Unsolved Mystery
#> 3 10.1371/journal.pbio.1000146/references             Unsolved Mystery
#> 4 10.1371/journal.pbio.1000146/body                   Unsolved Mystery
#> 5 10.1371/journal.pbio.1000146/supporting_information Unsolved Mystery
```

### Faceted search

Facet on multiple fields


```r
facetplos(q='alcohol', facet.field=c('journal','subject'), facet.limit=5)
#> $facet_queries
#> NULL
#> 
#> $facet_fields
#> $facet_fields$journal
#> # A tibble: 5 x 2
#>   term                             value
#>   <chr>                            <chr>
#> 1 plos one                         31656
#> 2 plos genetics                    702  
#> 3 plos medicine                    699  
#> 4 plos neglected tropical diseases 637  
#> 5 plos pathogens                   459  
#> 
#> $facet_fields$subject
#> # A tibble: 5 x 2
#>   term                          value
#>   <chr>                         <chr>
#> 1 biology and life sciences     33275
#> 2 medicine and health sciences  30309
#> 3 research and analysis methods 17340
#> 4 biochemistry                  14587
#> 5 physical sciences             12502
#> 
#> 
#> $facet_pivot
#> NULL
#> 
#> $facet_dates
#> NULL
#> 
#> $facet_ranges
#> NULL
```

Range faceting


```r
facetplos(q='*:*', url=url, facet.range='counter_total_all',
 facet.range.start=5, facet.range.end=100, facet.range.gap=10)
#> $facet_queries
#> NULL
#> 
#> $facet_fields
#> NULL
#> 
#> $facet_pivot
#> NULL
#> 
#> $facet_dates
#> NULL
#> 
#> $facet_ranges
#> $facet_ranges$counter_total_all
#> # A tibble: 10 x 2
#>    term  value
#>    <chr> <chr>
#>  1 5     39   
#>  2 15    12   
#>  3 25    33   
#>  4 35    47   
#>  5 45    47   
#>  6 55    44   
#>  7 65    59   
#>  8 75    80   
#>  9 85    114  
#> 10 95    136
```

### Highlight searches

Search for and highlight the term _alcohol_ in the abstract field only


```r
(out <- highplos(q='alcohol', hl.fl = 'abstract', rows=3))
#> $`10.1371/journal.pone.0218147`
#> $`10.1371/journal.pone.0218147`$abstract
#> [1] "Background: Binge drinking, an increasingly common form of <em>alcohol</em> use disorder, is associated"
#> 
#> 
#> $`10.1371/journal.pone.0138021`
#> $`10.1371/journal.pone.0138021`$abstract
#> [1] "Background and Aim: Harmful <em>alcohol</em> consumption has long been recognized as being the major"
#> 
#> 
#> $`10.1371/journal.pone.0201042`
#> $`10.1371/journal.pone.0201042`$abstract
#> [1] "\nAcute <em>alcohol</em> administration can lead to a loss of control over drinking. Several models argue"
```

And you can browse the results in your default browser


```r
highbrow(out)
```

![highbrow](man/figures/highbrow.png)

### Full text urls

Simple function to get full text urls for a DOI


```r
full_text_urls(doi='10.1371/journal.pone.0086169')
#> [1] "http://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0086169&type=manuscript"
```

### Full text xml given a DOI


```r
(out <- plos_fulltext(doi='10.1371/journal.pone.0086169'))
#> 1 full-text articles retrieved 
#> Min. Length: 110717 - Max. Length: 110717 
#> DOIs: 10.1371/journal.pone.0086169 ... 
#> 
#> NOTE: extract xml strings like output['<doi>']
```

Then parse the XML any way you like, here getting the abstract


```r
library("XML")
xpathSApply(xmlParse(out$`10.1371/journal.pone.0086169`), "//abstract", xmlValue)
#> [1] "Mammalian females pay high energetic costs for reproduction, the greatest of which is imposed by lactation. The synthesis of milk requires, in part, the mobilization of bodily reserves to nourish developing young. Numerous hypotheses have been advanced to predict how mothers will differentially invest in sons and daughters, however few studies have addressed sex-biased milk synthesis. Here we leverage the dairy cow model to investigate such phenomena. Using 2.39 million lactation records from 1.49 million dairy cows, we demonstrate that the sex of the fetus influences the capacity of the mammary gland to synthesize milk during lactation. Cows favor daughters, producing significantly more milk for daughters than for sons across lactation. Using a sub-sample of this dataset (N = 113,750 subjects) we further demonstrate that the effects of fetal sex interact dynamically across parities, whereby the sex of the fetus being gestated can enhance or diminish the production of milk during an established lactation. Moreover the sex of the fetus gestated on the first parity has persistent consequences for milk synthesis on the subsequent parity. Specifically, gestation of a daughter on the first parity increases milk production by ∼445 kg over the first two lactations. Our results identify a dramatic and sustained programming of mammary function by offspring in utero. Nutritional and endocrine conditions in utero are known to have pronounced and long-term effects on progeny, but the ways in which the progeny has sustained physiological effects on the dam have received little attention to date."
```

### Search within a field

There are a series of convience functions for searching within sections of articles.

* `plosauthor()`
* `plosabstract()`
* `plosfigtabcaps()`
* `plostitle()`
* `plossubject()`

For example:


```r
plossubject(q='marine ecology',  fl = c('id','journal'), limit = 10)
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1     2287     0
#> 
#> $data
#> # A tibble: 10 x 2
#>    id                                                  journal 
#>    <chr>                                               <chr>   
#>  1 10.1371/journal.pone.0092590                        PLoS ONE
#>  2 10.1371/journal.pone.0092590/title                  PLoS ONE
#>  3 10.1371/journal.pone.0092590/abstract               PLoS ONE
#>  4 10.1371/journal.pone.0092590/references             PLoS ONE
#>  5 10.1371/journal.pone.0092590/body                   PLoS ONE
#>  6 10.1371/journal.pone.0092590/introduction           PLoS ONE
#>  7 10.1371/journal.pone.0092590/results_and_discussion PLoS ONE
#>  8 10.1371/journal.pone.0092590/materials_and_methods  PLoS ONE
#>  9 10.1371/journal.pone.0092590/conclusions            PLoS ONE
#> 10 10.1371/journal.pone.0149852                        PLOS ONE
```

However, you can always just do this in `searchplos()` like `searchplos(q = "subject:science")`. See also the `fq` parameter. The above convenience functions are simply wrappers around `searchplos`, so take all the same parameters.

### Search by article views

Search with term _marine ecology_, by field _subject_, and limit to 5 results


```r
plosviews(search='marine ecology', byfield='subject', limit=5)
#>                             id counter_total_all
#> 5 10.1371/journal.pone.0012946              5309
#> 3 10.1371/journal.pone.0167128              5423
#> 1 10.1371/journal.pone.0092590             12839
#> 2 10.1371/journal.pone.0149852             18954
#> 4 10.1371/journal.pone.0011372             28244
```

### Visualize

Visualize word use across articles


```r
plosword(list('monkey','Helianthus','sunflower','protein','whale'), vis = 'TRUE')
#> $table
#>   No_Articles       Term
#> 1       14528     monkey
#> 2         646 Helianthus
#> 3        1876  sunflower
#> 4      164537    protein
#> 5        2142      whale
#> 
#> $plot
```

![wordusage](man/figures/unnamed-chunk-21-1.png)

### progress bars


```r
res <- searchplos(q='*:*', limit = 2000, progress = httr::progress())
#>  |=====================================| 100%
#>  |=====================================| 100%
#>  |=====================================| 100%
#>  |=====================================| 100%
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rplos/issues).
* License: MIT
* Get citation information for `rplos` in R doing `citation(package = 'rplos')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

---

This package is part of a richer suite called [fulltext](https://github.com/ropensci/fulltext), along with several other packages, that provides the ability to search for and retrieve full text of open access scholarly articles. We recommend using `fulltext` as the primary R interface to `rplos` unless your needs are limited to this single source.

---
rplos 1.0
=========

## MINOR IMPROVEMENTS

* fix a broken test on cran (#128)


rplos 0.9.0
===========

## MINOR IMPROVEMENTS

* functions that use solrium under the hood now have a `progress` parameter that you can pass `htt::progress()` to get progress information; especially useful for long running queries  (#124)
* move readme images to man/figures (#127)
* replace `dplyr::data_frame` with `dplyr::tibble` (#126)


rplos 0.8.6
============

## MINOR IMPROVEMENTS

* use `preserve_exact_body_bytes` for tests for plosabstract and plosfigtabcaps to avoid non-ascii text problems on debian clang devel (#125)


rplos 0.8.4
============

## MINOR IMPROVEMENTS

* update docs for `searchplos()` and all wrapper fxns to explain that internal pagination is used, but that users can do their own pagination if they like (#122)

## BUG FIXES

* fix to pagination in `searchplos()` and all wrapper fxns. large numbers were being passed as scientific notation, fixed now (#123)


rplos 0.8.2
===========

## NEW FEATURES

* Integration with `vcr` and `webmockr` packages for unit test stubbing

## BUG FIXES

* for `highbrow()` open pages with `https://doi.org` instead of `http://dx.doi.org` (#117)
* remove message `"Looping - printing progress ..."` from `searchplos()` (#120)
* fix internal pagination for `searchplos()`: were accidentally dropping `fq` statements if more than 1, woopsy  (#121)


rplos 0.8.0
===========

## NEW FEATURES

* Now using `solrium` for under the hood Solr interaction instead
of `solr` package (#106)
* Along with above change, the following: `facetplos`, `searchplos`,
and `highplos` lose parameter `verbose`, and gain parameters
`error` and `proxy` for changing how verbose error reporting is, and
for setting proxy details, respectively.
* Now using `crul` instead of `httr` for HTTP requests (#110)

## MINOR IMPROVEMENTS

* Fix to placement of images for README requested by CRAN (#114)
* Replaced `XML` with `xml2` (#112)
* `citations` function for PLOS rich citations is defunct as the
service is gone (#113)
* package `tm` dropped from Enhances (#111)
* added code of conduct, issue and pull request templates


rplos 0.6.4
===========

## MINOR IMPROVEMENTS

* URLs to full text XML have been changed - old URLs were working
but were going through 2 302 redirects to get there. Updated URLs.
(#107)

## BUG FIXES

* Fixed `content-type` check for `plos_fulltext()` function. XML
can be either `application/xml` or `text/xml` (#108)


rplos 0.6.0
===========

## MINOR IMPROVEMENTS

* Added notes to documentation for relavant functions for how to do
phrase searching. (#96) (#97) thanks @poldham
* Removed parameter `random` parameter from `citations()` function as it's
no longer available in the API (#103)
* Swapped out all uses of `dplyr::rbind_all()` for `dplyr::bind_rows()` (#105)
* `full_text_urls()` now gives back `NA` when DOIs for annotations are
given, which can be easily removed.

## BUG FIXES

* Fixed `full_text_urls()` function to create full text URLs for PLOS
Clinical Trials correctly (#104)

rplos 0.5.6
============

## MINOR IMPROVEMENTS

* move `ggplot2` from _Depends_ to _Imports_, and using `@importFrom` for
`ggplot2` functions, now all imports are using `@importFrom` (#99)
* Fixes for `httr::content()` to parse manually, and use explicit
encoding of `UTF-8` (#102)

rplos 0.5.4
===========

## MINOR IMPROVEMENTS

* Change `solr` dependency to require version `v0.1.6` or less (#94)

rplos 0.5.2
===========

## MINOR IMPROVEMENTS

* More tests added (#94)

## BUG FIXES

* Fix encoding in parsing of XML data in `plos_fulltext()` to
avoid unicode problems (#93)

rplos 0.5.0
===========

## MINOR IMPROVEMENTS

* Now importing non-Base R functions from `utils`, `stats`, and `methods` packages (#90)

## BUG FIXES

* Fixes for `httr` `v1` that broke `rplos` when length 0 list passed to `query` parameter (#89)


rplos 0.4.7
===========

## NEW FEATURES

* New function `citations()` for querying the PLOS Rich Citations API (http://api.richcitations.org/) (#88)

## BUG FIXES

* Added `vignettes/figure` to `.Rbuildignore` as requested by CRAN admin (#87)

rplos 0.4.6
===========

## NEW FEATURES

* API key no longer required (#86)

## MINOR IMPROVEMENTS

* `searchplos()` now returns a list of length two, `meta` and `data`, and `meta` is a data.frame of metadata for the search.
* Switched from CC0 to MIT license.
* No longer importing libraries `RCurl`, `data.table`, `googleVis`, `assertthat`, `RJSONIO`, and `stringr` (#79) (#82) (#84)
* Now importing `dplyr`.
* Moved `jsonlite` from Suggests to Imports. Replaces use of `RJSONIO`. (#80)
* `crossref()` now defunct. See package `rcrossref` https://github.com/ropensci/rcrossref. (#83)
* `highplos()` now uses `solr::solr_highlight()` to do highlight searches.
* `searchplos()`, `plosabstract()`, and other functions that wrap `searchplos()` now use `...` to pass in curl options to `httr::GET()`. You'll now get an error on using `callopts` parameter.
* Added manual file entry for the dataset `isocodes`.
* Reworked both `plosword()` and `plot_throughtime()` to have far less code, uses `httr` now instead of `RCurl`, but to the user, everything should be the same.
* Made documentation more clear on discrepancy between PLOS website behavior and `rplos` behavior, and how to make them match, or match more closely (#76)
* Added package level man file to allow `?rplos` to go to help page.

## BUG FIXES

* Removed some examples from `searchplos()` that are now not working for some unknown reason. (#81)
* Previously when user set `limit=0`, we still gave back data, this is fixed, and now the `meta` slot given back, and the `data` slot gives an `NA` (#85)

rplos 0.4.1
===========

## BUG FIXES

* Fixed some broken tests.

rplos 0.4.0
===========

## NEW FEATURES

* Errors from the data provider are reported now. At least we attempt to do so when they are given, for example if you specify `asc` or `desc` incorrectly with the `sort` parameter. See the `check_response()` function https://github.com/ropensci/rplos/blob/master/tests/testthat/test-check_response.R for examples.
* New functions `facetplos()` and `highplos()` using the `solr` R wrapper to the Solr indexing engine. The PLOS API just exposes the Solr endpoints, so we can use the general Solr wrapper package `solr` to allow more flexible Solr searching.
* New function `highbrow()` to visualize highlighting results in a browser.
* New function `plos_fulltext()` to get full text xml of PLOS articles. Helper function `full_text_urls()` constructs the URL's for full text xml.

## BUG FIXES

* Fixed bug in tests where we forgot to give a key. No key is required per se, but PLOS encourages it so we prevent a call from happening without at least a dumby key.
* Added function `check_response()` to check responses from the PLOS API, deals with capturing server error messages, and checking for correct content type, etc.

## IMPROVEMENTS

* Removed function `crossref_r()` as we are working on a package for the CrossRef API.
* Parameter arguments in `searchplos()`, `plosauthor()`, `plosfigtabcaps()`, `plossubject()`, and `plostitle()` were changed to match closer the Solr parameter names. `terms` to `q`. `fields` to `fl`. `toquery` to `fq`.
* Multiple values passed to `fields `
* `returndf` parameter is gone from `searchplos()`, `plosauthor()`, `plosfigtabcaps()`, `plossubject()`, and `plostitle()`. You can easily get raw JSON, etc. data using the `solr` package.
* Now using `httr` instead of `RCurl` in `plosviews()` function.

rplos 0.3.6
===========

## NEW FEATURES

* All search functions (searchplos(), plosabstract(), plosauthor(), plosfigtabcaps(), plossubject(), and plostitle()) gain highlighting argument, setting to TRUE (default=FALSE) returns matching sentence fragments that were matched. NOTE that if highlighting=TRUE the output can be a list of data.frame's if returndf=TRUE, with two named elements 'data' and 'highlighting', or a list of lists if returndf=FALSE.
* All search functions (searchplos(), plosabstract(), plosauthor(), plosfigtabcaps(), plossubject(), and plostitle()) gain sort argument. You can pass a field to sort by even if you don't return that field in the output, e.g., sort='counter_total_month desc'.
* A tiny function parsehighlight() added to parse out html code from highlighting output.

## BUG FIXES

* Some examples in docs didn't work - fixed them.
* Fixed bug in searchplos() that was causing elements of a return field to cause failure because they were longer than 1 (e.g., authors). Now concatenating elements of length > 1.
* Fixed bug in searchplos() that was causing elements of length 0 to cause failure. Now removing elements of length 0.
* Fixed parsehighlight function to return NA if highlighting return of length 0.
* Fixed broken test for plosauthor(), plosabstract(), and plot_throughtime().

rplos 0.3.0
===========

## NEW FEATURES

* Added httr::stop_for_status() calls to a few functions to give informative http status errors when they happen

## BUG FIXES

* Fixed bug in plot_throughtime() that was throwing errors and preventing fxn from working, thanks to Ben Bolker for the fix.
* Simplified code in many functions to have cleaner and simpler code.
* ... parameter in many functions changed to callopts=list(), which passes in curl options to a call to either RCurl::getForm() or httr::GET()
* Fixed bug in function plosviews() that caused errors in some calls. Now forces full document searches, so that you get views data back for full papers only, not sections of papers. See package alm (https://github.com/ropensci/alm) for more in depth PLOS article-level metrics.


rplos 0.2.0
===========

## NEW FEATURES

* All functions for interacting with the PLOS ALM (altmetrics) API have been removed, and are now in a separate package called alm (http://github.com/ropensci/alm).
* Convenience functions `plosabstract`, `plosauthor`, `plosfigtabcaps`, `plossubject`, and `plostitle`, that search specifically within those sections of papers now wrap `searchplos`, so they should behave the same way.
* ldfast() fxn added as an attempt to do ldply faster
* performance improvements in searchplos

## BUG FIXES

* Dependency on assertthat removed since it's not on CRAN.
* Fixed namespace conflicts by importing only functions needed from some packages.
* searchplos() now removes leading, trailing, and internal whitespace from character strings


rplos 0.1.1
===========

* remove alm*() functions so that this package now only wraps the PLoS search API.


rplos 0.1.0
===========

* The `almdateupdated` function has been deprecated - use `almupdated` instead.

* The `articlelength` function has been deprecated - didn't see the usefulness any longer.

* In general simplified and prettified code.

* Changed from using RCurl to httr in many functions, but not all.

* Added more examples for many functions.

* Added three internal functions: `concat_todf`, `addmissing`, and `getkey`.

* Added Karthik Ram as a package author.

## BUG FIXES

* All `url` arguments in functions put inside functions as they are not likely to change that often.

* Fixed `crossref` function, and added more examples.

## NEW FEATURES

* The `alm` function (previously `almplosallviews`) gains many ### new features: now allows up to 50 DOIs per call; you can specify the source you want to get alm data from as an argument; you can specify the year you want to get alm data from as an argument.

* Added the plosfields data file to get all the possible fields to use in function calls.


## NEW FUNCTIONS

* `almplosallviews` changed to `alm`.

* `almplotallviews` changed to `almplot`.

* `almevents` added to specifically search and get detailed events data for a specific source or N sources.

* `crossref_r` gets 20 random DOIs from Crossref.org.

* Added package startup message.

* `journalnamekey` function to get the short name keys for each PLoS Journal.


rplos 0.0-7
===========

## IMPROVEMENTS AND ### BUG FIXES

* ALM functions (any functions starting with alm) received updated arguments/parameters according to the ALM API version 3.0 changes.

* ### Bug fixes in general across library.

* Added tests.

* `almplosallviews` now outputs different output - two data.frames, one total metrics (summed across time), and history (for metrics for each time period specified in the search)

* `crossref` function returns R's native bibtype format.  See examples in `crossref` function documentation

rplos 0.0-5
===========

## IMPROVEMENTS AND BUG FIXES

* `almpub` changed to `almdatepub`

* changed help file `rplos` to `help` - use help('rplos') in R

* changed URL from http://ropensci.org/ to https://github.com/ropensci/rplos

* added sleep argument to `plosallviews` function to allow pauses between API calls when running `plosallviews` in a loop - this is an attempt to limit hitting the PLoS API too hard

* various other fixed to functions

* more examples added to some functions

## NEW FUNCTIONS

* added function `journalnamekey` to get short keys for journals to use in searching for specific journals


rplos 0.0-1
===============

## NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.0.4 Patched
* ubuntu 16.04 (on github actions), R 4.0.4
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 downstream dependency, with no errors. Summary at <https://github.com/ropensci/rplos/blob/master/revdep/README.md>

----------

This version fixes a broken test and future proofs all tests so they are skipped if the remote service is offline.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rplos/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rplos.git`
* Make sure to track progress upstream (i.e., on our version of `rplos` at `ropensci/rplos`) by doing `git remote add upstream https://github.com/ropensci/rplos.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/rplos`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.6.3 Patched (2020-03-11 r78147) |
|os       |macOS Catalina 10.15.4                      |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-04-07                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|rplos   |0.8.6 |0.9.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*