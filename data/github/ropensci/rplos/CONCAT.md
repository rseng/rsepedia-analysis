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

*Wow, no problems at all. :)**Wow, no problems at all. :)*rplos
=====

```{r echo=FALSE}
knitr::opts_chunk$set(
  fig.path = "man/figures/",
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.cap = ""
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rplos)](https://cranchecks.info/pkgs/rplos)
[![R-check](https://github.com/ropensci/rplos/workflows/R-check/badge.svg)](https://github.com/ropensci/rplos/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/rplos/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rplos?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rplos)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rplos)](https://cran.r-project.org/package=rplos)

## Install

You can get this package at CRAN [here](https://cran.r-project.org/package=rplos), or install it within R by doing

```{r eval=FALSE}
install.packages("rplos")
```

Or install the development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/rplos")
```

```{r}
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

```{r}
searchplos('ecology', 'id,publication_date', limit = 5)
```

Get DOIs for full article in PLoS One

```{r}
searchplos(q="*:*", fl='id', fq=list('journal_key:PLoSONE',
   'doc_type:full'), limit=5)
```

Query to get some PLOS article-level metrics, notice difference between two outputs

```{r}
out <- searchplos(q="*:*", fl=c('id','counter_total_all','alm_twitterCount'), fq='doc_type:full')
out_sorted <- searchplos(q="*:*", fl=c('id','counter_total_all','alm_twitterCount'),
   fq='doc_type:full', sort='counter_total_all desc')
head(out$data)
head(out_sorted$data)
```

A list of articles about social networks that are popular on a social network

```{r}
searchplos(q="*:*",fl=c('id','alm_twitterCount'),
   fq=list('doc_type:full','subject:"Social networks"','alm_twitterCount:[100 TO 10000]'),
   sort='counter_total_month desc')
```

Show all articles that have these two words less then about 15 words apart

```{r}
searchplos(q='everything:"sports alcohol"~15', fl='title', fq='doc_type:full', limit=3)
```

Narrow results to 7 words apart, changing the ~15 to ~7

```{r}
searchplos(q='everything:"sports alcohol"~7', fl='title', fq='doc_type:full', limit=3)
```

Remove DOIs for annotations (i.e., corrections) and Viewpoints articles

```{r}
searchplos(q='*:*', fl=c('id','article_type'),
   fq=list('-article_type:correction','-article_type:viewpoints'), limit=5)
```

### Faceted search

Facet on multiple fields

```{r}
facetplos(q='alcohol', facet.field=c('journal','subject'), facet.limit=5)
```

Range faceting

```{r}
facetplos(q='*:*', url=url, facet.range='counter_total_all',
 facet.range.start=5, facet.range.end=100, facet.range.gap=10)
```

### Highlight searches

Search for and highlight the term _alcohol_ in the abstract field only

```{r}
(out <- highplos(q='alcohol', hl.fl = 'abstract', rows=3))
```

And you can browse the results in your default browser

```{r eval=FALSE}
highbrow(out)
```

![highbrow](man/figures/highbrow.png)

### Full text urls

Simple function to get full text urls for a DOI

```{r}
full_text_urls(doi='10.1371/journal.pone.0086169')
```

### Full text xml given a DOI

```{r}
(out <- plos_fulltext(doi='10.1371/journal.pone.0086169'))
```

Then parse the XML any way you like, here getting the abstract

```{r}
library("XML")
xpathSApply(xmlParse(out$`10.1371/journal.pone.0086169`), "//abstract", xmlValue)
```

### Search within a field

There are a series of convience functions for searching within sections of articles.

* `plosauthor()`
* `plosabstract()`
* `plosfigtabcaps()`
* `plostitle()`
* `plossubject()`

For example:

```{r}
plossubject(q='marine ecology',  fl = c('id','journal'), limit = 10)
```

However, you can always just do this in `searchplos()` like `searchplos(q = "subject:science")`. See also the `fq` parameter. The above convenience functions are simply wrappers around `searchplos`, so take all the same parameters.

### Search by article views

Search with term _marine ecology_, by field _subject_, and limit to 5 results

```{r}
plosviews(search='marine ecology', byfield='subject', limit=5)
```

### Visualize

Visualize word use across articles

```{r fig.cap="wordusage"}
plosword(list('monkey','Helianthus','sunflower','protein','whale'), vis = 'TRUE')
```

### progress bars

```{r eval = FALSE}
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
---
title: Introduction to rplos
author: Scott Chamberlain
date: "2021-02-23"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Introduction to rplos}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



The `rplos` package interacts with the API services of [PLoS](https://plos.org/)
(Public Library of Science) Journals. You used to need an API key to work with
this package - that is no longer needed!

This tutorial will go through three use cases to demonstrate the kinds
of things possible in `rplos`.

* Search across PLoS papers in various sections of papers
* Search for terms and visualize results as a histogram OR as a plot through
time
* Text mining of scientific literature

### Load package from CRAN


```r
install.packages("rplos")
```


```r
library('rplos')
```

### Search across PLoS papers in various sections of papers

`searchplos` is a general search, and in this case searches for the term
**Helianthus** and returns the DOI's of matching papers


```r
searchplos(q = "Helianthus", fl = "id", limit = 5)
```

```
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1      646     0
#> 
#> $data
#> # A tibble: 5 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pone.0198869
#> 2 10.1371/journal.pone.0213065
#> 3 10.1371/journal.pone.0148280
#> 4 10.1371/journal.pone.0111982
#> 5 10.1371/journal.pone.0212371
```

Get only full article DOIs


```r
searchplos(q = "*:*", fl = 'id', fq = 'doc_type:full', start = 0, limit = 5)
```

```
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1   292789     0
#> 
#> $data
#> # A tibble: 5 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pcbi.0020071
#> 2 10.1371/journal.pbio.1000152
#> 3 10.1371/journal.pbio.1000153
#> 4 10.1371/journal.pbio.1000159
#> 5 10.1371/journal.pbio.1000165
```

Get DOIs for only PLoS One articles


```r
searchplos(q = "*:*", fl = 'id', fq = 'journal_key:PLoSONE',
           start = 0, limit = 5)
```

```
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1  2125942     0
#> 
#> $data
#> # A tibble: 5 x 1
#>   id                                       
#>   <chr>                                    
#> 1 10.1371/journal.pone.0002397/title       
#> 2 10.1371/journal.pone.0002397/abstract    
#> 3 10.1371/journal.pone.0002397/references  
#> 4 10.1371/journal.pone.0002397/body        
#> 5 10.1371/journal.pone.0002397/introduction
```

Get DOIs for full article in PLoS One


```r
searchplos(q = "*:*", fl = 'id',
   fq = list('journal_key:PLoSONE', 'doc_type:full'),
   start = 0, limit = 5)
```

```
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

Search for many terms


```r
q <- c('ecology','evolution','science')
lapply(q, function(x) searchplos(x, limit = 2))
```

```
#> [[1]]
#> [[1]]$meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1    55873     0
#> 
#> [[1]]$data
#> # A tibble: 2 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pone.0001248
#> 2 10.1371/journal.pone.0059813
#> 
#> 
#> [[2]]
#> [[2]]$meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1    82168     0
#> 
#> [[2]]$data
#> # A tibble: 2 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pbio.2002255
#> 2 10.1371/journal.pone.0205798
#> 
#> 
#> [[3]]
#> [[3]]$meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1   260220     0
#> 
#> [[3]]$data
#> # A tibble: 2 x 1
#>   id                          
#>   <chr>                       
#> 1 10.1371/journal.pone.0229237
#> 2 10.1371/journal.pone.0202320
```

### Search on specific sections

A suite of functions were created as light wrappers around `searchplos` as
a shorthand to search specific sections of a paper.

* `plosauthor` searchers in authors
* `plosabstract` searches in abstracts
* `plostitle` searches in titles
* `plosfigtabcaps` searches in figure and table captions
* `plossubject` searches in subject areas

`plosauthor` searches across authors, and in this case returns the authors of
the matching papers. the fl parameter determines what is returned


```r
plosauthor(q = "Eisen", fl = "author", limit = 5)
```

```
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1     1107     0
#> 
#> $data
#> # A tibble: 5 x 1
#>   author                                                                       
#>   <chr>                                                                        
#> 1 Myungsun Kang,Timothy J Eisen,Ellen A Eisen,Arup K Chakraborty,Herman N Eisen
#> 2 Myungsun Kang,Timothy J Eisen,Ellen A Eisen,Arup K Chakraborty,Herman N Eisen
#> 3 Myungsun Kang,Timothy J Eisen,Ellen A Eisen,Arup K Chakraborty,Herman N Eisen
#> 4 Myungsun Kang,Timothy J Eisen,Ellen A Eisen,Arup K Chakraborty,Herman N Eisen
#> 5 Myungsun Kang,Timothy J Eisen,Ellen A Eisen,Arup K Chakraborty,Herman N Eisen
```

`plosabstract` searches across abstracts, and in this case returns the id and
title of the matching papers


```r
plosabstract(q = 'drosophila', fl = 'id,title', limit = 5)
```

```
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1     3944     0
#> 
#> $data
#> # A tibble: 5 x 2
#>   id                     title                                                  
#>   <chr>                  <chr>                                                  
#> 1 10.1371/journal.pone.… Host Range and Specificity of the Drosophila C Virus   
#> 2 10.1371/journal.pgen.… Drosophila Myc restores immune homeostasis of Imd path…
#> 3 10.1371/journal.pone.… Exogenous expression of Drp1 plays neuroprotective rol…
#> 4 10.1371/journal.pone.… A Drosophila model for developmental nicotine exposure 
#> 5 10.1371/image.pbio.v0… PLoS Biology Issue Image | Vol. 6(5) May 2008
```

`plostitle` searches across titles, and in this case returns the title and
journal of the matching papers


```r
plostitle(q = 'drosophila', fl = 'title,journal', limit = 5)
```

```
#> $meta
#> # A tibble: 1 x 2
#>   numFound start
#>      <int> <int>
#> 1     2480     0
#> 
#> $data
#> # A tibble: 5 x 2
#>   journal  title                                                                
#>   <chr>    <chr>                                                                
#> 1 PLOS ONE Tandem Duplications and the Limits of Natural Selection in Drosophil…
#> 2 PLoS ONE A DNA Virus of Drosophila                                            
#> 3 PLOS ONE Nematocytes: Discovery and characterization of a novel anculeate hem…
#> 4 PLOS ONE Peptidergic control in a fruit crop pest: The spotted-wing drosophil…
#> 5 PLoS ONE In Vivo RNAi Rescue in Drosophila melanogaster with Genomic Transgen…
```

### Search terms & visualize results as a histogram OR as a plot through time

`plosword` allows you to search for 1 to K words and visualize the results
as a histogram, comparing number of matching papers for each word


```r
out <- plosword(list("monkey", "Helianthus", "sunflower", "protein", "whale"),
    vis = "TRUE")
out$table
```

```
#>   No_Articles       Term
#> 1       14528     monkey
#> 2         646 Helianthus
#> 3        1876  sunflower
#> 4      164537    protein
#> 5        2142      whale
```


```r
out$plot
```

![plot of chunk plosword1plot](../man/figures/plosword1plot-1.png)

You can also pass in curl options, in this case get verbose information on the
curl call.


```r
plosword('Helianthus', callopts = list(verbose = TRUE))
```

```
#> Number of articles with search term 
#>                                 646
```

### Visualize terms

`plot_throughtime` allows you to search for up to 2 words and visualize the
results as a line plot through time, comparing number of articles matching
through time. Visualize with the ggplot2 package, only up to two terms for now.


```r
library("ggplot2")
plot_throughtime(terms = "phylogeny", limit = 200) +
  geom_line(size = 2, color = 'black')
```

![plot of chunk throughtime1](../man/figures/throughtime1-1.png)

### More

See the _Faceted and highlighted searches_ and _Full text_ vignettes for
more `rplos` help.
---
title: Faceted and highlighted searches
author: Scott Chamberlain
date: "2021-02-23"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Faceted and highlighted searches}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



In addition to `searchplos()` and related searching functions, there are a 
few slightly different ways to search: faceting and highlighted searches. 
Faceting allows you to ask, e.g., how many articles are published in each of 
the PLOS journals. Highlighting allows you to ask, e.g., highlight terms that 
I search for in the text results given back, which can make downstream 
processing easier, and help visualize search results (see `highbrow()` below). 

### Load package from CRAN


```r
install.packages("rplos")
```


```r
library('rplos')
```

### Faceted search

Facet by journal


```r
facetplos(q='*:*', facet.field='journal')
```

```
#> $facet_queries
#> NULL
#> 
#> $facet_fields
#> $facet_fields$journal
#> # A tibble: 9 x 2
#>   term                             value  
#>   <chr>                            <chr>  
#> 1 plos one                         2073767
#> 2 plos genetics                    77191  
#> 3 plos neglected tropical diseases 72937  
#> 4 plos pathogens                   71789  
#> 5 plos computational biology       66404  
#> 6 plos biology                     44593  
#> 7 plos medicine                    32612  
#> 8 plos clinical trials             521    
#> 9 plos medicin                     9      
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

Using `facet.query` to get counts


```r
facetplos(q='*:*', facet.field='journal', facet.query='cell,bird')
```

```
#> $facet_queries
#> # A tibble: 1 x 2
#>   term      value
#>   <chr>     <int>
#> 1 cell,bird 12052
#> 
#> $facet_fields
#> $facet_fields$journal
#> # A tibble: 9 x 2
#>   term                             value  
#>   <chr>                            <chr>  
#> 1 plos one                         2073767
#> 2 plos genetics                    77191  
#> 3 plos neglected tropical diseases 72937  
#> 4 plos pathogens                   71789  
#> 5 plos computational biology       66404  
#> 6 plos biology                     44593  
#> 7 plos medicine                    32612  
#> 8 plos clinical trials             521    
#> 9 plos medicin                     9      
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

Date faceting


```r
facetplos(q='*:*', url=url, facet.date='publication_date',
  facet.date.start='NOW/DAY-5DAYS', facet.date.end='NOW', 
  facet.date.gap='+1DAY')
```

```
#> list()
```

### Highlighted search

Search for the term _alcohol_ in the abstracts of articles, return only 
10 results


```r
highplos(q='alcohol', hl.fl = 'abstract', rows=2)
```

```
#> $`10.1371/journal.pone.0218147`
#> $`10.1371/journal.pone.0218147`$abstract
#> [1] "Background: Binge drinking, an increasingly common form of <em>alcohol</em> use disorder, is associated"
#> 
#> 
#> $`10.1371/journal.pone.0138021`
#> $`10.1371/journal.pone.0138021`$abstract
#> [1] "Background and Aim: Harmful <em>alcohol</em> consumption has long been recognized as being the major"
```

Search for the term _alcohol_ in the abstracts of articles, and return 
fragment size of 20 characters, return only 5 results


```r
highplos(q='alcohol', hl.fl='abstract', hl.fragsize=20, rows=2)
```

```
#> $`10.1371/journal.pone.0218147`
#> $`10.1371/journal.pone.0218147`$abstract
#> [1] " common form of <em>alcohol</em>"
#> 
#> 
#> $`10.1371/journal.pone.0138021`
#> $`10.1371/journal.pone.0138021`$abstract
#> [1] ": Harmful <em>alcohol</em>"
```

Search for the term _experiment_ across all sections of an article, return 
id (DOI) and title fl only, search in full articles only 
(via `fq='doc_type:full'`), and return only 10 results


```r
highplos(q='everything:"experiment"', fl='id,title', fq='doc_type:full',
   rows=2)
```

```
#> $`10.1371/journal.pone.0154334`
#> $`10.1371/journal.pone.0154334`$everything
#> [1] " and designed the <em>experiments</em>: RJ CM AOC. Performed the <em>experiments</em>: RJ AOC. Analyzed the data: RJ. Contributed"
#> 
#> 
#> $`10.1371/journal.pone.0039681`
#> $`10.1371/journal.pone.0039681`$everything
#> [1] " Selection of Transcriptomics <em>Experiments</em> Improves Guilt-by-Association Analyses Transcriptomics <em>Experiment</em>"
```

### Visualize highligted searches

Browse highlighted fragments in your default browser

This first examle, we only looko at 10 results


```r
out <- highplos(q='alcohol', hl.fl = 'abstract', rows=10)
highbrow(out)
```

![highbrow1](../man/figures/highbrow.png)

But it works quickly with lots of results too


```r
out <- highplos(q='alcohol', hl.fl = 'abstract', rows=1200)
highbrow(out)
```

![highbrow2](../man/figures/highbrow_big.png)
---
title: Full text
author: Scott Chamberlain
date: "2021-02-23"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{Full text}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



Search functions in `rplos` can be used to get back full text in addition to 
any section of an article. However, if you prefer XML, this vignette is 
for you.

### Load package from CRAN


```r
install.packages("rplos")
```


```r
library('rplos')
```

### Get full text URLs

Create urls for full text articles in PLOS journals

Here's the URL for XML full text for the DOI `10.1371/journal.pone.0086169`


```r
full_text_urls(doi = '10.1371/journal.pone.0086169')
```

```
#> [1] "http://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0086169&type=manuscript"
```

And for the DOI `10.1371/journal.pbio.1001845`


```r
full_text_urls(doi = '10.1371/journal.pbio.1001845')
```

```
#> [1] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1001845&type=manuscript"
```

The function is vectorized, so you can pass in many DOIs


```r
full_text_urls(doi = c('10.1371/journal.pone.0086169', 
                       '10.1371/journal.pbio.1001845'))
```

```
#> [1] "http://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0086169&type=manuscript"    
#> [2] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1001845&type=manuscript"
```

Use `searchplos()` to get a lot of DOIs, then get the URLs for full text XML


```r
dois <- searchplos(q = "*:*", fq = 'doc_type:full', limit = 20)$data$id
full_text_urls(dois)
```

```
#>  [1] "http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.0020071&type=manuscript"
#>  [2] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000152&type=manuscript" 
#>  [3] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000153&type=manuscript" 
#>  [4] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000159&type=manuscript" 
#>  [5] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000165&type=manuscript" 
#>  [6] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000166&type=manuscript" 
#>  [7] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000167&type=manuscript" 
#>  [8] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000173&type=manuscript" 
#>  [9] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000175&type=manuscript" 
#> [10] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000176&type=manuscript" 
#> [11] "http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.0030121&type=manuscript"
#> [12] "http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.0030124&type=manuscript"
#> [13] "http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.0030125&type=manuscript"
#> [14] "http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.0030127&type=manuscript"
#> [15] "http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.0030130&type=manuscript"
#> [16] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000409&type=manuscript" 
#> [17] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000410&type=manuscript" 
#> [18] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000411&type=manuscript" 
#> [19] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000414&type=manuscript" 
#> [20] "http://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.1000421&type=manuscript"
```

### Get XML

Get full text XML of PLOS papers given a DOI


```r
plos_fulltext(doi = '10.1371/journal.pone.0086169')
```

```
#> 1 full-text articles retrieved 
#> Min. Length: 110717 - Max. Length: 110717 
#> DOIs: 10.1371/journal.pone.0086169 ... 
#> 
#> NOTE: extract xml strings like output['<doi>']
```

`plos_fulltext()` is vectorized, so you can pass in more than one DOI


```r
plos_fulltext(c('10.1371/journal.pone.0086169','10.1371/journal.pbio.1001845'))
```

```
#> 2 full-text articles retrieved 
#> Min. Length: 110717 - Max. Length: 143442 
#> DOIs: 10.1371/journal.pone.0086169 10.1371/journal.pbio.1001845 ... 
#> 
#> NOTE: extract xml strings like output['<doi>']
```

Get many DOIs, then index to get the full XML of the one you want 
(output not shown)


```r
dois <- searchplos(q = "*:*", fq = 'doc_type:full', limit = 3)$data$id
out <- plos_fulltext(dois)
xml <- out[dois[1]][[1]]
```

Extract the abstract from the XML


```r
if (requireNamespace("xml2")) {
  library("xml2")
  xml_text(xml_find_all(read_xml(xml), "//abstract"))
}
```

```
#> [1] "With the burgeoning immunological data in the scientific literature, scientists must increasingly rely on Internet resources to inform and enhance their work. Here we provide a brief overview of the adaptive immune response and summaries of immunoinformatics resources, emphasizing those with Web interfaces. These resources include searchable databases of epitopes and immune-related molecules, and analysis tools for T cell and B cell epitope prediction, vaccine design, and protein structure comparisons. There is an agreeable synergy between the growing collections in immune-related databases and the growing sophistication of analysis software; the databases provide the foundation for developing predictive computational tools, which in turn enable more rapid identification of immune responses to populate the databases. Collectively, these resources contribute to improved understanding of immune responses and escape, and evolution of pathogens under immune pressure. The public health implications are vast, including designing vaccines, understanding autoimmune diseases, and defining the correlates of immune protection."
```

Extract reference lists, just give back first few for brevity


```r
if (requireNamespace("xml2")) {
  library("xml2")
  xml_find_all(read_xml(out[[3]]), "//ref-list/ref")[1:2]
}
```

```
#> {xml_nodeset (2)}
#> [1] <ref id="pbio.1000153-Fetz1">\n  <label>1</label>\n  <element-citation pu ...
#> [2] <ref id="pbio.1000153-Wessberg1">\n  <label>2</label>\n  <element-citatio ...
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{insertnones}
\alias{insertnones}
\title{Function to insert "none" character strings where NULL values found to
faciliate combining results}
\usage{
insertnones(x)
}
\description{
Function to insert "none" character strings where NULL values found to
faciliate combining results
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fulltext.R
\name{plos_fulltext}
\alias{plos_fulltext}
\alias{print.plosft}
\title{Get full text xml of PLOS papers given a DOI}
\usage{
plos_fulltext(doi, ...)

\method{print}{plosft}(x, ...)
}
\arguments{
\item{doi}{One or more DOIs}

\item{...}{Curl options passed on to \code{\link[crul]{HttpClient}}}

\item{x}{Input to print method}
}
\value{
Character string of XML.
}
\description{
Get full text xml of PLOS papers given a DOI
}
\examples{
\dontrun{
plos_fulltext(doi='10.1371/journal.pone.0086169')
plos_fulltext(c('10.1371/journal.pone.0086169',
  '10.1371/journal.pbio.1001845'))
dois <- searchplos(q = "*:*", 
  fq = list('doc_type:full', 'article_type:"Research Article"'), 
  limit = 3)$data$id
out <- plos_fulltext(dois)
out[dois[1]]
out[1:2]

# Extract text from the XML strings - xml2 package required
if (requireNamespace("xml2")) {
  library("xml2")
  lapply(out, function(x){
    tmp <- xml2::read_xml(x)
    xml2::xml_find_all(tmp, "//ref-list//ref")
  })
}

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plostitle.R
\name{plostitle}
\alias{plostitle}
\title{Search PLoS Journals titles.}
\usage{
plostitle(
  q = NULL,
  fl = "id",
  fq = NULL,
  sort = NULL,
  start = 0,
  limit = 10,
  sleep = 6,
  errors = "simple",
  proxy = NULL,
  callopts = NULL,
  progress = NULL,
  ...
)
}
\arguments{
\item{q}{Search terms (character). You can search on specific fields by
doing 'field:your query'. For example, a real query on a specific field would
be 'author:Smith'.}

\item{fl}{Fields to return from search (character) [e.g., 'id,title'],
any combination of search fields (see the dataset \code{plosfields})}

\item{fq}{List specific fields to filter the query on (if NA, all queried).
The options for this parameter are the same as those for the fl parameter.
Note that using this parameter doesn't influence the actual query, but is used
to filter the results to a subset of those you want returned. For example,
if you want full articles only, you can do \code{'doc_type:full'}. In another example,
if you want only results from the journal PLOS One, you can do
\code{'journal_key:PLoSONE'}. See \code{\link{journalnamekey}} for journal
abbreviations.}

\item{sort}{Sort results according to a particular field, and specify ascending (asc)
or descending (desc) after a space; see examples. For example, to sort the
counter_total_all field in descending fashion, do sort='counter_total_all desc'}

\item{start}{Record to start at (used in combination with limit when
you need to cycle through more results than the max allowed=1000). See Pagination
below}

\item{limit}{Number of results to return (integer). Setting \code{limit=0} returns only
metadata. See Pagination below}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call to searchplos. However, if you are using searchplos in a loop or
lapply type call, do sleep parameter is used to prevent your IP address from being
blocked. You can only do 10 requests per minute, so one request every 6 seconds is
about right.}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{(list) optional curl options passed to \code{\link[crul]{HttpClient}}}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}.
only supports \code{httr::progress()}}

\item{...}{Additional Solr arguments}
}
\value{
Titles, in addition to any other fields requested in a data.frame.
}
\description{
Search PLoS Journals titles.
}
\details{
Details:
}
\section{Faceting}{

Read more about faceting here: url{http://wiki.apache.org/solr/SimpleFacetParameters}
}

\section{Website vs. API behavior}{

Don't be surprised if queries you perform in a scripting language, like using \code{rplos}
in R, give different results than when searching for articles on the PLOS website. I am
not sure what exact defaults they use on their website. There are a few things to consider.
You can tweak which types of articles are returned: Try using the \code{article_type}
filter in the \code{fq} parameter. For which journal to search, e.g., do
\code{'journal_key:PLoSONE'}. See \code{journalnamekey()} for journal
abbreviations.
}

\section{Phrase searching}{

To search phrases, e.g., \strong{synthetic biology} as a single item, rather than
separate occurrences of \strong{synthetic} and \strong{biology}, simply put double
quotes around the phrase. For example, to search for cases of \strong{synthetic biology},
do \code{searchplos(q = '"synthetic biology"')}.

You can modify phrase searches as well. For example,
\code{searchplos(q = '"synthetic biology" ~ 10')} asks for cases of
\strong{synthetic biology} within 10 words of each other. See examples.
}

\section{Pagination}{

The \code{searchplos} function and the many functions that are wrappers around 
\code{searchplos} all do paginatino internally for you. That is, if you request for 
example, 2000 results, the max you can get in any one request is 1000, so we'll do 
two requests for you. And so on for larger requests. 

You can always do your own paginatino by doing a lapply type call or a for loop 
to cycle through pages of results.
}

\examples{
\dontrun{
plostitle(q='drosophila', fl='title', limit=99)
plostitle(q='drosophila', fl=c('title','journal'), limit=10)
plostitle(q='drosophila',  limit = 5)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/journalnamekey.R
\name{journalnamekey}
\alias{journalnamekey}
\title{Get short keys for journals to use in searching specific journals.}
\usage{
journalnamekey(...)
}
\arguments{
\item{...}{optional curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
(character) journal name keys
}
\description{
Get short keys for journals to use in searching specific journals.
}
\examples{
\dontrun{
journalnamekey()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plossubject.R
\name{plossubject}
\alias{plossubject}
\title{Search PLoS Journals subjects.}
\usage{
plossubject(
  q = NULL,
  fl = "id",
  fq = NULL,
  sort = NULL,
  start = 0,
  limit = 10,
  sleep = 6,
  errors = "simple",
  proxy = NULL,
  callopts = NULL,
  progress = NULL,
  ...
)
}
\arguments{
\item{q}{Search terms (character). You can search on specific fields by
doing 'field:your query'. For example, a real query on a specific field would
be 'author:Smith'.}

\item{fl}{Fields to return from search (character) [e.g., 'id,title'],
any combination of search fields (see the dataset \code{plosfields})}

\item{fq}{List specific fields to filter the query on (if NA, all queried).
The options for this parameter are the same as those for the fl parameter.
Note that using this parameter doesn't influence the actual query, but is used
to filter the results to a subset of those you want returned. For example,
if you want full articles only, you can do \code{'doc_type:full'}. In another example,
if you want only results from the journal PLOS One, you can do
\code{'journal_key:PLoSONE'}. See \code{\link{journalnamekey}} for journal
abbreviations.}

\item{sort}{Sort results according to a particular field, and specify ascending (asc)
or descending (desc) after a space; see examples. For example, to sort the
counter_total_all field in descending fashion, do sort='counter_total_all desc'}

\item{start}{Record to start at (used in combination with limit when
you need to cycle through more results than the max allowed=1000). See Pagination
below}

\item{limit}{Number of results to return (integer). Setting \code{limit=0} returns only
metadata. See Pagination below}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call to searchplos. However, if you are using searchplos in a loop or
lapply type call, do sleep parameter is used to prevent your IP address from being
blocked. You can only do 10 requests per minute, so one request every 6 seconds is
about right.}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{(list) optional curl options passed to \code{\link[crul]{HttpClient}}}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}.
only supports \code{httr::progress()}}

\item{...}{Additional Solr arguments}
}
\value{
Subject content, in addition to any other fields requested in a
   data.frame.
}
\description{
Search PLoS Journals subjects.
}
\details{
Details:

See \url{https://journals.plos.org/plosone/browse} for subject areas.
}
\section{Faceting}{

Read more about faceting here: url{http://wiki.apache.org/solr/SimpleFacetParameters}
}

\section{Website vs. API behavior}{

Don't be surprised if queries you perform in a scripting language, like using \code{rplos}
in R, give different results than when searching for articles on the PLOS website. I am
not sure what exact defaults they use on their website. There are a few things to consider.
You can tweak which types of articles are returned: Try using the \code{article_type}
filter in the \code{fq} parameter. For which journal to search, e.g., do
\code{'journal_key:PLoSONE'}. See \code{journalnamekey()} for journal
abbreviations.
}

\section{Phrase searching}{

To search phrases, e.g., \strong{synthetic biology} as a single item, rather than
separate occurrences of \strong{synthetic} and \strong{biology}, simply put double
quotes around the phrase. For example, to search for cases of \strong{synthetic biology},
do \code{searchplos(q = '"synthetic biology"')}.

You can modify phrase searches as well. For example,
\code{searchplos(q = '"synthetic biology" ~ 10')} asks for cases of
\strong{synthetic biology} within 10 words of each other. See examples.
}

\section{Pagination}{

The \code{searchplos} function and the many functions that are wrappers around 
\code{searchplos} all do paginatino internally for you. That is, if you request for 
example, 2000 results, the max you can get in any one request is 1000, so we'll do 
two requests for you. And so on for larger requests. 

You can always do your own paginatino by doing a lapply type call or a for loop 
to cycle through pages of results.
}

\examples{
\dontrun{
plossubject('marine ecology', limit = 5)
plossubject(q='marine ecology',  fl = c('id','journal','title'), limit = 20)
plossubject(q='marine ecology', fl = c('id','journal'),
   fq='doc_type:full', limit = 9)
plossubject(q='marine ecology', fl = c('id','journal'),
   fq=list('doc_type:full','!article_type_facet:"Issue\%20Image"'),
   limit = 9)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_throughtime.R
\name{plot_throughtime}
\alias{plot_throughtime}
\title{Plot results through time for serach results from PLoS Journals.}
\usage{
plot_throughtime(terms, limit = NA, ...)
}
\arguments{
\item{terms}{search terms (character)}

\item{limit}{number of results to return (integer)}

\item{...}{optional curl options passed to \code{\link[crul]{HttpClient}}}
}
\value{
Number of search results (vis = FALSE), or number of search in 
a table and a histogram of results (vis = TRUE).
}
\description{
Plot results through time for serach results from PLoS Journals.
}
\examples{
\dontrun{
plot_throughtime(terms='phylogeny', limit=300)
plot_throughtime(list('drosophila','monkey'), 100)
plot_throughtime(list('drosophila','flower','dolphin','cell','cloud'), 100)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highbrow.r
\name{highbrow}
\alias{highbrow}
\title{Browse highlighted fragments in your default browser.}
\usage{
highbrow(input = NULL, output = NULL, browse = TRUE)
}
\arguments{
\item{input}{Input, usually output from a call to \code{\link[rplos]{highplos}}}

\item{output}{Path and file name for output file. If NULL, a temp file is used.}

\item{browse}{Browse file in your default browse immediately after file creation.
If FALSE, the file is written, but not opened.}
}
\description{
Browse highlighted fragments in your default browser.
}
\examples{
\dontrun{
out <- highplos(q='alcohol', hl.fl = 'abstract', rows=10)
highbrow(out)

out <- highplos(q='alcohol', hl.fl = 'abstract', rows=100)
highbrow(out)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plosabstract.R
\name{plosabstract}
\alias{plosabstract}
\title{Search PLoS Journals abstracts.}
\usage{
plosabstract(
  q = NULL,
  fl = "id",
  fq = NULL,
  sort = NULL,
  start = 0,
  limit = 10,
  sleep = 6,
  errors = "simple",
  proxy = NULL,
  callopts = NULL,
  progress = NULL,
  ...
)
}
\arguments{
\item{q}{Search terms (character). You can search on specific fields by
doing 'field:your query'. For example, a real query on a specific field would
be 'author:Smith'.}

\item{fl}{Fields to return from search (character) [e.g., 'id,title'],
any combination of search fields (see the dataset \code{plosfields})}

\item{fq}{List specific fields to filter the query on (if NA, all queried).
The options for this parameter are the same as those for the fl parameter.
Note that using this parameter doesn't influence the actual query, but is used
to filter the results to a subset of those you want returned. For example,
if you want full articles only, you can do \code{'doc_type:full'}. In another example,
if you want only results from the journal PLOS One, you can do
\code{'journal_key:PLoSONE'}. See \code{\link{journalnamekey}} for journal
abbreviations.}

\item{sort}{Sort results according to a particular field, and specify ascending (asc)
or descending (desc) after a space; see examples. For example, to sort the
counter_total_all field in descending fashion, do sort='counter_total_all desc'}

\item{start}{Record to start at (used in combination with limit when
you need to cycle through more results than the max allowed=1000). See Pagination
below}

\item{limit}{Number of results to return (integer). Setting \code{limit=0} returns only
metadata. See Pagination below}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call to searchplos. However, if you are using searchplos in a loop or
lapply type call, do sleep parameter is used to prevent your IP address from being
blocked. You can only do 10 requests per minute, so one request every 6 seconds is
about right.}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{(list) optional curl options passed to \code{\link[crul]{HttpClient}}}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}.
only supports \code{httr::progress()}}

\item{...}{Additional Solr arguments}
}
\value{
Abstract content, in addition to any other fields requested in a list.
}
\description{
Search PLoS Journals abstracts.
}
\details{
Details:
}
\section{Faceting}{

Read more about faceting here: url{http://wiki.apache.org/solr/SimpleFacetParameters}
}

\section{Website vs. API behavior}{

Don't be surprised if queries you perform in a scripting language, like using \code{rplos}
in R, give different results than when searching for articles on the PLOS website. I am
not sure what exact defaults they use on their website. There are a few things to consider.
You can tweak which types of articles are returned: Try using the \code{article_type}
filter in the \code{fq} parameter. For which journal to search, e.g., do
\code{'journal_key:PLoSONE'}. See \code{journalnamekey()} for journal
abbreviations.
}

\section{Phrase searching}{

To search phrases, e.g., \strong{synthetic biology} as a single item, rather than
separate occurrences of \strong{synthetic} and \strong{biology}, simply put double
quotes around the phrase. For example, to search for cases of \strong{synthetic biology},
do \code{searchplos(q = '"synthetic biology"')}.

You can modify phrase searches as well. For example,
\code{searchplos(q = '"synthetic biology" ~ 10')} asks for cases of
\strong{synthetic biology} within 10 words of each other. See examples.
}

\section{Pagination}{

The \code{searchplos} function and the many functions that are wrappers around 
\code{searchplos} all do paginatino internally for you. That is, if you request for 
example, 2000 results, the max you can get in any one request is 1000, so we'll do 
two requests for you. And so on for larger requests. 

You can always do your own paginatino by doing a lapply type call or a for loop 
to cycle through pages of results.
}

\examples{
\dontrun{
plosabstract(q = 'drosophila', fl='abstract', limit=10)
plosabstract(q = 'drosophila', fl=c('id','author'), limit = 5)
plosabstract(q = 'drosophila', fl='author', limit = 5)
plosabstract(q = 'drosophila', fl=c('id','author','title'), limit = 5)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rplos-package.R
\name{rplos-defunct}
\alias{rplos-defunct}
\title{Defunct functions in rplos}
\description{
\itemize{
 \item \code{\link{crossref}}: service no longer provided -
 see the package \code{rcrossref}
 \item \code{\link{citations}}: service no longer available
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{concat_todf}
\alias{concat_todf}
\title{Concatenate author names, if present, used in other functions.}
\usage{
concat_todf(x)
}
\arguments{
\item{x}{a single list element with PLoS API returned nested elements}
}
\value{
data.frame of results, with authors concatenated to single vector.
}
\description{
Concatenate author names, if present, used in other functions.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plosfields.R
\docType{data}
\name{isocodes}
\alias{isocodes}
\title{Country names and FIPS codes}
\description{
Country names and FIPS codes
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citations.R
\name{citations}
\alias{citations}
\title{This function is defunct.}
\usage{
citations(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plosfields.R
\docType{data}
\name{plosfields}
\alias{plosfields}
\title{PLoS API fields to use for searching/retreiving data.}
\description{
PLoS API fields to use for searching/retreiving data.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highplos.r
\name{highplos}
\alias{highplos}
\title{Do highlighted searches on PLOS Journals full-text content}
\usage{
highplos(
  q,
  fl = NULL,
  fq = NULL,
  hl.fl = NULL,
  hl.snippets = NULL,
  hl.fragsize = NULL,
  hl.q = NULL,
  hl.mergeContiguous = NULL,
  hl.requireFieldMatch = NULL,
  hl.maxAnalyzedChars = NULL,
  hl.alternateField = NULL,
  hl.maxAlternateFieldLength = NULL,
  hl.preserveMulti = NULL,
  hl.maxMultiValuedToExamine = NULL,
  hl.maxMultiValuedToMatch = NULL,
  hl.formatter = NULL,
  hl.simple.pre = NULL,
  hl.simple.post = NULL,
  hl.fragmenter = NULL,
  hl.fragListBuilder = NULL,
  hl.fragmentsBuilder = NULL,
  hl.boundaryScanner = NULL,
  hl.bs.maxScan = NULL,
  hl.bs.chars = NULL,
  hl.bs.type = NULL,
  hl.bs.language = NULL,
  hl.bs.country = NULL,
  hl.useFastVectorHighlighter = NULL,
  hl.usePhraseHighlighter = NULL,
  hl.highlightMultiTerm = NULL,
  hl.regex.slop = NULL,
  hl.regex.pattern = NULL,
  hl.regex.maxAnalyzedChars = NULL,
  start = 0,
  rows = NULL,
  errors = "simple",
  proxy = NULL,
  callopts = list(),
  sleep = 6,
  ...
)
}
\arguments{
\item{q}{Search terms (character). You can search on specific fields by
doing 'field:your query'. For example, a real query on a specific field would
be 'author:Smith'.}

\item{fl}{Fields to return from search (character) [e.g., 'id,title'],
any combination of search fields [type 'data(plosfields)', then
'plosfields'].}

\item{fq}{List specific fields to filter the query on (if NA, all queried).
The options for this parameter are the same as those for the fl parameter.
Note that using this parameter doesn't influence the actual query, but is used
to filter the resuls to a subset of those you want returned. For example,
if you want full articles only, you can do 'doc_type:full'. In another example,
if you want only results from the journal PLOS One, you can do
'journal_key:PLoSONE'. See journalnamekey() for journal
abbreviations.}

\item{hl.fl}{A comma-separated list of fields for which to generate highlighted snippets.
If left blank, the fields highlighted for the LuceneQParser are the defaultSearchField
(or the df param if used) and for the DisMax parser the qf fields are used. A '*' can
be used to match field globs, e.g. 'text_*' or even '*' to highlight on all fields where
highlighting is possible. When using '*', consider adding hl.requireFieldMatch=TRUE.}

\item{hl.snippets}{Max no. of highlighted snippets to generate per field. Note:
it is possible for any number of snippets from zero to this value to be generated.
This parameter accepts per-field overrides. Default: 1.}

\item{hl.fragsize}{The size, in characters, of the snippets (aka fragments) created by
the highlighter. In the original Highlighter, "0" indicates that the whole field value
should be used with no fragmenting. See
https://cwiki.apache.org/confluence/display/solr/HighlightingParameters for more info}

\item{hl.q}{Set a query request to be highlighted. It overrides q parameter for
highlighting. Solr query syntax is acceptable for this parameter.}

\item{hl.mergeContiguous}{Collapse contiguous fragments into a single fragment. "true"
indicates contiguous fragments will be collapsed into single fragment. This parameter
accepts per-field overrides. This parameter makes sense for the original Highlighter
only. Default: FALSE.}

\item{hl.requireFieldMatch}{If TRUE, then a field will only be highlighted if the
query matched in this particular field (normally, terms are highlighted in all
requested fields regardless of which field matched the query). This only takes effect
if "hl.usePhraseHighlighter" is TRUE. Default: FALSE.}

\item{hl.maxAnalyzedChars}{How many characters into a document to look for suitable
snippets. This parameter makes sense for the original Highlighter only. Default: 51200.
You can assign a large value to this parameter and use hl.fragsize=0 to return
highlighting in large fields that have size greater than 51200 characters.}

\item{hl.alternateField}{If a snippet cannot be generated (due to no terms matching),
you can specify a field to use as the fallback. This parameter accepts per-field overrides.}

\item{hl.maxAlternateFieldLength}{If hl.alternateField is specified, this parameter
specifies the maximum number of characters of the field to return. Any value less than or
equal to 0 means unlimited. Default: unlimited.}

\item{hl.preserveMulti}{Preserve order of values in a multiValued list. Default: FALSE.}

\item{hl.maxMultiValuedToExamine}{When highlighting a multiValued field, stop examining
the individual entries after looking at this many of them. Will potentially return 0
snippets if this limit is reached before any snippets are found. If maxMultiValuedToMatch
is also specified, whichever limit is hit first will terminate looking for more.
Default: Integer.MAX_VALUE}

\item{hl.maxMultiValuedToMatch}{When highlighting a multiValued field, stop examining
the individual entries after looking at this many matches are found. If
maxMultiValuedToExamine is also specified, whichever limit is hit first will terminate
looking for more. Default: Integer.MAX_VALUE}

\item{hl.formatter}{Specify a formatter for the highlight output. Currently the only
legal value is "simple", which surrounds a highlighted term with a customizable pre- and
post text snippet. This parameter accepts per-field overrides. This parameter makes
sense for the original Highlighter only.}

\item{hl.simple.pre}{The text which appears before and after a highlighted term when using
the simple formatter. This parameter accepts per-field overrides. The default values are
"<em>" and "</em>" This parameter makes sense for the original Highlighter only. Use
hl.tag.pre and hl.tag.post for FastVectorHighlighter (see example under hl.fragmentsBuilder)}

\item{hl.simple.post}{The text which appears before and after a highlighted term when using
the simple formatter. This parameter accepts per-field overrides. The default values are
"<em>" and "</em>" This parameter makes sense for the original Highlighter only. Use
hl.tag.pre and hl.tag.post for FastVectorHighlighter (see example under hl.fragmentsBuilder)}

\item{hl.fragmenter}{Specify a text snippet generator for highlighted text. The standard
fragmenter is gap (which is so called because it creates fixed-sized fragments with gaps
for multi-valued fields). Another option is regex, which tries to create fragments that
"look like" a certain regular expression. This parameter accepts per-field overrides.
Default: "gap"}

\item{hl.fragListBuilder}{Specify the name of SolrFragListBuilder.  This parameter
makes sense for FastVectorHighlighter only. To create a fragSize=0 with the
FastVectorHighlighter, use the SingleFragListBuilder. This field supports per-field
overrides.}

\item{hl.fragmentsBuilder}{Specify the name of SolrFragmentsBuilder. This parameter makes
sense for FastVectorHighlighter only.}

\item{hl.boundaryScanner}{Configures how the boundaries of fragments are determined. By
default, boundaries will split at the character level, creating a fragment such as "uick
brown fox jumps over the la". Valid entries are breakIterator or simple, with breakIterator
being the most commonly used. This parameter makes sense for FastVectorHighlighter only.}

\item{hl.bs.maxScan}{Specify the length of characters to be scanned by SimpleBoundaryScanner.
Default: 10.  This parameter makes sense for FastVectorHighlighter only.}

\item{hl.bs.chars}{Specify the boundary characters, used by SimpleBoundaryScanner.
This parameter makes sense for FastVectorHighlighter only.}

\item{hl.bs.type}{Specify one of CHARACTER, WORD, SENTENCE and LINE, used by
BreakIteratorBoundaryScanner. Default: WORD. This parameter makes sense for
FastVectorHighlighter only.}

\item{hl.bs.language}{Specify the language for Locale that is used by
BreakIteratorBoundaryScanner. This parameter makes sense for FastVectorHighlighter only.
Valid entries take the form of ISO 639-1 strings.}

\item{hl.bs.country}{Specify the country for Locale that is used by
BreakIteratorBoundaryScanner. This parameter makes sense for FastVectorHighlighter only.
Valid entries take the form of ISO 3166-1 alpha-2 strings.}

\item{hl.useFastVectorHighlighter}{Use FastVectorHighlighter. FastVectorHighlighter
requires the field is termVectors=on, termPositions=on and termOffsets=on. This
parameter accepts per-field overrides. Default: FALSE}

\item{hl.usePhraseHighlighter}{Use SpanScorer to highlight phrase terms only when
they appear within the query phrase in the document. Default: TRUE.}

\item{hl.highlightMultiTerm}{If the SpanScorer is also being used, enables highlighting
for range/wildcard/fuzzy/prefix queries. Default: FALSE. This parameter makes sense
for the original Highlighter only.}

\item{hl.regex.slop}{Factor by which the regex fragmenter can stray from the ideal
fragment size (given by hl.fragsize) to accomodate the regular expression. For
instance, a slop of 0.2 with fragsize of 100 should yield fragments between 80
and 120 characters in length. It is usually good to provide a slightly smaller
fragsize when using the regex fragmenter. Default: .6. This parameter makes sense
for the original Highlighter only.}

\item{hl.regex.pattern}{The regular expression for fragmenting. This could be
used to extract sentences (see example solrconfig.xml) This parameter makes sense
for the original Highlighter only.}

\item{hl.regex.maxAnalyzedChars}{Only analyze this many characters from a field
when using the regex fragmenter (after which, the fragmenter produces fixed-sized
fragments). Applying a complicated regex to a huge field is expensive.
Default: 10000. This parameter makes sense for the original Highlighter only.}

\item{start}{Record to start at (used in combination with limit when
you need to cycle through more results than the max allowed=1000)}

\item{rows}{Number of results to return (integer)}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{Optional additional curl options passed to
\code{\link[crul]{HttpClient}}}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call to searchplos. However, if you are using searchplos in a loop or
lapply type call, do sleep parameter is used to prevent your IP address from being
blocked. You can only do 10 requests per minute, so one request every 6 seconds is
about right.}

\item{...}{Further arguments passed on to solr_highlight}
}
\value{
A list.
}
\description{
Do highlighted searches on PLOS Journals full-text content
}
\examples{
\dontrun{
highplos(q='alcohol', hl.fl = 'abstract', rows=10)
highplos(q='alcohol', hl.fl = c('abstract','title'), rows=10)
highplos(q='everything:"sports alcohol"~7', hl.fl='everything')
highplos(q='alcohol', hl.fl='abstract', hl.fragsize=20, rows=5)
highplos(q='alcohol', hl.fl='abstract', hl.snippets=5, rows=5)
highplos(q='alcohol', hl.fl='abstract', hl.snippets=5,
  hl.mergeContiguous='true', rows=5)
highplos(q='alcohol', hl.fl='abstract', hl.useFastVectorHighlighter='true',
  rows=5)
highplos(q='everything:"experiment"', fq='doc_type:full', rows=100,
  hl.fl = 'title')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plosfigtabcaps.R
\name{plosfigtabcaps}
\alias{plosfigtabcaps}
\title{Search PLoS Journals figure and table captions.}
\usage{
plosfigtabcaps(
  q = NULL,
  fl = "id",
  fq = NULL,
  sort = NULL,
  start = 0,
  limit = 10,
  sleep = 6,
  errors = "simple",
  proxy = NULL,
  callopts = NULL,
  progress = NULL,
  ...
)
}
\arguments{
\item{q}{Search terms (character). You can search on specific fields by
doing 'field:your query'. For example, a real query on a specific field would
be 'author:Smith'.}

\item{fl}{Fields to return from search (character) [e.g., 'id,title'],
any combination of search fields (see the dataset \code{plosfields})}

\item{fq}{List specific fields to filter the query on (if NA, all queried).
The options for this parameter are the same as those for the fl parameter.
Note that using this parameter doesn't influence the actual query, but is used
to filter the results to a subset of those you want returned. For example,
if you want full articles only, you can do \code{'doc_type:full'}. In another example,
if you want only results from the journal PLOS One, you can do
\code{'journal_key:PLoSONE'}. See \code{\link{journalnamekey}} for journal
abbreviations.}

\item{sort}{Sort results according to a particular field, and specify ascending (asc)
or descending (desc) after a space; see examples. For example, to sort the
counter_total_all field in descending fashion, do sort='counter_total_all desc'}

\item{start}{Record to start at (used in combination with limit when
you need to cycle through more results than the max allowed=1000). See Pagination
below}

\item{limit}{Number of results to return (integer). Setting \code{limit=0} returns only
metadata. See Pagination below}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call to searchplos. However, if you are using searchplos in a loop or
lapply type call, do sleep parameter is used to prevent your IP address from being
blocked. You can only do 10 requests per minute, so one request every 6 seconds is
about right.}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{(list) optional curl options passed to \code{\link[crul]{HttpClient}}}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}.
only supports \code{httr::progress()}}

\item{...}{Additional Solr arguments}
}
\value{
fields that you specify to return in a data.frame, along with the
		DOI's found.
}
\description{
Search PLoS Journals figure and table captions.
}
\details{
Details:
}
\section{Faceting}{

Read more about faceting here: url{http://wiki.apache.org/solr/SimpleFacetParameters}
}

\section{Website vs. API behavior}{

Don't be surprised if queries you perform in a scripting language, like using \code{rplos}
in R, give different results than when searching for articles on the PLOS website. I am
not sure what exact defaults they use on their website. There are a few things to consider.
You can tweak which types of articles are returned: Try using the \code{article_type}
filter in the \code{fq} parameter. For which journal to search, e.g., do
\code{'journal_key:PLoSONE'}. See \code{journalnamekey()} for journal
abbreviations.
}

\section{Phrase searching}{

To search phrases, e.g., \strong{synthetic biology} as a single item, rather than
separate occurrences of \strong{synthetic} and \strong{biology}, simply put double
quotes around the phrase. For example, to search for cases of \strong{synthetic biology},
do \code{searchplos(q = '"synthetic biology"')}.

You can modify phrase searches as well. For example,
\code{searchplos(q = '"synthetic biology" ~ 10')} asks for cases of
\strong{synthetic biology} within 10 words of each other. See examples.
}

\section{Pagination}{

The \code{searchplos} function and the many functions that are wrappers around 
\code{searchplos} all do paginatino internally for you. That is, if you request for 
example, 2000 results, the max you can get in any one request is 1000, so we'll do 
two requests for you. And so on for larger requests. 

You can always do your own paginatino by doing a lapply type call or a for loop 
to cycle through pages of results.
}

\examples{
\dontrun{
plosfigtabcaps('ecology', 'id', limit=100)
plosfigtabcaps(q='ecology', fl='figure_table_caption', limit=10)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plosword.R
\name{plosword}
\alias{plosword}
\title{Search results on a keyword over all fields in PLoS Journals.}
\usage{
plosword(terms, vis = FALSE, ...)
}
\arguments{
\item{terms}{search terms (character)}

\item{vis}{visualize results in bar plot or not (TRUE or FALSE)}

\item{...}{Optional additional curl options passed to
\code{\link[crul]{HttpClient}}}
}
\value{
Number of search results (vis = FALSE), or number of search in a
table and a histogram of results (vis = TRUE).
}
\description{
Search results on a keyword over all fields in PLoS Journals.
}
\examples{
\dontrun{
plosword('Helianthus')
plosword(list('monkey','replication','design','sunflower','whale'),
   vis = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatarticleurl.R
\name{formatarticleurl}
\alias{formatarticleurl}
\title{Format a URL for a specific article in a specific PLoS journal.}
\usage{
formatarticleurl(doi, journal)
}
\arguments{
\item{doi}{digital object identifier for an article in PLoS Journals}

\item{journal}{quoted journal name (character)}
}
\value{
Get url for the article to use in your browser, etc.
}
\description{
Format a URL for a specific article in a specific PLoS journal.
}
\details{
Choose the journal abbreviation, and the appropriate base URL
   will be inclued in the built URL.  Options are:
   PLoSBiology, PLoSGenetics, PLoSComputationalBiology, PLoSMedicine,
  PLoSONE, PLoSNeglectedTropicalDiseases, or PLoSPathogens.
}
\examples{
\dontrun{
formatarticleurl(doi="10.1371/journal.pone.0004045", journal='PLoSONE')
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/facetplos.r
\name{facetplos}
\alias{facetplos}
\title{Do faceted searches on PLOS Journals full-text content}
\usage{
facetplos(
  q = "*:*",
  facet.query = NA,
  facet.field = NA,
  facet.prefix = NA,
  facet.sort = NA,
  facet.limit = NA,
  facet.offset = NA,
  facet.mincount = NA,
  facet.missing = NA,
  facet.method = NA,
  facet.enum.cache.minDf = NA,
  facet.threads = NA,
  facet.date = NA,
  facet.date.start = NA,
  facet.date.end = NA,
  facet.date.gap = NA,
  facet.date.hardend = NA,
  facet.date.other = NA,
  facet.date.include = NA,
  facet.range = NA,
  facet.range.start = NA,
  facet.range.end = NA,
  facet.range.gap = NA,
  facet.range.hardend = NA,
  facet.range.other = NA,
  facet.range.include = NA,
  start = NA,
  rows = NA,
  url = NA,
  sleep = 6,
  errors = "simple",
  proxy = NULL,
  callopts = list(),
  ...
)
}
\arguments{
\item{q}{Query terms.}

\item{facet.query}{This param allows you to specify an arbitrary query in the
Lucene default syntax to generate a facet count. By default, faceting returns
a count of the unique terms for a "field", while facet.query allows you to
determine counts for arbitrary terms or expressions. This parameter can be
specified multiple times to indicate that multiple queries should be used as
separate facet constraints. It can be particularly useful for numeric range
based facets, or prefix based facets -- see example below (i.e. price:[* TO 500]
and  price:[501 TO *]).}

\item{facet.field}{This param allows you to specify a field which should be
treated as a facet. It will iterate over each Term in the field and generate a
facet count using that Term as the constraint. This parameter can be specified
multiple times to indicate multiple facet fields. None of the other params in
this section will have any effect without specifying at least one field name
using this param.}

\item{facet.prefix}{Limits the terms on which to facet to those starting with
the given string prefix. Note that unlike fq, this does not change the search
results -- it merely reduces the facet values returned to those beginning with
the specified prefix. This parameter can be specified on a per field basis.}

\item{facet.sort}{See \code{\link[solrium]{solr_facet}}.}

\item{facet.limit}{This param indicates the maximum number of constraint counts
that should be returned for the facet fields. A negative value means unlimited.
Default: 100. Can be specified on a per field basis.}

\item{facet.offset}{This param indicates an offset into the list of constraints
to allow paging. Default: 0. This parameter can be specified on a per field basis.}

\item{facet.mincount}{This param indicates the minimum counts for facet fields
should be included in the response. Default: 0. This parameter can be specified
on a per field basis.}

\item{facet.missing}{Set to "true" this param indicates that in addition to the
Term based constraints of a facet field, a count of all matching results which
have no value for the field should be computed. Default: FALSE. This parameter
can be specified on a per field basis.}

\item{facet.method}{See \code{\link[solrium]{solr_facet}}.}

\item{facet.enum.cache.minDf}{This param indicates the minimum document frequency
(number of documents matching a term) for which the filterCache should be used
when determining the constraint count for that term. This is only used when
facet.method=enum method of faceting. A value greater than zero will decrease
memory usage of the filterCache, but increase the query time. When faceting on
a field with a very large number of terms, and you wish to decrease memory usage,
try a low value of 25 to 50 first. Default: 0, causing the filterCache to be used
for all terms in the field. This parameter can be specified on a per field basis.}

\item{facet.threads}{This param will cause loading the underlying fields used in
faceting to be executed in parallel with the number of threads specified. Specify
as facet.threads=# where # is the maximum number of threads used. Omitting this
parameter or specifying the thread count as 0 will not spawn any threads just as
before. Specifying a negative number of threads will spin up to Integer.MAX_VALUE
threads. Currently this is limited to the fields, range and query facets are not
yet supported. In at least one case this has reduced warmup times from 20 seconds
to under 5 seconds.}

\item{facet.date}{Specify names of fields (of type DateField) which should be
treated as date facets. Can be specified multiple times to indicate multiple
date facet fields.}

\item{facet.date.start}{The lower bound for the first date range for all Date
Faceting on this field. This should be a single date expression which may use
the DateMathParser syntax. Can be specified on a per field basis.}

\item{facet.date.end}{The minimum upper bound for the last date range for all
Date Faceting on this field (see facet.date.hardend for an explanation of what
the actual end value may be greater). This should be a single date expression
which may use the DateMathParser syntax. Can be specified on a per field basis.}

\item{facet.date.gap}{The size of each date range expressed as an interval to
be added to the lower bound using the DateMathParser syntax. Eg:
facet.date.gap=+1DAY. Can be specified on a per field basis.}

\item{facet.date.hardend}{A Boolean parameter instructing Solr what to do in the
event that facet.date.gap does not divide evenly between facet.date.start and
facet.date.end. If this is true, the last date range constraint will have an
upper bound of facet.date.end; if false, the last date range will have the smallest
possible upper bound greater then facet.date.end such that the range is exactly
facet.date.gap wide. Default: FALSE. This parameter can be specified on a per
field basis.}

\item{facet.date.other}{See \code{\link[solrium]{solr_facet}}.}

\item{facet.date.include}{See \code{\link[solrium]{solr_facet}}.}

\item{facet.range}{Indicates what field to create range facets for. Example:
facet.range=price&facet.range=age}

\item{facet.range.start}{The lower bound of the ranges. Can be specified on a
per field basis. Example: f.price.facet.range.start=0.0&f.age.facet.range.start=10}

\item{facet.range.end}{The upper bound of the ranges. Can be specified on a per
field basis. Example: f.price.facet.range.end=1000.0&f.age.facet.range.start=99}

\item{facet.range.gap}{The size of each range expressed as a value to be added
to the lower bound. For date fields, this should be expressed using the
DateMathParser syntax. (ie: facet.range.gap=+1DAY). Can be specified
on a per field basis. Example: f.price.facet.range.gap=100&f.age.facet.range.gap=10}

\item{facet.range.hardend}{A Boolean parameter instructing Solr what to do in the
event that facet.range.gap does not divide evenly between facet.range.start and
facet.range.end. If this is true, the last range constraint will have an upper
bound of facet.range.end; if false, the last range will have the smallest possible
upper bound greater then facet.range.end such that the range is exactly
facet.range.gap wide. Default: FALSE. This parameter can be specified on a
per field basis.}

\item{facet.range.other}{See \code{\link[solrium]{solr_facet}}.}

\item{facet.range.include}{See \code{\link[solrium]{solr_facet}}.}

\item{start}{Record to start at, default to beginning.}

\item{rows}{Number of records to return.}

\item{url}{URL endpoint}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call. However, if you are doing many calls in a loop or lapply type call,
sleep parameter is used to prevent your IP address from being blocked. You can only
do 10 requests per minute, so one request every 6 seconds is about right.}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{Further args passed on to \code{\link[crul]{HttpClient}}}

\item{...}{Further args to \code{\link[solrium]{solr_facet}}}
}
\value{
A list
}
\description{
Do faceted searches on PLOS Journals full-text content
}
\examples{
\dontrun{
# Facet on a single field
facetplos(q='*:*', facet.field='journal')
facetplos(q='alcohol', facet.field='article_type')

# Facet on multiple fields
facetplos(q='alcohol', facet.field=c('journal','subject'))

# Using mincount
facetplos(q='alcohol', facet.field='journal', facet.mincount='500')

# Using facet.query to get counts
## A single facet.query term
facetplos(q='*:*', facet.field='journal', facet.query='cell')
## Many facet.query terms
facetplos(q='*:*', facet.field='journal', facet.query='cell,bird')

# Range faceting
facetplos(q='*:*', url=url, facet.range='counter_total_all',
   facet.range.start=5, facet.range.end=1000, facet.range.gap=10)
facetplos(q='alcohol', facet.range='alm_facebookCount', facet.range.start=1000,
   facet.range.end=5000, facet.range.gap = 100)

# Range faceting with > 1 field, same settings
facetplos(q='*:*', url=url, facet.range=c('counter_total_all','alm_twitterCount'),
 facet.range.start=5, facet.range.end=1000, facet.range.gap=10)

# Range faceting with > 1 field, different settings
facetplos(q='*:*', url=url, facet.range=c('counter_total_all','alm_twitterCount'),
 f.counter_total_all.facet.range.start=5, f.counter_total_all.facet.range.end=1000,
 f.counter_total_all.facet.range.gap=10, f.alm_twitterCount.facet.range.start=5,
 f.alm_twitterCount.facet.range.end=1000, f.alm_twitterCount.facet.range.gap=10)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fulltext.R
\name{full_text_urls}
\alias{full_text_urls}
\title{Create urls for full text articles in PLOS journals.}
\usage{
full_text_urls(doi)
}
\arguments{
\item{doi}{One or more doi's}
}
\value{
One or more urls, same length as input vector of dois
}
\description{
Create urls for full text articles in PLOS journals.
}
\details{
We give \strong{NA} for DOIs that are for annotations. Those can 
easily be removed like \code{Filter(Negate(is.na), res)}
}
\examples{
\dontrun{
full_text_urls(doi='10.1371/journal.pone.0086169')
full_text_urls(doi='10.1371/journal.pbio.1001845')
full_text_urls(doi=c('10.1371/journal.pone.0086169',
  '10.1371/journal.pbio.1001845'))

# contains some annotation DOIs
dois <- searchplos(q = "*:*", fq='doc_type:full', limit=20)$data$id
full_text_urls(dois)

# contains no annotation DOIs
dois <- searchplos(q = "*:*", 
  fq=list('doc_type:full', 'article_type:"Research Article"'), 
limit=20)$data$id
full_text_urls(dois)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossref.R
\name{crossref}
\alias{crossref}
\title{Lookup article info via CrossRef with DOI.}
\usage{
crossref(...)
}
\description{
Lookup article info via CrossRef with DOI.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchplos.R
\name{searchplos}
\alias{searchplos}
\title{Base function to search PLoS Journals}
\usage{
searchplos(
  q = NULL,
  fl = "id",
  fq = NULL,
  sort = NULL,
  start = 0,
  limit = 10,
  sleep = 6,
  errors = "simple",
  proxy = NULL,
  callopts = list(),
  progress = NULL,
  ...
)
}
\arguments{
\item{q}{Search terms (character). You can search on specific fields by
doing 'field:your query'. For example, a real query on a specific field would
be 'author:Smith'.}

\item{fl}{Fields to return from search (character) [e.g., 'id,title'],
any combination of search fields (see the dataset \code{plosfields})}

\item{fq}{List specific fields to filter the query on (if NA, all queried).
The options for this parameter are the same as those for the fl parameter.
Note that using this parameter doesn't influence the actual query, but is used
to filter the results to a subset of those you want returned. For example,
if you want full articles only, you can do \code{'doc_type:full'}. In another example,
if you want only results from the journal PLOS One, you can do
\code{'journal_key:PLoSONE'}. See \code{\link{journalnamekey}} for journal
abbreviations.}

\item{sort}{Sort results according to a particular field, and specify ascending (asc)
or descending (desc) after a space; see examples. For example, to sort the
counter_total_all field in descending fashion, do sort='counter_total_all desc'}

\item{start}{Record to start at (used in combination with limit when
you need to cycle through more results than the max allowed=1000). See Pagination
below}

\item{limit}{Number of results to return (integer). Setting \code{limit=0} returns only
metadata. See Pagination below}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call to searchplos. However, if you are using searchplos in a loop or
lapply type call, do sleep parameter is used to prevent your IP address from being
blocked. You can only do 10 requests per minute, so one request every 6 seconds is
about right.}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{(list) optional curl options passed to \code{\link[crul]{HttpClient}}}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}.
only supports \code{httr::progress()}}

\item{...}{Additional Solr arguments}
}
\value{
An object of class "plos", with a list of length two, each
element being a list itself.
}
\description{
Base function to search PLoS Journals
}
\details{
Details:
}
\section{Faceting}{

Read more about faceting here: url{http://wiki.apache.org/solr/SimpleFacetParameters}
}

\section{Website vs. API behavior}{

Don't be surprised if queries you perform in a scripting language, like using \code{rplos}
in R, give different results than when searching for articles on the PLOS website. I am
not sure what exact defaults they use on their website. There are a few things to consider.
You can tweak which types of articles are returned: Try using the \code{article_type}
filter in the \code{fq} parameter. For which journal to search, e.g., do
\code{'journal_key:PLoSONE'}. See \code{journalnamekey()} for journal
abbreviations.
}

\section{Phrase searching}{

To search phrases, e.g., \strong{synthetic biology} as a single item, rather than
separate occurrences of \strong{synthetic} and \strong{biology}, simply put double
quotes around the phrase. For example, to search for cases of \strong{synthetic biology},
do \code{searchplos(q = '"synthetic biology"')}.

You can modify phrase searches as well. For example,
\code{searchplos(q = '"synthetic biology" ~ 10')} asks for cases of
\strong{synthetic biology} within 10 words of each other. See examples.
}

\section{Pagination}{

The \code{searchplos} function and the many functions that are wrappers around 
\code{searchplos} all do paginatino internally for you. That is, if you request for 
example, 2000 results, the max you can get in any one request is 1000, so we'll do 
two requests for you. And so on for larger requests. 

You can always do your own paginatino by doing a lapply type call or a for loop 
to cycle through pages of results.
}

\examples{
\dontrun{
searchplos(q='ecology', fl=c('id','publication_date'), limit = 2)
searchplos('ecology', fl=c('id','publication_date'), limit = 2)
searchplos('ecology', c('id','title'), limit = 2)

# Get only full article DOIs
out <- searchplos(q="*:*", fl='id', fq='doc_type:full', start=0, limit=250)
head(out$data)

# Get DOIs for only PLoS One articles
out <- searchplos(q="*:*", fl='id', fq='journal_key:PLoSONE', start=0, limit=15)
out$data

# Get DOIs for full article in PLoS One
out <- searchplos(q="*:*", fl='id', fq=list('journal_key:PLoSONE',
   'doc_type:full'), limit=50)
out$data

# Serch for many q
q <- c('ecology','evolution','science')
lapply(q, function(x) searchplos(x, limit=2))

# Query to get some PLOS article-level metrics, notice difference between two outputs
out <- searchplos(q="*:*", fl=c('id','counter_total_all','alm_twitterCount'),fq='doc_type:full')
out_sorted <- searchplos(q="*:*", fl=c('id','counter_total_all','alm_twitterCount'),
   fq='doc_type:full', sort='counter_total_all desc')
out$data
out_sorted$data

# Show me all articles that have these two words less then about 15 words apart.
searchplos(q='everything:"sports alcohol"~15', fl='title', fq='doc_type:full')

# Now let's try to narrow our results to 7 words apart. Here I'm changing the ~15 to ~7
searchplos(q='everything:"sports alcohol"~7', fl='title', fq='doc_type:full')

# A list of articles about social networks that are popular on a social network
searchplos(q="*:*",fl=c('id','alm_twitterCount'),
   fq=list('doc_type:full','subject:"Social networks"','alm_twitterCount:[100 TO 10000]'),
   sort='counter_total_month desc')

# Now, lets also only look at articles that have seen some activity on twitter.
# Add "fq=alm_twitterCount:[1 TO *]" as a parameter within the fq argument.
searchplos(q='everything:"sports alcohol"~7', fl=c('alm_twitterCount','title'),
   fq=list('doc_type:full','alm_twitterCount:[1 TO *]'))
searchplos(q='everything:"sports alcohol"~7', fl=c('alm_twitterCount','title'),
   fq=list('doc_type:full','alm_twitterCount:[1 TO *]'),
   sort='counter_total_month desc')

# Return partial doc parts
## Return Abstracts only
out <- searchplos(q='*:*', fl=c('doc_partial_body','doc_partial_parent_id'),
   fq=list('doc_type:partial', 'doc_partial_type:Abstract'), limit=3)
## Return Title's only
out <- searchplos(q='*:*', fl=c('doc_partial_body','doc_partial_parent_id'),
   fq=list('doc_type:partial', 'doc_partial_type:Title'), limit=3)

# Remove DOIs for annotations (i.e., corrections)
searchplos(q='*:*', fl=c('id','article_type'),
   fq='-article_type:correction', limit=100)

# Remove DOIs for annotations (i.e., corrections) and Viewpoints articles
searchplos(q='*:*', fl=c('id','article_type'),
   fq=list('-article_type:correction','-article_type:viewpoints'), limit=100)

# Get eissn codes
searchplos(q='*:*', fl=c('id','journal','eissn','cross_published_journal_eissn'),
   fq="doc_type:full", limit = 60)

searchplos(q='*:*', fl=c('id','journal','eissn','cross_published_journal_eissn'),
   limit = 2000)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rplos-package.R
\docType{package}
\name{rplos}
\alias{rplos}
\alias{rplos-package}
\title{Connect with PLoS API data}
\description{
\code{rplos} provides an R interface to the PLoS Search API. More
information about each function can be found in its help documentation.
}
\section{rplos functions}{


\pkg{rplos} functions make HTTP requests using the \pkg{crul} package,
and parse json using the \pkg{jsonlite} package.
}

\section{PLoS API key}{


You used to need an API key to use this package - no longer needed
}

\section{Tutorials}{


See the rOpenSci website for a tutorial:
https://ropensci.org/tutorials/rplos_tutorial.html
}

\section{Throttling}{

Beware, PLOS recently has started throttling requests. That is,
they will give error messages like "(503) Service Unavailable -
The server cannot process the request due to a high load", which
probably means you've done too many requests in a certain time period.

Here's what they say (http://api.plos.org/solr/faq/#solr_api_recommended_usage)
on the matter:

"Please limit your API requests to 7200 requests a day, 300 per hour, 10 per
minute and allow 5 seconds for your search to return results. If you exceed this
threshold, we will lock out your IP address. If you're a high-volume user of
the PLOS Search API and need more API requests a day, please contact us at
api@plos.org to discuss your options. We currently limit API users to no more
than five concurrent connections from a single IP address.""
}

\examples{
\dontrun{
searchplos(q='ecology', fl=c('id','publication_date'), limit = 2)

# Get only full article DOIs
out <- searchplos(q="*:*", fl='id', fq='doc_type:full', start=0, limit=250)
head(out$data)

# Get DOIs for only PLoS One articles
out <- searchplos(q="*:*", fl='id', fq='journal_key:PLoSONE',
  start=0, limit=15)
head(out$data)
}

}
\author{
Scott Chamberlain

Carl Boettiger \email{cboettig@gmail.com}

Karthik Ram \email{karthik.ram@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plosauthor.R
\name{plosauthor}
\alias{plosauthor}
\title{Search PLoS Journals authors.}
\usage{
plosauthor(
  q = NULL,
  fl = "id",
  fq = NULL,
  sort = NULL,
  start = 0,
  limit = 10,
  sleep = 6,
  errors = "simple",
  proxy = NULL,
  callopts = NULL,
  progress = NULL,
  ...
)
}
\arguments{
\item{q}{Search terms (character). You can search on specific fields by
doing 'field:your query'. For example, a real query on a specific field would
be 'author:Smith'.}

\item{fl}{Fields to return from search (character) [e.g., 'id,title'],
any combination of search fields (see the dataset \code{plosfields})}

\item{fq}{List specific fields to filter the query on (if NA, all queried).
The options for this parameter are the same as those for the fl parameter.
Note that using this parameter doesn't influence the actual query, but is used
to filter the results to a subset of those you want returned. For example,
if you want full articles only, you can do \code{'doc_type:full'}. In another example,
if you want only results from the journal PLOS One, you can do
\code{'journal_key:PLoSONE'}. See \code{\link{journalnamekey}} for journal
abbreviations.}

\item{sort}{Sort results according to a particular field, and specify ascending (asc)
or descending (desc) after a space; see examples. For example, to sort the
counter_total_all field in descending fashion, do sort='counter_total_all desc'}

\item{start}{Record to start at (used in combination with limit when
you need to cycle through more results than the max allowed=1000). See Pagination
below}

\item{limit}{Number of results to return (integer). Setting \code{limit=0} returns only
metadata. See Pagination below}

\item{sleep}{Number of seconds to wait between requests. No need to use this for
a single call to searchplos. However, if you are using searchplos in a loop or
lapply type call, do sleep parameter is used to prevent your IP address from being
blocked. You can only do 10 requests per minute, so one request every 6 seconds is
about right.}

\item{errors}{(character) One of simple or complete. Simple gives http code and
error message on an error, while complete gives both http code and error message,
and stack trace, if available.}

\item{proxy}{List of arguments for a proxy connection, including one or more of:
url, port, username, password, and auth. See \code{\link[crul]{proxy}} for
help, which is used to construct the proxy connection.}

\item{callopts}{(list) optional curl options passed to \code{\link[crul]{HttpClient}}}

\item{progress}{a function with logic for printing a progress
bar for an HTTP request, ultimately passed down to \pkg{curl}.
only supports \code{httr::progress()}}

\item{...}{Additional Solr arguments}
}
\value{
Author names, in addition to any other fields requested in a
   data.frame.
}
\description{
Search PLoS Journals authors.
}
\details{
Details:
}
\section{Faceting}{

Read more about faceting here: url{http://wiki.apache.org/solr/SimpleFacetParameters}
}

\section{Website vs. API behavior}{

Don't be surprised if queries you perform in a scripting language, like using \code{rplos}
in R, give different results than when searching for articles on the PLOS website. I am
not sure what exact defaults they use on their website. There are a few things to consider.
You can tweak which types of articles are returned: Try using the \code{article_type}
filter in the \code{fq} parameter. For which journal to search, e.g., do
\code{'journal_key:PLoSONE'}. See \code{journalnamekey()} for journal
abbreviations.
}

\section{Phrase searching}{

To search phrases, e.g., \strong{synthetic biology} as a single item, rather than
separate occurrences of \strong{synthetic} and \strong{biology}, simply put double
quotes around the phrase. For example, to search for cases of \strong{synthetic biology},
do \code{searchplos(q = '"synthetic biology"')}.

You can modify phrase searches as well. For example,
\code{searchplos(q = '"synthetic biology" ~ 10')} asks for cases of
\strong{synthetic biology} within 10 words of each other. See examples.
}

\section{Pagination}{

The \code{searchplos} function and the many functions that are wrappers around 
\code{searchplos} all do paginatino internally for you. That is, if you request for 
example, 2000 results, the max you can get in any one request is 1000, so we'll do 
two requests for you. And so on for larger requests. 

You can always do your own paginatino by doing a lapply type call or a for loop 
to cycle through pages of results.
}

\examples{
\dontrun{
plosauthor('Smith', 'id', limit=50)
plosauthor(q='Smith', fl=c('id','author'), limit=10)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{addmissing}
\alias{addmissing}
\title{Adds elements in a list that are missing because they were not returned
in the PLoS API call.}
\usage{
addmissing(x)
}
\arguments{
\item{x}{A list with PLoS API returned nested elements}
}
\value{
A list with the missing element added with an
		"na", if it is missing.
}
\description{
Adds elements in a list that are missing because they were not returned
in the PLoS API call.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plosviews.R
\name{plosviews}
\alias{plosviews}
\title{Search PLoS Journals by article views.}
\usage{
plosviews(search, byfield = NULL, views = "alltime", limit = NULL, ...)
}
\arguments{
\item{search}{search terms (character)}

\item{byfield}{field to search by, e.g., subject, author, etc. (character)}

\item{views}{views all time (alltime) or views last 30 days (last30)
(character)}

\item{limit}{number of results to return (integer)}

\item{...}{Optional additional curl options passed to 
\code{\link[crul]{HttpClient}}}
}
\description{
Search PLoS Journals by article views.
}
\examples{
\dontrun{
plosviews('10.1371/journal.pone.0002154', 'id', 'alltime')
plosviews('10.1371/journal.pone.0002154', 'id', 'last30')
plosviews('10.1371/journal.pone.0002154', 'id', 'alltime,last30')
plosviews(search='marine ecology', byfield='subject', limit=50)
plosviews(search='evolution', views = 'alltime', limit = 99)
plosviews('bird', views = 'alltime', limit = 99)
}
}
