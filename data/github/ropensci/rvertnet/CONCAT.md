rvertnet
=======



[![cran checks](https://cranchecks.info/badges/worst/rvertnet)](https://cranchecks.info/pkgs/rvertnet)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/rvertnet/workflows/R-check/badge.svg)](https://github.com/ropensci/rvertnet/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/rvertnet/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rvertnet?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rvertnet)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rvertnet)](https://cran.r-project.org/package=rvertnet)


`rvertnet` is a client for interacting with [VertNet.org](http://vertnet.org/).

VertNet.org API docs: <https://github.com/VertNet/webapp/wiki/The-API-search-function>

## Installation

Stable CRAN version


```r
install.packages("rvertnet")
```

Or development version from GitHub


```r
remotes::install_github("ropensci/rvertnet")
```


```r
library('rvertnet')
```

## Search by term

Search for _Aves_ in the state of _California_, limit to 10 records


```r
res <- searchbyterm(class = "Aves", stateprovince = "California",
    limit = 10, messages = FALSE)
```

Inspect metadata


```r
res$meta
#> $request_date
#> [1] "2021-05-12T23:08:19.676358"
#> 
#> $response_records
#> [1] 10
#> 
#> $submitted_query
#> [1] "class:Aves stateprovince:California"
#> 
#> $request_origin
#> [1] "45.505106,-122.675026"
#> 
#> $limit
#> [1] 10
#> 
#> $last_cursor
#> [1] "False:Cq8FCooDCtwC9wAAABn_____jIGJmo2LkZqL0o-QjYuek96WkZuah9LNz87M0s_H0s_H_wAA_3RtoKCZi4ygoP8AAP9dno-PmpGYlpGa_wAA_3N0bZaRm5qH_wAA_12biJz_AAD_c3Rtm5CcoJab_wAA_12cnozQkI2R0IqNkdKcnouek5CY0pyejNKQjZHSzs_Pzs7_AAD_c3-cnozQkI2R0IqNkdKcnouek5CY0pyejNKQjZHSzs_Pzs7_AAD__wD-__6MgYmajYuRmovSj5CNi56T3paRm5qH0s3PzszSz8fSz8f_AHRtoKCZi4ygoP8AXZ6Pj5qRmJaRmv8Ac3RtlpGbmof_AF2biJz_AHN0bZuQnKCWm_8AXZyejNCQjZHQio2R0pyei56TkJjSnJ6M0pCNkdLOz8_Ozv8Ac3-cnozQkI2R0IqNkdKcnouek5CY0pyejNKQjZHSzs_Pzs7_AP_-EAohBN0EkB08Gxk5AAAAAOb___9IClAAWgsJqvR3VAbTv8MQA2Cl4NG4AhINRG9jdW1lbnRJbmRleBruAShBTkQgKElTICJjdXN0b21lcl9uYW1lIiAiYXBwZW5naW5lIikgKElTICJncm91cF9uYW1lIiAic352ZXJ0bmV0LXBvcnRhbCIpIChJUyAibmFtZXNwYWNlIiAiaW5kZXgtMjAxMy0wOC0wOCIpIChJUyAiaW5kZXhfbmFtZSIgImR3YyIpIChBTkQgKE9SIChRVCAiQXZlcyIgInJ0ZXh0X2NsYXNzIikgKElTICJyYXRvbV9jbGFzcyIgImF2ZXMiKSkgKFFUICJDYWxpZm9ybmlhIiAicnRleHRfc3RhdGVwcm92aW5jZSIpKSk6GQoMKE4gb3JkZXJfaWQpEAEZAAAAAAAA8P9KBQgAQOgH"
#> 
#> $query_version
#> [1] "search.py 2016-08-15T16:43+02:00"
#> 
#> $matching_records
#> [1] ">10000"
#> 
#> $api_version
#> [1] "api.py 2017-11-24T12:16-03:00"
```

Inspect data. A `dplyr` data.frame is given back, so you get a nice brief data summary:


```r
res$data[,1:5]
#> # A tibble: 10 x 5
#>    higherclassification       stateprovince basisofrecord month decimallongitude
#>    <chr>                      <chr>         <chr>         <chr> <chr>           
#>  1 Animalia | Chordata |  | … California    PreservedSpe… 2     -121.7833       
#>  2 Animalia | Chordata |  | … California    PreservedSpe… 6     -122.15         
#>  3 Animalia | Chordata |  | … California    PreservedSpe… 5     -120.9014       
#>  4 Animalia; Chordata; Aves;… California    PreservedSpe… 1     -121.93300      
#>  5 Animalia; Chordata; Aves;… California    PreservedSpe… 1     -121.93300      
#>  6 Animalia; Chordata; Aves;… California    PreservedSpe… 7     -121.85760      
#>  7 Animalia; Chordata; Aves;… California    PreservedSpe… 7     -121.85760      
#>  8 Animalia; Chordata; Aves;… California    PreservedSpe… 7     -121.85760      
#>  9 Animalia; Chordata; Aves;… California    PreservedSpe… 7     -121.85760      
#> 10 Animalia; Chordata; Aves;… California    PreservedSpe… 6     -121.85760
```

Search for _Mustela nigripes_ in the states of _Wyoming_ or _South Dakota_, limit to 20 records


```r
res <- searchbyterm(specificepithet = "nigripes",
    stateprovince = "(wyoming OR south dakota)", 
    limit = 20, messages = FALSE)
res$data[,1:5]
#> # A tibble: 20 x 5
#>    month decimallongitude startdayofyear accessrights                    kingdom
#>    <chr> <chr>            <chr>          <chr>                           <chr>  
#>  1 12    -100.8276541162  336            http://vertnet.org/resources/n… Animal…
#>  2 03    -100.9827        64             http://vertnet.org/resources/n… Animal…
#>  3 1     -100.759483      1              http://vertnet.org/resources/n… Animal…
#>  4 3     -100.7373        67             http://biodiversity.ku.edu/res… Animal…
#>  5 11    <NA>             305            http://vertnet.org/resources/n… Animal…
#>  6 10    <NA>             282            <NA>                            Animal…
#>  7 8     <NA>             234            <NA>                            Animal…
#>  8 12    <NA>             342            <NA>                            Animal…
#>  9 12    <NA>             358            http://www.vertnet.org/resourc… Animal…
#> 10 1     <NA>             1              http://vertnet.org/resources/n… Animal…
#> 11 1     <NA>             1              http://vertnet.org/resources/n… Animal…
#> 12 11    <NA>             313            <NA>                            Animal…
#> 13 9     <NA>             272            <NA>                            Animal…
#> 14 12    <NA>             335            <NA>                            Animal…
#> 15 9     <NA>             259            <NA>                            Animal…
#> 16 10    <NA>             297            <NA>                            Animal…
#> 17 12    <NA>             339            <NA>                            Animal…
#> 18 11    <NA>             305            <NA>                            Animal…
#> 19 11    <NA>             315            <NA>                            Animal…
#> 20 <NA>  <NA>             <NA>           http://vertnet.org/resources/n… Animal…
```

### dplyr downstream

You can pass the data object directly on to `dplyr` functions. Here, we get a table of record counts by species in descending order.


```r
library("dplyr")
out <- searchbyterm(genus = "Ochotona", limit = 800)
out$data %>%
  group_by(scientificname) %>%
  summarise(count = length(scientificname)) %>%
  arrange(desc(count))
#> # A tibble: 20 x 2
#>    scientificname                  count
#>    <chr>                           <int>
#>  1 Ochotona princeps                 450
#>  2 Ochotona pallasi                  129
#>  3 Ochotona princeps saxatilis       103
#>  4 Ochotona hyperborea                30
#>  5 Ochotona dauurica                  21
#>  6 Ochotona collaris                  15
#>  7 Ochotona princeps figginsi         14
#>  8 Ochotona princeps taylori           8
#>  9 Ochotona princeps schisticeps       6
#> 10 Ochotona alpina                     4
#> 11 Ochotona princeps muiri             4
#> 12 Ochotona hyperborea mantchurica     3
#> 13 Ochotona princeps incana            3
#> 14 Ochotona princeps princeps          3
#> 15 Ochotona princeps murri             2
#> 16 Ochotona princeps brunnescens       1
#> 17 Ochotona princeps jewetti           1
#> 18 Ochotona princeps tutelata          1
#> 19 Ochotona princeps uinta             1
#> 20 Ochotona princeps ventorum          1
```


## Big data

Specifies a termwise search (like `searchbyterm()`), but requests that all available records be made available for download as a tab-delimited text file.


```r
bigsearch(genus = "ochotona", rfile = "pikaRecords", email = "big@@search.luv")
#> Processing request...
#>
#> Download of records file 'mydata' requested for 'you@gmail.com'
#>
#> Query/URL: "http://api.vertnet-portal.appspot.com/api/download?q=%7B%22q%22:%22genus:ochotona%22,%22n%22:%22mydata%22,%22e%22:%22you@gmail.com%22%7D"
#>
#> Thank you! Download instructions will be sent by email.
```

## Spatial search


```r
res <- spatialsearch(lat = 33.529, long = -105.694, radius = 2000,
    limit = 10, messages = FALSE)
res$data[,1:5]
#> # A tibble: 10 x 5
#>    month decimallongitude startdayofyear minimumelevationin… accessrights       
#>    <chr> <chr>            <chr>          <chr>               <chr>              
#>  1 07    -105.68633       193            2182.368            http://vertnet.org…
#>  2 07    -105.705479      196            2023.872            http://vertnet.org…
#>  3 07    -105.705479      196            2023.872            http://vertnet.org…
#>  4 07    -105.705479      196            2023.872            http://vertnet.org…
#>  5 07    -105.705479      196            2023.872            http://vertnet.org…
#>  6 07    -105.705479      196            2023.872            http://vertnet.org…
#>  7 07    -105.705479      196            2023.872            http://vertnet.org…
#>  8 07    -105.705479      196            2023.872            http://vertnet.org…
#>  9 07    -105.705479      196            2023.872            http://vertnet.org…
#> 10 07    -105.705479      196            2023.872            http://vertnet.org…
```

## Contributors

* Scott Chamberlain [@sckott](https://github.com/sckott)
* Chris Ray [@Pika8tona](https://github.com/Pika8tona)
* Vijay Barve [@vijaybarve](https://github.com/vijaybarve)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rvertnet/issues).
* License: MIT
* Get citation information for `rvertnet` in R doing `citation(package = 'rvertnet')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rvernet 0.8.2
===============

### MINOR IMPROVEMENTS

* vignette fix

rvertnet 0.8.0
===============

### NEW FEATURES

* `searchbyterm()` and `bigsearch()` reworked: both functions now have the first parameter as `...`, which accepts any valid query parameter. There were so many query parameters for these functions it was a bit overwhelming. See `?searchbyterm` docs for details  (#66)

### MINOR IMPROVEMENTS

* decode the request URL before printing to the R console so users can more easily see what request they have done (#67)
* vignette title fix (#68)

### BUG FIXES 

* `searchbyterm()` fix: booleans need to be converted to VertNet's expected `0/1` instead of `true/false` (#66)


rvertnet 0.7.0
===============

### BUG FIXES 

* add month and day params to `searchbyterm` (#64)


rvertnet 0.6.2
===============

### BUG FIXES 

* A small data source used in one function was on the web, and
was moved - that data source now within the pkg as quite small, and 
now pkg won't break when the file is moved again (#61) (#62)


rvertnet 0.6.0
===============

Added Code of Conduct.

### NEW FEATURES

* Now using `crul` package for HTTP requests instead of `httr` (#57)
* Note that `verbose` parameter has been replaced with `messages` throughout
the package.
* Now with function for search for trait data: `traitsearch()` (#55)

### DEFUNCT AND DEPRECATED

* All `dump` functions are now defunct. Those functions tried to help users
work with bulk Vertnet data - the setup has gotten too complex (#56)

### MINOR IMPROVEMENTS

* Improvements to documentation for `traitsearch()` function on what 
fields have given data (#58) thanks @gaurav
* `vertsearch()` and `searchbyterm()` gain new parameter `only_dwc`, which 
allows to optionally only return Darwin Core fields

### BUG FIXES 

* Small fix to `vertsummary()` (#59)


rvertnet 0.5.0
===============

### NEW FEATURES

* `searchbyterm()` gains new parameter `query` to allow full text search, 
much like `vertsearch()`, but with the ability to also use all the parameters
available in `searchbyterm()` (#53)

### MINOR IMPROVEMENTS

* Use `dplyr::bind_rows` instead of the deprecated `dplyr::rbind_all` (#51)
* remove personal email address from tests (#52)
* Namespace base R pkg fxn calls (`methods`/`stats`/`utils`), and removed 
some package dependencies that we didn't really need (`plyr`) (#54)

rvertnet 0.4.4
===============

### MINOR IMPROVEMENTS

* Updated docs to better indicate how to use the cursor feature (#49)
* Now using explicit encoding specification when using `httr::content()` (#47)

### BUG FIXES

* Fixed `externalptr` error in the internal `vert_GET()` function (#48)

rvertnet 0.4.1
===============

### BUG FIXES

* Fixed a bug in `bigsearch()` in which we had forgotten to do
internal conversion of logical input to 0/1 needed by the web
API (#46)

rvertnet 0.4.0
===============

### NEW FEATURES

* New set of functions to make working with VertNet data dumps
easier. `dump_links()` gives you links to various data dump
resources; `dump_init()` initialized a SQLite database connection;
`dump_tbl()` creates a `dplyr::tbl` object, which can then be used
in a `dplyr` query. This setup requires that the user manually
download data dumps uncompress, and load into SQLite. We hope to
make this process easier in the future. (#36)

### MINOR IMPROVEMENTS

* Fixes to `vertmap()` for new `ggplot2` version (#43)
* Added note to docs for `bigsearch()` for how to read in data
after obtaining the data (#44)

### BUG FIXES

* Fix to the `searchbyterm()` function. When the parameter `stateprovince`
was used, lead to error, as that param requires different handling than
other params. (#45)

rvertnet 0.3.4
===============

### NEW FEATURES

* New function `vert_id()` to get occurrence records by occurenceid,
that is, single occurrence ids. (#40)

### MINOR IMPROVEMENTS

* Explicitly import non-base R functions (#39)

### BUG FIXES

* Lowercase `occurenceID` to `occurrenceid` to simplify life (#41)

rvertnet 0.3.0
===============

### NEW FEATURES

* `searchbyterm()` and `bigsearch()` have some parameters that accept multiple values.
Fixed to allow this (#37)
* Internals of `searchbyterm()`, `spatialsearch()`, and `vertsearch()` reworked to
use cursor so we internally do paging for you for bigger result sets. (#25)

### MINOR IMPROVEMENTS

* Replaced `data.table` import with `dplyr`
* Using `skip_on_cran()` (#38)
* Minor vignette updates (#35)
* Metadata now returned in data requests (#33)

rvertnet 0.2.2
===============

Package completely reworked for the new VertNet API.

### NEW FEATURES

* The functions `vertavailablemaps()`, `vertlocations()`,
`vertoccurrence()`, `vertoccurrencecount()`, `vertproviders()`,
`verttaxa()` are now defunct. You can call these functions, but
they print an error message, saying they are defunct.
* Gained new functions `bigsearch()`, `searchbyterm()`,
`spatialsearch()`, and `vertsummary()`.
* Gained new author: Chris Ray

### MINOR IMPROVEMENTS

* `RJSONIO` replaced with `jsonlite`
* Changed from CC0 to MIT license

rvertnet 0.0-5
------------

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.0.5 Patched
* ubuntu 14.04 (on GitHub Actions), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 downstream dependency, with no problems found. Summary at https://github.com/ropensci/rvertnet/tree/master/revdep

--------

This version fixes a vignette issue, and remove LazyData in DESCRIPTION.

Sincerely,
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

* Submit an issue on the [Issues page](https://github.com/ropensci/rvertnet/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rvertnet.git`
* Make sure to track progress upstream (i.e., on our version of `rvertnet` at `ropensci/rvertnet`) by doing `git remote add upstream https://github.com/ropensci/rvertnet.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rvertnet`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rvertnet introduction}
%\VignetteEncoding{UTF-8}
-->



rvertnet introduction
=====================

`rvertnet` is a client for interacting with [VertNet.org](http://vertnet.org/).

## Installation

You can install the stable version from CRAN:


```r
install.packages("rvertnet")
```

Or the development version from GitHub using the `devtools` package:


```r
install.packages("devtools")
devtools::install_github("ropensci/rvertnet")
```


```r
library('rvertnet')
```

## Search by term

Search for _Aves_ in the state of _California_, limit to 10 records


```r
res <- searchbyterm(class = "Aves", state = "California", limit = 10, messages = FALSE)
```

All major functions (`searchbyterm()`, `spatialsearch()`, `vertsearch()`) give back a `meta` (for metadata, in a list) and `data` (for data, in a data.frame) slot. The metadata:


```r
res$meta
#> $request_date
#> [1] "2020-01-29T02:22:59.517328"
#> 
#> $response_records
#> [1] 10
#> 
#> $submitted_query
#> [1] "class:Aves"
#> 
#> $request_origin
#> [1] "45.505106,-122.675026"
#> 
#> $limit
#> [1] 10
#> 
#> $last_cursor
#> [1] "False:CskFCtIDCqQD9wAAABn_____jIGJmo2LkZqL0o-QjYuek96WkZuah9LNz87M0s_H0s_H_wAA_3RtoKCZi4ygoP8AAP9dno-PmpGYlpGa_wAA_3N0bZaRm5qH_wAA_12biJz_AAD_c3Rtm5CcoJab_wAA_12ektCQjZGWi5eQk5CYhtDPz87GmcjGztLGzcbJ0svPzprSxsfGy9Kdycqcz53NnMfGnpv_AAD_c3-ektCQjZGWi5eQk5CYhtDPz87GmcjGztLGzcbJ0svPzprSxsfGy9Kdycqcz53NnMfGnpv_AAD__wD-__6MgYmajYuRmovSj5CNi56T3paRm5qH0s3PzszSz8fSz8f_AHRtoKCZi4ygoP8AXZ6Pj5qRmJaRmv8Ac3RtlpGbmof_AF2biJz_AHN0bZuQnKCWm_8AXZ6S0JCNkZaLl5CTkJiG0M_PzsaZyMbO0sbNxsnSy8_OmtLGx8bL0p3JypzPnc2cx8aem_8Ac3-ektCQjZGWi5eQk5CYhtDPz87GmcjGztLGzcbJ0svPzprSxsfGy9Kdycqcz53NnMfGnpv_AP_-EAohBN0EkB08Gxk5AAAAAOb___9IClAAWgsJJ5RN7qiPsmIQA2Cb0daDBxINRG9jdW1lbnRJbmRleBrAAShBTkQgKElTICJjdXN0b21lcl9uYW1lIiAiYXBwZW5naW5lIikgKElTICJncm91cF9uYW1lIiAic352ZXJ0bmV0LXBvcnRhbCIpIChJUyAibmFtZXNwYWNlIiAiaW5kZXgtMjAxMy0wOC0wOCIpIChJUyAiaW5kZXhfbmFtZSIgImR3YyIpIChPUiAoUVQgIkF2ZXMiICJydGV4dF9jbGFzcyIpIChJUyAicmF0b21fY2xhc3MiICJhdmVzIikpKToZCgwoTiBvcmRlcl9pZCkQARkAAAAAAADw_0oFCABA6Ac"
#> 
#> $query_version
#> [1] "search.py 2016-08-15T16:43+02:00"
#> 
#> $matching_records
#> [1] ">10000"
#> 
#> $api_version
#> [1] "api.py 2017-11-24T12:16-03:00"
```

The data


```r
res$data
#> # A tibble: 10 x 46
#>    kingdom recordedby higherclassific… stateprovince basisofrecord month decimallongitude phylum references year 
#>    <chr>   <chr>      <chr>            <chr>         <chr>         <chr> <chr>            <chr>  <chr>      <chr>
#>  1 Animal… NSW STATE… Animalia | Chor… New South Wa… PreservedSpe… 11    153.316          Chord… http://po… 1866 
#>  2 Animal… TARONGA Z… Animalia | Chor… New South Wa… PreservedSpe… 6     151.33535        Chord… http://po… 2011 
#>  3 Animal… MCLENNAN,… Animalia | Chor… Queensland    PreservedSpe… 4     145.35           Chord… http://po… 1915 
#>  4 Animal… WATTS, MI… Animalia | Chor… New South Wa… PreservedSpe… 1     146.016          Chord… http://po… 1983 
#>  5 Animal… HOLCOMBE,… Animalia | Chor… New South Wa… PreservedSpe… 7     151.766          Chord… http://po… 1973 
#>  6 Animal… RHODES, C… Animalia | Chor… New South Wa… PreservedSpe… 10    150.866          Chord… http://po… 1931 
#>  7 Animal… C. HEDLEY… Animalia | Chor… Queensland    PreservedSpe… 10    142.6            Chord… http://po… 1907 
#>  8 Animal… ROBINSON,… Animalia | Chor… New South Wa… PreservedSpe… 10    149.266          Chord… http://po… 1897 
#>  9 Animal… SHARP, GE… Animalia | Chor… Queensland    PreservedSpe… 11    145.85           Chord… http://po… 1908 
#> 10 Animal… J.A. KEAS… Animalia | Chor… New South Wa… PreservedSpe… 3     146.266          Chord… http://po… 1958 
#> # … with 36 more variables: startdayofyear <chr>, taxonrank <chr>, specificepithet <chr>,
#> #   bibliographiccitation <chr>, family <chr>, countrycode <chr>, geodeticdatum <chr>,
#> #   coordinateuncertaintyinmeters <chr>, highergeography <chr>, accessrights <chr>, verbatimlocality <chr>,
#> #   verbatimeventdate <chr>, day <chr>, eventid <chr>, collectioncode <chr>, occurrencestatus <chr>,
#> #   locationremarks <chr>, coordinateprecision <chr>, institutioncode <chr>, scientificname <chr>, class <chr>,
#> #   decimallatitude <chr>, occurrenceid <chr>, language <chr>, license <chr>, country <chr>,
#> #   georeferenceverificationstatus <chr>, modified <chr>, eventdate <chr>, nomenclaturalcode <chr>, continent <chr>,
#> #   genus <chr>, order <chr>, catalognumber <chr>, enddayofyear <chr>, vernacularname <chr>
```

Search for _Mustela nigripes_ in the states of _Wyoming_ or _South Dakota_, limit to 20 records


```r
res <- searchbyterm(specificepithet = "nigripes", genus = "Mustela", state = "(wyoming OR south dakota)", limit = 20, messages = FALSE)
res$data
#> # A tibble: 20 x 76
#>    month decimallongitude startdayofyear accessrights kingdom verbatimcoordin… day   identificationv… occurrenceid
#>    <chr> <chr>            <chr>          <chr>        <chr>   <chr>            <chr> <chr>            <chr>       
#>  1 1     -88.305352       1              http://vert… Animal… decimal degrees  1     legacy           http://arct…
#>  2 03    -104.77472       74             http://vert… Animal… <NA>             15    legacy           http://arct…
#>  3 02    -103.731861      52             http://vert… Animal… <NA>             21    legacy           http://arct…
#>  4 12    -105.0137067407  349            http://vert… Animal… decimal degrees  15    field            http://arct…
#>  5 2     -103.067931      32             http://vert… Animal… <NA>             1     legacy           http://arct…
#>  6 1     -103.067931      1              http://vert… Animal… <NA>             1     legacy           http://arct…
#>  7 02    -103.067931      40             http://vert… Animal… <NA>             09    field            http://arct…
#>  8 05    -104.926320116   126            http://vert… Animal… <NA>             05    legacy           http://arct…
#>  9 02    -104.79742       42             http://vert… Animal… <NA>             11    field            http://arct…
#> 10 04    -106.1329632593  108            http://vert… Animal… decimal degrees  18    field            http://arct…
#> 11 10    -105.064706      304            http://vert… Animal… decimal degrees  31    legacy           http://arct…
#> 12 4     -106.3467709375  92             http://vert… Animal… decimal degrees  1     field            http://arct…
#> 13 05    -104.225829      133            http://vert… Animal… <NA>             13    field            http://arct…
#> 14 09    -105.873904      258            http://vert… Animal… <NA>             15    field            http://arct…
#> 15 12    -105.298898      362            http://vert… Animal… <NA>             28    legacy           http://arct…
#> 16 06    -105.376986      152            http://vert… Animal… <NA>             01    student          http://arct…
#> 17 11    -104.3831505257  305            http://vert… Animal… decimal degrees  01    student          http://arct…
#> 18 11    -104.7714765     314            http://vert… Animal… UTM              10    legacy           http://arct…
#> 19 09    -106.9094        267            http://vert… Animal… decimal degrees  23    legacy           http://arct…
#> 20 08    -107.5579841     234            http://vert… Animal… UTM              21    legacy           http://arct…
#> # … with 67 more variables: identificationqualifier <chr>, georeferenceddate <chr>, verbatimeventdate <chr>,
#> #   coordinateuncertaintyinmeters <chr>, higherclassification <chr>, sex <chr>, year <chr>, specificepithet <chr>,
#> #   basisofrecord <chr>, geodeticdatum <chr>, occurrenceremarks <chr>, highergeography <chr>, continent <chr>,
#> #   scientificname <chr>, language <chr>, institutionid <chr>, country <chr>, genus <chr>,
#> #   georeferenceprotocol <chr>, family <chr>, stateprovince <chr>, county <chr>, phylum <chr>, references <chr>,
#> #   georeferencedby <chr>, taxonrank <chr>, verbatimlocality <chr>, institutioncode <chr>, eventremarks <chr>,
#> #   organismid <chr>, eventtime <chr>, preparations <chr>, license <chr>, dynamicproperties <chr>,
#> #   georeferenceverificationstatus <chr>, modified <chr>, eventdate <chr>, individualcount <chr>,
#> #   bibliographiccitation <chr>, verbatimcoordinates <chr>, georeferencesources <chr>, nomenclaturalcode <chr>,
#> #   catalognumber <chr>, locality <chr>, informationwithheld <chr>, collectioncode <chr>, collectionid <chr>,
#> #   class <chr>, previousidentifications <chr>, identificationremarks <chr>, decimallatitude <chr>,
#> #   locationaccordingto <chr>, othercatalognumbers <chr>, identifiedby <chr>, associatedmedia <chr>, order <chr>,
#> #   enddayofyear <chr>, typestatus <chr>, recordedby <chr>, dateidentified <chr>, locationremarks <chr>,
#> #   associatedsequences <chr>, recordnumber <chr>, minimumelevationinmeters <chr>, maximumelevationinmeters <chr>,
#> #   lifestage <chr>, establishmentmeans <chr>
```

Search for class _Aves_, in the state of _Nevada_, with a coordinate uncertainty range (in meters) of less than 25 meters


```r
res <- searchbyterm(class = "Aves", stateprovince = "Nevada", error = "<25", messages = FALSE)
res$data
#> # A tibble: 1,000 x 91
#>    georeferencepro… higherclassific… stateprovince lifestage month decimallongitude phylum verbatimlongitu… year 
#>    <chr>            <chr>            <chr>         <chr>     <chr> <chr>            <chr>  <chr>            <chr>
#>  1 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad.     3     -117.73567       Chord… -117.7356796°    1886 
#>  2 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad      12    -119.19015       Chord… -119.1901537°    1912 
#>  3 MaNIS/HerpNet/O… Animalia; Chord… Nevada        Nestling  5     -119.55677       Chord… -119.5567779°    1918 
#>  4 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad.     9     -119.94328       Chord… -119.9432804°    1924 
#>  5 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad.     9     -119.94328       Chord… -119.9432804°    1924 
#>  6 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad.     9     -119.94328       Chord… -119.9432804°    1924 
#>  7 MaNIS/HerpNet/O… Animalia; Chord… Nevada        Downy     6     -119.93907       Chord… -119.9390724°    1936 
#>  8 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad.     6     -119.93907       Chord… -119.9390724°    1936 
#>  9 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad.     6     -119.93907       Chord… -119.9390724°    1936 
#> 10 MaNIS/HerpNet/O… Animalia; Chord… Nevada        U-Ad.     6     -119.93907       Chord… -119.9390724°    1936 
#> # … with 990 more rows, and 82 more variables: specificepithet <chr>, bibliographiccitation <chr>,
#> #   verbatimlatitude <chr>, family <chr>, locality <chr>, geodeticdatum <chr>, coordinateuncertaintyinmeters <chr>,
#> #   highergeography <chr>, continent <chr>, scientificnameauthorship <chr>, day <chr>, kingdom <chr>,
#> #   institutioncode <chr>, scientificname <chr>, preparations <chr>, sex <chr>, class <chr>, county <chr>,
#> #   decimallatitude <chr>, occurrenceid <chr>, language <chr>, license <chr>, basisofrecord <chr>, country <chr>,
#> #   collectioncode <chr>, modified <chr>, eventdate <chr>, verbatimeventdate <chr>, references <chr>, genus <chr>,
#> #   order <chr>, catalognumber <chr>, georeferencesources <chr>, recordedby <chr>, occurrenceremarks <chr>,
#> #   infraspecificepithet <chr>, georeferenceremarks <chr>, startdayofyear <chr>, dynamicproperties <chr>,
#> #   enddayofyear <chr>, recordnumber <chr>, minimumelevationinmeters <chr>, othercatalognumbers <chr>,
#> #   georeferenceddate <chr>, accessrights <chr>, institutionid <chr>, georeferencedby <chr>, taxonrank <chr>,
#> #   verbatimlocality <chr>, countrycode <chr>, georeferenceverificationstatus <chr>, occurrencestatus <chr>,
#> #   vernacularname <chr>, nomenclaturalcode <chr>, datasetname <chr>, collectionid <chr>,
#> #   verbatimcoordinatesystem <chr>, identificationverificationstatus <chr>, identificationqualifier <chr>,
#> #   locationremarks <chr>, organismid <chr>, individualcount <chr>, previousidentifications <chr>,
#> #   locationaccordingto <chr>, verbatimcoordinates <chr>, identifiedby <chr>, fieldnumber <chr>, disposition <chr>,
#> #   ownerinstitutioncode <chr>, rightsholder <chr>, identificationremarks <chr>, maximumelevationinmeters <chr>,
#> #   verbatimelevation <chr>, typestatus <chr>, habitat <chr>, dateidentified <chr>, eventtime <chr>,
#> #   samplingprotocol <chr>, associatedmedia <chr>, island <chr>, eventremarks <chr>, associatedoccurrences <chr>
```

## Spatial search

Spatial search service allows only to search on a point defined by latitude and longitude pair, with a radius (meters) from that point. All three parameters are required. 


```r
res <- spatialsearch(lat = 33.529, lon = -105.694, radius = 2000, limit = 10, messages = FALSE)
res$data
#> # A tibble: 10 x 62
#>    month decimallongitude startdayofyear minimumelevatio… accessrights kingdom day   identificationv… occurrenceid
#>    <chr> <chr>            <chr>          <chr>            <chr>        <chr>   <chr> <chr>            <chr>       
#>  1 07    -105.68633       193            2182.368         http://vert… Animal… 12    legacy           http://arct…
#>  2 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#>  3 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#>  4 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#>  5 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#>  6 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#>  7 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#>  8 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#>  9 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#> 10 07    -105.705479      196            2023.872         http://vert… Animal… 14    legacy           http://arct…
#> # … with 53 more variables: identificationqualifier <chr>, georeferenceddate <chr>, verbatimeventdate <chr>,
#> #   coordinateuncertaintyinmeters <chr>, higherclassification <chr>, sex <chr>, year <chr>, specificepithet <chr>,
#> #   basisofrecord <chr>, geodeticdatum <chr>, occurrenceremarks <chr>, highergeography <chr>, continent <chr>,
#> #   scientificname <chr>, language <chr>, institutionid <chr>, country <chr>, genus <chr>,
#> #   georeferenceprotocol <chr>, family <chr>, stateprovince <chr>, county <chr>, phylum <chr>, references <chr>,
#> #   georeferencedby <chr>, taxonrank <chr>, verbatimlocality <chr>, institutioncode <chr>, organismid <chr>,
#> #   maximumelevationinmeters <chr>, preparations <chr>, recordedby <chr>, license <chr>, dynamicproperties <chr>,
#> #   georeferenceverificationstatus <chr>, modified <chr>, eventdate <chr>, individualcount <chr>,
#> #   bibliographiccitation <chr>, georeferencesources <chr>, catalognumber <chr>, locality <chr>, recordnumber <chr>,
#> #   collectioncode <chr>, class <chr>, previousidentifications <chr>, decimallatitude <chr>,
#> #   locationaccordingto <chr>, othercatalognumbers <chr>, identifiedby <chr>, nomenclaturalcode <chr>, order <chr>,
#> #   enddayofyear <chr>
```

## Global full text search

`vertsearch()` provides a simple full text search against all fields. For more info see [the docs](https://github.com/VertNet/webapp/wiki/The-API-search-function#global-full-text-search). An example:


```r
res <- vertsearch(taxon = "aves", state = "california", limit = 10)
res$data
#> # A tibble: 10 x 57
#>    higherclassific… stateprovince basisofrecord month decimallongitude phylum references year  startdayofyear
#>    <chr>            <chr>         <chr>         <chr> <chr>            <chr>  <chr>      <chr> <chr>         
#>  1 Animalia | Chor… California    PreservedSpe… 2     -121.7833        Chord… http://po… 1974  54            
#>  2 Animalia | Chor… California    PreservedSpe… 6     -122.15          Chord… http://po… 1973  155           
#>  3 Animalia | Chor… California    PreservedSpe… 5     -120.9014        Chord… http://po… 1919  142           
#>  4 Animalia; Chord… South Caroli… PreservedSpe… 2     -79.86151        Chord… http://po… 1904  <NA>          
#>  5 Animalia; Chord… California    PreservedSpe… 1     -121.93300       Chord… http://po… 1908  <NA>          
#>  6 Animalia; Chord… California    PreservedSpe… 1     -121.93300       Chord… http://po… 1908  <NA>          
#>  7 Animalia; Chord… California    PreservedSpe… 7     -121.85760       Chord… http://po… 1907  <NA>          
#>  8 Animalia; Chord… California    PreservedSpe… 7     -121.85760       Chord… http://po… 1907  <NA>          
#>  9 Animalia; Chord… California    PreservedSpe… 7     -121.85760       Chord… http://po… 1907  <NA>          
#> 10 Animalia; Chord… California    PreservedSpe… 7     -121.85760       Chord… http://po… 1907  <NA>          
#> # … with 48 more variables: taxonrank <chr>, specificepithet <chr>, bibliographiccitation <chr>, family <chr>,
#> #   countrycode <chr>, geodeticdatum <chr>, coordinateuncertaintyinmeters <chr>, highergeography <chr>,
#> #   continent <chr>, verbatimlocality <chr>, day <chr>, kingdom <chr>, collectioncode <chr>, occurrencestatus <chr>,
#> #   coordinateprecision <chr>, institutioncode <chr>, scientificname <chr>, locality <chr>, class <chr>,
#> #   vernacularname <chr>, county <chr>, decimallatitude <chr>, occurrenceid <chr>, language <chr>, license <chr>,
#> #   country <chr>, georeferenceverificationstatus <chr>, modified <chr>, eventdate <chr>, nomenclaturalcode <chr>,
#> #   verbatimeventdate <chr>, genus <chr>, order <chr>, catalognumber <chr>, enddayofyear <chr>,
#> #   locationremarks <chr>, infraspecificepithet <chr>, georeferenceprotocol <chr>, lifestage <chr>, recordedby <chr>,
#> #   verbatimlatitude <chr>, scientificnameauthorship <chr>, preparations <chr>, georeferenceremarks <chr>, sex <chr>,
#> #   verbatimlongitude <chr>, georeferencesources <chr>, dynamicproperties <chr>
```

Limit the number of records returned (under 1000)


```r
res <- vertsearch("(kansas state OR KSU)", limit = 200)
res$data
#> # A tibble: 200 x 78
#>    individualcount georeferencepro… recordedby bibliographicci… stateprovince basisofrecord month decimallongitude
#>    <chr>           <chr>            <chr>      <chr>            <chr>         <chr>         <chr> <chr>           
#>  1 8               GEOLocate (Rios… H. W. Rob… Academy of Natu… Oklahoma      PreservedSpe… 10    -94.707552      
#>  2 11              GEOLocate (Rios… H. W. Rob… Academy of Natu… Oklahoma      PreservedSpe… 10    -94.707552      
#>  3 3               GEOLocate (Rios… H. W. Rob… Academy of Natu… Oklahoma      PreservedSpe… 10    -94.707552      
#>  4 <NA>            <NA>             <NA>       California Acad… Kansas        PreservedSpe… 11    -95.4569444444  
#>  5 <NA>            <NA>             <NA>       California Acad… Kansas        PreservedSpe… 5     -96.7475194444  
#>  6 <NA>            <NA>             <NA>       California Acad… Kansas        PreservedSpe… 8     -101.0889700000 
#>  7 1               VertNet Georefe… MCCOY, C … Carnegie Museum… Sonora        PreservedSpe… 8     -112.57         
#>  8 1               VertNet Georefe… MCCOY, C … Carnegie Museum… Sonora        PreservedSpe… 8     -111.37         
#>  9 1               VertNet Georefe… MCCOY, C … Carnegie Museum… Sonora        PreservedSpe… 8     -111.37         
#> 10 1               VertNet Georefe… MCCOY, C … Carnegie Museum… Oklahoma      PreservedSpe… 6     -100.49         
#> # … with 190 more rows, and 70 more variables: phylum <chr>, references <chr>, georeferencedby <chr>, year <chr>,
#> #   taxonrank <chr>, specificepithet <chr>, family <chr>, countrycode <chr>, locality <chr>, geodeticdatum <chr>,
#> #   coordinateuncertaintyinmeters <chr>, highergeography <chr>, continent <chr>, day <chr>, kingdom <chr>,
#> #   georeferenceddate <chr>, footprintwkt <chr>, institutioncode <chr>, scientificname <chr>, preparations <chr>,
#> #   disposition <chr>, class <chr>, identificationremarks <chr>, county <chr>, decimallatitude <chr>,
#> #   occurrenceid <chr>, language <chr>, license <chr>, country <chr>, georeferenceverificationstatus <chr>,
#> #   othercatalognumbers <chr>, infraspecificepithet <chr>, eventdate <chr>, identifiedby <chr>,
#> #   nomenclaturalcode <chr>, fieldnumber <chr>, verbatimeventdate <chr>, genus <chr>, order <chr>,
#> #   catalognumber <chr>, collectioncode <chr>, higherclassification <chr>, lifestage <chr>, startdayofyear <chr>,
#> #   occurrenceremarks <chr>, verbatimlocality <chr>, georeferencesources <chr>, verbatimcoordinatesystem <chr>,
#> #   institutionid <chr>, modified <chr>, dateidentified <chr>, enddayofyear <chr>, georeferenceremarks <chr>,
#> #   accessrights <chr>, occurrencestatus <chr>, sex <chr>, establishmentmeans <chr>, recordnumber <chr>,
#> #   collectionid <chr>, dynamicproperties <chr>, reproductivecondition <chr>, minimumelevationinmeters <chr>,
#> #   identificationverificationstatus <chr>, identificationqualifier <chr>, organismid <chr>,
#> #   maximumelevationinmeters <chr>, previousidentifications <chr>, locationaccordingto <chr>,
#> #   verbatimcoordinates <chr>, datasetname <chr>
```

Pass output of `vertsearch()` to a map


```r
out <- vertsearch(tax = "(mustela nivalis OR mustela erminea)")
vertmap(out)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

## Lots of data

For `searchbyterm()`, `spatialsearch()`, and `vertsearch()`, you can request more than 1000 records. VertNet limits each request to 1000 records, but internally in this package, if you request more than 1000 records, we'll continue to send requests to get all the records you want. See the [VertNet docs](https://github.com/VertNet/webapp/wiki/The-API-search-function#retrieving-large-result-sets) for more information on this.

## Email dump of data

`bigsearch()` specifies a termwise search (like `searchbyterm()`), but requests that all available records be made available for download as a tab-delimited text file.


```r
bigsearch(genus = "ochotona", rfile = "mydata", email = "you@gmail.com")
#> Processing request...
#> 
#> Download of records file 'mydata' requested for 'you@gmail.com'
#> 
#> Query/URL: "http://api.vertnet-portal.appspot.com/api/download?q=%7B%22q%22:%22genus:ochotona%22,%22n%22:%22mydata%22,%22e%22:%22you@gmail.com%22%7D"
#> 
#> Thank you! Download instructions will be sent by email.
```

## Messages

In the previous examples, we've suppressed messages for more concise output, but you can set `messages=TRUE` to get helpful messages - `messages=TRUE` is also the default setting so if you don't specify that parameter messages will be printed to the console. 


```r
res <- searchbyterm(class = "Aves", state = "California", limit = 10, messages = TRUE)
```
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.6.2 Patched (2020-01-25 r77715) |
|os       |macOS Mojave 10.14.6                        |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-01-29                                  |

# Dependencies

|package  |old   |new   |Δ  |
|:--------|:-----|:-----|:--|
|rvertnet |0.7.0 |0.8.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*