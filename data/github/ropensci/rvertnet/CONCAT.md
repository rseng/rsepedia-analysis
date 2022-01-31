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

*Wow, no problems at all. :)**Wow, no problems at all. :)*rvertnet
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

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

```{r eval=FALSE}
install.packages("rvertnet")
```

Or development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/rvertnet")
```

```{r}
library('rvertnet')
```

## Search by term

Search for _Aves_ in the state of _California_, limit to 10 records

```{r}
res <- searchbyterm(class = "Aves", stateprovince = "California",
    limit = 10, messages = FALSE)
```

Inspect metadata

```{r}
res$meta
```

Inspect data. A `dplyr` data.frame is given back, so you get a nice brief data summary:

```{r}
res$data[,1:5]
```

Search for _Mustela nigripes_ in the states of _Wyoming_ or _South Dakota_, limit to 20 records

```{r}
res <- searchbyterm(specificepithet = "nigripes",
    stateprovince = "(wyoming OR south dakota)", 
    limit = 20, messages = FALSE)
res$data[,1:5]
```

### dplyr downstream

You can pass the data object directly on to `dplyr` functions. Here, we get a table of record counts by species in descending order.

```{r eval=FALSE}
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

```{r eval=FALSE}
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

```{r}
res <- spatialsearch(lat = 33.529, long = -105.694, radius = 2000,
    limit = 10, messages = FALSE)
res$data[,1:5]
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
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rvertnet introduction}
%\VignetteEncoding{UTF-8}
-->

```{r, eval=TRUE, echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

rvertnet introduction
=====================

`rvertnet` is a client for interacting with [VertNet.org](http://vertnet.org/).

## Installation

You can install the stable version from CRAN:

```{r eval=FALSE}
install.packages("rvertnet")
```

Or the development version from GitHub using the `devtools` package:

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/rvertnet")
```

```{r}
library('rvertnet')
```

## Search by term

Search for _Aves_ in the state of _California_, limit to 10 records

```{r}
res <- searchbyterm(class = "Aves", state = "California", limit = 10, messages = FALSE)
```

All major functions (`searchbyterm()`, `spatialsearch()`, `vertsearch()`) give back a `meta` (for metadata, in a list) and `data` (for data, in a data.frame) slot. The metadata:

```{r}
res$meta
```

The data

```{r}
res$data
```

Search for _Mustela nigripes_ in the states of _Wyoming_ or _South Dakota_, limit to 20 records

```{r}
res <- searchbyterm(specificepithet = "nigripes", genus = "Mustela", state = "(wyoming OR south dakota)", limit = 20, messages = FALSE)
res$data
```

Search for class _Aves_, in the state of _Nevada_, with a coordinate uncertainty range (in meters) of less than 25 meters

```{r}
res <- searchbyterm(class = "Aves", stateprovince = "Nevada", error = "<25", messages = FALSE)
res$data
```

## Spatial search

Spatial search service allows only to search on a point defined by latitude and longitude pair, with a radius (meters) from that point. All three parameters are required. 

```{r}
res <- spatialsearch(lat = 33.529, lon = -105.694, radius = 2000, limit = 10, messages = FALSE)
res$data
```

## Global full text search

`vertsearch()` provides a simple full text search against all fields. For more info see [the docs](https://github.com/VertNet/webapp/wiki/The-API-search-function#global-full-text-search). An example:

```{r}
res <- vertsearch(taxon = "aves", state = "california", limit = 10)
res$data
```

Limit the number of records returned (under 1000)

```{r}
res <- vertsearch("(kansas state OR KSU)", limit = 200)
res$data
```

Pass output of `vertsearch()` to a map

```{r fig.width=8, fig.height=4}
out <- vertsearch(tax = "(mustela nivalis OR mustela erminea)")
vertmap(out)
```

## Lots of data

For `searchbyterm()`, `spatialsearch()`, and `vertsearch()`, you can request more than 1000 records. VertNet limits each request to 1000 records, but internally in this package, if you request more than 1000 records, we'll continue to send requests to get all the records you want. See the [VertNet docs](https://github.com/VertNet/webapp/wiki/The-API-search-function#retrieving-large-result-sets) for more information on this.

## Email dump of data

`bigsearch()` specifies a termwise search (like `searchbyterm()`), but requests that all available records be made available for download as a tab-delimited text file.

```{r eval=FALSE}
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

```{r}
res <- searchbyterm(class = "Aves", state = "California", limit = 10, messages = TRUE)
```
---
title: Introduction to rvertnet
author: Scott Chamberlain
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

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

![plot of chunk unnamed-chunk-13](../man/figures/unnamed-chunk-13-1.png)

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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vert_id.R
\name{vert_id}
\alias{vert_id}
\title{Search by Vertnet occurrence ID}
\usage{
vert_id(ids, compact = TRUE, messages = TRUE, ...)
}
\arguments{
\item{ids}{(character) VertNet IDs, one or more. Required.}

\item{compact}{(logical) Return a compact data frame. That is, remove
empty columns. Default: \code{TRUE}}

\item{messages}{(logical) Print progress and information messages.
Default: \code{TRUE}}

\item{...}{Curl arguments passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list, with data frame of search results, and list of metadata
}
\description{
Search by Vertnet occurrence ID
}
\details{
VertNet IDs can be a variety of things, some URIs
(i.e., with http://...), while others start with \code{urn}.

Internally in this function we filter data to darwin core terms only. To
see what terms we use, do
\code{readLines(
system.file("extdata", "simple_dwc_terms.txt", package = "rvertnet"))}.
Get in touch with us if these terms need correcting/are out of date. The
terms are from
https://github.com/tdwg/dwc/blob/master/dist/simple_dwc_horizontal.csv
}
\examples{
\dontrun{
vert_id(ids = "urn:catalog:CM:Herps:116520")
ids <- c("http://arctos.database.museum/guid/MSB:Mamm:56979?seid=1643089", 
         "urn:catalog:CM:Herps:116520",
         "urn:catalog:AUM:Fish:13271")
res <- vert_id(ids)
res$data$occurrenceid

out <- vertsearch(taxon = "aves", state = "california", limit = 5)
(ids <- out$data$occurrenceid)
res <- vert_id(ids)
identical(sort(res$data$occurrenceid), sort(ids))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vertsummary.R
\name{vertsummary}
\alias{vertsummary}
\title{Summarize a set of records downloaded from VertNet.}
\usage{
vertsummary(input, verbose = TRUE)
}
\arguments{
\item{input}{Output from \code{\link{vertsearch}},
\code{\link{searchbyterm}}, or \code{\link{spatialsearch}}. Required.}

\item{verbose}{Print progress and information messages. Default: TRUE}
}
\value{
A list of summary statistics
}
\description{
Creates a simple summary of data returned by a VertNet search.
}
\details{
\code{\link{vertsummary}} provides information on the sources,
types and extent of data returned by a VertNet search.
}
\examples{
\dontrun{
# get occurrence records
recs <- vertsearch("Junco hyemalis", limit = 10)

# summarize occurrence records
vertsummary(recs)

vertsummary(vertsearch("Oncorhynchus clarki henshawi"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{vertavailablemaps}
\alias{vertavailablemaps}
\title{This function is defunct.}
\usage{
vertavailablemaps(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{vertlocations}
\alias{vertlocations}
\title{This function is defunct.}
\usage{
vertlocations(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigsearch.R
\name{bigsearch}
\alias{bigsearch}
\title{Request to download a large number of VertNet records.}
\usage{
bigsearch(..., rfile, email, messages = TRUE, callopts = list())
}
\arguments{
\item{...}{arguments, must be named, see \code{\link[=searchbyterm]{searchbyterm()}} for details}

\item{rfile}{A name for the results file that you will download (character).
Required.}

\item{email}{An email address where you can be contacted when your records
are ready for download (character). Required.}

\item{messages}{(logical) Print progress and information messages.
Default: \code{TRUE}}

\item{callopts}{(named list) Curl arguments passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
Prints messages on progress, but returns NULL
}
\description{
Specifies a term-wise search (like \code{\link{searchbyterm}}) and requests
that all available  records be made available for download as a
tab-delimited text file.
}
\details{
\code{\link{bigsearch}} allows you to request records as a
tab-delimited text file. This is the best way to access a large number of
records, such as when your search results indicate that >1000 records are
available. You will be notified by email when your records are ready
for download.
}
\section{Reading data}{

We suggest reading data in with \code{data.table::fread()} - as it's very
fast for the sometimes large datasets
you will get from using this function, and is usually robust to
formatting issues.
}

\examples{
\dontrun{
# replace "big@search.luv" with your own email address
bigsearch(genus = "ochotona", rfile = "pikaRecords", email = "big@search.luv")

# Pass in curl options for curl debugging
bigsearch(genus = "ochotona", rfile = "pikaRecords",
  email = "big@search.luv", verbose = TRUE)

# Use more than one year query
bigsearch(class = "aves", year = c(">=1976", "<=1986"),
          rfile = "test-bigsearch1", email = "big@search.luv")
}
}
\references{
https://github.com/VertNet/webapp/wiki/The-API-search-function
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vertmap.R
\name{vertmap}
\alias{vertmap}
\title{Make a simple map to visualize VertNet data.}
\usage{
vertmap(
  input = NULL,
  mapdatabase = "world",
  region = ".",
  geom = geom_point,
  jitter = NULL
)
}
\arguments{
\item{input}{Output from \code{\link{vertsearch}},
\code{\link{searchbyterm}}, or \code{\link{spatialsearch}}. Must
include columns "decimallatitude" and "decimallongitude"}

\item{mapdatabase}{The base map on which your data are displayed; what you
choose here determines what you can choose in the region parameter; one
of: county, state, usa, world, world2, france, italy, or nz}

\item{region}{The region in which your data are displayed; to see region names
for the "world" database layer, run
\code{sort(unique(map_data("world")$region))} after loading packages maps
and ggplot2; to see region names for the US "state" layer, run
\code{sort(unique(map_data("state")$region))}}

\item{geom}{Specifies the type of object being plotted; one of: \code{geom_point} or
\code{geom_jitter} (do not use quotes)}

\item{jitter}{If \code{geom = geom_jitter}, the amount by which to jitter points in
width, height, or both. Default}
}
\value{
Map of record locations displayed on the selected base map
}
\description{
Plots record locations on a world or regional map using latitude/longitude
data returned by a VertNet search.
}
\details{
\code{vertmap} uses decimal latitude and longitude data in records generated by
an rvertnet search to display returned records on a specified base map. Taxa
are color-coded by scientific name, if available. Adapt the vertmap code to
construct maps according to your own specifications.
}
\examples{
\dontrun{
out <- vertsearch("Junco hyemalis") # get occurrence records
vertmap(out)                        # map occurrence records

# Records are color coded by dwc term "scientificname" - sometimes unavailble
out <- vertsearch("mustela nigripes")
vertmap(input = out, mapdatabase = "state")

# Use searchbyterm() to match records with mapped region
spec <- searchbyterm(genus = "ochotona", specificepithet = "princeps", state = "california",
limit = 200)
vertmap(input = spec, mapdatabase = "state", region = "california")

# Many species
splist <- c("Accipiter erythronemius", "Aix sponsa", "Haliaeetus leucocephalus",
		"Corvus corone", "Threskiornis molucca", "Merops malimbicus")
out <- lapply(splist, function(x) vertsearch(t=x, lim=100))
library("plyr")
out <- ldply(lapply(out, "[[", "data"))
vertmap(out)
## jitter points
library("ggplot2")
vertmap(out, geom = geom_jitter, jitter = position_jitter(1, 6))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traitsearch.R
\name{traitsearch}
\alias{traitsearch}
\title{Trait focused search}
\usage{
traitsearch(
  taxon = NULL,
  has_mass = FALSE,
  has_length = FALSE,
  has_sex = FALSE,
  has_lifestage = FALSE,
  length_type = NULL,
  length = NULL,
  mass = NULL,
  limit = 1000,
  compact = TRUE,
  messages = TRUE,
  callopts = list(),
  ...
)
}
\arguments{
\item{taxon}{(character) Taxonomic identifier or other text to search for}

\item{has_mass}{(logical) limit to records that have mass data (stored in
\code{massing}). Default: \code{FALSE}}

\item{has_length}{(logical) limit to records that have length data (stored
in \code{lengthinmm}). Default: \code{FALSE}}

\item{has_sex}{(logical) limit to records that have sex data (stored in
\code{sex}). Default: \code{FALSE}}

\item{has_lifestage}{(logical) limit to records that have lifestage data
(stored in \code{lifestage}). Default: \code{FALSE}}

\item{length_type}{(character) length type, one of 'total length',
'standard length', 'snout-vent length', 'head-body length', 'fork length',
'total length range', 'standard length range', 'snout-vent length range',
'head-body length range', 'fork length range'. (stored in \code{lengthtype})
Default: \code{NULL}}

\item{length}{(list) list of query terms for length, e.g., "< 100"}

\item{mass}{(list) list of query terms for mass, e.g., "< 100"}

\item{limit}{(numeric) Limit on the number of records returned. If >1000
results, we use a cursor internally, but you should still get up to the
results you asked for. See also
\code{\link{bigsearch}} to get larger result sets in a text file via email.}

\item{compact}{Return a compact data frame (boolean)}

\item{messages}{Print progress and information messages. Default: TRUE}

\item{callopts}{curl options in a list passed on to
\code{\link[crul]{HttpClient}}, see examples}

\item{...}{(character) Additional search terms. These must be unnamed}
}
\value{
a list, same as returned by \code{\link{vertsearch}}, with data
in the \code{data} slot
}
\description{
Trait focused search
}
\details{
Wraps \code{\link{vertsearch}}, with some of the same parameters,
but with additional parameters added to make querying for traits easy.
}
\examples{
\dontrun{
traitsearch(has_mass = TRUE, limit = 3)
traitsearch(has_lifestage = TRUE)
traitsearch(has_mass = TRUE, has_length = TRUE)
res <- traitsearch(length_type = "total length",
  length = list(">= 300", "<= 1000"))
summary(as.numeric(res$data$lengthinmm))
res <- traitsearch(has_mass = TRUE, mass = list(">= 20", "<= 500"))
summary(as.numeric(res$data$massing))

traitsearch(taxon = "aves", has_mass = TRUE, limit = 100)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{vertoccurrence}
\alias{vertoccurrence}
\title{This function is defunct.}
\usage{
vertoccurrence(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{rvertnet-defunct}
\alias{rvertnet-defunct}
\title{Defunct functions in rvertnet}
\description{
\itemize{
\item \code{\link{vertavailablemaps}}: Function is now defunct, i.e., not available anymore.
\item \code{\link{vertlocations}}: Function is now defunct, i.e., not available anymore.
\item \code{\link{vertoccurrence}}: Function is now defunct, i.e., not available anymore.
\item \code{\link{vertoccurrencecount}}: Function is now defunct, i.e., not available anymore.
\item \code{\link{vertproviders}}: Function is now defunct, i.e., not available anymore.
\item \code{\link{verttaxa}}: Function is now defunct, i.e., not available anymore.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{vertproviders}
\alias{vertproviders}
\title{This function is defunct.}
\usage{
vertproviders(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvertnet-package.R
\docType{package}
\name{rvertnet-package}
\alias{rvertnet-package}
\alias{rvertnet}
\title{Search VertNet archives using R}
\description{
There are a variety of ways to search VertNet
}
\section{Search by term}{


Search for \emph{Aves} in the state of \emph{California}, limit to 10 records, e.g.:

\code{searchbyterm(class = "Aves", state = "California", lim = 10, 
verbose = FALSE)}

Search for \emph{Mustela nigripes} in the states of \emph{Wyoming} or \emph{South Dakota},
limit to 20 records, e.g.:

\code{searchbyterm(genus = "Mustela", specificepithet = "nigripes", 
   state = "(wyoming OR south dakota)", limit = 20, verbose=FALSE)}
}

\section{Big data}{

Specifies a termwise search (like \code{searchbyterm()}), but requests that all
available records  be made available for download as a tab-delimited
text file.

\code{bigsearch(genus = "ochotona", rf = "pikaRecords", 
email = "big@search.luv")}
}

\section{Spatial search}{

\code{spatialsearch(lat = 33.529, lon = -105.694, radius = 2000, limit = 10, 
verbose = FALSE)}
}

\section{Full text search}{

Find records using a global full-text search of VertNet archives.

\code{vertsearch(taxon = "aves", state = "california")}
}

\section{No results?}{


It's possible to get no results when requesting data from VertNet,
then run the same function again 10 seconds later, and you do get a result.
I'm not sure why this is, something having to do with Vertnet's
infrastucture that I'm not aware of. Point is, if you are sure
you haven't made any mistakes with the parameters, etc., then
simply run the function call again.
}

\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchbyterm.R
\name{searchbyterm}
\alias{searchbyterm}
\title{Search by term}
\usage{
searchbyterm(
  ...,
  limit = 1000,
  compact = TRUE,
  messages = TRUE,
  only_dwc = TRUE,
  callopts = list()
)
}
\arguments{
\item{...}{arguments, must be named, see section \code{Parameters} for details.
Multiple inputs to a single parameter are supported, but you have to
construct that string yourself with \code{AND} or \code{OR} operators; see
examples below.}

\item{limit}{(numeric) Limit on the number of records returned. If >1000
results, we use a cursor internally, but you should still get up to the
results you asked for. See also \code{\link[=bigsearch]{bigsearch()}} to get larger
result sets in a text file via email.}

\item{compact}{(logical) Return a compact data frame}

\item{messages}{(logical) Print progress and information messages.
Default: \code{TRUE}}

\item{only_dwc}{(logical) whether or not to return only Darwin Core term
fields. Default: \code{TRUE}}

\item{callopts}{(named list) Curl arguments passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
A list with two slots:
\itemize{
\item meta: a named list of metadata for the search results
\item data: a data frame of search results, columns vary
}
}
\description{
Flexible search for records using keywords/terms
}
\details{
\code{searchbyterm()} builds a query from input parameters
based on Darwin Core (dwc) terms (for the full list of terms, see
https://code.google.com/p/darwincore/wiki/DarwinCoreTerms).
}
\section{Parameters}{


All these parameters can be passed in to \code{searchbyterm()}. All others
will be silently dropped.

See https://github.com/VertNet/webapp/wiki/The-API-search-function
for more details

\strong{taxon}
\itemize{
\item kingdom (character) Taxonomic kingdom
\item phylum (character) Taxonomic phylum
\item class (character) Taxonomic class
\item order (character) Taxonomic order
\item family (character) Taxonomic family
\item genus (character) Taxonomic genus
\item specificepithet (character) Taxonomic specific epithet, e.g. (sapiens
in Homo sapiens)
\item infraspecificepithet (character) Taxonomic infraspecific epithet
\item scientificname (character) scientific name
\item vernacularname (character) a verncular name
}

\strong{event}
\itemize{
\item year (numeric) Year or range of years designated by comparison
operators "<", ">", "<=" or ">=". You can pass in more than one of these
queries, in a vector. See example below
\item month (numeric) month or range of months designated by comparison
operators "<", ">", "<=" or ">=". You can pass in more than one of these
queries, in a vector. See example below
\item day (numeric) day or range of days designated by comparison
operators "<", ">", "<=" or ">=". You can pass in more than one of these
queries, in a vector. See example below
\item eventdate Event date associated with this occurrence record; yyyy-mm-dd
or the range yyyy-mm-dd/yyyy-mm-dd (character)
}

\strong{record level}
\itemize{
\item institutioncode (character) Code name for the provider/institution
of record
\item occurrenceid (character) Provider's unique identifier for this
occurrence record
\item catalognumber (character) Provider's catalog number or other ID for
this record
\item collectioncode (character) collection code
\item license (dcterms:license) license string
\item iptlicense (eml:intellectualRights) license string
\item basisofrecord (character) one of PreservedSpecimen, FossilSpecimen,
MaterialSample, Occurrence, MachineObservation, HumanObservation
\item hasmedia (logical) Record also references associated media, such as a
film or video
\item isfossil (logical) \code{dwc:basisOfRecord} is FossilSpecimen or collection
is a paleo collection
\item haslicense (logical) \code{dcterms:license} or \code{eml:intellectualRights} has a
license designated
}

\strong{identification}
\itemize{
\item typestatus (character) a type status
\item hastypestatus (logical) type status known or not
}

\strong{occurrence}
\itemize{
\item iptrecordid (character) (same as \code{dwc:occurrenceID})
\item recordedby (character) Collector name
\item recordnumber (character) record number
\item fieldnumber (character) field number
\item establishmentmeans (character) establishment means
\item wascaptive (logical) (\code{dwc:establishmentMeans} or occurrenceRemarks suggests
it was captive)
\item wasinvasive (logical) (was the organism recorded to be invasive where
and when it occurred)
\item sex (character) standardized sex from original sex field or extracted
from elsewhere in the record
\item lifestage (character) lifeStage from original sex field or extracted
from elsewhere in the record
\item preparations (not sure what this means)
\item hastissue (logical) Record is likely to reference tissues
\item reproductivecondition (not sure what this means)
}

\strong{location}
\itemize{
\item continent (character) Continent to search for occurrence
\item country (character) Country to search for occurrence
\item stateprovince (character) State or province to search for occurrence
\item county (character) County to search for occurrence
\item island (character) Island to search for occurrence
\item igroup (character) Island group to search for occurrence
\item municipality (character)
\item waterbody (character)
\item geodeticdatum (character)
\item georeferencedby (character)
\item georeferenceverificationstatus (character)
\item location a Google GeoField of the \code{dwc:decimalLatitude},
\code{dwc:decimalLongitude}
\item mappable (logical) Record includes valid coordinates in decimal latitude
and decimal longitude
}

\strong{geological context}
\itemize{
\item bed (character) geological bed
\item formation (character) geological formation
\item group (character) geological group
\item member (character) geological member
}

\strong{traits}
\itemize{
\item haslength (logical) (was a value for length extracted?)
\item hasmass (logical) (was a value for mass extraccted?)
\item hassex (logical) (does the record have sex?)
\item haslifestage (logical) (does the record have life stage?)
\item lengthtype (character) type of length measurement extracted from the
record, can refer to a number or to a range) ('total length',
'standard length', 'snout-vent length','head-body length', 'fork length',
'total length range', 'standard length range', 'snout-vent length range',
'head-body length range', 'fork length range'
\item lengthinmm (numeric) length measurement extracted from the record
\item massing (numeric) mass measurement extracted from the record (For
detailed information about trait extraction and aggregation and querying
via the VertNet portal, see http://vertnet.org/resources/traitsguide.html
}

\strong{data set}
\itemize{
\item gbifdatasetid (character) GBIF identifier for the data set
\item gbifpublisherid (character) GBIF identifier for the data publishing
organization
\item lastindexed (character) date (YYYY-MM-DD) the record was most recently
indexed into VertNet
\item networks (character) one of MaNIS, ORNIS, HerpNET, FishNet, VertNet,
Arctos, Paleo
\item migrator (character) the version of the migrator used to process the
data set, a date of form (YYYY-MM-DD)
\item orgcountry (character) the country where the organization is located
\item orgstateprovince (character) the first-level administrative unit where
the organization is located
}

\strong{index}
\itemize{
\item rank (character) a higher number means the record is more complete
with respect to georeferences, scientific names, and event dates
\item vntype (character) Type of record; "specimen" or "observation"
\item hashid (integer) a value to distribute records in 10k bins; 0-9998
}

\strong{other}
\itemize{
\item coordinateuncertaintyinmeters (character) Coordinate uncertainty
in meters (numeric) or range of uncertainty values designated by
comparison operators "<", ">", "<=", or ">="
}
}

\section{No results?}{


It's possible to get no results with a call to \code{searchbyterm()},
then run it again 10 seconds later, and you do get a result.
I'm not sure why this is, something having to do with Vertnet's
infrastucture that I'm not aware of. Point is, if you are sure
you haven't made any mistakes with the parameters, etc., then
simply run the function call again.
}

\examples{
\dontrun{
# Find multiple species
out <- searchbyterm(genus = "ochotona",
  specificepithet = "(princeps OR collaris)", limit=10)

# iptrecordid
searchbyterm(iptrecordid = "7108667e-1483-4d04-b204-6a44a73a5219")

# you can pass more than one, as above, in a single string in parens
records <- "(7108667e-1483-4d04-b204-6a44a73a5219 OR 1efe900e-bde2-45e7-9747-2b2c3e5f36c3)"
searchbyterm(iptrecordid = records, callopts = list(verbose = TRUE))

# Specifying a range (in meters) for uncertainty in spatial location
# (use quotes)
out <- searchbyterm(class = "aves", stateprovince = "nevada", 
  coordinateuncertaintyinmeters = "<25")
out <- searchbyterm(class = "aves", stateprovince = "california", year = 1976,
  coordinateuncertaintyinmeters = "<=1000")

# Specifying records by event date (use quotes)
out <- searchbyterm(class = "aves", stateprovince = "california",
  eventdate = "2009-03-25")
# ...but specifying a date range may not work
out <- searchbyterm(specificepithet = "nigripes",
  eventdate = "1935-09-01/1935-09-30")

# Pass in curl options for curl debugging
out <- searchbyterm(class = "aves", limit = 10,
 callopts = list(verbose = TRUE))

# Use more than one year query
searchbyterm(genus = "mustela", specificepithet = "nigripes",
   year = c('>=1900', '<=1940'))

searchbyterm(sex  = "male", limit = 30)$data$sex
searchbyterm(lifestage  = "juvenile", limit = 30)$data$lifestage
}
}
\references{
https://github.com/VertNet/webapp/wiki/The-API-search-function
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{verttaxa}
\alias{verttaxa}
\title{This function is defunct.}
\usage{
verttaxa(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vertsearch.R
\name{vertsearch}
\alias{vertsearch}
\title{Find records using a global full-text search of VertNet archives.}
\usage{
vertsearch(
  taxon = NULL,
  ...,
  limit = 1000,
  compact = TRUE,
  messages = TRUE,
  only_dwc = TRUE,
  callopts = list()
)
}
\arguments{
\item{taxon}{(character) Taxonomic identifier or other text to search for}

\item{...}{(character) Additional search terms. These must be unnamed}

\item{limit}{(numeric) Limit on the number of records returned. If >1000
results, we use a cursor internally, but you should still get up to the
results you asked for. See also
\code{\link{bigsearch}} to get larger result sets in a text file via email.}

\item{compact}{Return a compact data frame (boolean)}

\item{messages}{Print progress and information messages. Default: TRUE}

\item{only_dwc}{(logical) whether or not to return only Darwin Core term
fields. Default: \code{TRUE}}

\item{callopts}{curl options in a list passed on to
\code{\link[crul]{HttpClient}}, see examples}
}
\value{
A data frame of search results
}
\description{
Returns any record containing your target text in any field of the record.
}
\details{
\code{\link{vertsearch}} performs a nonspecific search for your
input within every record and field of the VertNet archives. For a more
specific search, try \code{\link{searchbyterm}}
}
\examples{
\dontrun{
out <- vertsearch(taxon = "aves", "california", limit=3)

# Limit the number of records returned (under 1000)
out <- vertsearch("(kansas state OR KSU)", limit = 200)
# Use bigsearch() to retrieve >1000 records

# Find multiple species using searchbyterm():
# a) returns a specific result
out <- searchbyterm(genus = "mustela", species = "(nivalis OR erminea)")
vertmap(out)

# b) returns a non-specific result
out <- vertsearch(taxon = "(mustela nivalis OR mustela erminea)")
vertmap(out)

# c) returns a non-specific result
splist <- c("mustela nivalis", "mustela erminea")
out <- lapply(splist, function(x) vertsearch(taxon = x, lim = 500))
library("plyr")
out <- ldply(lapply(out, "[[", "data"))
vertmap(out)

# curl options
vertsearch(taxon = "Aves", limit = 10, callopts = list(verbose = TRUE))
# vertsearch(taxon = "Aves", limit = 10, callopts = list(timeout_ms = 10))
}
}
\references{
\url{https://github.com/VertNet/webapp/wiki/The-API-search-function}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vertdump.R
\name{dump-defunct}
\alias{dump-defunct}
\alias{dump_init}
\alias{dump_tbl}
\alias{dump_links}
\title{These functions are defunct}
\usage{
dump_init(...)

dump_tbl(...)

dump_links(...)
}
\description{
These functions are defunct
}
\details{
The Vertnet dumps are too complicated to setup. Please get in
touch with Vertnet or maintainers of this package if you want help using
bulk data.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{vertoccurrencecount}
\alias{vertoccurrencecount}
\title{This function is defunct.}
\usage{
vertoccurrencecount(...)
}
\description{
This function is defunct.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialsearch.R
\name{spatialsearch}
\alias{spatialsearch}
\title{Find records within some distance of a point given latitude and longitude.}
\usage{
spatialsearch(
  lat,
  long,
  radius,
  limit = 1000,
  compact = TRUE,
  messages = TRUE,
  ...
)
}
\arguments{
\item{lat}{(numeric) Latitude of the central point, in decimal degrees
required.}

\item{long}{(numeric) Longitude of the central point, in decimal degrees
required.}

\item{radius}{(numeric) Radius to search, in meters. There is no default
value for this parameter. required.}

\item{limit}{(integer) Limit on the number of records returned. If >1000
results, we use a cursor internally, but you should still get up to the
results you asked for. See also \code{\link[=bigsearch]{bigsearch()}} to get larger result
sets in a text file via email.}

\item{compact}{(logical) Return a compact data frame. default: \code{TRUE}}

\item{messages}{(logical) Print progress and information messages.
Default: \code{TRUE}}

\item{...}{Curl arguments passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list with two slots:
\itemize{
\item meta: a named list of metadata for the search results
\item data: a data frame of search results, columns vary
}
}
\description{
Searches by decimal latitude and longitude to return any occurrence record
within the input distance (radius) of the input point.
}
\details{
\code{\link[=spatialsearch]{spatialsearch()}} finds all records of any taxa having
decimal lat/long coordinates within a given radius (in meters) of
your coordinates.
}
\examples{
\dontrun{
res <- spatialsearch(lat = 33.529, long = -105.694, radius = 2000,
  limit = 10)

# Pass in curl options for curl debugging
out <- spatialsearch(lat = 33.529, long = -105.694, radius = 2000,
  limit = 10, verbose = TRUE)
}
}
\references{
https://github.com/VertNet/webapp/wiki/The-API-search-function
}
