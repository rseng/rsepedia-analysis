pangaear
========



[![cran checks](https://cranchecks.info/badges/worst/pangaear)](https://cranchecks.info/pkgs/pangaear)
[![R-check](https://github.com/ropensci/pangaear/workflows/R-check/badge.svg)](https://github.com/ropensci/pangaear/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/pangaear/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/pangaear)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/pangaear)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/pangaear)](https://cran.r-project.org/package=pangaear)

`pangaear` is a data retrieval interface for the World Data Center PANGAEA (https://www.pangaea.de/). PANGAEA archieves published Earth & Environmental Science data under the following subjects: agriculture, atmosphere, biological classification, biosphere, chemistry, cryosphere, ecology, fisheries, geophysics, human dimensions, lakes & rives, land surface, lithosphere, oceans, and paleontology.

This package offers tools to interact with the PANGAEA Database, including functions for searching for data, fetching datasets by dataset ID, and working with the PANGAEA OAI-PMH service.

## Info

* Pangaea [website](https://www.pangaea.de/).
* Pangaea [OAI-PMH docs](https://wiki.pangaea.de/wiki/OAI-PMH).
* [OAI-PMH Spec](http://www.openarchives.org/OAI/openarchivesprotocol.html)

## Package API

 - pg_data
 - pg_list_metadata_formats
 - pg_identify
 - pg_list_records
 - pg_list_sets
 - pg_list_identifiers
 - pg_search
 - pg_get_record
 - pg_cache
 - pg_search_es
 - pg_cache_list
 - pg_cache_clear

## Installation

Stable version


```r
install.packages("pangaear")
```

Dev version


```r
remotes::install_github('ropensci/pangaear')
```


```r
library('pangaear')
```

## Search for data

This is a thin wrapper around the GUI search interface on the page <https://www.pangaea.de/>. Everything you can do there, you can do here.


```r
pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
#> # A tibble: 3 x 6
#>   score doi       size size_measure citation             supplement_to          
#>   <dbl> <chr>    <dbl> <chr>        <chr>                <chr>                  
#> 1 13.3  10.1594…     4 datasets     Krylova, EM; Sahlin… Krylova, EM; Sahling, …
#> 2 13.1  10.1594…     2 datasets     Simonyan, AV; Dultz… Simonyan, AV; Dultz, S…
#> 3  8.93 10.1594…   598 data points  WOCE Surface Veloci… <NA>
```

## Get data


```r
res <- pg_data(doi = '10.1594/PANGAEA.807580')
res[[1]]
#> <Pangaea data> 10.1594/PANGAEA.807580
#>   parent doi: 10.1594/PANGAEA.807580
#>   url:        https://doi.org/10.1594/PANGAEA.807580
#>   citation:   Schiebel, Ralf; Waniek, Joanna J; Bork, Matthias; Hemleben, Christoph (2001): Physical oceanography during METEOR cruise M36/6. PANGAEA, https://doi.org/10.1594/PANGAEA.807580, In supplement to: Schiebel, R et al. (2001): Planktic foraminiferal production stimulated by chlorophyll redistribution and entrainment of nutrients. Deep Sea Research Part I: Oceanographic Research Papers, 48(3), 721-740, https://doi.org/10.1016/S0967-0637(00)00065-0
#>   path:       /Users/sckott/Library/Caches/R/pangaear/10_1594_PANGAEA_807580.txt
#>   data:
#> # A tibble: 32,179 x 13
#>    Event       `Date/Time`   Latitude Longitude `Elevation [m]` `Depth water [m…
#>    <chr>       <chr>            <dbl>     <dbl>           <int>            <dbl>
#>  1 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             0   
#>  2 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             0.99
#>  3 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             1.98
#>  4 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             2.97
#>  5 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             3.96
#>  6 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             4.96
#>  7 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             5.95
#>  8 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             6.94
#>  9 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             7.93
#> 10 M36/6-CTD-… 1996-10-14T1…     49.0     -16.5           -4802             8.92
#> # … with 32,169 more rows, and 7 more variables: Press [dbar] <int>,
#> #   Temp [°C] <dbl>, Sal <dbl>, Tpot [°C] <dbl>, Sigma-theta [kg/m**3] <dbl>,
#> #   Sigma in situ [kg/m**3] <dbl>, Cond [mS/cm] <dbl>
```

Search for data then pass DOI to data function.


```r
res <- pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
pg_data(res$doi[3])[1:3]
#> [[1]]
#> <Pangaea data> 10.1594/PANGAEA.406110
#>   parent doi: 10.1594/PANGAEA.406110
#>   url:        https://doi.org/10.1594/PANGAEA.406110
#>   citation:   WOCE Surface Velocity Program, SVP (2006): Water temperature and current velocity from surface drifter SVP_9616641. PANGAEA, https://doi.org/10.1594/PANGAEA.406110
#>   path:       /Users/sckott/Library/Caches/R/pangaear/10_1594_PANGAEA_406110.txt
#>   data:
#> # A tibble: 101 x 10
#>    `Date/Time`  Latitude Longitude `Depth water [m… `Temp [°C]` `Cur vel U [cm/…
#>    <chr>           <dbl>     <dbl>            <int>       <dbl>            <dbl>
#>  1 1996-11-11T…     41.6     -125.                0        12.7            NA   
#>  2 1996-11-11T…     41.6     -125.                0        12.5            11.3 
#>  3 1996-11-12T…     41.6     -124.                0        12.4             2.91
#>  4 1996-11-12T…     41.7     -124.                0        12.3             3.64
#>  5 1996-11-12T…     41.7     -124.                0        11.9            23.4 
#>  6 1996-11-12T…     41.7     -124.                0        11.4            21.4 
#>  7 1996-11-13T…     41.8     -124.                0        11.1             0.21
#>  8 1996-11-13T…     41.8     -124.                0        11.2            -0.86
#>  9 1996-11-13T…     41.8     -124.                0        11.1             1.51
#> 10 1996-11-13T…     41.9     -124.                0        11.0            -5.58
#> # … with 91 more rows, and 4 more variables: Cur vel V [cm/s] <dbl>,
#> #   Latitude e <dbl>, Longitude e <dbl>, Code <int>
#> 
#> [[2]]
#> NULL
#> 
#> [[3]]
#> NULL
```

## OAI-PMH metadata


```r
# Identify the service
pg_identify()

# List metadata formats
pg_list_metadata_formats()

# List identifiers
pg_list_identifiers(from = Sys.Date() - 2, until = Sys.Date())

# List sets
pg_list_sets()

# List records
pg_list_records(from = Sys.Date() - 1, until = Sys.Date())

# Get a record
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")
```

## Contributors (reverse alphabetical)

* Naupaka Zimmerman
* Kara Woo
* Gavin Simpson
* Andrew MacDonald
* Scott Chamberlain

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pangaear/issues).
* License: MIT
* Get citation information for `pangaear` in R doing `citation(package = 'pangaear')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
pangaear 1.1.0
==============

## MINOR IMPROVEMENTS

* add markdown pkg to Suggests for vignette


pangaear 1.0.0
==============

## MINOR IMPROVEMENTS

* vcr caching for pg_data tests, local only as files are too big (#74)
* vignette title fix (#73)
* `pg_search` fix for searching with a bounding box that crosses 180/-180 longitude: not a fix in this package, but the remote data source fixed a problem reported from a user of this package (#71)

pangaear 0.8.2
==============

## BUG FIXES

* detected from cran checks, only failing on debian clang r-devel: change `pg_search` tests to use `preserve_exact_body_bytes = TRUE` in the `vcr::use_cassette()` call so `yaml` package doesn't fail on loading it (#72)


pangaear 0.8.0
==============

## NEW FEATURES

* now using package `hoardr` for managing caching, replacing package `rappdirs`; new object `pg_cache`, an R6 class, with methods for managing where you cache files, deleting them, listing them, etc. Importantly, you can set the full cache path now, not just the folder within a set directory (#66) (#69)

## MINOR IMPROVEMENTS

* most tests using `vcr` now for http reqeust caching (#70)
* output of `pg_data()` now includes a `metadata` slot with parsed metadata from text files (only included when the file is a txt file); `parameters` slot that's part of the metadata is partially parsed into an unnamed list (#67)

## BUG FIXES

* fix in `pg_data()`: Pangaea changed how links are organized for datasets, fixed now (#65)
* fix in `pg_data()`: response content type header changed - a space was added, breaking the content type check; now not sensitive to the space (#68)


pangaear 0.6.0
==============

## MINOR IMPROVEMENTS

* Added a vignette (#62)
* replaced `httr` pkg with `crul` for HTTP requests (#55)
* Some datasets require login on the Pangaea platform - we now detect this in `pg_data` and skip the file download with a message saying so. Eventually we hope to fix this to allow uesrs to input credentials to get the file(s) (#59)


pangaear 0.3.0
==============

## NEW FEATURES

* New function `pg_search_es` - an interface to Pangaea's Elasticsearch
query interface.

## MINOR IMPROVEMENTS

* added more tests (#57)
* now using markdown in docs 
* tidy code to 80 line width

## BUG FIXES

* fixed bug in oai functions due to changed base url for the 
Pangaea OAI server (#53)
* Fix to `pg_search` as search portal now has offset param - so if
more than 500 results need to page through them with the 
`offset` parameter (#56)


pangaear 0.2.4
==============

## MINOR IMPROVEMENTS

* Improved examples in the OAI methods to make sure they work 
regardless of when the user runs them

## BUG FIXES

* Fixes to `pg_search()` needed due to changes in the Pangaea 
website. Nearly identical functionality, but one parameter switch
(`env` param is now `topic`). More parameters may be added in the 
future. New fields are added to output, added to docs for the 
function. Now importing `jsonlite` as a result of these changes (#50)
* Fixes to `pg_data()` needed due to changes in the Pangaea 
website. Nothing should be different for the user. (#51)


pangaear 0.2.0
==============

## MINOR IMPROVEMENTS

* Dropped `XML`, using `xml2` now (#42)
* Using `rappdirs` package now for determing caching path on user's machine (#46)
* using `tibble` now for compact data.frame representations (#47)
* Dropped `methods` dependency (#48)
* Added test suite (#45)

## BUG FIXES

* Fixes to `pg_search()` - introduced due to changes in Pangaea website (#44)

pangaear 0.1.0
==============

* Released to CRAN.
## Test environments

* local macOS install, R 4.0.5 patched
* ubuntu 14.04 (on GitHub Actions), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version fixes the rmarkdown/markdown dependency issue for vignettes that Kurt emailed maintainers about.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/pangaear/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/pangaear.git`
* Make sure to track progress upstream (i.e., on our version of `pangaear` at `ropensci/pangaear`) by doing `git remote add upstream https://github.com/ropensci/pangaear.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/pangaear`
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Introduction to pangaear}
%\VignetteEncoding{UTF-8}
-->



Introduction to pangaear
========================

`pangaear` is a data retrieval interface for the World Data Center PANGAEA (https://www.pangaea.de/). PANGAEA archieves published Earth & Environmental Science data under the following subjects: agriculture, atmosphere, biological classification, biosphere, chemistry, cryosphere, ecology, fisheries, geophysics, human dimensions, lakes & rives, land surface, lithosphere, oceans, and paleontology.

## Installation

If you've not installed it yet, install from CRAN:


```r
install.packages("pangaear")
```

Or the development version:


```r
devtools::install_github("ropensci/pangaear")
```

## Load pangaear


```r
library("pangaear")
```

## Search for data

`pg_search` is a thin wrapper around the GUI search interface on the page <https://www.pangaea.de/>. Everything you can do there, you can do here.

For example, query for the term 'water', with a bounding box, and return only three results.


```r
pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
```

```
#> # A tibble: 3 x 6
#>   score doi            size size_measure citation                                                 supplement_to                                                                      
#>   <dbl> <chr>         <dbl> <chr>        <chr>                                                    <chr>                                                                              
#> 1 13.0  10.1594/PANG…     4 datasets     Krylova, EM; Sahling, H; Janssen, R (2010): A new genus… Krylova, EM; Sahling, H; Janssen, R (2010): Abyssogena: a new genus of the family …
#> 2 12.8  10.1594/PANG…     2 datasets     Simonyan, AV; Dultz, S; Behrens, H (2012): Diffusion tr… Simonyan, AV; Dultz, S; Behrens, H (2012): Diffusive transport of water in porous …
#> 3  8.78 10.1594/PANG…  1148 data points  WOCE Surface Velocity Program, SVP (2006): Water temper… <NA>
```

The resulting `data.frame` has details about different studies, and you can use the DOIs (Digital Object Identifiers) to get data and metadata for any studies you're interested in.

### Another search option

There's another search option with the `pg_search_es` function. It is an interface to the Pangaea Elasticsearch interface. This provides a very flexible interface for search Pangaea data - though it is different from what you're used to with the Pangaea website. 


```r
(res <- pg_search_es())
```

```
#> # A tibble: 10 x 46
#>    `_index` `_type` `_id` `_score` `_source.intern… `_source.parent… `_source.minEle… `_source.sf-aut… `_source.parent… `_source.techKe… `_source.geocod… `_source.sp-log…
#>  * <chr>    <chr>   <chr>    <dbl> <chr>            <chr>                       <dbl> <chr>                       <int> <list>           <list>                      <int>
#>  1 pangaea… panmd   9018…        1 2020-01-17T15:1… https://doi.org…             -0.4 Anhaus Philipp#…           901247 <chr [654]>      <chr [5]>                       1
#>  2 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              2   Anhaus Philipp#…           901247 <chr [652]>      <chr [5]>                       1
#>  3 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              2   Anhaus Philipp#…           901247 <chr [652]>      <chr [5]>                       1
#>  4 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              2   Anhaus Philipp#…           901247 <chr [652]>      <chr [5]>                       1
#>  5 pangaea… panmd   9017…        1 2020-01-17T15:1… <NA>                          2   Vuilleumier Lau…               NA <chr [51]>       <chr [6]>                       3
#>  6 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.pan…              0.3 Peeken Ilka#Mur…           901742 <chr [101]>      <chr [6]>                       3
#>  7 pangaea… panmd   9017…        1 2020-01-17T15:1… <NA>                         NA   Schreuder Laura…               NA <chr [37]>       <chr [3]>                       1
#>  8 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              0   Schreuder Laura…           901739 <chr [39]>       <chr [7]>                       1
#>  9 pangaea… panmd   9016…        1 2020-01-17T15:1… <NA>                         10   Augustine John$…               NA <chr [25]>       <chr [6]>                       3
#> 10 pangaea… panmd   9016…        1 2020-01-17T15:1… <NA>                          2   Augustine John$…               NA <chr [29]>       <chr [6]>                       3
#> # … with 34 more variables: `_source.agg-campaign` <list>, `_source.agg-author` <list>, `_source.eastBoundLongitude` <dbl>, `_source.URI` <chr>, `_source.agg-pubYear` <int>,
#> #   `_source.minDateTime` <chr>, `_source.agg-geometry` <chr>, `_source.xml-thumb` <chr>, `_source.xml` <chr>, `_source.sf-idDataSet` <int>, `_source.elevationGeocode` <chr>,
#> #   `_source.agg-method` <list>, `_source.maxDateTime` <chr>, `_source.xml-sitemap` <chr>, `_source.westBoundLongitude` <dbl>, `_source.northBoundLatitude` <dbl>,
#> #   `_source.sp-dataStatus` <int>, `_source.nDataPoints` <int>, `_source.sp-hidden` <lgl>, `_source.agg-location` <list>, `_source.internal-source` <chr>,
#> #   `_source.agg-basis` <chr>, `_source.southBoundLatitude` <dbl>, `_source.boost` <dbl>, `_source.oaiSet` <list>, `_source.maxElevation` <dbl>, `_source.agg-mainTopic` <list>,
#> #   `_source.agg-topic` <list>, `_source.agg-project` <list>, `_source.meanPosition.lat` <dbl>, `_source.meanPosition.lon` <dbl>, `_source.geoCoverage.type` <chr>,
#> #   `_source.geoCoverage.coordinates` <list>, `_source.geoCoverage.geometries` <list>
```

The returned data.frame has a lot of columns. You can limit columns returned with the `source` parameter.

There are attributes on the data.frame that give you the total number of results found as well as the max score found. 


```r
attributes(res)
```

```
#> $names
#>  [1] "_index"                          "_type"                           "_id"                             "_score"                          "_source.internal-datestamp"     
#>  [6] "_source.parentURI"               "_source.minElevation"            "_source.sf-authortitle"          "_source.parentIdDataSet"         "_source.techKeyword"            
#> [11] "_source.geocodes"                "_source.sp-loginOption"          "_source.agg-campaign"            "_source.agg-author"              "_source.eastBoundLongitude"     
#> [16] "_source.URI"                     "_source.agg-pubYear"             "_source.minDateTime"             "_source.agg-geometry"            "_source.xml-thumb"              
#> [21] "_source.xml"                     "_source.sf-idDataSet"            "_source.elevationGeocode"        "_source.agg-method"              "_source.maxDateTime"            
#> [26] "_source.xml-sitemap"             "_source.westBoundLongitude"      "_source.northBoundLatitude"      "_source.sp-dataStatus"           "_source.nDataPoints"            
#> [31] "_source.sp-hidden"               "_source.agg-location"            "_source.internal-source"         "_source.agg-basis"               "_source.southBoundLatitude"     
#> [36] "_source.boost"                   "_source.oaiSet"                  "_source.maxElevation"            "_source.agg-mainTopic"           "_source.agg-topic"              
#> [41] "_source.agg-project"             "_source.meanPosition.lat"        "_source.meanPosition.lon"        "_source.geoCoverage.type"        "_source.geoCoverage.coordinates"
#> [46] "_source.geoCoverage.geometries" 
#> 
#> $row.names
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $class
#> [1] "tbl_df"     "tbl"        "data.frame"
#> 
#> $total
#> [1] 390620
#> 
#> $max_score
#> [1] 1
```

```r
attr(res, "total")
```

```
#> [1] 390620
```

```r
attr(res, "max_score")
```

```
#> [1] 1
```

To get to the DOIs for each study, use 


```r
gsub("https://doi.org/", "", res$`_source.URI`)
```

```
#>  [1] "10.1594/PANGAEA.901810"                        "10.1594/PANGAEA.901736"                        "10.1594/PANGAEA.901733"                       
#>  [4] "10.1594/PANGAEA.901732"                        "https://doi.pangaea.de/10.1594/PANGAEA.901710" "https://doi.pangaea.de/10.1594/PANGAEA.901709"
#>  [7] "10.1594/PANGAEA.901739"                        "10.1594/PANGAEA.901738"                        "https://doi.pangaea.de/10.1594/PANGAEA.901697"
#> [10] "https://doi.pangaea.de/10.1594/PANGAEA.901695"
```


## Get data

The function `pg_data` fetches datasets for studies by their DOIs.


```r
res <- pg_data(doi = '10.1594/PANGAEA.807580')
res[[1]]
```

```
#> <Pangaea data> 10.1594/PANGAEA.807580
#>   parent doi: 10.1594/PANGAEA.807580
#>   url:        https://doi.org/10.1594/PANGAEA.807580
#>   citation:   Schiebel, Ralf; Waniek, Joanna J; Bork, Matthias; Hemleben, Christoph (2001): Physical oceanography during METEOR cruise M36/6. PANGAEA, https://doi.org/10.1594/PANGAEA.807580, In supplement to: Schiebel, R et al. (2001): Planktic foraminiferal production stimulated by chlorophyll redistribution and entrainment of nutrients. Deep Sea Research Part I: Oceanographic Research Papers, 48(3), 721-740, https://doi.org/10.1016/S0967-0637(00)00065-0
#>   path:       /Users/sckott/Library/Caches/R/pangaear/10_1594_PANGAEA_807580.txt
#>   data:
#> # A tibble: 32,179 x 13
#>    Event       `Date/Time`    Latitude Longitude `Elevation [m]` `Depth water [m… `Press [dbar]` `Temp [°C]`   Sal `Tpot [°C]` `Sigma-theta [kg/m… `Sigma in situ [kg… `Cond [mS/cm]`
#>    <chr>       <chr>             <dbl>     <dbl>           <int>            <dbl>          <int>       <dbl> <dbl>       <dbl>               <dbl>               <dbl>          <dbl>
#>  1 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             0                 0        15.7  35.7        15.7                26.4                26.4           44.4
#>  2 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             0.99              1        15.7  35.7        15.7                26.4                26.4           44.4
#>  3 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             1.98              2        15.7  35.7        15.7                26.4                26.4           44.4
#>  4 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             2.97              3        15.7  35.7        15.7                26.4                26.4           44.4
#>  5 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             3.96              4        15.7  35.7        15.7                26.4                26.4           44.4
#>  6 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             4.96              5        15.7  35.7        15.7                26.4                26.4           44.4
#>  7 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             5.95              6        15.7  35.7        15.7                26.4                26.4           44.4
#>  8 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             6.94              7        15.7  35.7        15.7                26.4                26.4           44.4
#>  9 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             7.93              8        15.7  35.7        15.7                26.4                26.4           44.4
#> 10 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             8.92              9        15.7  35.7        15.7                26.4                26.4           44.4
#> # … with 32,169 more rows
```

Search for data then pass one or more DOIs to the `pg_data` function.


```r
res <- pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
pg_data(res$doi[3])[1:3]
```

```
#> [[1]]
#> <Pangaea data> 10.1594/PANGAEA.405695
#>   parent doi: 10.1594/PANGAEA.405695
#>   url:        https://doi.org/10.1594/PANGAEA.405695
#>   citation:   WOCE Surface Velocity Program, SVP (2006): Water temperature and current velocity from surface drifter SVP_9524470. PANGAEA, https://doi.org/10.1594/PANGAEA.405695
#>   path:       /Users/sckott/Library/Caches/R/pangaear/10_1594_PANGAEA_405695.txt
#>   data:
#> # A tibble: 192 x 10
#>    `Date/Time`      Latitude Longitude `Depth water [m]` `Temp [°C]` `UC [cm/s]` `VC [cm/s]` `Lat e` `Lon e`  Code
#>    <chr>               <dbl>     <dbl>             <int>       <dbl>       <dbl>       <dbl>   <dbl>   <dbl> <int>
#>  1 1995-12-21T18:00     42.9     -125.                 0        12.4       NA          NA     0.0001  0.0001     1
#>  2 1995-12-22T00:00     42.9     -125.                 0        12.4        4.24      -35.2   0       0          1
#>  3 1995-12-22T06:00     42.8     -125.                 0        12.4       -7.8       -38.7   0.0002  0.0001     1
#>  4 1995-12-22T12:00     42.7     -125.                 0        12.4      -12.7       -23.1   0.0001  0.0001     1
#>  5 1995-12-22T18:00     42.7     -125.                 0        12.4      -15.1       -13.7   0       0          1
#>  6 1995-12-23T00:00     42.7     -125.                 0        12.3      -24.1        -9.15  0.0001  0.0001     1
#>  7 1995-12-23T06:00     42.7     -125.                 0        12.3      -38.4         0.65  0.0001  0.0001     1
#>  8 1995-12-23T12:00     42.7     -125.                 0        12.3      -37.2        15.8   0       0          1
#>  9 1995-12-23T18:00     42.7     -125.                 0        12.3      -25.4        29.6   0.0001  0.0001     1
#> 10 1995-12-24T00:00     42.8     -125.                 0        12.4      -18.5        35.1   0.0002  0.0002     1
#> # … with 182 more rows
#> 
#> [[2]]
#> NULL
#> 
#> [[3]]
#> NULL
```


## OAI-PMH metadata

[OAI-PMH](https://wiki.pangaea.de/wiki/OAI-PMH) is a standard protocol for serving metadata around objects, in this case datasets. If you are already familiar with OAI-PMH you are in luck as you can can use what you know here. If not familiar, it's relatively straight-forward. 

Note that you can't get data through these functions, rather only metadata about datasets.

### Identify the service


```r
pg_identify()
```

```
#> <Pangaea>
#>   repositoryName: PANGAEA - Data Publisher for Earth & Environmental Science
#>   baseURL: https://ws.pangaea.de/oai/provider
#>   protocolVersion: 2.0
#>   adminEmail: tech@pangaea.de
#>   adminEmail: tech@pangaea.de
#>   earliestDatestamp: 2015-01-01T00:00:00Z
#>   deletedRecord: transient
#>   granularity: YYYY-MM-DDThh:mm:ssZ
#>   compression: gzip
#>   description: oaipangaea.de:oai:pangaea.de:doi:10.1594/PANGAEA.999999
```

### List metadata formats


```r
pg_list_metadata_formats()
```

```
#>   metadataPrefix                                                  schema                           metadataNamespace
#> 1         oai_dc          http://www.openarchives.org/OAI/2.0/oai_dc.xsd http://www.openarchives.org/OAI/2.0/oai_dc/
#> 2         pan_md       http://ws.pangaea.de/schemas/pangaea/MetaData.xsd              http://www.pangaea.de/MetaData
#> 3            dif  http://gcmd.gsfc.nasa.gov/Aboutus/xml/dif/dif_v9.4.xsd  http://gcmd.gsfc.nasa.gov/Aboutus/xml/dif/
#> 4       iso19139                http://www.isotc211.org/2005/gmd/gmd.xsd            http://www.isotc211.org/2005/gmd
#> 5  iso19139.iodp                http://www.isotc211.org/2005/gmd/gmd.xsd            http://www.isotc211.org/2005/gmd
#> 6      datacite3   http://schema.datacite.org/meta/kernel-3/metadata.xsd         http://datacite.org/schema/kernel-3
#> 7      datacite4 http://schema.datacite.org/meta/kernel-4.1/metadata.xsd         http://datacite.org/schema/kernel-4
```

### List identifiers


```r
pg_list_identifiers(from = Sys.Date() - 2, until = Sys.Date())
```

### List sets


```r
pg_list_sets()
```

```
#> # A tibble: 282 x 2
#>    setSpec   setName                                        
#>    <chr>     <chr>                                          
#>  1 ACD       PANGAEA set / keyword 'ACD' (2 data sets)      
#>  2 ASPS      PANGAEA set / keyword 'ASPS' (59 data sets)    
#>  3 AWIXRFraw PANGAEA set / keyword 'AWIXRFraw' (1 data sets)
#>  4 BAH1960   PANGAEA set / keyword 'BAH1960' (2 data sets)  
#>  5 BAH1961   PANGAEA set / keyword 'BAH1961' (2 data sets)  
#>  6 BAH1962   PANGAEA set / keyword 'BAH1962' (7 data sets)  
#>  7 BAH1963   PANGAEA set / keyword 'BAH1963' (7 data sets)  
#>  8 BAH1964   PANGAEA set / keyword 'BAH1964' (7 data sets)  
#>  9 BAH1965   PANGAEA set / keyword 'BAH1965' (7 data sets)  
#> 10 BAH1966   PANGAEA set / keyword 'BAH1966' (6 data sets)  
#> # … with 272 more rows
```

### List records


```r
pg_list_records(from = Sys.Date() - 1, until = Sys.Date())
```

### Get a record


```r
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")
```

```
#> $`oai:pangaea.de:doi:10.1594/PANGAEA.788382`
#> $`oai:pangaea.de:doi:10.1594/PANGAEA.788382`$header
#> # A tibble: 1 x 3
#>   identifier                                datestamp            setSpec                                           
#>   <chr>                                     <chr>                <chr>                                             
#> 1 oai:pangaea.de:doi:10.1594/PANGAEA.788382 2020-01-18T03:11:42Z citable;supplement;topicChemistry;topicLithosphere
#> 
#> $`oai:pangaea.de:doi:10.1594/PANGAEA.788382`$metadata
#> # A tibble: 1 x 13
#>   title            creator     source                publisher date   type   format  identifier       description               language rights       coverage               subject 
#>   <chr>            <chr>       <chr>                 <chr>     <chr>  <chr>  <chr>   <chr>            <chr>                     <chr>    <chr>        <chr>                  <chr>   
#> 1 Trace metals in… Demina, Ly… P.P. Shirshov Instit… PANGAEA   2012-… Datas… applic… https://doi.pan… Bioaccumulation of trace… en       CC-BY-3.0: … MEDIAN LATITUDE: 29.1… Archive…
```
pangaear
========

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![cran checks](https://cranchecks.info/badges/worst/pangaear)](https://cranchecks.info/pkgs/pangaear)
[![R-check](https://github.com/ropensci/pangaear/workflows/R-check/badge.svg)](https://github.com/ropensci/pangaear/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/pangaear/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/pangaear)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/pangaear)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/pangaear)](https://cran.r-project.org/package=pangaear)

`pangaear` is a data retrieval interface for the World Data Center PANGAEA (https://www.pangaea.de/). PANGAEA archieves published Earth & Environmental Science data under the following subjects: agriculture, atmosphere, biological classification, biosphere, chemistry, cryosphere, ecology, fisheries, geophysics, human dimensions, lakes & rives, land surface, lithosphere, oceans, and paleontology.

This package offers tools to interact with the PANGAEA Database, including functions for searching for data, fetching datasets by dataset ID, and working with the PANGAEA OAI-PMH service.

## Info

* Pangaea [website](https://www.pangaea.de/).
* Pangaea [OAI-PMH docs](https://wiki.pangaea.de/wiki/OAI-PMH).
* [OAI-PMH Spec](http://www.openarchives.org/OAI/openarchivesprotocol.html)

## Package API

```{r echo=FALSE, comment=NA, results='asis'}
cat(paste(" -", paste(getNamespaceExports("pangaear"), collapse = "\n - ")))
```

## Installation

Stable version

```{r eval=FALSE}
install.packages("pangaear")
```

Dev version

```{r eval=FALSE}
remotes::install_github('ropensci/pangaear')
```

```{r}
library('pangaear')
```

## Search for data

This is a thin wrapper around the GUI search interface on the page <https://www.pangaea.de/>. Everything you can do there, you can do here.

```{r}
pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
```

## Get data

```{r}
res <- pg_data(doi = '10.1594/PANGAEA.807580')
res[[1]]
```

Search for data then pass DOI to data function.

```{r}
res <- pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
pg_data(res$doi[3])[1:3]
```

## OAI-PMH metadata

```{r eval=FALSE}
# Identify the service
pg_identify()

# List metadata formats
pg_list_metadata_formats()

# List identifiers
pg_list_identifiers(from = Sys.Date() - 2, until = Sys.Date())

# List sets
pg_list_sets()

# List records
pg_list_records(from = Sys.Date() - 1, until = Sys.Date())

# Get a record
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")
```

## Contributors (reverse alphabetical)

* Naupaka Zimmerman
* Kara Woo
* Gavin Simpson
* Andrew MacDonald
* Scott Chamberlain

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pangaear/issues).
* License: MIT
* Get citation information for `pangaear` in R doing `citation(package = 'pangaear')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Introduction to pangaear}
%\VignetteEncoding{UTF-8}
-->

```{r echo=FALSE}
knitr::opts_chunk$set(
  fig.path = "img/",
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

Introduction to pangaear
========================

`pangaear` is a data retrieval interface for the World Data Center PANGAEA (https://www.pangaea.de/). PANGAEA archieves published Earth & Environmental Science data under the following subjects: agriculture, atmosphere, biological classification, biosphere, chemistry, cryosphere, ecology, fisheries, geophysics, human dimensions, lakes & rives, land surface, lithosphere, oceans, and paleontology.

## Installation

If you've not installed it yet, install from CRAN:

```{r eval=FALSE}
install.packages("pangaear")
```

Or the development version:

```{r eval=FALSE}
devtools::install_github("ropensci/pangaear")
```

## Load pangaear

```{r}
library("pangaear")
```

## Search for data

`pg_search` is a thin wrapper around the GUI search interface on the page <https://www.pangaea.de/>. Everything you can do there, you can do here.

For example, query for the term 'water', with a bounding box, and return only three results.

```{r}
pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
```

The resulting `data.frame` has details about different studies, and you can use the DOIs (Digital Object Identifiers) to get data and metadata for any studies you're interested in.

### Another search option

There's another search option with the `pg_search_es` function. It is an interface to the Pangaea Elasticsearch interface. This provides a very flexible interface for search Pangaea data - though it is different from what you're used to with the Pangaea website. 

```{r}
(res <- pg_search_es())
```

The returned data.frame has a lot of columns. You can limit columns returned with the `source` parameter.

There are attributes on the data.frame that give you the total number of results found as well as the max score found. 

```{r}
attributes(res)
attr(res, "total")
attr(res, "max_score")
```

To get to the DOIs for each study, use 

```{r}
gsub("https://doi.org/", "", res$`_source.URI`)
```


## Get data

The function `pg_data` fetches datasets for studies by their DOIs.

```{r}
res <- pg_data(doi = '10.1594/PANGAEA.807580')
res[[1]]
```

Search for data then pass one or more DOIs to the `pg_data` function.

```{r}
res <- pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
pg_data(res$doi[3])[1:3]
```


## OAI-PMH metadata

[OAI-PMH](https://wiki.pangaea.de/wiki/OAI-PMH) is a standard protocol for serving metadata around objects, in this case datasets. If you are already familiar with OAI-PMH you are in luck as you can can use what you know here. If not familiar, it's relatively straight-forward. 

Note that you can't get data through these functions, rather only metadata about datasets.

### Identify the service

```{r}
pg_identify()
```

### List metadata formats

```{r}
pg_list_metadata_formats()
```

### List identifiers

```{r eval=FALSE}
pg_list_identifiers(from = Sys.Date() - 2, until = Sys.Date())
```

### List sets

```{r}
pg_list_sets()
```

### List records

```{r eval=FALSE}
pg_list_records(from = Sys.Date() - 1, until = Sys.Date())
```

### Get a record

```{r}
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")
```
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Introduction to pangaear}
%\VignetteEncoding{UTF-8}
-->



Introduction to pangaear
========================

`pangaear` is a data retrieval interface for the World Data Center PANGAEA (https://www.pangaea.de/). PANGAEA archieves published Earth & Environmental Science data under the following subjects: agriculture, atmosphere, biological classification, biosphere, chemistry, cryosphere, ecology, fisheries, geophysics, human dimensions, lakes & rives, land surface, lithosphere, oceans, and paleontology.

## Installation

If you've not installed it yet, install from CRAN:


```r
install.packages("pangaear")
```

Or the development version:


```r
devtools::install_github("ropensci/pangaear")
```

## Load pangaear


```r
library("pangaear")
```

## Search for data

`pg_search` is a thin wrapper around the GUI search interface on the page <https://www.pangaea.de/>. Everything you can do there, you can do here.

For example, query for the term 'water', with a bounding box, and return only three results.


```r
pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
```

```
#> # A tibble: 3 x 6
#>   score doi            size size_measure citation                                                 supplement_to                                                                      
#>   <dbl> <chr>         <dbl> <chr>        <chr>                                                    <chr>                                                                              
#> 1 13.0  10.1594/PANG…     4 datasets     Krylova, EM; Sahling, H; Janssen, R (2010): A new genus… Krylova, EM; Sahling, H; Janssen, R (2010): Abyssogena: a new genus of the family …
#> 2 12.8  10.1594/PANG…     2 datasets     Simonyan, AV; Dultz, S; Behrens, H (2012): Diffusion tr… Simonyan, AV; Dultz, S; Behrens, H (2012): Diffusive transport of water in porous …
#> 3  8.78 10.1594/PANG…  1148 data points  WOCE Surface Velocity Program, SVP (2006): Water temper… <NA>
```

The resulting `data.frame` has details about different studies, and you can use the DOIs (Digital Object Identifiers) to get data and metadata for any studies you're interested in.

### Another search option

There's another search option with the `pg_search_es` function. It is an interface to the Pangaea Elasticsearch interface. This provides a very flexible interface for search Pangaea data - though it is different from what you're used to with the Pangaea website. 


```r
(res <- pg_search_es())
```

```
#> # A tibble: 10 x 46
#>    `_index` `_type` `_id` `_score` `_source.intern… `_source.parent… `_source.minEle… `_source.sf-aut… `_source.parent… `_source.techKe… `_source.geocod… `_source.sp-log…
#>  * <chr>    <chr>   <chr>    <dbl> <chr>            <chr>                       <dbl> <chr>                       <int> <list>           <list>                      <int>
#>  1 pangaea… panmd   9018…        1 2020-01-17T15:1… https://doi.org…             -0.4 Anhaus Philipp#…           901247 <chr [654]>      <chr [5]>                       1
#>  2 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              2   Anhaus Philipp#…           901247 <chr [652]>      <chr [5]>                       1
#>  3 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              2   Anhaus Philipp#…           901247 <chr [652]>      <chr [5]>                       1
#>  4 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              2   Anhaus Philipp#…           901247 <chr [652]>      <chr [5]>                       1
#>  5 pangaea… panmd   9017…        1 2020-01-17T15:1… <NA>                          2   Vuilleumier Lau…               NA <chr [51]>       <chr [6]>                       3
#>  6 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.pan…              0.3 Peeken Ilka#Mur…           901742 <chr [101]>      <chr [6]>                       3
#>  7 pangaea… panmd   9017…        1 2020-01-17T15:1… <NA>                         NA   Schreuder Laura…               NA <chr [37]>       <chr [3]>                       1
#>  8 pangaea… panmd   9017…        1 2020-01-17T15:1… https://doi.org…              0   Schreuder Laura…           901739 <chr [39]>       <chr [7]>                       1
#>  9 pangaea… panmd   9016…        1 2020-01-17T15:1… <NA>                         10   Augustine John$…               NA <chr [25]>       <chr [6]>                       3
#> 10 pangaea… panmd   9016…        1 2020-01-17T15:1… <NA>                          2   Augustine John$…               NA <chr [29]>       <chr [6]>                       3
#> # … with 34 more variables: `_source.agg-campaign` <list>, `_source.agg-author` <list>, `_source.eastBoundLongitude` <dbl>, `_source.URI` <chr>, `_source.agg-pubYear` <int>,
#> #   `_source.minDateTime` <chr>, `_source.agg-geometry` <chr>, `_source.xml-thumb` <chr>, `_source.xml` <chr>, `_source.sf-idDataSet` <int>, `_source.elevationGeocode` <chr>,
#> #   `_source.agg-method` <list>, `_source.maxDateTime` <chr>, `_source.xml-sitemap` <chr>, `_source.westBoundLongitude` <dbl>, `_source.northBoundLatitude` <dbl>,
#> #   `_source.sp-dataStatus` <int>, `_source.nDataPoints` <int>, `_source.sp-hidden` <lgl>, `_source.agg-location` <list>, `_source.internal-source` <chr>,
#> #   `_source.agg-basis` <chr>, `_source.southBoundLatitude` <dbl>, `_source.boost` <dbl>, `_source.oaiSet` <list>, `_source.maxElevation` <dbl>, `_source.agg-mainTopic` <list>,
#> #   `_source.agg-topic` <list>, `_source.agg-project` <list>, `_source.meanPosition.lat` <dbl>, `_source.meanPosition.lon` <dbl>, `_source.geoCoverage.type` <chr>,
#> #   `_source.geoCoverage.coordinates` <list>, `_source.geoCoverage.geometries` <list>
```

The returned data.frame has a lot of columns. You can limit columns returned with the `source` parameter.

There are attributes on the data.frame that give you the total number of results found as well as the max score found. 


```r
attributes(res)
```

```
#> $names
#>  [1] "_index"                          "_type"                           "_id"                             "_score"                          "_source.internal-datestamp"     
#>  [6] "_source.parentURI"               "_source.minElevation"            "_source.sf-authortitle"          "_source.parentIdDataSet"         "_source.techKeyword"            
#> [11] "_source.geocodes"                "_source.sp-loginOption"          "_source.agg-campaign"            "_source.agg-author"              "_source.eastBoundLongitude"     
#> [16] "_source.URI"                     "_source.agg-pubYear"             "_source.minDateTime"             "_source.agg-geometry"            "_source.xml-thumb"              
#> [21] "_source.xml"                     "_source.sf-idDataSet"            "_source.elevationGeocode"        "_source.agg-method"              "_source.maxDateTime"            
#> [26] "_source.xml-sitemap"             "_source.westBoundLongitude"      "_source.northBoundLatitude"      "_source.sp-dataStatus"           "_source.nDataPoints"            
#> [31] "_source.sp-hidden"               "_source.agg-location"            "_source.internal-source"         "_source.agg-basis"               "_source.southBoundLatitude"     
#> [36] "_source.boost"                   "_source.oaiSet"                  "_source.maxElevation"            "_source.agg-mainTopic"           "_source.agg-topic"              
#> [41] "_source.agg-project"             "_source.meanPosition.lat"        "_source.meanPosition.lon"        "_source.geoCoverage.type"        "_source.geoCoverage.coordinates"
#> [46] "_source.geoCoverage.geometries" 
#> 
#> $row.names
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $class
#> [1] "tbl_df"     "tbl"        "data.frame"
#> 
#> $total
#> [1] 390620
#> 
#> $max_score
#> [1] 1
```

```r
attr(res, "total")
```

```
#> [1] 390620
```

```r
attr(res, "max_score")
```

```
#> [1] 1
```

To get to the DOIs for each study, use 


```r
gsub("https://doi.org/", "", res$`_source.URI`)
```

```
#>  [1] "10.1594/PANGAEA.901810"                        "10.1594/PANGAEA.901736"                        "10.1594/PANGAEA.901733"                       
#>  [4] "10.1594/PANGAEA.901732"                        "https://doi.pangaea.de/10.1594/PANGAEA.901710" "https://doi.pangaea.de/10.1594/PANGAEA.901709"
#>  [7] "10.1594/PANGAEA.901739"                        "10.1594/PANGAEA.901738"                        "https://doi.pangaea.de/10.1594/PANGAEA.901697"
#> [10] "https://doi.pangaea.de/10.1594/PANGAEA.901695"
```


## Get data

The function `pg_data` fetches datasets for studies by their DOIs.


```r
res <- pg_data(doi = '10.1594/PANGAEA.807580')
res[[1]]
```

```
#> <Pangaea data> 10.1594/PANGAEA.807580
#>   parent doi: 10.1594/PANGAEA.807580
#>   url:        https://doi.org/10.1594/PANGAEA.807580
#>   citation:   Schiebel, Ralf; Waniek, Joanna J; Bork, Matthias; Hemleben, Christoph (2001): Physical oceanography during METEOR cruise M36/6. PANGAEA, https://doi.org/10.1594/PANGAEA.807580, In supplement to: Schiebel, R et al. (2001): Planktic foraminiferal production stimulated by chlorophyll redistribution and entrainment of nutrients. Deep Sea Research Part I: Oceanographic Research Papers, 48(3), 721-740, https://doi.org/10.1016/S0967-0637(00)00065-0
#>   path:       /Users/sckott/Library/Caches/R/pangaear/10_1594_PANGAEA_807580.txt
#>   data:
#> # A tibble: 32,179 x 13
#>    Event       `Date/Time`    Latitude Longitude `Elevation [m]` `Depth water [m… `Press [dbar]` `Temp [°C]`   Sal `Tpot [°C]` `Sigma-theta [kg/m… `Sigma in situ [kg… `Cond [mS/cm]`
#>    <chr>       <chr>             <dbl>     <dbl>           <int>            <dbl>          <int>       <dbl> <dbl>       <dbl>               <dbl>               <dbl>          <dbl>
#>  1 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             0                 0        15.7  35.7        15.7                26.4                26.4           44.4
#>  2 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             0.99              1        15.7  35.7        15.7                26.4                26.4           44.4
#>  3 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             1.98              2        15.7  35.7        15.7                26.4                26.4           44.4
#>  4 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             2.97              3        15.7  35.7        15.7                26.4                26.4           44.4
#>  5 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             3.96              4        15.7  35.7        15.7                26.4                26.4           44.4
#>  6 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             4.96              5        15.7  35.7        15.7                26.4                26.4           44.4
#>  7 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             5.95              6        15.7  35.7        15.7                26.4                26.4           44.4
#>  8 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             6.94              7        15.7  35.7        15.7                26.4                26.4           44.4
#>  9 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             7.93              8        15.7  35.7        15.7                26.4                26.4           44.4
#> 10 M36/6-CTD-… 1996-10-14T12…     49.0     -16.5           -4802             8.92              9        15.7  35.7        15.7                26.4                26.4           44.4
#> # … with 32,169 more rows
```

Search for data then pass one or more DOIs to the `pg_data` function.


```r
res <- pg_search(query = 'water', bbox = c(-124.2, 41.8, -116.8, 46.1), count = 3)
pg_data(res$doi[3])[1:3]
```

```
#> [[1]]
#> <Pangaea data> 10.1594/PANGAEA.405695
#>   parent doi: 10.1594/PANGAEA.405695
#>   url:        https://doi.org/10.1594/PANGAEA.405695
#>   citation:   WOCE Surface Velocity Program, SVP (2006): Water temperature and current velocity from surface drifter SVP_9524470. PANGAEA, https://doi.org/10.1594/PANGAEA.405695
#>   path:       /Users/sckott/Library/Caches/R/pangaear/10_1594_PANGAEA_405695.txt
#>   data:
#> # A tibble: 192 x 10
#>    `Date/Time`      Latitude Longitude `Depth water [m]` `Temp [°C]` `UC [cm/s]` `VC [cm/s]` `Lat e` `Lon e`  Code
#>    <chr>               <dbl>     <dbl>             <int>       <dbl>       <dbl>       <dbl>   <dbl>   <dbl> <int>
#>  1 1995-12-21T18:00     42.9     -125.                 0        12.4       NA          NA     0.0001  0.0001     1
#>  2 1995-12-22T00:00     42.9     -125.                 0        12.4        4.24      -35.2   0       0          1
#>  3 1995-12-22T06:00     42.8     -125.                 0        12.4       -7.8       -38.7   0.0002  0.0001     1
#>  4 1995-12-22T12:00     42.7     -125.                 0        12.4      -12.7       -23.1   0.0001  0.0001     1
#>  5 1995-12-22T18:00     42.7     -125.                 0        12.4      -15.1       -13.7   0       0          1
#>  6 1995-12-23T00:00     42.7     -125.                 0        12.3      -24.1        -9.15  0.0001  0.0001     1
#>  7 1995-12-23T06:00     42.7     -125.                 0        12.3      -38.4         0.65  0.0001  0.0001     1
#>  8 1995-12-23T12:00     42.7     -125.                 0        12.3      -37.2        15.8   0       0          1
#>  9 1995-12-23T18:00     42.7     -125.                 0        12.3      -25.4        29.6   0.0001  0.0001     1
#> 10 1995-12-24T00:00     42.8     -125.                 0        12.4      -18.5        35.1   0.0002  0.0002     1
#> # … with 182 more rows
#> 
#> [[2]]
#> NULL
#> 
#> [[3]]
#> NULL
```


## OAI-PMH metadata

[OAI-PMH](https://wiki.pangaea.de/wiki/OAI-PMH) is a standard protocol for serving metadata around objects, in this case datasets. If you are already familiar with OAI-PMH you are in luck as you can can use what you know here. If not familiar, it's relatively straight-forward. 

Note that you can't get data through these functions, rather only metadata about datasets.

### Identify the service


```r
pg_identify()
```

```
#> <Pangaea>
#>   repositoryName: PANGAEA - Data Publisher for Earth & Environmental Science
#>   baseURL: https://ws.pangaea.de/oai/provider
#>   protocolVersion: 2.0
#>   adminEmail: tech@pangaea.de
#>   adminEmail: tech@pangaea.de
#>   earliestDatestamp: 2015-01-01T00:00:00Z
#>   deletedRecord: transient
#>   granularity: YYYY-MM-DDThh:mm:ssZ
#>   compression: gzip
#>   description: oaipangaea.de:oai:pangaea.de:doi:10.1594/PANGAEA.999999
```

### List metadata formats


```r
pg_list_metadata_formats()
```

```
#>   metadataPrefix                                                  schema                           metadataNamespace
#> 1         oai_dc          http://www.openarchives.org/OAI/2.0/oai_dc.xsd http://www.openarchives.org/OAI/2.0/oai_dc/
#> 2         pan_md       http://ws.pangaea.de/schemas/pangaea/MetaData.xsd              http://www.pangaea.de/MetaData
#> 3            dif  http://gcmd.gsfc.nasa.gov/Aboutus/xml/dif/dif_v9.4.xsd  http://gcmd.gsfc.nasa.gov/Aboutus/xml/dif/
#> 4       iso19139                http://www.isotc211.org/2005/gmd/gmd.xsd            http://www.isotc211.org/2005/gmd
#> 5  iso19139.iodp                http://www.isotc211.org/2005/gmd/gmd.xsd            http://www.isotc211.org/2005/gmd
#> 6      datacite3   http://schema.datacite.org/meta/kernel-3/metadata.xsd         http://datacite.org/schema/kernel-3
#> 7      datacite4 http://schema.datacite.org/meta/kernel-4.1/metadata.xsd         http://datacite.org/schema/kernel-4
```

### List identifiers


```r
pg_list_identifiers(from = Sys.Date() - 2, until = Sys.Date())
```

### List sets


```r
pg_list_sets()
```

```
#> # A tibble: 282 x 2
#>    setSpec   setName                                        
#>    <chr>     <chr>                                          
#>  1 ACD       PANGAEA set / keyword 'ACD' (2 data sets)      
#>  2 ASPS      PANGAEA set / keyword 'ASPS' (59 data sets)    
#>  3 AWIXRFraw PANGAEA set / keyword 'AWIXRFraw' (1 data sets)
#>  4 BAH1960   PANGAEA set / keyword 'BAH1960' (2 data sets)  
#>  5 BAH1961   PANGAEA set / keyword 'BAH1961' (2 data sets)  
#>  6 BAH1962   PANGAEA set / keyword 'BAH1962' (7 data sets)  
#>  7 BAH1963   PANGAEA set / keyword 'BAH1963' (7 data sets)  
#>  8 BAH1964   PANGAEA set / keyword 'BAH1964' (7 data sets)  
#>  9 BAH1965   PANGAEA set / keyword 'BAH1965' (7 data sets)  
#> 10 BAH1966   PANGAEA set / keyword 'BAH1966' (6 data sets)  
#> # … with 272 more rows
```

### List records


```r
pg_list_records(from = Sys.Date() - 1, until = Sys.Date())
```

### Get a record


```r
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")
```

```
#> $`oai:pangaea.de:doi:10.1594/PANGAEA.788382`
#> $`oai:pangaea.de:doi:10.1594/PANGAEA.788382`$header
#> # A tibble: 1 x 3
#>   identifier                                datestamp            setSpec                                           
#>   <chr>                                     <chr>                <chr>                                             
#> 1 oai:pangaea.de:doi:10.1594/PANGAEA.788382 2020-01-18T03:11:42Z citable;supplement;topicChemistry;topicLithosphere
#> 
#> $`oai:pangaea.de:doi:10.1594/PANGAEA.788382`$metadata
#> # A tibble: 1 x 13
#>   title            creator     source                publisher date   type   format  identifier       description               language rights       coverage               subject 
#>   <chr>            <chr>       <chr>                 <chr>     <chr>  <chr>  <chr>   <chr>            <chr>                     <chr>    <chr>        <chr>                  <chr>   
#> 1 Trace metals in… Demina, Ly… P.P. Shirshov Instit… PANGAEA   2012-… Datas… applic… https://doi.pan… Bioaccumulation of trace… en       CC-BY-3.0: … MEDIAN LATITUDE: 29.1… Archive…
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caching.R
\name{pg_cache}
\alias{pg_cache}
\title{Caching}
\description{
Manage cached \code{pangaear} files with \pkg{hoardr}
}
\details{
The dafault cache directory is
\code{paste0(rappdirs::user_cache_dir(), "/R/pangaear")}, but you can set
your own path using \code{cache_path_set()}

\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}
\section{Useful user functions}{

\itemize{
\item \code{pg_cache$cache_path_get()} get cache path
\item \code{pg_cache$cache_path_set()} set cache path
\item \code{pg_cache$list()} returns a character vector of full path file names
\item \code{pg_cache$files()} returns file objects with metadata
\item \code{pg_cache$details()} returns files with details
\item \code{pg_cache$delete()} delete specific files
\item \code{pg_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
pg_cache

# list files in cache
pg_cache$list()

# delete certain database files
# pg_cache$delete("file path")
# pg_cache$list()

# delete all files in cache
# pg_cache$delete_all()
# pg_cache$list()

# set a different cache path from the default
# pg_cache$cache_path_set(full_path = "/Foo/Bar")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_identify.R
\name{pg_identify}
\alias{pg_identify}
\title{Identify information about the Pangaea repository}
\usage{
pg_identify(...)
}
\arguments{
\item{...}{Curl debugging options passed on to \code{\link[oai:id]{oai::id()}}}
}
\value{
list
}
\description{
Identify information about the Pangaea repository
}
\examples{
\dontrun{
pg_identify()
}
}
\references{
\href{https://www.openarchives.org/pmh/}{OAI-PMH documentation}
}
\seealso{
wraps \code{\link[oai:id]{oai::id()}}

Other oai methods: 
\code{\link{pg_get_record}()},
\code{\link{pg_list_identifiers}()},
\code{\link{pg_list_metadata_formats}()},
\code{\link{pg_list_records}()},
\code{\link{pg_list_sets}()}
}
\concept{oai methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_list_identifiers.R
\name{pg_list_identifiers}
\alias{pg_list_identifiers}
\title{List identifiers of the Pangaea repository}
\usage{
pg_list_identifiers(
  prefix = "oai_dc",
  from = NULL,
  until = NULL,
  set = NULL,
  token = NULL,
  as = "df",
  ...
)
}
\arguments{
\item{prefix}{A character string to specify the metadata format in OAI-PMH
requests issued to the repository. The default (\code{oai_dc}) corresponds
to the mandatory OAI unqualified Dublin Core metadata schema.}

\item{from}{Character string giving datestamp to be used as lower bound for
datestamp-based selective harvesting (i.e., only harvest records with
datestamps in the given range). Dates and times must be encoded using
ISO 8601. The trailing Z must be used when including time. OAI-PMH implies
UTC for data/time specifications.}

\item{until}{Character string giving a datestamp to be used as an upper
bound, for datestamp-based selective harvesting (i.e., only harvest records
with datestamps in the given range).}

\item{set}{A character string giving a set to be used for selective
harvesting (i.e., only harvest records in the given set).}

\item{token}{(character) a token previously provided by the server to
resume a request where it last left off. 50 is max number of records
returned. We will loop for you internally to get all the records you
asked for.}

\item{as}{(character) What to return. One of "df" (for data.frame; default),
"list", or "raw" (raw text)}

\item{...}{Curl debugging options passed on to \code{\link[oai:list_identifiers]{oai::list_identifiers()}}}
}
\value{
XML character string, data.frame, or list, depending on what
requested with the \code{as} parameter
}
\description{
List identifiers of the Pangaea repository
}
\examples{
\dontrun{
pg_list_identifiers(
  from = paste0(Sys.Date() - 4, "T00:00:00Z"),
  until = paste0(Sys.Date() - 3, "T18:00:00Z")
)
pg_list_identifiers(set="geocode1", from=Sys.Date()-1, until=Sys.Date())
pg_list_identifiers(prefix="iso19139", from=Sys.Date()-1, until=Sys.Date())
pg_list_identifiers(prefix="dif",
  from = paste0(Sys.Date() - 2, "T00:00:00Z"),
  until = paste0(Sys.Date() - 1, "T18:00:00Z")
)
}
}
\references{
\href{https://www.openarchives.org/pmh/}{OAI-PMH documentation}
}
\seealso{
wraps \code{\link[oai:list_identifiers]{oai::list_identifiers()}}

Other oai methods: 
\code{\link{pg_get_record}()},
\code{\link{pg_identify}()},
\code{\link{pg_list_metadata_formats}()},
\code{\link{pg_list_records}()},
\code{\link{pg_list_sets}()}
}
\concept{oai methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_list_metadata_formats.R
\name{pg_list_metadata_formats}
\alias{pg_list_metadata_formats}
\title{Get metadata formats from the Pangaea repository}
\usage{
pg_list_metadata_formats(...)
}
\arguments{
\item{...}{Curl debugging options passed on to \code{\link[oai:list_metadataformats]{oai::list_metadataformats()}}}
}
\value{
data.frame
}
\description{
Get metadata formats from the Pangaea repository
}
\examples{
\dontrun{
pg_list_metadata_formats()
}
}
\references{
\href{https://www.openarchives.org/pmh/}{OAI-PMH documentation}
}
\seealso{
wraps \code{\link[oai:list_metadataformats]{oai::list_metadataformats()}}

Other oai methods: 
\code{\link{pg_get_record}()},
\code{\link{pg_identify}()},
\code{\link{pg_list_identifiers}()},
\code{\link{pg_list_records}()},
\code{\link{pg_list_sets}()}
}
\concept{oai methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_cache.R
\name{pg_cache_list}
\alias{pg_cache_list}
\title{cache list}
\usage{
pg_cache_list(...)
}
\arguments{
\item{...}{ignored}
}
\description{
cache list
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_get_record.R
\name{pg_get_record}
\alias{pg_get_record}
\title{Get record from the Pangaea repository}
\usage{
pg_get_record(identifier, prefix = "oai_dc", as = "df", ...)
}
\arguments{
\item{identifier}{Dataset identifier. See Examples.}

\item{prefix}{A character string to specify the metadata format in OAI-PMH
requests issued to the repository. The default (\code{oai_dc}) corresponds
to the mandatory OAI unqualified Dublin Core metadata schema.}

\item{as}{(character) What to return. One of "df" (for data.frame; default),
"list", or "raw" (raw text)}

\item{...}{Curl debugging options passed on to \code{\link[oai:get_records]{oai::get_records()}}}
}
\value{
XML character string, data.frame, or list, depending on what
requested with the \code{as} parameter
}
\description{
Get record from the Pangaea repository
}
\examples{
\dontrun{
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.269656",
prefix="iso19139")
pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.269656",
prefix="dif")

# invalid record id
# pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.11111")
# pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.11111",
#   prefix="adfadf")
}
}
\references{
\href{https://www.openarchives.org/pmh/}{OAI-PMH documentation}
}
\seealso{
wraps \code{\link[oai:get_records]{oai::get_records()}}

Other oai methods: 
\code{\link{pg_identify}()},
\code{\link{pg_list_identifiers}()},
\code{\link{pg_list_metadata_formats}()},
\code{\link{pg_list_records}()},
\code{\link{pg_list_sets}()}
}
\concept{oai methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_search.R
\name{pg_search}
\alias{pg_search}
\title{Search the Pangaea database}
\usage{
pg_search(
  query,
  count = 10,
  offset = 0,
  topic = NULL,
  bbox = NULL,
  mindate = NULL,
  maxdate = NULL,
  ...
)
}
\arguments{
\item{query}{(character) Query terms. You can refine a search by prefixing
the term(s) with a category, one of citation, reference, parameter, event,
project, campaign, or basis. See examples.}

\item{count}{(integer) Number of items to return. Default: 10. Maximum: 500.
Use \code{offset} parameter to page through results - see examples}

\item{offset}{(integer) Record number to start at. Default: 0}

\item{topic}{(character) topic area: one of NULL (all areas), "Agriculture",
"Atomosphere", "Biological Classification", "Biospshere", "Chemistry",
"Cryosphere", "Ecology", "Fisheries", "Geophysics", "Human Dimensions",
"Lakes & Rivers", "Land Surface", "Lithosphere", "Oceans", "Paleontology"}

\item{bbox}{(numeric) A bounding box, of the form: minlon, minlat, maxlon,
maxlat}

\item{mindate, maxdate}{(character) Dates to search for, of the form
"2014-10-28"}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
tibble/data.frame with the structure:
\itemize{
\item score: match score, higher is a better match
\item doi: the DOI for the data package
\item size: size number
\item size_measure: size measure, one of "data points" or "datasets"
\item citation: citation for the data package
\item supplement_to: citation for what the data package is a supplement to
}
}
\description{
Search the Pangaea database
}
\details{
This is a thin wrapper around the GUI search interface on the page
\url{https://www.pangaea.de}. Everything you can do there, you can do here.
}
\examples{
\dontrun{
pg_search(query='water')
pg_search(query='water', count=2)
pg_search(query='water', count=20)
pg_search(query='water', mindate="2013-06-01", maxdate="2013-07-01")
pg_search(query='water', bbox=c(-124.2, 41.8, -116.8, 46.1))
pg_search(query='reference:Archer')
pg_search(query='parameter:"carbon dioxide"')
pg_search(query='event:M2-track')
pg_search(query='event:TT011_2-CTD31')
pg_search(query='project:Joint Global Ocean Flux Study')
pg_search(query='campaign:M2')
pg_search(query='basis:Meteor')

# paging with count and offset
# max is 500 records per request - if you need > 500, use offset and count
res1 <- pg_search(query = "florisphaera", count = 500, offset = 0)
res2 <- pg_search(query = "florisphaera", count = 500, offset = 500)
res3 <- pg_search(query = "florisphaera", count = 500, offset = 1000)
do.call("rbind.data.frame", list(res1, res2, res3))

# get attributes: maxScore, totalCount, and offset
res <- pg_search(query='water', bbox=c(-124.2, 41.8, -116.8, 46.1))
attributes(res)
attr(res, "maxScore")
attr(res, "totalCount")
attr(res, "offset")

# curl options
pg_search(query='citation:Archer', verbose = TRUE)
}
}
\seealso{
\code{\link[=pg_search_es]{pg_search_es()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_list_sets.R
\name{pg_list_sets}
\alias{pg_list_sets}
\title{List the set structure of the Pangaea repository}
\usage{
pg_list_sets(token = NULL, as = "df", ...)
}
\arguments{
\item{token}{(character) a token previously provided by the server to
resume a request where it last left off. 50 is max number of records
returned. We will loop for you internally to get all the records you
asked for.}

\item{as}{(character) What to return. One of "df" (for data.frame; default),
"list", or "raw" (raw text)}

\item{...}{Curl debugging options passed on to \code{\link[oai:list_sets]{oai::list_sets()}}}
}
\value{
XML character string, data.frame, or list, depending on what
requested with the \code{as} parameter
}
\description{
List the set structure of the Pangaea repository
}
\examples{
\dontrun{
pg_list_sets()
pg_list_sets(as = "list")
pg_list_sets(as = "raw")
}
}
\references{
\href{https://www.openarchives.org/pmh/}{OAI-PMH documentation}
}
\seealso{
wraps \code{\link[oai:list_sets]{oai::list_sets()}}

Other oai methods: 
\code{\link{pg_get_record}()},
\code{\link{pg_identify}()},
\code{\link{pg_list_identifiers}()},
\code{\link{pg_list_metadata_formats}()},
\code{\link{pg_list_records}()}
}
\concept{oai methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_search_es.R
\name{pg_search_es}
\alias{pg_search_es}
\title{Search the Pangaea database with Elasticsearch}
\usage{
pg_search_es(
  query = NULL,
  size = 10,
  from = NULL,
  source = NULL,
  df = NULL,
  analyzer = NULL,
  default_operator = NULL,
  explain = NULL,
  sort = NULL,
  track_scores = NULL,
  timeout = NULL,
  terminate_after = NULL,
  search_type = NULL,
  lowercase_expanded_terms = NULL,
  analyze_wildcard = NULL,
  version = FALSE,
  ...
)
}
\arguments{
\item{query}{(character) Query terms..}

\item{size}{(character) The number of hits to return. Pass in as a
character string to avoid problems with large number conversion to
scientific notation. Default: 10. The default maximum is 10,000 - however,
you can change this default maximum by changing the
\code{index.max_result_window} index level parameter.}

\item{from}{(character) The starting from index of the hits to return.
Pass in as a character string to avoid problems with large number
conversion to scientific notation. Default: 0}

\item{source}{(character) character vector of fields to return}

\item{df}{(character) The default field to use when no field prefix is
defined within the query.}

\item{analyzer}{(character) The analyzer name to be used when analyzing the
query string.}

\item{default_operator}{(character) The default operator to be used, can be
\code{AND} or \code{OR}. Default: \code{OR}}

\item{explain}{(logical) For each hit, contain an explanation of how
scoring of the hits was computed. Default: \code{FALSE}}

\item{sort}{(character) Sorting to perform. Can either be in the form of
fieldName, or \code{fieldName:asc}/\code{fieldName:desc}. The fieldName
can either be an actual field within the document, or the special
\verb{_score} name to indicate sorting based on scores. There can be several
sort parameters (order is important).}

\item{track_scores}{(logical) When sorting, set to \code{TRUE} in order to
still track scores and return them as part of each hit.}

\item{timeout}{(numeric) A search timeout, bounding the search request to
be executed within the specified time value and bail with the hits
accumulated up to that point when expired. Default: no timeout.}

\item{terminate_after}{(numeric) The maximum number of documents to collect
for each shard, upon reaching which the query execution will terminate
early. If set, the response will have a boolean field terminated_early to
indicate whether the query execution has actually terminated_early.
Default: no terminate_after}

\item{search_type}{(character) The type of the search operation to perform.
Can be \code{query_then_fetch} (default) or \code{dfs_query_then_fetch}.
Types \code{scan} and \code{count} are deprecated.}

\item{lowercase_expanded_terms}{(logical) Should terms be automatically
lowercased or not. Default: \code{TRUE}.}

\item{analyze_wildcard}{(logical) Should wildcard and prefix queries be
analyzed or not. Default: \code{FALSE}}

\item{version}{(logical) Print the document version with each document.}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
tibble/data.frame, empty if no results
}
\description{
Search the Pangaea database with Elasticsearch
}
\details{
An interface to Pangaea's Elasticsearch query interface.
You can also just use \href{https://github.com/ropensci/elastic}{elastic}
package to interact with it. The base URL is
https://ws.pangaea.de/es/pangaea/panmd/_search
}
\examples{
\dontrun{
(res <- pg_search_es())
attributes(res)
attr(res, "total")
attr(res, "max_score")

pg_search_es(query = 'water', source = c('parentURI', 'minElevation'))
pg_search_es(query = 'water', size = 3)
pg_search_es(query = 'water', size = 3, from = 10)

pg_search_es(query = 'water sky', default_operator = "OR")
pg_search_es(query = 'water sky', default_operator = "AND")

pg_search_es(query = 'water', sort = "minElevation")
pg_search_es(query = 'water', sort = "minElevation:desc")
}
}
\seealso{
\code{\link[=pg_search]{pg_search()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_data.R
\name{pg_data}
\alias{pg_data}
\title{Download data from Pangaea.}
\usage{
pg_data(doi, overwrite = TRUE, mssgs = TRUE, ...)
}
\arguments{
\item{doi}{DOI of Pangaeae single dataset, or of a collection of datasets.
Expects either just a DOI of the form \code{10.1594/PANGAEA.746398}, or with
the URL part in front, like
\url{https://doi.pangaea.de/10.1594/PANGAEA.746398}}

\item{overwrite}{(logical) Ovewrite a file if one is found with the same name}

\item{mssgs}{(logical) print information messages. Default: \code{TRUE}}

\item{...}{Curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
One or more items of class pangaea, each with the doi, parent doi
(if many dois within a parent doi), url, citation, path, and data object.
Data object depends on what kind of file it is. For tabular data, we print
the first 10 columns or so; for a zip file we list the files in the zip
(but leave it up to the user to dig unzip and get files from the zip file);
for png files, we point the user to read the file in with \code{\link[png:readPNG]{png::readPNG()}}
}
\description{
Grabs data as a dataframe or list of dataframes from a Pangaea data
repository URI; see: \url{https://www.pangaea.de/}
}
\details{
Data files are stored in an operating system appropriate location.
Run \code{pg_cache$cache_path_get()} to get the storage location
on your machine. See \link{pg_cache} for more information, including how to
set a different base path for downloaded files.

Some files/datasets require the user to be logged in. For now we
just pass on these - that is, give back nothing other than metadata.
}
\examples{
\dontrun{
# a single file
(res <- pg_data(doi='10.1594/PANGAEA.807580'))
res[[1]]$doi
res[[1]]$citation
res[[1]]$data
res[[1]]$metadata

# another single file
pg_data(doi='10.1594/PANGAEA.807584')

# Many files
(res <- pg_data(doi='10.1594/PANGAEA.761032'))
res[[1]]
res[[2]]

# Manipulating the cache
## list files in the cache
pg_cache$list()

## clear all data
# pg_cache$delete_all()
pg_cache$list()

## clear a single dataset by DOI
pg_data(doi='10.1594/PANGAEA.812093')
pg_cache$list()
path <- grep("PANGAEA.812093", pg_cache$list(), value = TRUE)
pg_cache$delete(path)
pg_cache$list()

# search for datasets, then pass in DOIs
(searchres <- pg_search(query = 'birds', count = 20))
pg_data(searchres$doi[1])

# png file
pg_data(doi = "10.1594/PANGAEA.825428")

# zip file
pg_data(doi = "10.1594/PANGAEA.860500")

# login required
## we skip file download
pg_data("10.1594/PANGAEA.788547")
}
}
\references{
\url{https://www.pangaea.de}
}
\author{
Naupaka Zimmerman, Scott Chamberlain
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_cache.R
\name{pg_cache_clear}
\alias{pg_cache_clear}
\title{cache path clear}
\usage{
pg_cache_clear(...)
}
\arguments{
\item{...}{ignored}
}
\description{
cache path clear
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pangaear-package.R
\docType{package}
\name{pangaear-package}
\alias{pangaear-package}
\alias{pangaear}
\title{Client for the Pangaea Database}
\description{
\href{https://www.pangaea.de/}{Pangaea database}
}
\details{
Package includes tools to interact with the Pangaea Database,
including functions for searching for data, fetching datasets by
dataset ID, working with the Pangaea OAI-PMH service, and
Elasticsearch service.
}
\section{Getting data}{

The main workhorse function for getting data is \code{\link[=pg_data]{pg_data()}}.
One thing you may want to do is set a different path for caching
the data you download: see \link{pg_cache} for details
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pg_list_records.R
\name{pg_list_records}
\alias{pg_list_records}
\title{List records from Pangaea}
\usage{
pg_list_records(
  prefix = "oai_dc",
  from = NULL,
  until = NULL,
  set = NULL,
  token = NULL,
  as = "df",
  ...
)
}
\arguments{
\item{prefix}{A character string to specify the metadata format in OAI-PMH
requests issued to the repository. The default (\code{oai_dc}) corresponds
to the mandatory OAI unqualified Dublin Core metadata schema.}

\item{from}{Character string giving datestamp to be used as lower bound for
datestamp-based selective harvesting (i.e., only harvest records with
datestamps in the given range). Dates and times must be encoded using
ISO 8601. The trailing Z must be used when including time. OAI-PMH implies
UTC for data/time specifications.}

\item{until}{Character string giving a datestamp to be used as an upper
bound, for datestamp-based selective harvesting (i.e., only harvest records
with datestamps in the given range).}

\item{set}{A character string giving a set to be used for selective
harvesting (i.e., only harvest records in the given set).}

\item{token}{(character) a token previously provided by the server to
resume a request where it last left off. 50 is max number of records
returned. We will loop for you internally to get all the records you
asked for.}

\item{as}{(character) What to return. One of "df" (for data.frame; default),
"list", or "raw" (raw text)}

\item{...}{Curl debugging options passed on to \code{\link[oai:list_records]{oai::list_records()}}}
}
\value{
XML character string, data.frame, or list, depending on what
requested witht the \code{as} parameter
}
\description{
List records from Pangaea
}
\examples{
\dontrun{
pg_list_records(set='citable', from=Sys.Date()-1, until=Sys.Date())

# When no results found > "'noRecordsMatch'"
# pg_list_records(set='geomound', from='2015-01-01', until='2015-01-01')

pg_list_records(prefix="iso19139", set='citable', from=Sys.Date()-1,
  until=Sys.Date())

## FIXME - below are broken
# pg_list_records(prefix="dif", set='citable', from=Sys.Date()-4,
#   until=Sys.Date())
# pg_list_records(prefix="dif", set='project4094', from=Sys.Date()-4,
#   until=Sys.Date())
}
}
\references{
\href{https://www.openarchives.org/pmh/}{OAI-PMH documentation}
}
\seealso{
wraps \code{\link[oai:list_records]{oai::list_records()}}

Other oai methods: 
\code{\link{pg_get_record}()},
\code{\link{pg_identify}()},
\code{\link{pg_list_identifiers}()},
\code{\link{pg_list_metadata_formats}()},
\code{\link{pg_list_sets}()}
}
\concept{oai methods}
