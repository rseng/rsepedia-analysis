rbace
=====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rbace)](https://cranchecks.info/pkgs/rbace)
[![R-check](https://github.com/ropensci/rbace/workflows/R-check/badge.svg)](https://github.com/ropensci/rbace/actions?query=workflow%3AR-check)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rbace?color=C9A115)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rbace)](https://cran.r-project.org/package=rbace)


Client for interacting with the Bielefeld Academic Search Engine API.

Docs: https://docs.ropensci.org/rbace/

BASE API docs: https://www.base-search.net/about/download/base_interface.pdf

Access: The BASE API is IP address AND user-agent (see note below) restricted. The user agent is set correctly if you use this package, but you still need to get your IP address(es) white-listed by BASE. Request access at: https://www.base-search.net/about/en/contact.php - Note: the BASE website has a search portal you can use from anywhere; it's just the API that is IP and user-agent restricted.

Terminology:

- an IP address is the numeric label identifying a computer or server. the IP address for a computer can change, e.g., if you connect to a VPN
- a user-agent is a string of text that identifies the software requesting data from a server (in this case BASE's API).

Data from BASE (Bielefeld Academic Search Engine) https://www.base-search.net

[<img src="man/figures/BASE_search_engine_logo.svg.png" width="300">](https://www.base-search.net)

## Install


```r
install.packages("rbace")
```

or the dev version


```r
remotes::install_github("ropensci/rbace")
# OR the below should install the same thing
install.packages("rbace", repos = "https://dev.ropensci.org")
```


```r
library("rbace")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rbace/issues).
* License: MIT
* Get citation information for `rbace` in R doing `citation(package = 'rbace')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
rbace 0.2.2
===========

### MINOR IMPROVEMENTS

* improved some tests


rbace 0.2.0
===========

### MINOR IMPROVEMENTS

First version to CRAN.

* replace `dplyr` with `data.table::rbindlist` (#6)
* two new functions: `bs_profile()`, `bs_repositories()` (#9)
* add faceting to `bs_search()` (#10)
* `bs_search()` now uses RETRY/GET http requests (#14)
* `bs_search()` gains new parameter `filter` that's passed to `fq` param for the Solr server (#15)
* `filter` parameter in `bs_search()` can be `character` and `AsIs` in case you need to prevent html escaping (#16)
## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

---

This version improves some tests.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rbace/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rbace.git`
* Make sure to track progress upstream (i.e., on our version of `rbace` at `ropensci/rbace`) by doing `git remote add upstream https://github.com/ropensci/rbace.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rbace`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
rbace
=====

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/rbace)](https://cranchecks.info/pkgs/rbace)
[![R-check](https://github.com/ropensci/rbace/workflows/R-check/badge.svg)](https://github.com/ropensci/rbace/actions?query=workflow%3AR-check)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rbace?color=C9A115)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rbace)](https://cran.r-project.org/package=rbace)


Client for interacting with the Bielefeld Academic Search Engine API.

Docs: https://docs.ropensci.org/rbace/

BASE API docs: https://www.base-search.net/about/download/base_interface.pdf

Access: The BASE API is IP address AND user-agent (see note below) restricted. The user agent is set correctly if you use this package, but you still need to get your IP address(es) white-listed by BASE. Request access at: https://www.base-search.net/about/en/contact.php - Note: the BASE website has a search portal you can use from anywhere; it's just the API that is IP and user-agent restricted.

Terminology:

- an IP address is the numeric label identifying a computer or server. the IP address for a computer can change, e.g., if you connect to a VPN
- a user-agent is a string of text that identifies the software requesting data from a server (in this case BASE's API).

Data from BASE (Bielefeld Academic Search Engine) https://www.base-search.net

[<img src="man/figures/BASE_search_engine_logo.svg.png" width="300">](https://www.base-search.net)

## Install

```{r eval=FALSE}
install.packages("rbace")
```

or the dev version

```{r eval=FALSE}
remotes::install_github("ropensci/rbace")
# OR the below should install the same thing
install.packages("rbace", repos = "https://dev.ropensci.org")
```

```{r}
library("rbace")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rbace/issues).
* License: MIT
* Get citation information for `rbace` in R doing `citation(package = 'rbace')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "rbace"
author: "Scott Chamberlain"
date: "2020-10-12"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to rbace}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Client for interacting with the Bielefeld Academic Search Engine API.

Docs: https://docs.ropensci.org/rbace/

## Install


```r
install.packages("rbace", repos = "https://dev.ropensci.org")
```


```r
library("rbace")
```

## Get the profile for a repository



```r
bs_profile(target = "ftjhin")
#> # A tibble: 9 x 2
#>   name               value                                           
#>   <chr>              <chr>                                           
#> 1 activation_date    2019-12-05                                      
#> 2 country            de                                              
#> 3 name               HiN - Alexander von Humboldt im Netz (E-Journal)
#> 4 name_en            HiN - Alexander von Humboldt im Netz (E-Journal)
#> 5 num_non_oa_records 0                                               
#> 6 num_oa_cc_records  279                                             
#> 7 num_oa_pd_records  0                                               
#> 8 num_oa_records     279                                             
#> 9 num_records        279
```

## List repositories for a collection



```r
bs_repositories(coll = "ceu")
#> # A tibble: 3,056 x 2
#>    name                                                            internal_name
#>    <chr>                                                           <chr>        
#>  1 Gornye nauki i tekhnologii (E-Journal)                          ftjmst       
#>  2 Hanser (Carl Hanser Verlag - via Crossref)                      crhanserverl…
#>  3 Università della Calabria: Archivio Institutionale delle Tesi … ftunivcalabr…
#>  4 hogrefe (via Crossref)                                          crhogrefe    
#>  5 Iberoamericana Vervuert (via Crossref)                          crvervuert   
#>  6 Manchester University Press (via Crossref)                      crmanchestupr
#>  7 Bezmialem Vakıf Üniversitesi Kurumsal Akademik Arşiv            ftbezmialem  
#>  8 ForAP - Forschungsergebnisse von Absolventen und Promovierende… ftjforap     
#>  9 Movement and Nutrition in Health and Disease (E-Journal)        ftjmnhd      
#> 10 Radcliffe Group (via Crossref)                                  crradcliffe  
#> # … with 3,046 more rows
```

## Search

perform a search


```r
(res <- bs_search(coll = 'it', query = 'dccreator:manghi', boost = TRUE))
#> $docs
#> # A tibble: 10 x 32
#>    dchdate dcdocid dccontinent dccountry dccollection dcprovider dctitle
#>    <chr>   <chr>   <chr>       <chr>     <chr>        <chr>      <chr>  
#>  1 2017-0… 90f58a… ceu         it        ftpuma       PUMAlab (… Sfide …
#>  2 2017-0… 2c669c… ceu         it        ftpuma       PUMAlab (… OpenAI…
#>  3 2017-0… af1126… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  4 2017-0… fcaa8f… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  5 2017-0… b4843e… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  6 2020-0… 338c5e… ceu         it        ftunivmodena Archivio … Multi-…
#>  7 2017-0… 9be017… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  8 2017-0… 412c04… ceu         it        ftpuma       PUMAlab (… EFG191…
#>  9 2017-0… 967c7f… ceu         it        ftpuma       PUMAlab (… OpenAI…
#> 10 2017-0… 3d820d… ceu         it        ftpuma       PUMAlab (… DRIVER…
#> # … with 25 more variables: dccreator <chr>, dcperson <chr>, dcsubject <chr>,
#> #   dcdescription <chr>, dcpublisher <chr>, dcdate <chr>, dcyear <chr>,
#> #   dctype <chr>, dctypenorm <chr>, dcformat <chr>, dccontenttype <chr>,
#> #   dcidentifier <chr>, dclink <chr>, dcsource <chr>, dclanguage <chr>,
#> #   dcrelation <chr>, dcrights <chr>, dcoa <chr>, dclang <chr>,
#> #   dcautoclasscode <chr>, dcdeweyfull <chr>, dcdeweyhuns <chr>,
#> #   dcdeweytens <chr>, dcdeweyones <chr>, dcdoi <chr>
#> 
#> $facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "87"
#> attr(,"q")
#> [1] "creator:manghi"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"fq")
#> [1] " country:it -collection:(ftjmethode OR ftpenamultimedia)"
#> attr(,"bq")
#> [1] "oa:1^2"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 4871
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "6.5124907"
```

get the search metadata


```r
bs_meta(res)
#> $query
#> # A tibble: 1 x 4
#>   q         fl                                  fq                         start
#>   <chr>     <chr>                               <chr>                      <dbl>
#> 1 creator:… dccollection,dccontenttype,dcconti… " country:it -collection:…     0
#> 
#> $response
#> # A tibble: 1 x 2
#>   status num_found
#>    <dbl>     <dbl>
#> 1      0      4871
```

repository "ftubbiepub" containing the terms "lossau" and "summann" (search in the whole document)


```r
res <- bs_search(target = 'ftubbiepub', query = 'lossau summann', hits = 3)
res$docs$dcsubject
#> [1] "Elektronische Bibliothek; Suchmaschine; Information Retrieval; Bielefeld Academic Search Engine; ddc:020"                                                                                                                                              
#> [2] "Informatik; Massenmedien; Book; Bielefeld Academic Search Engine; Library; Documentation; Communication; Informatics; Libro; Biblioteca; Kommunikation; Documentazione; Comunicazione; Dokumentation; Informatica; Buch- und Bibliothekswesen; ddc:004"
res$docs$dccreator
#> [1] "Summann, Friedrich; Lossau, Norbert" "Summann, Friedrich; Lossau, Norbert"
res$docs$dcidentifier
#> [1] "https://nbn-resolving.org/urn:nbn:de:0070-pub-16809824; https://pub.uni-bielefeld.de/record/1680982; https://pub.uni-bielefeld.de/download/1680982/2315314"
#> [2] "https://nbn-resolving.org/urn:nbn:de:0070-pub-25166732; https://pub.uni-bielefeld.de/record/2516673; https://pub.uni-bielefeld.de/download/2516673/2516679"
attributes(res)
#> $names
#> [1] "docs"   "facets"
#> 
#> $status
#> [1] 0
#> 
#> $QTime
#> [1] "9"
#> 
#> $q
#> [1] "lossau summann"
#> 
#> $fl
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> 
#> $fq
#> [1] "collection:ftubbiepub"
#> 
#> $rows
#> [1] "3"
#> 
#> $name
#> [1] "response"
#> 
#> $numFound
#> [1] 2
#> 
#> $start
#> [1] 0
#> 
#> $maxScore
#> [1] "2.395762"
bs_meta(res)
#> $query
#> # A tibble: 1 x 4
#>   q           fl                                             fq            start
#>   <chr>       <chr>                                          <chr>         <dbl>
#> 1 lossau sum… dccollection,dccontenttype,dccontinent,dccoun… collection:f…     0
#> 
#> $response
#> # A tibble: 1 x 2
#>   status num_found
#>    <dbl>     <dbl>
#> 1      0         2
```

Italian repositories containing the term "manghi'
in the "dccreator" field (author).  The flag "boost" pushes open
access documents upwards in the result list


```r
(res <- bs_search(coll = 'it', query = 'dccreator:manghi', boost_oa = TRUE))
#> $docs
#> # A tibble: 10 x 31
#>    dchdate dcdocid dccontinent dccountry dccollection dcprovider dctitle
#>    <chr>   <chr>   <chr>       <chr>     <chr>        <chr>      <chr>  
#>  1 2017-0… 90f58a… ceu         it        ftpuma       PUMAlab (… Sfide …
#>  2 2017-0… 2c669c… ceu         it        ftpuma       PUMAlab (… OpenAI…
#>  3 2017-0… 9be017… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  4 2017-0… 412c04… ceu         it        ftpuma       PUMAlab (… EFG191…
#>  5 2017-0… 967c7f… ceu         it        ftpuma       PUMAlab (… OpenAI…
#>  6 2017-0… 3d820d… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  7 2017-0… 242b25… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  8 2017-0… af1126… ceu         it        ftpuma       PUMAlab (… DRIVER…
#>  9 2017-0… fcaa8f… ceu         it        ftpuma       PUMAlab (… DRIVER…
#> 10 2017-0… b4843e… ceu         it        ftpuma       PUMAlab (… DRIVER…
#> # … with 24 more variables: dccreator <chr>, dcperson <chr>, dcsubject <chr>,
#> #   dcdescription <chr>, dcpublisher <chr>, dcdate <chr>, dcyear <chr>,
#> #   dctype <chr>, dctypenorm <chr>, dcformat <chr>, dccontenttype <chr>,
#> #   dcidentifier <chr>, dclink <chr>, dcsource <chr>, dclanguage <chr>,
#> #   dcrelation <chr>, dcrights <chr>, dcoa <chr>, dclang <chr>,
#> #   dcautoclasscode <chr>, dcdeweyfull <chr>, dcdeweyhuns <chr>,
#> #   dcdeweytens <chr>, dcdeweyones <chr>
#> 
#> $facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "10"
#> attr(,"q")
#> [1] "creator:manghi"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"fq")
#> [1] " country:it -collection:(ftjmethode OR ftpenamultimedia)"
#> attr(,"bq")
#> [1] "oa:1^2"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 4871
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "6.5124907"
```

terms "schmidt" in dccreator field (author) and "biology" in dctitle.
The response starts after record 5 (offset=5) and contains max.
5 hits with the fields dctitle, dccreator and dcyear


```r
(res <- bs_search(query = 'dccreator:schmidt dctitle:biology',
  hits = 5, offset = 5, fields = c('dctitle', 'dccreator', 'dcyear')))
#> $docs
#> # A tibble: 5 x 3
#>   dctitle                            dccreator               dcyear
#>   <chr>                              <chr>                   <chr> 
#> 1 Biology: A degenerative affliction Schmidt, Charles        2016  
#> 2 Biology of Aging                   Schmidt, Barbara        2014  
#> 3 Root systems biology               Schmidt, Wolfgang       2014  
#> 4 Tumor Biology and Oncology         Osieka, R.; Schmidt, C. 1986  
#> 5 Molecular and Cellular Biology     Schmidt, Maxine G       2008  
#> 
#> $facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "112"
#> attr(,"q")
#> [1] "creator:schmidt title:biology"
#> attr(,"fl")
#> [1] "dctitle,dccreator,dcyear"
#> attr(,"start")
#> [1] 5
#> attr(,"fq")
#> [1] "-collection:(ftpakistanrl OR ftubheidelojs OR ftunivzuliaojs OR ftsaludcuba OR ftunivsriwijaya OR ftunivmanitobao2 OR ftjba OR ftunivcasablanca OR ftulbbonndc OR ftjmethode OR ftjqe OR ftjajiks OR ftunivmendozaojs OR ftunisafricaojs OR ftunivpelojs OR ftpenamultimedia OR fttoobunivet OR ftjpjmd OR fthamzanwadiuniv OR ftjedu OR ftjberumpun OR ftjavacient OR ftiaingorontalo OR ftstmikpelitanus OR ftunivfssparaojs OR crunivchicagopr OR crjohnbenjaminsp OR ftstikesbanyuwan OR ftjofg)"
#> attr(,"rows")
#> [1] "5"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 504
#> attr(,"maxScore")
#> [1] "6.60967"
```

term "unix" and published between 1983 and 2009, sorted by year
of publication (dcyear) in descending order


```r
(res <- bs_search(query = 'unix dcyear:[1983 TO 2009]',
  sortby = 'dcyear desc'))
#> $docs
#> # A tibble: 10 x 36
#>    dchdate dcdocid dccontinent dccountry dccollection dcprovider dctitle
#>    <chr>   <chr>   <chr>       <chr>     <chr>        <chr>      <chr>  
#>  1 2019-0… bd4d83… ceu         ch        ftmdpi       MDPI Open… An Imm…
#>  2 2018-1… 8ebef1… ceu         cz        ftunivzlin   Univerzit… Archiv…
#>  3 2020-0… 166b3b… ceu         by        ftbelarusia… Белорусск… Концеп…
#>  4 2019-0… 4f0a72… ceu         eu        fteudl       EUDL Euro… Using …
#>  5 2017-0… f98996… ceu         pl        ftunivwrocl… Bibliotek… Automa…
#>  6 2019-0… 3ad36f… cna         us        ftnasantrs   NASA Tech… Deskto…
#>  7 2015-0… 683c71… cna         us        fthighwire   HighWire … NAViGa…
#>  8 2020-0… ff0126… ceu         by        fttunivgomel Гомельски… Операц…
#>  9 2019-0… 9cc65f… ceu         de        ftmonarchch… TU Chemni… UNIX-S…
#> 10 2018-1… 4ad3ba… ceu         de        ftdara       da|ra - R… Sacram…
#> # … with 29 more variables: dccreator <chr>, dcperson <chr>, dcsubject <chr>,
#> #   dcdescription <chr>, dcpublisher <chr>, dcdate <chr>, dcyear <chr>,
#> #   dctype <chr>, dctypenorm <chr>, dcformat <chr>, dccontenttype <chr>,
#> #   dcidentifier <chr>, dclink <chr>, dcsource <chr>, dclanguage <chr>,
#> #   dcrelation <chr>, dcrights <chr>, dcrightsnorm <chr>,
#> #   dcautoclasscode <chr>, dcdeweyfull <chr>, dcdeweyhuns <chr>,
#> #   dcdeweytens <chr>, dcdeweyones <chr>, dcdoi <chr>, dcoa <chr>,
#> #   dclang <chr>, dccontributor <chr>, dccoverage <chr>, dcclasscode <chr>
#> 
#> $facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "414"
#> attr(,"q")
#> [1] "unix year:[1983 TO 2009]"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"sort")
#> [1] "dcyear_sort desc"
#> attr(,"fq")
#> [1] "-collection:(ftpakistanrl OR ftubheidelojs OR ftunivzuliaojs OR ftsaludcuba OR ftunivsriwijaya OR ftunivmanitobao2 OR ftjba OR ftunivcasablanca OR ftulbbonndc OR ftjmethode OR ftjqe OR ftjajiks OR ftunivmendozaojs OR ftunisafricaojs OR ftunivpelojs OR ftpenamultimedia OR fttoobunivet OR ftjpjmd OR fthamzanwadiuniv OR ftjedu OR ftjberumpun OR ftjavacient OR ftiaingorontalo OR ftstmikpelitanus OR ftunivfssparaojs OR crunivchicagopr OR crjohnbenjaminsp OR ftstikesbanyuwan OR ftjofg)"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 11753
#> attr(,"start")
#> [1] 0
```

raw XML output


```r
bs_search(target = 'ftubbiepub', query = 'lossau summann', raw = TRUE)
#> [1] "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<response>\n<lst name=\"responseHeader\"><bool name=\"zkConnected\">true</bool><int name=\"status\">0</int><int name=\"QTime\">10</int><lst name=\"params\"><str name=\"q\">lossau summann</str><str name=\"fl\">dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm</str><str name=\"fq\">collection:ftubbiepub</str></lst></lst><result name=\"response\" numFound=\"2\" start=\"0\" maxScore=\"2.395762\"><doc><date name=\"dchdate\">2020-08-05T10:49:47Z</date><str name=\"dcdocid\">c190469010e21e37ce99aa75d6d0d423b20d535c2b51242c7dac1eb5fa729691</str><str name=\"dccontinent\">ceu</str><str name=\"dccountry\">de</str><str name=\"dccollection\">ftubbiepub</str><str name=\"dcprovider\">PUB - Publikationen an der Universität Bielefeld</str><str name=\"dctitle\">Search engine technology and digital libraries - Moving from theory to practice</str><arr name=\"dccreator\"><str>Summann, Friedrich</str><str>Lossau, Norbert</str></arr><arr name=\"dcperson\"><str>Summann, Friedrich</str><str>Lossau, Norbert</str></arr><arr name=\"dcsubject\"><str>Elektronische Bibliothek</str><str>Suchmaschine</str><str>Information Retrieval</str><str>Bielefeld Academic Search Engine</str><str>ddc:020</str></arr><str name=\"dcdescription\">Summann F, Lossau N. Search engine technology and digital libraries - Moving from theory to practice. D-lib magazine . 2004;10(9). ; This article describes the development of a modern search engine-based information retrieval setting from its envisioning through to its technological realization. The author takes up the thread of an earlier article on this subject (ZfBB 51 (2004), 5/6), this time from a technical viewpoint. After presenting the conceptual considerations of the initial stages,this article deals principally with the technological aspects of the project.</str><arr name=\"dcpublisher\"><str>CNRI Acct</str></arr><str name=\"dcdate\">2004</str><int name=\"dcyear\">2004</int><arr name=\"dctype\"><str>info:eu-repo/semantics/article</str><str>doc-type:article</str><str>text</str></arr><arr name=\"dctypenorm\"><str>121</str></arr><arr name=\"dcidentifier\"><str>https://nbn-resolving.org/urn:nbn:de:0070-pub-16809824</str><str>https://pub.uni-bielefeld.de/record/1680982</str><str>https://pub.uni-bielefeld.de/download/1680982/2315314</str></arr><str name=\"dclink\">https://nbn-resolving.org/urn:nbn:de:0070-pub-16809824</str><arr name=\"dclanguage\"><str>eng</str></arr><arr name=\"dcrelation\"><str>info:eu-repo/semantics/altIdentifier/doi/10.1045/september2004-lossau</str><str>info:eu-repo/semantics/altIdentifier/issn/1082-9873</str><str>https://nbn-resolving.org/urn:nbn:de:0070-pub-16809824</str><str>https://pub.uni-bielefeld.de/record/1680982</str><str>https://pub.uni-bielefeld.de/download/1680982/2315314</str></arr><str name=\"dcrights\">info:eu-repo/semantics/openAccess ; https://rightsstatements.org/vocab/InC/1.0/</str><arr name=\"dcauthorid\"><str>Summann, Friedrich | orcid:0000-0002-6297-3348</str></arr><arr name=\"dcorcid\"><str>0000-0002-6297-3348</str></arr><arr name=\"dcdeweyfull\"><str>020</str></arr><arr name=\"dcdeweyhuns\"><str>0</str></arr><arr name=\"dcdeweytens\"><str>02</str></arr><arr name=\"dcdeweyones\"><str>020</str></arr><arr name=\"dcclasscode\"><str>020</str></arr><arr name=\"dcdoi\"><str>10.1045/september2004-lossau</str></arr><int name=\"dcoa\">1</int><arr name=\"dclang\"><str>eng</str></arr></doc><doc><date name=\"dchdate\">2020-08-05T10:53:05Z</date><str name=\"dcdocid\">88f9f2b6473596dca40499c60c37157b8c5b7a7e6e431d1ab7b20067032d8d3f</str><str name=\"dccontinent\">ceu</str><str name=\"dccountry\">de</str><str name=\"dccollection\">ftubbiepub</str><str name=\"dcprovider\">PUB - Publikationen an der Universität Bielefeld</str><str name=\"dctitle\">Suchmaschinentechnologie und Digitale Bibliotheken: Von der Theorie zur Praxis</str><arr name=\"dccreator\"><str>Summann, Friedrich</str><str>Lossau, Norbert</str></arr><arr name=\"dcperson\"><str>Summann, Friedrich</str><str>Lossau, Norbert</str></arr><arr name=\"dcsubject\"><str>Informatik</str><str>Massenmedien</str><str>Book</str><str>Bielefeld Academic Search Engine</str><str>Library</str><str>Documentation</str><str>Communication</str><str>Informatics</str><str>Libro</str><str>Biblioteca</str><str>Kommunikation</str><str>Documentazione</str><str>Comunicazione</str><str>Dokumentation</str><str>Informatica</str><str>Buch- und Bibliothekswesen</str><str>ddc:004</str></arr><str name=\"dcdescription\">Summann F, Lossau N. Suchmaschinentechnologie und Digitale Bibliotheken: Von der Theorie zur Praxis. Zeitschrift für Bibliothekswesen und Bibliographie (ZfBB) . 2005;52(1):13-19. ; Der folgende Aufsatz beschreibt aus technischer Sicht den Weg von der Konzeption und Vision einer modernen suchmaschinenbasierten Suchumgebung zu ihrer technologischen Umsetzung. Er nimmt den Faden, der im ersten Teil ( ZfBB 51 (2004), 5/6) beschrieben wurde, unter technischen Gesichtspunkten wieder auf. Dabei werden neben den konzeptionellen Ausgangsüberlegungen schwerpunktmäßig die technologischen Aspekte beleuchtet. ; This article describes the development of a modern search enginebased information retrieval setting from its envisioning through to its technological realization. The author takes up the thread of an earlier article on this subject (ZfBB 51 (2004), 5/6), this time from a technical viewpoint. After presenting the conceptual considerations of the initial stages, this article deals principally with the technological aspects of the project.</str><arr name=\"dcpublisher\"><str>Vittorio Klostermann GMBH</str></arr><str name=\"dcdate\">2005</str><int name=\"dcyear\">2005</int><arr name=\"dctype\"><str>info:eu-repo/semantics/article</str><str>doc-type:article</str><str>text</str></arr><arr name=\"dctypenorm\"><str>121</str></arr><arr name=\"dcidentifier\"><str>https://nbn-resolving.org/urn:nbn:de:0070-pub-25166732</str><str>https://pub.uni-bielefeld.de/record/2516673</str><str>https://pub.uni-bielefeld.de/download/2516673/2516679</str></arr><str name=\"dclink\">https://nbn-resolving.org/urn:nbn:de:0070-pub-25166732</str><arr name=\"dclanguage\"><str>deu</str></arr><arr name=\"dcrelation\"><str>info:eu-repo/semantics/altIdentifier/issn/0044-2380</str><str>https://nbn-resolving.org/urn:nbn:de:0070-pub-25166732</str><str>https://pub.uni-bielefeld.de/record/2516673</str><str>https://pub.uni-bielefeld.de/download/2516673/2516679</str></arr><str name=\"dcrights\">info:eu-repo/semantics/openAccess ; https://rightsstatements.org/vocab/InC/1.0/</str><arr name=\"dcautoclasscode\"><str>020</str></arr><arr name=\"dcauthorid\"><str>Summann, Friedrich | orcid:0000-0002-6297-3348</str></arr><arr name=\"dcorcid\"><str>0000-0002-6297-3348</str></arr><arr name=\"dcdeweyfull\"><str>004</str><str>020</str></arr><arr name=\"dcdeweyhuns\"><str>0</str><str>0</str></arr><arr name=\"dcdeweytens\"><str>00</str><str>02</str></arr><arr name=\"dcdeweyones\"><str>004</str><str>020</str></arr><arr name=\"dcclasscode\"><str>004</str></arr><int name=\"dcoa\">1</int><arr name=\"dclang\"><str>ger</str></arr></doc></result>\n</response>\n"
```

list output


```r
bs_search(target = 'ftubbiepub', query = 'lossau summann', parse = "list")
#> $docs
#> $docs[[1]]
#> $docs[[1]]$dchdate
#> [1] "2020-08-05T10:49:47Z"
#> 
#> $docs[[1]]$dcdocid
#> [1] "c190469010e21e37ce99aa75d6d0d423b20d535c2b51242c7dac1eb5fa729691"
#> 
#> $docs[[1]]$dccontinent
#> [1] "ceu"
#> 
#> $docs[[1]]$dccountry
#> [1] "de"
#> 
#> $docs[[1]]$dccollection
#> [1] "ftubbiepub"
#> 
#> $docs[[1]]$dcprovider
#> [1] "PUB - Publikationen an der Universität Bielefeld"
#> 
#> $docs[[1]]$dctitle
#> [1] "Search engine technology and digital libraries - Moving from theory to practice"
#> 
#> $docs[[1]]$dccreator
#> [1] "Summann, Friedrich; Lossau, Norbert"
#> 
#> $docs[[1]]$dcperson
#> [1] "Summann, Friedrich; Lossau, Norbert"
#> 
#> $docs[[1]]$dcsubject
#> [1] "Elektronische Bibliothek; Suchmaschine; Information Retrieval; Bielefeld Academic Search Engine; ddc:020"
#> 
#> $docs[[1]]$dcdescription
#> [1] "Summann F, Lossau N. Search engine technology and digital libraries - Moving from theory to practice. D-lib magazine . 2004;10(9). ; This article describes the development of a modern search engine-based information retrieval setting from its envisioning through to its technological realization. The author takes up the thread of an earlier article on this subject (ZfBB 51 (2004), 5/6), this time from a technical viewpoint. After presenting the conceptual considerations of the initial stages,this article deals principally with the technological aspects of the project."
#> 
#> $docs[[1]]$dcpublisher
#> [1] "CNRI Acct"
#> 
#> $docs[[1]]$dcdate
#> [1] "2004"
#> 
#> $docs[[1]]$dcyear
#> [1] "2004"
#> 
#> $docs[[1]]$dctype
#> [1] "info:eu-repo/semantics/article; doc-type:article; text"
#> 
#> $docs[[1]]$dctypenorm
#> [1] "121"
#> 
#> $docs[[1]]$dcidentifier
#> [1] "https://nbn-resolving.org/urn:nbn:de:0070-pub-16809824; https://pub.uni-bielefeld.de/record/1680982; https://pub.uni-bielefeld.de/download/1680982/2315314"
#> 
#> $docs[[1]]$dclink
#> [1] "https://nbn-resolving.org/urn:nbn:de:0070-pub-16809824"
#> 
#> $docs[[1]]$dclanguage
#> [1] "eng"
#> 
#> $docs[[1]]$dcrelation
#> [1] "info:eu-repo/semantics/altIdentifier/doi/10.1045/september2004-lossau; info:eu-repo/semantics/altIdentifier/issn/1082-9873; https://nbn-resolving.org/urn:nbn:de:0070-pub-16809824; https://pub.uni-bielefeld.de/record/1680982; https://pub.uni-bielefeld.de/download/1680982/2315314"
#> 
#> $docs[[1]]$dcrights
#> [1] "info:eu-repo/semantics/openAccess ; https://rightsstatements.org/vocab/InC/1.0/"
#> 
#> $docs[[1]]$dcauthorid
#> [1] "Summann, Friedrich | orcid:0000-0002-6297-3348"
#> 
#> $docs[[1]]$dcorcid
#> [1] "0000-0002-6297-3348"
#> 
#> $docs[[1]]$dcdeweyfull
#> [1] "020"
#> 
#> $docs[[1]]$dcdeweyhuns
#> [1] "0"
#> 
#> $docs[[1]]$dcdeweytens
#> [1] "02"
#> 
#> $docs[[1]]$dcdeweyones
#> [1] "020"
#> 
#> $docs[[1]]$dcclasscode
#> [1] "020"
#> 
#> $docs[[1]]$dcdoi
#> [1] "10.1045/september2004-lossau"
#> 
#> $docs[[1]]$dcoa
#> [1] "1"
#> 
#> $docs[[1]]$dclang
#> [1] "eng"
#> 
#> 
#> $docs[[2]]
#> $docs[[2]]$dchdate
#> [1] "2020-08-05T10:53:05Z"
#> 
#> $docs[[2]]$dcdocid
#> [1] "88f9f2b6473596dca40499c60c37157b8c5b7a7e6e431d1ab7b20067032d8d3f"
#> 
#> $docs[[2]]$dccontinent
#> [1] "ceu"
#> 
#> $docs[[2]]$dccountry
#> [1] "de"
#> 
#> $docs[[2]]$dccollection
#> [1] "ftubbiepub"
#> 
#> $docs[[2]]$dcprovider
#> [1] "PUB - Publikationen an der Universität Bielefeld"
#> 
#> $docs[[2]]$dctitle
#> [1] "Suchmaschinentechnologie und Digitale Bibliotheken: Von der Theorie zur Praxis"
#> 
#> $docs[[2]]$dccreator
#> [1] "Summann, Friedrich; Lossau, Norbert"
#> 
#> $docs[[2]]$dcperson
#> [1] "Summann, Friedrich; Lossau, Norbert"
#> 
#> $docs[[2]]$dcsubject
#> [1] "Informatik; Massenmedien; Book; Bielefeld Academic Search Engine; Library; Documentation; Communication; Informatics; Libro; Biblioteca; Kommunikation; Documentazione; Comunicazione; Dokumentation; Informatica; Buch- und Bibliothekswesen; ddc:004"
#> 
#> $docs[[2]]$dcdescription
#> [1] "Summann F, Lossau N. Suchmaschinentechnologie und Digitale Bibliotheken: Von der Theorie zur Praxis. Zeitschrift für Bibliothekswesen und Bibliographie (ZfBB) . 2005;52(1):13-19. ; Der folgende Aufsatz beschreibt aus technischer Sicht den Weg von der Konzeption und Vision einer modernen suchmaschinenbasierten Suchumgebung zu ihrer technologischen Umsetzung. Er nimmt den Faden, der im ersten Teil ( ZfBB 51 (2004), 5/6) beschrieben wurde, unter technischen Gesichtspunkten wieder auf. Dabei werden neben den konzeptionellen Ausgangsüberlegungen schwerpunktmäßig die technologischen Aspekte beleuchtet. ; This article describes the development of a modern search enginebased information retrieval setting from its envisioning through to its technological realization. The author takes up the thread of an earlier article on this subject (ZfBB 51 (2004), 5/6), this time from a technical viewpoint. After presenting the conceptual considerations of the initial stages, this article deals principally with the technological aspects of the project."
#> 
#> $docs[[2]]$dcpublisher
#> [1] "Vittorio Klostermann GMBH"
#> 
#> $docs[[2]]$dcdate
#> [1] "2005"
#> 
#> $docs[[2]]$dcyear
#> [1] "2005"
#> 
#> $docs[[2]]$dctype
#> [1] "info:eu-repo/semantics/article; doc-type:article; text"
#> 
#> $docs[[2]]$dctypenorm
#> [1] "121"
#> 
#> $docs[[2]]$dcidentifier
#> [1] "https://nbn-resolving.org/urn:nbn:de:0070-pub-25166732; https://pub.uni-bielefeld.de/record/2516673; https://pub.uni-bielefeld.de/download/2516673/2516679"
#> 
#> $docs[[2]]$dclink
#> [1] "https://nbn-resolving.org/urn:nbn:de:0070-pub-25166732"
#> 
#> $docs[[2]]$dclanguage
#> [1] "deu"
#> 
#> $docs[[2]]$dcrelation
#> [1] "info:eu-repo/semantics/altIdentifier/issn/0044-2380; https://nbn-resolving.org/urn:nbn:de:0070-pub-25166732; https://pub.uni-bielefeld.de/record/2516673; https://pub.uni-bielefeld.de/download/2516673/2516679"
#> 
#> $docs[[2]]$dcrights
#> [1] "info:eu-repo/semantics/openAccess ; https://rightsstatements.org/vocab/InC/1.0/"
#> 
#> $docs[[2]]$dcautoclasscode
#> [1] "020"
#> 
#> $docs[[2]]$dcauthorid
#> [1] "Summann, Friedrich | orcid:0000-0002-6297-3348"
#> 
#> $docs[[2]]$dcorcid
#> [1] "0000-0002-6297-3348"
#> 
#> $docs[[2]]$dcdeweyfull
#> [1] "004; 020"
#> 
#> $docs[[2]]$dcdeweyhuns
#> [1] "0; 0"
#> 
#> $docs[[2]]$dcdeweytens
#> [1] "00; 02"
#> 
#> $docs[[2]]$dcdeweyones
#> [1] "004; 020"
#> 
#> $docs[[2]]$dcclasscode
#> [1] "004"
#> 
#> $docs[[2]]$dcoa
#> [1] "1"
#> 
#> $docs[[2]]$dclang
#> [1] "ger"
#> 
#> 
#> 
#> $facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "9"
#> attr(,"q")
#> [1] "lossau summann"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"fq")
#> [1] "collection:ftubbiepub"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 2
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "2.395762"
```


```r
out <- list()
system.time(
  for (i in 1:3) {
    out[[i]] <- bs_search(target = 'ftubbiepub', query = 'lossau summann',
      hits = 1)
  }
)
#>    user  system elapsed 
#>   0.077   0.004   5.134
out
#> [[1]]
#> [[1]]$docs
#> # A tibble: 1 x 31
#>   dchdate dcdocid dccontinent dccountry dccollection dcprovider dctitle
#>   <chr>   <chr>   <chr>       <chr>     <chr>        <chr>      <chr>  
#> 1 2020-0… c19046… ceu         de        ftubbiepub   PUB - Pub… Search…
#> # … with 24 more variables: dccreator <chr>, dcperson <chr>, dcsubject <chr>,
#> #   dcdescription <chr>, dcpublisher <chr>, dcdate <chr>, dcyear <chr>,
#> #   dctype <chr>, dctypenorm <chr>, dcidentifier <chr>, dclink <chr>,
#> #   dclanguage <chr>, dcrelation <chr>, dcrights <chr>, dcauthorid <chr>,
#> #   dcorcid <chr>, dcdeweyfull <chr>, dcdeweyhuns <chr>, dcdeweytens <chr>,
#> #   dcdeweyones <chr>, dcclasscode <chr>, dcdoi <chr>, dcoa <chr>, dclang <chr>
#> 
#> [[1]]$facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "8"
#> attr(,"q")
#> [1] "lossau summann"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"fq")
#> [1] "collection:ftubbiepub"
#> attr(,"rows")
#> [1] "1"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 2
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "2.395762"
#> 
#> [[2]]
#> [[2]]$docs
#> # A tibble: 1 x 31
#>   dchdate dcdocid dccontinent dccountry dccollection dcprovider dctitle
#>   <chr>   <chr>   <chr>       <chr>     <chr>        <chr>      <chr>  
#> 1 2020-0… c19046… ceu         de        ftubbiepub   PUB - Pub… Search…
#> # … with 24 more variables: dccreator <chr>, dcperson <chr>, dcsubject <chr>,
#> #   dcdescription <chr>, dcpublisher <chr>, dcdate <chr>, dcyear <chr>,
#> #   dctype <chr>, dctypenorm <chr>, dcidentifier <chr>, dclink <chr>,
#> #   dclanguage <chr>, dcrelation <chr>, dcrights <chr>, dcauthorid <chr>,
#> #   dcorcid <chr>, dcdeweyfull <chr>, dcdeweyhuns <chr>, dcdeweytens <chr>,
#> #   dcdeweyones <chr>, dcclasscode <chr>, dcdoi <chr>, dcoa <chr>, dclang <chr>
#> 
#> [[2]]$facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "11"
#> attr(,"q")
#> [1] "lossau summann"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"fq")
#> [1] "collection:ftubbiepub"
#> attr(,"rows")
#> [1] "1"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 2
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "2.395762"
#> 
#> [[3]]
#> [[3]]$docs
#> # A tibble: 1 x 31
#>   dchdate dcdocid dccontinent dccountry dccollection dcprovider dctitle
#>   <chr>   <chr>   <chr>       <chr>     <chr>        <chr>      <chr>  
#> 1 2020-0… c19046… ceu         de        ftubbiepub   PUB - Pub… Search…
#> # … with 24 more variables: dccreator <chr>, dcperson <chr>, dcsubject <chr>,
#> #   dcdescription <chr>, dcpublisher <chr>, dcdate <chr>, dcyear <chr>,
#> #   dctype <chr>, dctypenorm <chr>, dcidentifier <chr>, dclink <chr>,
#> #   dclanguage <chr>, dcrelation <chr>, dcrights <chr>, dcauthorid <chr>,
#> #   dcorcid <chr>, dcdeweyfull <chr>, dcdeweyhuns <chr>, dcdeweytens <chr>,
#> #   dcdeweyones <chr>, dcclasscode <chr>, dcdoi <chr>, dcoa <chr>, dclang <chr>
#> 
#> [[3]]$facets
#> list()
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "11"
#> attr(,"q")
#> [1] "lossau summann"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"fq")
#> [1] "collection:ftubbiepub"
#> attr(,"rows")
#> [1] "1"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 2
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "2.395762"
```

### Faceting


```r
bs_search(query = "unix", facets = c("dcsubject", "dcyear"),
  facet_limit = 10)
#> Warning: `data_frame()` is deprecated as of tibble 1.1.0.
#> Please use `tibble()` instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_warnings()` to see where this warning was generated.
#> $docs
#> # A tibble: 10 x 15
#>    dcdocid dchdate dccollection dcprovider dctitle dcdescription dctype
#>    <chr>   <chr>   <chr>        <chr>      <chr>   <chr>         <chr> 
#>  1 cc7b9e… 2020-1… ftwikibooks  WikiBooks… Wikibo… "UNIX Basics… Book  
#>  2 89ddd7… 2020-1… ftwikibooks  WikiBooks… Wikibo… "Unix on lai… Book  
#>  3 211eab… 2020-1… ftwikibooks  WikiBooks… Wikibo… ""            Book  
#>  4 884007… 2020-1… ftwikibooks  WikiBooks… Wikibo… "This is the… Book  
#>  5 b8b4e1… 2020-1… ftwikibooks  WikiBooks… Wikibo… "О чём эта к… Book  
#>  6 820d1a… 2020-1… ftwikibooks  WikiBooks… Wikibo… "TODO Talk a… Book  
#>  7 8c8c49… 2020-1… ftwikibooks  WikiBooks… Wikibo… "Mac os X es… Book  
#>  8 a5a3ad… 2020-1… ftwikibooks  WikiBooks… Wikibo… "¿Qué signif… Book  
#>  9 0d1778… 2020-1… ftwikibooks  WikiBooks… Wikibo… "Unix/Linux"  Book  
#> 10 670021… 2020-1… ftwikibooks  WikiBooks… Wikibo… "stub Ubuntu… Book  
#> # … with 8 more variables: dctypenorm <chr>, dclanguage <chr>, dclang <chr>,
#> #   dcidentifier <chr>, dclink <chr>, dccontinent <chr>, dccountry <chr>,
#> #   dcoa <chr>
#> 
#> $facets
#> $facets$dcsubject
#> # A tibble: 10 x 2
#>    name                                      value
#>    <chr>                                     <chr>
#>  1 Computer Programming and Software         578  
#>  2 UNIX                                      316  
#>  3 Software                                  278  
#>  4 Unix                                      277  
#>  5 Computing and Computers                   217  
#>  6 Computing                                 191  
#>  7 COMPUTING                                 190  
#>  8 And Information Science                   189  
#>  9 Linux                                     188  
#> 10 99 General And Miscellaneous//Mathematics 187  
#> 
#> $facets$dcyear
#> # A tibble: 10 x 2
#>    name  value
#>    <chr> <chr>
#>  1 1994  786  
#>  2 1995  688  
#>  3 1993  673  
#>  4 2018  599  
#>  5 1996  586  
#>  6 1997  573  
#>  7 1998  556  
#>  8 1992  544  
#>  9 2017  521  
#> 10 1999  517  
#> 
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "177"
#> attr(,"q")
#> [1] "unix"
#> attr(,"facet.limit")
#> [1] "10"
#> attr(,"facet.field")
#> [1] "f_dcsubjectf_dcyear"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"facet.mincount")
#> [1] "1"
#> attr(,"fq")
#> [1] "-collection:(ftpakistanrl OR ftubheidelojs OR ftunivzuliaojs OR ftsaludcuba OR ftunivsriwijaya OR ftunivmanitobao2 OR ftjba OR ftunivcasablanca OR ftulbbonndc OR ftjmethode OR ftjqe OR ftjajiks OR ftunivmendozaojs OR ftunisafricaojs OR ftunivpelojs OR ftpenamultimedia OR fttoobunivet OR ftjpjmd OR fthamzanwadiuniv OR ftjedu OR ftjberumpun OR ftjavacient OR ftiaingorontalo OR ftstmikpelitanus OR ftunivfssparaojs OR crunivchicagopr OR crjohnbenjaminsp OR ftstikesbanyuwan OR ftjofg)"
#> attr(,"facet")
#> [1] "true"
#> attr(,"facet.sort")
#> [1] "count"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 21343
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "3.6168933"
bs_search(query = "unix", facets = c("dcsubject", "dcyear"),
  f_dcsubject = '"computer science"',
  facet_limit = 10, verbose = TRUE)
#> $docs
#> # A tibble: 10 x 33
#>    dchdate dcdocid dccontinent dccountry dccollection dcprovider dctitle
#>    <chr>   <chr>   <chr>       <chr>     <chr>        <chr>      <chr>  
#>  1 2019-1… 8631c6… cww         org       ftdatacite   DataCite … GFT: A…
#>  2 2020-0… e428bc… cww         org       ftdatacite   DataCite … GFT: A…
#>  3 2020-0… 6c572c… ceu         ch        ftethz       ETH Züric… GFT: A…
#>  4 2019-1… 124a35… cww         org       ftdatacite   DataCite … Shared…
#>  5 2019-1… 8b4ec3… ceu         eu        ftopengrey   Open Grey… Standa…
#>  6 2020-0… cb3c89… cww         org       ftdatacite   DataCite … Shared…
#>  7 2020-0… e981b2… cww         org       ftdatacite   DataCite … HP-Obe…
#>  8 2019-1… c9d33a… cww         org       ftdatacite   DataCite … HP-Obe…
#>  9 2019-1… 7aede2… ceu         eu        ftopengrey   Open Grey… Depend…
#> 10 2020-0… a6467e… ceu         ch        ftethz       ETH Züric… HP-Obe…
#> # … with 26 more variables: dccreator <chr>, dcperson <chr>, dcsubject <chr>,
#> #   dcdescription <chr>, dcpublisher <chr>, dcdate <chr>, dcyear <chr>,
#> #   dctype <chr>, dctypenorm <chr>, dcformat <chr>, dcidentifier <chr>,
#> #   dclink <chr>, dclanguage <chr>, dcrelation <chr>, dcrights <chr>,
#> #   dcdoi <chr>, dcoa <chr>, dclang <chr>, dccontenttype <chr>, dcsource <chr>,
#> #   dcdeweyfull <chr>, dcdeweyhuns <chr>, dcdeweytens <chr>, dcdeweyones <chr>,
#> #   dcclasscode <chr>, dccontributor <chr>
#> 
#> $facets
#> $facets$dcsubject
#> # A tibble: 10 x 2
#>    name                                                                    value
#>    <chr>                                                                   <chr>
#>  1 computer science                                                        37   
#>  2 Data processing                                                         14   
#>  3 UNIX (BETRIEBSSYSTEME)                                                  14   
#>  4 UNIX (OPERATING SYSTEMS)                                                14   
#>  5 info:eu-repo/classification/ddc/004                                     14   
#>  6 technical report                                                        8    
#>  7 090 - Electronics and electrical engineering                            6    
#>  8 DISTRIBUTED APPLICATIONS + CLOUD COMPUTING + GRID COMPUTING (COMPUTER … 6    
#>  9 MINICOMPUTER + WORKSTATIONS (COMPUTERSYSTEME)                           6    
#> 10 MINICOMPUTERS + WORKSTATIONS (COMPUTER SYSTEMS)                         6    
#> 
#> $facets$dcyear
#> # A tibble: 10 x 2
#>    name  value
#>    <chr> <chr>
#>  1 1994  8    
#>  2 1990  3    
#>  3 1993  3    
#>  4 2003  3    
#>  5 1981  2    
#>  6 1983  2    
#>  7 1995  2    
#>  8 1984  1    
#>  9 1985  1    
#> 10 1986  1    
#> 
#> 
#> attr(,"status")
#> [1] 0
#> attr(,"QTime")
#> [1] "136"
#> attr(,"q")
#> [1] "unix"
#> attr(,"facet.limit")
#> [1] "10"
#> attr(,"facet.field")
#> [1] "f_dcsubjectf_dcyear"
#> attr(,"fl")
#> [1] "dccollection,dccontenttype,dccontinent,dccountry,dccreator,dcauthorid,dcdate,dcdescription,dcdocid,dcdoi,dcformat,dcidentifier,dclang,dclanguage,dclink,dcorcid,dcperson,dcpublisher,dcrights,dcsource,dcsubject,dctitle,dcyear,dctype,dcclasscode,dctypenorm,dcdeweyfull,dcdeweyhuns,dcdeweytens,dcdeweyones,dcautoclasscode,dcrelation,dccontributor,dccoverage,dchdate,dcoa,dcrightsnorm"
#> attr(,"facet.mincount")
#> [1] "1"
#> attr(,"fq")
#> [1] "f_dcsubject:\"computer science\"-collection:(ftpakistanrl OR ftubheidelojs OR ftunivzuliaojs OR ftsaludcuba OR ftunivsriwijaya OR ftunivmanitobao2 OR ftjba OR ftunivcasablanca OR ftulbbonndc OR ftjmethode OR ftjqe OR ftjajiks OR ftunivmendozaojs OR ftunisafricaojs OR ftunivpelojs OR ftpenamultimedia OR fttoobunivet OR ftjpjmd OR fthamzanwadiuniv OR ftjedu OR ftjberumpun OR ftjavacient OR ftiaingorontalo OR ftstmikpelitanus OR ftunivfssparaojs OR crunivchicagopr OR crjohnbenjaminsp OR ftstikesbanyuwan OR ftjofg)"
#> attr(,"facet")
#> [1] "true"
#> attr(,"facet.sort")
#> [1] "count"
#> attr(,"name")
#> [1] "response"
#> attr(,"numFound")
#> [1] 37
#> attr(,"start")
#> [1] 0
#> attr(,"maxScore")
#> [1] "1.6781588"
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_search.R
\name{bs_retry_options}
\alias{bs_retry_options}
\title{bs_search retry options}
\usage{
bs_retry_options(
  pause_base = 1,
  pause_cap = 60,
  pause_min = 1,
  times = 3,
  terminate_on = NULL,
  retry_only_on = NULL,
  onwait = NULL
)
}
\arguments{
\item{pause_base, pause_cap, pause_min}{basis, maximum, and minimum for
calculating wait time for retry.}

\item{times}{the maximum number of times to retry.}

\item{terminate_on, retry_only_on}{a vector of HTTP status codes.}

\item{onwait}{a callback function if the request will be retried and
a wait time is being applied.}
}
\value{
a named list with the parameters given to this function
}
\description{
bs_search retry options
}
\details{
see \link[crul:HttpClient]{crul::HttpClient} for more detailed explanation of these
parameters
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_profile.R
\name{bs_profile}
\alias{bs_profile}
\title{Get the profile for a repository}
\usage{
bs_profile(target, ...)
}
\arguments{
\item{target}{(character) Internal name of a single repository as
delivered in \code{\link[=bs_repositories]{bs_repositories()}}}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a data.frame, of two columns: "name", "value".
"name" holds the "value" description. you can pivot the
data.frame to wide by e.g., \code{tidyr::pivot_wider(x)}
}
\description{
Get the profile for a repository
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_search.R
\name{bs_search}
\alias{bs_search}
\alias{bs_meta}
\title{Search BASE}
\usage{
bs_search(
  query = NULL,
  target = NULL,
  coll = NULL,
  boost_oa = FALSE,
  hits = NULL,
  offset = NULL,
  fields = NULL,
  sortby = NULL,
  facets = NULL,
  facet_limit = 100,
  facet_sort = NULL,
  filter = NULL,
  raw = FALSE,
  parse = "df",
  retry = bs_retry_options(),
  ...
)

bs_meta(x)
}
\arguments{
\item{query}{(character) query string. For syntax details see Appendix,
section "Query syntax"}

\item{target}{(character) Internal name of a single repository as
delivered in \code{\link[=bs_repositories]{bs_repositories()}}}

\item{coll}{(character) collection code. For existing, pre-defined
collections see Appendix, section "Collection-related queries"}

\item{boost_oa}{(logical) Push open access documents upwards in the
result list. Default: \code{FALSE}}

\item{hits}{(integer) number of results to return. Default: 10. Max: 100}

\item{offset}{(integer) record to start at. Default: 0. Max: 1000}

\item{fields}{(character) Fields to return. This doesn't appear to be
working though. The result records only contain fields listed in the
comma-separated field list. For existing, pre-defined fields see Appendix,
section "Fields"}

\item{sortby}{(character) field to sort by. A sort ordering must include
a single field name (see Appendix, section "Fields", table column
"Sorting"), followed by a whitespace (escaped as + or \%20 in URL strings),
followed by sort direction (asc or desc). Default: sorts by relevance}

\item{facets}{(character) The response contains an extra section
"facet_counts/facet_fields" with fields from the comma-separated facets
list. This section provides a breakdown or summary of the results. From the
user's perspective, faceted search breaks up search results into multiple
categories, typically showing counts for each, and allows the user to
"drill down" or further restrict their search results based on those facets.
Use of faceting does not affect the results section of a search response.
For existing, pre-defined facet fields see Appendix, section "Fields",
table column "Facet".}

\item{facet_limit}{(numeric) Maximum number of constraint counts that
should be returned for the facet fields. Default: 100; min:1; max: 500}

\item{facet_sort}{(character) Ordering of the facet field constraints:
count - sort by count (highest count first);  index - alphabetical sorting.
Default: count}

\item{filter}{(character) a string with the value to be used. html escaping
will be automatically done; embed string in \code{I()} to avoid html escaping.
This parameter gets used by \code{fq} solr parameter on the server}

\item{raw}{(logical) If \code{TRUE} returns raw XML, default: \code{FALSE}}

\item{parse}{(character) One of 'list' or 'df'}

\item{retry}{(list) use \code{\link[=bs_retry_options]{bs_retry_options()}} to make a named list of
retry options to pass on to the HTTP request. default values are passed
for you, but you can change them by setting an option in
\code{bs_retry_options()}}

\item{...}{Facet field based query options (See Facet below) or curl
options passed on to \link[crul:verb-GET]{crul::verb-GET}}

\item{x}{input to \code{bs_meta}}
}
\value{
XML as character string if \code{parse = FALSE} or data.frame
}
\description{
Search BASE
}
\details{
BASE asks that requests are not more frequent than 1 per second,
so we enforce the rate limit internally. if you do a single request not
in a for loop/lapply type situation, this won't be inoked, but will
if doing a for loop/lapply call, and there's no sleep invoked
}
\section{Facet}{

You can optionally pass in search term for specific facet fields.
See example. For existing, pre-defined facet fields see Appendix at
https://www.base-search.net/about/download/base_interface.pdf,
section "Fields", table column "Facet"
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbace-package.R
\docType{package}
\name{rbace-package}
\alias{rbace-package}
\alias{rbace}
\title{rbace}
\description{
Bielefeld Academic Search Engine Client
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_repositories.R
\name{bs_repositories}
\alias{bs_repositories}
\title{List repositories for a collection}
\usage{
bs_repositories(coll, ...)
}
\arguments{
\item{coll}{(character) collection code. For existing, pre-defined
collections see Appendix, section "Collection-related queries" in
the Appendix of
https://www.base-search.net/about/download/base_interface.pdf}

\item{...}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\value{
a data.frame of two columns: "name", "internal_name"
}
\description{
List repositories for a collection
}
