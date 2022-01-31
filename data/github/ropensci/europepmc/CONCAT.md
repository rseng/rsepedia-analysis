# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
europepmc - R Interface to Europe PMC RESTful Web Service
=== 





[![R build status](https://github.com/ropensci/europepmc/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/europepmc/actions)
[![Build status](https://ci.appveyor.com/api/projects/status/f8xtpvhhr074lk44?svg=true)](https://ci.appveyor.com/project/sckott/europepmc)
[![codecov.io](https://codecov.io/github/ropensci/europepmc/coverage.svg?branch=master)](https://codecov.io/github/ropensci/europepmc?branch=master)
[![cran version](https://www.r-pkg.org/badges/version/europepmc)](https://cran.r-project.org/package=europepmc)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/europepmc)](https://github.com/r-hub/cranlogs.app)
[![](https://badges.ropensci.org/29_status.svg)](https://github.com/ropensci/software-review/issues/29)

europepmc facilitates access to the [Europe PMC RESTful Web
Service](https://europepmc.org/RestfulWebService). The client furthermore supports the [Europe PMC Annotations API](https://europepmc.org/AnnotationsApi) to retrieve text-mined concepts and terms per article.

[Europe PMC](https://europepmc.org/) covers life science literature and
gives access to open access full texts. Europe
PMC ingests all PubMed content and extends its index with other literature and patent sources.

For more infos on Europe PMC, see:

<https://europepmc.org/About>

Levchenko, M., Gou, Y., Graef, F., Hamelers, A., Huang, Z., Ide-Smith, M., … McEntyre, J. (2017). Europe PMC in 2017. Nucleic Acids Research, 46(D1), D1254–D1260. <https://doi.org/10.1093/nar/gkx1005>

## Implemented API methods

This client supports the following API methods from the [Articles RESTful API](https://europepmc.org/RestfulWebService):

|API-Method     |Description                                                                                  |R functions                                |
|:--------------|:--------------------------------------------------------------------------------------------|:------------------------------------------|
|search         |Search Europe PMC and get detailed metadata                                                  |`epmc_search()`, `epmc_details()`, `epmc_search_by_doi()`          |
|profile        |Obtain a summary of hit counts for several Europe PMC databases                              |`epmc_profile()`                           |
|citations      |Load metadata representing citing articles for a given publication                           |`epmc_citations()`                         |
|references     |Retrieve the reference section of a publication                                               |`epmc_refs()`                              |
|databaseLinks  |Get links to biological databases such as UniProt or ENA                                     |`epmc_db()`, `epmc_db_count()`             |
|labslinks      |Access links to Europe PMC provided by third parties                                         |`epmc_lablinks()`, `epmc_lablinks_count()` |
|fullTextXML    |Fetch full-texts deposited in PMC                                                            |`epmc_ftxt()`                              |
|bookXML        |retrieve book XML formatted full text for the Open Access subset of the Europe PMC bookshelf |`epmc_ftxt_book()`                         |

From the [Europe PMC Annotations API](https://europepmc.org/AnnotationsApi):

|API-Method     |Description |R functions |
|:-----------|:-------------|:-------------|
annotationsByArticleIds | Get the annotations contained in the list of articles specified | `epmc_annotations_by_id()` |

## Installation

From CRAN

```r
install.packages("europepmc")
```

The latest development version can be installed using the
[remotes](https://github.com/r-lib/remotes/) package:


```r
require(remotes)
install_github("ropensci/europepmc")
```

Loading into R


```r
library(europepmc)
```

## Search Europe PMC

The search covers both metadata (e.g. abstracts or title) and full texts. To
build your query, please refer to the comprehensive guidance on how to search
Europe PMC: <https://europepmc.org/help>. Provide your query in the Europe
PMC search syntax to `epmc_search()`. 


```r
europepmc::epmc_search(query = '"2019-nCoV" OR "2019nCoV"')
#> # A tibble: 100 × 29
#>    id       source pmid     doi   title authorString journalTitle issue journalVolume
#>    <chr>    <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr> <chr>        
#>  1 33406042 MED    33406042 10.1… 2019… Xiao M, Liu… IEEE/ACM Tr… 4     18           
#>  2 34059225 MED    34059225 10.1… Livi… Santillan-G… Med Intensi… 5     45           
#>  3 34181072 MED    34181072 10.1… Self… Varghese JJ… Support Car… <NA>  <NA>         
#>  4 34108756 MED    34108756 10.2… COVI… Gabarron E,… Bull World … 6     99           
#>  5 33197230 MED    33197230 10.2… Sear… Lazarus JV,… J Med Inter… 11    22           
#>  6 33181701 MED    33181701 10.1… The … Kim YJ, Qia… Medicine (B… 46    99           
#>  7 34291001 MED    34291001 10.4… How … Moradi G, G… Med J Islam… <NA>  35           
#>  8 33521188 MED    33521188 10.1… Tota… Chen AZ, Sh… Arthroplast… <NA>  8            
#>  9 33009914 MED    33009914 10.1… A da… Zhu Z, Meng… Database (O… <NA>  2020         
#> 10 32341597 MED    32341597 10.1… Obes… Carretero G… Rev Clin Es… 6     220          
#> # … with 90 more rows, and 20 more variables: pubYear <chr>, journalIssn <chr>,
#> #   pageInfo <chr>, pubType <chr>, isOpenAccess <chr>, inEPMC <chr>,
#> #   inPMC <chr>, hasPDF <chr>, hasBook <chr>, hasSuppl <chr>,
#> #   citedByCount <int>, hasReferences <chr>, hasTextMinedTerms <chr>,
#> #   hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, pmcid <chr>, versionNumber <int>
```

Be aware that Europe PMC expands queries with MeSH synonyms by default. You can turn this behavior off using the `synonym = FALSE` parameter.

By default, `epmc_search()` returns 100 records. To adjust the limit, simply use
the `limit` parameter.

See vignette [Introducing europepmc, an R interface to Europe PMC RESTful API](https://docs.ropensci.org/europepmc/articles/introducing-europepmc.html) for a long-form documentation about how to search Europe PMC with this client.

## Creating proper review graphs with `epmc_hits_trend()`

There is also a nice function allowing you to easily create review graphs like described in Maëlle
Salmon's [blog post](https://masalmon.eu/2017/05/14/evergreenreviewgraph/):


```r
tt_oa <- europepmc::epmc_hits_trend("Malaria", period = 1995:2019, synonym = FALSE)
tt_oa
#> # A tibble: 25 × 3
#>     year all_hits query_hits
#>    <int>    <dbl>      <dbl>
#>  1  1995   449064       1495
#>  2  1996   458526       1572
#>  3  1997   456744       1873
#>  4  1998   474613       1762
#>  5  1999   493745       1947
#>  6  2000   532019       2092
#>  7  2001   545674       2187
#>  8  2002   561425       2378
#>  9  2003   588572       2612
#> 10  2004   628141       2845
#> # … with 15 more rows
# we use ggplot2 for plotting the graph
library(ggplot2)
ggplot(tt_oa, aes(year, query_hits / all_hits)) + 
  geom_point() + 
  geom_line() +
  xlab("Year published") + 
  ylab("Proportion of articles on Malaria in Europe PMC")
```

![plot of chunk unnamed-chunk-4](man/figures/unnamed-chunk-4-1.png)

For more info, read the vignette about creating literature review graphs:

<https://docs.ropensci.org/europepmc/articles/evergreenreviewgraphs.html>

## Re-use of europepmc

Check out the tidypmc package

<https://github.com/ropensci/tidypmc>

The package maintainer, Chris Stubben (@cstubben), has also created an Shiny App that allows you to search and browse Europe PMC:

<https://github.com/cstubben/euPMC>



## Other ways to access Europe PubMed Central

### Other APIs

- Data dumps: <https://europepmc.org/FtpSite>
- OAI service: <https://europepmc.org/OaiService>
- SOAP web service: <https://europepmc.org/SoapWebServices>
- Grants RESTful (Grist) API: <https://europepmc.org/GristAPI>

### Other R clients

- use rOpenSci's `oai` to get metadata and full text via Europe PMC's OAI interface: <https://github.com/ropensci/oai>
- use rOpenSci's `rentrez` to interact with [NCBI databases](https://www.ncbi.nlm.nih.gov/) such as PubMed: <https://github.com/ropensci/rentrez>
- rOpenSci's `fulltext` package gives access to supplementary material of open access life-science publications in Europe PMC: <https://github.com/ropensci/fulltext>

## Meta

Please note that this project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/europepmc/CONDUCT.html). By participating in this project you agree to abide by its terms.

License: GPL-3

Please use the issue tracker for bug reporting and feature requests.

---

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# europepmc 0.4.1

- Implement API version 6.6 support

Minor changes:

- Allow full-text retrieval only by PMCID, thanks @ESPoppelaars
- Bug fix: `epmc_search()` limit can now be larger than the expected number of results

# europepmc 0.4

- Support API version 6.3
- Europe PMC Annotations ID integration

New functionalities:

- `epmc_annotations_by_id()` get text-mined concepts and full-texts from Europe PMC indexed full-texts;  `epmc_tm()`and `epmc_tm_count()` deprecated, use `epmc_annotations_by_id()` instead
- `epmc_search_by_doi()` get PubMed metadata by DOI names

Minor changes:

- improved retry in case API call fails accidentially
- preprint records collection added
- retrieve full-text by PMID
- improved long-form documentation
- improved testing
- new official docs site
- use *tibble instead of deprecated *data_frame functions

# europepmc 0.3

- Implement API version 6.0

## Minor changes

- improved feedback when calling the API
- link to most current paper from the Europe PMC team

# europepmc 0.2

- Move to HTTPS
- new `epmc_hits_trends()` function to obtain data for review graphs (thanks @maelle)
- new vignette "Making proper trend graphs" and updated search documentation

## Minor changes

- fix sort param
- rename `jsonlite::rbind_pages()` function
- improve `europepmc::epmc_tm()` output

# europepmc 0.1.4

- fixed example in vignette which lead to warnings
- synonym search is operational again

# europepmc 0.1.3

## Minor changes

- [removed explicit API versioning, so that the client now always supports the most recent API version #13](https://github.com/ropensci/europepmc/issues/13)
- set user agent to "ropensci/europepmc"
- `epmc_db()`, `epmc_db_count()`: add PRIDE archive as external database

# europepmc 0.1.2

- cache HTTP 500 errors which sometimes occur and re-try up to five times. It is based on [googlesheet's approach](https://github.com/jennybc/googlesheets/commit/a91403ecb8ab5d8059bf14a9f9878ab68a829f0a)
- new function `epmc_profile()` to get an overview of hit counts for several databases or publication types
- update imported packages in DESCRIPTION

# europepmc 0.1.1

Implement [RESTful API v4.5.3](https://europepmc.org/docs/Europe_PMC_RESTful_Release_Notes.pdf)

## Major changes

- `epmc_search()`: implement cursorMark to paginate through results
- `epmc_search()`: added sort parameter
- `epmc_search()`: [support of `raw` output file #7](https://github.com/ropensci/europepmc/issues/7)

## Minor changes

- `epmc_search()` and other functions return non-nested data.frames as tibbles to better support the tidyverse
- `epmc_search()` improve error handling when nothing was found
- `epmc_details()` [added MeSH qualifer #8]((https://github.com/ropensci/europepmc/issues/8)
- remove NBK` as data source for `epmc_details()`, use PMIDs (`MED`) instead
- fix warnings regarding vignettes and imported dependencies reported by CRAN

# europepmc 0.1

Initial submission to CRAN

## Major changes

Support of the following  Europe PMC RESTful API methods:

- search
- citations
- references
- databaseLinks
- labsLinks
- textMinedTerms
- fullTextXML
- bookXML

Changes made during the ropensci onboarding review by @toph-allen <https://github.com/ropensci/software-review/issues/29>

Answering to @cstubben reports and suggestions:

- [search returns data.frame of lists #1](https://github.com/ropensci/europepmc/issues/1)
- [removed hit_count from epmc_search results #4](https://github.com/ropensci/europepmc/issues/4)
- [set default batch_size to 1000 and limit to 25? #5](https://github.com/ropensci/europepmc/issues/4)

# europepmc v 0.4.1

## Test environments

* local OS X install (Platform: aarch64-apple-darwin20 (64-bit)), R version 4.1.0 (2021-05-18)
* GitHub Actions initialized with `usethis::use_github_action("check-standard")`: windows-latest (release), macOS-latest (release), ubuntu-20.04 (release, devel)
* Win-Builder 

## R CMD check results

On local machine (OS):

Status: OK

Win-Builder:

Status: OK


## Reverse dependencies

* I have run R CMD check on downstream dependencies using revdepcheck::revdep_check() and found no problems related to this new version.

---

This submission implements the most recent API changes. This will fix current R checks problems. It will also prevent warnings and errors on CRAN when the internet resource should fail.

Thanks!

Najko Jahn
europepmc - R Interface to Europe PMC RESTful Web Service
=== 


```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```


[![R build status](https://github.com/ropensci/europepmc/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/europepmc/actions)
[![Build status](https://ci.appveyor.com/api/projects/status/f8xtpvhhr074lk44?svg=true)](https://ci.appveyor.com/project/sckott/europepmc)
[![codecov.io](https://codecov.io/github/ropensci/europepmc/coverage.svg?branch=master)](https://codecov.io/github/ropensci/europepmc?branch=master)
[![cran version](https://www.r-pkg.org/badges/version/europepmc)](https://cran.r-project.org/package=europepmc)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/europepmc)](https://github.com/r-hub/cranlogs.app)
[![](https://badges.ropensci.org/29_status.svg)](https://github.com/ropensci/software-review/issues/29)

europepmc facilitates access to the [Europe PMC RESTful Web
Service](https://europepmc.org/RestfulWebService). The client furthermore supports the [Europe PMC Annotations API](https://europepmc.org/AnnotationsApi) to retrieve text-mined concepts and terms per article.

[Europe PMC](https://europepmc.org/) covers life science literature and
gives access to open access full texts. Europe
PMC ingests all PubMed content and extends its index with other literature and patent sources.

For more infos on Europe PMC, see:

<https://europepmc.org/About>

Levchenko, M., Gou, Y., Graef, F., Hamelers, A., Huang, Z., Ide-Smith, M., … McEntyre, J. (2017). Europe PMC in 2017. Nucleic Acids Research, 46(D1), D1254–D1260. <https://doi.org/10.1093/nar/gkx1005>

## Implemented API methods

This client supports the following API methods from the [Articles RESTful API](https://europepmc.org/RestfulWebService):

|API-Method     |Description                                                                                  |R functions                                |
|:--------------|:--------------------------------------------------------------------------------------------|:------------------------------------------|
|search         |Search Europe PMC and get detailed metadata                                                  |`epmc_search()`, `epmc_details()`, `epmc_search_by_doi()`          |
|profile        |Obtain a summary of hit counts for several Europe PMC databases                              |`epmc_profile()`                           |
|citations      |Load metadata representing citing articles for a given publication                           |`epmc_citations()`                         |
|references     |Retrieve the reference section of a publication                                               |`epmc_refs()`                              |
|databaseLinks  |Get links to biological databases such as UniProt or ENA                                     |`epmc_db()`, `epmc_db_count()`             |
|labslinks      |Access links to Europe PMC provided by third parties                                         |`epmc_lablinks()`, `epmc_lablinks_count()` |
|fullTextXML    |Fetch full-texts deposited in PMC                                                            |`epmc_ftxt()`                              |
|bookXML        |retrieve book XML formatted full text for the Open Access subset of the Europe PMC bookshelf |`epmc_ftxt_book()`                         |

From the [Europe PMC Annotations API](https://europepmc.org/AnnotationsApi):

|API-Method     |Description |R functions |
|:-----------|:-------------|:-------------|
annotationsByArticleIds | Get the annotations contained in the list of articles specified | `epmc_annotations_by_id()` |

## Installation

From CRAN

```r
install.packages("europepmc")
```

The latest development version can be installed using the
[remotes](https://github.com/r-lib/remotes/) package:


```r
require(remotes)
install_github("ropensci/europepmc")
```

Loading into R

```{r}
library(europepmc)
```

## Search Europe PMC

The search covers both metadata (e.g. abstracts or title) and full texts. To
build your query, please refer to the comprehensive guidance on how to search
Europe PMC: <https://europepmc.org/help>. Provide your query in the Europe
PMC search syntax to `epmc_search()`. 

```{r}
europepmc::epmc_search(query = '"2019-nCoV" OR "2019nCoV"')
```

Be aware that Europe PMC expands queries with MeSH synonyms by default. You can turn this behavior off using the `synonym = FALSE` parameter.

By default, `epmc_search()` returns 100 records. To adjust the limit, simply use
the `limit` parameter.

See vignette [Introducing europepmc, an R interface to Europe PMC RESTful API](https://docs.ropensci.org/europepmc/articles/introducing-europepmc.html) for a long-form documentation about how to search Europe PMC with this client.

## Creating proper review graphs with `epmc_hits_trend()`

There is also a nice function allowing you to easily create review graphs like described in Maëlle
Salmon's [blog post](https://masalmon.eu/2017/05/14/evergreenreviewgraph/):

```{r, fig.path="man/figures/"}
tt_oa <- europepmc::epmc_hits_trend("Malaria", period = 1995:2019, synonym = FALSE)
tt_oa
# we use ggplot2 for plotting the graph
library(ggplot2)
ggplot(tt_oa, aes(year, query_hits / all_hits)) + 
  geom_point() + 
  geom_line() +
  xlab("Year published") + 
  ylab("Proportion of articles on Malaria in Europe PMC")
```

For more info, read the vignette about creating literature review graphs:

<https://docs.ropensci.org/europepmc/articles/evergreenreviewgraphs.html>

## Re-use of europepmc

Check out the tidypmc package

<https://github.com/ropensci/tidypmc>

The package maintainer, Chris Stubben (@cstubben), has also created an Shiny App that allows you to search and browse Europe PMC:

<https://github.com/cstubben/euPMC>



## Other ways to access Europe PubMed Central

### Other APIs

- Data dumps: <https://europepmc.org/FtpSite>
- OAI service: <https://europepmc.org/OaiService>
- SOAP web service: <https://europepmc.org/SoapWebServices>
- Grants RESTful (Grist) API: <https://europepmc.org/GristAPI>

### Other R clients

- use rOpenSci's `oai` to get metadata and full text via Europe PMC's OAI interface: <https://github.com/ropensci/oai>
- use rOpenSci's `rentrez` to interact with [NCBI databases](https://www.ncbi.nlm.nih.gov/) such as PubMed: <https://github.com/ropensci/rentrez>
- rOpenSci's `fulltext` package gives access to supplementary material of open access life-science publications in Europe PMC: <https://github.com/ropensci/fulltext>

## Meta

Please note that this project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/europepmc/CONDUCT.html). By participating in this project you agree to abide by its terms.

License: GPL-3

Please use the issue tracker for bug reporting and feature requests.

---

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "Making trend graphs"
author: "Najko Jahn"
date: "2021-08-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Making trend graphs}
  \usepackage[utf8]{inputenc}
---



Trend graphs in literature reviews show the development of concepts in scholarly communication. Some trend graphs, however, don't acknowledge that the number of scholarly publications is growing each year, but simply display the absolute number of hits they have found for a given concept. Noam Ross called these misleading graphs evergreen review graphs because of their enduring popularity in  review papers. Examples can be found on Twitter under the Hashtag [#evergreenreviewgraph](https://twitter.com/hashtag/evergreenreviewgraph).

This vignette guides you how to make proper trend graphs when reviewing Europe PMC literature. In these graphs, the number of hits found is divided by the total number of records indexed in Europe PMC for a given search query. 

## Preparing proper review graphs with `epmc_hits_trend()`

We use `epmc_hits_trend()` function, which was firstly introduced in Maëlle Salmon's blog post about "How not to make an evergreen review graph"[^1]. The function takes a query in the Europe PMC search syntax[^2] and the period of years over which to perform the search as arguments, and returns a data-frame with year, total number of hits (`all_hits`) and number of hits for the query (`query_hits`).   



```r
library(europepmc)
europepmc::epmc_hits_trend(query = "aspirin", period = 2010:2016)
#> # A tibble: 7 × 3
#>    year all_hits query_hits
#>   <int>    <dbl>      <dbl>
#> 1  2010   851021       7219
#> 2  2011   904618       7888
#> 3  2012   945899       9054
#> 4  2013  1003845      10127
#> 5  2014  1055435      10895
#> 6  2015  1095859      11750
#> 7  2016  1116157      12099
```

By default, synonym search is disabled and only Medline/PubMed index is searched.

[^1]: <https://masalmon.eu/2017/05/14/evergreenreviewgraph/>

[^2]: Europe PMC Search Syntax: <https://europepmc.org/Help#mostofsearch>

## Use Cases 

### Use Case: Growth of Open Access Literature

There is a growing interest in knowing the proportion of open access to scholarly literature. Europe PMC allows searching for open access content with the [`OPEN_ACCESS:Y` parameter](https://europepmc.org/search?query=OPEN_ACCESS:Y&page=1&sortby=Relevance). At the moment, Europe PMC contains 3,740,002 open access full-texts. Let's see how they are relatively distributed over the period 1995 - 2016.


```r
tt_oa <- europepmc::epmc_hits_trend("OPEN_ACCESS:Y", period = 1995:2016, synonym = FALSE)
tt_oa
#> # A tibble: 22 × 3
#>     year all_hits query_hits
#>    <int>    <dbl>      <dbl>
#>  1  1995   449064       3337
#>  2  1996   458526       3508
#>  3  1997   456744       3665
#>  4  1998   474613       3836
#>  5  1999   493745       3918
#>  6  2000   532019       4328
#>  7  2001   545674       5479
#>  8  2002   561426       5894
#>  9  2003   588572       7148
#> 10  2004   628141       9795
#> # … with 12 more rows
# we use ggplot2 for plotting the graph
library(ggplot2)
ggplot(tt_oa, aes(year, query_hits / all_hits)) + 
  geom_point() + 
  geom_line() +
  xlab("Year published") + 
  ylab("Proportion of OA full-texts in Europe PMC")
```

<img src="../vignettes/oa_pmc-1.png" title="oa in europe pmc" alt="oa in europe pmc" style="display: block; margin: auto;" />

Be careful with the interpretation of the slower growth in the last years because there are several ways how open access content is added to Europe PMC including the digitalization of back issues.[^3]

[^3]: See section "Content Growth" in: McEntyre JR, Ananiadou S, Andrews S, et al. UKPMC: a full text article resource
    for the life sciences. *Nucleic Acids Research*. 2011;39(Database):D58–D65. <https://doi.org/10.1093/nar/gkq1063>.

### Use Case: Cited open source software in scholarly publications

Another nice use case for trend graphs is to study how code and software repositories are cited in scientific literature. In recent years, it has become a good practice not only to re-use openly available software, but also to cite them. The FORCE11 Software Citation Working Group states:

> In general, we believe that software should be cited on the same basis as any other research product such as a paper or book; that is, authors should cite the appropriate set of software products just as they cite the appropriate set of papers. [(doi:10.7717/peerj-cs.86)](https://doi.org/10.7717/peerj-cs.86)

So let's see whether we can find evidence for this evolving practice by creating a proper review graph. As a start, we examine these four general purpose hosting services for version-controlled code:

- [code.google.com](https://code.google.com/)
- [github.com](https://github.com/)
- [sourceforge.net](https://sourceforge.net/)
- [bitbucket.org](https://bitbucket.org/)

and, of course, [CRAN](https://cran.r-project.org/), the R archive network.

#### How to query Europe PMC?

We only want to search reference lists. Because Europe PMC does not index references for its complete collection, we use `has_reflist:y` to restrict our search to those publications with reference lists. These literature sections can be searched with the `REF:` parameter.

Let's prepare the queries for links to the above mentioned code hosting services:


```r
dvcs <- c("code.google.com", "github.com", 
          "sourceforge.net", "bitbucket.org", "cran.r-project.org")
# make queries including reference section
dvcs_query <- paste0('REF:"', dvcs, '"')
```

and get publications for which Europe PMC gives access to reference lists for normalizing the review graph.


```r
library(dplyr)
my_df <- purrr::map_df(dvcs_query, function(x) {
  # get number of publications with indexed reference lists
  refs_hits <- 
    europepmc::epmc_hits_trend("has_reflist:y", period = 2009:2016, synonym = FALSE)$query_hits
  # get hit count querying for code repositories 
  europepmc::epmc_hits_trend(x, period = 2009:2016, synonym = FALSE) %>% 
    dplyr::mutate(query_id = x) %>%
    dplyr::mutate(refs_hits = refs_hits) %>%
    dplyr::select(year, all_hits, refs_hits, query_hits, query_id)
}) 
my_df
#> # A tibble: 40 × 5
#>     year all_hits refs_hits query_hits query_id                 
#>    <int>    <dbl>     <dbl>      <dbl> <chr>                    
#>  1  2009   793068    555477         13 "REF:\"code.google.com\""
#>  2  2010   851021    540514         40 "REF:\"code.google.com\""
#>  3  2011   904618    603060         65 "REF:\"code.google.com\""
#>  4  2012   945899    635448         92 "REF:\"code.google.com\""
#>  5  2013  1003845    761512        135 "REF:\"code.google.com\""
#>  6  2014  1055435    797039        140 "REF:\"code.google.com\""
#>  7  2015  1095859    779562        117 "REF:\"code.google.com\""
#>  8  2016  1116157    782630         65 "REF:\"code.google.com\""
#>  9  2009   793068    555477          2 "REF:\"github.com\""     
#> 10  2010   851021    540514         10 "REF:\"github.com\""     
#> # … with 30 more rows

### total
hits_summary <- my_df %>% 
  group_by(query_id) %>% 
  summarise(all = sum(query_hits)) %>% 
  arrange(desc(all))
hits_summary
#> # A tibble: 5 × 2
#>   query_id                       all
#>   <chr>                        <dbl>
#> 1 "REF:\"cran.r-project.org\""  8221
#> 2 "REF:\"github.com\""          1609
#> 3 "REF:\"code.google.com\""      667
#> 4 "REF:\"sourceforge.net\""      643
#> 5 "REF:\"bitbucket.org\""         94
```

The proportion of papers where Europe PMC was able to make the cited literature available was 70 for the period 2009-2016. There also seems to be a time-lag between indexing reference lists because the absolute number of publication was decreasing over the years. This is presumably because Europe PMC also includes delayed open access content, i.e. content which is not added immediately with the original publication.[^4]

[^4]: Ebd.

Now, let's make a proper review graph normalizing our query results with the number of publications with indexed references.


```r
library(ggplot2)
ggplot(my_df, aes(factor(year), query_hits / refs_hits, group = query_id, 
                  color = query_id)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 2) +
  scale_color_brewer(name = "Query", palette = "Set1")+
  xlab("Year published") +
  ylab("Proportion of articles in Europe PMC")
```

<img src="../vignettes/software_lit-1.png" title="literature links to software in europe pmc" alt="literature links to software in europe pmc" style="display: block; margin: auto;" />

#### Discussion and Conclusion

Although this figure illustrates the relative popularity of citing code hosted by CRAN and GitHub in recent years, there are some limits that needs to be discussed. As said before, Europe PMC does not extract reference lists from every indexed publication. It furthermore remains open whether and to what extent software is cited outside the reference section, i.e. as footnote or in the acknowledgements. 

Another problem of our query approach is that we did not consider that DOIs can also be used to cite software, a best-practice implemented by [Zenodo and GitHub](https://guides.github.com/activities/citable-code/) or the [The Journal of Open Source Software](https://joss.theoj.org/). 

Lastly, it actually remains unclear, which and what kind of software is cited how often. We could also not control if authors just cited the homepages and not a particular source code repository. One paper can also cite more than one code repository, which is also not represented in the trend graph. 

To conclude, a proper trend graph on the extent of software citation can only be the start for a more sophisticated approach that mines links to software repositories from scientific literature and fetches metadata about these code repositories from the hosting facilities.

## Conclusion

This vignette presented first steps on how to make trend graphs with `europepmc`. As our use-cases suggest, please carefully consider how you queried Europe PMC in the interpretation of your graph. Although trend graphs are a nice way to illustrate the development of certain concepts in scientific literature or recent trends in scholarly communication, they must be put in context in order to become meaningful. 

## Acknowledgements

Big thanks to Maëlle Salmon for getting me started to write this vignette.
---
title: "Overview"
author: "Najko Jahn"
date: "2021-08-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Overview}
  \usepackage[utf8]{inputenc}
---




## What is searched?

[Europe PMC](https://europepmc.org/) is a repository of life science literature. Europe PMC ingests all PubMed content and extends its index with other literature and patent sources. 

For more background on Europe PMC, see:

<https://europepmc.org/About>

Levchenko, M., Gou, Y., Graef, F., Hamelers, A., Huang, Z., Ide-Smith, M., … McEntyre, J. (2017). Europe PMC in 2017. Nucleic Acids Research, 46(D1), D1254–D1260. <https://doi.org/10.1093/nar/gkx1005>

## How to search Europe PMC with R?

This client supports the [Europe PMC search syntax](https://europepmc.org/Help#SSR). If you are unfamiliar with searching Europe PMC, check out the [Europe PMC query builder](https://europepmc.org/advancesearch), a very nice tool that helps you to build queries. To make use of Europe PMC queries in R, copy & paste the search string to the search functions of this package. 

In the following, some examples demonstrate how to search Europe PMC with R.

### Free search

`empc_search()` is the main function to query Europe PMC. It searches both metadata and fulltexts. 



```r
library(europepmc)
europepmc::epmc_search('malaria')
#> # A tibble: 100 × 29
#>    id       source pmid     doi   title authorString journalTitle issue journalVolume
#>    <chr>    <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr> <chr>        
#>  1 34100426 MED    34100426 10.4… New … Lima MN, Ba… Neural Rege… 1     17           
#>  2 33341138 MED    33341138 10.1… Trip… Wang J, Xu … Lancet       10267 396          
#>  3 33341139 MED    33341139 10.1… Trip… van der Plu… Lancet       10267 396          
#>  4 33535760 MED    33535760 10.3… THE … Damiani E, … Acta Med Hi… 2     18           
#>  5 33530764 MED    33530764 10.1… Disc… Hoarau M, V… J Enzyme In… 1     36           
#>  6 33372863 MED    33372863 10.1… ATP2… Lamy A, Mac… Emerg Micro… 1     10           
#>  7 33594960 MED    33594960 10.1… Mana… Kambale-Kom… Hematology   1     26           
#>  8 34283002 MED    34283002 10.1… <i>P… Alhassan AM… Pharm Biol   1     59           
#>  9 34184352 MED    34184352 10.1… Stru… Chhibber-Go… Protein Sci  9     30           
#> 10 34419123 MED    34419123 10.1… Burd… Dao F, Djon… Parasit Vec… 1     14           
#> # … with 90 more rows, and 20 more variables: pubYear <chr>, journalIssn <chr>,
#> #   pageInfo <chr>, pubType <chr>, isOpenAccess <chr>, inEPMC <chr>,
#> #   inPMC <chr>, hasPDF <chr>, hasBook <chr>, hasSuppl <chr>,
#> #   citedByCount <int>, hasReferences <chr>, hasTextMinedTerms <chr>,
#> #   hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, pmcid <chr>, versionNumber <int>
```

It is worth noting that Europe PMC expands queries with MeSH synonyms by default, a behavior which can be turned off with the `synonym` parameter. 


```r
europepmc::epmc_search('malaria', synonym = FALSE)
#> # A tibble: 100 × 29
#>    id        source pmid     doi   title authorString journalTitle issue journalVolume
#>    <chr>     <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr> <chr>        
#>  1 33341139  MED    33341139 10.1… Trip… van der Plu… Lancet       10267 396          
#>  2 33341138  MED    33341138 10.1… Trip… Wang J, Xu … Lancet       10267 396          
#>  3 34100426  MED    34100426 10.4… New … Lima MN, Ba… Neural Rege… 1     17           
#>  4 34184352  MED    34184352 10.1… Stru… Chhibber-Go… Protein Sci  9     30           
#>  5 34380494  MED    34380494 10.1… Publ… Heuschen AK… Malar J      1     20           
#>  6 33530764  MED    33530764 10.1… Disc… Hoarau M, V… J Enzyme In… 1     36           
#>  7 34399767  MED    34399767 10.1… Inve… Njau J, Sil… Malar J      1     20           
#>  8 PPR385006 PPR    <NA>     10.2… Temp… Ingholt MM,… <NA>         <NA>  <NA>         
#>  9 34419123  MED    34419123 10.1… Burd… Dao F, Djon… Parasit Vec… 1     14           
#> 10 34376219  MED    34376219 10.1… An a… Wanzira H, … BMC Health … 1     21           
#> # … with 90 more rows, and 20 more variables: pubYear <chr>, journalIssn <chr>,
#> #   pageInfo <chr>, pubType <chr>, isOpenAccess <chr>, inEPMC <chr>,
#> #   inPMC <chr>, hasPDF <chr>, hasBook <chr>, hasSuppl <chr>,
#> #   citedByCount <int>, hasReferences <chr>, hasTextMinedTerms <chr>,
#> #   hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, pmcid <chr>, versionNumber <int>
```

To get an exact match, use quotes as in the following example:


```r
europepmc::epmc_search('"Human malaria parasites"')
#> # A tibble: 100 × 29
#>    id        source pmid     doi   title authorString journalTitle pubYear journalIssn
#>    <chr>     <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr>   <chr>      
#>  1 34415329  MED    34415329 10.1… Func… Kimata-Arig… J Biochem    2021    "0021-924x…
#>  2 34087264  MED    34087264 10.1… Dive… Goh XT, Lim… Mol Biochem… 2021    "0166-6851…
#>  3 34400833  MED    34400833 10.1… A he… Tintó-Font … Nat Microbi… 2021    "2058-5276"
#>  4 33789941  MED    33789941 10.1… Addi… Kwon H, Sim… mSphere      2021    "2379-5042"
#>  5 34211355  MED    34211355 <NA>  An E… Clark NF, T… Yale J Biol… 2021    "0044-0086…
#>  6 34362867  MED    34362867 10.4… High… Lai MY, Raf… Trop Biomed  2021    "0127-5720…
#>  7 33693917  MED    33693917 10.1… Non-… Antinori S,… J Travel Med 2021    "1195-1982…
#>  8 32470136  MED    32470136 10.1… C-te… Kimata-Arig… J Biochem    2020    "0021-924x…
#>  9 PPR353209 PPR    <NA>     10.1… 5-me… Liu M, Guo … <NA>         2021     <NA>      
#> 10 33797521  MED    33797521 10.4… Comp… Mat Salleh … Trop Biomed  2021    "0127-5720…
#> # … with 90 more rows, and 20 more variables: pubType <chr>,
#> #   isOpenAccess <chr>, inEPMC <chr>, inPMC <chr>, hasPDF <chr>, hasBook <chr>,
#> #   hasSuppl <chr>, citedByCount <int>, hasReferences <chr>,
#> #   hasTextMinedTerms <chr>, hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, journalVolume <chr>, pageInfo <chr>,
#> #   issue <chr>, pmcid <chr>, versionNumber <int>
```

### Managing search results

By default, 100 records are returned, but the number of results can be expanded or limited with the `limit` parameter. 



```r
europepmc::epmc_search('"Human malaria parasites"', limit = 10)
#> # A tibble: 10 × 28
#>    id        source pmid     doi   title authorString journalTitle pubYear journalIssn
#>    <chr>     <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr>   <chr>      
#>  1 34415329  MED    34415329 10.1… Func… Kimata-Arig… J Biochem    2021    "0021-924x…
#>  2 34087264  MED    34087264 10.1… Dive… Goh XT, Lim… Mol Biochem… 2021    "0166-6851…
#>  3 34400833  MED    34400833 10.1… A he… Tintó-Font … Nat Microbi… 2021    "2058-5276"
#>  4 33789941  MED    33789941 10.1… Addi… Kwon H, Sim… mSphere      2021    "2379-5042"
#>  5 34211355  MED    34211355 <NA>  An E… Clark NF, T… Yale J Biol… 2021    "0044-0086…
#>  6 34362867  MED    34362867 10.4… High… Lai MY, Raf… Trop Biomed  2021    "0127-5720…
#>  7 33693917  MED    33693917 10.1… Non-… Antinori S,… J Travel Med 2021    "1195-1982…
#>  8 32470136  MED    32470136 10.1… C-te… Kimata-Arig… J Biochem    2020    "0021-924x…
#>  9 PPR353209 PPR    <NA>     10.1… 5-me… Liu M, Guo … <NA>         2021     <NA>      
#> 10 33797521  MED    33797521 10.4… Comp… Mat Salleh … Trop Biomed  2021    "0127-5720…
#> # … with 19 more variables: pubType <chr>, isOpenAccess <chr>, inEPMC <chr>,
#> #   inPMC <chr>, hasPDF <chr>, hasBook <chr>, hasSuppl <chr>,
#> #   citedByCount <int>, hasReferences <chr>, hasTextMinedTerms <chr>,
#> #   hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, journalVolume <chr>, pageInfo <chr>,
#> #   issue <chr>, pmcid <chr>
```

Results are sorted by relevance. Other options via the `sort` parameter are 

- `sort = 'cited'` by the number of citation, descending from the most cited publication
- `sort = 'date'` by date published starting with the most recent publication

### Search by DOIs

Sometimes, you would like to check, if articles are indexed in Europe PMC using DOI names, a widely used identifier for scholarly articles. Use `epmc_search_by_doi()` for this purpose.


```r
my_dois <- c(
  "10.1159/000479962",
  "10.1002/sctm.17-0081",
  "10.1161/strokeaha.117.018077",
  "10.1007/s12017-017-8447-9"
  )
europepmc::epmc_search_by_doi(doi = my_dois)
#> # A tibble: 4 × 28
#>   id       source pmid     doi   title authorString journalTitle issue journalVolume
#>   <chr>    <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr> <chr>        
#> 1 28957815 MED    28957815 10.1… Clin… Schnieder M… Eur Neurol   5-6   78           
#> 2 28941317 MED    28941317 10.1… Conc… Doeppner TR… Stem Cells … 11    6            
#> 3 29018132 MED    29018132 10.1… One-… Psychogios … Stroke       11    48           
#> 4 28623611 MED    28623611 10.1… Defe… Carboni E, … Neuromolecu… 2-3   19           
#> # … with 19 more variables: pubYear <chr>, journalIssn <chr>, pageInfo <chr>,
#> #   pubType <chr>, isOpenAccess <chr>, inEPMC <chr>, inPMC <chr>, hasPDF <chr>,
#> #   hasBook <chr>, hasSuppl <chr>, citedByCount <int>, hasReferences <chr>,
#> #   hasTextMinedTerms <chr>, hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, pmcid <chr>
```

### Output options

By default, a non-nested data frame printed as tibble is returned. 
Other formats are `output = "id_list"` returning a list of IDs and sources, 
and output = "'raw'"" for getting full metadata as list. 
Please be aware that these lists can become very large.

### More advanced options to search Europe PMC

#### Author search

Use the Europe PMC query syntax to search by author names:


```r
europepmc::epmc_search('AUTH:"Salmon Maelle"')
#> # A tibble: 10 × 28
#>    id       source pmid     doi   title authorString journalTitle issue journalVolume
#>    <chr>    <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr> <chr>        
#>  1 30378432 MED    30378432 10.1… When… Milà C, Sal… Environ Sci… 22    52           
#>  2 29778830 MED    29778830 10.1… Wear… Salmon M, M… Environ Int  <NA>  117          
#>  3 29751338 MED    29751338 10.1… Use … Kumar MK, S… Environ Pol… <NA>  239          
#>  4 29330030 MED    29330030 10.1… Heal… Mueller N, … Prev Med     <NA>  109          
#>  5 29626773 MED    29626773 10.1… Deve… Sanchez M, … Sci Total E… <NA>  634          
#>  6 29088243 MED    29088243 10.1… Time… Schumacher … PLoS One     10    12           
#>  7 28606699 MED    28606699 10.1… Inte… Tonne C, Sa… Int J Hyg E… 6     220          
#>  8 28708095 MED    28708095 10.3… Pred… Sanchez M, … Int J Envir… 7     14           
#>  9 27063588 MED    27063588 10.2… A sy… Salmon M, S… Euro Survei… 13    21           
#> 10 26250543 MED    26250543 10.1… Baye… Salmon M, S… Biom J       6     57           
#> # … with 19 more variables: pubYear <chr>, journalIssn <chr>, pageInfo <chr>,
#> #   pubType <chr>, isOpenAccess <chr>, inEPMC <chr>, inPMC <chr>, hasPDF <chr>,
#> #   hasBook <chr>, hasSuppl <chr>, citedByCount <int>, hasReferences <chr>,
#> #   hasTextMinedTerms <chr>, hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, pmcid <chr>
```

[Europe PMC Advanced Search](https://europepmc.org/advancesearch) has a auto-suggest field for author names if you feel unsure how the name you are searching for is indexed in Europe PMC. Using the Boolean `OR` operator allows searching for more than one spelling variant:


```r
q <- 'AUTH:"PÜHLER Alfred" OR AUTH:"Pühler Alfred Prof. Dr." OR AUTH:"Puhler A"'
europepmc::epmc_search(q, limit = 1000)
#> # A tibble: 590 × 29
#>    id        source pmid     pmcid doi   title authorString journalTitle journalVolume
#>    <chr>     <chr>  <chr>    <chr> <chr> <chr> <chr>        <chr>        <chr>        
#>  1 34367203  MED    34367203 PMC8… 10.3… ExoS… Geiger O, S… Front Plant… 12           
#>  2 34361893  MED    34361893 PMC8… 10.3… Indi… Hassa J, Kl… Microorgani… 9            
#>  3 34040261  MED    34040261 PMC8… 10.1… Swar… Warnat-Herr… Nature       594          
#>  4 33589928  MED    33589928 <NA>  10.1… Impl… Mayer G, Mü… Brief Bioin… <NA>         
#>  5 33643369  MED    33643369 PMC7… 10.3… Exop… Castellani … Front Plant… 12           
#>  6 33441124  MED    33441124 PMC7… 10.1… Dise… Aschenbrenn… Genome Med   13           
#>  7 PPR264825 PPR    <NA>     <NA>  10.2… The … Droste J, O… <NA>         <NA>         
#>  8 33220679  MED    33220679 <NA>  10.1… Glob… Nilsson JF,… FEMS Microb… 97           
#>  9 33348776  MED    33348776 PMC7… 10.3… The … Maus I, Tub… Microorgani… 8            
#> 10 33296687  MED    33296687 PMC7… 10.1… Long… Bernardes J… Immunity     53           
#> # … with 580 more rows, and 20 more variables: pubYear <chr>,
#> #   journalIssn <chr>, pageInfo <chr>, pubType <chr>, isOpenAccess <chr>,
#> #   inEPMC <chr>, inPMC <chr>, hasPDF <chr>, hasBook <chr>, hasSuppl <chr>,
#> #   citedByCount <int>, hasReferences <chr>, hasTextMinedTerms <chr>,
#> #   hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, issue <chr>, versionNumber <int>
```

There is a considerable overlap between common names. The integration of ORCID, a persistent author identifier, allows unambiguous search for personal publications in Europe PMC. For example, here's how to search for publications written by Bernd Weisshaar (ORCID: <https://orcid.org/0000-0002-7635-3473>) sorted by the number of times cited in descending order:


```r
europepmc::epmc_search('AUTHORID:"0000-0002-7635-3473"', limit = 200, sort = "cited")
#> # A tibble: 150 × 28
#>    id       source pmid     doi   title authorString journalTitle issue journalVolume
#>    <chr>    <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr> <chr>        
#>  1 21873998 MED    21873998 10.1… The … Wang X, Wan… Nat Genet    10    43           
#>  2 20674465 MED    20674465 10.1… MYB … Dubos C, St… Trends Plan… 10    15           
#>  3 11597504 MED    11597504 10.1… The … Stracke R, … Curr Opin P… 5     4            
#>  4 11906833 MED    11906833 10.1… bZIP… Jakoby M, W… Trends Plan… 3     7            
#>  5 14756321 MED    14756321 10.1… An A… Rosso MG, L… Plant Mol B… 1-2   53           
#>  6 12679534 MED    12679534 10.1… The … Heim MA, Ja… Mol Biol Ev… 5     20           
#>  7 11080161 MED    11080161 10.1… Tran… Jin H, Comi… EMBO J       22    19           
#>  8 15361138 MED    15361138 10.1… Comp… Zimmermann … Plant J      1     40           
#>  9 15255866 MED    15255866 10.1… TT2,… Baudry A, H… Plant J      3     39           
#> 10 17419845 MED    17419845 10.1… Diff… Stracke R, … Plant J      4     50           
#> # … with 140 more rows, and 19 more variables: pubYear <chr>,
#> #   journalIssn <chr>, pageInfo <chr>, pubType <chr>, isOpenAccess <chr>,
#> #   inEPMC <chr>, inPMC <chr>, hasPDF <chr>, hasBook <chr>, hasSuppl <chr>,
#> #   citedByCount <int>, hasReferences <chr>, hasTextMinedTerms <chr>,
#> #   hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, pmcid <chr>
```

#### Annotations 

Europe PMC provides text-mined annotations contained in abstracts and open access full-text articles.

These automatically identified concepts and term can be retrieved at the article-level:


```r
europepmc::epmc_annotations_by_id(c("MED:28585529", "PMC:PMC1664601"))
#> # A tibble: 774 × 13
#>    source ext_id   pmcid      prefix exact postfix name  uri   id    type  section
#>    <chr>  <chr>    <chr>      <chr>  <chr> <chr>   <chr> <chr> <chr> <chr> <chr>  
#>  1 MED    28585529 PMC5467160 "tive… Beta… " allo… Beta… http… http… Clin… Title …
#>  2 MED    28585529 PMC5467160 "nomi… genes ".\nRa… gene  http… http… Sequ… Title …
#>  3 MED    28585529 PMC5467160 "nomi… genes " is o… gene  http… http… Sequ… Abstra…
#>  4 MED    28585529 PMC5467160 " One… genes " are … gene  http… http… Sequ… Abstra…
#>  5 MED    28585529 PMC5467160 " ide… beet  " (Bet… Beta… http… http… Clin… Abstra…
#>  6 MED    28585529 PMC5467160 "ify … Beta… " ssp.… Beta… http… http… Clin… Abstra…
#>  7 MED    28585529 PMC5467160 "ulga… gene  " Rz2 … gene  http… http… Sequ… Abstra…
#>  8 MED    28585529 PMC5467160 "e ge… geno… " sequ… geno… http… http… Sequ… Abstra…
#>  9 MED    28585529 PMC5467160 "eque… beet  ". Our… Beta… http… http… Clin… Abstra…
#> 10 MED    28585529 PMC5467160 "disc… genes " rele… gene  http… http… Sequ… Abstra…
#> # … with 764 more rows, and 2 more variables: provider <chr>, subType <chr>
```

To obtain a list of articles where Europe PMC has text-minded annotations, either subset the resulting data.frame 


```r
tt <- epmc_search("malaria")
tt[tt$hasTextMinedTerms == "Y" | tt$hasTMAccessionNumbers == "Y",]
#> # A tibble: 94 × 29
#>    id        source pmid     doi   title authorString journalTitle issue journalVolume
#>    <chr>     <chr>  <chr>    <chr> <chr> <chr>        <chr>        <chr> <chr>        
#>  1 34100426  MED    34100426 10.4… New … Lima MN, Ba… Neural Rege… 1     17           
#>  2 33535760  MED    33535760 10.3… THE … Damiani E, … Acta Med Hi… 2     18           
#>  3 33530764  MED    33530764 10.1… Disc… Hoarau M, V… J Enzyme In… 1     36           
#>  4 33372863  MED    33372863 10.1… ATP2… Lamy A, Mac… Emerg Micro… 1     10           
#>  5 33594960  MED    33594960 10.1… Mana… Kambale-Kom… Hematology   1     26           
#>  6 34283002  MED    34283002 10.1… <i>P… Alhassan AM… Pharm Biol   1     59           
#>  7 34184352  MED    34184352 10.1… Stru… Chhibber-Go… Protein Sci  9     30           
#>  8 34362867  MED    34362867 10.4… High… Lai MY, Raf… Trop Biomed  3     38           
#>  9 34399767  MED    34399767 10.1… Inve… Njau J, Sil… Malar J      1     20           
#> 10 PPR385006 PPR    <NA>     10.2… Temp… Ingholt MM,… <NA>         <NA>  <NA>         
#> # … with 84 more rows, and 20 more variables: pubYear <chr>, journalIssn <chr>,
#> #   pageInfo <chr>, pubType <chr>, isOpenAccess <chr>, inEPMC <chr>,
#> #   inPMC <chr>, hasPDF <chr>, hasBook <chr>, hasSuppl <chr>,
#> #   citedByCount <int>, hasReferences <chr>, hasTextMinedTerms <chr>,
#> #   hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>, pmcid <chr>, versionNumber <int>
```

or expand the query choosing an annotation type or provider from the [Europe PMC Advanced Search](https://europepmc.org/advancesearch) query builder.


```r
epmc_search('malaria AND (ANNOTATION_TYPE:"Cell") AND (ANNOTATION_PROVIDER:"Europe PMC")')
#> # A tibble: 100 × 28
#>    id       source pmid     pmcid  doi   title  authorString  journalTitle issue
#>    <chr>    <chr>  <chr>    <chr>  <chr> <chr>  <chr>         <chr>        <chr>
#>  1 31782768 MED    31782768 PMC79… 10.1… Incre… Jongo SA, Ch… Clin Infect… 11   
#>  2 31808816 MED    31808816 PMC76… 10.1… Retin… Villaverde C… J Pediatric… 5    
#>  3 30989220 MED    30989220 PMC73… 10.1… Clini… Enane LA, Su… J Pediatric… 3    
#>  4 31300826 MED    31300826 PMC72… 10.1… Black… Opoka RO, Wa… Clin Infect… 11   
#>  5 31807752 MED    31807752 <NA>   10.1… Malar… Marcombe S, … J Med Entom… 3    
#>  6 31505001 MED    31505001 <NA>   10.1… Acute… Oshomah-Bell… J Trop Pedi… 2    
#>  7 31687768 MED    31687768 <NA>   10.1… Evalu… Ferdinand DY… Trans R Soc… 3    
#>  8 31693130 MED    31693130 PMC71… 10.1… Reduc… Kingston HWF… J Infect Dis 9    
#>  9 31679146 MED    31679146 <NA>   10.1… A Sys… Thiengsusuk … Eur J Drug … 2    
#> 10 30852586 MED    30852586 <NA>   10.1… An Ex… Woodford J, … J Infect Dis 6    
#> # … with 90 more rows, and 19 more variables: journalVolume <chr>,
#> #   pubYear <chr>, journalIssn <chr>, pageInfo <chr>, pubType <chr>,
#> #   isOpenAccess <chr>, inEPMC <chr>, inPMC <chr>, hasPDF <chr>, hasBook <chr>,
#> #   hasSuppl <chr>, citedByCount <int>, hasReferences <chr>,
#> #   hasTextMinedTerms <chr>, hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>
```

#### Data integrations

Another nice feature of Europe PMC is to search for cross-references between Europe PMC to other databases. For instance, to get publications cited by
entries in the [Protein Data bank in Europe](https://www.ebi.ac.uk/pdbe/node/1) published 2016:


```r
europepmc::epmc_search('(HAS_PDB:y) AND FIRST_PDATE:2016')
#> # A tibble: 100 × 28
#>    id       source pmid     pmcid  doi   title  authorString  journalTitle issue
#>    <chr>    <chr>  <chr>    <chr>  <chr> <chr>  <chr>         <chr>        <chr>
#>  1 27989121 MED    27989121 PMC58… 10.1… Short… Lin J, Pozha… Biochemistry 2    
#>  2 27815281 MED    27815281 PMC52… 10.1… Struc… Wakamatsu T,… Appl Enviro… 2    
#>  3 28035004 MED    28035004 PMC53… 10.1… Struc… Waz S, Nakam… J Biol Chem  7    
#>  4 28030602 MED    28030602 PMC51… 10.1… Struc… Christensen … PLoS One     12   
#>  5 28066558 MED    28066558 PMC51… 10.1… Struc… Gai Z, Wang … Cell Discov  <NA> 
#>  6 28024149 MED    28024149 PMC53… 10.1… Cryst… Kuk AC, Mash… Nat Struct … 2    
#>  7 28031486 MED    28031486 PMC52… 10.1… Struc… Sevrioukova … Proc Natl A… 3    
#>  8 28011634 MED    28011634 PMC53… 10.1… Struc… Levdikov VM,… J Biol Chem  7    
#>  9 28009010 MED    28009010 PMC51… 10.1… Struc… Zhao H, Wei … Sci Rep      <NA> 
#> 10 28197319 MED    28197319 PMC53… 10.1… Struc… Johannes JW,… ACS Med Che… 2    
#> # … with 90 more rows, and 19 more variables: journalVolume <chr>,
#> #   pubYear <chr>, journalIssn <chr>, pageInfo <chr>, pubType <chr>,
#> #   isOpenAccess <chr>, inEPMC <chr>, inPMC <chr>, hasPDF <chr>, hasBook <chr>,
#> #   hasSuppl <chr>, citedByCount <int>, hasReferences <chr>,
#> #   hasTextMinedTerms <chr>, hasDbCrossReferences <chr>, hasLabsLinks <chr>,
#> #   hasTMAccessionNumbers <chr>, firstIndexDate <chr>,
#> #   firstPublicationDate <chr>
```

The following sources are supported

- **CHEBI** a database and ontology of chemical entities of biological interest <https://www.ebi.ac.uk/chebi/>
- **CHEMBL** a database of bioactive drug-like small molecules <https://www.ebi.ac.uk/chembldb/>
- **EMBL** now ENA, provides a comprehensive record of the world's nucleotide sequencing information <https://www.ebi.ac.uk/ena/browser/>
- **INTACT** provides a freely available, open source database system and analysis tools for molecular interaction data <https://www.ebi.ac.uk/intact/>
- **INTERPRO** provides functional analysis of proteins by classifying them into families and predicting domains and important sites <https://www.ebi.ac.uk/interpro/>
- **OMIM** a comprehensive and authoritative compendium of human genes and genetic phenotypes <https://www.omim.org/about>
- **PDB** European resource for the collection, organisation and dissemination of data on biological macromolecular structures <https://www.ebi.ac.uk/pdbe/>
- **UNIPROT** comprehensive and freely accessible resource of protein sequence and functional information <https://www.uniprot.org/>
- **PRIDE** PRIDE Archive - proteomics data repository <https://www.ebi.ac.uk/pride/archive/>

To retrieve metadata about these external database links, use `europepmc_epmc_db()`. 

#### Citations and reference sections

Europe PMC let us also obtain citation metadata and reference sections. For retrieving citation metadata per article, use


```r
europepmc::epmc_citations("9338777", limit = 500)
#> # A tibble: 233 × 11
#>    id     source citationType title authorString journalAbbrevia… pubYear volume
#>    <chr>  <chr>  <chr>        <chr> <chr>        <chr>              <int> <chr> 
#>  1 33353… MED    review-arti… Xeno… Galow AM, G… Int J Mol Sci       2020 21    
#>  2 31565… MED    research-ar… Regu… Chung HC, N… J Vet Sci           2019 20    
#>  3 30230… MED    research su… Bioe… Legallais C… Adv Healthc Mat…    2018 7     
#>  4 30264… MED    research su… Porc… Fiebig U, F… Xenotransplanta…    2018 25    
#>  5 29756… MED    historical … Infe… Weiss RA.    Xenotransplanta…    2018 25    
#>  6 29642… MED    research su… Trac… Kawasaki J,… Viruses             2018 10    
#>  7 28768… MED    research su… Pres… Kawasaki J,… J Virol             2017 91    
#>  8 28437… MED    research su… Thre… Colon-Moran… Virology            2017 507   
#>  9 28054… MED    research su… Anti… Inoue Y, Yo… Ann Biomed Eng      2017 45    
#> 10 27832… MED    research-ar… Tran… Kim N, Choi… PLoS One            2016 11    
#> # … with 223 more rows, and 3 more variables: issue <chr>, citedByCount <int>,
#> #   pageInfo <chr>
```

For reference section from an article:


```r
europepmc::epmc_refs("28632490", limit = 200)
#> # A tibble: 169 × 19
#>    id       source citationType title authorString journalAbbrevia… issue pubYear
#>    <chr>    <chr>  <chr>        <chr> <chr>        <chr>            <chr>   <int>
#>  1 12002480 MED    JOURNAL ART… Tric… Adolfsson-E… Chemosphere      9-10     2002
#>  2 18795164 MED    JOURNAL ART… In v… Ahn KC, Zha… Environ Health … 9        2008
#>  3 18556606 MED    JOURNAL ART… Effe… Aiello AE, … Am J Public Hea… 8        2008
#>  4 17683018 MED    JOURNAL ART… Cons… Aiello AE, … Clin Infect Dis  <NA>     2007
#>  5 15273108 MED    JOURNAL ART… Rela… Aiello AE, … Antimicrob Agen… 8        2004
#>  6 18207219 MED    JOURNAL ART… The … Allmyr M, H… Sci Total Envir… 1        2008
#>  7 17007908 MED    JOURNAL ART… Tric… Allmyr M, A… Sci Total Envir… 1        2006
#>  8 26948762 MED    JOURNAL ART… Pres… Alvarez-Riv… J Chromatogr A   <NA>     2016
#>  9 23192912 MED    JOURNAL ART… Expo… Anderson SE… Toxicol Sci      1        2012
#> 10 25837385 MED    JOURNAL ART… Obse… Vladar EK, … Methods Cell Bi… <NA>     2015
#> # … with 159 more rows, and 11 more variables: volume <chr>, pageInfo <chr>,
#> #   citedOrder <int>, match <chr>, essn <chr>, issn <chr>,
#> #   publicationTitle <chr>, publisherLoc <chr>, publisherName <chr>,
#> #   externalLink <chr>, doi <chr>
```

#### Fulltext access

Europe PMC gives not only access to metadata, but also to full-texts. Adding `AND (OPEN_ACCESS:y)` to your search query, returns only those articles where Europe PMC has also the fulltext.

Fulltext as xml document can accessed via the PMID or the PubMed Central ID (PMCID):


```r
europepmc::epmc_ftxt("PMC3257301")
#> {xml_document}
#> <article article-type="research-article" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:mml="http://www.w3.org/1998/Math/MathML">
#> [1] <front>\n  <journal-meta>\n    <journal-id journal-id-type="nlm-ta">PLoS  ...
#> [2] <body>\n  <sec id="s1">\n    <title>Introduction</title>\n    <p>Atmosphe ...
#> [3] <back>\n  <ack>\n    <p>We would like to thank Dr. C. Gourlay and Dr. T.  ...
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_search_by_doi.r
\name{epmc_search_by_doi_}
\alias{epmc_search_by_doi_}
\title{Search Europe PMC by a DOI name}
\usage{
epmc_search_by_doi_(doi, .pb = NULL, output = NULL)
}
\arguments{
\item{doi}{character vector containing DOI names.}

\item{.pb}{progress bar object}

\item{output}{character, what kind of output should be returned. One of
'parsed', 'id_list' or 'raw' As default, parsed key metadata will be
returned as data.frame. 'id_list' returns a list of IDs and sources. Use
'raw' to get full metadata as list. Please be aware that these lists can
become very large.}
}
\description{
Please use \code{\link{epmc_search_by_doi}} instead. It calls this
method, returning open access status information from all your requests.
}
\examples{
\dontrun{
  epmc_search_by_doi_("10.1159/000479962")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_citations.r
\name{epmc_citations}
\alias{epmc_citations}
\title{Get citations for a given publication}
\usage{
epmc_citations(ext_id = NULL, data_src = "med", limit = 100, verbose = TRUE)
}
\arguments{
\item{ext_id}{character, publication identifier}

\item{data_src}{character, data source, by default Pubmed/MedLine index will
  be searched.


  The following three letter codes represent the sources
  Europe PubMed Central supports:
  \describe{
    \item{agr}{Agricola is a bibliographic database of citations to the
    agricultural literature created by the US National Agricultural Library
    and its co-operators.}
    \item{cba}{Chinese Biological Abstracts}
    \item{ctx}{CiteXplore}
    \item{eth}{EthOs Theses, i.e. PhD theses (British Library)}
    \item{hir}{NHS Evidence}
    \item{med}{PubMed/Medline NLM}
    \item{nbk}{Europe PMC Book metadata}
    \item{pat}{Biological Patents}
    \item{pmc}{PubMed Central}
    }}

\item{limit}{integer, number of results. By default, this function
returns 100 records.}

\item{verbose}{logical, print some information on what is going on.}
}
\value{
Metadata of citing documents as data.frame
}
\description{
Finds works that cite a given publication.
}
\examples{
\dontrun{
epmc_citations("PMC3166943", data_src = "pmc")
epmc_citations("9338777")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_hits.r
\name{epmc_hits}
\alias{epmc_hits}
\title{Get search result count}
\usage{
epmc_hits(query = NULL, ...)
}
\arguments{
\item{query}{query in the Europe PMC syntax}

\item{...}{add query parameters from `epmc_search()`, e.g. synonym=true}
}
\description{
Search over Europe PMC and retrieve the number of results found
}
\examples{
 \dontrun{
 epmc_hits('abstract:"burkholderia pseudomallei"')
 epmc_hits('AUTHORID:"0000-0002-7635-3473"')
 }
}
\seealso{
\code{\link{epmc_search}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_details.r
\name{epmc_details}
\alias{epmc_details}
\title{Get details for individual records}
\usage{
epmc_details(ext_id = NULL, data_src = "med")
}
\arguments{
\item{ext_id}{character, publication identifier}

\item{data_src}{character, data source, by default Pubmed/MedLine index will
be searched.
Other sources Europe PubMed Central supports are:
\describe{
  \item{agr}{Agricola is a bibliographic database of citations to the
  agricultural literature created by the US National Agricultural Library
  and its co-operators.}
  \item{cba}{Chinese Biological Abstracts}
  \item{ctx}{CiteXplore}
  \item{eth}{EthOs Theses, i.e. PhD theses (British Library)}
  \item{hir}{NHS Evidence}
  \item{med}{PubMed/Medline NLM}
  \item{pat}{Biological Patents}
  \item{pmc}{PubMed Central}
  \item{ppr}{Preprint records}
  }}
}
\value{
list of data frames
}
\description{
This function returns parsed metadata for a given publication ID
including abstract, full text links, author details including ORCID and affiliation,
MeSH terms, chemicals, grants.
}
\examples{
\dontrun{
epmc_details(ext_id = "26980001")
epmc_details(ext_id = "24270414")

# PMC record
epmc_details(ext_id = "PMC4747116", data_src = "pmc")

# Other sources:
# Agricolo
epmc_details("IND43783977", data_src = "agr")
# Biological Patents
epmc_details("EP2412369", data_src = "pat")
# Chinese Biological Abstracts
epmc_details("583843", data_src = "cba")
# CiteXplore
epmc_details("C6802", data_src = "ctx")
# NHS Evidence
epmc_details("338638", data_src = "hir")
# Theses
epmc_details("409323", data_src = "eth")
# Preprint
epmc_details("PPR158112", data_src = "ppr")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_search.r
\name{epmc_search_}
\alias{epmc_search_}
\title{Get one page of results when searching Europe PubMed Central}
\usage{
epmc_search_(
  query = NULL,
  limit = 100,
  output = "parsed",
  page_token = NULL,
  ...
)
}
\arguments{
\item{query}{character, search query. For more information on how to
build a search query, see \url{http://europepmc.org/Help}}

\item{limit}{integer, limit the number of records you wish to retrieve.
By default, 25 are returned.}

\item{output}{character, what kind of output should be returned. One of 'parsed', 'id_list'
or 'raw' As default, parsed key metadata will be returned as data.frame.
'id_list returns a list of IDs and sources.
Use 'raw' to get full metadata as list. Please be aware that these lists
can become very large.}

\item{page_token}{cursor marking the page}

\item{...}{further params from \code{\link{epmc_search}}}
}
\description{
In general, use \code{\link{epmc_search}} instead. It calls this function, calling all
pages within the defined limit.
}
\seealso{
\link{epmc_search}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_hits_trend.R
\name{epmc_hits_trend}
\alias{epmc_hits_trend}
\title{Get the yearly number of hits for a query and the total
yearly number of hits for a given period}
\usage{
epmc_hits_trend(query, synonym = TRUE, data_src = "med", period = 1975:2016)
}
\arguments{
\item{query}{query in the Europe PMC syntax}

\item{synonym}{logical, synonym search. If TRUE, synonym terms from MeSH
terminology and the UniProt synonym list are queried, too. Disabled by
default.}

\item{data_src}{character, data source, by default Pubmed/MedLine index  (\code{med})
will be searched.
The following three letter codes represent the sources, which are currently
supported
\describe{
  \item{agr}{Agricola is a bibliographic database of citations to the
  agricultural literature created by the US National Agricultural Library
  and its co-operators.}
  \item{cba}{Chinese Biological Abstracts}
  \item{ctx}{CiteXplore}
  \item{eth}{EthOs Theses, i.e. PhD theses (British Library)}
  \item{hir}{NHS Evidence}
  \item{med}{PubMed/Medline NLM}
  \item{nbk}{Europe PMC Book metadata}
  \item{pat}{Biological Patents}
  \item{pmc}{PubMed Central}
  \item{ppr}{Preprint records}
  }}

\item{period}{a vector of years (numeric) over which to perform the search}
}
\value{
a data.frame (dplyr tbl_df) with year, total number of hits
 (all_hits) and number of hits for the query (query_hits)
}
\description{
Get the yearly number of hits for a query and the total
yearly number of hits for a given period
}
\details{
A similar function was used in
 \url{https://masalmon.eu/2017/05/14/evergreenreviewgraph/} where
 it was advised to not plot no. of hits over time for a query,
 but to normalize it by the total no. of hits.
}
\examples{
\dontrun{
# aspirin as query
epmc_hits_trend('aspirin', period = 2006:2016, synonym = FALSE)
# link to cran packages in reference lists
epmc_hits_trend('REF:"cran.r-project.org*"', period = 2006:2016, synonym = FALSE)
# more complex with publication type review
epmc_hits_trend('(REF:"cran.r-project.org*") AND (PUB_TYPE:"Review" OR PUB_TYPE:"review-article")',
period = 2006:2016, synonym = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/europepmc.r
\docType{package}
\name{europepmc}
\alias{europepmc}
\title{europepmc - an R client for the Europe PMC RESTful article API}
\description{
What is europepmc?:

europepmc facilitates access to Europe PMC RESTful Web Service. Europe
PMC covers life science literature and gives access to open access full
texts. Coverage is not only restricted to Europe, but articles and
abstracts are indexed from all over the world. Europe PMC ingests all PubMed
content and extends its index with other sources, including Agricola, a
bibliographic database of citations to the agricultural literature, or
Biological Patents.

Besides searching abstracts and full text, europepmc can be used to
retrieve reference sections and citations, text-mined terms or cross-links
to other databases hosted by the European Bioinformatics Institute (EBI).

For more information about Europe PMC, see their current paper:
Levchenko, M., Gou, Y., Graef, F., Hamelers, A., Huang, Z., Ide-Smith, M.,
 … McEntyre, J. (2017). Europe PMC in 2017. Nucleic Acids Research, 46(D1),
 D1254–D1260. \doi{10.1093/nar/gkx1005}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_refs.r
\name{epmc_refs}
\alias{epmc_refs}
\title{Get references for a given publication}
\usage{
epmc_refs(ext_id = NULL, data_src = "med", limit = 100, verbose = TRUE)
}
\arguments{
\item{ext_id}{character, publication identifier}

\item{data_src}{character, data source, by default Pubmed/MedLine index will
  be searched.


  The following three letter codes represent the sources
  Europe PubMed Central supports:
  \describe{
    \item{agr}{Agricola is a bibliographic database of citations to the
    agricultural literature created by the US National Agricultural Library
    and its co-operators.}
    \item{cba}{Chinese Biological Abstracts}
    \item{ctx}{CiteXplore}
    \item{eth}{EthOs Theses, i.e. PhD theses (British Library)}
    \item{hir}{NHS Evidence}
    \item{med}{PubMed/Medline NLM}
    \item{nbk}{Europe PMC Book metadata}
    \item{pat}{Biological Patents}
    \item{pmc}{PubMed Central}
    }}

\item{limit}{integer, number of results. By default, this function
returns 100 records.}

\item{verbose}{logical, print some information on what is going on.}
}
\value{
returns reference section as tibble
}
\description{
This function retrieves all the works listed in the bibliography of a given
 article.
}
\examples{
\dontrun{
epmc_refs("PMC3166943", data_src = "pmc")
epmc_refs("25378340")
epmc_refs("21753913")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_search_by_doi.r
\name{epmc_search_by_doi}
\alias{epmc_search_by_doi}
\title{Search Europe PMC by DOIs}
\usage{
epmc_search_by_doi(doi = NULL, output = "parsed")
}
\arguments{
\item{doi, }{character vector containing DOI names.}

\item{output}{character, what kind of output should be returned. One of
'parsed', 'id_list' or 'raw' As default, parsed key metadata will be
returned as data.frame. 'id_list' returns a list of IDs and sources. Use
'raw' to get full metadata as list. Please be aware that these lists can
become very large.}
}
\description{
Look up DOIs indexed in Europe PMC and get metadata back.
}
\examples{
\dontrun{
# single DOI name
epmc_search_by_doi(doi = "10.1161/strokeaha.117.018077")
# multiple DOIname in a vector
my_dois <- c(
  "10.1159/000479962",
  "10.1002/sctm.17-0081",
  "10.1161/strokeaha.117.018077",
  "10.1007/s12017-017-8447-9")
epmc_search_by_doi(doi = my_dois)
# full metadata
epmc_search_by_doi(doi = my_dois, output = "raw")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_db.r
\name{epmc_db}
\alias{epmc_db}
\title{Retrieve external database entities referenced in a given publication}
\usage{
epmc_db(
  ext_id = NULL,
  data_src = "med",
  db = NULL,
  limit = 100,
  verbose = TRUE
)
}
\arguments{
\item{ext_id}{character, publication identifier}

\item{data_src}{character, data source, by default Pubmed/MedLine index will
  be searched.


  The following three letter codes represent the sources
  Europe PubMed Central supports:
  \describe{
    \item{agr}{Agricola is a bibliographic database of citations to the
    agricultural literature created by the US National Agricultural Library
    and its co-operators.}
    \item{cba}{Chinese Biological Abstracts}
    \item{ctx}{CiteXplore}
    \item{eth}{EthOs Theses, i.e. PhD theses (British Library)}
    \item{hir}{NHS Evidence}
    \item{med}{PubMed/Medline NLM}
    \item{nbk}{Europe PMC Book metadata}
    \item{pat}{Biological Patents}
    \item{pmc}{PubMed Central}
    }}

\item{db}{character, specify database:
\describe{
\item{'ARXPR'}{Array Express, a database of functional genomics experiments}
\item{'CHEBI'}{a database and ontology of chemical entities of biological
    interest}
 \item{'CHEMBL'}{a database of bioactive drug-like small molecules}
 \item{'EMBL'}{now ENA, provides a comprehensive record of the world's
 nucleotide sequencing information}
 \item{'INTACT'}{provides a freely available, open
    source database system and analysis tools for molecular interaction data}
 \item{'INTERPRO'}{provides functional analysis of proteins by classifying
    them into families and predicting domains and important sites}
 \item{'OMIM'}{a comprehensive and authoritative compendium of human genes and
    genetic phenotypes}
 \item{'PDB'}{European resource for the collection,
    organisation and dissemination of data on biological macromolecular
    structures}
 \item{'UNIPROT'}{comprehensive and freely accessible
    resource of protein sequence and functional information}
 \item{'PRIDE'}{PRIDE Archive - proteomics data repository}}}

\item{limit}{integer, number of results. By default, this function
returns 100 records.}

\item{verbose}{logical, print some information on what is going on.}
}
\value{
Cross-references as data.frame
}
\description{
This function returns EBI database entities referenced in a publication from
Europe PMC RESTful Web Service.
}
\examples{
  \dontrun{
  epmc_db("12368864", db = "uniprot", limit = 150)
  epmc_db("25249410", db = "embl")
  epmc_db("14756321", db = "uniprot")
  epmc_db("11805837", db = "pride")
  }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_lablinks_count.r
\name{epmc_lablinks_count}
\alias{epmc_lablinks_count}
\title{Summarise links to external sources}
\usage{
epmc_lablinks_count(ext_id = NULL, data_src = "med")
}
\arguments{
\item{ext_id}{publication identifier}

\item{data_src}{data source, by default Pubmed/MedLine index will be searched.
The following three letter codes represents the sources
Europe PubMed Central supports:
\describe{
  \item{agr}{Agricola is a bibliographic database of citations to the
  agricultural literature created by the US National Agricultural Library
  and its co-operators.}
  \item{cba}{Chinese Biological Abstracts}
  \item{ctx}{CiteXplore}
  \item{eth}{EthOs Theses, i.e. PhD theses (British Library)}
  \item{hir}{NHS Evidence}
  \item{med}{PubMed/Medline NLM}
  \item{nbk}{Europe PMC Book metadata}
  \item{pat}{Biological Patents}
  \item{pmc}{PubMed Central}
  }}
}
\value{
data.frame with counts for each database
}
\description{
With the External Link services, Europe PMC allows third parties to publish
links from Europe PMC to other webpages or tools. Current External Link
providers, which can be selected through Europe PMC's advanced search,
include Wikipedia, Dryad Digital Repository or the institutional repo of
Bielefeld University. For more information, see
\url{http://europepmc.org/labslink}.
}
\examples{
   \dontrun{
   epmc_lablinks_count("24023770")
   epmc_lablinks_count("PMC3986813", data_src = "pmc")
   }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_search.r
\name{epmc_search}
\alias{epmc_search}
\title{Search Europe PMC publication database}
\usage{
epmc_search(
  query = NULL,
  output = "parsed",
  synonym = TRUE,
  verbose = TRUE,
  limit = 100,
  sort = NULL
)
}
\arguments{
\item{query}{character, search query. For more information on how to build a
search query, see \url{http://europepmc.org/Help}}

\item{output}{character, what kind of output should be returned. One of
'parsed', 'id_list' or 'raw' As default, parsed key metadata will be
returned as data.frame. 'id_list' returns a list of IDs and sources. Use
'raw' to get full metadata as list. Please be aware that these lists can
become very large.}

\item{synonym}{logical, synonym search. If TRUE, synonym terms from MeSH
terminology and the UniProt synonym list are queried, too.
In order to replicate results from the website, with the Rest API
you need to turn synonyms ON!}

\item{verbose}{logical, print progress bar. Activated by default.}

\item{limit}{integer, limit the number of records you wish to retrieve. By
default, 100 are returned.}

\item{sort}{character, relevance ranking is used by default. Use
\code{sort = 'cited'} for sorting by the number of citations, or
\code{sort = 'date'} by the most recent publications.}
}
\value{
tibble
}
\description{
This is the main function to search Europe PMC RESTful Web
  Service (\url{http://europepmc.org/RestfulWebService}). It fully supports
  the comprehensive Europe PMC query language. Simply copy & paste your query
  terms to R. To get familiar with the Europe PMC query syntax, check the
  Advanced Search Query Builder \url{https://europepmc.org/advancesearch}.
}
\examples{
\dontrun{
#Search articles for 'Gabi-Kat'
my.data <- epmc_search(query='Gabi-Kat')

#Get article metadata by DOI
my.data <- epmc_search(query = 'DOI:10.1007/bf00197367')

#Get article metadata by PubMed ID (PMID)
my.data <- epmc_search(query = 'EXT_ID:22246381')

#Get only PLOS Genetics article with EMBL database references
my.data <- epmc_search(query = 'ISSN:1553-7404 HAS_EMBL:y')
#Limit search to 250 PLOS Genetics articles
my.data <- epmc_search(query = 'ISSN:1553-7404', limit = 250)

# exclude MeSH synonyms in search
my.data <- epmc_search(query = 'aspirin', synonym = FALSE)

# get 100 most cited atricles from PLOS ONE publsihed in 2014
epmc_search(query = '(ISSN:1932-6203) AND FIRST_PDATE:2014', sort = 'cited')

# print number of records found
attr(my.data, "hit_count")

# change output

}
}
\seealso{
\url{http://europepmc.org/Help}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_lablinks.r
\name{epmc_lablinks}
\alias{epmc_lablinks}
\title{Get links to external sources}
\usage{
epmc_lablinks(
  ext_id = NULL,
  data_src = "med",
  lab_id = NULL,
  limit = 100,
  verbose = TRUE
)
}
\arguments{
\item{ext_id}{publication identifier}

\item{data_src}{data source, by default Pubmed/MedLine index will be searched.
The following three letter codes represents the sources
Europe PubMed Central supports:
\describe{
  \item{agr}{Agricola is a bibliographic database of citations to the
  agricultural literature created by the US National Agricultural Library
  and its co-operators.}
  \item{cba}{Chinese Biological Abstracts}
  \item{ctx}{CiteXplore}
  \item{eth}{EthOs Theses, i.e. PhD theses (British Library)}
  \item{hir}{NHS Evidence}
  \item{med}{PubMed/Medline NLM}
  \item{nbk}{Europe PMC Book metadata}
  \item{pat}{Biological Patents}
  \item{pmc}{PubMed Central}
  }}

\item{lab_id}{character vector, identifiers of the external link service.
Use Europe PMC's advanced search form to find ids.}

\item{limit}{Number of records to be returned. By default, this function
returns 100 records.}

\item{verbose}{print information about what's going on}
}
\value{
Links found as nested data_frame
}
\description{
With the External Link services, Europe PMC allows third parties to publish
links from Europe PMC to other webpages or tools. Current External Link
providers, which can be selected through Europe PMC's advanced search,
include Wikipedia, Dryad Digital Repository or other open services.
For more information, see
\url{http://europepmc.org/labslink}.
}
\examples{
  \dontrun{
  # Fetch links
  epmc_lablinks("24007304")
  # Link to Altmetric (lab_id = "1562")
  epmc_lablinks("25389392", lab_id = "1562")

  # Links to Wikipedia
  epmc_lablinks("24007304", lab_id = "1507")

  # Link to full text copy archived through the institutional repo of
  Bielefeld University
  epmc_lablinks("12736239", lab_id = "1056")
  }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_profile.r
\name{epmc_profile}
\alias{epmc_profile}
\title{Obtain a summary of hit counts}
\usage{
epmc_profile(query = NULL, synonym = TRUE)
}
\arguments{
\item{query}{character, search query. For more information on how to
build a search query, see \url{http://europepmc.org/Help}}

\item{synonym}{logical, synonym search. If TRUE, synonym terms from MeSH
terminology and the UniProt synonym list are queried, too. Enabled by
default.}
}
\description{
This functions returns the number of results found for your query,
  and breaks it down to the various publication types, data sources, and
  subsets Europe PMC provides.
}
\examples{
\dontrun{
  epmc_profile('malaria')
  # use field search, e.g. query materials and reference section for
  # mentions of "ropensci"
  epmc_profile('(METHODS:"ropensci")')
 }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_ftxt.r
\name{epmc_ftxt}
\alias{epmc_ftxt}
\title{Fetch Europe PMC full texts}
\usage{
epmc_ftxt(ext_id = NULL)
}
\arguments{
\item{ext_id}{character, PMCID. 
All full text publications have external IDs starting 'PMC_'}
}
\value{
xml_document
}
\description{
This function loads full texts into R. Full texts are in XML format and are
only provided for the Open Access subset of Europe PMC.
}
\examples{
  \dontrun{
  epmc_ftxt("PMC3257301")
  epmc_ftxt("PMC3639880")
  }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_db_count.r
\name{epmc_db_count}
\alias{epmc_db_count}
\title{Retrieve the number of database links from Europe PMC publication database}
\usage{
epmc_db_count(ext_id = NULL, data_src = "med")
}
\arguments{
\item{ext_id}{character, publication identifier}

\item{data_src}{character, data source, by default Pubmed/MedLine index will
be searched.}
}
\value{
data.frame with counts for each database
}
\description{
This function returns the number of EBI database links associated with a
publication.
}
\details{
Europe PMC supports cross-references between literature and the
  following databases:
 \describe{
 \item{'ARXPR'}{Array Express, a database of functional genomics experiments}
 \item{'CHEBI'}{a database and ontology of chemical entities of biological
     interest}
  \item{'CHEMBL'}{a database of bioactive drug-like small molecules}
  \item{'EMBL'}{now ENA, provides a comprehensive record of the world's
  nucleotide sequencing information}
  \item{'INTACT'}{provides a freely available, open
     source database system and analysis tools for molecular interaction data}
  \item{'INTERPRO'}{provides functional analysis of proteins by classifying
     them into families and predicting domains and important sites}
  \item{'OMIM'}{a comprehensive and authoritative compendium of human genes and
     genetic phenotypes}
  \item{'PDB'}{European resource for the collection,
     organisation and dissemination of data on biological macromolecular
     structures}
  \item{'UNIPROT'}{comprehensive and freely accessible
     resource of protein sequence and functional information}
  \item{'PRIDE'}{PRIDE Archive - proteomics data repository}}
}
\examples{
  \dontrun{
  epmc_db_count(ext_id = "10779411")
  epmc_db_count(ext_id = "PMC3245140", data_src = "PMC")
  }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations_by_id.R
\name{epmc_annotations_by_id}
\alias{epmc_annotations_by_id}
\title{Get annotations by article}
\usage{
epmc_annotations_by_id(ids = NULL)
}
\arguments{
\item{ids, }{character vector with publication identifiers
following the structure "source:ext_id", e.g. `"MED:28585529"`}
}
\value{
returns text-mined annotations in a tidy format with the following
variables

\describe{
  \item{source}{Publication data source}
  \item{ext_id}{Article Identifier}
  \item{pmcid}{PMCID that locates full-text in Pubmed Central}
  \item{prefix}{Text snipped found before the annotation}
  \item{exact}{Annotated entity}
  \item{postfix}{Text snipped found after the annotation}
  \item{name}{Targeted entity}
  \item{uri}{Uniform link dictionary entry for targeted entity}
  \item{id}{URL to full-text occurence of the annotation}
  \item{type}{Type of annotation like Chemicals}
  \item{section}{Article section mentioning the annotation like Methods}
  \item{provider}{Annotation data provider}
  \item{subtype}{Sub-data provider}
}
}
\description{
Retrieve text-mined annotations contained in abstracts and open access
full-text articles.
}
\examples{
\dontrun{
  annotations_by_id("MED:28585529")
  # multiple ids
  annotations_by_id(c("MED:28585529", "PMC:PMC1664601"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epmc_ftxt_book.r
\name{epmc_ftxt_book}
\alias{epmc_ftxt_book}
\title{Fetch Europe PMC books}
\usage{
epmc_ftxt_book(ext_id = NULL)
}
\arguments{
\item{ext_id}{character, publication identifier. All book full texts are accessible
either by the PMID or the 'NBK' book number.}
}
\value{
xml_document
}
\description{
Use this function to retrieve book XML formatted full text for the Open
Access subset of the Europe PMC bookshelf.
}
\examples{
  \dontrun{
  epmc_ftxt_book("NBK32884")
  }
}
