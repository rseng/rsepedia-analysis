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
