wikitaxa
========



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/wikitaxa)](https://cranchecks.info/pkgs/wikitaxa)
[![R-check](https://github.com/ropensci/wikitaxa/workflows/R-check/badge.svg)](https://github.com/ropensci/wikitaxa/actions/)
[![codecov](https://codecov.io/gh/ropensci/wikitaxa/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/wikitaxa)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/wikitaxa)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/wikitaxa)](https://cran.r-project.org/package=wikitaxa)

`wikitaxa` - taxonomy data from Wikipedia/Wikidata/Wikispecies

`wikitaxa` docs: https://docs.ropensci.org/wikitaxa/

See also the taxize book: https://books.ropensci.org/taxize/


### Low level API

The low level API is meant for power users and gives you more control,
but requires more knowledge.

* `wt_wiki_page()`
* `wt_wiki_page_parse()`
* `wt_wiki_url_build()`
* `wt_wiki_url_parse()`
* `wt_wikispecies_parse()`
* `wt_wikicommons_parse()`
* `wt_wikipedia_parse()`

### High level API

The high level API is meant to be easier and faster to use.

* `wt_data()`
* `wt_data_id()`
* `wt_wikispecies()`
* `wt_wikicommons()`
* `wt_wikipedia()`

Search functions:

* `wt_wikicommons_search()`
* `wt_wikispecies_search()`
* `wt_wikipedia_search()`

## Installation

CRAN version


```r
install.packages("wikitaxa")
```

Dev version


```r
remotes::install_github("ropensci/wikitaxa")
```


```r
library('wikitaxa')
```

## Contributors

* [Ethan Welty](https://github.com/ezwelty)
* [Scott Chamberlain](https://github.com/sckott)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/wikitaxa/issues).
* License: MIT
* Get citation information for `wikitaxa` in R doing `citation(package = 'wikitaxa')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/wikitaxa/blob/master/CODE_OF_CONDUCT.md
wikitaxa 0.4.0
==============

### MINOR IMPROVEMENTS

* remove docs link to httr that's not used in package (#19)
* remove egs from readme, change vignette name


wikitaxa 0.3.0
==============

### MINOR IMPROVEMENTS

* integration with vcr for test caching for all HTTP requests (#17) (#18)
* link to `taxize` book and `wikitaxa` vignette in readme (#16)

### BUG FIXES

* fix to `wt_wikipedia()` to separate `<br>` tags appropriately (#15)


wikitaxa 0.2.0
==============

### BUG FIXES

* `wt_wikicommons()` fails better now when a page does not exist, and is now consitent with the rest of package (#14)
* `wt_wikicommons()` fixed - classification objects were not working correctly as the data used is a hot mess - tried to improve parsing of that text (#13)
* `wt_data()` fix - was failing due to i think a change in the internal pkg `WikidataR` (#12)


wikitaxa 0.1.4
==============

### NEW FEATURES

* `wt_wikipedia()` and `wt_wikipedia_search()` gain parameter `wiki`
to give the wiki language, which defaults to `en` (#9)

### MINOR IMPROVEMENTS

* move some examples to dontrun (#11)


wikitaxa 0.1.0
==============

### NEW FEATURES

* Released to CRAN
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
## Test environments

* local OS X install, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 2 reverse dependencies.
  (Summary at <https://github.com/ropensci/wikitaxa/blob/master/revdep/README.md>). No problems were found.

---

This version fixes a documentation link to a package that is not in Depends,  Imports, or Suggests in this package.

This is a re-submission of this version, removing Remotes from DESCRIPTION.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/wikitaxa/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/wikitaxa.git`
* Make sure to track progress upstream (i.e., on our version of `wikitaxa` at `ropensci/wikitaxa`) by doing `git remote add upstream https://github.com/ropensci/wikitaxa.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/wikitaxa`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

```

```
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.2 (2020-06-22) |
|os       |macOS Catalina 10.15.5       |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |US/Pacific                   |
|date     |2020-06-28                   |

# Dependencies

|package  |old   |new   |Δ  |
|:--------|:-----|:-----|:--|
|wikitaxa |0.3.0 |0.4.0 |*  |
|openssl  |NA    |1.4.1 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*wikitaxa
========

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/wikitaxa)](https://cranchecks.info/pkgs/wikitaxa)
[![R-check](https://github.com/ropensci/wikitaxa/workflows/R-check/badge.svg)](https://github.com/ropensci/wikitaxa/actions/)
[![codecov](https://codecov.io/gh/ropensci/wikitaxa/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/wikitaxa)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/wikitaxa)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/wikitaxa)](https://cran.r-project.org/package=wikitaxa)

`wikitaxa` - taxonomy data from Wikipedia/Wikidata/Wikispecies

`wikitaxa` docs: https://docs.ropensci.org/wikitaxa/

See also the taxize book: https://books.ropensci.org/taxize/


### Low level API

The low level API is meant for power users and gives you more control,
but requires more knowledge.

* `wt_wiki_page()`
* `wt_wiki_page_parse()`
* `wt_wiki_url_build()`
* `wt_wiki_url_parse()`
* `wt_wikispecies_parse()`
* `wt_wikicommons_parse()`
* `wt_wikipedia_parse()`

### High level API

The high level API is meant to be easier and faster to use.

* `wt_data()`
* `wt_data_id()`
* `wt_wikispecies()`
* `wt_wikicommons()`
* `wt_wikipedia()`

Search functions:

* `wt_wikicommons_search()`
* `wt_wikispecies_search()`
* `wt_wikipedia_search()`

## Installation

CRAN version

```{r eval=FALSE}
install.packages("wikitaxa")
```

Dev version

```{r eval=FALSE}
remotes::install_github("ropensci/wikitaxa")
```

```{r}
library('wikitaxa')
```

## Contributors

* [Ethan Welty](https://github.com/ezwelty)
* [Scott Chamberlain](https://github.com/sckott)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/wikitaxa/issues).
* License: MIT
* Get citation information for `wikitaxa` in R doing `citation(package = 'wikitaxa')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/wikitaxa/blob/master/CODE_OF_CONDUCT.md
---
title: "Introduction to the wikitaxa package"
author: "Scott Chamberlain"
date: "2020-06-28"
output: 
    html_document:
        toc: true
        toc_float: true
        theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to the wikitaxa package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



`wikitaxa` - Taxonomy data from Wikipedia

The goal of `wikitaxa` is to allow search and taxonomic data retrieval from
across many Wikimedia sites, including: Wikipedia, Wikicommons, and
Wikispecies.

There are lower level and higher level parts to the package API:

### Low level API

The low level API is meant for power users and gives you more control,
but requires more knowledge.

* `wt_wiki_page()`
* `wt_wiki_page_parse()`
* `wt_wiki_url_build()`
* `wt_wiki_url_parse()`
* `wt_wikispecies_parse()`
* `wt_wikicommons_parse()`
* `wt_wikipedia_parse()`

### High level API

The high level API is meant to be easier and faster to use.

* `wt_data()`
* `wt_data_id()`
* `wt_wikispecies()`
* `wt_wikicommons()`
* `wt_wikipedia()`

Search functions:

* `wt_wikicommons_search()`
* `wt_wikispecies_search()`
* `wt_wikipedia_search()`

## Installation

CRAN version


```r
install.packages("wikitaxa")
```

Dev version


```r
remotes::install_github("ropensci/wikitaxa")
```


```r
library("wikitaxa")
```

## wiki data


```r
z <- wt_data("Poa annua")
names(z)
#> [1] "labels"       "descriptions" "aliases"      "sitelinks"    "claims"
head(z$labels)
#>   language            value
#> 1       pt        Poa annua
#> 2       is   Varpasveifgras
#> 3       pl Wiechlina roczna
#> 4       fr   Pâturin annuel
#> 5       es        Poa annua
#> 6       en        Poa annua
```

Get a Wikidata ID


```r
wt_data_id("Mimulus foliatus")
#> [1] "Q6495130"
#> attr(,"class")
#> [1] "wiki_id"
```

## wikipedia

lower level


```r
pg <- wt_wiki_page("https://en.wikipedia.org/wiki/Malus_domestica")
res <- wt_wiki_page_parse(pg)
res$iwlinks
#> [1] "https://commons.wikimedia.org/wiki/Category:Apples"         
#> [2] "https://commons.wikimedia.org/wiki/Category:Apple_cultivars"
#> [3] "https://www.wikidata.org/wiki/Q158657"                      
#> [4] "https://www.wikidata.org/wiki/Q18674606"                    
#> [5] "https://species.wikimedia.org/wiki/Malus_pumila"            
#> [6] "https://species.wikimedia.org/wiki/Malus_domestica"
```

higher level


```r
res <- wt_wikipedia("Malus domestica")
res$common_names
#> # A tibble: 1 x 2
#>   name  language
#>   <chr> <chr>   
#> 1 Apple en
res$classification
#> # A tibble: 3 x 2
#>   rank       name             
#>   <chr>      <chr>            
#> 1 plainlinks ""               
#> 2 binomial   "Malus domestica"
#> 3 <NA>       ""
```

choose a wikipedia language


```r
# French
wt_wikipedia(name = "Malus domestica", wiki = "fr")
#> $langlinks
#> # A tibble: 60 x 5
#>    lang   url                                   langname   autonym    `*`       
#>    <chr>  <chr>                                 <chr>      <chr>      <chr>     
#>  1 als    https://als.wikipedia.org/wiki/Kultu… Alemannis… Alemannis… Kulturapf…
#>  2 am     https://am.wikipedia.org/wiki/%E1%89… amharique  አማርኛ       ቱፋሕ       
#>  3 ast    https://ast.wikipedia.org/wiki/Malus… asturien   asturianu  Malus dom…
#>  4 az     https://az.wikipedia.org/wiki/M%C9%9… azéri      azərbayca… Mədəni al…
#>  5 bat-s… https://bat-smg.wikipedia.org/wiki/V… Samogitian žemaitėška Vuobelės  
#>  6 bg     https://bg.wikipedia.org/wiki/%D0%94… bulgare    български  Домашна я…
#>  7 bpy    https://bpy.wikipedia.org/wiki/%E0%A… bishnupri… বিষ্ণুপ্র…    আপেল      
#>  8 ca     https://ca.wikipedia.org/wiki/Pomera… catalan    català     Pomera co…
#>  9 cs     https://cs.wikipedia.org/wiki/Jablo%… tchèque    čeština    Jabloň do…
#> 10 csb    https://csb.wikipedia.org/wiki/Dom%C… kachoube   kaszëbsczi Domôcô ja…
#> # … with 50 more rows
#> 
#> $externallinks
#>  [1] "http://www.cabi-publishing.org/pdf/Books/0851995926/0851995926_Chap01.pdf"                 
#>  [2] "http://www.umass.edu/fruitadvisor/fruitnotes/ontheorigin.pdf"                              
#>  [3] "http://www.applegenome.org"                                                                
#>  [4] "http://societeradio-canada.info/emissions/les_annees_lumiere/2010-2011/"                   
#>  [5] "http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.654.html"                           
#>  [6] "https://gallica.bnf.fr/ark:/12148/bpt6k28582v"                                             
#>  [7] "http://worldcat.org/issn/1471-2229&lang=fr"                                                
#>  [8] "https://www.ncbi.nlm.nih.gov/pubmed/26924309"                                              
#>  [9] "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4770685"                                      
#> [10] "https://dx.doi.org/10.1186%2Fs12870-016-0739-y"                                            
#> [11] "https://doi.org/10.1186/s12870-016-0739-y"                                                 
#> [12] "http://www.hamblenne.be/LISTE_RGF.pdf"                                                     
#> [13] "http://www.arcticapples.com/blog/john/demystifying-arctic-apples#.UaeX-Jzjmw5"             
#> [14] "http://www.cctec.cornell.edu/plants/GENEVA-Apple-Rootstocks-Comparison-Chart-120911.pdf"   
#> [15] "https://commons.wikimedia.org/wiki/Category:Malus_domestica?uselang=fr"                    
#> [16] "http://www.tela-botanica.org/page:eflore"                                                  
#> [17] "http://www.tela-botanica.org/bdtfx-nn-40744"                                               
#> [18] "http://www.cbif.gc.ca/acp/fra/siti/regarder?tsn=516655"                                    
#> [19] "http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=516655"      
#> [20] "http://www.cbif.gc.ca/acp/fra/siti/regarder?tsn=25262"                                     
#> [21] "http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=25262"       
#> [22] "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lin=s&p=has_linkout&id=3750"       
#> [23] "http://www.ars-grin.gov/"                                                                  
#> [24] "https://npgsweb.ars-grin.gov/gringlobal/taxonomydetail.aspx?104681"                        
#> [25] "https://www.biolib.cz/en/"                                                                 
#> [26] "https://www.biolib.cz/en/taxon/id39552/"                                                   
#> [27] "http://inpn.mnhn.fr/isb/espece/cd_nom/107207"                                              
#> [28] "http://www.bmlisieux.com/normandie/roblet.htm"                                             
#> [29] "http://site.voila.fr/babadubonsai/docum/docmal.html"                                       
#> [30] "http://www.fruitiers.net"                                                                  
#> [31] "http://www.inra.fr/hyppz/CULTURES/3c---003.htm"                                            
#> [32] "http://cat.inist.fr/?aModele=afficheN&cpsidt=15506238"                                     
#> [33] "http://www.omafra.gov.on.ca/french/crops/facts/98-014.htm"                                 
#> [34] "http://www.gardenaction.co.uk/fruit_veg_diary/fruit_veg_mini_project_september_2_apple.asp"
#> 
#> $common_names
#> # A tibble: 1 x 2
#>   name               language
#>   <chr>              <chr>   
#> 1 Pommier domestique fr      
#> 
#> $classification
#> # A tibble: 0 x 0
#> 
#> $synonyms
#> list()
# Slovak
wt_wikipedia(name = "Malus domestica", wiki = "sk")
#> $langlinks
#> # A tibble: 60 x 5
#>    lang   url                               langname        autonym    `*`      
#>    <chr>  <chr>                             <chr>           <chr>      <chr>    
#>  1 als    https://als.wikipedia.org/wiki/K… Alemannisch     Alemannis… Kulturap…
#>  2 am     https://am.wikipedia.org/wiki/%E… amharčina       አማርኛ       ቱፋሕ      
#>  3 ast    https://ast.wikipedia.org/wiki/M… astúrčina       asturianu  Malus do…
#>  4 az     https://az.wikipedia.org/wiki/M%… azerbajdžančina azərbayca… Mədəni a…
#>  5 bat-s… https://bat-smg.wikipedia.org/wi… Samogitian      žemaitėška Vuobelės 
#>  6 bg     https://bg.wikipedia.org/wiki/%D… bulharčina      български  Домашна …
#>  7 bpy    https://bpy.wikipedia.org/wiki/%… bišnuprijskoma… বিষ্ণুপ্র…    আপেল     
#>  8 ca     https://ca.wikipedia.org/wiki/Po… katalánčina     català     Pomera c…
#>  9 cs     https://cs.wikipedia.org/wiki/Ja… čeština         čeština    Jabloň d…
#> 10 csb    https://csb.wikipedia.org/wiki/D… kašubčina       kaszëbsczi Domôcô j…
#> # … with 50 more rows
#> 
#> $externallinks
#> list()
#> 
#> $common_names
#> # A tibble: 1 x 2
#>   name          language
#>   <chr>         <chr>   
#> 1 Jabloň domáca sk      
#> 
#> $classification
#> # A tibble: 0 x 0
#> 
#> $synonyms
#> list()
# Vietnamese
wt_wikipedia(name = "Malus domestica", wiki = "vi")
#> $langlinks
#> # A tibble: 60 x 5
#>    lang     url                                 langname    autonym   `*`       
#>    <chr>    <chr>                               <chr>       <chr>     <chr>     
#>  1 als      https://als.wikipedia.org/wiki/Kul… Alemannisch Alemanni… Kulturapf…
#>  2 am       https://am.wikipedia.org/wiki/%E1%… Tiếng Amha… አማርኛ      ቱፋሕ       
#>  3 ast      https://ast.wikipedia.org/wiki/Mal… Tiếng Astu… asturianu Malus dom…
#>  4 az       https://az.wikipedia.org/wiki/M%C9… Tiếng Azer… azərbayc… Mədəni al…
#>  5 zh-min-… https://zh-min-nan.wikipedia.org/w… Chinese (M… Bân-lâm-… Phōng-kó-…
#>  6 bg       https://bg.wikipedia.org/wiki/%D0%… Tiếng Bulg… български Домашна я…
#>  7 ca       https://ca.wikipedia.org/wiki/Pome… Tiếng Cata… català    Pomera co…
#>  8 cs       https://cs.wikipedia.org/wiki/Jabl… Tiếng Séc   čeština   Jabloň do…
#>  9 da       https://da.wikipedia.org/wiki/Almi… Tiếng Đan … dansk     Almindeli…
#> 10 de       https://de.wikipedia.org/wiki/Kult… Tiếng Đức   Deutsch   Kulturapf…
#> # … with 50 more rows
#> 
#> $externallinks
#>  [1] "http://biology.umaine.edu/Amelanchier/Rosaceae_2007.pdf"                                       
#>  [2] "//dx.doi.org/10.1007%2Fs00606-007-0539-9"                                                      
#>  [3] "https://npgsweb.ars-grin.gov/gringlobal/taxonomydetail.aspx?410495"                            
#>  [4] "https://npgsweb.ars-grin.gov/gringlobal/taxonomydetail.aspx?30530"                             
#>  [5] "http://www.uga.edu/fruit/apple.html"                                                           
#>  [6] "//dx.doi.org/10.3732%2Fajb.93.3.357"                                                           
#>  [7] "http://www.plosgenetics.org/article/info:doi%2F10.1371%2Fjournal.pgen.1002703"                 
#>  [8] "//www.ncbi.nlm.nih.gov/pmc/articles/PMC3349737"                                                
#>  [9] "//www.ncbi.nlm.nih.gov/pubmed/22589740"                                                        
#> [10] "//dx.doi.org/10.1371%2Fjournal.pgen.1002703"                                                   
#> [11] "http://news.sciencemag.org/sciencenow/2012/05/scienceshot-the-secret-history-o.html"           
#> [12] "http://www.plantpress.com/wildlife/o523-apple.php"                                             
#> [13] "http://cahnrsnews.wsu.edu/2010/08/29/apple-cup-rivals-contribute-to-apple-genome-sequencing/"  
#> [14] "http://www.nature.com/ng/journal/v42/n10/full/ng.654.html"                                     
#> [15] "http://www.alphagalileo.org/ViewItem.aspx?ItemId=83717&CultureCode=en"                         
#> [16] "http://www.ornl.gov/sci/techresources/Human_Genome/project/info.shtml"                         
#> [17] "https://commons.wikimedia.org/wiki/Apple?uselang=vi"                                           
#> [18] "https://commons.wikimedia.org/wiki/Category:Malus_domestica?uselang=vi"                        
#> [19] "http://bachkhoatoanthu.vass.gov.vn/noidung/tudien/Lists/GiaiNghia/View_Detail.aspx?ItemID=5007"
#> [20] "http://www.eol.org/pages/629094"                                                               
#> [21] "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3750"                               
#> [22] "http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=516655"          
#> [23] "http://www.catalogueoflife.org/col/details/species/id/19538828/synonym/19539435"               
#> [24] "https://www.biolib.cz/cz/taxon/id39552"                                                        
#> [25] "https://gd.eppo.int/taxon/MABSD"                                                               
#> [26] "https://npgsweb.ars-grin.gov/gringlobal/taxonomydetail.aspx?id=104681"                         
#> [27] "http://www.ipni.org/ipni/idPlantNameSearch.do?id=726282-1"                                     
#> [28] "https://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=516655"         
#> [29] "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=3750"                    
#> [30] "http://www.nzor.org.nz/names/14d024a2-d821-48e3-95d8-f0dd206c70a0"                             
#> [31] "http://www.pfaf.org/user/Plant.aspx?LatinName=Malus+domestica"                                 
#> [32] "http://www.theplantlist.org/tpl1.1/record/rjp-454"                                             
#> [33] "http://www.plantsoftheworldonline.org/taxon/urn:lsid:ipni.org:names:726282-1"                  
#> [34] "http://legacy.tropicos.org/Name/27804420"                                                      
#> [35] "https://vicflora.rbg.vic.gov.au/flora/taxon/e41b929d-b709-4f4c-8dbe-2a9241e2342b"              
#> [36] "http://www.ipni.org/ipni/idPlantNameSearch.do?id=60476301-2"                                   
#> [37] "http://www.plantsoftheworldonline.org/taxon/urn:lsid:ipni.org:names:60476301-2"                
#> [38] "http://legacy.tropicos.org/Name/100473089"                                                     
#> 
#> $common_names
#> # A tibble: 1 x 2
#>   name            language
#>   <chr>           <chr>   
#> 1 Malus domestica vi      
#> 
#> $classification
#> # A tibble: 0 x 0
#> 
#> $synonyms
#> list()
```

search


```r
wt_wikipedia_search(query = "Pinus")
#> $batchcomplete
#> [1] ""
#> 
#> $continue
#> $continue$sroffset
#> [1] 10
#> 
#> $continue$continue
#> [1] "-||"
#> 
#> 
#> $query
#> $query$searchinfo
#> $query$searchinfo$totalhits
#> [1] 3374
#> 
#> $query$searchinfo$suggestion
#> [1] "penis"
#> 
#> $query$searchinfo$suggestionsnippet
#> [1] "penis"
#> 
#> 
#> $query$search
#> # A tibble: 10 x 7
#>       ns title     pageid  size wordcount snippet                    timestamp  
#>    <int> <chr>      <int> <int>     <int> <chr>                      <chr>      
#>  1     0 Pine       39389 36555      4058 "A pine is any conifer in… 2020-06-26…
#>  2     0 Pinus po… 532941 31087      3069 "misidentified it as <spa… 2020-06-19…
#>  3     0 Pinus co… 507717 20486      2343 "all pines (member specie… 2020-05-25…
#>  4     0 Pinus je… 463015  9130      1008 "long, with a large (15 t… 2019-12-22…
#>  5     0 Pinus st… 464301 31478      3815 "3 ft) tall &amp; wide. M… 2020-06-22…
#>  6     0 Pinus re… 507802  7501       783 "&quot;<span class=\"sear… 2020-05-08…
#>  7     0 Pinus lo… 649634 15408      1741 "sometimes form dense for… 2020-04-14…
#>  8     0 Pinus la… 459402 11464      1338 "Fire affected this speci… 2020-01-15…
#>  9     0 Pinus ni… 438963 11947      1421 "hypodermal cells. P. nig… 2020-04-03…
#> 10     0 Pinus mu… 438946 11964       901 "encyclopedia) is still r… 2020-06-17…
```

search supports languages


```r
wt_wikipedia_search(query = "Pinus", wiki = "fr")
#> $batchcomplete
#> [1] ""
#> 
#> $continue
#> $continue$sroffset
#> [1] 10
#> 
#> $continue$continue
#> [1] "-||"
#> 
#> 
#> $query
#> $query$searchinfo
#> $query$searchinfo$totalhits
#> [1] 990
#> 
#> 
#> $query$search
#> # A tibble: 10 x 7
#>       ns title     pageid  size wordcount snippet                    timestamp  
#>    <int> <chr>      <int> <int>     <int> <chr>                      <chr>      
#>  1     0 Pin (pl…   89798 83647      9325 "<span class=\"searchmatc… 2020-05-23…
#>  2     0 Pinus p…  121544 31274      3892 "<span class=\"searchmatc… 2020-04-30…
#>  3     0 Pinus c…   98421  8237       959 "<span class=\"searchmatc… 2019-05-30…
#>  4     0 Pinus n…  950330 26623      3013 "recycler}}. <span class=… 2020-03-30…
#>  5     0 Pin syl…  121562 13725      1611 "<span class=\"searchmatc… 2020-03-22…
#>  6     0 Pinus h…  117280 22257      2671 "<span class=\"searchmatc… 2020-05-07…
#>  7     0 Pin par…  138378  8763       916 "<span class=\"searchmatc… 2020-04-26…
#>  8     0 Pinus s…  776950 11662      1628 "les articles homonymes, … 2020-04-14…
#>  9     0 Pinus m… 2480854 21747      2310 "<span class=\"searchmatc… 2019-02-25…
#> 10     0 Pinus u… 3208429  6316       720 "significations, voir Pin… 2020-03-09…
```


## wikicommons

lower level


```r
pg <- wt_wiki_page("https://commons.wikimedia.org/wiki/Abelmoschus")
res <- wt_wikicommons_parse(pg)
res$common_names[1:3]
#> [[1]]
#> [[1]]$name
#> [1] "okra"
#> 
#> [[1]]$language
#> [1] "en"
#> 
#> 
#> [[2]]
#> [[2]]$name
#> [1] "مسكي"
#> 
#> [[2]]$language
#> [1] "ar"
#> 
#> 
#> [[3]]
#> [[3]]$name
#> [1] "Abelmoş"
#> 
#> [[3]]$language
#> [1] "az"
```

higher level


```r
res <- wt_wikicommons("Abelmoschus")
res$classification
#> # A tibble: 15 x 2
#>    rank       name            
#>    <chr>      <chr>           
#>  1 Domain     "Eukaryota"     
#>  2 unranked   "Archaeplastida"
#>  3 Regnum     "Plantae"       
#>  4 Cladus     "angiosperms"   
#>  5 Cladus     "eudicots"      
#>  6 Cladus     "core eudicots" 
#>  7 Cladus     "superrosids"   
#>  8 Cladus     "rosids"        
#>  9 Cladus     "eurosids II"   
#> 10 Ordo       "Malvales"      
#> 11 Familia    "Malvaceae"     
#> 12 Subfamilia "Malvoideae"    
#> 13 Tribus     "Hibisceae"     
#> 14 Genus      "Abelmoschus"   
#> 15 Authority  " Medik. (1787)"
res$common_names
#> # A tibble: 19 x 2
#>    name             language
#>    <chr>            <chr>   
#>  1 okra             en      
#>  2 مسكي             ar      
#>  3 Abelmoş          az      
#>  4 Bamja            bs      
#>  5 Ibiškovec        cs      
#>  6 Bisameibisch     de      
#>  7 Okrat            fi      
#>  8 Abelmosco        gl      
#>  9 Abelmošus        hr      
#> 10 Ybiškė           lt      
#> 11 അബെൽമോസ്കസ്        ml      
#> 12 Абельмош         mrj     
#> 13 Abelmoskusslekta nn      
#> 14 Piżmian          pl      
#> 15 Абельмош         ru      
#> 16 Okrasläktet      sv      
#> 17 Абельмош         udm     
#> 18 Chi Vông vang    vi      
#> 19 黄葵属           zh
```

search


```r
wt_wikicommons_search(query = "Pinus")
#> $batchcomplete
#> [1] ""
#> 
#> $continue
#> $continue$sroffset
#> [1] 10
#> 
#> $continue$continue
#> [1] "-||"
#> 
#> 
#> $query
#> $query$searchinfo
#> $query$searchinfo$totalhits
#> [1] 270
#> 
#> 
#> $query$search
#> # A tibble: 10 x 7
#>       ns title                    pageid size  wordcount snippet timestamp      
#>    <int> <chr>                     <int> <lgl>     <int> <chr>   <chr>          
#>  1     0 Pinus sylvestris           9066 NA            0 ""      2020-05-13T19:…
#>  2     0 Pinus ponderosa          250435 NA            0 ""      2020-04-18T15:…
#>  3     0 Pinus nigra               64703 NA            0 ""      2018-03-06T10:…
#>  4     0 Pinus                     82071 NA            0 ""      2017-05-28T10:…
#>  5     0 Pinus mugo               132442 NA            0 ""      2019-07-26T09:…
#>  6     0 Rogów Arboretum        10563490 NA            0 ""      2020-01-01T13:…
#>  7     0 Pinus contorta           186918 NA            0 ""      2020-01-19T19:…
#>  8     0 Anacortes Community F…  2989013 NA            0 ""      2014-12-10T15:…
#>  9     0 Pinus halepensis         172181 NA            0 ""      2018-05-05T10:…
#> 10     0 Pinus brutia             139389 NA            0 ""      2014-11-23T11:…
```


## wikispecies

lower level


```r
pg <- wt_wiki_page("https://species.wikimedia.org/wiki/Malus_domestica")
res <- wt_wikispecies_parse(pg, types = "common_names")
res$common_names[1:3]
#> [[1]]
#> [[1]]$name
#> [1] "Ябълка"
#> 
#> [[1]]$language
#> [1] "български"
#> 
#> 
#> [[2]]
#> [[2]]$name
#> [1] "Poma, pomera"
#> 
#> [[2]]$language
#> [1] "català"
#> 
#> 
#> [[3]]
#> [[3]]$name
#> [1] "jabloň domácí"
#> 
#> [[3]]$language
#> [1] "čeština"
```

higher level


```r
res <- wt_wikispecies("Malus domestica")
res$classification
#> # A tibble: 8 x 2
#>   rank        name         
#>   <chr>       <chr>        
#> 1 Superregnum Eukaryota    
#> 2 Regnum      Plantae      
#> 3 Cladus      Angiosperms  
#> 4 Cladus      Eudicots     
#> 5 Cladus      Core eudicots
#> 6 Cladus      Rosids       
#> 7 Cladus      Eurosids I   
#> 8 Ordo        Rosales
res$common_names
#> # A tibble: 22 x 2
#>    name          language  
#>    <chr>         <chr>     
#>  1 Ябълка        български 
#>  2 Poma, pomera  català    
#>  3 jabloň domácí čeština   
#>  4 Apfel         Deutsch   
#>  5 Μηλιά         Ελληνικά  
#>  6 Apple         English   
#>  7 Manzano       español   
#>  8 Aed-õunapuu   eesti     
#>  9 Tarhaomenapuu suomi     
#> 10 Aapel         Nordfriisk
#> # … with 12 more rows
```

search


```r
wt_wikispecies_search(query = "Pinus")
#> $batchcomplete
#> [1] ""
#> 
#> $continue
#> $continue$sroffset
#> [1] 10
#> 
#> $continue$continue
#> [1] "-||"
#> 
#> 
#> $query
#> $query$searchinfo
#> $query$searchinfo$totalhits
#> [1] 515
#> 
#> 
#> $query$search
#> # A tibble: 10 x 7
#>       ns title         pageid  size wordcount snippet                timestamp  
#>    <int> <chr>          <int> <int>     <int> <chr>                  <chr>      
#>  1     0 Pinus         1.74e4  5737       784 "Familia: Pinaceae Ge… 2020-06-06…
#>  2     0 Pinus halepe… 4.51e4  4047       580 "Pinaceae Genus: <spa… 2019-12-20…
#>  3     0 Pinus pinea   4.51e4  1949       406 "Familia: Pinaceae Ge… 2019-10-19…
#>  4     0 Pinus veitch… 1.34e6  1450       181 "Familia: Pinaceae Ge… 2019-07-19…
#>  5     0 Pinus pumila  7.35e4  1395       189 "Pinaceae Genus: <spa… 2019-07-14…
#>  6     0 Pinus subg. … 3.01e5   358        27 "Pinaceae Genus: <spa… 2019-11-24…
#>  7     0 Pinus clausa  4.50e4  1552       208 "Pinaceae Genus: <spa… 2019-08-15…
#>  8     0 Pinus pseudo… 1.48e6  2114       310 "Genus: <span class=\… 2020-05-21…
#>  9     0 Pinus pinast… 1.32e6  2764       379 "Pinaceae Genus: <spa… 2019-12-20…
#> 10     0 Pinus nigra … 3.27e5  1799       138 "Genus: <span class=\… 2020-03-02…
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikipages.R
\name{wt_wiki_url_parse}
\alias{wt_wiki_url_parse}
\title{Parse MediaWiki Page URL}
\usage{
wt_wiki_url_parse(url)
}
\arguments{
\item{url}{(character) MediaWiki page url.}
}
\value{
a list with elements:
\itemize{
\item wiki - wiki language
\item type - wikipedia type
\item page - page name
}
}
\description{
Parse a MediaWiki page url into its component parts (wiki name, wiki type,
and page title). Supports both static page urls and their equivalent API
calls.
}
\examples{
wt_wiki_url_parse(url="https://en.wikipedia.org/wiki/Malus_domestica")
wt_wiki_url_parse("https://en.wikipedia.org/w/api.php?page=Malus_domestica")
}
\seealso{
Other MediaWiki functions: 
\code{\link{wt_wiki_page_parse}()},
\code{\link{wt_wiki_page}()},
\code{\link{wt_wiki_url_build}()}
}
\concept{MediaWiki functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikipages.R
\name{wt_wiki_page_parse}
\alias{wt_wiki_page_parse}
\title{Parse MediaWiki Page}
\usage{
wt_wiki_page_parse(
  page,
  types = c("langlinks", "iwlinks", "externallinks"),
  tidy = FALSE
)
}
\arguments{
\item{page}{(\link[crul:HttpResponse]{crul::HttpResponse}) Result of \code{\link[=wt_wiki_page]{wt_wiki_page()}}}

\item{types}{(character) List of properties to parse.}

\item{tidy}{(logical). tidy output to data.frames when possible.
Default: \code{FALSE}}
}
\value{
a list
}
\description{
Parses common properties from the result of a MediaWiki API page call.
}
\details{
Available properties currently not parsed:
title, displaytitle, pageid, revid, redirects, text, categories,
links, templates, images, sections, properties, ...
}
\examples{
\dontrun{
pg <- wt_wiki_page("https://en.wikipedia.org/wiki/Malus_domestica")
wt_wiki_page_parse(pg)
}
}
\seealso{
Other MediaWiki functions: 
\code{\link{wt_wiki_page}()},
\code{\link{wt_wiki_url_build}()},
\code{\link{wt_wiki_url_parse}()}
}
\concept{MediaWiki functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikicommons.R
\name{wt_wikicommons}
\alias{wt_wikicommons}
\alias{wt_wikicommons_parse}
\alias{wt_wikicommons_search}
\title{WikiCommons}
\usage{
wt_wikicommons(name, utf8 = TRUE, ...)

wt_wikicommons_parse(
  page,
  types = c("langlinks", "iwlinks", "externallinks", "common_names", "classification"),
  tidy = FALSE
)

wt_wikicommons_search(query, limit = 10, offset = 0, utf8 = TRUE, ...)
}
\arguments{
\item{name}{(character) Wiki name - as a page title, must be length 1}

\item{utf8}{(logical) If \code{TRUE}, encodes most (but not all) non-ASCII
characters as UTF-8 instead of replacing them with hexadecimal escape
sequences. Default: \code{TRUE}}

\item{...}{curl options, passed on to \code{\link[httr:GET]{httr::GET()}}}

\item{page}{(\code{\link[httr:response]{httr::response()}}) Result of \code{\link[=wt_wiki_page]{wt_wiki_page()}}}

\item{types}{(character) List of properties to parse}

\item{tidy}{(logical). tidy output to data.frame's if possible.
Default: \code{FALSE}}

\item{query}{(character) query terms}

\item{limit}{(integer) number of results to return. Default: 10}

\item{offset}{(integer) record to start at. Default: 0}
}
\value{
\code{wt_wikicommons} returns a list, with slots:
\itemize{
\item langlinks - language page links
\item externallinks - external links
\item common_names - a data.frame with \code{name} and \code{language} columns
\item classification - a data.frame with \code{rank} and \code{name} columns
}

\code{wt_wikicommons_parse} returns a list

\code{wt_wikicommons_search} returns a list with slots for \code{continue} and
\code{query}, where \code{query} holds the results, with \code{query$search} slot with
the search results
}
\description{
WikiCommons
}
\examples{
\dontrun{
# high level
wt_wikicommons(name = "Malus domestica")
wt_wikicommons(name = "Pinus contorta")
wt_wikicommons(name = "Ursus americanus")
wt_wikicommons(name = "Balaenoptera musculus")

wt_wikicommons(name = "Category:Poeae")
wt_wikicommons(name = "Category:Pinaceae")

# low level
pg <- wt_wiki_page("https://commons.wikimedia.org/wiki/Malus_domestica")
wt_wikicommons_parse(pg)

# search wikicommons
# FIXME: utf=FALSE for now until curl::curl_escape fix 
# https://github.com/jeroen/curl/issues/228
wt_wikicommons_search(query = "Pinus", utf8 = FALSE)

## use search results to dig into pages
res <- wt_wikicommons_search(query = "Pinus", utf8 = FALSE)
lapply(res$query$search$title[1:3], wt_wikicommons)
}
}
\references{
\url{https://www.mediawiki.org/wiki/API:Search} for help on search
}
\concept{Wikicommons functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikipages.R
\name{wt_wiki_url_build}
\alias{wt_wiki_url_build}
\title{Build MediaWiki Page URL}
\usage{
wt_wiki_url_build(
  wiki,
  type = NULL,
  page = NULL,
  api = FALSE,
  action = "parse",
  redirects = TRUE,
  format = "json",
  utf8 = TRUE,
  prop = c("text", "langlinks", "categories", "links", "templates", "images",
    "externallinks", "sections", "revid", "displaytitle", "iwlinks", "properties")
)
}
\arguments{
\item{wiki}{(character | list) Either the wiki name or a list with
\verb{$wiki}, \verb{$type}, and \verb{$page} (the output of \code{\link[=wt_wiki_url_parse]{wt_wiki_url_parse()}}).}

\item{type}{(character) Wiki type.}

\item{page}{(character) Wiki page title.}

\item{api}{(boolean) Whether to return an API call or a static page url
(default). If \code{FALSE}, all following (API-only) arguments are ignored.}

\item{action}{(character) See \url{https://en.wikipedia.org/w/api.php}
for supported actions. This function currently only supports "parse".}

\item{redirects}{(boolean) If the requested page is set to a redirect,
resolve it.}

\item{format}{(character) See \url{https://en.wikipedia.org/w/api.php}
for supported output formats.}

\item{utf8}{(boolean) If \code{TRUE}, encodes most (but not all) non-ASCII
characters as UTF-8 instead of replacing them with hexadecimal escape
sequences.}

\item{prop}{(character) Properties to retrieve, either as a character vector
or pipe-delimited string. See
\url{https://en.wikipedia.org/w/api.php?action=help&modules=parse} for
supported properties.}
}
\value{
a URL (character)
}
\description{
Builds a MediaWiki page url from its component parts (wiki name, wiki type,
and page title). Supports both static page urls and their equivalent API
calls.
}
\examples{
wt_wiki_url_build(wiki = "en", type = "wikipedia", page = "Malus domestica")
wt_wiki_url_build(
  wt_wiki_url_parse("https://en.wikipedia.org/wiki/Malus_domestica"))
wt_wiki_url_build("en", "wikipedia", "Malus domestica", api = TRUE)
}
\seealso{
Other MediaWiki functions: 
\code{\link{wt_wiki_page_parse}()},
\code{\link{wt_wiki_page}()},
\code{\link{wt_wiki_url_parse}()}
}
\concept{MediaWiki functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikipedia.R
\name{wt_wikipedia}
\alias{wt_wikipedia}
\alias{wt_wikipedia_parse}
\alias{wt_wikipedia_search}
\title{Wikipedia}
\usage{
wt_wikipedia(name, wiki = "en", utf8 = TRUE, ...)

wt_wikipedia_parse(
  page,
  types = c("langlinks", "iwlinks", "externallinks", "common_names", "classification"),
  tidy = FALSE
)

wt_wikipedia_search(
  query,
  wiki = "en",
  limit = 10,
  offset = 0,
  utf8 = TRUE,
  ...
)
}
\arguments{
\item{name}{(character) Wiki name - as a page title, must be length 1}

\item{wiki}{(character) wiki language. default: en. See \link{wikipedias} for
language codes.}

\item{utf8}{(logical) If \code{TRUE}, encodes most (but not all) non-ASCII
characters as UTF-8 instead of replacing them with hexadecimal escape
sequences. Default: \code{TRUE}}

\item{...}{curl options, passed on to \code{\link[httr:GET]{httr::GET()}}}

\item{page}{(\code{\link[httr:response]{httr::response()}}) Result of \code{\link[=wt_wiki_page]{wt_wiki_page()}}}

\item{types}{(character) List of properties to parse}

\item{tidy}{(logical). tidy output to data.frame's if possible.
Default: \code{FALSE}}

\item{query}{(character) query terms}

\item{limit}{(integer) number of results to return. Default: 10}

\item{offset}{(integer) record to start at. Default: 0}
}
\value{
\code{wt_wikipedia} returns a list, with slots:
\itemize{
\item langlinks - language page links
\item externallinks - external links
\item common_names - a data.frame with \code{name} and \code{language} columns
\item classification - a data.frame with \code{rank} and \code{name} columns
\item synonyms - a character vector with taxonomic names
}

\code{wt_wikipedia_parse} returns a list with same slots determined by
the \code{types} parmeter

\code{wt_wikipedia_search} returns a list with slots for \code{continue} and
\code{query}, where \code{query} holds the results, with \code{query$search} slot with
the search results
}
\description{
Wikipedia
}
\examples{
\dontrun{
# high level
wt_wikipedia(name = "Malus domestica")
wt_wikipedia(name = "Malus domestica", wiki = "fr")
wt_wikipedia(name = "Malus domestica", wiki = "da")

# low level
pg <- wt_wiki_page("https://en.wikipedia.org/wiki/Malus_domestica")
wt_wikipedia_parse(pg)
wt_wikipedia_parse(pg, tidy = TRUE)

# search wikipedia
# FIXME: utf=FALSE for now until curl::curl_escape fix 
# https://github.com/jeroen/curl/issues/228
wt_wikipedia_search(query = "Pinus", utf8=FALSE)
wt_wikipedia_search(query = "Pinus", wiki = "fr", utf8=FALSE)
wt_wikipedia_search(query = "Pinus", wiki = "br", utf8=FALSE)

## curl options
# wt_wikipedia_search(query = "Pinus", verbose = TRUE, utf8=FALSE)

## use search results to dig into pages
res <- wt_wikipedia_search(query = "Pinus", utf8=FALSE)
lapply(res$query$search$title[1:3], wt_wikipedia)
}
}
\references{
\url{https://www.mediawiki.org/wiki/API:Search} for help on search
}
\concept{Wikipedia functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wiki.R
\name{wt_data}
\alias{wt_data}
\alias{wt_data_id}
\title{Wikidata taxonomy data}
\usage{
wt_data(x, property = NULL, ...)

wt_data_id(x, language = "en", limit = 10, ...)
}
\arguments{
\item{x}{(character) a taxonomic name}

\item{property}{(character) a property id, e.g., P486}

\item{...}{curl options passed on to \code{httr::GET()}}

\item{language}{(character) two letter language code}

\item{limit}{(integer) records to return. Default: 10}
}
\value{
\code{wt_data} searches Wikidata, and returns a list with elements:
\itemize{
\item labels - data.frame with columns: language, value
\item descriptions - data.frame with columns: language, value
\item aliases - data.frame with columns: language, value
\item sitelinks - data.frame with columns: site, title
\item claims - data.frame with columns: claims, property_value,
property_description, value (comma separted values in string)
}

\code{wt_data_id} gets the Wikidata ID for the searched term, and
returns the ID as character
}
\description{
Wikidata taxonomy data
}
\details{
Note that \code{wt_data} can take a while to run since when fetching
claims it has to do so one at a time for each claim

You can search things other than taxonomic names with \code{wt_data} if you
like
}
\examples{
\dontrun{
# search by taxon name
# wt_data("Mimulus alsinoides")

# choose which properties to return
wt_data(x="Mimulus foliatus", property = c("P846", "P815"))

# get a taxonomic identifier
wt_data_id("Mimulus foliatus")
# the id can be passed directly to wt_data()
# wt_data(wt_data_id("Mimulus foliatus"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikipages.R
\name{wt_wiki_page}
\alias{wt_wiki_page}
\title{Get MediaWiki Page from API}
\usage{
wt_wiki_page(url, ...)
}
\arguments{
\item{url}{(character) MediaWiki page url.}

\item{...}{Arguments passed to \code{\link[=wt_wiki_url_build]{wt_wiki_url_build()}} if \code{url}
is a static page url.}
}
\value{
an \code{HttpResponse} response object from \pkg{crul}
}
\description{
Supports both static page urls and their equivalent API calls.
}
\details{
If the URL given is for a human readable html page,
we convert it to equivalent API call - if URL is already an API call,
we just use that.
}
\examples{
\dontrun{
wt_wiki_page("https://en.wikipedia.org/wiki/Malus_domestica")
}
}
\seealso{
Other MediaWiki functions: 
\code{\link{wt_wiki_page_parse}()},
\code{\link{wt_wiki_url_build}()},
\code{\link{wt_wiki_url_parse}()}
}
\concept{MediaWiki functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikitaxa-package.R
\docType{package}
\name{wikitaxa-package}
\alias{wikitaxa-package}
\alias{wikitaxa}
\title{wikitaxa}
\description{
Taxonomic Information from Wikipedia
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}

Ethan Welty
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikitaxa-package.R
\docType{data}
\name{wikipedias}
\alias{wikipedias}
\title{List of Wikipedias}
\description{
data.frame of 295 rows, with 3 columns:
\itemize{
\item language - language
\item language_local - language in local name
\item wiki - langugae code for the wiki
}
}
\details{
From \url{https://meta.wikimedia.org/wiki/List_of_Wikipedias}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikispecies.R
\name{wt_wikispecies}
\alias{wt_wikispecies}
\alias{wt_wikispecies_parse}
\alias{wt_wikispecies_search}
\title{WikiSpecies}
\usage{
wt_wikispecies(name, utf8 = TRUE, ...)

wt_wikispecies_parse(
  page,
  types = c("langlinks", "iwlinks", "externallinks", "common_names", "classification"),
  tidy = FALSE
)

wt_wikispecies_search(query, limit = 10, offset = 0, utf8 = TRUE, ...)
}
\arguments{
\item{name}{(character) Wiki name - as a page title, must be length 1}

\item{utf8}{(logical) If \code{TRUE}, encodes most (but not all) non-ASCII
characters as UTF-8 instead of replacing them with hexadecimal escape
sequences. Default: \code{TRUE}}

\item{...}{curl options, passed on to \code{\link[httr:GET]{httr::GET()}}}

\item{page}{(\code{\link[httr:response]{httr::response()}}) Result of \code{\link[=wt_wiki_page]{wt_wiki_page()}}}

\item{types}{(character) List of properties to parse}

\item{tidy}{(logical). tidy output to data.frame's if possible.
Default: \code{FALSE}}

\item{query}{(character) query terms}

\item{limit}{(integer) number of results to return. Default: 10}

\item{offset}{(integer) record to start at. Default: 0}
}
\value{
\code{wt_wikispecies} returns a list, with slots:
\itemize{
\item langlinks - language page links
\item externallinks - external links
\item common_names - a data.frame with \code{name} and \code{language} columns
\item classification - a data.frame with \code{rank} and \code{name} columns
}

\code{wt_wikispecies_parse} returns a list

\code{wt_wikispecies_search} returns a list with slots for \code{continue} and
\code{query}, where \code{query} holds the results, with \code{query$search} slot with
the search results
}
\description{
WikiSpecies
}
\examples{
\dontrun{
# high level
wt_wikispecies(name = "Malus domestica")
wt_wikispecies(name = "Pinus contorta")
wt_wikispecies(name = "Ursus americanus")
wt_wikispecies(name = "Balaenoptera musculus")

# low level
pg <- wt_wiki_page("https://species.wikimedia.org/wiki/Abelmoschus")
wt_wikispecies_parse(pg)

# search wikispecies
# FIXME: utf=FALSE for now until curl::curl_escape fix 
# https://github.com/jeroen/curl/issues/228
wt_wikispecies_search(query = "pine tree", utf8=FALSE)

## use search results to dig into pages
res <- wt_wikispecies_search(query = "pine tree", utf8=FALSE)
lapply(res$query$search$title[1:3], wt_wikispecies)
}
}
\references{
\url{https://www.mediawiki.org/wiki/API:Search} for help on search
}
\concept{Wikispecies functions}
