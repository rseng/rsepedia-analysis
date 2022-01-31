rbhl
====



[![cran checks](https://cranchecks.info/badges/worst/rbhl)](https://cranchecks.info/pkgs/rbhl)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/rbhl/workflows/R-check/badge.svg)](https://github.com/ropensci/rbhl/actions?query=workflow%3AR-check)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rbhl?color=2ECC71)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rbhl)](https://cran.r-project.org/package=rbhl)
[![codecov.io](https://codecov.io/github/ropensci/rbhl/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rbhl?branch=master)

`rbhl` is an R interface to the Biodiversity Heritage Library API.

## Info

Authentication:

* Get your Biodiversity Heritage Library API key at <https://www.biodiversitylibrary.org/getapikey.aspx>
* Put your API key either in your .Renviron file (or similar) with the name `BHL_KEY`, or your .Rprofile file with the name `bhl_key`, and we'll use your API key without you having to enter it every time you run a search.

Documentation:

* Biodiversity Heritage Library API documentation at <https://www.biodiversitylibrary.org/docs/api3.html>
* Biodiversity Heritage Library OpenURL documentation at <https://www.biodiversitylibrary.org/openurlhelp.aspx>

## Installation

Stable version


```r
install.packages("rbhl")
```

Development version


```r
remotes::install_github("ropensci/rbhl")
```


```r
library("rbhl")
```

## Output formats

You can output various formats using the `as` parameter, setting to `table`, `list`, `json` or `xml`.

The default is usually `table`:


```r
bhl_authorsearch(name='dimmock')
#> # A tibble: 6 x 4
#>   AuthorID Name                 CreatorUrl                               Dates  
#>   <chr>    <chr>                <chr>                                    <chr>  
#> 1 189035   "Dimmock, Anna Kath… https://www.biodiversitylibrary.org/cre… <NA>   
#> 2 59023    "Dimmock, G"         https://www.biodiversitylibrary.org/cre… <NA>   
#> 3 189042   "Dimmock, Geo "      https://www.biodiversitylibrary.org/cre… <NA>   
#> 4 189021   "Dimmock, George"    https://www.biodiversitylibrary.org/cre… <NA>   
#> 5 1970     "Dimmock, George,"   https://www.biodiversitylibrary.org/cre… 1852-  
#> 6 8126     "Dimmock, George,"   https://www.biodiversitylibrary.org/cre… 1852-1…
```

list output


```r
bhl_authorsearch(name='dimmock', as='list')$Result[[1]]
#> $AuthorID
#> [1] "189035"
#> 
#> $Name
#> [1] "Dimmock, Anna Katherina"
#> 
#> $CreatorUrl
#> [1] "https://www.biodiversitylibrary.org/creator/189035"
```

XML output


```r
bhl_authorsearch(name='dimmock', as='xml')
#> [1] "﻿<?xml version=\"1.0\" encoding=\"utf-8\"?><Response xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"><Status>ok</Status><ErrorMessage /><Result><Author><AuthorID>189035</AuthorID><Name>Dimmock, Anna Katherina</Name><CreatorUrl>https://www.biodiversitylibrary.org/creator/189035</CreatorUrl></Author><Author><AuthorID>59023</AuthorID><Name>Dimmock, G</Name><CreatorUrl>https://www.biodiversitylibrary.org/creator/59023</CreatorUrl></Author><Author><AuthorID>189042</AuthorID><Name>Dimmock, Geo </Name><CreatorUrl>https://www.biodiversitylibrary.org/creator/189042</CreatorUrl></Author><Author><AuthorID>189021</AuthorID><Name>Dimmock, George</Name><CreatorUrl>https://www.biodiversitylibrary.org/creator/189021</CreatorUrl></Author><Author><AuthorID>1970</AuthorID><Name>Dimmock, George,</Name><Dates>1852-</Dates><CreatorUrl>https://www.biodiversitylibrary.org/creator/1970</CreatorUrl></Author><Author><AuthorID>8126</AuthorID><Name>Dimmock, George,</Name><Dates>1852-1930</Dates><CreatorUrl>https://www.biodiversitylibrary.org/creator/8126</CreatorUrl></Author></Result></Response>"
```

JSON output


```r
bhl_authorsearch(name='dimmock', as='json')
#> [1] "{\"Status\":\"ok\",\"ErrorMessage\":\"\",\"Result\":[{\"AuthorID\":\"189035\",\"Name\":\"Dimmock, Anna Katherina\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189035\"},{\"AuthorID\":\"59023\",\"Name\":\"Dimmock, G\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/59023\"},{\"AuthorID\":\"189042\",\"Name\":\"Dimmock, Geo \",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189042\"},{\"AuthorID\":\"189021\",\"Name\":\"Dimmock, George\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189021\"},{\"AuthorID\":\"1970\",\"Name\":\"Dimmock, George,\",\"Dates\":\"1852-\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/1970\"},{\"AuthorID\":\"8126\",\"Name\":\"Dimmock, George,\",\"Dates\":\"1852-1930\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/8126\"}]}"
```

## Get title metadata


```r
bhl_gettitlemetadata(1726, items = TRUE, as = "list")$Result[[1]]$Items
#> [[1]]
#> [[1]]$ItemID
#> [1] 16800
#> 
#> [[1]]$ItemUrl
#> [1] "https://www.biodiversitylibrary.org/item/16800"
```

## Book search


```r
bhl_publicationsearch('cocos island costa rica birds')
#> # A tibble: 2 x 17
#>   BHLType FoundIn  Volume Authors  PartUrl  PartID Genre Title   ContainerTitle 
#>   <chr>   <chr>    <chr>  <list>   <chr>    <chr>  <chr> <chr>   <chr>          
#> 1 Part    Metadata 2      <df[,1]… https:/… 69838  Arti… Field … Proceedings of…
#> 2 Item    Metadata <NA>   <df[,1]… <NA>     <NA>   Book  The bi… <NA>           
#> # … with 8 more variables: Series <chr>, Date <chr>, PageRange <chr>,
#> #   ItemID <chr>, TitleID <chr>, ExternalUrl <chr>, ItemUrl <chr>,
#> #   TitleUrl <chr>
```

## Get languages


```r
bhl_getlanguages()
#> # A tibble: 80 x 2
#>    LanguageCode LanguageName        
#>    <chr>        <chr>               
#>  1 AFR          Afrikaans           
#>  2 ARA          Arabic              
#>  3 ARC          Aramaic             
#>  4 MAP          Austronesian (Other)
#>  5 BUL          Bulgarian           
#>  6 BUR          Burmese             
#>  7 CAR          Carib               
#>  8 CAT          Catalan             
#>  9 CEL          Celtic (Other)      
#> 10 CHI          Chinese             
#> # … with 70 more rows
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rbhl/issues).
* License: MIT
* Get citation information for `rbhl` in R doing `citation(package = 'rbhl')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
rbhl 0.9.2
===============

### MINOR IMPROVEMENTS

* vignette fix


rbhl 0.9.0
===============

### NEW FEATURES

* Re-worked package functions to work with BHL API version 3. Some functions are defunct, see `?rbhl-defunct` after installation. Some BHL API methods were removed and others folded into other API methods; see https://www.biodiversitylibrary.org/docs/api3.html for all the details (#25) (#27)

### MINOR IMPROVEMENTS

* fix missing title in vignette (#26)
* fix instructions for setting a BHL API key (#24)


rbhl 0.8.0
===============

### MINOR IMPROVEMENTS

* moved to markdown docs and tidied docs to 80 line width (#22)
* replaced `httr` with `crul` for http client (#21)
* added many more tests (#2)
* `bhl_getpartendnote()` and `bhl_gettitleendnote` methods removed
due to their removal in the BHL API

rbhl 0.3.0
===============

### NEW FEATURES

* Added code of conduct
* New function `bhl_getinstitutions()` (#16)
* the compact `data.frame` outputs created when `as = "table"` are now
replaced by using `tibble`. in addition, the data.frame used to be
printed as a data.frame, but was actually inside a list accessible by
`x$data` to get to the data.frame. You no longer need to do that. Now
you run a function, and the output is a data.frame (and a `tbl_df`) (#18)
* Now using `xml2` instead of `XML`

### BUG FIXES

* `bhl_getpartnames` was down but it back up now (#17)


rbhl 0.2.0
===============

### NEW FEATURES

* Gains new funtions `bhl_getauthoparts()`, `bhl_getitemparts()`,
`bhl_getpartbibtex()`, `bhl_getpartbyidentifier()`, `bhl_getpartendnote()`,
`bhl_getpartmetadata()`, `bhl_getpartnames()`, `bhl_getsubjectparts()`,
`bhl_getunpublishedparts()`, and `bhl_partsearch()` (#6)
* New function `getpages()` to get multiple OCR pages at once, given
an itemid (#4)

### BUG FIXES

* Import all non-base R functions (#13)
* Fix to README, requested by CRAN, to an image that didn't have an alt
attribute (#8) (#11) (#12)
* Removed arguments in `bhl_getpageocrtext()` for names and ocr that
weren't actually available in the API method. (#3)
* Some functions were missing parameters, and have now been added in where
appropriate. (#5)

rbhl 0.1.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local macOS install, R 4.0.5 Patched
* ubuntu 14.04 (on GitHub Actions), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version fixes a vignette issue.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/rbhl/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rbhl.git`
* Make sure to track progress upstream (i.e., on our version of `rbhl` at `ropensci/rbhl`) by doing `git remote add upstream https://github.com/ropensci/rbhl.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rbhl`

### Forum

Our [discussion forum](https://discuss.ropensci.org) is another good place to get help.
<!-- Do not share your api key in this issue -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rbhl vignette}
%\VignetteEncoding{UTF-8}
-->



rbhl Introduction
=================

`rbhl` is an R interface to the Biodiversity Heritage Library.

## Installation

Stable version from CRAN


```r
install.packages("rbhl")
```

Or, development version from GitHub


```r
remotes::install_github("ropensci/rbhl")
```

Then load library


```r
library("rbhl")
```

## Output formats

You can output various formats using the `as` parameter, setting to 'json', 'xml', 'list', or 'table'. Note that not all options are available in all functions due to varying returned data objects.

Raw output, in xml format


```r
bhl_authorsearch(name='dimmock', as='xml')
```

json format


```r
bhl_authorsearch(name='dimmock', as='json')
```

```
#> [1] "{\"Status\":\"ok\",\"ErrorMessage\":\"\",\"Result\":[{\"AuthorID\":\"189035\",\"Name\":\"Dimmock, Anna Katherina\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189035\"},{\"AuthorID\":\"59023\",\"Name\":\"Dimmock, G\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/59023\"},{\"AuthorID\":\"189042\",\"Name\":\"Dimmock, Geo \",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189042\"},{\"AuthorID\":\"189021\",\"Name\":\"Dimmock, George\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189021\"},{\"AuthorID\":\"1970\",\"Name\":\"Dimmock, George,\",\"Dates\":\"1852-\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/1970\"},{\"AuthorID\":\"8126\",\"Name\":\"Dimmock, George,\",\"Dates\":\"1852-1930\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/8126\"}]}"
```

Or get a list


```r
bhl_authorsearch(name='dimmock', as='list')
```

```
#> $Status
#> [1] "ok"
#> 
#> $ErrorMessage
#> [1] ""
#> 
#> $Result
#> $Result[[1]]
#> $Result[[1]]$AuthorID
#> [1] "189035"
...
```

The default option in most functions is to parse to a data.frame


```r
bhl_authorsearch(name='dimmock')
```

```
#> # A tibble: 6 x 4
#>   AuthorID Name                      CreatorUrl                                         Dates    
#>   <chr>    <chr>                     <chr>                                              <chr>    
#> 1 189035   "Dimmock, Anna Katherina" https://www.biodiversitylibrary.org/creator/189035 <NA>     
#> 2 59023    "Dimmock, G"              https://www.biodiversitylibrary.org/creator/59023  <NA>     
#> 3 189042   "Dimmock, Geo "           https://www.biodiversitylibrary.org/creator/189042 <NA>     
#> 4 189021   "Dimmock, George"         https://www.biodiversitylibrary.org/creator/189021 <NA>     
#> 5 1970     "Dimmock, George,"        https://www.biodiversitylibrary.org/creator/1970   1852-    
#> 6 8126     "Dimmock, George,"        https://www.biodiversitylibrary.org/creator/8126   1852-1930
```

## Some examples of function calls

### Get title metadata


```r
bhl_gettitlemetadata(1726, items = TRUE)
```

```
#> $Status
#> [1] "ok"
#> 
#> $ErrorMessage
#> [1] ""
#> 
#> $Result
#> $Result[[1]]
#> $Result[[1]]$TitleID
#> [1] 1726
...
```

### Book search


```r
bhl_publicationsearch('cocos island costa rica birds')
```

```
#> # A tibble: 2 x 17
#>   BHLType FoundIn  Volume Authors   PartUrl       PartID Genre  Title            ContainerTitle       Series Date  PageRange ItemID TitleID ExternalUrl   ItemUrl      TitleUrl      
#>   <chr>   <chr>    <chr>  <list>    <chr>         <chr>  <chr>  <chr>            <chr>                <chr>  <chr> <chr>     <chr>  <chr>   <chr>         <chr>        <chr>         
#> 1 Part    Metadata 2      <df[,1] … https://www.… 69838  Artic… Field notes on … Proceedings of the … 4      1919  189--258  <NA>   <NA>    <NA>          <NA>         <NA>          
#> 2 Item    Metadata <NA>   <df[,1] … <NA>          <NA>   Book   The birds of Co… <NA>                 <NA>   <NA>  <NA>      167380 89765   http://hdl.h… https://www… https://www.b…
```

### Get languages

This function gets a list of languages in which materials in BHL have been written. 


```r
bhl_getlanguages()
```

```
#> # A tibble: 75 x 2
#>    LanguageCode LanguageName        
#>    <chr>        <chr>               
#>  1 AFR          Afrikaans           
#>  2 ARA          Arabic              
#>  3 ARC          Aramaic             
#>  4 MAP          Austronesian (Other)
#>  5 BUL          Bulgarian           
#>  6 BUR          Burmese             
#>  7 CAR          Carib               
#>  8 CAT          Catalan             
#>  9 CEL          Celtic (Other)      
#> 10 CHI          Chinese             
#> # … with 65 more rows
```

## Use case for BHL data

Here, we search for keywords, then plot the data using `ggplot2`


```r
x <- bhl_publicationsearch('Selborne', year=1825, collection=4, language='eng')
pages <- bhl_pages(itemid = x$ItemID[1])
```

Examine some of the pages


```r
lapply(pages[7:10], "[[", "OcrText")
```

```
#> $`25934962`
#> [1] "^\\ \n\n\n\nC^nX \n\n\n\nN^ \n\n\n\nACCOUNT \n\n\n\nOF THE \n\n\n\nREV. GILBERT WHITE. \n\n\n\n\" Gilbert White was the eldest son of John \n\n White, of Selborne, Esq., and of Anne, the \n\n daughter of Thomas Holt, rector of Streatham \n\n in Surrey. He was born at Selborne on July 18, \n\n 1720, and received his school education at Bas- \n\n ingstoke, under the Rev. Thomas Warton, vicar \n\n of that place, and father of those two distin- \n\n guished literary characters, Dr. Joseph Warton, \n\n master of Winchester school, and Mr. Thomas \n\n Warton, poetry professor at Oxford. He was \n\n admitted at Oriel College, Oxford, in December, \n\n 1739, and took his degree of Bachelor of Arts \n\n in June, 1743. In March, 1744, he was elected \n\n Fellow of his college. He became Master of \n\n Arts in October, 1746, and was admitted one of \n\n the senior proctors of the University in April, \n\n 1752. Being of an unambitious temper, and \n\n strongly attached to the charms of rural scenery, \n\n he early fixed his residence in his native village, \n\n\n\nflnromr raiL»£4 fi» 6 \n\n\n\nN. C. State College \n\n\n\n"
#> 
#> $`25934961`
#> [1] "Vlll \n\n\n\nwhere he spent the greater part of his life in lit- \n\n erary occupations, and especially in the study of \n\n Nature. This he followed with patient assiduity, \n\n and a mind ever open to the lessons of piety and \n\n benevolence, which such a study is so well cal- \n\n culated to afford. Though several occasions of- \n\n fered of settling upon a college living, he could \n\n never persuade himself to quit the beloved spot, \n\n which was, indeed, a peculiarly happy situation \n\n for an observer. Thus his days passed tranquil \n\n and serene, with scarcely any other vicissitudes \n\n than those of the seasons, till they closed at a \n\n mature age on June 26, 1793.\" \n\n\n\na «a .k i* \n\n\n\n"
#> 
#> $`25934960`
#> [1] "INDEX. \n\n\n\nPage \n\n\n\nAffections of Birds . . 181 \n\n AylesHolt 39 \n\n\n\nBarometers 315 \n\n\n\nBat 49, 86, 116 \n\n\n\nBat, Large 98 \n\n\n\nBirds, Language of . . 266 \n\n\n\nBirds of Passage 68, 145, 197 \n\n\n\nBlack Act 32 \n\n\n\nBlackcap 129 \n\n\n\nBotany 257 \n\n\n\nBotany of Selborne . . 259 \n\n\n\nBullfinch 63 \n\n\n\nBunting 57 \n\n\n\nBurning of Heath ... 34 \n\n\n\nButcher-bird .... 77 \n\n\n\nBuzzards, Honey . . . 137 \n\n\n\n... 53 \n\n\n\n... 240 \n\n\n\n55, 63, 166 \n\n\n\n... 308 \n\n\n\n... 273 \n\n\n\n... 294 \n\n\n\n236 \n\n\n\n176 \n\n\n\n21 \n\n\n\n277 \n\n\n\n280 \n\n\n\n284 \n\n\n\nCanary Birds \n\n\n\nCat and Leveret \n\n\n\nChaffinches . \n\n\n\nChinese Dog \n\n\n\nCliff, Fall of a . . . \n\n\n\nCoccus \n\n\n\nCondensation by Trees \n\n Congregation of Birds \n\n Cornua Ammonis . . \n\n Crickets, Field . . . \n\n Crickets, House . . \n\n Crickets, Mole . . . \n\n\n\nCrossbeaks 164 \n\n\n\nCrossbills 48 \n\n\n\nCuckoo .... 155, 157 \n\n Curious Fossil Shell . . 21 \n\n Curlew, Stone 65, 79, 112, 311 \n\n\n\nDeer . . . \n\n\n\nDeer, Moose \n\n Deer-stealers \n\n Dogs . . . \n\n Dove, Ring . \n\n Dove, Stock \n\n\n\nEagle-owl \n\n Echoes . . \n\n Elm, Broad-leaved \n\n\n\n30,60 \n\n 101, 106 \n\n . 33 \n\n . 309 \n\n . 142 \n\n . 139 \n\n\n\n. 97 \n\n\n\n251,314 \n\n\n\n. 17 \n\n\n\n96 \n\n\n\nEmployments . \n\n\n\nFalco .... \n\n Falcon, Peregrine \n\n Fieldfares . . \n\n Fish .... \n\n Fishes, Gold and Silver \n\n Flight of Birds . \n\n Fly-catcher . . \n\n Forest of Wolmer \n\n Forest-fly . . . \n\n Fossil Shell . . \n\n Fossil Wood \n\n Fowls, Language of \n\n Freestone . . \n\n\n\nGame .... \n\n Garden . . . \n\n German Silktail \n\n Gipsies . . . \n\n Goat-sucker . . \n\n Gossamer . . \n\n Grosbeaks . . \n\n\n\nHarvest Mouse . \n\n Harvest Bug \n\n Hawk and Hens \n\n Hawk, Sparrow \n\n Haws .... \n\n Hedgehogs . . \n\n Heliotropes . . \n\n Herons . . . \n\n Himantopus . . \n\n Hoopoes . . . \n\n House Cat . . \n\n\n\nIdiot Boy, Propens: \n\n Instinct \n\n\n\nJackdaws \n\n\n\nLakes .... \n\n Lamperns . . \n\n Land-springs \n\n Lark, Willow . \n\n Lark, Grasshopper \n\n Leprosy . . . \n\n\n\nPage \n\n\n\n27 \n\n\n\n47 \n\n\n\n307 \n\n\n\n100, 167 \n\n\n\n48 \n\n\n\n297 \n\n\n\n263 \n\n\n\n44,63 \n\n\n\n27 \n\n 188 \n\n\n\n21 \n\n 311 \n\n 266 \n\n\n\n22 \n\n\n\n29 \n\n 249 \n\n\n\n52 \n\n 229 \n\n\n\n85 \n\n 224 \n\n\n\n48 \n\n\n\nity \n\n\n\n59 \n\n 113 \n\n 269 \n\n 138 \n\n\n\n52 \n\n\n\n99 \n\n 271 \n\n\n\n84 \n\n 285 \n\n\n\n47 \n\n 105, 240 \n\n\n\nof an 234 \n\n 302 \n\n\n\n81 \n\n\n\n37 \n\n 74 \n\n\n\n207 \n\n 76 \n\n 67 \n\n\n\n246 \n\n\n\n"
#> 
#> $`25934959`
#> [1] "X \n\n\n\nINDEX. \n\n\n\nLinnets 56 \n\n\n\nLizard 79, 87 \n\n\n\nLizard, Green .... 83 \n\n Loaches ...... 74 \n\n\n\nManor of Selborne . . 25 \n\n Maps of Scotland . . . 136 \n\n\n\nMartins 119 \n\n\n\nMartins, House 119,194,291,300 \n\n Martins, Sand .... 208 \n\n Martins, Black .... 213 \n\n\n\nMice 51 \n\n\n\nMigrating Birds . . 88, 169 \n\n Migrations of Grallae . . 175 \n\n Missel-thrush .... 222 \n\n\n\nNotes of Owls and Cuckoos 174 \n\n Nuthatch 69 \n\n\n\nOaks 16, 18 \n\n\n\nOrnithology of Selborne . 121 \n\n\n\nOtter 105 \n\n\n\nOwl, Fern . . . 117,163 \n\n\n\nOwls 49 \n\n\n\nOwls, White .... 184 \n\n\n\nPairing of Birds . . . 104 \n\n\n\nPeacocks 115 \n\n\n\nPettichaps 305 \n\n\n\nPonds on Chalk-hills . . 237 \n\n\n\nPoor 26 \n\n\n\nPopulation 26 \n\n\n\nPulveratrices .... 163 \n\n\n\nRain 26, 313 \n\n\n\nRaven, Tree .... 19 \n\n\n\nRedwings 167 \n\n\n\nReed-sparrow . . . . 161 \n\n Ringousels, 78, 82, 89, 90, 107* \n\n 120, 164, 198 \n\n Rooks 184,312 \n\n\n\nRooks, White \n\n Rush Candles \n\n\n\nSalads . . . \n\n Salicaria . . \n\n Sandpiper \n\n Sandstone \n\n\n\n62 \n\n 231 \n\n\n\n250 \n\n 94 \n\n 77 \n\n 24 \n\n\n\nScopoh's Annus Primus 109,162 \n\n Sheep 196 \n\n\n\nPago \n\n\n\nSinging Birds . . 149, 157 \n\n Singing Birds, Silence of 159 \n\n\n\nSmother-fly 296 \n\n\n\nSnipes 69 \n\n\n\nSociality of Brutes . . 227 \n\n Soft-billed Birds ... 132 \n\n\n\nSoils 16 \n\n\n\nSpiracula of Animals . . 60 \n\n Sticklebacks .... 74 \n\n\n\nStreams 15 \n\n\n\nSummer Birds of Passage 68 \n\n 145, 197 \n\n Summer Birds, Return of 222 \n\n Summer Evening Walk . 91 \n\n Sussex Downs .... 195 \n\n Swallows 42, 53, 88, 119, 178 \n\n 187, 200, 305 \n\n Swallows, Torpidity of . 244 \n\n Swifts . . 98, 213, 256, 292 \n\n\n\nTeals 184 \n\n\n\nThe Holt 40 \n\n\n\nThe Plestor 17 \n\n\n\nTitlark 153 \n\n\n\nTitmouse 134 \n\n\n\nToads 70, 83 \n\n\n\nTortoise .... 179,288 \n\n Tortoise, Land .... 165 \n\n Turnip-fly 114 \n\n\n\nVernal and Autumnal Cro- \n\n cus 262 \n\n\n\nVillage of Selborne . . 13 \n\n Vipers 72,239 \n\n\n\nWater-rats .... 45, 97 \n\n\n\nWater-newt 72 \n\n\n\nWater-eft 75 \n\n\n\nWeather 331 \n\n\n\nWheatear 57 \n\n\n\nWinter Birds of Passage 148 \n\n Wolmer Pond .... 37 \n\n Woodcocks . . . 167, 171 \n\n Wood-pigeons .... 139 \n\n\n\nWorms 242 \n\n\n\nWren, Golden-crowned . 69 \n\n Wren, Willow .... 66 \n\n\n\nYellow-hammer . . . 153 \n\n\n\n"
```

With this text data, you can do text mining of the content to answer cool questions you have. 
rbhl
====

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![cran checks](https://cranchecks.info/badges/worst/rbhl)](https://cranchecks.info/pkgs/rbhl)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/rbhl/workflows/R-check/badge.svg)](https://github.com/ropensci/rbhl/actions?query=workflow%3AR-check)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rbhl?color=2ECC71)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/rbhl)](https://cran.r-project.org/package=rbhl)
[![codecov.io](https://codecov.io/github/ropensci/rbhl/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rbhl?branch=master)

`rbhl` is an R interface to the Biodiversity Heritage Library API.

## Info

Authentication:

* Get your Biodiversity Heritage Library API key at <https://www.biodiversitylibrary.org/getapikey.aspx>
* Put your API key either in your .Renviron file (or similar) with the name `BHL_KEY`, or your .Rprofile file with the name `bhl_key`, and we'll use your API key without you having to enter it every time you run a search.

Documentation:

* Biodiversity Heritage Library API documentation at <https://www.biodiversitylibrary.org/docs/api3.html>
* Biodiversity Heritage Library OpenURL documentation at <https://www.biodiversitylibrary.org/openurlhelp.aspx>

## Installation

Stable version

```{r eval=FALSE}
install.packages("rbhl")
```

Development version

```{r eval=FALSE}
remotes::install_github("ropensci/rbhl")
```

```{r}
library("rbhl")
```

## Output formats

You can output various formats using the `as` parameter, setting to `table`, `list`, `json` or `xml`.

The default is usually `table`:

```{r}
bhl_authorsearch(name='dimmock')
```

list output

```{r}
bhl_authorsearch(name='dimmock', as='list')$Result[[1]]
```

XML output

```{r}
bhl_authorsearch(name='dimmock', as='xml')
```

JSON output

```{r}
bhl_authorsearch(name='dimmock', as='json')
```

## Get title metadata

```{r}
bhl_gettitlemetadata(1726, items = TRUE, as = "list")$Result[[1]]$Items
```

## Book search

```{r}
bhl_publicationsearch('cocos island costa rica birds')
```

## Get languages

```{r}
bhl_getlanguages()
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rbhl/issues).
* License: MIT
* Get citation information for `rbhl` in R doing `citation(package = 'rbhl')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{rbhl vignette}
%\VignetteEncoding{UTF-8}
-->

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
# options(encoding = "native.enc")
knitr::opts_chunk$set(
  comment = "#>", 
  warning = FALSE, 
  message = FALSE,
  cache = FALSE
)

```

rbhl Introduction
=================

`rbhl` is an R interface to the Biodiversity Heritage Library.

## Installation

Stable version from CRAN

```{r installrbhlcran, eval=FALSE}
install.packages("rbhl")
```

Or, development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/rbhl")
```

Then load library

```{r loadrbhl}
library("rbhl")
```

## Output formats

You can output various formats using the `as` parameter, setting to 'json', 'xml', 'list', or 'table'. Note that not all options are available in all functions due to varying returned data objects.

Raw output, in xml format

```{r authorsearch1, messages=FALSE, warning=FALSE, eval=FALSE}
bhl_authorsearch(name='dimmock', as='xml')
```

json format

```{r authorsearch2, messages=FALSE, warning=FALSE}
bhl_authorsearch(name='dimmock', as='json')
```

Or get a list

```{r output.lines=10}
bhl_authorsearch(name='dimmock', as='list')
``` 

The default option in most functions is to parse to a data.frame

```{r authorsearch4, messages=FALSE, warning=FALSE}
bhl_authorsearch(name='dimmock')
``` 

## Some examples of function calls

### Get title metadata

```{r output.lines=10}
bhl_gettitlemetadata(1726, items = TRUE)
```

### Book search

```{r booksearch, messages=FALSE, warning=FALSE}
bhl_publicationsearch('cocos island costa rica birds')
```

### Get languages

This function gets a list of languages in which materials in BHL have been written. 

```{r}
bhl_getlanguages()
```

## Use case for BHL data

Here, we search for keywords, then plot the data using `ggplot2`

```{r cache=TRUE}
x <- bhl_publicationsearch('Selborne', year=1825, collection=4, language='eng')
pages <- bhl_pages(itemid = x$ItemID[1])
```

Examine some of the pages

```{r}
lapply(pages[7:10], "[[", "OcrText")
```

With this text data, you can do text mining of the content to answer cool questions you have. 
---
title: Introduction to rbhl
author: Scott Chamberlain
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

`rbhl` is an R interface to the Biodiversity Heritage Library.

## Installation

Stable version from CRAN


```r
install.packages("rbhl")
```

Or, development version from GitHub


```r
remotes::install_github("ropensci/rbhl")
```

Then load library


```r
library("rbhl")
```

## Output formats

You can output various formats using the `as` parameter, setting to 'json', 'xml', 'list', or 'table'. Note that not all options are available in all functions due to varying returned data objects.

Raw output, in xml format


```r
bhl_authorsearch(name='dimmock', as='xml')
```

json format


```r
bhl_authorsearch(name='dimmock', as='json')
```

```
#> [1] "{\"Status\":\"ok\",\"ErrorMessage\":\"\",\"Result\":[{\"AuthorID\":\"189035\",\"Name\":\"Dimmock, Anna Katherina\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189035\"},{\"AuthorID\":\"59023\",\"Name\":\"Dimmock, G\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/59023\"},{\"AuthorID\":\"189042\",\"Name\":\"Dimmock, Geo \",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189042\"},{\"AuthorID\":\"189021\",\"Name\":\"Dimmock, George\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/189021\"},{\"AuthorID\":\"1970\",\"Name\":\"Dimmock, George,\",\"Dates\":\"1852-\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/1970\"},{\"AuthorID\":\"8126\",\"Name\":\"Dimmock, George,\",\"Dates\":\"1852-1930\",\"CreatorUrl\":\"https://www.biodiversitylibrary.org/creator/8126\"}]}"
```

Or get a list


```r
bhl_authorsearch(name='dimmock', as='list')
```

```
#> $Status
#> [1] "ok"
#> 
#> $ErrorMessage
#> [1] ""
#> 
#> $Result
#> $Result[[1]]
#> $Result[[1]]$AuthorID
#> [1] "189035"
...
```

The default option in most functions is to parse to a data.frame


```r
bhl_authorsearch(name='dimmock')
```

```
#> # A tibble: 6 x 4
#>   AuthorID Name                      CreatorUrl                                         Dates    
#>   <chr>    <chr>                     <chr>                                              <chr>    
#> 1 189035   "Dimmock, Anna Katherina" https://www.biodiversitylibrary.org/creator/189035 <NA>     
#> 2 59023    "Dimmock, G"              https://www.biodiversitylibrary.org/creator/59023  <NA>     
#> 3 189042   "Dimmock, Geo "           https://www.biodiversitylibrary.org/creator/189042 <NA>     
#> 4 189021   "Dimmock, George"         https://www.biodiversitylibrary.org/creator/189021 <NA>     
#> 5 1970     "Dimmock, George,"        https://www.biodiversitylibrary.org/creator/1970   1852-    
#> 6 8126     "Dimmock, George,"        https://www.biodiversitylibrary.org/creator/8126   1852-1930
```

## Some examples of function calls

### Get title metadata


```r
bhl_gettitlemetadata(1726, items = TRUE)
```

```
#> $Status
#> [1] "ok"
#> 
#> $ErrorMessage
#> [1] ""
#> 
#> $Result
#> $Result[[1]]
#> $Result[[1]]$TitleID
#> [1] 1726
...
```

### Book search


```r
bhl_publicationsearch('cocos island costa rica birds')
```

```
#> # A tibble: 2 x 17
#>   BHLType FoundIn  Volume Authors   PartUrl       PartID Genre  Title            ContainerTitle       Series Date  PageRange ItemID TitleID ExternalUrl   ItemUrl      TitleUrl      
#>   <chr>   <chr>    <chr>  <list>    <chr>         <chr>  <chr>  <chr>            <chr>                <chr>  <chr> <chr>     <chr>  <chr>   <chr>         <chr>        <chr>         
#> 1 Part    Metadata 2      <df[,1] … https://www.… 69838  Artic… Field notes on … Proceedings of the … 4      1919  189--258  <NA>   <NA>    <NA>          <NA>         <NA>          
#> 2 Item    Metadata <NA>   <df[,1] … <NA>          <NA>   Book   The birds of Co… <NA>                 <NA>   <NA>  <NA>      167380 89765   http://hdl.h… https://www… https://www.b…
```

### Get languages

This function gets a list of languages in which materials in BHL have been written. 


```r
bhl_getlanguages()
```

```
#> # A tibble: 75 x 2
#>    LanguageCode LanguageName        
#>    <chr>        <chr>               
#>  1 AFR          Afrikaans           
#>  2 ARA          Arabic              
#>  3 ARC          Aramaic             
#>  4 MAP          Austronesian (Other)
#>  5 BUL          Bulgarian           
#>  6 BUR          Burmese             
#>  7 CAR          Carib               
#>  8 CAT          Catalan             
#>  9 CEL          Celtic (Other)      
#> 10 CHI          Chinese             
#> # … with 65 more rows
```

## Use case for BHL data

Here, we search for keywords, then plot the data using `ggplot2`


```r
x <- bhl_publicationsearch('Selborne', year=1825, collection=4, language='eng')
pages <- bhl_pages(itemid = x$ItemID[1])
```

Examine some of the pages


```r
lapply(pages[7:10], "[[", "OcrText")
```

```
#> $`25934962`
#> [1] "^\\ \n\n\n\nC^nX \n\n\n\nN^ \n\n\n\nACCOUNT \n\n\n\nOF THE \n\n\n\nREV. GILBERT WHITE. \n\n\n\n\" Gilbert White was the eldest son of John \n\n White, of Selborne, Esq., and of Anne, the \n\n daughter of Thomas Holt, rector of Streatham \n\n in Surrey. He was born at Selborne on July 18, \n\n 1720, and received his school education at Bas- \n\n ingstoke, under the Rev. Thomas Warton, vicar \n\n of that place, and father of those two distin- \n\n guished literary characters, Dr. Joseph Warton, \n\n master of Winchester school, and Mr. Thomas \n\n Warton, poetry professor at Oxford. He was \n\n admitted at Oriel College, Oxford, in December, \n\n 1739, and took his degree of Bachelor of Arts \n\n in June, 1743. In March, 1744, he was elected \n\n Fellow of his college. He became Master of \n\n Arts in October, 1746, and was admitted one of \n\n the senior proctors of the University in April, \n\n 1752. Being of an unambitious temper, and \n\n strongly attached to the charms of rural scenery, \n\n he early fixed his residence in his native village, \n\n\n\nflnromr raiL»£4 fi» 6 \n\n\n\nN. C. State College \n\n\n\n"
#> 
#> $`25934961`
#> [1] "Vlll \n\n\n\nwhere he spent the greater part of his life in lit- \n\n erary occupations, and especially in the study of \n\n Nature. This he followed with patient assiduity, \n\n and a mind ever open to the lessons of piety and \n\n benevolence, which such a study is so well cal- \n\n culated to afford. Though several occasions of- \n\n fered of settling upon a college living, he could \n\n never persuade himself to quit the beloved spot, \n\n which was, indeed, a peculiarly happy situation \n\n for an observer. Thus his days passed tranquil \n\n and serene, with scarcely any other vicissitudes \n\n than those of the seasons, till they closed at a \n\n mature age on June 26, 1793.\" \n\n\n\na «a .k i* \n\n\n\n"
#> 
#> $`25934960`
#> [1] "INDEX. \n\n\n\nPage \n\n\n\nAffections of Birds . . 181 \n\n AylesHolt 39 \n\n\n\nBarometers 315 \n\n\n\nBat 49, 86, 116 \n\n\n\nBat, Large 98 \n\n\n\nBirds, Language of . . 266 \n\n\n\nBirds of Passage 68, 145, 197 \n\n\n\nBlack Act 32 \n\n\n\nBlackcap 129 \n\n\n\nBotany 257 \n\n\n\nBotany of Selborne . . 259 \n\n\n\nBullfinch 63 \n\n\n\nBunting 57 \n\n\n\nBurning of Heath ... 34 \n\n\n\nButcher-bird .... 77 \n\n\n\nBuzzards, Honey . . . 137 \n\n\n\n... 53 \n\n\n\n... 240 \n\n\n\n55, 63, 166 \n\n\n\n... 308 \n\n\n\n... 273 \n\n\n\n... 294 \n\n\n\n236 \n\n\n\n176 \n\n\n\n21 \n\n\n\n277 \n\n\n\n280 \n\n\n\n284 \n\n\n\nCanary Birds \n\n\n\nCat and Leveret \n\n\n\nChaffinches . \n\n\n\nChinese Dog \n\n\n\nCliff, Fall of a . . . \n\n\n\nCoccus \n\n\n\nCondensation by Trees \n\n Congregation of Birds \n\n Cornua Ammonis . . \n\n Crickets, Field . . . \n\n Crickets, House . . \n\n Crickets, Mole . . . \n\n\n\nCrossbeaks 164 \n\n\n\nCrossbills 48 \n\n\n\nCuckoo .... 155, 157 \n\n Curious Fossil Shell . . 21 \n\n Curlew, Stone 65, 79, 112, 311 \n\n\n\nDeer . . . \n\n\n\nDeer, Moose \n\n Deer-stealers \n\n Dogs . . . \n\n Dove, Ring . \n\n Dove, Stock \n\n\n\nEagle-owl \n\n Echoes . . \n\n Elm, Broad-leaved \n\n\n\n30,60 \n\n 101, 106 \n\n . 33 \n\n . 309 \n\n . 142 \n\n . 139 \n\n\n\n. 97 \n\n\n\n251,314 \n\n\n\n. 17 \n\n\n\n96 \n\n\n\nEmployments . \n\n\n\nFalco .... \n\n Falcon, Peregrine \n\n Fieldfares . . \n\n Fish .... \n\n Fishes, Gold and Silver \n\n Flight of Birds . \n\n Fly-catcher . . \n\n Forest of Wolmer \n\n Forest-fly . . . \n\n Fossil Shell . . \n\n Fossil Wood \n\n Fowls, Language of \n\n Freestone . . \n\n\n\nGame .... \n\n Garden . . . \n\n German Silktail \n\n Gipsies . . . \n\n Goat-sucker . . \n\n Gossamer . . \n\n Grosbeaks . . \n\n\n\nHarvest Mouse . \n\n Harvest Bug \n\n Hawk and Hens \n\n Hawk, Sparrow \n\n Haws .... \n\n Hedgehogs . . \n\n Heliotropes . . \n\n Herons . . . \n\n Himantopus . . \n\n Hoopoes . . . \n\n House Cat . . \n\n\n\nIdiot Boy, Propens: \n\n Instinct \n\n\n\nJackdaws \n\n\n\nLakes .... \n\n Lamperns . . \n\n Land-springs \n\n Lark, Willow . \n\n Lark, Grasshopper \n\n Leprosy . . . \n\n\n\nPage \n\n\n\n27 \n\n\n\n47 \n\n\n\n307 \n\n\n\n100, 167 \n\n\n\n48 \n\n\n\n297 \n\n\n\n263 \n\n\n\n44,63 \n\n\n\n27 \n\n 188 \n\n\n\n21 \n\n 311 \n\n 266 \n\n\n\n22 \n\n\n\n29 \n\n 249 \n\n\n\n52 \n\n 229 \n\n\n\n85 \n\n 224 \n\n\n\n48 \n\n\n\nity \n\n\n\n59 \n\n 113 \n\n 269 \n\n 138 \n\n\n\n52 \n\n\n\n99 \n\n 271 \n\n\n\n84 \n\n 285 \n\n\n\n47 \n\n 105, 240 \n\n\n\nof an 234 \n\n 302 \n\n\n\n81 \n\n\n\n37 \n\n 74 \n\n\n\n207 \n\n 76 \n\n 67 \n\n\n\n246 \n\n\n\n"
#> 
#> $`25934959`
#> [1] "X \n\n\n\nINDEX. \n\n\n\nLinnets 56 \n\n\n\nLizard 79, 87 \n\n\n\nLizard, Green .... 83 \n\n Loaches ...... 74 \n\n\n\nManor of Selborne . . 25 \n\n Maps of Scotland . . . 136 \n\n\n\nMartins 119 \n\n\n\nMartins, House 119,194,291,300 \n\n Martins, Sand .... 208 \n\n Martins, Black .... 213 \n\n\n\nMice 51 \n\n\n\nMigrating Birds . . 88, 169 \n\n Migrations of Grallae . . 175 \n\n Missel-thrush .... 222 \n\n\n\nNotes of Owls and Cuckoos 174 \n\n Nuthatch 69 \n\n\n\nOaks 16, 18 \n\n\n\nOrnithology of Selborne . 121 \n\n\n\nOtter 105 \n\n\n\nOwl, Fern . . . 117,163 \n\n\n\nOwls 49 \n\n\n\nOwls, White .... 184 \n\n\n\nPairing of Birds . . . 104 \n\n\n\nPeacocks 115 \n\n\n\nPettichaps 305 \n\n\n\nPonds on Chalk-hills . . 237 \n\n\n\nPoor 26 \n\n\n\nPopulation 26 \n\n\n\nPulveratrices .... 163 \n\n\n\nRain 26, 313 \n\n\n\nRaven, Tree .... 19 \n\n\n\nRedwings 167 \n\n\n\nReed-sparrow . . . . 161 \n\n Ringousels, 78, 82, 89, 90, 107* \n\n 120, 164, 198 \n\n Rooks 184,312 \n\n\n\nRooks, White \n\n Rush Candles \n\n\n\nSalads . . . \n\n Salicaria . . \n\n Sandpiper \n\n Sandstone \n\n\n\n62 \n\n 231 \n\n\n\n250 \n\n 94 \n\n 77 \n\n 24 \n\n\n\nScopoh's Annus Primus 109,162 \n\n Sheep 196 \n\n\n\nPago \n\n\n\nSinging Birds . . 149, 157 \n\n Singing Birds, Silence of 159 \n\n\n\nSmother-fly 296 \n\n\n\nSnipes 69 \n\n\n\nSociality of Brutes . . 227 \n\n Soft-billed Birds ... 132 \n\n\n\nSoils 16 \n\n\n\nSpiracula of Animals . . 60 \n\n Sticklebacks .... 74 \n\n\n\nStreams 15 \n\n\n\nSummer Birds of Passage 68 \n\n 145, 197 \n\n Summer Birds, Return of 222 \n\n Summer Evening Walk . 91 \n\n Sussex Downs .... 195 \n\n Swallows 42, 53, 88, 119, 178 \n\n 187, 200, 305 \n\n Swallows, Torpidity of . 244 \n\n Swifts . . 98, 213, 256, 292 \n\n\n\nTeals 184 \n\n\n\nThe Holt 40 \n\n\n\nThe Plestor 17 \n\n\n\nTitlark 153 \n\n\n\nTitmouse 134 \n\n\n\nToads 70, 83 \n\n\n\nTortoise .... 179,288 \n\n Tortoise, Land .... 165 \n\n Turnip-fly 114 \n\n\n\nVernal and Autumnal Cro- \n\n cus 262 \n\n\n\nVillage of Selborne . . 13 \n\n Vipers 72,239 \n\n\n\nWater-rats .... 45, 97 \n\n\n\nWater-newt 72 \n\n\n\nWater-eft 75 \n\n\n\nWeather 331 \n\n\n\nWheatear 57 \n\n\n\nWinter Birds of Passage 148 \n\n Wolmer Pond .... 37 \n\n Woodcocks . . . 167, 171 \n\n Wood-pigeons .... 139 \n\n\n\nWorms 242 \n\n\n\nWren, Golden-crowned . 69 \n\n Wren, Willow .... 66 \n\n\n\nYellow-hammer . . . 153 \n\n\n\n"
```

With this text data, you can do text mining of the content to answer cool questions you have. 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_namelist.R
\name{bhl_namelist}
\alias{bhl_namelist}
\title{List the unique names}
\usage{
bhl_namelist(...)
}
\description{
List the unique names
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getauthormetadata.R
\name{bhl_getauthormetadata}
\alias{bhl_getauthormetadata}
\title{Get metadata about an author}
\usage{
bhl_getauthormetadata(
  id,
  id_type = "bhl",
  pubs = FALSE,
  as = "table",
  key = NULL,
  ...
)
}
\arguments{
\item{id}{(numeric) the identifier of an individual author}

\item{id_type}{(character) the type of identifier (bhl, biostor, viaf).
Default: "bhl"}

\item{pubs}{(logical) \code{TRUE} to return the subject's publications.
Default: \code{FALSE}}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Get metadata about an author
}
\details{
You may choose to include a list of the author's publications.
The BHLType element identifies the type of each publication (Title or Part).
}
\examples{
\dontrun{
bhl_getauthormetadata(87509)
bhl_getauthormetadata(87509, pubs = TRUE, verbose = TRUE)
bhl_getauthormetadata(87509, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_subjectsearch.R
\name{bhl_subjectsearch}
\alias{bhl_subjectsearch}
\title{Return a list of subjects that match (fully or partially) the specified
search string.}
\usage{
bhl_subjectsearch(subject, as = "table", key = NULL, ...)
}
\arguments{
\item{subject}{the full or partial subject for which to search (character)}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Return a list of subjects that match (fully or partially) the specified
search string.
}
\examples{
\dontrun{
bhl_subjectsearch('diptera')
bhl_subjectsearch('diptera', "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getlanguages.R
\name{bhl_getlanguages}
\alias{bhl_getlanguages}
\title{Get a list of languages in which materials in BHL have been written.}
\usage{
bhl_getlanguages(as = "table", key = NULL, ...)
}
\arguments{
\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Get a list of languages in which materials in BHL have been written.
}
\examples{
\dontrun{
bhl_getlanguages()
bhl_getlanguages('json')
bhl_getlanguages('xml')
bhl_getlanguages('list')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getinstitutions.R
\name{bhl_getinstitutions}
\alias{bhl_getinstitutions}
\title{GetInstitutions data}
\usage{
bhl_getinstitutions(as = "table", key = NULL, ...)
}
\arguments{
\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
GetInstitutions data
}
\details{
Get a list of institutions which have contributed to BHL in some
way. These institutions may fill roles such as contributors, rights
holders, and scanning institutions.
}
\examples{
\dontrun{
bhl_getinstitutions()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_gettitleitems.R
\name{bhl_gettitleitems}
\alias{bhl_gettitleitems}
\title{Return a list of a title's items (books)}
\usage{
bhl_gettitleitems(...)
}
\description{
Return a list of a title's items (books)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getunpublisheditems.R
\name{bhl_getunpublisheditems}
\alias{bhl_getunpublisheditems}
\title{Return a list of the identifiers of all unpublished items}
\usage{
bhl_getunpublisheditems(...)
}
\description{
Return a list of the identifiers of all unpublished items
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getpages.R
\name{bhl_pages}
\alias{bhl_pages}
\title{Get many OCR-generated pages given a single item id}
\usage{
bhl_pages(itemid, key = NULL, ...)
}
\arguments{
\item{itemid}{the item id (character). Required}

\item{key}{Your BHL API key, either enter, or loads from your
\code{.Renviron} as \code{BHL_KEY} or from \code{.Rprofile} as
\code{bhl_key}}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Get many OCR-generated pages given a single item id
}
\examples{
\dontrun{
x <- bhl_publicationsearch('Selborne', year=1825, collection=4,
  language='eng')
bhl_pages(itemid = x$ItemID[1])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getitembyidentifier.R
\name{bhl_getitembyidentifier}
\alias{bhl_getitembyidentifier}
\title{Get item by identifier}
\usage{
bhl_getitembyidentifier(...)
}
\description{
Get item by identifier
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbhl-package.R
\docType{package}
\name{rbhl-package}
\alias{rbhl-package}
\alias{rbhl}
\title{rbhl}
\description{
R interface to the Biodiversity Heritage Library API.

You need an API key to use the Biodiversity Heritage Library API. Get your
BHL API key at \url{http://www.biodiversitylibrary.org/getapikey.aspx}.
Put your API key either in your .Renviron file (or similar) with the name
\code{BHL_KEY}, or in your .Rprofile file with the name \code{bhl_key}, and the
we'll use your API key without you having to enter
it every time you run a search.

See \link{rbhl-defunct} for defunct functions.
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_namesearch.R
\name{bhl_namesearch}
\alias{bhl_namesearch}
\title{Search for a particular name.}
\usage{
bhl_namesearch(name = NULL, as = "table", key = NULL, ...)
}
\arguments{
\item{name}{species name (character)}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Names both with and without NameBank identifiers are returned.
}
\examples{
\dontrun{
bhl_namesearch('poa annua')
bhl_namesearch(name='helianthus annuus')
bhl_namesearch(name='helianthus annuus', as='xml')
bhl_namesearch(name='helianthus annuus', as='json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_booksearch.R
\name{bhl_booksearch}
\alias{bhl_booksearch}
\title{Book search}
\usage{
bhl_booksearch(...)
}
\description{
Book search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getitempages.R
\name{bhl_getitempages}
\alias{bhl_getitempages}
\title{Return a list of an item's pages.}
\usage{
bhl_getitempages(...)
}
\description{
Return a list of an item's pages.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_openurl.R
\name{bhl_openurl}
\alias{bhl_openurl}
\title{Not sure how this differs from their other API...}
\usage{
bhl_openurl(
  genre = NULL,
  title = NULL,
  aufirst = NULL,
  aulast = NULL,
  date = NULL,
  spage = NULL,
  issue = NULL,
  version = 0.1,
  as = "list",
  key = NULL,
  ...
)
}
\arguments{
\item{genre}{Book genre}

\item{title}{Book title}

\item{aufirst}{First author}

\item{aulast}{Last author}

\item{date}{Date of publication}

\item{spage}{Start page}

\item{issue}{Issue number}

\item{version}{One of 0.1 or 1.0}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Not sure how this differs from their other API...
}
\examples{
\dontrun{
bhl_openurl(
	genre="book",
	title="Manual+of+North+American+Diptera",
	aufirst="Samuel Wendell",
	aulast="Williston",
	date=1908,
	spage=16)

bhl_openurl(genre="book", title="Manual+of+North+American+Diptera",
   aufirst="Samuel Wendell", aulast="Williston", date=1908, spage=16)

bhl_openurl(genre="book", title="Manual+of+North+American+Diptera",
   aufirst="Samuel Wendell", aulast="Williston", date=1908, spage=16,
   as='xml')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbhl-package.R
\name{rbhl-defunct}
\alias{rbhl-defunct}
\title{Defunct functions in rbhl}
\description{
\itemize{
\item \code{bhl_getpartendnote()}: BHL removed this API method.
\item \code{bhl_gettitleendnote()}: BHL removed this API method.
\item \code{\link[=bhl_booksearch]{bhl_booksearch()}}: method removed in API v3, see
\code{\link[=bhl_publicationsearch]{bhl_publicationsearch()}}
\item \code{\link[=bhl_partsearch]{bhl_partsearch()}}: method removed in API v3, see
\code{\link[=bhl_publicationsearch]{bhl_publicationsearch()}}
\item \code{\link[=bhl_titlesearchsimple]{bhl_titlesearchsimple()}}: method removed in API v3, see
\code{\link[=bhl_publicationsearch]{bhl_publicationsearch()}}
\item \code{\link[=bhl_getauthorparts]{bhl_getauthorparts()}}: method removed in API v3, see
\code{\link[=bhl_getauthormetadata]{bhl_getauthormetadata()}}
\item \code{\link[=bhl_getauthortitles]{bhl_getauthortitles()}}: method removed in API v3, see
\code{\link[=bhl_getauthormetadata]{bhl_getauthormetadata()}}
\item \code{\link[=bhl_getitembyidentifier]{bhl_getitembyidentifier()}}: method removed in API v3, see
\code{\link[=bhl_getitemmetadata]{bhl_getitemmetadata()}}
\item \code{\link[=bhl_getitempages]{bhl_getitempages()}}: method removed in API v3, see
\code{\link[=bhl_getitemmetadata]{bhl_getitemmetadata()}}
\item \code{\link[=bhl_getitemparts]{bhl_getitemparts()}}: method removed in API v3, see
\code{\link[=bhl_getitemmetadata]{bhl_getitemmetadata()}}
\item \code{\link[=bhl_getpagenames]{bhl_getpagenames()}}: method removed in API v3, see
\code{\link[=bhl_getpagemetadata]{bhl_getpagemetadata()}}
\item \code{\link[=bhl_getpageocrtext]{bhl_getpageocrtext()}}: method removed in API v3, see
\code{\link[=bhl_getpagemetadata]{bhl_getpagemetadata()}}
\item \code{\link[=bhl_getpartbyidentifier]{bhl_getpartbyidentifier()}}: method removed in API v3, see
\code{\link[=bhl_getpartmetadata]{bhl_getpartmetadata()}}
\item \code{\link[=bhl_getpartnames]{bhl_getpartnames()}}: method removed in API v3, see
\code{\link[=bhl_getpartmetadata]{bhl_getpartmetadata()}}
\item \code{\link[=bhl_getsubjectparts]{bhl_getsubjectparts()}}: method removed in API v3, see
\code{\link[=bhl_subjectmetadata]{bhl_subjectmetadata()}}
\item \code{\link[=bhl_getsubjecttitles]{bhl_getsubjecttitles()}}: method removed in API v3, see
\code{\link[=bhl_subjectmetadata]{bhl_subjectmetadata()}}
\item \code{\link[=bhl_getpartbibtex]{bhl_getpartbibtex()}}: method removed in API v3
\item \code{\link[=bhl_gettitlebibTex]{bhl_gettitlebibTex()}}: method removed in API v3
\item \code{\link[=bhl_gettitlebyidentifier]{bhl_gettitlebyidentifier()}}: API method merged in API v3, see
\code{\link[=bhl_gettitlemetadata]{bhl_gettitlemetadata()}}
\item \code{\link[=bhl_gettitleitems]{bhl_gettitleitems()}}: API method merged in API v3, see
\code{\link[=bhl_gettitlemetadata]{bhl_gettitlemetadata()}}
\item \code{\link[=bhl_getunpublisheditems]{bhl_getunpublisheditems()}}: method removed in API v3
\item \code{\link[=bhl_getunpublishedparts]{bhl_getunpublishedparts()}}: method removed in API v3
\item \code{\link[=bhl_getunpublishedtitles]{bhl_getunpublishedtitles()}}: method removed in API v3
\item \code{\link[=bhl_namecount]{bhl_namecount()}}: method removed in API v3
\item \code{\link[=bhl_namelist]{bhl_namelist()}}: method removed in API v3
\item \code{bhl_namegetdetail()}: method renamed - see \code{\link[=bhl_namemetadata]{bhl_namemetadata()}}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_gettitlebibtex.R
\name{bhl_gettitlebibTex}
\alias{bhl_gettitlebibTex}
\title{Return a citation for a title, using the BibTeX format}
\usage{
bhl_gettitlebibTex(...)
}
\description{
Return a citation for a title, using the BibTeX format
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getpartbyidentifier.R
\name{bhl_getpartbyidentifier}
\alias{bhl_getpartbyidentifier}
\title{Return a list of the identifiers of all unpublished items}
\usage{
bhl_getpartbyidentifier(...)
}
\description{
Return a list of the identifiers of all unpublished items
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getsubjecttitles.R
\name{bhl_getsubjecttitles}
\alias{bhl_getsubjecttitles}
\title{Return a list of titles associated with a subject}
\usage{
bhl_getsubjecttitles(...)
}
\description{
Return a list of titles associated with a subject
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_gettitlebyidentifier.R
\name{bhl_gettitlebyidentifier}
\alias{bhl_gettitlebyidentifier}
\title{Find and return metadata about a title or titles that match a specific
identifier}
\usage{
bhl_gettitlebyidentifier(...)
}
\description{
Find and return metadata about a title or titles that match a specific
identifier
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getauthorparts.R
\name{bhl_getauthorparts}
\alias{bhl_getauthorparts}
\title{Author parts}
\usage{
bhl_getauthorparts(...)
}
\description{
Author parts
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getunpublishedparts.R
\name{bhl_getunpublishedparts}
\alias{bhl_getunpublishedparts}
\title{Return a list of the identifiers of all unpublished parts
(articles, chapters, etc)}
\usage{
bhl_getunpublishedparts(...)
}
\description{
Return a list of the identifiers of all unpublished parts
(articles, chapters, etc)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getpageocrtext.R
\name{bhl_getpageocrtext}
\alias{bhl_getpageocrtext}
\title{Return the OCR-generated text of a page.}
\usage{
bhl_getpageocrtext(...)
}
\description{
Return the OCR-generated text of a page.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getpagenames.R
\name{bhl_getpagenames}
\alias{bhl_getpagenames}
\title{Return a list of names that appear on a page.}
\usage{
bhl_getpagenames(...)
}
\description{
Return a list of names that appear on a page.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getpartmetadata.R
\name{bhl_getpartmetadata}
\alias{bhl_getpartmetadata}
\title{Return a list of an item's pages.}
\usage{
bhl_getpartmetadata(partid, key = NULL, ...)
}
\arguments{
\item{partid}{The identifier of an individual part (article, chapter, etc)
(numeric)}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Return a list of an item's pages.
}
\examples{
\dontrun{
bhl_getpartmetadata(10409)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getitemmetadata.R
\name{bhl_getitemmetadata}
\alias{bhl_getitemmetadata}
\title{Return metadata about an item.}
\usage{
bhl_getitemmetadata(
  itemid = NULL,
  pages = TRUE,
  ocr = FALSE,
  parts = FALSE,
  as = "table",
  key = NULL,
  ...
)
}
\arguments{
\item{itemid}{item id (character)}

\item{pages}{return the items pages (\code{TRUE}/\code{FALSE})}

\item{ocr}{(logical) \code{TRUE} to return the ocr for the item's pages.
Setting this to \code{TRUE} apparently doesn't return any actual ocr text,
but leaving parameter here for now.}

\item{parts}{(logical) \code{TRUE} to return the item's parts. Setting this
to \code{TRUE} apparently doesn't return any parts text, but leaving
parameter here for now.}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
You may choose to include a list of the item's pages.
}
\examples{
\dontrun{
bhl_getitemmetadata('16800', TRUE)
bhl_getitemmetadata('16800', TRUE, as='xml')
bhl_getitemmetadata('16800', TRUE, as='json')
bhl_getitemmetadata('16800', TRUE, as='list')
bhl_getitemmetadata(20419, pages=FALSE, parts=TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_gettitlemetadata.R
\name{bhl_gettitlemetadata}
\alias{bhl_gettitlemetadata}
\title{Get title metadata}
\usage{
bhl_gettitlemetadata(titleid = NA, items = FALSE, as = "list", key = NULL, ...)
}
\arguments{
\item{titleid}{the identifier of an individual title (numeric)}

\item{items}{(logical) TRUE of FALSE (default) to inclue items}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Return metadata about a title. You may choose to include a list of
the items (books) associated  with the title.
}
\examples{
\dontrun{
bhl_gettitlemetadata(1726, items = TRUE)
bhl_gettitlemetadata(1726, as='list')
bhl_gettitlemetadata(1726, as='xml')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_namecount.R
\name{bhl_namecount}
\alias{bhl_namecount}
\title{Return the number of unique names found on pages in BHL}
\usage{
bhl_namecount(...)
}
\description{
Return the number of unique names found on pages in BHL
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_subjectmetadata.R
\name{bhl_subjectmetadata}
\alias{bhl_subjectmetadata}
\title{Get metadata about a subject}
\usage{
bhl_subjectmetadata(subject, pubs = FALSE, as = "table", key = NULL, ...)
}
\arguments{
\item{subject}{(character) he subject for which to return metadata}

\item{pubs}{(logical) \code{TRUE} to return the subject's publications.
Default: \code{FALSE}}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Get metadata about a subject
}
\details{
You may choose to include a list of the subject's publications.
The BHLType element identifies the type of each publication (Title or Part).
}
\examples{
\dontrun{
bhl_subjectmetadata(subject = "water")
x <- bhl_subjectmetadata(subject = "water", pubs = TRUE)
head(x$Publications[[1]])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_publicationsearch.R
\name{bhl_publicationsearchadv}
\alias{bhl_publicationsearchadv}
\title{Publication search advanced}
\usage{
bhl_publicationsearchadv(
  title = NULL,
  titleop = NULL,
  authorname = NULL,
  year = NULL,
  subject = NULL,
  language = NULL,
  collection = NULL,
  notes = NULL,
  notesop = NULL,
  text = NULL,
  textop = NULL,
  page = NULL,
  as = "table",
  key = NULL,
  ...
)
}
\arguments{
\item{title}{(character) a title for which to search}

\item{titleop}{(character) 'all' to search for all specified words in the
title fields; 'phrase' to search for the exact phrase specified}

\item{authorname}{(character) an author name for which to search}

\item{year}{(integer) a four-digit publication year for which to search}

\item{subject}{(character) a subject for which to search}

\item{language}{(character) a language code; search will only return
publications in the specified language}

\item{collection}{(character) a collection id; search will only return
publications from the specfied collection}

\item{notes}{(character) one or more words for which to search in the
publication notes}

\item{notesop}{(character) 'all' to search for all specified words in the
notes field; 'phrase' to search for the exact phrase specified}

\item{text}{(character) one or more words for which to search in the text
of the publications}

\item{textop}{(character) 'all' to search for all specified words in the
text field; 'phrase' to search for the exact phrase specified}

\item{page}{(integer) 1 for the first page of results, 2 for the second,
and so on}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Publication search advanced
}
\examples{
\dontrun{
bhl_publicationsearchadv(title = "cocos island", titleop = "all",
  authorname = "gifford", year = 1919)
bhl_publicationsearchadv(title = "cocos island", authorname = "gifford")
}
}
\seealso{
Other publication-search: 
\code{\link{bhl_publicationsearch}()}
}
\concept{publication-search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getauthortitles.R
\name{bhl_getauthortitles}
\alias{bhl_getauthortitles}
\title{Author titles}
\usage{
bhl_getauthortitles(...)
}
\description{
Author titles
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getpartbibtex.R
\name{bhl_getpartbibtex}
\alias{bhl_getpartbibtex}
\title{Return a citation for a part, using the BibTeX format.}
\usage{
bhl_getpartbibtex(...)
}
\description{
Return a citation for a part, using the BibTeX format.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_authorsearch.R
\name{bhl_authorsearch}
\alias{bhl_authorsearch}
\title{Author search}
\usage{
bhl_authorsearch(name = NULL, as = "table", key = NULL, ...)
}
\arguments{
\item{name}{full or partial name of the author for which to search
(last name listed first, i.e. 'Darwin, Charles') (character)}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Return a list of authors that match (fully or partially)
the specified search string.
}
\details{
The namessearched are those contained in MARC 100a, 110a,
111a, 700a, 710a, and 711a library records.
}
\examples{
\dontrun{
bhl_authorsearch(name='dimmock')
bhl_authorsearch(name='Jones')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getitemparts.R
\name{bhl_getitemparts}
\alias{bhl_getitemparts}
\title{Return a list of an item's parts}
\usage{
bhl_getitemparts(...)
}
\description{
Return a list of an item's parts
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_publicationsearch.R
\name{bhl_publicationsearch}
\alias{bhl_publicationsearch}
\title{Publication search}
\usage{
bhl_publicationsearch(
  q = NULL,
  full_text = FALSE,
  page = NULL,
  as = "table",
  key = NULL,
  ...
)
}
\arguments{
\item{q}{(character) the text for which to search}

\item{full_text}{(character) \code{TRUE} for a catalog plus full-text search;
\code{FALSE} for a catalog-only search. default: \code{FALSE}}

\item{page}{(integer) page to return}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Publication search
}
\examples{
\dontrun{
bhl_publicationsearch('cocos island costa rica birds')
bhl_publicationsearch('poa annua')
bhl_publicationsearch('helianthus annuus')
}
}
\seealso{
Other publication-search: 
\code{\link{bhl_publicationsearchadv}()}
}
\concept{publication-search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getcollections.R
\name{bhl_getcollections}
\alias{bhl_getcollections}
\title{Get a list of collections which are used to group titles and items. A single
collection may contain either titles or items, but not both.}
\usage{
bhl_getcollections(as = "table", key = NULL, ...)
}
\arguments{
\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Get a list of collections which are used to group titles and items. A single
collection may contain either titles or items, but not both.
}
\examples{
\dontrun{
bhl_getcollections()
bhl_getcollections(as = 'list')
bhl_getcollections(as = 'json')
bhl_getcollections(as = 'xml')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_pagesearch.R
\name{bhl_pagesearch}
\alias{bhl_pagesearch}
\title{Search an item for pages containing the specified text}
\usage{
bhl_pagesearch(id, text, as = "table", key = NULL, ...)
}
\arguments{
\item{id}{(integer) BHL identifier of the item to be searched}

\item{text}{(character) the text for which to search}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Search an item for pages containing the specified text
}
\examples{
\dontrun{
bhl_pagesearch(22004, "dog")
bhl_pagesearch(22004, "dog", as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_namegetdetail.R
\name{bhl_namemetadata}
\alias{bhl_namemetadata}
\title{Get basic title, item, and page metadata for each page on which the
specified name appears.}
\usage{
bhl_namemetadata(
  namebankid = NULL,
  name = NULL,
  id_type = NULL,
  as = "list",
  key = NULL,
  ...
)
}
\arguments{
\item{namebankid}{(numeric) (not used if 'name' specified) NameBank
identifier for a name}

\item{name}{(character) (not used if 'namebankid' specified) a
complete name string}

\item{id_type}{(character) the type of identifier (namebank, eol, gni,
ion, col, gbif, itis, ipni, worms). Not used if name is specified}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Get basic title, item, and page metadata for each page on which the
specified name appears.
}
\examples{
\dontrun{
res <- bhl_namemetadata(name = 'poa annua linnaeus')
str(res)
# bhl_namemetadata(namebankid = 3501464)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getunpublishedtitles.R
\name{bhl_getunpublishedtitles}
\alias{bhl_getunpublishedtitles}
\title{Return a list of the identifiers of all unpublished titles}
\usage{
bhl_getunpublishedtitles(...)
}
\description{
Return a list of the identifiers of all unpublished titles
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbhl-package.R
\docType{data}
\name{rbhlmethods}
\alias{rbhlmethods}
\title{Data.frame of all the BHL API methods from the BHL website.}
\description{
Data.frame of all the BHL API methods from the BHL website.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getpartnames.R
\name{bhl_getpartnames}
\alias{bhl_getpartnames}
\title{Return a list of scientific names associated with a part}
\usage{
bhl_getpartnames(...)
}
\description{
Return a list of scientific names associated with a part
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getsubjectparts.R
\name{bhl_getsubjectparts}
\alias{bhl_getsubjectparts}
\title{Return a list of parts (articles, chapters, etc) associated with a subject}
\usage{
bhl_getsubjectparts(...)
}
\description{
Return a list of parts (articles, chapters, etc) associated with a subject
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_getpagemetadata.R
\name{bhl_getpagemetadata}
\alias{bhl_getpagemetadata}
\title{Return metadata about a page.}
\usage{
bhl_getpagemetadata(
  page = NULL,
  ocr = FALSE,
  names = FALSE,
  as = "table",
  key = NULL,
  ...
)
}
\arguments{
\item{page}{page number to get}

\item{ocr}{return ocr text of the page (\code{TRUE}/\code{FALSE})}

\item{names}{return the names that appear on the page (\code{TRUE}/\code{FALSE})}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
You may choose to include the OCR text and a list of names found on the page
}
\examples{
\dontrun{
bhl_getpagemetadata(page=1328690, ocr=TRUE)
bhl_getpagemetadata(page=1328690, ocr=TRUE, as='json')
bhl_getpagemetadata(page=1328690, ocr=TRUE, as='xml')
bhl_getpagemetadata(page=1328690, ocr=TRUE, as='list')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_titlesearchsimple.R
\name{bhl_titlesearchsimple}
\alias{bhl_titlesearchsimple}
\title{Perform a simple title search}
\usage{
bhl_titlesearchsimple(...)
}
\description{
Perform a simple title search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_partsearch.R
\name{bhl_partsearch}
\alias{bhl_partsearch}
\title{Search for parts of books}
\usage{
bhl_partsearch(...)
}
\description{
Search for parts of books
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhl_bioherlib.R
\name{bhl_bioherlib}
\alias{bhl_bioherlib}
\title{Search BHL across many API methods.}
\usage{
bhl_bioherlib(
  method = "GetPageMetadata",
  pageid = NULL,
  ocr = FALSE,
  names = FALSE,
  as = "table",
  key = NULL,
  ...
)
}
\arguments{
\item{method}{The API method to use.}

\item{pageid}{The identifier of an individual page in a scanned book.}

\item{ocr}{return ocr text of the page (logical). Default: \code{FALSE}}

\item{names}{return the names that appear on the page (logical).
Default: \code{FALSE}}

\item{as}{(character) Return a list ("list"), json ("json"), xml ("xml"),
or parsed table ("table", default). Note that \code{as="table"} can give
different data format back depending on the function - for example,
sometimes a data.frame and sometimes a character vector.}

\item{key}{Your BHL API key, either enter, or loads from your \code{.Renviron}
as \code{BHL_KEY}
or from \code{.Rprofile} as \code{bhl_key}.}

\item{...}{Curl options passed on to \code{\link[crul:HttpClient]{crul::HttpClient()}}}
}
\description{
Search BHL across many API methods.
}
\examples{
\dontrun{
bhl_bioherlib(method='GetPageMetadata', pageid=1328690, ocr=TRUE,
  names=TRUE)
bhl_bioherlib(method='GetPageMetadata', pageid=1328690, ocr=TRUE,
  names=TRUE, as="xml")
bhl_bioherlib(method='GetPageMetadata', pageid=1328690, ocr=TRUE,
  names=TRUE, as="list")
}
}
