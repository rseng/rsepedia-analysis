

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/oai/actions/workflows/R-check.yml/badge.svg)](https://github.com/ropensci/oai/actions/workflows/R-check.yml)
[![cran checks](https://cranchecks.info/badges/worst/oai)](https://cranchecks.info/pkgs/oai)
[![codecov.io](https://codecov.io/github/ropensci/oai/coverage.svg?branch=master)](https://codecov.io/github/ropensci/oai?branch=master) 
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/oai?color=2ED968)](https://github.com/r-hub/cranlogs.app) 
[![cran version](https://www.r-pkg.org/badges/version/oai)](https://cran.r-project.org/package=oai) 
[![](https://badges.ropensci.org/19_status.svg)](https://github.com/ropensci/software-review/issues/19)

`oai` is an R client to work with OAI-PMH (Open Archives Initiative Protocol for Metadata Harvesting) services, a protocol developed by the Open Archives Initiative (https://en.wikipedia.org/wiki/Open_Archives_Initiative). OAI-PMH uses XML data format transported over HTTP.

OAI-PMH Info:

* Wikipedia (https://en.wikipedia.org/wiki/Open_Archives_Initiative_Protocol_for_Metadata_Harvesting)
* OAI V2 specification (http://www.openarchives.org/OAI/openarchivesprotocol.html)

`oai` is built on `xml2` and `httr`. In addition, we give back data.frame's whenever possible to make data comprehension, manipulation, and visualization easier. We also have functions to fetch a large directory of OAI-PMH services - it isn't exhaustive, but does contain a lot.

OAI-PMH instead of paging with e.g., `page` and `per_page` parameters, uses (optionally) `resumptionTokens`, optionally with an expiration date. These tokens can be used to continue on to the next chunk of data, if the first request did not get to the end. Often, OAI-PMH services limit each request to 50 records, but this may vary by provider, I don't know for sure. The API of this package is such that we `while` loop for you internally until we get all records. We may in the future expose e.g., a `limit` parameter so you can say how many records you want, but we haven't done this yet.

## Install

Install from CRAN


```r
install.packages("oai")
```

Development version


```r
devtools::install_github("ropensci/oai")
```


```r
library("oai")
```

## Identify


```r
id("http://oai.datacite.org/oai")
#>   repositoryName                      baseURL protocolVersion
#> 1       DataCite https://oai.datacite.org/oai             2.0
#>             adminEmail    earliestDatestamp deletedRecord          granularity
#> 1 support@datacite.org 2011-01-01T00:00:00Z    persistent YYYY-MM-DDThh:mm:ssZ
#>   compression compression.1                                    description
#> 1        gzip       deflate oaioai.datacite.org:oai:oai.datacite.org:12425
```

## ListIdentifiers


```r
list_identifiers(from = '2018-05-01T', until = '2018-06-01T')
#> # A tibble: 85 x 5
#>    identifier         datestamp     setSpec             setSpec.1      setSpec.2
#>    <chr>              <chr>         <chr>               <chr>          <chr>    
#>  1 cf7fbc99-de82-41a… 2018-05-31T1… installation:791e3… dataset_type:… country:…
#>  2 09d5405e-ca86-45f… 2018-05-30T1… installation:804b8… dataset_type:… country:…
#>  3 4b64d1f2-31c2-40c… 2018-05-30T1… installation:804b8… dataset_type:… country:…
#>  4 884378d6-d591-476… 2018-05-29T1… installation:a1650… dataset_type:… country:…
#>  5 18799ce9-1a66-40f… 2018-05-14T1… installation:d1b0a… dataset_type:… country:…
#>  6 7e91aacb-c994-41e… 2018-05-21T1… installation:d5b61… dataset_type:… country:…
#>  7 f83746ee-4cf2-4e6… 2018-05-08T0… installation:c4195… dataset_type:… country:…
#>  8 a3533a61-6f88-443… 2018-05-08T1… installation:06d75… dataset_type:… country:…
#>  9 ba9b66a3-2d11-419… 2018-05-05T2… installation:d1b0a… dataset_type:… country:…
#> 10 78b696d9-8f0d-41a… 2018-05-05T2… installation:d1b0a… dataset_type:… country:…
#> # … with 75 more rows
```

## Count Identifiers


```r
count_identifiers()
#>                            url   count
#> 1 http://export.arxiv.org/oai2 1586724
```

## ListRecords


```r
list_records(from = '2018-05-01T', until = '2018-05-15T')
#> # A tibble: 42 x 26
#>    identifier datestamp setSpec setSpec.1 setSpec.2 title publisher identifier.1
#>    <chr>      <chr>     <chr>   <chr>     <chr>     <chr> <chr>     <chr>       
#>  1 18799ce9-… 2018-05-… instal… dataset_… country:… Bird… Sokoine … https://www…
#>  2 f83746ee-… 2018-05-… instal… dataset_… country:… NDFF… Dutch Na… https://www…
#>  3 a3533a61-… 2018-05-… instal… dataset_… country:… EDP … EDP - En… https://www…
#>  4 ba9b66a3-… 2018-05-… instal… dataset_… country:… Ende… Sokoine … https://www…
#>  5 78b696d9-… 2018-05-… instal… dataset_… country:… Ende… Sokoine … https://www…
#>  6 c791b255-… 2018-05-… instal… dataset_… country:… Ende… Sokoine … https://www…
#>  7 b929ccda-… 2018-05-… instal… dataset_… country:… List… Sokoine … https://www…
#>  8 da285c2a-… 2018-05-… instal… dataset_… country:… Moni… Corporac… https://www…
#>  9 87372877-… 2018-05-… instal… dataset_… country:… Moni… Corporac… https://www…
#> 10 ed7d4c25-… 2018-05-… instal… dataset_… country:… Samo… Ministry… https://www…
#> # … with 32 more rows, and 18 more variables: subject <chr>, source <chr>,
#> #   description <chr>, description.1 <chr>, type <chr>, creator <chr>,
#> #   date <chr>, language <chr>, coverage <chr>, coverage.1 <chr>, format <chr>,
#> #   source.1 <chr>, subject.1 <chr>, creator.1 <chr>, coverage.2 <chr>,
#> #   description.2 <chr>, creator.2 <chr>, subject.2 <chr>
```

## GetRecords


```r
ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", "d981c07d-bc43-40a2-be1f-e786e25106ac")
get_records(ids)
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$header
#> # A tibble: 1 x 3
#>   identifier              datestamp      setSpec                                
#>   <chr>                   <chr>          <chr>                                  
#> 1 87832186-00ea-44dd-a6b… 2018-06-29T12… installation:729a7375-b120-4e4f-bb81-a…
#> 
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$metadata
#> # A tibble: 0 x 0
#> 
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$header
#> # A tibble: 1 x 3
#>   identifier              datestamp      setSpec                                
#>   <chr>                   <chr>          <chr>                                  
#> 1 d981c07d-bc43-40a2-be1… 2018-01-21T21… installation:804b8dd0-07ac-4a30-bf92-3…
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$metadata
#> # A tibble: 1 x 12
#>   title  publisher  identifier  subject  source  description type  creator date 
#>   <chr>  <chr>      <chr>       <chr>    <chr>   <chr>       <chr> <chr>   <chr>
#> 1 Peces… Instituto… https://ww… peces, … http:/… Caracteriz… Data… Fernan… 2018…
#> # … with 3 more variables: language <chr>, coverage <chr>, format <chr>
```

## List MetadataFormats


```r
list_metadataformats(id = "87832186-00ea-44dd-a6bf-c2896c4d09b4")
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#>   metadataPrefix                                                   schema
#> 1         oai_dc           http://www.openarchives.org/OAI/2.0/oai_dc.xsd
#> 2            eml http://rs.gbif.org/schema/eml-gbif-profile/1.0.2/eml.xsd
#>                             metadataNamespace
#> 1 http://www.openarchives.org/OAI/2.0/oai_dc/
#> 2          eml://ecoinformatics.org/eml-2.1.1
```

## List Sets


```r
list_sets("http://api.gbif.org/v1/oai-pmh/registry")
#> # A tibble: 597 x 2
#>    setSpec                     setName         
#>    <chr>                       <chr>           
#>  1 dataset_type                per dataset type
#>  2 dataset_type:OCCURRENCE     occurrence      
#>  3 dataset_type:CHECKLIST      checklist       
#>  4 dataset_type:METADATA       metadata        
#>  5 dataset_type:SAMPLING_EVENT sampling_event  
#>  6 country                     per country     
#>  7 country:AD                  Andorra         
#>  8 country:AM                  Armenia         
#>  9 country:AO                  Angola          
#> 10 country:AR                  Argentina       
#> # … with 587 more rows
```

## Examples of other OAI providers

### Biodiversity Heritage Library

Identify


```r
id("http://www.biodiversitylibrary.org/oai")
#>                                 repositoryName
#> 1 Biodiversity Heritage Library OAI Repository
#>                                   baseURL protocolVersion
#> 1 https://www.biodiversitylibrary.org/oai             2.0
#>                    adminEmail earliestDatestamp deletedRecord granularity
#> 1 oai@biodiversitylibrary.org        2006-01-01            no  YYYY-MM-DD
#>                                                        description
#> 1 oaibiodiversitylibrary.org:oai:biodiversitylibrary.org:item/1000
```

Get records


```r
get_records(c("oai:biodiversitylibrary.org:item/7", "oai:biodiversitylibrary.org:item/9"),
            url = "http://www.biodiversitylibrary.org/oai")
#> $`oai:biodiversitylibrary.org:item/7`
#> $`oai:biodiversitylibrary.org:item/7`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/7 2016-01-26T06:05:19Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/7`$metadata
#> # A tibble: 1 x 10
#>   title   creator  subject  description  publisher contributor type  identifier 
#>   <chr>   <chr>    <chr>    <chr>        <chr>     <chr>       <chr> <chr>      
#> 1 Die Mu… Fleisch… Bogor;I… pt.5:v.1 (1… Leiden :… Missouri B… text… https://ww…
#> # … with 2 more variables: language <chr>, rights <chr>
#> 
#> 
#> $`oai:biodiversitylibrary.org:item/9`
#> $`oai:biodiversitylibrary.org:item/9`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/9 2016-01-26T06:05:19Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/9`$metadata
#> # A tibble: 1 x 10
#>   title   creator  subject  description  publisher contributor type  identifier 
#>   <chr>   <chr>    <chr>    <chr>        <chr>     <chr>       <chr> <chr>      
#> 1 Die Mu… Fleisch… Bogor;I… pt.5:v.3 (1… Leiden :… Missouri B… text… https://ww…
#> # … with 2 more variables: language <chr>, rights <chr>
```


## Acknowledgements

Michal Bojanowski thanks National Science Centre for support through grant 2012/07/D/HS6/01971.


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/oai/issues).
* License: MIT
* Get citation information for `oai` in R doing `citation(package = 'oai')`
* Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.
oai 0.3.2
=========

### MINOR IMPROVEMENTS

* vignette fix, use markdown in suggests
* update readme, use ropensci coc

oai 0.3.0
=========

### NEW FEATURES

* `id()` gains `as` parameter so user can ask for different outputs (parsed list/data.frame, or raw xml/text) (#54)
* most OAI functions changed their default url to `http://api.gbif.org/v1/oai-pmh/registry`, while `count_identifiers()` changed it's default url to `http://export.arxiv.org/oai2`. the previous default url for Datacite was too unreliable (was often unresponsive)

### MINOR IMPROVEMENTS

* add grant information for one author (#48) (#49)
* code of conduct urls fixed
* now using markdown supported documentation (#56)
* replace `tibble::as_data_frame` with `tibble::as_tibble` throughout package

### BUG FIXES

* fix to `update_providers()`; html page that we scrape had changed (#57)


oai 0.2.2
=========

### NEW FEATURES

* Added new parsers in `get_records()` specific to different OAI prefixes. Currently has
parsers for `oai_dc` and `oai_datacite`. For prefixes we don't have
parsers for we return raw XML so you can parse it yourself. (#45)

### MINOR IMPROVEMENTS

* Replace `xml2::xml_find_one()` with `xml2::xml_find_first()` (#39)
* Update URLs in `DESCRIPTION` file (#43)
* Using `tibble` now for compact data.frame instead of internal
methods for the same (#44)
* `as` parameter in `get_records()` now has options `parsed` or `raw`, which
replaces `df` `list`, or `raw`


oai 0.2.0
=========

### NEW FEATURES

* A set of new functions for dealing with larger data results:
`dump_raw_to_txt()`, `dump_to_rds()`, and `dump_raw_to_db()`.
They can be used with `oai` functions `list_identifiers()`, `list_sets()`,
and `list_records()` (#9) (#15) (#21) thanks @mbojan

### MINOR IMPROVEMENTS

* Sped up some tests (#19)
* Better description of OAI protocol in the `DESCRIPTION` file (#28)
* Commented about where some internal functions come from (#31)
* Import `plyr` for `rbind.fill()` (#32)
* Including now some examples using OAI-PMH with GBIF and BHL (#33)
* Better error handling! (#10) (#12) (#22) (#27) thanks @mbojan

### BUG FIXES

* Fixed bug where `list_identifiers()` threw error when no result was found (#13)
* Dealing better with bad inputs to `as` parameter - stop with informative message
now instead of failing without anything returned (#34)

oai 0.1.0
=========

* Released to CRAN.
## Test environments

* local macOS install, R 4.0.5 patched
* ubuntu 14.04 (on GitHub Actions), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 5 downstream dependencies, with 
no problems related to this package. Results are at
<https://github.com/ropensci/oai/blob/master/revdep/README.md>

------

This version fixes the rmarkdown/markdown dependency issue for vignettes that Kurt emailed maintainers about.

Thanks!
Scott Chamberlain
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{oai introduction}
%\VignetteEncoding{UTF-8}
-->



oai introduction
================

A general purpose client to work with any 'OAI-PMH' service. 

The 'OAI-PMH' protocol is described at <http://www.openarchives.org/OAI/openarchivesprotocol.html>.

The main functions follow the OAI-PMH verbs:

* `GetRecord`
* `Identify`
* `ListIdentifiers`
* `ListMetadataFormats`
* `ListRecords`
* `ListSets`

## Get oai

Install from CRAN


```r
install.packages("oai")
```

Or install the development version from GitHub


```r
devtools::install_github("ropensci/oai")
```

Load `oai`


```r
library("oai")
```

## Identify


```r
id("http://oai.datacite.org/oai")
#>   repositoryName                      baseURL protocolVersion
#> 1   DataCite MDS https://oai.datacite.org/oai             2.0
#>             adminEmail    earliestDatestamp deletedRecord
#> 1 support@datacite.org 2011-01-01T00:00:00Z    persistent
#>            granularity compression compression.1
#> 1 YYYY-MM-DDThh:mm:ssZ        gzip       deflate
#>                                      description
#> 1 oaioai.datacite.org:oai:oai.datacite.org:12425
```

## ListIdentifiers


```r
list_identifiers(from = '2018-05-01T', until = '2018-09-01T')
#> # A tibble: 2,904 x 5
#>    identifier        datestamp    setSpec           setSpec.1     setSpec.2
#>    <chr>             <chr>        <chr>             <chr>         <chr>    
#>  1 cbe4cdda-2045-41… 2018-08-30T… installation:842… dataset_type… country:…
#>  2 3fb7ddd8-07c0-49… 2018-08-28T… installation:842… dataset_type… country:…
#>  3 34585f24-1ffe-47… 2018-08-28T… installation:842… dataset_type… country:…
#>  4 8ac24ec8-1e7b-40… 2018-08-28T… installation:842… dataset_type… country:…
#>  5 e381b970-9b62-46… 2018-08-28T… installation:842… dataset_type… country:…
#>  6 4496b1b1-1ea9-4e… 2018-08-28T… installation:842… dataset_type… country:…
#>  7 fc5d68a6-b511-4c… 2018-08-27T… installation:73e… dataset_type… country:…
#>  8 5e0073bc-de2a-4b… 2018-08-27T… installation:688… dataset_type… country:…
#>  9 ab2684cf-62e5-4f… 2018-08-27T… installation:804… dataset_type… country:…
#> 10 34fbcf59-d9bb-47… 2018-08-27T… installation:7c5… dataset_type… country:…
#> # … with 2,894 more rows
```

## Count Identifiers


```r
count_identifiers()
```

## ListRecords


```r
list_records(from = '2018-05-01T', until = '2018-05-15T')
#> # A tibble: 44 x 26
#>    identifier datestamp setSpec setSpec.1 setSpec.2 title publisher
#>    <chr>      <chr>     <chr>   <chr>     <chr>     <chr> <chr>    
#>  1 18799ce9-… 2018-05-… instal… dataset_… country:… Bird… Sokoine …
#>  2 79f51633-… 2018-05-… instal… dataset_… country:… Impl… Aïgos SAS
#>  3 f83746ee-… 2018-05-… instal… dataset_… country:… NDFF… Dutch Na…
#>  4 a3533a61-… 2018-05-… instal… dataset_… country:… EDP … EDP - En…
#>  5 ba9b66a3-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  6 78b696d9-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  7 c791b255-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  8 b929ccda-… 2018-05-… instal… dataset_… country:… List… Sokoine …
#>  9 da285c2a-… 2018-05-… instal… dataset_… country:… Moni… Corporac…
#> 10 87372877-… 2018-05-… instal… dataset_… country:… Moni… Corporac…
#> # … with 34 more rows, and 19 more variables: identifier.1 <chr>,
#> #   subject <chr>, source <chr>, description <chr>, description.1 <chr>,
#> #   type <chr>, creator <chr>, date <chr>, language <chr>, coverage <chr>,
#> #   coverage.1 <chr>, format <chr>, source.1 <chr>, subject.1 <chr>,
#> #   coverage.2 <chr>, creator.1 <chr>, description.2 <chr>,
#> #   creator.2 <chr>, subject.2 <chr>
```

## GetRecords


```r
ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", "d981c07d-bc43-40a2-be1f-e786e25106ac")
get_records(ids)
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$header
#> # A tibble: 1 x 3
#>   identifier             datestamp      setSpec                            
#>   <chr>                  <chr>          <chr>                              
#> 1 87832186-00ea-44dd-a6… 2018-06-29T12… installation:729a7375-b120-4e4f-bb…
#> 
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$metadata
#> # A tibble: 0 x 0
#> 
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$header
#> # A tibble: 1 x 3
#>   identifier             datestamp      setSpec                            
#>   <chr>                  <chr>          <chr>                              
#> 1 d981c07d-bc43-40a2-be… 2018-01-21T21… installation:804b8dd0-07ac-4a30-bf…
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$metadata
#> # A tibble: 1 x 12
#>   title publisher identifier subject source description type  creator date 
#>   <chr> <chr>     <chr>      <chr>   <chr>  <chr>       <chr> <chr>   <chr>
#> 1 Pece… Institut… https://w… peces,… http:… Caracteriz… Data… Fernan… 2018…
#> # … with 3 more variables: language <chr>, coverage <chr>, format <chr>
```

## List MetadataFormats


```r
list_metadataformats(id = "87832186-00ea-44dd-a6bf-c2896c4d09b4")
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#>   metadataPrefix                                                   schema
#> 1         oai_dc           http://www.openarchives.org/OAI/2.0/oai_dc.xsd
#> 2            eml http://rs.gbif.org/schema/eml-gbif-profile/1.0.2/eml.xsd
#>                             metadataNamespace
#> 1 http://www.openarchives.org/OAI/2.0/oai_dc/
#> 2          eml://ecoinformatics.org/eml-2.1.1
```

## List Sets


```r
list_sets("http://api.gbif.org/v1/oai-pmh/registry")
#> # A tibble: 572 x 2
#>    setSpec                     setName         
#>    <chr>                       <chr>           
#>  1 dataset_type                per dataset type
#>  2 dataset_type:OCCURRENCE     occurrence      
#>  3 dataset_type:CHECKLIST      checklist       
#>  4 dataset_type:METADATA       metadata        
#>  5 dataset_type:SAMPLING_EVENT sampling_event  
#>  6 country                     per country     
#>  7 country:AD                  Andorra         
#>  8 country:AO                  Angola          
#>  9 country:AR                  Argentina       
#> 10 country:AT                  Austria         
#> # … with 562 more rows
```

### Biodiversity Heritage Library

Identify


```r
id("http://www.biodiversitylibrary.org/oai")
#>                                 repositoryName
#> 1 Biodiversity Heritage Library OAI Repository
#>                                   baseURL protocolVersion
#> 1 https://www.biodiversitylibrary.org/oai             2.0
#>                    adminEmail earliestDatestamp deletedRecord granularity
#> 1 oai@biodiversitylibrary.org        2006-01-01            no  YYYY-MM-DD
#>                                                        description
#> 1 oaibiodiversitylibrary.org:oai:biodiversitylibrary.org:item/1000
```

Get records


```r
get_records(c("oai:biodiversitylibrary.org:item/7", "oai:biodiversitylibrary.org:item/9"),
            url = "http://www.biodiversitylibrary.org/oai")
#> $`oai:biodiversitylibrary.org:item/7`
#> $`oai:biodiversitylibrary.org:item/7`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/7 2016-07-13T09:13:41Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/7`$metadata
#> # A tibble: 1 x 11
#>   title creator subject description publisher contributor date  type 
#>   <chr> <chr>   <chr>   <chr>       <chr>     <chr>       <chr> <chr>
#> 1 Die … Fleisc… Bogor;… pt.5:v.1 (… Leiden :… Missouri B… 1900… text…
#> # … with 3 more variables: identifier <chr>, language <chr>, rights <chr>
#> 
#> 
#> $`oai:biodiversitylibrary.org:item/9`
#> $`oai:biodiversitylibrary.org:item/9`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/9 2016-07-13T09:13:41Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/9`$metadata
#> # A tibble: 1 x 11
#>   title creator subject description publisher contributor date  type 
#>   <chr> <chr>   <chr>   <chr>       <chr>     <chr>       <chr> <chr>
#> 1 Die … Fleisc… Bogor;… pt.5:v.3 (… Leiden :… Missouri B… 1906… text…
#> # … with 3 more variables: identifier <chr>, language <chr>, rights <chr>
```
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 3.6.1 Patched (2019-09-05 r77154) |
|os       |macOS Mojave 10.14.6                        |
|system   |x86_64, darwin15.6.0                        |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2019-09-06                                  |

# Dependencies

|package  |old   |new        |Δ  |
|:--------|:-----|:----------|:--|
|oai      |0.2.2 |0.2.9.9910 |*  |
|ellipsis |NA    |0.2.0.9000 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>",
  cache.path = "inst/cache/"
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/oai/actions/workflows/R-check.yml/badge.svg)](https://github.com/ropensci/oai/actions/workflows/R-check.yml)
[![cran checks](https://cranchecks.info/badges/worst/oai)](https://cranchecks.info/pkgs/oai)
[![codecov.io](https://codecov.io/github/ropensci/oai/coverage.svg?branch=master)](https://codecov.io/github/ropensci/oai?branch=master) 
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/oai?color=2ED968)](https://github.com/r-hub/cranlogs.app) 
[![cran version](https://www.r-pkg.org/badges/version/oai)](https://cran.r-project.org/package=oai) 
[![](https://badges.ropensci.org/19_status.svg)](https://github.com/ropensci/software-review/issues/19)

`oai` is an R client to work with OAI-PMH (Open Archives Initiative Protocol for Metadata Harvesting) services, a protocol developed by the Open Archives Initiative (https://en.wikipedia.org/wiki/Open_Archives_Initiative). OAI-PMH uses XML data format transported over HTTP.

OAI-PMH Info:

* Wikipedia (https://en.wikipedia.org/wiki/Open_Archives_Initiative_Protocol_for_Metadata_Harvesting)
* OAI V2 specification (http://www.openarchives.org/OAI/openarchivesprotocol.html)

`oai` is built on `xml2` and `httr`. In addition, we give back data.frame's whenever possible to make data comprehension, manipulation, and visualization easier. We also have functions to fetch a large directory of OAI-PMH services - it isn't exhaustive, but does contain a lot.

OAI-PMH instead of paging with e.g., `page` and `per_page` parameters, uses (optionally) `resumptionTokens`, optionally with an expiration date. These tokens can be used to continue on to the next chunk of data, if the first request did not get to the end. Often, OAI-PMH services limit each request to 50 records, but this may vary by provider, I don't know for sure. The API of this package is such that we `while` loop for you internally until we get all records. We may in the future expose e.g., a `limit` parameter so you can say how many records you want, but we haven't done this yet.

## Install

Install from CRAN

```{r eval=FALSE}
install.packages("oai")
```

Development version

```{r eval=FALSE}
devtools::install_github("ropensci/oai")
```

```{r}
library("oai")
```

## Identify

```{r}
id("http://oai.datacite.org/oai")
```

## ListIdentifiers

```{r}
list_identifiers(from = '2018-05-01T', until = '2018-06-01T')
```

## Count Identifiers

```{r cache=TRUE}
count_identifiers()
```

## ListRecords

```{r}
list_records(from = '2018-05-01T', until = '2018-05-15T')
```

## GetRecords

```{r}
ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", "d981c07d-bc43-40a2-be1f-e786e25106ac")
get_records(ids)
```

## List MetadataFormats

```{r}
list_metadataformats(id = "87832186-00ea-44dd-a6bf-c2896c4d09b4")
```

## List Sets

```{r}
list_sets("http://api.gbif.org/v1/oai-pmh/registry")
```

## Examples of other OAI providers

### Biodiversity Heritage Library

Identify

```{r}
id("http://www.biodiversitylibrary.org/oai")
```

Get records

```{r}
get_records(c("oai:biodiversitylibrary.org:item/7", "oai:biodiversitylibrary.org:item/9"),
            url = "http://www.biodiversitylibrary.org/oai")
```


## Acknowledgements

Michal Bojanowski thanks National Science Centre for support through grant 2012/07/D/HS6/01971.


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/oai/issues).
* License: MIT
* Get citation information for `oai` in R doing `citation(package = 'oai')`
* Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{oai introduction}
%\VignetteEncoding{UTF-8}
-->

```{r echo=FALSE}
knitr::opts_chunk$set(
	comment = "#>",
	collapse = TRUE,
	warning = FALSE,
	message = FALSE
)
```

oai introduction
================

A general purpose client to work with any 'OAI-PMH' service. 

The 'OAI-PMH' protocol is described at <http://www.openarchives.org/OAI/openarchivesprotocol.html>.

The main functions follow the OAI-PMH verbs:

* `GetRecord`
* `Identify`
* `ListIdentifiers`
* `ListMetadataFormats`
* `ListRecords`
* `ListSets`

## Get oai

Install from CRAN

```{r install, eval=FALSE}
install.packages("oai")
```

Or install the development version from GitHub

```{r installgh, eval=FALSE}
devtools::install_github("ropensci/oai")
```

Load `oai`

```{r load}
library("oai")
```

## Identify

```{r}
id("http://oai.datacite.org/oai")
```

## ListIdentifiers

```{r}
list_identifiers(from = '2018-05-01T', until = '2018-09-01T')
```

## Count Identifiers

```{r eval=FALSE}
count_identifiers()
```

## ListRecords

```{r}
list_records(from = '2018-05-01T', until = '2018-05-15T')
```

## GetRecords

```{r}
ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", "d981c07d-bc43-40a2-be1f-e786e25106ac")
get_records(ids)
```

## List MetadataFormats

```{r}
list_metadataformats(id = "87832186-00ea-44dd-a6bf-c2896c4d09b4")
```

## List Sets

```{r}
list_sets("http://api.gbif.org/v1/oai-pmh/registry")
```

### Biodiversity Heritage Library

Identify

```{r}
id("http://www.biodiversitylibrary.org/oai")
```

Get records

```{r}
get_records(c("oai:biodiversitylibrary.org:item/7", "oai:biodiversitylibrary.org:item/9"),
            url = "http://www.biodiversitylibrary.org/oai")
```
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{oai introduction}
%\VignetteEncoding{UTF-8}
-->



oai introduction
================

A general purpose client to work with any 'OAI-PMH' service. 

The 'OAI-PMH' protocol is described at <http://www.openarchives.org/OAI/openarchivesprotocol.html>.

The main functions follow the OAI-PMH verbs:

* `GetRecord`
* `Identify`
* `ListIdentifiers`
* `ListMetadataFormats`
* `ListRecords`
* `ListSets`

## Get oai

Install from CRAN


```r
install.packages("oai")
```

Or install the development version from GitHub


```r
devtools::install_github("ropensci/oai")
```

Load `oai`


```r
library("oai")
```

## Identify


```r
id("http://oai.datacite.org/oai")
#>   repositoryName                      baseURL protocolVersion
#> 1   DataCite MDS https://oai.datacite.org/oai             2.0
#>             adminEmail    earliestDatestamp deletedRecord
#> 1 support@datacite.org 2011-01-01T00:00:00Z    persistent
#>            granularity compression compression.1
#> 1 YYYY-MM-DDThh:mm:ssZ        gzip       deflate
#>                                      description
#> 1 oaioai.datacite.org:oai:oai.datacite.org:12425
```

## ListIdentifiers


```r
list_identifiers(from = '2018-05-01T', until = '2018-09-01T')
#> # A tibble: 2,904 x 5
#>    identifier        datestamp    setSpec           setSpec.1     setSpec.2
#>    <chr>             <chr>        <chr>             <chr>         <chr>    
#>  1 cbe4cdda-2045-41… 2018-08-30T… installation:842… dataset_type… country:…
#>  2 3fb7ddd8-07c0-49… 2018-08-28T… installation:842… dataset_type… country:…
#>  3 34585f24-1ffe-47… 2018-08-28T… installation:842… dataset_type… country:…
#>  4 8ac24ec8-1e7b-40… 2018-08-28T… installation:842… dataset_type… country:…
#>  5 e381b970-9b62-46… 2018-08-28T… installation:842… dataset_type… country:…
#>  6 4496b1b1-1ea9-4e… 2018-08-28T… installation:842… dataset_type… country:…
#>  7 fc5d68a6-b511-4c… 2018-08-27T… installation:73e… dataset_type… country:…
#>  8 5e0073bc-de2a-4b… 2018-08-27T… installation:688… dataset_type… country:…
#>  9 ab2684cf-62e5-4f… 2018-08-27T… installation:804… dataset_type… country:…
#> 10 34fbcf59-d9bb-47… 2018-08-27T… installation:7c5… dataset_type… country:…
#> # … with 2,894 more rows
```

## Count Identifiers


```r
count_identifiers()
```

## ListRecords


```r
list_records(from = '2018-05-01T', until = '2018-05-15T')
#> # A tibble: 44 x 26
#>    identifier datestamp setSpec setSpec.1 setSpec.2 title publisher
#>    <chr>      <chr>     <chr>   <chr>     <chr>     <chr> <chr>    
#>  1 18799ce9-… 2018-05-… instal… dataset_… country:… Bird… Sokoine …
#>  2 79f51633-… 2018-05-… instal… dataset_… country:… Impl… Aïgos SAS
#>  3 f83746ee-… 2018-05-… instal… dataset_… country:… NDFF… Dutch Na…
#>  4 a3533a61-… 2018-05-… instal… dataset_… country:… EDP … EDP - En…
#>  5 ba9b66a3-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  6 78b696d9-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  7 c791b255-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  8 b929ccda-… 2018-05-… instal… dataset_… country:… List… Sokoine …
#>  9 da285c2a-… 2018-05-… instal… dataset_… country:… Moni… Corporac…
#> 10 87372877-… 2018-05-… instal… dataset_… country:… Moni… Corporac…
#> # … with 34 more rows, and 19 more variables: identifier.1 <chr>,
#> #   subject <chr>, source <chr>, description <chr>, description.1 <chr>,
#> #   type <chr>, creator <chr>, date <chr>, language <chr>, coverage <chr>,
#> #   coverage.1 <chr>, format <chr>, source.1 <chr>, subject.1 <chr>,
#> #   coverage.2 <chr>, creator.1 <chr>, description.2 <chr>,
#> #   creator.2 <chr>, subject.2 <chr>
```

## GetRecords


```r
ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", "d981c07d-bc43-40a2-be1f-e786e25106ac")
get_records(ids)
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$header
#> # A tibble: 1 x 3
#>   identifier             datestamp      setSpec                            
#>   <chr>                  <chr>          <chr>                              
#> 1 87832186-00ea-44dd-a6… 2018-06-29T12… installation:729a7375-b120-4e4f-bb…
#> 
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$metadata
#> # A tibble: 0 x 0
#> 
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$header
#> # A tibble: 1 x 3
#>   identifier             datestamp      setSpec                            
#>   <chr>                  <chr>          <chr>                              
#> 1 d981c07d-bc43-40a2-be… 2018-01-21T21… installation:804b8dd0-07ac-4a30-bf…
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$metadata
#> # A tibble: 1 x 12
#>   title publisher identifier subject source description type  creator date 
#>   <chr> <chr>     <chr>      <chr>   <chr>  <chr>       <chr> <chr>   <chr>
#> 1 Pece… Institut… https://w… peces,… http:… Caracteriz… Data… Fernan… 2018…
#> # … with 3 more variables: language <chr>, coverage <chr>, format <chr>
```

## List MetadataFormats


```r
list_metadataformats(id = "87832186-00ea-44dd-a6bf-c2896c4d09b4")
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#>   metadataPrefix                                                   schema
#> 1         oai_dc           http://www.openarchives.org/OAI/2.0/oai_dc.xsd
#> 2            eml http://rs.gbif.org/schema/eml-gbif-profile/1.0.2/eml.xsd
#>                             metadataNamespace
#> 1 http://www.openarchives.org/OAI/2.0/oai_dc/
#> 2          eml://ecoinformatics.org/eml-2.1.1
```

## List Sets


```r
list_sets("http://api.gbif.org/v1/oai-pmh/registry")
#> # A tibble: 572 x 2
#>    setSpec                     setName         
#>    <chr>                       <chr>           
#>  1 dataset_type                per dataset type
#>  2 dataset_type:OCCURRENCE     occurrence      
#>  3 dataset_type:CHECKLIST      checklist       
#>  4 dataset_type:METADATA       metadata        
#>  5 dataset_type:SAMPLING_EVENT sampling_event  
#>  6 country                     per country     
#>  7 country:AD                  Andorra         
#>  8 country:AO                  Angola          
#>  9 country:AR                  Argentina       
#> 10 country:AT                  Austria         
#> # … with 562 more rows
```

### Biodiversity Heritage Library

Identify


```r
id("http://www.biodiversitylibrary.org/oai")
#>                                 repositoryName
#> 1 Biodiversity Heritage Library OAI Repository
#>                                   baseURL protocolVersion
#> 1 https://www.biodiversitylibrary.org/oai             2.0
#>                    adminEmail earliestDatestamp deletedRecord granularity
#> 1 oai@biodiversitylibrary.org        2006-01-01            no  YYYY-MM-DD
#>                                                        description
#> 1 oaibiodiversitylibrary.org:oai:biodiversitylibrary.org:item/1000
```

Get records


```r
get_records(c("oai:biodiversitylibrary.org:item/7", "oai:biodiversitylibrary.org:item/9"),
            url = "http://www.biodiversitylibrary.org/oai")
#> $`oai:biodiversitylibrary.org:item/7`
#> $`oai:biodiversitylibrary.org:item/7`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/7 2016-07-13T09:13:41Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/7`$metadata
#> # A tibble: 1 x 11
#>   title creator subject description publisher contributor date  type 
#>   <chr> <chr>   <chr>   <chr>       <chr>     <chr>       <chr> <chr>
#> 1 Die … Fleisc… Bogor;… pt.5:v.1 (… Leiden :… Missouri B… 1900… text…
#> # … with 3 more variables: identifier <chr>, language <chr>, rights <chr>
#> 
#> 
#> $`oai:biodiversitylibrary.org:item/9`
#> $`oai:biodiversitylibrary.org:item/9`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/9 2016-07-13T09:13:41Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/9`$metadata
#> # A tibble: 1 x 11
#>   title creator subject description publisher contributor date  type 
#>   <chr> <chr>   <chr>   <chr>       <chr>     <chr>       <chr> <chr>
#> 1 Die … Fleisc… Bogor;… pt.5:v.3 (… Leiden :… Missouri B… 1906… text…
#> # … with 3 more variables: identifier <chr>, language <chr>, rights <chr>
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_metadataformats.R
\name{list_metadataformats}
\alias{list_metadataformats}
\title{List available metadata formats from various providers.}
\usage{
list_metadataformats(
  url = "http://api.gbif.org/v1/oai-pmh/registry",
  id = NULL,
  ...
)
}
\arguments{
\item{url}{(character) OAI-PMH base url. Defaults to the URL for
arXiv's OAI-PMH server (http://export.arxiv.org/oai2)
or GBIF's OAI-PMH server (http://api.gbif.org/v1/oai-pmh/registry)}

\item{id}{The OAI-PMH identifier for the record. Optional.}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
List available metadata formats from various providers.
}
\examples{
\dontrun{
list_metadataformats()

# no metadatformats for an identifier
list_metadataformats(id = "9da8a65a-1b9b-487c-a564-d184a91a2705")

# metadatformats available for an identifier
list_metadataformats(id = "ad7295e0-3261-4028-8308-b2047d51d408")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_providers.R
\name{load_providers}
\alias{load_providers}
\title{Load an updated cache}
\usage{
load_providers(path = NULL, envir = .GlobalEnv)
}
\arguments{
\item{path}{location where cache is located. Leaving to NULL loads
the version in the installed package}

\item{envir}{R environment to load data in to.}
}
\value{
loads the object providers into the working space.
}
\description{
Load an updated cache
}
\details{
Loads the data object providers into the global workspace.
}
\examples{
\dontrun{
# By default the new providers table goes to directory ".", so just
# load from there
update_providers()
load_providers(path=".")

# Loads the version in the package
load_providers()
}
}
\seealso{
\code{\link[=update_providers]{update_providers()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_records.R
\name{list_records}
\alias{list_records}
\title{List records}
\usage{
list_records(
  url = "http://api.gbif.org/v1/oai-pmh/registry",
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
\item{url}{(character) OAI-PMH base url. Defaults to the URL for
arXiv's OAI-PMH server (http://export.arxiv.org/oai2)
or GBIF's OAI-PMH server (http://api.gbif.org/v1/oai-pmh/registry)}

\item{prefix}{specifies the metadata format that the records will be
returned in. Default: \code{oai_dc}}

\item{from}{specifies that records returned must have been
created/update/deleted on or after this date.}

\item{until}{specifies that records returned must have been
created/update/deleted on or before this date.}

\item{set}{specifies the set that returned records must belong to.}

\item{token}{(character) a token previously provided by the server to
resume a request where it last left off. 50 is max number of records
returned. We will loop for you internally to get all the records you
asked for.}

\item{as}{(character) What to return. One of "df" (for data.frame; default),
"list", or "raw" (raw text)}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
List records
}
\examples{
\dontrun{
# By default you get back a single data.frame
list_records(from = '2018-05-01T00:00:00Z', until = '2018-05-03T00:00:00Z')
list_records(from = '2018-05-01T', until = '2018-05-04T')

# Get a list
list_records(from = '2018-05-01T', until = '2018-05-04T', as = "list")

# Get raw text
list_records(from = '2018-05-01T', until = '2018-05-04T', as = "raw")
list_records(from = '2018-05-01T', until = '2018-05-04T', as = "raw")

# Use a resumption token
# list_records(token =
#  "1443799900201,2015-09-01T00:00:00Z,2015-10-01T23:59:59Z,50,null,oai_dc")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_record.R
\name{get_records}
\alias{get_records}
\title{Get records}
\usage{
get_records(
  ids,
  prefix = "oai_dc",
  url = "http://api.gbif.org/v1/oai-pmh/registry",
  as = "parsed",
  ...
)
}
\arguments{
\item{ids}{The OAI-PMH identifier for the record. One or more. Required.}

\item{prefix}{specifies the metadata format that the records will be
returned in. Default: \code{oai_dc}}

\item{url}{(character) OAI-PMH base url. Defaults to the URL for
arXiv's OAI-PMH server (http://export.arxiv.org/oai2)
or GBIF's OAI-PMH server (http://api.gbif.org/v1/oai-pmh/registry)}

\item{as}{(character) What to return. One of "parsed" (default),
or "raw" (raw text)}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\value{
a named list of data.frame's, or lists, or raw text
}
\description{
Get records
}
\details{
There are some finite set of results based on the OAI prefix.
We will provide parsers as we have time, and as users express interest.
For prefix types we have parsers for we return a list of data.frame's,
for each identifier, one data.frame for the \code{header} bits of data, and
one data.frame for the \code{metadata} bits of data.

For prefixes we don't have parsers for, we fall back to returning raw
XML, so you can at least parse the XML yourself.

Because some XML nodes are duplicated, we join values together of
duplicated node names, separated by a semicolon (\verb{;}) with no
spaces. You can seprarate them yourself easily.
}
\examples{
\dontrun{
get_records("87832186-00ea-44dd-a6bf-c2896c4d09b4")

ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", 
  "d981c07d-bc43-40a2-be1f-e786e25106ac")
(res <- get_records(ids))
lapply(res, "[[", "header")
lapply(res, "[[", "metadata")
do.call(rbind, lapply(res, "[[", "header"))
do.call(rbind, lapply(res, "[[", "metadata"))

# Get raw text
get_records("d981c07d-bc43-40a2-be1f-e786e25106ac", as = "raw")

# from arxiv.org
get_records("oai:arXiv.org:0704.0001", url = "http://export.arxiv.org/oai2")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oai-package.R
\docType{data}
\name{providers}
\alias{providers}
\title{Metadata providers data.frame.}
\value{
A data.frame of three columns:
\itemize{
\item repo_name - Name of the OAI repository
\item base_url - Base URL of the OAI repository
\item oai_identifier - OAI identifier for the OAI repository
}
}
\description{
Metadata providers data.frame.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oai-package.R
\docType{package}
\name{oai-package}
\alias{oai-package}
\alias{oai}
\title{OAI-PMH Client}
\description{
oai is an R client to work with OAI-PMH (Open Archives
Initiative Protocol for Metadata Harvesting) services, a protocol
developed by the Open Archives Initiative
(https://en.wikipedia.org/wiki/Open_Archives_Initiative).
OAI-PMH uses XML data format transported over HTTP.
}
\section{OAI-PMH Info}{

See the OAI-PMH V2 specification at
\url{http://www.openarchives.org/OAI/openarchivesprotocol.html}
}

\section{Implementation details}{

oai is built on \pkg{xml2} and \pkg{httr}. In addition, we give back
data.frame's whenever possible to make data comprehension, manipulation,
and visualization easier. We also have functions to fetch a large directory
of OAI-PMH services - it isn't exhaustive, but does contain a lot.
}

\section{Paging}{

Instead of paging with e.g., \code{page} and \code{per_page} parameters,
OAI-PMH uses (optionally) \code{resumptionTokens}, with an optional
expiration date. These tokens can be used to continue on to the next chunk
of data, if the first request did not get to the end. Often, OAI-PMH
services limit each request to 50 records, but this may vary by provider,
I don't know for sure. The API of this package is such that we \code{while}
loop for you internally until we get all records. We may in the future
expose e.g., a \code{limit} parameter so you can say how many records
you want, but we haven't done this yet.
}

\section{Acknowledgements}{

Michal Bojanowski contributions were supported by (Polish) National Science
Center (NCN) through grant 2012/07/D/HS6/01971.
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}

Michal Bojanowski \email{michal2992@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dumpers.R
\name{dumpers}
\alias{dumpers}
\alias{dump_raw_to_txt}
\alias{dump_to_rds}
\alias{dump_raw_to_db}
\title{Result dumpers}
\usage{
dump_raw_to_txt(
  res,
  args,
  as,
  file_pattern = "oaidump",
  file_dir = ".",
  file_ext = ".xml"
)

dump_to_rds(
  res,
  args,
  as,
  file_pattern = "oaidump",
  file_dir = ".",
  file_ext = ".rds"
)

dump_raw_to_db(res, args, as, dbcon, table_name, field_name, ...)
}
\arguments{
\item{res}{results, depends on \code{as}, not to be specified by the user}

\item{args}{list, query arguments, not to be specified by the user}

\item{as}{character, type of result to return, not to be specified by the
user}

\item{file_pattern, file_dir, file_ext}{character respectively: initial part of
the file name, directory name, and file extension used to create file
names. These arguments are passed to \code{\link[=tempfile]{tempfile()}} arguments
\code{pattern}, \code{tmpdir}, and \code{fileext} respectively.}

\item{dbcon}{\pkg{DBI}-compliant database connection}

\item{table_name}{character, name of the database table to write into}

\item{field_name}{character, name of the field in database table to write
into}

\item{...}{arguments passed to/from other functions}
}
\value{
Dumpers should return \code{NULL} or a value that will be collected
and returned by the function using the dumper.

\code{dump_raw_to_txt} returns the name of the created file.

\code{dump_to_rds} returns the name of the created file.

\code{dump_xml_to_db} returns \code{NULL}
}
\description{
Result dumpers are functions allowing to handle the chunks of results from
OAI-PMH service "on the fly". Handling can include processing, writing to
files, databases etc.
}
\details{
Often the result of a request to a OAI-PMH service are so large that it is
split into chunks that need to be requested separately using
\code{resumptionToken}. By default functions like
\code{\link[=list_identifiers]{list_identifiers()}} or \code{\link[=list_records]{list_records()}} request these
chunks under the hood and return all concatenated in a single R object. It
is convenient but insufficient when dealing with large result sets that
might not fit into RAM. A result dumper is a function that is called on
each result chunk. Dumper functions can write chunks to files or databases,
include initial pre-processing or extraction, and so on.

A result dumper needs to be function that accepts at least the arguments:
\code{res}, \code{args}, \code{as}. They will get values by the enclosing
function internally. There may be additional arguments, including \code{...}.
Dumpers should return \code{NULL} or a value that will
be collected and returned by the function calling the dumper (e.g.
\code{\link[=list_records]{list_records()}}).

Currently result dumpers can be used with functions:
\code{\link[=list_identifiers]{list_identifiers()}}, \code{\link[=list_records]{list_records()}}, and \code{\link[=list_sets]{list_sets()}}.
To use a dumper with one of these functions you need to:
\itemize{
\item Pass it as an additional argument \code{dumper}
\item Pass optional addtional arguments to the dumper function in a list
as the \code{dumper_args} argument
}

See Examples. Below we provide more details on the dumpers currently
implemented.

\code{dump_raw_to_txt} writes raw XML to text files. It requires
\code{as=="raw"}. File names are created using \code{\link[=tempfile]{tempfile()}}. By
default they are written in the current working directory and have a format
\code{oaidump*.xml} where \code{*} is a random string in hex.

\code{dump_to_rds} saves results in an \code{.rds} file via \code{\link[=saveRDS]{saveRDS()}}.
Type of object being saved is determined by the \code{as} argument. File names
are generated in the same way as by \code{dump_raw_to_txt}, but with default
extension \code{.rds}

\code{dump_xml_to_db} writes raw XML to a single text column of a table in a
database. Requires \code{as == "raw"}. Database connection \code{dbcon}
should be a connection object as created by \code{\link[DBI:dbConnect]{DBI::dbConnect()}} from
package \pkg{DBI}. As such, it can connect to any database supported by
\pkg{DBI}. The records are written to a field \code{field_name} in a table
\code{table_name} using \code{\link[DBI:dbWriteTable]{DBI::dbWriteTable()}}. If the table does not
exist, it is created. If it does, the records are appended. Any additional
arguments are passed to \code{\link[DBI:dbWriteTable]{DBI::dbWriteTable()}}
}
\examples{
\dontrun{

### Dumping raw XML to text files

# This will write a set of XML files to a temporary directory
fnames <- list_identifiers(from="2018-06-01T",
                           until="2018-06-14T",
                           as="raw",
                           dumper=dump_raw_to_txt,
                           dumper_args=list(file_dir=tempdir()))
# vector of file names created
str(fnames)
all( file.exists(fnames) )
# clean-up
unlink(fnames)


### Dumping raw XML to a database

# Connect to in-memory SQLite database
con <- DBI::dbConnect(RSQLite::SQLite(), dbname=":memory:")
# Harvest and dump the results into field "bar" of table "foo"
list_identifiers(from="2018-06-01T",
                 until="2018-06-14T",
                 as="raw",
                 dumper=dump_raw_to_db,
                 dumper_args=list(dbcon=con,
                                  table_name="foo",
                                  field_name="bar") )
# Count records, should be 101
DBI::dbGetQuery(con, "SELECT count(*) as no_records FROM foo")

DBI::dbDisconnect(con)




}
}
\references{
OAI-PMH specification
\url{https://www.openarchives.org/OAI/openarchivesprotocol.html}
}
\seealso{
Functions supporting the dumpers:
\code{\link[=list_identifiers]{list_identifiers()}}, \code{\link[=list_sets]{list_sets()}}, and \code{\link[=list_records]{list_records()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_identifiers.R
\name{count_identifiers}
\alias{count_identifiers}
\title{Count OAI-PMH identifiers for a data provider.}
\usage{
count_identifiers(url = "http://export.arxiv.org/oai2", prefix = "oai_dc", ...)
}
\arguments{
\item{url}{(character) OAI-PMH base url. Defaults to the URL for
arXiv's OAI-PMH server (http://export.arxiv.org/oai2)
or GBIF's OAI-PMH server (http://api.gbif.org/v1/oai-pmh/registry)}

\item{prefix}{Specifies the metadata format that the records will be
returned in}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Count OAI-PMH identifiers for a data provider.
}
\details{
Note that some OAI providers do not include the entry
\code{completeListSize}
(\url{http://www.openarchives.org/OAI/openarchivesprotocol.html#FlowControl})
in which case we return an NA - which does not mean 0, but rather we don't
know.
}
\examples{
\dontrun{
count_identifiers()

# curl options
# library("httr")
# count_identifiers(config = verbose())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oai_available.R
\name{oai_available}
\alias{oai_available}
\title{Test of OAI-PMH service is available}
\usage{
oai_available(u, ...)
}
\arguments{
\item{u}{base URL to OAI-PMH service}

\item{...}{other arguments passed to \code{\link[=id]{id()}}}
}
\value{
\code{TRUE} or \code{FALSE} if the service is available.
}
\description{
Silently test if OAI-PMH service is available under the URL provided.
}
\examples{
\dontrun{
url_list <- list(
  archivesic="http://archivesic.ccsd.cnrs.fr/oai/oai.php",
  datacite = "http://oai.datacite.org/oai",

  # No OAI-PMH here
  google = "http://google.com"
)

sapply(url_list, oai_available)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_identifiers.R
\name{list_identifiers}
\alias{list_identifiers}
\title{List OAI-PMH identifiers}
\usage{
list_identifiers(
  url = "http://api.gbif.org/v1/oai-pmh/registry",
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
\item{url}{(character) OAI-PMH base url. Defaults to the URL for
arXiv's OAI-PMH server (http://export.arxiv.org/oai2)
or GBIF's OAI-PMH server (http://api.gbif.org/v1/oai-pmh/registry)}

\item{prefix}{Specifies the metadata format that the records will be
returned in.}

\item{from}{specifies that records returned must have been
created/update/deleted on or after this date.}

\item{until}{specifies that records returned must have been
created/update/deleted on or before this date.}

\item{set}{specifies the set that returned records must belong to.}

\item{token}{a token previously provided by the server to resume a request
where it last left off.}

\item{as}{(character) What to return. One of "df" (for data.frame; default),
"list", or "raw" (raw text)}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
List OAI-PMH identifiers
}
\examples{
\dontrun{
# from
recently <- format(Sys.Date() - 1, "\%Y-\%m-\%d")
list_identifiers(from = recently)

# from and until
list_identifiers(from = '2018-06-01T', until = '2018-06-14T')

# set parameter - here, using ANDS - Australian National Data Service
list_identifiers(from = '2018-09-01T', until = '2018-09-05T',
  set = "dataset_type:CHECKLIST")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_providers.R
\name{update_providers}
\alias{update_providers}
\title{Update the locally stored OAI-PMH data providers table.}
\usage{
update_providers(path = ".", ...)
}
\arguments{
\item{path}{Path to put data in.}

\item{...}{Curl options passed on to \code{\link[httr:GET]{httr::GET()}}}
}
\description{
Data comes from
\url{http://www.openarchives.org/Register/BrowseSites}. It includes the
oai-identifier (if they have one) and the base URL. The website has
the name of the data provider too, but not provided in the data pulled
down here, but you can grab the name using the example below.
}
\details{
This table is scraped from
\url{http://www.openarchives.org/Register/BrowseSites}.
I would get it from \url{http://www.openarchives.org/pmh/registry/ListFriends},
but it does not include repository names.

This function updates the table for you. Does take a while though, so
go get a coffee.
}
\examples{
\dontrun{
update_providers()
load_providers()
}
}
\seealso{
\code{\link[=load_providers]{load_providers()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/id.R
\name{id}
\alias{id}
\title{Identify the OAI-PMH service for each data provider.}
\usage{
id(url, as = "parsed", ...)
}
\arguments{
\item{url}{(character) OAI-PMH base url. Defaults to the URL for
arXiv's OAI-PMH server (http://export.arxiv.org/oai2)
or GBIF's OAI-PMH server (http://api.gbif.org/v1/oai-pmh/registry)}

\item{as}{(character) What to return. One of "parsed" (default),
or "raw" (raw text)}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
Identify the OAI-PMH service for each data provider.
}
\examples{
\dontrun{
# arxiv
id("http://export.arxiv.org/oai2")

# GBIF - http://www.gbif.org/
id("http://api.gbif.org/v1/oai-pmh/registry")

# get back text instead of parsed
id("http://export.arxiv.org/oai2", as = "raw")
id("http://api.gbif.org/v1/oai-pmh/registry", as = "raw")

# curl options
library("httr")
id("http://export.arxiv.org/oai2", config = verbose())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_sets.R
\name{list_sets}
\alias{list_sets}
\title{List sets}
\usage{
list_sets(
  url = "http://api.gbif.org/v1/oai-pmh/registry",
  token = NULL,
  as = "df",
  ...
)
}
\arguments{
\item{url}{(character) OAI-PMH base url. Defaults to the URL for
arXiv's OAI-PMH server (http://export.arxiv.org/oai2)
or GBIF's OAI-PMH server (http://api.gbif.org/v1/oai-pmh/registry)}

\item{token}{(character) a token previously provided by the server to
resume a request where it last left off}

\item{as}{(character) What to return. One of "df" (for data.frame; default),
"list", or "raw" (raw text)}

\item{...}{Curl options passed on to \code{\link[httr]{GET}}}
}
\description{
List sets
}
\examples{
\dontrun{
# Get back a data.frame
list_sets()

# Get back a list
list_sets(as = "list")

# Get back raw text
list_sets(as = "raw")

# curl options
library("httr")
list_sets(config = verbose())
}
}
