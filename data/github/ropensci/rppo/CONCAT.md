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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# rppo

[![Build
Status](https://travis-ci.org/ropensci/rppo.svg?branch=master)](https://travis-ci.org/ropensci/rppo)
[![codecov.io](https://codecov.io/github/r-lib/covr/coverage.svg?branch=master)](https://codecov.io/github/r-lib/covr?branch=master)
[![](https://badges.ropensci.org/207_status.svg)](https://github.com/ropensci/software-review/issues/207)

The global plant phenology data portal, is an aggregation of plant
phenological observations from
[USA-NPN](https://www.usanpn.org/usa-national-phenology-network),
[NEON](https://www.neonscience.org/), and
[PEP725](https://www.pep725.eu/) representing 20 million phenological
observations from across North America and Europe. The PPO data portal
utilizes the [Plant Phenology
Ontology](https://github.com/PlantPhenoOntology/ppo/) (PPO) to align
phenological terms and measurements from the various databases. The rppo
R package enables programmatic access to all data contained in the PPO
data portal incuding selected classes contained in the PPO itself.

For information on how data is assembled for the PPO data portal, visit
the [ppo-data-pipeline git
repository](https://github.com/biocodellc/ppo-data-pipeline).

## Installation

The production version of rppo is accessible on CRAN:

``` r
install.packages("rppo")  
#> installing the source package 'rppo'
library(rppo)
```

You can install the development version of rppo from github with:

``` r
install.packages("devtools")
devtools::install_github("ropensci/rppo")
library(rppo)
```

## Examples

Following are a couple of brief examples to illustrate how to get
started with rppo. We recommend visiting the [rppo
vignette](https://htmlpreview.github.io/?https://github.com/ropensci/rppo/blob/master/vignettes/rppo-vignette.html)
for a more complete set of examples on using the rppo package, as well
as viewing man pages for rppo functions in the R environment, using
`?ppo_data` and
`?ppo_terms`.

``` r
# query all results from day 1 through 100 in a particular bounding box, 
# limited to 2 records
r <- ppo_data(fromDay = 1, toDay = 100, bbox="37,-120,38,-119", limit=2, timeLimit=4)
#> sending request for data ...
#> https://www.plantphenology.org/api/v2/download/?q=%2Blatitude:>=37+AND+%2Blatitude:<=38+AND+%2Blongitude:>=-120+AND+%2Blongitude:<=-119+AND+%2BdayOfYear:>=1+AND+%2BdayOfYear:<=100+AND+source:USA-NPN,NEON&source=latitude,longitude,year,dayOfYear,termID&limit=2

# view the data returned
print(r$data)
#>   dayOfYear year   genus specificEpithet latitude longitude
#> 1        33 2017 Quercus       douglasii 37.11144 -119.7315
#> 2        96 2017  Bromus        diandrus 37.11144 -119.7315
#>                                                                                                                                                                                                                                            termID
#> 1                                                                                                                                                                                                 obo:PPO_0002610,obo:PPO_0002013,obo:PPO_0002000
#> 2 obo:PPO_0002601,obo:PPO_0002610,obo:PPO_0002005,obo:PPO_0002604,obo:PPO_0002605,obo:PPO_0002013,obo:PPO_0002003,obo:PPO_0002000,obo:PPO_0002602,obo:PPO_0002006,obo:PPO_0002007,obo:PPO_0002004,obo:PPO_0002008,obo:PPO_0002603,obo:PPO_0002600
#>   source
#> 1   NEON
#> 2   NEON
#>                                                              eventId
#> 1 https://n2t.net/ark:/21547/Amn2cd982ca2-6147-4a63-a864-f4e556420562
#> 2 https://n2t.net/ark:/21547/Amn2d1a3e6de-7885-404f-828f-9ebf63248d68

# view the number of possible records returned
print(r$number_possible)
#> [1] 7251

# return a data frame of present
presentTerms <- ppo_terms(present = TRUE, timeLimit=3)
#> sending request for terms ...

# print the 2nd present term returned
print(presentTerms[2,])
#>            termID                            label
#> 2 obo:PPO_0002358 abscised fruits or seeds present
#>                                                                                                                                                                                                                                                                                                                                  definition
#> 2 An 'abscised fruit or seed presence' (PPO:0002059) trait that is a 'quality of' (RO:0000080) a 'whole plant' (PO:0000003) from which at least one 'ripe fruit' (PPO:0001045) has been abscised or removed by an herbivore or that has at least one 'ripe fruit' (PPO:0001045) that has abscised at least one 'mature seed' (PPO:0001024).
#>                                          uri
#> 2 https://purl.obolibrary.org/obo/PPO_0002358
```

## Citation

To cite the ‘rppo’ R package in publications
use:

``` 
   'John Deck, Brian Stucky, Ramona Walls, Kjell Bolmgren, Ellen Denny, Robert Guralnick' (2018). rppo: An interface to the Plant Phenology Ontology and associated data store.  R package version 1.0
   https://github.com/ropensci/rppo
```

## Code of Conduct

View our [code of conduct](https://github.com/ropensci/rppo/blob/master/CONDUCT.md)

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# rppo 1.0.2
## New Features
 * updated release to gracefully check on server status

## Resubmission
This is a resubmission. In this version I have:

* updated code to check if a server is responding or not.  Previous versions
* checked for non-200 status messages.  This release checks if server is even 
* responding at all.  Currently the server responding to all API calls has
* completely failed and will be rebuilt by June 7th.  Meanwhile, we are re-releasing
* the R package here so it fails gracefully.
* Please NOTE that the WIN package provides a note about the spelling of the word "phenology", which is actually correct.

## Test environments
* local OS X install, R 3.5.0
* ubuntu 16.04 (on travis-ci), R 4.0.2

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies
There are currently no downstream dependencies for this package

---
title: "rppo vignette"
author: "John Deck"
date: "2018-05-23"
output:
 html_document:
    keep_md: yes
vignette: |
  %\VignetteIndexEntry{rppo Vignette} 
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---




The rppo package contains just two functions.  One to query terms from the Plant Phenology Ontology (PPO) and another to query the data global plant phenology data portal (PPO data portal).  Following are three examples which illustrate use of these functions: the first two sections illustrate the `ppo_data` and `ppo_terms` functions and the third section illustrates how to use the functions together.

### ppo_terms function
It is frequently useful to look through the list of present and absent terms contained in the PPO.   The `ppo_terms` function returns present terms, absent terms, or both, with columns containing a termID, label, definition and full URI for each term.  Use the termIDs returned from this function to query terms in the `ppo_data` function.  The following example returns the present terms into a "present_terms" data frame and a sample slice from the dataframe.


```r
present_terms <- ppo_terms(present = TRUE, timeLimit = 3)
# print the first five rows, with just the termIDs and labels
print(present_terms[1:5,c("termID","label")])
#>            termID                                              label
#> 1 obo:PPO_0002300                           plant structures present
#> 2 obo:PPO_0002301                           new shoot system present
#> 3 obo:PPO_0002302 new above-ground shoot-borne shoot systems present
#> 4 obo:PPO_0002311                         breaking leaf buds present
#> 5 obo:PPO_0002303     new shoot systems emerging from ground present
```

### ppo_data function
The `ppo_data` function queries the PPO Data Portal, passing values to the database and extracting matching results. The results of the `ppo_data` function are returned as a list with five elements: 1) a data frame containing data, 2) a readme string containing usage information and some statistics about the query itself, 3) a citation string containing information about proper citation, 4) a number_possible integer indicating the total number of results if a limit has been specified, and 5) a status code returned from the service. 

The "df" variable below is populated with results from the data element in the results list, with an example slice of data showing the first record.

```r
results <- ppo_data(genus = "Quercus", fromYear = 2013, toYear = 2013, fromDay = 100, toDay = 110, termID = 'obo:PPO_0002313', limit = 10, timeLimit = 2)
df <- results$data
print(df[1:1,])
#>   dayOfYear year   genus specificEpithet latitude longitude
#> 1       106 2013 Quercus          lobata 34.67545 -120.0407
#>                                                                                                                                                                                                                                            termID
#> 1 obo:PPO_0002316,obo:PPO_0002322,obo:PPO_0002018,obo:PPO_0002318,obo:PPO_0002022,obo:PPO_0002312,obo:PPO_0002313,obo:PPO_0002014,obo:PPO_0002024,obo:PPO_0002015,obo:PPO_0002000,obo:PPO_0002017,obo:PPO_0002020,obo:PPO_0002320,obo:PPO_0002315
#>    source                               eventId
#> 1 USA-NPN http://n2t.net/ark:/21547/Amg22054478
```

The readme and citation files returned by the list of results can be accessed by calling the readme and citation elements.  Note that the the file "citation_and_data_use_policies.txt" that is referred to in the readme file can be accessed using `cat(results$citation)` 

```r
cat(results$readme)
#> The following contains information about your download from the Global Plant 
#> Phenology Database.  Please refer to the citation_and_data_use_policies.txt 
#> file for important information about data usage policies, licensing, and 
#> citation protocols for each dataset.  This file contains summary information 
#> about the query that was run.  
#> 
#> data file = data.csv
#> date query ran = Wed May 23 2018 20:22:34 GMT-0400 (EDT)
#> query = +genus:Quercus AND +plantStructurePresenceTypes:"http://purl.obolibrary.org/obo/PPO_0002313" AND +year:>=2013 AND +year:<=2013 AND +dayOfYear:>=100 AND +dayOfYear:<=110 AND source:USA-NPN,NEON
#> fields returned = dayOfYear,year,genus,specificEpithet,latitude,longitude,source,eventId
#> user specified limit = 10
#> total results possible = 518
#> total results returned = 0
```

The results lists also shows the number of possible results in the results set, which is useful if the submitted query had a limit.  For example, in the query above, the limit is set to 10 but we want to know how many records were possible if the limit was not set.

```r
cat(results$number_possible)
#> 518
```

### working with terms and data together
Here we will generate a data frame showing the frequency of "present" and "absent" terms for a particular query.  The query is for genus = "Quercus" and latitude > 47.  For each row in the returned data frame `ppo_data` will typically return multiple terms in the termID field, corresponding to phenological stages as defined by the PPO.  For our example, we will generate a frequency table of the number of times "present" or "absent" term occur in the entire returned dataset.  Note that the termID field returned by `ppo_data` will return "presence" terms in addition to "present" and "absent" terms, while the `ppo_terms` function only returns "present" and "absent" terms.  Thus, our frequency distribution only counts the number of "present" and "absent" terms [For an in-depth discussion of the difference between "presence", "present", and "absent", see https://www.frontiersin.org/articles/10.3389/fpls.2018.00517/full].  Finally, since termIDs are returned as URI identifiers and not easily readable text, this example maps termIDs to labels. The resulting data frame shows two columns: 1) a column of term labels, and 2) a frequency of the number of times this label appeared in the result set. 


```r
###############################################################################
# Generate a frequency data frame showing the number of times each termID
# is populated for genus equals "Quercus" above latitude of 47
# Note that all latitude/longitude queries need to be in the format of a
# bounding box
###############################################################################
df <- ppo_data(
  genus = "Quercus", 
  bbox="47,-180,90,180", timeLimit = 2)
#> sending request for data ...
#> https://www.plantphenology.org/api/v2/download/?q=%2Bgenus:Quercus+AND+%2Blatitude:>=47+AND+%2Blatitude:<=90+AND+%2Blongitude:>=-180+AND+%2Blongitude:<=180+AND+source:USA-NPN,NEON&source=latitude,longitude,year,dayOfYear,termID
# return just the termID column
t1 <- df$data[,c('termID')]
# paste each cell into one string
t2<-paste(t1, collapse = ",")
# split strings at ,
t3<-strsplit(t2, ",")
# create a frequency table as a data frame
freqFrame <- as.data.frame(table(t3))

# create a new data frame that we want to populate
resultFrame <- data.frame(
  label = character(), 
  frequency = integer(), 
  stringsAsFactors = FALSE)

###############################################################################
# Replace termIDs with labels in frequency frame
###############################################################################
# fetch "present" and "absent" terms using `ppo_terms`
termList <- ppo_terms(absent = TRUE, present = TRUE, timeLimit =2);
#> sending request for terms ...

# loop all "present"" and "absent" terms
for (term in 1:nrow(termList)) {
  termListTermID<-termList[term,'termID'];
  termListLabel<-termList[term,'label'];
  # loop all rows that have a frequency generated
  for (row in 1:nrow(freqFrame)) {
    freqFrameTermID = freqFrame[row,'t3']
    freqFrameFrequency = freqFrame[row,'Freq']
    # Populate resultFrame with matching "present" or "absent" labels.
    # In this step, we will ignore "presence" terms
    # found in the frequency frame since the ppo_terms only returns
    # "present" and "absent" terms. 
    if (freqFrameTermID == termListTermID) {
      resultFrame[nrow(resultFrame)+1,] <- c(termListLabel,freqFrameFrequency)
    }
  }
}

# print results, showing term labels and a frequency count
print(resultFrame)
#>                                                 label frequency
#> 1                            new shoot system present        32
#> 2  new above-ground shoot-borne shoot systems present        32
#> 3                          breaking leaf buds present        32
#> 4                                   leaf buds present        32
#> 5                       non-dormant leaf buds present        32
#> 6                             vascular leaves present       101
#> 7                                 true leaves present       101
#> 8                        unfolded true leaves present        69
#> 9          non-senescing unfolded true leaves present        22
#> 10              immature unfolded true leaves present        22
#> 11             expanding unfolded true leaves present        22
#> 12                      senescing true leaves present         6
#> 13                      expanding true leaves present        54
#> 14                      unfolding true leaves present        32
#> 15                    reproductive structures present        28
#> 16                          floral structures present        16
#> 17             non-senesced floral structures present         9
#> 18                                    flowers present         7
#> 19                       non-senesced flowers present         7
#> 20                               open flowers present         7
#> 21                   pollen-releasing flowers present         4
#> 22                                     fruits present        12
#> 23                            ripening fruits present        12
#> 24                            abscised leaves present         4
#> 25                          breaking leaf buds absent       159
#> 26                       senescing true leaves absent       342
#> 27                        unfolded true leaves absent       159
#> 28                          mature true leaves absent       159
#> 29          non-senescing unfolded true leaves absent       159
#> 30              expanding unfolded true leaves absent       159
#> 31               immature unfolded true leaves absent       159
#> 32               expanded immature true leaves absent       159
#> 33                  unopened floral structures absent       175
#> 34              non-senesced floral structures absent       175
#> 35          pollen-releasing floral structures absent       533
#> 36                      open floral structures absent       175
#> 37                    pollen-releasing flowers absent       358
#> 38                                open flowers absent       181
#> 39               pollen-releasing flower heads absent       358
#> 40                           open flower heads absent       181
#> 41                               unripe fruits absent       176
#> 42                             ripening fruits absent       176
#> 43                                 ripe fruits absent       362
#> 44                             abscised leaves absent       365
#> 45                    abscised fruits or seeds absent       365
#> 46                     abscised cones or seeds absent       365
```
---
output: github_document
always_allow_html: yes
editor_options: 
  chunk_output_type: console
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# rppo
[![Build Status](https://travis-ci.org/ropensci/rppo.svg?branch=master)](https://travis-ci.org/ropensci/rppo)
[![codecov.io](https://codecov.io/github/r-lib/covr/coverage.svg?branch=master)](https://codecov.io/github/r-lib/covr?branch=master)
[![](https://badges.ropensci.org/207_status.svg)](https://github.com/ropensci/onboarding/issues/207)



The global plant phenology data portal, or PPO data portal, is an aggregation of 
plant phenological observations from [USA-NPN](https://www.usanpn.org/usa-national-phenology-network), [NEON](https://www.neonscience.org/), and [PEP725](http://www.pep725.eu/) representing 
20 million phenological observations from across North America and Europe.  The PPO data portal utilizes the [Plant
Phenology Ontology](https://github.com/PlantPhenoOntology/ppo/) (PPO) to align phenological
terms and measurements from the various databases. The rppo R package enables programmatic access to all data contained in the PPO data portal incuding selected classes contained in the PPO itself.  

For information on how data is assembled for the PPO data portal, visit the 
[ppo-data-pipeline git repository](https://github.com/biocodellc/ppo-data-pipeline).

## Installation

The production version of rppo is accessible on CRAN:
```{r}
install.packages("rppo")  
library(rppo)
```

You can install the development version of rppo from github with:

```{r gh-installation, eval = FALSE}
install.packages("devtools")
devtools::install_github("ropensci/rppo")
library(rppo)
```

## Examples

Following are a couple of brief examples to illustrate how to get started with rppo.  We recommend visiting the [rppo vignette](http://htmlpreview.github.io/?https:/x/github.com/ropensci/rppo/blob/master/vignettes/rppo-vignette.html) for a more complete set of examples on using the rppo package, as well as viewing man pages for rppo functions in the R environment, using `?ppo_data` and `?ppo_terms`.  

```{r examples}
# query all results from day 1 through 100 in a particular bounding box, 
# limited to 2 records
r <- ppo_data(fromDay = 1, toDay = 100, bbox="37,-120,38,-119", limit=2, timeLimit=5)

# view the data returned
print(r$data)

# view the number of possible records returned
print(r$number_possible)

# return a data frame of present
presentTerms <- ppo_terms(present = TRUE, timeLimit=3)

# print the 2nd present term returned
print(presentTerms[2,])
```

## Citation
To cite the 'rppo' R package in publications use:

```
   'John Deck, Brian Stucky, Ramona Walls, Kjell Bolmgren, Ellen Denny, Robert Guralnick' (2018). rppo: An interface to the Plant Phenology Ontology and associated data store.  R package version 1.0
   https://github.com/ropensci/rppo
```

## Code of Conduct
View our [code of conduct](CONDUCT.md)

[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "rppo vignette"
author: "John Deck"
date: "`r Sys.Date()`"
output:
 html_document:
    keep_md: yes
vignette: |
  %\VignetteIndexEntry{rppo Vignette} 
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---


```{r setup, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "rppo-vignette-"
)
library(rppo)
```

The rppo package contains just two functions.  One to query terms from the Plant Phenology Ontology (PPO) and another to query the data global plant phenology data portal (PPO data portal).  Following are three examples which illustrate use of these functions: the first two sections illustrate the `ppo_data` and `ppo_terms` functions and the third section illustrates how to use the functions together.

### ppo_terms function
It is frequently useful to look through the list of present and absent terms contained in the PPO.   The `ppo_terms` function returns present terms, absent terms, or both, with columns containing a termID, label, definition and full URI for each term.  Use the termIDs returned from this function to query terms in the `ppo_data` function.  The following example returns the present terms into a "present_terms" data frame and a sample slice from the dataframe.

```{r term example, echo=TRUE, message=FALSE, warning=FALSE}
present_terms <- ppo_terms(present = TRUE)
# print the first five rows, with just the termIDs and labels
print(present_terms[1:5,c("termID","label")])
```

### ppo_data function
The `ppo_data` function queries the PPO Data Portal, passing values to the database and extracting matching results. The results of the `ppo_data` function are returned as a list with five elements: 1) a data frame containing data, 2) a readme string containing usage information and some statistics about the query itself, 3) a citation string containing information about proper citation, 4) a number_possible integer indicating the total number of results if a limit has been specified, and 5) a status code returned from the service. 

The "df" variable below is populated with results from the data element in the results list, with an example slice of data showing the first record.
```{r data example, echo=TRUE, message=FALSE, warning=FALSE, paged.print=TRUE}
results <- ppo_data(genus = "Quercus", fromYear = 2013, toYear = 2013, fromDay = 100, toDay = 110, termID = 'obo:PPO_0002313', limit = 10, timeLimit = 2)
df <- results$data
print(df[1:1,])
```

The readme and citation files returned by the list of results can be accessed by calling the readme and citation elements.  Note that the the file "citation_and_data_use_policies.txt" that is referred to in the readme file can be accessed using `cat(results$citation)` 
```{r readme results example}
cat(results$readme)
```

The results lists also shows the number of possible results in the results set, which is useful if the submitted query had a limit.  For example, in the query above, the limit is set to 10 but we want to know how many records were possible if the limit was not set.
```{r readme possible results example}
cat(results$number_possible)
```

### working with terms and data together
Here we will generate a data frame showing the frequency of "present" and "absent" terms for a particular query.  The query is for genus = "Quercus" and latitude > 47.  For each row in the returned data frame `ppo_data` will typically return multiple terms in the termID field, corresponding to phenological stages as defined by the PPO.  For our example, we will generate a frequency table of the number of times "present" or "absent" term occur in the entire returned dataset.  Note that the termID field returned by `ppo_data` will return "presence" terms in addition to "present" and "absent" terms, while the `ppo_terms` function only returns "present" and "absent" terms.  Thus, our frequency distribution only counts the number of "present" and "absent" terms [For an in-depth discussion of the difference between "presence", "present", and "absent", see https://www.frontiersin.org/articles/10.3389/fpls.2018.00517/full].  Finally, since termIDs are returned as URI identifiers and not easily readable text, this example maps termIDs to labels. The resulting data frame shows two columns: 1) a column of term labels, and 2) a frequency of the number of times this label appeared in the result set. 

```{r workting with terms and data example}
###############################################################################
# Generate a frequency data frame showing the number of times each termID
# is populated for genus equals "Quercus" above latitude of 47
# Note that all latitude/longitude queries need to be in the format of a
# bounding box
###############################################################################
df <- ppo_data(
  genus = "Quercus", 
  limit="10", timeLimit = 4)
# return just the termID column
t1 <- df$data[,c('termID')]
# paste each cell into one string
t2<-paste(t1, collapse = ",")
# split strings at ,
t3<-strsplit(t2, ",")
# create a frequency table as a data frame
freqFrame <- as.data.frame(table(t3))

# create a new data frame that we want to populate
resultFrame <- data.frame(
  label = character(), 
  frequency = integer(), 
  stringsAsFactors = FALSE)

###############################################################################
# Replace termIDs with labels in frequency frame
###############################################################################
# fetch "present" and "absent" terms using `ppo_terms`
termList <- ppo_terms(absent = TRUE, present = TRUE, timeLimit = 2);

# loop all "present"" and "absent" terms
if (!is.null(termList)) {
  for (term in 1:nrow(termList)) {
    termListTermID<-termList[term,'termID'];
    termListLabel<-termList[term,'label'];
    # loop all rows that have a frequency generated
    for (row in 1:nrow(freqFrame)) {
      freqFrameTermID = freqFrame[row,'t3']
      freqFrameFrequency = freqFrame[row,'Freq']
      # Populate resultFrame with matching "present" or "absent" labels.
      # In this step, we will ignore "presence" terms
      # found in the frequency frame since the ppo_terms only returns
      # "present" and "absent" terms. 
      if (freqFrameTermID == termListTermID) {
        resultFrame[nrow(resultFrame)+1,] <- c(termListLabel,freqFrameFrequency)
      }
    }
  }
} else {
  message("termList is null, likely due to a server response issue.  Try increasing the timeLimit or try again later. If the problem persists email the authors.")
}


# print results, showing term labels and a frequency count
print(resultFrame)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppo_terms.R
\name{ppo_terms}
\alias{ppo_terms}
\title{Access Terms From the Plant Phenology Ontology}
\usage{
ppo_terms(present = FALSE, absent = FALSE, timeLimit = 4)
}
\arguments{
\item{present}{(boolean) If TRUE then return all "present" phenological stages}

\item{absent}{(boolean) IF TRUE then return all "absent" phenological stages.}

\item{timeLimit}{(integer) set the limit ofthe amount of time to wait for a response}
}
\value{
data.frame
}
\description{
Access present and absent terms from the Plant Phenology Ontology
}
\details{
The ppo_terms function returns terms from the Plant Phenology Ontology (PPO).
The function only accepts parameters for "present" or "absent" terms.
The response populates a data frame with: termID, label, description, and
URI.  Use the termID values in submitting termID values to the
\code{\link{ppo_data}} function.  The label and description fields are
extracted from the Plant Phenology Ontology and are useful in
determining the proper term to query on.  The URI field contains a link to
the term itself which is useful for determining superclass and subclass
relationships for each term.
For more information on the PPO ontology itself, we suggest loading the PPO
\url{https://github.com/PlantPhenoOntology/ppo} with
protege \url{https://protege.stanford.edu/}
}
\examples{
presentTerms <- ppo_terms(present = TRUE, timeLimit = 1)

absentTerms <- ppo_terms(absent = TRUE, timeLimit = 1)
}
\keyword{lookup}
\keyword{trait}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppo_data.R
\name{ppo_data}
\alias{ppo_data}
\title{Access Data From the Global Plant Phenology Data Portal}
\usage{
ppo_data(
  genus = NULL,
  specificEpithet = NULL,
  termID = NULL,
  fromYear = NULL,
  toYear = NULL,
  fromDay = NULL,
  toDay = NULL,
  bbox = NULL,
  limit = NULL,
  timeLimit = 60
)
}
\arguments{
\item{genus}{(string) a plant genus name}

\item{specificEpithet}{(string) a plant specific epithet}

\item{termID}{(string) A single termID from the plant phenology ontology.
See the \code{\link{ppo_terms}} function for more information.}

\item{fromYear}{(integer) return data from the specified year}

\item{toYear}{(integer) return data up to and including the specified year}

\item{fromDay}{(integer) return data starting from the specified day}

\item{toDay}{(integer) return data up to and including the specified day}

\item{bbox}{(string) return data within a bounding box. Format is
\code{lat,long,lat,long} and is structured as a string.  Use this website:
http://boundingbox.klokantech.com/ to quickly grab a bbox (set format on
bottom left to csv and be sure to switch the order from
long, lat, long, lat to lat, long, lat, long).}

\item{limit}{(integer) limit returned data to a specified number of records}

\item{timeLimit}{(integer) set the limit ofthe amount of time to wait for a response}
}
\value{
Return value containing a list with the following components:
\itemize{
 \item {`data`: A data frame containing data}
 \item {`readme`: A string with information about the return package}
 \item {`citation`: A string with citation information}
 \item {`number_possible`: An integer with total possible results}
 \item {`status_code`: An integer with status code returned from server}
}
}
\description{
Access data from the global plant phenology data portal
(PPO data portal)
}
\details{
The ppo_data function returns a list containing the following information:
a readme file, citation information, a data frame with data, an integer with
the number of records returned and a status code.  The function is called with
parameters that correspond to values contained in the data itself which act
as a filter on the returned record set.
}
\examples{

r1 <- ppo_data(genus = "Quercus", termID='obo:PPO_0002313', limit=10, timeLimit = 4)

r2 <- ppo_data(fromDay = 1, toDay = 100, bbox="37,-120,38,-119", limit=10, timeLimit = 4)

my_data_frame <- r2$data
}
\keyword{data}
\keyword{download}
\keyword{phenology}
\keyword{plant}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rppo.R
\docType{package}
\name{rppo-package}
\alias{rppo}
\alias{rppo-package}
\title{Access the Global Plant Phenology Data Portal}
\description{
Access data from the global plant phenology data portal
(PPO data portal) and phenology terms from the Plant Phenology Ontology
(PPO)
}
\details{
\pkg{rppo} enables users to query the global plant phenology
data portal (PPO data portal). The PPO data portal is an aggregation
of phenology data from several different data sources.  Currently it
contains USA-NPN, NEON, and PEP725 data sources.  The PPO data portal
harvests data using the ppo-data-pipeline, with code available at
\url{https://github.com/biocodellc/ppo-data-pipeline/}.  All phenological
terms in the data portal are aligned using the Plant Phenology Ontology
(PPO), available at \url{https://github.com/PlantPhenoOntology/ppo}.


Two functions  are contained in the \pkg{rppo}:
\code{\link{ppo_terms}} allows users to discover present and absent
phenological stages while \code{\link{ppo_data}} enables
users to query the PPO data portal. The \pkg{rppo} package source code is
available at \url{https://github.com/ropensci/rppo/}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/rppo/}
  \item \url{https://github.com/ropensci/rppo/}
  \item Report bugs at \url{https://github.com/ropensci/rppo/issues}
}

}
\author{
\strong{Maintainer}: John Deck \email{jdeck88@gmail.com}

Authors:
\itemize{
  \item Brian Stucky \email{stuckyb@flmnh.ufl.edu}
  \item Ramona Walls \email{rwalls@cyverse.org}
  \item Kjell Bolmgren \email{Kjell.Bolmgren@slu.se}
  \item Ellen Denny \email{ellen@usanpn.org}
  \item Robert Guralnick \email{rguralnick@flmnh.ufl.edu}
}

}
