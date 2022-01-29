
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![codecov](https://codecov.io/gh/ropensci/fingertipsR/branch/master/graph/badge.svg?token=MpVheRqaRo)](https://codecov.io/gh/ropensci/fingertipsR)
[![](https://badges.ropensci.org/168_status.svg)](https://github.com/ropensci/software-review/issues/168)
[![R build
status](https://github.com/ropensci/fingertipsR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/fingertipsR/actions)

# fingertipsR

This is an R package to interact with Public Health England’s
[Fingertips](http://fingertips.phe.org.uk/) data tool. Fingertips is a
major public repository of population and public health indicators for
England. The site presents the information in many ways to improve
accessibility for a wide range of audiences ranging from public health
professionals and researchers to the general public. The information
presented is a mixture of data available from other public sources, and
those that are available through user access agreements with other
organisations. The source of each indicator presented is available using
the `indicator_metadata()` function.

This package can be used to load data from the Fingertips API into R for
further use.

## Installation

### rOpenSci

Get the latest released, stable version from rOpenSci:

``` r
install.packages("fingertipsR", repos = "https://dev.ropensci.org")
```

### With remotes

You can install the latest development version from github using
[remotes](https://github.com/r-lib/remotes):

``` r
# install.packages("remotes")
remotes::install_github("rOpenSci/fingertipsR",
                        build_vignettes = TRUE,
                        dependencies = "suggests")
```

## Example

This is an example of a workflow for downloading data for the indicator
on *Healthy Life Expectancy at Birth* from the *Public Health Outcomes
Framework* profile.

The `profiles()` function presents all of the available profiles:

``` r
library(fingertipsR)
profs <- profiles()
profs <- profs[grepl("Public Health Outcomes Framework", profs$ProfileName),]
head(profs)
#> # A tibble: 6 x 4
#>   ProfileID ProfileName                        DomainID DomainName              
#>       <int> <chr>                                 <int> <chr>                   
#> 1        19 Public Health Outcomes Framework    1000049 A. Overarching indicato~
#> 2        19 Public Health Outcomes Framework    1000041 B. Wider determinants o~
#> 3        19 Public Health Outcomes Framework    1000042 C. Health improvement   
#> 4        19 Public Health Outcomes Framework    1000043 D. Health protection    
#> 5        19 Public Health Outcomes Framework    1000044 E. Healthcare and prema~
#> 6        19 Public Health Outcomes Framework 1938132983 Supporting information
```

This table shows that the `ProfileID` for the Public Health Outcomes
Framework is 19. This can be used as an input for the `indicators()`
function:

``` r
profid <- 19
inds <- indicators(ProfileID = profid)
print(inds[grepl("Healthy", inds$IndicatorName), c("IndicatorID", "IndicatorName")])
#> # A tibble: 2 x 2
#>   IndicatorID IndicatorName                          
#>         <int> <fct>                                  
#> 1       90362 A01a - Healthy life expectancy at birth
#> 2       93505 A01a - Healthy life expectancy at 65
```

Healthy Life Expectancy at Birth has the `IndicatorID` equal to 90362.

Finally, the data can be extracted using the `fingertips_data()`
function using that `IndicatorID`:

``` r
indid <- 90362
df <- fingertips_data(IndicatorID = indid, AreaTypeID = 202)
head(df)
#>   IndicatorID                    IndicatorName ParentCode ParentName  AreaCode
#> 1       90362 Healthy life expectancy at birth       <NA>       <NA> E92000001
#> 2       90362 Healthy life expectancy at birth       <NA>       <NA> E92000001
#> 3       90362 Healthy life expectancy at birth  E92000001    England E12000001
#> 4       90362 Healthy life expectancy at birth  E92000001    England E12000002
#> 5       90362 Healthy life expectancy at birth  E92000001    England E12000003
#> 6       90362 Healthy life expectancy at birth  E92000001    England E12000004
#>                          AreaName AreaType    Sex      Age CategoryType
#> 1                         England  England   Male All ages         <NA>
#> 2                         England  England Female All ages         <NA>
#> 3               North East region   Region   Male All ages         <NA>
#> 4               North West region   Region   Male All ages         <NA>
#> 5 Yorkshire and the Humber region   Region   Male All ages         <NA>
#> 6            East Midlands region   Region   Male All ages         <NA>
#>   Category Timeperiod    Value LowerCI95.0limit UpperCI95.0limit
#> 1     <NA>  2009 - 11 63.02647         62.87787         63.17508
#> 2     <NA>  2009 - 11 64.03794         63.88135         64.19453
#> 3     <NA>  2009 - 11 59.71114         59.19049         60.23179
#> 4     <NA>  2009 - 11 60.76212         60.39880         61.12544
#> 5     <NA>  2009 - 11 60.84033         60.38649         61.29417
#> 6     <NA>  2009 - 11 62.60207         62.07083         63.13332
#>   LowerCI99.8limit UpperCI99.8limit Count Denominator Valuenote RecentTrend
#> 1               NA               NA    NA          NA      <NA>        <NA>
#> 2               NA               NA    NA          NA      <NA>        <NA>
#> 3               NA               NA    NA          NA      <NA>        <NA>
#> 4               NA               NA    NA          NA      <NA>        <NA>
#> 5               NA               NA    NA          NA      <NA>        <NA>
#> 6               NA               NA    NA          NA      <NA>        <NA>
#>   ComparedtoEnglandvalueorpercentiles ComparedtoRegionvalueorpercentiles
#> 1                        Not compared                       Not compared
#> 2                        Not compared                       Not compared
#> 3                               Worse                       Not compared
#> 4                               Worse                       Not compared
#> 5                               Worse                       Not compared
#> 6                             Similar                       Not compared
#>   TimeperiodSortable Newdata Comparedtogoal
#> 1           20090000    <NA>           <NA>
#> 2           20090000    <NA>           <NA>
#> 3           20090000    <NA>           <NA>
#> 4           20090000    <NA>           <NA>
#> 5           20090000    <NA>           <NA>
#> 6           20090000    <NA>           <NA>
```

## Use

Please see the vignettes for information on use.

``` r
browseVignettes("fingertipsR")
```

## More information

-   Please note that the ‘fingertipsR’ project is released with a
    [Contributor Code of
    Conduct](https://github.com/ropensci/fingertipsR/blob/master/CODE_OF_CONDUCT.md).
    By contributing to this project, you agree to abide by its terms.
-   License: [GPL-3](https://opensource.org/licenses/GPL-3.0)

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# fingertipsR 1.0.8

* Fixed issue introduced to `fingertips_data()` where `AreaTypeID = "All"`
* More `AreaTypeID`s available for `nearest_neighbour()`
* `nearest_neighbour_areatypeids()` added

# fingertipsR 1.0.7

* No change to previous version

# fingertipsR 1.0.6 (2021-05-17)

* fixed warning message for multiple IndicatorIDs passed to indicator_metadata()
* fixed bug in indicator_metadata(IndicatorID = "All") - thanks Luke Bradley for pointing it out
* indicator_metadata now includes `Impact of COVID19 field`

# fingertipsR 1.0.5 (2020-09-16)

* `indicator_metadata()` accepts `IndicatorID = "All"` 
* GitHub actions added

# fingertipsR 1.0.4 (2020-06-06)

* no change for users

# fingertipsR 1.0.3 (2020-04-16)

* fixed `fingertips_data()` bug introduced in 1.0.2 where multiple IndicatorIDs and ProfileIDs are provided

# fingertipsR 1.0.2 (2020-02-09)

* Bug fix for `fingertips_data()` where `AreaTypeID = "All"

# fingertipsR 1.0.1 (2020-01-28)

* Improved flexibility of `deprivation_decile()` so new deprivation data can be drawn from website without being added to the package
* url_only argument added to `fingertips_data()` function, allowing user to retrieve the API url(s) used to download their data

# fingertipsR 1.0.0 (2020-01-08)

* default value for `AreaTypeID` in `fingertips_data()` and `deprivation_decile()` removed (along with a helpful error message if `AreaTypeID` is missing)
* progress bar added when AreaTypeID = "All" in `fingertips_data()`
* fixed issue where GP deprivation_decile returned data with 0 records
* added more allowable AreaTypes to `deprivation_decile()` function
* added 2019 year to `deprivation_decile()`
* `AreaTypeID = "All"` in `fingertips_data()` function now faster and accurate for specific profile

 
# fingertipsR 0.2.9 (2019-09-25)

* fixed issue caused by v0.2.8 for users behind organisational firewall

# fingertipsR 0.2.8 (2019-09-09)

* functions will provide informative message when they fail because of no response from the API

# fingertipsR 0.2.7 (2019-07-18)
 * Added option for `AreaTypeID = "All"` in `fingertips_data()`
 * Added `ProfileID` argument to `area_types()` function
 * `category_types()` returns field called `CategoryType` which is joinable to `fingertips_data()` output when `categorytype = TRUE`
 * Added a retry function to handle occasions when the API times out when it shouldn't

# fingertipsR 0.2.6 (2019-06-07)

* fixed `indicator_order()` function

* `nearest_neighbours()` now allows AreaTypeID 154 (CCG unchanged plus new 2018), but no longer allows AreaTypeID 153

* `nearest_neighbours()` function has been fixed

# fingertipsR 0.2.4 (2019-05-19)

* no changes from a user perspective to previous versions

# fingertipsR 0.2.3 (2019-05-14)

* no changes from a user perspective to previous versions

# fingertipsR 0.2.2 (2019-04-08)

* no changes from a user perspective to previous versions

# fingertipsR 0.2.1 (2019-03-07)

 * `deprivation_decile()` for `AreaTypeID = 7` (General Practice) now only contains 2015 deprivation deciles. 2010 to 2012 have been removed

* bug fix around entering vectors of `IndicatorID`s and `ProfileID`s into `fingertips_data()` function

# fingertipsR 0.2.0 (2018-11-12)

* `select_indicators()` fixed for selecting more than one indicator

* `deprivation_decile()` now includes MSOA (AreaTypeID = 3)

* removed deprecated `inequalities` argument from `fingertips_data()`

* increased speed of `indicators()` function (and therefore `select_indicators()` and other functions that rely on it)

* IndicatorName field in `indicators()` table returns short name rather than long name

# fingertipsR 0.1.9 (2018-08-31)

* New field "Compared to goal" in `fingertips_data()`

* New field name "Compared to [AreaType]" in `fingertips_data()`

* Improved tests

# fingertipsR 0.1.8 (2018-07-05)

* fingertips_data function adapted for "New data" field in API

# fingertipsR 0.1.7 (27/05/2018)

* nearest_neighbours() modified to include `measure` parameter - so user can get 15 CIPFA nn for upper tier local authorities

* data retrieval functions now use local proxy settings

# fingertipsR 0.1.6 (03/05/2018)

* nearest_neighbours() function added

* indicator_order() function added

# fingertipsR 0.1.5 (06/02/2018)

* corrected fingertips_stats to give accurate stats

# fingertipsR 0.1.4 (03/02/2018)

* modifications to the fingertipsR paper

* badges added to README

* package approved by ropensci

* fingertips_stats function added to give high level statistics of indicators in Fingertips

* indicator_areatypes now links to API rather than built in dataset

* indicators_unique function provides unique table of indicators

# fingertipsR 0.1.3 (5/10/2017)

* API structure updated to include 99.8 and 95 confidence intervals. Reflected in the outputs of `fingertips_data`. **NOTE** earlier versions of the package will not work anymore because of the underlying change in the API structure

# fingertipsR 0.1.2 (27/9/2017)

* fixed issue with rank and some `fingertips_data` queries

* removed dependency on tidyjson as a result of its removal from CRAN

# fingertipsR 0.1.1 (7/9/2017)

* `select_indicators()` allows user to point and click to select indicators

* stringsAsFactors parameter available in `fingertips_data()`

* automatically filter for `CategoryType = FALSE` in `fingertips_data()` - this can be set to `TRUE` if needed

* rank of area and polarity of indicator returned from `fingertips_data()` where `rank = TRUE` (polarity can also be found in `indicator_metadata()`)

* `fingertips_redred` highlights which areas are statistically different to comparator *and* trending in the wrong direction

* `category_types()` lookup function to support ordering where categories exist (eg, deprivation decile)

* `areatypes_by_indicators()` to help users determine which indicators are available for each area type (and vice versa)

* A new vignette demonstrating how some of the new functions can be used

# fingertipsR version 0.1.0 (17/6/2017)

This package allows the user to retrieve tables of:

* indicators, domains and profiles and their relationships to each other
* data related to indicators for geographies where the data are already available on the Fingertips website
* indicator metadata
* deprivation data for geographies that are available on the Fingertips website
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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
# CONTRIBUTING #

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the `fingertipsR` project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
This is a new submission.

This package provides data from the Fingertips API owned by Public Health England.

It was archived from CRAN on 2020-08-26 because of some tests failing. These tests were nothing to do with package working as it was meant to, but rather because data would be removed from the API. The tests that now require data from the API make use of the `skip_on_cran()` function, and the package uses GitHub actions to test these functions each month once the data are updated on the API.

## Test 

* local Windows 10 install, R 4.0.2
* used Travis to check on Linux (2020-09-08)

## R CMD check results

There were no NOTEs, WARNINGs or ERRORs.

## Downstream dependencies

No errors with downstream dependencies
Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://stackoverflow.com/> or <https://community.rstudio.com/>.

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](https://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

---

Brief description of the problem

```r
# insert reprex here
```
---
title: "Plotting healthy life expectancy and life expectancy by deprivation for English local authorities"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Life expectancy by deprivation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This worked example attempts to document a common workflow a user might follow when using the `fingertipsR` package.

`fingertipsR` provides users the ability to import data from the [Fingertips](http://fingertips.phe.org.uk/) website. Fingertips is a major repository of public health indicators in England. The site is structured in the following way:

* Profiles - these contain indicators related to a broad theme (such as a risk factor or a disease topic etc)
* Domains - these are subcategories of profiles, and break the profiles down to themes within the broader theme
* Indicators - this is the lowest level of the structure and sit within the domains. Indicators are presented at different time periods, geographies, sexes, ageband and categories.

This example demonstrates how you can plot healthy life expectancy and life expectancy by geographical regions for a given year of data that fingertips contains. So, *where to start*?

## Where to start

There is one function in the `fingertipsR` package that extracts data from the Fingertips API: `fingertips_data()`. This function has the following inputs:

* IndicatorID
* AreaCode
* DomainID
* ProfileID
* AreaTypeID
* ParentAreaTypeID1

At least one of *IndicatorID*, *DomainID* or *ProfileID* must be complete. These fields relate to each other as described in the introduction. *AreaTypeID* is also required, and determines the geography for which data is extracted. In this case we want County and Unitary Authority level. *AreaCode* needs completing if you are extracting data for a particular area or group of areas only. *ParentAreaTypeID* requires an area type code that the *AreaTypeID* maps to at a higher level of geography. For example, County and Unitary Authorities map to a higher level of geography called Government Office Regions. These mappings can be identified using the `area_types()` function. If ignored, a *ParentAreaTypeID* will be chosen automatically.

Therefore, the inputs to the `fingertips_data` function that we need to find out are the ID codes for:

* IndicatorID 
* AreaTypeID
* ParentAreaTypeID

We need to begin by calling the `fingertipsR` package: 

```r
library(fingertipsR)
```

## IndicatorID

There are two indicators we are interested in for this exercise. Without consulting the [Fingertips website](https://fingertips.phe.org.uk/  "Fingertips"), we know approximately what they are called:

* Healthy life expectancy
* Life expectancy

We can use the `indicators()` function to return a list of all the indicators within Fingertips. We can then filter the name field for the term *life expectancy* (note, the IndicatorName field has been converted to lower case in the following code chunk to ensure matches will not be overlooked as a result of upper case letters).


```r
inds <- indicators_unique()
life_expectancy <- inds[grepl("life expectancy", tolower(inds$IndicatorName)),]
```


| IndicatorID|IndicatorName                                                                       |
|-----------:|:-----------------------------------------------------------------------------------|
|       90362|Healthy life expectancy at birth                                                    |
|       90366|Life expectancy at birth                                                            |
|       90825|Inequality in healthy life expectancy at birth ENGLAND                              |
|       91102|Life expectancy at 65                                                               |
|       92031|Inequality in healthy life expectancy at birth LA                                   |
|       92901|Inequality in life expectancy at birth                                              |
|       93190|Inequality in life expectancy at 65                                                 |
|       93505|Healthy life expectancy at 65                                                       |
|       93523|Disability-free life expectancy at 65                                               |
|       93562|Disability-free life expectancy at birth                                            |
|         650|Life expectancy - MSOA based                                                        |
|       93249|Disability free life expectancy, (Upper age band 85+)                               |
|       93283|Life expectancy at birth, (upper age band 90+)                                      |
|       93285|Life expectancy at birth, (upper age band 85+)                                      |
|       93298|Healthy life expectancy, (upper age band 85+)                                       |
|       92641|Life expectancy at 75 (SPOT: NHSOD 1b)                                              |
|       90365|Gap in life expectancy at birth between each local authority and England as a whole |

The two indicators we are interested in from this table are:

* 90362
* 90366

## AreaTypeID

We can work out what the *AreaTypeID* codes we need using the function `area_types()`. We've decided that we want to produce the graph at County and Unitary Authority level. From the section [Where to start] we need codes for *AreaTypeID* and *ParentAreaTypeID.*


```r
areaTypes <- area_types()
```













