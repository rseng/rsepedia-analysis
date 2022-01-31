
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













---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

[![codecov](https://codecov.io/gh/ropensci/fingertipsR/branch/master/graph/badge.svg?token=MpVheRqaRo)](https://codecov.io/gh/ropensci/fingertipsR)
[![](https://badges.ropensci.org/168_status.svg)](https://github.com/ropensci/software-review/issues/168)
[![R build status](https://github.com/ropensci/fingertipsR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/fingertipsR/actions)

# fingertipsR

This is an R package to interact with Public Health England's [Fingertips](http://fingertips.phe.org.uk/) data tool. Fingertips is a major public repository of population and public health indicators for England. The site presents the information in many ways to improve accessibility for a wide range of audiences ranging from public health professionals and researchers to the general public. The information presented is a mixture of data available from other public sources, and those that are available through user access agreements with other organisations. The source of each indicator presented is available using the `indicator_metadata()` function.

This package can be used to load data from the Fingertips API into R for further use. 

## Installation

### rOpenSci

Get the latest released, stable version from rOpenSci:

```{r CRAN-install, eval=FALSE}
install.packages("fingertipsR", repos = "https://dev.ropensci.org")
```

### With remotes

You can install the latest development version from github using [remotes](https://github.com/r-lib/remotes):

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("rOpenSci/fingertipsR",
                        build_vignettes = TRUE,
                        dependencies = "suggests")
```


## Example

This is an example of a workflow for downloading data for the indicator on *Healthy Life Expectancy at Birth* from the *Public Health Outcomes Framework* profile.

The `profiles()` function presents all of the available profiles:

```{r profiles example}
library(fingertipsR)
profs <- profiles()
profs <- profs[grepl("Public Health Outcomes Framework", profs$ProfileName),]
head(profs)
```

This table shows that the `ProfileID` for the Public Health Outcomes Framework is 19. This can be used as an input for the `indicators()` function:

```{r indicators example}
profid <- 19
inds <- indicators(ProfileID = profid)
print(inds[grepl("Healthy", inds$IndicatorName), c("IndicatorID", "IndicatorName")])
```

Healthy Life Expectancy at Birth has the `IndicatorID` equal to 90362.

Finally, the data can be extracted using the `fingertips_data()` function using that `IndicatorID`:

```{r fingertips_data example}
indid <- 90362
df <- fingertips_data(IndicatorID = indid, AreaTypeID = 202)
head(df)
```

## Use

Please see the vignettes for information on use.

```{r use, eval=FALSE}
browseVignettes("fingertipsR")
```

## More information

* Please note that the 'fingertipsR' project is released with a
[Contributor Code of Conduct](https://github.com/ropensci/fingertipsR/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
* License: [GPL-3](https://opensource.org/licenses/GPL-3.0)

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Plotting healthy life expectancy and life expectancy by deprivation for English local authorities"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Plotting healthy life expectancy and life expectancy by deprivation for English local authorities}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```r
knitr::opts_chunk$set(fig.path = "charts-")
```

This worked example attempts to document a common workflow a user might follow when using the `fingertipsR` package.

`fingertipsR` provides users the ability to import data from the [Fingertips](http://fingertips.phe.org.uk/) website. Fingertips is a major repository of public health indicators in England. The site is structured in the following way:

* Profiles - these contain indicators related to a broad theme (such as a risk factor or a disease topic etc)
* Domains - these are subcategories of profiles, and break the profiles down to themes within the broader theme
* Indicators - this is the lowest level of the structure and sit within the domains. Indicators are presented at different time periods, geographies, sexes, ageband and categories.

This example demonstrates how you can plot healthy life expectancy and life expectancy by geographical regions for a given year of data that fingertips contains. So, *where to start*?

## Where to start

There is one function in the `fingertipsR` package that imports data from the Fingertips API; `fingertips_data()`. This function has the following inputs:

* IndicatorID
* AreaCode
* DomainID
* ProfileID
* AreaTypeID
* ParentAreaTypeID

At least one of *IndicatorID*, *DomainID* or *ProfileID* must be complete. These fields relate to each other as described in the introduction. *AreaTypeID* is also required, and determines the geography for which data is imported. In this case we want County and Unitary Authority level. *AreaCode* needs completing if you are extracting data for a particular area or group of areas only. *ParentAreaTypeID* requires an area type code that the *AreaTypeID* maps to at a higher level of geography. For example, County and Unitary Authorities map to a higher level of geography called Government Office Regions. These mappings can be identified using the `area_types()` function. If ignored, a *ParentAreaTypeID* will be chosen automatically.

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

![plot of chunk display_area_types](charts-display_area_types-1.png)

The table shows that the *AreaTypeID* for County and Unitary Authority level is 202. The third column, *ParentAreaTypeID*, shows the IDs of the area types that these map to. In the case of County and Unitary Authorities, these are:


| AreaTypeID|AreaTypeName                               | ParentAreaTypeID|ParentAreaTypeName                        |
|----------:|:------------------------------------------|----------------:|:-----------------------------------------|
|        202|Upper tier local authorities (4/19 - 3/20) |                6|Government Office Region                  |
|        202|Upper tier local authorities (4/19 - 3/20) |              104|PHEC 2015 new plus PHEC 2013 unchanged    |
|        202|Upper tier local authorities (4/19 - 3/20) |            10105|Depriv. decile (IMD2015, 4/19 boundaries) |
|        202|Upper tier local authorities (4/19 - 3/20) |            10113|Depriv. deciles (IMD2019)                 |
|        202|Upper tier local authorities (4/19 - 3/20) |              126|Combined authorities                      |

*ParentAreaTypeID* is 6 by default for the `fingertips_data()` function for `AreaTypeID` of 202 (this value changes if different `AreaTypeID`s are entered), so we can stick with that in this example. Use the `area_types()` function to understand more about how areas map to each other.

## Extracting the data

Finally, we can use the `fingertips_data()` function with the inputs we have determined previously.


```r
indicators <- c(90362, 90366)
data <- fingertips_data(IndicatorID = indicators,
                        AreaTypeID = 202)
```


```
## 
## 
## |  &nbsp;  | IndicatorID |      IndicatorName       | ParentCode |
## |:--------:|:-----------:|:------------------------:|:----------:|
## | **8023** |    90366    | Life expectancy at birth | E12000005  |
## | **8024** |    90366    | Life expectancy at birth | E12000006  |
## | **8025** |    90366    | Life expectancy at birth | E12000008  |
## | **8026** |    90366    | Life expectancy at birth | E12000005  |
## | **8027** |    90366    | Life expectancy at birth | E12000008  |
## | **8028** |    90366    | Life expectancy at birth | E12000005  |
## 
## Table: Table continues below
## 
##  
## 
## |  &nbsp;  |       ParentName       | AreaCode  |    AreaName    |
## |:--------:|:----------------------:|:---------:|:--------------:|
## | **8023** |  West Midlands region  | E10000028 | Staffordshire  |
## | **8024** | East of England region | E10000029 |    Suffolk     |
## | **8025** |   South East region    | E10000030 |     Surrey     |
## | **8026** |  West Midlands region  | E10000031 |  Warwickshire  |
## | **8027** |   South East region    | E10000032 |  West Sussex   |
## | **8028** |  West Midlands region  | E10000034 | Worcestershire |
## 
## Table: Table continues below
## 
##  
## 
## |  &nbsp;  |        AreaType         |  Sex   |   Age    | CategoryType | Category |
## |:--------:|:-----------------------:|:------:|:--------:|:------------:|:--------:|
## | **8023** | County & UA (4/19-3/20) | Female | All ages |      NA      |    NA    |
## | **8024** | County & UA (4/19-3/20) | Female | All ages |      NA      |    NA    |
## | **8025** | County & UA (4/19-3/20) | Female | All ages |      NA      |    NA    |
## | **8026** | County & UA (4/19-3/20) | Female | All ages |      NA      |    NA    |
## | **8027** | County & UA (4/19-3/20) | Female | All ages |      NA      |    NA    |
## | **8028** | County & UA (4/19-3/20) | Female | All ages |      NA      |    NA    |
## 
## Table: Table continues below
## 
##  
## 
## |  &nbsp;  | Timeperiod | Value | LowerCI95.0limit | UpperCI95.0limit |
## |:--------:|:----------:|:-----:|:----------------:|:----------------:|
## | **8023** | 2017 - 19  | 83.45 |      83.24       |      83.66       |
## | **8024** | 2017 - 19  | 84.25 |      84.02       |      84.47       |
## | **8025** | 2017 - 19  | 85.33 |      85.15       |      85.51       |
## | **8026** | 2017 - 19  | 83.85 |      83.58       |      84.11       |
## | **8027** | 2017 - 19  | 84.21 |        84        |      84.42       |
## | **8028** | 2017 - 19  | 83.8  |      83.54       |      84.06       |
## 
## Table: Table continues below
## 
##  
## 
## |  &nbsp;  | LowerCI99.8limit | UpperCI99.8limit | Count | Denominator | Valuenote |
## |:--------:|:----------------:|:----------------:|:-----:|:-----------:|:---------:|
## | **8023** |        NA        |        NA        |  NA   |     NA      |    NA     |
## | **8024** |        NA        |        NA        |  NA   |     NA      |    NA     |
## | **8025** |        NA        |        NA        |  NA   |     NA      |    NA     |
## | **8026** |        NA        |        NA        |  NA   |     NA      |    NA     |
## | **8027** |        NA        |        NA        |  NA   |     NA      |    NA     |
## | **8028** |        NA        |        NA        |  NA   |     NA      |    NA     |
## 
## Table: Table continues below
## 
##  
## 
## |  &nbsp;  |     RecentTrend      | ComparedtoEnglandvalueorpercentiles |
## |:--------:|:--------------------:|:-----------------------------------:|
## | **8023** | Cannot be calculated |               Similar               |
## | **8024** | Cannot be calculated |               Better                |
## | **8025** | Cannot be calculated |               Better                |
## | **8026** | Cannot be calculated |               Better                |
## | **8027** | Cannot be calculated |               Better                |
## | **8028** | Cannot be calculated |               Better                |
## 
## Table: Table continues below
## 
##  
## 
## |  &nbsp;  | ComparedtoRegionvalueorpercentiles | TimeperiodSortable | Newdata |
## |:--------:|:----------------------------------:|:------------------:|:-------:|
## | **8023** |               Better               |      20170000      |   NA    |
## | **8024** |               Better               |      20170000      |   NA    |
## | **8025** |               Better               |      20170000      |   NA    |
## | **8026** |               Better               |      20170000      |   NA    |
## | **8027** |              Similar               |      20170000      |   NA    |
## | **8028** |               Better               |      20170000      |   NA    |
## 
## Table: Table continues below
## 
##  
## 
## |  &nbsp;  | Comparedtogoal |
## |:--------:|:--------------:|
## | **8023** |       NA       |
## | **8024** |       NA       |
## | **8025** |       NA       |
## | **8026** |       NA       |
## | **8027** |       NA       |
## | **8028** |       NA       |
```

The data frame returned by `fingertips_data()` contains 26 variables.  For this exercise, we are only interested in a few of them and for the time period 2012-14:

* IndicatorID
* AreaCode
* ParentAreaName
* Sex
* Timeperiod
* Value

The data frame also contains data for the parent area, and for England, so we want to filter it to remove these too.


```r
cols <- c("IndicatorID", "AreaCode", "ParentName", "Sex", "Timeperiod", "Value")

area_type_name <- table(data$AreaType) # tally each group in the AreaType field

area_type_name <- area_type_name[area_type_name == max(area_type_name)] # pick the group with the highest frequency
area_type_name <- names(area_type_name) # retrieve the name

data <- data[data$AreaType == area_type_name &
               data$Timeperiod == "2012 - 14", cols]
```

## Plotting outputs

Using `ggplot2` it is possible to plot the outputs.


```r
library(ggplot2)
ggplot(data, aes(x = reorder(ParentName, Value, median), y = Value, col = factor(IndicatorID))) +
        geom_boxplot(data = data[data$IndicatorID == 90366, ]) +
        geom_boxplot(data = data[data$IndicatorID == 90362, ]) +
        facet_wrap(~ Sex) +
        scale_colour_manual(name = "Indicator",
                            breaks = c("90366", "90362"),
                            labels = c("Life expectancy", "Healthy life expectancy"),
                            values = c("#128c4a", "#88c857")) +
        labs(x = "Region",
             y = "Age",
             title = "Life expectancy and healthy life expectancy at birth \nfor Upper Tier Local Authorities within England regions (2012 - 2014)") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45,
                                         hjust = 1))
```

![plot of chunk life-expectancy](charts-life-expectancy-1.png)

## Other useful functions

The plot above makes use of the fields that are within the dataset by default when using the `fingertips_data()` function. There is also a `deprivation_decile()` function, which provides an indicator of deprivation for each geographical area (see `?deprivation_decile`).

Not all indicators are available for every geography. To understand how indicators are mapped to different gegoraphies, there is a function `indicator_areatypes()`.

To understand more about what comprises each indicator, there is the `indicator_metadata()` function, which provides the information on the definitions page of the Fingertips website.

Finally, the `nearest_neighbours()` function provides groups of statistically similar area for some of the geographies that are available. The geographies these are available for, and their sources, are documented within the function documentation (`?nearest_neighbours`).
---
title: "Interactively selecting indicators"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Interactively selecting indicators and identifying poorly performing areas}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```r
knitr::opts_chunk$set(fig.path = "charts-")
```

This short example introduces two new features included in the 0.1.1 version of the package:

* interactively selecting indicators
* highlighting the areas for each indicator that are deteriorating and are already statistically significantly worse than the England (or parent) value

The following libraries are needed for this vignette:

```r
library(fingertipsR)
library(ggplot2)
```

To begin with, we want to select the indicators that we want to analyse in this example. The `select_indicators()` function helps us do this:


```r
inds <- select_indicators()
```

After running the above code a browser window opens and will spend some time loading while it accesses the available indicators. Once loaded, the search bar in the top right corner allows the user to type in any searches to help locate the indicators of interest. This example will use the following indicators, selected at random, which belong to the Public Health Outcomes Framework profile. Note, the indicator IDs of the selected indicators are displayed on the left-hand side. The user can click an indicator for a second time to deselect the indicator if it has been selected unnecessarily.


```
##  [1] 90630 10101 10301 92313 10401 11401 10501 92314 11502 10601 20101 20201 20601 20301 20602 90284 90832
## [18] 90285 22001 22002 90244
```

The second function to be highlighted in this vignette is `fingertips_redred()`. This will return a data frame of the data for all of the areas that are performing significantly worse than the chosen benchmark and are deteriorating.


```r
df <- fingertips_redred(inds, AreaTypeID = 202, Comparator = "England")
```

The `geom_tile()` function from `ggplot2` can be used to visualise the poorly performing areas:


```r
df$IndicatorName <- sapply(strwrap(df$IndicatorName, 60, 
                                   simplify = FALSE), 
                           paste, collapse= "\n")
p <- ggplot(df, aes(IndicatorName, AreaName)) + 
        geom_point(colour = "darkred") + 
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45,
                                         hjust = 1,
                                         size = rel(0.85)),
              axis.text.y = element_text(size = rel(0.9))) +
        labs(y = "Upper Tier Local Authority",
             x = "Indicator")
print(p)
```

![plot of chunk redred](charts-redred-1.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_indicators.R
\name{select_indicators}
\alias{select_indicators}
\title{Select indicator}
\usage{
select_indicators()
}
\value{
A numeric vector of indicator IDs
}
\description{
Point and click method of selecting indicators and assigning them to object.
Note, this function can take up to a few minutes to run (depending on
internet connection speeds).
}
\examples{
\dontrun{
# Opens a browser window allowing the user to select indicators by their name, domain and profile
inds <- select_indicators()}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/area_types.R
\name{area_types}
\alias{area_types}
\title{Area types}
\usage{
area_types(AreaTypeName = NULL, AreaTypeID = NULL, ProfileID = NULL, path)
}
\arguments{
\item{AreaTypeName}{Character vector, description of the area type; default
is NULL}

\item{AreaTypeID}{Numeric vector, the Fingertips ID for the area type;
default is NULL}

\item{ProfileID}{Numeric vector, id of profiles of interest}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of area type ids and their descriptions
}
\description{
Outputs a data frame of area type ids, their descriptions, and how they map
to parent area types. To understand more on mappings of areas, see the Where
to start section of the Life Expectancy vignette.
}
\examples{
\dontrun{
# Returns a data frame with all levels of area and how they map to one another
area_types()

# Returns a data frame of county and unitary authority mappings
 area_types("counties")

# Returns a data frame of both counties, district
# and unitary authorities and their respective mappings
areas <- c("counties","district")
area_types(areas)

# Uses AreaTypeID to filter area types
area_types(AreaTypeID = 152)}
}
\seealso{
\code{\link{indicators}} for indicator lookups,
  \code{\link{profiles}} for profile lookups,
  \code{\link{deprivation_decile}} for deprivation decile lookups,
  \code{\link{category_types}} for category lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups,
  \code{\link{indicators_unique}} for unique indicatorids and their names,
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an area and
  \code{\link{indicator_order}} for the order indicators are presented on the
  Fingertips website within a Domain

Other lookup functions: 
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_fingertips_api}
\alias{get_fingertips_api}
\title{Retrieve data from a given Fingertips API url}
\usage{
get_fingertips_api(api_path)
}
\arguments{
\item{api_path}{string; the API url to retrieve data from}
}
\description{
Retrieve data from a given Fingertips API url
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{indicators}
\alias{indicators}
\title{Live indicators and the profiles and domains they belong to}
\usage{
indicators(ProfileID = NULL, DomainID = NULL, path)
}
\arguments{
\item{ProfileID}{Numeric vector, id of profiles of interest}

\item{DomainID}{Numeric vector, id of domains of interest}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of indicators within a profile or domain.
}
\description{
Outputs a data frame of indicators within a profile or domain
}
\examples{
\dontrun{
# Returns a complete data frame of indicators and their domains and profiles
indicators()

# Returns a data frame of all of the indicators in the Public Health Outcomes Framework
indicators(ProfileID = 19)}
}
\seealso{
\code{\link{area_types}} for area type  and their parent mappings,
  \code{\link{indicator_metadata}} for indicator metadata,
  \code{\link{profiles}} for profile lookups,
  \code{\link{deprivation_decile}} for deprivation decile lookups,
  \code{\link{category_types}} for category lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups,
  \code{\link{indicators_unique}} for unique indicatorids and their names,
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an area and
  \code{\link{indicator_order}} for the order indicators are presented on the
  Fingertips website within a Domain

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fingertips_data.R
\name{fingertips_data}
\alias{fingertips_data}
\title{Fingertips data}
\usage{
fingertips_data(
  IndicatorID = NULL,
  AreaCode = NULL,
  DomainID = NULL,
  ProfileID = NULL,
  AreaTypeID,
  ParentAreaTypeID = NULL,
  categorytype = FALSE,
  rank = FALSE,
  url_only = FALSE,
  path
)
}
\arguments{
\item{IndicatorID}{Numeric vector, id of the indicator of interest}

\item{AreaCode}{Character vector, ONS area code of area of interest}

\item{DomainID}{Numeric vector, id of domains of interest}

\item{ProfileID}{Numeric vector, id of profiles of interest. Indicator
polarity can vary between profiles therefore if using one of the comparison
fields it is recommended to complete this field as well as IndicatorID. If
IndicatorID is populated, ProfileID can be ignored or must be the same
length as IndicatorID (but can contain NAs).}

\item{AreaTypeID}{Numeric vector, the Fingertips ID for the area type. This
argument accepts "All", which returns data for all available area types for
the indicator(s), though this can take a long time to run}

\item{ParentAreaTypeID}{Numeric vector, the comparator area type for the data
extracted; if NULL the function will use the first record for the specified
`AreaTypeID` from the area_types() function}

\item{categorytype}{TRUE or FALSE, determines whether the final table
includes categorytype data where it exists. Default to FALSE}

\item{rank}{TRUE or FALSE, the rank of the area compared to other areas for
that combination of indicator, sex, age, categorytype and category along
with the indicator's polarity. 1 is lowest NAs will be bottom and ties will
return the average position. The total count of areas with a non-NA value
are returned also in AreaValuesCount}

\item{url_only}{TRUE or FALSE, return only the url of the api call as a
character vector}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of data extracted from the Fingertips API
}
\description{
Outputs a data frame of data from
\href{https://fingertips.phe.org.uk/}{Fingertips}. Note, this function can
take up to a few minutes to run (depending on internet connection speeds and
parameter selection).
}
\details{
Note, polarity of an indicator is not automatically returned (eg,
  whether a low value is good, bad or neither). Use the rank field for this
  to be returned (though it adds a lot of time to the query)
}
\examples{
\dontrun{
# Returns data for the two selected domains at county and unitary authority geography
doms <- c(1000049,1938132983)
fingdata <- fingertips_data(DomainID = doms, AreaTypeID = 202)

# Returns data at local authority district geography (AreaTypeID = 101)
# for the indicator with the id 22401
fingdata <- fingertips_data(22401, AreaTypeID = 101)

# Returns same indicator with different comparisons due to indicator polarity
# differences between profiles on the website
# It is recommended to check the website to ensure consistency between your
# data extract here and the polarity required
fingdata <- fingertips_data(rep(90282,2),
                            ProfileID = c(19,93),
                            AreaTypeID = 202,
                            AreaCode = "E06000008")
fingdata <- fingdata[order(fingdata$TimeperiodSortable, fingdata$Sex),]

# Returns data for all available area types for an indicator
fingdata <- fingertips_data(10101, AreaTypeID = "All")}
}
\seealso{
Other data extract functions: 
\code{\link{fingertips_redred}()}
}
\concept{data extract functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{indicators_unique}
\alias{indicators_unique}
\title{Live indicators}
\usage{
indicators_unique(ProfileID = NULL, DomainID = NULL, path)
}
\arguments{
\item{ProfileID}{Numeric vector, id of profiles of interest}

\item{DomainID}{Numeric vector, id of domains of interest}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of indicator ids and names
}
\description{
Outputs a data frame of indicators (their id and name only). Note, this
function can take up to a few minutes to run (depending on internet
connection speeds)
}
\examples{
\dontrun{
indicators_unique(ProfileID = 21)}
}
\seealso{
\code{\link{indicators}} for indicators and their parent domains and
  profiles, \code{\link{area_types}} for area type  and their parent
  mappings, \code{\link{indicator_metadata}} for indicator metadata and
  \code{\link{profiles}} for profile lookups and
  \code{\link{deprivation_decile}} for deprivation decile lookups and
  \code{\link{category_types}} for category lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups and
  \code{\link{indicator_order}} for the order indicators are presented on the
  Fingertips website within a Domain

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profiles.R
\name{profiles}
\alias{profiles}
\title{Live profiles}
\usage{
profiles(ProfileID = NULL, ProfileName = NULL, path)
}
\arguments{
\item{ProfileID}{Numeric vector, id of profiles of interest}

\item{ProfileName}{Character vector, full name of profile(s)}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of live profile ids and names along with their domain
  names and ids.
}
\description{
Outputs a data frame of live profiles that data are available for in
Fingertips \url{http://fingertips.phe.org.uk/}
}
\examples{
\dontrun{
# Returns a complete data frame of domains and their profiles
profiles()

# Returns a data frame of all of the domains in the Public Health Outcomes Framework
profiles(ProfileName = "Public Health Outcomes Framework")}
}
\seealso{
\code{\link{area_types}} for area type  and their parent mappings,
  \code{\link{indicators}} for indicator lookups,
  \code{\link{indicator_metadata}} for indicator metadata,
  \code{\link{deprivation_decile}} for deprivation decile lookups,
  \code{\link{category_types}} for category lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups,
  \code{\link{indicators_unique}} for unique indicatorids and their names,
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an area and
  \code{\link{indicator_order}} for the order indicators are presented on the
  Fingertips website within a Domain

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fingertipsR.R
\docType{package}
\name{fingertipsR}
\alias{fingertipsR}
\title{fingertipsR: A package for extracting the data behind the Fingertips website
(\url{https://fingertips.phe.org.uk/})}
\description{
The fingertipsR package provides two categories of important functions: lookup
and data extract.
}
\section{Lookup functions}{
 The lookup functions are to provide users the
  ability to understand the ID inputs for the data extract functions.
}

\section{Data extract functions}{
 Using ID codes as inputs, the data extract
  functions allow the user to extract data from the Fingertips API.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicator_metadata.R
\name{indicator_metadata}
\alias{indicator_metadata}
\title{Indicator metadata}
\usage{
indicator_metadata(IndicatorID = NULL, DomainID = NULL, ProfileID = NULL, path)
}
\arguments{
\item{IndicatorID}{Numeric vector, id of the indicator of interest. Also accepts "All".}

\item{DomainID}{Numeric vector, id of domains of interest}

\item{ProfileID}{Numeric vector, id of profiles of interest. Indicator
polarity can vary between profiles therefore if using one of the comparison
fields it is recommended to complete this field as well as IndicatorID. If
IndicatorID is populated, ProfileID can be ignored or must be the same
length as IndicatorID (but can contain NAs).}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
The metadata associated with each indicator/domain/profile identified
}
\description{
Outputs a data frame containing the metadata for selected indicators. Note, this
function can take up to a few minutes to run (depending on internet
connection speeds)
}
\examples{
\dontrun{
# Returns metadata for indicator ID 90362 and 1107
indicatorIDs <- c(90362, 1107)
indicator_metadata(indicatorIDs)

# Returns metadata for the indicators within the domain 1000101
indicator_metadata(DomainID = 1000101)

# Returns metadata for the indicators within the profile with the ID 129
indicator_metadata(ProfileID = 129)}
}
\seealso{
\code{\link{indicators}} for indicator lookups,
  \code{\link{profiles}} for profile lookups,
  \code{\link{deprivation_decile}} for deprivation lookups,
  \code{\link{area_types}} for area types and their parent mappings,
  \code{\link{category_types}} for category lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups,
  \code{\link{indicators_unique}} for unique indicatorids and their names,
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an area and
  \code{\link{indicator_order}} for the order indicators are presented on the
  Fingertips website within a Domain

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fingertips_ensure_api_available}
\alias{fingertips_ensure_api_available}
\title{Check if the given Fingertips API endpoint is available}
\usage{
fingertips_ensure_api_available(endpoint = fingertips_endpoint())
}
\arguments{
\item{endpoint}{string, the API base URL to check}
}
\value{
\code{TRUE} if the API is available, otherwise \code{stop()} is called.
}
\description{
Check if the given Fingertips API endpoint is available
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhancements.R
\name{fingertips_redred}
\alias{fingertips_redred}
\title{Red significance and red trend}
\usage{
fingertips_redred(Comparator = "England", ...)
}
\arguments{
\item{Comparator}{String, either "England" or "Parent" to determine which
field to compare the spot value significance to}

\item{...}{Parameters provided to fingertips_data()}
}
\value{
A data frame of data extracted from
the Fingertips API
}
\description{
Filters data returned by the fingertips_data function for values for areas
that are trending statistically significantly worse and the spot value is
significantly worse than the comparator (England or Parent) value in the
latest year of that indicator
}
\examples{
\dontrun{
# Returns data for the two selected domains at county and unitary authority geography
reddata <- fingertips_redred(ProfileID = 26, AreaTypeID = 102)}
}
\seealso{
Other data extract functions: 
\code{\link{fingertips_data}()}
}
\concept{data extract functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/area_types.R
\name{nearest_neighbours}
\alias{nearest_neighbours}
\title{Nearest neighbours}
\usage{
nearest_neighbours(AreaCode, AreaTypeID, measure, path)
}
\arguments{
\item{AreaCode}{Character vector, ONS area code of area of interest}

\item{AreaTypeID}{AreaTypeID of the nearest neighbours (see
\code{\link{nearest_neighbour_areatypeids}}) for available IDs}

\item{measure}{deprecated. Previously a string; when AreaTypeID = 102 measure
must be either "CIPFA" for CIPFA local authority nearest neighbours or
"CSSN" for Children's services statistical neighbours}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A character vector of area codes
}
\description{
Outputs a character vector of similar areas for given area. Currently returns
similar areas for Clinical Commissioning Groups (old and new) based on
\href{https://www.england.nhs.uk/publication/similar-10-ccg-explorer-tool/}{NHS
England's similar CCG explorer tool} or lower and upper tier local
authorities based on
\href{https://www.cipfastats.net/resources/nearestneighbours/}{CIPFA's
Nearest Neighbours Model} or upper tier local authorities based on
\href{https://www.gov.uk/government/publications/local-authority-interactive-tool-lait}{Children's
services statistical neighbour benchmarking tool}
}
\examples{
\dontrun{
nearest_neighbours(AreaCode = "E38000002", AreaTypeID = 154)}
}
\seealso{
\code{\link{nearest_neighbour_areatypeids}} for the AreaTypeIDs
  available for this function

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fingertips_endpoint}
\alias{fingertips_endpoint}
\title{Get the default fingertips API endpoint}
\usage{
fingertips_endpoint()
}
\value{
A character string with the HTTP URL of the Fingertips API
}
\description{
Get the default fingertips API endpoint
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/area_types.R
\name{category_types}
\alias{category_types}
\title{Category types}
\usage{
category_types(path)
}
\arguments{
\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of category type ids and their descriptions
}
\description{
Outputs a data frame of category type ids, their name (along with a short name)
}
\examples{
\dontrun{
# Returns the deprivation category types
cats <- category_types()
cats[cats$CategoryTypeId == 1,]}
}
\seealso{
\code{\link{indicators}} for indicator lookups,
  \code{\link{profiles}} for profile lookups,
  \code{\link{deprivation_decile}} for deprivation decile lookups,
  \code{\link{area_types}} for area type lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups,
  \code{\link{indicators_unique}} for unique indicatorids and their names,
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an area and
  \code{\link{indicator_order}} for the order indicators are presented on the
  Fingertips website within a Domain

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/area_types.R
\name{indicator_areatypes}
\alias{indicator_areatypes}
\title{Area types by indicator}
\usage{
indicator_areatypes(IndicatorID, AreaTypeID, path)
}
\arguments{
\item{IndicatorID}{integer; the Indicator ID (can be ignored or of length 1).
Takes priority over AreaTypeID if both are entered}

\item{AreaTypeID}{integer; the Area Type ID (can be ignored or of length 1)}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of indicator ids and area type ids
}
\description{
Outputs a data frame of indicator ids and the area type ids that exist for
that indicator
}
\examples{
\dontrun{
indicator_areatypes(IndicatorID = 10101)}
}
\seealso{
\code{\link{indicators}} for indicator lookups,
  \code{\link{profiles}} for profile lookups,
  \code{\link{deprivation_decile}} for deprivation decile lookups,
  \code{\link{area_types}} for area type lookups,
  \code{\link{category_types}} for category type lookups,
  \code{\link{indicators_unique}} for unique indicatorids and their names,
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an area and
  \code{\link{indicator_order}} for the order indicators are presented on the
  Fingertips website within a Domain

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhancements.R
\name{fingertips_stats}
\alias{fingertips_stats}
\title{High level statistics on Fingertips data}
\usage{
fingertips_stats()
}
\value{
A string that summarises the high level statistics of indicators and
  profiles in Fingertips
}
\description{
A sentence that summarises the number of indicators, unique indicators and
profiles
}
\examples{
\dontrun{
# Returns a sentence describing number of indicators and profiles in Fingertips
fingertips_stats()}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{add_timestamp}
\alias{add_timestamp}
\title{Add timestamp onto end of api url to prevent caching issues}
\usage{
add_timestamp(api_path)
}
\arguments{
\item{api_path}{string; the API url to retrieve data from}
}
\description{
Add timestamp onto end of api url to prevent caching issues
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{indicator_order}
\alias{indicator_order}
\title{Indicator order number}
\usage{
indicator_order(DomainID, AreaTypeID, ParentAreaTypeID, path)
}
\arguments{
\item{DomainID}{Numeric vector, id of domains of interest}

\item{AreaTypeID}{Numeric vector, the Fingertips ID for the area type. This
argument accepts "All", which returns data for all available area types for
the indicator(s), though this can take a long time to run}

\item{ParentAreaTypeID}{Numeric vector, the comparator area type for the data
extracted; if NULL the function will use the first record for the specified
`AreaTypeID` from the area_types() function}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A data frame of indicator ids and sequence number
}
\description{
Outputs a tibble of indicator ids and their sequence number for the provided
domain and area type. This enables the user to order the indicators as they
are ordered on the Fingertips website.
}
\examples{
\dontrun{
indicator_order(DomainID = 1938133161, AreaTypeID = 102, ParentAreaTypeID = 6)}
}
\seealso{
\code{\link{indicators}} for indicators and their parent domains and profiles,
  \code{\link{area_types}} for area type and their parent mappings,
  \code{\link{indicator_metadata}} for indicator metadata,
  \code{\link{profiles}} for profile lookups,
  \code{\link{deprivation_decile}} for deprivation decile lookups,
  \code{\link{category_types}} for category lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups and
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an area

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{deprivation_decile}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/area_types.R
\name{nearest_neighbour_areatypeids}
\alias{nearest_neighbour_areatypeids}
\title{Nearest neighbours area type ids}
\usage{
nearest_neighbour_areatypeids()
}
\value{
table of AreaTypeIDs
}
\description{
Outputs a table of AreaTypeIDs available for the nearest_neighbour function
}
\examples{
\dontrun{
nearest_neighbour_areatypeids()}
}
\seealso{
\code{\link{nearest_neighbours}} to access the geogaphy codes of the
  nearest neighbours for a locality
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fingertips_deframe}
\alias{fingertips_deframe}
\title{fingertips_deframe}
\usage{
fingertips_deframe(data)
}
\arguments{
\item{data}{list whose first item is a vector of names, and second item is a
list. Items 1 and 2 must be equal length}
}
\description{
mimic tibble::deframe() without needing to import the function
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprivation_decile.R
\name{deprivation_decile}
\alias{deprivation_decile}
\title{Deprivation deciles}
\usage{
deprivation_decile(AreaTypeID, Year = 2019, path)
}
\arguments{
\item{AreaTypeID}{Integer value; this function uses the IndicatorIDs 91872,
93275 and 93553, please use the \code{indicator_areatypes()} function to
see what AreaTypeIDs are available}

\item{Year}{Integer value, representing the year of IMD release to be
applied, limited to 2015 or 2019}

\item{path}{String; Fingertips API address. Function will default to the
correct address}
}
\value{
A lookup table providing deprivation decile and area code
}
\description{
Outputs a data frame allocating deprivation decile to  area code based on the
Indices of Multiple Deprivation (IMD) produced by Department of Communities
and Local Government
}
\details{
This function uses the fingertips_data function to filter for the
  Index of multiple deprivation score for the year and area supplied, and
  returns the area code, along with the score and the deprivation decile,
  which is calculated using the ntile function from dplyr
}
\examples{
\dontrun{
# Return 2019 deprivation scores for Sustainability and Transformation Footprints
deprivation_decile(120, 2019)}
}
\seealso{
\code{\link{indicators}} for indicator lookups,
  \code{\link{profiles}} for profile lookups,
  \code{\link{indicator_metadata}} for the metadata for each indicator,
  \code{\link{area_types}} for area types and their parent mappings,
  \code{\link{category_types}} for category lookups,
  \code{\link{indicator_areatypes}} for indicators by area types lookups,
  \code{\link{indicators_unique}} for unique indicatorids and their names,
  \code{\link{nearest_neighbours}} for a vector of nearest neighbours for an
  area and \code{\link{indicator_order}} for the order indicators are
  presented on the Fingertips website within a Domain

Other lookup functions: 
\code{\link{area_types}()},
\code{\link{category_types}()},
\code{\link{indicator_areatypes}()},
\code{\link{indicator_metadata}()},
\code{\link{indicator_order}()},
\code{\link{indicators_unique}()},
\code{\link{indicators}()},
\code{\link{nearest_neighbours}()},
\code{\link{profiles}()}
}
\concept{lookup functions}
