
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cde <img src="https://docs.ropensci.org/cde/reference/figures/logo.png" align="right" height=140/>

[![R-CMD-check](https://github.com/ropensci/cde/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cde/actions)
[![codecov](https://codecov.io/gh/ropensci/cde/branch/master/graph/badge.svg?token=F4R6nEywTx)](https://codecov.io/gh/ropensci/cde)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://badges.ropensci.org/284_status.svg)](https://github.com/ropensci/onboarding/issues/284)
[![DOI](https://zenodo.org/badge/92712854.svg)](https://zenodo.org/badge/latestdoi/92712854)
[![status](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589/status.svg)](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589)
[![CRAN
Version](http://www.r-pkg.org/badges/version/cde)](http://www.r-pkg.org/pkg/cde)
[![](http://cranlogs.r-pkg.org/badges/cde)](http://cran.rstudio.com/web/packages/cde/index.html)

## Introduction

**NOTE: The EA has recently changed the format of their API, so none of
the download functions in the package are currently working. This is
being worked on and a new release will be made once fixed.**

Within Europe, the [Water Framework
Directive](http://ec.europa.eu/environment/water/water-framework/index_en.html)
(WFD) sets EU-wide standards for how the quality of surface- and
ground-waters across Europe is assessed and classified. Assessment of
quality using the WFD is based on a range of elements that vary
depending on the type of water being assessed and are combined to give
an overall classification of waterbodies into five classes (High, Good,
Moderate, Poor and Bad) for surface waters and two classes (Good or
Poor) for groundwaters.

In the UK the Environment Agency (EA) is the competent authority
responsible for monitoring and assessment of water quality within
England. The EA have made the reporting data relating to the
requirements of the WFD available via the Catchment Data Explorer (CDE)
website, <https://environment.data.gov.uk/catchment-planning/>.

`cde` is a package for R which facilitates searching and download of the
WFD reporting data for all waterbodies from the EA CDE website. The
ability to access these data from within the R environment allows for
efficient collation and interrogation of data and reproducible analysis
of trends or patterns in water quality and pressures on waterbodies
across England. There are also some inconsistencies in the way in which
the data are structured within the original CDE website; `cde` provides
consistently named and structured output which facilitates further
analysis.

The types of data that can be downloaded are: WFD status classification
data, Reasons for Not Achieving Good (RNAG) status, objectives set for
waterbodies, measures put in place to improve water quality and details
of associated protected areas.

The CDE data are made available under the [Open Government Licence
v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/)
and use of the data accessed by and contained within this package
implies acceptance of these licence conditions.

## Installation

You can install the current development version from github with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/cde")
```

## Basic usage

See the [Get started](https://docs.ropensci.org/cde/articles/cde.html)
vignette or
[Reference](https://docs.ropensci.org/cde/reference/index.html) sections
above for details of the different functions.

## Contributing

For details of how to contribute to this package, see
[here](https://docs.ropensci.org/cde/CONTRIBUTING.html).

<!-- README.md is generated from README.Rmd. Please edit that file -->

# cde <img src="https://docs.ropensci.org/cde/reference/figures/logo.png" align="right" height=140/>

[![R-CMD-check](https://github.com/ropensci/cde/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cde/actions)
[![codecov](https://codecov.io/gh/ropensci/cde/branch/master/graph/badge.svg?token=F4R6nEywTx)](https://codecov.io/gh/ropensci/cde)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://badges.ropensci.org/284_status.svg)](https://github.com/ropensci/onboarding/issues/284)
[![DOI](https://zenodo.org/badge/92712854.svg)](https://zenodo.org/badge/latestdoi/92712854)
[![status](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589/status.svg)](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589)
[![CRAN
Version](http://www.r-pkg.org/badges/version/cde)](http://www.r-pkg.org/pkg/cde)
[![](http://cranlogs.r-pkg.org/badges/cde)](http://cran.rstudio.com/web/packages/cde/index.html)

## Introduction

**NOTE: The EA has recently changed the format of their API, so none of
the download functions in the package are currently working. This is
being worked on and a new release will be made once fixed.**

Within Europe, the [Water Framework
Directive](http://ec.europa.eu/environment/water/water-framework/index_en.html)
(WFD) sets EU-wide standards for how the quality of surface- and
ground-waters across Europe is assessed and classified. Assessment of
quality using the WFD is based on a range of elements that vary
depending on the type of water being assessed and are combined to give
an overall classification of waterbodies into five classes (High, Good,
Moderate, Poor and Bad) for surface waters and two classes (Good or
Poor) for groundwaters.

In the UK the Environment Agency (EA) is the competent authority
responsible for monitoring and assessment of water quality within
England. The EA have made the reporting data relating to the
requirements of the WFD available via the Catchment Data Explorer (CDE)
website, <https://environment.data.gov.uk/catchment-planning/>.

`cde` is a package for R which facilitates searching and download of the
WFD reporting data for all waterbodies from the EA CDE website. The
ability to access these data from within the R environment allows for
efficient collation and interrogation of data and reproducible analysis
of trends or patterns in water quality and pressures on waterbodies
across England. There are also some inconsistencies in the way in which
the data are structured within the original CDE website; `cde` provides
consistently named and structured output which facilitates further
analysis.

The types of data that can be downloaded are: WFD status classification
data, Reasons for Not Achieving Good (RNAG) status, objectives set for
waterbodies, measures put in place to improve water quality and details
of associated protected areas.

The CDE data are made available under the [Open Government Licence
v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/)
and use of the data accessed by and contained within this package
implies acceptance of these licence conditions.

## Installation

You can install the current development version from github with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/cde")
```

## Basic usage

Details of how to use the package can be found at
<https://docs.ropensci.org/cde>.

## Contributing

For details of how to contribute to this package, see
[here](https://docs.ropensci.org/cde/CONTRIBUTING.html).

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# cde (development version)

# cde 0.4.1.9000

# cde 0.4.1

### BUG FIX

  * Fixed bug in year subsetting caused by accidental deletion of endyear code

cde 0.4.0 (2019-05-17)
=========================

### NEW FEATURES

  * Adding `plot` and `print` methods for all output following rOpensci review
  * Adding `cde_df` class for output data
  
### MINOR IMPROVEMENTS

  * Updating all documentation following review
  * Added reference table of all columns returned from CDE site to vignette

cde v 0.4.1
================
Rob Briers
2019-09-03

## Test environments

  - Local Win 7 Enterprise, R 3.6.0 (via R CMD check –as-cran)
  - Local Windows 10, R 3.6.0 (via R CMD check –as-cran)
  - ubuntu 14.04.5, R: release (travis-ci)
  - ubuntu 14.04.5, R: old-rel (travis-ci)
  - ubuntu 14.04.5, R: devel (travis-ci)
  - macOS High Sierra 10.13.3, R: release (travis-ci)
  - macOS High Sierra 10.13.3, R: old-rel (travis-ci)
  - Fedora Linux, R-devel, clang, gfortran (rhub)
  - Ubuntu Linux 16.04 LTS, R-release, GCC (rhub)
  - Windows Server 2008 R2 SP1, R-devel, 32⁄64 bit (rhub)

## R CMD check results

There were no ERRORs or WARNINGs.

There is one NOTE:

CDE (4:53) RNAG (11:64) WFD (3:44, 9:3, 10:49) cde (8:18) rOpenSci
(17:42) waterbodies (9:28, 12:22) Possibly mis-spelled words in
DESCRIPTION:

These are all correct and are mostly abbreviations explained in
supporting docs.

## Downstream dependencies

There aren’t any.

## Resubmission notes

The title of the package has been reduced to less than 65 characters and
the additional LICENCE file and reference to this in the DESCRIPTION
have been removed.
---
title: 'cde - R package to retrieve data from the Environment Agency Catchment Data Explorer site'
tags:
  - R
  - water quality
  - Water Framework Directive
authors:
  - name: Robert A Briers
    orcid: 0000-0003-0341-1203
    affiliation: 1
affiliations:
 - name: School of Applied Sciences, Edinburgh Napier University
   index: 1
date: 17 May 2019
bibliography: paper.bib
---

# Summary

Globally, issues around water quality and quantity are expected to increase in coming decades, set against a background of already widespread degradation [@WWAPUnitedNationsWorldWaterAssessmentProgramme2018]. Within Europe, the [Water Framework Directive](http://ec.europa.eu/environment/water/water-framework/) (WFD) set EU-wide standards for how the quality of surface waters and groundwater across Europe is assessed and classified [@Communities2000]. Assessment of quality under the WFD is based on a range of elements that vary depending on the type of water being assessed. The elements cover biological, chemical and hydromorphological/quantitative components with a hierarchical structure. These are combined to give an overall classification of waterbodies into five classes (High, Good, Moderate, Poor and Bad) for surface waters and two classes (Good or Poor) for groundwaters. The overall aim of the WFD is that all European surface waters and groundwater will reach at least Good status by 2027, although there have been a number of issues with implementation [@Voulvoulis2017].

In the UK, the Environment Agency (EA) is the competent authority responsible for monitoring and assessment of water quality within England. The EA have made all of the reporting data relating to the requirements of the WFD available via the Catchment Data Explorer (CDE) website [https://environment.data.gov.uk/catchment-planning/](https://environment.data.gov.uk/catchment-planning/). The ``cde`` package for R [@RCoreTeam2019] provides functions that facilitate the querying, download and plotting of the data available on the CDE site. The ability to access these data from within the R environment allows for efficient collation and interrogation of data and reproducible analysis of trends or patterns in water quality and pressures on waterbodies across England. There are also some inconsistencies in the way in which the data are structured within the original CDE website; ``cde`` provides consistently named and structured output which facilitates further analysis.

The package allows users to search for waterbodies, catchments or River Basin Districts that match given search strings. Having identified the relevant sites, the following types of data can be downloaded:

*Status classification*: either for overall waterbody classification or at a range of more detailed levels relating to specific quality elements (e.g., ecological, chemical or priority substances). These can be downloaded for specified year ranges and for specific waterbody types (such as lakes).

*Reasons for Not Achieving Good status*: for catchments or River Basin Districts where there are waterbodies that have not achieved Good status, the package provides the functionality to download a summary of the Reasons for Not Achieving Good (RNAG) data. This gives a range of information regarding the relevant pressures identified as contributing to the current status, classified according to a standard hierarchy given on the CDE website.

*Objectives for waterbodies*: where less than Good status has been achieved, data on the objectives that have been set in terms of status aimed for in the longer term can be downloaded, for specific target years and for specified levels of classification.

*Measures to achieve objectives*: details of actions that have been put in place or are proposed to achieve the objectives set (currently in relation to the target objective set for 2021). Only data linked to the achievement of a specific outcome in terms of status are included.

*Protected areas*: a summary of associated protected areas (such as Special Areas of Conservation or Nitrate Vulnerable Zones), again at a range of levels from individual waterbodies to whole River Basin Districts.

For each of the types of data that can be downloaded, summary plots can also be produced. These differ depending on the type of data, but an example showing the percentage of water bodies in each status class (derived from the ``get_status`` function) is shown in Figure 1.

![A plot of status classification data for the Lark Operational Catchment between 2013 and 2015](lark plot-1.png){ width=80% }

# Acknowledgements

Thanks to Matt Starr of the EA and Dave Reynolds of Epimorphics Ltd for useful discussions about the CDE API and providing a full site listing to help development.

# References
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

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

Please note that the `cde` project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>

<!-- README.md is generated from README.Rmd. Please edit that file -->
cde internal data
-----------------

For the purpose of constructing API calls, `cde` makes calls to an dataframe (`ea_wbids`) contained within the `data` folder. This consists of a table of the details (name and index number for each site/catchment) of all waterbodies that can be downloaded from <https://environment.data.gov.uk/catchment-planning/>. The dataframe is documented in `?ea_wbids`. The file here (`ea_wbids.csv`) contains a copy of the internal data content. The data are provided under the terms of the [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).
---
title: "cde: a run through"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{cde: Get Started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Introduction

Within Europe, the [Water Framework Directive](http://ec.europa.eu/environment/water/water-framework/index_en.html) (WFD) sets EU-wide standards for how the quality of surface- and ground-waters across Europe is assessed and classified. Assessment of quality using the WFD is based on a range of elements that vary depending on the type of water being assessed and are combined to give an overall classification of waterbodies into five classes (High, Good, Moderate, Poor and Bad) for surface waters and two classes (Good or Poor) for groundwaters.

In the UK the Environment Agency (EA) is the competent authority responsible for monitoring and assessment of water quality within England. The EA have made the reporting data relating to the requirements of the WFD available via the Catchment Data Explorer (CDE) website, [https://environment.data.gov.uk/catchment-planning/](https://environment.data.gov.uk/catchment-planning/). 

`cde` is a package for R which facilitates searching and download of the WFD reporting data for all waterbodies from the EA CDE website.

The types of data that can be downloaded are: WFD status classification data, Reasons for Not Achieving Good (RNAG) status, objectives set for waterbodies, measures put in place to improve water quality and details of associated protected areas.

The CDE data are made available under the [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) and use of the data accessed by and contained within this package implies acceptance of these licence conditions.

## Installation

You can install the stable version of `cde` from CRAN with:


```r
install.packages("cde")
```

Or you can install the current development version from github with:


```r
# if you have not done so already
# install.packages("remotes")
remotes::install_github("ropensci/cde")
```

## Searching for sites

The `search_sites` function allows you to search for waterbodies, Operational or Management Catchments or River Basin Districts that contain a match or partial match for a specified search string (which is case-sensitive). There is a hierarchical relationship between waterbodies, catchments and River Basin Districts (RBD) as shown [here](https://environment.data.gov.uk/catchment-planning/help#help-catchment-hierarchy). As an example, we will search for waterbodies containing the name "Lark".


```r
# load the package
library(cde)

# search for waterbodies containing the name "Lark"
lark_wb<-search_names(string="Lark", column="name")
```

The dataframe returned contains details of all the waterbodies containing the string "Lark" in their name. The details returned include waterbody id codes (WBID), type of waterbody, Operational and Management Catchment names and River Basin District.


```r
# show the top 6 rows of the 'name' column
head(lark_wb$name)
#> [1] "Lark (US Hawstead)"                   
#> [2] "Lark downstream of Mill Street Bridge"
#> [3] "Lark (Hawstead to Abbey Gardens)"     
#> [4] "Lark (Abbey Gardens to Mildenhall)"   
#> [5] "Lark"                                 
#> [6] "Lark - Fynn (d/s confluence)"
```

To search for Operational Catchments containing the same string we would use the following code.


```r
lark_oc<-search_names(string="Lark", column="OC")
```

## Retrieving quality status classification data

Having located a waterbody, catchment or River Basin District that we want to retrieve data for, we can use the `get_status` function to retrieve the status classification information from the CDE website. We can extract the data for a specific year, or a range of years. For Operational/Management Catchment or River Basin District level downloads, we can also extract information just for a specific waterbody type (such as rivers) or for all waterbody types. In addition it is possible to extract classification data relating to a specific element of the classification.

The overall classification is made up of a number of different elements in a hierarchy. Details of the hierarchy of classification levels can be found [here](https://environment.data.gov.uk/catchment-planning/help#help-classification-hierarchy). By default it retrieves the "Overall Water Body"" classification status, but by specifying the `level`, information on a specific level of classification can be retrieved. The possible values are:

Level 1 | Level 2 | Level 4
--- | --- | ---
Ecological | Biological quality elements | Overall Water Body
Chemical | Chemical Status element | -
Quantitative | Hydromorphological Supporting Elements | -
 - | Other Substances | -
 - | Physico-chemical quality elements | -
 - | Priority hazardous substances | -
 - | Priority substances | -
 - | Quantitative Status element | - 
 - | Specific pollutants | -
 - | Supporting elements | -

The function returns an object of class `cde_df` (basically a dataframe with custom print and plot methods) containing the status (and other details) for the specified combination of column, value, level and dates. Note that during 2013 and 2014 waterbodies were classified under both Cycle 1 and Cycle 2 methodologies. The status information extracted for these years is just for the Cycle 2 classification, to avoid double counting. There was also a change in some of the environmental standards applied to chemical aspects of status assessment between cycles, so there may be some noticeable changes in status between these years. See [here](https://environment.data.gov.uk/catchment-planning/help#help-surface-water-chemical-classification) for more details.

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).


```r
# extract overall waterbody status classification data for a single 
# waterbody in all years

# first decide which waterbody, we can use one from the first search 
# above (need the WBID information)
head(lark_wb)
#>                WBID                                  name  type    OC
#> 1911 GB105033042920                    Lark (US Hawstead) River  Lark
#> 1912 GB105033043052 Lark downstream of Mill Street Bridge River  Lark
#> 1914 GB105033042940      Lark (Hawstead to Abbey Gardens) River  Lark
#> 1918 GB105033043051    Lark (Abbey Gardens to Mildenhall) River  Lark
#> 2197 GB105035040360                                  Lark River Deben
#> 2200 GB105035040300          Lark - Fynn (d/s confluence) River Deben
#>                    MC     RBD
#> 1911 Cam and Ely Ouse Anglian
#> 1912 Cam and Ely Ouse Anglian
#> 1914 Cam and Ely Ouse Anglian
#> 1918 Cam and Ely Ouse Anglian
#> 2197     East Suffolk Anglian
#> 2200     East Suffolk Anglian

# we will get data for the first waterbody here (WBID: GB105033042920, 
# name: Lark (US Hawstead))
lark_hawstead<-get_status(ea_name="GB105033042920", column="WBID")

# the dataframe returned contains all of the data for this site in all 
# years (we did not specify year/year range).
lark_hawstead
#>  river_basin_district management_catchment operational_catchment
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#> With an additional 18 columns of data. 
#> Row values may be truncated to fit console.

# just a quick look at the actual status data
table(lark_hawstead$status)
#> 
#>     Good Moderate 
#>        1        7
```

An example of a higher level download, specifying a year range and type (in this case Rivers).


```r
# download status data for a given year range and type of waterbody
lark_OC_rivers<-get_status(ea_name="Lark", column="OC", startyr=2013, endyr=2015, type="River")
# print out the results
lark_OC_rivers
#>  river_basin_district management_catchment operational_catchment
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#>               Anglian     Cam and Ely Ouse                  Lark
#> With an additional 26 rows and 15 columns of data. 
#> Row values may be truncated to fit console.
```

To get information about status classification in relation to a specific level in the classification, we can specify `level` as well (see table above for options and [here](https://environment.data.gov.uk/catchment-planning/help#help-classification-hierarchy) for more details on the classification levels used).


```r
# download Chemical status for rivers in all years
lark_OC_rivers_chem<-get_status(ea_name="Lark", column="OC", type="River", level="Chemical")
```

## Plotting quality status classification data

The `get_status` function, along with other `get_...` functions, has a `plot` method which provides quick overview plots of status classes, giving a plot of percentages of waterbodies in different status classes for the combination of criteria specified. Plotting is only possible for Operational/Management Catchment or River Basin District downloads.


```r
# get overall waterbody status information for the Lark OC between 2013 and 2015
lark_OC_2013_15 <- get_status(ea_name="Lark", column="OC", startyr=2013, endyr=2015)
# plot the data
plot(lark_OC_2013_15)
```

<img src="figure/lark plot-1.png" title="plot of chunk lark plot" alt="plot of chunk lark plot" style="display: block; margin: auto;" />

For plots, the colour scheme used is based on the `viridis` palette. For `get_status` and `get_objectives` an alternative colour scheme, based on the WFD-defined status class colours, can be used instead by setting `scheme="wfd"` within a `plot` call. Also if a single year is specified, a standard (as opposed to stacked) barplot is produced as shown below.


```r
# get the overall waterbody status information for rivers in the Lark OC in 2015
lark_OC_rivers_2015 <- get_status(ea_name="Lark", column="OC", startyr=2015, type="River")
# plot these data, using WFD colour scheme
plot(lark_OC_rivers_2015, scheme="wfd")
```

<img src="figure/lark riverplot wfd-1.png" title="plot of chunk lark riverplot wfd" alt="plot of chunk lark riverplot wfd" style="display: block; margin: auto;" />

## Reasons for Not Achieving Good status

Not all waterbodies in the Lark Operational Catchment example above have achieved Good status. The `get_rnag` function downloads Reasons for Not Achieving Good (RNAG) data, which allow us to find out more detail on the pressures on the waterbodies that have been assessed to be driving the failure. RNAG data are only available from 2013 onwards. The RNAG data can be extracted for specific years, and also for specific classification levels, as per the status data above.

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).


```r
# what are the RNAG for the Lark OC between 2013 and 2015
lark_OC_RNAG_2013_15<-get_rnag(ea_name="Lark", column="OC", startyr=2013, endyr=2015)
```

Plots of RNAG data are given as frequency histograms of the occurence of information in the `pressure_tier_3` column. For details of this, see the [reference list](../articles/cde-output-reference.html).


```r
# plot RNAG data for the Lark OC, between 2013 and 2015
plot(lark_OC_RNAG_2013_15)
```

<img src="figure/lark RNAG plot-1.png" title="plot of chunk lark RNAG plot" alt="plot of chunk lark RNAG plot" style="display: block; margin: auto;" />

## Objectives set for waterbodies

For those waterbodies that are at less than Good status, objectives are set to indicate what status is aimed for in the longer term. The objectives are set in relation to what is determined to be achievable in the given timescale. Therefore objectives have been set in relation to the 6-year cycle of assessment (so years 2015, 2021 and 2027, then also 2040 and 2050 for long-term objectives). Using the `get_objectives` function, we can download objectives for waterbodies, catchments or River Basin Districts. Objectives can be downloaded for a specific year (2015, 2021, 2027, 2040 or 2050), level of classification and waterbody type as per the `get_status` function. Note however that not all waterbodies have objectives set for all years, levels or types. If no objectives are set for the criteria specified, a message is given.

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).


```r
# download the objectives set for 2015 for the Lark Operational Catchment
lark_OC_obj_2015<-get_objectives(ea_name="Lark", column="OC", year=2015)
```

Plotting of objectives is similar to that of `get_status` data, except the status classes represent the target objectives predicted to be achieved by the date specified.

```r
# plot the objectives for the Lark OC in 2015
plot(lark_OC_obj_2015)
```

<img src="figure/lark obj plot-1.png" title="plot of chunk lark obj plot" alt="plot of chunk lark obj plot" style="display: block; margin: auto;" />

## Protected Areas

The `get_pa` function downloads details of the protected areas associated with a waterbody, catchment or River Basin District. The protected areas listed include those designated under conservation reasons, such as SACs (Habitats and Species Directive), pollution reduction, such as Nitrate Vulnerable Zones (Nitrates Directive) or human use (Bathing Water Directive).

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).


```r
# get details of the protected areas within the Lark Operational Catchment
lark_OC_pa<-get_pa(ea_name="Lark", column="OC")
```

Plotting the output of `get_pa` produces a frequency histogram of the `protected_area_type` column within the area specified.


```r
plot(lark_OC_pa)
```

<img src="figure/lark pa plot-1.png" title="plot of chunk lark pa plot" alt="plot of chunk lark pa plot" style="display: block; margin: auto;" />

## Measures put in place to improve status

Measures are the planned actions that are intended to achieve the objectives set for given waterbodies/catchments etc. The `get_measures` function downloads the details of the measures in place or proposed. These data are very patchy, so this will quite often return no data. Only the name and column need to be specified for this - measures are not specified in more detail than this. Measures can be plotted in the same way as RNAG or Protected Areas data (frequency histogram).

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).


```r
# what measures are there for the Lark Operational Catchment?
lark_OC_meas<-get_measures(ea_name="Lark", column="OC")
```

Plotting the output of `get_measures` produces a frequency histogram of the `measure_category_1` column within the area specified.


```r
plot(lark_OC_meas)
```

<img src="figure/lark measures plot-1.png" title="plot of chunk lark measures plot" alt="plot of chunk lark measures plot" style="display: block; margin: auto;" />
