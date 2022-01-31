
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
---
title: "cde v 0.4.1"
author: "Rob Briers"
date: "`r Sys.Date()`"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Test environments

* Local Win 7 Enterprise, R 3.6.0 (via R CMD check --as-cran)
* Local Windows 10, R 3.6.0 (via R CMD check --as-cran)
* ubuntu 14.04.5, R: release (travis-ci)
* ubuntu 14.04.5, R: old-rel (travis-ci)
* ubuntu 14.04.5, R: devel (travis-ci)
* macOS High Sierra 10.13.3, R: release (travis-ci)
* macOS High Sierra 10.13.3, R: old-rel (travis-ci)
* Fedora Linux, R-devel, clang, gfortran (rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (rhub)
* Windows Server 2008 R2 SP1, R-devel, 32⁄64 bit (rhub)

## R CMD check results

There were no ERRORs or WARNINGs. 

There is one NOTE:

CDE (4:53)
     RNAG (11:64)
     WFD (3:44, 9:3, 10:49)
     cde (8:18)
     rOpenSci (17:42)
     waterbodies (9:28, 12:22)
     Possibly mis-spelled words in DESCRIPTION:

  
These are all correct and are mostly abbreviations explained in supporting docs.

## Downstream dependencies

There aren't any.

## Resubmission notes

The title of the package has been reduced to less than 65 characters 
and the additional LICENCE file and reference to this in the DESCRIPTION
have been removed.
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

# cde <img src="https://docs.ropensci.org/cde/reference/figures/logo.png" align="right" height=140/>
[![R-CMD-check](https://github.com/ropensci/cde/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cde/actions) [![codecov](https://codecov.io/gh/ropensci/cde/branch/master/graph/badge.svg?token=F4R6nEywTx)](https://codecov.io/gh/ropensci/cde) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![](https://badges.ropensci.org/284_status.svg)](https://github.com/ropensci/onboarding/issues/284) [![DOI](https://zenodo.org/badge/92712854.svg)](https://zenodo.org/badge/latestdoi/92712854)
[![status](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589/status.svg)](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589)
[![CRAN Version](http://www.r-pkg.org/badges/version/cde)](http://www.r-pkg.org/pkg/cde)
[![](http://cranlogs.r-pkg.org/badges/cde)](http://cran.rstudio.com/web/packages/cde/index.html)

## Introduction

**NOTE: The EA has recently changed the format of their API, so none of the download functions in the package are currently working. This is being worked on and a new release will be made once fixed.**

Within Europe, the [Water Framework Directive](http://ec.europa.eu/environment/water/water-framework/index_en.html) (WFD) sets EU-wide standards for how the quality of surface- and ground-waters across Europe is assessed and classified. Assessment of quality using the WFD is based on a range of elements that vary depending on the type of water being assessed and are combined to give an overall classification of waterbodies into five classes (High, Good, Moderate, Poor and Bad) for surface waters and two classes (Good or Poor) for groundwaters.

In the UK the Environment Agency (EA) is the competent authority responsible for monitoring and assessment of water quality within England. The EA have made the reporting data relating to the requirements of the WFD available via the Catchment Data Explorer (CDE) website, [https://environment.data.gov.uk/catchment-planning/](https://environment.data.gov.uk/catchment-planning/). 

`cde` is a package for R which facilitates searching and download of the WFD reporting data for all waterbodies from the EA CDE website. The ability to access these data from within the R environment allows for efficient collation and interrogation of data and reproducible analysis of trends or patterns in water quality and pressures on waterbodies across England. There are also some inconsistencies in the way in which the data are structured within the original CDE website; ``cde`` provides consistently named and structured output which facilitates further analysis.

The types of data that can be downloaded are: WFD status classification data, Reasons for Not Achieving Good (RNAG) status, objectives set for waterbodies, measures put in place to improve water quality and details of associated protected areas.

The CDE data are made available under the [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) and use of the data accessed by and contained within this package implies acceptance of these licence conditions.

## Installation

You can install the current development version from github with:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/cde")
```

## Basic usage

Details of how to use the package can be found at [https://docs.ropensci.org/cde](https://docs.ropensci.org/cde).

## Contributing

For details of how to contribute to this package, see [here](https://docs.ropensci.org/cde/CONTRIBUTING.html).

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)


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

# cde <img src="https://docs.ropensci.org/cde/reference/figures/logo.png" align="right" height=140/>
[![R-CMD-check](https://github.com/ropensci/cde/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/cde/actions) [![codecov](https://codecov.io/gh/ropensci/cde/branch/master/graph/badge.svg?token=F4R6nEywTx)](https://codecov.io/gh/ropensci/cde) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![](https://badges.ropensci.org/284_status.svg)](https://github.com/ropensci/onboarding/issues/284) [![DOI](https://zenodo.org/badge/92712854.svg)](https://zenodo.org/badge/latestdoi/92712854)
[![status](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589/status.svg)](http://joss.theoj.org/papers/0d35f75e861fcf47556d70571e226589)
[![CRAN Version](http://www.r-pkg.org/badges/version/cde)](http://www.r-pkg.org/pkg/cde)
[![](http://cranlogs.r-pkg.org/badges/cde)](http://cran.rstudio.com/web/packages/cde/index.html)

## Introduction

**NOTE: The EA has recently changed the format of their API, so none of the download functions in the package are currently working. This is being worked on and a new release will be made once fixed.**

Within Europe, the [Water Framework Directive](http://ec.europa.eu/environment/water/water-framework/index_en.html) (WFD) sets EU-wide standards for how the quality of surface- and ground-waters across Europe is assessed and classified. Assessment of quality using the WFD is based on a range of elements that vary depending on the type of water being assessed and are combined to give an overall classification of waterbodies into five classes (High, Good, Moderate, Poor and Bad) for surface waters and two classes (Good or Poor) for groundwaters.

In the UK the Environment Agency (EA) is the competent authority responsible for monitoring and assessment of water quality within England. The EA have made the reporting data relating to the requirements of the WFD available via the Catchment Data Explorer (CDE) website, [https://environment.data.gov.uk/catchment-planning/](https://environment.data.gov.uk/catchment-planning/). 

`cde` is a package for R which facilitates searching and download of the WFD reporting data for all waterbodies from the EA CDE website. The ability to access these data from within the R environment allows for efficient collation and interrogation of data and reproducible analysis of trends or patterns in water quality and pressures on waterbodies across England. There are also some inconsistencies in the way in which the data are structured within the original CDE website; ``cde`` provides consistently named and structured output which facilitates further analysis.

The types of data that can be downloaded are: WFD status classification data, Reasons for Not Achieving Good (RNAG) status, objectives set for waterbodies, measures put in place to improve water quality and details of associated protected areas.

The CDE data are made available under the [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) and use of the data accessed by and contained within this package implies acceptance of these licence conditions.

## Installation

You can install the current development version from github with:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/cde")
```

## Basic usage

See the [Get started](https://docs.ropensci.org/cde/articles/cde.html) vignette or [Reference](https://docs.ropensci.org/cde/reference/index.html) sections above for details of the different functions.

## Contributing

For details of how to contribute to this package, see [here](https://docs.ropensci.org/cde/CONTRIBUTING.html).
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
## cde internal data

For the purpose of constructing API calls, `cde` makes calls to an dataframe (`ea_wbids`) contained within the `data` folder. This consists of a table of the details (name and index number for each site/catchment) of all waterbodies that can be downloaded from [https://environment.data.gov.uk/catchment-planning/](https://environment.data.gov.uk/catchment-planning/). The dataframe is documented in `?ea_wbids`.  The file here (`ea_wbids.csv`) contains a copy of the internal data content. The data are provided under the terms of the [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).
---
title: "cde: a run through"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{cde: Get Started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

Within Europe, the [Water Framework Directive](http://ec.europa.eu/environment/water/water-framework/index_en.html) (WFD) sets EU-wide standards for how the quality of surface- and ground-waters across Europe is assessed and classified. Assessment of quality using the WFD is based on a range of elements that vary depending on the type of water being assessed and are combined to give an overall classification of waterbodies into five classes (High, Good, Moderate, Poor and Bad) for surface waters and two classes (Good or Poor) for groundwaters.

In the UK the Environment Agency (EA) is the competent authority responsible for monitoring and assessment of water quality within England. The EA have made the reporting data relating to the requirements of the WFD available via the Catchment Data Explorer (CDE) website, [https://environment.data.gov.uk/catchment-planning/](https://environment.data.gov.uk/catchment-planning/). 

`cde` is a package for R which facilitates searching and download of the WFD reporting data for all waterbodies from the EA CDE website.

The types of data that can be downloaded are: WFD status classification data, Reasons for Not Achieving Good (RNAG) status, objectives set for waterbodies, measures put in place to improve water quality and details of associated protected areas.

The CDE data are made available under the [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) and use of the data accessed by and contained within this package implies acceptance of these licence conditions.

## Installation

You can install the current development version from github with:

```{r gh-installation, eval = FALSE}
# if you have not done so already
# install.packages("remotes")
remotes::install_github("ropensci/cde")
```

## Searching for sites

The `search_sites` function allows you to search for waterbodies, Operational or Management Catchments or River Basin Districts that contain a match or partial match for a specified search string (which is case-sensitive). There is a hierarchical relationship between waterbodies, catchments and River Basin Districts (RBD) as shown [here](https://environment.data.gov.uk/catchment-planning/help#help-catchment-hierarchy). As an example, we will search for waterbodies containing the name "Lark".

```{r load and search, eval = TRUE}
# load the package
library(cde)

# search for waterbodies containing the name "Lark"
lark_wb<-search_names(string="Lark", column="name")
```

The dataframe returned contains details of all the waterbodies containing the string "Lark" in their name. The details returned include waterbody id codes (WBID), type of waterbody, Operational and Management Catchment names and River Basin District.

```{r show lark_wb content,eval=TRUE}
# show the top 6 rows of the 'name' column
head(lark_wb$name)
```

To search for Operational Catchments containing the same string we would use the following code.

```{r lark search, eval=TRUE}
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

```{r status data, eval=TRUE}
# extract overall waterbody status classification data for a single 
# waterbody in all years

# first decide which waterbody, we can use one from the first search 
# above (need the WBID information)
head(lark_wb)

# we will get data for the first waterbody here (WBID: GB105033042920, 
# name: Lark (US Hawstead))
lark_hawstead<-get_status(ea_name="GB105033042920", column="WBID")

# the dataframe returned contains all of the data for this site in all 
# years (we did not specify year/year range).
lark_hawstead

# just a quick look at the actual status data
table(lark_hawstead$status)
```

An example of a higher level download, specifying a year range and type (in this case Rivers).

```{r lark river,eval=TRUE}
# download status data for a given year range and type of waterbody
lark_OC_rivers<-get_status(ea_name="Lark", column="OC", startyr=2013, endyr=2015, type="River")
# print out the results
lark_OC_rivers
```

To get information about status classification in relation to a specific level in the classification, we can specify `level` as well (see table above for options and [here](https://environment.data.gov.uk/catchment-planning/help#help-classification-hierarchy) for more details on the classification levels used).

```{r lark rivers chem, eval=TRUE}
# download Chemical status for rivers in all years
lark_OC_rivers_chem<-get_status(ea_name="Lark", column="OC", type="River", level="Chemical")
```

## Plotting quality status classification data

The `get_status` function, along with other `get_...` functions, has a `plot` method which provides quick overview plots of status classes, giving a plot of percentages of waterbodies in different status classes for the combination of criteria specified. Plotting is only possible for Operational/Management Catchment or River Basin District downloads.

```{r lark plot, fig.height=4, fig.width=6.5, fig.align="center", eval=TRUE}
# get overall waterbody status information for the Lark OC between 2013 and 2015
lark_OC_2013_15 <- get_status(ea_name="Lark", column="OC", startyr=2013, endyr=2015)
# plot the data
plot(lark_OC_2013_15)
```

For plots, the colour scheme used is based on the `viridis` palette. For `get_status` and `get_objectives` an alternative colour scheme, based on the WFD-defined status class colours, can be used instead by setting `scheme="wfd"` within a `plot` call. Also if a single year is specified, a standard (as opposed to stacked) barplot is produced as shown below.

```{r lark riverplot wfd,fig.height=4, fig.width=6.5, fig.align="center", eval=TRUE}
# get the overall waterbody status information for rivers in the Lark OC in 2015
lark_OC_rivers_2015 <- get_status(ea_name="Lark", column="OC", startyr=2015, type="River")
# plot these data, using WFD colour scheme
plot(lark_OC_rivers_2015, scheme="wfd")
```

## Reasons for Not Achieving Good status

Not all waterbodies in the Lark Operational Catchment example above have achieved Good status. The `get_rnag` function downloads Reasons for Not Achieving Good (RNAG) data, which allow us to find out more detail on the pressures on the waterbodies that have been assessed to be driving the failure. The RNAG data can be extracted for specific classification levels, as per the status data above.

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).

```{r RNAG in Lark, eval=TRUE}
# what are the RNAG for the Lark OC
lark_OC_RNAG <- get_rnag(ea_name="Lark", column="OC")
```

Plots of RNAG data are given as frequency histograms of the occurence of information in the `pressure_tier_3` column. For details of this, see the [reference list](../articles/cde-output-reference.html).

```{r lark RNAG plot,fig.height=4, fig.width=6.5, fig.align="center", eval=TRUE}
# plot RNAG data for the Lark OC
plot(lark_OC_RNAG)
```

## Objectives set for waterbodies

For those waterbodies that are at less than Good status, objectives are set to indicate what status is aimed for in the longer term. The objectives are set in relation to what is determined to be achievable in the given timescale. Therefore objectives have been set in relation to the 6-year cycle of assessment (so years 2015, 2021 and 2027, then also 2040 and 2050 for long-term objectives). Using the `get_objectives` function, we can download objectives for waterbodies, catchments or River Basin Districts. Objectives can be downloaded for a specific year (2015, 2021, 2027, 2040 or 2050), level of classification and waterbody type as per the `get_status` function. Note however that not all waterbodies have objectives set for all years, levels or types. If no objectives are set for the criteria specified, a message is given.

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).

```{r lark obj, eval=TRUE}
# download the objectives set for 2015 for the Lark Operational Catchment
lark_OC_obj_2015<-get_objectives(ea_name="Lark", column="OC", year=2015)
```

Plotting of objectives is similar to that of `get_status` data, except the status classes represent the target objectives predicted to be achieved by the date specified.
```{r lark obj plot,fig.height=4, fig.width=6.5, fig.align="center", eval=TRUE}
# plot the objectives for the Lark OC in 2015
plot(lark_OC_obj_2015)
```

## Protected Areas

The `get_pa` function downloads details of the protected areas associated with a waterbody, catchment or River Basin District. The protected areas listed include those designated under conservation reasons, such as SACs (Habitats and Species Directive), pollution reduction, such as Nitrate Vulnerable Zones (Nitrates Directive) or human use (Bathing Water Directive).

For details of the meaning of the the different columns returned, see the [output reference list](../articles/cde-output-reference.html).

```{r lark PA, eval=TRUE}
# get details of the protected areas within the Lark Operational Catchment
lark_OC_pa<-get_pa(ea_name="Lark", column="OC")
```

Plotting the output of `get_pa` produces a frequency histogram of the `protected_area_type` column within the area specified.

```{r lark pa plot,fig.height=4, fig.width=6.5, fig.align="center", eval=TRUE}
plot(lark_OC_pa)
```
---
title: "cde: Reference list of all output columns"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{cde-output-reference}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
The table below gives information on all the current columns of data retrieved by the different `get_...` functions of `cde`; names, what they mean and which functions return them. The EA do modify the format of their data from time to time, so this will be updated as necessary.

Column name | Description | Output from which functions
--- | --- | ---
river_basin_district | Name of the River Basin District. This is the top level in the [catchment hierarchy](https://environment.data.gov.uk/catchment-planning/help#help-catchment-hierarchy). | All get_... functions
management_catchment | Name of the Management Catchment, the next level down the hierarchy, see link above. | All get_... functions
operational_catchment | name of the Operational Catchment, the next level down from Management Catchments. | All get_... functions
waterbody_id | Waterbodies are the smallest unit in the catchment hierarchy, representing all or part of a river system, lake, estuary etc. Each has a unique ID (WBID). | All get_... functions
status | WFD status classification (High, Good, Moderate, Poor or Bad). For objectives this represents the target (aimed for) status class. | All get_... functions except get_pa
classification_item | Broad category of classification level (ecological, chemical or quantitative). | get_objectives, get_rnag & get_status
classification_level | The level within the range of classification elements that the status class refers to (see details [here](https://environment.data.gov.uk/catchment-planning/help#help-surface-waters-classification-hierarchy)). | get_objectives, get_rnag & get_status
cycle | Which cycle of assessment the data come from. The first cycle. The first cycle ran from 2009 to 2015 and cycle 2 to 2015-2021. There were changes in some environmental standards between cycles and for a two year period (2013-2014) where both assessment methods were employed, `cde` only returns cycle 2 assessments. | get_objectives, get_rnag & get_status
water_body | Name of the waterbody. | get_objectives, get_rnag & get_status
water_body_type | Type of waterbody. Values are: River, Lake, CoastalWater, TransitionalWater, GroundWaterBody | get_objectives, get_rnag & get_status
year | Year to which the data refer to. | get_objectives, get_rnag & get_status
easting | British National Grid Reference easting, i.e. x coordinate of location. | get_objectives & get_status
hydromorphological_designation | Hydromorphological designation of the waterbody (e.g. whether the body is artificial or heavily modified). | get_objectives & get_status
ngr | Ordnance Survey Landranger grid reference of the location. | get_objectives & get_status
northing | British National Grid Reference northing, i.e. y coordinate of location.  | get_objectives & get_status
activity | Classification of the activity linked to the RNAG identified, for example "Forestry", "Septic Tanks". | get_rnag
activity_certainty | Categorical assessment of the certainty with which the RNAG identified is linked to the activity specifed.  | get_rnag
business_sector | Categorisation of the business sector linked to the activity causing the RNAG, for example "Transport", "Retail sector", "Forestry". | get_rnag
category | Broad categorisation of the RNAG, for example "Industry", "Agriculture and rural land management" or "Water Industry". | get_rnag
category_certainty | Categorical assessment of the certainty with which the RNAG identified is linked to the category specifed. | get_rnag
certainty | Categorical assessment of the certainty of the status class assignment. | get_status
classification_element | Specific quality element or metric linked to the RNAG. | get_rnag
confidence | Confidence in the status class assignment (0-1 scale). | get_status
estimated_start_date | Estimated date at which measure specified will begin.  | get_measures
funding_stream | Source of funding for measure specified. | get_measures
id | Unique numerical code of RNAG (may be linked to measures specified).  | get_rnag
investigation_outcome | Outcome of the investigation into the RNAG. | get_rnag
lead_organisation | Organisation leading the work on the measure specified. | get_measures
linked_rnags | Numerical code(s) of RNAGS linked to measures specified (see 'id' above). | get_measures
measure_category_1 | Top level classification of measure specified. | get_measures
measure_category_2 | Next level down (more detailed) classification of measure specified. | get_measures
measure_category_3 | Most detailed classification of measure specified. | get_measures
measure_reference_code | Numerical code of measure type. | get_measures
measure_type | Classification of type of measure employed. | get_measures
national_swmi_header | SWMI (Significant Water Management Issue) classification at the national level, for example "Pollution from waste water" or "Physical modifications". | get_rnag
objectivetype | Type of objective set, either "Predicted" or "Objective" | get_objectives
pressure_tier_1 | Detailed classification of the pressure related to the RNAG identified, for example "Irrigation", "Priority Substances (pesticides)". | get_rnag
pressure_tier_2 | Specific substance or pressure associated with RNAG identiifed, for example "Phosphate", "High temperature". | get_rnag
pressure_tier_3 | High level classification of the pressure related to the RNAG identifed, for example "Physical modification", "Dissolved oxygen (DO)". | get_rnag
protected_area_code | Unique code of each protected area. | get_pa
protected_area_label | Name or code of protected area. | get_pa
protected_area_type | Designation of protected area, for example "Nitrates Directive" ,  "Bathing Water Directive". | get_pa
reason_type | Type of RNAG, either "RFF" - reason for failure or "RFD" - reason for deterioration. | get_rnag
reasons_for_alternative_objectives | Detailed text as to why alternative objectives i.e. a status less than good, have been set. | get_objectives
rnag_url |  URL to page on CDE site giving specifics of RNAG. | get_rnag
sector_of_lead_organisation | Sectoral classification of organisation leading measures work. | get_measures
see_also | Additional details about the protected area specified. | get_pa
swmi | SWMI (Significant Water Management Issue) classification of the RNAG, for example "Diffuse source", "Invasive non-native species". | get_rnag
swmi_certainty | Categorical assessment of the certainy with which the RNAG is linked to the SWMI classification specified. | get_rnag
title | Title of the measures project specified. | get_measures

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cde-class_and_methods.r
\name{print.cde_df}
\alias{print.cde_df}
\title{Print method for cde_df}
\usage{
\method{print}{cde_df}(x, ...)
}
\arguments{
\item{x}{An object of class \code{cde_df}.}

\item{...}{Other arguments passed on to individual methods. None 
implemented at present.}
}
\description{
Custom \code{print} method for objects of class \code{cde_df}.
Formats output to fit current width of console, keeping full column names 
but truncating row values as required.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_objectives.r
\name{get_objectives}
\alias{get_objectives}
\title{Retrieve Objectives set for waterbodies}
\usage{
get_objectives(
  ea_name = NULL,
  column = NULL,
  level = "Overall Water Body",
  year = NULL,
  type = NULL
)
}
\arguments{
\item{ea_name}{A string representing the description (\code{name} for 
\code{OC}, \code{MC} or \code{RBD} level downloads or \code{WBID} for 
individual waterbodies) of the features to be extracted. For example 
to extract data for the whole of the Humber RBD, this would be "Humber"; 
also see examples. Must be an exact match to the values used in the EA 
database. Use the \code{\link{search_names}} function to search for 
specific values.}

\item{column}{The column to be searched. Possible options are
\code{WBID} (waterbody id), \code{OC} (Operational Catchment), \code{MC}
(Management Catchment) and \code{RBD} (River Basin District)}

\item{level}{The level within the WFD quality classification elements that 
objectives have been set at. For full details of the hierarchy of elements 
within the classification used, see \url{https://environment.data.gov.uk/catchment-planning/help#help-classification-hierarchy}.

Defaults to 'Overall Water Body'. Possible values for the different levels 
retrieved by the function are shown below.

\tabular{ccc}{
 \strong{Level 1} \tab \strong{Level 2} \tab \strong{Level 4}\cr
Ecological \tab Biological quality elements \tab Overall Water Body \cr
Chemical \tab Chemical Status element \tab -\cr
Quantitative \tab Hydromorphological Supporting Elements \tab -\cr
 - \tab Other Substances \tab -\cr
 - \tab Physico-chemical quality elements \tab -\cr
 - \tab Priority hazardous substances \tab -\cr
 - \tab Priority substances \tab -\cr
 - \tab Quantitative Status element \tab - \cr
 - \tab Specific pollutants \tab -\cr
 - \tab Supporting elements \tab -\cr
}}

\item{year}{The year that objectives are set for, either 2015, 
2021, 2027, 2040 or 2050. If not given then objectives for all years 
are returned. Note that objectives may not be set for all years.}

\item{type}{Type of waterbody to be extracted. For Operational/Management
catchment level or RBD level queries, the data can also be subset by
waterbody type. Possible values are \code{River}, \code{Lake},
\code{GroundWaterBody}, \code{TransitionalWater} or \code{CoastalWater}.}
}
\value{
An object of class \code{cde_df} containing the details 
of the objectives set for the specified set of waterbodies.
For details of the meaning of the the different columns returned, 
see \url{https://docs.ropensci.org/cde/articles/cde-output-reference.html}.
}
\description{
Retrieves details of objectives set for waterbodies in terms
of predicted classification from EA Catchment Data Explorer site.
Data can be retrieved by specifying waterbody id
(\code{WBID}), Management Catchment (\code{MC}), Operational
Catchment (\code{OC}) or River Basin District (\code{RBD}).
Start year (\code{startyr}) and end year (\code{endyr}) allow
specific timeranges to be downloaded.
For Management Catchment (\code{MC}), Operational
Catchment (\code{OC}) or River Basin District (\code{RBD}) level
downloads, waterbody \code{type} can also be specified to allow
extraction of specific waterbody types (River, Lake etc).
}
\examples{
# get all objectives set for waterbody GB112071065700
get_objectives(ea_name="GB112071065700", column="WBID")

# get the objectives set for Lakes in the Humber RBD, for the year 2021
get_objectives(ea_name="Humber", column="RBD", year=2021, type="Lake")

# get the objectives set for Rivers in the Avon Warwickshire
# Management Catchment in relation to Chemical status
get_objectives(ea_name="Avon Warwickshire", column="MC", level="Chemical", type="River")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cde-class_and_methods.r
\name{plot.cde_df}
\alias{plot.cde_df}
\title{Plot method for \code{cde_df} output}
\usage{
\method{plot}{cde_df}(x, ...)
}
\arguments{
\item{x}{An object of class \code{cde_df} to be plotted.}

\item{...}{Other arguments passed on to individual methods. The only other
argument implemented at present is \code{scheme}. For \code{status} and 
\code{objectives} data this defines which colour scheme to use with plots. It 
defaults to a viridis-based scheme (\code{scheme="vir"}). Alternatively, the 
colours specified in the WFD document can be used by specifying 
\code{scheme="wfd"}.}
}
\description{
Default plots of the output main \code{get_} functions.
Details of the plots for different data are given below.
 
 For \code{status} and \code{objectives} produces a (stacked) 
 percentage barplot of waterbody observed or predicted (objective) 
 status information for a given set of data. 

For \code{rnag}, \code{measures} or \code{pa} produces a frequency 
histogram. The columns plotted for each data type are given below:

\itemize{
  \item{rnag} 
  (pressure_tier_3)
  \item {measures} 
  (measure_category_1)
  \item {pa} 
  (protected_area_type)
}
The full detail of the different data being plotted can be found in 
the EA Catchment Data Explorer API reference:
\url{https://environment.data.gov.uk/catchment-planning/ui/reference}
     
Plotting is only possible for MC, OC or RBD downloads.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_pa.r
\name{get_pa}
\alias{get_pa}
\title{Retrieve Protected Area Information}
\usage{
get_pa(ea_name = NULL, column = NULL)
}
\arguments{
\item{ea_name}{A string representing the description (\code{name} for 
\code{OC}, \code{MC} or \code{RBD} level downloads or \code{WBID} for 
individual waterbodies) of the features to be extracted. For example 
to extract data for the whole of the Humber RBD, this would be "Humber"; 
also see examples. Must be an exact match to the values used in the 
EA database. Use the \code{\link{search_names}} function to search 
for specific values.}

\item{column}{The column to be searched. Possible options are
\code{WBID} (waterbody id), \code{OC} (Operational Catchment), \code{MC}
(Management Catchment) and \code{RBD} (River Basin District)}
}
\value{
An object of class \code{cde_df} containing the details of the 
Protected Areas associated with the waterbodies.
For details of the meaning of the the different columns returned, 
see \url{https://docs.ropensci.org/cde/articles/cde-output-reference.html}.
}
\description{
Retrieves details of Protected Areas associated with 
waterbodies, catchments or River Basin Districts from the EA 
Catchment Data Explorer site.
Data can be retrieved by specifying waterbody id
(\code{WBID}), Management Catchment (\code{MC}), Operational
Catchment (\code{OC}) or River Basin District (\code{RBD}).
}
\examples{
# get protected areas associated with waterbody GB112071065700
get_pa(ea_name="GB112071065700", column="WBID")

# get the protected areas associated with the Humber RBD
get_pa(ea_name="Humber", column="RBD")

# get the protected areas associated with the Avon Warwickshire
# Management Catchment
get_pa(ea_name="Avon Warwickshire", column="MC")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_status.r
\name{get_status}
\alias{get_status}
\title{Retrieve WFD Status Classification Data}
\usage{
get_status(
  ea_name = NULL,
  column = NULL,
  level = "Overall Water Body",
  startyr = NULL,
  endyr = NULL,
  type = NULL
)
}
\arguments{
\item{ea_name}{A string representing the description (\code{name} for 
\code{OC}, \code{MC} or \code{RBD} level downloads or \code{WBID} for 
individual waterbodies) of the features to be extracted. For example 
to extract data for the whole of the Humber RBD, this would be "Humber"; 
also see examples. Must be an exact match to the values used in the EA 
database. Use the \code{\link{search_names}} function to search for 
specific values.}

\item{column}{The column to be searched. Possible options are
\code{WBID} (waterbody id), \code{OC} (Operational Catchment), \code{MC}
(Management Catchment) and \code{RBD} (River Basin District)}

\item{level}{The level within the WFD quality classification elements that 
objectives have been set at. For full details of the hierarchy of elements 
within the classification used, see \url{https://environment.data.gov.uk/catchment-planning/help/usage#catchment-hierarchy}.

Defaults to 'Overall Water Body'. Possible values for the different levels 
retrieved by the function are shown below.
\tabular{ccc}{
\strong{Level 1} \tab \strong{Level 2} \tab \strong{Level 4}\cr
Ecological \tab Biological quality elements \tab Overall Water Body\cr
Chemical \tab Chemical Status element \tab -\cr
Quantitative \tab Hydromorphological Supporting Elements \tab -\cr
- \tab Other Substances \tab -\cr
- \tab Physico-chemical quality elements \tab -\cr
- \tab Priority hazardous substances \tab -\cr
- \tab Priority substances \tab -\cr
- \tab Quantitative Status element \tab - \cr
- \tab Specific pollutants \tab -\cr
- \tab Supporting elements \tab -\cr
}}

\item{startyr}{The data can be extracted for specific years using the
\code{startyr} and \code{endyr} arguments. If only \code{startyr} is
specified this extracts for a particular year. If no years are specified
all years are returned.}

\item{endyr}{The data can be extracted for specific years using the
\code{startyr} and \code{endyr} arguments. The \code{endyr} should
only be specified if \code{startyr} is also included, otherwise an 
error is returned.}

\item{type}{Type of waterbody to be extracted. For Operational/Management
catchment level or RBD level queries, the data can also be subset by
waterbody type. Possible values are \code{River}, \code{Lake},
\code{GroundWaterBody}, \code{TransitionalWater} or \code{CoastalWater}.}
}
\value{
An object of class \code{cde_df} containing the classification 
details for the specified combination of criteria.
For details of the meaning of the the different columns returned, 
see \url{https://docs.ropensci.org/cde/articles/cde-output-reference.html}.
}
\description{
Retrieves WFD Status classification data from EA Catchment 
Data Explorer site. Data can be retrieved by specifying waterbody id
(\code{WBID}), Management Catchment (\code{MC}), Operational
Catchment (\code{OC}) or River Basin District (\code{RBD}).
Start year (\code{startyr}) and end year (\code{endyr}) allow
specific timeranges to be downloaded.
For Management Catchment (\code{MC}), Operational
Catchment (\code{OC}) or River Basin District (\code{RBD}) level
downloads, waterbody \code{type} can also be specified to allow
extraction of specific waterbody types (River, Lake etc).
}
\examples{
# get Overall Water Body status classification for waterbody GB520804714300
get_status(ea_name="GB520804714300", column="WBID")

# get status class based on Priority substances for waterbody GB520804714300
get_status(ea_name="GB520804714300", column="WBID", level="Priority substances")

# get the Overall Water Body status of Lakes in the Humber RBD, between
# 2012 and 2014
get_status(ea_name="Humber", column="RBD", startyr=2012, endyr=2014, type="Lake")

# get the Overall Water Body status for Rivers in the Avon Warwickshire
# Operational Catchment in 2011
get_status(ea_name="Avon Warwickshire", column="MC", startyr=2011, type="River")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cde.r
\docType{package}
\name{cde}
\alias{cde}
\title{cde: Download Water Framework Directive (WFD) data from the 
Environment Agency Catchment Data Explorer (CDE) website.}
\description{
Facilitates searching and download of the WFD-related
data for all waterbodies within the Environment Agency area (i.e. England).
The types of data that can be downloaded are: WFD status classification 
data, Reasons for Not Achieving Good (RNAG) status, objectives set for 
waterbodies, and details of associated protected areas. Default plots can 
also be produced from thedata downloaded (form of plot depends on data type).
}
\details{
The website that is accessed is: 
\url{https://environment.data.gov.uk/catchment-planning/}.
The data accessed by and included within the package are made available 
under the Open Government Licence v3.0 
\url{https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_names.r
\name{search_names}
\alias{search_names}
\title{Search database of site names}
\usage{
search_names(string = NULL, column = NULL)
}
\arguments{
\item{string}{The search string to be matched (case-sensitive). Will match 
whole or partial strings in the column values.}

\item{column}{The column to be searched. Possible options are
\code{WBID}, \code{name}, \code{OC} (Operational Catchment), \code{MC}
(Management Catchment) and \code{RBD} (River Basin District)}
}
\value{
A data frame containing the details of all the sites that match
the search string (full or partial matches) in the column specified.
Columns returned are defined in \code{\link{ea_wbids}}.
}
\description{
Searches the listing of EA monitoring sites to find rows
that contain the string provided. Can search by WBID (\code{WBID}), name 
(\code{name}), Management Catchment (\code{MC}), Operational Catchment 
(\code{OC}) or River Basin District (\code{RBD}). There is a hierarchical 
relationship between these levels as shown at \url{https://environment.data.gov.uk/catchment-planning/help#help-catchment-hierarchy}.

The search is done on a local copy of the waterbody listing contained in 
the \code{\link{ea_wbids}} object rather than connecting to the 
EA site.
}
\examples{
# search for sites containing "Tadnoll" in the name
search_names(string="Tadnoll", column="name")

# search for Operational Catchments containing "Cornwall"
search_names(string="Cornwall", column="OC")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ea_wbids.r
\docType{data}
\name{ea_wbids}
\alias{ea_wbids}
\title{Details of name and index of all sites/catchments.}
\format{
A data frame with 5237 rows and 9 variables:
\describe{
  \item{WBID}{identifier for individual waterbodies}
  \item{name}{detailed name of the site/catchment}
  \item{type}{type of waterbody (River, Lake etc.)}
  \item{OC}{Operational Catchment name}
  \item{OC_num}{Index number of the Operational Catchment}
  \item{MC}{Management Catchment name}
  \item{MC_num}{Index number of the Management Catchment}
  \item{RBD}{River Basin District name}
  \item{RBD_num}{Index number of the River Basin District}
}
For details of the hierarchy of the different catchment types, see
\url{https://environment.data.gov.uk/catchment-planning/help#help-catchment-hierarchy}
}
\source{
\url{https://environment.data.gov.uk/catchment-planning/}
}
\usage{
ea_wbids
}
\description{
Dataframe used by `cde` to construct API calls.
The data included are made available under the Open Government Licence v3.0 
\url{https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/}.
Use of the data accessed by and contained within this package implies 
acceptance of these licence conditions.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_rnag.r
\name{get_rnag}
\alias{get_rnag}
\title{Retrieve Reasons for Not Achieving Good Status}
\usage{
get_rnag(ea_name = NULL, column = NULL, type = NULL)
}
\arguments{
\item{ea_name}{A string representing the description (\code{name} for 
\code{OC}, \code{MC} or \code{RBD} level downloads or \code{WBID} for 
individual waterbodies) of the features to be extracted. For example 
to extract data for the whole of the Humber RBD, this would be "Humber"; 
also see examples. Must be an exact match to the values used in the EA 
database. Use the \code{\link{search_names}} function to search for 
specific values.}

\item{column}{The column to be searched. Possible options are
\code{WBID} (waterbody id), \code{OC} (Operational Catchment), \code{MC}
(Management Catchment) and \code{RBD} (River Basin District)}

\item{type}{Type of waterbody to be extracted. For Operational/Management
catchment level or RBD level queries, the data can also be subset by
waterbody type. Possible values are \code{River}, \code{Lake},
\code{GroundWaterBody}, \code{TransitionalWater} or \code{CoastalWater}.}
}
\value{
An object of class \code{cde_df} containing the details of the 
Reasons for Not Achieving Good Status for the specified combination 
of criteria.
For details of the meaning of the the different columns returned, 
see \url{https://docs.ropensci.org/cde/articles/cde-output-reference.html}.
}
\description{
Retrieves details of Reasons for Not Achieving Good (RNAG)
status and Reasons For Failure (RFF) from EA Catchment Data Explorer site.
Data can be retrieved by specifying waterbody id
(\code{WBID}), Management Catchment (\code{MC}), Operational
Catchment (\code{OC}) or River Basin District (\code{RBD}).
Start year (\code{startyr}) and end year (\code{endyr}) allow
specific timeranges to be downloaded.
For Management Catchment (\code{MC}), Operational
Catchment (\code{OC}) or River Basin District (\code{RBD}) level
downloads, waterbody \code{type} can also be specified to allow
extraction of specific waterbody types (River, Lake etc).
Data are presented at the level of individual elements that are the
reasons for not achieving good status.
}
\examples{
# get all RNAG issues identified for waterbody GB112071065700
get_rnag("GB112071065700", "WBID")

# get the RNAG issues for Lakes in the Humber RBD, between
# 2013 and 2014
get_rnag(ea_name="Humber", column="RBD", type="Lake")

# get the RNAG issues for Rivers in the Avon Warwickshire
# Management Catchment
get_rnag(ea_name="Avon Warwickshire", column="MC", type="River")

}
