
# tbeptools

[![R-CMD-check](https://github.com/tbep-tech/tbeptools/workflows/R-CMD-check/badge.svg)](https://github.com/tbep-tech/tbeptools/actions)
[![pkgdown](https://github.com/tbep-tech/tbeptools/workflows/pkgdown/badge.svg)](https://github.com/tbep-tech/tbeptools/actions)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03485/status.svg)](https://doi.org/10.21105/joss.03485)
[![Codecov test coverage](https://codecov.io/gh/tbep-tech/tbeptools/branch/master/graph/badge.svg)](https://codecov.io/gh/tbep-tech/tbeptools?branch=master)
[![DOI](https://zenodo.org/badge/184627857.svg)](https://zenodo.org/badge/latestdoi/184627857)

R package for Tampa Bay Estuary Program functions. Please see the [vignette](https://tbep-tech.github.io/tbeptools/articles/intro.html) for a full description.

<img src="man/figures/logo.png" align="center" width="125"/>

Please cite this package as follows: 

Beck, M.W., Schrandt, M.N., Wessel, M.R., Sherwood, E.T., Raulerson, G.E., Budihal Prasad, A.A., Best, B.D., (2021). tbeptools: An R package for synthesizing estuarine data for environmental research. Journal of Open Source Software, 6(65), 3485, https://doi.org/10.21105/joss.03485

# Installation

The package can be installed from [r-universe](https://tbep-tech.r-universe.dev).  The source code is available on the tbep-tech GitHub group web page: <https://github.com/tbep-tech/tbeptools>.  Note that tbeptools only needs to be installed once, but it needs to be loaded every new R session (i.e., `library(tbeptools)`).

```r
# enable repos
options(repos = c(
    tbeptech = 'https://tbep-tech.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

# install tbeptools
install.packages('tbeptools')

# load tbeptools
library(tbeptools)
```

After the package is loaded, you can view the help files for each function by typing a question mark followed by the function name, e.g., `?read_importwq`, on the console.  The help files provide a brief description of what each function does and the required arguments that are needed to run the function.

# Package vignettes

The vignettes are organized by topic and are an excellent place to start for understanding how to use the package. Currently, there are five vignettes available for tbeptools:

* [Intro to TBEP tools](https://tbep-tech.github.io/tbeptools/articles/intro.html): A general overview of the package with specific examples of functions for working with the water quality report card
* [Tampa Bay Nekton Index](https://tbep-tech.github.io/tbeptools/articles/tbni.html): Overview of functions to import, analyze, and plot results for the Tampa Bay Nekton Index
* [Tampa Bay Benthic Index](https://tbep-tech.github.io/tbeptools/articles/tbbi.html): Overview of functions to import data for Tampa Bay Benthic Index, under development
* [Tidal Creeks Assessment](https://tbep-tech.github.io/tbeptools/articles/tidalcreeks.html): Overview of functions to import, analyze, and plot results for the assessment of tidal creeks in southwest Florida
* [Seagrass Transect Data](https://tbep-tech.github.io/tbeptools/articles/seagrasstransect.html): Overview of functions to import, analyze, and plot results for the seagrass transect data collected in Tampa Bay

# Usage

Functions in tbeptools fall in three categories depending on mode of use.  Each function is named using a prefix for the mode of use, followed by what the function does. The prefixes are:

* `read`: Import current data from the main ftp site.

* `anlz`: Analyze or summarize the imported data. 

* `show`: Create a plot of the analyzed data.

The functions can be easily found in RStudio after loading the package and typing the prefix at the command line.  An autofill dialog box will pop up showing all functions that apply for the prefix. This eliminates the need for searching for individual functions if all you know is the category of function you need (e.g., `read`, `anlz`, or `show`).

Each function also includes a semi-descriptive suffix that generally describes what category it applies to (e.g, water quality, seagrass) and what it does (e.g., imports, formats).  These follow a loose convention that attempts to strike a balance between description and brevity.  The optimal balance is often hard to achieve.  To aid in understanding, we provide a brief description of suffixes that are used more than once.  

Suffix descriptions:

* `attain`: Analyze functions that summarize data relative to attainment categories specific to bay segments
* `ave`, `med`: Analyze functions that summarize data into averages or medians
* `benthic`: Applies to benthic monitoring data used for the Tampa Bay Benthic Index
* `fim`: Applies to data from the Fisheries Independent Monitoring program used for the Tampa Bay Nekton Index
* `form`: An intermediate function for formatting imported data for downstream analysis
* `import`: A function used to import data from a source external to the package
* `indic`: A function that analyzes or plots individual tidal creek indicator values, as opposed to integrated creek scores
* `iwr`: Functions or data that apply to the Impaired Waters Rule (IWR) data maintained by the Florida Department of Environmental Protection used as source data for the tidal creek functions
* `matrix`: A plotting function that creates a report card style matrix
* `met`: A function that analyses or plots individual metrics for integrated indices, e.g., TBBI, TBNI
* `phyto`: Applies to phytoplankton data from the Hillsborough County Environmental Protection Commission
* `plotly`: A plotting function that returns an interactive plotly object
* `scr`: A function that analyses or plots summary scores for integrated indices, e.g., TBBI, TBNI
* `seg`, `site`: Functions that analyze or plot results relative to bay segments or individual monitoring sites 
* `tbbi`: Applies to the Tampa Bay Benthic Index (TBBI)
* `tbni`: Applies to the Tampa Bay Nekton Index (TBNI)
* `tdlcrk`: Applies to tidal creeks
* `transect`: Applies to seagrass transect data
* `wq`: Applies to water quality

The function [reference page](https://tbep-tech.github.io/tbeptools/reference/index.html) can also be viewed for a complete list of functions organized by category, a description of what they do, and links to the help files. 

The following example demonstrates use of a subset of the functions for water quality data to read a file from the Hillsborough County Environmental Protection Commission long-term monitoring dataset (available from <https://www.tampabay.wateratlas.usf.edu/>), analyze monthly and annual averages by major bay segments of Tampa Bay, and plot an annual time series for one of the bay segments.

```r
# load the package
library(tbeptools)

# read current data
wqdat <- read_importwq(xlsx = "wqdata.xlsx", download_latest = TRUE)
wqdat
```

```
## # A tibble: 26,611 x 22
##   bay_segment epchc_station SampleTime             yr    mo
##   <chr>               <dbl> <dttm>              <dbl> <dbl>
## 1 HB                      6 2021-06-08 10:59:00  2021     6
## 2 HB                      7 2021-06-08 11:13:00  2021     6
## 3 HB                      8 2021-06-08 14:15:00  2021     6
## 4 MTB                     9 2021-06-08 13:14:00  2021     6
## 5 MTB                    11 2021-06-08 11:30:00  2021     6
## # ... with 26,606 more rows, and 17 more variables:
## #   Latitude <dbl>, Longitude <dbl>, Total_Depth_m <dbl>,
## #   Sample_Depth_m <dbl>, tn <dbl>, tn_q <chr>, sd_m <dbl>,
## #   sd_raw_m <dbl>, sd_q <chr>, chla <dbl>, chla_q <chr>,
## #   Sal_Top_ppth <dbl>, Sal_Mid_ppth <dbl>,
## #   Sal_Bottom_ppth <dbl>, Temp_Water_Top_degC <dbl>,
## #   Temp_Water_Mid_degC <dbl>, ...
```

```r
# analyze monthly and annual means by bay segment
avedat <- anlz_avedat(wqdat)
avedat
```

```
## $ann
## # A tibble: 584 x 4
##      yr bay_segment var         val
##   <dbl> <chr>       <chr>     <dbl>
## 1  1974 HB          mean_chla 22.4 
## 2  1974 LTB         mean_chla  4.24
## 3  1974 MTB         mean_chla  9.66
## 4  1974 OTB         mean_chla 10.2 
## 5  1975 HB          mean_chla 27.9 
## # ... with 579 more rows
## 
## $mos
## # A tibble: 4,484 x 5
##   bay_segment    yr    mo var         val
##   <chr>       <dbl> <dbl> <chr>     <dbl>
## 1 HB           1974     1 mean_chla 36.2 
## 2 LTB          1974     1 mean_chla  1.75
## 3 MTB          1974     1 mean_chla 11.5 
## 4 OTB          1974     1 mean_chla  4.4 
## 5 HB           1974     2 mean_chla 42.4 
## # ... with 4,479 more rows
```

```r
# show annual time series of chlorophyll for Hillsborough bay segment
show_thrplot(wqdat, bay_segment = "HB", yrrng = c(1975, 2020))
```

<img src="man/figures/thrplotex-1.jpeg" align="center"/>

Functions in `tbeptools` also support the creation of content for interactive, online dashboards that can facilitate more informed decisions without requiring an intimate understanding of the R programming language or the methods for analysis.  These dashboards include assessments for [water quality](https://shiny.tbeptech.org/wq-dash/), [seagrasses](https://shiny.tbep.org/seagrasstransect-dash/), [nekton communities](http://shiny.tbeptech.org/nekton-dash), and [tidal creeks](https://shiny.tbep.org/tidalcreek-dash/).

# Issues and suggestions

Please report any issues and suggestions on the [issues link](https://github.com/tbep-tech/tbeptools/issues) for the repository.  A guide to posting issues can be found [here](.github/ISSUE_TEMPLATE.md).

# Contributing

Please view our [contributing](.github/CONTRIBUTING.md) guidelines for any changes or pull requests.
---
title: 'tbeptools: An R package for synthesizing estuarine data for environmental research'
tags:
  - R
  - estuary
  - Tampa Bay
  - water quality
  - reporting
authors:
  - name: Marcus W. Beck^[Corresponding author]
    orcid: 0000-0002-4996-0059
    affiliation: 1 
  - name: Meagan N. Schrandt
    orcid: 0000-0002-0482-5072
    affiliation: 2
  - name: Michael R. Wessel
    affiliation: 3
  - name: Edward T. Sherwood
    orcid: 0000-0001-5330-302X
    affiliation: 1
  - name: Gary E. Raulerson
    orcid: 0000-0002-5920-5743
    affiliation: 1  
  - name: Adhokshaja Achar Budihal Prasad
    orcid: 0000-0002-6485-6858
    affiliation: 4
  - name: Benjamin D. Best
    orcid: 0000-0002-2686-0784
    affiliation: 5
affiliations:
  - name: Tampa Bay Estuary Program, St. Petersburg, Florida, USA
    index: 1
  - name: Fish and Wildlife Research Institute, Florida Fish and Wildlife Conservation Commission, St. Petersburg, Florida, USA
    index: 2
  - name: Janicki Environmental, Inc., St. Petersburg, Florida, USA
    index: 3
  - name: University of South Florida, Tampa, Florida, USA
    index: 4
  - name: EcoQuants, LLC, Santa Barbara, California, USA
    index: 5
date: 27 April 2021
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---



# Summary

Many environmental programs report on the status and trends of natural resources to inform management decisions for protecting or restoring environmental condition.  The National Estuary Program (NEP) in the United States is one example of a resource management institution focused on "estuaries of national significance" that provides place-based solutions to managing coastal resources.  There are 28 NEPs in the United States, each with similar but location-specific programmatic goals to address environmental challenges related to water quality, alteration of hydrologic flows, invasive species, climate change,  declines in fish and wildlife populations, pathogens and other contaminants, and stormwater management.  A critical need of each NEP is the synthesis of data from disparate sources that can inform management response to address these environmental challenges. 

The Tampa Bay Estuary Program (TBEP) in Florida, USA is responsible for developing and implementing a place-based plan to sustain historical and future progress in the restoration of Tampa Bay [@tbep1017].  The needs of TBEP for reporting on indicators of environmental condition are similar to other environmental organizations.  Multiple local and regional partners collect data that are used for different reporting products.  Without data synthesis tools that are transparent, accessible, and reproducible, NEP staff and colleagues waste time and resources compiling information by hand.  The `tbeptools` R software package can be used for routine development of reporting products, allowing for more efficient use of limited resources and a more effective approach to communicate research to environmental decision-makers. Functions in `tbeptools` also support the creation of content for interactive, online dashboards that can facilitate more informed decisions without requiring an intimate understanding of the R programming language or the methods for analysis.       

The `tbeptools` package also addresses challenges associated with data retrieval and assessment of environmental data relative to important policy or management targets.  For each environmental indicator, functions are included to import required data directly from sources, removing the need to manually obtain information prior to reporting.  These functions are integrated into summary report cards that are generated automatically through continuous integration services (i.e., [GitHub Actions](https://github.com/features/actions)) that free the analyst from external downloads, analysis, and copying of results that can introduce errors in reporting.  Similar packages provide seamless access to database services [e.g., `dataRetrieval`, @DeCicco21], but few packages link these data sources directly to analysis and reporting as in `tbeptools` [but see `wqindex`, @Thorley18].  Management targets and regulatory thresholds based on summary assessments of data, either direct from sources or included as supplementary data in the package, are also hard-coded into the functions. 

# Statement of need

The `tbeptools` R package was developed to automate data synthesis and analysis for many of the environmental indicators for Tampa Bay, with more general application to commonly available datasets for estuaries.  The functions in the package were developed to extract methods from existing technical documents and to make them available in an open source programming environment.  By making these tools available as an R package, routine assessments are now accomplished more quickly and other researchers can use the tools to develop more specific analysis pipelines.  

Most of the NEPs do not have analysis software to operationalize data import, analysis, and plotting for reporting.  Recently, a similar software package, `peptools` [@Beck21a], was developed for the Peconic Estuary Partnership (New York, USA) using many of the functions in `tbeptools` to develop reporting products for a new water quality monitoring program.  This successful technology transfer demonstrates the added value of presenting these methods in an open source environment available for discovery and reuse by others.  We expect other NEPs to begin using these tools as their application becomes more widespread among estuarine researchers.

Beyond the NEPs, `tbeptools` is an effective example of an R package for implementing technical methods in existing literature and reports that can be used to support environmental monitoring and assessment needs for science-based decisions.  To this end, the `tbeptools` package was also created to support the development of online dashboards created in R Shiny [@Chang21].  Dashboards are powerful tools to increase accessibility for end users to engage with scientific products without the need to understand technical details in their creation.  However, providing the underlying methods as source code in an R package increases transparency and reproducibility of reporting products if users require a more detailed understanding of how the content was created.  Currently, the `tbeptools` package supports dashboards created by TBEP for the assessment of [water quality](https://shiny.tbeptech.org/wq-dash/) [Figure \ref{fig:dashex}, @Beck20a], [seagrasses](https://shiny.tbep.org/seagrasstransect-dash/) [@Beck20b], [nekton communities](http://shiny.tbeptech.org/nekton-dash) [@Beck20c], and [tidal creeks](https://shiny.tbep.org/tidalcreek-dash/) [@Beck20d].  Resource management agencies or similar institutions could follow this approach to facilitate development of front-end products for more informed decision-making. 

![The TBEP water quality dashboard, demonstrating use of the `tbeptools` R package to generate summary plots for specific bay segments.\label{fig:dashex}](paper_files/figure-latex/dashex.JPG){width='100%'}

# Example usage

The function names were chosen with a typical analysis workflow in mind, where functions are available to `read` data from a source (typically from an online repository or stable URL), `anlz` to analyze the imported data using methods in existing technical documents or published papers, and to `show` the results as a summary graphic for use by environmental managers.  The functions are used to report on water quality [@Beck21b], fisheries [@Schrandt21], benthic condition [@Karlen20], tidal creeks [@Wessel21], and seagrass transect data [@Sherwood17]. The [vignettes](https://tbep-tech.github.io/tbeptools/articles/intro.html) for the package are topically organized to describe the functions that apply to each of the indicators.

The following example demonstrates use of a subset of the functions for water quality data to read a file from the Hillsborough County Environmental Protection Commission long-term monitoring dataset (available from <https://www.tampabay.wateratlas.usf.edu/>), analyze monthly and annual averages by major bay segments of Tampa Bay, and plot an annual time series for one of the bay segments.


```r
# load the package
library(tbeptools)

# read current data
wqdat <- read_importwq(xlsx = 'wqdata.xlsx', download_latest = TRUE)
wqdat
```

```
## # A tibble: 26,611 x 22
##   bay_segment epchc_station SampleTime             yr    mo
##   <chr>               <dbl> <dttm>              <dbl> <dbl>
## 1 HB                      6 2021-06-08 10:59:00  2021     6
## 2 HB                      7 2021-06-08 11:13:00  2021     6
## 3 HB                      8 2021-06-08 14:15:00  2021     6
## 4 MTB                     9 2021-06-08 13:14:00  2021     6
## 5 MTB                    11 2021-06-08 11:30:00  2021     6
## # ... with 26,606 more rows, and 17 more variables:
## #   Latitude <dbl>, Longitude <dbl>, Total_Depth_m <dbl>,
## #   Sample_Depth_m <dbl>, tn <dbl>, tn_q <chr>, sd_m <dbl>,
## #   sd_raw_m <dbl>, sd_q <chr>, chla <dbl>, chla_q <chr>,
## #   Sal_Top_ppth <dbl>, Sal_Mid_ppth <dbl>,
## #   Sal_Bottom_ppth <dbl>, Temp_Water_Top_degC <dbl>,
## #   Temp_Water_Mid_degC <dbl>, ...
```

```r
# analyze monthly and annual means by bay segment
avedat <- anlz_avedat(wqdat)
avedat
```

```
## $ann
## # A tibble: 584 x 4
##      yr bay_segment var         val
##   <dbl> <chr>       <chr>     <dbl>
## 1  1974 HB          mean_chla 22.4 
## 2  1974 LTB         mean_chla  4.24
## 3  1974 MTB         mean_chla  9.66
## 4  1974 OTB         mean_chla 10.2 
## 5  1975 HB          mean_chla 27.9 
## # ... with 579 more rows
## 
## $mos
## # A tibble: 4,484 x 5
##   bay_segment    yr    mo var         val
##   <chr>       <dbl> <dbl> <chr>     <dbl>
## 1 HB           1974     1 mean_chla 36.2 
## 2 LTB          1974     1 mean_chla  1.75
## 3 MTB          1974     1 mean_chla 11.5 
## 4 OTB          1974     1 mean_chla  4.4 
## 5 HB           1974     2 mean_chla 42.4 
## # ... with 4,479 more rows
```

```r
# show annual time series of chlorophyll for Hillsborough bay segment
show_thrplot(wqdat, bay_segment = 'HB', yrrng = c(1975, 2020))
```

![](paper_files/figure-latex/thrplotex-1.jpeg)<!-- --> 


# Acknowledgements

We acknowledge our many local and regional partners for their continuing collaborative efforts in working towards a healthy Tampa Bay, in particular the [Tampa Bay Nitrogen Management Consortium](https://tbep.org/our-work/boards-committees/nitrogen-management-consortium/). The `tbeptools` software would not be possible without data provided by our partners.  We also thank two reviewers for providing useful feedback that improved the manuscript.

# References
# Issues reporting guide

Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://stackoverflow.com/> or <https://community.rstudio.com/>.

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](https://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

---

Brief description of the problem

```r
# insert reprex here
```
# Contributing to tbeptools

This outlines how to propose a change to tbeptools. For more detailed
info about contributing to this, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

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
*  New code should follow the tidyverse [style guide](https://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the tbeptools project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See tidyverse [development contributing guide](https://rstd.io/tidy-contrib)
for further details.
---
title: 'tbeptools: An R package for synthesizing estuarine data for environmental research'
tags:
  - R
  - estuary
  - Tampa Bay
  - water quality
  - reporting
authors:
  - name: Marcus W. Beck^[Corresponding author]
    orcid: 0000-0002-4996-0059
    affiliation: 1 
  - name: Meagan N. Schrandt
    orcid: 0000-0002-0482-5072
    affiliation: 2
  - name: Michael R. Wessel
    affiliation: 3
  - name: Edward T. Sherwood
    orcid: 0000-0001-5330-302X
    affiliation: 1
  - name: Gary E. Raulerson
    orcid: 0000-0002-5920-5743
    affiliation: 1  
  - name: Adhokshaja Achar Budihal Prasad
    orcid: 0000-0002-6485-6858
    affiliation: 4
  - name: Benjamin D. Best
    orcid: 0000-0002-2686-0784
    affiliation: 5
affiliations:
  - name: Tampa Bay Estuary Program, St. Petersburg, Florida, USA
    index: 1
  - name: Fish and Wildlife Research Institute, Florida Fish and Wildlife Conservation Commission, St. Petersburg, Florida, USA
    index: 2
  - name: Janicki Environmental, Inc., St. Petersburg, Florida, USA
    index: 3
  - name: University of South Florida, Tampa, Florida, USA
    index: 4
  - name: EcoQuants, LLC, Santa Barbara, California, USA
    index: 5
date: 27 April 2021
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

```{r setup, warning = F, message = F, echo = F}
knitr::opts_chunk$set(message = F, warning = F, echo = T, cache = F, dev.args = list(family = 'serif'), dpi = 400, dev = 'jpeg')

options(tibble.width = 60, tibble.print_max = 5, tibble.print_min = 5)
```

# Summary

Many environmental programs report on the status and trends of natural resources to inform management decisions for protecting or restoring environmental condition.  The National Estuary Program (NEP) in the United States is one example of a resource management institution focused on "estuaries of national significance" that provides place-based solutions to managing coastal resources.  There are 28 NEPs in the United States, each with similar but location-specific programmatic goals to address environmental challenges related to water quality, alteration of hydrologic flows, invasive species, climate change,  declines in fish and wildlife populations, pathogens and other contaminants, and stormwater management.  A critical need of each NEP is the synthesis of data from disparate sources that can inform management response to address these environmental challenges. 

The Tampa Bay Estuary Program (TBEP) in Florida, USA is responsible for developing and implementing a place-based plan to sustain historical and future progress in the restoration of Tampa Bay [@tbep1017].  The needs of TBEP for reporting on indicators of environmental condition are similar to other environmental organizations.  Multiple local and regional partners collect data that are used for different reporting products.  Without data synthesis tools that are transparent, accessible, and reproducible, NEP staff and colleagues waste time and resources compiling information by hand.  The `tbeptools` R software package can be used for routine development of reporting products, allowing for more efficient use of limited resources and a more effective approach to communicate research to environmental decision-makers. Functions in `tbeptools` also support the creation of content for interactive, online dashboards that can facilitate more informed decisions without requiring an intimate understanding of the R programming language or the methods for analysis.       

The `tbeptools` package also addresses challenges associated with data retrieval and assessment of environmental data relative to important policy or management targets.  For each environmental indicator, functions are included to import required data directly from sources, removing the need to manually obtain information prior to reporting.  These functions are integrated into summary report cards that are generated automatically through continuous integration services (i.e., [GitHub Actions](https://github.com/features/actions)) that free the analyst from external downloads, analysis, and copying of results that can introduce errors in reporting.  Similar packages provide seamless access to database services [e.g., `dataRetrieval`, @DeCicco21], but few packages link these data sources directly to analysis and reporting as in `tbeptools` [but see `wqindex`, @Thorley18].  Management targets and regulatory thresholds based on summary assessments of data, either direct from sources or included as supplementary data in the package, are also hard-coded into the functions. 

# Statement of need

The `tbeptools` R package was developed to automate data synthesis and analysis for many of the environmental indicators for Tampa Bay, with more general application to commonly available datasets for estuaries.  The functions in the package were developed to extract methods from existing technical documents and to make them available in an open source programming environment.  By making these tools available as an R package, routine assessments are now accomplished more quickly and other researchers can use the tools to develop more specific analysis pipelines.  

Most of the NEPs do not have analysis software to operationalize data import, analysis, and plotting for reporting.  Recently, a similar software package, `peptools` [@Beck21a], was developed for the Peconic Estuary Partnership (New York, USA) using many of the functions in `tbeptools` to develop reporting products for a new water quality monitoring program.  This successful technology transfer demonstrates the added value of presenting these methods in an open source environment available for discovery and reuse by others.  We expect other NEPs to begin using these tools as their application becomes more widespread among estuarine researchers.

Beyond the NEPs, `tbeptools` is an effective example of an R package for implementing technical methods in existing literature and reports that can be used to support environmental monitoring and assessment needs for science-based decisions.  To this end, the `tbeptools` package was also created to support the development of online dashboards created in R Shiny [@Chang21].  Dashboards are powerful tools to increase accessibility for end users to engage with scientific products without the need to understand technical details in their creation.  However, providing the underlying methods as source code in an R package increases transparency and reproducibility of reporting products if users require a more detailed understanding of how the content was created.  Currently, the `tbeptools` package supports dashboards created by TBEP for the assessment of [water quality](https://shiny.tbeptech.org/wq-dash/) [Figure \ref{fig:dashex}, @Beck20a], [seagrasses](https://shiny.tbep.org/seagrasstransect-dash/) [@Beck20b], [nekton communities](http://shiny.tbeptech.org/nekton-dash) [@Beck20c], and [tidal creeks](https://shiny.tbep.org/tidalcreek-dash/) [@Beck20d].  Resource management agencies or similar institutions could follow this approach to facilitate development of front-end products for more informed decision-making. 

![The TBEP water quality dashboard, demonstrating use of the `tbeptools` R package to generate summary plots for specific bay segments.\label{fig:dashex}](paper_files/figure-latex/dashex.JPG){width='100%'}

# Usage

The function names were chosen with a typical analysis workflow in mind, where functions are available to `read` data from a source (typically from an online repository or stable URL), `anlz` to analyze the imported data using methods in existing technical documents or published papers, and to `show` the results as a summary graphic for use by environmental managers.  The functions also include a suffix that is semi-descriptive of the type of data for which they apply (e.g., water quality, seagrass) and what they do (e.g., import, format).  The package [landing page](https://tbep-tech.github.io/tbeptools/) describes more of these suffixes.  Overall, the functions are used to report on water quality [@Beck21b], fisheries [@Schrandt21], benthic condition [@Karlen20], tidal creeks [@Wessel21], and seagrass transect data [@Sherwood17]. The [vignettes](https://tbep-tech.github.io/tbeptools/articles/intro.html) for the package are topically organized to describe the functions that apply to each of the indicators.

The following example demonstrates use of a subset of the functions for water quality data to read a file from the Hillsborough County Environmental Protection Commission long-term monitoring dataset (available from <https://www.tampabay.wateratlas.usf.edu/>), analyze monthly and annual averages by major bay segments of Tampa Bay, and plot an annual time series for one of the bay segments.

```{r thrplotex, fig.height = 3.5, fig.width = 7, dev = 'jpeg'}
# load the package
library(tbeptools)

# read current data
wqdat <- read_importwq(xlsx = 'wqdata.xlsx', download_latest = TRUE)
wqdat

# analyze monthly and annual means by bay segment
avedat <- anlz_avedat(wqdat)
avedat

# show annual time series of chlorophyll for Hillsborough bay segment
show_thrplot(wqdat, bay_segment = 'HB', yrrng = c(1975, 2020))
```
```{r, echo = F, results = 'hide'}
file.remove('wqdata.xlsx')
file.copy('paper_files/figure-latex/thrplotex-1.jpeg', '../man/figures/')
```

# Acknowledgements

We acknowledge our many local and regional partners for their continuing collaborative efforts in working towards a healthy Tampa Bay, in particular the [Tampa Bay Nitrogen Management Consortium](https://tbep.org/our-work/boards-committees/nitrogen-management-consortium/). The `tbeptools` software would not be possible without data provided by our partners.  We also thank two reviewers for providing useful feedback that improved the manuscript.

# References
---
title: "Intro to TBEP Tools"
csl: stylefile.csl
bibliography: refs.bib
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Intro to TBEP Tools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = F, warning = F, 
  fig.align = 'center'
)

# libraries
library(tbeptools) 
library(bookdown)
library(ggplot2)
library(dplyr)
library(knitr)

# spelling::spell_check_files("vignettes/intro.Rmd")
```

## Background

Dashboard: https://shiny.tbep.org/wq-dash

This vignette provides an overview of the functions in tbeptools that can be used to work with water quality data in Tampa Bay.  View the other vignettes for topical introductions to other reporting products (e.g., seagrasess, tidal creeks, etc.).

The environmental recovery of Tampa Bay is an exceptional success story for coastal water quality management. Nitrogen loads in the mid 1970s have been estimated at 8.2 million kg/yr, with approximately 5.5 million kg/yr entering the upper Bay alone [@Poe05,@Greening06].  Reduced water clarity associated with phytoplankton biomass contributed to a dramatic reduction in the areal coverage of seagrass [@Tomasko05] and development of hypoxic events, causing a decline in benthic faunal production [@Santos80].  Extensive efforts to reduce nutrient loads to the Bay occurred by the late 1970s, with the most notable being improvements in infrastructure for wastewater treatment in 1979.  Improvements in water clarity and decreases in chlorophyll concentrations were observed Bay-wide in the 1980s, with conditions generally remaining constant to present day [@Beck15].

Tracking changes in environmental condition from the past to present day would not have been possible without a long-term monitoring dataset. Data have been collected monthly by the Environmental Protection Commission of Hillsborough County since 1974 [@Sherwood16;@TBEP17].  Samples are taken at forty-five stations by water collection or monitoring sonde at bottom, mid- or surface depths, depending on parameter.  The locations of monitoring stations are fixed and cover the entire Bay from the uppermost mesohaline sections to the lowermost euhaline portions that have direct interaction with the Gulf of Mexico.  Up to 515 observations are available for different parameters at each station, e.g., nitrogen, chlorophyll-a, and secchi depth. 

Data collected from the monitoring program are processed and maintained in a spreadsheet titled `RWMDataSpreadsheet_ThroughCurrentReportMonth.xlsx` at <ftp://ftp.epchc.org/EPC_ERM_FTP/WQM_Reports/>.  These data include observations at all stations and for all parameters throughout the period of record.  To date, there have been no systematic tools for importing, analyzing, and reporting information from these data. The **tbeptools** package provides was developed to address this need.

```{r tbmap, out.width = '80%', echo = F, fig.cap = 'Locations of long-term monitoring stations in Tampa Bay. The Bay is separated into four segments defined by chemical, physical, and geopolitical boundaries.'}
knitr::include_graphics('tb_map.png')
```

## Read

The main function for importing water quality data is `read_importwq()`.  This function downloads the latest file if one is not already available at the location specified by the `xlsx` input argument.

First, create a character path for the location of the file.  If one does not exist, specify a desired location and name for the downloaded file.  Here, we want to put the file in the vignettes folder and name is current_results.xls.  Note that this file path is relative to the root working directly for the current R session.  You can view the working directory with `getwd()`.

```{r}
xlsx <- 'vignettes/current_results.xls'
```

Now we pass this `xlsx` object to the `read_importwq()` function. 

```{r, eval = F}
ecpdata <- read_importwq(xlsx)
```
```
#> Error in read_importwq("empty") : file.exists(xlsx) is not TRUE
```

We get an error message from the function indicating that the file is not found. This makes sense because the file doesn't exist yet, so we need to tell the function to download the latest file.  This is done by changing the `download_latest` argument to `TRUE` (the default is `FALSE`). 

```{r, eval = F}
ecpdata <- read_importwq(xlsx, download_latest = TRUE)
```
```
#> File vignettes/current_results.xls does not exist, replacing with downloaded file...

#> trying URL 'ftp://ftp.epchc.org/EPC_ERM_FTP/WQM_Reports/RWMDataSpreadsheet_ThroughCurrentReportMonth.xlsx'
 length 24562051 bytes (23.4 MB)
```
Now we get the same message, but with an indication that file on the server is being downloaded. We'll have the data downloaded and saved to the `epcdata` object after it finishes downloading. 

If we try to run the function again after downloading the data from the server, we get the following message.  This check is done to make sure that the data are not unnecessarily downloaded if the current matches the file on the server.

```{r, eval = F}
ecpdata <- read_importwq(xlsx, download_latest = TRUE)
```
```
#> File is current...
```

Every time that tbeptools is used to work with the monitoring data, `read_importwq()` should be used to import the data. You will always receive the message `File is current...` if your local file matches the one on the server.  However, new data are regularly collected and posted on the server.  If `download_latest = TRUE` and your local file is out of date, you will receive the following message:

```
#> Replacing local file with current...
```

The final argument `na` indicates which fields in the downloaded spreadsheet are treated as blank values and assigned to `NA`. Any number of strings can be added to this function to replace fields with `NA` values.  

After the data are successfully imported, you can view them from the assigned object: 

```{r}
epcdata
```

These data include the bay segment name, station number, sample time, year, month, latitude, longitude, station depth, sample depth, secchi depth, and chlorophyll.  Note that the monitoring data include additional parameters.  Chlorophyll and secchi depth are currently the only parameters returned by `read_importwq()` given the reporting indicators used below. 

An import function is also available to download and format phytoplankton cell count data.  The `read_importphyto()` function works similarly as the import function for the water quality data.  Start by specifying a path where the data should be downloaded and set `download_latest` to `TRUE`.  This function will download and summarize data from the file `PlanktonDataList_ThroughCurrentReportMonth.xlsx` on the EPC website.

```{r eval = F}
xlsx <- 'phyto_data.xlsx'
phytodata <- read_importphyto(xlsx, download_latest = T)
```
```
#> File vignettes/phyto_data.xlsx does not exist, replacing with downloaded file...

#> trying URL 'ftp://ftp.epchc.org/EPC_ERM_FTP/WQM_Reports/PlanktonDataList_ThroughCurrentReportMonth.xlsx'
 length 12319508 bytes (11.7 MB)
```

After the phytoplankton data are successfully imported, you can view them from the assigned object: 

```{r phyto_import, echo = F, message = F, include = F}
# local file path
xlsx <- 'phyto_data.xlsx'

# load data and some light formatting
phytodata <- read_importphyto(xlsx, download_latest = T)
```
```{r}
phytodata
```

These data are highly summarized from the raw data file available online.  Cell counts (as number of cells per 0.1mL) for selected taxa are summed for each station by quarters (i.e., Jan/Feb/Mar, Apr/May/Jun, etc.).  The quarter is indicated in the `yrqrt` column specified by the starting date of each quarter (e.g., `1975-07-01` is the quarter Jul/Aug/Sep for 1975).  These data are primarily used to support analyses in the water quality dashboard: <https://shiny.tbep.org/wq-dash/>

## Analyze {.tabset}

The functions `anlz_avedat()` and `anlz_avedatsite()` summarize the station data by bay segments or by sites, respectively.  Both functions return annual means for chlorophyll and light attenuation (based on Secchi depth measurements) and monthly means by year for chlorophyll and light attenuation.  These summaries are then used to determine if bay segment targets for water quality are met using the `anlz_attain()` and `anlz_attainsite()` function.

Here we use `anlz_avedat()` to summarize the data by bay segment to estimate annual and monthly means for chlorophyll and light attenuation.  The output is a two-element list for the annual (`ann`) and monthly (`mos`) means by segment.

```{r}
avedat <- anlz_avedat(epcdata)
avedat
```

This output can then be further analyzed with `anlz_attain()` to determine if the bay segment outcomes are met in each year.  The results are used by the plotting functions described below.  In short, the `chl_la` column indicates the categorical outcome for chlorophyll and light attenuation for each segment.  The outcomes are integer values from zero to three.  The relative exceedances of water quality thresholds for each segment, both in duration and magnitude, are indicated by higher integer values.  

```{r}
anlz_attain(avedat)
```

Similar information can be obtained for individual sites using `anlz_avedatsite()` and `anlz_attainsite()`.  The main difference is that a yes/no column `met`is added that indicates only if the target was above or below the segment threshold for each site.

```{r}
anlz_avedatsite(epcdata) %>% anlz_attainsite
```

## Show

External package libraries in R can be used to plot the time series data.  Here's an example using the popular [ggplot2](https://ggplot2.tidyverse.org/) package.  Some data wrangling with the [dplyr](https://dplyr.tidyverse.org/) is done first to filter the data we want to plot.

```{r, fig.height = 3, fig.width = 8}
toplo <- epcdata %>% 
  filter(epchc_station == '52')

ggplot(toplo, aes(x = SampleTime, y = chla)) + 
  geom_line() + 
  geom_point() + 
  scale_y_log10() + 
  labs(
    y = 'Chlorophyll-a concentration (ug/L)', 
    x = NULL, 
    title = 'Chlorophyll trends',
    subtitle = 'Hillsborough Bay station 52, all dates'
    ) + 
  theme_bw()
```

The `show_thrplot()` function provides a more descriptive assessment of annual trends for a chosen bay segment relative to defined targets or thresholds. In this plot we show the annual averages across stations Old Tampa bay (`bay_segment = "OTB"`) for chlorophyll (`thr = "chla"`).  The red line shows annual trends and the horizontal blue lines indicate the thresholds and targets for chlorophyll-a that are specific to Old Tampa Bay.  The dashed and dotted blue lines indicate +1 and +2 standard errors for the management target shown by the filled line.  The target and standard errors are considered when identifying the annual segment outcome for chlorophyll.

```{r, fig.height = 5, fig.width = 8}
show_thrplot(epcdata, bay_segment = "OTB", thr = "chla")
```

We can show the same plot but for light attenuation by changing the `thr = "chla"` to `thr = "la"`.  Note the change in the horizontal reference lines for the light attenuation target.

```{r, fig.height = 5, fig.width = 8}
show_thrplot(epcdata, bay_segment = "OTB", thr = "la")
```

The year range to plot can also be specified using the `yrrng` argument, where the default is the year range from `epcdata`.

```{r, fig.height = 5, fig.width = 8}
show_thrplot(epcdata, bay_segment = "OTB", thr = "la", yrrng = c(2000, 2018))
```

The `show_thrplot()` function uses results from the `anlz_avedat()` function.  For example, you can retrieve the values from the above plot as follows: 

```{r}
epcdata %>% 
  anlz_avedat %>% 
  .[['ann']] %>% 
  filter(bay_segment == 'OTB') %>% 
  filter(var == 'mean_la') %>% 
  filter(yr >= 2000 & yr <= 2018)
```

Similarly, the `show_boxplot()` function provides an assessment of seasonal changes in chlorophyll or light attenuation values by bay segment.  The most recent year is highlighted in red by default. This allows a simple evaluation of how the most recent year compared to historical averages.  The large exceedance value is shown in blue text and as the dotted line.  This corresponds to a "large" magnitude change of +2 standard errors above the bay segment threshold and is the same dotted line shown in `show_thrplot()`.    

```{r, fig.height = 5, fig.width = 8}
show_boxplot(epcdata, param = 'chla', bay_segment = "OTB")
show_boxplot(epcdata, param = 'la', bay_segment = "HB")
```

A different subset of years and selected year of interest can also be viewed by changing the `yrrng` and `yrsel` arguments.  Here we show 1980 compared to monthly averages from 2008 to 2018. 

```{r, fig.height = 5, fig.width = 8}
show_boxplot(epcdata, param = 'chla', bay_segment = "OTB", yrrng = c(2008, 2018), yrsel = 1980)
```

The `show_thrplot()` function is useful to understand annual variation in chlorophyll and light attenuation relative to management targets for each bay segment.  The information from these plots can provide an understanding of how the annual reporting outcomes are determined.  As noted above, an outcome integer from zero to three is assigned to each bay segment for each annual estimate of chlorophyll and light attenuation.  These outcomes are based on both the exceedance of the annual estimate above the threshold or target (blue lines in `show_thrplot()`) and duration of the exceedance for the years prior.  The following graphic describes this logic [@tbep0400]. 

```{r, echo = F, fig.cap = 'Outcomes for annual estimates of water quality are assigned an integer value from zero to three depending on both magnitude and duration of the exceedence.', out.width = '80%'}
knitr::include_graphics('outints.PNG')
```

These outcomes are assigned for both chlorophyll and light attenuation. The duration criteria are determined based on whether the exceedance was observed for years prior to the current year. The exceedance criteria for chlorophyll and light-attenuation are specific to each segment.  The tbeptools package contains a `targets` data file that is a reference for determining annual outcomes.  This file is loaded automatically with the package and can be viewed from the command line.

```{r}
targets
```

The final plotting function is `show_matrix()`, which creates an annual reporting matrix that reflects the combined outcomes for chlorophyll and light attenuation. Tracking the attainment of bay segment specific targets for these indicators provides the framework from which bay management actions are developed and initiated.  For each year and segment, a color-coded management action is assigned:

<span style="color:#33FF3B; text-shadow: 0 0 3px #333;">__Stay the Course__</span>: Continue planned projects. Report data via annual progress reports and Baywide Environmental Monitoring Report. 

<span style="color:#F9FF33; text-shadow: 0 0 3px #333;">__Caution__</span>: Review monitoring data and nitrogen loading estimates. Begin/continue TAC and Management Board development of specific management recommendations.

<span style="color:#FF3333; text-shadow: 0 0 3px #333;">__On Alert__</span>: Finalize development and implement appropriate management actions to get back on track.

The management category or action is based on the combination of outcomes for chlorophyll and light attenuation [@tbep0400].

```{r, echo = F, fig.cap = 'Management action categories assigned to each bay segment and year based on chlorophyll and light attenuation outcomes.', out.width = '80%'}
knitr::include_graphics('matrixcats.PNG')
```

The results can be viewed with `show_matrix()`.

```{r, fig.height = 8, fig.width = 3}
show_matrix(epcdata)
```

The matrix is also a `ggplot` object and its layout can be changed using `ggplot` elements. Note the use of `txtsz = NULL` to remove the color labels. 

```{r, fig.height = 1.5, fig.width = 8}
show_matrix(epcdata, txtsz = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = sort(unique(epcdata$yr))) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
```

If preferred, the matrix can also be returned in an HTML table that can be sorted and scrolled. Only the first ten rows are shown by defaul.  The default number of rows (10) can be changed with the \code{nrows} argument.  Use a very large number to show all rows.

```{r}
show_matrix(epcdata, asreact = TRUE)
```

A plotly (interactive, dynamic plot) can be returned by setting the `plotly` argument to `TRUE`. 

```{r, fig.height = 8, fig.width = 3}
show_matrix(epcdata, plotly = TRUE)
```

Results can also be obtained for a selected year. Outcomes can be returned in tabular format with `anlz_yrattain()`.  This table also shows segment averages for chlorophyll and light attenuation, including the associated targets. 

```{r}
anlz_yrattain(epcdata, yrsel = 2018)
```

A map showing if individual sites achieved chlorophyll targets can be obtained with `show_sitemap()`.  The station averages for chlorophyll for the selected year are shown next to each point.  Stations in red failed to meet the segment target.    

```{r, fig.height = 7, fig.width = 5}
show_sitemap(epcdata, yrsel = 2018)
```

The `show_sitemap()` function also includes an argument to specify a particular monthly range for the selected year.  If this option is chosen, averages are shown as continuous values at each station. 

```{r, fig.height = 7, fig.width = 6}
show_sitemap(epcdata, yrsel = 2018, mosel = c(7, 9))
```

Bay segment exceedances can also be viewed in a matrix using `show_wqmatrix()`.  The thresholds for these values correspond to the Florida DEP criteria (or a large exceedance defined as +2 standard errors above the segment target).  

```{r, fig.height = 8, fig.width = 3}
show_wqmatrix(epcdata)
```

By default, the `show_wqmatrix()` function returns chlorophyll exceedances by segment.  Light attenuation exceedances can be viewed by changing the `param` argument.  

```{r, fig.height = 8, fig.width = 3}
show_wqmatrix(epcdata, param = 'la')
```

The results from `show_matrix()` and `show_wqmatrix()` can be combined for an individual segment using the `show_segmatrix()` function.  This is useful to understand which water quality parameter is driving the management outcome for a given year. The plot shows the light attenuation and chlorophyll outcomes from `show_wqmatrix()` next to the segment management outcomes from `show_matrix()`.  Only one segment can be plotted for each function call. 

```{r, fig.height = 8, fig.width = 2.5}
show_segmatrix(epcdata, bay_segment = 'OTB')
```

Finally, all segment plots can be shown together using the `show_segplotly()` function that combines chlorophyll and secchi data for a given segment.  This function combines outputs from `show_thrplot()` and `show_segmatrix()`. The final plot is interactive and can be zoomed by dragging the mouse pointer over a section of the plot. Information about each cell or value can be seen by hovering over a location in the plot.

```{r, out.width = '100%', fig.height = 6, fig.width = 11}
show_segplotly(epcdata, width = 1000, height = 600)
```

From these plots, we can quickly view a summary of the environmental history of water quality in Tampa Bay.  Degraded conditions were common early in the period of record, particularly for Old Tampa Bay and Hillsborough Bay.  Conditions began to improve by the late 1980s and early 1990s, with good conditions persisting to present day. However, recent trends in Old Tampa Bay have shown conditions changing from "stay the course" to "caution".  

# References

---
title: "Tampa Bay Nekton Index"
csl: stylefile.csl
bibliography: tbnirefs.bib
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Tampa Bay Nekton Index}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = F, warning = F, 
  fig.align = 'center'
)

# libraries
library(tbeptools) 
library(bookdown)
library(patchwork)
library(dplyr)
library(knitr)
library(mapview)
library(ggplot2)
library(sf)

st_crs(fimstations) <- 4326

# spelling::spell_check_files("vignettes/tbni.Rmd")
```

## Background

The Tampa Bay Nekton Index (TBNI) [@tbep0418;@tbep1219]) is a multimetric assessment method that quantifies the ecological health of the nekton community in Tampa Bay.  The index provides a complementary approach to evaluating environmental condition that is supported by other assessment methods currently available for Tampa Bay (e.g., water quality report card, Benthic index, etc.). The tbeptools package includes several functions described below to import data required for the index, analyze the data to calculate metrics and index scores, and plot the results to view trends over time.  Each of the functions are described in detail below. 

The TBNI uses catch data from the Florida Fish and Wildlife Conservation Commission (FWC) Fish and Wildlife Research Institute’s (FWRI) [Fisheries-Independent Monitoring (FIM) program](https://myfwc.com/research/saltwater/reef-fish/monitoring/fim-stratified-random-sampling/).  Catch results from a center-bag seine have the longest and most consistent record in the FIM database and were used to develop the TBNI.  These include counts and taxa identification for individuals caught in near shore areas, generally as early recruits, juveniles, and smaller-bodied nekton. All fish and selected invertebrates are identified to the lowest possible taxon (usually species), counted, and a subset are measured. Current protocols were established in 1998 and TBNI estimates are unavailable prior to this date.  

## Data import and included datasets

Data required for calculating TBNI scores can be imported into the current R session using the `read_importfim()` function.  This function downloads the latest FIM file from an FTP site if the data have not already been downloaded to the location specified by the input arguments.

To download the data, first create a character path for the location of the file. If one does not exist, specify a desired location and name for the downloaded file. Here, we want to put the file on the desktop in our home directory and name it `fimdata.csv`. 

```{r, eval = F}
csv <- '~/Desktop/fimdata.csv'
fimdata <- read_importfim(csv)
```

Running the above code will return the following error: 

```{r, eval = F}
#> Error in read_importfim(csv) : file.exists(csv) is not TRUE
```

We get an error message from the function indicating that the file is not found. This makes sense because the file doesn’t exist yet, so we need to tell the function to download the latest file. This is done by changing the `download_latest` argument to `TRUE` (the default is `FALSE`).

```{r, eval = F}
fimdata <- read_importfim(csv, download_latest = T)
```
```{r, eval = F}
#> File ~/Desktop/fimdata.csv does not exist, replacing with downloaded file...

#> trying URL 'ftp://ftp.floridamarine.org/users/fim/tmac/NektonIndex/TampaBay_NektonIndexData.csv' length 11083878 bytes (10.6 MB)
```

Now we get an indication that the file on the server is being downloaded. When the download is complete, we’ll have the data downloaded and saved to the `fimdata` object in the current R session.

If we try to run the function again after downloading the data from the server, we get the following message. This check is done to make sure that the data are not unnecessarily downloaded if the current matches the file on the server.

```{r, eval = F}
fimdata <- read_importfim(csv, download_latest = T)
```
```{r, eval = F}
#> File is current..
```

Every time that tbeptools is used to work with the FIM data, `read_importfim()` should be used to import the data. You will always receive the message `File is current...` if your local file matches the one on the server. However, new data are regularly collected and posted on the server. If `download_latest = TRUE` and your local file is out of date, you will receive the following message:

```{r, eval = F}
#> Replacing local file with current...
```

After the data are successfully imported, you can view them from the assigned object:

```{r}
head(fimdata)
```

The imported data are formatted for calculating the TBNI.  The columns include a `Reference` for the FIM sampling site, the sampling date, sampling `Zone`, sampling `Grid`, `NODCCODE` as a unique identifier for species, sample year, sample month, total catch as `Total_N`, scientific name, a column indicating if the species is included in the index, and several columns indicating species-specific information required for the metrics. For the final columns, a separate lookup table is provided in the package that is merged with the imported FIM data.  This file, `tbnispp`, can be viewed anytime the package is loaded: 

```{r}
head(tbnispp)
```

The `read_importfim()` function can also return a simple features object of sampled stations in the raw FIM data by setting \code{locs = TRUE}.  These data are matched to the appropriate bay segments for tabulating TBNI scores.  The resulting dataset indicates where sampling has occurred and can be mapped with the `mapview()` function.  For ease of use, a dataset named `fimstations` is included in tbeptools.      

```{r, eval = F}
fimstations <- read_importfim(csv, download_latest = TRUE, locs = TRUE)
mapview(fimstations, zcol = 'bay_segment')
```
```{r, echo = F, out.width = '100%'}
mapview(fimstations, zcol = 'bay_segment')
```

The `read_importfim()` function processes the observed data as needed for the TBNI, including merging the rows with the `tbnispp` and `fimstations` data. Once imported, the metrics and scores can be calculated.  

## Calculating metrics and TBNI scores

Metrics and scores for the Tampa Bay Nekton Index can be calculated using two functions.  The `anlz_tbnimet()` function calculates all raw metrics and the `anlz_tbniscr()` function calculates scored metrics and the final TBNI score.  Both functions use the imported and formatted FIM data as input.  

The TBNI includes five metrics that were sensitive to stressor gradients and provide unique information about Nekton community response to environmental conditions.  The metrics include: 

* `NumTaxa`: Species richness

* `BenthicTaxa`: Species richness for benthic taxa

* `TaxaSelect`: Number of "selected" species (i.e., commercially and/or recreationally important)

* `NumGuilds`: Number of trophic guilds

* `Shannon`: Shannon Diversity (H)

Raw metrics are first calculated from the observed data and then scaled to a standard score from 0 - 10 by accounting for expected relationships to environmental gradients and 5th/95th percentiles of the distributions.  The final TBNI score is the summed average of the scores ranging from 0 - 100.  

The raw metrics, scored metrics, and final TBNI score is returned with the `anlz_tbniscr()` function. 

```{r}
tbniscr <- anlz_tbniscr(fimdata)
head(tbniscr)
```

The five metrics chosen for the TBNI were appropriate for the Tampa Bay dataset and were selected from a larger pool of candidate metrics.  All potential metrics can be calculated using the `anlz_tbnimet()` function.  These metrics can be used in standalone assessments or for developing a Nekton index outside of Tampa Bay.  The argument `all = TRUE` must be used to return all metrics, otherwise only the selected five for the TBNI are returned.

```{r}
tbnimet <- anlz_tbnimet(fimdata, all = T)
head(tbnimet)
```

## Plotting results

The TBNI scores can be viewed as annual averages using the `show_tbniscr()`, `show_tbniscrall()` and `show_tbnimatrix()` functions.  The `show_tbniscr()` creates a line graph of values over time for each bay segment, whereas the `show_tbniscrall()` function plots an overall average across bay segments over time.  The `show_tbnimatrix()` plots the annual bay segment averages as categorical values in a conventional "stoplight" graphic.  The input to each function is the output from the `anlz_tbniscr()` function.  

```{r, fig.height = 5, fig.width = 7}
show_tbniscr(tbniscr)
```

```{r, fig.height = 5, fig.width = 7}
show_tbniscrall(tbniscr)
```

```{r, fig.height = 7, fig.width = 3}
show_tbnimatrix(tbniscr)
```

Each of the plots can also be produced as [plotly](https://plotly.com/r/) interactive plots by setting `plotly = TRUE` inside each function. 

```{r, fig.height = 5, fig.width = 7}
show_tbniscr(tbniscr, plotly = T)
```

```{r, fig.height = 5, fig.width = 7}
show_tbniscrall(tbniscr, plotly = T)
```

```{r, fig.height = 7, fig.width = 3}
show_tbnimatrix(tbniscr, plotly = T)
```

The breakpoints for the categorical outcomes of the TBNI scores shown by the colors in each graph are based on the 33rd and 50th percentiles of the distribution of all TBNI scores calculated for Tampa Bay.  This plotting option is provided for consistency with existing TBEP reporting tools, e.g., the [water quality report card](https://tbeptech.org/TBEP_TECH_PUBS/2020/TBEP_01_20_2019_Decision_Matrix_Update.pdf) returned by `show_matrix()`.  The categorical outcomes serve as management guidelines each year for activities to support environmental resources of the Bay: <span style="color:#33FF3B; text-shadow: 0 0 3px #333;">__Stay the Course__</span>, <span style="color:#F9FF33; text-shadow: 0 0 3px #333;">__Caution__</span>, and <span style="color:#FF3333; text-shadow: 0 0 3px #333;">__On Alert__</span> [@tbep0105].  

The graphs returned by the plotting functions are `ggplot` objects that can be further modified.  They can be combined below using [patchwork](https://patchwork.data-imaginist.com/) in a single graphic showing the trends over time as both categorical outcomes in the matrix and continuous scores in the bottom plot. 

```{r, fig.height = 7, fig.width = 7}
p1 <- show_tbnimatrix(tbniscr, txtsz = NULL, rev = TRUE, position = 'bottom') +
  scale_y_continuous(expand = c(0,0), breaks = c(1998:2020)) +
  coord_flip() +
  theme(axis.text.x = element_blank())

p2 <- show_tbniscr(tbniscr)

p1 + p2 + plot_layout(ncol = 1, heights = c(0.3, 1))
```

# References
---
title: "Tampa Bay Benthic Index"
csl: stylefile.csl
bibliography: tbbirefs.bib
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Tampa Bay Benthic Index}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = F, warning = F, 
  fig.align = 'center'
)

# libraries
library(tbeptools)
library(bookdown)
library(patchwork)
library(dplyr)
library(knitr)
library(mapview)
library(ggplot2)
library(sf)
library(tibble)

# spelling::spell_check_files("vignettes/tbbi.Rmd")
```

## Background

The Tampa Bay Benthic Index (TBBI) [@tbep0620;@tbep0106]) is an assessment method that quantifies the ecological health of the benthic community in Tampa Bay.  The index provides a complementary approach to evaluating environmental condition that is supported by other assessment methods currently available for the region (e.g., water quality report card, nekton index, etc.). The tbeptools package includes several functions described below to import data required for the index and plot the results to view trends over time.  Each of the functions are described in detail below. 

The TBBI uses data from the Tampa Bay Benthic Monitoring Program as part of the Hillsborough Country Environmental Protection Commission (EPC).  The data are updated annually on a FTP site maintained by EPC, typically in December after Summer/Fall sampling.  This is the same website that hosts water quality data used for the [water quality report card](https://tbep-tech.github.io/tbeptools/articles/intro.html).  The required data for the TBBI are more extensive than the water quality report card and the data are made available as a zipped folder of csv files, available [here](ftp://ftp.epchc.org/EPC_ERM_FTP/Benthic_Monitoring/DataDeliverables/BaseTables.zip).  The process for downloading and working with the data are similar as for the other functions in tbeptools.

## Data import and included datasets

Data for calculating TBBI scores can be imported into the current R session using the `read_importbenthic()` function.  This function downloads the zipped folder of base tables used for the TBBI from the EPC FTP site if the data have not already been downloaded to the location specified by the input arguments.  

To download the data with teptools, first create a character path for the location where you want to download the zipped files. If one does not exist, specify a location and name for the downloaded file. Here, we name the folder `benthic.zip` and download it on our desktop.

```{r, eval = F}
path <- '~/Desktop/benthic.zip'
benthicdata <- read_importbenthos(path)
```

Running the above code will return the following error: 

```{r, eval = F}
#> Error in read_importbenthic() : File at path does not exist, use download_latest = TRUE
```

We get an error message from the function indicating that the file is not found. This makes sense because the file does not exist yet, so we need to tell the function to download the latest file. This is done by changing the `download_latest` argument to `TRUE` (the default is `FALSE`).

```{r, eval = F}
benthicdata <- read_importbenthic(path, download_latest = T)
```
```{r, eval = F}
#> File ~/Desktop/benthic.zip does not exist, replacing with downloaded file...

#> trying URL 'ftp://ftp.epchc.org/EPC_ERM_FTP/Benthic_Monitoring/DataDeliverables/BaseTables.zip' length 37122877 bytes (35.4 MB)
```

Now we get an indication that the file on the server is being downloaded. When the download is complete, we’ll have the data downloaded and saved to the `benthicdata` object in the current R session.

If we try to run the function again after downloading the data, we get the following message. This check is done to make sure that the data are not unnecessarily downloaded if the current file matches the file on the server.

```{r, eval = F}
benthicdata <- read_importbenthic(path, download_latest = T)
```
```{r, eval = F}
#> File is current..
```

Every time that tbeptools is used to work with the benthic data, `read_importbenthic()` should be used to import the data. You will always receive the message `File is current...` if your local file matches the one on the server. However, data are periodically updated and posted on the server. If `download_latest = TRUE` and your local file is out of date, you will receive the following message:

```{r, eval = F}
#> Replacing local file with current...
```

## Calculating TBBI scores

After the data are imported, you can view them from the assigned object.  The data are provided as a nested tibble that includes three different datasets: station information, field sample data (salinity), and detailed taxa information.  

```{r}
benthicdata
```

The individual datasets can be viewed by extracting them from the parent object using the `deframe()` function from the tibble package. 

```{r}
# see all
deframe(benthicdata)

# get only station dat
deframe(benthicdata)[['stations']]
```

The `anlz_tbbiscr()` function uses the nested `benthicdata` to estimate the TBBI scores at each site. The TBBI scores typically range from 0 to 100 and are grouped into categories that describe the general condition of the benthic community.  Scores less than 73 are considered "degraded", scores between 73 and 87 are "intermediate", and scores greater than 87 are "healthy".  Locations that were sampled but no organisms were found are assigned a score of zero and a category of "empty sample".  The total abundance (`TotalAbundance`, organisms/m2), species richness (`SpeciesRichness`) and bottom salinity (`Salinity`, psu) are also provided.  Some metrics for the TBBI are corrected for salinity and bottom measurements taken at the time of sampling are required for accurate calculation of the TBBI.

```{r}
tbbiscr <- anlz_tbbiscr(benthicdata)
tbbiscr
```

## Plotting results

The TBBI scores can be viewed as annual averages for each bay segment using the `show_tbbimatrix()` function.  The `show_tbbimatrix()` plots the annual bay segment averages as categorical values in a conventional "stoplight" graphic. A baywide estimate is also returned, one based on all samples across all locations ("All") and another weighted by the relative surface areas of each bay segment ("All (wt)"). The input to `show_tbbimatrix()` function is the output from the `anlz_tbbiscr()` function.  

```{r, fig.height = 7, fig.width = 3}
show_tbbimatrix(tbbiscr)
```

The matrix can also be produced as a [plotly](https://plotly.com/r/) interactive plot by setting `plotly = TRUE` inside the function. 

```{r, fig.height = 7, fig.width = 4}
show_tbbimatrix(tbbiscr, plotly = T)
```

# References
---
title: "Seagrass Transect Data"
csl: stylefile.csl
bibliography: seagrasstransect.bib
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Seagrass Transect Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = F, warning = F, 
  fig.align = 'center'
)

library(tbeptools)
library(mapview)

sf::st_crs(trnpts) <- 4326
sf::st_crs(trnlns) <- 4326

# spelling::spell_check_files("vignettes/seagrasstransect.Rmd")
```

Each year, TBEP partners collect seagrass transect data at fixed locations in Tampa Bay.  Data have been collected since the mid 1990s and are hosted online at the [Tampa Bay Water Atlas](https://www.tampabay.wateratlas.usf.edu/) by the University of South Florida Water Institute.  Functions are available in tbeptools for downloading, analyzing, and plotting these data. 

## Data import and included datasets

There are two datasets included in tbeptools that show the actively monitored transect locations in Tampa Bay.  The `trnpts` dataset is a point object for the starting location of each transect and the `trnlns` dataset is a line object showing the approximate direction and length of each transect beginning at each point in `trnpts`.  Each dataset also includes the `MonAgency` column that indicates which monitoring agency collects the data at each transect. 

```{r}
trnpts
trnlns
```

The two datasets are `sf()` (simple features) objects and are easily mapped with `mapview()` to view their locations.

```{r, out.width = '100%'}
cols <- c("#E16A86", "#CB7F2F", "#9F9400", "#50A315", "#00AC79", "#00AAB7", "#009ADE", "#A87BE4", "#DA65C3")

mapview(trnpts, zcol = 'MonAgency', lwd = 0, legend = F, homebutton = F, col.regions = cols) + 
  mapview(trnlns, zcol = 'MonAgency', homebutton = F, layer.name = 'Monitoring Agency', lwd = 4, color = cols)
```

The transect data can be downloaded from the Water Atlas using the `read_transect()` function.  The only required argument for this function is `training`, which indicates if you want to download training data or the complete dataset, i.e., `training = TRUE` or `training = FALSE` (default).  In the former case, a small dataset is downloaded that includes only data collected during an annual training event. These are primarily used internally by TBEP staff to assess precision among different training crews. The data are downloaded as a JSON object and formatted internally using the `read_formtransect()` function.  Shoot density is reported as number of shoots per square meter and is corrected for the quadrat size entered in the raw data. Abundance is reported as a numeric value from 0 -5 for Braun-Blanquet coverage estimates and blade length is in cm.  

```{r}
# import training data
traindat <- read_transect(training = TRUE)

# view the data
traindat
```

Change the `training` argument to `FALSE` to download the entire transect database.  This may take a few seconds.

```{r}
# import entire transct dataset as JSON
transect <- read_transect(training = FALSE)

# view the data
transect
```

The columns in the complete transect database describe the crew (`Crew`), the monitoring agency (`MonitoringAgency`), sample date (`Date`), transect name (`Transect`), the meter location for the quadrat along the transect (`Site`, m), depth at the site (`Depth`, cm), Seagrass species (`Savspecies`), distance of the seagrass edge on the transect (`SeagrassEdge`, m), the seagrass variable (`var`), average value of the variable (`aveval`), and standard deviation of the variable if appropriate (`sdval`).  

If the raw, unformatted transect data are preferred, use the `raw = TRUE` argument for `read_transect()`. 

```{r}
# raw transect data
transectraw <- read_transect(training = FALSE, raw = TRUE)

# view the data
transectraw
```

## Calculating seagrass frequency occurrence

The rest of the seagrass functions in tbeptools were developed to work with the complete database.  Only the `show_complot()` function (see below) was developed for the training data.  The rest of the functions can be used to estimate and plot frequency occurrence data.  

The `anlz_transectocc()` function summarizes frequency occurrence for all transects and dates by collapsing species results across quadrats within each transect.  Abundance and frequency occurrence are estimated as in Sherwood et al. 2017, equations 1 and 2 [@tbep0917]. In short, frequency occurrence is estimated as the number of instances a species was observed along a transect divided by the number of placements along a transect and average abundance was estimated as the sum of species-specific Braun-Blanquet scores divided by the number of placements along a transect.  The estimates are obtained for all seagrass species including Caulerpa spp., whereas all attached and drift algae species are aggregated.

```{r}
transectocc <- anlz_transectocc(transect)
transectocc
```

The second function, `anlz_transectave()`, takes the results from `anlz_transectocc()` and estimates annual results across major bay segments for all seagrass species by averaging frequency occurrence across transects.  This function is used internally within the `show_transectmatrix()` function to create summary plots. The frequency occurrence estimates are also binned into categories for simple trend assessments, e.g., red < 25%, orange 25-50%, yellow 50-75%, and green > 75%.  Results for specific bay segments and annual ranges can be filtered with the `bay_segment` and `yrrng` arguments.

```{r}
transectave <- anlz_transectave(transectocc)
transectave
```

The third function, `anlz_transectavespp()`, takes the results from `anlz_transectocc()` and estimates annual averages across major bay segments as in the last function, but results are retained for individual species.  This function is used internally within the `show_transectavespp()` function to create summary plots.  All summaries are aggregated across the selected bay segments, i.e., the default is to average by species/year across all segments.  Results for an individual bay segment can be returned with the appropriate argument, e.g., by using `bay_segment = 'OTB'` to select only Old Tampa Bay.  Results can also be filtered by specific species using the `species` argument, where the default is to return all.  *Caulerpa spp.* are also included. 

```{r}
transectavespp <- anlz_transectavespp(transectocc)
transectavespp
```

Results for individual bay segments from `anlz_transectavespp()` can be retained by setting the `by_seg` argument to `TRUE`.  Note that totals are not returned in this case. 

```{r}
transectavespp <- anlz_transectavespp(transectocc, by_seg = TRUE)
transectavespp
```

## Plotting results

There is one plotting function for the training data.  The `show_compplot()` function is used to compare training data between crews for a selected species (`species` argument) and variable (`varplo` argument). 

```{r, fig.height = 5, fig.width = 7}
show_compplot(traindat, yr = 2021, site = '1', species = 'Halodule', varplo = 'Abundance', base_size = 14)
```

The rest of the plotting functions work with the complete transect data. Data for an individual transect can be viewed with the `show_transect()` function by entering the transect (site) number, species (one to many), and variable to plot.  The plot shows relative values for the selected species and variable by distance along the transect (x-axis) and year of sampling (y-axis).  The plots provide an overall summary of temporal and spatial changes in the selected seagrass metric for an individual location.  

```{r, fig.height = 8, fig.width = 9}
show_transect(transect, site = 'S3T10', species = 'Halodule', varplo = 'Abundance')
```

The plot can also be produced as a [plotly](https://plotly.com/r/) interactive plot by setting `plotly = TRUE` inside the function Note that the size legend is removed in this option, but sizes can be viewed on mouseover of each point. 

```{r, fig.height = 8, fig.width = 8}
show_transect(transect, site = 'S3T10', species = 'Halodule', varplo = 'Abundance', plotly = T)
```

The `show_transect()` function can also be used to plot multiple species.  One to many species can be provided to the `species` argument. 

```{r, fig.height = 8, fig.width = 8}
show_transect(transect, site = 'S3T10', species = c('Halodule', 'Syringodium', 'Thalassia'), varplo = 'Abundance')
```

The plots can also be separated into facets for each species using `facet = TRUE`. This is useful to reduce overplotting of multiple species found at the same location.

```{r, fig.height = 8, fig.width = 8}
show_transect(transect, site = 'S3T10', species = c('Halodule', 'Syringodium', 'Thalassia'), varplo = 'Abundance', facet = TRUE)
```

The `show_transectsum()` function provides an alternative summary of data at an individual transect. This plot provides a quick visual assessment of how frequency occurrence or abundance for multiple species has changed over time at a selected transect.  Unlike `show_transect()`, the plot shows aggregated results across quadrats along the transect and uses summarized data from the `anlz_transectocc()` function as input.  

```{r, out.width = '100%'}
show_transectsum(transectocc, site = 'S3T10')
```

A summary matrix of frequency occurrence estimates across all species can be plotted with `show_transectmatrix()`.  This uses results from the `anlz_transectocc()` and `anlz_transectave()` functions to estimate annual averages by bay segment.  The continuous frequency occurrence estimates are binned into color categories described above, as in Table 1 in [@tbep0816]. 

```{r, fig.height = 7, fig.width = 4}
show_transectmatrix(transectocc)
```

The default color scheme is based on arbitrary breaks at 25, 50, and 75 percent frequency occurrence.  These don't necessarily translate to any ecological breakpoints. Use `neutral = TRUE` to use a neutral and continuous color palette.

```{r, fig.height = 7, fig.width = 4}
show_transectmatrix(transectocc, neutral = T)
```

The matrix can also be produced as a [plotly](https://plotly.com/r/) interactive plot by setting `plotly = TRUE` inside the function. 

```{r, fig.height = 7, fig.width = 5}
show_transectmatrix(transectocc, plotly = T)
```

Time series plots of annual averages of frequency occurrence estimates by each species can be shown with the `show_transectavespp()` function.  By default, all estimates are averaged across all bay segments for each species.  The plot is a representation of Figure 2 in [@tbep0816].

```{r, fig.height = 6, fig.width = 7}
show_transectavespp(transectocc)
```

Results for individual segments and species can be returned with the `bay_segment` and `species` arguments.  Use the argument `total = FALSE` to omit the total frequency occurrence from the plot.

```{r, fig.height = 6, fig.width = 7}
show_transectavespp(transectocc, bay_segment = 'LTB', species = c('Syringodium', 'Thalassia'), total = FALSE)
```

The plot can also be produced as a [plotly](https://plotly.com/r/) interactive plot by setting `plotly = TRUE` inside the function. 

```{r, fig.height = 6, fig.width = 7}
show_transectavespp(transectocc, bay_segment = 'LTB', species = c('Syringodium', 'Thalassia'), plotly = T)
```

As an alternative to plotting the species averages over time with `show_transectavespp()`, a table can be created by setting `asreact = TRUE`.  Filtering options that apply to the plot also apply to the table, e.g., filtering by the four major bay segments and specific year ranges.  Also note that the totals are not returned in the table. 

```{r}
show_transectavespp(transectocc, asreact = T, bay_segment = c('HB', 'OTB', 'MTB', 'LTB'), yrrng = c(2006, 2012))
```

# References
---
title: "Tidal Creek Assessment"
csl: stylefile.csl
bibliography: tidalrefs.bib
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Tidal Creek Assessment}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = F, warning = F, 
  fig.align = 'center', 
  cache = F
)

# libraries
library(tbeptools) 
library(bookdown)
library(dplyr)
library(knitr)
library(mapview)
library(plotly)

sf::st_crs(tidalcreeks) <- 4326

# spelling::spell_check_files("vignettes/tidalcreeks.Rmd")
```

## Background

Dashboard: https://shiny.tbep.org/tidalcreek-dash/

Tidal creeks (aka tributaries) are essential habitats in the Tampa Bay Estuary and serve as important focal points for understanding watershed inputs that affect water quality. A fundamental goal of the Tampa Bay Estuary Program is to develop effective nutrient management strategies to support the ecological function of tidal creeks. In partnership with Sarasota Bay NEP, Coastal & Heartland NEP, and local government and agency stakeholders, an open science framework has been developed for assessing the tidal creek condition based on a host of commonly collected water quality data [@tbep0216;@tbep0220;@tbep1121]. These assessments can support tracking of water quality management goals and can help refine restoration and management plans in priority tributaries, including those in need of hydrologic restoration that can support critical nursery habitats for sportfishes.

The tbeptools package includes a [simple features](https://r-spatial.github.io/sf/articles/sf1.html) spatial data object of the population of tidal creeks in southwest Florida, called `tidalcreeks()`. This includes `r nrow(tidalcreeks)` polyline features designated by a water body ID (`WBID`), creek id (`JEI`), and [FDEP class](https://floridadep.gov/dear/water-quality-standards/content/surface-water-quality-standards-classes-uses-criteria) (`class`, 1 for potable water, 2 for shellfish harvesting or propagation, and 3F/3M for freshwater/marine fish consumption, recreation, propagation and maintenance of a healthy, well-balanced population of fish and wildlife).

```{r, eval = T, cache = F}
mapview(tidalcreeks, homebutton = F, legend = F)
```

The tidal creek assessment framework was established based on data from the FDEP [Impaired Waters Rule](https://www.flrules.org/gateway/ChapterHome.asp?Chapter=62-303) database run 56 available [here](http://publicfiles.dep.state.fl.us/DEAR/IWR/) which includes data collected through January 10th 2019. However, the this framework intends to link to future IWR databases to refresh the site with new data as it becomes available.  Raw data from the IWR database required for assessment is provided in the tbeptools package in the `iwrraw()` data object.  

## Assessment

The tidal creek assessment framework includes both a "report card" and "indicators" assessment which are provided as separate tabs in the dashboard. The report card provides an assessment of total nitrogen concentrations (the limiting nutrient in these creeks) based on annual geometric average concentrations relative to standards developed for contributing freshwater streams. The indicators are based a several water quality metrics derived as outcomes of our study to describe tidal creek condition and provide insights into site specific attributes of these creeks that may govern overall creek health.

# Report Card

The report card is similar to the TBEP [water quality report card](https://shiny.tbep.org/wq-dash/) in that tidal creeks are assigned to categories within an assessment framework intended to serve as both a mechanism for evaluating data relative to the need for management action, and to identify stewardship goals that, if properly pursued, may preclude the need for any regulatory actions. These categories were established based principally on fish as a biological response indicator. Tidal creeks are assigned to one of five categories:

<span style="color:#ADD8E6; text-shadow: 0 0 3px #333;">__No Data__</span>: Data are unavailable for evaluation.

<span style="color:#33FF3B; text-shadow: 0 0 3px #333;">__Monitor__</span>: Creek is at or below nitrogen concentrations that protect individual creek types within the larger population of creeks.

<span style="color:#F9FF33; text-shadow: 0 0 3px #333;">__Caution__</span>: Creek nutrients showing signs of elevated nutrient concentrations that may increase risk of eutrophic condition.

<span style="color:#FFA500; text-shadow: 0 0 3px #333;">__Investigate__</span>: Creek nutrient concentrations above margin of safety to protect creek from potential impairment.

<span style="color:#FF7F50; text-shadow: 0 0 3px #333;">__Prioritize__</span>: Creek nutrient concentrations have exceeded regulatory standard for associated freshwater portion of tributary indicating that actions are needed to identify remediative measures to reduce nutrients to the creek.

Conceptually, these thresholds appear in the figure below.

```{r, echo = F, fig.cap = 'Scoring rubrik for tidal creeks based on nitrogen thresholds.', out.width = '80%'}
knitr::include_graphics('tidalcreekreport.PNG')
```

The Prioritize category was defined based on Florida's freshwater stream numeric nutrient criteria (NNC).Two different freshwater stream NNC are applicable to our region; the West Central NNC of 1.65 mg/l and  Peninsular region NNC of 1.54 mg/l. The histograms in the above figure represent a range of annual geometric mean (AGM) nitrogen concentrations associated with the Prioritize and Investigate categories which are based on the NNC. In the example above, the maximum expected distribution of AGMs not to exceed of 1.65 mg/l with a 1:3 exceedence probability (as defined in F.A.C. 62-303) was generated using monte carlo simulation and the highest observed standard deviation from data collected during the first creeks study. The Investigate category was then defined as an explicit margin of safety by adjusting the distribution to find the grand geometric average that would result in a 1:20 chance of exceeding 1.65 mg/l. Assignment of a creek into the Caution category depended on a creek length adjustment as described below to protect smaller creeks from elevated nutrient concentrations. 

The `tidaltargets()` data object included in tbeptools includes these thresholds.  Note that the "Caution" category is a function of creek length. 

```{r}
tidaltargets
```

A scoring algorithm was derived to define the final report card outcome for each creek using the entire ten year record of available data based on the following criteria. A single exceedance of the Prioritize and Investigate categories in any year of the ten year record would result in a classification of that creek into the respective category unless at least three other years of data were below the threshold level for that category. Creeks were assigned the next lower category if only one AGM for TN was above a given level while multiple other years (i.e., more than two) were below the given levels defining the cutoff points for each category. For example, a creek with at least 4 years of data and only a single exceedance of the Prioritize threshold would be assigned the Investigate category. Outcomes are exemplified below.

## Report Card Functions

The two primary functions for the tidal creek assessments are `anlz_tdlcrk()` to obtain the scores and `show_tdlcrk()` to view an interactive map of the results.  The `anlz_tdlcrk()` function uses the included `tidalcreeks()` and `iwrraw()` datasets to estimate the scores: 

```{r}
results <- anlz_tdlcrk(tidalcreeks, iwrraw)
results
```

The results include a unique creek identifier (`id`, based on the `wbid` and `JEI` fields), the waterbody id (`wbid`), the creek ID (`JEI`), the FDEP class (`class`), and results from the assessment in the remaining columns.  The columns `monitor`, `caution`, `investigate`, and `prioritize` indicate the number of years from 2011 to 2021 that the nitrogen values were within the ranges appropriate for the creek type as specified within `tidaltargets()`.  The `score` column indicates the overall category assigned to the creek for the period of record.  Note that many creeks are assigned a `No Data` value if sufficient data were unavailable.  A summation of the four component columns (`monitor`, `caution`, `investigate`, and `prioritize`) provides the number of years for which data were available at a creek. 

The `show_tdlcrk()` function can be used with the output of `anlz_tdlcrk()` to view an interactive map of the results. Creeks are color-coded by the exceedance categories, with "No Data" creeks shown in light blue.  

```{r, eval = T, cache = F}
show_tdlcrk(results)
```

A report card style matrix can be plotted using the `show_tdlcrkmatrix()` function that shows the overall creek score and the number of years of data that were used to estimate the overall score. The plot shows a matrix with rows for individual creeks and columns for overall creek score.  The columns show an overall creek score and the number of years in the prior ten years that nitrogen values at a creek were assigned to each of the four score categories.  Number of years is mapped to cell transparency.  By default, the plot shows creeks with a marine WBID (water body identifier) designation as `3M` or `2`.  This can be changed with the `class` argument (i.e., `class c('3M', '2', '3F', '1')` for marine and freshwater WBIDs). 

```{r, fig.height = 10, fig.width = 6}
show_tdlcrkmatrix(results)
```

## Indicator Functions
Water quality Indicators were developed to provide context for interpreting the report card outcomes as described in detail in Wessel et al. 2021 and include thresholds for total nitrogen (>1.1 mg/l), chlorophyll a (>11 ug/l), dissolved oxygen (< 42 % saturation), a trophic state index score (>55), the chlorophyll/nitrogen ratio (>15) and a ratio of the nitrates in the tidal and freshwater portion of the creek (>1) (if data are available).  The results for each indicator relative to the established thresholds are calculated on an annual basis and then synthesized for the 10 year period by calculating the percentage of annual outcomes exceeding the identified threshold indicator values out of the total number of years with available data. An integrative summary for all indicators is presented using a standardized polar coordinate system and Radar Charts to provide a single multi-metric summary plot of the results across indicators.  

The `anlz_tdlcrkindic()` function generates these annual outcomes for each `wbid`/`JEI` combination.

```{r}
results <- anlz_tdlcrkindic(tidalcreeks, iwrraw)
head(results)
```

Individual creek indicators are summarized using a multivariate response plot called a "radar plot" that indicates the percentage of years where each indicator exceeded its respective threshold value. These plots are created by using the `radar = TRUE` argument with `anlz_tdlcrkindic()` function and then using those results with the `show_tdlcrkradar()` function. The radar plots only apply to the marine WBIDs of the tidal creeks (Florida DEP class 2, 3M). Indicators without data for the creek do not have a point on the plot.

```{r, fig.width = 5, fig.height = 5}
cntdat <- anlz_tdlcrkindic(tidalcreeks, iwrraw, yr = 2021, radar = T)
show_tdlcrkradar(id = 495, cntdat = cntdat)
```


General descriptive plots of the annual outcomes are provided with interactive [plotly](https://plot.ly/r/) graphics using the `show_tdlcrkindic()` and `show_tdlcrkindiccdf()` functions. 

The `show_tdlcrindic()` function produces bar plots of annual outcomes at the selected creek. The creek to plot is selected with the `id` argument as an integer that is used to filter results from the `anlz_tdlcrkindic()` function, where the latter is passed to the `cntdat` argument.  The `thrsel` argument plots dotted red lines based on the threshold values.  Each year has its own unique color.  

```{r}
cntdat <- anlz_tdlcrkindic(tidalcreeks, iwrraw, yr = 2021)
show_tdlcrkindic(id = 495, cntdat = cntdat, thrsel = TRUE)
```

The `show_tdlcrkindiccdf()` function is similar except that empirical cumulative distribution functions (CDF) are plotted to evaluate outcomes for a specific creek relative to the entire distribution of creeks in southwest Florida. Each indicator and each year for the selected creek are plotted on the CDF curves.  Location of the points indicate both a comparison to the population and the trajectory of indicators over time (i.e., brown are older observations and blue are more recent).  Holding the mouse cursor over a point shows the year and holding the cursor over the line shows the percentile value from the CDF. 

```{r}
show_tdlcrkindiccdf(id = 495, cntdat = cntdat, thrsel = TRUE)
```


# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_reactable.R
\name{show_reactable}
\alias{show_reactable}
\title{Create reactable table from matrix data}
\usage{
show_reactable(totab, colfun, nrows = 10)
}
\arguments{
\item{totab}{A data frame in wide format of summarized results}

\item{colfun}{Function specifying how colors are treated in cell background}

\item{nrows}{numeric specifying number of rows in the table}
}
\value{
A \code{\link[reactable]{reactable}} table
}
\description{
Create reactable table from matrix data
}
\details{
This function is used internally within \code{\link{show_matrix}} and \code{\link{show_wqmatrix}}
}
\examples{
data(targets)
data(epcdata)

library(tidyr)
library(dplyr)

# data
totab <- anlz_avedat(epcdata) \%>\%
  .$ann \%>\%
  filter(var \%in\% 'mean_chla') \%>\%
  left_join(targets, by = 'bay_segment') \%>\%
  select(bay_segment, yr, val, chla_thresh) \%>\%
  mutate(
    bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB')),
     outcome = case_when(
      val < chla_thresh ~ 'green',
      val >= chla_thresh ~ 'red'
    )
  ) \%>\%
  select(bay_segment, yr, outcome) \%>\%
  spread(bay_segment, outcome)

# color function
colfun <- function(x){

  out <- case_when(
    x == 'red' ~ '#FF3333',
    x == 'green' ~ '#33FF3B'
  )

  return(out)

}

show_reactable(totab, colfun)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_thrplot.R
\name{show_thrplot}
\alias{show_thrplot}
\title{Plot annual water quality values, targets, and thresholds for a segment}
\usage{
show_thrplot(
  epcdata,
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  thr = c("chla", "la"),
  trgs = NULL,
  yrrng = c(1975, 2021),
  family = NA,
  labelexp = TRUE,
  txtlab = TRUE,
  thrs = FALSE,
  partialyr = FALSE
)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{bay_segment}{chr string for the bay segment, one of "OTB", "HB", "MTB", "LTB"}

\item{thr}{chr string indicating which water quality value and appropriate target/threshold to plot, one of "chl" for chlorophyll and "la" for light availability}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets/thresholds, defaults to \code{\link{targets}}}

\item{yrrng}{numeric vector indicating min, max years to include}

\item{family}{optional chr string indicating font family for text labels}

\item{labelexp}{logical indicating if y axis and target labels are plotted as expressions, default \code{TRUE}}

\item{txtlab}{logical indicating if a text label for the target value is shown in the plot}

\item{thrs}{logical indicating if reference lines are shown only for the regulatory threshold}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}
}
\value{
A \code{\link[ggplot2]{ggplot}} object
}
\description{
Plot annual water quality values, targets, and thresholds for a bay segment
}
\examples{
show_thrplot(epcdata, bay_segment = 'OTB', thr = 'chl')
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_tdlcrk.R
\name{anlz_tdlcrk}
\alias{anlz_tdlcrk}
\title{Estimate tidal creek report card scores}
\usage{
anlz_tdlcrk(tidalcreeks, iwrraw, tidtrgs = NULL, yr = 2021)
}
\arguments{
\item{tidalcreeks}{\code{\link[sf]{sf}} object for population of tidal creeks}

\item{iwrraw}{FDEP impaired waters data base as \code{\link{data.frame}}}

\item{tidtrgs}{optional \code{data.frame} for tidal creek nitrogen targets, defaults to \code{\link{tidaltargets}}}

\item{yr}{numeric for reference year to evaluate, scores are based on the planning period beginning ten years prior to this date}
}
\value{
A \code{\link{data.frame}} with the report card scores for each creek, as prioritize, investigate, caution, monitor, or no data
}
\description{
Estimate tidal creek report card scores
}
\examples{
anlz_tdlcrk(tidalcreeks, iwrraw, yr = 2021)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_matrixplotly.R
\name{show_matrixplotly}
\alias{show_matrixplotly}
\title{Creates a plotly matrix from any matrix function input}
\usage{
show_matrixplotly(
  mat,
  family = NA,
  tooltip = "Result",
  width = NULL,
  height = NULL
)
}
\arguments{
\item{mat}{input matrix as output from \code{\link{show_matrix}}, \code{\link{show_segmatrix}}, \code{\link{show_wqmatrix}}, or \code{\link{show_tbnimatrix}}}

\item{family}{optional chr string indicating font family for text labels}

\item{tooltip}{chr string indicating the column name for tooltip}

\item{width}{numeric for width of the plot in pixels}

\item{height}{numeric for height of the plot in pixels}
}
\value{
A \code{\link[plotly]{plotly}} data object
}
\description{
Creates a plotly matrix from any matrix function input
}
\examples{
mat <- show_wqmatrix(epcdata)
show_matrixplotly(mat)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_formfim.R
\name{read_formfim}
\alias{read_formfim}
\title{Format FIM data for the Tampa Bay Nekton Index}
\usage{
read_formfim(datin, locs = FALSE)
}
\arguments{
\item{datin}{input \code{data.frame} loaded from \code{\link{read_importfim}}}

\item{locs}{logical indicating if a spatial features object is returned with locations of each FIM sampling station}
}
\value{
A formatted \code{data.frame} with FIM data if \code{locs = FALSE}, otherwise a simple features object if \code{locs = TRUE}
}
\description{
Format FIM data for the Tampa Bay Nekton Index
}
\details{
Function is used internally within \code{\link{read_importfim}}
}
\examples{
\dontrun{
# file path
csv <- '~/Desktop/fimraw.csv'

# load and assign to object
fimdata <- read_importfim(csv)
}
}
\seealso{
\code{\link{read_importfim}}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tdlcrkradar.R
\name{show_tdlcrkradar}
\alias{show_tdlcrkradar}
\title{Radar plots for tidal creek indicators}
\usage{
show_tdlcrkradar(
  id,
  cntdat,
  col = "#338080E6",
  ptsz = 1,
  lbsz = 0.8,
  valsz = 1,
  brdwd = 5
)
}
\arguments{
\item{id}{numeric indicating the \code{id} number of the tidal creek to plot}

\item{cntdat}{output from \code{\link{anlz_tdlcrkindic}}}

\item{col}{color input for polygon and line portions}

\item{ptsz}{numeric size of points}

\item{lbsz}{numeric for size of text labels}

\item{valsz}{numeric for size of numeric value labels}

\item{brdwd}{numeric for polygon border width}
}
\value{
A radar plot
}
\description{
Radar plots for tidal creek indicators
}
\details{
See details in \code{\link{anlz_tdlcrkindic}} for an explanation of the indicators

Internal code borrowed heavily from \code{\link[fmsb]{radarchart}}.
}
\examples{
cntdat <- anlz_tdlcrkindic(tidalcreeks, iwrraw, yr = 2021, radar = TRUE)
show_tdlcrkradar(495, cntdat)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_transect.R
\name{read_transect}
\alias{read_transect}
\title{Import JSON seagrass transect data from Water Atlas}
\usage{
read_transect(training = FALSE, raw = FALSE)
}
\arguments{
\item{training}{logical if training data are imported or the complete database}

\item{raw}{logical indicating if raw, unformatted data are returned, see details}
}
\value{
data frame
}
\description{
Import JSON seagrass transect data from Water Atlas
}
\details{
The function imports a JSON file from the USF Water Atlas.  If \code{training = TRUE}, a dataset from the TBEP training survey is imported from \url{http://dev.seagrass.wateratlas.usf.edu/api/assessments/training}.  If \code{training = FALSE}, the entire transect survey database is imported from \url{http://dev.seagrass.wateratlas.usf.edu/api/assessments/all__use-with-care}.

Abundance is reported as a numeric value from 0 -5 for Braun-Blanquet coverage estimates, blade length is in cm, and short shoot density is number of shoots per square meter.  The short density is corrected for quadrat size included in the raw data.

If \code{raw = TRUE}, the unformatted data are returned.  The default is to use formatting that allows the raw data to be used with the downstream functions. The raw data may have extra information that may be of use outside of the plotting functions in this package.
}
\examples{
\dontrun{
# get training data
transect <- read_transect(training = TRUE)

# import all transect data
transect <- read_transect()
}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_yrattain.R
\name{anlz_yrattain}
\alias{anlz_yrattain}
\title{Get attainment categories for a chosen year}
\usage{
anlz_yrattain(epcdata, yrsel, partialyr = FALSE)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{yrsel}{numeric indicating chosen year}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}
}
\value{
A \code{data.frame} for the chosen year and all bay segments showing the bay segment averages for chloropyll concentration, light attenuations, segment targets, and attainment categories.
}
\description{
Get attainment categories for a chosen year
}
\examples{

# defaults to current year
anlz_yrattain(epcdata, yrsel = 2021)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_formwq.R
\name{read_formwq}
\alias{read_formwq}
\title{Format water quality data}
\usage{
read_formwq(datin)
}
\arguments{
\item{datin}{input \code{data.frame} loaded from \code{\link{read_importwq}}}
}
\value{
A lightly formatted \code{data.frame} with chloropyll and secchi observations
}
\description{
Format water quality data
}
\details{
Secchi data VOB depths or secchis < 0.5 ft from bottom are assigned \code{NA}, function is used internally within \code{\link{read_importwq}}
}
\examples{
\dontrun{
# file path
xlsx <- '~/Desktop/2018_Results_Updated.xls'

# load and assign to object
epcdata <- read_importwq(xlsx, download_latest = T)
}
}
\seealso{
\code{\link{read_importwq}}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_tbnimet.R
\name{anlz_tbnimet}
\alias{anlz_tbnimet}
\title{Get all raw metrics for Tampa Bay Nekton Index}
\usage{
anlz_tbnimet(fimdata, all = FALSE)
}
\arguments{
\item{fimdata}{\code{data.frame} formatted from \code{read_importfim}}

\item{all}{logical indicating if only TBNI metrics are returned (default), otherwise all are calcualted}
}
\value{
A data frame of raw metrics in wide fomat.  If \code{all = TRUE}, all metrics are returned, otherwise only \code{NumTaxa}, \code{BenthicTaxa}, \code{TaxaSelect}, \code{NumGuilds}, and \code{Shannon} are returned.
}
\description{
Get all raw metrics for Tampa Bay Nekton Index
}
\details{
All raw metrics are returned in addition to those required for the TBNI.  Each row shows metric values for a station, year, and month where fish catch data were available.
}
\examples{
anlz_tbnimet(fimdata)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_formtransect.R
\name{read_formtransect}
\alias{read_formtransect}
\title{Format seagrass transect data from Water Atlas}
\usage{
read_formtransect(jsn, training = FALSE, raw = FALSE)
}
\arguments{
\item{jsn}{A data frame returned from \code{\link[jsonlite]{fromJSON}}}

\item{training}{logical if input are transect training data or complete database}

\item{raw}{logical indicating if raw, unformatted data are returned, see details}
}
\value{
data frame in long format
}
\description{
Format seagrass transect data from Water Atlas
}
\details{
Shoot density is reported as number of shoots per square meter and is corrected for the quadrat size entered in the raw data.  Shoot density and blade height (cm) are based on averages across random observations at each transect point that are entered separately in the data form. Abundance is reported as a numeric value from 0 - 5 for Braun-Blanquet coverage estimates.

If \code{raw = TRUE}, the unformatted data are returned.  The default is to use formatting that allows the raw data to be used with the downstream functions. The raw data may have extra information that may be of use outside of the plotting functions in this package.
}
\examples{
library(jsonlite)

\dontrun{
# all transect data
url <- 'http://dev.seagrass.wateratlas.usf.edu/api/assessments/all__use-with-care'
jsn <- fromJSON(url)
trndat <- read_formtransect(jsn)
}

# training transect data
url <- 'http://dev.seagrass.wateratlas.usf.edu/api/assessments/training'
jsn <- fromJSON(url)
trndat <- read_formtransect(jsn, training = TRUE)
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbnispp.R
\docType{data}
\name{tbnispp}
\alias{tbnispp}
\title{Reference table for Tampa Bay Nekton Index species classifications}
\format{
A data frame with 196 rows and 9 variables:
\describe{
  \item{NODCCODE}{chr}
  \item{ScientificName}{chr}
  \item{Include_TB_Index}{chr}
  \item{Hab_Cat}{chr}
  \item{Est_Cat}{chr}
  \item{Est_Use}{chr}
  \item{Feeding_Cat}{chr}
  \item{Feeding_Guild}{chr}
  \item{Selected_Taxa}{chr}
}
}
\usage{
tbnispp
}
\description{
Reference table for Tampa Bay Nekton Index species classifications
}
\examples{
\dontrun{
library(dplyr)

# import and clean
tbnispp <- read.csv('../tbni-proc/data/TBIndex_spp_codes.csv',
    header = TRUE, stringsAsFactors = FALSE) \%>\%
  mutate(
    NODCCODE = as.character(NODCCODE),
    NODCCODE = case_when(NODCCODE == "9.998e+09" ~ "9998000000",
                             TRUE ~ NODCCODE)
  )

save(tbnispp, file = 'data/tbnispp.RData', compress = 'xz')
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_wqmatrix.R
\name{show_wqmatrix}
\alias{show_wqmatrix}
\title{Create a colorized table for chlorophyll or light attenuation exceedances}
\usage{
show_wqmatrix(
  epcdata,
  param = c("chla", "la"),
  txtsz = 3,
  trgs = NULL,
  yrrng = c(1975, 2021),
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  asreact = FALSE,
  nrows = 10,
  abbrev = FALSE,
  family = NA,
  plotly = FALSE,
  partialyr = FALSE,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{param}{chr string for which parameter to plot, one of \code{"chla"} for chlorophyll or \code{"la"} for light attenuation}

\item{txtsz}{numeric for size of text in the plot, applies only if \code{tab = FALSE}}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets, defaults to \code{\link{targets}}}

\item{yrrng}{numeric vector indicating min, max years to include}

\item{bay_segment}{chr string for bay segments to include, one to all of "OTB", "HB", "MTB", "LTB"}

\item{asreact}{logical indicating if a \code{\link[reactable]{reactable}} object is returned}

\item{nrows}{if \code{asreact = TRUE}, a numeric specifying number of rows in the table}

\item{abbrev}{logical indicating if text labels in the plot are abbreviated as the first letter}

\item{family}{optional chr string indicating font family for text labels}

\item{plotly}{logical if matrix is created using plotly}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}
}
\value{
A static \code{\link[ggplot2]{ggplot}} object is returned if \code{asreact = FALSE}, otherwise a \code{\link[reactable]{reactable}} table is returned
}
\description{
Create a colorized table for chlorophyll or light attenuation exceedances
}
\examples{
show_wqmatrix(epcdata)
}
\seealso{
\code{\link{show_matrix}}, \code{\link{show_segmatrix}}
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fimdata.R
\docType{data}
\name{fimdata}
\alias{fimdata}
\title{FIM data for Tampa Bay Nekton Index current as of 07262021}
\format{
A data frame with 45945 rows and 19 variables:
\describe{
  \item{Reference}{chr}
  \item{Sampling_Date}{Date}
  \item{Latitude}{num}
  \item{Longitude}{num}
  \item{Zone}{chr}
  \item{Grid}{int}
  \item{NODCCODE}{chr}
  \item{Year}{num}
  \item{Month}{num}
  \item{Total_N}{num}
  \item{ScientificName}{chr}
  \item{Include_TB_Index}{chr}
  \item{Hab_Cat}{chr}
  \item{Est_Cat}{chr}
  \item{Est_Use}{chr}
  \item{Feeding_Cat}{chr}
  \item{Feeding_Guild}{chr}
  \item{Selected_Taxa}{chr}
  \item{bay_segment}{chr}
  }
}
\usage{
fimdata
}
\description{
FIM data for Tampa Bay Nekton Index current as of 07262021
}
\examples{
\dontrun{
csv <- '~/Desktop/fimdata.csv'

fimdata <- read_importfim(csv, download_latest = TRUE)

save(fimdata, file = 'data/fimdata.RData', compress = 'xz')

}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidaltargets.R
\docType{data}
\name{tidaltargets}
\alias{tidaltargets}
\title{Tidal creek nitrogen targets}
\format{
A data frame with 2 rows and 4 variables:
\describe{
  \item{region}{chr}
  \item{prioritize}{num}
  \item{investigate}{num}
  \item{caution}{num}
}
}
\usage{
tidaltargets
}
\description{
Tidal creek nitrogen targets
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tdlcrkmatrix.R
\name{show_tdlcrkmatrix}
\alias{show_tdlcrkmatrix}
\title{Plot the tidal creek report card matrix}
\usage{
show_tdlcrkmatrix(
  dat,
  class = c("3M", "2"),
  score = c("Prioritize", "Investigate", "Caution", "Monitor"),
  family = NA,
  size = 11
)
}
\arguments{
\item{dat}{input creek score data returned from \code{\link{anlz_tdlcrk}}}

\item{class}{character vector indicating which creek classes to show, one to many of \code{'3M'}, \code{'2'}, \code{'3F'}, and \code{'1'}.  Defaults to marine only (\code{'3M', '2'}).}

\item{score}{character vector of score categories to include, one to many of \code{'Prioritize'}, \code{'Investigate'}, \code{'Caution'}, and \code{'Monitor'}. Defaults to all.}

\item{family}{optional chr string indicating font family for text labels}

\item{size}{numeric for text and line scaling}
}
\value{
A static \code{\link[ggplot2]{ggplot}} object is returned.
}
\description{
Plot the tidal creek report card matrix
}
\details{
The plot shows a matrix with rows for individual creeks and columns for overall creek score.  The columns show an overall creek score and the number of years in the prior ten years that nitrogen values at a creek were assigned to each of the four score categories.  Number of years is mapped to cell transparency.
}
\examples{
dat <- anlz_tdlcrk(tidalcreeks, iwrraw, yr = 2021)
show_tdlcrkmatrix(dat)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/targets.R
\docType{data}
\name{targets}
\alias{targets}
\title{Bay segment targets}
\format{
A data frame with 4 rows and 8 variables:
\describe{
  \item{bay_segment}{chr}
  \item{name}{chr}
  \item{chla_target}{num}
  \item{chla_smallex}{num}
  \item{chla_thresh}{num}
  \item{la_target}{num}
  \item{la_smallex}{num}
  \item{la_thresh}{num}
}
}
\usage{
targets
}
\description{
Bay segment specific management targets including low and high magnitude exceedance thresholds
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_transectsum.R
\name{show_transectsum}
\alias{show_transectsum}
\title{Plot frequency occurrence for a seagrass transect by time for all species}
\usage{
show_transectsum(
  transectocc,
  site,
  species = c("Halodule", "Syringodium", "Thalassia", "Halophila", "Ruppia",
    "Caulerpa"),
  yrrng = c(1998, 2021),
  abund = FALSE,
  sppcol = NULL
)
}
\arguments{
\item{transectocc}{data frame returned by \code{\link{anlz_transectocc}}}

\item{site}{chr string indicating site results to plot}

\item{species}{chr string indicating which species to plot}

\item{yrrng}{numeric indicating year ranges to evaluate}

\item{abund}{logical indicating if abundance averages are plotted instead of frequency occurrence}

\item{sppcol}{character vector of alternative colors to use for each species, must have length of six}
}
\value{
A \code{\link[plotly]{plotly}} object
}
\description{
Plot frequency occurrence for a seagrass transect by time for all species
}
\details{
This plot provides a quick visual assessment of how frequency occurrence or abundance for multiple species has changed over time at a selected transect.
}
\examples{
\dontrun{
transect <- read_transect()
}
transectocc <- anlz_transectocc(transect)
show_transectsum(transectocc, site = 'S3T10')
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_transectmatrix.R
\name{show_transectmatrix}
\alias{show_transectmatrix}
\title{Show matrix of seagrass frequency occurrence by bay segments and year}
\usage{
show_transectmatrix(
  transectocc,
  bay_segment = c("OTB", "HB", "MTB", "LTB", "BCB"),
  total = TRUE,
  neutral = FALSE,
  yrrng = c(1998, 2021),
  alph = 1,
  txtsz = 3,
  family = NA,
  rev = FALSE,
  position = "top",
  plotly = FALSE,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{transectocc}{data frame returned by \code{\link{anlz_transectocc}}}

\item{bay_segment}{chr string for the bay segment, one to many of "HB", "OTB", "MTB", "LTB", "TCB", "MR", "BCB"}

\item{total}{logical indicating if average frequency occurrence is calculated for the entire bay across segments}

\item{neutral}{logical indicating if a neutral and continuous color scheme is used}

\item{yrrng}{numeric indicating year ranges to evaluate}

\item{alph}{numeric indicating alpha value for score category colors}

\item{txtsz}{numeric for size of text in the plot}

\item{family}{optional chr string indicating font family for text labels}

\item{rev}{logical if factor levels for bay segments are reversed}

\item{position}{chr string of location for bay segment labels, default on top, passed to \code{\link[ggplot2]{scale_x_discrete}}}

\item{plotly}{logical if matrix is created using plotly}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}
}
\value{
A \code{\link[ggplot2]{ggplot}} object showing trends over time for each bay segment if \code{plotly = FALSE}, otherwise a \code{\link[plotly]{plotly}} object
}
\description{
Show matrix of seagrass frequency occurrence by bay segments and year
}
\details{
Results are based on averages across species by date and transect in each bay segment

The color scheme is based on arbitrary breaks at 25, 50, and 75 percent frequency occurrence.  These don't necessarily translate to any ecological breakpoints.  Use \code{neutral = TRUE} to use a neutral and continuous color palette.
}
\examples{
\dontrun{
transect <- read_transect()
}
transectocc <- anlz_transectocc(transect)
show_transectmatrix(transectocc)
}
\references{
This plot is a representation of Table 1 in R. Johansson (2016) Seagrass Transect Monitoring in Tampa Bay: A Summary of Findings from 1997 through 2015, Technical report #08-16, Tampa Bay Estuary Program, St. Petersburg, Florida.
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_segmatrix.R
\name{show_segmatrix}
\alias{show_segmatrix}
\title{Create a colorized table for water quality outcomes and exceedances by segment}
\usage{
show_segmatrix(
  epcdata,
  txtsz = 3,
  trgs = NULL,
  yrrng = c(1975, 2021),
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  abbrev = FALSE,
  family = NA,
  historic = FALSE,
  plotly = FALSE,
  partialyr = FALSE,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{txtsz}{numeric for size of text in the plot, applies only if \code{tab = FALSE}}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets, defaults to \code{\link{targets}}}

\item{yrrng}{numeric vector indicating min, max years to include}

\item{bay_segment}{chr string for bay segments to include, only one of "OTB", "HB", "MTB", "LTB"}

\item{abbrev}{logical indicating if text labels in the plot are abbreviated as the first letter}

\item{family}{optional chr string indicating font family for text labels}

\item{historic}{logical if historic data are used from 2005 and earlier}

\item{plotly}{logical if matrix is created using plotly}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}
}
\value{
A static \code{\link[ggplot2]{ggplot}} object is returned
}
\description{
Create a colorized table for water quality outcomes by segment that includes the management action and chlorophyll, and light attenuation exceedances
}
\details{
This function provides a combined output for the \code{\link{show_wqmatrix}} and \code{\link{show_matrix}} functions. Only one bay segment can be plotted for each function call.
}
\examples{
show_segmatrix(epcdata, bay_segment = 'OTB')
}
\seealso{
\code{\link{show_wqmatrix}}, \code{\link{show_matrix}}
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_attain.R
\name{anlz_attain}
\alias{anlz_attain}
\title{Get attainment categories}
\usage{
anlz_attain(avedat, magdurout = FALSE, trgs = NULL)
}
\arguments{
\item{avedat}{result returned from \code{\link{anlz_avedat}}}

\item{magdurout}{logical indicating if the separate magnitude and duration estimates are returned}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets, defaults to \code{\link{targets}}}
}
\value{
A \code{data.frame} for each year and bay segment showing the attainment category
}
\description{
Get attainment categories for each year and bay segment using chlorophyll and light attenuation
}
\examples{
avedat <- anlz_avedat(epcdata)
anlz_attain(avedat)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_transectavespp.R
\name{show_transectavespp}
\alias{show_transectavespp}
\title{Show annual averages of seagrass frequency occurrence by bay segments, year, and species}
\usage{
show_transectavespp(
  transectocc,
  bay_segment = c("OTB", "HB", "MTB", "LTB", "BCB"),
  yrrng = c(1998, 2021),
  species = c("Halodule", "Syringodium", "Thalassia", "Halophila", "Ruppia",
    "Caulerpa"),
  total = TRUE,
  alph = 1,
  family = NA,
  plotly = FALSE,
  asreact = FALSE,
  width = NULL,
  height = NULL,
  sppcol = NULL
)
}
\arguments{
\item{transectocc}{data frame returned by \code{\link{anlz_transectocc}}}

\item{bay_segment}{chr string for the bay segment, one to many of "OTB", "HB", "MTB", "LTB", "BCB"}

\item{yrrng}{numeric indicating year ranges to evaluate}

\item{species}{chr string of species to summarize, one to many of "Halodule", "Syringodium", "Thalassia", "Ruppia", "Halophila", "Caulerpa"}

\item{total}{logical indicating if total frequency occurrence for all species is also returned, only applies if \code{asreact = FALSE}}

\item{alph}{numeric indicating alpha value for score category colors}

\item{family}{optional chr string indicating font family for text labels}

\item{plotly}{logical if matrix is created using plotly}

\item{asreact}{logical if a reactable table is returned instead of a plot}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{sppcol}{character vector of alternative colors to use for each species, must have length of six}
}
\value{
If \code{asreact = F}, a \code{\link[ggplot2]{ggplot}} or \code{\link[plotly]{plotly}} (if \code{plotly = T}) object is returned showing trends over time by species for selected bay segments.  If \code{asreact = T}, a \code{\link[reactable]{reactable}} table showing results by year, segment, and species is returned.
}
\description{
Show annual averages of seagrass frequency occurrence by bay segments, year, and species
}
\details{
Results are based on averages across species by date and transect in each bay segment
}
\examples{
\dontrun{
transect <- read_transect()
}
transectocc <- anlz_transectocc(transect)
show_transectavespp(transectocc)
}
\references{
The plot is a representation of figure 2 in Johansson, R. (2016) Seagrass Transect Monitoring in Tampa Bay: A Summary of Findings from 1997 through 2015, Technical report #08-16, Tampa Bay Estuary Program, St. Petersburg, Florida.

The table is a representation of table 2, p. 163 in Yarbro, L. A., and P. R. Carlson, Jr., eds. 2016. Seagrass Integrated Mapping and Monitoring Program: Mapping and Monitoring Report No. 2. Fish and Wildlife Research Institute Technical Report TR-17 version 2. vi + 281 p.
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_tdlcrkindic.R
\name{anlz_tdlcrkindic}
\alias{anlz_tdlcrkindic}
\title{Analyze tidal creek water quality indicators}
\usage{
anlz_tdlcrkindic(tidalcreeks, iwrraw, yr = 2021, radar = FALSE)
}
\arguments{
\item{tidalcreeks}{\code{\link[sf]{sf}} object for population of tidal creeks}

\item{iwrraw}{FDEP impaired waters rule data base as \code{\link{data.frame}}}

\item{yr}{numeric for reference year to evaluate, scores are based on the planning period beginning ten years prior to this date}

\item{radar}{logical indicating if output is for \code{\link{show_tdlcrkradar}}, see details}
}
\value{
A \code{\link{data.frame}} with the indicator values for each tidal creek
}
\description{
Estimate tidal creek water quality indicators to support report card scores
}
\details{
Annual geometric means for additional water quality data available at each wbid, JEI combination.  Florida trophic state index values are also estimated where data are available.

Nitrogen ratios are estimated for JEIs that cover source (upstream, freshwater) and tidal (downstream) WBIDs, defined as the ratio of concentrations between the two (i.e., ratios > 1 mean source has higher concentrations).  Nitrogen ratios for a given year reflect the ratio of the median nitrogen concentrations when they were measured in both a source and tidal segment during the same day.  Note that a ratio of one can be obtained if both the source and tidal segments are at minimum detection.

Indicators for years where more than 10\% of observations exceed DO saturation criteria are also estimated.  The \code{do_bnml} and \code{do_prop} columns show a 1 or 0 for a given year to indicate if more than ten percent of observations were below DO percent saturation of 42.  The first column is based on a binomial probability exceedance that takes into account the number of observations in the year and the second column is based on a simple proportional estimate from the raw data.

If \code{radar = TRUE}, output is returned in a format for use with \code{\link{show_tdlcrkradar}}  Specifically, results are calculated as the percentage of years where an indicator exceeds a relevant threshold.  This only applies to the marine WBIDs of the tidal creeks (Florida DEP class 2, 3M).  Six indicators are returned with percentage exceedances based on total nitrogen (\code{tn_ind}) greater than 1.1 mg/L, chlorophyll (\code{chla_ind}) greater than 11 ug/L, trophic state index (\code{tsi_ind}) greater than 55 (out of 100), nitrate/nitrite ratio between marine and upstream segments (\code{nox_ind}) greater than one, chlorophyll and total nitrogen ratios > 15, and percentage of years more where than ten percent of observations were below DO percent saturation of 42.
}
\examples{
dat <- anlz_tdlcrkindic(tidalcreeks, iwrraw, yr = 2021)
head(dat)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_transectave.R
\name{anlz_transectave}
\alias{anlz_transectave}
\title{Get annual averages of seagrass frequency occurrence by bay segments and year}
\usage{
anlz_transectave(
  transectocc,
  bay_segment = c("OTB", "HB", "MTB", "LTB", "BCB"),
  total = TRUE,
  yrrng = c(1998, 2021),
  rev = FALSE
)
}
\arguments{
\item{transectocc}{data frame returned by \code{\link{anlz_transectocc}}}

\item{bay_segment}{chr string for the bay segment, one to many of "OTB", "HB", "MTB", "LTB", "BCB"}

\item{total}{logical indicating if average frequency occurrence is calculated for the entire bay across segments}

\item{yrrng}{numeric indicating year ranges to evaluate}

\item{rev}{logical if factor levels for bay segments are reversed}
}
\value{
A data frame of annual averages by bay segment
}
\description{
Get annual averages of seagrass frequency occurrence by bay segments and year
}
\details{
The \code{focat} column returned in the results shows a color category based on arbitrary breaks of the frequency occurrence estimates (\code{foest}) at 25, 50, and 75 percent.  These don't necessarily translate to any ecological breakpoints.
}
\examples{
\dontrun{
transect <- read_transect()
}
transectocc <- anlz_transectocc(transect)
anlz_transectave(transectocc)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgseg.R
\docType{data}
\name{sgseg}
\alias{sgseg}
\title{Seagrass segment reporting boundaries for southwest Florida}
\format{
A simple features \code{\link[sf]{sf}} object (POLYGON) with 22 features and 1 field, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
\describe{
  \item{segment}{chr}
}
}
\usage{
sgseg
}
\description{
Seagrass segment reporting boundaries for southwest Florida
}
\details{
These polygons are used by Southwest Florida Water Management District for summarizing seagrass coverage estimates by major coastal and estuarine boundaries.
}
\examples{
\dontrun{
library(sf)
library(tidyverse)
library(tools)

# create sf object of boundaries
# make sure projection does not change
sgseg <- st_read(
  dsn = '~/Desktop/TBEP/GISboundaries/Seagrass_Segment_Boundaries/Seagrass_Segment_Boundaries.shp',
  drivers = 'ESRI Shapefile'
  ) \%>\%
  select(segment = SEAGRASSSE) \%>\%
  mutate(
     segment = tolower(segment),
     segment = case_when(
        segment == 'terra ciea bay' ~ 'Terra Ceia Bay',
        T ~ segment
     ),
     segment = toTitleCase(segment)
  )

# save
save(sgseg, file = 'data/sgseg.RData', compress = 'xz')

}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_avedatsite.R
\name{anlz_avedatsite}
\alias{anlz_avedatsite}
\title{Estimate annual means by site}
\usage{
anlz_avedatsite(epcdata, partialyr = FALSE)
}
\arguments{
\item{epcdata}{\code{data.frame} formatted from \code{read_importwq}}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}
}
\value{
Mean estimates for chlorophyll and secchi
}
\description{
Estimate annual means by site for chlorophyll and secchi data
}
\examples{
# view average estimates
anlz_avedatsite(epcdata)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbshed.R
\docType{data}
\name{tbshed}
\alias{tbshed}
\title{Spatial data object of Tampa Bay watershed}
\format{
A simple features \code{\link[sf]{sf}} object (POLYGON) with 1 feature and 0 fields, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
}
\usage{
tbshed
}
\description{
Spatial data object of Tampa Bay watershed, includes the bay proper
}
\examples{
library(sf)
plot(st_geometry(tbshed))
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_boxplot.R
\name{show_boxplot}
\alias{show_boxplot}
\title{Plot monthly chlorophyll or light attenuation values for a segment}
\usage{
show_boxplot(
  epcdata,
  param = c("chla", "la"),
  yrsel = NULL,
  yrrng = c(1975, 2021),
  ptsz = 0.5,
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  trgs = NULL,
  family = NA,
  labelexp = TRUE,
  txtlab = TRUE,
  partialyr = FALSE
)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{param}{chr string for which parameter to plot, one of \code{"chla"} for chlorophyll or \code{"la"} for light attenuation}

\item{yrsel}{numeric for year to emphasize, shown as separate red points on the plot}

\item{yrrng}{numeric vector indicating min, max years to include}

\item{ptsz}{numeric indicating point size of observations not in \code{yrsel}}

\item{bay_segment}{chr string for the bay segment, one of "OTB", "HB", "MTB", "LTB"}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets, defaults to \code{\link{targets}}}

\item{family}{optional chr string indicating font family for text labels}

\item{labelexp}{logical indicating if y axis and target labels are plotted as expressions, default \code{TRUE}}

\item{txtlab}{logical indicating if a text label for the target value is shown in the plot}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}
}
\value{
A \code{\link[ggplot2]{ggplot}} object
}
\description{
Plot monthly chlorophyll or light attenuation values for a bay segment
}
\details{
Points not included in \code{yrsel} are plotted over the box plots using \code{\link[ggplot2]{position_jitter}}. Use \code{ptsz = -1} to suppress.  The dotted line in the plot shows the large exceedance value.
}
\examples{
show_boxplot(epcdata, bay_segment = 'OTB')
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_matrix.R
\name{show_matrix}
\alias{show_matrix}
\title{Create a colorized table for indicator reporting}
\usage{
show_matrix(
  epcdata,
  txtsz = 3,
  trgs = NULL,
  yrrng = NULL,
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  asreact = FALSE,
  nrows = 10,
  abbrev = FALSE,
  family = NA,
  historic = FALSE,
  plotly = FALSE,
  partialyr = FALSE,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{txtsz}{numeric for size of text in the plot, applies only if \code{tab = FALSE}}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets, defaults to \code{\link{targets}}}

\item{yrrng}{numeric vector indicating min, max years to include, defaults to range of years in \code{epcdata}}

\item{bay_segment}{chr string for bay segments to include, one to all of "OTB", "HB", "MTB", "LTB"}

\item{asreact}{logical indicating if a \code{\link[reactable]{reactable}} object is returned}

\item{nrows}{if \code{asreact = TRUE}, a numeric specifying number of rows in the table}

\item{abbrev}{logical indicating if text labels in the plot are abbreviated as the first letter}

\item{family}{optional chr string indicating font family for text labels}

\item{historic}{logical if historic data are used from 2005 and earlier}

\item{plotly}{logical if matrix is created using plotly}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}
}
\value{
A static \code{\link[ggplot2]{ggplot}} object is returned if \code{asreact = FALSE}, otherwise a \code{\link[reactable]{reactable}} table is returned
}
\description{
Create a colorized table for indicator reporting
}
\examples{
show_matrix(epcdata)
}
\seealso{
\code{\link{show_wqmatrix}}, \code{\link{show_segmatrix}}
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_attainsite.R
\name{anlz_attainsite}
\alias{anlz_attainsite}
\title{Get site attainments}
\usage{
anlz_attainsite(
  avedatsite,
  thr = c("chla", "la"),
  trgs = NULL,
  yrrng = NULL,
  thrs = FALSE
)
}
\arguments{
\item{avedatsite}{result returned from \code{\link{anlz_avedatsite}}}

\item{thr}{chr string indicating with water quality value and appropriate threshold to to plot, one of "chl" for chlorophyll and "la" for light availability}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets, defaults to \code{\link{targets}}}

\item{yrrng}{optional numeric value for year to return, defaults to all}

\item{thrs}{logical indicating if attainment category is relative to targets (default) or thresholds}
}
\value{
a \code{data.frame} for each year and site showing the attainment category
}
\description{
Get site attainment categories for chlorophyll or light attenuation
}
\details{
This function is a simplication of the attainment categories returned by \code{\link{anlz_attain}}.  Sites are only compared to the targets/thresholds that apply separately for chlorophyll or light attenuation.
}
\examples{
avedatsite <- anlz_avedatsite(epcdata)
anlz_attainsite(avedatsite)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbniref.R
\docType{data}
\name{tbniref}
\alias{tbniref}
\title{Reference conditions for Tampa Bay Nekton Index metrics}
\format{
A data frame with 16 rows and 12 variables:
\describe{
  \item{bay_segment}{chr}
  \item{Season}{chr}
  \item{NumTaxa_P5}{num}
  \item{NumTaxa_P95}{num}
  \item{BenthicTaxa_P5}{num}
  \item{BenthicTaxa_P95}{num}
  \item{TaxaSelect_P5}{num}
  \item{TaxaSelect_P95}{num}
  \item{NumGuilds_P5}{num}
  \item{NumGuilds_P95}{num}
  \item{Shannon_P5}{num}
  \item{Shannon_P95}{num}
}
}
\usage{
tbniref
}
\description{
Reference conditions for Tampa Bay Nekton Index metrics
}
\examples{
\dontrun{

library(tbeptools)

tbniref <- anlz_tbnimet(fimdata) \%>\%
  dplyr::filter(between(Year, 1998, 2015)) \%>\%
  dplyr::select(Season, bay_segment, NumTaxa, BenthicTaxa, TaxaSelect, NumGuilds, Shannon) \%>\%
  dplyr::group_by(bay_segment, Season) \%>\%
  dplyr::summarize(NumTaxa_P5 = round(quantile(NumTaxa, probs = 0.05)),
                   NumTaxa_P95 = round(quantile(NumTaxa, probs = 0.95)),
                   BenthicTaxa_P5 = round(quantile(BenthicTaxa, probs = 0.05)),
                   BenthicTaxa_P95 = round(quantile(BenthicTaxa, probs = 0.95)),
                   TaxaSelect_P5 = round(quantile(TaxaSelect, probs = 0.05)),
                   TaxaSelect_P95 = round(quantile(TaxaSelect, probs = 0.95)),
                   NumGuilds_P5 = round(quantile(NumGuilds, probs = 0.05)),
                   NumGuilds_P95 = round(quantile(NumGuilds, probs = 0.95)),
                   Shannon_P5 = quantile(Shannon, probs = 0.05),
                   Shannon_P95 = quantile(Shannon, probs = 0.95))

save(tbniref, file = 'data/tbniref.RData', compress = 'xz')

}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_tbniave.R
\name{anlz_tbniave}
\alias{anlz_tbniave}
\title{Get annual averages of Tampa Bay Nekton Index scores by bay segment}
\usage{
anlz_tbniave(
  tbniscr,
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  rev = FALSE,
  perc = c(32, 46)
)
}
\arguments{
\item{tbniscr}{input data frame as returned by \code{\link{anlz_tbniscr}}}

\item{bay_segment}{chr string for the bay segment, one to many of "OTB", "HB", "MTB", "LTB"}

\item{rev}{logical if factor levels for bay segments are reversed}

\item{perc}{numeric values indicating break points for score categories}
}
\value{
A data frame of annual averages by bay segment
}
\description{
Get annual averages of Tampa Bay Nekton Index scores by bay segment
}
\examples{
tbniscr <- anlz_tbniscr(fimdata)
anlz_tbniave(tbniscr)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_tbniscr.R
\name{anlz_tbniscr}
\alias{anlz_tbniscr}
\title{Get Tampa Bay Nekton Index scores}
\usage{
anlz_tbniscr(fimdata, raw = TRUE)
}
\arguments{
\item{fimdata}{\code{data.frame} formatted from \code{\link{read_importfim}}}

\item{raw}{logical indicating if raw metric values are also returned}
}
\value{
A data frame of metrics and TBNI scores in wide format.
}
\description{
Get Tampa Bay Nekton Index scores
}
\details{
This function calculates raw and scored metrics for the TBNI, including \code{NumTaxa}, \code{BenthicTaxa}, \code{TaxaSelect}, \code{NumGuilds}, and \code{Shannon}.  The total TBNI score is returned as \code{TBNI_Score}.
}
\examples{
anlz_tbniscr(fimdata)
}
\seealso{
\code{\link{anlz_tbnimet}}
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_tbbimed.R
\name{anlz_tbbimed}
\alias{anlz_tbbimed}
\title{Get annual medians of Tampa Bay Benthic Index scores by bay segment}
\usage{
anlz_tbbimed(
  tbbiscr,
  bay_segment = c("HB", "OTB", "MTB", "LTB", "TCB", "MR", "BCB", "All", "All (wt)"),
  rev = FALSE,
  yrrng = c(1993, 2019)
)
}
\arguments{
\item{tbbiscr}{input data frame as returned by \code{\link{anlz_tbbiscr}}}

\item{bay_segment}{chr string for the bay segment, one to many of "HB", "OTB", "MTB", "LTB", "TCB", "MR", "BCB", "All", "All (wt)"}

\item{rev}{logical if factor levels for bay segments are reversed}

\item{yrrng}{numeric indicating year ranges to evaluate}
}
\value{
A data frame of annual medians by bay segment
}
\description{
Get annual medians of Tampa Bay Benthic Index scores by bay segment
}
\details{
Additional summaries are provided for the entire bay, as a summary across categories ("All") and a summary weighted across the relative sizes of each bay segment ("All (wt)").

Only sampling funded by TBEP and as part of the routine EPC benthic monitoring program are included in the final categories.
}
\examples{
tbbiscr <- anlz_tbbiscr(benthicdata)
anlz_tbbimed(tbbiscr)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_formphyto.R
\name{read_formphyto}
\alias{read_formphyto}
\title{Format phytoplankton data}
\usage{
read_formphyto(datin)
}
\arguments{
\item{datin}{input \code{data.frame} loaded from \code{\link{read_importphyto}}}
}
\value{
A formatted \code{data.frame} with phytoplankton count data
}
\description{
Format phytoplankton data
}
\details{
Only seven taxonomic groups are summarized. Pyrodinium bahamense, Karenia brevis, Tripos hircus, Pseudo-nitzschia sp., and Pseudo-nitzschia pungens are retained at the species level.  Bacillariophyta and Cyanobacteria are retained at the phylum level.  All other taxa are grouped into an "other" category.
}
\examples{
\dontrun{
# file path
xlsx <- '~/Desktop/phyto_data.xlsx'

# load and assign to object
phytodata <- read_importphyto(xlsx)
}
}
\seealso{
\code{\link{read_importphyto}}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tbniscrall.R
\name{show_tbniscrall}
\alias{show_tbniscrall}
\title{Plot Tampa Bay Nekton Index scores over time as average across bay segments}
\usage{
show_tbniscrall(
  tbniscr,
  perc = c(32, 46),
  alph = 0.3,
  ylim = c(22, 58),
  rev = FALSE,
  plotly = FALSE
)
}
\arguments{
\item{tbniscr}{input dat frame as returned by \code{\link{anlz_tbniscr}}}

\item{perc}{numeric values indicating break points for score categories}

\item{alph}{numeric indicating alpha value for score category colors}

\item{ylim}{numeric for y axis limits}

\item{rev}{logical if factor levels for bay segments are reversed}

\item{plotly}{logical if matrix is created using plotly}
}
\value{
A \code{\link[ggplot2]{ggplot}} object showing trends over time in TBNI scores for each bay segment or a \code{\link[plotly]{plotly}} object if \code{plotly = TRUE}
}
\description{
Plot Tampa Bay Nekton Index scores over time as average across bay segments
}
\examples{
tbniscr <- anlz_tbniscr(fimdata)
show_tbniscrall(tbniscr)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trnlns.R
\docType{data}
\name{trnlns}
\alias{trnlns}
\title{Seagrass transect locations}
\format{
A \code{sf} LINESTRING object
}
\usage{
trnlns
}
\description{
Seagrass transect locations
}
\examples{
\dontrun{
library(sf)
library(dplyr)

trnlns <- st_read('T:/05_GIS/SEAGRASS_TRANSECTS/transect_routes.shp') \%>\%
   st_transform(crs = 4326) \%>\%
   dplyr::filter(!as.character(Site) \%in\% c('S8T1', 'S8T2', 'S8T3', 'S3T2')) \%>\%
   dplyr::mutate_if(is.factor, as.character) \%>\%
   dplyr::filter(Site \%in\% trnpts$TRAN_ID)

# add bearing, positive counter-clockwise from east
bearing <- lapply(trnlns$geometry, function(x) geosphere::bearing(x[, c(1:2)])[[1]]) \%>\%
  unlist()

trnlns$bearing <- bearing

save(trnlns, file = 'data/trnlns.RData', compress = 'xz')
#' }
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_formbenthic.R
\name{read_formbenthic}
\alias{read_formbenthic}
\title{Format benthic data for the Tampa Bay Benthic Index}
\usage{
read_formbenthic(pathin)
}
\arguments{
\item{pathin}{A path to unzipped csv files with base tables used to calculate benthic index}
}
\value{
A nested \code{\link[tibble]{tibble}} of station, field sample, and taxa data
}
\description{
Format benthic data for the Tampa Bay Benthic Index
}
\details{
Function is used internally within \code{\link{read_importbenthic}}
}
\examples{
\dontrun{

# location to download data
path <- '~/Desktop/benthic.zip'

# download
urlin <- 'ftp://ftp.epchc.org/EPC_ERM_FTP/Benthic_Monitoring/DataDeliverables/BaseTables.zip'
read_dlcurrent(path, download_latest = TRUE, urlin = urlin)

# unzip
tmppth <- tempfile()
utils::unzip(path, exdir = tmppth, overwrite = TRUE)

# format benthic data
read_formbenthic(pathin = tmppth)

# remove temporary path
unlink(tmppth, recursive = TRUE)

}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_segplotly.R
\name{show_segplotly}
\alias{show_segplotly}
\title{Plot chlorophyll and secchi data together with matrix outcomes}
\usage{
show_segplotly(
  epcdata,
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  yrrng = c(1975, 2021),
  family = NULL,
  partialyr = FALSE,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{bay_segment}{chr string for the bay segment, one of "OTB", "HB", "MTB", "LTB"}

\item{yrrng}{numeric for year range to plot}

\item{family}{optional chr string indicating font family for text labels}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}

\item{width}{numeric for width of the plot in pixels}

\item{height}{numeric for height of the plot in pixels}
}
\value{
An interactive plotly object
}
\description{
Plot chlorophyll and secchi data together with matrix outcomes
}
\details{
This function combines outputs from \code{\link{show_thrplot}} and \code{\link{show_segmatrix}} for a selected bay segment. The plot is interactive and can be zoomed by dragging the mouse pointer over a section of the plot. Information about each cell or value can be seen by hovering over a location in the plot.
}
\examples{
show_segplotly(epcdata)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_tbbiscr.R
\name{anlz_tbbiscr}
\alias{anlz_tbbiscr}
\title{Get Tampa Bay Benthic Index scores}
\usage{
anlz_tbbiscr(benthicdata)
}
\arguments{
\item{benthicdata}{nested \code{\link[tibble]{tibble}} formatted from \code{\link{read_importbenthic}}}
}
\value{
A single data frame of TBBI scores for each site.
}
\description{
Get Tampa Bay Benthic Index scores
}
\details{
This function calculates scores for the TBBI based on station, taxa, and field sample data.  The total TBBI scores are returned as \code{TBBI} and \code{TBBICat}, where the latter is a categorical description of the scores.
}
\examples{
anlz_tbbiscr(benthicdata)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_dlcurrent.R
\name{read_dlcurrent}
\alias{read_dlcurrent}
\title{Download latest file from epchc.org}
\usage{
read_dlcurrent(xlsx, download_latest = TRUE, urlin)
}
\arguments{
\item{xlsx}{chr string path for local excel file, to overwrite it not current}

\item{download_latest}{logical to download latest file regardless of local copy}

\item{urlin}{url for file location}
}
\value{
The local copy specified in the path by \code{xlsx} is overwritten by the new file is not current or \code{download_latest = TRUE}.  The function does nothing if \code{download_latest = FALSE}.
}
\description{
Download latest file from epchc.org
}
\details{
The local copy is checked against a temporary file downloaded from the location specified by \code{urlin}.  The local file is replaced with the downloaded file if the MD5 hashes are different.
}
\examples{
\dontrun{
xlsx <- '~/Desktop/2018_Results_Updated.xls'
urlin <- 'ftp://ftp.epchc.org/EPC_ERM_FTP/WQM_Reports/'
urlin <- paste0(urlin, 'RWMDataSpreadsheet_ThroughCurrentReportMonth.xlsx')
read_dlcurrent(xlsx, urlin = urlin)
}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbseglines.R
\docType{data}
\name{tbseglines}
\alias{tbseglines}
\title{Spatial data object of lines defining major Tampa Bay segments}
\format{
A simple features \code{\link[sf]{sf}} object (LINESTRING) with 3 features and 1 field, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
}
\usage{
tbseglines
}
\description{
Spatial data object of lines defining major Tampa Bay segments
}
\examples{
library(sf)
plot(st_geometry(tbseglines))
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_transectavespp.R
\name{anlz_transectavespp}
\alias{anlz_transectavespp}
\title{Get annual averages of seagrass frequency occurrence by bay segments, year, and species}
\usage{
anlz_transectavespp(
  transectocc,
  bay_segment = c("OTB", "HB", "MTB", "LTB", "BCB"),
  yrrng = c(1998, 2021),
  species = c("Halodule", "Syringodium", "Thalassia", "Ruppia", "Halophila",
    "Caulerpa"),
  total = TRUE,
  by_seg = FALSE
)
}
\arguments{
\item{transectocc}{data frame returned by \code{\link{anlz_transectocc}}}

\item{bay_segment}{chr string for the bay segment, one to many of "OTB", "HB", "MTB", "LTB", "BCB"}

\item{yrrng}{numeric indicating year ranges to evaluate}

\item{species}{chr string of species to summarize, one to many of "Halodule", "Syringodium", "Thalassia", "Ruppia", "Halophila", "Caulerpa"}

\item{total}{logical indicating if total frequency occurrence for all species is also returned}

\item{by_seg}{logical indicating if separate results by bay segments are retained}
}
\value{
A data frame of annual averages by bay segment
}
\description{
Get annual averages of seagrass frequency occurrence by bay segments, year, and species
}
\details{
Frequency occurrence estimates are averaged across segments in \code{bay_segment} if \code{by_seg = F}, i.e., separate results by location are not returned.  Results are retained by bay segment if \code{by_seg = T}.  Also note that totals across species (\code{total = T}) are not returned if \code{by_seg = T}.
}
\examples{
\dontrun{
transect <- read_transect()
}
transectocc <- anlz_transectocc(transect)
anlz_transectavespp(transectocc)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tdlcrk.R
\name{show_tdlcrk}
\alias{show_tdlcrk}
\title{Make a map for tidal creek report card}
\usage{
show_tdlcrk(dat, weight = 1.5)
}
\arguments{
\item{dat}{input creek score data returned from \code{\link{anlz_tdlcrk}}}

\item{weight}{numeric for weight of polylines, passed to \code{\link[leaflet]{addPolylines}}}
}
\value{
A \code{\link[leaflet]{leaflet}} object
}
\description{
Make a map for tidal creek report card
}
\examples{
dat <- anlz_tdlcrk(tidalcreeks, iwrraw, yr = 2021)
show_tdlcrk(dat)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tdlcrkline.R
\name{show_tdlcrkline}
\alias{show_tdlcrkline}
\title{Add a line or annotation to a plotly graph}
\usage{
show_tdlcrkline(
  varin = c("CHLAC", "TN", "chla_tn_ratio", "DO", "tsi", "no23_ratio"),
  thrsel = FALSE,
  horiz = TRUE,
  annotate = FALSE
)
}
\arguments{
\item{varin}{chr string for the indicator}

\item{thrsel}{logical if something is returned, otherwise NULL, this is a hack for working with the plotly output}

\item{horiz}{logical indicating if output is horizontal or vertical}

\item{annotate}{logical indicating if output is line or annotation text}
}
\value{
A list object passed to the layout argument of plotly, either shapes or annotate depending on user input
}
\description{
Add a line or annotation to a plotly graph for the tidal creek indicators
}
\details{
This function is used internally within \code{\link{show_tdlcrkindic}} and \code{\link{show_tdlcrkindiccdf}}
}
\examples{
# code for vertical line output, chloropyll
show_tdlcrkline('CHLAC', thrsel = TRUE)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_compplot.R
\name{show_compplot}
\alias{show_compplot}
\title{Make a bar plot for transect training group comparisons}
\usage{
show_compplot(
  transect,
  yr,
  site,
  species = c("Halodule", "Syringodium", "Thalassia", "Halophila", "Ruppia"),
  varplo = c("Abundance", "Blade Length", "Short Shoot Density"),
  base_size = 18,
  xtxt = 10,
  size = 1
)
}
\arguments{
\item{transect}{data frame returned by \code{\link{read_transect}}}

\item{yr}{numeric for year of training data to plot}

\item{site}{chr string indicating site results to plot}

\item{species}{chr string indicating which species to plot}

\item{varplo}{chr string indicating which variable to plot}

\item{base_size}{numeric indicating text scaling size for plot}

\item{xtxt}{numeric indicating text size for x-axis labels}

\item{size}{numeric indicating line size}
}
\value{
A \code{\link[ggplot2]{ggplot}} object
}
\description{
Make a bar plot for transect training group comparisons
}
\examples{
transect <- read_transect(training = TRUE)
show_compplot(transect, yr = 2021, site = '1', species = 'Halodule', varplo = 'Abundance')
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidalcreeks.R
\docType{data}
\name{tidalcreeks}
\alias{tidalcreeks}
\title{Spatial data object of tidal creeks}
\format{
A simple features \code{\link[sf]{sf}} object (MULTILINESTRING) with 609 features and 6 fields, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
\describe{
  \item{id}{num}
  \item{wbid}{chr}
  \item{JEI}{chr}
  \item{class}{chr}
  \item{name}{chr}
  \item{Creek_Length_m}{num}
}
}
\usage{
tidalcreeks
}
\description{
Spatial data object of tidal creeks
}
\examples{
\dontrun{
library(sf)
library(tidyverse)

prj <- 4326

# create sf object of creek population, join with creek length data
tidalcreeks <- st_read(
  dsn = '../../02_DOCUMENTS/tidal_creeks/TidalCreek_ALL_Line_WBID61.shp',
  drivers = 'ESRI Shapefile'
  ) \%>\%
  st_transform(crs = prj) \%>\%
  mutate(
    id = 1:nrow(.)
  ) \%>\%
  select(id, name = Name, JEI = CreekID, wbid = WBID, class = CLASS, Creek_Length_m = Total_m)

# save
save(tidalcreeks, file = 'data/tidalcreeks.RData', compress = 'xz')

}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fimstations.R
\docType{data}
\name{fimstations}
\alias{fimstations}
\title{Spatial data object of FIM stations including Tampa Bay segments}
\format{
A simple features \code{\link[sf]{sf}} object (POINT) with 6796 features and 2 fields, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
\describe{
  \item{Reference}{num}
  \item{bay_segment}{chr}
}
}
\usage{
fimstations
}
\description{
Spatial data object of FIM stations including Tampa Bay segments
}
\examples{
\dontrun{
# file path
csv <- '~/Desktop/fimraw.csv'

# load and assign to object
fimstations <- read_importfim(csv, download_latest = TRUE, locs = TRUE)
save(fimstations, file = 'data/fimstations.RData', compress = 'xz')
}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tbniscrplotly.R
\name{show_tbniscrplotly}
\alias{show_tbniscrplotly}
\title{Creates a plotly object for TBNI score plots}
\usage{
show_tbniscrplotly(p, width = NULL, height = NULL)
}
\arguments{
\item{p}{\code{\link[ggplot2]{ggplot}} object as output from \code{\link{show_tbniscr}} or \code{\link{show_tbniscrall}}}

\item{width}{numeric for width of the plot in pixels}

\item{height}{numeric for height of the plot in pixels}
}
\value{
A \code{\link[plotly]{plotly}} data object
}
\description{
Creates a plotly object for TBNI score plots
}
\examples{
tbniscr <- anlz_tbniscr(fimdata)
p <- show_tbniscrall(tbniscr)
show_tbniscrplotly(p)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsmap.R
\docType{data}
\name{bsmap}
\alias{bsmap}
\title{Terrain basemap}
\format{
A \code{\link[ggmap]{ggmap}} object
}
\usage{
bsmap
}
\description{
Terrain basemap
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tdlcrkindiccdf.R
\name{show_tdlcrkindiccdf}
\alias{show_tdlcrkindiccdf}
\title{Plotly empirical CDF plots of tidal creek context indicators}
\usage{
show_tdlcrkindiccdf(
  id,
  cntdat,
  yr = 2021,
  thrsel = FALSE,
  pal = c("#5C4A42", "#427355", "#004F7E")
)
}
\arguments{
\item{id}{numeric indicating the \code{id} number of the tidal creek to plot}

\item{cntdat}{output from \code{\link{anlz_tdlcrkindic}}}

\item{yr}{numeric indicating reference year}

\item{thrsel}{logical if threshold lines and annotations are shown on the plots}

\item{pal}{vector of colors for the palette}
}
\value{
A plotly object
}
\description{
Plotly empirical CDF plots of tidal creek context indicators
}
\details{
This function returns several empirical cumulative distribution plots for the tidal creek context indicators.  Points on the plot indicate the observed values and percentiles for the creek specified by \code{id}. The percentiles and CDF values are defined by the "population" of creeks in \code{cntdat}.  Points in the plots are color-coded by sample year to evaluate temporal trends, if any.
}
\examples{
cntdat <- anlz_tdlcrkindic(tidalcreeks, iwrraw, yr = 2021)
show_tdlcrkindiccdf(495, cntdat, thrsel = TRUE)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tdlcrkindic.R
\name{show_tdlcrkindic}
\alias{show_tdlcrkindic}
\title{Plotly barplots of tidal creek context indicators}
\usage{
show_tdlcrkindic(
  id,
  cntdat,
  yr = 2021,
  thrsel = FALSE,
  pal = c("#5C4A42", "#427355", "#004F7E")
)
}
\arguments{
\item{id}{numeric indicating the \code{id} number of the tidal creek to plot}

\item{cntdat}{output from \code{\link{anlz_tdlcrkindic}}}

\item{yr}{numeric indicating reference year}

\item{thrsel}{logical if threshold lines and annotations are shown on the plots}

\item{pal}{vector of colors for the palette}
}
\value{
A plotly object
}
\description{
Plotly barplots of tidal creek context indicators
}
\examples{
cntdat <- anlz_tdlcrkindic(tidalcreeks, iwrraw, yr = 2021)
show_tdlcrkindic(495, cntdat, thrsel = TRUE)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seagrass.R
\docType{data}
\name{seagrass}
\alias{seagrass}
\title{Seagrass coverage by year}
\format{
A data frame used to create the flagship seagrass coverage graphic:
\describe{
  \item{Year}{int}
  \item{Acres}{num}
  \item{Hectares}{num}
}
}
\usage{
seagrass
}
\description{
Seagrass coverage by year
}
\examples{
\dontrun{

seagrass <- structure(list(
   Year = c(1950L, 1982L, 1988L, 1990L, 1992L, 1994L, 1996L,
        1999L, 2001L, 2004L, 2006L, 2008L, 2010L, 2012L, 2014L,
        2016L, 2018L, 2020L),
   Acres = c(40420, 21650, 23285, 25226, 25753, 26518, 26916,
        24841, 26078, 27021, 28299, 29647, 32897, 34642, 40294.71,
        41655.16, 40651.55, 34298),
   Hectares = c(16357.39, 8761.44, 9423.11, 10208.6, 10421.87,
        10731.45, 10892.52, 10052.8, 10553.39, 10935.01, 11452.2,
        11997.72, 13312.94, 14019.27, 16306.69, 16857.25, 16451.1,
        13880)
 ), class = "data.frame", row.names = c(NA, -18L))

save(seagrass, file = 'data/seagrass.RData', compress = 'xz')

}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_importbenthic.R
\name{read_importbenthic}
\alias{read_importbenthic}
\title{Download and import benthic data for Tampa Bay}
\usage{
read_importbenthic(path, download_latest = FALSE, remove = FALSE)
}
\arguments{
\item{path}{chr string for local path where the zipped folder will be downloaded, must include .zip extension}

\item{download_latest}{logical to download latest if a more recent dataset is available}

\item{remove}{logical if the downloaded folder is removed after unzipping}
}
\value{
A nested \code{tibble} of station, taxa, and field sample data
}
\description{
Download and import benthic data for Tampa Bay
}
\details{
This function downloads and unzips a folder of base tables used to calculate the benthic index from \url{ftp://ftp.epchc.org/EPC_ERM_FTP/Benthic_Monitoring/DataDeliverables/BaseTables.zip}.
}
\examples{
\dontrun{
# location to download data
path <- '~/Desktop/benthic.zip'

# load and assign to object
benthicdata <- read_importbenthic(path, download_latest = TRUE)

}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stations.R
\docType{data}
\name{stations}
\alias{stations}
\title{Bay stations by segment}
\format{
A data frame with 45 rows and 4 variables:
\describe{
  \item{bay_segment}{chr}
  \item{epchc_station}{num}
  \item{Latitude}{num}
  \item{Longitude}{num}
}
}
\usage{
stations
}
\description{
Bay stations by segment
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tbbimatrix.R
\name{show_tbbimatrix}
\alias{show_tbbimatrix}
\title{Plot a matrix of Tampa Bay Benthic Index scores over time by bay segment}
\usage{
show_tbbimatrix(
  tbbiscr,
  bay_segment = c("HB", "OTB", "MTB", "LTB", "TCB", "MR", "BCB", "All", "All (wt)"),
  yrrng = c(1993, 2019),
  alph = 1,
  txtsz = 3,
  family = NA,
  rev = FALSE,
  position = "top",
  plotly = FALSE,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{tbbiscr}{input data frame as returned by \code{\link{anlz_tbbiscr}}}

\item{bay_segment}{chr string for the bay segment, one to many of "HB", "OTB", "MTB", "LTB", "TCB", "MR", "BCB", "All", "All (wt)"}

\item{yrrng}{numeric indicating year ranges to evaluate}

\item{alph}{numeric indicating alpha value for score category colors}

\item{txtsz}{numeric for size of text in the plot}

\item{family}{optional chr string indicating font family for text labels}

\item{rev}{logical if factor levels for bay segments are reversed}

\item{position}{chr string of location for bay segment labels, default on top, passed to \code{\link[ggplot2]{scale_x_discrete}}}

\item{plotly}{logical if matrix is created using plotly}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}
}
\value{
A \code{\link[ggplot2]{ggplot}} object showing trends over time in TBBI scores for each bay segment if \code{plotly = FALSE}, otherwise a \code{\link[plotly]{plotly}} object
}
\description{
Plot a matrix of Tampa Bay Benthic Index scores over time by bay segment
}
\details{
Additional summaries are provided for the entire bay, as a summary across categories ("All") and a summary weighted across the relative sizes of each bay segment ("All (wt)").
}
\examples{
tbbiscr <- anlz_tbbiscr(benthicdata)
show_tbbimatrix(tbbiscr)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trnpts.R
\docType{data}
\name{trnpts}
\alias{trnpts}
\title{Seagrass transect starting locations}
\format{
A \code{sf} POINT object
}
\usage{
trnpts
}
\description{
Seagrass transect starting locations
}
\examples{
\dontrun{
library(sf)
library(dplyr)
library(tbeptools)

trnpts <- st_read('T:/05_GIS/SEAGRASS_TRANSECTS/TransectBasics2019.shp') \%>\%
   st_transform(crs = 4326) \%>\%
   dplyr::rename(MonAgency = 'MON_AGENCY') \%>\%
   dplyr::filter(!as.character(TRAN_ID) \%in\% c('S8T1', 'S8T2', 'S8T3', 'S3T2')) \%>\%
   sf::st_intersection(sf::st_make_valid(tbsegshed)) \%>\%
   dplyr::select(-long_name) \%>\%
   dplyr::mutate_if(is.factor, as.character)

save(trnpts, file = 'data/trnpts.RData', compress = 'xz')
}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tbniscr.R
\name{show_tbniscr}
\alias{show_tbniscr}
\title{Plot Tampa Bay Nekton Index scores over time by bay segment}
\usage{
show_tbniscr(
  tbniscr,
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  perc = c(32, 46),
  alph = 0.3,
  ylim = c(22, 58),
  rev = FALSE,
  plotly = FALSE,
  family = NA,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{tbniscr}{input dat frame as returned by \code{\link{anlz_tbniscr}}}

\item{bay_segment}{chr string for the bay segment, one to many of "OTB", "HB", "MTB", "LTB"}

\item{perc}{numeric values indicating break points for score categories}

\item{alph}{numeric indicating alpha value for score category colors}

\item{ylim}{numeric for y axis limits}

\item{rev}{logical if factor levels for bay segments are reversed}

\item{plotly}{logical if matrix is created using plotly}

\item{family}{optional chr string indicating font family for text labels}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}
}
\value{
A \code{\link[ggplot2]{ggplot}} object showing trends over time in TBNI scores for each bay segment or a \code{\link[plotly]{plotly}} object if \code{plotly = TRUE}
}
\description{
Plot Tampa Bay Nekton Index scores over time by bay segment
}
\examples{
tbniscr <- anlz_tbniscr(fimdata)
show_tbniscr(tbniscr)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_tbnimatrix.R
\name{show_tbnimatrix}
\alias{show_tbnimatrix}
\title{Plot a matrix of Tampa Bay Nekton Index scores over time by bay segment}
\usage{
show_tbnimatrix(
  tbniscr,
  bay_segment = c("OTB", "HB", "MTB", "LTB"),
  perc = c(32, 46),
  alph = 0.3,
  txtsz = 3,
  family = NA,
  rev = FALSE,
  position = "top",
  plotly = FALSE,
  width = NULL,
  height = NULL
)
}
\arguments{
\item{tbniscr}{input data frame as returned by \code{\link{anlz_tbniscr}}}

\item{bay_segment}{chr string for the bay segment, one to many of "OTB", "HB", "MTB", "LTB"}

\item{perc}{numeric values indicating break points for score categories}

\item{alph}{numeric indicating alpha value for score category colors}

\item{txtsz}{numeric for size of text in the plot}

\item{family}{optional chr string indicating font family for text labels}

\item{rev}{logical if factor levels for bay segments are reversed}

\item{position}{chr string of location for bay segment labels, default on top, passed to \code{\link[ggplot2]{scale_x_discrete}}}

\item{plotly}{logical if matrix is created using plotly}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}
}
\value{
A \code{\link[ggplot2]{ggplot}} object showing trends over time in TBNI scores for each bay segment
}
\description{
Plot a matrix of Tampa Bay Nekton Index scores over time by bay segment
}
\examples{
tbniscr <- anlz_tbniscr(fimdata)
show_tbnimatrix(tbniscr)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epcdata.R
\docType{data}
\name{epcdata}
\alias{epcdata}
\title{All bay data as of 01242022}
\format{
A data frame with 26791 rows and 26 variables:
\describe{
  \item{bay_segment}{chr}
  \item{epchc_station}{num}
  \item{SampleTime}{POSIXct}
  \item{yr}{num}
  \item{mo}{num}
  \item{Latitude}{num}
  \item{Longitude}{num}
  \item{Total_Depth_m}{num}
  \item{Sample_Depth_m}{num}
  \item{tn}{num}
  \item{tn_q}{chr}
  \item{sd_m}{num}
  \item{sd_raw_m}{num}
  \item{sd_q}{chr}
  \item{chla}{num}
  \item{chla_q}{chr}
  \item{Sal_Top_ppth}{num}
  \item{Sal_Mid_ppth}{num}
  \item{Sal_Bottom_ppth}{num}
  \item{Temp_Water_Top_degC}{num}
  \item{Temp_Water_Mid_degC}{num}
  \item{Temp_Water_Bottom_degC}{num}
  \item{Turbidity_JTU-NTU}{num}
  \item{Turbidity_Q}{num}
  \item{Color_345_F45_PCU}{num}
  \item{Color_345_F45_Q}{num}
  }
}
\usage{
epcdata
}
\description{
All bay data as of 01242021
}
\examples{
\dontrun{
xlsx <- '~/Desktop/epcdata.xls'
epcdata <- read_importwq(xlsx, download_latest = TRUE)

nrow(epcdata)
ncol(epcdata)

save(epcdata, file = 'data/epcdata.RData', compress = 'xz')
}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbsegshed.R
\docType{data}
\name{tbsegshed}
\alias{tbsegshed}
\title{Spatial data object of Tampa Bay segments plus watersheds}
\format{
A simple features \code{\link[sf]{sf}} object (POLYGON) with 7 features and 2 fields, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
\describe{
  \item{long_name}{chr}
  \item{bay_segment}{chr}
}
}
\usage{
tbsegshed
}
\description{
Spatial data object of Tampa Bay segments plus waterhseds
}
\examples{
library(sf)
plot(st_geometry(tbsegshed))
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_transectocc.R
\name{anlz_transectocc}
\alias{anlz_transectocc}
\title{Get seagrass average abundance and occurrence across transects}
\usage{
anlz_transectocc(transect)
}
\arguments{
\item{transect}{data frame returned by \code{\link{read_transect}}}
}
\value{
A data frame with abundance and frequency occurrence estimates aggregated by species, transect, and date.  The nsites column is the total number of placements that were sampled along a transect for a particular date.
}
\description{
Get seagrass average abundance and occurrence across transects
}
\details{
Abundance and frequency occurrence are estimated as in Sherwood et al. 2017, equations 1 and 2.  In short, frequency occurrence is estimated as the number of instances a species was observed along a transect divided by the number of placements along a transect and average abundance was estimated as the sum of species-specific Braun-Blanquet scores divided by the number of placements along a transect.  The estimates are obtained for all seagrass species including Caulerpa, whereas all attached and drift algae species are aggregated.
}
\examples{
\dontrun{
transect <- read_transect()
}
anlz_transectocc(transect)
}
\references{
Sherwood, E.T., Greening, H.S., Johansson, J.O.R., Kaufman, K., Raulerson, G.E. 2017. Tampa Bay (Florida, USA): Documenting seagrass recovery since the 1980's and reviewing the benefits. Southeastern Geographer. 57(3):294-319.
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_chkdate.R
\name{read_chkdate}
\alias{read_chkdate}
\title{Compare date of local xlsx file with the same file on the server}
\usage{
read_chkdate(urlin, xlsx)
}
\arguments{
\item{urlin}{chr string of full path to file on the server}

\item{xlsx}{chr string of full path to local file}
}
\value{
A logical vector indicating if the local file is current
}
\description{
Compare date of local xlsx file with the same file on the server
}
\examples{
\dontrun{
urlin <- 'ftp://ftp.epchc.org/EPC_ERM_FTP/WQM_Reports/'
urlin <- paste0(urlin, 'RWMDataSpreadsheet_ThroughCurrentReportMonth.xlsx')
xlsx <- '~/Desktop/2018_Results_Updated.xls'
read_chkdate(urlin, xlsx)
}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_avedat.R
\name{anlz_avedat}
\alias{anlz_avedat}
\title{Estimate annual means}
\usage{
anlz_avedat(epcdata, partialyr = FALSE)
}
\arguments{
\item{epcdata}{\code{data.frame} formatted from \code{read_importwq}}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter}
}
\value{
Mean estimates for chlorophyll and secchi
}
\description{
Estimate annual means for chlorophyll and secchi data
}
\examples{
# view average estimates
anlz_avedat(epcdata)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_hydroload.R
\name{anlz_hydroload}
\alias{anlz_hydroload}
\title{Estimate hydrological estimates and adjustment factors for bay segments}
\usage{
anlz_hydroload(yrs, noaa_key = NULL, trace = FALSE)
}
\arguments{
\item{yrs}{numeric vector indicating years to return}

\item{noaa_key}{user-supplied NOAA key, see details}

\item{trace}{logical indicating if function progress is printed in the consol}
}
\value{
A data frame with hydrological load estimates by bay segments for the requested years
}
\description{
Estimate hydrological estimates and adjustment factors for bay segments
}
\details{
This function uses rainfall and streamflow data from NOAA and USGS and requires an API key.  See the "Authentication" section under the help file for \code{\link[rnoaa]{ncdc}}.  This key can be added to the R environment file and called for later use, see the examples.

These estimates are used in annual compliance assessment reports produced by the Tampa Bay Nitrogen Management Consortium. Load estimates and adjustment factors are based on regression models in https://drive.google.com/file/d/11NT0NQ2WbPO6pVZaD7P7Z6qjcwO1jxHw/view?usp=drivesdk
}
\examples{
\dontrun{
# this function requires an API key
# save it to the R environment file (only once)
# save the key, do only once
cat("NOAA_KEY=XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n",
   file=file.path(normalizePath("~/"), ".Renviron"),
   append=TRUE)

# retrieve the key after saving, may need to restart R
noaa_key <- Sys.getenv('NOAA_key')

# get estimates for 2021
anlz_hydroload(2021, noaa_key)

}
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_importwq.R
\name{read_importwq}
\alias{read_importwq}
\title{Load local water quality file}
\usage{
read_importwq(xlsx, download_latest = FALSE, na = "")
}
\arguments{
\item{xlsx}{chr string path for local excel file, to overwrite if not current}

\item{download_latest}{logical passed to \code{\link{read_dlcurrent}} to download raw data and compare with existing in \code{xlsx} if available}

\item{na}{chr vector of strings to interpret as \code{NA}, passed to \code{\link[readxl]{read_xlsx}}}
}
\value{
A \code{data.frame} of formatted water quality data.
}
\description{
Load local water quality file
}
\details{
Loads the "RWMDataSpreadsheet" worksheet from the file located at \code{xlsx}.  The file is downloaded from \url{ftp://ftp.epchc.org/EPC_ERM_FTP/WQM_Reports/RWMDataSpreadsheet_ThroughCurrentReportMonth.xlsx}.
}
\examples{
\dontrun{
# file path
xlsx <- '~/Desktop/2018_Results_Updated.xls'

# load and assign to object
epcdata <- read_importwq(xlsx, download_latest = T)
}
}
\seealso{
\code{\link{read_formwq}}, \code{\link{read_importphyto}}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sitemap.R
\name{show_sitemap}
\alias{show_sitemap}
\title{Map site attainment categories for a selected year}
\usage{
show_sitemap(
  epcdata,
  yrsel,
  mosel = c(1, 12),
  param = c("chla", "la"),
  trgs = NULL,
  thrs = FALSE,
  partialyr = FALSE
)
}
\arguments{
\item{epcdata}{data frame of epc data returned by \code{\link{read_importwq}}}

\item{yrsel}{numeric for year to plot}

\item{mosel}{optional numeric of length one or two for mapping results for a specific month or month range in a given year, default full year}

\item{param}{chr string for which parameter to plot, one of \code{"chla"} for chlorophyll or \code{"la"} for light attenuation}

\item{trgs}{optional \code{data.frame} for annual bay segment water quality targets, defaults to \code{\link{targets}}, only applies if \code{mosel = c(1, 12)}}

\item{thrs}{logical indicating if attainment category is relative to targets (default) or thresholds, passed to \code{\link{anlz_attainsite}}, only applies if \code{mosel = c(1, 12)}}

\item{partialyr}{logical indicating if incomplete annual data for the most recent year are approximated by five year monthly averages for each parameter, only applies if \code{mosel = c(1, 12)}}
}
\value{
A static \code{ggplot} object is returned
}
\description{
Map site attainment categories for a selected year
}
\examples{
show_sitemap(epcdata, yrsel = 2021)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_importfim.R
\name{read_importfim}
\alias{read_importfim}
\title{Load local FIM data for the Tampa Bay Nekton Index}
\usage{
read_importfim(csv, download_latest = FALSE, locs = FALSE)
}
\arguments{
\item{csv}{chr string path for local csv file, to overwrite if not current}

\item{download_latest}{logical passed to \code{\link{read_dlcurrent}} to download raw data and compare with existing in \code{csv} if available}

\item{locs}{logical indicating if a spatial features object is returned with locations of each FIM sampling station}
}
\value{
A formatted \code{data.frame} with FIM data if \code{locs = FALSE}, otherwise a simple features object if \code{locs = TRUE}
}
\description{
Load local FIM data for the Tampa Bay Nekton Index
}
\details{
Data downloaded from \url{'ftp://ftp.floridamarine.org/users/fim/tmac/NektonIndex/TampaBay_NektonIndexData.csv'}.
}
\examples{
\dontrun{
# file path
csv <- '~/Desktop/fimraw.csv'

# load and assign to object
fimdata <- read_importfim(csv, download_latest = TRUE)
}
}
\seealso{
\code{\link{read_formwq}}, \code{\link{read_importphyto}}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_transect.R
\name{show_transect}
\alias{show_transect}
\title{Plot results for a seagrass transect by time and location}
\usage{
show_transect(
  transect,
  site,
  species = c("Halodule", "Syringodium", "Thalassia", "Halophila", "Ruppia",
    "Caulerpa"),
  yrrng = c(1998, 2021),
  varplo = c("Abundance", "Blade Length", "Short Shoot Density"),
  base_size = 12,
  facet = FALSE,
  ncol = NULL,
  plotly = FALSE,
  width = NULL,
  height = NULL,
  sppcol = NULL
)
}
\arguments{
\item{transect}{data frame returned by \code{\link{read_transect}}}

\item{site}{chr string indicating site results to plot}

\item{species}{chr string indicating one to many of which species to plot}

\item{yrrng}{numeric indicating year ranges to evaluate}

\item{varplo}{chr string indicating which variable to plot}

\item{base_size}{numeric indicating text scaling size for plot}

\item{facet}{logical indicating if plots are separated into facets by species}

\item{ncol}{numeric indicating number of columns if \code{facet = TRUE}}

\item{plotly}{logical if plot is created using plotly}

\item{width}{numeric for width of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{height}{numeric for height of the plot in pixels, only applies of \code{plotly = TRUE}}

\item{sppcol}{character vector of alternative colors to use for each species, must have length of six}
}
\value{
A \code{\link[ggplot2]{ggplot}} object
}
\description{
Plot results for a seagrass transect by time and location
}
\details{
All sites along a transect that were surveyed are shown in the plot, including those where the selected species were not found.  The latter is colored in grey hollow points.  Species options include Halodule, Syringodium, Thalassia, Halophila, Ruppia, and/or Caulerpa.

Note that if \code{plotly = TRUE}, the size legend is not shown.
}
\examples{
\dontrun{
transect <- read_transect()
}

# one species
show_transect(transect, site = 'S3T10', species = 'Halodule', varplo = 'Abundance')

# multiple species, one plot
show_transect(transect, site = 'S3T10',
  species = c('Halodule', 'Syringodium', 'Thalassia', 'Halophila', 'Ruppia', 'Caulerpa'),
  varplo = 'Abundance')

# multiple species, multiple plots
show_transect(transect, site = 'S3T10',
  species = c('Halodule', 'Syringodium', 'Thalassia', 'Halophila', 'Ruppia', 'Caulerpa'),
  varplo = 'Abundance', facet = TRUE)
}
\concept{show}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transect.R
\docType{data}
\name{transect}
\alias{transect}
\title{Seagrass transect data for Tampa Bay current as of 01212022}
\format{
A data frame with 140460 rows and 11 variables:
\describe{
  \item{Crew}{chr}
  \item{MonitoringAgency}{chr}
  \item{Date}{Date}
  \item{Transect}{chr}
  \item{Site}{chr}
  \item{Depth}{int}
  \item{Savspecies}{chr}
  \item{SeagrassEdge}{num}
  \item{var}{chr}
  \item{aveval}{num}
  \item{sdval}{num}
  }
}
\usage{
transect
}
\description{
Seagrass transect data for Tampa Bay current as of 01212022
}
\examples{
\dontrun{

transect <- read_transect()

save(transect, file = 'data/transect.RData', compress = 'xz')

}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_importphyto.R
\name{read_importphyto}
\alias{read_importphyto}
\title{Load local phytoplankton cell count file}
\usage{
read_importphyto(xlsx, download_latest = FALSE, na = "")
}
\arguments{
\item{xlsx}{chr string path for local excel file, to overwrite if not current}

\item{download_latest}{logical passed to \code{\link{read_dlcurrent}} to download raw data and compare with existing in \code{xlsx} if available}

\item{na}{chr vector of strings to interpret as \code{NA}, passed to \code{\link[readxl]{read_xlsx}}}
}
\value{
A \code{data.frame} of formatted water quality data.
}
\description{
Load local phytoplankton cell count file
}
\details{
Pytoplankton cell count data downloaded from ftp://ftp.epchc.org/EPC_ERM_FTP/WQM_Reports/, file PlanktonDataList_ThroughCurrentReportMonth.xlsx
}
\examples{
\dontrun{
# file path
xlsx <- '~/Desktop/phyto_data.xlsx'

# load and assign to object
phytodata <- read_importphyto(xlsx)
}
}
\seealso{
\code{\link{read_importwq}}
}
\concept{read}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iwrraw.R
\docType{data}
\name{iwrraw}
\alias{iwrraw}
\title{FDEP IWR run 61}
\format{
A data frame 405682 rows and 15 variables
}
\usage{
iwrraw
}
\description{
Florida Department of Environmental Protection, Impaired Waters Rule, Run 61
}
\examples{
\dontrun{
library(dplyr)

load(file = '../../02_DOCUMENTS/tidal_creeks/iwrraw_run61.RData')
iwrraw <- sf::st_set_geometry(iwrraw, NULL) \%>\%
  rename(JEI = jei)
save(iwrraw, file = 'data/iwrraw.RData', compress = 'xz')
}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benthicdata.R
\docType{data}
\name{benthicdata}
\alias{benthicdata}
\title{Benthic data for the Tampa Bay Benthic Index current as of 08012021}
\format{
A nested \code{\link[tibble]{tibble}} with 3 rows and 2 variables:
\describe{
  \item{name}{chr}
  \item{value}{list}
  }
}
\usage{
benthicdata
}
\description{
Benthic data for the Tampa Bay Benthic Index current as of 08012021
}
\examples{
\dontrun{
# location to download data
path <- '~/Desktop/benthic.zip'

# load and assign to object
benthicdata <- read_importbenthic(path, download_latest = TRUE, remove = TRUE)

save(benthicdata, file = 'data/benthicdata.RData', compress = 'xz')

}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgmanagement.R
\docType{data}
\name{sgmanagement}
\alias{sgmanagement}
\title{Seagrass management areas for Tampa Bay}
\format{
A simple features \code{\link[sf]{sf}} object (MULTIPOLYGON) with 30 features and 1 field, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
\describe{
  \item{areas}{int}
}
}
\usage{
sgmanagement
}
\description{
Seagrass management areas for Tampa Bay
}
\details{
These polygons are seagrass management areas for Tampa Bay that provide a finer division of areas within major segments (\code{\link{tbseg}}) having relevance for locations of seagrass beds.
}
\examples{
\dontrun{
library(sf)
library(tidyverse)
library(tools)

# NAD83(HARN) / Florida West (ftUS)
# same as sgseg
prj <- 2882

# create sf object of boundaries
sgmanagement <- st_read(
  dsn = '~/Desktop/TBEP/GISboundaries/Seagrass_Management_Areas/TBEP_SG_MA_FINAL_Projectfix.shp',
  drivers = 'ESRI Shapefile'
  ) \%>\%
  select(areas = TBEP_SG_MA) \%>\%
  st_zm() \%>\%
  st_transform(prj)

# save
save(sgmanagement, file = 'data/sgmanagement.RData', compress = 'xz')

}
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_refs.R
\name{anlz_refs}
\alias{anlz_refs}
\title{Convert references csv to bib}
\usage{
anlz_refs(path)
}
\arguments{
\item{path}{chr string of path to reference csv file or data frame object}
}
\value{
A data frame with references formatted as bib entries
}
\description{
Convert references csv to bib
}
\examples{

# input and format
path <- 'https://raw.githubusercontent.com/tbep-tech/tbep-refs/master/tbep-refs.csv'
bibs <- anlz_refs(path)

\dontrun{
# save output
 writeLines(bibs, 'formatted.bib')
}
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anlz_iwrraw.R
\name{anlz_iwrraw}
\alias{anlz_iwrraw}
\title{Format raw IWR data}
\usage{
anlz_iwrraw(iwrraw, tidalcreeks, yr = 2021)
}
\arguments{
\item{iwrraw}{FDEP impaired waters rule data base as \code{\link{data.frame}}}

\item{tidalcreeks}{\code{\link[sf]{sf}} object for population of tidal creeks}

\item{yr}{numeric for reference year to evaluate, scores are based on the planning period beginning ten years prior to this date}
}
\value{
A \code{\link{data.frame}} with the formatted data
}
\description{
Format raw IWR data
}
\details{
The function subsets the raw IWR data for the selected value in \code{yr} and the ten years prior to \code{yr} and subsets by the creek population in \code{\link{tidalcreeks}}. Select water quality parameters in \code{masterCode} are filtered and some of the names are combined for continuity.
}
\examples{
anlz_iwrraw(iwrraw, tidalcreeks, yr = 2021)
}
\concept{analyze}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tbseg.R
\docType{data}
\name{tbseg}
\alias{tbseg}
\title{Spatial data object of Tampa Bay segments}
\format{
A simple features \code{\link[sf]{sf}} object (POLYGON) with 4 features and 2 fields, +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
\describe{
  \item{long_name}{chr}
  \item{bay_segment}{chr}
}
}
\usage{
tbseg
}
\description{
Note that these boundaries are not used for formal analysis and are only used as visual aids in mapping.
}
\details{
Spatial data object of Tampa Bay segments
}
\examples{
library(sf)
plot(st_geometry(tbseg))
}
\concept{data}
\keyword{datasets}
