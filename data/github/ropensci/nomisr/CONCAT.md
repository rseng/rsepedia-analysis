
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nomisr

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/nomisr)](https://cran.r-project.org/package=nomisr)
[![GitHub
tag](https://img.shields.io/github/tag/ropensci/nomisr.svg)](https://github.com/ropensci/nomisr)
[![](https://cranlogs.r-pkg.org/badges/grand-total/nomisr)](https://cran.r-project.org/package=nomisr)
[![R build
status](https://github.com/ropensci/nomisr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/nomisr/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/nomisr/master.svg)](https://codecov.io/github/ropensci/nomisr?branch=master)
[![ropensci](https://badges.ropensci.org/190_status.svg)](https://github.com/ropensci/software-review/issues/190)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1157908.svg)](https://doi.org/10.5281/zenodo.1157908)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.00859/status.svg)](https://doi.org/10.21105/joss.00859)
<!-- badges: end -->

`nomisr` is for accessing UK official statistics from the
[Nomis](https://www.nomisweb.co.uk/) database through R. Nomis contains
data from the Census, the Labour Force Survey, DWP benefit statistics
and other economic and demographic data, and is maintained on behalf of
the Office for National Statistics by the University of Durham.

The `nomisr` package provides functions to find what data is available,
the variables and query options for different datasets and a function
for downloading data. `nomisr` returns data in
[`tibble`](https://cran.r-project.org/package=tibble) format. Most of
the data available through `nomisr` is based around statistical
geographies, with a handful of exceptions.

The package is for demographers, economists, geographers, public health
researchers and any other researchers who are interested in geographic
factors. The package aims to aid reproducibility, reduce the need to
manually download area profiles, and allow easy linking of different
datasets covering the same geographic area.

## Installation

`nomisr` is available on CRAN:

``` r
install.packages("nomisr")
```

You can install the development version `nomisr` from github with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/nomisr")
```

## Using `nomisr`

`nomisr` contains functions to search for datasets, identify the query
options for different datasets and retrieve data from queries, all done
with [`tibbles`](https://tibble.tidyverse.org/), to take advantage of
how `tibble` manages list-columns. The use of metadata queries, rather
than simply downloading all available data, is useful to avoid
overwhelming the rate limits of the API.

There are `nomisr` vignette
[introduction](https://docs.evanodell.com/nomisr/articles/introduction.html)
has details on all available functions and basic demonstrations of their
use.

The example below demonstrates a workflow to retrieve the latest data on
Jobseeker’s Allowance with rates and proportions, on a national level,
with all male claimants and workforce.

``` r
 library(nomisr)
 jobseekers_search <- nomis_search(name = "*Jobseeker*")
 
 tibble::glimpse(jobseekers_search)
#> Rows: 17
#> Columns: 14
#> $ agencyid                             <chr> "NOMIS", "NOMIS", "NOMIS", "NOMIS…
#> $ id                                   <chr> "NM_1_1", "NM_4_1", "NM_8_1", "NM…
#> $ uri                                  <chr> "Nm-1d1", "Nm-4d1", "Nm-8d1", "Nm…
#> $ version                              <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
#> $ annotations.annotation               <list> [<data.frame[10 x 2]>], [<data.f…
#> $ components.attribute                 <list> [<data.frame[7 x 4]>], [<data.fr…
#> $ components.dimension                 <list> [<data.frame[5 x 3]>], [<data.fr…
#> $ components.primarymeasure.conceptref <chr> "OBS_VALUE", "OBS_VALUE", "OBS_VA…
#> $ components.timedimension.codelist    <chr> "CL_1_1_TIME", "CL_4_1_TIME", "CL…
#> $ components.timedimension.conceptref  <chr> "TIME", "TIME", "TIME", "TIME", "…
#> $ description.value                    <chr> "Records the number of people cla…
#> $ description.lang                     <chr> "en", "en", NA, "en", "en", "en",…
#> $ name.value                           <chr> "Jobseeker's Allowance with rates…
#> $ name.lang                            <chr> "en", "en", "en", "en", "en", "en…

 jobseekers_measures <- nomis_get_metadata("NM_1_1", "measures")
 
 tibble::glimpse(jobseekers_measures)
#> Rows: 4
#> Columns: 3
#> $ id             <chr> "20100", "20201", "20202", "20203"
#> $ label.en       <chr> "claimants", "workforce", "active", "residence"
#> $ description.en <chr> "claimants", "workforce", "active", "residence"
 
 jobseekers_geography <- nomis_get_metadata("NM_1_1", "geography", "TYPE")
 
 tail(jobseekers_geography)
#> # A tibble: 6 × 3
#>   id      label.en                                           description.en     
#>   <chr>   <chr>                                              <chr>              
#> 1 TYPE490 government office regions tec / lec based          government office …
#> 2 TYPE491 government office regions (former inc. Merseyside) government office …
#> 3 TYPE492 standard statistical regions                       standard statistic…
#> 4 TYPE496 pre-1996 local authority districts                 pre-1996 local aut…
#> 5 TYPE498 pre-1996 counties / scottish regions               pre-1996 counties …
#> 6 TYPE499 countries                                          countries
 
 jobseekers_sex <- nomis_get_metadata("NM_1_1", "sex", "TYPE")
 
 tibble::glimpse(jobseekers_sex)
#> Rows: 3
#> Columns: 4
#> $ id             <chr> "5", "6", "7"
#> $ parentCode     <chr> "7", "7", NA
#> $ label.en       <chr> "Male", "Female", "Total"
#> $ description.en <chr> "Male", "Female", "Total"
 
 z <- nomis_get_data(id = "NM_1_1", time = "latest", geography = "TYPE499",
                     measures=c(20100, 20201), sex=5)
#> No encoding supplied: defaulting to UTF-8.
 
 tibble::glimpse(z)
#> Rows: 70
#> Columns: 34
#> $ DATE                <chr> "2021-12", "2021-12", "2021-12", "2021-12", "2021-…
#> $ DATE_NAME           <chr> "December 2021", "December 2021", "December 2021",…
#> $ DATE_CODE           <chr> "2021-12", "2021-12", "2021-12", "2021-12", "2021-…
#> $ DATE_TYPE           <chr> "date", "date", "date", "date", "date", "date", "d…
#> $ DATE_TYPECODE       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ DATE_SORTORDER      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ GEOGRAPHY           <dbl> 2092957697, 2092957697, 2092957697, 2092957697, 20…
#> $ GEOGRAPHY_NAME      <chr> "United Kingdom", "United Kingdom", "United Kingdo…
#> $ GEOGRAPHY_CODE      <chr> "K02000001", "K02000001", "K02000001", "K02000001"…
#> $ GEOGRAPHY_TYPE      <chr> "countries", "countries", "countries", "countries"…
#> $ GEOGRAPHY_TYPECODE  <dbl> 499, 499, 499, 499, 499, 499, 499, 499, 499, 499, …
#> $ GEOGRAPHY_SORTORDER <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,…
#> $ SEX                 <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,…
#> $ SEX_NAME            <chr> "Male", "Male", "Male", "Male", "Male", "Male", "M…
#> $ SEX_CODE            <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,…
#> $ SEX_TYPE            <chr> "sex", "sex", "sex", "sex", "sex", "sex", "sex", "…
#> $ SEX_TYPECODE        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ SEX_SORTORDER       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ ITEM                <dbl> 1, 1, 2, 2, 3, 3, 4, 4, 9, 9, 1, 1, 2, 2, 3, 3, 4,…
#> $ ITEM_NAME           <chr> "Total claimants", "Total claimants", "Students on…
#> $ ITEM_CODE           <dbl> 1, 1, 2, 2, 3, 3, 4, 4, 9, 9, 1, 1, 2, 2, 3, 3, 4,…
#> $ ITEM_TYPE           <chr> "item", "item", "item", "item", "item", "item", "i…
#> $ ITEM_TYPECODE       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ ITEM_SORTORDER      <dbl> 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 1, 1, 2, 2, 3,…
#> $ MEASURES            <dbl> 20100, 20201, 20100, 20201, 20100, 20201, 20100, 2…
#> $ MEASURES_NAME       <chr> "Persons claiming JSA", "Workplace-based estimates…
#> $ OBS_VALUE           <dbl> 73931.0, 0.3, NA, NA, NA, NA, NA, NA, NA, NA, 6791…
#> $ OBS_STATUS          <chr> "A", "A", "Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q", …
#> $ OBS_STATUS_NAME     <chr> "Normal Value", "Normal Value", "These figures are…
#> $ OBS_CONF            <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, F…
#> $ OBS_CONF_NAME       <chr> "Free (free for publication)", "Free (free for pub…
#> $ URN                 <chr> "Nm-1d1d32348e0d2092957697d5d1d20100", "Nm-1d1d323…
#> $ RECORD_OFFSET       <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, …
#> $ RECORD_COUNT        <dbl> 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70…
```

There is a lot of data available through Nomis, and there are some
limits to the amount of data that can be retrieved within a certain
period of time, although those are not published. For more details, see
the [full API documentation](https://www.nomisweb.co.uk/api/v01/help)
from Nomis. Full package documentation is available at
[docs.evanodell.com/nomisr](https://docs.evanodell.com/nomisr).

## Meta

Bug reports, suggestions, and code contributions are all welcome. Please
see
[CONTRIBUTING.md](https://github.com/ropensci/nomisr/blob/master/CONTRIBUTING.md)
for details.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/nomisr/blob/master/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.

Please note that this project is not affiliated with the Office for
National Statistics or the University of Durham (who run Nomis on behalf
of the Office for National Statistics).

Please use the reference below when citing `nomisr`, or use
`citation(package = 'nomisr')`.

Odell, (2018). nomisr: Access ‘Nomis’ UK Labour Market Data. *Journal of
Open Source Software*, 3(27), 859, doi:
[10.21105/joss.00859](https://doi.org/10.21105/joss.00859).

A BibTeX entry for LaTeX users is

    @Article{odell2018,
        title = {{nomisr}: Access Nomis UK Labour Market Data With R},
        volume = {3},
        doi = {10.21105/joss.00859},
        number = {27},
        journal = {The Journal of Open Source Software},
        year = {2018},
        month = {July},
        pages = {859},
        author = {Evan Odell},
      }

License:
[MIT](https://github.com/ropensci/nomisr/blob/master/LICENSE.md)

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


# nomisr 0.4.5

* Better error message when API returns empty data in some circumstances.

* Suppressed printing of column types when reading CSV files (#25, thanks jackobailey)

* Removed 2nd vignette due to errors 


# nomisr 0.4.4

* replaced `as.tibble` with `as_tibble`

* Includes second vignette "Work and Health Indicators with nomisr"

# nomisr 0.4.3

* Fix readme links

# nomisr 0.4.2

* Error handling improvements when using non-existent parameters, and clarifies
  error messages when no data is available for a given query.
  
* Removes redundant call to API (#19), thanks @Chrisjb

* New `tidy` parameter in `nomis_get_metadata()` to convert names to snake_case.

* Now using the `snakecase` package to implement name cleaning, 
  providing a broader range of naming styles.
  
* `nomis_get_metadata()` now makes existence of time concept explicit in the 
  tibble returned by `nomis_get_metadata({id})`.


# nomisr 0.4.1

* Adding `query_id` parameter to `nomis_get_data()`

* Changed documentation to use `roxygen` markdown

## Bug fixes

* Fixed bug where the `select` parameter in `nomis_get_data()` didn't work if 
  "OBS_VALUE" was not one of the variables. (@JimShady, #12)


# nomisr 0.4.0

* Version bump for CRAN

* Citation now refers to JOSS paper

* Some minor changes to internal code for easier maintenance

* Documentation updates to clarify difference between `time` and `date` 
parameters in `nomis_get_data()`


# nomisr 0.3.2 (non-CRAN release)

## New features and function changes

* The `tidy` parameter in `nomis_get_data()` now defaults to `FALSE` in order
to preserve existing workflows.

# nomisr 0.3.1 (non-CRAN release)

## New features and function changes

* `nomis_get_data()` now includes `tidy` and `tidy_style` parameters. 
`nomis_get_data()` now defaults to converting all variable names to 
"snake_case". `tidy_style` also accepts "camelCase" and "period.case" as name
style options.

## Bug fixes

* Dot queries to `nomis_get_data()` work with case-insensitive measurements 
more persistently. 

## Other updates

* Clarification of need to specified as `NULL` unused named parameters in
`nomis_get_data()` when using similarly named parameters in `...`.

# nomisr 0.3.0

## New features and function changes

* New `nomis_codelist()` function, which returns the internal coding for 
different concepts used by the NOMIS API in a `tibble`, given a dataset 
ID and a concept name. 

* The `additional_queries` parameter in `nomis_get_data()` and 
`nomis_get_metadata()` has been deprecated and will eventually be removed. 
Please use the `...` parameter for queries including concepts not available 
through the default parameters.

* The `sex` parameter in `nomis_get_data()` will also work with datasets that 
use "gender" instead of "sex".

## Internal changes and bug fixes

* Uses `rsdmx` to parse metadata, fixing #7.

# nomisr 0.2.0

* Improved API key handling (#5)

* Increased test coverage

* Adding rOpenSci reviewers to DESCRIPTION file.


# nomisr 0.1.0

* Moved to rOpenSci github repository

* Added API key functionality, which is not required by the API but is 
useful for large requests.

* In interactive sessions, users are asked if they want to continue when 
calling more than 15 pages of data at a time.

# nomisr 0.0.2

* Introduction of additional parameters to the `nomis_get_data()` and 
`nomis_codes()` functions, improvements to documentation.

# nomisr 0.0.1

* 1st release. Rudimentary functions for retrieving information on available 
datasets and downloading datasets from nomis, with some limited parameters.

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
(https:contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
# Contributing to `nomisr`

Thank you for any and all contributions! Following these guidelines will help streamline the process of contributing and make sure that we're all on the same page. While we ask that you read this guide and follow it to the best of your abilities, we welcome contributions from all, regardless of your level of experience.

By participating in this project, you agree to abide by the [code of conduct](https://github.com/ropensci/nomisr/blob/master/CONDUCT.md).

# Types of contributions 

Don't feel that you must be a computer whiz to make meaningful contributions. Feel free to:

- Identify areas for future development ([open an Issue](https://github.com/ropensci/nomisr/issues))
- Identify issues/bugs ([open an Issue](https://github.com/ropensci/nomisr/issues))
- Write tutorials/vignettes ([open a Pull Request](https://github.com/ropensci/nomisr/pulls) to contribute to the ones here, or make your own elsewhere and send us a link)
- Add functionality ([open a Pull Request](https://github.com/ropensci/nomisr/pulls))
- Fix bugs ([open a Pull Request](https://github.com/ropensci/nomisr/pulls))

# New to GitHub?

Getting ready to make your first contribution? Here are a couple of tutorials you may wish to check out:

- [Tutorial for first-timers](https://github.com/Roshanjossey/first-contributions)
- [How to contribute (in-depth lessons)](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github)
- [GitHub on setup](https://help.github.com/articles/set-up-git)
- [GitHub on pull requests](https://help.github.com/articles/using-pull-requests/).)


# How to contribute code

- Fork the repository
- Clone the repository from GitHub to your computer e.g,. `git clone https://github.com/ropensci/nomisr.git`
- Make sure to track progress upstream (i.e., on our version of `nomisr` at `ropensci/nomisr`)
  - `git remote add upstream https://github.com/ropensci/nomisr.git`
  - Before making changes make sure to pull changes in from upstream with `git pull upstream`
- Make your changes
  - For changes beyond minor typos, add an item to NEWS.md describing the changes and add yourself to the DESCRIPTION file as a contributor
- Push to your GitHub account
- Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/nomisr`

# Code formatting

- In general follow the convention of <https://r-pkgs.had.co.nz/r.html#style> (snake_case functions and argument names, etc.)
- Where there is conflict, default to the style of `nomisr`
- Use explicit package imports (i.e. package_name::package_function) and avoid @import if at all possible

## Release summary

This is an resubmission of the `nomisr` package, version number 0.4.5. It 
features improved error messages and makes output less verbose when reading 
CSV files. It also removes a vignette that was causing erros due to reliance on 
old package versions.

## Test environments
* local macOS install, R 4.1.2
* win-builder (devel and release)
* Windows Server 2019, release (on GitHub Actions)
* macOS Catalina 10.15, release (on GitHub Actions)
* ubuntu 20.04 (devel and release) (on GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 0 notes
---
title: 'nomisr: Access ''Nomis'' UK Labour Market Data'
authors:
- affiliation: '1'
  name: Evan Odell
  orcid: 0000-0003-1845-808X
date: "2018-07-27"
output:
  html_document:
    keep_md: yes
bibliography: paper.bib
tags:
- R
- geography
- official statistics
affiliations:
- index: 1
  name: Disability Rights UK
---



The University of Durham runs the Nomis database of labour market statistics on behalf of the UK's Office for National Statistics [-@ons1981]. As of publication, Nomis contains 1,249 datasets, almost all of which are based around differing statistical geographies. All the data is freely available and does not require users to create accounts to download data, and Nomis provides an interactive web platform for downloading data. However, like all GUI downloading systems, there is a risk that users will select the wrong option without realising their mistake, and downloading multiple datasets is tedious, repetitive work. The `nomisr` package provides functions to identify datasets available through Nomis, the  variables and query options for those datasets, and a function for downloading data, including combining large datasets into a single `tibble` [@muller2018]

`nomisr` is designed around a three stage workflow, with functions for each stage:

1. Identifying available datasets, using the `nomis_data_info()` without any parameters to return a tibble with the names and basic metadata of all available datasets, or the `nomis_search()` function to retrieve a tibble of datasets matching a given search term.

1. Identifying metadata, including "concepts", the name Nomis uses for variables that can be specified when retrieving data. This is done using the `nomis_get_metadata()` function.

1. Downloading data, using the `nomis_get_data()` function, which requires the ID of a given dataset, and accepts parameters specifying the geographic concept (either all geographies of a given area type or one or more specific geographic areas of a given type) and any other concepts used to specify queries. 

`nomisr` is able to return specific releases of a given dataset (to aid reproducible research) or the most recent available data. The `nomis_get_data()` function includes common parameters built into the function, and also accepts unquoted concepts and concept values using quasiquotation [@henry2018], in order to accomodate the wide range of concepts used by different Nomis datasets.

Data downloaded with `nomisr` is in a `tibble` that is ready for plotting (see Figure 1) or mapping with other R packages. Specifying geographies in `nomisr` requires using Nomis' internal geography coding, but all data downloads include standard ONS geography codes to aid mapping.

![](paper_files/figure-html/plot-demo-1.png)<!-- -->
Figure 1. Data retrieved with `nomisr` is ready for plotting.


`nomisr` is available on GitHub at <https://github.com/ropensci/nomisr>.

# Acknowledgements

I thank Paul Egeler for his contributions to `nomisr` and his code review comments, and Christophe Dervieux for his code review comments.


# References
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

# nomisr

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/nomisr)](https://cran.r-project.org/package=nomisr)
[![GitHub tag](https://img.shields.io/github/tag/ropensci/nomisr.svg)](https://github.com/ropensci/nomisr)
[![](https://cranlogs.r-pkg.org/badges/grand-total/nomisr)](https://cran.r-project.org/package=nomisr)
[![R build status](https://github.com/ropensci/nomisr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/nomisr/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/nomisr/master.svg)](https://codecov.io/github/ropensci/nomisr?branch=master)
[![ropensci](https://badges.ropensci.org/190_status.svg)](https://github.com/ropensci/software-review/issues/190)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1157908.svg)](https://doi.org/10.5281/zenodo.1157908)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.00859/status.svg)](https://doi.org/10.21105/joss.00859)
<!-- badges: end -->

`nomisr` is for accessing UK official statistics from the [Nomis](https://www.nomisweb.co.uk/) database through R. Nomis contains data from the Census, the Labour Force Survey, DWP benefit statistics and other economic and demographic data, and is maintained on behalf of the Office for National Statistics by the University of Durham.

The `nomisr` package provides functions to find what data is available, the  variables and query options for different datasets and a function for downloading data. `nomisr` returns data in [`tibble`](https://cran.r-project.org/package=tibble) format. Most of the data available through `nomisr` is based around statistical geographies, with a handful of exceptions.

The package is for demographers, economists, geographers, public health researchers and any other researchers who are interested in geographic factors. The package aims to aid reproducibility, reduce the need to manually download area profiles, and allow easy linking of different datasets covering the same geographic area.

## Installation

`nomisr` is available on CRAN:

```{r cran-installation, eval = FALSE}
install.packages("nomisr")
```

You can install the development version `nomisr` from github with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/nomisr")
```

## Using `nomisr`

`nomisr` contains functions to search for datasets, identify the query options for different datasets and retrieve data from queries, all done with [`tibbles`](https://tibble.tidyverse.org/), to take advantage of how `tibble` manages list-columns. The use of metadata queries, rather than simply downloading all available data, is useful to avoid overwhelming the rate limits of the API. 

There are `nomisr` vignette [introduction](https://docs.evanodell.com/nomisr/articles/introduction.html) has details on all available functions and basic demonstrations of their use.

The example below demonstrates a workflow to retrieve the latest data on Jobseeker's Allowance with rates and proportions, on a national level, with all male claimants and workforce.

```{r example}
 library(nomisr)
 jobseekers_search <- nomis_search(name = "*Jobseeker*")
 
 tibble::glimpse(jobseekers_search)

 jobseekers_measures <- nomis_get_metadata("NM_1_1", "measures")
 
 tibble::glimpse(jobseekers_measures)
 
 jobseekers_geography <- nomis_get_metadata("NM_1_1", "geography", "TYPE")
 
 tail(jobseekers_geography)
 
 jobseekers_sex <- nomis_get_metadata("NM_1_1", "sex", "TYPE")
 
 tibble::glimpse(jobseekers_sex)
 
 z <- nomis_get_data(id = "NM_1_1", time = "latest", geography = "TYPE499",
                     measures=c(20100, 20201), sex=5)
 
 tibble::glimpse(z)
```

There is a lot of data available through Nomis, and there are some limits to the amount of data that can be retrieved within a certain period of time, although those are not published. For more details, see the [full API documentation](https://www.nomisweb.co.uk/api/v01/help) from Nomis. Full package documentation is available at [docs.evanodell.com/nomisr](https://docs.evanodell.com/nomisr).


## Meta

Bug reports, suggestions, and code contributions are all welcome. Please see [CONTRIBUTING.md](https://github.com/ropensci/nomisr/blob/master/CONTRIBUTING.md) for details.

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/nomisr/blob/master/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

Please note that this project is not affiliated with the Office for National Statistics or the University of Durham (who run Nomis on behalf of the Office for National Statistics).

Please use the reference below when citing `nomisr`, or use `citation(package = 'nomisr')`.

Odell, (2018). nomisr: Access 'Nomis' UK Labour Market Data. _Journal of Open Source Software_, 3(27), 859, doi: [10.21105/joss.00859](https://doi.org/10.21105/joss.00859).

A BibTeX entry for LaTeX users is
```
@Article{odell2018,
    title = {{nomisr}: Access Nomis UK Labour Market Data With R},
    volume = {3},
    doi = {10.21105/joss.00859},
    number = {27},
    journal = {The Journal of Open Source Software},
    year = {2018},
    month = {July},
    pages = {859},
    author = {Evan Odell},
  }
```
License: [MIT](https://github.com/ropensci/nomisr/blob/master/LICENSE.md)

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: 'nomisr: Access ''Nomis'' UK Labour Market Data'
authors:
- affiliation: '1'
  name: Evan Odell
  orcid: 0000-0003-1845-808X
date: "`r Sys.Date()`"
output:
  html_document:
    keep_md: yes
bibliography: paper.bib
tags:
- R
- geography
- official statistics
affiliations:
- index: 1
  name: Disability Rights UK
---

```{r options, message=FALSE, warning=FALSE, include=FALSE}
library(nomisr)
library(ggplot2)
library(dplyr)

knitr::opts_chunk$set(cache = FALSE)

options(width = 90, tibble.max_extra_cols = 0)
```

The University of Durham runs the Nomis database of labour market statistics on behalf of the UK's Office for National Statistics [-@ons1981]. As of publication, Nomis contains 1,249 datasets, almost all of which are based around differing statistical geographies. All the data is freely available and does not require users to create accounts to download data, and Nomis provides an interactive web platform for downloading data. However, like all GUI downloading systems, there is a risk that users will select the wrong option without realising their mistake, and downloading multiple datasets is tedious, repetitive work. The `nomisr` package provides functions to identify datasets available through Nomis, the  variables and query options for those datasets, and a function for downloading data, including combining large datasets into a single `tibble` [@muller2018]

`nomisr` is designed around a three stage workflow, with functions for each stage:

1. Identifying available datasets, using the `nomis_data_info()` without any parameters to return a tibble with the names and basic metadata of all available datasets, or the `nomis_search()` function to retrieve a tibble of datasets matching a given search term.

1. Identifying metadata, including "concepts", the name Nomis uses for variables that can be specified when retrieving data. This is done using the `nomis_get_metadata()` function.

1. Downloading data, using the `nomis_get_data()` function, which requires the ID of a given dataset, and accepts parameters specifying the geographic concept (either all geographies of a given area type or one or more specific geographic areas of a given type) and any other concepts used to specify queries. 

`nomisr` is able to return specific releases of a given dataset (to aid reproducible research) or the most recent available data. The `nomis_get_data()` function includes common parameters built into the function, and also accepts unquoted concepts and concept values using quasiquotation [@henry2018], in order to accomodate the wide range of concepts used by different Nomis datasets.

Data downloaded with `nomisr` is in a `tibble` that is ready for plotting (see Figure 1) or mapping with other R packages. Specifying geographies in `nomisr` requires using Nomis' internal geography coding, but all data downloads include standard ONS geography codes to aid mapping.

```{r plot-demo, echo=FALSE, message=FALSE, warning=FALSE, dpi = 600}

mort_data <- nomis_get_data(id = "NM_161_1", date = "2016",
                            geography = "TYPE480", 
                            cause_of_death = "440", 
                            sex = 0, age = 0, MEASURE = 6)

mort_data <- mort_data %>%
  arrange(OBS_VALUE) %>%               # sort your dataframe
  mutate(GEOGRAPHY_NAME = factor(GEOGRAPHY_NAME, unique(GEOGRAPHY_NAME))) #

blood_cancer <- ggplot(mort_data, aes(x = GEOGRAPHY_NAME, 
                                      y = OBS_VALUE, fill = OBS_VALUE)) +
  geom_col() + 
  scale_fill_viridis_c(option = "plasma", name = "% of Deaths") + 
  labs(title = "Percentage of Total Deaths from Blood Cancers",
       subtitle = "Regions of England and Wales, 2016") + 
  xlab("Region") + 
  ylab("Percentage of Deaths") + 
  theme(legend.position="none",
        axis.text.x = element_text(angle=30, hjust=1))

blood_cancer
```
Figure 1. Data retrieved with `nomisr` is ready for plotting.


`nomisr` is available on GitHub at <https://github.com/ropensci/nomisr>.

# Acknowledgements

I thank Paul Egeler for his contributions to `nomisr` and his code review comments, and Christophe Dervieux for his code review comments.


# References
---
title: "Introduction to nomisr"
author: "Evan Odell"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        theme: flatly
vignette: >
  %\VignetteIndexEntry{Introduction to nomisr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction to the Nomis API

`nomisr` is for accessing [Nomis](https://www.nomisweb.co.uk/) data with R. The Nomis API is free to access without registration, and contains up-to-date official statistics, including data from the Census, the Labour Force Survey, DWP benefit statistics and other economic and demographic data. Nomis is maintained on behalf of the Office for National Statistics by the University of Durham.

There is a lot of data available through Nomis, and there are some limits to the amount of data that can be retrieved within a certain period of time, although those are not published. For more details, see the [full API documentation](https://www.nomisweb.co.uk/api/v01/help) from Nomis.

Nomis data is based around administrative and statistical geographies, and a particular geography should be specified when downloading data.

`nomisr` is designed around a pipeline of three key functions: `nomis_data_info()`, `nomis_get_metadata()` and `nomis_get_data()`. The `nomis_overview()`, `nomis_content_type()` and `nomis_search()` functions can assist with this.

## Querying data availability

The `nomis_data_info()` function is focused on the structure and coverage of the available datasets.

Use the `nomis_data_info()` function without any parameters to get a tibble with metadata for all available datasets:

```
x <- nomis_data_info()
head(x)
```


`nomis_data_info()` can also be used to query metadata from a specific dataset, using its ID. The example below uses the "LC4408EW - Tenure by number of persons per bedroom in household by household type" dataset from the 2011 census, which has the ID "NM_893_1". 

```{r specific-dataset, echo=TRUE}
library(nomisr)
y <- nomis_data_info("NM_893_1")

tibble::glimpse(y)
```

When a tibble with metadata for all datasets or a specific dataset is returned, three of the columns, `annotations.annotation`, `components.attribute` and `components.dimension`, are list-columns of data frames. `annotations.annotation` contains metadata on the dataset, including units and current status. `components.attribute` contains more detailed status metadata. `components.dimension` contains the grouping and summary variables available in the dataset, which vary between different datasets.

The example below shows how to access data stored in list columns returned from the Nomis API. In the case of requests for metadata from a single dataset, the three columns are all lists with a length of 1. If requesting all dataset information with `nomis_data_info()`, each row is a list of length 1. Each list contains a data.frame, of varrying dimensions depending on the column and dataset. You can unnest individual list-columns to display their data in the same row as data from the rest of the tibble. Due to the differing lengths of the list-columns returned by `nomis_data_info()`, only one list-column can be unnested at a time.

```{r specific-dataset-exam, echo=TRUE}
library(dplyr, warn.conflicts = F)

y$annotations.annotation %>% class()

y$annotations.annotation %>% length()

y$annotations.annotation[[1]] %>% class()

y %>% pull(annotations.annotation) %>% class()

y %>% pull(annotations.annotation) %>% .[[1]] %>% class()

y %>% pull(annotations.annotation) %>% purrr::pluck() %>% class()

## Unnesting list columns
y %>% tidyr::unnest(annotations.annotation) %>% glimpse()
```

### Searching for data

`nomisr` also contains the `nomis_search()` function to search for datasets on particular topics. `nomis_search()` can be used to search in one or more of dataset names, descriptions, keywords, content type and units. If using multiple parameters, `nomis_search()` will return information on all datasets that match one or more parameters. Character vectors of strings can be used in searches, and likewise `nomis_search()` will return information on datasets that match one or more queries. The * is used as a wildcard symbol. `nomis_search()` returns metadata in the same format as `nomis_data_info()`, including using list-columns. The `nomis_content_type()` function can assist in identifying content type IDs for `nomis_search()`.

```{r data-searching, echo=TRUE}
a <- nomis_search(name = '*jobseekers*', keywords = 'Claimants')

tibble::glimpse(a)

a %>% tidyr::unnest(components.attribute) %>% glimpse()

b <- nomis_search(keywords = c('Claimants', '*Year*'))

tibble::glimpse(b)

b %>% tidyr::unnest(components.attribute) %>% glimpse()

```


### Other ways to access metadata

`nomis_overview()` returns a tibble with a generalised overview of a given dataset.


```{r overview, echo=TRUE}
q <- nomis_overview("NM_1650_1")

q %>% tidyr::unnest(name) %>% glimpse()

```

`nomis_overview()` has a `select` parameter that can be used to select only particular elements of the overview to return.

```{r overview-select, echo=TRUE}
s <- nomis_overview("NM_1650_1", select = c("units", "keywords"))
 
s %>% tidyr::unnest(name) %>% glimpse()
```


## Querying data variables

Vast amounts of data are available through Nomis and so to avoid overwhelming the API, it is good practice to query what concepts are available, using `nomis_get_metadata()`. While the other metadata functions can return concept metadata, `nomis_get_metadata()` provides greater flexibility and specificity over the returned metadata than `nomis_overview()` and `nomis_data_info()`.

The example below queries some of the metadata available through the API for the "LC4408EW - Tenure bynumber of persons per bedroom in household by household type" dataset.

### Getting concepts

If provided with just a dataset ID, `nomis_get_metadata()` will return the concepts available for the given dataset. 

```{r get-metadata, echo=TRUE}
a <- nomis_get_metadata(id = "NM_893_1")

a
```

### Concept Values
  
If provided with a concept name it returns the available values for that concept. However, in some cases, espescially with the geography concept, there are multiple options available, which Nomis labels types. In that case `nomis_get_metadata()` returns the values of the lowest indexed type available.
  
```
b <- nomis_get_metadata(id = "NM_893_1", concept = "GEOGRAPHY")
```

We can now pass a generic "type" string to the `type` parameter in `nomis_get_metadata()`, which returns all available geography types for dataset "NM_893_1".

  
```
c <- nomis_get_metadata(id = "NM_893_1", concept = "geography", type = "type")
```

  
Passing a specific type to the `type` parameter, in this case "TYPE460" for all post-2010 parliamentary constituencies, returns a tibble with geographic codes for those specific constituencies, which can be used to filter queries.
  
```
d <- nomis_get_metadata(id = "NM_893_1", concept = "geography", type = "TYPE460")

```


The vast majority (98% as of February 2018) of Nomis datasets include a geographic variable.

## Downloading data

Using the information above, we can now query the latest data on bedroom occupancy per household type in different NHS clinical commissioning groups.

```
z <- nomis_get_data(id = "NM_893_1", time = "latest", geography = "TYPE266")
```

We can also query bedroom occupancy per household type in the Manchester, Gorton and Manchester, Withington parliamentary constituencies.

```
x <- nomis_get_data(id = "NM_893_1", time = "latest",
                    geography = c("1929380119", "1929380120"))
```


`nomisr` also allows for time series queries. The example below shows how to retrieve the percentage of the workforce claiming Jobseekers Allowance from January 2015 to January 2020, inclusive, for each region of the UK, divided by male and female claimants, with an accompanying graph.


```{r jsa-claimaints, eval=FALSE}
library(ggplot2)
library(dplyr)
library(nomisr)

jsa <- nomis_get_data(id = "NM_1_1", time = "2018-01-2021-10",
                      geography = "TYPE480", measures=20201,
                      sex=c(5,6), item = 1, tidy = TRUE)

jsa <- jsa %>% 
  mutate(date = as.Date(paste0(date, "-01")),
         obs_value = obs_value/100)

theme_set(theme_bw())

p_jsa <- ggplot(jsa, aes(x = date, y = obs_value, colour = sex_name)) + 
  geom_line(size = 1.15) +
  scale_colour_viridis_d(end = 0.75, begin = 0.1, name = "Gender") + 
  scale_x_date(breaks = "6 months", date_labels = "%b %Y") + 
  scale_y_continuous(labels = scales::percent) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom") + 
  labs(x = "Date", y= "JSA Claimants (Percentage of Workforce)") + 
  facet_wrap(~geography_name, scales = "free_y")

p_jsa
```

![](p_jsa.png)


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codelist.R
\name{nomis_codelist}
\alias{nomis_codelist}
\title{Nomis codelists}
\usage{
nomis_codelist(id, concept, search = NULL)
}
\arguments{
\item{id}{A string with the ID of the particular dataset. Must be specified.}

\item{concept}{A string with the variable concept to return options for. If
left empty, returns all the variables for the dataset specified by \code{id}.
Codes are not case sensitive and must be specified.}

\item{search}{Search for codes that contain a given string. The wildcard
character \code{*} can be added to the beginning and/or end of each
search string. Search strings are not case sensitive.
Defaults to \code{NULL}. Note that the search function is not very powerful
for some datasets.}
}
\value{
A tibble with the codes used to query specific concepts.
}
\description{
Nomis uses its own internal coding for query concepts. \code{nomis_codelist}
returns the codes for a given concept in a \code{tibble}, given a dataset
ID and a concept name.
Note that some codelists, particularly geography, can be very large.
}
\examples{
\donttest{
x <- nomis_codelist("NM_1_1", "item")


# Searching for codes ending with "london"
y <- nomis_codelist("NM_1_1", "geography", search = "*london")


z <- nomis_codelist("NM_161_1", "cause_of_death")
}

}
\seealso{
\code{\link[=nomis_data_info]{nomis_data_info()}}

\code{\link[=nomis_get_metadata]{nomis_get_metadata()}}

\code{\link[=nomis_overview]{nomis_overview()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nomisr-package.R
\docType{package}
\name{nomisr}
\alias{nomisr-package}
\title{nomisr: Access Nomis UK Labour Market Data with R}
\description{
Access UK official statistics from the Nomis database. Nomis
includes data from the Census, the Labour Force Survey, DWP benefit
statistics and other economic and demographic data from the Office for
National Statistics.
}
\details{
The package provides functions to find what data is available, metadata,
including the variables and query options for different datasets and
a function for downloading data.

The full API documentation and optional registration for an API key is
available at \url{https://www.nomisweb.co.uk/api/v01/help}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api-key.R
\name{nomis_api_key}
\alias{nomis_api_key}
\title{Nomis API Key}
\usage{
nomis_api_key(check_env = FALSE)
}
\arguments{
\item{check_env}{If TRUE, will check the environment variable
\code{NOMIS_API_KEY} first before asking for user input.}
}
\description{
Assign or reassign API key for Nomis.
}
\details{
The Nomis API has an optional key. Using the key means that 100,000
rows can be returned per call, which can speed up larger data requests and
reduce the chances of being rate limited or having requests timing out.

By default, \code{nomisr} will look for the environment variable
\code{NOMIS_API_KEY} when the package is loaded. If found, the API key will
be stored in the session option \code{nomisr.API.key}. If you would like to
reload the API key or would like to manually enter one in, this function
may be used.

You can sign up for an API key
\href{https://www.nomisweb.co.uk/myaccount/userjoin.asp}{here}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{nomis_search}
\alias{nomis_search}
\title{Search Nomis datasets}
\usage{
nomis_search(
  name = NULL,
  description = NULL,
  keywords = NULL,
  content_type = NULL,
  units = NULL,
  tidy = FALSE
)
}
\arguments{
\item{name}{A string or character vector of strings to search for in
available dataset names. Defaults to \code{NULL}.}

\item{description}{A string or character vector of strings to search for in
available dataset descriptions. Note that \code{description} looks for
complete matches, so wildcards should be used at the start and end of
each string. Defaults to \code{NULL}.}

\item{keywords}{A string or character vector of strings to search for in
available dataset keywords. Defaults to \code{NULL}.}

\item{content_type}{A string or character vector of strings to search for
in available dataset content types. \code{content_type} can include an
optional ID for that content type. Defaults to \code{NULL}.}

\item{units}{A string or character vector of strings to search for in
available dataset units. Defaults to \code{NULL}.}

\item{tidy}{If \code{TRUE}, converts tibble names to snakecase.}
}
\value{
A tibble with details on all datasets matching the search query.
}
\description{
A function to search for datasets on given topics. In the case of multiple
search parameters, returns metadata on all datasets matching  one or more
of the parameters. The wildcard character \code{*} can be added to the
beginning and/or end of each search string.
}
\examples{
\donttest{
x <- nomis_search(name = "*seekers*")

y <- nomis_search(keywords = "Claimants")

# Return metadata of all datasets with content_type "sources".
a <- nomis_search(content_type = "sources")


# Return metadata of all datasets with content_type "sources" and
# source ID "acses"
b <- nomis_search(content_type = "sources-acses")
}

}
\seealso{
\code{\link[=nomis_content_type]{nomis_content_type()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_download.R
\name{nomis_get_data}
\alias{nomis_get_data}
\title{Retrieve Nomis datasets}
\usage{
nomis_get_data(
  id,
  time = NULL,
  date = NULL,
  geography = NULL,
  sex = NULL,
  measures = NULL,
  additional_queries = NULL,
  exclude_missing = FALSE,
  select = NULL,
  tidy = FALSE,
  tidy_style = "snake_case",
  query_id = NULL,
  ...
)
}
\arguments{
\item{id}{A string containing the ID of the dataset to retrieve,
in \code{"nm_***_*"} format. The \code{id} parameter is not case sensitive.}

\item{time}{Parameter for selecting dates and date ranges. Accepts either a
single date value, or two date values and returns all data between the two
date values, There are two styles of values that can be used to query time.

The first is one or two of \code{"latest"} (returns the latest available
data), \code{"previous"} (the date prior to \code{"latest"}), \code{"prevyear"}
(the date one year prior to \code{"latest"}) or \code{"first"}
(the oldest available data for the dataset).

The second style is to use or a specific date or multiple dates, in the
style of the time variable codelist, which can be found using the
\code{\link[=nomis_get_metadata]{nomis_get_metadata()}} function.

Values for the \code{time} and \code{date} parameters should not be used
at the same time. If they are, the function will retrieve data based
on the the \code{date} parameter.

Defaults to \code{NULL}.}

\item{date}{Parameter for selecting specific dates. Accepts one or more date
values. If given multiple values, only data for the given dates will be
returned, but there is no limit to the number of data values. For example,
\code{date=c("latest, latestMINUS3, latestMINUS6")} will return the latest
data, data from three months prior to the latest data and six months prior
to the latest data. There are two styles of values that can be used to
query time.

The first is one or more of \code{"latest"} (returns the latest available
data), \code{"previous"} (the date prior to \code{"latest"}),
\code{"prevyear"} (the date one year prior to \code{"latest"}) or
\code{"first"} (the oldest available data for the dataset).

The second style is to use or a specific date or multiple dates, in the
style of the time variable codelist, which can be found using the
\code{\link[=nomis_get_metadata]{nomis_get_metadata()}} function.

Values for the \code{time} and \code{date} parameters should not be used at
the same time. If they are, the function will retrieve data based on the
the \code{date} parameter.

Defaults to \code{NULL}.}

\item{geography}{The code of the geographic area to return data for. If
\code{NULL}, returns data for all available geographic areas, subject to
other parameters. Defaults to \code{NULL}. In the rare instance that a
geographic variable does not exist, if not \code{NULL}, the function
will return an error.}

\item{sex}{The code for sexes/genders to include in the dataset.
Accepts a string or number, or a vector of strings or numbers.
\code{nomisr} automatically voids any queries for sex if it is not an
available code in the requested dataset. Defaults to \code{NULL} and returns
all available sex/gender data.

There are two different codings used for sex, depending on the dataset. For
datasets using \code{"SEX"}, \code{7} will return results for males and females, \code{6}
only females and \code{5} only males. Defaults to \code{NULL}, equivalent to
\code{c(5,6,7)} for datasets where sex is an option. For datasets using
\code{"C_SEX"}, \code{0} will return results for males and females, \code{1} only males
and \code{2} only females. Some datasets use \code{"GENDER"} with the same values
as \code{"SEX"}, which works with both \verb{sex = <code>} and \verb{gender = <code>}
as a dot parameter.}

\item{measures}{The code for the statistical measure(s) to include in the
data. Accepts a single string or number, or a list of strings or numbers.
If \code{NULL}, returns data for all available statistical measures subject
to other parameters. Defaults to \code{NULL}.}

\item{additional_queries}{Any other additional queries to pass to the API.
See \url{https://www.nomisweb.co.uk/api/v01/help} for instructions on
query structure. Defaults to \code{NULL}. Deprecated in package versions greater
than 0.2.0 and will eventually be removed in a future version.}

\item{exclude_missing}{If \code{TRUE}, excludes all missing values.
Defaults to \code{FALSE}.}

\item{select}{A character vector of one or more variables to include in
the returned data, excluding all others. \code{select} is not case sensitive.}

\item{tidy}{Logical parameter. If \code{TRUE}, converts variable names to
\code{snake_case}, or another style as specified by the
\code{tidy_style} parameter. Defaults to \code{FALSE}. The default variable name style
from the API is SCREAMING_SNAKE_CASE.}

\item{tidy_style}{The style to convert variable names to, if
\code{tidy = TRUE}. Accepts one of \code{"snake_case"}, \code{"camelCase"}
and \code{"period.case"}, or any of the \code{case} options accepted by
\code{\link[snakecase:to_any_case]{snakecase::to_any_case()}}. Defaults to \code{"snake_case"}.}

\item{query_id}{Results can be labelled as belonging to a certain query
made to the API. \code{query_id} accepts any value as a string, and will
be included in every row of the tibble returned by \code{nomis_get_data}
in a column labelled "QUERY_ID" in the default SCREAMING_SNAKE_CASE
used by the API. Defaults to \code{NULL}.}

\item{...}{Use to pass any other parameters to the API. Useful for passing
concepts that are not available through the default parameters. Only accepts
concepts identified in \code{\link[=nomis_get_metadata]{nomis_get_metadata()}} and concept values
identified in \code{\link[=nomis_codelist]{nomis_codelist()}}. Parameters can be quoted or
unquoted. Each parameter should have a name and a value. For example,
\code{CAUSE_OF_DEATH = 10300} when querying dataset \code{"NM_161_1"}.
Parameters are not case sensitive. Note that R using partial matching for
function variables, and so passing a parameter with the same opening
characters as one of the above-named parameters can cause an error unless
the value of the named parameter is specified, including as \code{NULL}.
See example below:}
}
\value{
A tibble containing the selected dataset. By default, all tibble
columns except for the \code{"OBS_VALUE"} column are parsed as characters.
}
\description{
To find the code options for a given dataset, use
\code{\link[=nomis_get_metadata]{nomis_get_metadata()}} for specific codes, and
\code{\link[=nomis_codelist]{nomis_codelist()}} for code values.

This can be a slow process if querying significant amounts of
data. Guest users are limited to 25,000 rows per query, although
\code{nomisr} identifies queries that will return more than 25,000 rows,
sending individual queries and combining the results of those queries into
a single tibble. In interactive sessions, \code{nomisr} will warn you if
guest users are requesting more than 350,000 rows of data, and if
registered users are requesting more than 1,500,000 rows.

Note the difference between the \code{time} and \code{date}
parameters. The \code{time} and \code{date} parameters should not be used at the same
time. If they are, the function will retrieve data based on the the
\code{date} parameter. If given more than one query, \code{time} will
return all data available between those queries, inclusively, while
\code{date} will only return data for the exact queries specified. So
\code{time = c("first", "latest")} will return all data, while
\code{date = c("first", "latest")} will return only the earliest and latest
data published.
}
\examples{
\donttest{

# Return data on Jobseekers Allowance for each country in the UK
jobseekers_country <- nomis_get_data(
  id = "NM_1_1", time = "latest",
  geography = "TYPE499",
  measures = c(20100, 20201), sex = 5
)

# Return data on Jobseekers Allowance for Wigan
jobseekers_wigan <- nomis_get_data(
  id = "NM_1_1", time = "latest",
  geography = "1879048226",
  measures = c(20100, 20201), sex = "5"
)

# annual population survey - regional - employment by occupation
emp_by_occupation <- nomis_get_data(
  id = "NM_168_1", time = "latest",
  geography = "2013265925", sex = "0",
  select = c(
    "geography_code",
    "C_OCCPUK11H_0_NAME", "obs_vAlUE"
  )
)

# Deaths in 2016 and 2015 by three specified causes,
# identified with nomis_get_metadata()
death <- nomis_get_data("NM_161_1",
  date = c("2016", "2015"),
  geography = "TYPE480",
  cause_of_death = c(10300, 102088, 270)
)

# All causes of death in London in 2016
london_death <- nomis_get_data("NM_161_1",
  date = c("2016"),
  geography = "2013265927", sex = 1, age = 0
)
}
\dontrun{
# Results in an error because `measure` is mistaken for `measures`
mort_data1 <- nomis_get_data(
  id = "NM_161_1", date = "2016",
  geography = "TYPE464", sex = 0, cause_of_death = "10381",
  age = 0, measure = 6
)

# Does not error because `measures` is specified
mort_data2 <- nomis_get_data(
  id = "NM_161_1", date = "2016",
  geography = "TYPE464", sex = 0, measures = NULL,
  cause_of_death = "10381", age = 0, measure = 6
)
}

}
\seealso{
\code{\link[=nomis_data_info]{nomis_data_info()}}

\code{\link[=nomis_get_metadata]{nomis_get_metadata()}}

\code{\link[=nomis_codelist]{nomis_codelist()}}

\code{\link[=nomis_overview]{nomis_overview()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/content_type.R
\name{nomis_content_type}
\alias{nomis_content_type}
\title{Nomis dataset content types}
\usage{
nomis_content_type(content_type, id = NULL)
}
\arguments{
\item{content_type}{A string with the content type to return metadata on.}

\item{id}{A string with an optional \code{content_type} id.}
}
\value{
A tibble with metadata on a given content type.
}
\description{
Nomis content type metadata is included in annotation tags, in the form of
\verb{contenttype/<contenttype>} in the \code{annotationtitle} column in
the \code{annotations.annotation} list-column returned from
\code{\link[=nomis_data_info]{nomis_data_info()}}. For example, the content types returned from
dataset "NM_1658_1", using \code{nomis_data_info("NM_1658_1")}, are
"geoglevel", "2001census" and "sources".
}
\examples{
\donttest{
a <- nomis_content_type("sources")

tibble::glimpse(a)

b <- nomis_content_type("sources", id = "census")

tibble::glimpse(b)
}

}
\seealso{
\code{\link[=nomis_search]{nomis_search()}}

\code{\link[=nomis_data_info]{nomis_data_info()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{nomis_get_metadata}
\alias{nomis_get_metadata}
\title{Nomis metadata concepts and types}
\usage{
nomis_get_metadata(
  id,
  concept = NULL,
  type = NULL,
  search = NULL,
  additional_queries = NULL,
  ...,
  tidy = FALSE
)
}
\arguments{
\item{id}{The ID of the particular dataset. Returns no data if not specified.}

\item{concept}{A string with the variable concept to return options for. If
left empty, returns all the variables for the dataset specified by \code{id}.
Codes are not case sensitive. Defaults to \code{NULL}.}

\item{type}{A string with options for a particular code value, to return
types of variables available for a given code. Defaults to \code{NULL}. If
\code{concept == NULL}, \code{type} will be ignored.}

\item{search}{A string or character vector of strings to search for in the
metadata. Defaults to \code{NULL}. As in \code{\link[=nomis_search]{nomis_search()}}, the
wildcard character \code{*} can be added to the beginning and/or end of each
search string.}

\item{additional_queries}{Any other additional queries to pass to the API.
See \url{https://www.nomisweb.co.uk/api/v01/help} for instructions on
query structure. Defaults to \code{NULL}. Deprecated in package
versions greater than 0.2.0 and will eventually be removed.}

\item{...}{Use to pass any other parameters to the API.}

\item{tidy}{If \code{TRUE}, converts tibble names to snakecase.}
}
\value{
A tibble with metadata options for queries using \code{\link[=nomis_get_data]{nomis_get_data()}}.
}
\description{
Retrieve all concept code options of all Nomis datasets,
concept code options for a given dataset, or the all the options for a given
concept variable from a particular dataset. Specifying \code{concept} will
return all the options for a given variable in a particular dataset.

If looking for a more detailed overview of all available
metadata for a given dataset, see \code{\link[=nomis_overview]{nomis_overview()}}.
}
\examples{
\donttest{
a <- nomis_get_metadata("NM_1_1")

print(a)

b <- nomis_get_metadata("NM_1_1", "geography")

tibble::glimpse(b)

# returns all types of geography
c <- nomis_get_metadata("NM_1_1", "geography", "TYPE")

tibble::glimpse(c)

# returns geography types available within Wigan
d <- nomis_get_metadata("NM_1_1", "geography", "1879048226")

tibble::glimpse(d)

e <- nomis_get_metadata("NM_1_1", "item", geography = 1879048226, sex = 5)

print(e)

f <- nomis_get_metadata("NM_1_1", "item", search = "*married*")

tibble::glimpse(f)
}

}
\seealso{
\code{\link[=nomis_data_info]{nomis_data_info()}}

\code{\link[=nomis_get_data]{nomis_get_data()}}

\code{\link[=nomis_overview]{nomis_overview()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset_info.R
\name{nomis_data_info}
\alias{nomis_data_info}
\title{Nomis data structures}
\usage{
nomis_data_info(id, tidy = FALSE)
}
\arguments{
\item{id}{Dataset ID. If empty, returns data on all available datasets.
If the ID of a dataset, returns metadata for that particular dataset.}

\item{tidy}{If \code{TRUE}, converts tibble names to snakecase.}
}
\value{
A tibble with all available datasets and their metadata.
}
\description{
Retrieve metadata on the structure and available variables for all available
data sets or the information available in a specific dataset based on its ID.
}
\examples{
\donttest{

# Get info on all datasets
x <- nomis_data_info()

tibble::glimpse(x)

# Get info on a particular dataset
y <- nomis_data_info("NM_1658_1")

tibble::glimpse(y)
}

}
\seealso{
\code{\link[=nomis_get_data]{nomis_get_data()}}

\code{\link[=nomis_get_metadata]{nomis_get_metadata()}}

\code{\link[=nomis_overview]{nomis_overview()}}

\code{\link[=nomis_codelist]{nomis_codelist()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overview.R
\name{nomis_overview}
\alias{nomis_overview}
\title{Nomis dataset overview}
\usage{
nomis_overview(id, select = NULL)
}
\arguments{
\item{id}{The ID of the particular dataset. Returns no data if not specified.}

\item{select}{A string or character vector of one or more overview parts to
select, excluding all others. \code{select} is not case sensitive. The
options for \code{select} are described below, and are taken from the
\href{https://www.nomisweb.co.uk/api/v01/help}{Nomis API help page}.}
}
\value{
A tibble with two columns, one a character vector with the name of
the metadata category, and the other a list column of values for each
category.
}
\description{
Returns an overview of available metadata for a given dataset.
}
\section{Overview part options}{


\describe{
\item{DatasetInfo}{General dataset information such as name, description,
sub-description, mnemonic, restricted access and status}
\item{Coverage}{Shows the geographic coverage of the main geography
dimension in this dataset (e.g. United Kingdom, England and Wales etc.)}
\item{Keywords}{The keywords allocated to the dataset}
\item{Units}{The units of measure supported by the dataset}
\item{ContentTypes}{The classifications allocated to this dataset}
\item{DateMetadata}{Information about the first release, last update and
next update}
\item{Contact}{Details for the point of contact for this dataset}
\item{Analyses}{Show the available analysis breakdowns of this dataset}
\item{Dimensions}{Individual dimension information (e.g. sex, geography,
date, etc.)}
\item{Dimension-concept}{Allows a specific dimension to be selected (e.g.
dimension-geography would allow information about geography dimension). This
is not used if "Dimensions" is specified too.}
\item{Codes}{Full list of selectable codes, excluding Geography, which as a
list of Types instead. (Requires "Dimensions" to be selected too)}
\item{Codes-concept}{Full list of selectable codes for a specific dimension,
excluding Geography, which as a list of Types instead. This is not used if
"Codes" is specified too (Requires "Dimensions" or equivalent to be
selected too)}
\item{DimensionMetadata}{Any available metadata attached at the dimensional
level (Requires "Dimensions" or equivalent to be selected too)}
\item{Make}{Information about whether user defined codes can be created with
the MAKE parameter when querying data (Requires "Dimensions" or equivalent
to be selected too)}
\item{DatasetMetadata}{Metadata attached at the dataset level}
}
}

\examples{
\donttest{
library(dplyr)

q <- nomis_overview("NM_1650_1")

q \%>\%
  tidyr::unnest(name) \%>\%
  glimpse()

s <- nomis_overview("NM_1650_1", select = c("Units", "Keywords"))

s \%>\%
  tidyr::unnest(name) \%>\%
  glimpse()
}

}
\seealso{
\code{\link[=nomis_data_info]{nomis_data_info()}}

\code{\link[=nomis_get_metadata]{nomis_get_metadata()}}
}
