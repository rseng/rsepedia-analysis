
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nomisr

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/nomisr)](https://cran.r-project.org/package=nomisr)
[![GitHub
tag](https://img.shields.io/github/tag/ropensci/nomisr.svg)](https://github.com/ropensci/nomisr)
[![](https://cranlogs.r-pkg.org/badges/grand-total/nomisr)](https://cran.r-project.org/package=nomisr)
[![R build
status](https://github.com/ropensci/nomisr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/nomisr/actions)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/evanodell/nomisr?branch=master&svg=true)](https://ci.appveyor.com/project/evanodell/nomisr)
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

There are two `nomisr` vignettes. The
[introduction](https://docs.evanodell.com/nomisr/articles/introduction.html)
has details on all available functions and basic demonstrations of their
use. The [Work and Health Indicators with
nomisr](https://docs.evanodell.com/nomisr/articles/introduction-to-work-and-health-nomis-indicators)
vignette shows a worked-through demonstration demonstrating the use of
key indicators, courtesy of [Nina Robery](https://github.com/ninarobery)
from Public Health England.

The example below demostrates a workflow to retrieve the latest data on
Jobseeker’s Allowance with rates and proportions, on a national level,
with all male claimants and workforce.

``` r
 library(nomisr)
 jobseekers_search <- nomis_search(name = "*Jobseeker*")
 
 tibble::glimpse(jobseekers_search)
#> Rows: 17
#> Columns: 14
#> $ agencyid                             <chr> "NOMIS", "NOMIS", "NOMIS", "NOMI…
#> $ id                                   <chr> "NM_1_1", "NM_4_1", "NM_8_1", "N…
#> $ uri                                  <chr> "Nm-1d1", "Nm-4d1", "Nm-8d1", "N…
#> $ version                              <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
#> $ annotations.annotation               <list> [<data.frame[10 x 2]>, <data.fr…
#> $ components.attribute                 <list> [<data.frame[7 x 4]>, <data.fra…
#> $ components.dimension                 <list> [<data.frame[5 x 3]>, <data.fra…
#> $ components.primarymeasure.conceptref <chr> "OBS_VALUE", "OBS_VALUE", "OBS_V…
#> $ components.timedimension.codelist    <chr> "CL_1_1_TIME", "CL_4_1_TIME", "C…
#> $ components.timedimension.conceptref  <chr> "TIME", "TIME", "TIME", "TIME", …
#> $ description.value                    <chr> "Records the number of people cl…
#> $ description.lang                     <chr> "en", "en", NA, "en", "en", "en"…
#> $ name.value                           <chr> "Jobseeker's Allowance with rate…
#> $ name.lang                            <chr> "en", "en", "en", "en", "en", "e…

 jobseekers_measures <- nomis_get_metadata("NM_1_1", "measures")
 
 tibble::glimpse(jobseekers_measures)
#> Rows: 4
#> Columns: 3
#> $ id             <chr> "20100", "20201", "20202", "20203"
#> $ label.en       <chr> "claimants", "workforce", "active", "residence"
#> $ description.en <chr> "claimants", "workforce", "active", "residence"
 
 jobseekers_geography <- nomis_get_metadata("NM_1_1", "geography", "TYPE")
 
 tail(jobseekers_geography)
#> # A tibble: 6 x 3
#>   id      label.en                          description.en                      
#>   <chr>   <chr>                             <chr>                               
#> 1 TYPE490 government office regions tec / … government office regions tec / lec…
#> 2 TYPE491 government office regions (forme… government office regions (former i…
#> 3 TYPE492 standard statistical regions      standard statistical regions        
#> 4 TYPE496 pre-1996 local authority distric… pre-1996 local authority districts  
#> 5 TYPE498 pre-1996 counties / scottish reg… pre-1996 counties / scottish regions
#> 6 TYPE499 countries                         countries
 
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
#> 
#> ── Column specification ────────────────────────────────────────────────────────
#> cols(
#>   .default = col_double(),
#>   DATE = col_character(),
#>   DATE_NAME = col_character(),
#>   DATE_CODE = col_character(),
#>   DATE_TYPE = col_character(),
#>   GEOGRAPHY_NAME = col_character(),
#>   GEOGRAPHY_CODE = col_character(),
#>   GEOGRAPHY_TYPE = col_character(),
#>   SEX_NAME = col_character(),
#>   SEX_TYPE = col_character(),
#>   ITEM_NAME = col_character(),
#>   ITEM_TYPE = col_character(),
#>   MEASURES_NAME = col_character(),
#>   OBS_STATUS = col_character(),
#>   OBS_STATUS_NAME = col_character(),
#>   OBS_CONF = col_logical(),
#>   OBS_CONF_NAME = col_character(),
#>   URN = col_character()
#> )
#> ℹ Use `spec()` for the full column specifications.
 
 tibble::glimpse(z)
#> Rows: 70
#> Columns: 34
#> $ DATE                <chr> "2020-11", "2020-11", "2020-11", "2020-11", "2020…
#> $ DATE_NAME           <chr> "November 2020", "November 2020", "November 2020"…
#> $ DATE_CODE           <chr> "2020-11", "2020-11", "2020-11", "2020-11", "2020…
#> $ DATE_TYPE           <chr> "date", "date", "date", "date", "date", "date", "…
#> $ DATE_TYPECODE       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ DATE_SORTORDER      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ GEOGRAPHY           <dbl> 2092957697, 2092957697, 2092957697, 2092957697, 2…
#> $ GEOGRAPHY_NAME      <chr> "United Kingdom", "United Kingdom", "United Kingd…
#> $ GEOGRAPHY_CODE      <chr> "K02000001", "K02000001", "K02000001", "K02000001…
#> $ GEOGRAPHY_TYPE      <chr> "countries", "countries", "countries", "countries…
#> $ GEOGRAPHY_TYPECODE  <dbl> 499, 499, 499, 499, 499, 499, 499, 499, 499, 499,…
#> $ GEOGRAPHY_SORTORDER <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1…
#> $ SEX                 <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
#> $ SEX_NAME            <chr> "Male", "Male", "Male", "Male", "Male", "Male", "…
#> $ SEX_CODE            <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
#> $ SEX_TYPE            <chr> "sex", "sex", "sex", "sex", "sex", "sex", "sex", …
#> $ SEX_TYPECODE        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ SEX_SORTORDER       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ ITEM                <dbl> 1, 1, 2, 2, 3, 3, 4, 4, 9, 9, 1, 1, 2, 2, 3, 3, 4…
#> $ ITEM_NAME           <chr> "Total claimants", "Total claimants", "Students o…
#> $ ITEM_CODE           <dbl> 1, 1, 2, 2, 3, 3, 4, 4, 9, 9, 1, 1, 2, 2, 3, 3, 4…
#> $ ITEM_TYPE           <chr> "item", "item", "item", "item", "item", "item", "…
#> $ ITEM_TYPECODE       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ ITEM_SORTORDER      <dbl> 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 1, 1, 2, 2, 3…
#> $ MEASURES            <dbl> 20100, 20201, 20100, 20201, 20100, 20201, 20100, …
#> $ MEASURES_NAME       <chr> "Persons claiming JSA", "Workplace-based estimate…
#> $ OBS_VALUE           <dbl> 179827.0, 0.9, NA, NA, NA, NA, NA, NA, NA, NA, 17…
#> $ OBS_STATUS          <chr> "A", "A", "Q", "Q", "Q", "Q", "Q", "Q", "Q", "Q",…
#> $ OBS_STATUS_NAME     <chr> "Normal Value", "Normal Value", "These figures ar…
#> $ OBS_CONF            <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, …
#> $ OBS_CONF_NAME       <chr> "Free (free for publication)", "Free (free for pu…
#> $ URN                 <chr> "Nm-1d1d32331e0d2092957697d5d1d20100", "Nm-1d1d32…
#> $ RECORD_OFFSET       <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,…
#> $ RECORD_COUNT        <dbl> 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 7…
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

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)




# nomisr 0.4.4.9000

* Better error message when API returns empty data in some circumstances.

* Suppressed printing of column types when reading CSV files (#25, thanks jackobailey)


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
CSV files.


## Test environments
* local macOS install, R 4.1.1
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
