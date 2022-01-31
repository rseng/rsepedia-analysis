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

# rdhs <img src="tools/logo.png" align="right" style="padding-left:10px;background-color:white;" />

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis-CI Build
Status](https://travis-ci.org/ropensci/rdhs.png?branch=master)](https://travis-ci.org/ropensci/rdhs)
[![codecov.io](https://codecov.io/github/ropensci/rdhs/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rdhs?branch=master)
[![Documentation via
pkgdown](https://github.com/ropensci/rdhs/raw/master/tools/pkgdownshield.png)](https://docs.ropensci.org/rdhs/)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/rdhs)](https://cran.r-project.org/package=rdhs)
[![Downloads from Rstudio
mirror](https://cranlogs.r-pkg.org/badges/grand-total/rdhs)](https://www.r-pkg.org:443/pkg/rdhs)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/rdhs)](https://cran.r-project.org/package=rdhs)
[![rOpenSci](https://badges.ropensci.org/238_status.svg)](https://github.com/ropensci/software-review/issues/238)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2423635.svg)](https://doi.org/10.5281/zenodo.2423635)

## Motivation

The Demographic and Health Surveys (DHS) Program has collected
population survey data from over 90 countries for over 30 years. In many
countries, DHS provide the key data that mark progress towards targets
such as the Sustainable Development Goals (SDGs) and inform health
policy. Though standard health indicators are routinely published in
survey final reports, much of the value of DHS is derived from the
ability to download and analyse standardized microdata datasets for
subgroup analysis, pooled multi-country analysis, and extended research
studies. The suite of tools within `rdhs` improves the accessibility of
these datasets for statistical analysis with R, with aim to support
reproducible global health research and simplify common analytical
pipelines.

> For questions regarding how to analyse DHS survey data, please read
> the DHS website’s data section first. If you have any questions after
> this then please create an
> [issue](https://github.com/ropensci/rdhs/issues) with your question.
> It is really likely that your question will help other people and so
> posting them publically as an issue may help others with similar
> questions.

-----

`rdhs` is a package for management and analysis of [Demographic and
Health Survey (DHS)](https://www.dhsprogram.com) data. This includes
functionality to:

1.  Access standard indicator data (i.e. [DHS
    STATcompiler](https://www.statcompiler.com/)) in R via the [DHS
    API](https://api.dhsprogram.com/).
2.  Identify surveys and datasets relevant to a particular analysis.
3.  Download survey datasets from the [DHS
    website](https://dhsprogram.com/data/available-datasets.cfm).
4.  Load datasets and associated metadata into R.
5.  Extract variables and combining datasets for pooled multi-survey
    analyses.

## Installation

You can install the latest version from
[`CRAN`](https://cran.r-project.org/package=rdhs) using:

``` r
install.packages("rdhs")
```

You can also install the development version of `rdhs` with the latest
patches from github with:

``` r
#install.packages("devtools")
devtools::install_github("ropensci/rdhs")
```

``` r
# Load the package
library(rdhs)
```

## Getting started

To be able to **download survey datasets from the DHS website**, you
will need to **set up an account with the DHS website**, which will
enable you to request access to the datasets. Instructions on how to do
this can be found
[here](https://dhsprogram.com/data/Access-Instructions.cfm). The email,
password, and project name that were used to create the account will
then need to be provided to `rdhs` when attempting to download datasets.

-----

  - Request dataset access from the DHS website
    [here](https://dhsprogram.com/data/Access-Instructions.cfm).

  - Full functionality is described in the tutorial
    [here](https://docs.ropensci.org/rdhs/articles/introduction.html).

  - An example workflow using `rdhs` to calculate trends in anemia
    prevalence is available
    [here](https://docs.ropensci.org/rdhs/articles/anemia.html).

## Basic Functionality

### Query the [DHS API](https://api.dhsprogram.com/).

Obtain survey estimates for Malaria prevalence among children from the
Democratic Republic of Congo and Tanzania in the last 5 years (since
2013) that included rapid diagnostic tests
(RDTs).

``` r
dhs_indicators(indicatorIds = "ML_PMAL_C_RDT", returnFields=c("IndicatorId", "ShortName"))
#>                             ShortName   IndicatorId
#> 1 Malaria prevalence according to RDT ML_PMAL_C_RDT

dhs_data(countryIds = c("CD","TZ"), indicatorIds = "ML_PMAL_C_RDT", surveyYearStart = 2013,
       returnFields=c("Indicator", "SurveyId", "Value", "SurveyYearLabel", "CountryName"))
#>                             Indicator  SurveyId SurveyYearLabel Value
#> 1 Malaria prevalence according to RDT CD2013DHS         2013-14  30.8
#> 2 Malaria prevalence according to RDT TZ2015DHS         2015-16  14.4
#> 3 Malaria prevalence according to RDT TZ2017MIS            2017   7.3
#>                 CountryName
#> 1 Congo Democratic Republic
#> 2                  Tanzania
#> 3                  Tanzania
```

### Identify survey datasets

Now, obtain survey microdatasets to analyze these same indicators. Query
the *surveyCharacteristics* endpoint to identify the survey
characteristic ID for malaria RDT testing.

``` r
## call with no arguments to return all characterstics
sc <- dhs_survey_characteristics()
sc[grepl("Malaria", sc$SurveyCharacteristicName), ]
#>    SurveyCharacteristicID SurveyCharacteristicName
#> 58                     96            Malaria - DBS
#> 59                     90     Malaria - Microscopy
#> 60                     89            Malaria - RDT
#> 61                     57 Malaria bednet inventory
```

Use `dhs_surveys()` identify surveys for the countries and years of
interest.

``` r
## what are the countryIds - we can find that using this API request
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))

## find all the surveys that match the search criteria
survs <- dhs_surveys(surveyCharacteristicIds = 89, countryIds = c("CD","TZ"), surveyYearStart = 2013)
```

Lastly, identify the datasets required for download. By default, the
recommended option is to download either the spss (.sav), `fileFormat =
"SV"`, or the flat file (.dat), `fileFormat = "FL"` datasets. The flat
is quicker, but there are still one or two very old datasets that don’t
read correctly, whereas the .sav files are slower to read in but so far
no datasets have been found that don’t read in correctly. The household
member recode (`PR`) reports the RDT status for children under
five.

``` r
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "FL", fileType = "PR")
str(datasets)
#> 'data.frame':    3 obs. of  13 variables:
#>  $ FileFormat          : chr  "Flat ASCII data (.dat)" "Flat ASCII data (.dat)" "Flat ASCII data (.dat)"
#>  $ FileSize            : int  6595349 6491292 2171918
#>  $ DatasetType         : chr  "Survey Datasets" "Survey Datasets" "Survey Datasets"
#>  $ SurveyNum           : int  421 485 529
#>  $ SurveyId            : chr  "CD2013DHS" "TZ2015DHS" "TZ2017MIS"
#>  $ FileType            : chr  "Household Member Recode" "Household Member Recode" "Household Member Recode"
#>  $ FileDateLastModified: chr  "September, 19 2016 09:58:23" "September, 28 2019 17:58:28" "June, 11 2019 15:38:22"
#>  $ SurveyType          : chr  "DHS" "DHS" "MIS"
#>  $ SurveyYearLabel     : chr  "2013-14" "2015-16" "2017"
#>  $ SurveyYear          : chr  "2013" "2015" "2017"
#>  $ DHS_CountryCode     : chr  "CD" "TZ" "TZ"
#>  $ FileName            : chr  "CDPR61FL.ZIP" "TZPR7BFL.ZIP" "TZPR7IFL.ZIP"
#>  $ CountryName         : chr  "Congo Democratic Republic" "Tanzania" "Tanzania"
```

### Download datasets

We can now go ahead and download our datasets. To be able to download
survey datasets from the DHS website, you will need to set up an account
with them to enable you to request access to the datasets. Instructions
on how to do this can be found
[here](https://dhsprogram.com/data/Access-Instructions.cfm). The email,
password, and project name that were used to create the account will
then need to be provided to `rdhs` when attempting to download datasets.

Once we have created an account, we need to set up our credentials using
the function `set_rdhs_config()`. This will require providing as
arguments your `email` and `project` for which you want to download
datasets from. You will then be prompted for your password.

You can also specify a directory for datasets and API calls to be cached
to using `cache_path`. In order to comply with CRAN, this function will
also ask you for your permission to write to files outside your
temporary directory, and you must type out the filename for the
`config_path` - “rdhs.json”. (See [introduction
vignette](https://docs.ropensci.org/rdhs/articles/introduction.html) for
specific format for config, or `?set_rdhs_config`).

``` r
## login
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "rdhs R package development",
                config_path = "rdhs.json",
                global = FALSE)
#> Writing your configuration to:
#>    -> rdhs.json
```

The path to your config is saved between sessions so you only have to
set this once. With your credentials set, all API requests will be
cached within the `cache_path` directory provided so that these can be
returned when working remotely or with a poor internet connection.

``` r
# the first time this will take a few seconds 
microbenchmark::microbenchmark(dhs_datasets(surveyYearStart = 1986),times = 1)
#> Unit: milliseconds
#>                                  expr     min      lq    mean  median      uq
#>  dhs_datasets(surveyYearStart = 1986) 46.3744 46.3744 46.3744 46.3744 46.3744
#>      max neval
#>  46.3744     1

# after caching, results will be available instantly
microbenchmark::microbenchmark(dhs_datasets(surveyYearStart = 1986),times = 1)
#> Unit: milliseconds
#>                                  expr      min       lq     mean   median
#>  dhs_datasets(surveyYearStart = 1986) 1.410894 1.410894 1.410894 1.410894
#>        uq      max neval
#>  1.410894 1.410894     1
```

Now download datasets by providing a list of desired dataset filenames.

``` r
# download datasets
downloads <- get_datasets(datasets$FileName)

str(downloads)
#> List of 3
#>  $ CDPR61FL: chr "/home/oj/.cache/rdhs/datasets/CDPR61FL.rds"
#>  $ TZPR7BFL: chr "/home/oj/.cache/rdhs/datasets/TZPR7BFL.rds"
#>  $ TZPR7IFL: chr "/home/oj/.cache/rdhs/datasets/TZPR7IFL.rds"
#>  - attr(*, "reformat")= logi FALSE
```

### Load datasets into R

The `get_datasets()` function returns a vector with a file path to the
saved location of the downloaded datasets. These are read using
`readRDS()`:

``` r
# read in first dataset
cdpr <- readRDS(downloads$CDPR61FL)
```

Value labels are stored as attributes to each of the columns of the data
frame using the `labelled` class (see `haven::labelled` or our
introduction vignette for more details). Variable labels are stored in
the `label` attribute.

### Extract variables and pool datasets

The client also caches all variable labels to quickly query variables in
each survey *without* loading the datasets.

``` r
# rapid diagnostic test search
vars <- search_variable_labels(datasets$FileName, search_terms = "malaria rapid test")
```

Then extract these variables from the datasets. Optionally, geographic
data may be added.

``` r
# and now extract the data
extract <- extract_dhs(vars, add_geo = FALSE)
#> Starting Survey 1 out of 3 surveys:CDPR61FL
#> Starting Survey 2 out of 3 surveys:TZPR7BFL
#> Starting Survey 3 out of 3 surveys:TZPR7IFL
```

The returned object is a list of extracted datasets.

Dataset extracts can alternate be specified by providing a vector of
surveys and vector of variable names:

``` r
# and grab the questions from this now utilising the survey variables
vars <- search_variables(datasets$FileName, variables = c("hv024","hml35"))

# and now extract the data
extract <- extract_dhs(vars, add_geo = FALSE)
#> Starting Survey 1 out of 3 surveys:CDPR61FL
#> Starting Survey 2 out of 3 surveys:TZPR7BFL
#> Starting Survey 3 out of 3 surveys:TZPR7IFL
```

Finally, the two datasets are pooled using the function
`rbind_labelled()`. This function works specifically with our lists of
labelled `data.frame`s. Labels are specified for each variable: for
`hv024` all labels are retained (concatenate) but for `hml35` labels
across both datasets to be “Neg” and “Pos”.

``` r
# now let's try our second extraction
extract <- rbind_labelled(extract,
                          labels = list("hv024" = "concatenate",
                                        "hml35" = c("Neg"=0, "Pos"=1)))
```

There is also an option to process downloaded datasets with labelled
variables coded as strings, rather than labelled variables. This is
specified by the argument `reformat=TRUE`.

``` r
# identify questions but specifying the reformat argument
questions <- search_variables(datasets$FileName, variables = c("hv024", "hml35"),
                                     reformat=TRUE)

# and now extract the data
extract <- extract_dhs(questions, add_geo = FALSE)
#> Starting Survey 1 out of 3 surveys:CDPR61FL
#> Starting Survey 2 out of 3 surveys:TZPR7BFL
#> Starting Survey 3 out of 3 surveys:TZPR7IFL

# group our results
extract <- rbind_labelled(extract)

# our hv024 variable is now just character strings, so you can decide when/how to factor/label it later
str(extract)
#> Classes 'dhs_dataset' and 'data.frame':  208595 obs. of  4 variables:
#>  $ hv024   : chr  "equateur" "equateur" "equateur" "equateur" ...
#>   ..- attr(*, "label")= chr "Province"
#>  $ hml35   : chr  NA NA NA NA ...
#>   ..- attr(*, "label")= chr "Result of malaria rapid test"
#>  $ SurveyId: chr  "CD2013DHS" "CD2013DHS" "CD2013DHS" "CD2013DHS" ...
#>  $ DATASET : chr  "CDPR61FL" "CDPR61FL" "CDPR61FL" "CDPR61FL" ...
```

-----

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
## rdhs 0.7.3

* Reference internal access to `model_datasets` with `rdhs::model_datasets` to avoid errors if `rdhs` namespace is not loaded.

## rdhs 0.7.2

* `available_datasets` patch (#115)

## rdhs 0.7.1

* Remove class `"dhs_dataset"` from downloaded micro data sets. This class is not 
  anywhere and it creates an error for dplyr_v1.0. 
  
  Cached datasets will need to be re-downloaded after updating to clear the 
  dhs_dataset clas.

* Replace `readLines()` with `brio::read_lines()` to make parsers robust to 
  Windows encoding issues (similar to https://stackoverflow.com/questions/18789330/r-on-windows-character-encoding-hell).

* Use `"sf"` as default download method for `download_boundaries(..., method = "sf")`.
  Add arguments `quiet_download` and `quiet_parse = TRUE` to 
  `download_boundaries()`. `quiet_download` (default `FALSE`) controls `download.file()` 
  messages. `quiet_parse` (default `TRUE`) controls messages from `sf::st_read()` when
  `method = "sf"`.

## rdhs 0.7.0

* Add CITATION info.
* New `download_boundaries` for downloading spatial boundaries using (#71)
* New `dhs_gps_data_format` for DHS GPS Information (#74)
* Tibbles can be specified correctly as data.frame format (#89)
* Config creation on Windows 10 fixed (#91)
* Typos and messaging fixed (#78, #84, #87, #92)
* `unzip_special` correctly detects 4Gb files (#43)

## rdhs 0.6.3

* Addresses CRAN fail on windows 
* New `delabel_df` for converting labelled data frames to characters (#54)

## rdhs 0.6.2

* Duplicate labels when parsing flat data files corrected (#79)

## rdhs 0.6.1

* `extraction(add_geo=TRUE)` correction for Kenya 2014 surveys (#67)
* Geospatial covariate data sets now supported correctly (#64)

## rdhs 0.6.0

* New `as_factor.labelled` for backward compatibility with `haven <2.0.0` 
`labelled` classes.

* `model_datasets` now internal and exported dataset (#60).

## rdhs 0.5.2

* New vignettes: `country_codes`

* Documentation typos corrected (#55)

## rdhs 0.5.0

* New `set_rdhs_config` for providing login credentials. This deprecates
`set_dhs_credentials`. 

* New `get_rdhs_config` shows the credentials currently used by `rdhs`

## rdhs 0.4.0

* Permissionn from user now required for file saving.

* API requests can now ignore any cached responses (`force = TRUE`) argument 
(#23):

```R
dat <- dhs_countries(force = TRUE)
```

* `get_datasets(clear_cache = TRUE)` will clear the cached available datasets,
enabling newly requested datasets to be downloaded (@kaisero, #29).

* geojson objects can now be requested from the API so that you can return
geojson objects for mapping purposes (#28) e.g. :

```R
d <- dhs_data(countryIds = "SN", surveyYearStart = 2014, 
              breakdown = "subnational", returnGeometry = TRUE,
              f = "geojson")
              
# convert to spatial object
sp <- geojsonio::as.json(d) %>% geojsonio::geojson_sp

```

## rdhs 0.3.0

* New `dhs_data()`, `dhs_countries()` and other API functions (`dhs_x()`). 

## rdhs 0.2.1

* `authenticate_dhs()` now works with short project names.

## rdhs 0.2.0

* New vignettes: `anemia`

* New `read_dhs_flat()` for reading flat datasets.

## rdhs 0.1.0

* Initial share on Feb, 24th 2018 to colleagues at UNC.
## Test environments
* windows-latest (release) on Github Actions
* macOS-latest (release) on Github Actions
* ubuntu-20.04 (release) on Github Actions
* ubuntu-20.04 (devel) on Github Actions

## R CMD check results
## rdhs 0.7.2

* Duration: 11m 27.9s
* 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
* R CMD check succeeded

## Downstream dependencies

There are 2 non-strong downstream dependencies.
---
title: 'rdhs: an R package to interact with The Demographic and Health Surveys (DHS) Program data sets'
tags:
  - R
  - DHS
  - API
  - survey analysis
authors:
  - name: Oliver J Watson
    orcid: 0000-0003-2374-0741
    affiliation: "1"
  - name: Jeffrey W Eaton
    orcid: 0000-0001-7728-728X
    affiliation: "1"
affiliations:
  - name: 1.	MRC Centre for Global Infectious Disease Analysis, Department of Infectious Disease Epidemiology, Imperial College London
    index: 1
date: 8th May 2018
bibliography: paper.bib
---

# Summary

The Demographic and Health Surveys (DHS) Program has collected and disseminated population survey data from over 90 countries for over 30 years. In many countries, DHS provide the key data that mark progress towards targets such as the Sustainable Development Goals (SDGs) and inform health policy such as detailing trends in child mortality [@Silva2012] and characterising the distribution of malaria control interventions in Africa in order to map the burden of malaria since the year 2000 [@Bhatt2015]. Though standard health indicators are routinely published in survey final reports, much of the value of DHS is derived from the ability to download and analyse standardized microdata datasets for subgroup analysis, pooled multi-country analysis, and extended research studies. 

The analysis of the microdata datasets, however, requires a 'clean' dataset that contains all the desired information. One of the main challenges when interacting with the raw DHS datasets is isolating the required dataset variables across different countries. Since the DHS Program started, there have been 7 'phases' of questionnaires used between 1984 - 2018. The data from each phase then recoded to consistency and comparability across surveys. However, new questions are often included or ammended between different phases of the DHS program, which results in variable names sometimes changing between different phases. As well as this, there are a number of country specific records that are not part of model questionnaires. As such, it can become increasingly difficult to identify which variables to use within your final 'clean' dataset. 

The rdhs package was designed to facilitate the management and processing of DHS survey data. This occurs through both functioning as an API client, allowing access to all data provided within the DHS API, and helping to download the raw datasets from the DHS website and read them into conventional R data structures. In overview, the package provides a suite of tools for the following:

 + 1. Access standard survey indicators through the DHS Program API. This data includes the same summarised health indicators that is available through the DHS web application STATcompiler, as well as additional data endpoints that provide metadata for the conducted surveys and raw datasets.  
 + 2. Identify all survey datasets that include a particular topic or indicator relevant to a particular analysis. 
 + 3. Directly download survey datasets from the DHS website. 
 + 4. Load datasets and data dictionaries into R.
 + 5. Extract variables and pool harmonized datasets for multi-survey analysis. 

The functionality provided represents the output of conversations with numerous research groups globally, and serves to simplify commonly required analytical pipelines. The end result aims to increase the end user accessibility to the raw data and create a tool that supports reproducible global health research. Furthermore, the package is hoped to enable researches in lower middle income countries, which constitute the majority of countries that are surveyed as part of the DHS program, to analyse their data without the need for proprietary software. 


# References
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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
*  Look at the Travis build status before and after making changes.
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

### Writing tests that require downloading datasets

In order to check your contribution passes all tests, you will needt be able to
run the full test suite. Instructions for how to do this can be found 
[here](https://ropensci.github.io/rdhs/articles/testing.html). 

Alternatively, you can use the package's CI with [travis](https://travis-ci.org/ropensci/rdhs),
to see if the test suite has passed. This test suite uses the encrypted 
"rdhs.json" config within the package root to download datasets. The datasets
that are available to download with this config can be viewed within the
"available_test_suite_datasets.csv". 

If your contributions need access to a dataset that is not available, e.g. a new 
dataset that affects the flat file parsers, then open a PR and they will be
requested from the DHS website.

### Code of Conduct

Please note that the rdhs project is released with a
[Contributor Code of Conduct](CONDUCT.md). By contributing to this
project you agree to abide by its terms. See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

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
# Location
* from: github.com/schloerke/leaflet-providers@urlProtocol

* Inspiration taken from https://github.com/leaflet-extras/leaflet-providers/commit/dea786a3219f9cc824b8e96903a17f46ca9a5afc to use the 'old' relative url protocols and to 'upgrade' them at js runtime.



# Notes...

* Copy/paste provider information into `providers.json`
```js
var providers = L.TileLayer.Provider.providers;
JSON.stringify(providers, null, "  ");
```
  * `./data-raw/providerNames.R` was re-ran to update to the latest providers

* Some providers had their protocols turned into '//'.
  * This allows browsers to pick the protocol
  * To stop files from the protocols staying as files, a ducktape patch was applied to `L.TileLayer.prototype.initialize` and `L.TileLayer.WMS.prototype.initialize`
---
output:
  rmarkdown::github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# rdhs <img src="tools/logo.png" align="right" style="padding-left:10px;background-color:white;" />

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/rdhs.png?branch=master)](https://travis-ci.org/ropensci/rdhs)
[![codecov.io](https://codecov.io/github/ropensci/rdhs/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rdhs?branch=master)
[![Documentation via pkgdown](https://github.com/ropensci/rdhs/raw/master/tools/pkgdownshield.png)](https://docs.ropensci.org/rdhs/)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/rdhs)](https://cran.r-project.org/package=rdhs)
[![Downloads from Rstudio mirror](https://cranlogs.r-pkg.org/badges/grand-total/rdhs)](https://www.r-pkg.org:443/pkg/rdhs)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/rdhs)](https://cran.r-project.org/package=rdhs)
[![rOpenSci](https://badges.ropensci.org/238_status.svg)](https://github.com/ropensci/software-review/issues/238)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2423635.svg)](https://doi.org/10.5281/zenodo.2423635)

## Motivation

The Demographic and Health Surveys (DHS) Program has collected population survey data from over 90 countries for over 30 years. In many countries, DHS provide the key data that mark progress towards targets such as the Sustainable Development Goals (SDGs) and inform health policy. Though standard health indicators are routinely published in survey final reports, much of the value of DHS is derived from the ability to download and analyse standardized microdata datasets for subgroup analysis, pooled multi-country analysis, and extended research studies. The suite of tools within `rdhs` improves the accessibility of these datasets for statistical analysis with R, with aim to support reproducible global health research and simplify common analytical pipelines.

 > For questions regarding how to analyse DHS survey data, please read the DHS website's data section first. If you have any questions after this then please create an [issue](https://github.com/ropensci/rdhs/issues) with your question. It is really likely that your question will help other people and so posting them publically as an issue may help others with similar questions.

---

`rdhs` is a package for management and analysis of [Demographic and Health Survey (DHS)](https://www.dhsprogram.com) data. This includes functionality to:

1. Access standard indicator data (i.e. [DHS STATcompiler](https://www.statcompiler.com/)) in R via the [DHS API](https://api.dhsprogram.com/).
1. Identify surveys and datasets relevant to a particular analysis.
1. Download survey datasets from the [DHS website](https://dhsprogram.com/data/available-datasets.cfm).
1. Load datasets and associated metadata into R.
1. Extract variables and combining datasets for pooled multi-survey analyses.

## Installation

You can install the latest version from [`CRAN`](https://cran.r-project.org/package=rdhs) using:

```{r cran-installation, message=FALSE, eval = FALSE}
install.packages("rdhs")
```

You can also install the development version of `rdhs` with the latest patches from github with:

```{r gh_installation, message=FALSE, eval = FALSE}
#install.packages("devtools")
devtools::install_github("ropensci/rdhs")
```

```{r}
# Load the package
library(rdhs)
```

## Getting started

To be able to **download survey datasets from the DHS website**, you will need to **set up an account with the DHS website**, which will enable you to request access to the datasets. Instructions on how to do this can be found [here](https://dhsprogram.com/data/Access-Instructions.cfm). The email, password, and project name that were used to create the account will then need to be provided to `rdhs` when attempting to download datasets. 

---

* Request dataset access from the DHS website [here](https://dhsprogram.com/data/Access-Instructions.cfm).

* Full functionality is described in the tutorial [here](https://docs.ropensci.org/rdhs/articles/introduction.html).

* An example workflow using `rdhs` to calculate trends in anemia prevalence is available [here](https://docs.ropensci.org/rdhs/articles/anemia.html).

## Basic Functionality

### Query the [DHS API](https://api.dhsprogram.com/).

Obtain survey estimates for Malaria prevalence among children from the Democratic Republic of Congo and Tanzania in the last 5 years (since 2013) that included rapid diagnostic tests (RDTs).

```{r api, message = FALSE}
dhs_indicators(indicatorIds = "ML_PMAL_C_RDT", returnFields=c("IndicatorId", "ShortName"))

dhs_data(countryIds = c("CD","TZ"), indicatorIds = "ML_PMAL_C_RDT", surveyYearStart = 2013,
       returnFields=c("Indicator", "SurveyId", "Value", "SurveyYearLabel", "CountryName"))
```

### Identify survey datasets

Now, obtain survey microdatasets to analyze these same indicators. Query the *surveyCharacteristics* endpoint to identify the survey characteristic ID for malaria RDT testing.

```{r sc}
## call with no arguments to return all characterstics
sc <- dhs_survey_characteristics()
sc[grepl("Malaria", sc$SurveyCharacteristicName), ]
```

Use `dhs_surveys()` identify surveys for the countries and years of interest.

```{r surv}
## what are the countryIds - we can find that using this API request
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))

## find all the surveys that match the search criteria
survs <- dhs_surveys(surveyCharacteristicIds = 89, countryIds = c("CD","TZ"), surveyYearStart = 2013)
```

Lastly, identify the datasets required for download. By default, the recommended option is to download either the spss (.sav), `fileFormat = "SV"`, or the flat file (.dat), `fileFormat = "FL"` datasets. The flat is quicker, but there are still one or two very old datasets that don't read correctly, whereas the .sav files are slower to read in but so far no datasets have been found that don't read in correctly. The household member recode (`PR`) reports the RDT status for children under five.

```{r }
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "FL", fileType = "PR")
str(datasets)
```

### Download datasets

We can now go ahead and download our datasets. To be able to download survey datasets from the DHS website, you will need to set up an account with them to enable you to request access to the datasets. Instructions on how to do this can be found [here](https://dhsprogram.com/data/Access-Instructions.cfm). The email, password, and project name that were used to create the account will then need to be provided to `rdhs` when attempting to download datasets. 

Once we have created an account, we need to set up our credentials using the function `set_rdhs_config()`. This will require providing as arguments your `email` and `project` for which you want to download datasets from. You will then be prompted for your password.

You can also specify a directory for datasets and API calls to be cached to using `cache_path`. In order to comply with CRAN, this function will also ask you for your permission to write to files outside your temporary directory, and you must type out the filename for the `config_path` - "rdhs.json". (See [introduction vignette](https://docs.ropensci.org/rdhs/articles/introduction.html) for specific format for config, or `?set_rdhs_config`). 

```{r client , R.options = list("rappdir_permission" = TRUE)}
## login
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "rdhs R package development",
                config_path = "rdhs.json",
                global = FALSE)
```

The path to your config is saved between sessions so you only have to set this once. With your credentials set, all API requests will be cached within the `cache_path` directory provided so that these can be returned when working remotely or with a poor internet connection.

```{r client_api_cache}
# the first time this will take a few seconds 
microbenchmark::microbenchmark(dhs_datasets(surveyYearStart = 1986),times = 1)

# after caching, results will be available instantly
microbenchmark::microbenchmark(dhs_datasets(surveyYearStart = 1986),times = 1)
```

Now download datasets by providing a list of desired dataset filenames.

```{r download, message=FALSE}
# download datasets
downloads <- get_datasets(datasets$FileName)

str(downloads)
```

### Load datasets into R

The `get_datasets()` function returns a vector with a file path to the saved location of the downloaded datasets. These are read using `readRDS()`:

```{r read a dataset}
# read in first dataset
cdpr <- readRDS(downloads$CDPR61FL)
```

Value labels are stored as attributes to each of the columns of the data frame using the `labelled` class (see `haven::labelled` or our introduction vignette for more details). Variable labels are stored in the `label` attribute.

### Extract variables and pool datasets

The client also caches all variable labels to quickly query variables in each survey *without* loading the datasets.

```{r questions}
# rapid diagnostic test search
vars <- search_variable_labels(datasets$FileName, search_terms = "malaria rapid test")
```

Then extract these variables from the datasets. Optionally, geographic data may be added.

```{r extract_questions}
# and now extract the data
extract <- extract_dhs(vars, add_geo = FALSE)
```

The returned object is a list of extracted datasets.

Dataset extracts can alternate be specified by providing a vector of surveys and vector of variable names:

```{r extract_variables}
# and grab the questions from this now utilising the survey variables
vars <- search_variables(datasets$FileName, variables = c("hv024","hml35"))

# and now extract the data
extract <- extract_dhs(vars, add_geo = FALSE)
```

Finally, the two datasets are pooled using the function `rbind_labelled()`. This function works specifically with our lists of labelled `data.frame`s. Labels are specified for each variable: for `hv024` all labels are retained (concatenate) but for `hml35` labels across both datasets to be "Neg" and "Pos".

```{r rbind_labelled}
# now let's try our second extraction
extract <- rbind_labelled(extract,
                          labels = list("hv024" = "concatenate",
                                        "hml35" = c("Neg"=0, "Pos"=1)))
```


There is also an option to process downloaded datasets with labelled variables coded as strings, rather than labelled variables. This is specified by the argument `reformat=TRUE`.

```{r reformat}
# identify questions but specifying the reformat argument
questions <- search_variables(datasets$FileName, variables = c("hv024", "hml35"),
                                     reformat=TRUE)

# and now extract the data
extract <- extract_dhs(questions, add_geo = FALSE)

# group our results
extract <- rbind_labelled(extract)

# our hv024 variable is now just character strings, so you can decide when/how to factor/label it later
str(extract)
```

---

 [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Interacting with geojson API returns"
author: "OJ Watson"
date: "2018-09-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Interacting with geojson API results}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

The API can return geojson objects, which can be very useful for quickly creating maps of the 
DHS indicator data. These are large objects though, and can cause the API to return slowly if
too much data is returned. 

In this demonstration we will be using 2 other packages, `leaflet` and `geojson`.

---


```r
# load our package
library(rdhs)

# install other packages
# install.packages("geojson")
# install.packages("leaflet")
```

```r
# make request
d <- dhs_data(countryIds = "SN",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2014,
              breakdown = "subnational",
              returnGeometry = TRUE,
              f = "geojson")

# convert to sp
m <- geojsonio::as.json(d)
nc <- geojsonio::geojson_sp(m) 

# plot using leaflet
pal <- leaflet::colorNumeric("viridis", NULL)

leaflet::leaflet(nc[nc$IndicatorId=="FE_FRTR_W_A15",]) %>%
  leaflet::addTiles() %>%
  leaflet::addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal(log10(Value)),
              label = ~paste0(CharacteristicLabel, ": ", formatC(Value, big.mark = ","))) %>%
  leaflet::addLegend(pal = pal, values = ~log10(Value), opacity = 1.0,
            labFormat = leaflet::labelFormat(transform = function(x) round(10^x)),title = ~Indicator[1])
```

<!--html_preserve--><div id="htmlwidget-c2831947ec24dd72452f" style="width:576px;height:576px;" class="leaflet html-widget"></div>
<script type="application/json" data-for="htmlwidget-c2831947ec24dd72452f">{"x":{"options":{"crs":{"crsClass":"L.CRS.EPSG3857","code":null,"proj4def":null,"projectedBounds":null,"options":{}}},"calls":[{"method":"addTiles","args":["//{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",null,null,{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false,"attribution":"&copy; <a href=\"http://openstreetmap.org\">OpenStreetMap<\/a> contributors, <a href=\"http://creativecommons.org/licenses/by-sa/2.0/\">CC-BY-SA<\/a>"}]},{"method":"addPolygons","args":[[[[{"lng":[-16.5334,-16.5351,-16.5355,-16.5346,-16.5321,-16.5313,-16.5305,-16.5297,-16.5296,-16.5288,-16.5288,-16.5255,-16.5255,-16.5246,-16.5247,-16.523,-16.523,-16.5221,-16.5213,-16.5204,-16.52,-16.518,-16.518,-16.518,-16.5188,-16.5188,-16.5196,-16.5196,-16.5188,-16.5197,-16.5197,-16.523,-16.523,-16.5238,-16.5238,-16.5246,-16.5246,-16.5263,-16.5271,-16.528,-16.528,-16.5288,-16.5296,-16.5299,-16.5305,-16.5313,-16.5334],"lat":[15.7959,15.7959,15.7962,15.8038,15.8096,15.8105,15.817,15.8179,15.8205,15.8213,15.8246,15.8329,15.8371,15.8379,15.8405,15.8429,15.8463,15.8471,15.8587,15.8596,15.8617,15.8604,15.858,15.8504,15.8496,15.8479,15.8471,15.8404,15.8388,15.8379,15.8362,15.8288,15.8254,15.8246,15.8229,15.8221,15.8188,15.8171,15.8129,15.8121,15.8087,15.8079,15.8037,15.8035,15.8029,15.7987,15.7959]},{"lng":[-16.5163,-16.5155,-16.5155,-16.5168,-16.5175,-16.5213,-16.5213,-16.5192,-16.5163],"lat":[15.8662,15.8679,15.8729,15.875,15.875,15.8712,15.8688,15.8659,15.8662]},{"lng":[-16.4872,-16.4871,-16.4879,-16.4883,-16.4896,-16.4896,-16.4897,-16.4884,-16.4872],"lat":[15.9162,15.9213,15.9232,15.9242,15.9237,15.9219,15.9171,15.9158,15.9162]},{"lng":[-16.4959,-16.4946,-16.4946,-16.4946,-16.4946,-16.4955,-16.4946,-16.4955,-16.4955,-16.4955,-16.4955,-16.4963,-16.4963,-16.4946,-16.4951,-16.4967,-16.4975,-16.4996,-16.4996,-16.4971,-16.4988,-16.498,-16.4988,-16.498,-16.4997,-16.498,-16.498,-16.4988,-16.4988,-16.4975,-16.4959],"lat":[15.9359,15.9371,15.9401,15.9437,15.9454,15.9463,15.9488,15.9496,15.9538,15.9552,15.9596,15.9604,15.9638,15.9654,15.9658,15.9658,15.965,15.9629,15.9604,15.9579,15.9554,15.9538,15.9496,15.9479,15.9454,15.9438,15.9396,15.9388,15.9371,15.9359,15.9359]},{"lng":[-16.4988,-16.4988,-16.4988,-16.4988,-16.4963,-16.4964,-16.4976,-16.5001,-16.5021,-16.5022,-16.503,-16.503,-16.503,-16.5021,-16.5021,-16.5009,-16.4988],"lat":[15.9896,15.9911,15.9944,15.9963,15.9996,16.0062,16.0083,16.0083,16.0062,16.002,16.0012,15.9996,15.9979,15.9971,15.9904,15.9884,15.9896]},{"lng":[-16.503,-16.503,-16.503,-16.5022,-16.5021,-16.5005,-16.5013,-16.5035,-16.5064,-16.5064,-16.5064,-16.5063,-16.5043,-16.503],"lat":[16.0204,16.0236,16.0254,16.0263,16.0304,16.0346,16.0387,16.0392,16.0345,16.0233,16.0216,16.0196,16.0192,16.0204]},{"lng":[-16.4772,-16.4772,-16.4755,-16.4755,-16.4738,-16.4739,-16.473,-16.473,-16.4743,-16.4768,-16.4785,-16.4797,-16.4806,-16.483,-16.483,-16.4839,-16.4839,-16.4851,-16.488,-16.4889,-16.4876,-16.4868,-16.4793,-16.4772],"lat":[16.0421,16.0454,16.0471,16.0479,16.0513,16.0571,16.0579,16.0596,16.0609,16.0608,16.06,16.0587,16.0554,16.0529,16.0521,16.0512,16.0479,16.0475,16.0438,16.0388,16.0383,16.0392,16.0392,16.0421]},{"lng":[-16.4968,-16.4975,-16.498,-16.4982,-16.4986,-16.4989,-16.4992,-16.5006,-16.5034,-16.5091,-16.509,-16.509,-16.5092,-16.5098,-16.5107,-16.5105,-16.5113,-16.5113,-16.513,-16.513,-16.5138,-16.5138,-16.5138,-16.5147,-16.5146,-16.5155,-16.5155,-16.5163,-16.5163,-16.5171,-16.5171,-16.518,-16.518,-16.518,-16.5188,-16.5188,-16.5196,-16.5196,-16.5205,-16.5213,-16.52,-16.5188,-16.5171,-16.5163,-16.5163,-16.5155,-16.5155,-16.5138,-16.5138,-16.5146,-16.5146,-16.5155,-16.5146,-16.5138,-16.5138,-16.513,-16.5122,-16.5096,-16.5096,-16.5088,-16.5088,-16.5097,-16.5096,-16.5105,-16.5096,-16.5096,-16.5088,-16.5088,-16.5088,-16.508,-16.5096,-16.5084,-16.5071,-16.5088,-16.5088,-16.5067,-16.5055,-16.5055,-16.508,-16.508,-16.5072,-16.5072,-16.508,-16.508,-16.508,-16.508,-16.508,-16.5072,-16.5072,-16.5056,-16.5055,-16.5064,-16.5064,-16.5047,-16.5047,-16.5047,-16.5047,-16.5043,-16.5031,-16.503,-16.5022,-16.5021,-16.5014,-16.5014,-16.4993,-16.4972,-16.4969,-16.4968,-16.4955,-16.4955,-16.4943,-16.4914,-16.4903,-16.4899,-16.4892,-16.4888,-16.4897,-16.4901,-16.4901,-16.4901,-16.49,-16.4898,-16.4891,-16.489,-16.4891,-16.4891,-16.4891,-16.4893,-16.4896,-16.4896,-16.4895,-16.4895,-16.4896,-16.4897,-16.4899,-16.4901,-16.4902,-16.4904,-16.4905,-16.4908,-16.491,-16.4912,-16.4914,-16.4917,-16.4917,-16.4917,-16.4917,-16.4919,-16.4921,-16.4923,-16.4928,-16.4931,-16.4933,-16.4937,-16.4938,-16.4938,-16.4938,-16.4937,-16.4938,-16.4938,-16.494,-16.4941,-16.4946,-16.4949,-16.4958,-16.4961,-16.4964,-16.4967,-16.4968],"lat":[16.0732,16.0736,16.0737,16.0736,16.0732,16.0727,16.0717,16.0696,16.0664,16.0598,16.0507,16.0415,16.0354,16.0211,16.0122,16.0121,15.9996,15.9988,15.9946,15.9862,15.9854,15.9567,15.9546,15.9538,15.9488,15.948,15.9321,15.9312,15.9237,15.9229,15.9187,15.9179,15.902,15.9004,15.8996,15.8954,15.8946,15.8904,15.8896,15.8779,15.8775,15.8788,15.8879,15.8888,15.8912,15.892,15.8938,15.8971,15.9038,15.9046,15.9088,15.9096,15.9171,15.9179,15.9221,15.9229,15.9287,15.9329,15.9354,15.9362,15.9413,15.9421,15.9463,15.9471,15.9479,15.9546,15.9554,15.9596,15.9729,15.9737,15.9813,15.9816,15.9838,15.9854,15.9871,15.9875,15.9887,15.993,15.9954,16.0013,16.0021,16.0079,16.0087,16.0137,16.0204,16.0223,16.0354,16.0363,16.0379,16.0404,16.0445,16.0454,16.0471,16.0504,16.0527,16.0596,16.0646,16.065,16.0637,16.0595,16.0587,16.0538,16.0529,16.0512,16.0483,16.0496,16.0574,16.0592,16.0571,16.0495,16.0492,16.052,16.0543,16.0555,16.0567,16.0572,16.0577,16.058,16.0581,16.0593,16.0596,16.0599,16.0605,16.0607,16.0613,16.0616,16.0619,16.0622,16.0632,16.0635,16.0639,16.0643,16.0646,16.0648,16.0662,16.0664,16.067,16.067,16.0671,16.0677,16.068,16.068,16.0683,16.0684,16.0686,16.0689,16.0693,16.0693,16.0692,16.0694,16.0697,16.07,16.0703,16.0706,16.0707,16.0708,16.0712,16.0715,16.0718,16.0721,16.0722,16.0723,16.0724,16.0725,16.0728,16.073,16.073,16.073,16.0732]},{"lng":[-14.9459,-14.9503,-14.9528,-14.9545,-14.9556,-14.9564,-14.9573,-14.9575,-14.9584,-14.9609,-14.9631,-14.9689,-14.9739,-14.9803,-14.9861,-14.9903,-14.9945,-14.9981,-14.9992,-15,-15.0025,-15.0031,-15.005,-15.007,-15.0089,-15.0117,-15.0145,-15.017,-15.0206,-15.0248,-15.0306,-15.0353,-15.0411,-15.0475,-15.0534,-15.057,-15.0606,-15.0639,-15.0661,-15.0686,-15.0709,-15.0734,-15.077,-15.0809,-15.0845,-15.0903,-15.097,-15.1017,-15.105,-15.1092,-15.1128,-15.1153,-15.1159,-15.1142,-15.1125,-15.1095,-15.105,-15.1,-15.095,-15.0898,-15.0873,-15.0864,-15.0875,-15.0903,-15.0936,-15.097,-15.1006,-15.1039,-15.1073,-15.1125,-15.1189,-15.1253,-15.1317,-15.1375,-15.1434,-15.1506,-15.157,-15.1636,-15.1684,-15.1734,-15.1798,-15.1848,-15.1903,-15.1959,-15.2017,-15.2075,-15.2123,-15.2173,-15.222,-15.2273,-15.232,-15.2384,-15.245,-15.2506,-15.2559,-15.2617,-15.2667,-15.272,-15.2778,-15.2828,-15.2878,-15.2936,-15.3003,-15.3067,-15.3095,-15.3159,-15.3217,-15.3267,-15.3314,-15.337,-15.342,-15.3478,-15.3542,-15.36,-15.3642,-15.3692,-15.3748,-15.3798,-15.3848,-15.3889,-15.3945,-15.3995,-15.4059,-15.4125,-15.417,-15.4206,-15.4242,-15.4281,-15.4303,-15.4328,-15.4345,-15.4375,-15.4406,-15.445,-15.45,-15.4564,-15.462,-15.4681,-15.4728,-15.4761,-15.4778,-15.4825,-15.4884,-15.4942,-15.4989,-15.4992,-15.5039,-15.5073,-15.5092,-15.5098,-15.5103,-15.5114,-15.5125,-15.5142,-15.5178,-15.5223,-15.5264,-15.532,-15.5367,-15.5425,-15.5489,-15.5561,-15.5628,-15.5689,-15.5756,-15.5814,-15.587,-15.5923,-15.5973,-15.6025,-15.6075,-15.6142,-15.6206,-15.6239,-15.6253,-15.6264,-15.6275,-15.6289,-15.6289,-15.6295,-15.6336,-15.6378,-15.6428,-15.6484,-15.6534,-15.6592,-15.6639,-15.6695,-15.6759,-15.6825,-15.6898,-15.6961,-15.7017,-15.7078,-15.7136,-15.7186,-15.7245,-15.7286,-15.7345,-15.7411,-15.7475,-15.7534,-15.7603,-15.7664,-15.7723,-15.7781,-15.7836,-15.7886,-15.7945,-15.7998,-15.8006,-15.8019,-15.8022,-15.8023,-15.8048,-15.8092,-15.8156,-15.8223,-15.8255,-15.8258,-15.827,-15.832,-15.837,-15.842,-15.8489,-15.855,-15.8592,-15.8614,-15.8616,-15.8617,-15.8636,-15.8675,-15.872,-15.8764,-15.8814,-15.8873,-15.8936,-15.8995,-15.9067,-15.9123,-15.9175,-15.9231,-15.9278,-15.9328,-15.9378,-15.9382,-15.9383,-15.9425,-15.9467,-15.9509,-15.9559,-15.96,-15.9673,-15.9731,-15.9789,-15.9848,-15.9892,-15.9948,-15.9992,-16.0006,-16.0061,-16.0128,-16.0178,-16.0225,-16.0275,-16.0331,-16.0375,-16.0417,-16.0484,-16.0481,-16.0517,-16.0561,-16.0598,-16.0636,-16.0667,-16.0689,-16.0711,-16.0725,-16.0726,-16.0756,-16.0789,-16.0817,-16.0848,-16.0889,-16.0931,-16.0967,-16.0989,-16.0986,-16.0995,-16.102,-16.105,-16.1089,-16.1125,-16.1175,-16.1234,-16.1306,-16.1375,-16.145,-16.1506,-16.1548,-16.1589,-16.1623,-16.1656,-16.1698,-16.1734,-16.1767,-16.1773,-16.1853,-16.1957,-16.208,-16.2223,-16.2406,-16.249,-16.2585,-16.2693,-16.2717,-16.2741,-16.2796,-16.2808,-16.2836,-16.2888,-16.2964,-16.3011,-16.3043,-16.3051,-16.3047,-16.3079,-16.3127,-16.3175,-16.3218,-16.3262,-16.329,-16.3318,-16.3318,-16.3306,-16.3318,-16.3346,-16.3378,-16.3534,-16.3565,-16.3604,-16.3632,-16.3616,-16.3612,-16.3577,-16.3517,-16.3497,-16.3493,-16.3513,-16.3608,-16.3676,-16.3716,-16.3748,-16.3807,-16.3887,-16.3903,-16.3935,-16.4006,-16.4094,-16.4146,-16.4209,-16.4313,-16.4377,-16.4417,-16.4424,-16.4464,-16.4484,-16.4488,-16.448,-16.4484,-16.4496,-16.45,-16.4484,-16.448,-16.4508,-16.4512,-16.4512,-16.448,-16.4444,-16.4436,-16.444,-16.4448,-16.4468,-16.4504,-16.4576,-16.4639,-16.4666,-16.4666,-16.468,-16.4692,-16.468,-16.467,-16.468,-16.468,-16.4651,-16.4635,-16.4618,-16.4593,-16.4585,-16.456,-16.4543,-16.4531,-16.4535,-16.4551,-16.456,-16.4601,-16.4634,-16.4651,-16.466,-16.4676,-16.4697,-16.4697,-16.4705,-16.4705,-16.4722,-16.4755,-16.4764,-16.4739,-16.4739,-16.473,-16.473,-16.4743,-16.4763,-16.4772,-16.481,-16.4868,-16.4876,-16.4901,-16.4918,-16.493,-16.4935,-16.4976,-16.5005,-16.5005,-16.5013,-16.5013,-16.4963,-16.4947,-16.4946,-16.4955,-16.498,-16.4972,-16.498,-16.498,-16.4963,-16.4963,-16.498,-16.498,-16.4967,-16.4958,-16.4921,-16.4913,-16.4913,-16.4892,-16.4859,-16.4834,-16.4817,-16.4796,-16.4788,-16.4788,-16.4775,-16.4771,-16.4788,-16.4788,-16.4775,-16.4742,-16.4734,-16.4717,-16.4693,-16.468,-16.4671,-16.4684,-16.47,-16.4709,-16.4742,-16.4759,-16.4775,-16.4809,-16.4855,-16.4855,-16.483,-16.483,-16.4838,-16.4838,-16.4867,-16.4892,-16.4909,-16.4926,-16.4955,-16.4946,-16.4946,-16.4938,-16.4938,-16.493,-16.4938,-16.4938,-16.4921,-16.4938,-16.4922,-16.4921,-16.4913,-16.4921,-16.4913,-16.4913,-16.4905,-16.4905,-16.4892,-16.4851,-16.4842,-16.4838,-16.4846,-16.4846,-16.4859,-16.4871,-16.4838,-16.4838,-16.4855,-16.4855,-16.4859,-16.4871,-16.4867,-16.4865,-16.4863,-16.4876,-16.4888,-16.4896,-16.4905,-16.4905,-16.4918,-16.4919,-16.4921,-16.4904,-16.4905,-16.4913,-16.4913,-16.4951,-16.4988,-16.4988,-16.5005,-16.5005,-16.5005,-16.4996,-16.4996,-16.5013,-16.5005,-16.5005,-16.4988,-16.4988,-16.4996,-16.4988,-16.5,-16.5026,-16.5038,-16.503,-16.503,-16.5046,-16.5038,-16.5047,-16.503,-16.5046,-16.5038,-16.503,-16.5038,-16.5038,-16.503,-16.5038,-16.5038,-16.5046,-16.5046,-16.5055,-16.5055,-16.5046,-16.5046,-16.5071,-16.5071,-16.5063,-16.5063,-16.5071,-16.5071,-16.5088,-16.5088,-16.508,-16.508,-16.5088,-16.5088,-16.5097,-16.5105,-16.5113,-16.5121,-16.513,-16.5129,-16.5138,-16.5146,-16.5155,-16.5155,-16.5142,-16.5138,-16.5128,-16.5113,-16.5113,-16.5096,-16.5096,-16.508,-16.5071,-16.5071,-16.508,-16.5088,-16.5104,-16.5105,-16.513,-16.513,-16.5155,-16.5163,-16.5171,-16.5171,-16.518,-16.5196,-16.5205,-16.5213,-16.5238,-16.5238,-16.5247,-16.5255,-16.528,-16.5297,-16.533,-16.533,-16.5355,-16.5363,-16.5371,-16.538,-16.5413,-16.543,-16.5438,-16.5471,-16.5471,-16.5488,-16.5488,-16.5505,-16.5505,-16.5538,-16.5546,-16.5563,-16.5563,-16.558,-16.5588,-16.5605,-16.5605,-16.5638,-16.5671,-16.5671,-16.5721,-16.573,-16.5755,-16.5763,-16.5788,-16.5813,-16.588,-16.5888,-16.5946,-16.5946,-16.5955,-16.6022,-16.6071,-16.6096,-16.613,-16.6138,-16.6155,-16.6163,-16.618,-16.6213,-16.6271,-16.6271,-16.633,-16.6355,-16.6438,-16.6463,-16.6513,-16.6513,-16.6538,-16.6546,-16.6571,-16.6571,-16.6605,-16.6613,-16.6655,-16.6655,-16.6721,-16.6721,-16.6755,-16.678,-16.6805,-16.6805,-16.6847,-16.6863,-16.6905,-16.6905,-16.6955,-16.6963,-16.6988,-16.6996,-16.7021,-16.7038,-16.7071,-16.7105,-16.7163,-16.7163,-16.7246,-16.7263,-16.7271,-16.7305,-16.7307,-16.7102,-16.7012,-16.6958,-16.6911,-16.6901,-16.6892,-16.6886,-16.6888,-16.6883,-16.6887,-16.6889,-16.6908,-16.6973,-16.6995,-16.7031,-16.7033,-16.7075,-16.7127,-16.7155,-16.7187,-16.7209,-16.7231,-16.7247,-16.7253,-16.7274,-16.7277,-16.7285,-16.7337,-16.7338,-16.7373,-16.738,-16.7401,-16.7401,-16.7423,-16.7432,-16.7415,-16.7405,-16.7381,-16.7366,-16.7333,-16.7319,-16.7305,-16.7233,-16.718,-16.7167,-16.7155,-16.7136,-16.7094,-16.7081,-16.7073,-16.7048,-16.7022,-16.701,-16.6961,-16.6922,-16.6907,-16.6867,-16.6844,-16.6828,-16.6811,-16.676,-16.6745,-16.673,-16.6701,-16.6687,-16.6669,-16.6635,-16.6626,-16.6608,-16.6557,-16.6536,-16.6476,-16.6428,-16.6423,-16.6412,-16.6384,-16.6365,-16.6338,-16.6312,-16.6307,-16.6306,-16.6285,-16.6245,-16.6207,-16.6165,-16.6133,-16.6121,-16.6121,-16.6123,-16.6129,-16.6134,-16.614,-16.6141,-16.6128,-16.6092,-16.6072,-16.6045,-16.6037,-16.6012,-16.5993,-16.5958,-16.5907,-16.5836,-16.5824,-16.581,-16.5722,-16.5696,-16.5671,-16.5647,-16.5596,-16.5571,-16.554,-16.5519,-16.5471,-16.5451,-16.5418,-16.54,-16.5381,-16.535,-16.5332,-16.5307,-16.5279,-16.5251,-16.5237,-16.5225,-16.5202,-16.5186,-16.5173,-16.5158,-16.5145,-16.5135,-16.5132,-16.5129,-16.5124,-16.5123,-16.5123,-16.5135,-16.5143,-16.5156,-16.5164,-16.5162,-16.5158,-16.5146,-16.5137,-16.5124,-16.5111,-16.5096,-16.5077,-16.506,-16.5044,-16.503,-16.502,-16.5007,-16.499,-16.4974,-16.4954,-16.4941,-16.4915,-16.4901,-16.4889,-16.4873,-16.4861,-16.4845,-16.4832,-16.4817,-16.4804,-16.4791,-16.4779,-16.4759,-16.4751,-16.4736,-16.4709,-16.4691,-16.4659,-16.4641,-16.4618,-16.4594,-16.4564,-16.4546,-16.4538,-16.4526,-16.4515,-16.4509,-16.4483,-16.4458,-16.4436,-16.4414,-16.4385,-16.4362,-16.4338,-16.4311,-16.428,-16.4255,-16.4228,-16.4189,-16.4153,-16.4107,-16.4039,-16.3979,-16.3932,-16.3878,-16.3839,-16.3796,-16.3773,-16.3752,-16.3732,-16.3719,-16.3693,-16.366,-16.363,-16.3587,-16.3547,-16.3515,-16.3496,-16.3468,-16.3448,-16.3421,-16.3411,-16.3393,-16.3371,-16.3356,-16.3346,-16.3344,-16.3354,-16.3371,-16.3388,-16.3416,-16.3446,-16.3454,-16.3457,-16.345,-16.3436,-16.3415,-16.338,-16.335,-16.3332,-16.3305,-16.3274,-16.3237,-16.3202,-16.3171,-16.3142,-16.3112,-16.3102,-16.3093,-16.3088,-16.3084,-16.308,-16.3064,-16.3045,-16.3017,-16.3,-16.2969,-16.295,-16.2931,-16.2908,-16.2885,-16.2868,-16.2858,-16.2849,-16.2832,-16.2806,-16.278,-16.276,-16.2742,-16.2723,-16.2708,-16.27,-16.2692,-16.2685,-16.268,-16.2681,-16.2685,-16.2689,-16.2692,-16.2697,-16.2699,-16.2693,-16.2688,-16.2683,-16.2672,-16.2661,-16.2647,-16.2632,-16.2616,-16.2601,-16.2587,-16.2572,-16.2551,-16.2526,-16.2499,-16.2475,-16.245,-16.2438,-16.2426,-16.2416,-16.2392,-16.2376,-16.2373,-16.2372,-16.2373,-16.2377,-16.2374,-16.2371,-16.2364,-16.2345,-16.2321,-16.2304,-16.227,-16.2249,-16.2188,-16.2136,-16.2116,-16.2086,-16.205,-16.2041,-16.2012,-16.1989,-16.1968,-16.1962,-16.1951,-16.1947,-16.1947,-16.1947,-16.1947,-16.1947,-16.1942,-16.1943,-16.1944,-16.1931,-16.1921,-16.1904,-16.1891,-16.1886,-16.1873,-16.1858,-16.1832,-16.1824,-16.1808,-16.18,-16.1796,-16.1789,-16.1784,-16.1782,-16.1783,-16.1787,-16.1782,-16.1781,-16.1779,-16.1776,-16.1771,-16.1767,-16.1764,-16.1759,-16.1755,-16.1761,-16.1766,-16.1771,-16.1771,-16.1772,-16.1771,-16.177,-16.177,-16.176,-16.1735,-16.1705,-16.1674,-16.1638,-16.1605,-16.1577,-16.1544,-16.1522,-16.1486,-16.147,-16.1422,-16.1395,-16.1368,-16.1345,-16.1317,-16.1284,-16.1235,-16.1213,-16.1177,-16.1158,-16.1114,-16.1095,-16.1074,-16.1049,-16.1025,-16.0992,-16.0971,-16.0948,-16.0938,-16.0921,-16.0905,-16.0879,-16.0863,-16.0843,-16.0816,-16.0794,-16.0759,-16.0728,-16.0703,-16.0681,-16.0658,-16.0602,-16.0577,-16.0545,-16.0521,-16.0495,-16.0469,-16.0439,-16.0408,-16.0377,-16.0361,-16.0348,-16.0321,-16.0264,-16.0176,-16.0156,-16.0139,-16.0099,-16.0078,-16.0058,-16.004,-16.0018,-16.0001,-15.999,-15.9973,-15.9968,-15.9962,-15.9954,-15.9944,-15.9936,-15.9929,-15.9922,-15.9916,-15.9911,-15.9904,-15.9899,-15.9891,-15.9884,-15.9873,-15.9854,-15.9811,-15.9777,-15.9744,-15.972,-15.9688,-15.963,-15.9575,-15.9532,-15.9479,-15.9434,-15.9369,-15.9331,-15.9276,-15.9244,-15.9224,-15.9196,-15.9171,-15.9163,-15.9157,-15.9147,-15.9133,-15.9103,-15.9075,-15.9048,-15.9026,-15.9009,-15.8992,-15.8976,-15.8962,-15.8942,-15.8915,-15.8837,-15.877,-15.8686,-15.8635,-15.8626,-15.8623,-15.8626,-15.8633,-15.8644,-15.8651,-15.8646,-15.862,-15.8588,-15.855,-15.8515,-15.8467,-15.842,-15.8385,-15.8342,-15.8305,-15.8264,-15.8234,-15.8201,-15.8161,-15.8104,-15.8034,-15.8028,-15.8023,-15.8018,-15.8002,-15.7989,-15.7962,-15.793,-15.7892,-15.7868,-15.7851,-15.7829,-15.7815,-15.7798,-15.7779,-15.7765,-15.775,-15.7737,-15.7729,-15.7727,-15.772,-15.7715,-15.7699,-15.7678,-15.7659,-15.7633,-15.7617,-15.7603,-15.7584,-15.7552,-15.7508,-15.7457,-15.7413,-15.7369,-15.7322,-15.7274,-15.7223,-15.7158,-15.7117,-15.7066,-15.7003,-15.6946,-15.687,-15.6791,-15.6722,-15.6664,-15.6628,-15.6554,-15.6498,-15.6435,-15.639,-15.6356,-15.6354,-15.6406,-15.6525,-15.6621,-15.6751,-15.6876,-15.692,-15.6903,-15.6841,-15.6773,-15.676,-15.6755,-15.6752,-15.6749,-15.6744,-15.6739,-15.6735,-15.673,-15.6726,-15.6711,-15.6689,-15.6666,-15.6649,-15.6624,-15.6608,-15.6587,-15.6565,-15.6544,-15.6509,-15.6484,-15.6457,-15.6434,-15.6408,-15.6379,-15.6351,-15.6322,-15.6307,-15.6272,-15.6251,-15.6218,-15.6185,-15.6148,-15.6122,-15.6099,-15.6078,-15.6064,-15.603,-15.6009,-15.5991,-15.5962,-15.5937,-15.5909,-15.5868,-15.5846,-15.5815,-15.5785,-15.5773,-15.5757,-15.5736,-15.5698,-15.5665,-15.563,-15.5597,-15.5562,-15.5523,-15.5496,-15.5452,-15.5432,-15.5408,-15.5386,-15.5364,-15.5348,-15.5313,-15.5276,-15.5238,-15.5204,-15.5175,-15.5148,-15.513,-15.5103,-15.5078,-15.5054,-15.5046,-15.503,-15.501,-15.499,-15.4973,-15.4953,-15.4923,-15.4904,-15.4884,-15.4864,-15.4844,-15.4817,-15.4786,-15.4759,-15.4731,-15.471,-15.4686,-15.4672,-15.4646,-15.4617,-15.4594,-15.4569,-15.4532,-15.4493,-15.4449,-15.4425,-15.4397,-15.4371,-15.4357,-15.4357,-15.4358,-15.4353,-15.4361,-15.436,-15.437,-15.4371,-15.4377,-15.4378,-15.4392,-15.4395,-15.4402,-15.4405,-15.4408,-15.4416,-15.4419,-15.4423,-15.443,-15.4435,-15.4446,-15.4454,-15.4458,-15.4185,-15.3912,-15.3639,-15.3366,-15.3093,-15.2807,-15.252,-15.2234,-15.1947,-15.1661,-15.1375,-15.1088,-15.0887,-15.0855,-15.0802,-15.0515,-15.0229,-15.0162,-14.9989,-14.9896,-14.9629,-14.9362,-14.9095,-14.8891,-14.8686,-14.8496,-14.836,-14.8203,-14.8111,-14.8036,-14.8028,-14.8027,-14.8002,-14.7981,-14.796,-14.7922,-14.7822,-14.7822,-14.7734,-14.7733,-14.7711,-14.7703,-14.7693,-14.7608,-14.7534,-14.7462,-14.7379,-14.7345,-14.7323,-14.7315,-14.7279,-14.7277,-14.7191,-14.7191,-14.7056,-14.7051,-14.6926,-14.6888,-14.683,-14.6642,-14.6634,-14.6628,-14.659,-14.6587,-14.6581,-14.6577,-14.6514,-14.6329,-14.6271,-14.6076,-14.6053,-14.5993,-14.5894,-14.5821,-14.5738,-14.5733,-14.5717,-14.5674,-14.5658,-14.5604,-14.5566,-14.5484,-14.5413,-14.5316,-14.5247,-14.519,-14.513,-14.5025,-14.4985,-14.4976,-14.494,-14.494,-14.4939,-14.4937,-14.4913,-14.4904,-14.4895,-14.4882,-14.4879,-14.4858,-14.4846,-14.4837,-14.4834,-14.4824,-14.4814,-14.4799,-14.4776,-14.4755,-14.4739,-14.4713,-14.4707,-14.4703,-14.4693,-14.4666,-14.456,-14.4472,-14.4356,-14.428,-14.4166,-14.408,-14.3998,-14.3799,-14.3683,-14.3681,-14.3515,-14.3513,-14.3302,-14.3301,-14.3301,-14.33,-14.3286,-14.3236,-14.2941,-14.2904,-14.2846,-14.2775,-14.2414,-14.2414,-14.2385,-14.2384,-14.1875,-14.1851,-14.1825,-14.1792,-14.1792,-14.1642,-14.1565,-14.1564,-14.1545,-14.1062,-14.0565,-14.0531,-14.034,-14.0253,-14.0016,-13.9989,-13.9948,-13.9737,-13.9645,-13.9482,-13.9227,-13.8972,-13.8717,-13.8462,-13.8194,-13.8163,-13.7886,-13.7874,-13.7868,-13.7864,-13.7834,-13.7564,-13.7265,-13.7162,-13.6966,-13.6667,-13.6368,-13.6069,-13.577,-13.5471,-13.5323,-13.5172,-13.5021,-13.4869,-13.4717,-13.4565,-13.4413,-13.4261,-13.4109,-13.3958,-13.3806,-13.3654,-13.3502,-13.335,-13.3198,-13.3046,-13.2894,-13.2742,-13.259,-13.2438,-13.2141,-13.2106,-13.1845,-13.1548,-13.1252,-13.0956,-13.0669,-13.0659,-13.0363,-13.0174,-12.9989,-12.9745,-12.9537,-12.9506,-12.9247,-12.8988,-12.8729,-12.847,-12.8211,-12.7952,-12.7701,-12.7449,-12.7198,-12.703,-12.6861,-12.6693,-12.6524,-12.651,-12.6429,-12.6334,-12.6239,-12.6147,-12.6144,-12.6049,-12.5954,-12.5942,-12.5859,-12.5866,-12.587,-12.5878,-12.5885,-12.5892,-12.5904,-12.591,-12.5914,-12.592,-12.5928,-12.5937,-12.594,-12.5948,-12.5953,-12.596,-12.5967,-12.597,-12.5987,-12.5992,-12.6,-12.6006,-12.6006,-12.6001,-12.6001,-12.6005,-12.6005,-12.6011,-12.6011,-12.6007,-12.6013,-12.6008,-12.6006,-12.6014,-12.6019,-12.602,-12.6022,-12.6022,-12.6026,-12.6031,-12.6032,-12.6038,-12.6042,-12.6042,-12.6049,-12.6053,-12.6051,-12.6058,-12.6066,-12.6074,-12.6087,-12.6098,-12.6107,-12.6118,-12.6128,-12.6129,-12.6126,-12.612,-12.611,-12.6102,-12.6117,-12.6126,-12.614,-12.6144,-12.6156,-12.6155,-12.6152,-12.615,-12.6149,-12.6141,-12.6141,-12.6146,-12.6155,-12.6163,-12.6177,-12.6192,-12.6205,-12.6213,-12.6209,-12.6212,-12.6211,-12.6213,-12.6211,-12.6219,-12.6228,-12.6235,-12.6237,-12.6247,-12.6444,-12.6641,-12.6861,-12.6911,-12.6961,-12.6998,-12.7036,-12.7073,-12.7103,-12.7134,-12.717,-12.72,-12.7236,-12.7281,-12.7317,-12.7361,-12.7403,-12.7448,-12.7498,-12.7542,-12.7592,-12.7636,-12.7686,-12.7803,-12.7811,-12.7825,-12.785,-12.7845,-12.7825,-12.7825,-12.7828,-12.7831,-12.7834,-12.785,-12.7873,-12.7909,-12.7939,-12.7981,-12.8075,-12.815,-12.8209,-12.8259,-12.8314,-12.837,-12.8436,-12.8489,-12.8525,-12.8553,-12.857,-12.8606,-12.8645,-12.8636,-12.8673,-12.8717,-12.8753,-12.8784,-12.882,-12.8856,-12.8886,-12.8903,-12.8925,-12.8898,-12.8856,-12.8792,-12.8734,-12.8686,-12.8628,-12.857,-12.8506,-12.8464,-12.8431,-12.8417,-12.8406,-12.8392,-12.8392,-12.8403,-12.8425,-12.845,-12.8486,-12.8514,-12.855,-12.8581,-12.862,-12.8664,-12.8706,-12.8742,-12.8786,-12.8836,-12.8881,-12.8934,-12.8975,-12.902,-12.9064,-12.9114,-12.9164,-12.9217,-12.9253,-12.9289,-12.9317,-12.9342,-12.9359,-12.9375,-12.9392,-12.9395,-12.9395,-12.9398,-12.9392,-12.9386,-12.9375,-12.9356,-12.9328,-12.9336,-12.9386,-12.9423,-12.9475,-12.9492,-12.9492,-12.9509,-12.9523,-12.9523,-12.9539,-12.9556,-12.955,-12.9553,-12.9548,-12.9561,-12.9566,-12.9595,-12.9639,-12.9709,-12.9764,-12.9814,-12.9864,-12.9928,-12.9984,-12.9992,-13.0034,-13.0084,-13.0148,-13.0211,-13.0275,-13.0325,-13.0373,-13.0417,-13.0467,-13.0523,-13.0592,-13.065,-13.0686,-13.0731,-13.0767,-13.082,-13.0864,-13.0909,-13.0931,-13.0932,-13.0936,-13.0967,-13.0992,-13.1031,-13.1036,-13.1031,-13.1017,-13.0989,-13.0964,-13.0931,-13.0903,-13.0889,-13.0886,-13.0895,-13.0898,-13.0906,-13.0914,-13.0959,-13.1011,-13.1053,-13.1098,-13.1148,-13.1192,-13.1234,-13.1278,-13.1314,-13.1367,-13.1403,-13.1453,-13.1506,-13.1564,-13.1614,-13.165,-13.1686,-13.1731,-13.1781,-13.1839,-13.1903,-13.1967,-13.2023,-13.2075,-13.2131,-13.2181,-13.225,-13.23,-13.2339,-13.2361,-13.2392,-13.2417,-13.2445,-13.2467,-13.2484,-13.25,-13.2503,-13.2489,-13.247,-13.245,-13.2417,-13.2361,-13.2295,-13.2239,-13.2211,-13.2192,-13.2186,-13.2217,-13.2253,-13.2289,-13.2325,-13.2359,-13.2395,-13.2439,-13.2475,-13.2511,-13.2548,-13.2592,-13.2636,-13.2673,-13.2711,-13.2742,-13.2764,-13.2778,-13.2809,-13.2831,-13.2848,-13.2859,-13.2859,-13.2839,-13.2789,-13.2748,-13.2714,-13.2728,-13.2764,-13.2806,-13.2842,-13.2878,-13.2906,-13.2945,-13.2967,-13.2989,-13.3006,-13.3025,-13.3039,-13.3048,-13.305,-13.3053,-13.3064,-13.3078,-13.3103,-13.3131,-13.3161,-13.3198,-13.3198,-13.3242,-13.3259,-13.3253,-13.3248,-13.3248,-13.3245,-13.3239,-13.327,-13.33,-13.3336,-13.3381,-13.3417,-13.347,-13.3511,-13.3556,-13.3592,-13.3625,-13.3648,-13.3678,-13.37,-13.3731,-13.3753,-13.3775,-13.3786,-13.3803,-13.3804,-13.3804,-13.3817,-13.3828,-13.3836,-13.3848,-13.3856,-13.3873,-13.3889,-13.3906,-13.392,-13.3945,-13.3973,-13.4017,-13.4061,-13.4114,-13.4156,-13.4217,-13.4264,-13.4309,-13.4359,-13.4411,-13.4448,-13.4492,-13.4523,-13.4559,-13.4598,-13.462,-13.465,-13.4681,-13.4714,-13.4748,-13.4778,-13.4811,-13.4839,-13.4873,-13.49,-13.4942,-13.4992,-13.5006,-13.505,-13.5081,-13.5111,-13.5134,-13.515,-13.5145,-13.5134,-13.5098,-13.5056,-13.5014,-13.4992,-13.4973,-13.4923,-13.4886,-13.4861,-13.4842,-13.485,-13.4881,-13.4925,-13.4945,-13.4992,-13.5003,-13.5042,-13.5081,-13.5106,-13.5148,-13.5198,-13.522,-13.5284,-13.5348,-13.5411,-13.5478,-13.5534,-13.56,-13.5659,-13.5723,-13.5792,-13.5859,-13.5909,-13.5964,-13.6014,-13.6056,-13.6106,-13.6148,-13.6203,-13.6253,-13.6325,-13.6389,-13.6445,-13.6489,-13.6523,-13.6556,-13.6592,-13.6639,-13.6703,-13.6761,-13.6786,-13.6823,-13.6853,-13.6889,-13.6928,-13.6964,-13.7009,-13.7053,-13.7089,-13.712,-13.712,-13.7109,-13.7089,-13.7061,-13.7036,-13.7009,-13.6981,-13.6979,-13.697,-13.697,-13.6995,-13.7039,-13.7103,-13.7153,-13.72,-13.7242,-13.7278,-13.732,-13.7361,-13.7398,-13.7439,-13.7464,-13.7506,-13.7548,-13.76,-13.7656,-13.7711,-13.7784,-13.7848,-13.7911,-13.7967,-13.8025,-13.8084,-13.8131,-13.8173,-13.8223,-13.8267,-13.8309,-13.8348,-13.8384,-13.8425,-13.8473,-13.8506,-13.8528,-13.855,-13.8573,-13.8617,-13.865,-13.8673,-13.8681,-13.8692,-13.8709,-13.8703,-13.8706,-13.8709,-13.8711,-13.872,-13.8742,-13.8778,-13.8845,-13.8895,-13.8953,-13.9011,-13.9061,-13.9092,-13.9114,-13.9153,-13.9181,-13.9225,-13.927,-13.9306,-13.935,-13.94,-13.9459,-13.95,-13.9531,-13.9595,-13.9653,-13.9689,-13.9725,-13.9728,-13.9717,-13.9698,-13.9684,-13.9664,-13.9639,-13.9639,-13.9661,-13.9706,-13.975,-13.9786,-13.9825,-13.9842,-13.9836,-13.9817,-13.9781,-13.9725,-13.9678,-13.9636,-13.9606,-13.9631,-13.9667,-13.9711,-13.975,-13.9786,-13.9823,-13.9861,-13.9898,-13.9936,-13.9981,-13.9992,-14.0036,-14.0103,-14.0173,-14.0239,-14.0289,-14.0325,-14.0356,-14.0392,-14.0417,-14.0434,-14.0464,-14.0486,-14.0531,-14.0573,-14.062,-14.0664,-14.07,-14.0736,-14.0773,-14.0811,-14.0853,-14.09,-14.0945,-14.0981,-14.1025,-14.1061,-14.11,-14.1114,-14.1148,-14.1175,-14.1211,-14.1256,-14.13,-14.1339,-14.1375,-14.1406,-14.1431,-14.1453,-14.1484,-14.1528,-14.1578,-14.1623,-14.1667,-14.1703,-14.1728,-14.1764,-14.1785,-14.1785,-14.1795,-14.1831,-14.1861,-14.1898,-14.1942,-14.1981,-14.2011,-14.2039,-14.2073,-14.21,-14.2125,-14.2161,-14.2206,-14.2242,-14.2295,-14.2339,-14.2375,-14.2409,-14.2453,-14.2495,-14.2528,-14.2564,-14.2592,-14.2625,-14.2684,-14.2728,-14.2759,-14.2775,-14.277,-14.2773,-14.2767,-14.2764,-14.2764,-14.2775,-14.2825,-14.2889,-14.2945,-14.2998,-14.3039,-14.3086,-14.3136,-14.3195,-14.325,-14.3303,-14.3348,-14.3378,-14.3389,-14.337,-14.3356,-14.3345,-14.3348,-14.3342,-14.3336,-14.3331,-14.3361,-14.34,-14.3445,-14.35,-14.3545,-14.355,-14.3625,-14.3695,-14.3767,-14.3834,-14.3903,-14.3953,-14.3992,-14.3992,-14.4006,-14.4056,-14.41,-14.4153,-14.4189,-14.4236,-14.4292,-14.4342,-14.4392,-14.4425,-14.4453,-14.4478,-14.4498,-14.4534,-14.4584,-14.4656,-14.472,-14.4792,-14.4856,-14.4911,-14.4978,-14.4992,-14.5034,-14.5106,-14.517,-14.5228,-14.5286,-14.5345,-14.5395,-14.5434,-14.5478,-14.5523,-14.5586,-14.5628,-14.5661,-14.5706,-14.5753,-14.5717,-14.5784,-14.5853,-14.5925,-14.5984,-14.6042,-14.6184,-14.6256,-14.6328,-14.6373,-14.6389,-14.6406,-14.6442,-14.6514,-14.6578,-14.6642,-14.6709,-14.6773,-14.6828,-14.6886,-14.6942,-14.7,-14.7064,-14.7123,-14.7186,-14.7239,-14.7281,-14.7306,-14.7342,-14.7398,-14.7448,-14.7473,-14.7495,-14.7525,-14.7556,-14.7595,-14.765,-14.7709,-14.7745,-14.7792,-14.7834,-14.7886,-14.7942,-14.7992,-14.8031,-14.8067,-14.8111,-14.817,-14.8228,-14.8267,-14.8309,-14.8361,-14.8417,-14.8481,-14.8539,-14.8595,-14.861,-14.8639,-14.8678,-14.8714,-14.8753,-14.8795,-14.8839,-14.8884,-14.8917,-14.8945,-14.8978,-14.9014,-14.9061,-14.912,-14.9142,-14.9167,-14.9234,-14.9281,-14.9345,-14.9403,-14.9409,-14.9459],"lat":[16.6424,16.6449,16.6496,16.6552,16.6613,16.6674,16.6721,16.6735,16.6796,16.6843,16.689,16.6902,16.6921,16.6927,16.691,16.6879,16.6852,16.6815,16.6795,16.6779,16.6738,16.669,16.664,16.659,16.6543,16.6499,16.6457,16.6415,16.6379,16.6352,16.6335,16.6313,16.6296,16.6299,16.6313,16.6346,16.6377,16.6418,16.6465,16.6513,16.656,16.6607,16.664,16.6674,16.6704,16.6718,16.6707,16.6685,16.6649,16.6621,16.6585,16.6543,16.6479,16.6424,16.6371,16.6329,16.6304,16.6288,16.6268,16.6249,16.6202,16.614,16.6077,16.6035,16.5999,16.5963,16.5927,16.5893,16.5857,16.5835,16.5824,16.5829,16.5832,16.5846,16.5857,16.5852,16.5857,16.5849,16.5824,16.5802,16.5793,16.5771,16.5754,16.5738,16.5721,16.5707,16.5685,16.566,16.5638,16.5615,16.5593,16.5585,16.5588,16.5599,16.5618,16.5629,16.5649,16.5665,16.5679,16.5696,16.5715,16.5727,16.5732,16.5721,16.5721,16.571,16.5696,16.5671,16.5649,16.5635,16.561,16.561,16.5613,16.5596,16.5582,16.556,16.5543,16.5521,16.5496,16.5468,16.5452,16.5443,16.5432,16.5438,16.5463,16.5496,16.5529,16.5563,16.561,16.5657,16.5713,16.5752,16.579,16.5818,16.5835,16.5827,16.581,16.5793,16.5771,16.5749,16.5749,16.5727,16.571,16.5693,16.5671,16.567,16.5649,16.5613,16.5563,16.5502,16.544,16.5385,16.5327,16.5279,16.5243,16.5213,16.5185,16.5168,16.5146,16.5129,16.5121,16.5118,16.5121,16.5127,16.5129,16.514,16.5152,16.5171,16.519,16.5207,16.5227,16.5229,16.5221,16.5185,16.5129,16.5074,16.5032,16.4997,16.4996,16.4982,16.4954,16.4924,16.4902,16.4885,16.4863,16.4846,16.4824,16.4807,16.4799,16.479,16.4785,16.479,16.4802,16.4813,16.4824,16.4843,16.4854,16.4865,16.4877,16.4882,16.4885,16.4896,16.4893,16.4904,16.4918,16.4929,16.494,16.4957,16.4968,16.4988,16.4991,16.4996,16.4997,16.4998,16.5007,16.5018,16.5024,16.5013,16.4998,16.4996,16.499,16.4968,16.4957,16.4949,16.4946,16.4957,16.4982,16.4996,16.4997,16.4997,16.501,16.504,16.5068,16.5093,16.5113,16.5124,16.5127,16.5124,16.5121,16.5104,16.5082,16.5065,16.5043,16.5021,16.4999,16.4997,16.4996,16.4974,16.4946,16.4915,16.4893,16.4863,16.486,16.4874,16.4885,16.4896,16.4921,16.4932,16.493,16.4929,16.4927,16.4918,16.4893,16.4871,16.4849,16.4832,16.4829,16.4815,16.4804,16.4799,16.4804,16.4829,16.4849,16.4882,16.4921,16.4954,16.4988,16.4996,16.4997,16.5015,16.504,16.5065,16.5093,16.5118,16.5152,16.5185,16.5232,16.5293,16.534,16.5388,16.5429,16.546,16.5493,16.5513,16.5524,16.5521,16.5518,16.5515,16.5499,16.5468,16.544,16.5404,16.5368,16.5338,16.5302,16.5265,16.5259,16.5231,16.5195,16.5183,16.5211,16.5327,16.5334,16.5311,16.5239,16.5205,16.5171,16.5052,16.4936,16.4797,16.4745,16.4698,16.4642,16.4566,16.4491,16.4109,16.4037,16.4021,16.4029,16.4029,16.3997,16.3941,16.375,16.3703,16.3627,16.3512,16.344,16.34,16.3307,16.3289,16.3241,16.3149,16.3014,16.3,16.2891,16.2807,16.2763,16.2723,16.2668,16.2584,16.2473,16.2361,16.227,16.221,16.2166,16.2161,16.215,16.2154,16.2174,16.2182,16.2178,16.2126,16.2079,16.2028,16.2019,16.1943,16.1864,16.1784,16.1673,16.1633,16.1537,16.1494,16.1438,16.1374,16.1287,16.1259,16.1179,16.1135,16.1119,16.1095,16.1028,16.0984,16.094,16.0892,16.0877,16.0845,16.0826,16.0825,16.0804,16.0796,16.0789,16.0776,16.0763,16.0713,16.0684,16.0684,16.07,16.07,16.0692,16.0692,16.0683,16.0671,16.0667,16.0667,16.0659,16.0658,16.0642,16.0642,16.0634,16.0634,16.0604,16.0579,16.0571,16.0471,16.0437,16.0404,16.0371,16.0346,16.0329,16.0321,16.0304,16.0292,16.0304,16.0337,16.0375,16.0375,16.0367,16.0358,16.0334,16.0329,16.0317,16.0308,16.0271,16.0245,16.0237,16.0179,16.0129,16.0104,15.9996,15.997,15.9954,15.9921,15.9913,15.9896,15.9879,15.9787,15.977,15.9763,15.975,15.975,15.9779,15.9796,15.9829,15.985,15.9867,15.9867,15.9875,15.9896,15.9912,15.9937,15.9942,15.9913,15.9896,15.9887,15.9867,15.9858,15.985,15.985,15.9883,15.9879,15.9829,15.9817,15.9817,15.9825,15.9825,15.9817,15.9842,15.9842,15.9813,15.9795,15.9762,15.9721,15.9712,15.9696,15.9666,15.9667,15.965,15.965,15.9629,15.9612,15.9538,15.9529,15.9504,15.9495,15.9487,15.9471,15.9454,15.9437,15.9421,15.9404,15.9396,15.938,15.9371,15.9321,15.9313,15.9296,15.9284,15.9284,15.9292,15.9279,15.927,15.9262,15.9258,15.9238,15.9187,15.9137,15.9105,15.9037,15.9034,15.9046,15.9095,15.911,15.9129,15.9134,15.9121,15.9071,15.9063,15.9046,15.9017,15.9063,15.9137,15.9171,15.9221,15.9229,15.9254,15.9258,15.9296,15.9329,15.9338,15.946,15.9487,15.9496,15.9571,15.9596,15.9605,15.9637,15.9663,15.9746,15.9754,15.9788,15.98,15.9792,15.9779,15.9763,15.9696,15.9671,15.9654,15.9638,15.9621,15.9588,15.9546,15.9538,15.9529,15.9512,15.9505,15.9496,15.9471,15.9462,15.9363,15.9355,15.9271,15.9262,15.9221,15.9179,15.9154,15.9146,15.9062,15.9054,15.9038,15.9021,15.8988,15.8979,15.8912,15.8904,15.8813,15.8804,15.8738,15.8729,15.8596,15.8588,15.8563,15.8554,15.8513,15.8504,15.8479,15.8475,15.8479,15.8502,15.8537,15.8563,15.8587,15.8621,15.8662,15.8662,15.8604,15.8596,15.8562,15.8545,15.8488,15.8463,15.8454,15.8429,15.8346,15.8338,15.8296,15.8287,15.8221,15.8213,15.8163,15.8112,15.8079,15.8071,15.8029,15.7979,15.7962,15.7904,15.7887,15.7863,15.7829,15.782,15.7771,15.7712,15.7696,15.7671,15.7621,15.7604,15.7579,15.7563,15.7537,15.7521,15.7471,15.7437,15.7412,15.7396,15.7371,15.7338,15.7321,15.7304,15.7263,15.7196,15.7179,15.7087,15.7054,15.7013,15.6971,15.6921,15.687,15.6787,15.6763,15.668,15.6662,15.6654,15.6521,15.6421,15.6387,15.6337,15.6313,15.6296,15.627,15.6254,15.6187,15.6129,15.6121,15.6045,15.5988,15.5854,15.5813,15.5754,15.5746,15.5721,15.5696,15.5671,15.5663,15.5621,15.5596,15.5545,15.5537,15.5471,15.5462,15.5412,15.5355,15.5329,15.532,15.5279,15.5246,15.5196,15.5187,15.5121,15.5095,15.5071,15.5046,15.5021,15.4979,15.4946,15.4887,15.4821,15.4813,15.4721,15.4687,15.4663,15.4621,15.4614,15.4525,15.4471,15.443,15.4379,15.4348,15.4333,15.4303,15.4278,15.4259,15.4226,15.422,15.4167,15.4086,15.4068,15.403,15.4021,15.3963,15.3905,15.3866,15.3804,15.3726,15.3693,15.3645,15.364,15.3577,15.3549,15.3523,15.3406,15.3395,15.3313,15.3241,15.3166,15.3152,15.31,15.3026,15.2994,15.2991,15.297,15.2967,15.2948,15.2949,15.2939,15.2917,15.291,15.2904,15.2903,15.2892,15.2835,15.282,15.2796,15.2768,15.2737,15.2722,15.2693,15.2685,15.2667,15.265,15.2629,15.2625,15.2603,15.2561,15.2541,15.2522,15.2486,15.2461,15.2453,15.2426,15.2409,15.2405,15.2367,15.2357,15.2304,15.2283,15.2275,15.2271,15.226,15.2247,15.2226,15.2178,15.2115,15.2086,15.2034,15.2002,15.1985,15.1957,15.1891,15.1824,15.1785,15.1746,15.1692,15.166,15.1621,15.1563,15.1457,15.1496,15.1526,15.1532,15.154,15.1548,15.1573,15.1649,15.1714,15.1783,15.1792,15.1797,15.1886,15.1924,15.1959,15.1987,15.2111,15.2141,15.2224,15.2254,15.2355,15.2389,15.2456,15.2486,15.2524,15.2585,15.2619,15.2662,15.2693,15.2687,15.2686,15.2688,15.2696,15.2703,15.2713,15.2726,15.2743,15.2766,15.2783,15.2795,15.281,15.2828,15.2836,15.2839,15.2844,15.2849,15.2861,15.2868,15.288,15.29,15.2917,15.2935,15.2953,15.2974,15.3001,15.3022,15.3045,15.3066,15.3083,15.31,15.3119,15.3133,15.3147,15.3158,15.3186,15.3208,15.3224,15.3237,15.325,15.3267,15.3281,15.3296,15.3314,15.3332,15.3341,15.3361,15.3378,15.3393,15.3417,15.3434,15.3461,15.3476,15.3499,15.3518,15.3545,15.356,15.3569,15.3582,15.3591,15.3591,15.3577,15.356,15.3536,15.3511,15.3482,15.3455,15.3431,15.3406,15.3368,15.3339,15.3311,15.3276,15.3258,15.3252,15.3257,15.3255,15.3245,15.3233,15.3212,15.3188,15.3171,15.3164,15.3174,15.3192,15.3213,15.3233,15.3239,15.3242,15.3235,15.3213,15.3196,15.3167,15.3136,15.3105,15.3088,15.3056,15.3029,15.2994,15.2957,15.2925,15.2878,15.2847,15.282,15.2793,15.2769,15.2758,15.2743,15.2715,15.2689,15.266,15.2636,15.2609,15.2579,15.2561,15.2537,15.2513,15.25,15.249,15.2488,15.2509,15.2521,15.2548,15.2596,15.2634,15.2662,15.2691,15.2706,15.2723,15.2725,15.2728,15.273,15.2725,15.2722,15.2717,15.2703,15.2696,15.269,15.2672,15.2646,15.2618,15.2596,15.257,15.2536,15.2498,15.2467,15.2425,15.2392,15.2357,15.2313,15.2288,15.2247,15.2232,15.2201,15.2176,15.2149,15.21,15.2049,15.2006,15.1971,15.1941,15.1907,15.1883,15.1855,15.1831,15.1807,15.1771,15.1741,15.1712,15.1692,15.1677,15.1658,15.163,15.1607,15.1563,15.1537,15.15,15.1443,15.1382,15.1345,15.131,15.1278,15.1247,15.1219,15.1194,15.118,15.1162,15.1147,15.1113,15.1079,15.1068,15.1046,15.1014,15.1003,15.0971,15.0954,15.0934,15.0912,15.0863,15.0805,15.0766,15.0711,15.0663,15.0617,15.0561,15.0491,15.0454,15.0417,15.037,15.0337,15.0306,15.0273,15.024,15.0219,15.0182,15.0161,15.0118,15.0115,15.0092,15.0078,15.0058,15.0015,15.0006,14.9953,14.9934,14.9908,14.9887,14.9865,14.9844,14.9827,14.9803,14.9784,14.9767,14.9752,14.9732,14.9712,14.9695,14.9682,14.9662,14.9648,14.9629,14.963,14.9633,14.9632,14.9636,14.9635,14.9626,14.9628,14.9631,14.9631,14.9626,14.9621,14.9617,14.9612,14.9608,14.9609,14.9603,14.9597,14.9585,14.9579,14.9574,14.957,14.9564,14.9564,14.9565,14.9565,14.9566,14.9573,14.9575,14.958,14.9581,14.9585,14.9588,14.9597,14.9602,14.9608,14.9607,14.9611,14.961,14.961,14.9606,14.9605,14.9604,14.9598,14.9599,14.9597,14.9594,14.9593,14.9595,14.959,14.9587,14.959,14.959,14.9589,14.9586,14.9577,14.9568,14.9567,14.9563,14.9559,14.9563,14.9562,14.9555,14.9557,14.9554,14.9553,14.9553,14.9557,14.9564,14.9583,14.9599,14.9618,14.9634,14.9653,14.967,14.9686,14.9701,14.9708,14.9716,14.9721,14.9723,14.9726,14.9721,14.9722,14.9722,14.9722,14.9726,14.9731,14.9736,14.9741,14.9746,14.9749,14.9756,14.9761,14.9766,14.9768,14.9771,14.9776,14.9785,14.9788,14.9786,14.9788,14.9792,14.9802,14.9816,14.9826,14.9835,14.9849,14.9867,14.9887,14.9908,14.9934,14.9964,15.0027,15.0064,15.0109,15.0165,15.02,15.025,15.0334,15.0378,15.0417,15.0455,15.0494,15.0512,15.0511,15.0485,15.0441,15.0389,15.0327,15.0287,15.0234,15.0183,15.0123,15.0071,15.0033,15.0009,14.9993,14.9979,14.9801,14.9688,14.968,14.9661,14.9635,14.9619,14.9611,14.9601,14.9588,14.9568,14.9547,14.9521,14.9497,14.9479,14.9459,14.9438,14.9412,14.9391,14.9358,14.933,14.9306,14.9282,14.9249,14.9222,14.9179,14.915,14.9119,14.9085,14.9042,14.9013,14.8984,14.8977,14.8982,14.8987,14.9,14.9023,14.9062,14.9086,14.9132,14.9208,14.9314,14.9416,14.9459,14.9455,14.9421,14.9305,14.9206,14.9151,14.9102,14.903,14.8965,14.8884,14.8794,14.8794,14.8832,14.8859,14.8815,14.8721,14.8622,14.8578,14.8546,14.8572,14.8597,14.8615,14.8626,14.8635,14.8649,14.866,14.8673,14.8695,14.8697,14.8691,14.8692,14.869,14.869,14.8688,14.8687,14.8685,14.8684,14.8688,14.8686,14.8676,14.8674,14.8674,14.8673,14.8671,14.8669,14.8669,14.8659,14.8656,14.865,14.8645,14.8641,14.8634,14.8632,14.8621,14.8611,14.8597,14.8591,14.8577,14.856,14.8551,14.8531,14.8505,14.8482,14.8469,14.8448,14.8441,14.843,14.8417,14.8395,14.8378,14.8352,14.8335,14.8316,14.8292,14.8273,14.825,14.8237,14.8237,14.8245,14.8253,14.8259,14.8272,14.8287,14.8304,14.8314,14.8322,14.8333,14.8339,14.8349,14.8358,14.8366,14.837,14.8377,14.8386,14.8389,14.8393,14.8395,14.8395,14.8398,14.8401,14.8401,14.8395,14.839,14.8384,14.8377,14.8371,14.8363,14.8358,14.8355,14.8346,14.8334,14.8322,14.8308,14.8293,14.8281,14.8259,14.8254,14.8249,14.8241,14.823,14.8222,14.8208,14.8177,14.8155,14.8127,14.8101,14.8078,14.8051,14.8003,14.7966,14.7929,14.7886,14.7842,14.781,14.7766,14.7726,14.7667,14.763,14.7582,14.7532,14.7476,14.7445,14.728,14.7115,14.695,14.6786,14.6621,14.6569,14.6516,14.6464,14.6412,14.6359,14.6307,14.6255,14.6218,14.6212,14.6202,14.615,14.6097,14.6348,14.6457,14.6516,14.6684,14.6851,14.7019,14.6785,14.6551,14.6551,14.6413,14.6233,14.6206,14.6183,14.6181,14.6179,14.6162,14.6151,14.614,14.612,14.6057,14.6057,14.6005,14.6005,14.5992,14.5986,14.598,14.5926,14.5877,14.5829,14.5829,14.583,14.5831,14.5827,14.5831,14.583,14.5827,14.5827,14.5818,14.5817,14.5807,14.5804,14.5793,14.5763,14.5761,14.5761,14.5758,14.5758,14.5755,14.5754,14.5745,14.5735,14.5738,14.5747,14.5746,14.5745,14.5742,14.5717,14.568,14.5678,14.5667,14.5638,14.5627,14.561,14.561,14.5604,14.562,14.5623,14.5628,14.5632,14.5627,14.5623,14.5621,14.5621,14.5619,14.5619,14.562,14.562,14.5626,14.563,14.5634,14.5639,14.564,14.5647,14.5651,14.5654,14.5655,14.5658,14.5665,14.5675,14.5702,14.5724,14.5743,14.577,14.5776,14.5779,14.5787,14.5787,14.579,14.5807,14.5802,14.5791,14.5799,14.58,14.5801,14.5812,14.5823,14.5823,14.5838,14.5838,14.5864,14.5864,14.5864,14.5864,14.586,14.5845,14.5758,14.5747,14.573,14.5709,14.5588,14.5588,14.5578,14.5578,14.5412,14.5404,14.5395,14.5385,14.5385,14.5335,14.531,14.531,14.5304,14.5145,14.4984,14.4973,14.493,14.4911,14.4864,14.4859,14.484,14.4741,14.4696,14.4617,14.4493,14.437,14.4246,14.4123,14.4117,14.4116,14.411,14.411,14.4109,14.4109,14.4109,14.4103,14.4096,14.4094,14.4089,14.4083,14.4076,14.4069,14.4062,14.4055,14.4052,14.4048,14.4163,14.4278,14.4393,14.4507,14.4622,14.4737,14.4852,14.4967,14.5081,14.5196,14.5311,14.5425,14.554,14.5655,14.577,14.5884,14.5999,14.6114,14.6069,14.6064,14.6025,14.5981,14.5936,14.5892,14.5849,14.5847,14.5803,14.5764,14.5726,14.5654,14.5592,14.5583,14.562,14.5657,14.5694,14.573,14.5767,14.5804,14.5915,14.6025,14.6135,14.6362,14.6589,14.6816,14.7043,14.7087,14.7329,14.7616,14.7902,14.818,14.8188,14.8474,14.8761,14.8798,14.9052,14.9057,14.9071,14.9079,14.9091,14.91,14.911,14.912,14.9129,14.9139,14.9151,14.9169,14.9179,14.9194,14.9206,14.922,14.9229,14.9245,14.9263,14.9279,14.9294,14.9308,14.9321,14.9333,14.9341,14.9356,14.9363,14.9375,14.9386,14.9401,14.9412,14.943,14.944,14.9454,14.9466,14.9478,14.9492,14.9502,14.9512,14.9521,14.953,14.9552,14.956,14.9568,14.958,14.9593,14.9608,14.9623,14.9628,14.9632,14.9634,14.964,14.9647,14.965,14.9656,14.9666,14.9675,14.9682,14.969,14.9703,14.9719,14.9724,14.9726,14.973,14.9742,14.9751,14.9759,14.9766,14.9775,14.9783,14.9794,14.98,14.9808,14.9813,14.9822,14.9824,14.9829,14.9837,14.985,14.9861,14.9872,14.9882,14.9905,14.9916,14.9933,14.9955,14.9973,15.0004,15.0273,15.0543,15.087,15.0857,15.0877,15.0907,15.0943,15.0974,15.1015,15.1057,15.109,15.1129,15.1163,15.119,15.1224,15.1249,15.1274,15.1302,15.1321,15.1346,15.1365,15.1393,15.141,15.1479,15.1485,15.1493,15.1568,15.1635,15.1704,15.1774,15.184,15.191,15.1979,15.2035,15.2082,15.2102,15.2132,15.2154,15.2165,15.2154,15.2104,15.2082,15.2068,15.2052,15.2057,15.2077,15.211,15.2149,15.2204,15.2224,15.2265,15.2271,15.2304,15.2332,15.2365,15.2404,15.2438,15.2471,15.251,15.2565,15.2613,15.2657,15.2685,15.2679,15.2668,15.2649,15.2635,15.2624,15.2632,15.2663,15.2699,15.2754,15.281,15.2865,15.2918,15.2982,15.3029,15.3077,15.311,15.3149,15.3182,15.3224,15.3257,15.3282,15.331,15.3343,15.3368,15.3388,15.3415,15.3435,15.346,15.3485,15.3513,15.3529,15.3549,15.3568,15.3602,15.3635,15.3677,15.3724,15.3779,15.3832,15.3888,15.3957,15.4027,15.4093,15.4157,15.4218,15.4274,15.4324,15.4365,15.4413,15.4432,15.4452,15.4471,15.4527,15.4579,15.4621,15.4663,15.4649,15.469,15.4743,15.4804,15.4874,15.4938,15.499,15.4996,15.5032,15.5057,15.5057,15.504,15.5018,15.4996,15.4988,15.4971,15.4967,15.4949,15.4927,15.4918,15.491,15.4899,15.4877,15.4857,15.4827,15.4804,15.479,15.4788,15.4799,15.4832,15.4857,15.489,15.491,15.4938,15.4963,15.4996,15.4997,15.5004,15.5043,15.5082,15.5193,15.5193,15.5254,15.531,15.5354,15.5396,15.5429,15.5474,15.5529,15.559,15.5652,15.5707,15.5768,15.5829,15.5857,15.5877,15.5902,15.5929,15.5946,15.5974,15.5999,15.6027,15.606,15.6079,15.6113,15.6129,15.6149,15.616,15.6179,15.6199,15.6232,15.626,15.6279,15.629,15.6282,15.6274,15.6257,15.6235,15.6218,15.6196,15.6193,15.6213,15.6246,15.6293,15.6335,15.6382,15.6424,15.6471,15.6524,15.6579,15.6649,15.6704,15.6752,15.6802,15.6838,15.6852,15.6849,15.6863,15.6904,15.6954,15.7015,15.7057,15.709,15.7124,15.7157,15.7196,15.7229,15.7257,15.729,15.7324,15.7354,15.7382,15.7407,15.744,15.7474,15.7515,15.7549,15.7563,15.7602,15.7649,15.7704,15.7765,15.7821,15.7871,15.7893,15.7921,15.7957,15.8013,15.8046,15.8079,15.811,15.8143,15.8185,15.8218,15.8265,15.8313,15.8368,15.8421,15.8477,15.8538,15.8607,15.8677,15.8738,15.8793,15.884,15.8865,15.8893,15.8913,15.8918,15.8946,15.8985,15.9049,15.911,15.9165,15.9227,15.929,15.9329,15.9371,15.9404,15.9429,15.9463,15.9482,15.9507,15.9535,15.9568,15.9607,15.9654,15.9696,15.9743,15.9782,15.9832,15.9879,15.994,15.9993,15.9996,15.9997,16.0049,16.011,16.0174,16.0235,16.0296,16.0352,16.0404,16.046,16.0513,16.0549,16.0588,16.0615,16.064,16.066,16.0685,16.0699,16.0715,16.0743,16.0763,16.0779,16.0813,16.084,16.0879,16.0913,16.0946,16.0993,16.1035,16.1074,16.104,16.1004,16.096,16.0927,16.0882,16.0849,16.0804,16.0777,16.0781,16.0782,16.0807,16.0854,16.0896,16.0943,16.0996,16.106,16.1115,16.1149,16.1179,16.1207,16.1223,16.1238,16.126,16.1296,16.1338,16.1385,16.1449,16.1488,16.1492,16.1493,16.148,16.1477,16.1449,16.1413,16.1371,16.134,16.1318,16.1318,16.131,16.1313,16.1318,16.1324,16.1335,16.134,16.1352,16.1357,16.1354,16.1346,16.1324,16.1307,16.1285,16.1254,16.1235,16.1204,16.1188,16.1165,16.1163,16.1154,16.1138,16.111,16.1074,16.104,16.1004,16.0982,16.0971,16.0985,16.1018,16.1052,16.109,16.1124,16.1157,16.119,16.1215,16.1243,16.1277,16.1315,16.1385,16.144,16.149,16.1532,16.1574,16.1615,16.1657,16.1667,16.1713,16.1768,16.1815,16.1843,16.1846,16.1824,16.1802,16.1774,16.1738,16.1707,16.1679,16.1643,16.1615,16.1574,16.1543,16.1515,16.1493,16.1471,16.1454,16.1452,16.1443,16.1435,16.1418,16.1402,16.1385,16.1363,16.1335,16.1313,16.1282,16.1254,16.1227,16.119,16.116,16.1138,16.1179,16.1227,16.1274,16.1321,16.1346,16.1388,16.1435,16.1496,16.1543,16.1599,16.166,16.1729,16.1799,16.1868,16.1929,16.1977,16.1996,16.2002,16.2018,16.2032,16.2029,16.2049,16.2088,16.2121,16.2154,16.2196,16.2221,16.2246,16.2279,16.2307,16.2327,16.2324,16.2307,16.2302,16.229,16.2302,16.2335,16.2368,16.2424,16.2479,16.2529,16.2585,16.2632,16.2674,16.2729,16.2777,16.2804,16.2829,16.2863,16.2896,16.2952,16.3013,16.3063,16.3096,16.3113,16.3135,16.3163,16.3207,16.3254,16.3288,16.3313,16.3346,16.3379,16.3413,16.3446,16.3479,16.3513,16.3538,16.354,16.3549,16.3554,16.3552,16.3557,16.3574,16.3607,16.3649,16.3682,16.3729,16.3782,16.3824,16.3857,16.3882,16.391,16.3935,16.3963,16.3993,16.4027,16.406,16.4093,16.4121,16.4146,16.4171,16.4204,16.4229,16.4263,16.4296,16.4352,16.439,16.4432,16.4465,16.449,16.4518,16.4549,16.4582,16.4624,16.4671,16.4718,16.4757,16.4785,16.4802,16.4829,16.4854,16.4888,16.4935,16.4968,16.4996,16.4997,16.501,16.504,16.5082,16.5115,16.514,16.5174,16.5215,16.5254,16.5293,16.5335,16.5382,16.5415,16.544,16.5474,16.5493,16.5504,16.551,16.551,16.5493,16.5465,16.5443,16.5407,16.5365,16.5329,16.534,16.5374,16.5415,16.5468,16.5532,16.5588,16.5649,16.571,16.5765,16.5827,16.5846,16.5838,16.5821,16.5799,16.5768,16.5746,16.5724,16.571,16.5721,16.5738,16.5765,16.5804,16.5865,16.5915,16.5971,16.6027,16.6096,16.6157,16.6218,16.6282,16.6321,16.6354,16.6382,16.6393,16.6404,16.6404,16.6402,16.6399,16.6396,16.6388,16.6385,16.6404,16.6424,16.6429,16.6435,16.6454,16.6482,16.6499,16.6532,16.6557,16.6571,16.6563,16.6538,16.6504,16.646,16.6418,16.6371,16.6335,16.6313,16.631,16.6299,16.6296,16.6288,16.6271,16.6263,16.6262,16.626,16.6257,16.6263,16.6274,16.6285,16.6296,16.6315,16.6349,16.6374,16.6402,16.639,16.6363,16.6327,16.6296,16.6288,16.6296,16.6288,16.6285,16.6282,16.6279,16.6279,16.6274,16.6271,16.6268,16.6293,16.6349,16.6402,16.6435,16.6432,16.6424,16.6427,16.6432,16.6424,16.6407,16.639,16.6374,16.6374,16.6377,16.6388,16.6393,16.6371,16.634,16.6299,16.6263,16.6246,16.6265,16.6313,16.636,16.6402,16.644,16.6474,16.6485,16.6468,16.6449,16.6424,16.6396,16.6374,16.6385,16.6404,16.6438,16.6468,16.6496,16.6507,16.649,16.6463,16.6432,16.641,16.6393,16.6385,16.6368,16.6379,16.6389,16.6407,16.644,16.6471,16.6504,16.6529,16.6557,16.6527,16.6485,16.6443,16.6407,16.6371,16.6349,16.636,16.6407,16.6454,16.646,16.6435,16.6427,16.6424,16.6418,16.6424]}]],[[{"lng":[-17.4742,-17.4734,-17.4726,-17.4713,-17.4713,-17.4734,-17.4736,-17.475,-17.4755,-17.4751,-17.4746,-17.4763,-17.4759,-17.4742],"lat":[14.6583,14.6575,14.6575,14.6562,14.6554,14.6534,14.6533,14.6533,14.6546,14.6553,14.6562,14.6571,14.6584,14.6583]},{"lng":[-17.4001,-17.398,-17.398,-17.4009,-17.4025,-17.4038,-17.403,-17.4038,-17.403,-17.4026,-17.4001],"lat":[14.6634,14.6654,14.668,14.6717,14.6725,14.6712,14.6696,14.6688,14.6679,14.6633,14.6634]},{"lng":[-17.5443,-17.5438,-17.5438,-17.5442,-17.545,-17.5455,-17.5455,-17.5451,-17.5443],"lat":[14.7417,14.7421,14.7446,14.745,14.745,14.7446,14.7421,14.7417,14.7417]},{"lng":[-17.5392,-17.5371,-17.5376,-17.5384,-17.5405,-17.5401,-17.5392],"lat":[14.7475,14.7496,14.75,14.75,14.7479,14.7475,14.7475]},{"lng":[-17.5134,-17.5142,-17.5159,-17.5167,-17.5184,-17.5197,-17.5196,-17.5182,-17.5159,-17.5142,-17.5129,-17.513,-17.5134],"lat":[14.7584,14.7591,14.7592,14.7583,14.7584,14.7571,14.7562,14.7561,14.7559,14.755,14.7562,14.758,14.7584]},{"lng":[-17.4775,-17.4771,-17.4771,-17.4771,-17.4776,-17.4792,-17.4796,-17.4796,-17.4797,-17.4792,-17.4775],"lat":[14.7692,14.7696,14.7711,14.7712,14.7716,14.7717,14.7713,14.7705,14.7695,14.7692,14.7692]},{"lng":[-16.7366,-16.7381,-16.7405,-16.7415,-16.7432,-16.7423,-16.7401,-16.7401,-16.738,-16.7373,-16.7338,-16.7337,-16.7285,-16.7277,-16.7274,-16.7253,-16.7247,-16.7231,-16.7209,-16.7187,-16.7155,-16.7127,-16.7075,-16.7033,-16.7031,-16.6995,-16.6973,-16.6908,-16.6889,-16.6887,-16.6883,-16.6888,-16.6886,-16.6892,-16.6901,-16.6911,-16.6958,-16.7012,-16.7102,-16.7307,-16.7313,-16.7347,-16.7363,-16.7405,-16.7421,-16.7455,-16.7455,-16.748,-16.7513,-16.7546,-16.7546,-16.7596,-16.7613,-16.763,-16.7638,-16.7655,-16.7663,-16.768,-16.7713,-16.773,-16.7746,-16.7746,-16.7763,-16.7763,-16.7805,-16.7821,-16.7838,-16.7863,-16.7872,-16.7913,-16.7988,-16.7988,-16.8055,-16.8063,-16.8088,-16.8155,-16.8188,-16.8213,-16.8238,-16.8263,-16.8296,-16.8321,-16.8338,-16.8363,-16.8363,-16.8388,-16.8397,-16.8413,-16.8446,-16.8463,-16.8463,-16.848,-16.8505,-16.8521,-16.8572,-16.8588,-16.8663,-16.8663,-16.8688,-16.8688,-16.873,-16.8738,-16.8763,-16.8771,-16.8838,-16.8838,-16.8871,-16.888,-16.893,-16.893,-16.9055,-16.9055,-16.9096,-16.9105,-16.913,-16.9196,-16.9205,-16.9238,-16.9305,-16.9313,-16.9388,-16.9405,-16.9496,-16.9497,-16.958,-16.9588,-16.9588,-16.9704,-16.9755,-16.9755,-16.9788,-16.9788,-16.988,-16.988,-16.9938,-16.9946,-16.9992,-17.0063,-17.0063,-17.0146,-17.0163,-17.0255,-17.0255,-17.0363,-17.0363,-17.0413,-17.0413,-17.0472,-17.0472,-17.0596,-17.0597,-17.0646,-17.0646,-17.068,-17.0738,-17.0738,-17.0763,-17.0772,-17.0909,-17.0917,-17.0946,-17.098,-17.1063,-17.108,-17.1105,-17.1138,-17.1155,-17.1163,-17.118,-17.1196,-17.1213,-17.123,-17.1238,-17.1255,-17.1255,-17.1288,-17.1301,-17.1309,-17.1325,-17.1375,-17.1417,-17.1517,-17.1551,-17.1634,-17.1684,-17.1785,-17.1817,-17.1951,-17.1975,-17.2017,-17.2084,-17.2125,-17.2217,-17.2275,-17.2334,-17.235,-17.2417,-17.2434,-17.2459,-17.2476,-17.2509,-17.2526,-17.2551,-17.2567,-17.2592,-17.2609,-17.2742,-17.2759,-17.2859,-17.2875,-17.2951,-17.2967,-17.2992,-17.3009,-17.3018,-17.3067,-17.3092,-17.3109,-17.3134,-17.3251,-17.3268,-17.3284,-17.3309,-17.3326,-17.3376,-17.3392,-17.3509,-17.3525,-17.3568,-17.3592,-17.3626,-17.365,-17.3676,-17.3692,-17.3742,-17.3767,-17.38,-17.3826,-17.3842,-17.3867,-17.39,-17.3951,-17.3967,-17.3992,-17.4034,-17.4042,-17.4067,-17.4076,-17.4109,-17.4117,-17.4142,-17.415,-17.4168,-17.4176,-17.4218,-17.4226,-17.4242,-17.4293,-17.4309,-17.4333,-17.4376,-17.4425,-17.4434,-17.4459,-17.4492,-17.4542,-17.455,-17.4584,-17.4608,-17.4626,-17.4627,-17.4651,-17.4684,-17.4692,-17.4718,-17.4725,-17.4759,-17.4792,-17.4809,-17.4842,-17.4896,-17.49,-17.4934,-17.4951,-17.4967,-17.4992,-17.5009,-17.5017,-17.5051,-17.5068,-17.5101,-17.5142,-17.5184,-17.5196,-17.5209,-17.5234,-17.5243,-17.5259,-17.5263,-17.5292,-17.5301,-17.5313,-17.5313,-17.5313,-17.533,-17.5225,-17.5217,-17.5201,-17.5184,-17.508,-17.508,-17.5071,-17.508,-17.5059,-17.5025,-17.5,-17.4976,-17.4951,-17.493,-17.493,-17.49,-17.488,-17.4888,-17.4875,-17.4859,-17.483,-17.4829,-17.4817,-17.48,-17.478,-17.478,-17.4805,-17.4805,-17.478,-17.4788,-17.4784,-17.4767,-17.4743,-17.473,-17.4735,-17.4738,-17.473,-17.473,-17.4747,-17.4746,-17.4726,-17.4717,-17.4701,-17.4688,-17.4684,-17.4658,-17.4646,-17.4642,-17.4634,-17.4601,-17.4584,-17.4542,-17.4534,-17.4517,-17.4488,-17.4488,-17.4463,-17.4463,-17.4463,-17.4442,-17.4396,-17.4396,-17.4386,-17.438,-17.4358,-17.4334,-17.4321,-17.4315,-17.4313,-17.4313,-17.4313,-17.4338,-17.4338,-17.4318,-17.4267,-17.4251,-17.4221,-17.4205,-17.4205,-17.4226,-17.425,-17.4259,-17.4268,-17.428,-17.4284,-17.43,-17.4317,-17.4325,-17.433,-17.4322,-17.4338,-17.433,-17.4355,-17.4355,-17.4342,-17.4334,-17.4317,-17.4301,-17.4293,-17.4267,-17.4251,-17.4238,-17.4246,-17.4234,-17.4225,-17.4205,-17.4205,-17.4221,-17.4221,-17.4234,-17.425,-17.4271,-17.4271,-17.4263,-17.4242,-17.4217,-17.418,-17.4179,-17.4192,-17.4225,-17.425,-17.4276,-17.4321,-17.433,-17.433,-17.4321,-17.4313,-17.4313,-17.4296,-17.4217,-17.4209,-17.4192,-17.4159,-17.4101,-17.4092,-17.4059,-17.405,-17.4017,-17.4009,-17.3967,-17.3959,-17.3867,-17.3859,-17.3684,-17.3675,-17.3625,-17.3617,-17.3584,-17.3576,-17.3542,-17.3534,-17.3475,-17.3392,-17.3258,-17.3242,-17.3125,-17.3075,-17.3059,-17.3017,-17.3,-17.2992,-17.2975,-17.2942,-17.2901,-17.2859,-17.2826,-17.2784,-17.2759,-17.2742,-17.2726,-17.2692,-17.2676,-17.2659,-17.2642,-17.2617,-17.2601,-17.2525,-17.2484,-17.2467,-17.2351,-17.2284,-17.2259,-17.2217,-17.2208,-17.2192,-17.2084,-17.2076,-17.205,-17.2026,-17.1971,-17.1971,-17.193,-17.193,-17.1867,-17.1859,-17.1755,-17.1755,-17.1696,-17.1663,-17.1663,-17.1613,-17.1605,-17.1588,-17.1571,-17.1546,-17.1521,-17.1521,-17.1505,-17.1496,-17.1496,-17.1471,-17.1459,-17.1446,-17.1446,-17.143,-17.143,-17.1413,-17.1372,-17.1346,-17.1334,-17.1267,-17.1242,-17.1159,-17.1113,-17.1113,-17.1088,-17.1088,-17.1063,-17.1063,-17.1038,-17.1013,-17.1005,-17.098,-17.0971,-17.0946,-17.0946,-17.0938,-17.093,-17.093,-17.0921,-17.0922,-17.0896,-17.0896,-17.0896,-17.088,-17.0805,-17.0805,-17.0747,-17.0755,-17.0726,-17.0712,-17.0709,-17.07,-17.0684,-17.0638,-17.0638,-17.0613,-17.0605,-17.0596,-17.0596,-17.0596,-17.0588,-17.0588,-17.0559,-17.0523,-17.0509,-17.0459,-17.0451,-17.0434,-17.0392,-17.0376,-17.0351,-17.0325,-17.0309,-17.025,-17.0217,-17.0209,-17.0167,-17.0151,-17.0142,-17.0117,-17.0059,-16.9992,-16.9972,-16.9971,-16.9867,-16.9834,-16.9805,-16.9751,-16.9734,-16.9709,-16.9692,-16.9655,-16.9613,-16.9596,-16.9551,-16.9525,-16.9513,-16.9505,-16.9496,-16.9471,-16.9471,-16.9463,-16.9463,-16.9396,-16.9355,-16.9355,-16.9338,-16.933,-16.9313,-16.9313,-16.9321,-16.933,-16.9338,-16.9338,-16.9346,-16.9346,-16.9355,-16.9355,-16.9346,-16.9355,-16.9346,-16.9338,-16.9322,-16.9321,-16.9313,-16.9309,-16.9308,-16.9302,-16.9275,-16.9259,-16.9226,-16.9201,-16.9175,-16.9105,-16.9071,-16.8988,-16.8942,-16.8901,-16.8863,-16.8855,-16.8821,-16.878,-16.8771,-16.875,-16.8738,-16.874,-16.8744,-16.8746,-16.8755,-16.875,-16.8738,-16.873,-16.873,-16.8721,-16.8721,-16.8713,-16.8713,-16.8705,-16.8705,-16.8713,-16.8713,-16.8671,-16.8671,-16.8678,-16.8688,-16.8688,-16.8675,-16.8662,-16.8659,-16.8592,-16.853,-16.8521,-16.8513,-16.8513,-16.8484,-16.8475,-16.8384,-16.8325,-16.8301,-16.8284,-16.8242,-16.8226,-16.8213,-16.8225,-16.8242,-16.8255,-16.8255,-16.8242,-16.8226,-16.8201,-16.8159,-16.815,-16.8134,-16.813,-16.8151,-16.8176,-16.8196,-16.8197,-16.8167,-16.8084,-16.8055,-16.8055,-16.8038,-16.7998,-16.7849,-16.7835,-16.7818,-16.7798,-16.778,-16.7766,-16.7763,-16.7737,-16.7725,-16.771,-16.7693,-16.7677,-16.7655,-16.7641,-16.7633,-16.7625,-16.7615,-16.7593,-16.7585,-16.7573,-16.7559,-16.7548,-16.7537,-16.7533,-16.7526,-16.7527,-16.7533,-16.7535,-16.7551,-16.7559,-16.7565,-16.7574,-16.7597,-16.7617,-16.763,-16.7634,-16.7635,-16.7624,-16.7612,-16.7601,-16.7588,-16.7574,-16.7557,-16.7539,-16.7523,-16.7507,-16.7493,-16.7479,-16.7468,-16.7456,-16.745,-16.7437,-16.7428,-16.7413,-16.7402,-16.7392,-16.7375,-16.7356,-16.7341,-16.7329,-16.7321,-16.7304,-16.7291,-16.7278,-16.7265,-16.7256,-16.7248,-16.7241,-16.7235,-16.7223,-16.7222,-16.7228,-16.7232,-16.7239,-16.7254,-16.7281,-16.7303,-16.7312,-16.7322,-16.7348,-16.7368,-16.7387,-16.7401,-16.7422,-16.7423,-16.7431,-16.7433,-16.7441,-16.7457,-16.746,-16.7476,-16.7474,-16.7486,-16.7488,-16.7488,-16.7484,-16.7475,-16.7472,-16.7463,-16.7452,-16.7437,-16.7425,-16.7409,-16.739,-16.7372,-16.7362,-16.7344,-16.7333,-16.7318,-16.7306,-16.729,-16.7277,-16.7264,-16.7246,-16.7232,-16.7219,-16.7208,-16.7198,-16.7187,-16.7178,-16.7166,-16.715,-16.7138,-16.7124,-16.7109,-16.7096,-16.7088,-16.7085,-16.7078,-16.7073,-16.7071,-16.7075,-16.7074,-16.7079,-16.7078,-16.7079,-16.7089,-16.7089,-16.7092,-16.7097,-16.7102,-16.7103,-16.7098,-16.71,-16.7102,-16.7105,-16.7104,-16.7091,-16.7088,-16.7079,-16.7079,-16.7083,-16.7079,-16.7077,-16.7072,-16.7068,-16.7065,-16.7069,-16.7064,-16.7061,-16.7053,-16.7054,-16.7052,-16.7038,-16.702,-16.7002,-16.6981,-16.6963,-16.6948,-16.6928,-16.6906,-16.6888,-16.6869,-16.685,-16.6832,-16.6811,-16.6795,-16.6777,-16.6755,-16.6733,-16.671,-16.6684,-16.6659,-16.6625,-16.6596,-16.6568,-16.6541,-16.6516,-16.6492,-16.6466,-16.6457,-16.6439,-16.6424,-16.6404,-16.6389,-16.6377,-16.6352,-16.6336,-16.6327,-16.6308,-16.6297,-16.628,-16.6263,-16.6259,-16.625,-16.6234,-16.6217,-16.6202,-16.6175,-16.6155,-16.614,-16.6119,-16.6106,-16.6089,-16.6073,-16.6057,-16.6037,-16.6017,-16.5997,-16.5972,-16.595,-16.592,-16.5906,-16.589,-16.5867,-16.5841,-16.5814,-16.5794,-16.5771,-16.5753,-16.5735,-16.5711,-16.5694,-16.5681,-16.5653,-16.5635,-16.5611,-16.559,-16.5574,-16.5545,-16.5526,-16.5512,-16.5506,-16.5495,-16.5493,-16.5484,-16.5474,-16.5459,-16.5443,-16.5435,-16.5426,-16.5412,-16.5405,-16.5394,-16.5387,-16.5382,-16.537,-16.5361,-16.5361,-16.5358,-16.5348,-16.5349,-16.5346,-16.5335,-16.533,-16.5322,-16.5317,-16.5314,-16.531,-16.5306,-16.5308,-16.5303,-16.5299,-16.5294,-16.5287,-16.5286,-16.5285,-16.5279,-16.5268,-16.5269,-16.5271,-16.5273,-16.5273,-16.5267,-16.5257,-16.5246,-16.5234,-16.5225,-16.5217,-16.5212,-16.5206,-16.5198,-16.519,-16.5186,-16.518,-16.5177,-16.5178,-16.5176,-16.5176,-16.5177,-16.5178,-16.518,-16.5181,-16.5184,-16.5189,-16.5193,-16.5199,-16.5207,-16.5215,-16.5223,-16.5229,-16.5234,-16.5241,-16.5241,-16.5259,-16.5278,-16.5285,-16.5296,-16.5303,-16.5311,-16.5322,-16.5339,-16.5348,-16.5362,-16.5374,-16.5386,-16.5403,-16.5421,-16.5436,-16.5452,-16.5467,-16.5481,-16.55,-16.5516,-16.5533,-16.5546,-16.5562,-16.5568,-16.5591,-16.5608,-16.562,-16.563,-16.5647,-16.5654,-16.5667,-16.5682,-16.5694,-16.5701,-16.5727,-16.5737,-16.575,-16.5768,-16.5776,-16.5791,-16.5802,-16.5812,-16.5819,-16.5829,-16.5845,-16.5849,-16.5862,-16.587,-16.5883,-16.5895,-16.5905,-16.5914,-16.5926,-16.594,-16.5948,-16.5958,-16.5969,-16.5982,-16.5994,-16.6011,-16.6022,-16.6043,-16.6056,-16.6073,-16.609,-16.6102,-16.6118,-16.6134,-16.615,-16.6163,-16.6171,-16.6192,-16.6206,-16.6209,-16.6214,-16.6211,-16.6212,-16.6208,-16.6206,-16.6201,-16.6196,-16.6185,-16.6169,-16.6156,-16.6142,-16.6123,-16.6108,-16.6104,-16.6088,-16.6075,-16.6064,-16.6056,-16.6053,-16.6047,-16.6045,-16.6047,-16.6046,-16.6044,-16.6047,-16.6047,-16.6051,-16.6055,-16.6051,-16.605,-16.6054,-16.6055,-16.6056,-16.6053,-16.6053,-16.605,-16.605,-16.6044,-16.6042,-16.6039,-16.6039,-16.6048,-16.6051,-16.6053,-16.6054,-16.6055,-16.6058,-16.6058,-16.6062,-16.6064,-16.6068,-16.607,-16.6075,-16.608,-16.6093,-16.6109,-16.613,-16.6139,-16.6158,-16.6222,-16.6244,-16.6264,-16.6285,-16.6306,-16.6331,-16.6356,-16.6383,-16.6408,-16.6429,-16.6451,-16.6484,-16.6486,-16.6474,-16.6468,-16.6453,-16.6445,-16.6437,-16.6428,-16.6411,-16.6395,-16.6374,-16.6353,-16.6339,-16.632,-16.6293,-16.6276,-16.6247,-16.6237,-16.6233,-16.6242,-16.6242,-16.6248,-16.6248,-16.6256,-16.6261,-16.6267,-16.6271,-16.6285,-16.6293,-16.6304,-16.631,-16.6319,-16.6331,-16.6346,-16.636,-16.6369,-16.6384,-16.6401,-16.6419,-16.6433,-16.6451,-16.6468,-16.6495,-16.6516,-16.6524,-16.6534,-16.6545,-16.6552,-16.656,-16.6568,-16.6578,-16.6583,-16.6594,-16.6608,-16.6615,-16.6633,-16.6648,-16.6652,-16.6663,-16.6678,-16.6695,-16.6707,-16.6717,-16.6729,-16.6741,-16.6769,-16.6793,-16.6813,-16.6833,-16.6855,-16.6841,-16.6802,-16.6773,-16.6734,-16.6685,-16.6652,-16.6628,-16.6615,-16.658,-16.6542,-16.6522,-16.6509,-16.6479,-16.6444,-16.6402,-16.6349,-16.6288,-16.6211,-16.6143,-16.607,-16.6015,-16.5963,-16.5913,-16.5871,-16.5799,-16.5738,-16.5645,-16.5608,-16.5576,-16.5519,-16.5438,-16.5377,-16.5323,-16.527,-16.5202,-16.5157,-16.5135,-16.5102,-16.5066,-16.5029,-16.5002,-16.4981,-16.4961,-16.4965,-16.4939,-16.4921,-16.4889,-16.4858,-16.4815,-16.4789,-16.4768,-16.4755,-16.4742,-16.4732,-16.4713,-16.4677,-16.4636,-16.4591,-16.4553,-16.4521,-16.4469,-16.4454,-16.4424,-16.4405,-16.438,-16.4353,-16.4334,-16.4314,-16.4291,-16.4263,-16.4232,-16.4198,-16.4169,-16.4138,-16.4107,-16.4084,-16.406,-16.4022,-16.3988,-16.394,-16.3905,-16.3901,-16.3879,-16.3858,-16.3833,-16.3788,-16.3764,-16.3736,-16.3709,-16.3683,-16.3632,-16.3592,-16.3556,-16.3512,-16.3485,-16.3462,-16.3427,-16.3395,-16.3377,-16.3346,-16.3319,-16.3301,-16.328,-16.3254,-16.3225,-16.3198,-16.3188,-16.313,-16.3099,-16.3081,-16.3049,-16.3029,-16.3008,-16.2988,-16.2976,-16.2939,-16.2911,-16.2885,-16.2867,-16.2841,-16.281,-16.2786,-16.2769,-16.2746,-16.2724,-16.2708,-16.2683,-16.2654,-16.2626,-16.2601,-16.2579,-16.255,-16.253,-16.2502,-16.2483,-16.2472,-16.2454,-16.2426,-16.2408,-16.2396,-16.237,-16.2357,-16.2345,-16.2329,-16.232,-16.2322,-16.2316,-16.2311,-16.2308,-16.2294,-16.2275,-16.2253,-16.2237,-16.2206,-16.2183,-16.2153,-16.2123,-16.2101,-16.2068,-16.2036,-16.2014,-16.1988,-16.1966,-16.1945,-16.1932,-16.1904,-16.188,-16.1854,-16.1817,-16.1782,-16.177,-16.177,-16.1771,-16.1772,-16.1771,-16.1771,-16.1766,-16.1761,-16.1755,-16.1759,-16.1764,-16.1767,-16.1771,-16.1776,-16.1779,-16.1781,-16.1782,-16.1787,-16.1783,-16.1782,-16.1784,-16.1789,-16.1796,-16.18,-16.1808,-16.1824,-16.1832,-16.1858,-16.1873,-16.1886,-16.1891,-16.1904,-16.1921,-16.1931,-16.1944,-16.1943,-16.1942,-16.1947,-16.1947,-16.1947,-16.1947,-16.1947,-16.1951,-16.1962,-16.1968,-16.1989,-16.2012,-16.2041,-16.205,-16.2086,-16.2116,-16.2136,-16.2188,-16.2249,-16.227,-16.2304,-16.2321,-16.2345,-16.2364,-16.2371,-16.2374,-16.2377,-16.2373,-16.2372,-16.2373,-16.2376,-16.2392,-16.2416,-16.2426,-16.2438,-16.245,-16.2475,-16.2499,-16.2526,-16.2551,-16.2572,-16.2587,-16.2601,-16.2616,-16.2632,-16.2647,-16.2661,-16.2672,-16.2683,-16.2688,-16.2693,-16.2699,-16.2697,-16.2692,-16.2689,-16.2685,-16.2681,-16.268,-16.2685,-16.2692,-16.27,-16.2708,-16.2723,-16.2742,-16.276,-16.278,-16.2806,-16.2832,-16.2849,-16.2858,-16.2868,-16.2885,-16.2908,-16.2931,-16.295,-16.2969,-16.3,-16.3017,-16.3045,-16.3064,-16.308,-16.3084,-16.3088,-16.3093,-16.3102,-16.3112,-16.3142,-16.3171,-16.3202,-16.3237,-16.3274,-16.3305,-16.3332,-16.335,-16.338,-16.3415,-16.3436,-16.345,-16.3457,-16.3454,-16.3446,-16.3416,-16.3388,-16.3371,-16.3354,-16.3344,-16.3346,-16.3356,-16.3371,-16.3393,-16.3411,-16.3421,-16.3448,-16.3468,-16.3496,-16.3515,-16.3547,-16.3587,-16.363,-16.366,-16.3693,-16.3719,-16.3732,-16.3752,-16.3773,-16.3796,-16.3839,-16.3878,-16.3932,-16.3979,-16.4039,-16.4107,-16.4153,-16.4189,-16.4228,-16.4255,-16.428,-16.4311,-16.4338,-16.4362,-16.4385,-16.4414,-16.4436,-16.4458,-16.4483,-16.4509,-16.4515,-16.4526,-16.4538,-16.4546,-16.4564,-16.4594,-16.4618,-16.4641,-16.4659,-16.4691,-16.4709,-16.4736,-16.4751,-16.4759,-16.4779,-16.4791,-16.4804,-16.4817,-16.4832,-16.4845,-16.4861,-16.4873,-16.4889,-16.4901,-16.4915,-16.4941,-16.4954,-16.4974,-16.499,-16.5007,-16.502,-16.503,-16.5044,-16.506,-16.5077,-16.5096,-16.5111,-16.5124,-16.5137,-16.5146,-16.5158,-16.5162,-16.5164,-16.5156,-16.5143,-16.5135,-16.5123,-16.5123,-16.5124,-16.5129,-16.5132,-16.5135,-16.5145,-16.5158,-16.5173,-16.5186,-16.5202,-16.5225,-16.5237,-16.5251,-16.5279,-16.5307,-16.5332,-16.535,-16.5381,-16.54,-16.5418,-16.5451,-16.5471,-16.5519,-16.554,-16.5571,-16.5596,-16.5647,-16.5671,-16.5696,-16.5722,-16.581,-16.5824,-16.5836,-16.5907,-16.5958,-16.5993,-16.6012,-16.6037,-16.6045,-16.6072,-16.6092,-16.6128,-16.6141,-16.614,-16.6134,-16.6129,-16.6123,-16.6121,-16.6121,-16.6133,-16.6165,-16.6207,-16.6245,-16.6285,-16.6306,-16.6307,-16.6312,-16.6338,-16.6365,-16.6384,-16.6412,-16.6423,-16.6428,-16.6476,-16.6536,-16.6557,-16.6608,-16.6626,-16.6635,-16.6669,-16.6687,-16.6701,-16.673,-16.6745,-16.676,-16.6811,-16.6828,-16.6844,-16.6867,-16.6907,-16.6922,-16.6961,-16.701,-16.7022,-16.7048,-16.7073,-16.7081,-16.7094,-16.7136,-16.7155,-16.7167,-16.718,-16.7233,-16.7305,-16.7319,-16.7333,-16.7366],"lat":[15.2967,15.297,15.2991,15.2994,15.3026,15.31,15.3152,15.3166,15.3241,15.3313,15.3395,15.3406,15.3523,15.3549,15.3577,15.364,15.3645,15.3693,15.3726,15.3804,15.3866,15.3905,15.3963,15.4021,15.403,15.4068,15.4086,15.4167,15.422,15.4226,15.4259,15.4278,15.4303,15.4333,15.4348,15.4379,15.443,15.4471,15.4525,15.4614,15.4596,15.4554,15.4512,15.4462,15.4429,15.4387,15.4379,15.4354,15.4296,15.4254,15.4238,15.4163,15.4129,15.4112,15.4088,15.4071,15.4046,15.4029,15.3979,15.392,15.3904,15.3888,15.3871,15.3854,15.3804,15.3771,15.3737,15.3712,15.3688,15.3638,15.3504,15.3496,15.3421,15.3396,15.3371,15.3238,15.3195,15.3138,15.3104,15.3046,15.3012,15.2971,15.2938,15.2913,15.2904,15.2879,15.2854,15.2837,15.2771,15.2754,15.2738,15.2713,15.2688,15.2646,15.2579,15.2546,15.2463,15.2454,15.2429,15.2421,15.2371,15.2346,15.2321,15.2295,15.2204,15.2196,15.2163,15.2137,15.2079,15.2071,15.1954,15.1946,15.1904,15.1888,15.1837,15.1771,15.1754,15.1696,15.1621,15.1596,15.1504,15.1471,15.1371,15.1363,15.1279,15.1262,15.1254,15.1129,15.1062,15.1054,15.1021,15.1012,15.0912,15.0904,15.0846,15.0821,15.0784,15.0712,15.0696,15.0612,15.0571,15.0487,15.0479,15.037,15.0363,15.0313,15.0304,15.0246,15.0238,15.0105,15.0096,15.0046,15.0037,14.9996,14.9946,14.9937,14.9912,14.9888,14.975,14.975,14.9721,14.9654,14.9571,14.9538,14.9512,14.9454,14.9446,14.9421,14.9404,14.9379,14.9313,14.9279,14.9237,14.9204,14.9187,14.9145,14.9125,14.9125,14.9109,14.9091,14.9042,14.8992,14.8983,14.8925,14.89,14.885,14.8834,14.8767,14.875,14.8742,14.8708,14.87,14.865,14.8608,14.8583,14.8567,14.8534,14.8533,14.8516,14.8517,14.85,14.85,14.8484,14.8483,14.8467,14.8467,14.8391,14.8392,14.8342,14.8342,14.83,14.83,14.8283,14.8283,14.8275,14.8267,14.825,14.825,14.8234,14.8175,14.8159,14.8158,14.8141,14.8142,14.8117,14.8116,14.8059,14.8058,14.8033,14.8034,14.8017,14.8017,14.8,14.8,14.7975,14.795,14.7933,14.7908,14.7908,14.7892,14.7884,14.7858,14.7859,14.7842,14.7833,14.7825,14.7825,14.7817,14.7816,14.7808,14.7809,14.78,14.78,14.7792,14.7784,14.7775,14.7775,14.775,14.775,14.7725,14.7709,14.7708,14.77,14.77,14.7683,14.7675,14.7667,14.7658,14.7642,14.7642,14.7641,14.7625,14.7625,14.7633,14.7633,14.7642,14.7642,14.7658,14.7675,14.7675,14.762,14.76,14.7592,14.7567,14.7575,14.7575,14.7592,14.7592,14.7559,14.7566,14.7567,14.7517,14.7525,14.7512,14.7492,14.7492,14.7501,14.7517,14.7517,14.7516,14.7496,14.7471,14.7469,14.7437,14.7413,14.7408,14.74,14.74,14.7392,14.7287,14.7271,14.7262,14.7204,14.7183,14.7183,14.715,14.715,14.7167,14.7138,14.7112,14.7075,14.7062,14.7054,14.7042,14.7042,14.702,14.6996,14.6984,14.6983,14.6971,14.6963,14.6938,14.6921,14.6904,14.6888,14.6875,14.6866,14.6867,14.6854,14.6834,14.6821,14.6813,14.6796,14.6779,14.6763,14.6742,14.6742,14.6759,14.6754,14.6742,14.6742,14.6761,14.6767,14.6767,14.6733,14.6725,14.6725,14.6717,14.6717,14.6687,14.6671,14.6637,14.6609,14.6587,14.6575,14.6538,14.6521,14.65,14.6488,14.6467,14.6467,14.6479,14.6492,14.6496,14.6501,14.6537,14.6587,14.6629,14.665,14.6675,14.6675,14.6696,14.6737,14.6771,14.6775,14.6758,14.6775,14.6767,14.6771,14.6792,14.6792,14.6775,14.6775,14.6796,14.6813,14.6821,14.6846,14.6854,14.6871,14.6884,14.6875,14.6883,14.6875,14.6892,14.685,14.6858,14.6846,14.6813,14.68,14.68,14.6821,14.6846,14.6871,14.6888,14.69,14.69,14.6921,14.6938,14.6954,14.6967,14.7,14.7012,14.7054,14.7067,14.7067,14.7092,14.7067,14.7113,14.7129,14.7162,14.7179,14.7187,14.7204,14.7237,14.7309,14.7308,14.7325,14.7342,14.735,14.7358,14.7359,14.7367,14.7367,14.7375,14.7375,14.7384,14.7392,14.74,14.74,14.7392,14.7392,14.7383,14.7383,14.7375,14.7375,14.7367,14.7367,14.7325,14.7258,14.7259,14.72,14.7175,14.7175,14.715,14.715,14.7142,14.7142,14.7125,14.7092,14.7092,14.7117,14.7117,14.71,14.71,14.7084,14.7067,14.7067,14.705,14.705,14.7034,14.7033,14.6992,14.6984,14.6975,14.6867,14.6833,14.6809,14.6792,14.6783,14.6784,14.6675,14.6675,14.665,14.6642,14.6588,14.6579,14.6537,14.6529,14.6467,14.6467,14.6362,14.6354,14.6296,14.6246,14.6238,14.6188,14.6162,14.6146,14.6113,14.6088,14.6029,14.6004,14.5979,14.5963,14.5913,14.5871,14.5867,14.5846,14.5829,14.5804,14.5788,14.5754,14.5712,14.5654,14.5642,14.5608,14.5583,14.5542,14.5496,14.5487,14.5446,14.5429,14.5396,14.5379,14.5354,14.5304,14.5237,14.5204,14.5171,14.5137,14.5121,14.5113,14.5096,14.5071,14.5062,14.5029,14.4979,14.4969,14.4945,14.4912,14.4846,14.4813,14.4763,14.4738,14.4725,14.4725,14.4725,14.4717,14.4717,14.4679,14.4663,14.4629,14.4588,14.4579,14.4571,14.4554,14.4546,14.4521,14.4492,14.4474,14.4467,14.4467,14.4475,14.4475,14.4459,14.4442,14.4442,14.4458,14.4458,14.4392,14.4392,14.44,14.4358,14.4358,14.4367,14.4366,14.4333,14.4259,14.4237,14.4229,14.4134,14.4142,14.4096,14.4041,14.4034,14.4033,14.4025,14.3988,14.3904,14.387,14.3825,14.3825,14.3813,14.3779,14.3771,14.3729,14.3704,14.3696,14.3637,14.3571,14.3504,14.3495,14.348,14.3437,14.3413,14.3313,14.3304,14.3204,14.3196,14.3163,14.3154,14.3079,14.3071,14.3038,14.3029,14.3013,14.3004,14.2962,14.2937,14.2904,14.2896,14.2859,14.285,14.2848,14.2842,14.2859,14.2859,14.2833,14.2834,14.2762,14.2704,14.2588,14.2542,14.2517,14.2479,14.2454,14.2413,14.2321,14.2279,14.2259,14.2263,14.2273,14.2306,14.2321,14.2329,14.2358,14.2312,14.2304,14.2271,14.2263,14.2221,14.2212,14.218,14.2171,14.2087,14.2079,14.2013,14.1946,14.1879,14.1873,14.1862,14.1846,14.1842,14.1848,14.185,14.1817,14.1763,14.1713,14.1704,14.1671,14.1641,14.1641,14.1575,14.1542,14.1542,14.155,14.1592,14.1592,14.1571,14.1558,14.1558,14.1546,14.1521,14.1508,14.1509,14.1492,14.1492,14.15,14.15,14.1487,14.1475,14.1475,14.1454,14.1429,14.14,14.1367,14.1346,14.1338,14.1321,14.124,14.1397,14.14,14.1403,14.1411,14.1414,14.1403,14.1394,14.1387,14.138,14.1379,14.138,14.1385,14.1395,14.1403,14.1412,14.142,14.1424,14.142,14.1412,14.1407,14.1401,14.1403,14.141,14.1415,14.1432,14.1448,14.1458,14.1468,14.1483,14.1502,14.1518,14.1532,14.1546,14.1562,14.1576,14.1594,14.1609,14.1625,14.1635,14.1655,14.1674,14.1691,14.1698,14.1705,14.1708,14.171,14.1715,14.1729,14.1739,14.1754,14.1766,14.1789,14.1806,14.182,14.1825,14.1827,14.1827,14.1823,14.1829,14.1848,14.1864,14.1879,14.1889,14.1895,14.1915,14.1932,14.1952,14.1971,14.1988,14.2011,14.2035,14.2048,14.2064,14.208,14.2091,14.2104,14.2113,14.2111,14.2112,14.2126,14.2144,14.2156,14.2173,14.2189,14.2201,14.2209,14.2222,14.2231,14.2252,14.2273,14.2294,14.2317,14.2348,14.2376,14.2407,14.2432,14.2454,14.2478,14.2504,14.2524,14.2539,14.2552,14.2574,14.2599,14.2623,14.2638,14.2658,14.2669,14.2684,14.2699,14.2713,14.2731,14.2752,14.2772,14.2791,14.2807,14.2823,14.284,14.2861,14.2876,14.2894,14.2918,14.294,14.2963,14.2984,14.3004,14.3024,14.3052,14.3069,14.3097,14.3117,14.314,14.3166,14.3186,14.321,14.3237,14.3262,14.3287,14.3306,14.3332,14.3363,14.339,14.3419,14.3449,14.3473,14.3499,14.3523,14.3551,14.3584,14.3608,14.363,14.3647,14.368,14.3707,14.373,14.3752,14.3778,14.3805,14.3835,14.3859,14.3891,14.3918,14.395,14.398,14.4004,14.4028,14.4046,14.4067,14.4081,14.4098,14.4121,14.4142,14.416,14.418,14.4199,14.4222,14.423,14.4239,14.4251,14.4259,14.4256,14.4259,14.4257,14.4257,14.426,14.4261,14.4258,14.4259,14.4266,14.4262,14.4263,14.4266,14.4267,14.4271,14.4275,14.428,14.4287,14.43,14.4314,14.4335,14.435,14.4366,14.4383,14.4391,14.4406,14.4421,14.4437,14.4451,14.4472,14.448,14.449,14.4501,14.4506,14.4515,14.4523,14.4525,14.4529,14.4536,14.4545,14.4554,14.4561,14.4586,14.4601,14.4611,14.4624,14.4646,14.4663,14.4677,14.4698,14.4712,14.4735,14.475,14.4765,14.4776,14.4797,14.4814,14.4834,14.4855,14.4871,14.4891,14.4909,14.4921,14.494,14.4955,14.4963,14.4992,14.5001,14.5015,14.503,14.5044,14.5061,14.508,14.5096,14.5111,14.5132,14.5151,14.517,14.5184,14.5197,14.5203,14.5224,14.5233,14.5252,14.5276,14.529,14.5303,14.5316,14.5328,14.535,14.5373,14.5386,14.541,14.5432,14.5455,14.5466,14.5484,14.5503,14.5521,14.5543,14.5575,14.5585,14.5599,14.5615,14.5624,14.5638,14.5655,14.5677,14.5691,14.5704,14.5715,14.5728,14.5741,14.5753,14.5766,14.5784,14.5795,14.5811,14.5831,14.5845,14.5864,14.588,14.5892,14.5904,14.5917,14.5934,14.5948,14.5964,14.5987,14.6006,14.6024,14.6037,14.6043,14.6052,14.606,14.6057,14.606,14.606,14.606,14.6061,14.6061,14.606,14.6059,14.6058,14.6057,14.6057,14.6059,14.6059,14.6056,14.6053,14.6054,14.6051,14.6047,14.6046,14.6046,14.6045,14.6046,14.6043,14.6042,14.6038,14.6042,14.6042,14.604,14.6038,14.6036,14.6038,14.604,14.6044,14.6044,14.605,14.6053,14.6057,14.606,14.6067,14.6071,14.6072,14.6072,14.6079,14.6089,14.6093,14.6099,14.6117,14.6122,14.613,14.6142,14.615,14.6161,14.6168,14.6176,14.6182,14.6195,14.6201,14.6217,14.6231,14.6252,14.6258,14.6273,14.6281,14.6302,14.6314,14.6322,14.6332,14.6341,14.6352,14.6363,14.6371,14.6382,14.6422,14.6441,14.646,14.6487,14.6513,14.6536,14.6564,14.659,14.6619,14.6646,14.6669,14.6682,14.67,14.6717,14.6737,14.6744,14.6777,14.6802,14.6831,14.6849,14.6865,14.6898,14.6921,14.6945,14.6964,14.6997,14.7019,14.704,14.7061,14.7091,14.7125,14.7149,14.7175,14.7193,14.7207,14.7228,14.7258,14.7282,14.7312,14.7339,14.7363,14.7395,14.7408,14.7429,14.7452,14.747,14.749,14.7519,14.7543,14.7558,14.758,14.7616,14.764,14.7661,14.7682,14.7699,14.773,14.7749,14.7768,14.7777,14.7793,14.7925,14.7936,14.7946,14.7961,14.7978,14.7991,14.8005,14.8027,14.8044,14.8057,14.8073,14.812,14.8148,14.8173,14.8187,14.8206,14.822,14.823,14.8247,14.8264,14.8281,14.83,14.8319,14.8329,14.8344,14.8363,14.838,14.8398,14.8411,14.8449,14.8468,14.8479,14.8491,14.8512,14.8529,14.8553,14.8591,14.8612,14.8632,14.8651,14.8664,14.8676,14.8688,14.8701,14.8718,14.8729,14.8742,14.8759,14.8772,14.8784,14.8796,14.8805,14.8822,14.8844,14.8858,14.8865,14.8874,14.8887,14.8905,14.8922,14.8938,14.895,14.8963,14.8988,14.9022,14.905,14.9053,14.9068,14.9082,14.9093,14.9108,14.9128,14.9143,14.9155,14.9166,14.918,14.9199,14.9219,14.9235,14.9253,14.9273,14.9292,14.9335,14.9362,14.9402,14.9451,14.948,14.9497,14.9506,14.9532,14.9553,14.9566,14.9574,14.9586,14.96,14.9618,14.9626,14.9642,14.9657,14.9673,14.9688,14.9704,14.9708,14.9724,14.9737,14.9757,14.9777,14.9796,14.9811,14.9824,14.9842,14.9862,14.9883,14.9897,14.9918,14.9945,14.9955,14.9967,14.997,14.9976,14.9978,14.9982,14.9984,14.9993,15.0008,15.0026,15.0036,15.0061,15.0084,15.0109,15.0125,15.0136,15.0143,15.0148,15.0153,15.0162,15.0176,15.0194,15.0213,15.022,15.0222,15.0235,15.0238,15.0243,15.0243,15.0243,15.0244,15.0245,15.0246,15.0246,15.0244,15.0242,15.0239,15.0236,15.0234,15.023,15.0228,15.0226,15.0222,15.0218,15.0221,15.0219,15.0212,15.021,15.021,15.0207,15.0199,15.0192,15.0189,15.018,15.0171,15.0165,15.0154,15.0144,15.0137,15.0127,15.0119,15.0107,15.0094,15.0086,15.0067,15.0057,15.0048,15.0039,15.0027,15.0018,15.001,15.0009,14.9982,14.9975,14.9969,14.9963,14.9957,14.9954,14.9951,14.9947,14.9938,14.9934,14.9925,14.9921,14.9913,14.9902,14.9903,14.9898,14.9888,14.9879,14.987,14.986,14.985,14.9844,14.9836,14.982,14.9806,14.9799,14.9783,14.9771,14.9761,14.9751,14.9732,14.9722,14.971,14.9694,14.9677,14.9668,14.965,14.9633,14.962,14.9598,14.9575,14.956,14.9546,14.953,14.9528,14.9525,14.9519,14.952,14.9527,14.9535,14.9534,14.9544,14.9558,14.9562,14.9574,14.9579,14.9589,14.9592,14.9603,14.961,14.962,14.9625,14.9628,14.9629,14.9648,14.9662,14.9682,14.9695,14.9712,14.9732,14.9752,14.9767,14.9784,14.9803,14.9827,14.9844,14.9865,14.9887,14.9908,14.9934,14.9953,15.0006,15.0015,15.0058,15.0078,15.0092,15.0115,15.0118,15.0161,15.0182,15.0219,15.024,15.0273,15.0306,15.0337,15.037,15.0417,15.0454,15.0491,15.0561,15.0617,15.0663,15.0711,15.0766,15.0805,15.0863,15.0912,15.0934,15.0954,15.0971,15.1003,15.1014,15.1046,15.1068,15.1079,15.1113,15.1147,15.1162,15.118,15.1194,15.1219,15.1247,15.1278,15.131,15.1345,15.1382,15.1443,15.15,15.1537,15.1563,15.1607,15.163,15.1658,15.1677,15.1692,15.1712,15.1741,15.1771,15.1807,15.1831,15.1855,15.1883,15.1907,15.1941,15.1971,15.2006,15.2049,15.21,15.2149,15.2176,15.2201,15.2232,15.2247,15.2288,15.2313,15.2357,15.2392,15.2425,15.2467,15.2498,15.2536,15.257,15.2596,15.2618,15.2646,15.2672,15.269,15.2696,15.2703,15.2717,15.2722,15.2725,15.273,15.2728,15.2725,15.2723,15.2706,15.2691,15.2662,15.2634,15.2596,15.2548,15.2521,15.2509,15.2488,15.249,15.25,15.2513,15.2537,15.2561,15.2579,15.2609,15.2636,15.266,15.2689,15.2715,15.2743,15.2758,15.2769,15.2793,15.282,15.2847,15.2878,15.2925,15.2957,15.2994,15.3029,15.3056,15.3088,15.3105,15.3136,15.3167,15.3196,15.3213,15.3235,15.3242,15.3239,15.3233,15.3213,15.3192,15.3174,15.3164,15.3171,15.3188,15.3212,15.3233,15.3245,15.3255,15.3257,15.3252,15.3258,15.3276,15.3311,15.3339,15.3368,15.3406,15.3431,15.3455,15.3482,15.3511,15.3536,15.356,15.3577,15.3591,15.3591,15.3582,15.3569,15.356,15.3545,15.3518,15.3499,15.3476,15.3461,15.3434,15.3417,15.3393,15.3378,15.3361,15.3341,15.3332,15.3314,15.3296,15.3281,15.3267,15.325,15.3237,15.3224,15.3208,15.3186,15.3158,15.3147,15.3133,15.3119,15.31,15.3083,15.3066,15.3045,15.3022,15.3001,15.2974,15.2953,15.2935,15.2917,15.29,15.288,15.2868,15.2861,15.2849,15.2844,15.2839,15.2836,15.2828,15.281,15.2795,15.2783,15.2766,15.2743,15.2726,15.2713,15.2703,15.2696,15.2688,15.2686,15.2687,15.2693,15.2662,15.2619,15.2585,15.2524,15.2486,15.2456,15.2389,15.2355,15.2254,15.2224,15.2141,15.2111,15.1987,15.1959,15.1924,15.1886,15.1797,15.1792,15.1783,15.1714,15.1649,15.1573,15.1548,15.154,15.1532,15.1526,15.1496,15.1457,15.1563,15.1621,15.166,15.1692,15.1746,15.1785,15.1824,15.1891,15.1957,15.1985,15.2002,15.2034,15.2086,15.2115,15.2178,15.2226,15.2247,15.226,15.2271,15.2275,15.2283,15.2304,15.2357,15.2367,15.2405,15.2409,15.2426,15.2453,15.2461,15.2486,15.2522,15.2541,15.2561,15.2603,15.2625,15.2629,15.265,15.2667,15.2685,15.2693,15.2722,15.2737,15.2768,15.2796,15.282,15.2835,15.2892,15.2903,15.2904,15.291,15.2917,15.2939,15.2949,15.2948,15.2967]}]],[[{"lng":[-16.5668,-16.5692,-16.5696,-16.5684,-16.5675,-16.5663,-16.5668],"lat":[13.6258,13.6258,13.6263,13.6275,13.6275,13.6262,13.6258]},{"lng":[-16.6396,-16.6417,-16.643,-16.643,-16.6409,-16.6359,-16.6346,-16.6346,-16.6359,-16.6367,-16.6392,-16.6396],"lat":[13.6571,13.6575,13.6562,13.6537,13.6517,13.6525,13.6537,13.6563,13.6567,13.655,13.655,13.6571]},{"lng":[-16.668,-16.6679,-16.6688,-16.6688,-16.6701,-16.6709,-16.6719,-16.6721,-16.6721,-16.6713,-16.6713,-16.6642,-16.6592,-16.6542,-16.6484,-16.6442,-16.6436,-16.6421,-16.6421,-16.6434,-16.645,-16.6484,-16.6509,-16.653,-16.6538,-16.6567,-16.6609,-16.6622,-16.6621,-16.6634,-16.6667,-16.668],"lat":[13.6554,13.6596,13.6604,13.6621,13.6633,13.6634,13.6623,13.662,13.6554,13.6546,13.6529,13.6467,13.645,13.6384,13.6359,13.6358,13.6365,13.6379,13.6404,13.6417,13.6417,13.6391,13.6392,13.6413,13.6437,13.6467,13.6475,13.6487,13.6529,13.6542,13.6542,13.6554]},{"lng":[-16.6546,-16.6559,-16.6567,-16.6608,-16.6621,-16.6576,-16.6546],"lat":[13.678,13.6792,13.6783,13.6783,13.6771,13.6767,13.678]},{"lng":[-16.6346,-16.6346,-16.6313,-16.6313,-16.6321,-16.6355,-16.6355,-16.6355,-16.6355,-16.6359,-16.6371,-16.6371,-16.638,-16.638,-16.6388,-16.6388,-16.6396,-16.6405,-16.6413,-16.6413,-16.6401,-16.6384,-16.6346],"lat":[13.6629,13.6663,13.6695,13.6721,13.6729,13.6779,13.6835,13.6861,13.6912,13.6917,13.6913,13.6837,13.6829,13.6779,13.6771,13.6721,13.6712,13.6654,13.6646,13.6612,13.66,13.66,13.6629]},{"lng":[-16.6759,-16.6722,-16.6734,-16.6768,-16.6792,-16.6808,-16.6813,-16.6784,-16.6759],"lat":[13.6934,13.6937,13.695,13.695,13.6934,13.6933,13.6921,13.6917,13.6934]},{"lng":[-16.6526,-16.6517,-16.6463,-16.6455,-16.6455,-16.6475,-16.6484,-16.6526,-16.6538,-16.653,-16.6538,-16.6549,-16.6588,-16.6588,-16.6605,-16.6605,-16.6584,-16.6526],"lat":[13.7775,13.7775,13.7829,13.7846,13.7855,13.7875,13.7875,13.7839,13.7829,13.7821,13.7804,13.7793,13.7754,13.7746,13.773,13.7721,13.7717,13.7775]},{"lng":[-16.6609,-16.6563,-16.6547,-16.6521,-16.6521,-16.6546,-16.6546,-16.6567,-16.658,-16.658,-16.6588,-16.6596,-16.6621,-16.6621,-16.6638,-16.6638,-16.6646,-16.6634,-16.6609],"lat":[13.7758,13.7804,13.7846,13.7871,13.7887,13.7904,13.7912,13.7934,13.7921,13.7904,13.7895,13.7845,13.7804,13.7795,13.7779,13.7771,13.7763,13.775,13.7758]},{"lng":[-16.6696,-16.6709,-16.6733,-16.6746,-16.6746,-16.6746,-16.6734,-16.6696],"lat":[13.7979,13.8025,13.8017,13.8004,13.7954,13.7945,13.7942,13.7979]},{"lng":[-16.65,-16.6472,-16.643,-16.6421,-16.643,-16.643,-16.645,-16.6463,-16.6555,-16.6555,-16.6563,-16.6525,-16.65],"lat":[13.7925,13.7962,13.8046,13.8055,13.807,13.8096,13.8108,13.8071,13.7979,13.7971,13.7962,13.7917,13.7925]},{"lng":[-16.638,-16.6363,-16.6346,-16.635,-16.6367,-16.6409,-16.6422,-16.6405,-16.6404,-16.6415,-16.6421,-16.6409,-16.638],"lat":[13.8054,13.8088,13.8104,13.8142,13.8134,13.8133,13.812,13.8087,13.8063,13.8041,13.8029,13.8025,13.8054]},{"lng":[-16.6117,-16.6076,-16.6051,-16.5992,-16.5967,-16.5913,-16.5905,-16.5905,-16.5951,-16.6051,-16.6084,-16.6109,-16.6134,-16.6176,-16.6209,-16.6234,-16.6276,-16.6313,-16.6321,-16.6338,-16.6363,-16.6384,-16.6418,-16.6446,-16.6446,-16.6463,-16.6463,-16.6426,-16.64,-16.6392,-16.6351,-16.6343,-16.6292,-16.6284,-16.6267,-16.6242,-16.6209,-16.6151,-16.6142,-16.6117],"lat":[13.8,13.8025,13.8034,13.81,13.8108,13.8171,13.8188,13.8229,13.8283,13.8283,13.8267,13.8267,13.825,13.8242,13.8225,13.8225,13.82,13.8154,13.8121,13.8096,13.8038,13.8017,13.8009,13.797,13.7946,13.7929,13.7888,13.785,13.785,13.7859,13.7867,13.7875,13.7875,13.7883,13.7884,13.7908,13.7925,13.7975,13.7975,13.8]},{"lng":[-16.7346,-16.7363,-16.7363,-16.7371,-16.7371,-16.7371,-16.7384,-16.7405,-16.7405,-16.7396,-16.738,-16.7351,-16.7346],"lat":[13.8413,13.8429,13.8462,13.8471,13.85,13.8554,13.8567,13.8554,13.8504,13.8496,13.8429,13.84,13.8413]},{"lng":[-16.7405,-16.7405,-16.7409,-16.7426,-16.743,-16.743,-16.7425,-16.7417,-16.7405],"lat":[13.9062,13.9071,13.9075,13.9075,13.9071,13.9054,13.905,13.905,13.9062]},{"lng":[-16.7463,-16.7463,-16.7476,-16.7497,-16.7496,-16.7505,-16.7505,-16.7492,-16.7463],"lat":[13.9196,13.9229,13.9242,13.9229,13.9187,13.9179,13.9154,13.915,13.9196]},{"lng":[-15.8028,-15.8034,-15.8104,-15.8161,-15.8201,-15.8234,-15.8264,-15.8305,-15.8342,-15.8385,-15.842,-15.8467,-15.8515,-15.855,-15.8588,-15.862,-15.8646,-15.8651,-15.8644,-15.8633,-15.8626,-15.8623,-15.8626,-15.8635,-15.8686,-15.877,-15.8837,-15.8915,-15.8942,-15.8962,-15.8976,-15.8992,-15.9009,-15.9026,-15.9048,-15.9075,-15.9103,-15.9133,-15.9147,-15.9157,-15.9163,-15.9171,-15.9196,-15.9224,-15.9244,-15.9276,-15.9331,-15.9369,-15.9434,-15.9479,-15.9532,-15.9575,-15.963,-15.9688,-15.972,-15.9744,-15.9777,-15.9811,-15.9854,-15.9873,-15.9884,-15.9891,-15.9899,-15.9904,-15.9911,-15.9916,-15.9922,-15.9929,-15.9936,-15.9944,-15.9954,-15.9962,-15.9968,-15.9973,-15.999,-16.0001,-16.0018,-16.004,-16.0058,-16.0078,-16.0099,-16.0139,-16.0156,-16.0176,-16.0264,-16.0321,-16.0348,-16.0361,-16.0377,-16.0408,-16.0439,-16.0469,-16.0495,-16.0521,-16.0545,-16.0577,-16.0602,-16.0658,-16.0681,-16.0703,-16.0728,-16.0759,-16.0794,-16.0816,-16.0843,-16.0863,-16.0879,-16.0905,-16.0921,-16.0938,-16.0948,-16.0971,-16.0992,-16.1025,-16.1049,-16.1074,-16.1095,-16.1114,-16.1158,-16.1177,-16.1213,-16.1235,-16.1284,-16.1317,-16.1345,-16.1368,-16.1395,-16.1422,-16.147,-16.1486,-16.1522,-16.1544,-16.1577,-16.1605,-16.1638,-16.1674,-16.1705,-16.1735,-16.176,-16.177,-16.1782,-16.1817,-16.1854,-16.188,-16.1904,-16.1932,-16.1945,-16.1966,-16.1988,-16.2014,-16.2036,-16.2068,-16.2101,-16.2123,-16.2153,-16.2183,-16.2206,-16.2237,-16.2253,-16.2275,-16.2294,-16.2308,-16.2311,-16.2316,-16.2322,-16.232,-16.2329,-16.2345,-16.2357,-16.237,-16.2396,-16.2408,-16.2426,-16.2454,-16.2472,-16.2483,-16.2502,-16.253,-16.255,-16.2579,-16.2601,-16.2626,-16.2654,-16.2683,-16.2708,-16.2724,-16.2746,-16.2769,-16.2786,-16.281,-16.2841,-16.2867,-16.2885,-16.2911,-16.2939,-16.2976,-16.2988,-16.3008,-16.3029,-16.3049,-16.3081,-16.3099,-16.313,-16.3188,-16.3198,-16.3225,-16.3254,-16.328,-16.3301,-16.3319,-16.3346,-16.3377,-16.3395,-16.3427,-16.3462,-16.3485,-16.3512,-16.3556,-16.3592,-16.3632,-16.3683,-16.3709,-16.3736,-16.3764,-16.3788,-16.3833,-16.3858,-16.3879,-16.3901,-16.3905,-16.394,-16.3988,-16.4022,-16.406,-16.4084,-16.4107,-16.4138,-16.4169,-16.4198,-16.4232,-16.4263,-16.4291,-16.4314,-16.4334,-16.4353,-16.438,-16.4405,-16.4424,-16.4454,-16.4469,-16.4521,-16.4553,-16.4591,-16.4636,-16.4677,-16.4713,-16.4732,-16.4742,-16.4755,-16.4768,-16.4789,-16.4815,-16.4858,-16.4889,-16.4921,-16.4939,-16.4965,-16.4961,-16.4981,-16.5002,-16.5029,-16.5066,-16.5102,-16.5135,-16.5157,-16.5202,-16.527,-16.5323,-16.5377,-16.5438,-16.5519,-16.5576,-16.5608,-16.5645,-16.5738,-16.5799,-16.5871,-16.5913,-16.5963,-16.6015,-16.607,-16.6143,-16.6211,-16.6288,-16.6349,-16.6402,-16.6444,-16.6479,-16.6509,-16.6522,-16.6542,-16.658,-16.6615,-16.6628,-16.6652,-16.6685,-16.6734,-16.6773,-16.6802,-16.6841,-16.6855,-16.6833,-16.6813,-16.6793,-16.6769,-16.6741,-16.6729,-16.6717,-16.6707,-16.6695,-16.6678,-16.6663,-16.6652,-16.6648,-16.6633,-16.6615,-16.6608,-16.6594,-16.6583,-16.6578,-16.6568,-16.656,-16.6552,-16.6545,-16.6534,-16.6524,-16.6516,-16.6495,-16.6468,-16.6451,-16.6433,-16.6419,-16.6401,-16.6384,-16.6369,-16.636,-16.6346,-16.6331,-16.6319,-16.631,-16.6304,-16.6293,-16.6285,-16.6271,-16.6267,-16.6261,-16.6256,-16.6248,-16.6248,-16.6242,-16.6242,-16.6233,-16.6237,-16.6247,-16.6276,-16.6293,-16.632,-16.6339,-16.6353,-16.6374,-16.6395,-16.6411,-16.6428,-16.6437,-16.6445,-16.6453,-16.6468,-16.6474,-16.6486,-16.6484,-16.6451,-16.6429,-16.6408,-16.6383,-16.6356,-16.6331,-16.6306,-16.6285,-16.6264,-16.6244,-16.6222,-16.6158,-16.6139,-16.613,-16.6109,-16.6093,-16.608,-16.6075,-16.607,-16.6068,-16.6064,-16.6062,-16.6058,-16.6058,-16.6055,-16.6054,-16.6053,-16.6051,-16.6048,-16.6039,-16.6039,-16.6042,-16.6044,-16.605,-16.605,-16.6053,-16.6053,-16.6056,-16.6055,-16.6054,-16.605,-16.6051,-16.6055,-16.6051,-16.6047,-16.6047,-16.6044,-16.6046,-16.6047,-16.6045,-16.6047,-16.6053,-16.6056,-16.6064,-16.6075,-16.6088,-16.6104,-16.6108,-16.6123,-16.6142,-16.6156,-16.6169,-16.6185,-16.6196,-16.6201,-16.6206,-16.6208,-16.6212,-16.6211,-16.6214,-16.6209,-16.6206,-16.6192,-16.6171,-16.6163,-16.615,-16.6134,-16.6118,-16.6102,-16.609,-16.6073,-16.6056,-16.6043,-16.6022,-16.6011,-16.5994,-16.5982,-16.5969,-16.5958,-16.5948,-16.594,-16.5926,-16.5914,-16.5905,-16.5895,-16.5883,-16.587,-16.5862,-16.5849,-16.5845,-16.5829,-16.5819,-16.5812,-16.5802,-16.5791,-16.5776,-16.5768,-16.575,-16.5737,-16.5727,-16.5701,-16.5694,-16.5682,-16.5667,-16.5654,-16.5647,-16.563,-16.562,-16.5608,-16.5591,-16.5568,-16.5562,-16.5546,-16.5533,-16.5516,-16.55,-16.5481,-16.5467,-16.5452,-16.5436,-16.5421,-16.5403,-16.5386,-16.5374,-16.5362,-16.5348,-16.5339,-16.5322,-16.5311,-16.5303,-16.5296,-16.5285,-16.5278,-16.5259,-16.5241,-16.5241,-16.5234,-16.5229,-16.5223,-16.5215,-16.5207,-16.5199,-16.5193,-16.5189,-16.5184,-16.5181,-16.518,-16.5178,-16.5177,-16.5176,-16.5176,-16.5178,-16.5177,-16.518,-16.5186,-16.519,-16.5198,-16.5206,-16.5212,-16.5217,-16.5225,-16.5234,-16.5246,-16.5257,-16.5267,-16.5273,-16.5273,-16.5271,-16.5269,-16.5268,-16.5279,-16.5285,-16.5286,-16.5287,-16.5294,-16.5299,-16.5303,-16.5308,-16.5306,-16.531,-16.5314,-16.5317,-16.5322,-16.533,-16.5335,-16.5346,-16.5349,-16.5348,-16.5358,-16.5361,-16.5361,-16.537,-16.5382,-16.5387,-16.5394,-16.5405,-16.5412,-16.5426,-16.5435,-16.5443,-16.5459,-16.5474,-16.5484,-16.5493,-16.5495,-16.5506,-16.5512,-16.5526,-16.5545,-16.5574,-16.559,-16.5611,-16.5635,-16.5653,-16.5681,-16.5694,-16.5711,-16.5735,-16.5753,-16.5771,-16.5794,-16.5814,-16.5841,-16.5867,-16.589,-16.5906,-16.592,-16.595,-16.5972,-16.5997,-16.6017,-16.6037,-16.6057,-16.6073,-16.6089,-16.6106,-16.6119,-16.614,-16.6155,-16.6175,-16.6202,-16.6217,-16.6234,-16.625,-16.6259,-16.6263,-16.628,-16.6297,-16.6308,-16.6327,-16.6336,-16.6352,-16.6377,-16.6389,-16.6404,-16.6424,-16.6439,-16.6457,-16.6466,-16.6492,-16.6516,-16.6541,-16.6568,-16.6596,-16.6625,-16.6659,-16.6684,-16.671,-16.6733,-16.6755,-16.6777,-16.6795,-16.6811,-16.6832,-16.685,-16.6869,-16.6888,-16.6906,-16.6928,-16.6948,-16.6963,-16.6981,-16.7002,-16.702,-16.7038,-16.7052,-16.7054,-16.7053,-16.7061,-16.7064,-16.7069,-16.7065,-16.7068,-16.7072,-16.7077,-16.7079,-16.7083,-16.7079,-16.7079,-16.7088,-16.7091,-16.7104,-16.7105,-16.7102,-16.71,-16.7098,-16.7103,-16.7102,-16.7097,-16.7092,-16.7089,-16.7089,-16.7079,-16.7078,-16.7079,-16.7074,-16.7075,-16.7071,-16.7073,-16.7078,-16.7085,-16.7088,-16.7096,-16.7109,-16.7124,-16.7138,-16.715,-16.7166,-16.7178,-16.7187,-16.7198,-16.7208,-16.7219,-16.7232,-16.7246,-16.7264,-16.7277,-16.729,-16.7306,-16.7318,-16.7333,-16.7344,-16.7362,-16.7372,-16.739,-16.7409,-16.7425,-16.7437,-16.7452,-16.7463,-16.7472,-16.7475,-16.7484,-16.7488,-16.7488,-16.7486,-16.7474,-16.7476,-16.746,-16.7457,-16.7441,-16.7433,-16.7431,-16.7423,-16.7422,-16.7401,-16.7387,-16.7368,-16.7348,-16.7322,-16.7312,-16.7303,-16.7281,-16.7254,-16.7239,-16.7232,-16.7228,-16.7222,-16.7223,-16.7235,-16.7241,-16.7248,-16.7256,-16.7265,-16.7278,-16.7291,-16.7304,-16.7321,-16.7329,-16.7341,-16.7356,-16.7375,-16.7392,-16.7402,-16.7413,-16.7428,-16.7437,-16.745,-16.7456,-16.7468,-16.7479,-16.7493,-16.7507,-16.7523,-16.7539,-16.7557,-16.7574,-16.7588,-16.7601,-16.7612,-16.7624,-16.7635,-16.7634,-16.763,-16.7617,-16.7597,-16.7574,-16.7565,-16.7559,-16.7551,-16.7535,-16.7533,-16.7527,-16.7526,-16.7533,-16.7537,-16.7548,-16.7559,-16.7573,-16.7585,-16.7593,-16.7615,-16.7625,-16.7633,-16.7641,-16.7655,-16.7677,-16.7693,-16.771,-16.7725,-16.7737,-16.7763,-16.7766,-16.778,-16.7798,-16.7818,-16.7835,-16.7849,-16.7998,-16.7996,-16.7996,-16.7946,-16.7938,-16.793,-16.793,-16.7913,-16.7913,-16.788,-16.7863,-16.7855,-16.7855,-16.7846,-16.7846,-16.7838,-16.7838,-16.7813,-16.778,-16.778,-16.778,-16.7771,-16.7772,-16.7763,-16.7763,-16.7755,-16.7755,-16.7746,-16.7746,-16.7738,-16.7738,-16.7746,-16.7746,-16.7738,-16.7738,-16.7721,-16.773,-16.7713,-16.7713,-16.7705,-16.7705,-16.7696,-16.7688,-16.7688,-16.768,-16.768,-16.7688,-16.7688,-16.768,-16.768,-16.7671,-16.7671,-16.768,-16.768,-16.7671,-16.768,-16.768,-16.7688,-16.768,-16.768,-16.7671,-16.7671,-16.7663,-16.7663,-16.7663,-16.7655,-16.7655,-16.7646,-16.7646,-16.7655,-16.7655,-16.7646,-16.7638,-16.763,-16.763,-16.7613,-16.7613,-16.76,-16.7584,-16.7563,-16.7563,-16.7571,-16.758,-16.7588,-16.7588,-16.758,-16.758,-16.7555,-16.7555,-16.7563,-16.7571,-16.7563,-16.7563,-16.7571,-16.7571,-16.7588,-16.7588,-16.7596,-16.763,-16.7629,-16.7638,-16.763,-16.763,-16.7621,-16.7613,-16.7605,-16.7605,-16.7596,-16.7596,-16.7605,-16.7604,-16.7613,-16.7613,-16.7588,-16.7588,-16.7571,-16.7563,-16.7546,-16.7538,-16.7509,-16.7484,-16.7451,-16.7425,-16.7418,-16.7384,-16.7317,-16.7284,-16.723,-16.7221,-16.7221,-16.723,-16.7238,-16.7238,-16.7238,-16.7255,-16.7226,-16.7192,-16.7176,-16.7155,-16.7155,-16.7146,-16.7125,-16.7109,-16.7105,-16.7092,-16.7059,-16.7025,-16.7,-16.6959,-16.6926,-16.69,-16.6884,-16.688,-16.6905,-16.6942,-16.7009,-16.7047,-16.7059,-16.7109,-16.7142,-16.7167,-16.7205,-16.7213,-16.7213,-16.7205,-16.7184,-16.7159,-16.7125,-16.7092,-16.7084,-16.7059,-16.7034,-16.6976,-16.6834,-16.6801,-16.6796,-16.6805,-16.6776,-16.6767,-16.6792,-16.6834,-16.6867,-16.69,-16.6926,-16.695,-16.7001,-16.7029,-16.7021,-16.7021,-16.7121,-16.7129,-16.7155,-16.7155,-16.7159,-16.7188,-16.7193,-16.7222,-16.7238,-16.7267,-16.7276,-16.7304,-16.7304,-16.7321,-16.743,-16.743,-16.7421,-16.7421,-16.7405,-16.738,-16.738,-16.7346,-16.7346,-16.7355,-16.7396,-16.7405,-16.7438,-16.7425,-16.7367,-16.7359,-16.7251,-16.7242,-16.7184,-16.7176,-16.7258,-16.7317,-16.7363,-16.7363,-16.7355,-16.7355,-16.7338,-16.7338,-16.733,-16.7338,-16.7338,-16.7355,-16.7355,-16.7325,-16.7292,-16.7284,-16.7251,-16.7242,-16.7221,-16.7221,-16.7234,-16.725,-16.728,-16.7271,-16.7271,-16.7271,-16.7288,-16.7288,-16.73,-16.7305,-16.7296,-16.7288,-16.7309,-16.7326,-16.7351,-16.7355,-16.7346,-16.7346,-16.7338,-16.7338,-16.7322,-16.7305,-16.7313,-16.7313,-16.7346,-16.7346,-16.7321,-16.7321,-16.7313,-16.7313,-16.728,-16.7246,-16.7213,-16.7167,-16.7142,-16.7117,-16.7101,-16.7063,-16.7013,-16.7005,-16.6963,-16.6934,-16.6917,-16.6892,-16.6884,-16.6851,-16.6834,-16.6813,-16.6805,-16.6763,-16.6767,-16.68,-16.6825,-16.6843,-16.6859,-16.6888,-16.6855,-16.6855,-16.6838,-16.6834,-16.6826,-16.6801,-16.6784,-16.6771,-16.6763,-16.6771,-16.6772,-16.6788,-16.6759,-16.6734,-16.6717,-16.6684,-16.6679,-16.6688,-16.6688,-16.668,-16.6675,-16.6659,-16.6621,-16.6621,-16.6567,-16.655,-16.6517,-16.6492,-16.6467,-16.6455,-16.6455,-16.648,-16.648,-16.6505,-16.6513,-16.6521,-16.6521,-16.6509,-16.6488,-16.6488,-16.648,-16.6471,-16.6442,-16.6434,-16.6421,-16.6418,-16.6405,-16.64,-16.6384,-16.6367,-16.6359,-16.6342,-16.6309,-16.6275,-16.6226,-16.6201,-16.6142,-16.6134,-16.6117,-16.6092,-16.6059,-16.6034,-16.6009,-16.5959,-16.5951,-16.5934,-16.5909,-16.5876,-16.5851,-16.5825,-16.5809,-16.5767,-16.5742,-16.5692,-16.568,-16.568,-16.5688,-16.5646,-16.5638,-16.5621,-16.5621,-16.5538,-16.553,-16.5509,-16.5488,-16.5488,-16.5472,-16.5471,-16.5413,-16.5421,-16.5421,-16.5446,-16.5522,-16.5521,-16.5492,-16.5475,-16.5467,-16.54,-16.5396,-16.5413,-16.5413,-16.5401,-16.5367,-16.5367,-16.5384,-16.5396,-16.5388,-16.5388,-16.5367,-16.5334,-16.5326,-16.53,-16.5292,-16.5267,-16.5259,-16.5217,-16.52,-16.5184,-16.5146,-16.513,-16.5168,-16.5201,-16.5259,-16.5267,-16.53,-16.5309,-16.5326,-16.5359,-16.5367,-16.5392,-16.5401,-16.5425,-16.5459,-16.5476,-16.5488,-16.5488,-16.5451,-16.5442,-16.5407,-16.5351,-16.5334,-16.5322,-16.5321,-16.5338,-16.5338,-16.5359,-16.5392,-16.5413,-16.5413,-16.5401,-16.5392,-16.5375,-16.5371,-16.5384,-16.54,-16.5409,-16.5433,-16.5467,-16.5476,-16.5534,-16.5542,-16.5567,-16.5575,-16.5609,-16.5617,-16.5651,-16.5659,-16.5676,-16.5692,-16.5725,-16.5746,-16.5746,-16.5717,-16.5705,-16.5767,-16.5817,-16.5851,-16.5867,-16.5938,-16.5942,-16.5984,-16.6013,-16.6013,-16.5996,-16.6009,-16.6017,-16.6051,-16.6067,-16.6113,-16.6155,-16.623,-16.6221,-16.6226,-16.6242,-16.6263,-16.6263,-16.6321,-16.6313,-16.6292,-16.6233,-16.62,-16.615,-16.6142,-16.6109,-16.61,-16.6059,-16.6051,-16.6043,-16.6084,-16.6093,-16.6109,-16.6142,-16.6168,-16.6175,-16.6259,-16.6267,-16.6301,-16.633,-16.6371,-16.6372,-16.633,-16.633,-16.6305,-16.6305,-16.6288,-16.6271,-16.6204,-16.6205,-16.6196,-16.6196,-16.6171,-16.6155,-16.6154,-16.6163,-16.6163,-16.6171,-16.6171,-16.6163,-16.6163,-16.6163,-16.6155,-16.6146,-16.6146,-16.6134,-16.61,-16.6076,-16.6063,-16.6063,-16.6068,-16.6092,-16.6109,-16.613,-16.613,-16.6101,-16.6092,-16.6076,-16.6042,-16.6009,-16.5976,-16.5875,-16.5867,-16.5851,-16.5742,-16.57,-16.5676,-16.5671,-16.5646,-16.5684,-16.5709,-16.5742,-16.5746,-16.5788,-16.5813,-16.5813,-16.5834,-16.5859,-16.5859,-16.5809,-16.5775,-16.5742,-16.5726,-16.5717,-16.5642,-16.5605,-16.5605,-16.5609,-16.5634,-16.5651,-16.5659,-16.5709,-16.5734,-16.575,-16.5771,-16.5763,-16.5763,-16.568,-16.568,-16.5638,-16.5592,-16.5591,-16.5584,-16.5576,-16.5567,-16.5558,-16.5527,-16.5511,-16.5492,-16.5477,-16.5464,-16.5457,-16.5453,-16.5462,-16.5463,-16.5462,-16.5455,-16.545,-16.5447,-16.5448,-16.5449,-16.5425,-16.5248,-16.5073,-16.4992,-16.4895,-16.4717,-16.4539,-16.4361,-16.4289,-16.4186,-16.4184,-16.4009,-16.3831,-16.365,-16.3475,-16.3298,-16.312,-16.2945,-16.2764,-16.2586,-16.2411,-16.2234,-16.2059,-16.1878,-16.17,-16.1523,-16.1348,-16.117,-16.1033,-16.0992,-16.0878,-16.0814,-16.0636,-16.0461,-16.0284,-16.0103,-16,-15.9992,-15.982,-15.9642,-15.9467,-15.9289,-15.9111,-15.8934,-15.8756,-15.8578,-15.8403,-15.8225,-15.8048,-15.787,-15.7692,-15.7514,-15.7339,-15.7159,-15.6981,-15.6806,-15.6628,-15.645,-15.6273,-15.6095,-15.6036,-15.5917,-15.5739,-15.5564,-15.5386,-15.5206,-15.5031,-15.4992,-15.4881,-15.4867,-15.4878,-15.4879,-15.4881,-15.4878,-15.4873,-15.4861,-15.4848,-15.4845,-15.4834,-15.4834,-15.4814,-15.4795,-15.4775,-15.4764,-15.4745,-15.4725,-15.4706,-15.4686,-15.4659,-15.4639,-15.4623,-15.4598,-15.4578,-15.4542,-15.4517,-15.4481,-15.4448,-15.4417,-15.4381,-15.4348,-15.4306,-15.4273,-15.4231,-15.4195,-15.4148,-15.4147,-15.4106,-15.4056,-15.4017,-15.397,-15.3928,-15.3886,-15.3853,-15.3811,-15.3761,-15.372,-15.3678,-15.3636,-15.3589,-15.3542,-15.3492,-15.3442,-15.3389,-15.3331,-15.3275,-15.3206,-15.3139,-15.3084,-15.3031,-15.2986,-15.2931,-15.2886,-15.285,-15.2814,-15.277,-15.2731,-15.2703,-15.2664,-15.2636,-15.2603,-15.2575,-15.2545,-15.2514,-15.2506,-15.2334,-15.2159,-15.1986,-15.1814,-15.1642,-15.1636,-15.1459,-15.1284,-15.1111,-15.0934,-15.0759,-15.0695,-15.0514,-15.0336,-15.0156,-14.9992,-14.9975,-14.9909,-14.9861,-14.9811,-14.9756,-14.97,-14.9642,-14.9586,-14.9523,-14.9453,-14.9389,-14.9323,-14.9259,-14.9195,-14.9136,-14.9081,-14.9023,-14.8973,-14.8923,-14.8873,-14.8823,-14.8778,-14.8734,-14.8684,-14.8631,-14.8595,-14.8556,-14.855,-14.8511,-14.847,-14.8431,-14.8407,-14.8234,-14.8213,-14.82,-14.8174,-14.8158,-14.8139,-14.8125,-14.8107,-14.8093,-14.8082,-14.8065,-14.8047,-14.8023,-14.8016,-14.7993,-14.7961,-14.7906,-14.7845,-14.7792,-14.7757,-14.7706,-14.7643,-14.7593,-14.754,-14.7497,-14.7441,-14.7362,-14.7303,-14.7239,-14.7201,-14.7186,-14.7175,-14.7161,-14.7145,-14.7135,-14.7121,-14.7105,-14.7095,-14.7092,-14.7069,-14.7042,-14.7016,-14.7009,-14.6999,-14.6973,-14.6943,-14.6919,-14.687,-14.682,-14.679,-14.675,-14.6709,-14.6673,-14.6631,-14.6587,-14.655,-14.6516,-14.6489,-14.6452,-14.6422,-14.6379,-14.6338,-14.6309,-14.6251,-14.6225,-14.6186,-14.6123,-14.609,-14.6053,-14.6043,-14.6033,-14.6013,-14.6005,-14.5992,-14.5985,-14.5972,-14.5959,-14.5945,-14.5934,-14.5922,-14.5907,-14.59,-14.5889,-14.5879,-14.5871,-14.5862,-14.5856,-14.5853,-14.5843,-14.5835,-14.5828,-14.5822,-14.5818,-14.5815,-14.5818,-14.5821,-14.5823,-14.5826,-14.5829,-14.5833,-14.5842,-14.5845,-14.5847,-14.5848,-14.5849,-14.5853,-14.5857,-14.5863,-14.5862,-14.5869,-14.5867,-14.5873,-14.5878,-14.5879,-14.5882,-14.5888,-14.5893,-14.5894,-14.5901,-14.5899,-14.5904,-14.5905,-14.5913,-14.5913,-14.5921,-14.5921,-14.5923,-14.5924,-14.5928,-14.5931,-14.5933,-14.5935,-14.593,-14.5926,-14.592,-14.592,-14.5914,-14.5912,-14.5907,-14.5907,-14.5899,-14.5893,-14.5894,-14.5889,-14.5884,-14.588,-14.5876,-14.587,-14.5867,-14.5864,-14.5863,-14.5859,-14.5852,-14.5848,-14.5845,-14.5846,-14.5847,-14.5845,-14.5846,-14.5848,-14.5854,-14.5855,-14.586,-14.5867,-14.5868,-14.5866,-14.5868,-14.5868,-14.5871,-14.5876,-14.5881,-14.5884,-14.5882,-14.5884,-14.5886,-14.5889,-14.5891,-14.5889,-14.5894,-14.5898,-14.5901,-14.5901,-14.5901,-14.5901,-14.59,-14.5903,-14.5908,-14.5909,-14.5915,-14.5912,-14.5918,-14.5919,-14.592,-14.5924,-14.5932,-14.5933,-14.5932,-14.5935,-14.5938,-14.5941,-14.5946,-14.5946,-14.595,-14.5953,-14.5954,-14.5954,-14.5961,-14.5961,-14.5975,-14.6003,-14.6029,-14.6055,-14.6085,-14.611,-14.6131,-14.6161,-14.6201,-14.6235,-14.6267,-14.6292,-14.6304,-14.6317,-14.6361,-14.6376,-14.6393,-14.6403,-14.6388,-14.6373,-14.6358,-14.6343,-14.6328,-14.6314,-14.6299,-14.6284,-14.6277,-14.6269,-14.6251,-14.6249,-14.6326,-14.6428,-14.6485,-14.6528,-14.6571,-14.6577,-14.6581,-14.6587,-14.659,-14.6628,-14.6634,-14.6642,-14.683,-14.6888,-14.6926,-14.7051,-14.7056,-14.7191,-14.7191,-14.7277,-14.7279,-14.7315,-14.7323,-14.7345,-14.7379,-14.7462,-14.7534,-14.7608,-14.7693,-14.7703,-14.7711,-14.7733,-14.7734,-14.7822,-14.7822,-14.7922,-14.796,-14.7981,-14.8002,-14.8027,-14.8028,-14.8036,-14.8111,-14.8203,-14.836,-14.8496,-14.8686,-14.8891,-14.9095,-14.9362,-14.9629,-14.9896,-14.9989,-15.0162,-15.0229,-15.0515,-15.0802,-15.0855,-15.0887,-15.1088,-15.1375,-15.1661,-15.1947,-15.2234,-15.252,-15.2807,-15.3093,-15.3366,-15.3639,-15.3912,-15.4185,-15.4458,-15.4454,-15.4446,-15.4435,-15.443,-15.4423,-15.4419,-15.4416,-15.4408,-15.4405,-15.4402,-15.4395,-15.4392,-15.4378,-15.4377,-15.4371,-15.437,-15.436,-15.4361,-15.4353,-15.4358,-15.4357,-15.4357,-15.4371,-15.4397,-15.4425,-15.4449,-15.4493,-15.4532,-15.4569,-15.4594,-15.4617,-15.4646,-15.4672,-15.4686,-15.471,-15.4731,-15.4759,-15.4786,-15.4817,-15.4844,-15.4864,-15.4884,-15.4904,-15.4923,-15.4953,-15.4973,-15.499,-15.501,-15.503,-15.5046,-15.5054,-15.5078,-15.5103,-15.513,-15.5148,-15.5175,-15.5204,-15.5238,-15.5276,-15.5313,-15.5348,-15.5364,-15.5386,-15.5408,-15.5432,-15.5452,-15.5496,-15.5523,-15.5562,-15.5597,-15.563,-15.5665,-15.5698,-15.5736,-15.5757,-15.5773,-15.5785,-15.5815,-15.5846,-15.5868,-15.5909,-15.5937,-15.5962,-15.5991,-15.6009,-15.603,-15.6064,-15.6078,-15.6099,-15.6122,-15.6148,-15.6185,-15.6218,-15.6251,-15.6272,-15.6307,-15.6322,-15.6351,-15.6379,-15.6408,-15.6434,-15.6457,-15.6484,-15.6509,-15.6544,-15.6565,-15.6587,-15.6608,-15.6624,-15.6649,-15.6666,-15.6689,-15.6711,-15.6726,-15.673,-15.6735,-15.6739,-15.6744,-15.6749,-15.6752,-15.6755,-15.676,-15.6773,-15.6841,-15.6903,-15.692,-15.6876,-15.6751,-15.6621,-15.6525,-15.6406,-15.6354,-15.6356,-15.639,-15.6435,-15.6498,-15.6554,-15.6628,-15.6664,-15.6722,-15.6791,-15.687,-15.6946,-15.7003,-15.7066,-15.7117,-15.7158,-15.7223,-15.7274,-15.7322,-15.7369,-15.7413,-15.7457,-15.7508,-15.7552,-15.7584,-15.7603,-15.7617,-15.7633,-15.7659,-15.7678,-15.7699,-15.7715,-15.772,-15.7727,-15.7729,-15.7737,-15.775,-15.7765,-15.7779,-15.7798,-15.7815,-15.7829,-15.7851,-15.7868,-15.7892,-15.793,-15.7962,-15.7989,-15.8002,-15.8018,-15.8023,-15.8028],"lat":[14.9801,14.9979,14.9993,15.0009,15.0033,15.0071,15.0123,15.0183,15.0234,15.0287,15.0327,15.0389,15.0441,15.0485,15.0511,15.0512,15.0494,15.0455,15.0417,15.0378,15.0334,15.025,15.02,15.0165,15.0109,15.0064,15.0027,14.9964,14.9934,14.9908,14.9887,14.9867,14.9849,14.9835,14.9826,14.9816,14.9802,14.9792,14.9788,14.9786,14.9788,14.9785,14.9776,14.9771,14.9768,14.9766,14.9761,14.9756,14.9749,14.9746,14.9741,14.9736,14.9731,14.9726,14.9722,14.9722,14.9722,14.9721,14.9726,14.9723,14.9721,14.9716,14.9708,14.9701,14.9686,14.967,14.9653,14.9634,14.9618,14.9599,14.9583,14.9564,14.9557,14.9553,14.9553,14.9554,14.9557,14.9555,14.9562,14.9563,14.9559,14.9563,14.9567,14.9568,14.9577,14.9586,14.9589,14.959,14.959,14.9587,14.959,14.9595,14.9593,14.9594,14.9597,14.9599,14.9598,14.9604,14.9605,14.9606,14.961,14.961,14.9611,14.9607,14.9608,14.9602,14.9597,14.9588,14.9585,14.9581,14.958,14.9575,14.9573,14.9566,14.9565,14.9565,14.9564,14.9564,14.957,14.9574,14.9579,14.9585,14.9597,14.9603,14.9609,14.9608,14.9612,14.9617,14.9621,14.9626,14.9631,14.9631,14.9628,14.9626,14.9635,14.9636,14.9632,14.9633,14.963,14.9629,14.9628,14.9625,14.962,14.961,14.9603,14.9592,14.9589,14.9579,14.9574,14.9562,14.9558,14.9544,14.9534,14.9535,14.9527,14.952,14.9519,14.9525,14.9528,14.953,14.9546,14.956,14.9575,14.9598,14.962,14.9633,14.965,14.9668,14.9677,14.9694,14.971,14.9722,14.9732,14.9751,14.9761,14.9771,14.9783,14.9799,14.9806,14.982,14.9836,14.9844,14.985,14.986,14.987,14.9879,14.9888,14.9898,14.9903,14.9902,14.9913,14.9921,14.9925,14.9934,14.9938,14.9947,14.9951,14.9954,14.9957,14.9963,14.9969,14.9975,14.9982,15.0009,15.001,15.0018,15.0027,15.0039,15.0048,15.0057,15.0067,15.0086,15.0094,15.0107,15.0119,15.0127,15.0137,15.0144,15.0154,15.0165,15.0171,15.018,15.0189,15.0192,15.0199,15.0207,15.021,15.021,15.0212,15.0219,15.0221,15.0218,15.0222,15.0226,15.0228,15.023,15.0234,15.0236,15.0239,15.0242,15.0244,15.0246,15.0246,15.0245,15.0244,15.0243,15.0243,15.0243,15.0238,15.0235,15.0222,15.022,15.0213,15.0194,15.0176,15.0162,15.0153,15.0148,15.0143,15.0136,15.0125,15.0109,15.0084,15.0061,15.0036,15.0026,15.0008,14.9993,14.9984,14.9982,14.9978,14.9976,14.997,14.9967,14.9955,14.9945,14.9918,14.9897,14.9883,14.9862,14.9842,14.9824,14.9811,14.9796,14.9777,14.9757,14.9737,14.9724,14.9708,14.9704,14.9688,14.9673,14.9657,14.9642,14.9626,14.9618,14.96,14.9586,14.9574,14.9566,14.9553,14.9532,14.9506,14.9497,14.948,14.9451,14.9402,14.9362,14.9335,14.9292,14.9273,14.9253,14.9235,14.9219,14.9199,14.918,14.9166,14.9155,14.9143,14.9128,14.9108,14.9093,14.9082,14.9068,14.9053,14.905,14.9022,14.8988,14.8963,14.895,14.8938,14.8922,14.8905,14.8887,14.8874,14.8865,14.8858,14.8844,14.8822,14.8805,14.8796,14.8784,14.8772,14.8759,14.8742,14.8729,14.8718,14.8701,14.8688,14.8676,14.8664,14.8651,14.8632,14.8612,14.8591,14.8553,14.8529,14.8512,14.8491,14.8479,14.8468,14.8449,14.8411,14.8398,14.838,14.8363,14.8344,14.8329,14.8319,14.83,14.8281,14.8264,14.8247,14.823,14.822,14.8206,14.8187,14.8173,14.8148,14.812,14.8073,14.8057,14.8044,14.8027,14.8005,14.7991,14.7978,14.7961,14.7946,14.7936,14.7925,14.7793,14.7777,14.7768,14.7749,14.773,14.7699,14.7682,14.7661,14.764,14.7616,14.758,14.7558,14.7543,14.7519,14.749,14.747,14.7452,14.7429,14.7408,14.7395,14.7363,14.7339,14.7312,14.7282,14.7258,14.7228,14.7207,14.7193,14.7175,14.7149,14.7125,14.7091,14.7061,14.704,14.7019,14.6997,14.6964,14.6945,14.6921,14.6898,14.6865,14.6849,14.6831,14.6802,14.6777,14.6744,14.6737,14.6717,14.67,14.6682,14.6669,14.6646,14.6619,14.659,14.6564,14.6536,14.6513,14.6487,14.646,14.6441,14.6422,14.6382,14.6371,14.6363,14.6352,14.6341,14.6332,14.6322,14.6314,14.6302,14.6281,14.6273,14.6258,14.6252,14.6231,14.6217,14.6201,14.6195,14.6182,14.6176,14.6168,14.6161,14.615,14.6142,14.613,14.6122,14.6117,14.6099,14.6093,14.6089,14.6079,14.6072,14.6072,14.6071,14.6067,14.606,14.6057,14.6053,14.605,14.6044,14.6044,14.604,14.6038,14.6036,14.6038,14.604,14.6042,14.6042,14.6038,14.6042,14.6043,14.6046,14.6045,14.6046,14.6046,14.6047,14.6051,14.6054,14.6053,14.6056,14.6059,14.6059,14.6057,14.6057,14.6058,14.6059,14.606,14.6061,14.6061,14.606,14.606,14.606,14.6057,14.606,14.6052,14.6043,14.6037,14.6024,14.6006,14.5987,14.5964,14.5948,14.5934,14.5917,14.5904,14.5892,14.588,14.5864,14.5845,14.5831,14.5811,14.5795,14.5784,14.5766,14.5753,14.5741,14.5728,14.5715,14.5704,14.5691,14.5677,14.5655,14.5638,14.5624,14.5615,14.5599,14.5585,14.5575,14.5543,14.5521,14.5503,14.5484,14.5466,14.5455,14.5432,14.541,14.5386,14.5373,14.535,14.5328,14.5316,14.5303,14.529,14.5276,14.5252,14.5233,14.5224,14.5203,14.5197,14.5184,14.517,14.5151,14.5132,14.5111,14.5096,14.508,14.5061,14.5044,14.503,14.5015,14.5001,14.4992,14.4963,14.4955,14.494,14.4921,14.4909,14.4891,14.4871,14.4855,14.4834,14.4814,14.4797,14.4776,14.4765,14.475,14.4735,14.4712,14.4698,14.4677,14.4663,14.4646,14.4624,14.4611,14.4601,14.4586,14.4561,14.4554,14.4545,14.4536,14.4529,14.4525,14.4523,14.4515,14.4506,14.4501,14.449,14.448,14.4472,14.4451,14.4437,14.4421,14.4406,14.4391,14.4383,14.4366,14.435,14.4335,14.4314,14.43,14.4287,14.428,14.4275,14.4271,14.4267,14.4266,14.4263,14.4262,14.4266,14.4259,14.4258,14.4261,14.426,14.4257,14.4257,14.4259,14.4256,14.4259,14.4251,14.4239,14.423,14.4222,14.4199,14.418,14.416,14.4142,14.4121,14.4098,14.4081,14.4067,14.4046,14.4028,14.4004,14.398,14.395,14.3918,14.3891,14.3859,14.3835,14.3805,14.3778,14.3752,14.373,14.3707,14.368,14.3647,14.363,14.3608,14.3584,14.3551,14.3523,14.3499,14.3473,14.3449,14.3419,14.339,14.3363,14.3332,14.3306,14.3287,14.3262,14.3237,14.321,14.3186,14.3166,14.314,14.3117,14.3097,14.3069,14.3052,14.3024,14.3004,14.2984,14.2963,14.294,14.2918,14.2894,14.2876,14.2861,14.284,14.2823,14.2807,14.2791,14.2772,14.2752,14.2731,14.2713,14.2699,14.2684,14.2669,14.2658,14.2638,14.2623,14.2599,14.2574,14.2552,14.2539,14.2524,14.2504,14.2478,14.2454,14.2432,14.2407,14.2376,14.2348,14.2317,14.2294,14.2273,14.2252,14.2231,14.2222,14.2209,14.2201,14.2189,14.2173,14.2156,14.2144,14.2126,14.2112,14.2111,14.2113,14.2104,14.2091,14.208,14.2064,14.2048,14.2035,14.2011,14.1988,14.1971,14.1952,14.1932,14.1915,14.1895,14.1889,14.1879,14.1864,14.1848,14.1829,14.1823,14.1827,14.1827,14.1825,14.182,14.1806,14.1789,14.1766,14.1754,14.1739,14.1729,14.1715,14.171,14.1708,14.1705,14.1698,14.1691,14.1674,14.1655,14.1635,14.1625,14.1609,14.1594,14.1576,14.1562,14.1546,14.1532,14.1518,14.1502,14.1483,14.1468,14.1458,14.1448,14.1432,14.1415,14.141,14.1403,14.1401,14.1407,14.1412,14.142,14.1424,14.142,14.1412,14.1403,14.1395,14.1385,14.138,14.1379,14.138,14.1387,14.1394,14.1403,14.1414,14.1411,14.1403,14.14,14.1397,14.124,14.1238,14.1221,14.1137,14.1104,14.1096,14.1079,14.1055,14.1037,14.0971,14.0913,14.0904,14.0871,14.0863,14.0746,14.0737,14.0721,14.0696,14.062,14.0587,14.0579,14.0571,14.0546,14.0537,14.0504,14.0496,14.0454,14.0446,14.0379,14.0371,14.0321,14.0312,14.0229,14.0221,14.0204,14.0188,14.0096,14.0071,14.0054,14.0046,13.9996,13.9913,13.9904,13.9863,13.9854,13.9779,13.9771,13.9738,13.9729,13.9671,13.9663,13.9579,13.9571,13.9513,13.9504,13.9496,13.9446,13.9437,13.9421,13.9346,13.9338,13.9213,13.9204,13.918,13.8971,13.8963,13.8896,13.8887,13.8771,13.8763,13.868,13.8671,13.8554,13.8546,13.8496,13.8471,13.8346,13.8333,13.8334,13.8354,13.8471,13.8479,13.8562,13.857,13.8588,13.8596,13.8621,13.8662,13.8688,13.8696,13.8729,13.8738,13.8887,13.8896,13.8921,13.8954,13.8988,13.9004,13.9038,13.9063,13.9071,13.9079,13.9129,13.9138,13.9221,13.9229,13.9354,13.9362,13.9446,13.9454,13.9487,13.9496,13.9579,13.9629,13.9646,13.9662,13.9696,13.9721,13.9762,13.9792,13.98,13.9825,13.9825,13.9817,13.9817,13.9891,13.99,13.9946,13.9963,13.9988,13.9996,14.0004,14.0174,14.0188,14.0204,14.0233,14.0233,14.0242,14.0263,14.0287,14.0304,14.0325,14.0325,14.0313,14.0309,14.0342,14.0359,14.0383,14.0384,14.04,14.0442,14.0442,14.0429,14.0379,14.0342,14.0317,14.0288,14.0267,14.0241,14.0209,14.02,14.0171,14.0154,14.0129,14.0113,14.0092,14.0092,14.0108,14.01,14.0109,14.0117,14.0142,14.0167,14.0167,14.0183,14.0171,14.0154,14.0125,14.0125,14.0062,14.0025,14.0025,14.0041,14.0067,14.0067,14.0042,14.0012,13.9996,13.9987,13.9887,13.9854,13.9829,13.9771,13.9767,13.9796,13.9809,13.978,13.9745,13.9717,13.9717,13.9688,13.9679,13.9654,13.9546,13.9521,13.9513,13.9471,13.9437,13.9404,13.9304,13.9279,13.9271,13.9254,13.9212,13.9179,13.9121,13.9108,13.9109,13.91,13.91,13.9092,13.9092,13.9083,13.9069,13.9059,13.9037,13.9012,13.9004,13.8979,13.8954,13.8871,13.8854,13.8846,13.8829,13.8813,13.8763,13.8734,13.8717,13.8733,13.8692,13.8692,13.8671,13.8663,13.8659,13.8675,13.8671,13.8662,13.8621,13.8571,13.8554,13.8538,13.8525,13.8546,13.8554,13.8604,13.8617,13.8592,13.8584,13.8554,13.8546,13.8504,13.8496,13.8471,13.8454,13.8412,13.8404,13.8379,13.8379,13.8363,13.8312,13.8262,13.8254,13.8238,13.8188,13.8155,13.8088,13.8041,13.8042,13.8025,13.8025,13.7988,13.7904,13.7854,13.7771,13.7742,13.7733,13.7733,13.7742,13.7742,13.775,13.7779,13.7804,13.7846,13.7858,13.7833,13.7825,13.7803,13.7784,13.7787,13.7812,13.7821,13.7838,13.785,13.785,13.7875,13.7875,13.7887,13.7971,13.7979,13.8029,13.8071,13.8075,13.805,13.805,13.8084,13.8079,13.8062,13.8029,13.8021,13.7991,13.7992,13.8029,13.8038,13.81,13.8108,13.8109,13.8125,13.8134,13.8146,13.8196,13.8229,13.8254,13.8288,13.8354,13.8363,13.8379,13.8384,13.8354,13.8321,13.8312,13.8271,13.8242,13.8242,13.8254,13.8284,13.8271,13.8242,13.8233,13.8234,13.8242,13.8242,13.8267,13.8275,13.83,13.8333,13.8325,13.8334,13.8333,13.835,13.8358,13.8375,13.8359,13.8358,13.8367,13.8367,13.8383,13.8392,13.8434,13.8442,13.8433,13.8475,13.8483,13.8534,13.8537,13.8571,13.8579,13.8612,13.8662,13.8679,13.8688,13.8763,13.8788,13.88,13.8829,13.8846,13.8871,13.892,13.8854,13.8846,13.8821,13.8771,13.8696,13.8662,13.8625,13.8625,13.8617,13.8625,13.8654,13.8671,13.8688,13.87,13.87,13.8692,13.8692,13.8679,13.8671,13.8663,13.8642,13.8642,13.865,13.865,13.8658,13.8658,13.8667,13.8683,13.8708,13.8708,13.8746,13.8788,13.845,13.8467,13.8467,13.8459,13.8459,13.845,13.845,13.8433,13.8425,13.8425,13.8417,13.8416,13.84,13.84,13.8388,13.8354,13.8334,13.8325,13.8325,13.8325,13.83,13.8296,13.8279,13.8262,13.8204,13.82,13.8234,13.8237,13.8271,13.8275,13.825,13.825,13.8287,13.83,13.83,13.8309,13.8309,13.8309,13.8317,13.8317,13.8308,13.8309,13.83,13.83,13.8292,13.8292,13.8284,13.8283,13.8267,13.825,13.8221,13.8171,13.8167,13.8146,13.8142,13.8083,13.8067,13.8067,13.7988,13.7975,13.795,13.7929,13.7913,13.7879,13.7875,13.7883,13.7883,13.7875,13.7821,13.7754,13.7687,13.7604,13.76,13.7617,13.7604,13.7555,13.7438,13.732,13.7308,13.7309,13.7325,13.7325,13.7333,13.7334,13.7342,13.7342,13.735,13.7312,13.7308,13.7292,13.73,13.7283,13.7283,13.7275,13.7267,13.7258,13.7258,13.7237,13.7146,13.6996,13.6912,13.6896,13.6846,13.6829,13.6795,13.6779,13.6646,13.6629,13.662,13.6579,13.6521,13.6521,13.6546,13.6554,13.6571,13.6579,13.6613,13.6621,13.664,13.6679,13.6696,13.6696,13.6637,13.6608,13.6609,13.6642,13.6637,13.6562,13.6558,13.6592,13.6592,13.6571,13.6546,13.65,13.65,13.6484,13.6475,13.645,13.6434,13.6433,13.6442,13.6442,13.6558,13.6575,13.6575,13.6596,13.6621,13.6534,13.6508,13.6492,13.6479,13.6437,13.6396,13.6379,13.6358,13.6359,13.6342,13.6342,13.6308,13.6292,13.6292,13.63,13.63,13.6263,13.6229,13.6225,13.625,13.625,13.6242,13.6241,13.6267,13.6267,13.6254,13.6238,13.6204,13.6121,13.6112,13.6071,13.6,13.6,13.5986,13.5974,13.5969,13.597,13.5996,13.6013,13.6031,13.6041,13.6045,13.6016,13.6003,13.6,13.5996,13.5982,13.5954,13.5935,13.5913,13.5902,13.5896,13.5896,13.5899,13.5899,13.59,13.5902,13.5902,13.5902,13.5904,13.5904,13.5904,13.5904,13.5907,13.5907,13.591,13.591,13.591,13.5913,13.5913,13.5915,13.5915,13.5915,13.5918,13.5918,13.5918,13.5921,13.5915,13.5915,13.5915,13.5918,13.5918,13.5918,13.5918,13.5918,13.5921,13.5921,13.5921,13.5918,13.5918,13.5921,13.5921,13.5921,13.5924,13.5915,13.5918,13.5918,13.5918,13.5913,13.5913,13.5913,13.5915,13.5907,13.591,13.591,13.591,13.5904,13.5904,13.5904,13.5907,13.5899,13.5902,13.5902,13.5902,13.5902,13.5896,13.5896,13.5896,13.5896,13.5895,13.589,13.5918,13.5979,13.6,13.6049,13.611,13.6171,13.6229,13.6285,13.6346,13.6402,13.6404,13.6457,13.6507,13.6557,13.6613,13.666,13.671,13.676,13.6807,13.6852,13.6899,13.6949,13.699,13.704,13.7077,13.7118,13.7154,13.719,13.7227,13.7263,13.7299,13.7327,13.7363,13.739,13.7429,13.7457,13.7458,13.7488,13.7515,13.7546,13.7577,13.7604,13.7635,13.7671,13.7699,13.7721,13.7752,13.7779,13.781,13.7832,13.7854,13.7877,13.7902,13.7915,13.7932,13.7949,13.7952,13.7949,13.7935,13.7918,13.789,13.7879,13.7854,13.7821,13.7788,13.7763,13.7729,13.7688,13.7657,13.7615,13.7577,13.7535,13.7496,13.7454,13.7443,13.7538,13.7635,13.7738,13.7835,13.794,13.7946,13.8002,13.8057,13.8118,13.8174,13.8238,13.826,13.819,13.8121,13.8046,13.7986,13.7979,13.7946,13.7954,13.7979,13.7993,13.801,13.8027,13.8043,13.8052,13.8054,13.8049,13.8046,13.804,13.8038,13.8027,13.8013,13.8002,13.7985,13.7965,13.7946,13.7927,13.7902,13.7877,13.7857,13.7838,13.7804,13.7782,13.7779,13.7754,13.7721,13.7693,13.7673,13.7829,13.7849,13.7881,13.7913,13.7935,13.7962,13.7979,13.8004,13.8029,13.8046,13.8064,13.8078,13.8094,13.8095,13.8098,13.8102,13.8107,13.8109,13.8109,13.8112,13.8115,13.8117,13.812,13.8126,13.813,13.814,13.8145,13.8152,13.8162,13.8174,13.8208,13.8232,13.827,13.8304,13.8331,13.8366,13.8398,13.8428,13.8436,13.8488,13.8549,13.8602,13.862,13.8651,13.8714,13.8791,13.8836,13.8858,13.8887,13.8905,13.8913,13.8939,13.8974,13.9,13.9024,13.9051,13.9077,13.9093,13.9111,13.9137,13.9166,13.9188,13.9211,13.9242,13.9263,13.9286,13.9312,13.9335,13.9366,13.9374,13.9382,13.9412,13.9436,13.9473,13.95,13.9522,13.9564,13.9599,13.963,13.9663,13.9698,13.9724,13.9761,13.9787,13.9811,13.9834,13.9855,13.9863,13.9887,13.991,13.994,13.9957,13.9984,14.0001,14.0007,14.0014,14.0022,14.0041,14.0066,14.0079,14.0097,14.0111,14.0127,14.0162,14.0192,14.0221,14.0249,14.0269,14.0307,14.0344,14.039,14.0408,14.0436,14.0462,14.0494,14.0517,14.0537,14.0572,14.0591,14.0617,14.0629,14.0676,14.0695,14.0723,14.0743,14.0773,14.0794,14.0824,14.0837,14.0854,14.0873,14.0892,14.0909,14.0948,14.0959,14.0976,14.0995,14.1021,14.1036,14.1058,14.1076,14.1102,14.1118,14.1143,14.1172,14.1181,14.1214,14.1234,14.1251,14.1278,14.13,14.1325,14.1351,14.1367,14.1378,14.1385,14.1396,14.1416,14.144,14.1465,14.1493,14.1527,14.1551,14.1581,14.161,14.162,14.1638,14.1663,14.1692,14.1721,14.1748,14.1768,14.1784,14.1812,14.1834,14.1859,14.1884,14.1909,14.1928,14.1944,14.1958,14.1971,14.1989,14.2012,14.2026,14.2057,14.2074,14.2091,14.213,14.2154,14.2173,14.2192,14.2213,14.2243,14.2253,14.2282,14.231,14.2348,14.2381,14.2402,14.2431,14.2457,14.2484,14.2514,14.2537,14.2552,14.2581,14.2592,14.2601,14.2623,14.2649,14.2676,14.2699,14.2721,14.274,14.2769,14.2812,14.2832,14.2868,14.2895,14.2906,14.2919,14.2958,14.2974,14.2986,14.301,14.3278,14.3546,14.3813,14.4081,14.4349,14.4617,14.4885,14.5153,14.5282,14.542,14.5623,14.5649,14.5672,14.5703,14.572,14.5733,14.5751,14.5754,14.5755,14.5758,14.5758,14.5761,14.5761,14.5763,14.5793,14.5804,14.5807,14.5817,14.5818,14.5827,14.5827,14.583,14.5831,14.5827,14.5831,14.583,14.5829,14.5829,14.5877,14.5926,14.598,14.5986,14.5992,14.6005,14.6005,14.6057,14.6057,14.612,14.614,14.6151,14.6162,14.6179,14.6181,14.6183,14.6206,14.6233,14.6413,14.6551,14.6551,14.6785,14.7019,14.6851,14.6684,14.6516,14.6457,14.6348,14.6097,14.615,14.6202,14.6212,14.6218,14.6255,14.6307,14.6359,14.6412,14.6464,14.6516,14.6569,14.6621,14.6786,14.695,14.7115,14.728,14.7445,14.7476,14.7532,14.7582,14.763,14.7667,14.7726,14.7766,14.781,14.7842,14.7886,14.7929,14.7966,14.8003,14.8051,14.8078,14.8101,14.8127,14.8155,14.8177,14.8208,14.8222,14.823,14.8241,14.8249,14.8254,14.8259,14.8281,14.8293,14.8308,14.8322,14.8334,14.8346,14.8355,14.8358,14.8363,14.8371,14.8377,14.8384,14.839,14.8395,14.8401,14.8401,14.8398,14.8395,14.8395,14.8393,14.8389,14.8386,14.8377,14.837,14.8366,14.8358,14.8349,14.8339,14.8333,14.8322,14.8314,14.8304,14.8287,14.8272,14.8259,14.8253,14.8245,14.8237,14.8237,14.825,14.8273,14.8292,14.8316,14.8335,14.8352,14.8378,14.8395,14.8417,14.843,14.8441,14.8448,14.8469,14.8482,14.8505,14.8531,14.8551,14.856,14.8577,14.8591,14.8597,14.8611,14.8621,14.8632,14.8634,14.8641,14.8645,14.865,14.8656,14.8659,14.8669,14.8669,14.8671,14.8673,14.8674,14.8674,14.8676,14.8686,14.8688,14.8684,14.8685,14.8687,14.8688,14.869,14.869,14.8692,14.8691,14.8697,14.8695,14.8673,14.866,14.8649,14.8635,14.8626,14.8615,14.8597,14.8572,14.8546,14.8578,14.8622,14.8721,14.8815,14.8859,14.8832,14.8794,14.8794,14.8884,14.8965,14.903,14.9102,14.9151,14.9206,14.9305,14.9421,14.9455,14.9459,14.9416,14.9314,14.9208,14.9132,14.9086,14.9062,14.9023,14.9,14.8987,14.8982,14.8977,14.8984,14.9013,14.9042,14.9085,14.9119,14.915,14.9179,14.9222,14.9249,14.9282,14.9306,14.933,14.9358,14.9391,14.9412,14.9438,14.9459,14.9479,14.9497,14.9521,14.9547,14.9568,14.9588,14.9601,14.9611,14.9619,14.9635,14.9661,14.968,14.9688,14.9801]}]],[[{"lng":[-16.7967,-16.7975,-16.798,-16.798,-16.7959,-16.7951,-16.7946,-16.7946,-16.7967],"lat":[12.4834,12.4834,12.4837,12.4846,12.4867,12.4867,12.4862,12.4854,12.4834]},{"lng":[-16.7959,-16.7947,-16.7946,-16.795,-16.7967,-16.798,-16.798,-16.7976,-16.7959],"lat":[12.4884,12.4896,12.4904,12.4908,12.4908,12.4896,12.4887,12.4884,12.4884]},{"lng":[-16.7922,-16.7921,-16.7925,-16.7942,-16.7955,-16.7942,-16.7926,-16.7922],"lat":[12.492,12.4938,12.4942,12.4942,12.4929,12.4917,12.4917,12.492]},{"lng":[-16.6746,-16.6759,-16.6825,-16.6838,-16.6825,-16.6759,-16.6746],"lat":[12.5038,12.505,12.505,12.5038,12.5025,12.5025,12.5038]},{"lng":[-16.515,-16.5138,-16.5151,-16.5192,-16.5196,-16.5176,-16.515],"lat":[12.5442,12.5471,12.5484,12.5483,12.5462,12.5441,12.5442]},{"lng":[-16.6988,-16.7001,-16.7017,-16.7059,-16.7075,-16.7092,-16.7109,-16.7122,-16.7121,-16.7193,-16.7217,-16.7238,-16.7238,-16.7209,-16.7209,-16.7284,-16.7317,-16.7358,-16.7371,-16.7372,-16.738,-16.738,-16.7388,-16.7388,-16.7396,-16.7396,-16.7371,-16.7355,-16.7346,-16.7267,-16.7242,-16.7208,-16.7168,-16.7151,-16.7142,-16.6951,-16.6942,-16.6909,-16.6901,-16.6817,-16.6788,-16.6788,-16.6796,-16.6796,-16.6805,-16.6805,-16.6796,-16.6796,-16.6788,-16.6788,-16.6817,-16.6859,-16.69,-16.693,-16.693,-16.6905,-16.6905,-16.6934,-16.6938,-16.693,-16.693,-16.6938,-16.6938,-16.6938,-16.6975,-16.6984,-16.7,-16.7034,-16.7059,-16.7071,-16.708,-16.7084,-16.7097,-16.7096,-16.7088,-16.7075,-16.7017,-16.7009,-16.6996,-16.6988,-16.6988],"lat":[12.5579,12.5592,12.5592,12.5567,12.5567,12.5592,12.5592,12.558,12.5571,12.5487,12.5458,12.5462,12.5479,12.5492,12.55,12.55,12.5525,12.5525,12.5512,12.5487,12.5479,12.5429,12.5421,12.5346,12.5337,12.5321,12.5279,12.5263,12.5229,12.515,12.5142,12.5108,12.5092,12.5092,12.5084,12.5084,12.5092,12.5092,12.51,12.51,12.5137,12.5163,12.517,12.5196,12.5204,12.5279,12.5288,12.5312,12.5321,12.5388,12.5408,12.5409,12.5434,12.5429,12.5404,12.5379,12.5371,12.535,12.5363,12.5371,12.5388,12.5396,12.5434,12.5454,12.5458,12.545,12.545,12.5475,12.5475,12.5463,12.5446,12.5408,12.5421,12.5463,12.5479,12.5492,12.5517,12.5534,12.5537,12.5555,12.5579]},{"lng":[-16.3234,-16.3226,-16.3221,-16.3221,-16.3226,-16.3242,-16.3255,-16.3242,-16.3234],"lat":[12.5892,12.5892,12.5896,12.5904,12.5908,12.5909,12.5896,12.5883,12.5892]},{"lng":[-16.3242,-16.3238,-16.3238,-16.323,-16.323,-16.3234,-16.3251,-16.3263,-16.3263,-16.3259,-16.3242],"lat":[12.5917,12.5921,12.5929,12.5937,12.5946,12.595,12.595,12.5938,12.592,12.5917,12.5917]},{"lng":[-16.3276,-16.3263,-16.3263,-16.3276,-16.33,-16.3313,-16.3313,-16.3292,-16.3276],"lat":[12.595,12.5962,12.5988,12.6,12.6,12.5987,12.5971,12.595,12.595]},{"lng":[-16.6967,-16.7009,-16.7021,-16.7021,-16.7021,-16.7009,-16.6993,-16.6984,-16.6925,-16.6896,-16.6896,-16.6918,-16.6934,-16.6967],"lat":[12.6042,12.6058,12.6046,12.5982,12.5954,12.5942,12.5942,12.595,12.595,12.5971,12.5979,12.6,12.6,12.6042]},{"lng":[-16.4663,-16.4663,-16.4684,-16.4705,-16.4705,-16.4684,-16.4663],"lat":[12.6538,12.6562,12.6592,12.6588,12.6554,12.6525,12.6538]},{"lng":[-16.7894,-16.7871,-16.7855,-16.7846,-16.7838,-16.7842,-16.7858,-16.788,-16.7888,-16.7892,-16.7896,-16.7896,-16.7894],"lat":[12.7504,12.7504,12.7538,12.7588,12.7596,12.7609,12.7608,12.7587,12.7554,12.7551,12.7546,12.7504,12.7504]},{"lng":[-12.3659,-12.3711,-12.3767,-12.3811,-12.3861,-12.3906,-12.395,-12.3986,-12.4023,-12.4059,-12.4089,-12.4125,-12.4148,-12.4184,-12.422,-12.4256,-12.43,-12.4336,-12.4375,-12.4411,-12.4448,-12.4478,-12.4481,-12.45,-12.4517,-12.4525,-12.4528,-12.4517,-12.4503,-12.4498,-12.4514,-12.4542,-12.4567,-12.4598,-12.462,-12.4642,-12.4664,-12.4681,-12.4698,-12.4714,-12.4727,-12.4728,-12.4753,-12.4781,-12.482,-12.4856,-12.4906,-12.4956,-12.4992,-12.5014,-12.5064,-12.5123,-12.5173,-12.5223,-12.5275,-12.5317,-12.5367,-12.542,-12.5467,-12.5511,-12.5564,-12.5606,-12.565,-12.5686,-12.5731,-12.5767,-12.5811,-12.5848,-12.5886,-12.5923,-12.5959,-12.5989,-12.6017,-12.6053,-12.6086,-12.6131,-12.6173,-12.6214,-12.6259,-12.6311,-12.6367,-12.6425,-12.6484,-12.6539,-12.6611,-12.6675,-12.6725,-12.6759,-12.6786,-12.682,-12.6853,-12.6856,-12.6861,-12.6641,-12.6444,-12.6247,-12.6237,-12.6235,-12.6228,-12.6219,-12.6211,-12.6213,-12.6211,-12.6212,-12.6209,-12.6213,-12.6205,-12.6192,-12.6177,-12.6163,-12.6155,-12.6146,-12.6141,-12.6141,-12.6149,-12.615,-12.6152,-12.6155,-12.6156,-12.6144,-12.614,-12.6126,-12.6117,-12.6102,-12.611,-12.612,-12.6126,-12.6129,-12.6128,-12.6118,-12.6107,-12.6098,-12.6087,-12.6074,-12.6066,-12.6058,-12.6051,-12.6053,-12.6049,-12.6042,-12.6042,-12.6038,-12.6032,-12.6031,-12.6026,-12.6022,-12.6022,-12.602,-12.6019,-12.6014,-12.6006,-12.6008,-12.6013,-12.6007,-12.6011,-12.6011,-12.6005,-12.6005,-12.6001,-12.6001,-12.6006,-12.6006,-12.6,-12.5992,-12.5987,-12.597,-12.5967,-12.596,-12.5953,-12.5948,-12.594,-12.5937,-12.5928,-12.592,-12.5914,-12.591,-12.5904,-12.5892,-12.5885,-12.5878,-12.587,-12.5866,-12.5859,-12.5942,-12.5954,-12.6049,-12.6144,-12.6147,-12.6239,-12.6334,-12.6429,-12.651,-12.6524,-12.6693,-12.6861,-12.703,-12.7198,-12.7449,-12.7701,-12.7952,-12.8211,-12.847,-12.8729,-12.8988,-12.9247,-12.9506,-12.9537,-12.9745,-12.9989,-13.0174,-13.0363,-13.0659,-13.0669,-13.0956,-13.1252,-13.1548,-13.1845,-13.2106,-13.2141,-13.2438,-13.259,-13.2742,-13.2894,-13.3046,-13.3198,-13.335,-13.3502,-13.3654,-13.3806,-13.3958,-13.4109,-13.4261,-13.4413,-13.4565,-13.4717,-13.4869,-13.5021,-13.5172,-13.5323,-13.5471,-13.577,-13.6069,-13.6368,-13.6667,-13.6966,-13.7162,-13.7265,-13.7564,-13.7834,-13.7864,-13.7868,-13.7874,-13.7886,-13.8163,-13.8194,-13.8462,-13.8717,-13.8972,-13.9227,-13.9482,-13.9645,-13.9737,-13.9948,-13.9989,-14.0016,-14.0253,-14.034,-14.0531,-14.0565,-14.1062,-14.1545,-14.1564,-14.1565,-14.1642,-14.1792,-14.1792,-14.1825,-14.1851,-14.1875,-14.2384,-14.2385,-14.2414,-14.2414,-14.2775,-14.2846,-14.2904,-14.2941,-14.3236,-14.3286,-14.33,-14.3301,-14.3301,-14.3302,-14.3513,-14.3515,-14.3681,-14.3683,-14.3799,-14.3998,-14.408,-14.4166,-14.428,-14.4356,-14.4472,-14.456,-14.4666,-14.4693,-14.4703,-14.4707,-14.4713,-14.4739,-14.4755,-14.4776,-14.4799,-14.4814,-14.4824,-14.4834,-14.4837,-14.4846,-14.4858,-14.4879,-14.4882,-14.4895,-14.4904,-14.4913,-14.4937,-14.4939,-14.494,-14.494,-14.4976,-14.4985,-14.5025,-14.513,-14.519,-14.5247,-14.5316,-14.5413,-14.5484,-14.5566,-14.5604,-14.5658,-14.5674,-14.5717,-14.5733,-14.5738,-14.5821,-14.5894,-14.5993,-14.6053,-14.6076,-14.6271,-14.6329,-14.6514,-14.6577,-14.6571,-14.6528,-14.6485,-14.6428,-14.6326,-14.6249,-14.6251,-14.6269,-14.6277,-14.6284,-14.6299,-14.6314,-14.6328,-14.6343,-14.6358,-14.6373,-14.6388,-14.6403,-14.6393,-14.6376,-14.6361,-14.6317,-14.6304,-14.6292,-14.6267,-14.6235,-14.6201,-14.6161,-14.6131,-14.611,-14.6085,-14.6055,-14.6029,-14.6003,-14.5975,-14.5961,-14.5961,-14.5954,-14.5954,-14.5953,-14.595,-14.5946,-14.5946,-14.5941,-14.5938,-14.5935,-14.5932,-14.5933,-14.5932,-14.5924,-14.592,-14.5919,-14.5918,-14.5912,-14.5915,-14.5909,-14.5908,-14.5903,-14.59,-14.5901,-14.5901,-14.5901,-14.5901,-14.5898,-14.5894,-14.5889,-14.5891,-14.5889,-14.5886,-14.5884,-14.5882,-14.5884,-14.5881,-14.5876,-14.5871,-14.5868,-14.5868,-14.5866,-14.5868,-14.5867,-14.586,-14.5855,-14.5854,-14.5848,-14.5846,-14.5845,-14.5847,-14.5846,-14.5845,-14.5848,-14.5852,-14.5859,-14.5863,-14.5864,-14.5867,-14.587,-14.5876,-14.588,-14.5884,-14.5889,-14.5894,-14.5893,-14.5899,-14.5907,-14.5907,-14.5912,-14.5914,-14.592,-14.592,-14.5926,-14.593,-14.5935,-14.5933,-14.5931,-14.5928,-14.5924,-14.5923,-14.5921,-14.5921,-14.5913,-14.5913,-14.5905,-14.5904,-14.5899,-14.5901,-14.5894,-14.5893,-14.5888,-14.5882,-14.5879,-14.5878,-14.5873,-14.5867,-14.5869,-14.5862,-14.5863,-14.5857,-14.5853,-14.5849,-14.5848,-14.5847,-14.5845,-14.5842,-14.5833,-14.5829,-14.5826,-14.5823,-14.5821,-14.5818,-14.5815,-14.5818,-14.5822,-14.5828,-14.5835,-14.5843,-14.5853,-14.5856,-14.5862,-14.5871,-14.5879,-14.5889,-14.59,-14.5907,-14.5922,-14.5934,-14.5945,-14.5959,-14.5972,-14.5985,-14.5992,-14.6005,-14.6013,-14.6033,-14.6043,-14.6053,-14.609,-14.6123,-14.6186,-14.6225,-14.6251,-14.6309,-14.6338,-14.6379,-14.6422,-14.6452,-14.6489,-14.6516,-14.655,-14.6587,-14.6631,-14.6673,-14.6709,-14.675,-14.679,-14.682,-14.687,-14.6919,-14.6943,-14.6973,-14.6999,-14.7009,-14.7016,-14.7042,-14.7069,-14.7092,-14.7095,-14.7105,-14.7121,-14.7135,-14.7145,-14.7161,-14.7175,-14.7186,-14.7201,-14.7239,-14.7303,-14.7362,-14.7441,-14.7497,-14.754,-14.7593,-14.7643,-14.7706,-14.7757,-14.7792,-14.7845,-14.7906,-14.7961,-14.7993,-14.8016,-14.8023,-14.8047,-14.8065,-14.8082,-14.8093,-14.8107,-14.8125,-14.8139,-14.8158,-14.8174,-14.82,-14.8213,-14.8234,-14.8407,-14.8392,-14.8353,-14.8317,-14.8281,-14.8245,-14.8209,-14.8184,-14.8153,-14.8131,-14.8106,-14.8086,-14.8061,-14.8039,-14.8023,-14.8006,-14.7986,-14.7973,-14.7964,-14.7945,-14.7928,-14.7925,-14.7917,-14.7914,-14.7911,-14.7731,-14.7548,-14.7367,-14.7186,-14.7178,-14.7136,-14.7103,-14.707,-14.7036,-14.7,-14.6959,-14.6925,-14.6884,-14.6842,-14.68,-14.6759,-14.6717,-14.6675,-14.6628,-14.6578,-14.6531,-14.6481,-14.6425,-14.637,-14.6311,-14.625,-14.6186,-14.6123,-14.605,-14.5981,-14.5925,-14.5867,-14.5809,-14.5753,-14.5695,-14.5645,-14.5586,-14.5536,-14.5492,-14.5445,-14.5392,-14.5348,-14.5306,-14.5267,-14.5225,-14.5184,-14.5153,-14.5109,-14.5073,-14.5042,-14.5006,-14.4992,-14.4967,-14.4945,-14.4917,-14.4892,-14.487,-14.4845,-14.4823,-14.4806,-14.4798,-14.4781,-14.477,-14.4761,-14.475,-14.4742,-14.4734,-14.4723,-14.472,-14.4661,-14.4606,-14.455,-14.4506,-14.4463,-14.4456,-14.4425,-14.44,-14.4373,-14.4342,-14.4306,-14.427,-14.4237,-14.4234,-14.4181,-14.4145,-14.41,-14.4056,-14.4014,-14.3975,-14.3934,-14.3892,-14.3848,-14.3803,-14.3759,-14.3714,-14.3667,-14.3623,-14.357,-14.3523,-14.347,-14.342,-14.337,-14.3311,-14.3267,-14.3217,-14.3161,-14.3106,-14.3048,-14.2992,-14.2936,-14.2886,-14.2839,-14.2789,-14.2725,-14.2664,-14.2606,-14.2542,-14.2473,-14.24,-14.2361,-14.2325,-14.23,-14.2273,-14.2239,-14.2206,-14.2205,-14.2203,-14.217,-14.2128,-14.2086,-14.2045,-14.1998,-14.1948,-14.1906,-14.1856,-14.18,-14.1753,-14.1689,-14.1634,-14.1575,-14.1514,-14.1456,-14.1386,-14.1323,-14.1259,-14.1231,-14.1225,-14.1198,-14.1164,-14.1175,-14.115,-14.1114,-14.1067,-14.1017,-14.097,-14.0923,-14.087,-14.0817,-14.0761,-14.0698,-14.0639,-14.0567,-14.0506,-14.0448,-14.0398,-14.0345,-14.0298,-14.0256,-14.0214,-14.0181,-14.0139,-14.0098,-14.0056,-14.0014,-13.9992,-13.9964,-13.9923,-13.9867,-13.9798,-13.9725,-13.9675,-13.9623,-13.9581,-13.9539,-13.95,-13.9436,-13.9373,-13.93,-13.9236,-13.9181,-13.9131,-13.9086,-13.9042,-13.8998,-13.8961,-13.8925,-13.8889,-13.8853,-13.8823,-13.88,-13.8775,-13.8761,-13.8745,-13.8728,-13.8698,-13.8686,-13.8678,-13.8589,-13.8501,-13.8499,-13.8492,-13.8403,-13.8334,-13.8267,-13.8211,-13.817,-13.8135,-13.8114,-13.8053,-13.8006,-13.7986,-13.7978,-13.7978,-13.7978,-13.8006,-13.8048,-13.8136,-13.8198,-13.8267,-13.8342,-13.8417,-13.8492,-13.8575,-13.8664,-13.8781,-13.8911,-13.8986,-13.9056,-13.9175,-13.93,-13.9425,-13.9539,-13.9631,-13.972,-13.9814,-13.992,-13.9992,-14.0042,-14.0173,-14.0295,-14.0406,-14.0481,-14.0584,-14.0686,-14.0823,-14.0948,-14.107,-14.115,-14.1198,-14.1253,-14.1331,-14.1417,-14.1509,-14.1603,-14.1714,-14.1817,-14.1858,-14.1934,-14.207,-14.2139,-14.2261,-14.2364,-14.2473,-14.257,-14.2673,-14.2789,-14.2911,-14.3036,-14.3106,-14.3242,-14.337,-14.3503,-14.3617,-14.372,-14.3831,-14.3928,-14.4,-14.4092,-14.4178,-14.427,-14.4336,-14.4406,-14.4459,-14.4517,-14.4611,-14.4748,-14.4881,-14.4992,-14.5003,-14.5111,-14.5203,-14.5306,-14.5381,-14.5456,-14.547,-14.5525,-14.5586,-14.5653,-14.5736,-14.5825,-14.592,-14.6025,-14.6148,-14.6264,-14.64,-14.647,-14.6598,-14.6723,-14.6834,-14.6942,-14.7045,-14.7139,-14.7228,-14.7311,-14.7381,-14.7436,-14.7489,-14.7545,-14.7592,-14.7642,-14.775,-14.7886,-14.8011,-14.8128,-14.8228,-14.8325,-14.8409,-14.847,-14.855,-14.8684,-14.8756,-14.8889,-14.9017,-14.9142,-14.9245,-14.9334,-14.9436,-14.9539,-14.9628,-14.9717,-14.98,-14.9889,-14.997,-14.9992,-14.9997,-15.0061,-15.0142,-15.0223,-15.0311,-15.0403,-15.0486,-15.0553,-15.0628,-15.0703,-15.0778,-15.0848,-15.0909,-15.0978,-15.1039,-15.1109,-15.1217,-15.1339,-15.145,-15.1548,-15.16,-15.165,-15.1703,-15.1742,-15.1767,-15.1839,-15.1923,-15.1998,-15.2017,-15.2025,-15.203,-15.203,-15.2034,-15.2034,-15.2048,-15.2039,-15.2067,-15.2081,-15.212,-15.217,-15.2239,-15.23,-15.2381,-15.2464,-15.2548,-15.2636,-15.2739,-15.2842,-15.295,-15.3075,-15.3198,-15.3334,-15.3377,-15.3459,-15.3595,-15.3717,-15.3834,-15.3936,-15.4056,-15.4173,-15.4286,-15.4398,-15.45,-15.4617,-15.4725,-15.4828,-15.4923,-15.4992,-15.5048,-15.5159,-15.5242,-15.5294,-15.5303,-15.5386,-15.5473,-15.5564,-15.567,-15.5803,-15.5939,-15.6009,-15.6084,-15.6206,-15.6323,-15.6448,-15.655,-15.6659,-15.6761,-15.6842,-15.6925,-15.7017,-15.712,-15.7239,-15.7364,-15.7503,-15.757,-15.7639,-15.7714,-15.7784,-15.7906,-15.8028,-15.8086,-15.8078,-15.8078,-15.8078,-15.8078,-15.8086,-15.8086,-15.8086,-15.8086,-15.8086,-15.8086,-15.8089,-15.8261,-15.8439,-15.8617,-15.8795,-15.8943,-15.8948,-15.8973,-15.9153,-15.9331,-15.9503,-15.9681,-15.9859,-15.9981,-15.9992,-16.0161,-16.0339,-16.0509,-16.0686,-16.0867,-16.1045,-16.1223,-16.14,-16.1573,-16.175,-16.1928,-16.2106,-16.2284,-16.2461,-16.2642,-16.2811,-16.2989,-16.3167,-16.3348,-16.3525,-16.3703,-16.3873,-16.4053,-16.4231,-16.4409,-16.4586,-16.4764,-16.4942,-16.4992,-16.5114,-16.5292,-16.547,-16.5648,-16.5828,-16.5836,-16.6006,-16.6175,-16.6353,-16.6534,-16.6711,-16.6889,-16.6998,-16.7048,-16.7075,-16.7075,-16.7142,-16.7203,-16.7286,-16.7328,-16.7306,-16.7324,-16.7315,-16.7285,-16.7285,-16.7291,-16.7302,-16.7309,-16.7315,-16.7326,-16.7342,-16.7359,-16.7373,-16.738,-16.7382,-16.7392,-16.74,-16.7412,-16.7446,-16.7459,-16.7461,-16.7471,-16.7469,-16.7469,-16.7465,-16.7454,-16.745,-16.7447,-16.7446,-16.743,-16.743,-16.7421,-16.7421,-16.7429,-16.743,-16.7421,-16.7421,-16.743,-16.7438,-16.7447,-16.7446,-16.7455,-16.7455,-16.7463,-16.7463,-16.7471,-16.7472,-16.748,-16.7488,-16.7488,-16.7504,-16.7513,-16.753,-16.753,-16.7555,-16.7547,-16.7554,-16.7555,-16.7563,-16.7563,-16.7571,-16.7571,-16.758,-16.758,-16.7588,-16.7588,-16.7596,-16.7605,-16.7613,-16.7613,-16.7621,-16.7621,-16.763,-16.763,-16.7663,-16.7671,-16.7705,-16.7713,-16.773,-16.7746,-16.7755,-16.7755,-16.7771,-16.7779,-16.7796,-16.7796,-16.7805,-16.7804,-16.7813,-16.7813,-16.7821,-16.783,-16.7837,-16.7838,-16.7846,-16.7846,-16.7855,-16.7855,-16.7863,-16.7863,-16.7871,-16.7871,-16.788,-16.7871,-16.788,-16.788,-16.7871,-16.7871,-16.7855,-16.7855,-16.7834,-16.7817,-16.7805,-16.7813,-16.7822,-16.7821,-16.7813,-16.7813,-16.7788,-16.7788,-16.778,-16.778,-16.7788,-16.7788,-16.778,-16.778,-16.7771,-16.7772,-16.7763,-16.7755,-16.7746,-16.7738,-16.7705,-16.7688,-16.768,-16.7671,-16.7671,-16.7663,-16.7663,-16.7655,-16.7663,-16.7663,-16.768,-16.768,-16.7688,-16.7721,-16.7738,-16.7738,-16.7709,-16.77,-16.7684,-16.7667,-16.7609,-16.7542,-16.7534,-16.75,-16.7484,-16.7425,-16.7384,-16.7317,-16.7313,-16.7358,-16.7417,-16.7434,-16.7442,-16.7467,-16.7484,-16.7501,-16.755,-16.7575,-16.7592,-16.7626,-16.7634,-16.7651,-16.7684,-16.7705,-16.7705,-16.7684,-16.7659,-16.7655,-16.7701,-16.773,-16.773,-16.7746,-16.7746,-16.7734,-16.7721,-16.7721,-16.7751,-16.7763,-16.7763,-16.7771,-16.7771,-16.7779,-16.7813,-16.783,-16.783,-16.7846,-16.7846,-16.7838,-16.7821,-16.7804,-16.7796,-16.7796,-16.7805,-16.7813,-16.7838,-16.7838,-16.7821,-16.7809,-16.7776,-16.7763,-16.7763,-16.7771,-16.7758,-16.7717,-16.7709,-16.7646,-16.7638,-16.7647,-16.7621,-16.7613,-16.7575,-16.7551,-16.7518,-16.7488,-16.7488,-16.7496,-16.7492,-16.7443,-16.7434,-16.7417,-16.7326,-16.7292,-16.7276,-16.7238,-16.7238,-16.7271,-16.7288,-16.7284,-16.7259,-16.7234,-16.7171,-16.7155,-16.7138,-16.7138,-16.7096,-16.7096,-16.7126,-16.7159,-16.7176,-16.7192,-16.7301,-16.7359,-16.7434,-16.7488,-16.7488,-16.7505,-16.7513,-16.7521,-16.7521,-16.7513,-16.7534,-16.7546,-16.7529,-16.753,-16.7534,-16.7542,-16.7567,-16.7605,-16.7621,-16.7638,-16.7646,-16.7663,-16.7663,-16.768,-16.768,-16.7696,-16.7705,-16.7705,-16.773,-16.7738,-16.7738,-16.7746,-16.7746,-16.7755,-16.7763,-16.7796,-16.7797,-16.775,-16.7717,-16.7701,-16.7676,-16.7642,-16.7625,-16.7601,-16.7584,-16.7571,-16.7621,-16.7634,-16.765,-16.7671,-16.7675,-16.7688,-16.7692,-16.7709,-16.773,-16.773,-16.7721,-16.773,-16.7713,-16.7696,-16.7705,-16.7696,-16.7696,-16.7688,-16.7688,-16.768,-16.7671,-16.7663,-16.7655,-16.7646,-16.7646,-16.7646,-16.7634,-16.7631,-16.7605,-16.7605,-16.7588,-16.7596,-16.7613,-16.7613,-16.7629,-16.7655,-16.7655,-16.7647,-16.7646,-16.7646,-16.763,-16.763,-16.7584,-16.7534,-16.7526,-16.745,-16.74,-16.738,-16.738,-16.7396,-16.7413,-16.7409,-16.7396,-16.7396,-16.7384,-16.7334,-16.7284,-16.7217,-16.7209,-16.7176,-16.7134,-16.7109,-16.7101,-16.705,-16.7042,-16.7017,-16.6992,-16.6967,-16.6951,-16.693,-16.693,-16.6942,-16.7009,-16.7017,-16.705,-16.7055,-16.7038,-16.7038,-16.7046,-16.7055,-16.7063,-16.7071,-16.7067,-16.7055,-16.7043,-16.7026,-16.7,-16.6992,-16.6925,-16.6859,-16.6784,-16.6709,-16.6683,-16.6671,-16.6663,-16.6646,-16.6646,-16.6675,-16.67,-16.6709,-16.6742,-16.6747,-16.6726,-16.6709,-16.6684,-16.6638,-16.6634,-16.6617,-16.6592,-16.6584,-16.6525,-16.6459,-16.6451,-16.64,-16.6392,-16.6329,-16.6329,-16.6359,-16.6368,-16.6371,-16.6359,-16.6342,-16.6317,-16.6305,-16.6305,-16.6313,-16.6313,-16.6338,-16.6338,-16.6346,-16.6355,-16.638,-16.638,-16.6388,-16.6388,-16.6396,-16.6396,-16.6394,-16.6371,-16.6351,-16.633,-16.6313,-16.6305,-16.6305,-16.6288,-16.6288,-16.6275,-16.6242,-16.6238,-16.6246,-16.6284,-16.6296,-16.6297,-16.6288,-16.628,-16.6226,-16.6209,-16.6167,-16.6159,-16.6117,-16.6109,-16.6083,-16.6075,-16.6042,-16.5988,-16.5988,-16.6013,-16.6009,-16.5975,-16.5942,-16.5917,-16.59,-16.588,-16.588,-16.5855,-16.5855,-16.5813,-16.5809,-16.5788,-16.5788,-16.5734,-16.5709,-16.5692,-16.565,-16.5588,-16.5588,-16.5576,-16.5572,-16.5563,-16.5563,-16.5555,-16.5546,-16.553,-16.5529,-16.5567,-16.5588,-16.5588,-16.5609,-16.5655,-16.5663,-16.5676,-16.5696,-16.5697,-16.5634,-16.5609,-16.5567,-16.5534,-16.5451,-16.5442,-16.5438,-16.5409,-16.5276,-16.5234,-16.5225,-16.52,-16.5184,-16.5146,-16.5146,-16.5075,-16.5063,-16.5063,-16.5013,-16.5013,-16.4976,-16.4934,-16.4921,-16.4922,-16.493,-16.4938,-16.4955,-16.4955,-16.4938,-16.4921,-16.4921,-16.493,-16.4917,-16.4876,-16.4847,-16.4855,-16.4871,-16.4871,-16.488,-16.488,-16.4901,-16.4922,-16.4905,-16.4905,-16.493,-16.493,-16.4905,-16.4904,-16.4892,-16.4876,-16.4846,-16.4846,-16.4863,-16.4859,-16.4775,-16.4746,-16.4746,-16.4759,-16.478,-16.4788,-16.4804,-16.4809,-16.4825,-16.485,-16.4855,-16.4846,-16.4834,-16.4817,-16.4792,-16.4788,-16.4788,-16.4805,-16.4805,-16.4826,-16.4851,-16.4868,-16.488,-16.488,-16.4817,-16.4776,-16.4768,-16.475,-16.473,-16.473,-16.4717,-16.47,-16.4675,-16.465,-16.4646,-16.4663,-16.468,-16.4671,-16.4651,-16.4642,-16.4621,-16.4613,-16.4592,-16.4584,-16.4571,-16.4596,-16.4596,-16.4638,-16.4638,-16.4659,-16.4701,-16.473,-16.473,-16.4738,-16.4738,-16.4775,-16.4792,-16.4805,-16.4805,-16.4796,-16.4797,-16.478,-16.478,-16.4805,-16.4813,-16.4834,-16.4884,-16.4905,-16.4905,-16.4884,-16.4842,-16.4813,-16.4796,-16.4776,-16.4717,-16.4642,-16.4592,-16.4567,-16.4542,-16.4525,-16.4513,-16.4522,-16.4538,-16.4538,-16.4546,-16.4542,-16.453,-16.4438,-16.4438,-16.4409,-16.4401,-16.4384,-16.4375,-16.4346,-16.4346,-16.4338,-16.4338,-16.4317,-16.4284,-16.4255,-16.4255,-16.4242,-16.4234,-16.4205,-16.4209,-16.4217,-16.423,-16.423,-16.428,-16.4271,-16.425,-16.4159,-16.415,-16.4134,-16.4126,-16.4059,-16.4017,-16.4009,-16.3934,-16.3918,-16.3892,-16.3792,-16.3776,-16.3767,-16.3717,-16.3709,-16.3651,-16.36,-16.3555,-16.3522,-16.3521,-16.3509,-16.3484,-16.3471,-16.3479,-16.3488,-16.3505,-16.353,-16.3521,-16.3521,-16.3505,-16.3496,-16.3497,-16.3467,-16.3434,-16.34,-16.3392,-16.335,-16.3292,-16.3284,-16.3267,-16.3259,-16.3234,-16.3225,-16.3134,-16.3125,-16.3068,-16.3042,-16.2992,-16.2959,-16.2942,-16.2909,-16.29,-16.2851,-16.2817,-16.2809,-16.275,-16.2684,-16.2663,-16.2663,-16.2688,-16.2688,-16.2676,-16.2642,-16.2617,-16.2584,-16.2576,-16.2542,-16.2509,-16.25,-16.2446,-16.2446,-16.2421,-16.2421,-16.2375,-16.2326,-16.2317,-16.2293,-16.2197,-16.2197,-16.2226,-16.2284,-16.2367,-16.24,-16.2417,-16.2467,-16.2484,-16.2501,-16.2542,-16.2592,-16.2675,-16.2684,-16.2709,-16.2751,-16.2801,-16.2826,-16.2851,-16.2867,-16.2909,-16.2918,-16.2938,-16.2955,-16.2955,-16.2976,-16.3009,-16.3017,-16.3051,-16.3067,-16.3159,-16.3176,-16.3192,-16.3225,-16.3233,-16.3259,-16.3284,-16.3309,-16.3326,-16.3359,-16.3421,-16.3429,-16.3463,-16.3463,-16.3471,-16.3471,-16.348,-16.348,-16.3505,-16.3567,-16.3609,-16.3642,-16.365,-16.3717,-16.3725,-16.3767,-16.3775,-16.3826,-16.3842,-16.3884,-16.3896,-16.3876,-16.3871,-16.3888,-16.3888,-16.3909,-16.3926,-16.3934,-16.3967,-16.3984,-16.4009,-16.4033,-16.4042,-16.4109,-16.4167,-16.4184,-16.4209,-16.425,-16.43,-16.4317,-16.4367,-16.4409,-16.4434,-16.4442,-16.4475,-16.4496,-16.4496,-16.4521,-16.4542,-16.4546,-16.4538,-16.4538,-16.4496,-16.4505,-16.4521,-16.4517,-16.4492,-16.4475,-16.4446,-16.4438,-16.4451,-16.4476,-16.4492,-16.4517,-16.4526,-16.453,-16.4555,-16.4555,-16.4584,-16.4625,-16.4634,-16.4709,-16.4717,-16.4809,-16.4817,-16.4842,-16.485,-16.4876,-16.4884,-16.4934,-16.4942,-16.4976,-16.5026,-16.5108,-16.5117,-16.5142,-16.5167,-16.5209,-16.5217,-16.5234,-16.5276,-16.5325,-16.5338,-16.5338,-16.533,-16.5321,-16.5309,-16.5293,-16.5246,-16.5259,-16.5301,-16.5305,-16.5296,-16.5296,-16.5288,-16.5288,-16.5242,-16.5184,-16.5175,-16.5158,-16.5117,-16.5076,-16.5067,-16.5038,-16.5105,-16.5113,-16.5134,-16.515,-16.5184,-16.52,-16.5217,-16.5247,-16.5251,-16.5263,-16.5263,-16.5254,-16.5255,-16.523,-16.523,-16.5251,-16.5275,-16.5301,-16.5317,-16.5321,-16.5321,-16.533,-16.533,-16.5346,-16.5347,-16.5355,-16.5355,-16.5371,-16.5371,-16.5396,-16.5396,-16.5413,-16.5413,-16.5384,-16.5371,-16.5363,-16.5347,-16.5346,-16.5372,-16.5372,-16.5492,-16.55,-16.5559,-16.5592,-16.5634,-16.5655,-16.5659,-16.5676,-16.5688,-16.5692,-16.5726,-16.575,-16.5775,-16.5784,-16.5817,-16.5851,-16.5867,-16.5913,-16.5922,-16.5921,-16.5913,-16.595,-16.5967,-16.6025,-16.6097,-16.6105,-16.6105,-16.6113,-16.6113,-16.6126,-16.615,-16.618,-16.6213,-16.6238,-16.6263,-16.6263,-16.6284,-16.633,-16.6329,-16.6397,-16.6397,-16.6413,-16.6413,-16.643,-16.6434,-16.6451,-16.6534,-16.6567,-16.6584,-16.6626,-16.6659,-16.6667,-16.67,-16.6755,-16.6755,-16.6788,-16.6788,-16.6771,-16.6771,-16.6738,-16.6738,-16.6746,-16.6755,-16.6763,-16.6763,-16.6705,-16.6696,-16.6663,-16.6613,-16.6646,-16.6634,-16.66,-16.6596,-16.6642,-16.6646,-16.6705,-16.6705,-16.6671,-16.6663,-16.665,-16.6647,-16.6667,-16.6701,-16.6718,-16.6776,-16.6834,-16.6847,-16.6792,-16.675,-16.6725,-16.6713,-16.6713,-16.6696,-16.6696,-16.668,-16.668,-16.6663,-16.6663,-16.6654,-16.6667,-16.6717,-16.6742,-16.6855,-16.6855,-16.6876,-16.6892,-16.6951,-16.6963,-16.6963,-16.6955,-16.6955,-16.6926,-16.6909,-16.6884,-16.6834,-16.6805,-16.6796,-16.6755,-16.6755,-16.6746,-16.6747,-16.6755,-16.6792,-16.6859,-16.69,-16.6917,-16.6921,-16.6901,-16.6875,-16.6863,-16.6917,-16.695,-16.6976,-16.7051,-16.7083,-16.7167,-16.7217,-16.7259,-16.7309,-16.7338,-16.7363,-16.7363,-16.7388,-16.7396,-16.7405,-16.7405,-16.7413,-16.7446,-16.7438,-16.7438,-16.7446,-16.7455,-16.748,-16.748,-16.7488,-16.75,-16.7517,-16.7534,-16.755,-16.7584,-16.7625,-16.768,-16.768,-16.7719,-16.778,-16.7805,-16.7805,-16.7822,-16.783,-16.7838,-16.7846,-16.7863,-16.7863,-16.788,-16.788,-16.7921,-16.7921,-16.793,-16.793,-16.7921,-16.7921,-16.7913,-16.7913,-16.7921,-16.7921,-16.7913,-16.7905,-16.7896,-16.7896,-16.7888,-16.7888,-16.7871,-16.7871,-16.7855,-16.7855,-16.7846,-16.7838,-16.783,-16.783,-16.7813,-16.7813,-16.7805,-16.7796,-16.7755,-16.7755,-16.7742,-16.7725,-16.7688,-16.768,-16.7679,-16.7667,-16.7634,-16.7597,-16.7584,-16.7521,-16.7513,-16.7513,-16.7501,-16.7459,-16.7442,-16.738,-16.7355,-16.7291,-16.7198,-16.7156,-16.7153,-16.7143,-16.7109,-16.7106,-16.685,-16.6692,-16.6692,-16.6692,-16.6678,-16.6675,-16.6662,-16.6655,-16.6639,-16.6656,-16.6665,-16.6666,-16.6665,-16.6445,-16.6378,-16.632,-16.6264,-16.6211,-16.617,-16.6106,-16.6036,-16.5967,-16.5903,-16.5831,-16.5767,-16.5706,-16.565,-16.5586,-16.5517,-16.5461,-16.5417,-16.5364,-16.5314,-16.5264,-16.5214,-16.5142,-16.5089,-16.5034,-16.4992,-16.4978,-16.4923,-16.4859,-16.48,-16.4748,-16.4692,-16.4628,-16.4578,-16.4517,-16.4453,-16.4384,-16.432,-16.4261,-16.4214,-16.4161,-16.4111,-16.4056,-16.3984,-16.3936,-16.3911,-16.3906,-16.3809,-16.3636,-16.3464,-16.3289,-16.312,-16.2945,-16.2773,-16.2603,-16.2428,-16.2328,-16.2259,-16.2084,-16.2011,-16.1834,-16.1653,-16.1634,-16.1592,-16.1542,-16.1503,-16.1461,-16.1406,-16.1356,-16.1309,-16.1261,-16.1206,-16.1142,-16.1078,-16.1014,-16.095,-16.0881,-16.0811,-16.0748,-16.0684,-16.062,-16.0564,-16.0506,-16.0448,-16.0392,-16.0342,-16.0289,-16.0234,-16.0184,-16.0142,-16.0098,-16.0053,-15.9998,-15.9992,-15.9945,-15.9895,-15.9845,-15.9795,-15.9745,-15.97,-15.965,-15.9614,-15.9564,-15.9392,-15.9217,-15.9039,-15.8878,-15.8864,-15.8792,-15.8614,-15.8439,-15.8261,-15.8086,-15.7906,-15.7728,-15.7553,-15.7375,-15.72,-15.702,-15.6845,-15.6831,-15.6656,-15.6484,-15.6311,-15.6139,-15.5967,-15.5792,-15.562,-15.5445,-15.5319,-15.5318,-15.5273,-15.51,-15.4992,-15.4925,-15.4753,-15.4578,-15.4406,-15.4236,-15.4186,-15.4061,-15.3945,-15.3775,-15.3606,-15.3587,-15.3434,-15.3336,-15.3164,-15.2986,-15.2814,-15.2684,-15.2514,-15.2342,-15.2173,-15.2117,-15.1939,-15.1761,-15.1586,-15.1409,-15.1364,-15.1231,-15.1056,-15.0878,-15.0703,-15.0525,-15.0348,-15.017,-14.9995,-14.9992,-14.9931,-14.9756,-14.9578,-14.9423,-14.94,-14.9225,-14.9048,-14.8873,-14.8736,-14.8692,-14.8517,-14.8342,-14.8164,-14.7986,-14.7809,-14.7634,-14.7459,-14.7281,-14.7103,-14.6925,-14.675,-14.6575,-14.6395,-14.622,-14.6045,-14.5867,-14.5689,-14.5511,-14.5336,-14.5161,-14.4992,-14.4981,-14.4806,-14.4631,-14.4453,-14.4275,-14.4098,-14.3923,-14.3745,-14.357,-14.3411,-14.3409,-14.3392,-14.3214,-14.3039,-14.2864,-14.2684,-14.2509,-14.2331,-14.2156,-14.1978,-14.18,-14.1625,-14.1589,-14.145,-14.1273,-14.1095,-14.0917,-14.0742,-14.0567,-14.0386,-14.0211,-14.0036,-13.9992,-13.9978,-13.98,-13.9625,-13.9448,-13.9273,-13.9095,-13.8906,-13.8567,-13.8389,-13.8211,-13.8034,-13.7859,-13.7684,-13.7503,-13.7328,-13.7153,-13.7123,-13.6975,-13.6798,-13.6623,-13.6445,-13.627,-13.6258,-13.6198,-13.602,-13.5845,-13.5667,-13.5517,-13.5489,-13.5314,-13.5136,-13.4992,-13.4961,-13.4781,-13.4606,-13.4468,-13.4428,-13.4253,-13.4075,-13.3898,-13.372,-13.3596,-13.3559,-13.3523,-13.345,-13.3406,-13.3361,-13.3278,-13.3225,-13.3186,-13.3148,-13.3092,-13.3069,-13.3036,-13.2984,-13.2906,-13.2861,-13.2817,-13.277,-13.2723,-13.2681,-13.2631,-13.2567,-13.2542,-13.2523,-13.2484,-13.2423,-13.235,-13.23,-13.2242,-13.2195,-13.2159,-13.2125,-13.2056,-13.2006,-13.197,-13.1925,-13.1898,-13.1859,-13.1817,-13.1759,-13.1734,-13.1692,-13.1636,-13.1578,-13.1509,-13.1453,-13.1414,-13.1386,-13.1353,-13.1317,-13.1292,-13.1261,-13.1253,-13.1225,-13.1175,-13.1117,-13.1053,-13.1006,-13.0948,-13.0884,-13.082,-13.0764,-13.07,-13.065,-13.0603,-13.0573,-13.0536,-13.0506,-13.0478,-13.0467,-13.0461,-13.045,-13.0428,-13.0411,-13.0409,-13.0423,-13.0434,-13.0461,-13.0481,-13.0506,-13.0525,-13.0548,-13.055,-13.057,-13.0606,-13.0639,-13.0673,-13.0686,-13.0684,-13.0678,-13.0656,-13.0653,-13.065,-13.0634,-13.0622,-13.0611,-13.0589,-13.0559,-13.0528,-13.0486,-13.0445,-13.04,-13.037,-13.03,-13.0245,-13.0181,-13.012,-13.0048,-12.9992,-12.9984,-12.9925,-12.9878,-12.9834,-12.9781,-12.9734,-12.9678,-12.9628,-12.9586,-12.9545,-12.9525,-12.952,-12.9509,-12.9503,-12.9501,-12.9492,-12.9486,-12.9475,-12.9461,-12.945,-12.9423,-12.9389,-12.9348,-12.9306,-12.9256,-12.9209,-12.9145,-12.9075,-12.9003,-12.8953,-12.8906,-12.887,-12.8836,-12.8809,-12.8778,-12.8742,-12.87,-12.8656,-12.8606,-12.855,-12.85,-12.845,-12.8392,-12.8328,-12.8292,-12.8298,-12.8319,-12.832,-12.8331,-12.8381,-12.8414,-12.8428,-12.8431,-12.8403,-12.8348,-12.8289,-12.8253,-12.8225,-12.8181,-12.8125,-12.8067,-12.7998,-12.7925,-12.7861,-12.7806,-12.7742,-12.7684,-12.7628,-12.7634,-12.7645,-12.7659,-12.7656,-12.7603,-12.7534,-12.747,-12.74,-12.7331,-12.7261,-12.7189,-12.7125,-12.7056,-12.6984,-12.692,-12.685,-12.6786,-12.6717,-12.6648,-12.6584,-12.6514,-12.645,-12.6378,-12.6323,-12.6286,-12.625,-12.622,-12.6192,-12.6173,-12.615,-12.6131,-12.6106,-12.6075,-12.6031,-12.597,-12.5911,-12.5875,-12.5873,-12.5895,-12.5906,-12.587,-12.5842,-12.582,-12.5781,-12.5717,-12.5675,-12.5634,-12.5592,-12.5548,-12.5498,-12.5442,-12.5392,-12.5345,-12.5303,-12.5253,-12.5203,-12.5156,-12.5106,-12.5073,-12.5017,-12.4992,-12.4961,-12.4903,-12.4842,-12.4784,-12.4725,-12.4678,-12.462,-12.4561,-12.4506,-12.445,-12.4398,-12.4356,-12.432,-12.427,-12.4228,-12.4164,-12.41,-12.4042,-12.3986,-12.3942,-12.3906,-12.3873,-12.3839,-12.3817,-12.3795,-12.3773,-12.375,-12.3728,-12.3698,-12.3661,-12.3617,-12.3575,-12.3548,-12.3523,-12.35,-12.347,-12.3448,-12.3434,-12.3386,-12.3339,-12.3289,-12.3248,-12.32,-12.315,-12.3103,-12.3048,-12.3009,-12.2998,-12.2948,-12.29,-12.2845,-12.2786,-12.2731,-12.2675,-12.2625,-12.2584,-12.2536,-12.2489,-12.2425,-12.2356,-12.2284,-12.2214,-12.2142,-12.2078,-12.2017,-12.1953,-12.1911,-12.1856,-12.1798,-12.175,-12.1695,-12.1645,-12.1595,-12.155,-12.1509,-12.1467,-12.1425,-12.1389,-12.1348,-12.1309,-12.1273,-12.1248,-12.1211,-12.1178,-12.1142,-12.1111,-12.107,-12.1034,-12.0978,-12.0909,-12.0845,-12.0773,-12.0703,-12.0631,-12.0567,-12.0511,-12.0461,-12.042,-12.0375,-12.0331,-12.0275,-12.0211,-12.0148,-12.0103,-12.0056,-12.0006,-11.9992,-11.9948,-11.9878,-11.9814,-11.975,-11.9681,-11.9625,-11.9584,-11.9542,-11.9498,-11.9453,-11.9411,-11.937,-11.932,-11.9278,-11.9228,-11.9189,-11.9148,-11.9106,-11.9073,-11.9042,-11.902,-11.8995,-11.8973,-11.8945,-11.8895,-11.8836,-11.8781,-11.8725,-11.8673,-11.8631,-11.8589,-11.8545,-11.8503,-11.845,-11.8409,-11.8348,-11.8289,-11.8234,-11.8178,-11.8123,-11.805,-11.7981,-11.7917,-11.7861,-11.7809,-11.7775,-11.7723,-11.7681,-11.7631,-11.7567,-11.7498,-11.7425,-11.7356,-11.73,-11.7242,-11.72,-11.7159,-11.7111,-11.7048,-11.6992,-11.6928,-11.6878,-11.6836,-11.6798,-11.6748,-11.6706,-11.6664,-11.6623,-11.6581,-11.6539,-11.6498,-11.6456,-11.6423,-11.6373,-11.6317,-11.6273,-11.6231,-11.6189,-11.6125,-11.6061,-11.6011,-11.5956,-11.59,-11.5839,-11.5775,-11.5711,-11.5648,-11.5586,-11.5514,-11.545,-11.5381,-11.5317,-11.5253,-11.5189,-11.512,-11.5056,-11.4992,-11.4986,-11.4923,-11.4853,-11.4789,-11.4725,-11.4661,-11.4603,-11.4548,-11.4492,-11.4436,-11.4378,-11.432,-11.4264,-11.4214,-11.4173,-11.4136,-11.4098,-11.4056,-11.4006,-11.395,-11.3898,-11.3842,-11.3792,-11.3731,-11.3723,-11.3714,-11.3725,-11.3734,-11.3742,-11.3745,-11.3739,-11.3728,-11.37,-11.3664,-11.3636,-11.3603,-11.3575,-11.3611,-11.3614,-11.3622,-11.3625,-11.3631,-11.3653,-11.3671,-11.3672,-11.3698,-11.3742,-11.3741,-11.3739,-11.3753,-11.3759,-11.3786,-11.3842,-11.3864,-11.3859,-11.3861,-11.3861,-11.3884,-11.3917,-11.3942,-11.395,-11.3967,-11.3995,-11.4036,-11.4081,-11.4123,-11.4123,-11.4106,-11.4106,-11.4136,-11.4189,-11.4225,-11.4256,-11.4284,-11.4298,-11.4309,-11.4331,-11.4373,-11.4436,-11.4478,-11.4509,-11.4489,-11.4439,-11.4392,-11.4336,-11.4286,-11.4236,-11.4195,-11.4161,-11.4156,-11.415,-11.4175,-11.42,-11.4234,-11.4256,-11.4284,-11.4292,-11.4273,-11.4253,-11.4217,-11.4216,-11.4206,-11.4195,-11.42,-11.4211,-11.4217,-11.4234,-11.4261,-11.4314,-11.4384,-11.4448,-11.4492,-11.4506,-11.4506,-11.4495,-11.4475,-11.4461,-11.4425,-11.4384,-11.4314,-11.4278,-11.4275,-11.4284,-11.4278,-11.4264,-11.4209,-11.4153,-11.4089,-11.4034,-11.3975,-11.3914,-11.3864,-11.3823,-11.3823,-11.3859,-11.3875,-11.387,-11.3886,-11.3923,-11.3964,-11.4009,-11.4009,-11.3959,-11.3925,-11.3906,-11.3892,-11.3895,-11.3911,-11.392,-11.3923,-11.3914,-11.3911,-11.3925,-11.3942,-11.3998,-11.4048,-11.4111,-11.4156,-11.417,-11.4156,-11.4131,-11.41,-11.4081,-11.4095,-11.412,-11.4136,-11.4123,-11.4103,-11.4061,-11.4034,-11.4009,-11.3986,-11.3989,-11.4011,-11.4048,-11.4084,-11.4114,-11.4109,-11.4067,-11.4017,-11.3967,-11.3923,-11.3867,-11.3803,-11.3753,-11.3711,-11.3673,-11.3667,-11.3673,-11.3684,-11.3692,-11.37,-11.3711,-11.3717,-11.3734,-11.375,-11.3773,-11.3809,-11.3856,-11.3878,-11.3884,-11.3881,-11.3892,-11.3914,-11.3939,-11.3992,-11.4034,-11.4084,-11.4139,-11.4189,-11.4225,-11.4206,-11.4175,-11.4164,-11.4175,-11.4203,-11.4217,-11.4214,-11.4185,-11.4184,-11.4175,-11.4173,-11.4131,-11.4081,-11.4056,-11.4075,-11.4125,-11.4184,-11.4214,-11.4228,-11.4236,-11.4253,-11.4289,-11.4313,-11.4331,-11.4403,-11.4425,-11.4428,-11.44,-11.4367,-11.4336,-11.4339,-11.437,-11.4414,-11.4478,-11.4548,-11.4598,-11.4656,-11.4703,-11.4739,-11.477,-11.4806,-11.485,-11.4881,-11.4914,-11.4959,-11.4992,-11.5028,-11.507,-11.5128,-11.5159,-11.5164,-11.5181,-11.5203,-11.5234,-11.527,-11.5292,-11.5306,-11.5286,-11.5253,-11.5234,-11.5223,-11.5214,-11.5217,-11.5217,-11.5248,-11.5289,-11.5353,-11.542,-11.5461,-11.5506,-11.5534,-11.555,-11.5567,-11.5573,-11.5561,-11.5542,-11.5514,-11.5486,-11.5461,-11.5434,-11.5406,-11.537,-11.5359,-11.5381,-11.5409,-11.5445,-11.547,-11.5492,-11.5514,-11.5545,-11.5595,-11.565,-11.5686,-11.5717,-11.5748,-11.5784,-11.5828,-11.5864,-11.587,-11.5881,-11.5895,-11.5925,-11.5953,-11.597,-11.597,-11.5967,-11.597,-11.597,-11.6011,-11.6039,-11.6073,-11.6131,-11.6189,-11.6245,-11.6289,-11.6311,-11.6298,-11.6292,-11.6292,-11.6317,-11.6359,-11.6423,-11.6473,-11.6509,-11.6542,-11.6606,-11.667,-11.6725,-11.6789,-11.6853,-11.6903,-11.6939,-11.6961,-11.6984,-11.7006,-11.7042,-11.7075,-11.7123,-11.7186,-11.725,-11.732,-11.7384,-11.742,-11.7434,-11.7439,-11.7459,-11.7478,-11.7503,-11.7531,-11.7559,-11.7581,-11.7592,-11.7598,-11.7586,-11.7592,-11.7634,-11.7689,-11.7756,-11.7806,-11.7839,-11.7864,-11.7886,-11.7911,-11.7939,-11.7973,-11.8014,-11.8056,-11.8106,-11.8156,-11.8192,-11.8223,-11.8239,-11.8261,-11.8275,-11.8278,-11.8289,-11.8289,-11.8325,-11.8384,-11.8431,-11.8489,-11.8539,-11.8584,-11.8611,-11.8642,-11.8673,-11.8714,-11.8786,-11.8842,-11.8873,-11.8873,-11.8873,-11.8853,-11.8834,-11.8811,-11.8792,-11.8789,-11.8784,-11.8764,-11.8731,-11.8695,-11.8667,-11.8653,-11.8656,-11.8673,-11.8673,-11.8692,-11.8711,-11.8742,-11.8786,-11.885,-11.892,-11.8975,-11.902,-11.9056,-11.9073,-11.9089,-11.9103,-11.9116,-11.9136,-11.9175,-11.9228,-11.9264,-11.93,-11.9342,-11.9375,-11.9409,-11.9459,-11.9517,-11.9586,-11.9639,-11.9681,-11.9711,-11.9734,-11.9734,-11.9761,-11.9792,-11.982,-11.9859,-11.99,-11.9931,-11.9961,-11.9989,-11.9992,-12.002,-12.0036,-12.005,-12.0067,-12.0098,-12.0125,-12.0167,-12.0211,-12.0256,-12.03,-12.0336,-12.035,-12.0367,-12.0348,-12.0334,-12.0336,-12.0359,-12.0389,-12.0425,-12.0461,-12.0498,-12.0534,-12.0564,-12.0586,-12.0609,-12.0617,-12.0611,-12.0614,-12.0636,-12.0673,-12.0717,-12.0739,-12.0753,-12.0758,-12.081,-12.0811,-12.0814,-12.0784,-12.0734,-12.0692,-12.0645,-12.0595,-12.0553,-12.0506,-12.0453,-12.0406,-12.0364,-12.0323,-12.0273,-12.0231,-12.0184,-12.0142,-12.0092,-12.005,-12.0009,-11.9992,-11.9967,-11.9925,-11.9889,-11.9848,-11.9814,-11.9781,-11.9739,-11.9703,-11.967,-11.9636,-11.9609,-11.9578,-11.9559,-11.9539,-11.952,-11.9509,-11.9498,-11.9484,-11.947,-11.9464,-11.9453,-11.9439,-11.9428,-11.9417,-11.9417,-11.9411,-11.9414,-11.9409,-11.9411,-11.9411,-11.9414,-11.9436,-11.9467,-11.9503,-11.9545,-11.9581,-11.962,-11.9656,-11.97,-11.9736,-11.9781,-11.9823,-11.9859,-11.99,-11.9939,-11.9975,-11.9992,-12.0003,-12.0042,-12.007,-12.0086,-12.0089,-12.0089,-12.0092,-12.0095,-12.01,-12.0103,-12.0125,-12.0142,-12.0136,-12.0123,-12.0098,-12.007,-12.005,-12.0053,-12.0061,-12.0056,-12.0036,-12.0017,-11.9992,-11.9975,-11.9939,-11.9903,-11.987,-11.985,-11.9845,-11.9856,-11.9878,-11.99,-11.9909,-11.99,-11.9881,-11.9861,-11.9842,-11.9817,-11.9798,-11.9792,-11.9786,-11.9786,-11.9798,-11.9814,-11.9836,-11.9859,-11.9881,-11.9898,-11.9911,-11.9928,-11.9959,-11.9986,-11.9992,-12.0009,-12.0034,-12.0056,-12.0078,-12.0095,-12.0117,-12.0134,-12.0148,-12.017,-12.0186,-12.0209,-12.0239,-12.027,-12.03,-12.0336,-12.04,-12.0456,-12.05,-12.057,-12.0628,-12.0686,-12.0728,-12.0778,-12.0823,-12.0867,-12.0909,-12.0945,-12.0975,-12.0989,-12.102,-12.105,-12.1086,-12.1123,-12.1117,-12.1061,-12.1003,-12.0961,-12.097,-12.0986,-12.0959,-12.0917,-12.0898,-12.0911,-12.0942,-12.0978,-12.1031,-12.1092,-12.115,-12.1214,-12.1245,-12.1267,-12.1295,-12.1353,-12.1411,-12.1453,-12.15,-12.1545,-12.1581,-12.1611,-12.1648,-12.1695,-12.1761,-12.1811,-12.1856,-12.1898,-12.1928,-12.2056,-12.2045,-12.2042,-12.2041,-12.2039,-12.2031,-12.2031,-12.2039,-12.204,-12.2042,-12.2036,-12.2017,-12.2003,-12.1998,-12.2014,-12.2045,-12.2061,-12.2067,-12.2081,-12.2084,-12.2086,-12.2117,-12.2159,-12.2198,-12.2225,-12.2242,-12.2249,-12.2249,-12.2253,-12.2245,-12.222,-12.2192,-12.2173,-12.2159,-12.2156,-12.2161,-12.2181,-12.2159,-12.2103,-12.2053,-12.1981,-12.1925,-12.1884,-12.1861,-12.185,-12.1867,-12.1878,-12.1878,-12.1856,-12.1823,-12.1789,-12.1784,-12.1778,-12.175,-12.1717,-12.1667,-12.1603,-12.1561,-12.1528,-12.1498,-12.1478,-12.1467,-12.147,-12.1486,-12.1514,-12.1559,-12.1609,-12.1678,-12.1736,-12.1778,-12.1803,-12.1814,-12.177,-12.1714,-12.1678,-12.1653,-12.1634,-12.1642,-12.1706,-12.1739,-12.1781,-12.1811,-12.182,-12.1873,-12.1928,-12.1984,-12.205,-12.2092,-12.2117,-12.2139,-12.2148,-12.2136,-12.2123,-12.2153,-12.22,-12.2253,-12.2311,-12.2361,-12.2398,-12.2434,-12.2442,-12.2442,-12.2441,-12.2439,-12.2439,-12.2475,-12.252,-12.2567,-12.2611,-12.2656,-12.2706,-12.2745,-12.2784,-12.2828,-12.2864,-12.2903,-12.2945,-12.2984,-12.302,-12.3056,-12.3092,-12.3125,-12.3162,-12.3189,-12.322,-12.3256,-12.33,-12.3348,-12.3392,-12.3445,-12.3503,-12.3559,-12.3603,-12.3659],"lat":[14.8396,14.8415,14.8427,14.8454,14.8474,14.8499,14.8527,14.856,14.8593,14.8627,14.8665,14.8699,14.8746,14.8779,14.8813,14.8849,14.8874,14.8907,14.894,14.8974,14.9007,14.9046,14.9052,14.9096,14.9149,14.9213,14.9279,14.9335,14.939,14.9454,14.9507,14.9549,14.9596,14.9638,14.9685,14.9732,14.9779,14.9835,14.9888,14.9943,14.9996,14.9999,15.0046,15.0085,15.0118,15.0152,15.0171,15.019,15.0199,15.0204,15.0221,15.0235,15.0254,15.0274,15.0293,15.0318,15.0338,15.0357,15.0377,15.0402,15.0421,15.0446,15.0474,15.0507,15.0532,15.0565,15.0593,15.0627,15.066,15.0693,15.0727,15.0768,15.0807,15.084,15.0882,15.0907,15.0935,15.096,15.0985,15.1004,15.1018,15.1029,15.104,15.1054,15.1052,15.1043,15.1021,15.0985,15.0943,15.0907,15.0877,15.0871,15.087,15.0543,15.0273,15.0004,14.9973,14.9955,14.9933,14.9916,14.9905,14.9882,14.9872,14.9861,14.985,14.9837,14.9829,14.9824,14.9822,14.9813,14.9808,14.98,14.9794,14.9783,14.9775,14.9766,14.9759,14.9751,14.9742,14.973,14.9726,14.9724,14.9719,14.9703,14.969,14.9682,14.9675,14.9666,14.9656,14.965,14.9647,14.964,14.9634,14.9632,14.9628,14.9623,14.9608,14.9593,14.958,14.9568,14.956,14.9552,14.953,14.9521,14.9512,14.9502,14.9492,14.9478,14.9466,14.9454,14.944,14.943,14.9412,14.9401,14.9386,14.9375,14.9363,14.9356,14.9341,14.9333,14.9321,14.9308,14.9294,14.9279,14.9263,14.9245,14.9229,14.922,14.9206,14.9194,14.9179,14.9169,14.9151,14.9139,14.9129,14.912,14.911,14.91,14.9091,14.9079,14.9071,14.9057,14.9052,14.8798,14.8761,14.8474,14.8188,14.818,14.7902,14.7616,14.7329,14.7087,14.7043,14.6816,14.6589,14.6362,14.6135,14.6025,14.5915,14.5804,14.5767,14.573,14.5694,14.5657,14.562,14.5583,14.5592,14.5654,14.5726,14.5764,14.5803,14.5847,14.5849,14.5892,14.5936,14.5981,14.6025,14.6064,14.6069,14.6114,14.5999,14.5884,14.577,14.5655,14.554,14.5425,14.5311,14.5196,14.5081,14.4967,14.4852,14.4737,14.4622,14.4507,14.4393,14.4278,14.4163,14.4048,14.4052,14.4055,14.4062,14.4069,14.4076,14.4083,14.4089,14.4094,14.4096,14.4103,14.4109,14.4109,14.4109,14.411,14.411,14.4116,14.4117,14.4123,14.4246,14.437,14.4493,14.4617,14.4696,14.4741,14.484,14.4859,14.4864,14.4911,14.493,14.4973,14.4984,14.5145,14.5304,14.531,14.531,14.5335,14.5385,14.5385,14.5395,14.5404,14.5412,14.5578,14.5578,14.5588,14.5588,14.5709,14.573,14.5747,14.5758,14.5845,14.586,14.5864,14.5864,14.5864,14.5864,14.5838,14.5838,14.5823,14.5823,14.5812,14.5801,14.58,14.5799,14.5791,14.5802,14.5807,14.579,14.5787,14.5787,14.5779,14.5776,14.577,14.5743,14.5724,14.5702,14.5675,14.5665,14.5658,14.5655,14.5654,14.5651,14.5647,14.564,14.5639,14.5634,14.563,14.5626,14.562,14.562,14.5619,14.5619,14.5621,14.5621,14.5623,14.5627,14.5632,14.5628,14.5623,14.562,14.5604,14.561,14.561,14.5627,14.5638,14.5667,14.5678,14.568,14.5717,14.5742,14.5745,14.5746,14.5747,14.5738,14.5735,14.5745,14.5754,14.5751,14.5733,14.572,14.5703,14.5672,14.5649,14.5623,14.542,14.5282,14.5153,14.4885,14.4617,14.4349,14.4081,14.3813,14.3546,14.3278,14.301,14.2986,14.2974,14.2958,14.2919,14.2906,14.2895,14.2868,14.2832,14.2812,14.2769,14.274,14.2721,14.2699,14.2676,14.2649,14.2623,14.2601,14.2592,14.2581,14.2552,14.2537,14.2514,14.2484,14.2457,14.2431,14.2402,14.2381,14.2348,14.231,14.2282,14.2253,14.2243,14.2213,14.2192,14.2173,14.2154,14.213,14.2091,14.2074,14.2057,14.2026,14.2012,14.1989,14.1971,14.1958,14.1944,14.1928,14.1909,14.1884,14.1859,14.1834,14.1812,14.1784,14.1768,14.1748,14.1721,14.1692,14.1663,14.1638,14.162,14.161,14.1581,14.1551,14.1527,14.1493,14.1465,14.144,14.1416,14.1396,14.1385,14.1378,14.1367,14.1351,14.1325,14.13,14.1278,14.1251,14.1234,14.1214,14.1181,14.1172,14.1143,14.1118,14.1102,14.1076,14.1058,14.1036,14.1021,14.0995,14.0976,14.0959,14.0948,14.0909,14.0892,14.0873,14.0854,14.0837,14.0824,14.0794,14.0773,14.0743,14.0723,14.0695,14.0676,14.0629,14.0617,14.0591,14.0572,14.0537,14.0517,14.0494,14.0462,14.0436,14.0408,14.039,14.0344,14.0307,14.0269,14.0249,14.0221,14.0192,14.0162,14.0127,14.0111,14.0097,14.0079,14.0066,14.0041,14.0022,14.0014,14.0007,14.0001,13.9984,13.9957,13.994,13.991,13.9887,13.9863,13.9855,13.9834,13.9811,13.9787,13.9761,13.9724,13.9698,13.9663,13.963,13.9599,13.9564,13.9522,13.95,13.9473,13.9436,13.9412,13.9382,13.9374,13.9366,13.9335,13.9312,13.9286,13.9263,13.9242,13.9211,13.9188,13.9166,13.9137,13.9111,13.9093,13.9077,13.9051,13.9024,13.9,13.8974,13.8939,13.8913,13.8905,13.8887,13.8858,13.8836,13.8791,13.8714,13.8651,13.862,13.8602,13.8549,13.8488,13.8436,13.8428,13.8398,13.8366,13.8331,13.8304,13.827,13.8232,13.8208,13.8174,13.8162,13.8152,13.8145,13.814,13.813,13.8126,13.812,13.8117,13.8115,13.8112,13.8109,13.8109,13.8107,13.8102,13.8098,13.8095,13.8094,13.8078,13.8064,13.8046,13.8029,13.8004,13.7979,13.7962,13.7935,13.7913,13.7881,13.7849,13.7829,13.7673,13.766,13.7629,13.7596,13.7563,13.7529,13.7496,13.7449,13.7407,13.736,13.7313,13.7265,13.7218,13.7171,13.7118,13.7063,13.7007,13.6954,13.6893,13.6838,13.6782,13.6713,13.6652,13.6585,13.6529,13.6432,13.6329,13.6227,13.6121,13.6115,13.6143,13.6179,13.6215,13.6252,13.6288,13.6318,13.6352,13.6382,13.641,13.644,13.6468,13.6499,13.6529,13.6552,13.6574,13.6596,13.6618,13.6635,13.6649,13.6665,13.6674,13.6685,13.6693,13.6696,13.6699,13.6688,13.6677,13.6665,13.6652,13.664,13.6624,13.661,13.6593,13.6565,13.6546,13.6529,13.6502,13.6477,13.6443,13.6418,13.639,13.6365,13.634,13.6307,13.6279,13.6246,13.6234,13.6213,13.6165,13.6127,13.6079,13.6032,13.5985,13.5935,13.5882,13.5821,13.5765,13.5704,13.5643,13.5579,13.5518,13.5457,13.5396,13.534,13.5329,13.5318,13.5304,13.5279,13.5263,13.526,13.5221,13.5174,13.5132,13.5093,13.506,13.5027,13.4997,13.4993,13.4974,13.494,13.4915,13.489,13.4863,13.4829,13.4804,13.4779,13.4752,13.4727,13.4699,13.4674,13.4654,13.4629,13.461,13.4593,13.4574,13.4554,13.4535,13.4538,13.456,13.4582,13.4599,13.4615,13.4629,13.4646,13.4663,13.4685,13.4707,13.4729,13.4738,13.4749,13.4763,13.4774,13.4777,13.4779,13.4807,13.4843,13.4885,13.4927,13.4963,13.4996,13.4997,13.4999,13.5035,13.5063,13.5093,13.5121,13.5143,13.5168,13.5196,13.5218,13.5235,13.5257,13.5265,13.5282,13.5299,13.5307,13.5324,13.5327,13.5335,13.5343,13.5352,13.5352,13.5379,13.5415,13.5396,13.5424,13.546,13.5488,13.551,13.554,13.5563,13.5585,13.5602,13.5618,13.5627,13.5629,13.5632,13.5627,13.5615,13.5596,13.5577,13.5557,13.5574,13.5602,13.5638,13.5665,13.5696,13.5727,13.5754,13.5764,13.5777,13.5804,13.5821,13.5824,13.5827,13.5807,13.5788,13.5763,13.5738,13.5704,13.5699,13.5693,13.5696,13.5693,13.5679,13.5663,13.5635,13.561,13.5582,13.5552,13.5518,13.5485,13.5452,13.541,13.5363,13.5315,13.526,13.5207,13.5152,13.511,13.5091,13.5077,13.5035,13.4997,13.4996,13.4993,13.494,13.4874,13.4802,13.4715,13.4621,13.4579,13.4554,13.4474,13.4382,13.426,13.4193,13.4127,13.406,13.3888,13.3793,13.3674,13.3593,13.3521,13.346,13.3393,13.3327,13.3274,13.3221,13.3188,13.3174,13.3174,13.3174,13.3188,13.3193,13.3215,13.3207,13.316,13.3107,13.306,13.3021,13.301,13.3002,13.2988,13.3002,13.3007,13.294,13.2902,13.2868,13.2854,13.2846,13.2868,13.2807,13.2713,13.2635,13.2568,13.2513,13.246,13.2421,13.2379,13.234,13.2331,13.2313,13.2307,13.2302,13.2321,13.2346,13.2379,13.2427,13.246,13.2435,13.2421,13.2393,13.2393,13.2379,13.2368,13.2374,13.2402,13.2427,13.246,13.2507,13.256,13.2613,13.2668,13.2721,13.2788,13.286,13.2949,13.3027,13.306,13.3049,13.3049,13.3066,13.3068,13.3093,13.3135,13.3215,13.3282,13.334,13.3354,13.3415,13.3488,13.356,13.3621,13.3588,13.354,13.3502,13.3474,13.3449,13.344,13.3435,13.3449,13.3468,13.3493,13.3521,13.356,13.3602,13.3654,13.3715,13.3788,13.3868,13.3954,13.404,13.414,13.4227,13.4254,13.424,13.4254,13.4274,13.4315,13.4349,13.4407,13.4488,13.4549,13.4549,13.4549,13.4549,13.456,13.4574,13.4602,13.464,13.4682,13.4715,13.4768,13.4821,13.4882,13.4927,13.4979,13.4993,13.4996,13.5035,13.5093,13.5146,13.5193,13.5246,13.5307,13.5374,13.544,13.5507,13.5574,13.5646,13.5721,13.5799,13.5879,13.5954,13.5935,13.5907,13.5874,13.5827,13.5746,13.5654,13.5574,13.5524,13.5493,13.5435,13.5374,13.5307,13.5188,13.5068,13.4998,13.4996,13.494,13.4807,13.4688,13.4615,13.4507,13.4388,13.4288,13.4202,13.4135,13.4054,13.3993,13.3935,13.3882,13.3827,13.3782,13.374,13.3707,13.3688,13.3668,13.3654,13.3656,13.366,13.366,13.3682,13.3693,13.3715,13.3735,13.3749,13.3774,13.3802,13.3835,13.3854,13.3888,13.3921,13.396,13.3968,13.3974,13.394,13.3888,13.3819,13.3807,13.3749,13.3693,13.364,13.3607,13.3593,13.3582,13.3574,13.3574,13.3582,13.3607,13.3621,13.3649,13.3682,13.3715,13.3674,13.3615,13.356,13.3521,13.3493,13.3482,13.3468,13.346,13.346,13.3454,13.3454,13.3468,13.3474,13.3307,13.3135,13.2968,13.2793,13.2621,13.2449,13.2446,13.2274,13.2102,13.1927,13.176,13.1593,13.1593,13.1602,13.1593,13.1602,13.1602,13.1602,13.1602,13.1607,13.1613,13.1613,13.1621,13.1621,13.1627,13.1627,13.1627,13.1627,13.1627,13.1621,13.1627,13.1627,13.1627,13.1627,13.1635,13.1635,13.1635,13.1627,13.1627,13.1635,13.1635,13.1635,13.1635,13.164,13.1635,13.1635,13.1635,13.164,13.164,13.164,13.1635,13.1635,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.1646,13.164,13.164,13.164,13.164,13.1607,13.1513,13.1407,13.1307,13.1232,13.1152,13.1093,13.1027,13.0921,13.0899,13.0883,13.0813,13.0806,13.0788,13.0774,13.0754,13.0738,13.0732,13.0727,13.0727,13.0726,13.073,13.0731,13.0735,13.0735,13.0734,13.0701,13.068,13.0657,13.0648,13.0636,13.0635,13.0617,13.0559,13.0527,13.0521,13.0429,13.0396,13.0346,13.0338,13.0213,13.0205,12.9996,12.9987,12.9871,12.9862,12.9796,12.9788,12.9746,12.9738,12.9712,12.9704,12.9662,12.9654,12.9496,12.9488,12.9471,12.9429,12.9404,12.9371,12.9354,12.9338,12.9287,12.9279,12.9271,12.9221,12.9212,12.9162,12.9154,12.9129,12.9121,12.9096,12.9087,12.9062,12.9054,12.9004,12.8996,12.8971,12.8962,12.8871,12.8862,12.8821,12.8754,12.8712,12.8646,12.8579,12.8546,12.8455,12.8446,12.8413,12.8379,12.8329,12.8296,12.8196,12.8188,12.8137,12.8129,12.8096,12.8088,12.8005,12.7997,12.7996,12.7927,12.7921,12.7912,12.7871,12.7862,12.7837,12.7829,12.7779,12.777,12.7763,12.7755,12.7721,12.7713,12.7696,12.7671,12.7654,12.7642,12.7642,12.7654,12.7688,12.7696,12.7738,12.7746,12.7762,12.7804,12.7829,12.7838,12.7896,12.7904,12.7938,12.7946,12.7996,12.8004,12.8054,12.8062,12.8129,12.8138,12.8221,12.8271,12.8296,12.8363,12.8371,12.8388,12.8387,12.8371,12.8363,12.8354,12.8313,12.8287,12.8263,12.8254,12.8121,12.8088,12.8021,12.7992,12.7992,12.8009,12.8009,12.805,12.805,12.8058,12.8059,12.8067,12.8116,12.8142,12.8142,12.8129,12.8117,12.8084,12.8066,12.8067,12.805,12.805,12.8033,12.8025,12.8009,12.8,12.8,12.7992,12.7992,12.7975,12.7954,12.7912,12.7892,12.7892,12.7879,12.7859,12.783,12.7821,12.7804,12.7771,12.7758,12.7754,12.7746,12.7742,12.7729,12.7629,12.7621,12.7588,12.7571,12.7537,12.7504,12.7487,12.7471,12.7429,12.7421,12.7388,12.7379,12.7362,12.7296,12.7287,12.7154,12.7096,12.7079,12.7046,12.7034,12.7034,12.7046,12.7079,12.7096,12.71,12.7142,12.7142,12.7204,12.7229,12.7238,12.7254,12.7279,12.7316,12.7325,12.735,12.7387,12.7404,12.7413,12.7458,12.7409,12.7409,12.7425,12.7425,12.74,12.7392,12.7396,12.7404,12.7437,12.7471,12.7492,12.7467,12.7458,12.7387,12.7354,12.7346,12.7338,12.7296,12.7279,12.7275,12.7292,12.7291,12.7309,12.7358,12.7392,12.7391,12.7338,12.7329,12.7313,12.7279,12.7271,12.7246,12.7237,12.7233,12.7221,12.7187,12.7154,12.715,12.7167,12.7175,12.7146,12.7104,12.7095,12.7021,12.6987,12.6971,12.6946,12.6929,12.6904,12.6888,12.6829,12.677,12.6754,12.672,12.6712,12.6646,12.6638,12.6612,12.6571,12.6529,12.6509,12.6508,12.6517,12.655,12.6558,12.6584,12.6583,12.6591,12.6579,12.6538,12.6508,12.6525,12.6513,12.6475,12.6479,12.6492,12.6492,12.6471,12.6454,12.6446,12.6379,12.6363,12.6329,12.6321,12.6313,12.6287,12.6279,12.6246,12.6238,12.6187,12.6179,12.6121,12.6112,12.5964,12.5938,12.5917,12.592,12.5954,12.5979,12.5979,12.5937,12.5913,12.5895,12.5879,12.5854,12.5771,12.5764,12.5763,12.5746,12.5729,12.5712,12.5667,12.5642,12.5634,12.5633,12.5675,12.5687,12.5704,12.5738,12.5754,12.58,12.5796,12.5754,12.5734,12.5733,12.5792,12.5783,12.5791,12.5792,12.5825,12.5825,12.5833,12.5842,12.585,12.585,12.5867,12.5884,12.5884,12.5904,12.5912,12.5925,12.5925,12.5933,12.5934,12.5938,12.5954,12.5996,12.6004,12.6045,12.6054,12.6121,12.6125,12.6113,12.6084,12.6075,12.6075,12.6067,12.6042,12.5975,12.5975,12.6033,12.6042,12.6054,12.6079,12.6096,12.6113,12.6142,12.6142,12.6133,12.6134,12.6137,12.6158,12.6167,12.6167,12.6129,12.6109,12.6109,12.6125,12.6125,12.6183,12.6217,12.6217,12.6267,12.6267,12.6321,12.6329,12.6358,12.6358,12.6371,12.6375,12.6358,12.6358,12.637,12.6446,12.6454,12.6496,12.6546,12.6571,12.6579,12.6621,12.6679,12.6713,12.6721,12.6821,12.6829,12.6938,12.6969,12.6988,12.6993,12.6976,12.6937,12.693,12.6904,12.6887,12.6862,12.685,12.685,12.6846,12.6829,12.6825,12.6812,12.6696,12.6687,12.6646,12.6583,12.6575,12.6575,12.6567,12.6559,12.655,12.655,12.6559,12.6567,12.6621,12.6629,12.6654,12.6675,12.6642,12.6625,12.6625,12.6634,12.6654,12.6662,12.6696,12.6738,12.6796,12.6809,12.6821,12.6829,12.6892,12.69,12.6917,12.6925,12.697,12.7012,12.7033,12.6971,12.6963,12.6929,12.6921,12.6879,12.6854,12.6838,12.6834,12.6863,12.6913,12.6917,12.6871,12.6845,12.6842,12.6821,12.6763,12.6691,12.6684,12.665,12.6641,12.6583,12.6584,12.6571,12.655,12.6467,12.6416,12.6417,12.6392,12.6392,12.6354,12.6346,12.6275,12.6271,12.6262,12.6213,12.6204,12.6166,12.6158,12.6171,12.6213,12.6221,12.6271,12.6304,12.6395,12.643,12.6446,12.6488,12.6496,12.6516,12.6533,12.6563,12.6613,12.6629,12.6754,12.6762,12.6779,12.6809,12.6812,12.6838,12.6896,12.6921,12.6938,12.6954,12.7046,12.7059,12.7059,12.7079,12.7096,12.7121,12.715,12.7158,12.7137,12.7113,12.7092,12.7087,12.7046,12.7038,12.7016,12.7025,12.7025,12.6995,12.6987,12.6942,12.6942,12.6959,12.6954,12.6937,12.6921,12.6896,12.6867,12.6866,12.6858,12.6837,12.6813,12.6758,12.6734,12.6725,12.6725,12.6705,12.6671,12.6658,12.6658,12.6684,12.6709,12.6696,12.6679,12.6638,12.6596,12.6567,12.6567,12.6587,12.6613,12.6634,12.6634,12.6621,12.6588,12.658,12.6529,12.6521,12.6509,12.6508,12.6546,12.6604,12.6613,12.6629,12.6675,12.6675,12.6663,12.6637,12.6629,12.6613,12.6596,12.6546,12.6504,12.6438,12.6417,12.64,12.6379,12.6321,12.63,12.63,12.6263,12.6221,12.62,12.6167,12.6167,12.6134,12.6125,12.6109,12.6109,12.6121,12.6171,12.6188,12.6204,12.6212,12.6234,12.6204,12.6104,12.6096,12.6067,12.6067,12.6042,12.6042,12.6013,12.5996,12.5987,12.5971,12.5942,12.5942,12.5987,12.6012,12.6017,12.6,12.5987,12.5975,12.5975,12.5963,12.5912,12.5855,12.5821,12.58,12.58,12.5792,12.5792,12.5784,12.5783,12.575,12.575,12.5675,12.5667,12.5667,12.5617,12.5617,12.5608,12.56,12.5592,12.5591,12.5617,12.5671,12.5738,12.5796,12.5809,12.5817,12.5837,12.5913,12.5921,12.5962,12.5996,12.6012,12.6046,12.6046,12.6029,12.6004,12.5983,12.5984,12.6,12.6,12.6042,12.6058,12.6067,12.6067,12.6075,12.6075,12.6067,12.6067,12.6075,12.6075,12.6058,12.605,12.6025,12.6025,12.5992,12.5992,12.595,12.5933,12.5925,12.5925,12.5983,12.5979,12.5971,12.5946,12.5912,12.5892,12.5884,12.59,12.59,12.5909,12.5917,12.5942,12.5942,12.5996,12.6004,12.6029,12.6038,12.6083,12.6117,12.6125,12.6108,12.6104,12.6096,12.6067,12.6033,12.595,12.5933,12.5933,12.5908,12.5892,12.5892,12.5867,12.5842,12.585,12.5859,12.5859,12.5875,12.5875,12.5859,12.5858,12.5875,12.5883,12.5892,12.5904,12.5938,12.5954,12.5975,12.5983,12.5992,12.5992,12.5984,12.5883,12.5875,12.5867,12.5867,12.5858,12.5859,12.5833,12.5834,12.585,12.585,12.5771,12.5746,12.5679,12.5662,12.5654,12.5621,12.5612,12.5596,12.5546,12.5483,12.5459,12.545,12.5442,12.5442,12.545,12.545,12.5442,12.5442,12.5459,12.5459,12.5471,12.5475,12.5487,12.5513,12.5529,12.555,12.555,12.5567,12.5584,12.5583,12.5609,12.5616,12.5625,12.5625,12.5592,12.5592,12.5625,12.5633,12.5659,12.5659,12.5684,12.5717,12.5717,12.5725,12.5725,12.5704,12.5679,12.5638,12.5625,12.5646,12.5654,12.5671,12.5721,12.5746,12.5762,12.5767,12.575,12.575,12.5771,12.5788,12.5792,12.5809,12.5809,12.5833,12.5833,12.5846,12.5862,12.5871,12.59,12.59,12.5892,12.5892,12.59,12.59,12.5909,12.5908,12.5917,12.5917,12.5925,12.5925,12.5934,12.5942,12.5966,12.5967,12.5959,12.5958,12.5942,12.5942,12.5934,12.5934,12.59,12.5875,12.5862,12.5829,12.5821,12.5771,12.5758,12.575,12.5746,12.5734,12.5734,12.5704,12.5696,12.5671,12.5662,12.5562,12.5517,12.5517,12.5525,12.5525,12.5558,12.5559,12.555,12.5546,12.5487,12.5454,12.5434,12.5425,12.5425,12.545,12.545,12.5404,12.5383,12.5404,12.5421,12.5429,12.5446,12.5462,12.5479,12.5492,12.55,12.5525,12.5517,12.5521,12.5563,12.5571,12.5587,12.5604,12.5621,12.5629,12.5663,12.5687,12.5704,12.5754,12.5913,12.5929,12.5938,12.5942,12.5963,12.6004,12.6029,12.6054,12.6096,12.6104,12.6225,12.6225,12.6284,12.63,12.63,12.6279,12.6267,12.6267,12.6279,12.6333,12.6342,12.6358,12.6358,12.6367,12.6367,12.635,12.635,12.6304,12.6288,12.6262,12.6246,12.6209,12.62,12.62,12.6129,12.6113,12.6079,12.6071,12.6037,12.6025,12.6017,12.5987,12.5921,12.588,12.5854,12.5845,12.5833,12.5788,12.5779,12.5713,12.5704,12.5687,12.5679,12.567,12.5659,12.565,12.565,12.5667,12.5667,12.5617,12.5617,12.5608,12.5591,12.5546,12.5538,12.5504,12.5479,12.5454,12.5438,12.5405,12.5371,12.5363,12.5237,12.5229,12.5162,12.5096,12.5054,12.4996,12.4996,12.4963,12.495,12.495,12.4938,12.4925,12.4913,12.4854,12.4846,12.4813,12.4787,12.4784,12.4771,12.4767,12.4792,12.4792,12.4733,12.4733,12.4712,12.4683,12.4642,12.4642,12.4629,12.4613,12.4596,12.4579,12.4563,12.4546,12.4529,12.4512,12.4504,12.4492,12.4542,12.4542,12.4654,12.4662,12.4683,12.4683,12.4625,12.4629,12.4671,12.4679,12.4696,12.4725,12.4725,12.4742,12.4767,12.4796,12.4821,12.4871,12.4904,12.4913,12.4962,12.4979,12.5,12.5,12.4975,12.4975,12.4979,12.5,12.5008,12.5038,12.5042,12.5025,12.5009,12.5008,12.5025,12.5058,12.5058,12.5092,12.5092,12.5137,12.5162,12.5171,12.5196,12.5221,12.5229,12.5271,12.5287,12.5295,12.5312,12.5345,12.5354,12.5379,12.5412,12.5454,12.5471,12.5484,12.5484,12.5467,12.5467,12.5434,12.5417,12.5354,12.5346,12.5303,12.5238,12.5188,12.5171,12.5154,12.5104,12.5096,12.5054,12.5046,12.5021,12.4987,12.4971,12.4887,12.4846,12.4837,12.4813,12.4804,12.4737,12.4729,12.4654,12.4646,12.4571,12.4562,12.452,12.4513,12.4487,12.4479,12.4454,12.443,12.4421,12.4404,12.4371,12.4362,12.4296,12.4287,12.4263,12.4246,12.4229,12.422,12.4163,12.4088,12.4062,12.405,12.405,12.4021,12.4004,12.3979,12.3967,12.3967,12.3929,12.3909,12.3871,12.3854,12.382,12.3809,12.3808,12.38,12.3712,12.3663,12.3562,12.3458,12.3395,12.3391,12.3388,12.3354,12.3352,12.3595,12.3595,12.3612,12.3616,12.3654,12.3659,12.3675,12.3679,12.3662,12.3638,12.3608,12.36,12.3595,12.3595,12.3582,12.3571,12.356,12.3543,12.3515,12.3513,12.3515,12.3518,12.3529,12.3532,12.354,12.3552,12.3568,12.3577,12.3579,12.3571,12.3543,12.3527,12.3507,12.3488,12.3471,12.3474,12.349,12.3507,12.352,12.3524,12.354,12.3549,12.3565,12.3582,12.3599,12.3607,12.3629,12.364,12.3649,12.3654,12.3649,12.3638,12.3621,12.3602,12.3582,12.3571,12.3577,12.3599,12.364,12.3688,12.3743,12.3827,12.3915,12.3999,12.409,12.4174,12.4265,12.4349,12.4438,12.4487,12.4521,12.4604,12.4574,12.4513,12.4452,12.4468,12.4496,12.4518,12.4549,12.4579,12.4593,12.4618,12.464,12.4663,12.4679,12.4688,12.4699,12.471,12.4718,12.4721,12.4724,12.4721,12.4715,12.4713,12.4702,12.469,12.4677,12.4668,12.4649,12.4629,12.4618,12.4602,12.4574,12.4549,12.4524,12.4513,12.4511,12.4493,12.4474,12.4457,12.4438,12.4418,12.4393,12.4374,12.4357,12.4365,12.4393,12.4421,12.4449,12.4477,12.4479,12.4468,12.4454,12.4435,12.4421,12.4402,12.4388,12.4377,12.4357,12.4343,12.4324,12.431,12.4296,12.4304,12.4388,12.4463,12.454,12.4624,12.4699,12.4782,12.486,12.4935,12.4996,12.4996,12.5018,12.5093,12.5145,12.5177,12.5254,12.5329,12.5413,12.5488,12.551,12.5565,12.5674,12.5827,12.5977,12.5992,12.6121,12.6154,12.6202,12.6252,12.6293,12.6407,12.6552,12.6696,12.6843,12.6838,12.6838,12.6838,12.6832,12.6832,12.6832,12.6832,12.6832,12.6832,12.6827,12.6827,12.6827,12.6827,12.6827,12.6826,12.6824,12.6824,12.6824,12.6824,12.6824,12.6824,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.681,12.681,12.681,12.681,12.681,12.681,12.681,12.681,12.6802,12.6802,12.6802,12.6802,12.6802,12.6802,12.6802,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6782,12.6777,12.6777,12.6774,12.6774,12.6774,12.6774,12.6774,12.6774,12.6765,12.6765,12.6765,12.6765,12.6763,12.6763,12.6763,12.6763,12.6758,12.6757,12.6757,12.6763,12.6763,12.6763,12.6763,12.6771,12.676,12.6765,12.6765,12.6765,12.6763,12.6763,12.6763,12.6768,12.6768,12.6768,12.6768,12.6765,12.6765,12.6765,12.6771,12.677,12.6765,12.676,12.6752,12.6743,12.6736,12.6735,12.6727,12.6721,12.6714,12.6713,12.6704,12.6696,12.669,12.6688,12.6682,12.6674,12.6665,12.6657,12.6651,12.6649,12.6615,12.6549,12.6524,12.6496,12.6445,12.6427,12.6393,12.6407,12.6493,12.6499,12.6507,12.6529,12.6546,12.6521,12.6493,12.6515,12.6538,12.6568,12.659,12.6585,12.6538,12.649,12.6457,12.6452,12.6454,12.6435,12.6424,12.6446,12.6482,12.6515,12.6518,12.6502,12.6465,12.644,12.6399,12.6368,12.634,12.6343,12.6385,12.6413,12.6429,12.6432,12.6435,12.6435,12.6402,12.6363,12.6396,12.6432,12.6474,12.6435,12.6427,12.6402,12.6382,12.6371,12.6365,12.6346,12.6335,12.6329,12.6338,12.6354,12.6363,12.6385,12.6393,12.6354,12.6321,12.6279,12.624,12.6177,12.6115,12.6054,12.6007,12.5952,12.5896,12.584,12.5785,12.5743,12.5696,12.5652,12.5604,12.5554,12.5493,12.5443,12.5407,12.5371,12.5338,12.5296,12.5282,12.5246,12.5199,12.5143,12.5077,12.5021,12.4997,12.4974,12.4927,12.4885,12.4846,12.4818,12.4793,12.4768,12.4727,12.4729,12.4743,12.4754,12.4763,12.4765,12.4761,12.476,12.4746,12.4729,12.4702,12.4682,12.4663,12.4679,12.4702,12.4732,12.476,12.4807,12.4871,12.4927,12.4988,12.4996,12.5043,12.5107,12.5163,12.5218,12.5271,12.5315,12.5352,12.5379,12.5407,12.5429,12.5452,12.5463,12.5463,12.5465,12.5446,12.5427,12.5393,12.5354,12.5313,12.5274,12.524,12.5213,12.5188,12.5168,12.5157,12.5138,12.5118,12.5107,12.5102,12.5068,12.5018,12.4997,12.4996,12.4985,12.4963,12.4927,12.4871,12.4835,12.4796,12.4782,12.4771,12.4738,12.4699,12.4671,12.466,12.4646,12.4649,12.4652,12.4646,12.4635,12.4629,12.4618,12.4604,12.4543,12.4488,12.4432,12.4363,12.4343,12.4346,12.4354,12.4357,12.436,12.4363,12.4365,12.4374,12.4377,12.4377,12.4371,12.4374,12.4371,12.4371,12.4374,12.4368,12.4371,12.4379,12.4382,12.4371,12.4338,12.4304,12.4263,12.4224,12.4168,12.4121,12.4074,12.4027,12.3985,12.396,12.3954,12.394,12.391,12.3854,12.3804,12.3749,12.3715,12.3677,12.3632,12.3593,12.3588,12.3618,12.3646,12.3677,12.3699,12.3721,12.3735,12.3757,12.3779,12.381,12.3829,12.3852,12.3874,12.3896,12.3932,12.3949,12.3955,12.3963,12.3965,12.396,12.3949,12.3935,12.3915,12.3904,12.3893,12.3879,12.3868,12.3849,12.3824,12.379,12.3771,12.3743,12.374,12.3735,12.3721,12.371,12.3682,12.3649,12.3615,12.3577,12.3529,12.3482,12.3435,12.3385,12.3338,12.3299,12.3265,12.3238,12.3213,12.3171,12.3124,12.3077,12.3035,12.3013,12.3024,12.3046,12.3068,12.309,12.3118,12.314,12.3163,12.3185,12.3199,12.3216,12.3221,12.3243,12.3265,12.3282,12.3296,12.3313,12.3327,12.3349,12.3379,12.3402,12.3424,12.3432,12.3435,12.3435,12.3438,12.344,12.3449,12.3443,12.3452,12.3468,12.3482,12.3499,12.3518,12.3535,12.3557,12.3579,12.3602,12.3629,12.366,12.3688,12.3724,12.3752,12.3788,12.3824,12.3865,12.3902,12.3943,12.3977,12.4021,12.4049,12.4085,12.4099,12.4102,12.411,12.4113,12.4115,12.4115,12.411,12.4099,12.4079,12.4054,12.4027,12.4002,12.3988,12.3982,12.3979,12.3952,12.3932,12.3913,12.391,12.3902,12.3904,12.3899,12.3907,12.391,12.3924,12.3952,12.3982,12.401,12.4032,12.406,12.409,12.4113,12.414,12.4163,12.419,12.4221,12.4249,12.4285,12.4271,12.4224,12.4177,12.4129,12.4088,12.4068,12.4057,12.4043,12.4032,12.4013,12.3985,12.396,12.3932,12.3907,12.3888,12.386,12.3871,12.3885,12.3902,12.3915,12.3932,12.3932,12.3935,12.3929,12.3918,12.3899,12.3865,12.3846,12.3818,12.3799,12.3793,12.3796,12.3799,12.3802,12.3815,12.3832,12.386,12.3888,12.391,12.3918,12.3907,12.3902,12.3924,12.3952,12.3979,12.4002,12.4032,12.406,12.4088,12.4118,12.4146,12.4174,12.4204,12.4238,12.426,12.4263,12.4235,12.421,12.4182,12.419,12.4199,12.4221,12.4238,12.4252,12.426,12.4268,12.4277,12.4288,12.4296,12.4296,12.4304,12.4307,12.4315,12.4324,12.434,12.4343,12.4349,12.4351,12.4352,12.436,12.4363,12.4357,12.4352,12.4346,12.4335,12.4321,12.431,12.4296,12.4285,12.4274,12.426,12.424,12.4215,12.4179,12.4146,12.4121,12.4102,12.409,12.4071,12.4057,12.4038,12.4046,12.4074,12.4102,12.4163,12.4227,12.4288,12.4357,12.4418,12.4474,12.4515,12.4552,12.4593,12.4629,12.4752,12.4763,12.4815,12.4856,12.4877,12.4938,12.4985,12.4996,12.4997,12.5013,12.501,12.4997,12.4949,12.4921,12.4888,12.4843,12.4829,12.4877,12.494,12.4996,12.5007,12.5054,12.509,12.5138,12.5199,12.5254,12.5293,12.5321,12.5346,12.5318,12.5263,12.5207,12.5154,12.5118,12.5124,12.5157,12.5199,12.5238,12.5293,12.5354,12.5402,12.5429,12.5435,12.546,12.5502,12.5538,12.5557,12.5565,12.5582,12.5604,12.5627,12.5654,12.569,12.5752,12.5813,12.5863,12.5902,12.5943,12.599,12.6029,12.6093,12.614,12.619,12.6232,12.6237,12.6288,12.6329,12.6377,12.644,12.6502,12.6557,12.6596,12.6615,12.6615,12.6618,12.6646,12.6674,12.6715,12.6757,12.6804,12.686,12.6896,12.691,12.6913,12.6949,12.701,12.7071,12.7135,12.719,12.7204,12.7193,12.7188,12.7188,12.7204,12.7213,12.7235,12.7263,12.7318,12.7352,12.7407,12.7468,12.7524,12.7557,12.7582,12.761,12.7665,12.7688,12.7721,12.7771,12.7827,12.7893,12.7949,12.8013,12.8079,12.8143,12.8204,12.826,12.8315,12.8313,12.829,12.8282,12.831,12.8363,12.8418,12.8446,12.8488,12.8538,12.8565,12.8613,12.8668,12.8724,12.8771,12.8802,12.8843,12.8885,12.8932,12.9002,12.9049,12.9082,12.9118,12.9157,12.9221,12.9249,12.9257,12.9238,12.921,12.9199,12.9207,12.9229,12.9257,12.9285,12.9349,12.941,12.9471,12.9532,12.9596,12.9657,12.9718,12.9774,12.9829,12.9877,12.991,12.9929,12.9879,12.9818,12.9749,12.9693,12.9646,12.9602,12.9582,12.9607,12.9627,12.9613,12.959,12.9624,12.9671,12.9715,12.9768,12.9832,12.9871,12.9927,12.9977,12.9996,12.9997,13.0003,13.0004,13.0032,13.0054,13.0096,13.0143,13.0163,13.0177,13.0215,13.0271,13.0332,13.0352,13.0388,13.0402,13.0413,13.0479,13.0527,13.0582,13.0624,13.066,13.0702,13.0757,13.0799,13.0824,13.0829,13.0829,13.0849,13.086,13.0879,13.0913,13.0954,13.0988,13.1013,13.1054,13.1088,13.1113,13.1113,13.1113,13.1082,13.1096,13.1138,13.1199,13.1252,13.1302,13.134,13.1374,13.1421,13.1477,13.1527,13.156,13.161,13.1665,13.1699,13.174,13.181,13.1852,13.1877,13.1882,13.1888,13.1913,13.194,13.1982,13.2035,13.209,13.2152,13.2207,13.2257,13.2299,13.234,13.2382,13.2424,13.2468,13.2502,13.2557,13.2604,13.2646,13.2679,13.2727,13.2777,13.2824,13.2863,13.2882,13.2896,13.2929,13.2968,13.301,13.3043,13.3071,13.3104,13.3165,13.3227,13.3282,13.3321,13.3363,13.3418,13.3474,13.3521,13.359,13.3643,13.3657,13.3615,13.3579,13.3565,13.3577,13.359,13.3615,13.3663,13.3718,13.3782,13.3849,13.3899,13.3924,13.3929,13.3907,13.3871,13.3838,13.384,13.3846,13.3832,13.3824,13.3829,13.3849,13.3882,13.3929,13.3977,13.4024,13.4057,13.4099,13.4118,13.4124,13.4113,13.4113,13.4104,13.4068,13.4013,13.3952,13.3902,13.3852,13.381,13.3768,13.3727,13.3677,13.3624,13.356,13.3499,13.3438,13.3407,13.3393,13.3399,13.3377,13.334,13.3299,13.3249,13.3202,13.3157,13.3124,13.3093,13.3065,13.3057,13.3077,13.311,13.3152,13.3204,13.3252,13.3307,13.3377,13.3438,13.3493,13.3527,13.354,13.3543,13.3557,13.3577,13.3602,13.3643,13.3682,13.3724,13.3752,13.3749,13.3746,13.3788,13.3822,13.3857,13.3904,13.3954,13.4004,13.4052,13.4113,13.4177,13.4224,13.426,13.4296,13.4338,13.4393,13.4449,13.4504,13.4571,13.4627,13.4674,13.4715,13.474,13.4746,13.4743,13.4757,13.4782,13.4815,13.4871,13.4927,13.4979,13.4996,13.5021,13.5049,13.5065,13.5102,13.5135,13.516,13.5202,13.5235,13.5254,13.5252,13.5249,13.5268,13.5296,13.5335,13.5371,13.5377,13.5418,13.5457,13.5499,13.5532,13.556,13.5599,13.564,13.5679,13.5683,13.5721,13.5774,13.5829,13.5885,13.5927,13.5965,13.5993,13.6018,13.6046,13.6071,13.6104,13.616,13.6213,13.6268,13.6327,13.6379,13.6427,13.6468,13.6502,13.6535,13.6568,13.6602,13.6643,13.669,13.6738,13.6785,13.6849,13.6904,13.6952,13.6985,13.701,13.7029,13.7043,13.7046,13.7077,13.7077,13.7079,13.709,13.7113,13.7143,13.7165,13.7188,13.7215,13.7238,13.726,13.7282,13.731,13.7338,13.736,13.739,13.7413,13.744,13.7463,13.749,13.7521,13.7532,13.7549,13.7577,13.7613,13.764,13.7677,13.7713,13.774,13.7777,13.7813,13.7849,13.789,13.7932,13.7979,13.8029,13.8079,13.8132,13.8188,13.8243,13.8299,13.8363,13.8418,13.8474,13.8529,13.8585,13.864,13.8702,13.8771,13.8832,13.8902,13.8971,13.904,13.9088,13.9127,13.916,13.9188,13.9221,13.9254,13.9288,13.9313,13.9349,13.9374,13.9399,13.9432,13.946,13.9493,13.9527,13.9552,13.9568,13.9602,13.964,13.9696,13.9752,13.9821,13.9888,13.9957,13.9996,14.0018,14.0068,14.0121,14.0185,14.024,14.0282,14.0324,14.0374,14.044,14.0504,14.0565,14.0613,14.0663,14.0679,14.069,14.0727,14.0763,14.0799,14.0846,14.091,14.0971,14.1018,14.1065,14.1113,14.1177,14.1224,14.1274,14.1324,14.1365,14.1413,14.1477,14.1538,14.1607,14.1668,14.1724,14.1771,14.1818,14.1865,14.1921,14.1977,14.2029,14.2071,14.211,14.2122,14.216,14.2207,14.2254,14.2302,14.2357,14.2404,14.246,14.2513,14.256,14.2615,14.2663,14.2704,14.2743,14.2785,14.2818,14.281,14.2821,14.2849,14.2846,14.2843,14.2857,14.2871,14.2888,14.2915,14.294,14.2968,14.3002,14.304,14.3096,14.3138,14.3177,14.321,14.3243,14.3279,14.3279,14.3268,14.3296,14.336,14.3413,14.3457,14.3485,14.3535,14.3588,14.3629,14.3663,14.3682,14.3671,14.3657,14.3663,14.3704,14.3752,14.379,14.3802,14.3788,14.376,14.3738,14.3763,14.3796,14.3838,14.3871,14.3877,14.3868,14.3888,14.3913,14.3938,14.3965,14.3949,14.3996,14.4004,14.4015,14.4027,14.406,14.4079,14.4143,14.4171,14.421,14.426,14.4307,14.4363,14.4427,14.4479,14.4521,14.4555,14.4568,14.4624,14.4693,14.476,14.4802,14.4829,14.486,14.4902,14.4957,14.4996,14.4997,14.5018,14.5079,14.5124,14.5165,14.5213,14.5268,14.5332,14.5393,14.5449,14.5496,14.5513,14.5507,14.5507,14.5524,14.5552,14.5602,14.5657,14.5713,14.5774,14.5829,14.5877,14.5913,14.5949,14.601,14.6071,14.6115,14.6149,14.6171,14.6182,14.621,14.6246,14.6288,14.6335,14.639,14.646,14.6515,14.6554,14.6582,14.6588,14.6585,14.6596,14.6624,14.6671,14.6732,14.6763,14.6777,14.6813,14.6854,14.6904,14.6965,14.6957,14.6921,14.6921,14.696,14.7021,14.704,14.704,14.7024,14.7029,14.7057,14.7104,14.7152,14.7213,14.7268,14.7324,14.7365,14.7385,14.7402,14.7415,14.7435,14.7468,14.7502,14.7563,14.7582,14.7605,14.7624,14.764,14.7652,14.7677,14.7696,14.7724,14.7749,14.7768,14.7802,14.7829,14.7854,14.7888,14.7921,14.7949,14.7982,14.8015,14.8049,14.8082,14.8121,14.8154,14.8196,14.8238,14.8271,14.8296,14.8315,14.8343,14.836,14.8374,14.8385,14.8385,14.8396]}]],[[{"lng":[-16.5334,-16.5351,-16.5355,-16.5346,-16.5321,-16.5313,-16.5305,-16.5297,-16.5296,-16.5288,-16.5288,-16.5255,-16.5255,-16.5246,-16.5247,-16.523,-16.523,-16.5221,-16.5213,-16.5204,-16.52,-16.518,-16.518,-16.518,-16.5188,-16.5188,-16.5196,-16.5196,-16.5188,-16.5197,-16.5197,-16.523,-16.523,-16.5238,-16.5238,-16.5246,-16.5246,-16.5263,-16.5271,-16.528,-16.528,-16.5288,-16.5296,-16.5299,-16.5305,-16.5313,-16.5334],"lat":[15.7959,15.7959,15.7962,15.8038,15.8096,15.8105,15.817,15.8179,15.8205,15.8213,15.8246,15.8329,15.8371,15.8379,15.8405,15.8429,15.8463,15.8471,15.8587,15.8596,15.8617,15.8604,15.858,15.8504,15.8496,15.8479,15.8471,15.8404,15.8388,15.8379,15.8362,15.8288,15.8254,15.8246,15.8229,15.8221,15.8188,15.8171,15.8129,15.8121,15.8087,15.8079,15.8037,15.8035,15.8029,15.7987,15.7959]},{"lng":[-16.5163,-16.5155,-16.5155,-16.5168,-16.5175,-16.5213,-16.5213,-16.5192,-16.5163],"lat":[15.8662,15.8679,15.8729,15.875,15.875,15.8712,15.8688,15.8659,15.8662]},{"lng":[-16.4872,-16.4871,-16.4879,-16.4883,-16.4896,-16.4896,-16.4897,-16.4884,-16.4872],"lat":[15.9162,15.9213,15.9232,15.9242,15.9237,15.9219,15.9171,15.9158,15.9162]},{"lng":[-16.4959,-16.4946,-16.4946,-16.4946,-16.4946,-16.4955,-16.4946,-16.4955,-16.4955,-16.4955,-16.4955,-16.4963,-16.4963,-16.4946,-16.4951,-16.4967,-16.4975,-16.4996,-16.4996,-16.4971,-16.4988,-16.498,-16.4988,-16.498,-16.4997,-16.498,-16.498,-16.4988,-16.4988,-16.4975,-16.4959],"lat":[15.9359,15.9371,15.9401,15.9437,15.9454,15.9463,15.9488,15.9496,15.9538,15.9552,15.9596,15.9604,15.9638,15.9654,15.9658,15.9658,15.965,15.9629,15.9604,15.9579,15.9554,15.9538,15.9496,15.9479,15.9454,15.9438,15.9396,15.9388,15.9371,15.9359,15.9359]},{"lng":[-16.4988,-16.4988,-16.4988,-16.4988,-16.4963,-16.4964,-16.4976,-16.5001,-16.5021,-16.5022,-16.503,-16.503,-16.503,-16.5021,-16.5021,-16.5009,-16.4988],"lat":[15.9896,15.9911,15.9944,15.9963,15.9996,16.0062,16.0083,16.0083,16.0062,16.002,16.0012,15.9996,15.9979,15.9971,15.9904,15.9884,15.9896]},{"lng":[-16.503,-16.503,-16.503,-16.5022,-16.5021,-16.5005,-16.5013,-16.5035,-16.5064,-16.5064,-16.5064,-16.5063,-16.5043,-16.503],"lat":[16.0204,16.0236,16.0254,16.0263,16.0304,16.0346,16.0387,16.0392,16.0345,16.0233,16.0216,16.0196,16.0192,16.0204]},{"lng":[-16.4772,-16.4772,-16.4755,-16.4755,-16.4738,-16.4739,-16.473,-16.473,-16.4743,-16.4768,-16.4785,-16.4797,-16.4806,-16.483,-16.483,-16.4839,-16.4839,-16.4851,-16.488,-16.4889,-16.4876,-16.4868,-16.4793,-16.4772],"lat":[16.0421,16.0454,16.0471,16.0479,16.0513,16.0571,16.0579,16.0596,16.0609,16.0608,16.06,16.0587,16.0554,16.0529,16.0521,16.0512,16.0479,16.0475,16.0438,16.0388,16.0383,16.0392,16.0392,16.0421]},{"lng":[-16.4968,-16.4975,-16.498,-16.4982,-16.4986,-16.4989,-16.4992,-16.5006,-16.5034,-16.5091,-16.509,-16.509,-16.5092,-16.5098,-16.5107,-16.5105,-16.5113,-16.5113,-16.513,-16.513,-16.5138,-16.5138,-16.5138,-16.5147,-16.5146,-16.5155,-16.5155,-16.5163,-16.5163,-16.5171,-16.5171,-16.518,-16.518,-16.518,-16.5188,-16.5188,-16.5196,-16.5196,-16.5205,-16.5213,-16.52,-16.5188,-16.5171,-16.5163,-16.5163,-16.5155,-16.5155,-16.5138,-16.5138,-16.5146,-16.5146,-16.5155,-16.5146,-16.5138,-16.5138,-16.513,-16.5122,-16.5096,-16.5096,-16.5088,-16.5088,-16.5097,-16.5096,-16.5105,-16.5096,-16.5096,-16.5088,-16.5088,-16.5088,-16.508,-16.5096,-16.5084,-16.5071,-16.5088,-16.5088,-16.5067,-16.5055,-16.5055,-16.508,-16.508,-16.5072,-16.5072,-16.508,-16.508,-16.508,-16.508,-16.508,-16.5072,-16.5072,-16.5056,-16.5055,-16.5064,-16.5064,-16.5047,-16.5047,-16.5047,-16.5047,-16.5043,-16.5031,-16.503,-16.5022,-16.5021,-16.5014,-16.5014,-16.4993,-16.4972,-16.4969,-16.4968,-16.4955,-16.4955,-16.4943,-16.4914,-16.4903,-16.4899,-16.4892,-16.4888,-16.4897,-16.4901,-16.4901,-16.4901,-16.49,-16.4898,-16.4891,-16.489,-16.4891,-16.4891,-16.4891,-16.4893,-16.4896,-16.4896,-16.4895,-16.4895,-16.4896,-16.4897,-16.4899,-16.4901,-16.4902,-16.4904,-16.4905,-16.4908,-16.491,-16.4912,-16.4914,-16.4917,-16.4917,-16.4917,-16.4917,-16.4919,-16.4921,-16.4923,-16.4928,-16.4931,-16.4933,-16.4937,-16.4938,-16.4938,-16.4938,-16.4937,-16.4938,-16.4938,-16.494,-16.4941,-16.4946,-16.4949,-16.4958,-16.4961,-16.4964,-16.4967,-16.4968],"lat":[16.0732,16.0736,16.0737,16.0736,16.0732,16.0727,16.0717,16.0696,16.0664,16.0598,16.0507,16.0415,16.0354,16.0211,16.0122,16.0121,15.9996,15.9988,15.9946,15.9862,15.9854,15.9567,15.9546,15.9538,15.9488,15.948,15.9321,15.9312,15.9237,15.9229,15.9187,15.9179,15.902,15.9004,15.8996,15.8954,15.8946,15.8904,15.8896,15.8779,15.8775,15.8788,15.8879,15.8888,15.8912,15.892,15.8938,15.8971,15.9038,15.9046,15.9088,15.9096,15.9171,15.9179,15.9221,15.9229,15.9287,15.9329,15.9354,15.9362,15.9413,15.9421,15.9463,15.9471,15.9479,15.9546,15.9554,15.9596,15.9729,15.9737,15.9813,15.9816,15.9838,15.9854,15.9871,15.9875,15.9887,15.993,15.9954,16.0013,16.0021,16.0079,16.0087,16.0137,16.0204,16.0223,16.0354,16.0363,16.0379,16.0404,16.0445,16.0454,16.0471,16.0504,16.0527,16.0596,16.0646,16.065,16.0637,16.0595,16.0587,16.0538,16.0529,16.0512,16.0483,16.0496,16.0574,16.0592,16.0571,16.0495,16.0492,16.052,16.0543,16.0555,16.0567,16.0572,16.0577,16.058,16.0581,16.0593,16.0596,16.0599,16.0605,16.0607,16.0613,16.0616,16.0619,16.0622,16.0632,16.0635,16.0639,16.0643,16.0646,16.0648,16.0662,16.0664,16.067,16.067,16.0671,16.0677,16.068,16.068,16.0683,16.0684,16.0686,16.0689,16.0693,16.0693,16.0692,16.0694,16.0697,16.07,16.0703,16.0706,16.0707,16.0708,16.0712,16.0715,16.0718,16.0721,16.0722,16.0723,16.0724,16.0725,16.0728,16.073,16.073,16.073,16.0732]},{"lng":[-14.9459,-14.9503,-14.9528,-14.9545,-14.9556,-14.9564,-14.9573,-14.9575,-14.9584,-14.9609,-14.9631,-14.9689,-14.9739,-14.9803,-14.9861,-14.9903,-14.9945,-14.9981,-14.9992,-15,-15.0025,-15.0031,-15.005,-15.007,-15.0089,-15.0117,-15.0145,-15.017,-15.0206,-15.0248,-15.0306,-15.0353,-15.0411,-15.0475,-15.0534,-15.057,-15.0606,-15.0639,-15.0661,-15.0686,-15.0709,-15.0734,-15.077,-15.0809,-15.0845,-15.0903,-15.097,-15.1017,-15.105,-15.1092,-15.1128,-15.1153,-15.1159,-15.1142,-15.1125,-15.1095,-15.105,-15.1,-15.095,-15.0898,-15.0873,-15.0864,-15.0875,-15.0903,-15.0936,-15.097,-15.1006,-15.1039,-15.1073,-15.1125,-15.1189,-15.1253,-15.1317,-15.1375,-15.1434,-15.1506,-15.157,-15.1636,-15.1684,-15.1734,-15.1798,-15.1848,-15.1903,-15.1959,-15.2017,-15.2075,-15.2123,-15.2173,-15.222,-15.2273,-15.232,-15.2384,-15.245,-15.2506,-15.2559,-15.2617,-15.2667,-15.272,-15.2778,-15.2828,-15.2878,-15.2936,-15.3003,-15.3067,-15.3095,-15.3159,-15.3217,-15.3267,-15.3314,-15.337,-15.342,-15.3478,-15.3542,-15.36,-15.3642,-15.3692,-15.3748,-15.3798,-15.3848,-15.3889,-15.3945,-15.3995,-15.4059,-15.4125,-15.417,-15.4206,-15.4242,-15.4281,-15.4303,-15.4328,-15.4345,-15.4375,-15.4406,-15.445,-15.45,-15.4564,-15.462,-15.4681,-15.4728,-15.4761,-15.4778,-15.4825,-15.4884,-15.4942,-15.4989,-15.4992,-15.5039,-15.5073,-15.5092,-15.5098,-15.5103,-15.5114,-15.5125,-15.5142,-15.5178,-15.5223,-15.5264,-15.532,-15.5367,-15.5425,-15.5489,-15.5561,-15.5628,-15.5689,-15.5756,-15.5814,-15.587,-15.5923,-15.5973,-15.6025,-15.6075,-15.6142,-15.6206,-15.6239,-15.6253,-15.6264,-15.6275,-15.6289,-15.6289,-15.6295,-15.6336,-15.6378,-15.6428,-15.6484,-15.6534,-15.6592,-15.6639,-15.6695,-15.6759,-15.6825,-15.6898,-15.6961,-15.7017,-15.7078,-15.7136,-15.7186,-15.7245,-15.7286,-15.7345,-15.7411,-15.7475,-15.7534,-15.7603,-15.7664,-15.7723,-15.7781,-15.7836,-15.7886,-15.7945,-15.7998,-15.8006,-15.8019,-15.8022,-15.8023,-15.8048,-15.8092,-15.8156,-15.8223,-15.8255,-15.8258,-15.827,-15.832,-15.837,-15.842,-15.8489,-15.855,-15.8592,-15.8614,-15.8616,-15.8617,-15.8636,-15.8675,-15.872,-15.8764,-15.8814,-15.8873,-15.8936,-15.8995,-15.9067,-15.9123,-15.9175,-15.9231,-15.9278,-15.9328,-15.9378,-15.9382,-15.9383,-15.9425,-15.9467,-15.9509,-15.9559,-15.96,-15.9673,-15.9731,-15.9789,-15.9848,-15.9892,-15.9948,-15.9992,-16.0006,-16.0061,-16.0128,-16.0178,-16.0225,-16.0275,-16.0331,-16.0375,-16.0417,-16.0484,-16.0481,-16.0517,-16.0561,-16.0598,-16.0636,-16.0667,-16.0689,-16.0711,-16.0725,-16.0726,-16.0756,-16.0789,-16.0817,-16.0848,-16.0889,-16.0931,-16.0967,-16.0989,-16.0986,-16.0995,-16.102,-16.105,-16.1089,-16.1125,-16.1175,-16.1234,-16.1306,-16.1375,-16.145,-16.1506,-16.1548,-16.1589,-16.1623,-16.1656,-16.1698,-16.1734,-16.1767,-16.1773,-16.1853,-16.1957,-16.208,-16.2223,-16.2406,-16.249,-16.2585,-16.2693,-16.2717,-16.2741,-16.2796,-16.2808,-16.2836,-16.2888,-16.2964,-16.3011,-16.3043,-16.3051,-16.3047,-16.3079,-16.3127,-16.3175,-16.3218,-16.3262,-16.329,-16.3318,-16.3318,-16.3306,-16.3318,-16.3346,-16.3378,-16.3534,-16.3565,-16.3604,-16.3632,-16.3616,-16.3612,-16.3577,-16.3517,-16.3497,-16.3493,-16.3513,-16.3608,-16.3676,-16.3716,-16.3748,-16.3807,-16.3887,-16.3903,-16.3935,-16.4006,-16.4094,-16.4146,-16.4209,-16.4313,-16.4377,-16.4417,-16.4424,-16.4464,-16.4484,-16.4488,-16.448,-16.4484,-16.4496,-16.45,-16.4484,-16.448,-16.4508,-16.4512,-16.4512,-16.448,-16.4444,-16.4436,-16.444,-16.4448,-16.4468,-16.4504,-16.4576,-16.4639,-16.4666,-16.4666,-16.468,-16.4692,-16.468,-16.467,-16.468,-16.468,-16.4651,-16.4635,-16.4618,-16.4593,-16.4585,-16.456,-16.4543,-16.4531,-16.4535,-16.4551,-16.456,-16.4601,-16.4634,-16.4651,-16.466,-16.4676,-16.4697,-16.4697,-16.4705,-16.4705,-16.4722,-16.4755,-16.4764,-16.4739,-16.4739,-16.473,-16.473,-16.4743,-16.4763,-16.4772,-16.481,-16.4868,-16.4876,-16.4901,-16.4918,-16.493,-16.4935,-16.4976,-16.5005,-16.5005,-16.5013,-16.5013,-16.4963,-16.4947,-16.4946,-16.4955,-16.498,-16.4972,-16.498,-16.498,-16.4963,-16.4963,-16.498,-16.498,-16.4967,-16.4958,-16.4921,-16.4913,-16.4913,-16.4892,-16.4859,-16.4834,-16.4817,-16.4796,-16.4788,-16.4788,-16.4775,-16.4771,-16.4788,-16.4788,-16.4775,-16.4742,-16.4734,-16.4717,-16.4693,-16.468,-16.4671,-16.4684,-16.47,-16.4709,-16.4742,-16.4759,-16.4775,-16.4809,-16.4855,-16.4855,-16.483,-16.483,-16.4838,-16.4838,-16.4867,-16.4892,-16.4909,-16.4926,-16.4955,-16.4946,-16.4946,-16.4938,-16.4938,-16.493,-16.4938,-16.4938,-16.4921,-16.4938,-16.4922,-16.4921,-16.4913,-16.4921,-16.4913,-16.4913,-16.4905,-16.4905,-16.4892,-16.4851,-16.4842,-16.4838,-16.4846,-16.4846,-16.4859,-16.4871,-16.4838,-16.4838,-16.4855,-16.4855,-16.4859,-16.4871,-16.4867,-16.4865,-16.4863,-16.4876,-16.4888,-16.4896,-16.4905,-16.4905,-16.4918,-16.4919,-16.4921,-16.4904,-16.4905,-16.4913,-16.4913,-16.4951,-16.4988,-16.4988,-16.5005,-16.5005,-16.5005,-16.4996,-16.4996,-16.5013,-16.5005,-16.5005,-16.4988,-16.4988,-16.4996,-16.4988,-16.5,-16.5026,-16.5038,-16.503,-16.503,-16.5046,-16.5038,-16.5047,-16.503,-16.5046,-16.5038,-16.503,-16.5038,-16.5038,-16.503,-16.5038,-16.5038,-16.5046,-16.5046,-16.5055,-16.5055,-16.5046,-16.5046,-16.5071,-16.5071,-16.5063,-16.5063,-16.5071,-16.5071,-16.5088,-16.5088,-16.508,-16.508,-16.5088,-16.5088,-16.5097,-16.5105,-16.5113,-16.5121,-16.513,-16.5129,-16.5138,-16.5146,-16.5155,-16.5155,-16.5142,-16.5138,-16.5128,-16.5113,-16.5113,-16.5096,-16.5096,-16.508,-16.5071,-16.5071,-16.508,-16.5088,-16.5104,-16.5105,-16.513,-16.513,-16.5155,-16.5163,-16.5171,-16.5171,-16.518,-16.5196,-16.5205,-16.5213,-16.5238,-16.5238,-16.5247,-16.5255,-16.528,-16.5297,-16.533,-16.533,-16.5355,-16.5363,-16.5371,-16.538,-16.5413,-16.543,-16.5438,-16.5471,-16.5471,-16.5488,-16.5488,-16.5505,-16.5505,-16.5538,-16.5546,-16.5563,-16.5563,-16.558,-16.5588,-16.5605,-16.5605,-16.5638,-16.5671,-16.5671,-16.5721,-16.573,-16.5755,-16.5763,-16.5788,-16.5813,-16.588,-16.5888,-16.5946,-16.5946,-16.5955,-16.6022,-16.6071,-16.6096,-16.613,-16.6138,-16.6155,-16.6163,-16.618,-16.6213,-16.6271,-16.6271,-16.633,-16.6355,-16.6438,-16.6463,-16.6513,-16.6513,-16.6538,-16.6546,-16.6571,-16.6571,-16.6605,-16.6613,-16.6655,-16.6655,-16.6721,-16.6721,-16.6755,-16.678,-16.6805,-16.6805,-16.6847,-16.6863,-16.6905,-16.6905,-16.6955,-16.6963,-16.6988,-16.6996,-16.7021,-16.7038,-16.7071,-16.7105,-16.7163,-16.7163,-16.7246,-16.7263,-16.7271,-16.7305,-16.7307,-16.7102,-16.7012,-16.6958,-16.6911,-16.6901,-16.6892,-16.6886,-16.6888,-16.6883,-16.6887,-16.6889,-16.6908,-16.6973,-16.6995,-16.7031,-16.7033,-16.7075,-16.7127,-16.7155,-16.7187,-16.7209,-16.7231,-16.7247,-16.7253,-16.7274,-16.7277,-16.7285,-16.7337,-16.7338,-16.7373,-16.738,-16.7401,-16.7401,-16.7423,-16.7432,-16.7415,-16.7405,-16.7381,-16.7366,-16.7333,-16.7319,-16.7305,-16.7233,-16.718,-16.7167,-16.7155,-16.7136,-16.7094,-16.7081,-16.7073,-16.7048,-16.7022,-16.701,-16.6961,-16.6922,-16.6907,-16.6867,-16.6844,-16.6828,-16.6811,-16.676,-16.6745,-16.673,-16.6701,-16.6687,-16.6669,-16.6635,-16.6626,-16.6608,-16.6557,-16.6536,-16.6476,-16.6428,-16.6423,-16.6412,-16.6384,-16.6365,-16.6338,-16.6312,-16.6307,-16.6306,-16.6285,-16.6245,-16.6207,-16.6165,-16.6133,-16.6121,-16.6121,-16.6123,-16.6129,-16.6134,-16.614,-16.6141,-16.6128,-16.6092,-16.6072,-16.6045,-16.6037,-16.6012,-16.5993,-16.5958,-16.5907,-16.5836,-16.5824,-16.581,-16.5722,-16.5696,-16.5671,-16.5647,-16.5596,-16.5571,-16.554,-16.5519,-16.5471,-16.5451,-16.5418,-16.54,-16.5381,-16.535,-16.5332,-16.5307,-16.5279,-16.5251,-16.5237,-16.5225,-16.5202,-16.5186,-16.5173,-16.5158,-16.5145,-16.5135,-16.5132,-16.5129,-16.5124,-16.5123,-16.5123,-16.5135,-16.5143,-16.5156,-16.5164,-16.5162,-16.5158,-16.5146,-16.5137,-16.5124,-16.5111,-16.5096,-16.5077,-16.506,-16.5044,-16.503,-16.502,-16.5007,-16.499,-16.4974,-16.4954,-16.4941,-16.4915,-16.4901,-16.4889,-16.4873,-16.4861,-16.4845,-16.4832,-16.4817,-16.4804,-16.4791,-16.4779,-16.4759,-16.4751,-16.4736,-16.4709,-16.4691,-16.4659,-16.4641,-16.4618,-16.4594,-16.4564,-16.4546,-16.4538,-16.4526,-16.4515,-16.4509,-16.4483,-16.4458,-16.4436,-16.4414,-16.4385,-16.4362,-16.4338,-16.4311,-16.428,-16.4255,-16.4228,-16.4189,-16.4153,-16.4107,-16.4039,-16.3979,-16.3932,-16.3878,-16.3839,-16.3796,-16.3773,-16.3752,-16.3732,-16.3719,-16.3693,-16.366,-16.363,-16.3587,-16.3547,-16.3515,-16.3496,-16.3468,-16.3448,-16.3421,-16.3411,-16.3393,-16.3371,-16.3356,-16.3346,-16.3344,-16.3354,-16.3371,-16.3388,-16.3416,-16.3446,-16.3454,-16.3457,-16.345,-16.3436,-16.3415,-16.338,-16.335,-16.3332,-16.3305,-16.3274,-16.3237,-16.3202,-16.3171,-16.3142,-16.3112,-16.3102,-16.3093,-16.3088,-16.3084,-16.308,-16.3064,-16.3045,-16.3017,-16.3,-16.2969,-16.295,-16.2931,-16.2908,-16.2885,-16.2868,-16.2858,-16.2849,-16.2832,-16.2806,-16.278,-16.276,-16.2742,-16.2723,-16.2708,-16.27,-16.2692,-16.2685,-16.268,-16.2681,-16.2685,-16.2689,-16.2692,-16.2697,-16.2699,-16.2693,-16.2688,-16.2683,-16.2672,-16.2661,-16.2647,-16.2632,-16.2616,-16.2601,-16.2587,-16.2572,-16.2551,-16.2526,-16.2499,-16.2475,-16.245,-16.2438,-16.2426,-16.2416,-16.2392,-16.2376,-16.2373,-16.2372,-16.2373,-16.2377,-16.2374,-16.2371,-16.2364,-16.2345,-16.2321,-16.2304,-16.227,-16.2249,-16.2188,-16.2136,-16.2116,-16.2086,-16.205,-16.2041,-16.2012,-16.1989,-16.1968,-16.1962,-16.1951,-16.1947,-16.1947,-16.1947,-16.1947,-16.1947,-16.1942,-16.1943,-16.1944,-16.1931,-16.1921,-16.1904,-16.1891,-16.1886,-16.1873,-16.1858,-16.1832,-16.1824,-16.1808,-16.18,-16.1796,-16.1789,-16.1784,-16.1782,-16.1783,-16.1787,-16.1782,-16.1781,-16.1779,-16.1776,-16.1771,-16.1767,-16.1764,-16.1759,-16.1755,-16.1761,-16.1766,-16.1771,-16.1771,-16.1772,-16.1771,-16.177,-16.177,-16.176,-16.1735,-16.1705,-16.1674,-16.1638,-16.1605,-16.1577,-16.1544,-16.1522,-16.1486,-16.147,-16.1422,-16.1395,-16.1368,-16.1345,-16.1317,-16.1284,-16.1235,-16.1213,-16.1177,-16.1158,-16.1114,-16.1095,-16.1074,-16.1049,-16.1025,-16.0992,-16.0971,-16.0948,-16.0938,-16.0921,-16.0905,-16.0879,-16.0863,-16.0843,-16.0816,-16.0794,-16.0759,-16.0728,-16.0703,-16.0681,-16.0658,-16.0602,-16.0577,-16.0545,-16.0521,-16.0495,-16.0469,-16.0439,-16.0408,-16.0377,-16.0361,-16.0348,-16.0321,-16.0264,-16.0176,-16.0156,-16.0139,-16.0099,-16.0078,-16.0058,-16.004,-16.0018,-16.0001,-15.999,-15.9973,-15.9968,-15.9962,-15.9954,-15.9944,-15.9936,-15.9929,-15.9922,-15.9916,-15.9911,-15.9904,-15.9899,-15.9891,-15.9884,-15.9873,-15.9854,-15.9811,-15.9777,-15.9744,-15.972,-15.9688,-15.963,-15.9575,-15.9532,-15.9479,-15.9434,-15.9369,-15.9331,-15.9276,-15.9244,-15.9224,-15.9196,-15.9171,-15.9163,-15.9157,-15.9147,-15.9133,-15.9103,-15.9075,-15.9048,-15.9026,-15.9009,-15.8992,-15.8976,-15.8962,-15.8942,-15.8915,-15.8837,-15.877,-15.8686,-15.8635,-15.8626,-15.8623,-15.8626,-15.8633,-15.8644,-15.8651,-15.8646,-15.862,-15.8588,-15.855,-15.8515,-15.8467,-15.842,-15.8385,-15.8342,-15.8305,-15.8264,-15.8234,-15.8201,-15.8161,-15.8104,-15.8034,-15.8028,-15.8023,-15.8018,-15.8002,-15.7989,-15.7962,-15.793,-15.7892,-15.7868,-15.7851,-15.7829,-15.7815,-15.7798,-15.7779,-15.7765,-15.775,-15.7737,-15.7729,-15.7727,-15.772,-15.7715,-15.7699,-15.7678,-15.7659,-15.7633,-15.7617,-15.7603,-15.7584,-15.7552,-15.7508,-15.7457,-15.7413,-15.7369,-15.7322,-15.7274,-15.7223,-15.7158,-15.7117,-15.7066,-15.7003,-15.6946,-15.687,-15.6791,-15.6722,-15.6664,-15.6628,-15.6554,-15.6498,-15.6435,-15.639,-15.6356,-15.6354,-15.6406,-15.6525,-15.6621,-15.6751,-15.6876,-15.692,-15.6903,-15.6841,-15.6773,-15.676,-15.6755,-15.6752,-15.6749,-15.6744,-15.6739,-15.6735,-15.673,-15.6726,-15.6711,-15.6689,-15.6666,-15.6649,-15.6624,-15.6608,-15.6587,-15.6565,-15.6544,-15.6509,-15.6484,-15.6457,-15.6434,-15.6408,-15.6379,-15.6351,-15.6322,-15.6307,-15.6272,-15.6251,-15.6218,-15.6185,-15.6148,-15.6122,-15.6099,-15.6078,-15.6064,-15.603,-15.6009,-15.5991,-15.5962,-15.5937,-15.5909,-15.5868,-15.5846,-15.5815,-15.5785,-15.5773,-15.5757,-15.5736,-15.5698,-15.5665,-15.563,-15.5597,-15.5562,-15.5523,-15.5496,-15.5452,-15.5432,-15.5408,-15.5386,-15.5364,-15.5348,-15.5313,-15.5276,-15.5238,-15.5204,-15.5175,-15.5148,-15.513,-15.5103,-15.5078,-15.5054,-15.5046,-15.503,-15.501,-15.499,-15.4973,-15.4953,-15.4923,-15.4904,-15.4884,-15.4864,-15.4844,-15.4817,-15.4786,-15.4759,-15.4731,-15.471,-15.4686,-15.4672,-15.4646,-15.4617,-15.4594,-15.4569,-15.4532,-15.4493,-15.4449,-15.4425,-15.4397,-15.4371,-15.4357,-15.4357,-15.4358,-15.4353,-15.4361,-15.436,-15.437,-15.4371,-15.4377,-15.4378,-15.4392,-15.4395,-15.4402,-15.4405,-15.4408,-15.4416,-15.4419,-15.4423,-15.443,-15.4435,-15.4446,-15.4454,-15.4458,-15.4185,-15.3912,-15.3639,-15.3366,-15.3093,-15.2807,-15.252,-15.2234,-15.1947,-15.1661,-15.1375,-15.1088,-15.0887,-15.0855,-15.0802,-15.0515,-15.0229,-15.0162,-14.9989,-14.9896,-14.9629,-14.9362,-14.9095,-14.8891,-14.8686,-14.8496,-14.836,-14.8203,-14.8111,-14.8036,-14.8028,-14.8027,-14.8002,-14.7981,-14.796,-14.7922,-14.7822,-14.7822,-14.7734,-14.7733,-14.7711,-14.7703,-14.7693,-14.7608,-14.7534,-14.7462,-14.7379,-14.7345,-14.7323,-14.7315,-14.7279,-14.7277,-14.7191,-14.7191,-14.7056,-14.7051,-14.6926,-14.6888,-14.683,-14.6642,-14.6634,-14.6628,-14.659,-14.6587,-14.6581,-14.6577,-14.6514,-14.6329,-14.6271,-14.6076,-14.6053,-14.5993,-14.5894,-14.5821,-14.5738,-14.5733,-14.5717,-14.5674,-14.5658,-14.5604,-14.5566,-14.5484,-14.5413,-14.5316,-14.5247,-14.519,-14.513,-14.5025,-14.4985,-14.4976,-14.494,-14.494,-14.4939,-14.4937,-14.4913,-14.4904,-14.4895,-14.4882,-14.4879,-14.4858,-14.4846,-14.4837,-14.4834,-14.4824,-14.4814,-14.4799,-14.4776,-14.4755,-14.4739,-14.4713,-14.4707,-14.4703,-14.4693,-14.4666,-14.456,-14.4472,-14.4356,-14.428,-14.4166,-14.408,-14.3998,-14.3799,-14.3683,-14.3681,-14.3515,-14.3513,-14.3302,-14.3301,-14.3301,-14.33,-14.3286,-14.3236,-14.2941,-14.2904,-14.2846,-14.2775,-14.2414,-14.2414,-14.2385,-14.2384,-14.1875,-14.1851,-14.1825,-14.1792,-14.1792,-14.1642,-14.1565,-14.1564,-14.1545,-14.1062,-14.0565,-14.0531,-14.034,-14.0253,-14.0016,-13.9989,-13.9948,-13.9737,-13.9645,-13.9482,-13.9227,-13.8972,-13.8717,-13.8462,-13.8194,-13.8163,-13.7886,-13.7874,-13.7868,-13.7864,-13.7834,-13.7564,-13.7265,-13.7162,-13.6966,-13.6667,-13.6368,-13.6069,-13.577,-13.5471,-13.5323,-13.5172,-13.5021,-13.4869,-13.4717,-13.4565,-13.4413,-13.4261,-13.4109,-13.3958,-13.3806,-13.3654,-13.3502,-13.335,-13.3198,-13.3046,-13.2894,-13.2742,-13.259,-13.2438,-13.2141,-13.2106,-13.1845,-13.1548,-13.1252,-13.0956,-13.0669,-13.0659,-13.0363,-13.0174,-12.9989,-12.9745,-12.9537,-12.9506,-12.9247,-12.8988,-12.8729,-12.847,-12.8211,-12.7952,-12.7701,-12.7449,-12.7198,-12.703,-12.6861,-12.6693,-12.6524,-12.651,-12.6429,-12.6334,-12.6239,-12.6147,-12.6144,-12.6049,-12.5954,-12.5942,-12.5859,-12.5866,-12.587,-12.5878,-12.5885,-12.5892,-12.5904,-12.591,-12.5914,-12.592,-12.5928,-12.5937,-12.594,-12.5948,-12.5953,-12.596,-12.5967,-12.597,-12.5987,-12.5992,-12.6,-12.6006,-12.6006,-12.6001,-12.6001,-12.6005,-12.6005,-12.6011,-12.6011,-12.6007,-12.6013,-12.6008,-12.6006,-12.6014,-12.6019,-12.602,-12.6022,-12.6022,-12.6026,-12.6031,-12.6032,-12.6038,-12.6042,-12.6042,-12.6049,-12.6053,-12.6051,-12.6058,-12.6066,-12.6074,-12.6087,-12.6098,-12.6107,-12.6118,-12.6128,-12.6129,-12.6126,-12.612,-12.611,-12.6102,-12.6117,-12.6126,-12.614,-12.6144,-12.6156,-12.6155,-12.6152,-12.615,-12.6149,-12.6141,-12.6141,-12.6146,-12.6155,-12.6163,-12.6177,-12.6192,-12.6205,-12.6213,-12.6209,-12.6212,-12.6211,-12.6213,-12.6211,-12.6219,-12.6228,-12.6235,-12.6237,-12.6247,-12.6444,-12.6641,-12.6861,-12.6911,-12.6961,-12.6998,-12.7036,-12.7073,-12.7103,-12.7134,-12.717,-12.72,-12.7236,-12.7281,-12.7317,-12.7361,-12.7403,-12.7448,-12.7498,-12.7542,-12.7592,-12.7636,-12.7686,-12.7803,-12.7811,-12.7825,-12.785,-12.7845,-12.7825,-12.7825,-12.7828,-12.7831,-12.7834,-12.785,-12.7873,-12.7909,-12.7939,-12.7981,-12.8075,-12.815,-12.8209,-12.8259,-12.8314,-12.837,-12.8436,-12.8489,-12.8525,-12.8553,-12.857,-12.8606,-12.8645,-12.8636,-12.8673,-12.8717,-12.8753,-12.8784,-12.882,-12.8856,-12.8886,-12.8903,-12.8925,-12.8898,-12.8856,-12.8792,-12.8734,-12.8686,-12.8628,-12.857,-12.8506,-12.8464,-12.8431,-12.8417,-12.8406,-12.8392,-12.8392,-12.8403,-12.8425,-12.845,-12.8486,-12.8514,-12.855,-12.8581,-12.862,-12.8664,-12.8706,-12.8742,-12.8786,-12.8836,-12.8881,-12.8934,-12.8975,-12.902,-12.9064,-12.9114,-12.9164,-12.9217,-12.9253,-12.9289,-12.9317,-12.9342,-12.9359,-12.9375,-12.9392,-12.9395,-12.9395,-12.9398,-12.9392,-12.9386,-12.9375,-12.9356,-12.9328,-12.9336,-12.9386,-12.9423,-12.9475,-12.9492,-12.9492,-12.9509,-12.9523,-12.9523,-12.9539,-12.9556,-12.955,-12.9553,-12.9548,-12.9561,-12.9566,-12.9595,-12.9639,-12.9709,-12.9764,-12.9814,-12.9864,-12.9928,-12.9984,-12.9992,-13.0034,-13.0084,-13.0148,-13.0211,-13.0275,-13.0325,-13.0373,-13.0417,-13.0467,-13.0523,-13.0592,-13.065,-13.0686,-13.0731,-13.0767,-13.082,-13.0864,-13.0909,-13.0931,-13.0932,-13.0936,-13.0967,-13.0992,-13.1031,-13.1036,-13.1031,-13.1017,-13.0989,-13.0964,-13.0931,-13.0903,-13.0889,-13.0886,-13.0895,-13.0898,-13.0906,-13.0914,-13.0959,-13.1011,-13.1053,-13.1098,-13.1148,-13.1192,-13.1234,-13.1278,-13.1314,-13.1367,-13.1403,-13.1453,-13.1506,-13.1564,-13.1614,-13.165,-13.1686,-13.1731,-13.1781,-13.1839,-13.1903,-13.1967,-13.2023,-13.2075,-13.2131,-13.2181,-13.225,-13.23,-13.2339,-13.2361,-13.2392,-13.2417,-13.2445,-13.2467,-13.2484,-13.25,-13.2503,-13.2489,-13.247,-13.245,-13.2417,-13.2361,-13.2295,-13.2239,-13.2211,-13.2192,-13.2186,-13.2217,-13.2253,-13.2289,-13.2325,-13.2359,-13.2395,-13.2439,-13.2475,-13.2511,-13.2548,-13.2592,-13.2636,-13.2673,-13.2711,-13.2742,-13.2764,-13.2778,-13.2809,-13.2831,-13.2848,-13.2859,-13.2859,-13.2839,-13.2789,-13.2748,-13.2714,-13.2728,-13.2764,-13.2806,-13.2842,-13.2878,-13.2906,-13.2945,-13.2967,-13.2989,-13.3006,-13.3025,-13.3039,-13.3048,-13.305,-13.3053,-13.3064,-13.3078,-13.3103,-13.3131,-13.3161,-13.3198,-13.3198,-13.3242,-13.3259,-13.3253,-13.3248,-13.3248,-13.3245,-13.3239,-13.327,-13.33,-13.3336,-13.3381,-13.3417,-13.347,-13.3511,-13.3556,-13.3592,-13.3625,-13.3648,-13.3678,-13.37,-13.3731,-13.3753,-13.3775,-13.3786,-13.3803,-13.3804,-13.3804,-13.3817,-13.3828,-13.3836,-13.3848,-13.3856,-13.3873,-13.3889,-13.3906,-13.392,-13.3945,-13.3973,-13.4017,-13.4061,-13.4114,-13.4156,-13.4217,-13.4264,-13.4309,-13.4359,-13.4411,-13.4448,-13.4492,-13.4523,-13.4559,-13.4598,-13.462,-13.465,-13.4681,-13.4714,-13.4748,-13.4778,-13.4811,-13.4839,-13.4873,-13.49,-13.4942,-13.4992,-13.5006,-13.505,-13.5081,-13.5111,-13.5134,-13.515,-13.5145,-13.5134,-13.5098,-13.5056,-13.5014,-13.4992,-13.4973,-13.4923,-13.4886,-13.4861,-13.4842,-13.485,-13.4881,-13.4925,-13.4945,-13.4992,-13.5003,-13.5042,-13.5081,-13.5106,-13.5148,-13.5198,-13.522,-13.5284,-13.5348,-13.5411,-13.5478,-13.5534,-13.56,-13.5659,-13.5723,-13.5792,-13.5859,-13.5909,-13.5964,-13.6014,-13.6056,-13.6106,-13.6148,-13.6203,-13.6253,-13.6325,-13.6389,-13.6445,-13.6489,-13.6523,-13.6556,-13.6592,-13.6639,-13.6703,-13.6761,-13.6786,-13.6823,-13.6853,-13.6889,-13.6928,-13.6964,-13.7009,-13.7053,-13.7089,-13.712,-13.712,-13.7109,-13.7089,-13.7061,-13.7036,-13.7009,-13.6981,-13.6979,-13.697,-13.697,-13.6995,-13.7039,-13.7103,-13.7153,-13.72,-13.7242,-13.7278,-13.732,-13.7361,-13.7398,-13.7439,-13.7464,-13.7506,-13.7548,-13.76,-13.7656,-13.7711,-13.7784,-13.7848,-13.7911,-13.7967,-13.8025,-13.8084,-13.8131,-13.8173,-13.8223,-13.8267,-13.8309,-13.8348,-13.8384,-13.8425,-13.8473,-13.8506,-13.8528,-13.855,-13.8573,-13.8617,-13.865,-13.8673,-13.8681,-13.8692,-13.8709,-13.8703,-13.8706,-13.8709,-13.8711,-13.872,-13.8742,-13.8778,-13.8845,-13.8895,-13.8953,-13.9011,-13.9061,-13.9092,-13.9114,-13.9153,-13.9181,-13.9225,-13.927,-13.9306,-13.935,-13.94,-13.9459,-13.95,-13.9531,-13.9595,-13.9653,-13.9689,-13.9725,-13.9728,-13.9717,-13.9698,-13.9684,-13.9664,-13.9639,-13.9639,-13.9661,-13.9706,-13.975,-13.9786,-13.9825,-13.9842,-13.9836,-13.9817,-13.9781,-13.9725,-13.9678,-13.9636,-13.9606,-13.9631,-13.9667,-13.9711,-13.975,-13.9786,-13.9823,-13.9861,-13.9898,-13.9936,-13.9981,-13.9992,-14.0036,-14.0103,-14.0173,-14.0239,-14.0289,-14.0325,-14.0356,-14.0392,-14.0417,-14.0434,-14.0464,-14.0486,-14.0531,-14.0573,-14.062,-14.0664,-14.07,-14.0736,-14.0773,-14.0811,-14.0853,-14.09,-14.0945,-14.0981,-14.1025,-14.1061,-14.11,-14.1114,-14.1148,-14.1175,-14.1211,-14.1256,-14.13,-14.1339,-14.1375,-14.1406,-14.1431,-14.1453,-14.1484,-14.1528,-14.1578,-14.1623,-14.1667,-14.1703,-14.1728,-14.1764,-14.1785,-14.1785,-14.1795,-14.1831,-14.1861,-14.1898,-14.1942,-14.1981,-14.2011,-14.2039,-14.2073,-14.21,-14.2125,-14.2161,-14.2206,-14.2242,-14.2295,-14.2339,-14.2375,-14.2409,-14.2453,-14.2495,-14.2528,-14.2564,-14.2592,-14.2625,-14.2684,-14.2728,-14.2759,-14.2775,-14.277,-14.2773,-14.2767,-14.2764,-14.2764,-14.2775,-14.2825,-14.2889,-14.2945,-14.2998,-14.3039,-14.3086,-14.3136,-14.3195,-14.325,-14.3303,-14.3348,-14.3378,-14.3389,-14.337,-14.3356,-14.3345,-14.3348,-14.3342,-14.3336,-14.3331,-14.3361,-14.34,-14.3445,-14.35,-14.3545,-14.355,-14.3625,-14.3695,-14.3767,-14.3834,-14.3903,-14.3953,-14.3992,-14.3992,-14.4006,-14.4056,-14.41,-14.4153,-14.4189,-14.4236,-14.4292,-14.4342,-14.4392,-14.4425,-14.4453,-14.4478,-14.4498,-14.4534,-14.4584,-14.4656,-14.472,-14.4792,-14.4856,-14.4911,-14.4978,-14.4992,-14.5034,-14.5106,-14.517,-14.5228,-14.5286,-14.5345,-14.5395,-14.5434,-14.5478,-14.5523,-14.5586,-14.5628,-14.5661,-14.5706,-14.5753,-14.5717,-14.5784,-14.5853,-14.5925,-14.5984,-14.6042,-14.6184,-14.6256,-14.6328,-14.6373,-14.6389,-14.6406,-14.6442,-14.6514,-14.6578,-14.6642,-14.6709,-14.6773,-14.6828,-14.6886,-14.6942,-14.7,-14.7064,-14.7123,-14.7186,-14.7239,-14.7281,-14.7306,-14.7342,-14.7398,-14.7448,-14.7473,-14.7495,-14.7525,-14.7556,-14.7595,-14.765,-14.7709,-14.7745,-14.7792,-14.7834,-14.7886,-14.7942,-14.7992,-14.8031,-14.8067,-14.8111,-14.817,-14.8228,-14.8267,-14.8309,-14.8361,-14.8417,-14.8481,-14.8539,-14.8595,-14.861,-14.8639,-14.8678,-14.8714,-14.8753,-14.8795,-14.8839,-14.8884,-14.8917,-14.8945,-14.8978,-14.9014,-14.9061,-14.912,-14.9142,-14.9167,-14.9234,-14.9281,-14.9345,-14.9403,-14.9409,-14.9459],"lat":[16.6424,16.6449,16.6496,16.6552,16.6613,16.6674,16.6721,16.6735,16.6796,16.6843,16.689,16.6902,16.6921,16.6927,16.691,16.6879,16.6852,16.6815,16.6795,16.6779,16.6738,16.669,16.664,16.659,16.6543,16.6499,16.6457,16.6415,16.6379,16.6352,16.6335,16.6313,16.6296,16.6299,16.6313,16.6346,16.6377,16.6418,16.6465,16.6513,16.656,16.6607,16.664,16.6674,16.6704,16.6718,16.6707,16.6685,16.6649,16.6621,16.6585,16.6543,16.6479,16.6424,16.6371,16.6329,16.6304,16.6288,16.6268,16.6249,16.6202,16.614,16.6077,16.6035,16.5999,16.5963,16.5927,16.5893,16.5857,16.5835,16.5824,16.5829,16.5832,16.5846,16.5857,16.5852,16.5857,16.5849,16.5824,16.5802,16.5793,16.5771,16.5754,16.5738,16.5721,16.5707,16.5685,16.566,16.5638,16.5615,16.5593,16.5585,16.5588,16.5599,16.5618,16.5629,16.5649,16.5665,16.5679,16.5696,16.5715,16.5727,16.5732,16.5721,16.5721,16.571,16.5696,16.5671,16.5649,16.5635,16.561,16.561,16.5613,16.5596,16.5582,16.556,16.5543,16.5521,16.5496,16.5468,16.5452,16.5443,16.5432,16.5438,16.5463,16.5496,16.5529,16.5563,16.561,16.5657,16.5713,16.5752,16.579,16.5818,16.5835,16.5827,16.581,16.5793,16.5771,16.5749,16.5749,16.5727,16.571,16.5693,16.5671,16.567,16.5649,16.5613,16.5563,16.5502,16.544,16.5385,16.5327,16.5279,16.5243,16.5213,16.5185,16.5168,16.5146,16.5129,16.5121,16.5118,16.5121,16.5127,16.5129,16.514,16.5152,16.5171,16.519,16.5207,16.5227,16.5229,16.5221,16.5185,16.5129,16.5074,16.5032,16.4997,16.4996,16.4982,16.4954,16.4924,16.4902,16.4885,16.4863,16.4846,16.4824,16.4807,16.4799,16.479,16.4785,16.479,16.4802,16.4813,16.4824,16.4843,16.4854,16.4865,16.4877,16.4882,16.4885,16.4896,16.4893,16.4904,16.4918,16.4929,16.494,16.4957,16.4968,16.4988,16.4991,16.4996,16.4997,16.4998,16.5007,16.5018,16.5024,16.5013,16.4998,16.4996,16.499,16.4968,16.4957,16.4949,16.4946,16.4957,16.4982,16.4996,16.4997,16.4997,16.501,16.504,16.5068,16.5093,16.5113,16.5124,16.5127,16.5124,16.5121,16.5104,16.5082,16.5065,16.5043,16.5021,16.4999,16.4997,16.4996,16.4974,16.4946,16.4915,16.4893,16.4863,16.486,16.4874,16.4885,16.4896,16.4921,16.4932,16.493,16.4929,16.4927,16.4918,16.4893,16.4871,16.4849,16.4832,16.4829,16.4815,16.4804,16.4799,16.4804,16.4829,16.4849,16.4882,16.4921,16.4954,16.4988,16.4996,16.4997,16.5015,16.504,16.5065,16.5093,16.5118,16.5152,16.5185,16.5232,16.5293,16.534,16.5388,16.5429,16.546,16.5493,16.5513,16.5524,16.5521,16.5518,16.5515,16.5499,16.5468,16.544,16.5404,16.5368,16.5338,16.5302,16.5265,16.5259,16.5231,16.5195,16.5183,16.5211,16.5327,16.5334,16.5311,16.5239,16.5205,16.5171,16.5052,16.4936,16.4797,16.4745,16.4698,16.4642,16.4566,16.4491,16.4109,16.4037,16.4021,16.4029,16.4029,16.3997,16.3941,16.375,16.3703,16.3627,16.3512,16.344,16.34,16.3307,16.3289,16.3241,16.3149,16.3014,16.3,16.2891,16.2807,16.2763,16.2723,16.2668,16.2584,16.2473,16.2361,16.227,16.221,16.2166,16.2161,16.215,16.2154,16.2174,16.2182,16.2178,16.2126,16.2079,16.2028,16.2019,16.1943,16.1864,16.1784,16.1673,16.1633,16.1537,16.1494,16.1438,16.1374,16.1287,16.1259,16.1179,16.1135,16.1119,16.1095,16.1028,16.0984,16.094,16.0892,16.0877,16.0845,16.0826,16.0825,16.0804,16.0796,16.0789,16.0776,16.0763,16.0713,16.0684,16.0684,16.07,16.07,16.0692,16.0692,16.0683,16.0671,16.0667,16.0667,16.0659,16.0658,16.0642,16.0642,16.0634,16.0634,16.0604,16.0579,16.0571,16.0471,16.0437,16.0404,16.0371,16.0346,16.0329,16.0321,16.0304,16.0292,16.0304,16.0337,16.0375,16.0375,16.0367,16.0358,16.0334,16.0329,16.0317,16.0308,16.0271,16.0245,16.0237,16.0179,16.0129,16.0104,15.9996,15.997,15.9954,15.9921,15.9913,15.9896,15.9879,15.9787,15.977,15.9763,15.975,15.975,15.9779,15.9796,15.9829,15.985,15.9867,15.9867,15.9875,15.9896,15.9912,15.9937,15.9942,15.9913,15.9896,15.9887,15.9867,15.9858,15.985,15.985,15.9883,15.9879,15.9829,15.9817,15.9817,15.9825,15.9825,15.9817,15.9842,15.9842,15.9813,15.9795,15.9762,15.9721,15.9712,15.9696,15.9666,15.9667,15.965,15.965,15.9629,15.9612,15.9538,15.9529,15.9504,15.9495,15.9487,15.9471,15.9454,15.9437,15.9421,15.9404,15.9396,15.938,15.9371,15.9321,15.9313,15.9296,15.9284,15.9284,15.9292,15.9279,15.927,15.9262,15.9258,15.9238,15.9187,15.9137,15.9105,15.9037,15.9034,15.9046,15.9095,15.911,15.9129,15.9134,15.9121,15.9071,15.9063,15.9046,15.9017,15.9063,15.9137,15.9171,15.9221,15.9229,15.9254,15.9258,15.9296,15.9329,15.9338,15.946,15.9487,15.9496,15.9571,15.9596,15.9605,15.9637,15.9663,15.9746,15.9754,15.9788,15.98,15.9792,15.9779,15.9763,15.9696,15.9671,15.9654,15.9638,15.9621,15.9588,15.9546,15.9538,15.9529,15.9512,15.9505,15.9496,15.9471,15.9462,15.9363,15.9355,15.9271,15.9262,15.9221,15.9179,15.9154,15.9146,15.9062,15.9054,15.9038,15.9021,15.8988,15.8979,15.8912,15.8904,15.8813,15.8804,15.8738,15.8729,15.8596,15.8588,15.8563,15.8554,15.8513,15.8504,15.8479,15.8475,15.8479,15.8502,15.8537,15.8563,15.8587,15.8621,15.8662,15.8662,15.8604,15.8596,15.8562,15.8545,15.8488,15.8463,15.8454,15.8429,15.8346,15.8338,15.8296,15.8287,15.8221,15.8213,15.8163,15.8112,15.8079,15.8071,15.8029,15.7979,15.7962,15.7904,15.7887,15.7863,15.7829,15.782,15.7771,15.7712,15.7696,15.7671,15.7621,15.7604,15.7579,15.7563,15.7537,15.7521,15.7471,15.7437,15.7412,15.7396,15.7371,15.7338,15.7321,15.7304,15.7263,15.7196,15.7179,15.7087,15.7054,15.7013,15.6971,15.6921,15.687,15.6787,15.6763,15.668,15.6662,15.6654,15.6521,15.6421,15.6387,15.6337,15.6313,15.6296,15.627,15.6254,15.6187,15.6129,15.6121,15.6045,15.5988,15.5854,15.5813,15.5754,15.5746,15.5721,15.5696,15.5671,15.5663,15.5621,15.5596,15.5545,15.5537,15.5471,15.5462,15.5412,15.5355,15.5329,15.532,15.5279,15.5246,15.5196,15.5187,15.5121,15.5095,15.5071,15.5046,15.5021,15.4979,15.4946,15.4887,15.4821,15.4813,15.4721,15.4687,15.4663,15.4621,15.4614,15.4525,15.4471,15.443,15.4379,15.4348,15.4333,15.4303,15.4278,15.4259,15.4226,15.422,15.4167,15.4086,15.4068,15.403,15.4021,15.3963,15.3905,15.3866,15.3804,15.3726,15.3693,15.3645,15.364,15.3577,15.3549,15.3523,15.3406,15.3395,15.3313,15.3241,15.3166,15.3152,15.31,15.3026,15.2994,15.2991,15.297,15.2967,15.2948,15.2949,15.2939,15.2917,15.291,15.2904,15.2903,15.2892,15.2835,15.282,15.2796,15.2768,15.2737,15.2722,15.2693,15.2685,15.2667,15.265,15.2629,15.2625,15.2603,15.2561,15.2541,15.2522,15.2486,15.2461,15.2453,15.2426,15.2409,15.2405,15.2367,15.2357,15.2304,15.2283,15.2275,15.2271,15.226,15.2247,15.2226,15.2178,15.2115,15.2086,15.2034,15.2002,15.1985,15.1957,15.1891,15.1824,15.1785,15.1746,15.1692,15.166,15.1621,15.1563,15.1457,15.1496,15.1526,15.1532,15.154,15.1548,15.1573,15.1649,15.1714,15.1783,15.1792,15.1797,15.1886,15.1924,15.1959,15.1987,15.2111,15.2141,15.2224,15.2254,15.2355,15.2389,15.2456,15.2486,15.2524,15.2585,15.2619,15.2662,15.2693,15.2687,15.2686,15.2688,15.2696,15.2703,15.2713,15.2726,15.2743,15.2766,15.2783,15.2795,15.281,15.2828,15.2836,15.2839,15.2844,15.2849,15.2861,15.2868,15.288,15.29,15.2917,15.2935,15.2953,15.2974,15.3001,15.3022,15.3045,15.3066,15.3083,15.31,15.3119,15.3133,15.3147,15.3158,15.3186,15.3208,15.3224,15.3237,15.325,15.3267,15.3281,15.3296,15.3314,15.3332,15.3341,15.3361,15.3378,15.3393,15.3417,15.3434,15.3461,15.3476,15.3499,15.3518,15.3545,15.356,15.3569,15.3582,15.3591,15.3591,15.3577,15.356,15.3536,15.3511,15.3482,15.3455,15.3431,15.3406,15.3368,15.3339,15.3311,15.3276,15.3258,15.3252,15.3257,15.3255,15.3245,15.3233,15.3212,15.3188,15.3171,15.3164,15.3174,15.3192,15.3213,15.3233,15.3239,15.3242,15.3235,15.3213,15.3196,15.3167,15.3136,15.3105,15.3088,15.3056,15.3029,15.2994,15.2957,15.2925,15.2878,15.2847,15.282,15.2793,15.2769,15.2758,15.2743,15.2715,15.2689,15.266,15.2636,15.2609,15.2579,15.2561,15.2537,15.2513,15.25,15.249,15.2488,15.2509,15.2521,15.2548,15.2596,15.2634,15.2662,15.2691,15.2706,15.2723,15.2725,15.2728,15.273,15.2725,15.2722,15.2717,15.2703,15.2696,15.269,15.2672,15.2646,15.2618,15.2596,15.257,15.2536,15.2498,15.2467,15.2425,15.2392,15.2357,15.2313,15.2288,15.2247,15.2232,15.2201,15.2176,15.2149,15.21,15.2049,15.2006,15.1971,15.1941,15.1907,15.1883,15.1855,15.1831,15.1807,15.1771,15.1741,15.1712,15.1692,15.1677,15.1658,15.163,15.1607,15.1563,15.1537,15.15,15.1443,15.1382,15.1345,15.131,15.1278,15.1247,15.1219,15.1194,15.118,15.1162,15.1147,15.1113,15.1079,15.1068,15.1046,15.1014,15.1003,15.0971,15.0954,15.0934,15.0912,15.0863,15.0805,15.0766,15.0711,15.0663,15.0617,15.0561,15.0491,15.0454,15.0417,15.037,15.0337,15.0306,15.0273,15.024,15.0219,15.0182,15.0161,15.0118,15.0115,15.0092,15.0078,15.0058,15.0015,15.0006,14.9953,14.9934,14.9908,14.9887,14.9865,14.9844,14.9827,14.9803,14.9784,14.9767,14.9752,14.9732,14.9712,14.9695,14.9682,14.9662,14.9648,14.9629,14.963,14.9633,14.9632,14.9636,14.9635,14.9626,14.9628,14.9631,14.9631,14.9626,14.9621,14.9617,14.9612,14.9608,14.9609,14.9603,14.9597,14.9585,14.9579,14.9574,14.957,14.9564,14.9564,14.9565,14.9565,14.9566,14.9573,14.9575,14.958,14.9581,14.9585,14.9588,14.9597,14.9602,14.9608,14.9607,14.9611,14.961,14.961,14.9606,14.9605,14.9604,14.9598,14.9599,14.9597,14.9594,14.9593,14.9595,14.959,14.9587,14.959,14.959,14.9589,14.9586,14.9577,14.9568,14.9567,14.9563,14.9559,14.9563,14.9562,14.9555,14.9557,14.9554,14.9553,14.9553,14.9557,14.9564,14.9583,14.9599,14.9618,14.9634,14.9653,14.967,14.9686,14.9701,14.9708,14.9716,14.9721,14.9723,14.9726,14.9721,14.9722,14.9722,14.9722,14.9726,14.9731,14.9736,14.9741,14.9746,14.9749,14.9756,14.9761,14.9766,14.9768,14.9771,14.9776,14.9785,14.9788,14.9786,14.9788,14.9792,14.9802,14.9816,14.9826,14.9835,14.9849,14.9867,14.9887,14.9908,14.9934,14.9964,15.0027,15.0064,15.0109,15.0165,15.02,15.025,15.0334,15.0378,15.0417,15.0455,15.0494,15.0512,15.0511,15.0485,15.0441,15.0389,15.0327,15.0287,15.0234,15.0183,15.0123,15.0071,15.0033,15.0009,14.9993,14.9979,14.9801,14.9688,14.968,14.9661,14.9635,14.9619,14.9611,14.9601,14.9588,14.9568,14.9547,14.9521,14.9497,14.9479,14.9459,14.9438,14.9412,14.9391,14.9358,14.933,14.9306,14.9282,14.9249,14.9222,14.9179,14.915,14.9119,14.9085,14.9042,14.9013,14.8984,14.8977,14.8982,14.8987,14.9,14.9023,14.9062,14.9086,14.9132,14.9208,14.9314,14.9416,14.9459,14.9455,14.9421,14.9305,14.9206,14.9151,14.9102,14.903,14.8965,14.8884,14.8794,14.8794,14.8832,14.8859,14.8815,14.8721,14.8622,14.8578,14.8546,14.8572,14.8597,14.8615,14.8626,14.8635,14.8649,14.866,14.8673,14.8695,14.8697,14.8691,14.8692,14.869,14.869,14.8688,14.8687,14.8685,14.8684,14.8688,14.8686,14.8676,14.8674,14.8674,14.8673,14.8671,14.8669,14.8669,14.8659,14.8656,14.865,14.8645,14.8641,14.8634,14.8632,14.8621,14.8611,14.8597,14.8591,14.8577,14.856,14.8551,14.8531,14.8505,14.8482,14.8469,14.8448,14.8441,14.843,14.8417,14.8395,14.8378,14.8352,14.8335,14.8316,14.8292,14.8273,14.825,14.8237,14.8237,14.8245,14.8253,14.8259,14.8272,14.8287,14.8304,14.8314,14.8322,14.8333,14.8339,14.8349,14.8358,14.8366,14.837,14.8377,14.8386,14.8389,14.8393,14.8395,14.8395,14.8398,14.8401,14.8401,14.8395,14.839,14.8384,14.8377,14.8371,14.8363,14.8358,14.8355,14.8346,14.8334,14.8322,14.8308,14.8293,14.8281,14.8259,14.8254,14.8249,14.8241,14.823,14.8222,14.8208,14.8177,14.8155,14.8127,14.8101,14.8078,14.8051,14.8003,14.7966,14.7929,14.7886,14.7842,14.781,14.7766,14.7726,14.7667,14.763,14.7582,14.7532,14.7476,14.7445,14.728,14.7115,14.695,14.6786,14.6621,14.6569,14.6516,14.6464,14.6412,14.6359,14.6307,14.6255,14.6218,14.6212,14.6202,14.615,14.6097,14.6348,14.6457,14.6516,14.6684,14.6851,14.7019,14.6785,14.6551,14.6551,14.6413,14.6233,14.6206,14.6183,14.6181,14.6179,14.6162,14.6151,14.614,14.612,14.6057,14.6057,14.6005,14.6005,14.5992,14.5986,14.598,14.5926,14.5877,14.5829,14.5829,14.583,14.5831,14.5827,14.5831,14.583,14.5827,14.5827,14.5818,14.5817,14.5807,14.5804,14.5793,14.5763,14.5761,14.5761,14.5758,14.5758,14.5755,14.5754,14.5745,14.5735,14.5738,14.5747,14.5746,14.5745,14.5742,14.5717,14.568,14.5678,14.5667,14.5638,14.5627,14.561,14.561,14.5604,14.562,14.5623,14.5628,14.5632,14.5627,14.5623,14.5621,14.5621,14.5619,14.5619,14.562,14.562,14.5626,14.563,14.5634,14.5639,14.564,14.5647,14.5651,14.5654,14.5655,14.5658,14.5665,14.5675,14.5702,14.5724,14.5743,14.577,14.5776,14.5779,14.5787,14.5787,14.579,14.5807,14.5802,14.5791,14.5799,14.58,14.5801,14.5812,14.5823,14.5823,14.5838,14.5838,14.5864,14.5864,14.5864,14.5864,14.586,14.5845,14.5758,14.5747,14.573,14.5709,14.5588,14.5588,14.5578,14.5578,14.5412,14.5404,14.5395,14.5385,14.5385,14.5335,14.531,14.531,14.5304,14.5145,14.4984,14.4973,14.493,14.4911,14.4864,14.4859,14.484,14.4741,14.4696,14.4617,14.4493,14.437,14.4246,14.4123,14.4117,14.4116,14.411,14.411,14.4109,14.4109,14.4109,14.4103,14.4096,14.4094,14.4089,14.4083,14.4076,14.4069,14.4062,14.4055,14.4052,14.4048,14.4163,14.4278,14.4393,14.4507,14.4622,14.4737,14.4852,14.4967,14.5081,14.5196,14.5311,14.5425,14.554,14.5655,14.577,14.5884,14.5999,14.6114,14.6069,14.6064,14.6025,14.5981,14.5936,14.5892,14.5849,14.5847,14.5803,14.5764,14.5726,14.5654,14.5592,14.5583,14.562,14.5657,14.5694,14.573,14.5767,14.5804,14.5915,14.6025,14.6135,14.6362,14.6589,14.6816,14.7043,14.7087,14.7329,14.7616,14.7902,14.818,14.8188,14.8474,14.8761,14.8798,14.9052,14.9057,14.9071,14.9079,14.9091,14.91,14.911,14.912,14.9129,14.9139,14.9151,14.9169,14.9179,14.9194,14.9206,14.922,14.9229,14.9245,14.9263,14.9279,14.9294,14.9308,14.9321,14.9333,14.9341,14.9356,14.9363,14.9375,14.9386,14.9401,14.9412,14.943,14.944,14.9454,14.9466,14.9478,14.9492,14.9502,14.9512,14.9521,14.953,14.9552,14.956,14.9568,14.958,14.9593,14.9608,14.9623,14.9628,14.9632,14.9634,14.964,14.9647,14.965,14.9656,14.9666,14.9675,14.9682,14.969,14.9703,14.9719,14.9724,14.9726,14.973,14.9742,14.9751,14.9759,14.9766,14.9775,14.9783,14.9794,14.98,14.9808,14.9813,14.9822,14.9824,14.9829,14.9837,14.985,14.9861,14.9872,14.9882,14.9905,14.9916,14.9933,14.9955,14.9973,15.0004,15.0273,15.0543,15.087,15.0857,15.0877,15.0907,15.0943,15.0974,15.1015,15.1057,15.109,15.1129,15.1163,15.119,15.1224,15.1249,15.1274,15.1302,15.1321,15.1346,15.1365,15.1393,15.141,15.1479,15.1485,15.1493,15.1568,15.1635,15.1704,15.1774,15.184,15.191,15.1979,15.2035,15.2082,15.2102,15.2132,15.2154,15.2165,15.2154,15.2104,15.2082,15.2068,15.2052,15.2057,15.2077,15.211,15.2149,15.2204,15.2224,15.2265,15.2271,15.2304,15.2332,15.2365,15.2404,15.2438,15.2471,15.251,15.2565,15.2613,15.2657,15.2685,15.2679,15.2668,15.2649,15.2635,15.2624,15.2632,15.2663,15.2699,15.2754,15.281,15.2865,15.2918,15.2982,15.3029,15.3077,15.311,15.3149,15.3182,15.3224,15.3257,15.3282,15.331,15.3343,15.3368,15.3388,15.3415,15.3435,15.346,15.3485,15.3513,15.3529,15.3549,15.3568,15.3602,15.3635,15.3677,15.3724,15.3779,15.3832,15.3888,15.3957,15.4027,15.4093,15.4157,15.4218,15.4274,15.4324,15.4365,15.4413,15.4432,15.4452,15.4471,15.4527,15.4579,15.4621,15.4663,15.4649,15.469,15.4743,15.4804,15.4874,15.4938,15.499,15.4996,15.5032,15.5057,15.5057,15.504,15.5018,15.4996,15.4988,15.4971,15.4967,15.4949,15.4927,15.4918,15.491,15.4899,15.4877,15.4857,15.4827,15.4804,15.479,15.4788,15.4799,15.4832,15.4857,15.489,15.491,15.4938,15.4963,15.4996,15.4997,15.5004,15.5043,15.5082,15.5193,15.5193,15.5254,15.531,15.5354,15.5396,15.5429,15.5474,15.5529,15.559,15.5652,15.5707,15.5768,15.5829,15.5857,15.5877,15.5902,15.5929,15.5946,15.5974,15.5999,15.6027,15.606,15.6079,15.6113,15.6129,15.6149,15.616,15.6179,15.6199,15.6232,15.626,15.6279,15.629,15.6282,15.6274,15.6257,15.6235,15.6218,15.6196,15.6193,15.6213,15.6246,15.6293,15.6335,15.6382,15.6424,15.6471,15.6524,15.6579,15.6649,15.6704,15.6752,15.6802,15.6838,15.6852,15.6849,15.6863,15.6904,15.6954,15.7015,15.7057,15.709,15.7124,15.7157,15.7196,15.7229,15.7257,15.729,15.7324,15.7354,15.7382,15.7407,15.744,15.7474,15.7515,15.7549,15.7563,15.7602,15.7649,15.7704,15.7765,15.7821,15.7871,15.7893,15.7921,15.7957,15.8013,15.8046,15.8079,15.811,15.8143,15.8185,15.8218,15.8265,15.8313,15.8368,15.8421,15.8477,15.8538,15.8607,15.8677,15.8738,15.8793,15.884,15.8865,15.8893,15.8913,15.8918,15.8946,15.8985,15.9049,15.911,15.9165,15.9227,15.929,15.9329,15.9371,15.9404,15.9429,15.9463,15.9482,15.9507,15.9535,15.9568,15.9607,15.9654,15.9696,15.9743,15.9782,15.9832,15.9879,15.994,15.9993,15.9996,15.9997,16.0049,16.011,16.0174,16.0235,16.0296,16.0352,16.0404,16.046,16.0513,16.0549,16.0588,16.0615,16.064,16.066,16.0685,16.0699,16.0715,16.0743,16.0763,16.0779,16.0813,16.084,16.0879,16.0913,16.0946,16.0993,16.1035,16.1074,16.104,16.1004,16.096,16.0927,16.0882,16.0849,16.0804,16.0777,16.0781,16.0782,16.0807,16.0854,16.0896,16.0943,16.0996,16.106,16.1115,16.1149,16.1179,16.1207,16.1223,16.1238,16.126,16.1296,16.1338,16.1385,16.1449,16.1488,16.1492,16.1493,16.148,16.1477,16.1449,16.1413,16.1371,16.134,16.1318,16.1318,16.131,16.1313,16.1318,16.1324,16.1335,16.134,16.1352,16.1357,16.1354,16.1346,16.1324,16.1307,16.1285,16.1254,16.1235,16.1204,16.1188,16.1165,16.1163,16.1154,16.1138,16.111,16.1074,16.104,16.1004,16.0982,16.0971,16.0985,16.1018,16.1052,16.109,16.1124,16.1157,16.119,16.1215,16.1243,16.1277,16.1315,16.1385,16.144,16.149,16.1532,16.1574,16.1615,16.1657,16.1667,16.1713,16.1768,16.1815,16.1843,16.1846,16.1824,16.1802,16.1774,16.1738,16.1707,16.1679,16.1643,16.1615,16.1574,16.1543,16.1515,16.1493,16.1471,16.1454,16.1452,16.1443,16.1435,16.1418,16.1402,16.1385,16.1363,16.1335,16.1313,16.1282,16.1254,16.1227,16.119,16.116,16.1138,16.1179,16.1227,16.1274,16.1321,16.1346,16.1388,16.1435,16.1496,16.1543,16.1599,16.166,16.1729,16.1799,16.1868,16.1929,16.1977,16.1996,16.2002,16.2018,16.2032,16.2029,16.2049,16.2088,16.2121,16.2154,16.2196,16.2221,16.2246,16.2279,16.2307,16.2327,16.2324,16.2307,16.2302,16.229,16.2302,16.2335,16.2368,16.2424,16.2479,16.2529,16.2585,16.2632,16.2674,16.2729,16.2777,16.2804,16.2829,16.2863,16.2896,16.2952,16.3013,16.3063,16.3096,16.3113,16.3135,16.3163,16.3207,16.3254,16.3288,16.3313,16.3346,16.3379,16.3413,16.3446,16.3479,16.3513,16.3538,16.354,16.3549,16.3554,16.3552,16.3557,16.3574,16.3607,16.3649,16.3682,16.3729,16.3782,16.3824,16.3857,16.3882,16.391,16.3935,16.3963,16.3993,16.4027,16.406,16.4093,16.4121,16.4146,16.4171,16.4204,16.4229,16.4263,16.4296,16.4352,16.439,16.4432,16.4465,16.449,16.4518,16.4549,16.4582,16.4624,16.4671,16.4718,16.4757,16.4785,16.4802,16.4829,16.4854,16.4888,16.4935,16.4968,16.4996,16.4997,16.501,16.504,16.5082,16.5115,16.514,16.5174,16.5215,16.5254,16.5293,16.5335,16.5382,16.5415,16.544,16.5474,16.5493,16.5504,16.551,16.551,16.5493,16.5465,16.5443,16.5407,16.5365,16.5329,16.534,16.5374,16.5415,16.5468,16.5532,16.5588,16.5649,16.571,16.5765,16.5827,16.5846,16.5838,16.5821,16.5799,16.5768,16.5746,16.5724,16.571,16.5721,16.5738,16.5765,16.5804,16.5865,16.5915,16.5971,16.6027,16.6096,16.6157,16.6218,16.6282,16.6321,16.6354,16.6382,16.6393,16.6404,16.6404,16.6402,16.6399,16.6396,16.6388,16.6385,16.6404,16.6424,16.6429,16.6435,16.6454,16.6482,16.6499,16.6532,16.6557,16.6571,16.6563,16.6538,16.6504,16.646,16.6418,16.6371,16.6335,16.6313,16.631,16.6299,16.6296,16.6288,16.6271,16.6263,16.6262,16.626,16.6257,16.6263,16.6274,16.6285,16.6296,16.6315,16.6349,16.6374,16.6402,16.639,16.6363,16.6327,16.6296,16.6288,16.6296,16.6288,16.6285,16.6282,16.6279,16.6279,16.6274,16.6271,16.6268,16.6293,16.6349,16.6402,16.6435,16.6432,16.6424,16.6427,16.6432,16.6424,16.6407,16.639,16.6374,16.6374,16.6377,16.6388,16.6393,16.6371,16.634,16.6299,16.6263,16.6246,16.6265,16.6313,16.636,16.6402,16.644,16.6474,16.6485,16.6468,16.6449,16.6424,16.6396,16.6374,16.6385,16.6404,16.6438,16.6468,16.6496,16.6507,16.649,16.6463,16.6432,16.641,16.6393,16.6385,16.6368,16.6379,16.6389,16.6407,16.644,16.6471,16.6504,16.6529,16.6557,16.6527,16.6485,16.6443,16.6407,16.6371,16.6349,16.636,16.6407,16.6454,16.646,16.6435,16.6427,16.6424,16.6418,16.6424]}]],[[{"lng":[-17.4742,-17.4734,-17.4726,-17.4713,-17.4713,-17.4734,-17.4736,-17.475,-17.4755,-17.4751,-17.4746,-17.4763,-17.4759,-17.4742],"lat":[14.6583,14.6575,14.6575,14.6562,14.6554,14.6534,14.6533,14.6533,14.6546,14.6553,14.6562,14.6571,14.6584,14.6583]},{"lng":[-17.4001,-17.398,-17.398,-17.4009,-17.4025,-17.4038,-17.403,-17.4038,-17.403,-17.4026,-17.4001],"lat":[14.6634,14.6654,14.668,14.6717,14.6725,14.6712,14.6696,14.6688,14.6679,14.6633,14.6634]},{"lng":[-17.5443,-17.5438,-17.5438,-17.5442,-17.545,-17.5455,-17.5455,-17.5451,-17.5443],"lat":[14.7417,14.7421,14.7446,14.745,14.745,14.7446,14.7421,14.7417,14.7417]},{"lng":[-17.5392,-17.5371,-17.5376,-17.5384,-17.5405,-17.5401,-17.5392],"lat":[14.7475,14.7496,14.75,14.75,14.7479,14.7475,14.7475]},{"lng":[-17.5134,-17.5142,-17.5159,-17.5167,-17.5184,-17.5197,-17.5196,-17.5182,-17.5159,-17.5142,-17.5129,-17.513,-17.5134],"lat":[14.7584,14.7591,14.7592,14.7583,14.7584,14.7571,14.7562,14.7561,14.7559,14.755,14.7562,14.758,14.7584]},{"lng":[-17.4775,-17.4771,-17.4771,-17.4771,-17.4776,-17.4792,-17.4796,-17.4796,-17.4797,-17.4792,-17.4775],"lat":[14.7692,14.7696,14.7711,14.7712,14.7716,14.7717,14.7713,14.7705,14.7695,14.7692,14.7692]},{"lng":[-16.7366,-16.7381,-16.7405,-16.7415,-16.7432,-16.7423,-16.7401,-16.7401,-16.738,-16.7373,-16.7338,-16.7337,-16.7285,-16.7277,-16.7274,-16.7253,-16.7247,-16.7231,-16.7209,-16.7187,-16.7155,-16.7127,-16.7075,-16.7033,-16.7031,-16.6995,-16.6973,-16.6908,-16.6889,-16.6887,-16.6883,-16.6888,-16.6886,-16.6892,-16.6901,-16.6911,-16.6958,-16.7012,-16.7102,-16.7307,-16.7313,-16.7347,-16.7363,-16.7405,-16.7421,-16.7455,-16.7455,-16.748,-16.7513,-16.7546,-16.7546,-16.7596,-16.7613,-16.763,-16.7638,-16.7655,-16.7663,-16.768,-16.7713,-16.773,-16.7746,-16.7746,-16.7763,-16.7763,-16.7805,-16.7821,-16.7838,-16.7863,-16.7872,-16.7913,-16.7988,-16.7988,-16.8055,-16.8063,-16.8088,-16.8155,-16.8188,-16.8213,-16.8238,-16.8263,-16.8296,-16.8321,-16.8338,-16.8363,-16.8363,-16.8388,-16.8397,-16.8413,-16.8446,-16.8463,-16.8463,-16.848,-16.8505,-16.8521,-16.8572,-16.8588,-16.8663,-16.8663,-16.8688,-16.8688,-16.873,-16.8738,-16.8763,-16.8771,-16.8838,-16.8838,-16.8871,-16.888,-16.893,-16.893,-16.9055,-16.9055,-16.9096,-16.9105,-16.913,-16.9196,-16.9205,-16.9238,-16.9305,-16.9313,-16.9388,-16.9405,-16.9496,-16.9497,-16.958,-16.9588,-16.9588,-16.9704,-16.9755,-16.9755,-16.9788,-16.9788,-16.988,-16.988,-16.9938,-16.9946,-16.9992,-17.0063,-17.0063,-17.0146,-17.0163,-17.0255,-17.0255,-17.0363,-17.0363,-17.0413,-17.0413,-17.0472,-17.0472,-17.0596,-17.0597,-17.0646,-17.0646,-17.068,-17.0738,-17.0738,-17.0763,-17.0772,-17.0909,-17.0917,-17.0946,-17.098,-17.1063,-17.108,-17.1105,-17.1138,-17.1155,-17.1163,-17.118,-17.1196,-17.1213,-17.123,-17.1238,-17.1255,-17.1255,-17.1288,-17.1301,-17.1309,-17.1325,-17.1375,-17.1417,-17.1517,-17.1551,-17.1634,-17.1684,-17.1785,-17.1817,-17.1951,-17.1975,-17.2017,-17.2084,-17.2125,-17.2217,-17.2275,-17.2334,-17.235,-17.2417,-17.2434,-17.2459,-17.2476,-17.2509,-17.2526,-17.2551,-17.2567,-17.2592,-17.2609,-17.2742,-17.2759,-17.2859,-17.2875,-17.2951,-17.2967,-17.2992,-17.3009,-17.3018,-17.3067,-17.3092,-17.3109,-17.3134,-17.3251,-17.3268,-17.3284,-17.3309,-17.3326,-17.3376,-17.3392,-17.3509,-17.3525,-17.3568,-17.3592,-17.3626,-17.365,-17.3676,-17.3692,-17.3742,-17.3767,-17.38,-17.3826,-17.3842,-17.3867,-17.39,-17.3951,-17.3967,-17.3992,-17.4034,-17.4042,-17.4067,-17.4076,-17.4109,-17.4117,-17.4142,-17.415,-17.4168,-17.4176,-17.4218,-17.4226,-17.4242,-17.4293,-17.4309,-17.4333,-17.4376,-17.4425,-17.4434,-17.4459,-17.4492,-17.4542,-17.455,-17.4584,-17.4608,-17.4626,-17.4627,-17.4651,-17.4684,-17.4692,-17.4718,-17.4725,-17.4759,-17.4792,-17.4809,-17.4842,-17.4896,-17.49,-17.4934,-17.4951,-17.4967,-17.4992,-17.5009,-17.5017,-17.5051,-17.5068,-17.5101,-17.5142,-17.5184,-17.5196,-17.5209,-17.5234,-17.5243,-17.5259,-17.5263,-17.5292,-17.5301,-17.5313,-17.5313,-17.5313,-17.533,-17.5225,-17.5217,-17.5201,-17.5184,-17.508,-17.508,-17.5071,-17.508,-17.5059,-17.5025,-17.5,-17.4976,-17.4951,-17.493,-17.493,-17.49,-17.488,-17.4888,-17.4875,-17.4859,-17.483,-17.4829,-17.4817,-17.48,-17.478,-17.478,-17.4805,-17.4805,-17.478,-17.4788,-17.4784,-17.4767,-17.4743,-17.473,-17.4735,-17.4738,-17.473,-17.473,-17.4747,-17.4746,-17.4726,-17.4717,-17.4701,-17.4688,-17.4684,-17.4658,-17.4646,-17.4642,-17.4634,-17.4601,-17.4584,-17.4542,-17.4534,-17.4517,-17.4488,-17.4488,-17.4463,-17.4463,-17.4463,-17.4442,-17.4396,-17.4396,-17.4386,-17.438,-17.4358,-17.4334,-17.4321,-17.4315,-17.4313,-17.4313,-17.4313,-17.4338,-17.4338,-17.4318,-17.4267,-17.4251,-17.4221,-17.4205,-17.4205,-17.4226,-17.425,-17.4259,-17.4268,-17.428,-17.4284,-17.43,-17.4317,-17.4325,-17.433,-17.4322,-17.4338,-17.433,-17.4355,-17.4355,-17.4342,-17.4334,-17.4317,-17.4301,-17.4293,-17.4267,-17.4251,-17.4238,-17.4246,-17.4234,-17.4225,-17.4205,-17.4205,-17.4221,-17.4221,-17.4234,-17.425,-17.4271,-17.4271,-17.4263,-17.4242,-17.4217,-17.418,-17.4179,-17.4192,-17.4225,-17.425,-17.4276,-17.4321,-17.433,-17.433,-17.4321,-17.4313,-17.4313,-17.4296,-17.4217,-17.4209,-17.4192,-17.4159,-17.4101,-17.4092,-17.4059,-17.405,-17.4017,-17.4009,-17.3967,-17.3959,-17.3867,-17.3859,-17.3684,-17.3675,-17.3625,-17.3617,-17.3584,-17.3576,-17.3542,-17.3534,-17.3475,-17.3392,-17.3258,-17.3242,-17.3125,-17.3075,-17.3059,-17.3017,-17.3,-17.2992,-17.2975,-17.2942,-17.2901,-17.2859,-17.2826,-17.2784,-17.2759,-17.2742,-17.2726,-17.2692,-17.2676,-17.2659,-17.2642,-17.2617,-17.2601,-17.2525,-17.2484,-17.2467,-17.2351,-17.2284,-17.2259,-17.2217,-17.2208,-17.2192,-17.2084,-17.2076,-17.205,-17.2026,-17.1971,-17.1971,-17.193,-17.193,-17.1867,-17.1859,-17.1755,-17.1755,-17.1696,-17.1663,-17.1663,-17.1613,-17.1605,-17.1588,-17.1571,-17.1546,-17.1521,-17.1521,-17.1505,-17.1496,-17.1496,-17.1471,-17.1459,-17.1446,-17.1446,-17.143,-17.143,-17.1413,-17.1372,-17.1346,-17.1334,-17.1267,-17.1242,-17.1159,-17.1113,-17.1113,-17.1088,-17.1088,-17.1063,-17.1063,-17.1038,-17.1013,-17.1005,-17.098,-17.0971,-17.0946,-17.0946,-17.0938,-17.093,-17.093,-17.0921,-17.0922,-17.0896,-17.0896,-17.0896,-17.088,-17.0805,-17.0805,-17.0747,-17.0755,-17.0726,-17.0712,-17.0709,-17.07,-17.0684,-17.0638,-17.0638,-17.0613,-17.0605,-17.0596,-17.0596,-17.0596,-17.0588,-17.0588,-17.0559,-17.0523,-17.0509,-17.0459,-17.0451,-17.0434,-17.0392,-17.0376,-17.0351,-17.0325,-17.0309,-17.025,-17.0217,-17.0209,-17.0167,-17.0151,-17.0142,-17.0117,-17.0059,-16.9992,-16.9972,-16.9971,-16.9867,-16.9834,-16.9805,-16.9751,-16.9734,-16.9709,-16.9692,-16.9655,-16.9613,-16.9596,-16.9551,-16.9525,-16.9513,-16.9505,-16.9496,-16.9471,-16.9471,-16.9463,-16.9463,-16.9396,-16.9355,-16.9355,-16.9338,-16.933,-16.9313,-16.9313,-16.9321,-16.933,-16.9338,-16.9338,-16.9346,-16.9346,-16.9355,-16.9355,-16.9346,-16.9355,-16.9346,-16.9338,-16.9322,-16.9321,-16.9313,-16.9309,-16.9308,-16.9302,-16.9275,-16.9259,-16.9226,-16.9201,-16.9175,-16.9105,-16.9071,-16.8988,-16.8942,-16.8901,-16.8863,-16.8855,-16.8821,-16.878,-16.8771,-16.875,-16.8738,-16.874,-16.8744,-16.8746,-16.8755,-16.875,-16.8738,-16.873,-16.873,-16.8721,-16.8721,-16.8713,-16.8713,-16.8705,-16.8705,-16.8713,-16.8713,-16.8671,-16.8671,-16.8678,-16.8688,-16.8688,-16.8675,-16.8662,-16.8659,-16.8592,-16.853,-16.8521,-16.8513,-16.8513,-16.8484,-16.8475,-16.8384,-16.8325,-16.8301,-16.8284,-16.8242,-16.8226,-16.8213,-16.8225,-16.8242,-16.8255,-16.8255,-16.8242,-16.8226,-16.8201,-16.8159,-16.815,-16.8134,-16.813,-16.8151,-16.8176,-16.8196,-16.8197,-16.8167,-16.8084,-16.8055,-16.8055,-16.8038,-16.7998,-16.7849,-16.7835,-16.7818,-16.7798,-16.778,-16.7766,-16.7763,-16.7737,-16.7725,-16.771,-16.7693,-16.7677,-16.7655,-16.7641,-16.7633,-16.7625,-16.7615,-16.7593,-16.7585,-16.7573,-16.7559,-16.7548,-16.7537,-16.7533,-16.7526,-16.7527,-16.7533,-16.7535,-16.7551,-16.7559,-16.7565,-16.7574,-16.7597,-16.7617,-16.763,-16.7634,-16.7635,-16.7624,-16.7612,-16.7601,-16.7588,-16.7574,-16.7557,-16.7539,-16.7523,-16.7507,-16.7493,-16.7479,-16.7468,-16.7456,-16.745,-16.7437,-16.7428,-16.7413,-16.7402,-16.7392,-16.7375,-16.7356,-16.7341,-16.7329,-16.7321,-16.7304,-16.7291,-16.7278,-16.7265,-16.7256,-16.7248,-16.7241,-16.7235,-16.7223,-16.7222,-16.7228,-16.7232,-16.7239,-16.7254,-16.7281,-16.7303,-16.7312,-16.7322,-16.7348,-16.7368,-16.7387,-16.7401,-16.7422,-16.7423,-16.7431,-16.7433,-16.7441,-16.7457,-16.746,-16.7476,-16.7474,-16.7486,-16.7488,-16.7488,-16.7484,-16.7475,-16.7472,-16.7463,-16.7452,-16.7437,-16.7425,-16.7409,-16.739,-16.7372,-16.7362,-16.7344,-16.7333,-16.7318,-16.7306,-16.729,-16.7277,-16.7264,-16.7246,-16.7232,-16.7219,-16.7208,-16.7198,-16.7187,-16.7178,-16.7166,-16.715,-16.7138,-16.7124,-16.7109,-16.7096,-16.7088,-16.7085,-16.7078,-16.7073,-16.7071,-16.7075,-16.7074,-16.7079,-16.7078,-16.7079,-16.7089,-16.7089,-16.7092,-16.7097,-16.7102,-16.7103,-16.7098,-16.71,-16.7102,-16.7105,-16.7104,-16.7091,-16.7088,-16.7079,-16.7079,-16.7083,-16.7079,-16.7077,-16.7072,-16.7068,-16.7065,-16.7069,-16.7064,-16.7061,-16.7053,-16.7054,-16.7052,-16.7038,-16.702,-16.7002,-16.6981,-16.6963,-16.6948,-16.6928,-16.6906,-16.6888,-16.6869,-16.685,-16.6832,-16.6811,-16.6795,-16.6777,-16.6755,-16.6733,-16.671,-16.6684,-16.6659,-16.6625,-16.6596,-16.6568,-16.6541,-16.6516,-16.6492,-16.6466,-16.6457,-16.6439,-16.6424,-16.6404,-16.6389,-16.6377,-16.6352,-16.6336,-16.6327,-16.6308,-16.6297,-16.628,-16.6263,-16.6259,-16.625,-16.6234,-16.6217,-16.6202,-16.6175,-16.6155,-16.614,-16.6119,-16.6106,-16.6089,-16.6073,-16.6057,-16.6037,-16.6017,-16.5997,-16.5972,-16.595,-16.592,-16.5906,-16.589,-16.5867,-16.5841,-16.5814,-16.5794,-16.5771,-16.5753,-16.5735,-16.5711,-16.5694,-16.5681,-16.5653,-16.5635,-16.5611,-16.559,-16.5574,-16.5545,-16.5526,-16.5512,-16.5506,-16.5495,-16.5493,-16.5484,-16.5474,-16.5459,-16.5443,-16.5435,-16.5426,-16.5412,-16.5405,-16.5394,-16.5387,-16.5382,-16.537,-16.5361,-16.5361,-16.5358,-16.5348,-16.5349,-16.5346,-16.5335,-16.533,-16.5322,-16.5317,-16.5314,-16.531,-16.5306,-16.5308,-16.5303,-16.5299,-16.5294,-16.5287,-16.5286,-16.5285,-16.5279,-16.5268,-16.5269,-16.5271,-16.5273,-16.5273,-16.5267,-16.5257,-16.5246,-16.5234,-16.5225,-16.5217,-16.5212,-16.5206,-16.5198,-16.519,-16.5186,-16.518,-16.5177,-16.5178,-16.5176,-16.5176,-16.5177,-16.5178,-16.518,-16.5181,-16.5184,-16.5189,-16.5193,-16.5199,-16.5207,-16.5215,-16.5223,-16.5229,-16.5234,-16.5241,-16.5241,-16.5259,-16.5278,-16.5285,-16.5296,-16.5303,-16.5311,-16.5322,-16.5339,-16.5348,-16.5362,-16.5374,-16.5386,-16.5403,-16.5421,-16.5436,-16.5452,-16.5467,-16.5481,-16.55,-16.5516,-16.5533,-16.5546,-16.5562,-16.5568,-16.5591,-16.5608,-16.562,-16.563,-16.5647,-16.5654,-16.5667,-16.5682,-16.5694,-16.5701,-16.5727,-16.5737,-16.575,-16.5768,-16.5776,-16.5791,-16.5802,-16.5812,-16.5819,-16.5829,-16.5845,-16.5849,-16.5862,-16.587,-16.5883,-16.5895,-16.5905,-16.5914,-16.5926,-16.594,-16.5948,-16.5958,-16.5969,-16.5982,-16.5994,-16.6011,-16.6022,-16.6043,-16.6056,-16.6073,-16.609,-16.6102,-16.6118,-16.6134,-16.615,-16.6163,-16.6171,-16.6192,-16.6206,-16.6209,-16.6214,-16.6211,-16.6212,-16.6208,-16.6206,-16.6201,-16.6196,-16.6185,-16.6169,-16.6156,-16.6142,-16.6123,-16.6108,-16.6104,-16.6088,-16.6075,-16.6064,-16.6056,-16.6053,-16.6047,-16.6045,-16.6047,-16.6046,-16.6044,-16.6047,-16.6047,-16.6051,-16.6055,-16.6051,-16.605,-16.6054,-16.6055,-16.6056,-16.6053,-16.6053,-16.605,-16.605,-16.6044,-16.6042,-16.6039,-16.6039,-16.6048,-16.6051,-16.6053,-16.6054,-16.6055,-16.6058,-16.6058,-16.6062,-16.6064,-16.6068,-16.607,-16.6075,-16.608,-16.6093,-16.6109,-16.613,-16.6139,-16.6158,-16.6222,-16.6244,-16.6264,-16.6285,-16.6306,-16.6331,-16.6356,-16.6383,-16.6408,-16.6429,-16.6451,-16.6484,-16.6486,-16.6474,-16.6468,-16.6453,-16.6445,-16.6437,-16.6428,-16.6411,-16.6395,-16.6374,-16.6353,-16.6339,-16.632,-16.6293,-16.6276,-16.6247,-16.6237,-16.6233,-16.6242,-16.6242,-16.6248,-16.6248,-16.6256,-16.6261,-16.6267,-16.6271,-16.6285,-16.6293,-16.6304,-16.631,-16.6319,-16.6331,-16.6346,-16.636,-16.6369,-16.6384,-16.6401,-16.6419,-16.6433,-16.6451,-16.6468,-16.6495,-16.6516,-16.6524,-16.6534,-16.6545,-16.6552,-16.656,-16.6568,-16.6578,-16.6583,-16.6594,-16.6608,-16.6615,-16.6633,-16.6648,-16.6652,-16.6663,-16.6678,-16.6695,-16.6707,-16.6717,-16.6729,-16.6741,-16.6769,-16.6793,-16.6813,-16.6833,-16.6855,-16.6841,-16.6802,-16.6773,-16.6734,-16.6685,-16.6652,-16.6628,-16.6615,-16.658,-16.6542,-16.6522,-16.6509,-16.6479,-16.6444,-16.6402,-16.6349,-16.6288,-16.6211,-16.6143,-16.607,-16.6015,-16.5963,-16.5913,-16.5871,-16.5799,-16.5738,-16.5645,-16.5608,-16.5576,-16.5519,-16.5438,-16.5377,-16.5323,-16.527,-16.5202,-16.5157,-16.5135,-16.5102,-16.5066,-16.5029,-16.5002,-16.4981,-16.4961,-16.4965,-16.4939,-16.4921,-16.4889,-16.4858,-16.4815,-16.4789,-16.4768,-16.4755,-16.4742,-16.4732,-16.4713,-16.4677,-16.4636,-16.4591,-16.4553,-16.4521,-16.4469,-16.4454,-16.4424,-16.4405,-16.438,-16.4353,-16.4334,-16.4314,-16.4291,-16.4263,-16.4232,-16.4198,-16.4169,-16.4138,-16.4107,-16.4084,-16.406,-16.4022,-16.3988,-16.394,-16.3905,-16.3901,-16.3879,-16.3858,-16.3833,-16.3788,-16.3764,-16.3736,-16.3709,-16.3683,-16.3632,-16.3592,-16.3556,-16.3512,-16.3485,-16.3462,-16.3427,-16.3395,-16.3377,-16.3346,-16.3319,-16.3301,-16.328,-16.3254,-16.3225,-16.3198,-16.3188,-16.313,-16.3099,-16.3081,-16.3049,-16.3029,-16.3008,-16.2988,-16.2976,-16.2939,-16.2911,-16.2885,-16.2867,-16.2841,-16.281,-16.2786,-16.2769,-16.2746,-16.2724,-16.2708,-16.2683,-16.2654,-16.2626,-16.2601,-16.2579,-16.255,-16.253,-16.2502,-16.2483,-16.2472,-16.2454,-16.2426,-16.2408,-16.2396,-16.237,-16.2357,-16.2345,-16.2329,-16.232,-16.2322,-16.2316,-16.2311,-16.2308,-16.2294,-16.2275,-16.2253,-16.2237,-16.2206,-16.2183,-16.2153,-16.2123,-16.2101,-16.2068,-16.2036,-16.2014,-16.1988,-16.1966,-16.1945,-16.1932,-16.1904,-16.188,-16.1854,-16.1817,-16.1782,-16.177,-16.177,-16.1771,-16.1772,-16.1771,-16.1771,-16.1766,-16.1761,-16.1755,-16.1759,-16.1764,-16.1767,-16.1771,-16.1776,-16.1779,-16.1781,-16.1782,-16.1787,-16.1783,-16.1782,-16.1784,-16.1789,-16.1796,-16.18,-16.1808,-16.1824,-16.1832,-16.1858,-16.1873,-16.1886,-16.1891,-16.1904,-16.1921,-16.1931,-16.1944,-16.1943,-16.1942,-16.1947,-16.1947,-16.1947,-16.1947,-16.1947,-16.1951,-16.1962,-16.1968,-16.1989,-16.2012,-16.2041,-16.205,-16.2086,-16.2116,-16.2136,-16.2188,-16.2249,-16.227,-16.2304,-16.2321,-16.2345,-16.2364,-16.2371,-16.2374,-16.2377,-16.2373,-16.2372,-16.2373,-16.2376,-16.2392,-16.2416,-16.2426,-16.2438,-16.245,-16.2475,-16.2499,-16.2526,-16.2551,-16.2572,-16.2587,-16.2601,-16.2616,-16.2632,-16.2647,-16.2661,-16.2672,-16.2683,-16.2688,-16.2693,-16.2699,-16.2697,-16.2692,-16.2689,-16.2685,-16.2681,-16.268,-16.2685,-16.2692,-16.27,-16.2708,-16.2723,-16.2742,-16.276,-16.278,-16.2806,-16.2832,-16.2849,-16.2858,-16.2868,-16.2885,-16.2908,-16.2931,-16.295,-16.2969,-16.3,-16.3017,-16.3045,-16.3064,-16.308,-16.3084,-16.3088,-16.3093,-16.3102,-16.3112,-16.3142,-16.3171,-16.3202,-16.3237,-16.3274,-16.3305,-16.3332,-16.335,-16.338,-16.3415,-16.3436,-16.345,-16.3457,-16.3454,-16.3446,-16.3416,-16.3388,-16.3371,-16.3354,-16.3344,-16.3346,-16.3356,-16.3371,-16.3393,-16.3411,-16.3421,-16.3448,-16.3468,-16.3496,-16.3515,-16.3547,-16.3587,-16.363,-16.366,-16.3693,-16.3719,-16.3732,-16.3752,-16.3773,-16.3796,-16.3839,-16.3878,-16.3932,-16.3979,-16.4039,-16.4107,-16.4153,-16.4189,-16.4228,-16.4255,-16.428,-16.4311,-16.4338,-16.4362,-16.4385,-16.4414,-16.4436,-16.4458,-16.4483,-16.4509,-16.4515,-16.4526,-16.4538,-16.4546,-16.4564,-16.4594,-16.4618,-16.4641,-16.4659,-16.4691,-16.4709,-16.4736,-16.4751,-16.4759,-16.4779,-16.4791,-16.4804,-16.4817,-16.4832,-16.4845,-16.4861,-16.4873,-16.4889,-16.4901,-16.4915,-16.4941,-16.4954,-16.4974,-16.499,-16.5007,-16.502,-16.503,-16.5044,-16.506,-16.5077,-16.5096,-16.5111,-16.5124,-16.5137,-16.5146,-16.5158,-16.5162,-16.5164,-16.5156,-16.5143,-16.5135,-16.5123,-16.5123,-16.5124,-16.5129,-16.5132,-16.5135,-16.5145,-16.5158,-16.5173,-16.5186,-16.5202,-16.5225,-16.5237,-16.5251,-16.5279,-16.5307,-16.5332,-16.535,-16.5381,-16.54,-16.5418,-16.5451,-16.5471,-16.5519,-16.554,-16.5571,-16.5596,-16.5647,-16.5671,-16.5696,-16.5722,-16.581,-16.5824,-16.5836,-16.5907,-16.5958,-16.5993,-16.6012,-16.6037,-16.6045,-16.6072,-16.6092,-16.6128,-16.6141,-16.614,-16.6134,-16.6129,-16.6123,-16.6121,-16.6121,-16.6133,-16.6165,-16.6207,-16.6245,-16.6285,-16.6306,-16.6307,-16.6312,-16.6338,-16.6365,-16.6384,-16.6412,-16.6423,-16.6428,-16.6476,-16.6536,-16.6557,-16.6608,-16.6626,-16.6635,-16.6669,-16.6687,-16.6701,-16.673,-16.6745,-16.676,-16.6811,-16.6828,-16.6844,-16.6867,-16.6907,-16.6922,-16.6961,-16.701,-16.7022,-16.7048,-16.7073,-16.7081,-16.7094,-16.7136,-16.7155,-16.7167,-16.718,-16.7233,-16.7305,-16.7319,-16.7333,-16.7366],"lat":[15.2967,15.297,15.2991,15.2994,15.3026,15.31,15.3152,15.3166,15.3241,15.3313,15.3395,15.3406,15.3523,15.3549,15.3577,15.364,15.3645,15.3693,15.3726,15.3804,15.3866,15.3905,15.3963,15.4021,15.403,15.4068,15.4086,15.4167,15.422,15.4226,15.4259,15.4278,15.4303,15.4333,15.4348,15.4379,15.443,15.4471,15.4525,15.4614,15.4596,15.4554,15.4512,15.4462,15.4429,15.4387,15.4379,15.4354,15.4296,15.4254,15.4238,15.4163,15.4129,15.4112,15.4088,15.4071,15.4046,15.4029,15.3979,15.392,15.3904,15.3888,15.3871,15.3854,15.3804,15.3771,15.3737,15.3712,15.3688,15.3638,15.3504,15.3496,15.3421,15.3396,15.3371,15.3238,15.3195,15.3138,15.3104,15.3046,15.3012,15.2971,15.2938,15.2913,15.2904,15.2879,15.2854,15.2837,15.2771,15.2754,15.2738,15.2713,15.2688,15.2646,15.2579,15.2546,15.2463,15.2454,15.2429,15.2421,15.2371,15.2346,15.2321,15.2295,15.2204,15.2196,15.2163,15.2137,15.2079,15.2071,15.1954,15.1946,15.1904,15.1888,15.1837,15.1771,15.1754,15.1696,15.1621,15.1596,15.1504,15.1471,15.1371,15.1363,15.1279,15.1262,15.1254,15.1129,15.1062,15.1054,15.1021,15.1012,15.0912,15.0904,15.0846,15.0821,15.0784,15.0712,15.0696,15.0612,15.0571,15.0487,15.0479,15.037,15.0363,15.0313,15.0304,15.0246,15.0238,15.0105,15.0096,15.0046,15.0037,14.9996,14.9946,14.9937,14.9912,14.9888,14.975,14.975,14.9721,14.9654,14.9571,14.9538,14.9512,14.9454,14.9446,14.9421,14.9404,14.9379,14.9313,14.9279,14.9237,14.9204,14.9187,14.9145,14.9125,14.9125,14.9109,14.9091,14.9042,14.8992,14.8983,14.8925,14.89,14.885,14.8834,14.8767,14.875,14.8742,14.8708,14.87,14.865,14.8608,14.8583,14.8567,14.8534,14.8533,14.8516,14.8517,14.85,14.85,14.8484,14.8483,14.8467,14.8467,14.8391,14.8392,14.8342,14.8342,14.83,14.83,14.8283,14.8283,14.8275,14.8267,14.825,14.825,14.8234,14.8175,14.8159,14.8158,14.8141,14.8142,14.8117,14.8116,14.8059,14.8058,14.8033,14.8034,14.8017,14.8017,14.8,14.8,14.7975,14.795,14.7933,14.7908,14.7908,14.7892,14.7884,14.7858,14.7859,14.7842,14.7833,14.7825,14.7825,14.7817,14.7816,14.7808,14.7809,14.78,14.78,14.7792,14.7784,14.7775,14.7775,14.775,14.775,14.7725,14.7709,14.7708,14.77,14.77,14.7683,14.7675,14.7667,14.7658,14.7642,14.7642,14.7641,14.7625,14.7625,14.7633,14.7633,14.7642,14.7642,14.7658,14.7675,14.7675,14.762,14.76,14.7592,14.7567,14.7575,14.7575,14.7592,14.7592,14.7559,14.7566,14.7567,14.7517,14.7525,14.7512,14.7492,14.7492,14.7501,14.7517,14.7517,14.7516,14.7496,14.7471,14.7469,14.7437,14.7413,14.7408,14.74,14.74,14.7392,14.7287,14.7271,14.7262,14.7204,14.7183,14.7183,14.715,14.715,14.7167,14.7138,14.7112,14.7075,14.7062,14.7054,14.7042,14.7042,14.702,14.6996,14.6984,14.6983,14.6971,14.6963,14.6938,14.6921,14.6904,14.6888,14.6875,14.6866,14.6867,14.6854,14.6834,14.6821,14.6813,14.6796,14.6779,14.6763,14.6742,14.6742,14.6759,14.6754,14.6742,14.6742,14.6761,14.6767,14.6767,14.6733,14.6725,14.6725,14.6717,14.6717,14.6687,14.6671,14.6637,14.6609,14.6587,14.6575,14.6538,14.6521,14.65,14.6488,14.6467,14.6467,14.6479,14.6492,14.6496,14.6501,14.6537,14.6587,14.6629,14.665,14.6675,14.6675,14.6696,14.6737,14.6771,14.6775,14.6758,14.6775,14.6767,14.6771,14.6792,14.6792,14.6775,14.6775,14.6796,14.6813,14.6821,14.6846,14.6854,14.6871,14.6884,14.6875,14.6883,14.6875,14.6892,14.685,14.6858,14.6846,14.6813,14.68,14.68,14.6821,14.6846,14.6871,14.6888,14.69,14.69,14.6921,14.6938,14.6954,14.6967,14.7,14.7012,14.7054,14.7067,14.7067,14.7092,14.7067,14.7113,14.7129,14.7162,14.7179,14.7187,14.7204,14.7237,14.7309,14.7308,14.7325,14.7342,14.735,14.7358,14.7359,14.7367,14.7367,14.7375,14.7375,14.7384,14.7392,14.74,14.74,14.7392,14.7392,14.7383,14.7383,14.7375,14.7375,14.7367,14.7367,14.7325,14.7258,14.7259,14.72,14.7175,14.7175,14.715,14.715,14.7142,14.7142,14.7125,14.7092,14.7092,14.7117,14.7117,14.71,14.71,14.7084,14.7067,14.7067,14.705,14.705,14.7034,14.7033,14.6992,14.6984,14.6975,14.6867,14.6833,14.6809,14.6792,14.6783,14.6784,14.6675,14.6675,14.665,14.6642,14.6588,14.6579,14.6537,14.6529,14.6467,14.6467,14.6362,14.6354,14.6296,14.6246,14.6238,14.6188,14.6162,14.6146,14.6113,14.6088,14.6029,14.6004,14.5979,14.5963,14.5913,14.5871,14.5867,14.5846,14.5829,14.5804,14.5788,14.5754,14.5712,14.5654,14.5642,14.5608,14.5583,14.5542,14.5496,14.5487,14.5446,14.5429,14.5396,14.5379,14.5354,14.5304,14.5237,14.5204,14.5171,14.5137,14.5121,14.5113,14.5096,14.5071,14.5062,14.5029,14.4979,14.4969,14.4945,14.4912,14.4846,14.4813,14.4763,14.4738,14.4725,14.4725,14.4725,14.4717,14.4717,14.4679,14.4663,14.4629,14.4588,14.4579,14.4571,14.4554,14.4546,14.4521,14.4492,14.4474,14.4467,14.4467,14.4475,14.4475,14.4459,14.4442,14.4442,14.4458,14.4458,14.4392,14.4392,14.44,14.4358,14.4358,14.4367,14.4366,14.4333,14.4259,14.4237,14.4229,14.4134,14.4142,14.4096,14.4041,14.4034,14.4033,14.4025,14.3988,14.3904,14.387,14.3825,14.3825,14.3813,14.3779,14.3771,14.3729,14.3704,14.3696,14.3637,14.3571,14.3504,14.3495,14.348,14.3437,14.3413,14.3313,14.3304,14.3204,14.3196,14.3163,14.3154,14.3079,14.3071,14.3038,14.3029,14.3013,14.3004,14.2962,14.2937,14.2904,14.2896,14.2859,14.285,14.2848,14.2842,14.2859,14.2859,14.2833,14.2834,14.2762,14.2704,14.2588,14.2542,14.2517,14.2479,14.2454,14.2413,14.2321,14.2279,14.2259,14.2263,14.2273,14.2306,14.2321,14.2329,14.2358,14.2312,14.2304,14.2271,14.2263,14.2221,14.2212,14.218,14.2171,14.2087,14.2079,14.2013,14.1946,14.1879,14.1873,14.1862,14.1846,14.1842,14.1848,14.185,14.1817,14.1763,14.1713,14.1704,14.1671,14.1641,14.1641,14.1575,14.1542,14.1542,14.155,14.1592,14.1592,14.1571,14.1558,14.1558,14.1546,14.1521,14.1508,14.1509,14.1492,14.1492,14.15,14.15,14.1487,14.1475,14.1475,14.1454,14.1429,14.14,14.1367,14.1346,14.1338,14.1321,14.124,14.1397,14.14,14.1403,14.1411,14.1414,14.1403,14.1394,14.1387,14.138,14.1379,14.138,14.1385,14.1395,14.1403,14.1412,14.142,14.1424,14.142,14.1412,14.1407,14.1401,14.1403,14.141,14.1415,14.1432,14.1448,14.1458,14.1468,14.1483,14.1502,14.1518,14.1532,14.1546,14.1562,14.1576,14.1594,14.1609,14.1625,14.1635,14.1655,14.1674,14.1691,14.1698,14.1705,14.1708,14.171,14.1715,14.1729,14.1739,14.1754,14.1766,14.1789,14.1806,14.182,14.1825,14.1827,14.1827,14.1823,14.1829,14.1848,14.1864,14.1879,14.1889,14.1895,14.1915,14.1932,14.1952,14.1971,14.1988,14.2011,14.2035,14.2048,14.2064,14.208,14.2091,14.2104,14.2113,14.2111,14.2112,14.2126,14.2144,14.2156,14.2173,14.2189,14.2201,14.2209,14.2222,14.2231,14.2252,14.2273,14.2294,14.2317,14.2348,14.2376,14.2407,14.2432,14.2454,14.2478,14.2504,14.2524,14.2539,14.2552,14.2574,14.2599,14.2623,14.2638,14.2658,14.2669,14.2684,14.2699,14.2713,14.2731,14.2752,14.2772,14.2791,14.2807,14.2823,14.284,14.2861,14.2876,14.2894,14.2918,14.294,14.2963,14.2984,14.3004,14.3024,14.3052,14.3069,14.3097,14.3117,14.314,14.3166,14.3186,14.321,14.3237,14.3262,14.3287,14.3306,14.3332,14.3363,14.339,14.3419,14.3449,14.3473,14.3499,14.3523,14.3551,14.3584,14.3608,14.363,14.3647,14.368,14.3707,14.373,14.3752,14.3778,14.3805,14.3835,14.3859,14.3891,14.3918,14.395,14.398,14.4004,14.4028,14.4046,14.4067,14.4081,14.4098,14.4121,14.4142,14.416,14.418,14.4199,14.4222,14.423,14.4239,14.4251,14.4259,14.4256,14.4259,14.4257,14.4257,14.426,14.4261,14.4258,14.4259,14.4266,14.4262,14.4263,14.4266,14.4267,14.4271,14.4275,14.428,14.4287,14.43,14.4314,14.4335,14.435,14.4366,14.4383,14.4391,14.4406,14.4421,14.4437,14.4451,14.4472,14.448,14.449,14.4501,14.4506,14.4515,14.4523,14.4525,14.4529,14.4536,14.4545,14.4554,14.4561,14.4586,14.4601,14.4611,14.4624,14.4646,14.4663,14.4677,14.4698,14.4712,14.4735,14.475,14.4765,14.4776,14.4797,14.4814,14.4834,14.4855,14.4871,14.4891,14.4909,14.4921,14.494,14.4955,14.4963,14.4992,14.5001,14.5015,14.503,14.5044,14.5061,14.508,14.5096,14.5111,14.5132,14.5151,14.517,14.5184,14.5197,14.5203,14.5224,14.5233,14.5252,14.5276,14.529,14.5303,14.5316,14.5328,14.535,14.5373,14.5386,14.541,14.5432,14.5455,14.5466,14.5484,14.5503,14.5521,14.5543,14.5575,14.5585,14.5599,14.5615,14.5624,14.5638,14.5655,14.5677,14.5691,14.5704,14.5715,14.5728,14.5741,14.5753,14.5766,14.5784,14.5795,14.5811,14.5831,14.5845,14.5864,14.588,14.5892,14.5904,14.5917,14.5934,14.5948,14.5964,14.5987,14.6006,14.6024,14.6037,14.6043,14.6052,14.606,14.6057,14.606,14.606,14.606,14.6061,14.6061,14.606,14.6059,14.6058,14.6057,14.6057,14.6059,14.6059,14.6056,14.6053,14.6054,14.6051,14.6047,14.6046,14.6046,14.6045,14.6046,14.6043,14.6042,14.6038,14.6042,14.6042,14.604,14.6038,14.6036,14.6038,14.604,14.6044,14.6044,14.605,14.6053,14.6057,14.606,14.6067,14.6071,14.6072,14.6072,14.6079,14.6089,14.6093,14.6099,14.6117,14.6122,14.613,14.6142,14.615,14.6161,14.6168,14.6176,14.6182,14.6195,14.6201,14.6217,14.6231,14.6252,14.6258,14.6273,14.6281,14.6302,14.6314,14.6322,14.6332,14.6341,14.6352,14.6363,14.6371,14.6382,14.6422,14.6441,14.646,14.6487,14.6513,14.6536,14.6564,14.659,14.6619,14.6646,14.6669,14.6682,14.67,14.6717,14.6737,14.6744,14.6777,14.6802,14.6831,14.6849,14.6865,14.6898,14.6921,14.6945,14.6964,14.6997,14.7019,14.704,14.7061,14.7091,14.7125,14.7149,14.7175,14.7193,14.7207,14.7228,14.7258,14.7282,14.7312,14.7339,14.7363,14.7395,14.7408,14.7429,14.7452,14.747,14.749,14.7519,14.7543,14.7558,14.758,14.7616,14.764,14.7661,14.7682,14.7699,14.773,14.7749,14.7768,14.7777,14.7793,14.7925,14.7936,14.7946,14.7961,14.7978,14.7991,14.8005,14.8027,14.8044,14.8057,14.8073,14.812,14.8148,14.8173,14.8187,14.8206,14.822,14.823,14.8247,14.8264,14.8281,14.83,14.8319,14.8329,14.8344,14.8363,14.838,14.8398,14.8411,14.8449,14.8468,14.8479,14.8491,14.8512,14.8529,14.8553,14.8591,14.8612,14.8632,14.8651,14.8664,14.8676,14.8688,14.8701,14.8718,14.8729,14.8742,14.8759,14.8772,14.8784,14.8796,14.8805,14.8822,14.8844,14.8858,14.8865,14.8874,14.8887,14.8905,14.8922,14.8938,14.895,14.8963,14.8988,14.9022,14.905,14.9053,14.9068,14.9082,14.9093,14.9108,14.9128,14.9143,14.9155,14.9166,14.918,14.9199,14.9219,14.9235,14.9253,14.9273,14.9292,14.9335,14.9362,14.9402,14.9451,14.948,14.9497,14.9506,14.9532,14.9553,14.9566,14.9574,14.9586,14.96,14.9618,14.9626,14.9642,14.9657,14.9673,14.9688,14.9704,14.9708,14.9724,14.9737,14.9757,14.9777,14.9796,14.9811,14.9824,14.9842,14.9862,14.9883,14.9897,14.9918,14.9945,14.9955,14.9967,14.997,14.9976,14.9978,14.9982,14.9984,14.9993,15.0008,15.0026,15.0036,15.0061,15.0084,15.0109,15.0125,15.0136,15.0143,15.0148,15.0153,15.0162,15.0176,15.0194,15.0213,15.022,15.0222,15.0235,15.0238,15.0243,15.0243,15.0243,15.0244,15.0245,15.0246,15.0246,15.0244,15.0242,15.0239,15.0236,15.0234,15.023,15.0228,15.0226,15.0222,15.0218,15.0221,15.0219,15.0212,15.021,15.021,15.0207,15.0199,15.0192,15.0189,15.018,15.0171,15.0165,15.0154,15.0144,15.0137,15.0127,15.0119,15.0107,15.0094,15.0086,15.0067,15.0057,15.0048,15.0039,15.0027,15.0018,15.001,15.0009,14.9982,14.9975,14.9969,14.9963,14.9957,14.9954,14.9951,14.9947,14.9938,14.9934,14.9925,14.9921,14.9913,14.9902,14.9903,14.9898,14.9888,14.9879,14.987,14.986,14.985,14.9844,14.9836,14.982,14.9806,14.9799,14.9783,14.9771,14.9761,14.9751,14.9732,14.9722,14.971,14.9694,14.9677,14.9668,14.965,14.9633,14.962,14.9598,14.9575,14.956,14.9546,14.953,14.9528,14.9525,14.9519,14.952,14.9527,14.9535,14.9534,14.9544,14.9558,14.9562,14.9574,14.9579,14.9589,14.9592,14.9603,14.961,14.962,14.9625,14.9628,14.9629,14.9648,14.9662,14.9682,14.9695,14.9712,14.9732,14.9752,14.9767,14.9784,14.9803,14.9827,14.9844,14.9865,14.9887,14.9908,14.9934,14.9953,15.0006,15.0015,15.0058,15.0078,15.0092,15.0115,15.0118,15.0161,15.0182,15.0219,15.024,15.0273,15.0306,15.0337,15.037,15.0417,15.0454,15.0491,15.0561,15.0617,15.0663,15.0711,15.0766,15.0805,15.0863,15.0912,15.0934,15.0954,15.0971,15.1003,15.1014,15.1046,15.1068,15.1079,15.1113,15.1147,15.1162,15.118,15.1194,15.1219,15.1247,15.1278,15.131,15.1345,15.1382,15.1443,15.15,15.1537,15.1563,15.1607,15.163,15.1658,15.1677,15.1692,15.1712,15.1741,15.1771,15.1807,15.1831,15.1855,15.1883,15.1907,15.1941,15.1971,15.2006,15.2049,15.21,15.2149,15.2176,15.2201,15.2232,15.2247,15.2288,15.2313,15.2357,15.2392,15.2425,15.2467,15.2498,15.2536,15.257,15.2596,15.2618,15.2646,15.2672,15.269,15.2696,15.2703,15.2717,15.2722,15.2725,15.273,15.2728,15.2725,15.2723,15.2706,15.2691,15.2662,15.2634,15.2596,15.2548,15.2521,15.2509,15.2488,15.249,15.25,15.2513,15.2537,15.2561,15.2579,15.2609,15.2636,15.266,15.2689,15.2715,15.2743,15.2758,15.2769,15.2793,15.282,15.2847,15.2878,15.2925,15.2957,15.2994,15.3029,15.3056,15.3088,15.3105,15.3136,15.3167,15.3196,15.3213,15.3235,15.3242,15.3239,15.3233,15.3213,15.3192,15.3174,15.3164,15.3171,15.3188,15.3212,15.3233,15.3245,15.3255,15.3257,15.3252,15.3258,15.3276,15.3311,15.3339,15.3368,15.3406,15.3431,15.3455,15.3482,15.3511,15.3536,15.356,15.3577,15.3591,15.3591,15.3582,15.3569,15.356,15.3545,15.3518,15.3499,15.3476,15.3461,15.3434,15.3417,15.3393,15.3378,15.3361,15.3341,15.3332,15.3314,15.3296,15.3281,15.3267,15.325,15.3237,15.3224,15.3208,15.3186,15.3158,15.3147,15.3133,15.3119,15.31,15.3083,15.3066,15.3045,15.3022,15.3001,15.2974,15.2953,15.2935,15.2917,15.29,15.288,15.2868,15.2861,15.2849,15.2844,15.2839,15.2836,15.2828,15.281,15.2795,15.2783,15.2766,15.2743,15.2726,15.2713,15.2703,15.2696,15.2688,15.2686,15.2687,15.2693,15.2662,15.2619,15.2585,15.2524,15.2486,15.2456,15.2389,15.2355,15.2254,15.2224,15.2141,15.2111,15.1987,15.1959,15.1924,15.1886,15.1797,15.1792,15.1783,15.1714,15.1649,15.1573,15.1548,15.154,15.1532,15.1526,15.1496,15.1457,15.1563,15.1621,15.166,15.1692,15.1746,15.1785,15.1824,15.1891,15.1957,15.1985,15.2002,15.2034,15.2086,15.2115,15.2178,15.2226,15.2247,15.226,15.2271,15.2275,15.2283,15.2304,15.2357,15.2367,15.2405,15.2409,15.2426,15.2453,15.2461,15.2486,15.2522,15.2541,15.2561,15.2603,15.2625,15.2629,15.265,15.2667,15.2685,15.2693,15.2722,15.2737,15.2768,15.2796,15.282,15.2835,15.2892,15.2903,15.2904,15.291,15.2917,15.2939,15.2949,15.2948,15.2967]}]],[[{"lng":[-16.5668,-16.5692,-16.5696,-16.5684,-16.5675,-16.5663,-16.5668],"lat":[13.6258,13.6258,13.6263,13.6275,13.6275,13.6262,13.6258]},{"lng":[-16.6396,-16.6417,-16.643,-16.643,-16.6409,-16.6359,-16.6346,-16.6346,-16.6359,-16.6367,-16.6392,-16.6396],"lat":[13.6571,13.6575,13.6562,13.6537,13.6517,13.6525,13.6537,13.6563,13.6567,13.655,13.655,13.6571]},{"lng":[-16.668,-16.6679,-16.6688,-16.6688,-16.6701,-16.6709,-16.6719,-16.6721,-16.6721,-16.6713,-16.6713,-16.6642,-16.6592,-16.6542,-16.6484,-16.6442,-16.6436,-16.6421,-16.6421,-16.6434,-16.645,-16.6484,-16.6509,-16.653,-16.6538,-16.6567,-16.6609,-16.6622,-16.6621,-16.6634,-16.6667,-16.668],"lat":[13.6554,13.6596,13.6604,13.6621,13.6633,13.6634,13.6623,13.662,13.6554,13.6546,13.6529,13.6467,13.645,13.6384,13.6359,13.6358,13.6365,13.6379,13.6404,13.6417,13.6417,13.6391,13.6392,13.6413,13.6437,13.6467,13.6475,13.6487,13.6529,13.6542,13.6542,13.6554]},{"lng":[-16.6546,-16.6559,-16.6567,-16.6608,-16.6621,-16.6576,-16.6546],"lat":[13.678,13.6792,13.6783,13.6783,13.6771,13.6767,13.678]},{"lng":[-16.6346,-16.6346,-16.6313,-16.6313,-16.6321,-16.6355,-16.6355,-16.6355,-16.6355,-16.6359,-16.6371,-16.6371,-16.638,-16.638,-16.6388,-16.6388,-16.6396,-16.6405,-16.6413,-16.6413,-16.6401,-16.6384,-16.6346],"lat":[13.6629,13.6663,13.6695,13.6721,13.6729,13.6779,13.6835,13.6861,13.6912,13.6917,13.6913,13.6837,13.6829,13.6779,13.6771,13.6721,13.6712,13.6654,13.6646,13.6612,13.66,13.66,13.6629]},{"lng":[-16.6759,-16.6722,-16.6734,-16.6768,-16.6792,-16.6808,-16.6813,-16.6784,-16.6759],"lat":[13.6934,13.6937,13.695,13.695,13.6934,13.6933,13.6921,13.6917,13.6934]},{"lng":[-16.6526,-16.6517,-16.6463,-16.6455,-16.6455,-16.6475,-16.6484,-16.6526,-16.6538,-16.653,-16.6538,-16.6549,-16.6588,-16.6588,-16.6605,-16.6605,-16.6584,-16.6526],"lat":[13.7775,13.7775,13.7829,13.7846,13.7855,13.7875,13.7875,13.7839,13.7829,13.7821,13.7804,13.7793,13.7754,13.7746,13.773,13.7721,13.7717,13.7775]},{"lng":[-16.6609,-16.6563,-16.6547,-16.6521,-16.6521,-16.6546,-16.6546,-16.6567,-16.658,-16.658,-16.6588,-16.6596,-16.6621,-16.6621,-16.6638,-16.6638,-16.6646,-16.6634,-16.6609],"lat":[13.7758,13.7804,13.7846,13.7871,13.7887,13.7904,13.7912,13.7934,13.7921,13.7904,13.7895,13.7845,13.7804,13.7795,13.7779,13.7771,13.7763,13.775,13.7758]},{"lng":[-16.6696,-16.6709,-16.6733,-16.6746,-16.6746,-16.6746,-16.6734,-16.6696],"lat":[13.7979,13.8025,13.8017,13.8004,13.7954,13.7945,13.7942,13.7979]},{"lng":[-16.65,-16.6472,-16.643,-16.6421,-16.643,-16.643,-16.645,-16.6463,-16.6555,-16.6555,-16.6563,-16.6525,-16.65],"lat":[13.7925,13.7962,13.8046,13.8055,13.807,13.8096,13.8108,13.8071,13.7979,13.7971,13.7962,13.7917,13.7925]},{"lng":[-16.638,-16.6363,-16.6346,-16.635,-16.6367,-16.6409,-16.6422,-16.6405,-16.6404,-16.6415,-16.6421,-16.6409,-16.638],"lat":[13.8054,13.8088,13.8104,13.8142,13.8134,13.8133,13.812,13.8087,13.8063,13.8041,13.8029,13.8025,13.8054]},{"lng":[-16.6117,-16.6076,-16.6051,-16.5992,-16.5967,-16.5913,-16.5905,-16.5905,-16.5951,-16.6051,-16.6084,-16.6109,-16.6134,-16.6176,-16.6209,-16.6234,-16.6276,-16.6313,-16.6321,-16.6338,-16.6363,-16.6384,-16.6418,-16.6446,-16.6446,-16.6463,-16.6463,-16.6426,-16.64,-16.6392,-16.6351,-16.6343,-16.6292,-16.6284,-16.6267,-16.6242,-16.6209,-16.6151,-16.6142,-16.6117],"lat":[13.8,13.8025,13.8034,13.81,13.8108,13.8171,13.8188,13.8229,13.8283,13.8283,13.8267,13.8267,13.825,13.8242,13.8225,13.8225,13.82,13.8154,13.8121,13.8096,13.8038,13.8017,13.8009,13.797,13.7946,13.7929,13.7888,13.785,13.785,13.7859,13.7867,13.7875,13.7875,13.7883,13.7884,13.7908,13.7925,13.7975,13.7975,13.8]},{"lng":[-16.7346,-16.7363,-16.7363,-16.7371,-16.7371,-16.7371,-16.7384,-16.7405,-16.7405,-16.7396,-16.738,-16.7351,-16.7346],"lat":[13.8413,13.8429,13.8462,13.8471,13.85,13.8554,13.8567,13.8554,13.8504,13.8496,13.8429,13.84,13.8413]},{"lng":[-16.7405,-16.7405,-16.7409,-16.7426,-16.743,-16.743,-16.7425,-16.7417,-16.7405],"lat":[13.9062,13.9071,13.9075,13.9075,13.9071,13.9054,13.905,13.905,13.9062]},{"lng":[-16.7463,-16.7463,-16.7476,-16.7497,-16.7496,-16.7505,-16.7505,-16.7492,-16.7463],"lat":[13.9196,13.9229,13.9242,13.9229,13.9187,13.9179,13.9154,13.915,13.9196]},{"lng":[-15.8028,-15.8034,-15.8104,-15.8161,-15.8201,-15.8234,-15.8264,-15.8305,-15.8342,-15.8385,-15.842,-15.8467,-15.8515,-15.855,-15.8588,-15.862,-15.8646,-15.8651,-15.8644,-15.8633,-15.8626,-15.8623,-15.8626,-15.8635,-15.8686,-15.877,-15.8837,-15.8915,-15.8942,-15.8962,-15.8976,-15.8992,-15.9009,-15.9026,-15.9048,-15.9075,-15.9103,-15.9133,-15.9147,-15.9157,-15.9163,-15.9171,-15.9196,-15.9224,-15.9244,-15.9276,-15.9331,-15.9369,-15.9434,-15.9479,-15.9532,-15.9575,-15.963,-15.9688,-15.972,-15.9744,-15.9777,-15.9811,-15.9854,-15.9873,-15.9884,-15.9891,-15.9899,-15.9904,-15.9911,-15.9916,-15.9922,-15.9929,-15.9936,-15.9944,-15.9954,-15.9962,-15.9968,-15.9973,-15.999,-16.0001,-16.0018,-16.004,-16.0058,-16.0078,-16.0099,-16.0139,-16.0156,-16.0176,-16.0264,-16.0321,-16.0348,-16.0361,-16.0377,-16.0408,-16.0439,-16.0469,-16.0495,-16.0521,-16.0545,-16.0577,-16.0602,-16.0658,-16.0681,-16.0703,-16.0728,-16.0759,-16.0794,-16.0816,-16.0843,-16.0863,-16.0879,-16.0905,-16.0921,-16.0938,-16.0948,-16.0971,-16.0992,-16.1025,-16.1049,-16.1074,-16.1095,-16.1114,-16.1158,-16.1177,-16.1213,-16.1235,-16.1284,-16.1317,-16.1345,-16.1368,-16.1395,-16.1422,-16.147,-16.1486,-16.1522,-16.1544,-16.1577,-16.1605,-16.1638,-16.1674,-16.1705,-16.1735,-16.176,-16.177,-16.1782,-16.1817,-16.1854,-16.188,-16.1904,-16.1932,-16.1945,-16.1966,-16.1988,-16.2014,-16.2036,-16.2068,-16.2101,-16.2123,-16.2153,-16.2183,-16.2206,-16.2237,-16.2253,-16.2275,-16.2294,-16.2308,-16.2311,-16.2316,-16.2322,-16.232,-16.2329,-16.2345,-16.2357,-16.237,-16.2396,-16.2408,-16.2426,-16.2454,-16.2472,-16.2483,-16.2502,-16.253,-16.255,-16.2579,-16.2601,-16.2626,-16.2654,-16.2683,-16.2708,-16.2724,-16.2746,-16.2769,-16.2786,-16.281,-16.2841,-16.2867,-16.2885,-16.2911,-16.2939,-16.2976,-16.2988,-16.3008,-16.3029,-16.3049,-16.3081,-16.3099,-16.313,-16.3188,-16.3198,-16.3225,-16.3254,-16.328,-16.3301,-16.3319,-16.3346,-16.3377,-16.3395,-16.3427,-16.3462,-16.3485,-16.3512,-16.3556,-16.3592,-16.3632,-16.3683,-16.3709,-16.3736,-16.3764,-16.3788,-16.3833,-16.3858,-16.3879,-16.3901,-16.3905,-16.394,-16.3988,-16.4022,-16.406,-16.4084,-16.4107,-16.4138,-16.4169,-16.4198,-16.4232,-16.4263,-16.4291,-16.4314,-16.4334,-16.4353,-16.438,-16.4405,-16.4424,-16.4454,-16.4469,-16.4521,-16.4553,-16.4591,-16.4636,-16.4677,-16.4713,-16.4732,-16.4742,-16.4755,-16.4768,-16.4789,-16.4815,-16.4858,-16.4889,-16.4921,-16.4939,-16.4965,-16.4961,-16.4981,-16.5002,-16.5029,-16.5066,-16.5102,-16.5135,-16.5157,-16.5202,-16.527,-16.5323,-16.5377,-16.5438,-16.5519,-16.5576,-16.5608,-16.5645,-16.5738,-16.5799,-16.5871,-16.5913,-16.5963,-16.6015,-16.607,-16.6143,-16.6211,-16.6288,-16.6349,-16.6402,-16.6444,-16.6479,-16.6509,-16.6522,-16.6542,-16.658,-16.6615,-16.6628,-16.6652,-16.6685,-16.6734,-16.6773,-16.6802,-16.6841,-16.6855,-16.6833,-16.6813,-16.6793,-16.6769,-16.6741,-16.6729,-16.6717,-16.6707,-16.6695,-16.6678,-16.6663,-16.6652,-16.6648,-16.6633,-16.6615,-16.6608,-16.6594,-16.6583,-16.6578,-16.6568,-16.656,-16.6552,-16.6545,-16.6534,-16.6524,-16.6516,-16.6495,-16.6468,-16.6451,-16.6433,-16.6419,-16.6401,-16.6384,-16.6369,-16.636,-16.6346,-16.6331,-16.6319,-16.631,-16.6304,-16.6293,-16.6285,-16.6271,-16.6267,-16.6261,-16.6256,-16.6248,-16.6248,-16.6242,-16.6242,-16.6233,-16.6237,-16.6247,-16.6276,-16.6293,-16.632,-16.6339,-16.6353,-16.6374,-16.6395,-16.6411,-16.6428,-16.6437,-16.6445,-16.6453,-16.6468,-16.6474,-16.6486,-16.6484,-16.6451,-16.6429,-16.6408,-16.6383,-16.6356,-16.6331,-16.6306,-16.6285,-16.6264,-16.6244,-16.6222,-16.6158,-16.6139,-16.613,-16.6109,-16.6093,-16.608,-16.6075,-16.607,-16.6068,-16.6064,-16.6062,-16.6058,-16.6058,-16.6055,-16.6054,-16.6053,-16.6051,-16.6048,-16.6039,-16.6039,-16.6042,-16.6044,-16.605,-16.605,-16.6053,-16.6053,-16.6056,-16.6055,-16.6054,-16.605,-16.6051,-16.6055,-16.6051,-16.6047,-16.6047,-16.6044,-16.6046,-16.6047,-16.6045,-16.6047,-16.6053,-16.6056,-16.6064,-16.6075,-16.6088,-16.6104,-16.6108,-16.6123,-16.6142,-16.6156,-16.6169,-16.6185,-16.6196,-16.6201,-16.6206,-16.6208,-16.6212,-16.6211,-16.6214,-16.6209,-16.6206,-16.6192,-16.6171,-16.6163,-16.615,-16.6134,-16.6118,-16.6102,-16.609,-16.6073,-16.6056,-16.6043,-16.6022,-16.6011,-16.5994,-16.5982,-16.5969,-16.5958,-16.5948,-16.594,-16.5926,-16.5914,-16.5905,-16.5895,-16.5883,-16.587,-16.5862,-16.5849,-16.5845,-16.5829,-16.5819,-16.5812,-16.5802,-16.5791,-16.5776,-16.5768,-16.575,-16.5737,-16.5727,-16.5701,-16.5694,-16.5682,-16.5667,-16.5654,-16.5647,-16.563,-16.562,-16.5608,-16.5591,-16.5568,-16.5562,-16.5546,-16.5533,-16.5516,-16.55,-16.5481,-16.5467,-16.5452,-16.5436,-16.5421,-16.5403,-16.5386,-16.5374,-16.5362,-16.5348,-16.5339,-16.5322,-16.5311,-16.5303,-16.5296,-16.5285,-16.5278,-16.5259,-16.5241,-16.5241,-16.5234,-16.5229,-16.5223,-16.5215,-16.5207,-16.5199,-16.5193,-16.5189,-16.5184,-16.5181,-16.518,-16.5178,-16.5177,-16.5176,-16.5176,-16.5178,-16.5177,-16.518,-16.5186,-16.519,-16.5198,-16.5206,-16.5212,-16.5217,-16.5225,-16.5234,-16.5246,-16.5257,-16.5267,-16.5273,-16.5273,-16.5271,-16.5269,-16.5268,-16.5279,-16.5285,-16.5286,-16.5287,-16.5294,-16.5299,-16.5303,-16.5308,-16.5306,-16.531,-16.5314,-16.5317,-16.5322,-16.533,-16.5335,-16.5346,-16.5349,-16.5348,-16.5358,-16.5361,-16.5361,-16.537,-16.5382,-16.5387,-16.5394,-16.5405,-16.5412,-16.5426,-16.5435,-16.5443,-16.5459,-16.5474,-16.5484,-16.5493,-16.5495,-16.5506,-16.5512,-16.5526,-16.5545,-16.5574,-16.559,-16.5611,-16.5635,-16.5653,-16.5681,-16.5694,-16.5711,-16.5735,-16.5753,-16.5771,-16.5794,-16.5814,-16.5841,-16.5867,-16.589,-16.5906,-16.592,-16.595,-16.5972,-16.5997,-16.6017,-16.6037,-16.6057,-16.6073,-16.6089,-16.6106,-16.6119,-16.614,-16.6155,-16.6175,-16.6202,-16.6217,-16.6234,-16.625,-16.6259,-16.6263,-16.628,-16.6297,-16.6308,-16.6327,-16.6336,-16.6352,-16.6377,-16.6389,-16.6404,-16.6424,-16.6439,-16.6457,-16.6466,-16.6492,-16.6516,-16.6541,-16.6568,-16.6596,-16.6625,-16.6659,-16.6684,-16.671,-16.6733,-16.6755,-16.6777,-16.6795,-16.6811,-16.6832,-16.685,-16.6869,-16.6888,-16.6906,-16.6928,-16.6948,-16.6963,-16.6981,-16.7002,-16.702,-16.7038,-16.7052,-16.7054,-16.7053,-16.7061,-16.7064,-16.7069,-16.7065,-16.7068,-16.7072,-16.7077,-16.7079,-16.7083,-16.7079,-16.7079,-16.7088,-16.7091,-16.7104,-16.7105,-16.7102,-16.71,-16.7098,-16.7103,-16.7102,-16.7097,-16.7092,-16.7089,-16.7089,-16.7079,-16.7078,-16.7079,-16.7074,-16.7075,-16.7071,-16.7073,-16.7078,-16.7085,-16.7088,-16.7096,-16.7109,-16.7124,-16.7138,-16.715,-16.7166,-16.7178,-16.7187,-16.7198,-16.7208,-16.7219,-16.7232,-16.7246,-16.7264,-16.7277,-16.729,-16.7306,-16.7318,-16.7333,-16.7344,-16.7362,-16.7372,-16.739,-16.7409,-16.7425,-16.7437,-16.7452,-16.7463,-16.7472,-16.7475,-16.7484,-16.7488,-16.7488,-16.7486,-16.7474,-16.7476,-16.746,-16.7457,-16.7441,-16.7433,-16.7431,-16.7423,-16.7422,-16.7401,-16.7387,-16.7368,-16.7348,-16.7322,-16.7312,-16.7303,-16.7281,-16.7254,-16.7239,-16.7232,-16.7228,-16.7222,-16.7223,-16.7235,-16.7241,-16.7248,-16.7256,-16.7265,-16.7278,-16.7291,-16.7304,-16.7321,-16.7329,-16.7341,-16.7356,-16.7375,-16.7392,-16.7402,-16.7413,-16.7428,-16.7437,-16.745,-16.7456,-16.7468,-16.7479,-16.7493,-16.7507,-16.7523,-16.7539,-16.7557,-16.7574,-16.7588,-16.7601,-16.7612,-16.7624,-16.7635,-16.7634,-16.763,-16.7617,-16.7597,-16.7574,-16.7565,-16.7559,-16.7551,-16.7535,-16.7533,-16.7527,-16.7526,-16.7533,-16.7537,-16.7548,-16.7559,-16.7573,-16.7585,-16.7593,-16.7615,-16.7625,-16.7633,-16.7641,-16.7655,-16.7677,-16.7693,-16.771,-16.7725,-16.7737,-16.7763,-16.7766,-16.778,-16.7798,-16.7818,-16.7835,-16.7849,-16.7998,-16.7996,-16.7996,-16.7946,-16.7938,-16.793,-16.793,-16.7913,-16.7913,-16.788,-16.7863,-16.7855,-16.7855,-16.7846,-16.7846,-16.7838,-16.7838,-16.7813,-16.778,-16.778,-16.778,-16.7771,-16.7772,-16.7763,-16.7763,-16.7755,-16.7755,-16.7746,-16.7746,-16.7738,-16.7738,-16.7746,-16.7746,-16.7738,-16.7738,-16.7721,-16.773,-16.7713,-16.7713,-16.7705,-16.7705,-16.7696,-16.7688,-16.7688,-16.768,-16.768,-16.7688,-16.7688,-16.768,-16.768,-16.7671,-16.7671,-16.768,-16.768,-16.7671,-16.768,-16.768,-16.7688,-16.768,-16.768,-16.7671,-16.7671,-16.7663,-16.7663,-16.7663,-16.7655,-16.7655,-16.7646,-16.7646,-16.7655,-16.7655,-16.7646,-16.7638,-16.763,-16.763,-16.7613,-16.7613,-16.76,-16.7584,-16.7563,-16.7563,-16.7571,-16.758,-16.7588,-16.7588,-16.758,-16.758,-16.7555,-16.7555,-16.7563,-16.7571,-16.7563,-16.7563,-16.7571,-16.7571,-16.7588,-16.7588,-16.7596,-16.763,-16.7629,-16.7638,-16.763,-16.763,-16.7621,-16.7613,-16.7605,-16.7605,-16.7596,-16.7596,-16.7605,-16.7604,-16.7613,-16.7613,-16.7588,-16.7588,-16.7571,-16.7563,-16.7546,-16.7538,-16.7509,-16.7484,-16.7451,-16.7425,-16.7418,-16.7384,-16.7317,-16.7284,-16.723,-16.7221,-16.7221,-16.723,-16.7238,-16.7238,-16.7238,-16.7255,-16.7226,-16.7192,-16.7176,-16.7155,-16.7155,-16.7146,-16.7125,-16.7109,-16.7105,-16.7092,-16.7059,-16.7025,-16.7,-16.6959,-16.6926,-16.69,-16.6884,-16.688,-16.6905,-16.6942,-16.7009,-16.7047,-16.7059,-16.7109,-16.7142,-16.7167,-16.7205,-16.7213,-16.7213,-16.7205,-16.7184,-16.7159,-16.7125,-16.7092,-16.7084,-16.7059,-16.7034,-16.6976,-16.6834,-16.6801,-16.6796,-16.6805,-16.6776,-16.6767,-16.6792,-16.6834,-16.6867,-16.69,-16.6926,-16.695,-16.7001,-16.7029,-16.7021,-16.7021,-16.7121,-16.7129,-16.7155,-16.7155,-16.7159,-16.7188,-16.7193,-16.7222,-16.7238,-16.7267,-16.7276,-16.7304,-16.7304,-16.7321,-16.743,-16.743,-16.7421,-16.7421,-16.7405,-16.738,-16.738,-16.7346,-16.7346,-16.7355,-16.7396,-16.7405,-16.7438,-16.7425,-16.7367,-16.7359,-16.7251,-16.7242,-16.7184,-16.7176,-16.7258,-16.7317,-16.7363,-16.7363,-16.7355,-16.7355,-16.7338,-16.7338,-16.733,-16.7338,-16.7338,-16.7355,-16.7355,-16.7325,-16.7292,-16.7284,-16.7251,-16.7242,-16.7221,-16.7221,-16.7234,-16.725,-16.728,-16.7271,-16.7271,-16.7271,-16.7288,-16.7288,-16.73,-16.7305,-16.7296,-16.7288,-16.7309,-16.7326,-16.7351,-16.7355,-16.7346,-16.7346,-16.7338,-16.7338,-16.7322,-16.7305,-16.7313,-16.7313,-16.7346,-16.7346,-16.7321,-16.7321,-16.7313,-16.7313,-16.728,-16.7246,-16.7213,-16.7167,-16.7142,-16.7117,-16.7101,-16.7063,-16.7013,-16.7005,-16.6963,-16.6934,-16.6917,-16.6892,-16.6884,-16.6851,-16.6834,-16.6813,-16.6805,-16.6763,-16.6767,-16.68,-16.6825,-16.6843,-16.6859,-16.6888,-16.6855,-16.6855,-16.6838,-16.6834,-16.6826,-16.6801,-16.6784,-16.6771,-16.6763,-16.6771,-16.6772,-16.6788,-16.6759,-16.6734,-16.6717,-16.6684,-16.6679,-16.6688,-16.6688,-16.668,-16.6675,-16.6659,-16.6621,-16.6621,-16.6567,-16.655,-16.6517,-16.6492,-16.6467,-16.6455,-16.6455,-16.648,-16.648,-16.6505,-16.6513,-16.6521,-16.6521,-16.6509,-16.6488,-16.6488,-16.648,-16.6471,-16.6442,-16.6434,-16.6421,-16.6418,-16.6405,-16.64,-16.6384,-16.6367,-16.6359,-16.6342,-16.6309,-16.6275,-16.6226,-16.6201,-16.6142,-16.6134,-16.6117,-16.6092,-16.6059,-16.6034,-16.6009,-16.5959,-16.5951,-16.5934,-16.5909,-16.5876,-16.5851,-16.5825,-16.5809,-16.5767,-16.5742,-16.5692,-16.568,-16.568,-16.5688,-16.5646,-16.5638,-16.5621,-16.5621,-16.5538,-16.553,-16.5509,-16.5488,-16.5488,-16.5472,-16.5471,-16.5413,-16.5421,-16.5421,-16.5446,-16.5522,-16.5521,-16.5492,-16.5475,-16.5467,-16.54,-16.5396,-16.5413,-16.5413,-16.5401,-16.5367,-16.5367,-16.5384,-16.5396,-16.5388,-16.5388,-16.5367,-16.5334,-16.5326,-16.53,-16.5292,-16.5267,-16.5259,-16.5217,-16.52,-16.5184,-16.5146,-16.513,-16.5168,-16.5201,-16.5259,-16.5267,-16.53,-16.5309,-16.5326,-16.5359,-16.5367,-16.5392,-16.5401,-16.5425,-16.5459,-16.5476,-16.5488,-16.5488,-16.5451,-16.5442,-16.5407,-16.5351,-16.5334,-16.5322,-16.5321,-16.5338,-16.5338,-16.5359,-16.5392,-16.5413,-16.5413,-16.5401,-16.5392,-16.5375,-16.5371,-16.5384,-16.54,-16.5409,-16.5433,-16.5467,-16.5476,-16.5534,-16.5542,-16.5567,-16.5575,-16.5609,-16.5617,-16.5651,-16.5659,-16.5676,-16.5692,-16.5725,-16.5746,-16.5746,-16.5717,-16.5705,-16.5767,-16.5817,-16.5851,-16.5867,-16.5938,-16.5942,-16.5984,-16.6013,-16.6013,-16.5996,-16.6009,-16.6017,-16.6051,-16.6067,-16.6113,-16.6155,-16.623,-16.6221,-16.6226,-16.6242,-16.6263,-16.6263,-16.6321,-16.6313,-16.6292,-16.6233,-16.62,-16.615,-16.6142,-16.6109,-16.61,-16.6059,-16.6051,-16.6043,-16.6084,-16.6093,-16.6109,-16.6142,-16.6168,-16.6175,-16.6259,-16.6267,-16.6301,-16.633,-16.6371,-16.6372,-16.633,-16.633,-16.6305,-16.6305,-16.6288,-16.6271,-16.6204,-16.6205,-16.6196,-16.6196,-16.6171,-16.6155,-16.6154,-16.6163,-16.6163,-16.6171,-16.6171,-16.6163,-16.6163,-16.6163,-16.6155,-16.6146,-16.6146,-16.6134,-16.61,-16.6076,-16.6063,-16.6063,-16.6068,-16.6092,-16.6109,-16.613,-16.613,-16.6101,-16.6092,-16.6076,-16.6042,-16.6009,-16.5976,-16.5875,-16.5867,-16.5851,-16.5742,-16.57,-16.5676,-16.5671,-16.5646,-16.5684,-16.5709,-16.5742,-16.5746,-16.5788,-16.5813,-16.5813,-16.5834,-16.5859,-16.5859,-16.5809,-16.5775,-16.5742,-16.5726,-16.5717,-16.5642,-16.5605,-16.5605,-16.5609,-16.5634,-16.5651,-16.5659,-16.5709,-16.5734,-16.575,-16.5771,-16.5763,-16.5763,-16.568,-16.568,-16.5638,-16.5592,-16.5591,-16.5584,-16.5576,-16.5567,-16.5558,-16.5527,-16.5511,-16.5492,-16.5477,-16.5464,-16.5457,-16.5453,-16.5462,-16.5463,-16.5462,-16.5455,-16.545,-16.5447,-16.5448,-16.5449,-16.5425,-16.5248,-16.5073,-16.4992,-16.4895,-16.4717,-16.4539,-16.4361,-16.4289,-16.4186,-16.4184,-16.4009,-16.3831,-16.365,-16.3475,-16.3298,-16.312,-16.2945,-16.2764,-16.2586,-16.2411,-16.2234,-16.2059,-16.1878,-16.17,-16.1523,-16.1348,-16.117,-16.1033,-16.0992,-16.0878,-16.0814,-16.0636,-16.0461,-16.0284,-16.0103,-16,-15.9992,-15.982,-15.9642,-15.9467,-15.9289,-15.9111,-15.8934,-15.8756,-15.8578,-15.8403,-15.8225,-15.8048,-15.787,-15.7692,-15.7514,-15.7339,-15.7159,-15.6981,-15.6806,-15.6628,-15.645,-15.6273,-15.6095,-15.6036,-15.5917,-15.5739,-15.5564,-15.5386,-15.5206,-15.5031,-15.4992,-15.4881,-15.4867,-15.4878,-15.4879,-15.4881,-15.4878,-15.4873,-15.4861,-15.4848,-15.4845,-15.4834,-15.4834,-15.4814,-15.4795,-15.4775,-15.4764,-15.4745,-15.4725,-15.4706,-15.4686,-15.4659,-15.4639,-15.4623,-15.4598,-15.4578,-15.4542,-15.4517,-15.4481,-15.4448,-15.4417,-15.4381,-15.4348,-15.4306,-15.4273,-15.4231,-15.4195,-15.4148,-15.4147,-15.4106,-15.4056,-15.4017,-15.397,-15.3928,-15.3886,-15.3853,-15.3811,-15.3761,-15.372,-15.3678,-15.3636,-15.3589,-15.3542,-15.3492,-15.3442,-15.3389,-15.3331,-15.3275,-15.3206,-15.3139,-15.3084,-15.3031,-15.2986,-15.2931,-15.2886,-15.285,-15.2814,-15.277,-15.2731,-15.2703,-15.2664,-15.2636,-15.2603,-15.2575,-15.2545,-15.2514,-15.2506,-15.2334,-15.2159,-15.1986,-15.1814,-15.1642,-15.1636,-15.1459,-15.1284,-15.1111,-15.0934,-15.0759,-15.0695,-15.0514,-15.0336,-15.0156,-14.9992,-14.9975,-14.9909,-14.9861,-14.9811,-14.9756,-14.97,-14.9642,-14.9586,-14.9523,-14.9453,-14.9389,-14.9323,-14.9259,-14.9195,-14.9136,-14.9081,-14.9023,-14.8973,-14.8923,-14.8873,-14.8823,-14.8778,-14.8734,-14.8684,-14.8631,-14.8595,-14.8556,-14.855,-14.8511,-14.847,-14.8431,-14.8407,-14.8234,-14.8213,-14.82,-14.8174,-14.8158,-14.8139,-14.8125,-14.8107,-14.8093,-14.8082,-14.8065,-14.8047,-14.8023,-14.8016,-14.7993,-14.7961,-14.7906,-14.7845,-14.7792,-14.7757,-14.7706,-14.7643,-14.7593,-14.754,-14.7497,-14.7441,-14.7362,-14.7303,-14.7239,-14.7201,-14.7186,-14.7175,-14.7161,-14.7145,-14.7135,-14.7121,-14.7105,-14.7095,-14.7092,-14.7069,-14.7042,-14.7016,-14.7009,-14.6999,-14.6973,-14.6943,-14.6919,-14.687,-14.682,-14.679,-14.675,-14.6709,-14.6673,-14.6631,-14.6587,-14.655,-14.6516,-14.6489,-14.6452,-14.6422,-14.6379,-14.6338,-14.6309,-14.6251,-14.6225,-14.6186,-14.6123,-14.609,-14.6053,-14.6043,-14.6033,-14.6013,-14.6005,-14.5992,-14.5985,-14.5972,-14.5959,-14.5945,-14.5934,-14.5922,-14.5907,-14.59,-14.5889,-14.5879,-14.5871,-14.5862,-14.5856,-14.5853,-14.5843,-14.5835,-14.5828,-14.5822,-14.5818,-14.5815,-14.5818,-14.5821,-14.5823,-14.5826,-14.5829,-14.5833,-14.5842,-14.5845,-14.5847,-14.5848,-14.5849,-14.5853,-14.5857,-14.5863,-14.5862,-14.5869,-14.5867,-14.5873,-14.5878,-14.5879,-14.5882,-14.5888,-14.5893,-14.5894,-14.5901,-14.5899,-14.5904,-14.5905,-14.5913,-14.5913,-14.5921,-14.5921,-14.5923,-14.5924,-14.5928,-14.5931,-14.5933,-14.5935,-14.593,-14.5926,-14.592,-14.592,-14.5914,-14.5912,-14.5907,-14.5907,-14.5899,-14.5893,-14.5894,-14.5889,-14.5884,-14.588,-14.5876,-14.587,-14.5867,-14.5864,-14.5863,-14.5859,-14.5852,-14.5848,-14.5845,-14.5846,-14.5847,-14.5845,-14.5846,-14.5848,-14.5854,-14.5855,-14.586,-14.5867,-14.5868,-14.5866,-14.5868,-14.5868,-14.5871,-14.5876,-14.5881,-14.5884,-14.5882,-14.5884,-14.5886,-14.5889,-14.5891,-14.5889,-14.5894,-14.5898,-14.5901,-14.5901,-14.5901,-14.5901,-14.59,-14.5903,-14.5908,-14.5909,-14.5915,-14.5912,-14.5918,-14.5919,-14.592,-14.5924,-14.5932,-14.5933,-14.5932,-14.5935,-14.5938,-14.5941,-14.5946,-14.5946,-14.595,-14.5953,-14.5954,-14.5954,-14.5961,-14.5961,-14.5975,-14.6003,-14.6029,-14.6055,-14.6085,-14.611,-14.6131,-14.6161,-14.6201,-14.6235,-14.6267,-14.6292,-14.6304,-14.6317,-14.6361,-14.6376,-14.6393,-14.6403,-14.6388,-14.6373,-14.6358,-14.6343,-14.6328,-14.6314,-14.6299,-14.6284,-14.6277,-14.6269,-14.6251,-14.6249,-14.6326,-14.6428,-14.6485,-14.6528,-14.6571,-14.6577,-14.6581,-14.6587,-14.659,-14.6628,-14.6634,-14.6642,-14.683,-14.6888,-14.6926,-14.7051,-14.7056,-14.7191,-14.7191,-14.7277,-14.7279,-14.7315,-14.7323,-14.7345,-14.7379,-14.7462,-14.7534,-14.7608,-14.7693,-14.7703,-14.7711,-14.7733,-14.7734,-14.7822,-14.7822,-14.7922,-14.796,-14.7981,-14.8002,-14.8027,-14.8028,-14.8036,-14.8111,-14.8203,-14.836,-14.8496,-14.8686,-14.8891,-14.9095,-14.9362,-14.9629,-14.9896,-14.9989,-15.0162,-15.0229,-15.0515,-15.0802,-15.0855,-15.0887,-15.1088,-15.1375,-15.1661,-15.1947,-15.2234,-15.252,-15.2807,-15.3093,-15.3366,-15.3639,-15.3912,-15.4185,-15.4458,-15.4454,-15.4446,-15.4435,-15.443,-15.4423,-15.4419,-15.4416,-15.4408,-15.4405,-15.4402,-15.4395,-15.4392,-15.4378,-15.4377,-15.4371,-15.437,-15.436,-15.4361,-15.4353,-15.4358,-15.4357,-15.4357,-15.4371,-15.4397,-15.4425,-15.4449,-15.4493,-15.4532,-15.4569,-15.4594,-15.4617,-15.4646,-15.4672,-15.4686,-15.471,-15.4731,-15.4759,-15.4786,-15.4817,-15.4844,-15.4864,-15.4884,-15.4904,-15.4923,-15.4953,-15.4973,-15.499,-15.501,-15.503,-15.5046,-15.5054,-15.5078,-15.5103,-15.513,-15.5148,-15.5175,-15.5204,-15.5238,-15.5276,-15.5313,-15.5348,-15.5364,-15.5386,-15.5408,-15.5432,-15.5452,-15.5496,-15.5523,-15.5562,-15.5597,-15.563,-15.5665,-15.5698,-15.5736,-15.5757,-15.5773,-15.5785,-15.5815,-15.5846,-15.5868,-15.5909,-15.5937,-15.5962,-15.5991,-15.6009,-15.603,-15.6064,-15.6078,-15.6099,-15.6122,-15.6148,-15.6185,-15.6218,-15.6251,-15.6272,-15.6307,-15.6322,-15.6351,-15.6379,-15.6408,-15.6434,-15.6457,-15.6484,-15.6509,-15.6544,-15.6565,-15.6587,-15.6608,-15.6624,-15.6649,-15.6666,-15.6689,-15.6711,-15.6726,-15.673,-15.6735,-15.6739,-15.6744,-15.6749,-15.6752,-15.6755,-15.676,-15.6773,-15.6841,-15.6903,-15.692,-15.6876,-15.6751,-15.6621,-15.6525,-15.6406,-15.6354,-15.6356,-15.639,-15.6435,-15.6498,-15.6554,-15.6628,-15.6664,-15.6722,-15.6791,-15.687,-15.6946,-15.7003,-15.7066,-15.7117,-15.7158,-15.7223,-15.7274,-15.7322,-15.7369,-15.7413,-15.7457,-15.7508,-15.7552,-15.7584,-15.7603,-15.7617,-15.7633,-15.7659,-15.7678,-15.7699,-15.7715,-15.772,-15.7727,-15.7729,-15.7737,-15.775,-15.7765,-15.7779,-15.7798,-15.7815,-15.7829,-15.7851,-15.7868,-15.7892,-15.793,-15.7962,-15.7989,-15.8002,-15.8018,-15.8023,-15.8028],"lat":[14.9801,14.9979,14.9993,15.0009,15.0033,15.0071,15.0123,15.0183,15.0234,15.0287,15.0327,15.0389,15.0441,15.0485,15.0511,15.0512,15.0494,15.0455,15.0417,15.0378,15.0334,15.025,15.02,15.0165,15.0109,15.0064,15.0027,14.9964,14.9934,14.9908,14.9887,14.9867,14.9849,14.9835,14.9826,14.9816,14.9802,14.9792,14.9788,14.9786,14.9788,14.9785,14.9776,14.9771,14.9768,14.9766,14.9761,14.9756,14.9749,14.9746,14.9741,14.9736,14.9731,14.9726,14.9722,14.9722,14.9722,14.9721,14.9726,14.9723,14.9721,14.9716,14.9708,14.9701,14.9686,14.967,14.9653,14.9634,14.9618,14.9599,14.9583,14.9564,14.9557,14.9553,14.9553,14.9554,14.9557,14.9555,14.9562,14.9563,14.9559,14.9563,14.9567,14.9568,14.9577,14.9586,14.9589,14.959,14.959,14.9587,14.959,14.9595,14.9593,14.9594,14.9597,14.9599,14.9598,14.9604,14.9605,14.9606,14.961,14.961,14.9611,14.9607,14.9608,14.9602,14.9597,14.9588,14.9585,14.9581,14.958,14.9575,14.9573,14.9566,14.9565,14.9565,14.9564,14.9564,14.957,14.9574,14.9579,14.9585,14.9597,14.9603,14.9609,14.9608,14.9612,14.9617,14.9621,14.9626,14.9631,14.9631,14.9628,14.9626,14.9635,14.9636,14.9632,14.9633,14.963,14.9629,14.9628,14.9625,14.962,14.961,14.9603,14.9592,14.9589,14.9579,14.9574,14.9562,14.9558,14.9544,14.9534,14.9535,14.9527,14.952,14.9519,14.9525,14.9528,14.953,14.9546,14.956,14.9575,14.9598,14.962,14.9633,14.965,14.9668,14.9677,14.9694,14.971,14.9722,14.9732,14.9751,14.9761,14.9771,14.9783,14.9799,14.9806,14.982,14.9836,14.9844,14.985,14.986,14.987,14.9879,14.9888,14.9898,14.9903,14.9902,14.9913,14.9921,14.9925,14.9934,14.9938,14.9947,14.9951,14.9954,14.9957,14.9963,14.9969,14.9975,14.9982,15.0009,15.001,15.0018,15.0027,15.0039,15.0048,15.0057,15.0067,15.0086,15.0094,15.0107,15.0119,15.0127,15.0137,15.0144,15.0154,15.0165,15.0171,15.018,15.0189,15.0192,15.0199,15.0207,15.021,15.021,15.0212,15.0219,15.0221,15.0218,15.0222,15.0226,15.0228,15.023,15.0234,15.0236,15.0239,15.0242,15.0244,15.0246,15.0246,15.0245,15.0244,15.0243,15.0243,15.0243,15.0238,15.0235,15.0222,15.022,15.0213,15.0194,15.0176,15.0162,15.0153,15.0148,15.0143,15.0136,15.0125,15.0109,15.0084,15.0061,15.0036,15.0026,15.0008,14.9993,14.9984,14.9982,14.9978,14.9976,14.997,14.9967,14.9955,14.9945,14.9918,14.9897,14.9883,14.9862,14.9842,14.9824,14.9811,14.9796,14.9777,14.9757,14.9737,14.9724,14.9708,14.9704,14.9688,14.9673,14.9657,14.9642,14.9626,14.9618,14.96,14.9586,14.9574,14.9566,14.9553,14.9532,14.9506,14.9497,14.948,14.9451,14.9402,14.9362,14.9335,14.9292,14.9273,14.9253,14.9235,14.9219,14.9199,14.918,14.9166,14.9155,14.9143,14.9128,14.9108,14.9093,14.9082,14.9068,14.9053,14.905,14.9022,14.8988,14.8963,14.895,14.8938,14.8922,14.8905,14.8887,14.8874,14.8865,14.8858,14.8844,14.8822,14.8805,14.8796,14.8784,14.8772,14.8759,14.8742,14.8729,14.8718,14.8701,14.8688,14.8676,14.8664,14.8651,14.8632,14.8612,14.8591,14.8553,14.8529,14.8512,14.8491,14.8479,14.8468,14.8449,14.8411,14.8398,14.838,14.8363,14.8344,14.8329,14.8319,14.83,14.8281,14.8264,14.8247,14.823,14.822,14.8206,14.8187,14.8173,14.8148,14.812,14.8073,14.8057,14.8044,14.8027,14.8005,14.7991,14.7978,14.7961,14.7946,14.7936,14.7925,14.7793,14.7777,14.7768,14.7749,14.773,14.7699,14.7682,14.7661,14.764,14.7616,14.758,14.7558,14.7543,14.7519,14.749,14.747,14.7452,14.7429,14.7408,14.7395,14.7363,14.7339,14.7312,14.7282,14.7258,14.7228,14.7207,14.7193,14.7175,14.7149,14.7125,14.7091,14.7061,14.704,14.7019,14.6997,14.6964,14.6945,14.6921,14.6898,14.6865,14.6849,14.6831,14.6802,14.6777,14.6744,14.6737,14.6717,14.67,14.6682,14.6669,14.6646,14.6619,14.659,14.6564,14.6536,14.6513,14.6487,14.646,14.6441,14.6422,14.6382,14.6371,14.6363,14.6352,14.6341,14.6332,14.6322,14.6314,14.6302,14.6281,14.6273,14.6258,14.6252,14.6231,14.6217,14.6201,14.6195,14.6182,14.6176,14.6168,14.6161,14.615,14.6142,14.613,14.6122,14.6117,14.6099,14.6093,14.6089,14.6079,14.6072,14.6072,14.6071,14.6067,14.606,14.6057,14.6053,14.605,14.6044,14.6044,14.604,14.6038,14.6036,14.6038,14.604,14.6042,14.6042,14.6038,14.6042,14.6043,14.6046,14.6045,14.6046,14.6046,14.6047,14.6051,14.6054,14.6053,14.6056,14.6059,14.6059,14.6057,14.6057,14.6058,14.6059,14.606,14.6061,14.6061,14.606,14.606,14.606,14.6057,14.606,14.6052,14.6043,14.6037,14.6024,14.6006,14.5987,14.5964,14.5948,14.5934,14.5917,14.5904,14.5892,14.588,14.5864,14.5845,14.5831,14.5811,14.5795,14.5784,14.5766,14.5753,14.5741,14.5728,14.5715,14.5704,14.5691,14.5677,14.5655,14.5638,14.5624,14.5615,14.5599,14.5585,14.5575,14.5543,14.5521,14.5503,14.5484,14.5466,14.5455,14.5432,14.541,14.5386,14.5373,14.535,14.5328,14.5316,14.5303,14.529,14.5276,14.5252,14.5233,14.5224,14.5203,14.5197,14.5184,14.517,14.5151,14.5132,14.5111,14.5096,14.508,14.5061,14.5044,14.503,14.5015,14.5001,14.4992,14.4963,14.4955,14.494,14.4921,14.4909,14.4891,14.4871,14.4855,14.4834,14.4814,14.4797,14.4776,14.4765,14.475,14.4735,14.4712,14.4698,14.4677,14.4663,14.4646,14.4624,14.4611,14.4601,14.4586,14.4561,14.4554,14.4545,14.4536,14.4529,14.4525,14.4523,14.4515,14.4506,14.4501,14.449,14.448,14.4472,14.4451,14.4437,14.4421,14.4406,14.4391,14.4383,14.4366,14.435,14.4335,14.4314,14.43,14.4287,14.428,14.4275,14.4271,14.4267,14.4266,14.4263,14.4262,14.4266,14.4259,14.4258,14.4261,14.426,14.4257,14.4257,14.4259,14.4256,14.4259,14.4251,14.4239,14.423,14.4222,14.4199,14.418,14.416,14.4142,14.4121,14.4098,14.4081,14.4067,14.4046,14.4028,14.4004,14.398,14.395,14.3918,14.3891,14.3859,14.3835,14.3805,14.3778,14.3752,14.373,14.3707,14.368,14.3647,14.363,14.3608,14.3584,14.3551,14.3523,14.3499,14.3473,14.3449,14.3419,14.339,14.3363,14.3332,14.3306,14.3287,14.3262,14.3237,14.321,14.3186,14.3166,14.314,14.3117,14.3097,14.3069,14.3052,14.3024,14.3004,14.2984,14.2963,14.294,14.2918,14.2894,14.2876,14.2861,14.284,14.2823,14.2807,14.2791,14.2772,14.2752,14.2731,14.2713,14.2699,14.2684,14.2669,14.2658,14.2638,14.2623,14.2599,14.2574,14.2552,14.2539,14.2524,14.2504,14.2478,14.2454,14.2432,14.2407,14.2376,14.2348,14.2317,14.2294,14.2273,14.2252,14.2231,14.2222,14.2209,14.2201,14.2189,14.2173,14.2156,14.2144,14.2126,14.2112,14.2111,14.2113,14.2104,14.2091,14.208,14.2064,14.2048,14.2035,14.2011,14.1988,14.1971,14.1952,14.1932,14.1915,14.1895,14.1889,14.1879,14.1864,14.1848,14.1829,14.1823,14.1827,14.1827,14.1825,14.182,14.1806,14.1789,14.1766,14.1754,14.1739,14.1729,14.1715,14.171,14.1708,14.1705,14.1698,14.1691,14.1674,14.1655,14.1635,14.1625,14.1609,14.1594,14.1576,14.1562,14.1546,14.1532,14.1518,14.1502,14.1483,14.1468,14.1458,14.1448,14.1432,14.1415,14.141,14.1403,14.1401,14.1407,14.1412,14.142,14.1424,14.142,14.1412,14.1403,14.1395,14.1385,14.138,14.1379,14.138,14.1387,14.1394,14.1403,14.1414,14.1411,14.1403,14.14,14.1397,14.124,14.1238,14.1221,14.1137,14.1104,14.1096,14.1079,14.1055,14.1037,14.0971,14.0913,14.0904,14.0871,14.0863,14.0746,14.0737,14.0721,14.0696,14.062,14.0587,14.0579,14.0571,14.0546,14.0537,14.0504,14.0496,14.0454,14.0446,14.0379,14.0371,14.0321,14.0312,14.0229,14.0221,14.0204,14.0188,14.0096,14.0071,14.0054,14.0046,13.9996,13.9913,13.9904,13.9863,13.9854,13.9779,13.9771,13.9738,13.9729,13.9671,13.9663,13.9579,13.9571,13.9513,13.9504,13.9496,13.9446,13.9437,13.9421,13.9346,13.9338,13.9213,13.9204,13.918,13.8971,13.8963,13.8896,13.8887,13.8771,13.8763,13.868,13.8671,13.8554,13.8546,13.8496,13.8471,13.8346,13.8333,13.8334,13.8354,13.8471,13.8479,13.8562,13.857,13.8588,13.8596,13.8621,13.8662,13.8688,13.8696,13.8729,13.8738,13.8887,13.8896,13.8921,13.8954,13.8988,13.9004,13.9038,13.9063,13.9071,13.9079,13.9129,13.9138,13.9221,13.9229,13.9354,13.9362,13.9446,13.9454,13.9487,13.9496,13.9579,13.9629,13.9646,13.9662,13.9696,13.9721,13.9762,13.9792,13.98,13.9825,13.9825,13.9817,13.9817,13.9891,13.99,13.9946,13.9963,13.9988,13.9996,14.0004,14.0174,14.0188,14.0204,14.0233,14.0233,14.0242,14.0263,14.0287,14.0304,14.0325,14.0325,14.0313,14.0309,14.0342,14.0359,14.0383,14.0384,14.04,14.0442,14.0442,14.0429,14.0379,14.0342,14.0317,14.0288,14.0267,14.0241,14.0209,14.02,14.0171,14.0154,14.0129,14.0113,14.0092,14.0092,14.0108,14.01,14.0109,14.0117,14.0142,14.0167,14.0167,14.0183,14.0171,14.0154,14.0125,14.0125,14.0062,14.0025,14.0025,14.0041,14.0067,14.0067,14.0042,14.0012,13.9996,13.9987,13.9887,13.9854,13.9829,13.9771,13.9767,13.9796,13.9809,13.978,13.9745,13.9717,13.9717,13.9688,13.9679,13.9654,13.9546,13.9521,13.9513,13.9471,13.9437,13.9404,13.9304,13.9279,13.9271,13.9254,13.9212,13.9179,13.9121,13.9108,13.9109,13.91,13.91,13.9092,13.9092,13.9083,13.9069,13.9059,13.9037,13.9012,13.9004,13.8979,13.8954,13.8871,13.8854,13.8846,13.8829,13.8813,13.8763,13.8734,13.8717,13.8733,13.8692,13.8692,13.8671,13.8663,13.8659,13.8675,13.8671,13.8662,13.8621,13.8571,13.8554,13.8538,13.8525,13.8546,13.8554,13.8604,13.8617,13.8592,13.8584,13.8554,13.8546,13.8504,13.8496,13.8471,13.8454,13.8412,13.8404,13.8379,13.8379,13.8363,13.8312,13.8262,13.8254,13.8238,13.8188,13.8155,13.8088,13.8041,13.8042,13.8025,13.8025,13.7988,13.7904,13.7854,13.7771,13.7742,13.7733,13.7733,13.7742,13.7742,13.775,13.7779,13.7804,13.7846,13.7858,13.7833,13.7825,13.7803,13.7784,13.7787,13.7812,13.7821,13.7838,13.785,13.785,13.7875,13.7875,13.7887,13.7971,13.7979,13.8029,13.8071,13.8075,13.805,13.805,13.8084,13.8079,13.8062,13.8029,13.8021,13.7991,13.7992,13.8029,13.8038,13.81,13.8108,13.8109,13.8125,13.8134,13.8146,13.8196,13.8229,13.8254,13.8288,13.8354,13.8363,13.8379,13.8384,13.8354,13.8321,13.8312,13.8271,13.8242,13.8242,13.8254,13.8284,13.8271,13.8242,13.8233,13.8234,13.8242,13.8242,13.8267,13.8275,13.83,13.8333,13.8325,13.8334,13.8333,13.835,13.8358,13.8375,13.8359,13.8358,13.8367,13.8367,13.8383,13.8392,13.8434,13.8442,13.8433,13.8475,13.8483,13.8534,13.8537,13.8571,13.8579,13.8612,13.8662,13.8679,13.8688,13.8763,13.8788,13.88,13.8829,13.8846,13.8871,13.892,13.8854,13.8846,13.8821,13.8771,13.8696,13.8662,13.8625,13.8625,13.8617,13.8625,13.8654,13.8671,13.8688,13.87,13.87,13.8692,13.8692,13.8679,13.8671,13.8663,13.8642,13.8642,13.865,13.865,13.8658,13.8658,13.8667,13.8683,13.8708,13.8708,13.8746,13.8788,13.845,13.8467,13.8467,13.8459,13.8459,13.845,13.845,13.8433,13.8425,13.8425,13.8417,13.8416,13.84,13.84,13.8388,13.8354,13.8334,13.8325,13.8325,13.8325,13.83,13.8296,13.8279,13.8262,13.8204,13.82,13.8234,13.8237,13.8271,13.8275,13.825,13.825,13.8287,13.83,13.83,13.8309,13.8309,13.8309,13.8317,13.8317,13.8308,13.8309,13.83,13.83,13.8292,13.8292,13.8284,13.8283,13.8267,13.825,13.8221,13.8171,13.8167,13.8146,13.8142,13.8083,13.8067,13.8067,13.7988,13.7975,13.795,13.7929,13.7913,13.7879,13.7875,13.7883,13.7883,13.7875,13.7821,13.7754,13.7687,13.7604,13.76,13.7617,13.7604,13.7555,13.7438,13.732,13.7308,13.7309,13.7325,13.7325,13.7333,13.7334,13.7342,13.7342,13.735,13.7312,13.7308,13.7292,13.73,13.7283,13.7283,13.7275,13.7267,13.7258,13.7258,13.7237,13.7146,13.6996,13.6912,13.6896,13.6846,13.6829,13.6795,13.6779,13.6646,13.6629,13.662,13.6579,13.6521,13.6521,13.6546,13.6554,13.6571,13.6579,13.6613,13.6621,13.664,13.6679,13.6696,13.6696,13.6637,13.6608,13.6609,13.6642,13.6637,13.6562,13.6558,13.6592,13.6592,13.6571,13.6546,13.65,13.65,13.6484,13.6475,13.645,13.6434,13.6433,13.6442,13.6442,13.6558,13.6575,13.6575,13.6596,13.6621,13.6534,13.6508,13.6492,13.6479,13.6437,13.6396,13.6379,13.6358,13.6359,13.6342,13.6342,13.6308,13.6292,13.6292,13.63,13.63,13.6263,13.6229,13.6225,13.625,13.625,13.6242,13.6241,13.6267,13.6267,13.6254,13.6238,13.6204,13.6121,13.6112,13.6071,13.6,13.6,13.5986,13.5974,13.5969,13.597,13.5996,13.6013,13.6031,13.6041,13.6045,13.6016,13.6003,13.6,13.5996,13.5982,13.5954,13.5935,13.5913,13.5902,13.5896,13.5896,13.5899,13.5899,13.59,13.5902,13.5902,13.5902,13.5904,13.5904,13.5904,13.5904,13.5907,13.5907,13.591,13.591,13.591,13.5913,13.5913,13.5915,13.5915,13.5915,13.5918,13.5918,13.5918,13.5921,13.5915,13.5915,13.5915,13.5918,13.5918,13.5918,13.5918,13.5918,13.5921,13.5921,13.5921,13.5918,13.5918,13.5921,13.5921,13.5921,13.5924,13.5915,13.5918,13.5918,13.5918,13.5913,13.5913,13.5913,13.5915,13.5907,13.591,13.591,13.591,13.5904,13.5904,13.5904,13.5907,13.5899,13.5902,13.5902,13.5902,13.5902,13.5896,13.5896,13.5896,13.5896,13.5895,13.589,13.5918,13.5979,13.6,13.6049,13.611,13.6171,13.6229,13.6285,13.6346,13.6402,13.6404,13.6457,13.6507,13.6557,13.6613,13.666,13.671,13.676,13.6807,13.6852,13.6899,13.6949,13.699,13.704,13.7077,13.7118,13.7154,13.719,13.7227,13.7263,13.7299,13.7327,13.7363,13.739,13.7429,13.7457,13.7458,13.7488,13.7515,13.7546,13.7577,13.7604,13.7635,13.7671,13.7699,13.7721,13.7752,13.7779,13.781,13.7832,13.7854,13.7877,13.7902,13.7915,13.7932,13.7949,13.7952,13.7949,13.7935,13.7918,13.789,13.7879,13.7854,13.7821,13.7788,13.7763,13.7729,13.7688,13.7657,13.7615,13.7577,13.7535,13.7496,13.7454,13.7443,13.7538,13.7635,13.7738,13.7835,13.794,13.7946,13.8002,13.8057,13.8118,13.8174,13.8238,13.826,13.819,13.8121,13.8046,13.7986,13.7979,13.7946,13.7954,13.7979,13.7993,13.801,13.8027,13.8043,13.8052,13.8054,13.8049,13.8046,13.804,13.8038,13.8027,13.8013,13.8002,13.7985,13.7965,13.7946,13.7927,13.7902,13.7877,13.7857,13.7838,13.7804,13.7782,13.7779,13.7754,13.7721,13.7693,13.7673,13.7829,13.7849,13.7881,13.7913,13.7935,13.7962,13.7979,13.8004,13.8029,13.8046,13.8064,13.8078,13.8094,13.8095,13.8098,13.8102,13.8107,13.8109,13.8109,13.8112,13.8115,13.8117,13.812,13.8126,13.813,13.814,13.8145,13.8152,13.8162,13.8174,13.8208,13.8232,13.827,13.8304,13.8331,13.8366,13.8398,13.8428,13.8436,13.8488,13.8549,13.8602,13.862,13.8651,13.8714,13.8791,13.8836,13.8858,13.8887,13.8905,13.8913,13.8939,13.8974,13.9,13.9024,13.9051,13.9077,13.9093,13.9111,13.9137,13.9166,13.9188,13.9211,13.9242,13.9263,13.9286,13.9312,13.9335,13.9366,13.9374,13.9382,13.9412,13.9436,13.9473,13.95,13.9522,13.9564,13.9599,13.963,13.9663,13.9698,13.9724,13.9761,13.9787,13.9811,13.9834,13.9855,13.9863,13.9887,13.991,13.994,13.9957,13.9984,14.0001,14.0007,14.0014,14.0022,14.0041,14.0066,14.0079,14.0097,14.0111,14.0127,14.0162,14.0192,14.0221,14.0249,14.0269,14.0307,14.0344,14.039,14.0408,14.0436,14.0462,14.0494,14.0517,14.0537,14.0572,14.0591,14.0617,14.0629,14.0676,14.0695,14.0723,14.0743,14.0773,14.0794,14.0824,14.0837,14.0854,14.0873,14.0892,14.0909,14.0948,14.0959,14.0976,14.0995,14.1021,14.1036,14.1058,14.1076,14.1102,14.1118,14.1143,14.1172,14.1181,14.1214,14.1234,14.1251,14.1278,14.13,14.1325,14.1351,14.1367,14.1378,14.1385,14.1396,14.1416,14.144,14.1465,14.1493,14.1527,14.1551,14.1581,14.161,14.162,14.1638,14.1663,14.1692,14.1721,14.1748,14.1768,14.1784,14.1812,14.1834,14.1859,14.1884,14.1909,14.1928,14.1944,14.1958,14.1971,14.1989,14.2012,14.2026,14.2057,14.2074,14.2091,14.213,14.2154,14.2173,14.2192,14.2213,14.2243,14.2253,14.2282,14.231,14.2348,14.2381,14.2402,14.2431,14.2457,14.2484,14.2514,14.2537,14.2552,14.2581,14.2592,14.2601,14.2623,14.2649,14.2676,14.2699,14.2721,14.274,14.2769,14.2812,14.2832,14.2868,14.2895,14.2906,14.2919,14.2958,14.2974,14.2986,14.301,14.3278,14.3546,14.3813,14.4081,14.4349,14.4617,14.4885,14.5153,14.5282,14.542,14.5623,14.5649,14.5672,14.5703,14.572,14.5733,14.5751,14.5754,14.5755,14.5758,14.5758,14.5761,14.5761,14.5763,14.5793,14.5804,14.5807,14.5817,14.5818,14.5827,14.5827,14.583,14.5831,14.5827,14.5831,14.583,14.5829,14.5829,14.5877,14.5926,14.598,14.5986,14.5992,14.6005,14.6005,14.6057,14.6057,14.612,14.614,14.6151,14.6162,14.6179,14.6181,14.6183,14.6206,14.6233,14.6413,14.6551,14.6551,14.6785,14.7019,14.6851,14.6684,14.6516,14.6457,14.6348,14.6097,14.615,14.6202,14.6212,14.6218,14.6255,14.6307,14.6359,14.6412,14.6464,14.6516,14.6569,14.6621,14.6786,14.695,14.7115,14.728,14.7445,14.7476,14.7532,14.7582,14.763,14.7667,14.7726,14.7766,14.781,14.7842,14.7886,14.7929,14.7966,14.8003,14.8051,14.8078,14.8101,14.8127,14.8155,14.8177,14.8208,14.8222,14.823,14.8241,14.8249,14.8254,14.8259,14.8281,14.8293,14.8308,14.8322,14.8334,14.8346,14.8355,14.8358,14.8363,14.8371,14.8377,14.8384,14.839,14.8395,14.8401,14.8401,14.8398,14.8395,14.8395,14.8393,14.8389,14.8386,14.8377,14.837,14.8366,14.8358,14.8349,14.8339,14.8333,14.8322,14.8314,14.8304,14.8287,14.8272,14.8259,14.8253,14.8245,14.8237,14.8237,14.825,14.8273,14.8292,14.8316,14.8335,14.8352,14.8378,14.8395,14.8417,14.843,14.8441,14.8448,14.8469,14.8482,14.8505,14.8531,14.8551,14.856,14.8577,14.8591,14.8597,14.8611,14.8621,14.8632,14.8634,14.8641,14.8645,14.865,14.8656,14.8659,14.8669,14.8669,14.8671,14.8673,14.8674,14.8674,14.8676,14.8686,14.8688,14.8684,14.8685,14.8687,14.8688,14.869,14.869,14.8692,14.8691,14.8697,14.8695,14.8673,14.866,14.8649,14.8635,14.8626,14.8615,14.8597,14.8572,14.8546,14.8578,14.8622,14.8721,14.8815,14.8859,14.8832,14.8794,14.8794,14.8884,14.8965,14.903,14.9102,14.9151,14.9206,14.9305,14.9421,14.9455,14.9459,14.9416,14.9314,14.9208,14.9132,14.9086,14.9062,14.9023,14.9,14.8987,14.8982,14.8977,14.8984,14.9013,14.9042,14.9085,14.9119,14.915,14.9179,14.9222,14.9249,14.9282,14.9306,14.933,14.9358,14.9391,14.9412,14.9438,14.9459,14.9479,14.9497,14.9521,14.9547,14.9568,14.9588,14.9601,14.9611,14.9619,14.9635,14.9661,14.968,14.9688,14.9801]}]],[[{"lng":[-16.7967,-16.7975,-16.798,-16.798,-16.7959,-16.7951,-16.7946,-16.7946,-16.7967],"lat":[12.4834,12.4834,12.4837,12.4846,12.4867,12.4867,12.4862,12.4854,12.4834]},{"lng":[-16.7959,-16.7947,-16.7946,-16.795,-16.7967,-16.798,-16.798,-16.7976,-16.7959],"lat":[12.4884,12.4896,12.4904,12.4908,12.4908,12.4896,12.4887,12.4884,12.4884]},{"lng":[-16.7922,-16.7921,-16.7925,-16.7942,-16.7955,-16.7942,-16.7926,-16.7922],"lat":[12.492,12.4938,12.4942,12.4942,12.4929,12.4917,12.4917,12.492]},{"lng":[-16.6746,-16.6759,-16.6825,-16.6838,-16.6825,-16.6759,-16.6746],"lat":[12.5038,12.505,12.505,12.5038,12.5025,12.5025,12.5038]},{"lng":[-16.515,-16.5138,-16.5151,-16.5192,-16.5196,-16.5176,-16.515],"lat":[12.5442,12.5471,12.5484,12.5483,12.5462,12.5441,12.5442]},{"lng":[-16.6988,-16.7001,-16.7017,-16.7059,-16.7075,-16.7092,-16.7109,-16.7122,-16.7121,-16.7193,-16.7217,-16.7238,-16.7238,-16.7209,-16.7209,-16.7284,-16.7317,-16.7358,-16.7371,-16.7372,-16.738,-16.738,-16.7388,-16.7388,-16.7396,-16.7396,-16.7371,-16.7355,-16.7346,-16.7267,-16.7242,-16.7208,-16.7168,-16.7151,-16.7142,-16.6951,-16.6942,-16.6909,-16.6901,-16.6817,-16.6788,-16.6788,-16.6796,-16.6796,-16.6805,-16.6805,-16.6796,-16.6796,-16.6788,-16.6788,-16.6817,-16.6859,-16.69,-16.693,-16.693,-16.6905,-16.6905,-16.6934,-16.6938,-16.693,-16.693,-16.6938,-16.6938,-16.6938,-16.6975,-16.6984,-16.7,-16.7034,-16.7059,-16.7071,-16.708,-16.7084,-16.7097,-16.7096,-16.7088,-16.7075,-16.7017,-16.7009,-16.6996,-16.6988,-16.6988],"lat":[12.5579,12.5592,12.5592,12.5567,12.5567,12.5592,12.5592,12.558,12.5571,12.5487,12.5458,12.5462,12.5479,12.5492,12.55,12.55,12.5525,12.5525,12.5512,12.5487,12.5479,12.5429,12.5421,12.5346,12.5337,12.5321,12.5279,12.5263,12.5229,12.515,12.5142,12.5108,12.5092,12.5092,12.5084,12.5084,12.5092,12.5092,12.51,12.51,12.5137,12.5163,12.517,12.5196,12.5204,12.5279,12.5288,12.5312,12.5321,12.5388,12.5408,12.5409,12.5434,12.5429,12.5404,12.5379,12.5371,12.535,12.5363,12.5371,12.5388,12.5396,12.5434,12.5454,12.5458,12.545,12.545,12.5475,12.5475,12.5463,12.5446,12.5408,12.5421,12.5463,12.5479,12.5492,12.5517,12.5534,12.5537,12.5555,12.5579]},{"lng":[-16.3234,-16.3226,-16.3221,-16.3221,-16.3226,-16.3242,-16.3255,-16.3242,-16.3234],"lat":[12.5892,12.5892,12.5896,12.5904,12.5908,12.5909,12.5896,12.5883,12.5892]},{"lng":[-16.3242,-16.3238,-16.3238,-16.323,-16.323,-16.3234,-16.3251,-16.3263,-16.3263,-16.3259,-16.3242],"lat":[12.5917,12.5921,12.5929,12.5937,12.5946,12.595,12.595,12.5938,12.592,12.5917,12.5917]},{"lng":[-16.3276,-16.3263,-16.3263,-16.3276,-16.33,-16.3313,-16.3313,-16.3292,-16.3276],"lat":[12.595,12.5962,12.5988,12.6,12.6,12.5987,12.5971,12.595,12.595]},{"lng":[-16.6967,-16.7009,-16.7021,-16.7021,-16.7021,-16.7009,-16.6993,-16.6984,-16.6925,-16.6896,-16.6896,-16.6918,-16.6934,-16.6967],"lat":[12.6042,12.6058,12.6046,12.5982,12.5954,12.5942,12.5942,12.595,12.595,12.5971,12.5979,12.6,12.6,12.6042]},{"lng":[-16.4663,-16.4663,-16.4684,-16.4705,-16.4705,-16.4684,-16.4663],"lat":[12.6538,12.6562,12.6592,12.6588,12.6554,12.6525,12.6538]},{"lng":[-16.7894,-16.7871,-16.7855,-16.7846,-16.7838,-16.7842,-16.7858,-16.788,-16.7888,-16.7892,-16.7896,-16.7896,-16.7894],"lat":[12.7504,12.7504,12.7538,12.7588,12.7596,12.7609,12.7608,12.7587,12.7554,12.7551,12.7546,12.7504,12.7504]},{"lng":[-12.3659,-12.3711,-12.3767,-12.3811,-12.3861,-12.3906,-12.395,-12.3986,-12.4023,-12.4059,-12.4089,-12.4125,-12.4148,-12.4184,-12.422,-12.4256,-12.43,-12.4336,-12.4375,-12.4411,-12.4448,-12.4478,-12.4481,-12.45,-12.4517,-12.4525,-12.4528,-12.4517,-12.4503,-12.4498,-12.4514,-12.4542,-12.4567,-12.4598,-12.462,-12.4642,-12.4664,-12.4681,-12.4698,-12.4714,-12.4727,-12.4728,-12.4753,-12.4781,-12.482,-12.4856,-12.4906,-12.4956,-12.4992,-12.5014,-12.5064,-12.5123,-12.5173,-12.5223,-12.5275,-12.5317,-12.5367,-12.542,-12.5467,-12.5511,-12.5564,-12.5606,-12.565,-12.5686,-12.5731,-12.5767,-12.5811,-12.5848,-12.5886,-12.5923,-12.5959,-12.5989,-12.6017,-12.6053,-12.6086,-12.6131,-12.6173,-12.6214,-12.6259,-12.6311,-12.6367,-12.6425,-12.6484,-12.6539,-12.6611,-12.6675,-12.6725,-12.6759,-12.6786,-12.682,-12.6853,-12.6856,-12.6861,-12.6641,-12.6444,-12.6247,-12.6237,-12.6235,-12.6228,-12.6219,-12.6211,-12.6213,-12.6211,-12.6212,-12.6209,-12.6213,-12.6205,-12.6192,-12.6177,-12.6163,-12.6155,-12.6146,-12.6141,-12.6141,-12.6149,-12.615,-12.6152,-12.6155,-12.6156,-12.6144,-12.614,-12.6126,-12.6117,-12.6102,-12.611,-12.612,-12.6126,-12.6129,-12.6128,-12.6118,-12.6107,-12.6098,-12.6087,-12.6074,-12.6066,-12.6058,-12.6051,-12.6053,-12.6049,-12.6042,-12.6042,-12.6038,-12.6032,-12.6031,-12.6026,-12.6022,-12.6022,-12.602,-12.6019,-12.6014,-12.6006,-12.6008,-12.6013,-12.6007,-12.6011,-12.6011,-12.6005,-12.6005,-12.6001,-12.6001,-12.6006,-12.6006,-12.6,-12.5992,-12.5987,-12.597,-12.5967,-12.596,-12.5953,-12.5948,-12.594,-12.5937,-12.5928,-12.592,-12.5914,-12.591,-12.5904,-12.5892,-12.5885,-12.5878,-12.587,-12.5866,-12.5859,-12.5942,-12.5954,-12.6049,-12.6144,-12.6147,-12.6239,-12.6334,-12.6429,-12.651,-12.6524,-12.6693,-12.6861,-12.703,-12.7198,-12.7449,-12.7701,-12.7952,-12.8211,-12.847,-12.8729,-12.8988,-12.9247,-12.9506,-12.9537,-12.9745,-12.9989,-13.0174,-13.0363,-13.0659,-13.0669,-13.0956,-13.1252,-13.1548,-13.1845,-13.2106,-13.2141,-13.2438,-13.259,-13.2742,-13.2894,-13.3046,-13.3198,-13.335,-13.3502,-13.3654,-13.3806,-13.3958,-13.4109,-13.4261,-13.4413,-13.4565,-13.4717,-13.4869,-13.5021,-13.5172,-13.5323,-13.5471,-13.577,-13.6069,-13.6368,-13.6667,-13.6966,-13.7162,-13.7265,-13.7564,-13.7834,-13.7864,-13.7868,-13.7874,-13.7886,-13.8163,-13.8194,-13.8462,-13.8717,-13.8972,-13.9227,-13.9482,-13.9645,-13.9737,-13.9948,-13.9989,-14.0016,-14.0253,-14.034,-14.0531,-14.0565,-14.1062,-14.1545,-14.1564,-14.1565,-14.1642,-14.1792,-14.1792,-14.1825,-14.1851,-14.1875,-14.2384,-14.2385,-14.2414,-14.2414,-14.2775,-14.2846,-14.2904,-14.2941,-14.3236,-14.3286,-14.33,-14.3301,-14.3301,-14.3302,-14.3513,-14.3515,-14.3681,-14.3683,-14.3799,-14.3998,-14.408,-14.4166,-14.428,-14.4356,-14.4472,-14.456,-14.4666,-14.4693,-14.4703,-14.4707,-14.4713,-14.4739,-14.4755,-14.4776,-14.4799,-14.4814,-14.4824,-14.4834,-14.4837,-14.4846,-14.4858,-14.4879,-14.4882,-14.4895,-14.4904,-14.4913,-14.4937,-14.4939,-14.494,-14.494,-14.4976,-14.4985,-14.5025,-14.513,-14.519,-14.5247,-14.5316,-14.5413,-14.5484,-14.5566,-14.5604,-14.5658,-14.5674,-14.5717,-14.5733,-14.5738,-14.5821,-14.5894,-14.5993,-14.6053,-14.6076,-14.6271,-14.6329,-14.6514,-14.6577,-14.6571,-14.6528,-14.6485,-14.6428,-14.6326,-14.6249,-14.6251,-14.6269,-14.6277,-14.6284,-14.6299,-14.6314,-14.6328,-14.6343,-14.6358,-14.6373,-14.6388,-14.6403,-14.6393,-14.6376,-14.6361,-14.6317,-14.6304,-14.6292,-14.6267,-14.6235,-14.6201,-14.6161,-14.6131,-14.611,-14.6085,-14.6055,-14.6029,-14.6003,-14.5975,-14.5961,-14.5961,-14.5954,-14.5954,-14.5953,-14.595,-14.5946,-14.5946,-14.5941,-14.5938,-14.5935,-14.5932,-14.5933,-14.5932,-14.5924,-14.592,-14.5919,-14.5918,-14.5912,-14.5915,-14.5909,-14.5908,-14.5903,-14.59,-14.5901,-14.5901,-14.5901,-14.5901,-14.5898,-14.5894,-14.5889,-14.5891,-14.5889,-14.5886,-14.5884,-14.5882,-14.5884,-14.5881,-14.5876,-14.5871,-14.5868,-14.5868,-14.5866,-14.5868,-14.5867,-14.586,-14.5855,-14.5854,-14.5848,-14.5846,-14.5845,-14.5847,-14.5846,-14.5845,-14.5848,-14.5852,-14.5859,-14.5863,-14.5864,-14.5867,-14.587,-14.5876,-14.588,-14.5884,-14.5889,-14.5894,-14.5893,-14.5899,-14.5907,-14.5907,-14.5912,-14.5914,-14.592,-14.592,-14.5926,-14.593,-14.5935,-14.5933,-14.5931,-14.5928,-14.5924,-14.5923,-14.5921,-14.5921,-14.5913,-14.5913,-14.5905,-14.5904,-14.5899,-14.5901,-14.5894,-14.5893,-14.5888,-14.5882,-14.5879,-14.5878,-14.5873,-14.5867,-14.5869,-14.5862,-14.5863,-14.5857,-14.5853,-14.5849,-14.5848,-14.5847,-14.5845,-14.5842,-14.5833,-14.5829,-14.5826,-14.5823,-14.5821,-14.5818,-14.5815,-14.5818,-14.5822,-14.5828,-14.5835,-14.5843,-14.5853,-14.5856,-14.5862,-14.5871,-14.5879,-14.5889,-14.59,-14.5907,-14.5922,-14.5934,-14.5945,-14.5959,-14.5972,-14.5985,-14.5992,-14.6005,-14.6013,-14.6033,-14.6043,-14.6053,-14.609,-14.6123,-14.6186,-14.6225,-14.6251,-14.6309,-14.6338,-14.6379,-14.6422,-14.6452,-14.6489,-14.6516,-14.655,-14.6587,-14.6631,-14.6673,-14.6709,-14.675,-14.679,-14.682,-14.687,-14.6919,-14.6943,-14.6973,-14.6999,-14.7009,-14.7016,-14.7042,-14.7069,-14.7092,-14.7095,-14.7105,-14.7121,-14.7135,-14.7145,-14.7161,-14.7175,-14.7186,-14.7201,-14.7239,-14.7303,-14.7362,-14.7441,-14.7497,-14.754,-14.7593,-14.7643,-14.7706,-14.7757,-14.7792,-14.7845,-14.7906,-14.7961,-14.7993,-14.8016,-14.8023,-14.8047,-14.8065,-14.8082,-14.8093,-14.8107,-14.8125,-14.8139,-14.8158,-14.8174,-14.82,-14.8213,-14.8234,-14.8407,-14.8392,-14.8353,-14.8317,-14.8281,-14.8245,-14.8209,-14.8184,-14.8153,-14.8131,-14.8106,-14.8086,-14.8061,-14.8039,-14.8023,-14.8006,-14.7986,-14.7973,-14.7964,-14.7945,-14.7928,-14.7925,-14.7917,-14.7914,-14.7911,-14.7731,-14.7548,-14.7367,-14.7186,-14.7178,-14.7136,-14.7103,-14.707,-14.7036,-14.7,-14.6959,-14.6925,-14.6884,-14.6842,-14.68,-14.6759,-14.6717,-14.6675,-14.6628,-14.6578,-14.6531,-14.6481,-14.6425,-14.637,-14.6311,-14.625,-14.6186,-14.6123,-14.605,-14.5981,-14.5925,-14.5867,-14.5809,-14.5753,-14.5695,-14.5645,-14.5586,-14.5536,-14.5492,-14.5445,-14.5392,-14.5348,-14.5306,-14.5267,-14.5225,-14.5184,-14.5153,-14.5109,-14.5073,-14.5042,-14.5006,-14.4992,-14.4967,-14.4945,-14.4917,-14.4892,-14.487,-14.4845,-14.4823,-14.4806,-14.4798,-14.4781,-14.477,-14.4761,-14.475,-14.4742,-14.4734,-14.4723,-14.472,-14.4661,-14.4606,-14.455,-14.4506,-14.4463,-14.4456,-14.4425,-14.44,-14.4373,-14.4342,-14.4306,-14.427,-14.4237,-14.4234,-14.4181,-14.4145,-14.41,-14.4056,-14.4014,-14.3975,-14.3934,-14.3892,-14.3848,-14.3803,-14.3759,-14.3714,-14.3667,-14.3623,-14.357,-14.3523,-14.347,-14.342,-14.337,-14.3311,-14.3267,-14.3217,-14.3161,-14.3106,-14.3048,-14.2992,-14.2936,-14.2886,-14.2839,-14.2789,-14.2725,-14.2664,-14.2606,-14.2542,-14.2473,-14.24,-14.2361,-14.2325,-14.23,-14.2273,-14.2239,-14.2206,-14.2205,-14.2203,-14.217,-14.2128,-14.2086,-14.2045,-14.1998,-14.1948,-14.1906,-14.1856,-14.18,-14.1753,-14.1689,-14.1634,-14.1575,-14.1514,-14.1456,-14.1386,-14.1323,-14.1259,-14.1231,-14.1225,-14.1198,-14.1164,-14.1175,-14.115,-14.1114,-14.1067,-14.1017,-14.097,-14.0923,-14.087,-14.0817,-14.0761,-14.0698,-14.0639,-14.0567,-14.0506,-14.0448,-14.0398,-14.0345,-14.0298,-14.0256,-14.0214,-14.0181,-14.0139,-14.0098,-14.0056,-14.0014,-13.9992,-13.9964,-13.9923,-13.9867,-13.9798,-13.9725,-13.9675,-13.9623,-13.9581,-13.9539,-13.95,-13.9436,-13.9373,-13.93,-13.9236,-13.9181,-13.9131,-13.9086,-13.9042,-13.8998,-13.8961,-13.8925,-13.8889,-13.8853,-13.8823,-13.88,-13.8775,-13.8761,-13.8745,-13.8728,-13.8698,-13.8686,-13.8678,-13.8589,-13.8501,-13.8499,-13.8492,-13.8403,-13.8334,-13.8267,-13.8211,-13.817,-13.8135,-13.8114,-13.8053,-13.8006,-13.7986,-13.7978,-13.7978,-13.7978,-13.8006,-13.8048,-13.8136,-13.8198,-13.8267,-13.8342,-13.8417,-13.8492,-13.8575,-13.8664,-13.8781,-13.8911,-13.8986,-13.9056,-13.9175,-13.93,-13.9425,-13.9539,-13.9631,-13.972,-13.9814,-13.992,-13.9992,-14.0042,-14.0173,-14.0295,-14.0406,-14.0481,-14.0584,-14.0686,-14.0823,-14.0948,-14.107,-14.115,-14.1198,-14.1253,-14.1331,-14.1417,-14.1509,-14.1603,-14.1714,-14.1817,-14.1858,-14.1934,-14.207,-14.2139,-14.2261,-14.2364,-14.2473,-14.257,-14.2673,-14.2789,-14.2911,-14.3036,-14.3106,-14.3242,-14.337,-14.3503,-14.3617,-14.372,-14.3831,-14.3928,-14.4,-14.4092,-14.4178,-14.427,-14.4336,-14.4406,-14.4459,-14.4517,-14.4611,-14.4748,-14.4881,-14.4992,-14.5003,-14.5111,-14.5203,-14.5306,-14.5381,-14.5456,-14.547,-14.5525,-14.5586,-14.5653,-14.5736,-14.5825,-14.592,-14.6025,-14.6148,-14.6264,-14.64,-14.647,-14.6598,-14.6723,-14.6834,-14.6942,-14.7045,-14.7139,-14.7228,-14.7311,-14.7381,-14.7436,-14.7489,-14.7545,-14.7592,-14.7642,-14.775,-14.7886,-14.8011,-14.8128,-14.8228,-14.8325,-14.8409,-14.847,-14.855,-14.8684,-14.8756,-14.8889,-14.9017,-14.9142,-14.9245,-14.9334,-14.9436,-14.9539,-14.9628,-14.9717,-14.98,-14.9889,-14.997,-14.9992,-14.9997,-15.0061,-15.0142,-15.0223,-15.0311,-15.0403,-15.0486,-15.0553,-15.0628,-15.0703,-15.0778,-15.0848,-15.0909,-15.0978,-15.1039,-15.1109,-15.1217,-15.1339,-15.145,-15.1548,-15.16,-15.165,-15.1703,-15.1742,-15.1767,-15.1839,-15.1923,-15.1998,-15.2017,-15.2025,-15.203,-15.203,-15.2034,-15.2034,-15.2048,-15.2039,-15.2067,-15.2081,-15.212,-15.217,-15.2239,-15.23,-15.2381,-15.2464,-15.2548,-15.2636,-15.2739,-15.2842,-15.295,-15.3075,-15.3198,-15.3334,-15.3377,-15.3459,-15.3595,-15.3717,-15.3834,-15.3936,-15.4056,-15.4173,-15.4286,-15.4398,-15.45,-15.4617,-15.4725,-15.4828,-15.4923,-15.4992,-15.5048,-15.5159,-15.5242,-15.5294,-15.5303,-15.5386,-15.5473,-15.5564,-15.567,-15.5803,-15.5939,-15.6009,-15.6084,-15.6206,-15.6323,-15.6448,-15.655,-15.6659,-15.6761,-15.6842,-15.6925,-15.7017,-15.712,-15.7239,-15.7364,-15.7503,-15.757,-15.7639,-15.7714,-15.7784,-15.7906,-15.8028,-15.8086,-15.8078,-15.8078,-15.8078,-15.8078,-15.8086,-15.8086,-15.8086,-15.8086,-15.8086,-15.8086,-15.8089,-15.8261,-15.8439,-15.8617,-15.8795,-15.8943,-15.8948,-15.8973,-15.9153,-15.9331,-15.9503,-15.9681,-15.9859,-15.9981,-15.9992,-16.0161,-16.0339,-16.0509,-16.0686,-16.0867,-16.1045,-16.1223,-16.14,-16.1573,-16.175,-16.1928,-16.2106,-16.2284,-16.2461,-16.2642,-16.2811,-16.2989,-16.3167,-16.3348,-16.3525,-16.3703,-16.3873,-16.4053,-16.4231,-16.4409,-16.4586,-16.4764,-16.4942,-16.4992,-16.5114,-16.5292,-16.547,-16.5648,-16.5828,-16.5836,-16.6006,-16.6175,-16.6353,-16.6534,-16.6711,-16.6889,-16.6998,-16.7048,-16.7075,-16.7075,-16.7142,-16.7203,-16.7286,-16.7328,-16.7306,-16.7324,-16.7315,-16.7285,-16.7285,-16.7291,-16.7302,-16.7309,-16.7315,-16.7326,-16.7342,-16.7359,-16.7373,-16.738,-16.7382,-16.7392,-16.74,-16.7412,-16.7446,-16.7459,-16.7461,-16.7471,-16.7469,-16.7469,-16.7465,-16.7454,-16.745,-16.7447,-16.7446,-16.743,-16.743,-16.7421,-16.7421,-16.7429,-16.743,-16.7421,-16.7421,-16.743,-16.7438,-16.7447,-16.7446,-16.7455,-16.7455,-16.7463,-16.7463,-16.7471,-16.7472,-16.748,-16.7488,-16.7488,-16.7504,-16.7513,-16.753,-16.753,-16.7555,-16.7547,-16.7554,-16.7555,-16.7563,-16.7563,-16.7571,-16.7571,-16.758,-16.758,-16.7588,-16.7588,-16.7596,-16.7605,-16.7613,-16.7613,-16.7621,-16.7621,-16.763,-16.763,-16.7663,-16.7671,-16.7705,-16.7713,-16.773,-16.7746,-16.7755,-16.7755,-16.7771,-16.7779,-16.7796,-16.7796,-16.7805,-16.7804,-16.7813,-16.7813,-16.7821,-16.783,-16.7837,-16.7838,-16.7846,-16.7846,-16.7855,-16.7855,-16.7863,-16.7863,-16.7871,-16.7871,-16.788,-16.7871,-16.788,-16.788,-16.7871,-16.7871,-16.7855,-16.7855,-16.7834,-16.7817,-16.7805,-16.7813,-16.7822,-16.7821,-16.7813,-16.7813,-16.7788,-16.7788,-16.778,-16.778,-16.7788,-16.7788,-16.778,-16.778,-16.7771,-16.7772,-16.7763,-16.7755,-16.7746,-16.7738,-16.7705,-16.7688,-16.768,-16.7671,-16.7671,-16.7663,-16.7663,-16.7655,-16.7663,-16.7663,-16.768,-16.768,-16.7688,-16.7721,-16.7738,-16.7738,-16.7709,-16.77,-16.7684,-16.7667,-16.7609,-16.7542,-16.7534,-16.75,-16.7484,-16.7425,-16.7384,-16.7317,-16.7313,-16.7358,-16.7417,-16.7434,-16.7442,-16.7467,-16.7484,-16.7501,-16.755,-16.7575,-16.7592,-16.7626,-16.7634,-16.7651,-16.7684,-16.7705,-16.7705,-16.7684,-16.7659,-16.7655,-16.7701,-16.773,-16.773,-16.7746,-16.7746,-16.7734,-16.7721,-16.7721,-16.7751,-16.7763,-16.7763,-16.7771,-16.7771,-16.7779,-16.7813,-16.783,-16.783,-16.7846,-16.7846,-16.7838,-16.7821,-16.7804,-16.7796,-16.7796,-16.7805,-16.7813,-16.7838,-16.7838,-16.7821,-16.7809,-16.7776,-16.7763,-16.7763,-16.7771,-16.7758,-16.7717,-16.7709,-16.7646,-16.7638,-16.7647,-16.7621,-16.7613,-16.7575,-16.7551,-16.7518,-16.7488,-16.7488,-16.7496,-16.7492,-16.7443,-16.7434,-16.7417,-16.7326,-16.7292,-16.7276,-16.7238,-16.7238,-16.7271,-16.7288,-16.7284,-16.7259,-16.7234,-16.7171,-16.7155,-16.7138,-16.7138,-16.7096,-16.7096,-16.7126,-16.7159,-16.7176,-16.7192,-16.7301,-16.7359,-16.7434,-16.7488,-16.7488,-16.7505,-16.7513,-16.7521,-16.7521,-16.7513,-16.7534,-16.7546,-16.7529,-16.753,-16.7534,-16.7542,-16.7567,-16.7605,-16.7621,-16.7638,-16.7646,-16.7663,-16.7663,-16.768,-16.768,-16.7696,-16.7705,-16.7705,-16.773,-16.7738,-16.7738,-16.7746,-16.7746,-16.7755,-16.7763,-16.7796,-16.7797,-16.775,-16.7717,-16.7701,-16.7676,-16.7642,-16.7625,-16.7601,-16.7584,-16.7571,-16.7621,-16.7634,-16.765,-16.7671,-16.7675,-16.7688,-16.7692,-16.7709,-16.773,-16.773,-16.7721,-16.773,-16.7713,-16.7696,-16.7705,-16.7696,-16.7696,-16.7688,-16.7688,-16.768,-16.7671,-16.7663,-16.7655,-16.7646,-16.7646,-16.7646,-16.7634,-16.7631,-16.7605,-16.7605,-16.7588,-16.7596,-16.7613,-16.7613,-16.7629,-16.7655,-16.7655,-16.7647,-16.7646,-16.7646,-16.763,-16.763,-16.7584,-16.7534,-16.7526,-16.745,-16.74,-16.738,-16.738,-16.7396,-16.7413,-16.7409,-16.7396,-16.7396,-16.7384,-16.7334,-16.7284,-16.7217,-16.7209,-16.7176,-16.7134,-16.7109,-16.7101,-16.705,-16.7042,-16.7017,-16.6992,-16.6967,-16.6951,-16.693,-16.693,-16.6942,-16.7009,-16.7017,-16.705,-16.7055,-16.7038,-16.7038,-16.7046,-16.7055,-16.7063,-16.7071,-16.7067,-16.7055,-16.7043,-16.7026,-16.7,-16.6992,-16.6925,-16.6859,-16.6784,-16.6709,-16.6683,-16.6671,-16.6663,-16.6646,-16.6646,-16.6675,-16.67,-16.6709,-16.6742,-16.6747,-16.6726,-16.6709,-16.6684,-16.6638,-16.6634,-16.6617,-16.6592,-16.6584,-16.6525,-16.6459,-16.6451,-16.64,-16.6392,-16.6329,-16.6329,-16.6359,-16.6368,-16.6371,-16.6359,-16.6342,-16.6317,-16.6305,-16.6305,-16.6313,-16.6313,-16.6338,-16.6338,-16.6346,-16.6355,-16.638,-16.638,-16.6388,-16.6388,-16.6396,-16.6396,-16.6394,-16.6371,-16.6351,-16.633,-16.6313,-16.6305,-16.6305,-16.6288,-16.6288,-16.6275,-16.6242,-16.6238,-16.6246,-16.6284,-16.6296,-16.6297,-16.6288,-16.628,-16.6226,-16.6209,-16.6167,-16.6159,-16.6117,-16.6109,-16.6083,-16.6075,-16.6042,-16.5988,-16.5988,-16.6013,-16.6009,-16.5975,-16.5942,-16.5917,-16.59,-16.588,-16.588,-16.5855,-16.5855,-16.5813,-16.5809,-16.5788,-16.5788,-16.5734,-16.5709,-16.5692,-16.565,-16.5588,-16.5588,-16.5576,-16.5572,-16.5563,-16.5563,-16.5555,-16.5546,-16.553,-16.5529,-16.5567,-16.5588,-16.5588,-16.5609,-16.5655,-16.5663,-16.5676,-16.5696,-16.5697,-16.5634,-16.5609,-16.5567,-16.5534,-16.5451,-16.5442,-16.5438,-16.5409,-16.5276,-16.5234,-16.5225,-16.52,-16.5184,-16.5146,-16.5146,-16.5075,-16.5063,-16.5063,-16.5013,-16.5013,-16.4976,-16.4934,-16.4921,-16.4922,-16.493,-16.4938,-16.4955,-16.4955,-16.4938,-16.4921,-16.4921,-16.493,-16.4917,-16.4876,-16.4847,-16.4855,-16.4871,-16.4871,-16.488,-16.488,-16.4901,-16.4922,-16.4905,-16.4905,-16.493,-16.493,-16.4905,-16.4904,-16.4892,-16.4876,-16.4846,-16.4846,-16.4863,-16.4859,-16.4775,-16.4746,-16.4746,-16.4759,-16.478,-16.4788,-16.4804,-16.4809,-16.4825,-16.485,-16.4855,-16.4846,-16.4834,-16.4817,-16.4792,-16.4788,-16.4788,-16.4805,-16.4805,-16.4826,-16.4851,-16.4868,-16.488,-16.488,-16.4817,-16.4776,-16.4768,-16.475,-16.473,-16.473,-16.4717,-16.47,-16.4675,-16.465,-16.4646,-16.4663,-16.468,-16.4671,-16.4651,-16.4642,-16.4621,-16.4613,-16.4592,-16.4584,-16.4571,-16.4596,-16.4596,-16.4638,-16.4638,-16.4659,-16.4701,-16.473,-16.473,-16.4738,-16.4738,-16.4775,-16.4792,-16.4805,-16.4805,-16.4796,-16.4797,-16.478,-16.478,-16.4805,-16.4813,-16.4834,-16.4884,-16.4905,-16.4905,-16.4884,-16.4842,-16.4813,-16.4796,-16.4776,-16.4717,-16.4642,-16.4592,-16.4567,-16.4542,-16.4525,-16.4513,-16.4522,-16.4538,-16.4538,-16.4546,-16.4542,-16.453,-16.4438,-16.4438,-16.4409,-16.4401,-16.4384,-16.4375,-16.4346,-16.4346,-16.4338,-16.4338,-16.4317,-16.4284,-16.4255,-16.4255,-16.4242,-16.4234,-16.4205,-16.4209,-16.4217,-16.423,-16.423,-16.428,-16.4271,-16.425,-16.4159,-16.415,-16.4134,-16.4126,-16.4059,-16.4017,-16.4009,-16.3934,-16.3918,-16.3892,-16.3792,-16.3776,-16.3767,-16.3717,-16.3709,-16.3651,-16.36,-16.3555,-16.3522,-16.3521,-16.3509,-16.3484,-16.3471,-16.3479,-16.3488,-16.3505,-16.353,-16.3521,-16.3521,-16.3505,-16.3496,-16.3497,-16.3467,-16.3434,-16.34,-16.3392,-16.335,-16.3292,-16.3284,-16.3267,-16.3259,-16.3234,-16.3225,-16.3134,-16.3125,-16.3068,-16.3042,-16.2992,-16.2959,-16.2942,-16.2909,-16.29,-16.2851,-16.2817,-16.2809,-16.275,-16.2684,-16.2663,-16.2663,-16.2688,-16.2688,-16.2676,-16.2642,-16.2617,-16.2584,-16.2576,-16.2542,-16.2509,-16.25,-16.2446,-16.2446,-16.2421,-16.2421,-16.2375,-16.2326,-16.2317,-16.2293,-16.2197,-16.2197,-16.2226,-16.2284,-16.2367,-16.24,-16.2417,-16.2467,-16.2484,-16.2501,-16.2542,-16.2592,-16.2675,-16.2684,-16.2709,-16.2751,-16.2801,-16.2826,-16.2851,-16.2867,-16.2909,-16.2918,-16.2938,-16.2955,-16.2955,-16.2976,-16.3009,-16.3017,-16.3051,-16.3067,-16.3159,-16.3176,-16.3192,-16.3225,-16.3233,-16.3259,-16.3284,-16.3309,-16.3326,-16.3359,-16.3421,-16.3429,-16.3463,-16.3463,-16.3471,-16.3471,-16.348,-16.348,-16.3505,-16.3567,-16.3609,-16.3642,-16.365,-16.3717,-16.3725,-16.3767,-16.3775,-16.3826,-16.3842,-16.3884,-16.3896,-16.3876,-16.3871,-16.3888,-16.3888,-16.3909,-16.3926,-16.3934,-16.3967,-16.3984,-16.4009,-16.4033,-16.4042,-16.4109,-16.4167,-16.4184,-16.4209,-16.425,-16.43,-16.4317,-16.4367,-16.4409,-16.4434,-16.4442,-16.4475,-16.4496,-16.4496,-16.4521,-16.4542,-16.4546,-16.4538,-16.4538,-16.4496,-16.4505,-16.4521,-16.4517,-16.4492,-16.4475,-16.4446,-16.4438,-16.4451,-16.4476,-16.4492,-16.4517,-16.4526,-16.453,-16.4555,-16.4555,-16.4584,-16.4625,-16.4634,-16.4709,-16.4717,-16.4809,-16.4817,-16.4842,-16.485,-16.4876,-16.4884,-16.4934,-16.4942,-16.4976,-16.5026,-16.5108,-16.5117,-16.5142,-16.5167,-16.5209,-16.5217,-16.5234,-16.5276,-16.5325,-16.5338,-16.5338,-16.533,-16.5321,-16.5309,-16.5293,-16.5246,-16.5259,-16.5301,-16.5305,-16.5296,-16.5296,-16.5288,-16.5288,-16.5242,-16.5184,-16.5175,-16.5158,-16.5117,-16.5076,-16.5067,-16.5038,-16.5105,-16.5113,-16.5134,-16.515,-16.5184,-16.52,-16.5217,-16.5247,-16.5251,-16.5263,-16.5263,-16.5254,-16.5255,-16.523,-16.523,-16.5251,-16.5275,-16.5301,-16.5317,-16.5321,-16.5321,-16.533,-16.533,-16.5346,-16.5347,-16.5355,-16.5355,-16.5371,-16.5371,-16.5396,-16.5396,-16.5413,-16.5413,-16.5384,-16.5371,-16.5363,-16.5347,-16.5346,-16.5372,-16.5372,-16.5492,-16.55,-16.5559,-16.5592,-16.5634,-16.5655,-16.5659,-16.5676,-16.5688,-16.5692,-16.5726,-16.575,-16.5775,-16.5784,-16.5817,-16.5851,-16.5867,-16.5913,-16.5922,-16.5921,-16.5913,-16.595,-16.5967,-16.6025,-16.6097,-16.6105,-16.6105,-16.6113,-16.6113,-16.6126,-16.615,-16.618,-16.6213,-16.6238,-16.6263,-16.6263,-16.6284,-16.633,-16.6329,-16.6397,-16.6397,-16.6413,-16.6413,-16.643,-16.6434,-16.6451,-16.6534,-16.6567,-16.6584,-16.6626,-16.6659,-16.6667,-16.67,-16.6755,-16.6755,-16.6788,-16.6788,-16.6771,-16.6771,-16.6738,-16.6738,-16.6746,-16.6755,-16.6763,-16.6763,-16.6705,-16.6696,-16.6663,-16.6613,-16.6646,-16.6634,-16.66,-16.6596,-16.6642,-16.6646,-16.6705,-16.6705,-16.6671,-16.6663,-16.665,-16.6647,-16.6667,-16.6701,-16.6718,-16.6776,-16.6834,-16.6847,-16.6792,-16.675,-16.6725,-16.6713,-16.6713,-16.6696,-16.6696,-16.668,-16.668,-16.6663,-16.6663,-16.6654,-16.6667,-16.6717,-16.6742,-16.6855,-16.6855,-16.6876,-16.6892,-16.6951,-16.6963,-16.6963,-16.6955,-16.6955,-16.6926,-16.6909,-16.6884,-16.6834,-16.6805,-16.6796,-16.6755,-16.6755,-16.6746,-16.6747,-16.6755,-16.6792,-16.6859,-16.69,-16.6917,-16.6921,-16.6901,-16.6875,-16.6863,-16.6917,-16.695,-16.6976,-16.7051,-16.7083,-16.7167,-16.7217,-16.7259,-16.7309,-16.7338,-16.7363,-16.7363,-16.7388,-16.7396,-16.7405,-16.7405,-16.7413,-16.7446,-16.7438,-16.7438,-16.7446,-16.7455,-16.748,-16.748,-16.7488,-16.75,-16.7517,-16.7534,-16.755,-16.7584,-16.7625,-16.768,-16.768,-16.7719,-16.778,-16.7805,-16.7805,-16.7822,-16.783,-16.7838,-16.7846,-16.7863,-16.7863,-16.788,-16.788,-16.7921,-16.7921,-16.793,-16.793,-16.7921,-16.7921,-16.7913,-16.7913,-16.7921,-16.7921,-16.7913,-16.7905,-16.7896,-16.7896,-16.7888,-16.7888,-16.7871,-16.7871,-16.7855,-16.7855,-16.7846,-16.7838,-16.783,-16.783,-16.7813,-16.7813,-16.7805,-16.7796,-16.7755,-16.7755,-16.7742,-16.7725,-16.7688,-16.768,-16.7679,-16.7667,-16.7634,-16.7597,-16.7584,-16.7521,-16.7513,-16.7513,-16.7501,-16.7459,-16.7442,-16.738,-16.7355,-16.7291,-16.7198,-16.7156,-16.7153,-16.7143,-16.7109,-16.7106,-16.685,-16.6692,-16.6692,-16.6692,-16.6678,-16.6675,-16.6662,-16.6655,-16.6639,-16.6656,-16.6665,-16.6666,-16.6665,-16.6445,-16.6378,-16.632,-16.6264,-16.6211,-16.617,-16.6106,-16.6036,-16.5967,-16.5903,-16.5831,-16.5767,-16.5706,-16.565,-16.5586,-16.5517,-16.5461,-16.5417,-16.5364,-16.5314,-16.5264,-16.5214,-16.5142,-16.5089,-16.5034,-16.4992,-16.4978,-16.4923,-16.4859,-16.48,-16.4748,-16.4692,-16.4628,-16.4578,-16.4517,-16.4453,-16.4384,-16.432,-16.4261,-16.4214,-16.4161,-16.4111,-16.4056,-16.3984,-16.3936,-16.3911,-16.3906,-16.3809,-16.3636,-16.3464,-16.3289,-16.312,-16.2945,-16.2773,-16.2603,-16.2428,-16.2328,-16.2259,-16.2084,-16.2011,-16.1834,-16.1653,-16.1634,-16.1592,-16.1542,-16.1503,-16.1461,-16.1406,-16.1356,-16.1309,-16.1261,-16.1206,-16.1142,-16.1078,-16.1014,-16.095,-16.0881,-16.0811,-16.0748,-16.0684,-16.062,-16.0564,-16.0506,-16.0448,-16.0392,-16.0342,-16.0289,-16.0234,-16.0184,-16.0142,-16.0098,-16.0053,-15.9998,-15.9992,-15.9945,-15.9895,-15.9845,-15.9795,-15.9745,-15.97,-15.965,-15.9614,-15.9564,-15.9392,-15.9217,-15.9039,-15.8878,-15.8864,-15.8792,-15.8614,-15.8439,-15.8261,-15.8086,-15.7906,-15.7728,-15.7553,-15.7375,-15.72,-15.702,-15.6845,-15.6831,-15.6656,-15.6484,-15.6311,-15.6139,-15.5967,-15.5792,-15.562,-15.5445,-15.5319,-15.5318,-15.5273,-15.51,-15.4992,-15.4925,-15.4753,-15.4578,-15.4406,-15.4236,-15.4186,-15.4061,-15.3945,-15.3775,-15.3606,-15.3587,-15.3434,-15.3336,-15.3164,-15.2986,-15.2814,-15.2684,-15.2514,-15.2342,-15.2173,-15.2117,-15.1939,-15.1761,-15.1586,-15.1409,-15.1364,-15.1231,-15.1056,-15.0878,-15.0703,-15.0525,-15.0348,-15.017,-14.9995,-14.9992,-14.9931,-14.9756,-14.9578,-14.9423,-14.94,-14.9225,-14.9048,-14.8873,-14.8736,-14.8692,-14.8517,-14.8342,-14.8164,-14.7986,-14.7809,-14.7634,-14.7459,-14.7281,-14.7103,-14.6925,-14.675,-14.6575,-14.6395,-14.622,-14.6045,-14.5867,-14.5689,-14.5511,-14.5336,-14.5161,-14.4992,-14.4981,-14.4806,-14.4631,-14.4453,-14.4275,-14.4098,-14.3923,-14.3745,-14.357,-14.3411,-14.3409,-14.3392,-14.3214,-14.3039,-14.2864,-14.2684,-14.2509,-14.2331,-14.2156,-14.1978,-14.18,-14.1625,-14.1589,-14.145,-14.1273,-14.1095,-14.0917,-14.0742,-14.0567,-14.0386,-14.0211,-14.0036,-13.9992,-13.9978,-13.98,-13.9625,-13.9448,-13.9273,-13.9095,-13.8906,-13.8567,-13.8389,-13.8211,-13.8034,-13.7859,-13.7684,-13.7503,-13.7328,-13.7153,-13.7123,-13.6975,-13.6798,-13.6623,-13.6445,-13.627,-13.6258,-13.6198,-13.602,-13.5845,-13.5667,-13.5517,-13.5489,-13.5314,-13.5136,-13.4992,-13.4961,-13.4781,-13.4606,-13.4468,-13.4428,-13.4253,-13.4075,-13.3898,-13.372,-13.3596,-13.3559,-13.3523,-13.345,-13.3406,-13.3361,-13.3278,-13.3225,-13.3186,-13.3148,-13.3092,-13.3069,-13.3036,-13.2984,-13.2906,-13.2861,-13.2817,-13.277,-13.2723,-13.2681,-13.2631,-13.2567,-13.2542,-13.2523,-13.2484,-13.2423,-13.235,-13.23,-13.2242,-13.2195,-13.2159,-13.2125,-13.2056,-13.2006,-13.197,-13.1925,-13.1898,-13.1859,-13.1817,-13.1759,-13.1734,-13.1692,-13.1636,-13.1578,-13.1509,-13.1453,-13.1414,-13.1386,-13.1353,-13.1317,-13.1292,-13.1261,-13.1253,-13.1225,-13.1175,-13.1117,-13.1053,-13.1006,-13.0948,-13.0884,-13.082,-13.0764,-13.07,-13.065,-13.0603,-13.0573,-13.0536,-13.0506,-13.0478,-13.0467,-13.0461,-13.045,-13.0428,-13.0411,-13.0409,-13.0423,-13.0434,-13.0461,-13.0481,-13.0506,-13.0525,-13.0548,-13.055,-13.057,-13.0606,-13.0639,-13.0673,-13.0686,-13.0684,-13.0678,-13.0656,-13.0653,-13.065,-13.0634,-13.0622,-13.0611,-13.0589,-13.0559,-13.0528,-13.0486,-13.0445,-13.04,-13.037,-13.03,-13.0245,-13.0181,-13.012,-13.0048,-12.9992,-12.9984,-12.9925,-12.9878,-12.9834,-12.9781,-12.9734,-12.9678,-12.9628,-12.9586,-12.9545,-12.9525,-12.952,-12.9509,-12.9503,-12.9501,-12.9492,-12.9486,-12.9475,-12.9461,-12.945,-12.9423,-12.9389,-12.9348,-12.9306,-12.9256,-12.9209,-12.9145,-12.9075,-12.9003,-12.8953,-12.8906,-12.887,-12.8836,-12.8809,-12.8778,-12.8742,-12.87,-12.8656,-12.8606,-12.855,-12.85,-12.845,-12.8392,-12.8328,-12.8292,-12.8298,-12.8319,-12.832,-12.8331,-12.8381,-12.8414,-12.8428,-12.8431,-12.8403,-12.8348,-12.8289,-12.8253,-12.8225,-12.8181,-12.8125,-12.8067,-12.7998,-12.7925,-12.7861,-12.7806,-12.7742,-12.7684,-12.7628,-12.7634,-12.7645,-12.7659,-12.7656,-12.7603,-12.7534,-12.747,-12.74,-12.7331,-12.7261,-12.7189,-12.7125,-12.7056,-12.6984,-12.692,-12.685,-12.6786,-12.6717,-12.6648,-12.6584,-12.6514,-12.645,-12.6378,-12.6323,-12.6286,-12.625,-12.622,-12.6192,-12.6173,-12.615,-12.6131,-12.6106,-12.6075,-12.6031,-12.597,-12.5911,-12.5875,-12.5873,-12.5895,-12.5906,-12.587,-12.5842,-12.582,-12.5781,-12.5717,-12.5675,-12.5634,-12.5592,-12.5548,-12.5498,-12.5442,-12.5392,-12.5345,-12.5303,-12.5253,-12.5203,-12.5156,-12.5106,-12.5073,-12.5017,-12.4992,-12.4961,-12.4903,-12.4842,-12.4784,-12.4725,-12.4678,-12.462,-12.4561,-12.4506,-12.445,-12.4398,-12.4356,-12.432,-12.427,-12.4228,-12.4164,-12.41,-12.4042,-12.3986,-12.3942,-12.3906,-12.3873,-12.3839,-12.3817,-12.3795,-12.3773,-12.375,-12.3728,-12.3698,-12.3661,-12.3617,-12.3575,-12.3548,-12.3523,-12.35,-12.347,-12.3448,-12.3434,-12.3386,-12.3339,-12.3289,-12.3248,-12.32,-12.315,-12.3103,-12.3048,-12.3009,-12.2998,-12.2948,-12.29,-12.2845,-12.2786,-12.2731,-12.2675,-12.2625,-12.2584,-12.2536,-12.2489,-12.2425,-12.2356,-12.2284,-12.2214,-12.2142,-12.2078,-12.2017,-12.1953,-12.1911,-12.1856,-12.1798,-12.175,-12.1695,-12.1645,-12.1595,-12.155,-12.1509,-12.1467,-12.1425,-12.1389,-12.1348,-12.1309,-12.1273,-12.1248,-12.1211,-12.1178,-12.1142,-12.1111,-12.107,-12.1034,-12.0978,-12.0909,-12.0845,-12.0773,-12.0703,-12.0631,-12.0567,-12.0511,-12.0461,-12.042,-12.0375,-12.0331,-12.0275,-12.0211,-12.0148,-12.0103,-12.0056,-12.0006,-11.9992,-11.9948,-11.9878,-11.9814,-11.975,-11.9681,-11.9625,-11.9584,-11.9542,-11.9498,-11.9453,-11.9411,-11.937,-11.932,-11.9278,-11.9228,-11.9189,-11.9148,-11.9106,-11.9073,-11.9042,-11.902,-11.8995,-11.8973,-11.8945,-11.8895,-11.8836,-11.8781,-11.8725,-11.8673,-11.8631,-11.8589,-11.8545,-11.8503,-11.845,-11.8409,-11.8348,-11.8289,-11.8234,-11.8178,-11.8123,-11.805,-11.7981,-11.7917,-11.7861,-11.7809,-11.7775,-11.7723,-11.7681,-11.7631,-11.7567,-11.7498,-11.7425,-11.7356,-11.73,-11.7242,-11.72,-11.7159,-11.7111,-11.7048,-11.6992,-11.6928,-11.6878,-11.6836,-11.6798,-11.6748,-11.6706,-11.6664,-11.6623,-11.6581,-11.6539,-11.6498,-11.6456,-11.6423,-11.6373,-11.6317,-11.6273,-11.6231,-11.6189,-11.6125,-11.6061,-11.6011,-11.5956,-11.59,-11.5839,-11.5775,-11.5711,-11.5648,-11.5586,-11.5514,-11.545,-11.5381,-11.5317,-11.5253,-11.5189,-11.512,-11.5056,-11.4992,-11.4986,-11.4923,-11.4853,-11.4789,-11.4725,-11.4661,-11.4603,-11.4548,-11.4492,-11.4436,-11.4378,-11.432,-11.4264,-11.4214,-11.4173,-11.4136,-11.4098,-11.4056,-11.4006,-11.395,-11.3898,-11.3842,-11.3792,-11.3731,-11.3723,-11.3714,-11.3725,-11.3734,-11.3742,-11.3745,-11.3739,-11.3728,-11.37,-11.3664,-11.3636,-11.3603,-11.3575,-11.3611,-11.3614,-11.3622,-11.3625,-11.3631,-11.3653,-11.3671,-11.3672,-11.3698,-11.3742,-11.3741,-11.3739,-11.3753,-11.3759,-11.3786,-11.3842,-11.3864,-11.3859,-11.3861,-11.3861,-11.3884,-11.3917,-11.3942,-11.395,-11.3967,-11.3995,-11.4036,-11.4081,-11.4123,-11.4123,-11.4106,-11.4106,-11.4136,-11.4189,-11.4225,-11.4256,-11.4284,-11.4298,-11.4309,-11.4331,-11.4373,-11.4436,-11.4478,-11.4509,-11.4489,-11.4439,-11.4392,-11.4336,-11.4286,-11.4236,-11.4195,-11.4161,-11.4156,-11.415,-11.4175,-11.42,-11.4234,-11.4256,-11.4284,-11.4292,-11.4273,-11.4253,-11.4217,-11.4216,-11.4206,-11.4195,-11.42,-11.4211,-11.4217,-11.4234,-11.4261,-11.4314,-11.4384,-11.4448,-11.4492,-11.4506,-11.4506,-11.4495,-11.4475,-11.4461,-11.4425,-11.4384,-11.4314,-11.4278,-11.4275,-11.4284,-11.4278,-11.4264,-11.4209,-11.4153,-11.4089,-11.4034,-11.3975,-11.3914,-11.3864,-11.3823,-11.3823,-11.3859,-11.3875,-11.387,-11.3886,-11.3923,-11.3964,-11.4009,-11.4009,-11.3959,-11.3925,-11.3906,-11.3892,-11.3895,-11.3911,-11.392,-11.3923,-11.3914,-11.3911,-11.3925,-11.3942,-11.3998,-11.4048,-11.4111,-11.4156,-11.417,-11.4156,-11.4131,-11.41,-11.4081,-11.4095,-11.412,-11.4136,-11.4123,-11.4103,-11.4061,-11.4034,-11.4009,-11.3986,-11.3989,-11.4011,-11.4048,-11.4084,-11.4114,-11.4109,-11.4067,-11.4017,-11.3967,-11.3923,-11.3867,-11.3803,-11.3753,-11.3711,-11.3673,-11.3667,-11.3673,-11.3684,-11.3692,-11.37,-11.3711,-11.3717,-11.3734,-11.375,-11.3773,-11.3809,-11.3856,-11.3878,-11.3884,-11.3881,-11.3892,-11.3914,-11.3939,-11.3992,-11.4034,-11.4084,-11.4139,-11.4189,-11.4225,-11.4206,-11.4175,-11.4164,-11.4175,-11.4203,-11.4217,-11.4214,-11.4185,-11.4184,-11.4175,-11.4173,-11.4131,-11.4081,-11.4056,-11.4075,-11.4125,-11.4184,-11.4214,-11.4228,-11.4236,-11.4253,-11.4289,-11.4313,-11.4331,-11.4403,-11.4425,-11.4428,-11.44,-11.4367,-11.4336,-11.4339,-11.437,-11.4414,-11.4478,-11.4548,-11.4598,-11.4656,-11.4703,-11.4739,-11.477,-11.4806,-11.485,-11.4881,-11.4914,-11.4959,-11.4992,-11.5028,-11.507,-11.5128,-11.5159,-11.5164,-11.5181,-11.5203,-11.5234,-11.527,-11.5292,-11.5306,-11.5286,-11.5253,-11.5234,-11.5223,-11.5214,-11.5217,-11.5217,-11.5248,-11.5289,-11.5353,-11.542,-11.5461,-11.5506,-11.5534,-11.555,-11.5567,-11.5573,-11.5561,-11.5542,-11.5514,-11.5486,-11.5461,-11.5434,-11.5406,-11.537,-11.5359,-11.5381,-11.5409,-11.5445,-11.547,-11.5492,-11.5514,-11.5545,-11.5595,-11.565,-11.5686,-11.5717,-11.5748,-11.5784,-11.5828,-11.5864,-11.587,-11.5881,-11.5895,-11.5925,-11.5953,-11.597,-11.597,-11.5967,-11.597,-11.597,-11.6011,-11.6039,-11.6073,-11.6131,-11.6189,-11.6245,-11.6289,-11.6311,-11.6298,-11.6292,-11.6292,-11.6317,-11.6359,-11.6423,-11.6473,-11.6509,-11.6542,-11.6606,-11.667,-11.6725,-11.6789,-11.6853,-11.6903,-11.6939,-11.6961,-11.6984,-11.7006,-11.7042,-11.7075,-11.7123,-11.7186,-11.725,-11.732,-11.7384,-11.742,-11.7434,-11.7439,-11.7459,-11.7478,-11.7503,-11.7531,-11.7559,-11.7581,-11.7592,-11.7598,-11.7586,-11.7592,-11.7634,-11.7689,-11.7756,-11.7806,-11.7839,-11.7864,-11.7886,-11.7911,-11.7939,-11.7973,-11.8014,-11.8056,-11.8106,-11.8156,-11.8192,-11.8223,-11.8239,-11.8261,-11.8275,-11.8278,-11.8289,-11.8289,-11.8325,-11.8384,-11.8431,-11.8489,-11.8539,-11.8584,-11.8611,-11.8642,-11.8673,-11.8714,-11.8786,-11.8842,-11.8873,-11.8873,-11.8873,-11.8853,-11.8834,-11.8811,-11.8792,-11.8789,-11.8784,-11.8764,-11.8731,-11.8695,-11.8667,-11.8653,-11.8656,-11.8673,-11.8673,-11.8692,-11.8711,-11.8742,-11.8786,-11.885,-11.892,-11.8975,-11.902,-11.9056,-11.9073,-11.9089,-11.9103,-11.9116,-11.9136,-11.9175,-11.9228,-11.9264,-11.93,-11.9342,-11.9375,-11.9409,-11.9459,-11.9517,-11.9586,-11.9639,-11.9681,-11.9711,-11.9734,-11.9734,-11.9761,-11.9792,-11.982,-11.9859,-11.99,-11.9931,-11.9961,-11.9989,-11.9992,-12.002,-12.0036,-12.005,-12.0067,-12.0098,-12.0125,-12.0167,-12.0211,-12.0256,-12.03,-12.0336,-12.035,-12.0367,-12.0348,-12.0334,-12.0336,-12.0359,-12.0389,-12.0425,-12.0461,-12.0498,-12.0534,-12.0564,-12.0586,-12.0609,-12.0617,-12.0611,-12.0614,-12.0636,-12.0673,-12.0717,-12.0739,-12.0753,-12.0758,-12.081,-12.0811,-12.0814,-12.0784,-12.0734,-12.0692,-12.0645,-12.0595,-12.0553,-12.0506,-12.0453,-12.0406,-12.0364,-12.0323,-12.0273,-12.0231,-12.0184,-12.0142,-12.0092,-12.005,-12.0009,-11.9992,-11.9967,-11.9925,-11.9889,-11.9848,-11.9814,-11.9781,-11.9739,-11.9703,-11.967,-11.9636,-11.9609,-11.9578,-11.9559,-11.9539,-11.952,-11.9509,-11.9498,-11.9484,-11.947,-11.9464,-11.9453,-11.9439,-11.9428,-11.9417,-11.9417,-11.9411,-11.9414,-11.9409,-11.9411,-11.9411,-11.9414,-11.9436,-11.9467,-11.9503,-11.9545,-11.9581,-11.962,-11.9656,-11.97,-11.9736,-11.9781,-11.9823,-11.9859,-11.99,-11.9939,-11.9975,-11.9992,-12.0003,-12.0042,-12.007,-12.0086,-12.0089,-12.0089,-12.0092,-12.0095,-12.01,-12.0103,-12.0125,-12.0142,-12.0136,-12.0123,-12.0098,-12.007,-12.005,-12.0053,-12.0061,-12.0056,-12.0036,-12.0017,-11.9992,-11.9975,-11.9939,-11.9903,-11.987,-11.985,-11.9845,-11.9856,-11.9878,-11.99,-11.9909,-11.99,-11.9881,-11.9861,-11.9842,-11.9817,-11.9798,-11.9792,-11.9786,-11.9786,-11.9798,-11.9814,-11.9836,-11.9859,-11.9881,-11.9898,-11.9911,-11.9928,-11.9959,-11.9986,-11.9992,-12.0009,-12.0034,-12.0056,-12.0078,-12.0095,-12.0117,-12.0134,-12.0148,-12.017,-12.0186,-12.0209,-12.0239,-12.027,-12.03,-12.0336,-12.04,-12.0456,-12.05,-12.057,-12.0628,-12.0686,-12.0728,-12.0778,-12.0823,-12.0867,-12.0909,-12.0945,-12.0975,-12.0989,-12.102,-12.105,-12.1086,-12.1123,-12.1117,-12.1061,-12.1003,-12.0961,-12.097,-12.0986,-12.0959,-12.0917,-12.0898,-12.0911,-12.0942,-12.0978,-12.1031,-12.1092,-12.115,-12.1214,-12.1245,-12.1267,-12.1295,-12.1353,-12.1411,-12.1453,-12.15,-12.1545,-12.1581,-12.1611,-12.1648,-12.1695,-12.1761,-12.1811,-12.1856,-12.1898,-12.1928,-12.2056,-12.2045,-12.2042,-12.2041,-12.2039,-12.2031,-12.2031,-12.2039,-12.204,-12.2042,-12.2036,-12.2017,-12.2003,-12.1998,-12.2014,-12.2045,-12.2061,-12.2067,-12.2081,-12.2084,-12.2086,-12.2117,-12.2159,-12.2198,-12.2225,-12.2242,-12.2249,-12.2249,-12.2253,-12.2245,-12.222,-12.2192,-12.2173,-12.2159,-12.2156,-12.2161,-12.2181,-12.2159,-12.2103,-12.2053,-12.1981,-12.1925,-12.1884,-12.1861,-12.185,-12.1867,-12.1878,-12.1878,-12.1856,-12.1823,-12.1789,-12.1784,-12.1778,-12.175,-12.1717,-12.1667,-12.1603,-12.1561,-12.1528,-12.1498,-12.1478,-12.1467,-12.147,-12.1486,-12.1514,-12.1559,-12.1609,-12.1678,-12.1736,-12.1778,-12.1803,-12.1814,-12.177,-12.1714,-12.1678,-12.1653,-12.1634,-12.1642,-12.1706,-12.1739,-12.1781,-12.1811,-12.182,-12.1873,-12.1928,-12.1984,-12.205,-12.2092,-12.2117,-12.2139,-12.2148,-12.2136,-12.2123,-12.2153,-12.22,-12.2253,-12.2311,-12.2361,-12.2398,-12.2434,-12.2442,-12.2442,-12.2441,-12.2439,-12.2439,-12.2475,-12.252,-12.2567,-12.2611,-12.2656,-12.2706,-12.2745,-12.2784,-12.2828,-12.2864,-12.2903,-12.2945,-12.2984,-12.302,-12.3056,-12.3092,-12.3125,-12.3162,-12.3189,-12.322,-12.3256,-12.33,-12.3348,-12.3392,-12.3445,-12.3503,-12.3559,-12.3603,-12.3659],"lat":[14.8396,14.8415,14.8427,14.8454,14.8474,14.8499,14.8527,14.856,14.8593,14.8627,14.8665,14.8699,14.8746,14.8779,14.8813,14.8849,14.8874,14.8907,14.894,14.8974,14.9007,14.9046,14.9052,14.9096,14.9149,14.9213,14.9279,14.9335,14.939,14.9454,14.9507,14.9549,14.9596,14.9638,14.9685,14.9732,14.9779,14.9835,14.9888,14.9943,14.9996,14.9999,15.0046,15.0085,15.0118,15.0152,15.0171,15.019,15.0199,15.0204,15.0221,15.0235,15.0254,15.0274,15.0293,15.0318,15.0338,15.0357,15.0377,15.0402,15.0421,15.0446,15.0474,15.0507,15.0532,15.0565,15.0593,15.0627,15.066,15.0693,15.0727,15.0768,15.0807,15.084,15.0882,15.0907,15.0935,15.096,15.0985,15.1004,15.1018,15.1029,15.104,15.1054,15.1052,15.1043,15.1021,15.0985,15.0943,15.0907,15.0877,15.0871,15.087,15.0543,15.0273,15.0004,14.9973,14.9955,14.9933,14.9916,14.9905,14.9882,14.9872,14.9861,14.985,14.9837,14.9829,14.9824,14.9822,14.9813,14.9808,14.98,14.9794,14.9783,14.9775,14.9766,14.9759,14.9751,14.9742,14.973,14.9726,14.9724,14.9719,14.9703,14.969,14.9682,14.9675,14.9666,14.9656,14.965,14.9647,14.964,14.9634,14.9632,14.9628,14.9623,14.9608,14.9593,14.958,14.9568,14.956,14.9552,14.953,14.9521,14.9512,14.9502,14.9492,14.9478,14.9466,14.9454,14.944,14.943,14.9412,14.9401,14.9386,14.9375,14.9363,14.9356,14.9341,14.9333,14.9321,14.9308,14.9294,14.9279,14.9263,14.9245,14.9229,14.922,14.9206,14.9194,14.9179,14.9169,14.9151,14.9139,14.9129,14.912,14.911,14.91,14.9091,14.9079,14.9071,14.9057,14.9052,14.8798,14.8761,14.8474,14.8188,14.818,14.7902,14.7616,14.7329,14.7087,14.7043,14.6816,14.6589,14.6362,14.6135,14.6025,14.5915,14.5804,14.5767,14.573,14.5694,14.5657,14.562,14.5583,14.5592,14.5654,14.5726,14.5764,14.5803,14.5847,14.5849,14.5892,14.5936,14.5981,14.6025,14.6064,14.6069,14.6114,14.5999,14.5884,14.577,14.5655,14.554,14.5425,14.5311,14.5196,14.5081,14.4967,14.4852,14.4737,14.4622,14.4507,14.4393,14.4278,14.4163,14.4048,14.4052,14.4055,14.4062,14.4069,14.4076,14.4083,14.4089,14.4094,14.4096,14.4103,14.4109,14.4109,14.4109,14.411,14.411,14.4116,14.4117,14.4123,14.4246,14.437,14.4493,14.4617,14.4696,14.4741,14.484,14.4859,14.4864,14.4911,14.493,14.4973,14.4984,14.5145,14.5304,14.531,14.531,14.5335,14.5385,14.5385,14.5395,14.5404,14.5412,14.5578,14.5578,14.5588,14.5588,14.5709,14.573,14.5747,14.5758,14.5845,14.586,14.5864,14.5864,14.5864,14.5864,14.5838,14.5838,14.5823,14.5823,14.5812,14.5801,14.58,14.5799,14.5791,14.5802,14.5807,14.579,14.5787,14.5787,14.5779,14.5776,14.577,14.5743,14.5724,14.5702,14.5675,14.5665,14.5658,14.5655,14.5654,14.5651,14.5647,14.564,14.5639,14.5634,14.563,14.5626,14.562,14.562,14.5619,14.5619,14.5621,14.5621,14.5623,14.5627,14.5632,14.5628,14.5623,14.562,14.5604,14.561,14.561,14.5627,14.5638,14.5667,14.5678,14.568,14.5717,14.5742,14.5745,14.5746,14.5747,14.5738,14.5735,14.5745,14.5754,14.5751,14.5733,14.572,14.5703,14.5672,14.5649,14.5623,14.542,14.5282,14.5153,14.4885,14.4617,14.4349,14.4081,14.3813,14.3546,14.3278,14.301,14.2986,14.2974,14.2958,14.2919,14.2906,14.2895,14.2868,14.2832,14.2812,14.2769,14.274,14.2721,14.2699,14.2676,14.2649,14.2623,14.2601,14.2592,14.2581,14.2552,14.2537,14.2514,14.2484,14.2457,14.2431,14.2402,14.2381,14.2348,14.231,14.2282,14.2253,14.2243,14.2213,14.2192,14.2173,14.2154,14.213,14.2091,14.2074,14.2057,14.2026,14.2012,14.1989,14.1971,14.1958,14.1944,14.1928,14.1909,14.1884,14.1859,14.1834,14.1812,14.1784,14.1768,14.1748,14.1721,14.1692,14.1663,14.1638,14.162,14.161,14.1581,14.1551,14.1527,14.1493,14.1465,14.144,14.1416,14.1396,14.1385,14.1378,14.1367,14.1351,14.1325,14.13,14.1278,14.1251,14.1234,14.1214,14.1181,14.1172,14.1143,14.1118,14.1102,14.1076,14.1058,14.1036,14.1021,14.0995,14.0976,14.0959,14.0948,14.0909,14.0892,14.0873,14.0854,14.0837,14.0824,14.0794,14.0773,14.0743,14.0723,14.0695,14.0676,14.0629,14.0617,14.0591,14.0572,14.0537,14.0517,14.0494,14.0462,14.0436,14.0408,14.039,14.0344,14.0307,14.0269,14.0249,14.0221,14.0192,14.0162,14.0127,14.0111,14.0097,14.0079,14.0066,14.0041,14.0022,14.0014,14.0007,14.0001,13.9984,13.9957,13.994,13.991,13.9887,13.9863,13.9855,13.9834,13.9811,13.9787,13.9761,13.9724,13.9698,13.9663,13.963,13.9599,13.9564,13.9522,13.95,13.9473,13.9436,13.9412,13.9382,13.9374,13.9366,13.9335,13.9312,13.9286,13.9263,13.9242,13.9211,13.9188,13.9166,13.9137,13.9111,13.9093,13.9077,13.9051,13.9024,13.9,13.8974,13.8939,13.8913,13.8905,13.8887,13.8858,13.8836,13.8791,13.8714,13.8651,13.862,13.8602,13.8549,13.8488,13.8436,13.8428,13.8398,13.8366,13.8331,13.8304,13.827,13.8232,13.8208,13.8174,13.8162,13.8152,13.8145,13.814,13.813,13.8126,13.812,13.8117,13.8115,13.8112,13.8109,13.8109,13.8107,13.8102,13.8098,13.8095,13.8094,13.8078,13.8064,13.8046,13.8029,13.8004,13.7979,13.7962,13.7935,13.7913,13.7881,13.7849,13.7829,13.7673,13.766,13.7629,13.7596,13.7563,13.7529,13.7496,13.7449,13.7407,13.736,13.7313,13.7265,13.7218,13.7171,13.7118,13.7063,13.7007,13.6954,13.6893,13.6838,13.6782,13.6713,13.6652,13.6585,13.6529,13.6432,13.6329,13.6227,13.6121,13.6115,13.6143,13.6179,13.6215,13.6252,13.6288,13.6318,13.6352,13.6382,13.641,13.644,13.6468,13.6499,13.6529,13.6552,13.6574,13.6596,13.6618,13.6635,13.6649,13.6665,13.6674,13.6685,13.6693,13.6696,13.6699,13.6688,13.6677,13.6665,13.6652,13.664,13.6624,13.661,13.6593,13.6565,13.6546,13.6529,13.6502,13.6477,13.6443,13.6418,13.639,13.6365,13.634,13.6307,13.6279,13.6246,13.6234,13.6213,13.6165,13.6127,13.6079,13.6032,13.5985,13.5935,13.5882,13.5821,13.5765,13.5704,13.5643,13.5579,13.5518,13.5457,13.5396,13.534,13.5329,13.5318,13.5304,13.5279,13.5263,13.526,13.5221,13.5174,13.5132,13.5093,13.506,13.5027,13.4997,13.4993,13.4974,13.494,13.4915,13.489,13.4863,13.4829,13.4804,13.4779,13.4752,13.4727,13.4699,13.4674,13.4654,13.4629,13.461,13.4593,13.4574,13.4554,13.4535,13.4538,13.456,13.4582,13.4599,13.4615,13.4629,13.4646,13.4663,13.4685,13.4707,13.4729,13.4738,13.4749,13.4763,13.4774,13.4777,13.4779,13.4807,13.4843,13.4885,13.4927,13.4963,13.4996,13.4997,13.4999,13.5035,13.5063,13.5093,13.5121,13.5143,13.5168,13.5196,13.5218,13.5235,13.5257,13.5265,13.5282,13.5299,13.5307,13.5324,13.5327,13.5335,13.5343,13.5352,13.5352,13.5379,13.5415,13.5396,13.5424,13.546,13.5488,13.551,13.554,13.5563,13.5585,13.5602,13.5618,13.5627,13.5629,13.5632,13.5627,13.5615,13.5596,13.5577,13.5557,13.5574,13.5602,13.5638,13.5665,13.5696,13.5727,13.5754,13.5764,13.5777,13.5804,13.5821,13.5824,13.5827,13.5807,13.5788,13.5763,13.5738,13.5704,13.5699,13.5693,13.5696,13.5693,13.5679,13.5663,13.5635,13.561,13.5582,13.5552,13.5518,13.5485,13.5452,13.541,13.5363,13.5315,13.526,13.5207,13.5152,13.511,13.5091,13.5077,13.5035,13.4997,13.4996,13.4993,13.494,13.4874,13.4802,13.4715,13.4621,13.4579,13.4554,13.4474,13.4382,13.426,13.4193,13.4127,13.406,13.3888,13.3793,13.3674,13.3593,13.3521,13.346,13.3393,13.3327,13.3274,13.3221,13.3188,13.3174,13.3174,13.3174,13.3188,13.3193,13.3215,13.3207,13.316,13.3107,13.306,13.3021,13.301,13.3002,13.2988,13.3002,13.3007,13.294,13.2902,13.2868,13.2854,13.2846,13.2868,13.2807,13.2713,13.2635,13.2568,13.2513,13.246,13.2421,13.2379,13.234,13.2331,13.2313,13.2307,13.2302,13.2321,13.2346,13.2379,13.2427,13.246,13.2435,13.2421,13.2393,13.2393,13.2379,13.2368,13.2374,13.2402,13.2427,13.246,13.2507,13.256,13.2613,13.2668,13.2721,13.2788,13.286,13.2949,13.3027,13.306,13.3049,13.3049,13.3066,13.3068,13.3093,13.3135,13.3215,13.3282,13.334,13.3354,13.3415,13.3488,13.356,13.3621,13.3588,13.354,13.3502,13.3474,13.3449,13.344,13.3435,13.3449,13.3468,13.3493,13.3521,13.356,13.3602,13.3654,13.3715,13.3788,13.3868,13.3954,13.404,13.414,13.4227,13.4254,13.424,13.4254,13.4274,13.4315,13.4349,13.4407,13.4488,13.4549,13.4549,13.4549,13.4549,13.456,13.4574,13.4602,13.464,13.4682,13.4715,13.4768,13.4821,13.4882,13.4927,13.4979,13.4993,13.4996,13.5035,13.5093,13.5146,13.5193,13.5246,13.5307,13.5374,13.544,13.5507,13.5574,13.5646,13.5721,13.5799,13.5879,13.5954,13.5935,13.5907,13.5874,13.5827,13.5746,13.5654,13.5574,13.5524,13.5493,13.5435,13.5374,13.5307,13.5188,13.5068,13.4998,13.4996,13.494,13.4807,13.4688,13.4615,13.4507,13.4388,13.4288,13.4202,13.4135,13.4054,13.3993,13.3935,13.3882,13.3827,13.3782,13.374,13.3707,13.3688,13.3668,13.3654,13.3656,13.366,13.366,13.3682,13.3693,13.3715,13.3735,13.3749,13.3774,13.3802,13.3835,13.3854,13.3888,13.3921,13.396,13.3968,13.3974,13.394,13.3888,13.3819,13.3807,13.3749,13.3693,13.364,13.3607,13.3593,13.3582,13.3574,13.3574,13.3582,13.3607,13.3621,13.3649,13.3682,13.3715,13.3674,13.3615,13.356,13.3521,13.3493,13.3482,13.3468,13.346,13.346,13.3454,13.3454,13.3468,13.3474,13.3307,13.3135,13.2968,13.2793,13.2621,13.2449,13.2446,13.2274,13.2102,13.1927,13.176,13.1593,13.1593,13.1602,13.1593,13.1602,13.1602,13.1602,13.1602,13.1607,13.1613,13.1613,13.1621,13.1621,13.1627,13.1627,13.1627,13.1627,13.1627,13.1621,13.1627,13.1627,13.1627,13.1627,13.1635,13.1635,13.1635,13.1627,13.1627,13.1635,13.1635,13.1635,13.1635,13.164,13.1635,13.1635,13.1635,13.164,13.164,13.164,13.1635,13.1635,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.1646,13.164,13.164,13.164,13.164,13.1607,13.1513,13.1407,13.1307,13.1232,13.1152,13.1093,13.1027,13.0921,13.0899,13.0883,13.0813,13.0806,13.0788,13.0774,13.0754,13.0738,13.0732,13.0727,13.0727,13.0726,13.073,13.0731,13.0735,13.0735,13.0734,13.0701,13.068,13.0657,13.0648,13.0636,13.0635,13.0617,13.0559,13.0527,13.0521,13.0429,13.0396,13.0346,13.0338,13.0213,13.0205,12.9996,12.9987,12.9871,12.9862,12.9796,12.9788,12.9746,12.9738,12.9712,12.9704,12.9662,12.9654,12.9496,12.9488,12.9471,12.9429,12.9404,12.9371,12.9354,12.9338,12.9287,12.9279,12.9271,12.9221,12.9212,12.9162,12.9154,12.9129,12.9121,12.9096,12.9087,12.9062,12.9054,12.9004,12.8996,12.8971,12.8962,12.8871,12.8862,12.8821,12.8754,12.8712,12.8646,12.8579,12.8546,12.8455,12.8446,12.8413,12.8379,12.8329,12.8296,12.8196,12.8188,12.8137,12.8129,12.8096,12.8088,12.8005,12.7997,12.7996,12.7927,12.7921,12.7912,12.7871,12.7862,12.7837,12.7829,12.7779,12.777,12.7763,12.7755,12.7721,12.7713,12.7696,12.7671,12.7654,12.7642,12.7642,12.7654,12.7688,12.7696,12.7738,12.7746,12.7762,12.7804,12.7829,12.7838,12.7896,12.7904,12.7938,12.7946,12.7996,12.8004,12.8054,12.8062,12.8129,12.8138,12.8221,12.8271,12.8296,12.8363,12.8371,12.8388,12.8387,12.8371,12.8363,12.8354,12.8313,12.8287,12.8263,12.8254,12.8121,12.8088,12.8021,12.7992,12.7992,12.8009,12.8009,12.805,12.805,12.8058,12.8059,12.8067,12.8116,12.8142,12.8142,12.8129,12.8117,12.8084,12.8066,12.8067,12.805,12.805,12.8033,12.8025,12.8009,12.8,12.8,12.7992,12.7992,12.7975,12.7954,12.7912,12.7892,12.7892,12.7879,12.7859,12.783,12.7821,12.7804,12.7771,12.7758,12.7754,12.7746,12.7742,12.7729,12.7629,12.7621,12.7588,12.7571,12.7537,12.7504,12.7487,12.7471,12.7429,12.7421,12.7388,12.7379,12.7362,12.7296,12.7287,12.7154,12.7096,12.7079,12.7046,12.7034,12.7034,12.7046,12.7079,12.7096,12.71,12.7142,12.7142,12.7204,12.7229,12.7238,12.7254,12.7279,12.7316,12.7325,12.735,12.7387,12.7404,12.7413,12.7458,12.7409,12.7409,12.7425,12.7425,12.74,12.7392,12.7396,12.7404,12.7437,12.7471,12.7492,12.7467,12.7458,12.7387,12.7354,12.7346,12.7338,12.7296,12.7279,12.7275,12.7292,12.7291,12.7309,12.7358,12.7392,12.7391,12.7338,12.7329,12.7313,12.7279,12.7271,12.7246,12.7237,12.7233,12.7221,12.7187,12.7154,12.715,12.7167,12.7175,12.7146,12.7104,12.7095,12.7021,12.6987,12.6971,12.6946,12.6929,12.6904,12.6888,12.6829,12.677,12.6754,12.672,12.6712,12.6646,12.6638,12.6612,12.6571,12.6529,12.6509,12.6508,12.6517,12.655,12.6558,12.6584,12.6583,12.6591,12.6579,12.6538,12.6508,12.6525,12.6513,12.6475,12.6479,12.6492,12.6492,12.6471,12.6454,12.6446,12.6379,12.6363,12.6329,12.6321,12.6313,12.6287,12.6279,12.6246,12.6238,12.6187,12.6179,12.6121,12.6112,12.5964,12.5938,12.5917,12.592,12.5954,12.5979,12.5979,12.5937,12.5913,12.5895,12.5879,12.5854,12.5771,12.5764,12.5763,12.5746,12.5729,12.5712,12.5667,12.5642,12.5634,12.5633,12.5675,12.5687,12.5704,12.5738,12.5754,12.58,12.5796,12.5754,12.5734,12.5733,12.5792,12.5783,12.5791,12.5792,12.5825,12.5825,12.5833,12.5842,12.585,12.585,12.5867,12.5884,12.5884,12.5904,12.5912,12.5925,12.5925,12.5933,12.5934,12.5938,12.5954,12.5996,12.6004,12.6045,12.6054,12.6121,12.6125,12.6113,12.6084,12.6075,12.6075,12.6067,12.6042,12.5975,12.5975,12.6033,12.6042,12.6054,12.6079,12.6096,12.6113,12.6142,12.6142,12.6133,12.6134,12.6137,12.6158,12.6167,12.6167,12.6129,12.6109,12.6109,12.6125,12.6125,12.6183,12.6217,12.6217,12.6267,12.6267,12.6321,12.6329,12.6358,12.6358,12.6371,12.6375,12.6358,12.6358,12.637,12.6446,12.6454,12.6496,12.6546,12.6571,12.6579,12.6621,12.6679,12.6713,12.6721,12.6821,12.6829,12.6938,12.6969,12.6988,12.6993,12.6976,12.6937,12.693,12.6904,12.6887,12.6862,12.685,12.685,12.6846,12.6829,12.6825,12.6812,12.6696,12.6687,12.6646,12.6583,12.6575,12.6575,12.6567,12.6559,12.655,12.655,12.6559,12.6567,12.6621,12.6629,12.6654,12.6675,12.6642,12.6625,12.6625,12.6634,12.6654,12.6662,12.6696,12.6738,12.6796,12.6809,12.6821,12.6829,12.6892,12.69,12.6917,12.6925,12.697,12.7012,12.7033,12.6971,12.6963,12.6929,12.6921,12.6879,12.6854,12.6838,12.6834,12.6863,12.6913,12.6917,12.6871,12.6845,12.6842,12.6821,12.6763,12.6691,12.6684,12.665,12.6641,12.6583,12.6584,12.6571,12.655,12.6467,12.6416,12.6417,12.6392,12.6392,12.6354,12.6346,12.6275,12.6271,12.6262,12.6213,12.6204,12.6166,12.6158,12.6171,12.6213,12.6221,12.6271,12.6304,12.6395,12.643,12.6446,12.6488,12.6496,12.6516,12.6533,12.6563,12.6613,12.6629,12.6754,12.6762,12.6779,12.6809,12.6812,12.6838,12.6896,12.6921,12.6938,12.6954,12.7046,12.7059,12.7059,12.7079,12.7096,12.7121,12.715,12.7158,12.7137,12.7113,12.7092,12.7087,12.7046,12.7038,12.7016,12.7025,12.7025,12.6995,12.6987,12.6942,12.6942,12.6959,12.6954,12.6937,12.6921,12.6896,12.6867,12.6866,12.6858,12.6837,12.6813,12.6758,12.6734,12.6725,12.6725,12.6705,12.6671,12.6658,12.6658,12.6684,12.6709,12.6696,12.6679,12.6638,12.6596,12.6567,12.6567,12.6587,12.6613,12.6634,12.6634,12.6621,12.6588,12.658,12.6529,12.6521,12.6509,12.6508,12.6546,12.6604,12.6613,12.6629,12.6675,12.6675,12.6663,12.6637,12.6629,12.6613,12.6596,12.6546,12.6504,12.6438,12.6417,12.64,12.6379,12.6321,12.63,12.63,12.6263,12.6221,12.62,12.6167,12.6167,12.6134,12.6125,12.6109,12.6109,12.6121,12.6171,12.6188,12.6204,12.6212,12.6234,12.6204,12.6104,12.6096,12.6067,12.6067,12.6042,12.6042,12.6013,12.5996,12.5987,12.5971,12.5942,12.5942,12.5987,12.6012,12.6017,12.6,12.5987,12.5975,12.5975,12.5963,12.5912,12.5855,12.5821,12.58,12.58,12.5792,12.5792,12.5784,12.5783,12.575,12.575,12.5675,12.5667,12.5667,12.5617,12.5617,12.5608,12.56,12.5592,12.5591,12.5617,12.5671,12.5738,12.5796,12.5809,12.5817,12.5837,12.5913,12.5921,12.5962,12.5996,12.6012,12.6046,12.6046,12.6029,12.6004,12.5983,12.5984,12.6,12.6,12.6042,12.6058,12.6067,12.6067,12.6075,12.6075,12.6067,12.6067,12.6075,12.6075,12.6058,12.605,12.6025,12.6025,12.5992,12.5992,12.595,12.5933,12.5925,12.5925,12.5983,12.5979,12.5971,12.5946,12.5912,12.5892,12.5884,12.59,12.59,12.5909,12.5917,12.5942,12.5942,12.5996,12.6004,12.6029,12.6038,12.6083,12.6117,12.6125,12.6108,12.6104,12.6096,12.6067,12.6033,12.595,12.5933,12.5933,12.5908,12.5892,12.5892,12.5867,12.5842,12.585,12.5859,12.5859,12.5875,12.5875,12.5859,12.5858,12.5875,12.5883,12.5892,12.5904,12.5938,12.5954,12.5975,12.5983,12.5992,12.5992,12.5984,12.5883,12.5875,12.5867,12.5867,12.5858,12.5859,12.5833,12.5834,12.585,12.585,12.5771,12.5746,12.5679,12.5662,12.5654,12.5621,12.5612,12.5596,12.5546,12.5483,12.5459,12.545,12.5442,12.5442,12.545,12.545,12.5442,12.5442,12.5459,12.5459,12.5471,12.5475,12.5487,12.5513,12.5529,12.555,12.555,12.5567,12.5584,12.5583,12.5609,12.5616,12.5625,12.5625,12.5592,12.5592,12.5625,12.5633,12.5659,12.5659,12.5684,12.5717,12.5717,12.5725,12.5725,12.5704,12.5679,12.5638,12.5625,12.5646,12.5654,12.5671,12.5721,12.5746,12.5762,12.5767,12.575,12.575,12.5771,12.5788,12.5792,12.5809,12.5809,12.5833,12.5833,12.5846,12.5862,12.5871,12.59,12.59,12.5892,12.5892,12.59,12.59,12.5909,12.5908,12.5917,12.5917,12.5925,12.5925,12.5934,12.5942,12.5966,12.5967,12.5959,12.5958,12.5942,12.5942,12.5934,12.5934,12.59,12.5875,12.5862,12.5829,12.5821,12.5771,12.5758,12.575,12.5746,12.5734,12.5734,12.5704,12.5696,12.5671,12.5662,12.5562,12.5517,12.5517,12.5525,12.5525,12.5558,12.5559,12.555,12.5546,12.5487,12.5454,12.5434,12.5425,12.5425,12.545,12.545,12.5404,12.5383,12.5404,12.5421,12.5429,12.5446,12.5462,12.5479,12.5492,12.55,12.5525,12.5517,12.5521,12.5563,12.5571,12.5587,12.5604,12.5621,12.5629,12.5663,12.5687,12.5704,12.5754,12.5913,12.5929,12.5938,12.5942,12.5963,12.6004,12.6029,12.6054,12.6096,12.6104,12.6225,12.6225,12.6284,12.63,12.63,12.6279,12.6267,12.6267,12.6279,12.6333,12.6342,12.6358,12.6358,12.6367,12.6367,12.635,12.635,12.6304,12.6288,12.6262,12.6246,12.6209,12.62,12.62,12.6129,12.6113,12.6079,12.6071,12.6037,12.6025,12.6017,12.5987,12.5921,12.588,12.5854,12.5845,12.5833,12.5788,12.5779,12.5713,12.5704,12.5687,12.5679,12.567,12.5659,12.565,12.565,12.5667,12.5667,12.5617,12.5617,12.5608,12.5591,12.5546,12.5538,12.5504,12.5479,12.5454,12.5438,12.5405,12.5371,12.5363,12.5237,12.5229,12.5162,12.5096,12.5054,12.4996,12.4996,12.4963,12.495,12.495,12.4938,12.4925,12.4913,12.4854,12.4846,12.4813,12.4787,12.4784,12.4771,12.4767,12.4792,12.4792,12.4733,12.4733,12.4712,12.4683,12.4642,12.4642,12.4629,12.4613,12.4596,12.4579,12.4563,12.4546,12.4529,12.4512,12.4504,12.4492,12.4542,12.4542,12.4654,12.4662,12.4683,12.4683,12.4625,12.4629,12.4671,12.4679,12.4696,12.4725,12.4725,12.4742,12.4767,12.4796,12.4821,12.4871,12.4904,12.4913,12.4962,12.4979,12.5,12.5,12.4975,12.4975,12.4979,12.5,12.5008,12.5038,12.5042,12.5025,12.5009,12.5008,12.5025,12.5058,12.5058,12.5092,12.5092,12.5137,12.5162,12.5171,12.5196,12.5221,12.5229,12.5271,12.5287,12.5295,12.5312,12.5345,12.5354,12.5379,12.5412,12.5454,12.5471,12.5484,12.5484,12.5467,12.5467,12.5434,12.5417,12.5354,12.5346,12.5303,12.5238,12.5188,12.5171,12.5154,12.5104,12.5096,12.5054,12.5046,12.5021,12.4987,12.4971,12.4887,12.4846,12.4837,12.4813,12.4804,12.4737,12.4729,12.4654,12.4646,12.4571,12.4562,12.452,12.4513,12.4487,12.4479,12.4454,12.443,12.4421,12.4404,12.4371,12.4362,12.4296,12.4287,12.4263,12.4246,12.4229,12.422,12.4163,12.4088,12.4062,12.405,12.405,12.4021,12.4004,12.3979,12.3967,12.3967,12.3929,12.3909,12.3871,12.3854,12.382,12.3809,12.3808,12.38,12.3712,12.3663,12.3562,12.3458,12.3395,12.3391,12.3388,12.3354,12.3352,12.3595,12.3595,12.3612,12.3616,12.3654,12.3659,12.3675,12.3679,12.3662,12.3638,12.3608,12.36,12.3595,12.3595,12.3582,12.3571,12.356,12.3543,12.3515,12.3513,12.3515,12.3518,12.3529,12.3532,12.354,12.3552,12.3568,12.3577,12.3579,12.3571,12.3543,12.3527,12.3507,12.3488,12.3471,12.3474,12.349,12.3507,12.352,12.3524,12.354,12.3549,12.3565,12.3582,12.3599,12.3607,12.3629,12.364,12.3649,12.3654,12.3649,12.3638,12.3621,12.3602,12.3582,12.3571,12.3577,12.3599,12.364,12.3688,12.3743,12.3827,12.3915,12.3999,12.409,12.4174,12.4265,12.4349,12.4438,12.4487,12.4521,12.4604,12.4574,12.4513,12.4452,12.4468,12.4496,12.4518,12.4549,12.4579,12.4593,12.4618,12.464,12.4663,12.4679,12.4688,12.4699,12.471,12.4718,12.4721,12.4724,12.4721,12.4715,12.4713,12.4702,12.469,12.4677,12.4668,12.4649,12.4629,12.4618,12.4602,12.4574,12.4549,12.4524,12.4513,12.4511,12.4493,12.4474,12.4457,12.4438,12.4418,12.4393,12.4374,12.4357,12.4365,12.4393,12.4421,12.4449,12.4477,12.4479,12.4468,12.4454,12.4435,12.4421,12.4402,12.4388,12.4377,12.4357,12.4343,12.4324,12.431,12.4296,12.4304,12.4388,12.4463,12.454,12.4624,12.4699,12.4782,12.486,12.4935,12.4996,12.4996,12.5018,12.5093,12.5145,12.5177,12.5254,12.5329,12.5413,12.5488,12.551,12.5565,12.5674,12.5827,12.5977,12.5992,12.6121,12.6154,12.6202,12.6252,12.6293,12.6407,12.6552,12.6696,12.6843,12.6838,12.6838,12.6838,12.6832,12.6832,12.6832,12.6832,12.6832,12.6832,12.6827,12.6827,12.6827,12.6827,12.6827,12.6826,12.6824,12.6824,12.6824,12.6824,12.6824,12.6824,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.681,12.681,12.681,12.681,12.681,12.681,12.681,12.681,12.6802,12.6802,12.6802,12.6802,12.6802,12.6802,12.6802,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6782,12.6777,12.6777,12.6774,12.6774,12.6774,12.6774,12.6774,12.6774,12.6765,12.6765,12.6765,12.6765,12.6763,12.6763,12.6763,12.6763,12.6758,12.6757,12.6757,12.6763,12.6763,12.6763,12.6763,12.6771,12.676,12.6765,12.6765,12.6765,12.6763,12.6763,12.6763,12.6768,12.6768,12.6768,12.6768,12.6765,12.6765,12.6765,12.6771,12.677,12.6765,12.676,12.6752,12.6743,12.6736,12.6735,12.6727,12.6721,12.6714,12.6713,12.6704,12.6696,12.669,12.6688,12.6682,12.6674,12.6665,12.6657,12.6651,12.6649,12.6615,12.6549,12.6524,12.6496,12.6445,12.6427,12.6393,12.6407,12.6493,12.6499,12.6507,12.6529,12.6546,12.6521,12.6493,12.6515,12.6538,12.6568,12.659,12.6585,12.6538,12.649,12.6457,12.6452,12.6454,12.6435,12.6424,12.6446,12.6482,12.6515,12.6518,12.6502,12.6465,12.644,12.6399,12.6368,12.634,12.6343,12.6385,12.6413,12.6429,12.6432,12.6435,12.6435,12.6402,12.6363,12.6396,12.6432,12.6474,12.6435,12.6427,12.6402,12.6382,12.6371,12.6365,12.6346,12.6335,12.6329,12.6338,12.6354,12.6363,12.6385,12.6393,12.6354,12.6321,12.6279,12.624,12.6177,12.6115,12.6054,12.6007,12.5952,12.5896,12.584,12.5785,12.5743,12.5696,12.5652,12.5604,12.5554,12.5493,12.5443,12.5407,12.5371,12.5338,12.5296,12.5282,12.5246,12.5199,12.5143,12.5077,12.5021,12.4997,12.4974,12.4927,12.4885,12.4846,12.4818,12.4793,12.4768,12.4727,12.4729,12.4743,12.4754,12.4763,12.4765,12.4761,12.476,12.4746,12.4729,12.4702,12.4682,12.4663,12.4679,12.4702,12.4732,12.476,12.4807,12.4871,12.4927,12.4988,12.4996,12.5043,12.5107,12.5163,12.5218,12.5271,12.5315,12.5352,12.5379,12.5407,12.5429,12.5452,12.5463,12.5463,12.5465,12.5446,12.5427,12.5393,12.5354,12.5313,12.5274,12.524,12.5213,12.5188,12.5168,12.5157,12.5138,12.5118,12.5107,12.5102,12.5068,12.5018,12.4997,12.4996,12.4985,12.4963,12.4927,12.4871,12.4835,12.4796,12.4782,12.4771,12.4738,12.4699,12.4671,12.466,12.4646,12.4649,12.4652,12.4646,12.4635,12.4629,12.4618,12.4604,12.4543,12.4488,12.4432,12.4363,12.4343,12.4346,12.4354,12.4357,12.436,12.4363,12.4365,12.4374,12.4377,12.4377,12.4371,12.4374,12.4371,12.4371,12.4374,12.4368,12.4371,12.4379,12.4382,12.4371,12.4338,12.4304,12.4263,12.4224,12.4168,12.4121,12.4074,12.4027,12.3985,12.396,12.3954,12.394,12.391,12.3854,12.3804,12.3749,12.3715,12.3677,12.3632,12.3593,12.3588,12.3618,12.3646,12.3677,12.3699,12.3721,12.3735,12.3757,12.3779,12.381,12.3829,12.3852,12.3874,12.3896,12.3932,12.3949,12.3955,12.3963,12.3965,12.396,12.3949,12.3935,12.3915,12.3904,12.3893,12.3879,12.3868,12.3849,12.3824,12.379,12.3771,12.3743,12.374,12.3735,12.3721,12.371,12.3682,12.3649,12.3615,12.3577,12.3529,12.3482,12.3435,12.3385,12.3338,12.3299,12.3265,12.3238,12.3213,12.3171,12.3124,12.3077,12.3035,12.3013,12.3024,12.3046,12.3068,12.309,12.3118,12.314,12.3163,12.3185,12.3199,12.3216,12.3221,12.3243,12.3265,12.3282,12.3296,12.3313,12.3327,12.3349,12.3379,12.3402,12.3424,12.3432,12.3435,12.3435,12.3438,12.344,12.3449,12.3443,12.3452,12.3468,12.3482,12.3499,12.3518,12.3535,12.3557,12.3579,12.3602,12.3629,12.366,12.3688,12.3724,12.3752,12.3788,12.3824,12.3865,12.3902,12.3943,12.3977,12.4021,12.4049,12.4085,12.4099,12.4102,12.411,12.4113,12.4115,12.4115,12.411,12.4099,12.4079,12.4054,12.4027,12.4002,12.3988,12.3982,12.3979,12.3952,12.3932,12.3913,12.391,12.3902,12.3904,12.3899,12.3907,12.391,12.3924,12.3952,12.3982,12.401,12.4032,12.406,12.409,12.4113,12.414,12.4163,12.419,12.4221,12.4249,12.4285,12.4271,12.4224,12.4177,12.4129,12.4088,12.4068,12.4057,12.4043,12.4032,12.4013,12.3985,12.396,12.3932,12.3907,12.3888,12.386,12.3871,12.3885,12.3902,12.3915,12.3932,12.3932,12.3935,12.3929,12.3918,12.3899,12.3865,12.3846,12.3818,12.3799,12.3793,12.3796,12.3799,12.3802,12.3815,12.3832,12.386,12.3888,12.391,12.3918,12.3907,12.3902,12.3924,12.3952,12.3979,12.4002,12.4032,12.406,12.4088,12.4118,12.4146,12.4174,12.4204,12.4238,12.426,12.4263,12.4235,12.421,12.4182,12.419,12.4199,12.4221,12.4238,12.4252,12.426,12.4268,12.4277,12.4288,12.4296,12.4296,12.4304,12.4307,12.4315,12.4324,12.434,12.4343,12.4349,12.4351,12.4352,12.436,12.4363,12.4357,12.4352,12.4346,12.4335,12.4321,12.431,12.4296,12.4285,12.4274,12.426,12.424,12.4215,12.4179,12.4146,12.4121,12.4102,12.409,12.4071,12.4057,12.4038,12.4046,12.4074,12.4102,12.4163,12.4227,12.4288,12.4357,12.4418,12.4474,12.4515,12.4552,12.4593,12.4629,12.4752,12.4763,12.4815,12.4856,12.4877,12.4938,12.4985,12.4996,12.4997,12.5013,12.501,12.4997,12.4949,12.4921,12.4888,12.4843,12.4829,12.4877,12.494,12.4996,12.5007,12.5054,12.509,12.5138,12.5199,12.5254,12.5293,12.5321,12.5346,12.5318,12.5263,12.5207,12.5154,12.5118,12.5124,12.5157,12.5199,12.5238,12.5293,12.5354,12.5402,12.5429,12.5435,12.546,12.5502,12.5538,12.5557,12.5565,12.5582,12.5604,12.5627,12.5654,12.569,12.5752,12.5813,12.5863,12.5902,12.5943,12.599,12.6029,12.6093,12.614,12.619,12.6232,12.6237,12.6288,12.6329,12.6377,12.644,12.6502,12.6557,12.6596,12.6615,12.6615,12.6618,12.6646,12.6674,12.6715,12.6757,12.6804,12.686,12.6896,12.691,12.6913,12.6949,12.701,12.7071,12.7135,12.719,12.7204,12.7193,12.7188,12.7188,12.7204,12.7213,12.7235,12.7263,12.7318,12.7352,12.7407,12.7468,12.7524,12.7557,12.7582,12.761,12.7665,12.7688,12.7721,12.7771,12.7827,12.7893,12.7949,12.8013,12.8079,12.8143,12.8204,12.826,12.8315,12.8313,12.829,12.8282,12.831,12.8363,12.8418,12.8446,12.8488,12.8538,12.8565,12.8613,12.8668,12.8724,12.8771,12.8802,12.8843,12.8885,12.8932,12.9002,12.9049,12.9082,12.9118,12.9157,12.9221,12.9249,12.9257,12.9238,12.921,12.9199,12.9207,12.9229,12.9257,12.9285,12.9349,12.941,12.9471,12.9532,12.9596,12.9657,12.9718,12.9774,12.9829,12.9877,12.991,12.9929,12.9879,12.9818,12.9749,12.9693,12.9646,12.9602,12.9582,12.9607,12.9627,12.9613,12.959,12.9624,12.9671,12.9715,12.9768,12.9832,12.9871,12.9927,12.9977,12.9996,12.9997,13.0003,13.0004,13.0032,13.0054,13.0096,13.0143,13.0163,13.0177,13.0215,13.0271,13.0332,13.0352,13.0388,13.0402,13.0413,13.0479,13.0527,13.0582,13.0624,13.066,13.0702,13.0757,13.0799,13.0824,13.0829,13.0829,13.0849,13.086,13.0879,13.0913,13.0954,13.0988,13.1013,13.1054,13.1088,13.1113,13.1113,13.1113,13.1082,13.1096,13.1138,13.1199,13.1252,13.1302,13.134,13.1374,13.1421,13.1477,13.1527,13.156,13.161,13.1665,13.1699,13.174,13.181,13.1852,13.1877,13.1882,13.1888,13.1913,13.194,13.1982,13.2035,13.209,13.2152,13.2207,13.2257,13.2299,13.234,13.2382,13.2424,13.2468,13.2502,13.2557,13.2604,13.2646,13.2679,13.2727,13.2777,13.2824,13.2863,13.2882,13.2896,13.2929,13.2968,13.301,13.3043,13.3071,13.3104,13.3165,13.3227,13.3282,13.3321,13.3363,13.3418,13.3474,13.3521,13.359,13.3643,13.3657,13.3615,13.3579,13.3565,13.3577,13.359,13.3615,13.3663,13.3718,13.3782,13.3849,13.3899,13.3924,13.3929,13.3907,13.3871,13.3838,13.384,13.3846,13.3832,13.3824,13.3829,13.3849,13.3882,13.3929,13.3977,13.4024,13.4057,13.4099,13.4118,13.4124,13.4113,13.4113,13.4104,13.4068,13.4013,13.3952,13.3902,13.3852,13.381,13.3768,13.3727,13.3677,13.3624,13.356,13.3499,13.3438,13.3407,13.3393,13.3399,13.3377,13.334,13.3299,13.3249,13.3202,13.3157,13.3124,13.3093,13.3065,13.3057,13.3077,13.311,13.3152,13.3204,13.3252,13.3307,13.3377,13.3438,13.3493,13.3527,13.354,13.3543,13.3557,13.3577,13.3602,13.3643,13.3682,13.3724,13.3752,13.3749,13.3746,13.3788,13.3822,13.3857,13.3904,13.3954,13.4004,13.4052,13.4113,13.4177,13.4224,13.426,13.4296,13.4338,13.4393,13.4449,13.4504,13.4571,13.4627,13.4674,13.4715,13.474,13.4746,13.4743,13.4757,13.4782,13.4815,13.4871,13.4927,13.4979,13.4996,13.5021,13.5049,13.5065,13.5102,13.5135,13.516,13.5202,13.5235,13.5254,13.5252,13.5249,13.5268,13.5296,13.5335,13.5371,13.5377,13.5418,13.5457,13.5499,13.5532,13.556,13.5599,13.564,13.5679,13.5683,13.5721,13.5774,13.5829,13.5885,13.5927,13.5965,13.5993,13.6018,13.6046,13.6071,13.6104,13.616,13.6213,13.6268,13.6327,13.6379,13.6427,13.6468,13.6502,13.6535,13.6568,13.6602,13.6643,13.669,13.6738,13.6785,13.6849,13.6904,13.6952,13.6985,13.701,13.7029,13.7043,13.7046,13.7077,13.7077,13.7079,13.709,13.7113,13.7143,13.7165,13.7188,13.7215,13.7238,13.726,13.7282,13.731,13.7338,13.736,13.739,13.7413,13.744,13.7463,13.749,13.7521,13.7532,13.7549,13.7577,13.7613,13.764,13.7677,13.7713,13.774,13.7777,13.7813,13.7849,13.789,13.7932,13.7979,13.8029,13.8079,13.8132,13.8188,13.8243,13.8299,13.8363,13.8418,13.8474,13.8529,13.8585,13.864,13.8702,13.8771,13.8832,13.8902,13.8971,13.904,13.9088,13.9127,13.916,13.9188,13.9221,13.9254,13.9288,13.9313,13.9349,13.9374,13.9399,13.9432,13.946,13.9493,13.9527,13.9552,13.9568,13.9602,13.964,13.9696,13.9752,13.9821,13.9888,13.9957,13.9996,14.0018,14.0068,14.0121,14.0185,14.024,14.0282,14.0324,14.0374,14.044,14.0504,14.0565,14.0613,14.0663,14.0679,14.069,14.0727,14.0763,14.0799,14.0846,14.091,14.0971,14.1018,14.1065,14.1113,14.1177,14.1224,14.1274,14.1324,14.1365,14.1413,14.1477,14.1538,14.1607,14.1668,14.1724,14.1771,14.1818,14.1865,14.1921,14.1977,14.2029,14.2071,14.211,14.2122,14.216,14.2207,14.2254,14.2302,14.2357,14.2404,14.246,14.2513,14.256,14.2615,14.2663,14.2704,14.2743,14.2785,14.2818,14.281,14.2821,14.2849,14.2846,14.2843,14.2857,14.2871,14.2888,14.2915,14.294,14.2968,14.3002,14.304,14.3096,14.3138,14.3177,14.321,14.3243,14.3279,14.3279,14.3268,14.3296,14.336,14.3413,14.3457,14.3485,14.3535,14.3588,14.3629,14.3663,14.3682,14.3671,14.3657,14.3663,14.3704,14.3752,14.379,14.3802,14.3788,14.376,14.3738,14.3763,14.3796,14.3838,14.3871,14.3877,14.3868,14.3888,14.3913,14.3938,14.3965,14.3949,14.3996,14.4004,14.4015,14.4027,14.406,14.4079,14.4143,14.4171,14.421,14.426,14.4307,14.4363,14.4427,14.4479,14.4521,14.4555,14.4568,14.4624,14.4693,14.476,14.4802,14.4829,14.486,14.4902,14.4957,14.4996,14.4997,14.5018,14.5079,14.5124,14.5165,14.5213,14.5268,14.5332,14.5393,14.5449,14.5496,14.5513,14.5507,14.5507,14.5524,14.5552,14.5602,14.5657,14.5713,14.5774,14.5829,14.5877,14.5913,14.5949,14.601,14.6071,14.6115,14.6149,14.6171,14.6182,14.621,14.6246,14.6288,14.6335,14.639,14.646,14.6515,14.6554,14.6582,14.6588,14.6585,14.6596,14.6624,14.6671,14.6732,14.6763,14.6777,14.6813,14.6854,14.6904,14.6965,14.6957,14.6921,14.6921,14.696,14.7021,14.704,14.704,14.7024,14.7029,14.7057,14.7104,14.7152,14.7213,14.7268,14.7324,14.7365,14.7385,14.7402,14.7415,14.7435,14.7468,14.7502,14.7563,14.7582,14.7605,14.7624,14.764,14.7652,14.7677,14.7696,14.7724,14.7749,14.7768,14.7802,14.7829,14.7854,14.7888,14.7921,14.7949,14.7982,14.8015,14.8049,14.8082,14.8121,14.8154,14.8196,14.8238,14.8271,14.8296,14.8315,14.8343,14.836,14.8374,14.8385,14.8385,14.8396]}]],[[{"lng":[-16.5334,-16.5351,-16.5355,-16.5346,-16.5321,-16.5313,-16.5305,-16.5297,-16.5296,-16.5288,-16.5288,-16.5255,-16.5255,-16.5246,-16.5247,-16.523,-16.523,-16.5221,-16.5213,-16.5204,-16.52,-16.518,-16.518,-16.518,-16.5188,-16.5188,-16.5196,-16.5196,-16.5188,-16.5197,-16.5197,-16.523,-16.523,-16.5238,-16.5238,-16.5246,-16.5246,-16.5263,-16.5271,-16.528,-16.528,-16.5288,-16.5296,-16.5299,-16.5305,-16.5313,-16.5334],"lat":[15.7959,15.7959,15.7962,15.8038,15.8096,15.8105,15.817,15.8179,15.8205,15.8213,15.8246,15.8329,15.8371,15.8379,15.8405,15.8429,15.8463,15.8471,15.8587,15.8596,15.8617,15.8604,15.858,15.8504,15.8496,15.8479,15.8471,15.8404,15.8388,15.8379,15.8362,15.8288,15.8254,15.8246,15.8229,15.8221,15.8188,15.8171,15.8129,15.8121,15.8087,15.8079,15.8037,15.8035,15.8029,15.7987,15.7959]},{"lng":[-16.5163,-16.5155,-16.5155,-16.5168,-16.5175,-16.5213,-16.5213,-16.5192,-16.5163],"lat":[15.8662,15.8679,15.8729,15.875,15.875,15.8712,15.8688,15.8659,15.8662]},{"lng":[-16.4872,-16.4871,-16.4879,-16.4883,-16.4896,-16.4896,-16.4897,-16.4884,-16.4872],"lat":[15.9162,15.9213,15.9232,15.9242,15.9237,15.9219,15.9171,15.9158,15.9162]},{"lng":[-16.4959,-16.4946,-16.4946,-16.4946,-16.4946,-16.4955,-16.4946,-16.4955,-16.4955,-16.4955,-16.4955,-16.4963,-16.4963,-16.4946,-16.4951,-16.4967,-16.4975,-16.4996,-16.4996,-16.4971,-16.4988,-16.498,-16.4988,-16.498,-16.4997,-16.498,-16.498,-16.4988,-16.4988,-16.4975,-16.4959],"lat":[15.9359,15.9371,15.9401,15.9437,15.9454,15.9463,15.9488,15.9496,15.9538,15.9552,15.9596,15.9604,15.9638,15.9654,15.9658,15.9658,15.965,15.9629,15.9604,15.9579,15.9554,15.9538,15.9496,15.9479,15.9454,15.9438,15.9396,15.9388,15.9371,15.9359,15.9359]},{"lng":[-16.4988,-16.4988,-16.4988,-16.4988,-16.4963,-16.4964,-16.4976,-16.5001,-16.5021,-16.5022,-16.503,-16.503,-16.503,-16.5021,-16.5021,-16.5009,-16.4988],"lat":[15.9896,15.9911,15.9944,15.9963,15.9996,16.0062,16.0083,16.0083,16.0062,16.002,16.0012,15.9996,15.9979,15.9971,15.9904,15.9884,15.9896]},{"lng":[-16.503,-16.503,-16.503,-16.5022,-16.5021,-16.5005,-16.5013,-16.5035,-16.5064,-16.5064,-16.5064,-16.5063,-16.5043,-16.503],"lat":[16.0204,16.0236,16.0254,16.0263,16.0304,16.0346,16.0387,16.0392,16.0345,16.0233,16.0216,16.0196,16.0192,16.0204]},{"lng":[-16.4772,-16.4772,-16.4755,-16.4755,-16.4738,-16.4739,-16.473,-16.473,-16.4743,-16.4768,-16.4785,-16.4797,-16.4806,-16.483,-16.483,-16.4839,-16.4839,-16.4851,-16.488,-16.4889,-16.4876,-16.4868,-16.4793,-16.4772],"lat":[16.0421,16.0454,16.0471,16.0479,16.0513,16.0571,16.0579,16.0596,16.0609,16.0608,16.06,16.0587,16.0554,16.0529,16.0521,16.0512,16.0479,16.0475,16.0438,16.0388,16.0383,16.0392,16.0392,16.0421]},{"lng":[-16.4968,-16.4975,-16.498,-16.4982,-16.4986,-16.4989,-16.4992,-16.5006,-16.5034,-16.5091,-16.509,-16.509,-16.5092,-16.5098,-16.5107,-16.5105,-16.5113,-16.5113,-16.513,-16.513,-16.5138,-16.5138,-16.5138,-16.5147,-16.5146,-16.5155,-16.5155,-16.5163,-16.5163,-16.5171,-16.5171,-16.518,-16.518,-16.518,-16.5188,-16.5188,-16.5196,-16.5196,-16.5205,-16.5213,-16.52,-16.5188,-16.5171,-16.5163,-16.5163,-16.5155,-16.5155,-16.5138,-16.5138,-16.5146,-16.5146,-16.5155,-16.5146,-16.5138,-16.5138,-16.513,-16.5122,-16.5096,-16.5096,-16.5088,-16.5088,-16.5097,-16.5096,-16.5105,-16.5096,-16.5096,-16.5088,-16.5088,-16.5088,-16.508,-16.5096,-16.5084,-16.5071,-16.5088,-16.5088,-16.5067,-16.5055,-16.5055,-16.508,-16.508,-16.5072,-16.5072,-16.508,-16.508,-16.508,-16.508,-16.508,-16.5072,-16.5072,-16.5056,-16.5055,-16.5064,-16.5064,-16.5047,-16.5047,-16.5047,-16.5047,-16.5043,-16.5031,-16.503,-16.5022,-16.5021,-16.5014,-16.5014,-16.4993,-16.4972,-16.4969,-16.4968,-16.4955,-16.4955,-16.4943,-16.4914,-16.4903,-16.4899,-16.4892,-16.4888,-16.4897,-16.4901,-16.4901,-16.4901,-16.49,-16.4898,-16.4891,-16.489,-16.4891,-16.4891,-16.4891,-16.4893,-16.4896,-16.4896,-16.4895,-16.4895,-16.4896,-16.4897,-16.4899,-16.4901,-16.4902,-16.4904,-16.4905,-16.4908,-16.491,-16.4912,-16.4914,-16.4917,-16.4917,-16.4917,-16.4917,-16.4919,-16.4921,-16.4923,-16.4928,-16.4931,-16.4933,-16.4937,-16.4938,-16.4938,-16.4938,-16.4937,-16.4938,-16.4938,-16.494,-16.4941,-16.4946,-16.4949,-16.4958,-16.4961,-16.4964,-16.4967,-16.4968],"lat":[16.0732,16.0736,16.0737,16.0736,16.0732,16.0727,16.0717,16.0696,16.0664,16.0598,16.0507,16.0415,16.0354,16.0211,16.0122,16.0121,15.9996,15.9988,15.9946,15.9862,15.9854,15.9567,15.9546,15.9538,15.9488,15.948,15.9321,15.9312,15.9237,15.9229,15.9187,15.9179,15.902,15.9004,15.8996,15.8954,15.8946,15.8904,15.8896,15.8779,15.8775,15.8788,15.8879,15.8888,15.8912,15.892,15.8938,15.8971,15.9038,15.9046,15.9088,15.9096,15.9171,15.9179,15.9221,15.9229,15.9287,15.9329,15.9354,15.9362,15.9413,15.9421,15.9463,15.9471,15.9479,15.9546,15.9554,15.9596,15.9729,15.9737,15.9813,15.9816,15.9838,15.9854,15.9871,15.9875,15.9887,15.993,15.9954,16.0013,16.0021,16.0079,16.0087,16.0137,16.0204,16.0223,16.0354,16.0363,16.0379,16.0404,16.0445,16.0454,16.0471,16.0504,16.0527,16.0596,16.0646,16.065,16.0637,16.0595,16.0587,16.0538,16.0529,16.0512,16.0483,16.0496,16.0574,16.0592,16.0571,16.0495,16.0492,16.052,16.0543,16.0555,16.0567,16.0572,16.0577,16.058,16.0581,16.0593,16.0596,16.0599,16.0605,16.0607,16.0613,16.0616,16.0619,16.0622,16.0632,16.0635,16.0639,16.0643,16.0646,16.0648,16.0662,16.0664,16.067,16.067,16.0671,16.0677,16.068,16.068,16.0683,16.0684,16.0686,16.0689,16.0693,16.0693,16.0692,16.0694,16.0697,16.07,16.0703,16.0706,16.0707,16.0708,16.0712,16.0715,16.0718,16.0721,16.0722,16.0723,16.0724,16.0725,16.0728,16.073,16.073,16.073,16.0732]},{"lng":[-14.9459,-14.9503,-14.9528,-14.9545,-14.9556,-14.9564,-14.9573,-14.9575,-14.9584,-14.9609,-14.9631,-14.9689,-14.9739,-14.9803,-14.9861,-14.9903,-14.9945,-14.9981,-14.9992,-15,-15.0025,-15.0031,-15.005,-15.007,-15.0089,-15.0117,-15.0145,-15.017,-15.0206,-15.0248,-15.0306,-15.0353,-15.0411,-15.0475,-15.0534,-15.057,-15.0606,-15.0639,-15.0661,-15.0686,-15.0709,-15.0734,-15.077,-15.0809,-15.0845,-15.0903,-15.097,-15.1017,-15.105,-15.1092,-15.1128,-15.1153,-15.1159,-15.1142,-15.1125,-15.1095,-15.105,-15.1,-15.095,-15.0898,-15.0873,-15.0864,-15.0875,-15.0903,-15.0936,-15.097,-15.1006,-15.1039,-15.1073,-15.1125,-15.1189,-15.1253,-15.1317,-15.1375,-15.1434,-15.1506,-15.157,-15.1636,-15.1684,-15.1734,-15.1798,-15.1848,-15.1903,-15.1959,-15.2017,-15.2075,-15.2123,-15.2173,-15.222,-15.2273,-15.232,-15.2384,-15.245,-15.2506,-15.2559,-15.2617,-15.2667,-15.272,-15.2778,-15.2828,-15.2878,-15.2936,-15.3003,-15.3067,-15.3095,-15.3159,-15.3217,-15.3267,-15.3314,-15.337,-15.342,-15.3478,-15.3542,-15.36,-15.3642,-15.3692,-15.3748,-15.3798,-15.3848,-15.3889,-15.3945,-15.3995,-15.4059,-15.4125,-15.417,-15.4206,-15.4242,-15.4281,-15.4303,-15.4328,-15.4345,-15.4375,-15.4406,-15.445,-15.45,-15.4564,-15.462,-15.4681,-15.4728,-15.4761,-15.4778,-15.4825,-15.4884,-15.4942,-15.4989,-15.4992,-15.5039,-15.5073,-15.5092,-15.5098,-15.5103,-15.5114,-15.5125,-15.5142,-15.5178,-15.5223,-15.5264,-15.532,-15.5367,-15.5425,-15.5489,-15.5561,-15.5628,-15.5689,-15.5756,-15.5814,-15.587,-15.5923,-15.5973,-15.6025,-15.6075,-15.6142,-15.6206,-15.6239,-15.6253,-15.6264,-15.6275,-15.6289,-15.6289,-15.6295,-15.6336,-15.6378,-15.6428,-15.6484,-15.6534,-15.6592,-15.6639,-15.6695,-15.6759,-15.6825,-15.6898,-15.6961,-15.7017,-15.7078,-15.7136,-15.7186,-15.7245,-15.7286,-15.7345,-15.7411,-15.7475,-15.7534,-15.7603,-15.7664,-15.7723,-15.7781,-15.7836,-15.7886,-15.7945,-15.7998,-15.8006,-15.8019,-15.8022,-15.8023,-15.8048,-15.8092,-15.8156,-15.8223,-15.8255,-15.8258,-15.827,-15.832,-15.837,-15.842,-15.8489,-15.855,-15.8592,-15.8614,-15.8616,-15.8617,-15.8636,-15.8675,-15.872,-15.8764,-15.8814,-15.8873,-15.8936,-15.8995,-15.9067,-15.9123,-15.9175,-15.9231,-15.9278,-15.9328,-15.9378,-15.9382,-15.9383,-15.9425,-15.9467,-15.9509,-15.9559,-15.96,-15.9673,-15.9731,-15.9789,-15.9848,-15.9892,-15.9948,-15.9992,-16.0006,-16.0061,-16.0128,-16.0178,-16.0225,-16.0275,-16.0331,-16.0375,-16.0417,-16.0484,-16.0481,-16.0517,-16.0561,-16.0598,-16.0636,-16.0667,-16.0689,-16.0711,-16.0725,-16.0726,-16.0756,-16.0789,-16.0817,-16.0848,-16.0889,-16.0931,-16.0967,-16.0989,-16.0986,-16.0995,-16.102,-16.105,-16.1089,-16.1125,-16.1175,-16.1234,-16.1306,-16.1375,-16.145,-16.1506,-16.1548,-16.1589,-16.1623,-16.1656,-16.1698,-16.1734,-16.1767,-16.1773,-16.1853,-16.1957,-16.208,-16.2223,-16.2406,-16.249,-16.2585,-16.2693,-16.2717,-16.2741,-16.2796,-16.2808,-16.2836,-16.2888,-16.2964,-16.3011,-16.3043,-16.3051,-16.3047,-16.3079,-16.3127,-16.3175,-16.3218,-16.3262,-16.329,-16.3318,-16.3318,-16.3306,-16.3318,-16.3346,-16.3378,-16.3534,-16.3565,-16.3604,-16.3632,-16.3616,-16.3612,-16.3577,-16.3517,-16.3497,-16.3493,-16.3513,-16.3608,-16.3676,-16.3716,-16.3748,-16.3807,-16.3887,-16.3903,-16.3935,-16.4006,-16.4094,-16.4146,-16.4209,-16.4313,-16.4377,-16.4417,-16.4424,-16.4464,-16.4484,-16.4488,-16.448,-16.4484,-16.4496,-16.45,-16.4484,-16.448,-16.4508,-16.4512,-16.4512,-16.448,-16.4444,-16.4436,-16.444,-16.4448,-16.4468,-16.4504,-16.4576,-16.4639,-16.4666,-16.4666,-16.468,-16.4692,-16.468,-16.467,-16.468,-16.468,-16.4651,-16.4635,-16.4618,-16.4593,-16.4585,-16.456,-16.4543,-16.4531,-16.4535,-16.4551,-16.456,-16.4601,-16.4634,-16.4651,-16.466,-16.4676,-16.4697,-16.4697,-16.4705,-16.4705,-16.4722,-16.4755,-16.4764,-16.4739,-16.4739,-16.473,-16.473,-16.4743,-16.4763,-16.4772,-16.481,-16.4868,-16.4876,-16.4901,-16.4918,-16.493,-16.4935,-16.4976,-16.5005,-16.5005,-16.5013,-16.5013,-16.4963,-16.4947,-16.4946,-16.4955,-16.498,-16.4972,-16.498,-16.498,-16.4963,-16.4963,-16.498,-16.498,-16.4967,-16.4958,-16.4921,-16.4913,-16.4913,-16.4892,-16.4859,-16.4834,-16.4817,-16.4796,-16.4788,-16.4788,-16.4775,-16.4771,-16.4788,-16.4788,-16.4775,-16.4742,-16.4734,-16.4717,-16.4693,-16.468,-16.4671,-16.4684,-16.47,-16.4709,-16.4742,-16.4759,-16.4775,-16.4809,-16.4855,-16.4855,-16.483,-16.483,-16.4838,-16.4838,-16.4867,-16.4892,-16.4909,-16.4926,-16.4955,-16.4946,-16.4946,-16.4938,-16.4938,-16.493,-16.4938,-16.4938,-16.4921,-16.4938,-16.4922,-16.4921,-16.4913,-16.4921,-16.4913,-16.4913,-16.4905,-16.4905,-16.4892,-16.4851,-16.4842,-16.4838,-16.4846,-16.4846,-16.4859,-16.4871,-16.4838,-16.4838,-16.4855,-16.4855,-16.4859,-16.4871,-16.4867,-16.4865,-16.4863,-16.4876,-16.4888,-16.4896,-16.4905,-16.4905,-16.4918,-16.4919,-16.4921,-16.4904,-16.4905,-16.4913,-16.4913,-16.4951,-16.4988,-16.4988,-16.5005,-16.5005,-16.5005,-16.4996,-16.4996,-16.5013,-16.5005,-16.5005,-16.4988,-16.4988,-16.4996,-16.4988,-16.5,-16.5026,-16.5038,-16.503,-16.503,-16.5046,-16.5038,-16.5047,-16.503,-16.5046,-16.5038,-16.503,-16.5038,-16.5038,-16.503,-16.5038,-16.5038,-16.5046,-16.5046,-16.5055,-16.5055,-16.5046,-16.5046,-16.5071,-16.5071,-16.5063,-16.5063,-16.5071,-16.5071,-16.5088,-16.5088,-16.508,-16.508,-16.5088,-16.5088,-16.5097,-16.5105,-16.5113,-16.5121,-16.513,-16.5129,-16.5138,-16.5146,-16.5155,-16.5155,-16.5142,-16.5138,-16.5128,-16.5113,-16.5113,-16.5096,-16.5096,-16.508,-16.5071,-16.5071,-16.508,-16.5088,-16.5104,-16.5105,-16.513,-16.513,-16.5155,-16.5163,-16.5171,-16.5171,-16.518,-16.5196,-16.5205,-16.5213,-16.5238,-16.5238,-16.5247,-16.5255,-16.528,-16.5297,-16.533,-16.533,-16.5355,-16.5363,-16.5371,-16.538,-16.5413,-16.543,-16.5438,-16.5471,-16.5471,-16.5488,-16.5488,-16.5505,-16.5505,-16.5538,-16.5546,-16.5563,-16.5563,-16.558,-16.5588,-16.5605,-16.5605,-16.5638,-16.5671,-16.5671,-16.5721,-16.573,-16.5755,-16.5763,-16.5788,-16.5813,-16.588,-16.5888,-16.5946,-16.5946,-16.5955,-16.6022,-16.6071,-16.6096,-16.613,-16.6138,-16.6155,-16.6163,-16.618,-16.6213,-16.6271,-16.6271,-16.633,-16.6355,-16.6438,-16.6463,-16.6513,-16.6513,-16.6538,-16.6546,-16.6571,-16.6571,-16.6605,-16.6613,-16.6655,-16.6655,-16.6721,-16.6721,-16.6755,-16.678,-16.6805,-16.6805,-16.6847,-16.6863,-16.6905,-16.6905,-16.6955,-16.6963,-16.6988,-16.6996,-16.7021,-16.7038,-16.7071,-16.7105,-16.7163,-16.7163,-16.7246,-16.7263,-16.7271,-16.7305,-16.7307,-16.7102,-16.7012,-16.6958,-16.6911,-16.6901,-16.6892,-16.6886,-16.6888,-16.6883,-16.6887,-16.6889,-16.6908,-16.6973,-16.6995,-16.7031,-16.7033,-16.7075,-16.7127,-16.7155,-16.7187,-16.7209,-16.7231,-16.7247,-16.7253,-16.7274,-16.7277,-16.7285,-16.7337,-16.7338,-16.7373,-16.738,-16.7401,-16.7401,-16.7423,-16.7432,-16.7415,-16.7405,-16.7381,-16.7366,-16.7333,-16.7319,-16.7305,-16.7233,-16.718,-16.7167,-16.7155,-16.7136,-16.7094,-16.7081,-16.7073,-16.7048,-16.7022,-16.701,-16.6961,-16.6922,-16.6907,-16.6867,-16.6844,-16.6828,-16.6811,-16.676,-16.6745,-16.673,-16.6701,-16.6687,-16.6669,-16.6635,-16.6626,-16.6608,-16.6557,-16.6536,-16.6476,-16.6428,-16.6423,-16.6412,-16.6384,-16.6365,-16.6338,-16.6312,-16.6307,-16.6306,-16.6285,-16.6245,-16.6207,-16.6165,-16.6133,-16.6121,-16.6121,-16.6123,-16.6129,-16.6134,-16.614,-16.6141,-16.6128,-16.6092,-16.6072,-16.6045,-16.6037,-16.6012,-16.5993,-16.5958,-16.5907,-16.5836,-16.5824,-16.581,-16.5722,-16.5696,-16.5671,-16.5647,-16.5596,-16.5571,-16.554,-16.5519,-16.5471,-16.5451,-16.5418,-16.54,-16.5381,-16.535,-16.5332,-16.5307,-16.5279,-16.5251,-16.5237,-16.5225,-16.5202,-16.5186,-16.5173,-16.5158,-16.5145,-16.5135,-16.5132,-16.5129,-16.5124,-16.5123,-16.5123,-16.5135,-16.5143,-16.5156,-16.5164,-16.5162,-16.5158,-16.5146,-16.5137,-16.5124,-16.5111,-16.5096,-16.5077,-16.506,-16.5044,-16.503,-16.502,-16.5007,-16.499,-16.4974,-16.4954,-16.4941,-16.4915,-16.4901,-16.4889,-16.4873,-16.4861,-16.4845,-16.4832,-16.4817,-16.4804,-16.4791,-16.4779,-16.4759,-16.4751,-16.4736,-16.4709,-16.4691,-16.4659,-16.4641,-16.4618,-16.4594,-16.4564,-16.4546,-16.4538,-16.4526,-16.4515,-16.4509,-16.4483,-16.4458,-16.4436,-16.4414,-16.4385,-16.4362,-16.4338,-16.4311,-16.428,-16.4255,-16.4228,-16.4189,-16.4153,-16.4107,-16.4039,-16.3979,-16.3932,-16.3878,-16.3839,-16.3796,-16.3773,-16.3752,-16.3732,-16.3719,-16.3693,-16.366,-16.363,-16.3587,-16.3547,-16.3515,-16.3496,-16.3468,-16.3448,-16.3421,-16.3411,-16.3393,-16.3371,-16.3356,-16.3346,-16.3344,-16.3354,-16.3371,-16.3388,-16.3416,-16.3446,-16.3454,-16.3457,-16.345,-16.3436,-16.3415,-16.338,-16.335,-16.3332,-16.3305,-16.3274,-16.3237,-16.3202,-16.3171,-16.3142,-16.3112,-16.3102,-16.3093,-16.3088,-16.3084,-16.308,-16.3064,-16.3045,-16.3017,-16.3,-16.2969,-16.295,-16.2931,-16.2908,-16.2885,-16.2868,-16.2858,-16.2849,-16.2832,-16.2806,-16.278,-16.276,-16.2742,-16.2723,-16.2708,-16.27,-16.2692,-16.2685,-16.268,-16.2681,-16.2685,-16.2689,-16.2692,-16.2697,-16.2699,-16.2693,-16.2688,-16.2683,-16.2672,-16.2661,-16.2647,-16.2632,-16.2616,-16.2601,-16.2587,-16.2572,-16.2551,-16.2526,-16.2499,-16.2475,-16.245,-16.2438,-16.2426,-16.2416,-16.2392,-16.2376,-16.2373,-16.2372,-16.2373,-16.2377,-16.2374,-16.2371,-16.2364,-16.2345,-16.2321,-16.2304,-16.227,-16.2249,-16.2188,-16.2136,-16.2116,-16.2086,-16.205,-16.2041,-16.2012,-16.1989,-16.1968,-16.1962,-16.1951,-16.1947,-16.1947,-16.1947,-16.1947,-16.1947,-16.1942,-16.1943,-16.1944,-16.1931,-16.1921,-16.1904,-16.1891,-16.1886,-16.1873,-16.1858,-16.1832,-16.1824,-16.1808,-16.18,-16.1796,-16.1789,-16.1784,-16.1782,-16.1783,-16.1787,-16.1782,-16.1781,-16.1779,-16.1776,-16.1771,-16.1767,-16.1764,-16.1759,-16.1755,-16.1761,-16.1766,-16.1771,-16.1771,-16.1772,-16.1771,-16.177,-16.177,-16.176,-16.1735,-16.1705,-16.1674,-16.1638,-16.1605,-16.1577,-16.1544,-16.1522,-16.1486,-16.147,-16.1422,-16.1395,-16.1368,-16.1345,-16.1317,-16.1284,-16.1235,-16.1213,-16.1177,-16.1158,-16.1114,-16.1095,-16.1074,-16.1049,-16.1025,-16.0992,-16.0971,-16.0948,-16.0938,-16.0921,-16.0905,-16.0879,-16.0863,-16.0843,-16.0816,-16.0794,-16.0759,-16.0728,-16.0703,-16.0681,-16.0658,-16.0602,-16.0577,-16.0545,-16.0521,-16.0495,-16.0469,-16.0439,-16.0408,-16.0377,-16.0361,-16.0348,-16.0321,-16.0264,-16.0176,-16.0156,-16.0139,-16.0099,-16.0078,-16.0058,-16.004,-16.0018,-16.0001,-15.999,-15.9973,-15.9968,-15.9962,-15.9954,-15.9944,-15.9936,-15.9929,-15.9922,-15.9916,-15.9911,-15.9904,-15.9899,-15.9891,-15.9884,-15.9873,-15.9854,-15.9811,-15.9777,-15.9744,-15.972,-15.9688,-15.963,-15.9575,-15.9532,-15.9479,-15.9434,-15.9369,-15.9331,-15.9276,-15.9244,-15.9224,-15.9196,-15.9171,-15.9163,-15.9157,-15.9147,-15.9133,-15.9103,-15.9075,-15.9048,-15.9026,-15.9009,-15.8992,-15.8976,-15.8962,-15.8942,-15.8915,-15.8837,-15.877,-15.8686,-15.8635,-15.8626,-15.8623,-15.8626,-15.8633,-15.8644,-15.8651,-15.8646,-15.862,-15.8588,-15.855,-15.8515,-15.8467,-15.842,-15.8385,-15.8342,-15.8305,-15.8264,-15.8234,-15.8201,-15.8161,-15.8104,-15.8034,-15.8028,-15.8023,-15.8018,-15.8002,-15.7989,-15.7962,-15.793,-15.7892,-15.7868,-15.7851,-15.7829,-15.7815,-15.7798,-15.7779,-15.7765,-15.775,-15.7737,-15.7729,-15.7727,-15.772,-15.7715,-15.7699,-15.7678,-15.7659,-15.7633,-15.7617,-15.7603,-15.7584,-15.7552,-15.7508,-15.7457,-15.7413,-15.7369,-15.7322,-15.7274,-15.7223,-15.7158,-15.7117,-15.7066,-15.7003,-15.6946,-15.687,-15.6791,-15.6722,-15.6664,-15.6628,-15.6554,-15.6498,-15.6435,-15.639,-15.6356,-15.6354,-15.6406,-15.6525,-15.6621,-15.6751,-15.6876,-15.692,-15.6903,-15.6841,-15.6773,-15.676,-15.6755,-15.6752,-15.6749,-15.6744,-15.6739,-15.6735,-15.673,-15.6726,-15.6711,-15.6689,-15.6666,-15.6649,-15.6624,-15.6608,-15.6587,-15.6565,-15.6544,-15.6509,-15.6484,-15.6457,-15.6434,-15.6408,-15.6379,-15.6351,-15.6322,-15.6307,-15.6272,-15.6251,-15.6218,-15.6185,-15.6148,-15.6122,-15.6099,-15.6078,-15.6064,-15.603,-15.6009,-15.5991,-15.5962,-15.5937,-15.5909,-15.5868,-15.5846,-15.5815,-15.5785,-15.5773,-15.5757,-15.5736,-15.5698,-15.5665,-15.563,-15.5597,-15.5562,-15.5523,-15.5496,-15.5452,-15.5432,-15.5408,-15.5386,-15.5364,-15.5348,-15.5313,-15.5276,-15.5238,-15.5204,-15.5175,-15.5148,-15.513,-15.5103,-15.5078,-15.5054,-15.5046,-15.503,-15.501,-15.499,-15.4973,-15.4953,-15.4923,-15.4904,-15.4884,-15.4864,-15.4844,-15.4817,-15.4786,-15.4759,-15.4731,-15.471,-15.4686,-15.4672,-15.4646,-15.4617,-15.4594,-15.4569,-15.4532,-15.4493,-15.4449,-15.4425,-15.4397,-15.4371,-15.4357,-15.4357,-15.4358,-15.4353,-15.4361,-15.436,-15.437,-15.4371,-15.4377,-15.4378,-15.4392,-15.4395,-15.4402,-15.4405,-15.4408,-15.4416,-15.4419,-15.4423,-15.443,-15.4435,-15.4446,-15.4454,-15.4458,-15.4185,-15.3912,-15.3639,-15.3366,-15.3093,-15.2807,-15.252,-15.2234,-15.1947,-15.1661,-15.1375,-15.1088,-15.0887,-15.0855,-15.0802,-15.0515,-15.0229,-15.0162,-14.9989,-14.9896,-14.9629,-14.9362,-14.9095,-14.8891,-14.8686,-14.8496,-14.836,-14.8203,-14.8111,-14.8036,-14.8028,-14.8027,-14.8002,-14.7981,-14.796,-14.7922,-14.7822,-14.7822,-14.7734,-14.7733,-14.7711,-14.7703,-14.7693,-14.7608,-14.7534,-14.7462,-14.7379,-14.7345,-14.7323,-14.7315,-14.7279,-14.7277,-14.7191,-14.7191,-14.7056,-14.7051,-14.6926,-14.6888,-14.683,-14.6642,-14.6634,-14.6628,-14.659,-14.6587,-14.6581,-14.6577,-14.6514,-14.6329,-14.6271,-14.6076,-14.6053,-14.5993,-14.5894,-14.5821,-14.5738,-14.5733,-14.5717,-14.5674,-14.5658,-14.5604,-14.5566,-14.5484,-14.5413,-14.5316,-14.5247,-14.519,-14.513,-14.5025,-14.4985,-14.4976,-14.494,-14.494,-14.4939,-14.4937,-14.4913,-14.4904,-14.4895,-14.4882,-14.4879,-14.4858,-14.4846,-14.4837,-14.4834,-14.4824,-14.4814,-14.4799,-14.4776,-14.4755,-14.4739,-14.4713,-14.4707,-14.4703,-14.4693,-14.4666,-14.456,-14.4472,-14.4356,-14.428,-14.4166,-14.408,-14.3998,-14.3799,-14.3683,-14.3681,-14.3515,-14.3513,-14.3302,-14.3301,-14.3301,-14.33,-14.3286,-14.3236,-14.2941,-14.2904,-14.2846,-14.2775,-14.2414,-14.2414,-14.2385,-14.2384,-14.1875,-14.1851,-14.1825,-14.1792,-14.1792,-14.1642,-14.1565,-14.1564,-14.1545,-14.1062,-14.0565,-14.0531,-14.034,-14.0253,-14.0016,-13.9989,-13.9948,-13.9737,-13.9645,-13.9482,-13.9227,-13.8972,-13.8717,-13.8462,-13.8194,-13.8163,-13.7886,-13.7874,-13.7868,-13.7864,-13.7834,-13.7564,-13.7265,-13.7162,-13.6966,-13.6667,-13.6368,-13.6069,-13.577,-13.5471,-13.5323,-13.5172,-13.5021,-13.4869,-13.4717,-13.4565,-13.4413,-13.4261,-13.4109,-13.3958,-13.3806,-13.3654,-13.3502,-13.335,-13.3198,-13.3046,-13.2894,-13.2742,-13.259,-13.2438,-13.2141,-13.2106,-13.1845,-13.1548,-13.1252,-13.0956,-13.0669,-13.0659,-13.0363,-13.0174,-12.9989,-12.9745,-12.9537,-12.9506,-12.9247,-12.8988,-12.8729,-12.847,-12.8211,-12.7952,-12.7701,-12.7449,-12.7198,-12.703,-12.6861,-12.6693,-12.6524,-12.651,-12.6429,-12.6334,-12.6239,-12.6147,-12.6144,-12.6049,-12.5954,-12.5942,-12.5859,-12.5866,-12.587,-12.5878,-12.5885,-12.5892,-12.5904,-12.591,-12.5914,-12.592,-12.5928,-12.5937,-12.594,-12.5948,-12.5953,-12.596,-12.5967,-12.597,-12.5987,-12.5992,-12.6,-12.6006,-12.6006,-12.6001,-12.6001,-12.6005,-12.6005,-12.6011,-12.6011,-12.6007,-12.6013,-12.6008,-12.6006,-12.6014,-12.6019,-12.602,-12.6022,-12.6022,-12.6026,-12.6031,-12.6032,-12.6038,-12.6042,-12.6042,-12.6049,-12.6053,-12.6051,-12.6058,-12.6066,-12.6074,-12.6087,-12.6098,-12.6107,-12.6118,-12.6128,-12.6129,-12.6126,-12.612,-12.611,-12.6102,-12.6117,-12.6126,-12.614,-12.6144,-12.6156,-12.6155,-12.6152,-12.615,-12.6149,-12.6141,-12.6141,-12.6146,-12.6155,-12.6163,-12.6177,-12.6192,-12.6205,-12.6213,-12.6209,-12.6212,-12.6211,-12.6213,-12.6211,-12.6219,-12.6228,-12.6235,-12.6237,-12.6247,-12.6444,-12.6641,-12.6861,-12.6911,-12.6961,-12.6998,-12.7036,-12.7073,-12.7103,-12.7134,-12.717,-12.72,-12.7236,-12.7281,-12.7317,-12.7361,-12.7403,-12.7448,-12.7498,-12.7542,-12.7592,-12.7636,-12.7686,-12.7803,-12.7811,-12.7825,-12.785,-12.7845,-12.7825,-12.7825,-12.7828,-12.7831,-12.7834,-12.785,-12.7873,-12.7909,-12.7939,-12.7981,-12.8075,-12.815,-12.8209,-12.8259,-12.8314,-12.837,-12.8436,-12.8489,-12.8525,-12.8553,-12.857,-12.8606,-12.8645,-12.8636,-12.8673,-12.8717,-12.8753,-12.8784,-12.882,-12.8856,-12.8886,-12.8903,-12.8925,-12.8898,-12.8856,-12.8792,-12.8734,-12.8686,-12.8628,-12.857,-12.8506,-12.8464,-12.8431,-12.8417,-12.8406,-12.8392,-12.8392,-12.8403,-12.8425,-12.845,-12.8486,-12.8514,-12.855,-12.8581,-12.862,-12.8664,-12.8706,-12.8742,-12.8786,-12.8836,-12.8881,-12.8934,-12.8975,-12.902,-12.9064,-12.9114,-12.9164,-12.9217,-12.9253,-12.9289,-12.9317,-12.9342,-12.9359,-12.9375,-12.9392,-12.9395,-12.9395,-12.9398,-12.9392,-12.9386,-12.9375,-12.9356,-12.9328,-12.9336,-12.9386,-12.9423,-12.9475,-12.9492,-12.9492,-12.9509,-12.9523,-12.9523,-12.9539,-12.9556,-12.955,-12.9553,-12.9548,-12.9561,-12.9566,-12.9595,-12.9639,-12.9709,-12.9764,-12.9814,-12.9864,-12.9928,-12.9984,-12.9992,-13.0034,-13.0084,-13.0148,-13.0211,-13.0275,-13.0325,-13.0373,-13.0417,-13.0467,-13.0523,-13.0592,-13.065,-13.0686,-13.0731,-13.0767,-13.082,-13.0864,-13.0909,-13.0931,-13.0932,-13.0936,-13.0967,-13.0992,-13.1031,-13.1036,-13.1031,-13.1017,-13.0989,-13.0964,-13.0931,-13.0903,-13.0889,-13.0886,-13.0895,-13.0898,-13.0906,-13.0914,-13.0959,-13.1011,-13.1053,-13.1098,-13.1148,-13.1192,-13.1234,-13.1278,-13.1314,-13.1367,-13.1403,-13.1453,-13.1506,-13.1564,-13.1614,-13.165,-13.1686,-13.1731,-13.1781,-13.1839,-13.1903,-13.1967,-13.2023,-13.2075,-13.2131,-13.2181,-13.225,-13.23,-13.2339,-13.2361,-13.2392,-13.2417,-13.2445,-13.2467,-13.2484,-13.25,-13.2503,-13.2489,-13.247,-13.245,-13.2417,-13.2361,-13.2295,-13.2239,-13.2211,-13.2192,-13.2186,-13.2217,-13.2253,-13.2289,-13.2325,-13.2359,-13.2395,-13.2439,-13.2475,-13.2511,-13.2548,-13.2592,-13.2636,-13.2673,-13.2711,-13.2742,-13.2764,-13.2778,-13.2809,-13.2831,-13.2848,-13.2859,-13.2859,-13.2839,-13.2789,-13.2748,-13.2714,-13.2728,-13.2764,-13.2806,-13.2842,-13.2878,-13.2906,-13.2945,-13.2967,-13.2989,-13.3006,-13.3025,-13.3039,-13.3048,-13.305,-13.3053,-13.3064,-13.3078,-13.3103,-13.3131,-13.3161,-13.3198,-13.3198,-13.3242,-13.3259,-13.3253,-13.3248,-13.3248,-13.3245,-13.3239,-13.327,-13.33,-13.3336,-13.3381,-13.3417,-13.347,-13.3511,-13.3556,-13.3592,-13.3625,-13.3648,-13.3678,-13.37,-13.3731,-13.3753,-13.3775,-13.3786,-13.3803,-13.3804,-13.3804,-13.3817,-13.3828,-13.3836,-13.3848,-13.3856,-13.3873,-13.3889,-13.3906,-13.392,-13.3945,-13.3973,-13.4017,-13.4061,-13.4114,-13.4156,-13.4217,-13.4264,-13.4309,-13.4359,-13.4411,-13.4448,-13.4492,-13.4523,-13.4559,-13.4598,-13.462,-13.465,-13.4681,-13.4714,-13.4748,-13.4778,-13.4811,-13.4839,-13.4873,-13.49,-13.4942,-13.4992,-13.5006,-13.505,-13.5081,-13.5111,-13.5134,-13.515,-13.5145,-13.5134,-13.5098,-13.5056,-13.5014,-13.4992,-13.4973,-13.4923,-13.4886,-13.4861,-13.4842,-13.485,-13.4881,-13.4925,-13.4945,-13.4992,-13.5003,-13.5042,-13.5081,-13.5106,-13.5148,-13.5198,-13.522,-13.5284,-13.5348,-13.5411,-13.5478,-13.5534,-13.56,-13.5659,-13.5723,-13.5792,-13.5859,-13.5909,-13.5964,-13.6014,-13.6056,-13.6106,-13.6148,-13.6203,-13.6253,-13.6325,-13.6389,-13.6445,-13.6489,-13.6523,-13.6556,-13.6592,-13.6639,-13.6703,-13.6761,-13.6786,-13.6823,-13.6853,-13.6889,-13.6928,-13.6964,-13.7009,-13.7053,-13.7089,-13.712,-13.712,-13.7109,-13.7089,-13.7061,-13.7036,-13.7009,-13.6981,-13.6979,-13.697,-13.697,-13.6995,-13.7039,-13.7103,-13.7153,-13.72,-13.7242,-13.7278,-13.732,-13.7361,-13.7398,-13.7439,-13.7464,-13.7506,-13.7548,-13.76,-13.7656,-13.7711,-13.7784,-13.7848,-13.7911,-13.7967,-13.8025,-13.8084,-13.8131,-13.8173,-13.8223,-13.8267,-13.8309,-13.8348,-13.8384,-13.8425,-13.8473,-13.8506,-13.8528,-13.855,-13.8573,-13.8617,-13.865,-13.8673,-13.8681,-13.8692,-13.8709,-13.8703,-13.8706,-13.8709,-13.8711,-13.872,-13.8742,-13.8778,-13.8845,-13.8895,-13.8953,-13.9011,-13.9061,-13.9092,-13.9114,-13.9153,-13.9181,-13.9225,-13.927,-13.9306,-13.935,-13.94,-13.9459,-13.95,-13.9531,-13.9595,-13.9653,-13.9689,-13.9725,-13.9728,-13.9717,-13.9698,-13.9684,-13.9664,-13.9639,-13.9639,-13.9661,-13.9706,-13.975,-13.9786,-13.9825,-13.9842,-13.9836,-13.9817,-13.9781,-13.9725,-13.9678,-13.9636,-13.9606,-13.9631,-13.9667,-13.9711,-13.975,-13.9786,-13.9823,-13.9861,-13.9898,-13.9936,-13.9981,-13.9992,-14.0036,-14.0103,-14.0173,-14.0239,-14.0289,-14.0325,-14.0356,-14.0392,-14.0417,-14.0434,-14.0464,-14.0486,-14.0531,-14.0573,-14.062,-14.0664,-14.07,-14.0736,-14.0773,-14.0811,-14.0853,-14.09,-14.0945,-14.0981,-14.1025,-14.1061,-14.11,-14.1114,-14.1148,-14.1175,-14.1211,-14.1256,-14.13,-14.1339,-14.1375,-14.1406,-14.1431,-14.1453,-14.1484,-14.1528,-14.1578,-14.1623,-14.1667,-14.1703,-14.1728,-14.1764,-14.1785,-14.1785,-14.1795,-14.1831,-14.1861,-14.1898,-14.1942,-14.1981,-14.2011,-14.2039,-14.2073,-14.21,-14.2125,-14.2161,-14.2206,-14.2242,-14.2295,-14.2339,-14.2375,-14.2409,-14.2453,-14.2495,-14.2528,-14.2564,-14.2592,-14.2625,-14.2684,-14.2728,-14.2759,-14.2775,-14.277,-14.2773,-14.2767,-14.2764,-14.2764,-14.2775,-14.2825,-14.2889,-14.2945,-14.2998,-14.3039,-14.3086,-14.3136,-14.3195,-14.325,-14.3303,-14.3348,-14.3378,-14.3389,-14.337,-14.3356,-14.3345,-14.3348,-14.3342,-14.3336,-14.3331,-14.3361,-14.34,-14.3445,-14.35,-14.3545,-14.355,-14.3625,-14.3695,-14.3767,-14.3834,-14.3903,-14.3953,-14.3992,-14.3992,-14.4006,-14.4056,-14.41,-14.4153,-14.4189,-14.4236,-14.4292,-14.4342,-14.4392,-14.4425,-14.4453,-14.4478,-14.4498,-14.4534,-14.4584,-14.4656,-14.472,-14.4792,-14.4856,-14.4911,-14.4978,-14.4992,-14.5034,-14.5106,-14.517,-14.5228,-14.5286,-14.5345,-14.5395,-14.5434,-14.5478,-14.5523,-14.5586,-14.5628,-14.5661,-14.5706,-14.5753,-14.5717,-14.5784,-14.5853,-14.5925,-14.5984,-14.6042,-14.6184,-14.6256,-14.6328,-14.6373,-14.6389,-14.6406,-14.6442,-14.6514,-14.6578,-14.6642,-14.6709,-14.6773,-14.6828,-14.6886,-14.6942,-14.7,-14.7064,-14.7123,-14.7186,-14.7239,-14.7281,-14.7306,-14.7342,-14.7398,-14.7448,-14.7473,-14.7495,-14.7525,-14.7556,-14.7595,-14.765,-14.7709,-14.7745,-14.7792,-14.7834,-14.7886,-14.7942,-14.7992,-14.8031,-14.8067,-14.8111,-14.817,-14.8228,-14.8267,-14.8309,-14.8361,-14.8417,-14.8481,-14.8539,-14.8595,-14.861,-14.8639,-14.8678,-14.8714,-14.8753,-14.8795,-14.8839,-14.8884,-14.8917,-14.8945,-14.8978,-14.9014,-14.9061,-14.912,-14.9142,-14.9167,-14.9234,-14.9281,-14.9345,-14.9403,-14.9409,-14.9459],"lat":[16.6424,16.6449,16.6496,16.6552,16.6613,16.6674,16.6721,16.6735,16.6796,16.6843,16.689,16.6902,16.6921,16.6927,16.691,16.6879,16.6852,16.6815,16.6795,16.6779,16.6738,16.669,16.664,16.659,16.6543,16.6499,16.6457,16.6415,16.6379,16.6352,16.6335,16.6313,16.6296,16.6299,16.6313,16.6346,16.6377,16.6418,16.6465,16.6513,16.656,16.6607,16.664,16.6674,16.6704,16.6718,16.6707,16.6685,16.6649,16.6621,16.6585,16.6543,16.6479,16.6424,16.6371,16.6329,16.6304,16.6288,16.6268,16.6249,16.6202,16.614,16.6077,16.6035,16.5999,16.5963,16.5927,16.5893,16.5857,16.5835,16.5824,16.5829,16.5832,16.5846,16.5857,16.5852,16.5857,16.5849,16.5824,16.5802,16.5793,16.5771,16.5754,16.5738,16.5721,16.5707,16.5685,16.566,16.5638,16.5615,16.5593,16.5585,16.5588,16.5599,16.5618,16.5629,16.5649,16.5665,16.5679,16.5696,16.5715,16.5727,16.5732,16.5721,16.5721,16.571,16.5696,16.5671,16.5649,16.5635,16.561,16.561,16.5613,16.5596,16.5582,16.556,16.5543,16.5521,16.5496,16.5468,16.5452,16.5443,16.5432,16.5438,16.5463,16.5496,16.5529,16.5563,16.561,16.5657,16.5713,16.5752,16.579,16.5818,16.5835,16.5827,16.581,16.5793,16.5771,16.5749,16.5749,16.5727,16.571,16.5693,16.5671,16.567,16.5649,16.5613,16.5563,16.5502,16.544,16.5385,16.5327,16.5279,16.5243,16.5213,16.5185,16.5168,16.5146,16.5129,16.5121,16.5118,16.5121,16.5127,16.5129,16.514,16.5152,16.5171,16.519,16.5207,16.5227,16.5229,16.5221,16.5185,16.5129,16.5074,16.5032,16.4997,16.4996,16.4982,16.4954,16.4924,16.4902,16.4885,16.4863,16.4846,16.4824,16.4807,16.4799,16.479,16.4785,16.479,16.4802,16.4813,16.4824,16.4843,16.4854,16.4865,16.4877,16.4882,16.4885,16.4896,16.4893,16.4904,16.4918,16.4929,16.494,16.4957,16.4968,16.4988,16.4991,16.4996,16.4997,16.4998,16.5007,16.5018,16.5024,16.5013,16.4998,16.4996,16.499,16.4968,16.4957,16.4949,16.4946,16.4957,16.4982,16.4996,16.4997,16.4997,16.501,16.504,16.5068,16.5093,16.5113,16.5124,16.5127,16.5124,16.5121,16.5104,16.5082,16.5065,16.5043,16.5021,16.4999,16.4997,16.4996,16.4974,16.4946,16.4915,16.4893,16.4863,16.486,16.4874,16.4885,16.4896,16.4921,16.4932,16.493,16.4929,16.4927,16.4918,16.4893,16.4871,16.4849,16.4832,16.4829,16.4815,16.4804,16.4799,16.4804,16.4829,16.4849,16.4882,16.4921,16.4954,16.4988,16.4996,16.4997,16.5015,16.504,16.5065,16.5093,16.5118,16.5152,16.5185,16.5232,16.5293,16.534,16.5388,16.5429,16.546,16.5493,16.5513,16.5524,16.5521,16.5518,16.5515,16.5499,16.5468,16.544,16.5404,16.5368,16.5338,16.5302,16.5265,16.5259,16.5231,16.5195,16.5183,16.5211,16.5327,16.5334,16.5311,16.5239,16.5205,16.5171,16.5052,16.4936,16.4797,16.4745,16.4698,16.4642,16.4566,16.4491,16.4109,16.4037,16.4021,16.4029,16.4029,16.3997,16.3941,16.375,16.3703,16.3627,16.3512,16.344,16.34,16.3307,16.3289,16.3241,16.3149,16.3014,16.3,16.2891,16.2807,16.2763,16.2723,16.2668,16.2584,16.2473,16.2361,16.227,16.221,16.2166,16.2161,16.215,16.2154,16.2174,16.2182,16.2178,16.2126,16.2079,16.2028,16.2019,16.1943,16.1864,16.1784,16.1673,16.1633,16.1537,16.1494,16.1438,16.1374,16.1287,16.1259,16.1179,16.1135,16.1119,16.1095,16.1028,16.0984,16.094,16.0892,16.0877,16.0845,16.0826,16.0825,16.0804,16.0796,16.0789,16.0776,16.0763,16.0713,16.0684,16.0684,16.07,16.07,16.0692,16.0692,16.0683,16.0671,16.0667,16.0667,16.0659,16.0658,16.0642,16.0642,16.0634,16.0634,16.0604,16.0579,16.0571,16.0471,16.0437,16.0404,16.0371,16.0346,16.0329,16.0321,16.0304,16.0292,16.0304,16.0337,16.0375,16.0375,16.0367,16.0358,16.0334,16.0329,16.0317,16.0308,16.0271,16.0245,16.0237,16.0179,16.0129,16.0104,15.9996,15.997,15.9954,15.9921,15.9913,15.9896,15.9879,15.9787,15.977,15.9763,15.975,15.975,15.9779,15.9796,15.9829,15.985,15.9867,15.9867,15.9875,15.9896,15.9912,15.9937,15.9942,15.9913,15.9896,15.9887,15.9867,15.9858,15.985,15.985,15.9883,15.9879,15.9829,15.9817,15.9817,15.9825,15.9825,15.9817,15.9842,15.9842,15.9813,15.9795,15.9762,15.9721,15.9712,15.9696,15.9666,15.9667,15.965,15.965,15.9629,15.9612,15.9538,15.9529,15.9504,15.9495,15.9487,15.9471,15.9454,15.9437,15.9421,15.9404,15.9396,15.938,15.9371,15.9321,15.9313,15.9296,15.9284,15.9284,15.9292,15.9279,15.927,15.9262,15.9258,15.9238,15.9187,15.9137,15.9105,15.9037,15.9034,15.9046,15.9095,15.911,15.9129,15.9134,15.9121,15.9071,15.9063,15.9046,15.9017,15.9063,15.9137,15.9171,15.9221,15.9229,15.9254,15.9258,15.9296,15.9329,15.9338,15.946,15.9487,15.9496,15.9571,15.9596,15.9605,15.9637,15.9663,15.9746,15.9754,15.9788,15.98,15.9792,15.9779,15.9763,15.9696,15.9671,15.9654,15.9638,15.9621,15.9588,15.9546,15.9538,15.9529,15.9512,15.9505,15.9496,15.9471,15.9462,15.9363,15.9355,15.9271,15.9262,15.9221,15.9179,15.9154,15.9146,15.9062,15.9054,15.9038,15.9021,15.8988,15.8979,15.8912,15.8904,15.8813,15.8804,15.8738,15.8729,15.8596,15.8588,15.8563,15.8554,15.8513,15.8504,15.8479,15.8475,15.8479,15.8502,15.8537,15.8563,15.8587,15.8621,15.8662,15.8662,15.8604,15.8596,15.8562,15.8545,15.8488,15.8463,15.8454,15.8429,15.8346,15.8338,15.8296,15.8287,15.8221,15.8213,15.8163,15.8112,15.8079,15.8071,15.8029,15.7979,15.7962,15.7904,15.7887,15.7863,15.7829,15.782,15.7771,15.7712,15.7696,15.7671,15.7621,15.7604,15.7579,15.7563,15.7537,15.7521,15.7471,15.7437,15.7412,15.7396,15.7371,15.7338,15.7321,15.7304,15.7263,15.7196,15.7179,15.7087,15.7054,15.7013,15.6971,15.6921,15.687,15.6787,15.6763,15.668,15.6662,15.6654,15.6521,15.6421,15.6387,15.6337,15.6313,15.6296,15.627,15.6254,15.6187,15.6129,15.6121,15.6045,15.5988,15.5854,15.5813,15.5754,15.5746,15.5721,15.5696,15.5671,15.5663,15.5621,15.5596,15.5545,15.5537,15.5471,15.5462,15.5412,15.5355,15.5329,15.532,15.5279,15.5246,15.5196,15.5187,15.5121,15.5095,15.5071,15.5046,15.5021,15.4979,15.4946,15.4887,15.4821,15.4813,15.4721,15.4687,15.4663,15.4621,15.4614,15.4525,15.4471,15.443,15.4379,15.4348,15.4333,15.4303,15.4278,15.4259,15.4226,15.422,15.4167,15.4086,15.4068,15.403,15.4021,15.3963,15.3905,15.3866,15.3804,15.3726,15.3693,15.3645,15.364,15.3577,15.3549,15.3523,15.3406,15.3395,15.3313,15.3241,15.3166,15.3152,15.31,15.3026,15.2994,15.2991,15.297,15.2967,15.2948,15.2949,15.2939,15.2917,15.291,15.2904,15.2903,15.2892,15.2835,15.282,15.2796,15.2768,15.2737,15.2722,15.2693,15.2685,15.2667,15.265,15.2629,15.2625,15.2603,15.2561,15.2541,15.2522,15.2486,15.2461,15.2453,15.2426,15.2409,15.2405,15.2367,15.2357,15.2304,15.2283,15.2275,15.2271,15.226,15.2247,15.2226,15.2178,15.2115,15.2086,15.2034,15.2002,15.1985,15.1957,15.1891,15.1824,15.1785,15.1746,15.1692,15.166,15.1621,15.1563,15.1457,15.1496,15.1526,15.1532,15.154,15.1548,15.1573,15.1649,15.1714,15.1783,15.1792,15.1797,15.1886,15.1924,15.1959,15.1987,15.2111,15.2141,15.2224,15.2254,15.2355,15.2389,15.2456,15.2486,15.2524,15.2585,15.2619,15.2662,15.2693,15.2687,15.2686,15.2688,15.2696,15.2703,15.2713,15.2726,15.2743,15.2766,15.2783,15.2795,15.281,15.2828,15.2836,15.2839,15.2844,15.2849,15.2861,15.2868,15.288,15.29,15.2917,15.2935,15.2953,15.2974,15.3001,15.3022,15.3045,15.3066,15.3083,15.31,15.3119,15.3133,15.3147,15.3158,15.3186,15.3208,15.3224,15.3237,15.325,15.3267,15.3281,15.3296,15.3314,15.3332,15.3341,15.3361,15.3378,15.3393,15.3417,15.3434,15.3461,15.3476,15.3499,15.3518,15.3545,15.356,15.3569,15.3582,15.3591,15.3591,15.3577,15.356,15.3536,15.3511,15.3482,15.3455,15.3431,15.3406,15.3368,15.3339,15.3311,15.3276,15.3258,15.3252,15.3257,15.3255,15.3245,15.3233,15.3212,15.3188,15.3171,15.3164,15.3174,15.3192,15.3213,15.3233,15.3239,15.3242,15.3235,15.3213,15.3196,15.3167,15.3136,15.3105,15.3088,15.3056,15.3029,15.2994,15.2957,15.2925,15.2878,15.2847,15.282,15.2793,15.2769,15.2758,15.2743,15.2715,15.2689,15.266,15.2636,15.2609,15.2579,15.2561,15.2537,15.2513,15.25,15.249,15.2488,15.2509,15.2521,15.2548,15.2596,15.2634,15.2662,15.2691,15.2706,15.2723,15.2725,15.2728,15.273,15.2725,15.2722,15.2717,15.2703,15.2696,15.269,15.2672,15.2646,15.2618,15.2596,15.257,15.2536,15.2498,15.2467,15.2425,15.2392,15.2357,15.2313,15.2288,15.2247,15.2232,15.2201,15.2176,15.2149,15.21,15.2049,15.2006,15.1971,15.1941,15.1907,15.1883,15.1855,15.1831,15.1807,15.1771,15.1741,15.1712,15.1692,15.1677,15.1658,15.163,15.1607,15.1563,15.1537,15.15,15.1443,15.1382,15.1345,15.131,15.1278,15.1247,15.1219,15.1194,15.118,15.1162,15.1147,15.1113,15.1079,15.1068,15.1046,15.1014,15.1003,15.0971,15.0954,15.0934,15.0912,15.0863,15.0805,15.0766,15.0711,15.0663,15.0617,15.0561,15.0491,15.0454,15.0417,15.037,15.0337,15.0306,15.0273,15.024,15.0219,15.0182,15.0161,15.0118,15.0115,15.0092,15.0078,15.0058,15.0015,15.0006,14.9953,14.9934,14.9908,14.9887,14.9865,14.9844,14.9827,14.9803,14.9784,14.9767,14.9752,14.9732,14.9712,14.9695,14.9682,14.9662,14.9648,14.9629,14.963,14.9633,14.9632,14.9636,14.9635,14.9626,14.9628,14.9631,14.9631,14.9626,14.9621,14.9617,14.9612,14.9608,14.9609,14.9603,14.9597,14.9585,14.9579,14.9574,14.957,14.9564,14.9564,14.9565,14.9565,14.9566,14.9573,14.9575,14.958,14.9581,14.9585,14.9588,14.9597,14.9602,14.9608,14.9607,14.9611,14.961,14.961,14.9606,14.9605,14.9604,14.9598,14.9599,14.9597,14.9594,14.9593,14.9595,14.959,14.9587,14.959,14.959,14.9589,14.9586,14.9577,14.9568,14.9567,14.9563,14.9559,14.9563,14.9562,14.9555,14.9557,14.9554,14.9553,14.9553,14.9557,14.9564,14.9583,14.9599,14.9618,14.9634,14.9653,14.967,14.9686,14.9701,14.9708,14.9716,14.9721,14.9723,14.9726,14.9721,14.9722,14.9722,14.9722,14.9726,14.9731,14.9736,14.9741,14.9746,14.9749,14.9756,14.9761,14.9766,14.9768,14.9771,14.9776,14.9785,14.9788,14.9786,14.9788,14.9792,14.9802,14.9816,14.9826,14.9835,14.9849,14.9867,14.9887,14.9908,14.9934,14.9964,15.0027,15.0064,15.0109,15.0165,15.02,15.025,15.0334,15.0378,15.0417,15.0455,15.0494,15.0512,15.0511,15.0485,15.0441,15.0389,15.0327,15.0287,15.0234,15.0183,15.0123,15.0071,15.0033,15.0009,14.9993,14.9979,14.9801,14.9688,14.968,14.9661,14.9635,14.9619,14.9611,14.9601,14.9588,14.9568,14.9547,14.9521,14.9497,14.9479,14.9459,14.9438,14.9412,14.9391,14.9358,14.933,14.9306,14.9282,14.9249,14.9222,14.9179,14.915,14.9119,14.9085,14.9042,14.9013,14.8984,14.8977,14.8982,14.8987,14.9,14.9023,14.9062,14.9086,14.9132,14.9208,14.9314,14.9416,14.9459,14.9455,14.9421,14.9305,14.9206,14.9151,14.9102,14.903,14.8965,14.8884,14.8794,14.8794,14.8832,14.8859,14.8815,14.8721,14.8622,14.8578,14.8546,14.8572,14.8597,14.8615,14.8626,14.8635,14.8649,14.866,14.8673,14.8695,14.8697,14.8691,14.8692,14.869,14.869,14.8688,14.8687,14.8685,14.8684,14.8688,14.8686,14.8676,14.8674,14.8674,14.8673,14.8671,14.8669,14.8669,14.8659,14.8656,14.865,14.8645,14.8641,14.8634,14.8632,14.8621,14.8611,14.8597,14.8591,14.8577,14.856,14.8551,14.8531,14.8505,14.8482,14.8469,14.8448,14.8441,14.843,14.8417,14.8395,14.8378,14.8352,14.8335,14.8316,14.8292,14.8273,14.825,14.8237,14.8237,14.8245,14.8253,14.8259,14.8272,14.8287,14.8304,14.8314,14.8322,14.8333,14.8339,14.8349,14.8358,14.8366,14.837,14.8377,14.8386,14.8389,14.8393,14.8395,14.8395,14.8398,14.8401,14.8401,14.8395,14.839,14.8384,14.8377,14.8371,14.8363,14.8358,14.8355,14.8346,14.8334,14.8322,14.8308,14.8293,14.8281,14.8259,14.8254,14.8249,14.8241,14.823,14.8222,14.8208,14.8177,14.8155,14.8127,14.8101,14.8078,14.8051,14.8003,14.7966,14.7929,14.7886,14.7842,14.781,14.7766,14.7726,14.7667,14.763,14.7582,14.7532,14.7476,14.7445,14.728,14.7115,14.695,14.6786,14.6621,14.6569,14.6516,14.6464,14.6412,14.6359,14.6307,14.6255,14.6218,14.6212,14.6202,14.615,14.6097,14.6348,14.6457,14.6516,14.6684,14.6851,14.7019,14.6785,14.6551,14.6551,14.6413,14.6233,14.6206,14.6183,14.6181,14.6179,14.6162,14.6151,14.614,14.612,14.6057,14.6057,14.6005,14.6005,14.5992,14.5986,14.598,14.5926,14.5877,14.5829,14.5829,14.583,14.5831,14.5827,14.5831,14.583,14.5827,14.5827,14.5818,14.5817,14.5807,14.5804,14.5793,14.5763,14.5761,14.5761,14.5758,14.5758,14.5755,14.5754,14.5745,14.5735,14.5738,14.5747,14.5746,14.5745,14.5742,14.5717,14.568,14.5678,14.5667,14.5638,14.5627,14.561,14.561,14.5604,14.562,14.5623,14.5628,14.5632,14.5627,14.5623,14.5621,14.5621,14.5619,14.5619,14.562,14.562,14.5626,14.563,14.5634,14.5639,14.564,14.5647,14.5651,14.5654,14.5655,14.5658,14.5665,14.5675,14.5702,14.5724,14.5743,14.577,14.5776,14.5779,14.5787,14.5787,14.579,14.5807,14.5802,14.5791,14.5799,14.58,14.5801,14.5812,14.5823,14.5823,14.5838,14.5838,14.5864,14.5864,14.5864,14.5864,14.586,14.5845,14.5758,14.5747,14.573,14.5709,14.5588,14.5588,14.5578,14.5578,14.5412,14.5404,14.5395,14.5385,14.5385,14.5335,14.531,14.531,14.5304,14.5145,14.4984,14.4973,14.493,14.4911,14.4864,14.4859,14.484,14.4741,14.4696,14.4617,14.4493,14.437,14.4246,14.4123,14.4117,14.4116,14.411,14.411,14.4109,14.4109,14.4109,14.4103,14.4096,14.4094,14.4089,14.4083,14.4076,14.4069,14.4062,14.4055,14.4052,14.4048,14.4163,14.4278,14.4393,14.4507,14.4622,14.4737,14.4852,14.4967,14.5081,14.5196,14.5311,14.5425,14.554,14.5655,14.577,14.5884,14.5999,14.6114,14.6069,14.6064,14.6025,14.5981,14.5936,14.5892,14.5849,14.5847,14.5803,14.5764,14.5726,14.5654,14.5592,14.5583,14.562,14.5657,14.5694,14.573,14.5767,14.5804,14.5915,14.6025,14.6135,14.6362,14.6589,14.6816,14.7043,14.7087,14.7329,14.7616,14.7902,14.818,14.8188,14.8474,14.8761,14.8798,14.9052,14.9057,14.9071,14.9079,14.9091,14.91,14.911,14.912,14.9129,14.9139,14.9151,14.9169,14.9179,14.9194,14.9206,14.922,14.9229,14.9245,14.9263,14.9279,14.9294,14.9308,14.9321,14.9333,14.9341,14.9356,14.9363,14.9375,14.9386,14.9401,14.9412,14.943,14.944,14.9454,14.9466,14.9478,14.9492,14.9502,14.9512,14.9521,14.953,14.9552,14.956,14.9568,14.958,14.9593,14.9608,14.9623,14.9628,14.9632,14.9634,14.964,14.9647,14.965,14.9656,14.9666,14.9675,14.9682,14.969,14.9703,14.9719,14.9724,14.9726,14.973,14.9742,14.9751,14.9759,14.9766,14.9775,14.9783,14.9794,14.98,14.9808,14.9813,14.9822,14.9824,14.9829,14.9837,14.985,14.9861,14.9872,14.9882,14.9905,14.9916,14.9933,14.9955,14.9973,15.0004,15.0273,15.0543,15.087,15.0857,15.0877,15.0907,15.0943,15.0974,15.1015,15.1057,15.109,15.1129,15.1163,15.119,15.1224,15.1249,15.1274,15.1302,15.1321,15.1346,15.1365,15.1393,15.141,15.1479,15.1485,15.1493,15.1568,15.1635,15.1704,15.1774,15.184,15.191,15.1979,15.2035,15.2082,15.2102,15.2132,15.2154,15.2165,15.2154,15.2104,15.2082,15.2068,15.2052,15.2057,15.2077,15.211,15.2149,15.2204,15.2224,15.2265,15.2271,15.2304,15.2332,15.2365,15.2404,15.2438,15.2471,15.251,15.2565,15.2613,15.2657,15.2685,15.2679,15.2668,15.2649,15.2635,15.2624,15.2632,15.2663,15.2699,15.2754,15.281,15.2865,15.2918,15.2982,15.3029,15.3077,15.311,15.3149,15.3182,15.3224,15.3257,15.3282,15.331,15.3343,15.3368,15.3388,15.3415,15.3435,15.346,15.3485,15.3513,15.3529,15.3549,15.3568,15.3602,15.3635,15.3677,15.3724,15.3779,15.3832,15.3888,15.3957,15.4027,15.4093,15.4157,15.4218,15.4274,15.4324,15.4365,15.4413,15.4432,15.4452,15.4471,15.4527,15.4579,15.4621,15.4663,15.4649,15.469,15.4743,15.4804,15.4874,15.4938,15.499,15.4996,15.5032,15.5057,15.5057,15.504,15.5018,15.4996,15.4988,15.4971,15.4967,15.4949,15.4927,15.4918,15.491,15.4899,15.4877,15.4857,15.4827,15.4804,15.479,15.4788,15.4799,15.4832,15.4857,15.489,15.491,15.4938,15.4963,15.4996,15.4997,15.5004,15.5043,15.5082,15.5193,15.5193,15.5254,15.531,15.5354,15.5396,15.5429,15.5474,15.5529,15.559,15.5652,15.5707,15.5768,15.5829,15.5857,15.5877,15.5902,15.5929,15.5946,15.5974,15.5999,15.6027,15.606,15.6079,15.6113,15.6129,15.6149,15.616,15.6179,15.6199,15.6232,15.626,15.6279,15.629,15.6282,15.6274,15.6257,15.6235,15.6218,15.6196,15.6193,15.6213,15.6246,15.6293,15.6335,15.6382,15.6424,15.6471,15.6524,15.6579,15.6649,15.6704,15.6752,15.6802,15.6838,15.6852,15.6849,15.6863,15.6904,15.6954,15.7015,15.7057,15.709,15.7124,15.7157,15.7196,15.7229,15.7257,15.729,15.7324,15.7354,15.7382,15.7407,15.744,15.7474,15.7515,15.7549,15.7563,15.7602,15.7649,15.7704,15.7765,15.7821,15.7871,15.7893,15.7921,15.7957,15.8013,15.8046,15.8079,15.811,15.8143,15.8185,15.8218,15.8265,15.8313,15.8368,15.8421,15.8477,15.8538,15.8607,15.8677,15.8738,15.8793,15.884,15.8865,15.8893,15.8913,15.8918,15.8946,15.8985,15.9049,15.911,15.9165,15.9227,15.929,15.9329,15.9371,15.9404,15.9429,15.9463,15.9482,15.9507,15.9535,15.9568,15.9607,15.9654,15.9696,15.9743,15.9782,15.9832,15.9879,15.994,15.9993,15.9996,15.9997,16.0049,16.011,16.0174,16.0235,16.0296,16.0352,16.0404,16.046,16.0513,16.0549,16.0588,16.0615,16.064,16.066,16.0685,16.0699,16.0715,16.0743,16.0763,16.0779,16.0813,16.084,16.0879,16.0913,16.0946,16.0993,16.1035,16.1074,16.104,16.1004,16.096,16.0927,16.0882,16.0849,16.0804,16.0777,16.0781,16.0782,16.0807,16.0854,16.0896,16.0943,16.0996,16.106,16.1115,16.1149,16.1179,16.1207,16.1223,16.1238,16.126,16.1296,16.1338,16.1385,16.1449,16.1488,16.1492,16.1493,16.148,16.1477,16.1449,16.1413,16.1371,16.134,16.1318,16.1318,16.131,16.1313,16.1318,16.1324,16.1335,16.134,16.1352,16.1357,16.1354,16.1346,16.1324,16.1307,16.1285,16.1254,16.1235,16.1204,16.1188,16.1165,16.1163,16.1154,16.1138,16.111,16.1074,16.104,16.1004,16.0982,16.0971,16.0985,16.1018,16.1052,16.109,16.1124,16.1157,16.119,16.1215,16.1243,16.1277,16.1315,16.1385,16.144,16.149,16.1532,16.1574,16.1615,16.1657,16.1667,16.1713,16.1768,16.1815,16.1843,16.1846,16.1824,16.1802,16.1774,16.1738,16.1707,16.1679,16.1643,16.1615,16.1574,16.1543,16.1515,16.1493,16.1471,16.1454,16.1452,16.1443,16.1435,16.1418,16.1402,16.1385,16.1363,16.1335,16.1313,16.1282,16.1254,16.1227,16.119,16.116,16.1138,16.1179,16.1227,16.1274,16.1321,16.1346,16.1388,16.1435,16.1496,16.1543,16.1599,16.166,16.1729,16.1799,16.1868,16.1929,16.1977,16.1996,16.2002,16.2018,16.2032,16.2029,16.2049,16.2088,16.2121,16.2154,16.2196,16.2221,16.2246,16.2279,16.2307,16.2327,16.2324,16.2307,16.2302,16.229,16.2302,16.2335,16.2368,16.2424,16.2479,16.2529,16.2585,16.2632,16.2674,16.2729,16.2777,16.2804,16.2829,16.2863,16.2896,16.2952,16.3013,16.3063,16.3096,16.3113,16.3135,16.3163,16.3207,16.3254,16.3288,16.3313,16.3346,16.3379,16.3413,16.3446,16.3479,16.3513,16.3538,16.354,16.3549,16.3554,16.3552,16.3557,16.3574,16.3607,16.3649,16.3682,16.3729,16.3782,16.3824,16.3857,16.3882,16.391,16.3935,16.3963,16.3993,16.4027,16.406,16.4093,16.4121,16.4146,16.4171,16.4204,16.4229,16.4263,16.4296,16.4352,16.439,16.4432,16.4465,16.449,16.4518,16.4549,16.4582,16.4624,16.4671,16.4718,16.4757,16.4785,16.4802,16.4829,16.4854,16.4888,16.4935,16.4968,16.4996,16.4997,16.501,16.504,16.5082,16.5115,16.514,16.5174,16.5215,16.5254,16.5293,16.5335,16.5382,16.5415,16.544,16.5474,16.5493,16.5504,16.551,16.551,16.5493,16.5465,16.5443,16.5407,16.5365,16.5329,16.534,16.5374,16.5415,16.5468,16.5532,16.5588,16.5649,16.571,16.5765,16.5827,16.5846,16.5838,16.5821,16.5799,16.5768,16.5746,16.5724,16.571,16.5721,16.5738,16.5765,16.5804,16.5865,16.5915,16.5971,16.6027,16.6096,16.6157,16.6218,16.6282,16.6321,16.6354,16.6382,16.6393,16.6404,16.6404,16.6402,16.6399,16.6396,16.6388,16.6385,16.6404,16.6424,16.6429,16.6435,16.6454,16.6482,16.6499,16.6532,16.6557,16.6571,16.6563,16.6538,16.6504,16.646,16.6418,16.6371,16.6335,16.6313,16.631,16.6299,16.6296,16.6288,16.6271,16.6263,16.6262,16.626,16.6257,16.6263,16.6274,16.6285,16.6296,16.6315,16.6349,16.6374,16.6402,16.639,16.6363,16.6327,16.6296,16.6288,16.6296,16.6288,16.6285,16.6282,16.6279,16.6279,16.6274,16.6271,16.6268,16.6293,16.6349,16.6402,16.6435,16.6432,16.6424,16.6427,16.6432,16.6424,16.6407,16.639,16.6374,16.6374,16.6377,16.6388,16.6393,16.6371,16.634,16.6299,16.6263,16.6246,16.6265,16.6313,16.636,16.6402,16.644,16.6474,16.6485,16.6468,16.6449,16.6424,16.6396,16.6374,16.6385,16.6404,16.6438,16.6468,16.6496,16.6507,16.649,16.6463,16.6432,16.641,16.6393,16.6385,16.6368,16.6379,16.6389,16.6407,16.644,16.6471,16.6504,16.6529,16.6557,16.6527,16.6485,16.6443,16.6407,16.6371,16.6349,16.636,16.6407,16.6454,16.646,16.6435,16.6427,16.6424,16.6418,16.6424]}]],[[{"lng":[-17.4742,-17.4734,-17.4726,-17.4713,-17.4713,-17.4734,-17.4736,-17.475,-17.4755,-17.4751,-17.4746,-17.4763,-17.4759,-17.4742],"lat":[14.6583,14.6575,14.6575,14.6562,14.6554,14.6534,14.6533,14.6533,14.6546,14.6553,14.6562,14.6571,14.6584,14.6583]},{"lng":[-17.4001,-17.398,-17.398,-17.4009,-17.4025,-17.4038,-17.403,-17.4038,-17.403,-17.4026,-17.4001],"lat":[14.6634,14.6654,14.668,14.6717,14.6725,14.6712,14.6696,14.6688,14.6679,14.6633,14.6634]},{"lng":[-17.5443,-17.5438,-17.5438,-17.5442,-17.545,-17.5455,-17.5455,-17.5451,-17.5443],"lat":[14.7417,14.7421,14.7446,14.745,14.745,14.7446,14.7421,14.7417,14.7417]},{"lng":[-17.5392,-17.5371,-17.5376,-17.5384,-17.5405,-17.5401,-17.5392],"lat":[14.7475,14.7496,14.75,14.75,14.7479,14.7475,14.7475]},{"lng":[-17.5134,-17.5142,-17.5159,-17.5167,-17.5184,-17.5197,-17.5196,-17.5182,-17.5159,-17.5142,-17.5129,-17.513,-17.5134],"lat":[14.7584,14.7591,14.7592,14.7583,14.7584,14.7571,14.7562,14.7561,14.7559,14.755,14.7562,14.758,14.7584]},{"lng":[-17.4775,-17.4771,-17.4771,-17.4771,-17.4776,-17.4792,-17.4796,-17.4796,-17.4797,-17.4792,-17.4775],"lat":[14.7692,14.7696,14.7711,14.7712,14.7716,14.7717,14.7713,14.7705,14.7695,14.7692,14.7692]},{"lng":[-16.7366,-16.7381,-16.7405,-16.7415,-16.7432,-16.7423,-16.7401,-16.7401,-16.738,-16.7373,-16.7338,-16.7337,-16.7285,-16.7277,-16.7274,-16.7253,-16.7247,-16.7231,-16.7209,-16.7187,-16.7155,-16.7127,-16.7075,-16.7033,-16.7031,-16.6995,-16.6973,-16.6908,-16.6889,-16.6887,-16.6883,-16.6888,-16.6886,-16.6892,-16.6901,-16.6911,-16.6958,-16.7012,-16.7102,-16.7307,-16.7313,-16.7347,-16.7363,-16.7405,-16.7421,-16.7455,-16.7455,-16.748,-16.7513,-16.7546,-16.7546,-16.7596,-16.7613,-16.763,-16.7638,-16.7655,-16.7663,-16.768,-16.7713,-16.773,-16.7746,-16.7746,-16.7763,-16.7763,-16.7805,-16.7821,-16.7838,-16.7863,-16.7872,-16.7913,-16.7988,-16.7988,-16.8055,-16.8063,-16.8088,-16.8155,-16.8188,-16.8213,-16.8238,-16.8263,-16.8296,-16.8321,-16.8338,-16.8363,-16.8363,-16.8388,-16.8397,-16.8413,-16.8446,-16.8463,-16.8463,-16.848,-16.8505,-16.8521,-16.8572,-16.8588,-16.8663,-16.8663,-16.8688,-16.8688,-16.873,-16.8738,-16.8763,-16.8771,-16.8838,-16.8838,-16.8871,-16.888,-16.893,-16.893,-16.9055,-16.9055,-16.9096,-16.9105,-16.913,-16.9196,-16.9205,-16.9238,-16.9305,-16.9313,-16.9388,-16.9405,-16.9496,-16.9497,-16.958,-16.9588,-16.9588,-16.9704,-16.9755,-16.9755,-16.9788,-16.9788,-16.988,-16.988,-16.9938,-16.9946,-16.9992,-17.0063,-17.0063,-17.0146,-17.0163,-17.0255,-17.0255,-17.0363,-17.0363,-17.0413,-17.0413,-17.0472,-17.0472,-17.0596,-17.0597,-17.0646,-17.0646,-17.068,-17.0738,-17.0738,-17.0763,-17.0772,-17.0909,-17.0917,-17.0946,-17.098,-17.1063,-17.108,-17.1105,-17.1138,-17.1155,-17.1163,-17.118,-17.1196,-17.1213,-17.123,-17.1238,-17.1255,-17.1255,-17.1288,-17.1301,-17.1309,-17.1325,-17.1375,-17.1417,-17.1517,-17.1551,-17.1634,-17.1684,-17.1785,-17.1817,-17.1951,-17.1975,-17.2017,-17.2084,-17.2125,-17.2217,-17.2275,-17.2334,-17.235,-17.2417,-17.2434,-17.2459,-17.2476,-17.2509,-17.2526,-17.2551,-17.2567,-17.2592,-17.2609,-17.2742,-17.2759,-17.2859,-17.2875,-17.2951,-17.2967,-17.2992,-17.3009,-17.3018,-17.3067,-17.3092,-17.3109,-17.3134,-17.3251,-17.3268,-17.3284,-17.3309,-17.3326,-17.3376,-17.3392,-17.3509,-17.3525,-17.3568,-17.3592,-17.3626,-17.365,-17.3676,-17.3692,-17.3742,-17.3767,-17.38,-17.3826,-17.3842,-17.3867,-17.39,-17.3951,-17.3967,-17.3992,-17.4034,-17.4042,-17.4067,-17.4076,-17.4109,-17.4117,-17.4142,-17.415,-17.4168,-17.4176,-17.4218,-17.4226,-17.4242,-17.4293,-17.4309,-17.4333,-17.4376,-17.4425,-17.4434,-17.4459,-17.4492,-17.4542,-17.455,-17.4584,-17.4608,-17.4626,-17.4627,-17.4651,-17.4684,-17.4692,-17.4718,-17.4725,-17.4759,-17.4792,-17.4809,-17.4842,-17.4896,-17.49,-17.4934,-17.4951,-17.4967,-17.4992,-17.5009,-17.5017,-17.5051,-17.5068,-17.5101,-17.5142,-17.5184,-17.5196,-17.5209,-17.5234,-17.5243,-17.5259,-17.5263,-17.5292,-17.5301,-17.5313,-17.5313,-17.5313,-17.533,-17.5225,-17.5217,-17.5201,-17.5184,-17.508,-17.508,-17.5071,-17.508,-17.5059,-17.5025,-17.5,-17.4976,-17.4951,-17.493,-17.493,-17.49,-17.488,-17.4888,-17.4875,-17.4859,-17.483,-17.4829,-17.4817,-17.48,-17.478,-17.478,-17.4805,-17.4805,-17.478,-17.4788,-17.4784,-17.4767,-17.4743,-17.473,-17.4735,-17.4738,-17.473,-17.473,-17.4747,-17.4746,-17.4726,-17.4717,-17.4701,-17.4688,-17.4684,-17.4658,-17.4646,-17.4642,-17.4634,-17.4601,-17.4584,-17.4542,-17.4534,-17.4517,-17.4488,-17.4488,-17.4463,-17.4463,-17.4463,-17.4442,-17.4396,-17.4396,-17.4386,-17.438,-17.4358,-17.4334,-17.4321,-17.4315,-17.4313,-17.4313,-17.4313,-17.4338,-17.4338,-17.4318,-17.4267,-17.4251,-17.4221,-17.4205,-17.4205,-17.4226,-17.425,-17.4259,-17.4268,-17.428,-17.4284,-17.43,-17.4317,-17.4325,-17.433,-17.4322,-17.4338,-17.433,-17.4355,-17.4355,-17.4342,-17.4334,-17.4317,-17.4301,-17.4293,-17.4267,-17.4251,-17.4238,-17.4246,-17.4234,-17.4225,-17.4205,-17.4205,-17.4221,-17.4221,-17.4234,-17.425,-17.4271,-17.4271,-17.4263,-17.4242,-17.4217,-17.418,-17.4179,-17.4192,-17.4225,-17.425,-17.4276,-17.4321,-17.433,-17.433,-17.4321,-17.4313,-17.4313,-17.4296,-17.4217,-17.4209,-17.4192,-17.4159,-17.4101,-17.4092,-17.4059,-17.405,-17.4017,-17.4009,-17.3967,-17.3959,-17.3867,-17.3859,-17.3684,-17.3675,-17.3625,-17.3617,-17.3584,-17.3576,-17.3542,-17.3534,-17.3475,-17.3392,-17.3258,-17.3242,-17.3125,-17.3075,-17.3059,-17.3017,-17.3,-17.2992,-17.2975,-17.2942,-17.2901,-17.2859,-17.2826,-17.2784,-17.2759,-17.2742,-17.2726,-17.2692,-17.2676,-17.2659,-17.2642,-17.2617,-17.2601,-17.2525,-17.2484,-17.2467,-17.2351,-17.2284,-17.2259,-17.2217,-17.2208,-17.2192,-17.2084,-17.2076,-17.205,-17.2026,-17.1971,-17.1971,-17.193,-17.193,-17.1867,-17.1859,-17.1755,-17.1755,-17.1696,-17.1663,-17.1663,-17.1613,-17.1605,-17.1588,-17.1571,-17.1546,-17.1521,-17.1521,-17.1505,-17.1496,-17.1496,-17.1471,-17.1459,-17.1446,-17.1446,-17.143,-17.143,-17.1413,-17.1372,-17.1346,-17.1334,-17.1267,-17.1242,-17.1159,-17.1113,-17.1113,-17.1088,-17.1088,-17.1063,-17.1063,-17.1038,-17.1013,-17.1005,-17.098,-17.0971,-17.0946,-17.0946,-17.0938,-17.093,-17.093,-17.0921,-17.0922,-17.0896,-17.0896,-17.0896,-17.088,-17.0805,-17.0805,-17.0747,-17.0755,-17.0726,-17.0712,-17.0709,-17.07,-17.0684,-17.0638,-17.0638,-17.0613,-17.0605,-17.0596,-17.0596,-17.0596,-17.0588,-17.0588,-17.0559,-17.0523,-17.0509,-17.0459,-17.0451,-17.0434,-17.0392,-17.0376,-17.0351,-17.0325,-17.0309,-17.025,-17.0217,-17.0209,-17.0167,-17.0151,-17.0142,-17.0117,-17.0059,-16.9992,-16.9972,-16.9971,-16.9867,-16.9834,-16.9805,-16.9751,-16.9734,-16.9709,-16.9692,-16.9655,-16.9613,-16.9596,-16.9551,-16.9525,-16.9513,-16.9505,-16.9496,-16.9471,-16.9471,-16.9463,-16.9463,-16.9396,-16.9355,-16.9355,-16.9338,-16.933,-16.9313,-16.9313,-16.9321,-16.933,-16.9338,-16.9338,-16.9346,-16.9346,-16.9355,-16.9355,-16.9346,-16.9355,-16.9346,-16.9338,-16.9322,-16.9321,-16.9313,-16.9309,-16.9308,-16.9302,-16.9275,-16.9259,-16.9226,-16.9201,-16.9175,-16.9105,-16.9071,-16.8988,-16.8942,-16.8901,-16.8863,-16.8855,-16.8821,-16.878,-16.8771,-16.875,-16.8738,-16.874,-16.8744,-16.8746,-16.8755,-16.875,-16.8738,-16.873,-16.873,-16.8721,-16.8721,-16.8713,-16.8713,-16.8705,-16.8705,-16.8713,-16.8713,-16.8671,-16.8671,-16.8678,-16.8688,-16.8688,-16.8675,-16.8662,-16.8659,-16.8592,-16.853,-16.8521,-16.8513,-16.8513,-16.8484,-16.8475,-16.8384,-16.8325,-16.8301,-16.8284,-16.8242,-16.8226,-16.8213,-16.8225,-16.8242,-16.8255,-16.8255,-16.8242,-16.8226,-16.8201,-16.8159,-16.815,-16.8134,-16.813,-16.8151,-16.8176,-16.8196,-16.8197,-16.8167,-16.8084,-16.8055,-16.8055,-16.8038,-16.7998,-16.7849,-16.7835,-16.7818,-16.7798,-16.778,-16.7766,-16.7763,-16.7737,-16.7725,-16.771,-16.7693,-16.7677,-16.7655,-16.7641,-16.7633,-16.7625,-16.7615,-16.7593,-16.7585,-16.7573,-16.7559,-16.7548,-16.7537,-16.7533,-16.7526,-16.7527,-16.7533,-16.7535,-16.7551,-16.7559,-16.7565,-16.7574,-16.7597,-16.7617,-16.763,-16.7634,-16.7635,-16.7624,-16.7612,-16.7601,-16.7588,-16.7574,-16.7557,-16.7539,-16.7523,-16.7507,-16.7493,-16.7479,-16.7468,-16.7456,-16.745,-16.7437,-16.7428,-16.7413,-16.7402,-16.7392,-16.7375,-16.7356,-16.7341,-16.7329,-16.7321,-16.7304,-16.7291,-16.7278,-16.7265,-16.7256,-16.7248,-16.7241,-16.7235,-16.7223,-16.7222,-16.7228,-16.7232,-16.7239,-16.7254,-16.7281,-16.7303,-16.7312,-16.7322,-16.7348,-16.7368,-16.7387,-16.7401,-16.7422,-16.7423,-16.7431,-16.7433,-16.7441,-16.7457,-16.746,-16.7476,-16.7474,-16.7486,-16.7488,-16.7488,-16.7484,-16.7475,-16.7472,-16.7463,-16.7452,-16.7437,-16.7425,-16.7409,-16.739,-16.7372,-16.7362,-16.7344,-16.7333,-16.7318,-16.7306,-16.729,-16.7277,-16.7264,-16.7246,-16.7232,-16.7219,-16.7208,-16.7198,-16.7187,-16.7178,-16.7166,-16.715,-16.7138,-16.7124,-16.7109,-16.7096,-16.7088,-16.7085,-16.7078,-16.7073,-16.7071,-16.7075,-16.7074,-16.7079,-16.7078,-16.7079,-16.7089,-16.7089,-16.7092,-16.7097,-16.7102,-16.7103,-16.7098,-16.71,-16.7102,-16.7105,-16.7104,-16.7091,-16.7088,-16.7079,-16.7079,-16.7083,-16.7079,-16.7077,-16.7072,-16.7068,-16.7065,-16.7069,-16.7064,-16.7061,-16.7053,-16.7054,-16.7052,-16.7038,-16.702,-16.7002,-16.6981,-16.6963,-16.6948,-16.6928,-16.6906,-16.6888,-16.6869,-16.685,-16.6832,-16.6811,-16.6795,-16.6777,-16.6755,-16.6733,-16.671,-16.6684,-16.6659,-16.6625,-16.6596,-16.6568,-16.6541,-16.6516,-16.6492,-16.6466,-16.6457,-16.6439,-16.6424,-16.6404,-16.6389,-16.6377,-16.6352,-16.6336,-16.6327,-16.6308,-16.6297,-16.628,-16.6263,-16.6259,-16.625,-16.6234,-16.6217,-16.6202,-16.6175,-16.6155,-16.614,-16.6119,-16.6106,-16.6089,-16.6073,-16.6057,-16.6037,-16.6017,-16.5997,-16.5972,-16.595,-16.592,-16.5906,-16.589,-16.5867,-16.5841,-16.5814,-16.5794,-16.5771,-16.5753,-16.5735,-16.5711,-16.5694,-16.5681,-16.5653,-16.5635,-16.5611,-16.559,-16.5574,-16.5545,-16.5526,-16.5512,-16.5506,-16.5495,-16.5493,-16.5484,-16.5474,-16.5459,-16.5443,-16.5435,-16.5426,-16.5412,-16.5405,-16.5394,-16.5387,-16.5382,-16.537,-16.5361,-16.5361,-16.5358,-16.5348,-16.5349,-16.5346,-16.5335,-16.533,-16.5322,-16.5317,-16.5314,-16.531,-16.5306,-16.5308,-16.5303,-16.5299,-16.5294,-16.5287,-16.5286,-16.5285,-16.5279,-16.5268,-16.5269,-16.5271,-16.5273,-16.5273,-16.5267,-16.5257,-16.5246,-16.5234,-16.5225,-16.5217,-16.5212,-16.5206,-16.5198,-16.519,-16.5186,-16.518,-16.5177,-16.5178,-16.5176,-16.5176,-16.5177,-16.5178,-16.518,-16.5181,-16.5184,-16.5189,-16.5193,-16.5199,-16.5207,-16.5215,-16.5223,-16.5229,-16.5234,-16.5241,-16.5241,-16.5259,-16.5278,-16.5285,-16.5296,-16.5303,-16.5311,-16.5322,-16.5339,-16.5348,-16.5362,-16.5374,-16.5386,-16.5403,-16.5421,-16.5436,-16.5452,-16.5467,-16.5481,-16.55,-16.5516,-16.5533,-16.5546,-16.5562,-16.5568,-16.5591,-16.5608,-16.562,-16.563,-16.5647,-16.5654,-16.5667,-16.5682,-16.5694,-16.5701,-16.5727,-16.5737,-16.575,-16.5768,-16.5776,-16.5791,-16.5802,-16.5812,-16.5819,-16.5829,-16.5845,-16.5849,-16.5862,-16.587,-16.5883,-16.5895,-16.5905,-16.5914,-16.5926,-16.594,-16.5948,-16.5958,-16.5969,-16.5982,-16.5994,-16.6011,-16.6022,-16.6043,-16.6056,-16.6073,-16.609,-16.6102,-16.6118,-16.6134,-16.615,-16.6163,-16.6171,-16.6192,-16.6206,-16.6209,-16.6214,-16.6211,-16.6212,-16.6208,-16.6206,-16.6201,-16.6196,-16.6185,-16.6169,-16.6156,-16.6142,-16.6123,-16.6108,-16.6104,-16.6088,-16.6075,-16.6064,-16.6056,-16.6053,-16.6047,-16.6045,-16.6047,-16.6046,-16.6044,-16.6047,-16.6047,-16.6051,-16.6055,-16.6051,-16.605,-16.6054,-16.6055,-16.6056,-16.6053,-16.6053,-16.605,-16.605,-16.6044,-16.6042,-16.6039,-16.6039,-16.6048,-16.6051,-16.6053,-16.6054,-16.6055,-16.6058,-16.6058,-16.6062,-16.6064,-16.6068,-16.607,-16.6075,-16.608,-16.6093,-16.6109,-16.613,-16.6139,-16.6158,-16.6222,-16.6244,-16.6264,-16.6285,-16.6306,-16.6331,-16.6356,-16.6383,-16.6408,-16.6429,-16.6451,-16.6484,-16.6486,-16.6474,-16.6468,-16.6453,-16.6445,-16.6437,-16.6428,-16.6411,-16.6395,-16.6374,-16.6353,-16.6339,-16.632,-16.6293,-16.6276,-16.6247,-16.6237,-16.6233,-16.6242,-16.6242,-16.6248,-16.6248,-16.6256,-16.6261,-16.6267,-16.6271,-16.6285,-16.6293,-16.6304,-16.631,-16.6319,-16.6331,-16.6346,-16.636,-16.6369,-16.6384,-16.6401,-16.6419,-16.6433,-16.6451,-16.6468,-16.6495,-16.6516,-16.6524,-16.6534,-16.6545,-16.6552,-16.656,-16.6568,-16.6578,-16.6583,-16.6594,-16.6608,-16.6615,-16.6633,-16.6648,-16.6652,-16.6663,-16.6678,-16.6695,-16.6707,-16.6717,-16.6729,-16.6741,-16.6769,-16.6793,-16.6813,-16.6833,-16.6855,-16.6841,-16.6802,-16.6773,-16.6734,-16.6685,-16.6652,-16.6628,-16.6615,-16.658,-16.6542,-16.6522,-16.6509,-16.6479,-16.6444,-16.6402,-16.6349,-16.6288,-16.6211,-16.6143,-16.607,-16.6015,-16.5963,-16.5913,-16.5871,-16.5799,-16.5738,-16.5645,-16.5608,-16.5576,-16.5519,-16.5438,-16.5377,-16.5323,-16.527,-16.5202,-16.5157,-16.5135,-16.5102,-16.5066,-16.5029,-16.5002,-16.4981,-16.4961,-16.4965,-16.4939,-16.4921,-16.4889,-16.4858,-16.4815,-16.4789,-16.4768,-16.4755,-16.4742,-16.4732,-16.4713,-16.4677,-16.4636,-16.4591,-16.4553,-16.4521,-16.4469,-16.4454,-16.4424,-16.4405,-16.438,-16.4353,-16.4334,-16.4314,-16.4291,-16.4263,-16.4232,-16.4198,-16.4169,-16.4138,-16.4107,-16.4084,-16.406,-16.4022,-16.3988,-16.394,-16.3905,-16.3901,-16.3879,-16.3858,-16.3833,-16.3788,-16.3764,-16.3736,-16.3709,-16.3683,-16.3632,-16.3592,-16.3556,-16.3512,-16.3485,-16.3462,-16.3427,-16.3395,-16.3377,-16.3346,-16.3319,-16.3301,-16.328,-16.3254,-16.3225,-16.3198,-16.3188,-16.313,-16.3099,-16.3081,-16.3049,-16.3029,-16.3008,-16.2988,-16.2976,-16.2939,-16.2911,-16.2885,-16.2867,-16.2841,-16.281,-16.2786,-16.2769,-16.2746,-16.2724,-16.2708,-16.2683,-16.2654,-16.2626,-16.2601,-16.2579,-16.255,-16.253,-16.2502,-16.2483,-16.2472,-16.2454,-16.2426,-16.2408,-16.2396,-16.237,-16.2357,-16.2345,-16.2329,-16.232,-16.2322,-16.2316,-16.2311,-16.2308,-16.2294,-16.2275,-16.2253,-16.2237,-16.2206,-16.2183,-16.2153,-16.2123,-16.2101,-16.2068,-16.2036,-16.2014,-16.1988,-16.1966,-16.1945,-16.1932,-16.1904,-16.188,-16.1854,-16.1817,-16.1782,-16.177,-16.177,-16.1771,-16.1772,-16.1771,-16.1771,-16.1766,-16.1761,-16.1755,-16.1759,-16.1764,-16.1767,-16.1771,-16.1776,-16.1779,-16.1781,-16.1782,-16.1787,-16.1783,-16.1782,-16.1784,-16.1789,-16.1796,-16.18,-16.1808,-16.1824,-16.1832,-16.1858,-16.1873,-16.1886,-16.1891,-16.1904,-16.1921,-16.1931,-16.1944,-16.1943,-16.1942,-16.1947,-16.1947,-16.1947,-16.1947,-16.1947,-16.1951,-16.1962,-16.1968,-16.1989,-16.2012,-16.2041,-16.205,-16.2086,-16.2116,-16.2136,-16.2188,-16.2249,-16.227,-16.2304,-16.2321,-16.2345,-16.2364,-16.2371,-16.2374,-16.2377,-16.2373,-16.2372,-16.2373,-16.2376,-16.2392,-16.2416,-16.2426,-16.2438,-16.245,-16.2475,-16.2499,-16.2526,-16.2551,-16.2572,-16.2587,-16.2601,-16.2616,-16.2632,-16.2647,-16.2661,-16.2672,-16.2683,-16.2688,-16.2693,-16.2699,-16.2697,-16.2692,-16.2689,-16.2685,-16.2681,-16.268,-16.2685,-16.2692,-16.27,-16.2708,-16.2723,-16.2742,-16.276,-16.278,-16.2806,-16.2832,-16.2849,-16.2858,-16.2868,-16.2885,-16.2908,-16.2931,-16.295,-16.2969,-16.3,-16.3017,-16.3045,-16.3064,-16.308,-16.3084,-16.3088,-16.3093,-16.3102,-16.3112,-16.3142,-16.3171,-16.3202,-16.3237,-16.3274,-16.3305,-16.3332,-16.335,-16.338,-16.3415,-16.3436,-16.345,-16.3457,-16.3454,-16.3446,-16.3416,-16.3388,-16.3371,-16.3354,-16.3344,-16.3346,-16.3356,-16.3371,-16.3393,-16.3411,-16.3421,-16.3448,-16.3468,-16.3496,-16.3515,-16.3547,-16.3587,-16.363,-16.366,-16.3693,-16.3719,-16.3732,-16.3752,-16.3773,-16.3796,-16.3839,-16.3878,-16.3932,-16.3979,-16.4039,-16.4107,-16.4153,-16.4189,-16.4228,-16.4255,-16.428,-16.4311,-16.4338,-16.4362,-16.4385,-16.4414,-16.4436,-16.4458,-16.4483,-16.4509,-16.4515,-16.4526,-16.4538,-16.4546,-16.4564,-16.4594,-16.4618,-16.4641,-16.4659,-16.4691,-16.4709,-16.4736,-16.4751,-16.4759,-16.4779,-16.4791,-16.4804,-16.4817,-16.4832,-16.4845,-16.4861,-16.4873,-16.4889,-16.4901,-16.4915,-16.4941,-16.4954,-16.4974,-16.499,-16.5007,-16.502,-16.503,-16.5044,-16.506,-16.5077,-16.5096,-16.5111,-16.5124,-16.5137,-16.5146,-16.5158,-16.5162,-16.5164,-16.5156,-16.5143,-16.5135,-16.5123,-16.5123,-16.5124,-16.5129,-16.5132,-16.5135,-16.5145,-16.5158,-16.5173,-16.5186,-16.5202,-16.5225,-16.5237,-16.5251,-16.5279,-16.5307,-16.5332,-16.535,-16.5381,-16.54,-16.5418,-16.5451,-16.5471,-16.5519,-16.554,-16.5571,-16.5596,-16.5647,-16.5671,-16.5696,-16.5722,-16.581,-16.5824,-16.5836,-16.5907,-16.5958,-16.5993,-16.6012,-16.6037,-16.6045,-16.6072,-16.6092,-16.6128,-16.6141,-16.614,-16.6134,-16.6129,-16.6123,-16.6121,-16.6121,-16.6133,-16.6165,-16.6207,-16.6245,-16.6285,-16.6306,-16.6307,-16.6312,-16.6338,-16.6365,-16.6384,-16.6412,-16.6423,-16.6428,-16.6476,-16.6536,-16.6557,-16.6608,-16.6626,-16.6635,-16.6669,-16.6687,-16.6701,-16.673,-16.6745,-16.676,-16.6811,-16.6828,-16.6844,-16.6867,-16.6907,-16.6922,-16.6961,-16.701,-16.7022,-16.7048,-16.7073,-16.7081,-16.7094,-16.7136,-16.7155,-16.7167,-16.718,-16.7233,-16.7305,-16.7319,-16.7333,-16.7366],"lat":[15.2967,15.297,15.2991,15.2994,15.3026,15.31,15.3152,15.3166,15.3241,15.3313,15.3395,15.3406,15.3523,15.3549,15.3577,15.364,15.3645,15.3693,15.3726,15.3804,15.3866,15.3905,15.3963,15.4021,15.403,15.4068,15.4086,15.4167,15.422,15.4226,15.4259,15.4278,15.4303,15.4333,15.4348,15.4379,15.443,15.4471,15.4525,15.4614,15.4596,15.4554,15.4512,15.4462,15.4429,15.4387,15.4379,15.4354,15.4296,15.4254,15.4238,15.4163,15.4129,15.4112,15.4088,15.4071,15.4046,15.4029,15.3979,15.392,15.3904,15.3888,15.3871,15.3854,15.3804,15.3771,15.3737,15.3712,15.3688,15.3638,15.3504,15.3496,15.3421,15.3396,15.3371,15.3238,15.3195,15.3138,15.3104,15.3046,15.3012,15.2971,15.2938,15.2913,15.2904,15.2879,15.2854,15.2837,15.2771,15.2754,15.2738,15.2713,15.2688,15.2646,15.2579,15.2546,15.2463,15.2454,15.2429,15.2421,15.2371,15.2346,15.2321,15.2295,15.2204,15.2196,15.2163,15.2137,15.2079,15.2071,15.1954,15.1946,15.1904,15.1888,15.1837,15.1771,15.1754,15.1696,15.1621,15.1596,15.1504,15.1471,15.1371,15.1363,15.1279,15.1262,15.1254,15.1129,15.1062,15.1054,15.1021,15.1012,15.0912,15.0904,15.0846,15.0821,15.0784,15.0712,15.0696,15.0612,15.0571,15.0487,15.0479,15.037,15.0363,15.0313,15.0304,15.0246,15.0238,15.0105,15.0096,15.0046,15.0037,14.9996,14.9946,14.9937,14.9912,14.9888,14.975,14.975,14.9721,14.9654,14.9571,14.9538,14.9512,14.9454,14.9446,14.9421,14.9404,14.9379,14.9313,14.9279,14.9237,14.9204,14.9187,14.9145,14.9125,14.9125,14.9109,14.9091,14.9042,14.8992,14.8983,14.8925,14.89,14.885,14.8834,14.8767,14.875,14.8742,14.8708,14.87,14.865,14.8608,14.8583,14.8567,14.8534,14.8533,14.8516,14.8517,14.85,14.85,14.8484,14.8483,14.8467,14.8467,14.8391,14.8392,14.8342,14.8342,14.83,14.83,14.8283,14.8283,14.8275,14.8267,14.825,14.825,14.8234,14.8175,14.8159,14.8158,14.8141,14.8142,14.8117,14.8116,14.8059,14.8058,14.8033,14.8034,14.8017,14.8017,14.8,14.8,14.7975,14.795,14.7933,14.7908,14.7908,14.7892,14.7884,14.7858,14.7859,14.7842,14.7833,14.7825,14.7825,14.7817,14.7816,14.7808,14.7809,14.78,14.78,14.7792,14.7784,14.7775,14.7775,14.775,14.775,14.7725,14.7709,14.7708,14.77,14.77,14.7683,14.7675,14.7667,14.7658,14.7642,14.7642,14.7641,14.7625,14.7625,14.7633,14.7633,14.7642,14.7642,14.7658,14.7675,14.7675,14.762,14.76,14.7592,14.7567,14.7575,14.7575,14.7592,14.7592,14.7559,14.7566,14.7567,14.7517,14.7525,14.7512,14.7492,14.7492,14.7501,14.7517,14.7517,14.7516,14.7496,14.7471,14.7469,14.7437,14.7413,14.7408,14.74,14.74,14.7392,14.7287,14.7271,14.7262,14.7204,14.7183,14.7183,14.715,14.715,14.7167,14.7138,14.7112,14.7075,14.7062,14.7054,14.7042,14.7042,14.702,14.6996,14.6984,14.6983,14.6971,14.6963,14.6938,14.6921,14.6904,14.6888,14.6875,14.6866,14.6867,14.6854,14.6834,14.6821,14.6813,14.6796,14.6779,14.6763,14.6742,14.6742,14.6759,14.6754,14.6742,14.6742,14.6761,14.6767,14.6767,14.6733,14.6725,14.6725,14.6717,14.6717,14.6687,14.6671,14.6637,14.6609,14.6587,14.6575,14.6538,14.6521,14.65,14.6488,14.6467,14.6467,14.6479,14.6492,14.6496,14.6501,14.6537,14.6587,14.6629,14.665,14.6675,14.6675,14.6696,14.6737,14.6771,14.6775,14.6758,14.6775,14.6767,14.6771,14.6792,14.6792,14.6775,14.6775,14.6796,14.6813,14.6821,14.6846,14.6854,14.6871,14.6884,14.6875,14.6883,14.6875,14.6892,14.685,14.6858,14.6846,14.6813,14.68,14.68,14.6821,14.6846,14.6871,14.6888,14.69,14.69,14.6921,14.6938,14.6954,14.6967,14.7,14.7012,14.7054,14.7067,14.7067,14.7092,14.7067,14.7113,14.7129,14.7162,14.7179,14.7187,14.7204,14.7237,14.7309,14.7308,14.7325,14.7342,14.735,14.7358,14.7359,14.7367,14.7367,14.7375,14.7375,14.7384,14.7392,14.74,14.74,14.7392,14.7392,14.7383,14.7383,14.7375,14.7375,14.7367,14.7367,14.7325,14.7258,14.7259,14.72,14.7175,14.7175,14.715,14.715,14.7142,14.7142,14.7125,14.7092,14.7092,14.7117,14.7117,14.71,14.71,14.7084,14.7067,14.7067,14.705,14.705,14.7034,14.7033,14.6992,14.6984,14.6975,14.6867,14.6833,14.6809,14.6792,14.6783,14.6784,14.6675,14.6675,14.665,14.6642,14.6588,14.6579,14.6537,14.6529,14.6467,14.6467,14.6362,14.6354,14.6296,14.6246,14.6238,14.6188,14.6162,14.6146,14.6113,14.6088,14.6029,14.6004,14.5979,14.5963,14.5913,14.5871,14.5867,14.5846,14.5829,14.5804,14.5788,14.5754,14.5712,14.5654,14.5642,14.5608,14.5583,14.5542,14.5496,14.5487,14.5446,14.5429,14.5396,14.5379,14.5354,14.5304,14.5237,14.5204,14.5171,14.5137,14.5121,14.5113,14.5096,14.5071,14.5062,14.5029,14.4979,14.4969,14.4945,14.4912,14.4846,14.4813,14.4763,14.4738,14.4725,14.4725,14.4725,14.4717,14.4717,14.4679,14.4663,14.4629,14.4588,14.4579,14.4571,14.4554,14.4546,14.4521,14.4492,14.4474,14.4467,14.4467,14.4475,14.4475,14.4459,14.4442,14.4442,14.4458,14.4458,14.4392,14.4392,14.44,14.4358,14.4358,14.4367,14.4366,14.4333,14.4259,14.4237,14.4229,14.4134,14.4142,14.4096,14.4041,14.4034,14.4033,14.4025,14.3988,14.3904,14.387,14.3825,14.3825,14.3813,14.3779,14.3771,14.3729,14.3704,14.3696,14.3637,14.3571,14.3504,14.3495,14.348,14.3437,14.3413,14.3313,14.3304,14.3204,14.3196,14.3163,14.3154,14.3079,14.3071,14.3038,14.3029,14.3013,14.3004,14.2962,14.2937,14.2904,14.2896,14.2859,14.285,14.2848,14.2842,14.2859,14.2859,14.2833,14.2834,14.2762,14.2704,14.2588,14.2542,14.2517,14.2479,14.2454,14.2413,14.2321,14.2279,14.2259,14.2263,14.2273,14.2306,14.2321,14.2329,14.2358,14.2312,14.2304,14.2271,14.2263,14.2221,14.2212,14.218,14.2171,14.2087,14.2079,14.2013,14.1946,14.1879,14.1873,14.1862,14.1846,14.1842,14.1848,14.185,14.1817,14.1763,14.1713,14.1704,14.1671,14.1641,14.1641,14.1575,14.1542,14.1542,14.155,14.1592,14.1592,14.1571,14.1558,14.1558,14.1546,14.1521,14.1508,14.1509,14.1492,14.1492,14.15,14.15,14.1487,14.1475,14.1475,14.1454,14.1429,14.14,14.1367,14.1346,14.1338,14.1321,14.124,14.1397,14.14,14.1403,14.1411,14.1414,14.1403,14.1394,14.1387,14.138,14.1379,14.138,14.1385,14.1395,14.1403,14.1412,14.142,14.1424,14.142,14.1412,14.1407,14.1401,14.1403,14.141,14.1415,14.1432,14.1448,14.1458,14.1468,14.1483,14.1502,14.1518,14.1532,14.1546,14.1562,14.1576,14.1594,14.1609,14.1625,14.1635,14.1655,14.1674,14.1691,14.1698,14.1705,14.1708,14.171,14.1715,14.1729,14.1739,14.1754,14.1766,14.1789,14.1806,14.182,14.1825,14.1827,14.1827,14.1823,14.1829,14.1848,14.1864,14.1879,14.1889,14.1895,14.1915,14.1932,14.1952,14.1971,14.1988,14.2011,14.2035,14.2048,14.2064,14.208,14.2091,14.2104,14.2113,14.2111,14.2112,14.2126,14.2144,14.2156,14.2173,14.2189,14.2201,14.2209,14.2222,14.2231,14.2252,14.2273,14.2294,14.2317,14.2348,14.2376,14.2407,14.2432,14.2454,14.2478,14.2504,14.2524,14.2539,14.2552,14.2574,14.2599,14.2623,14.2638,14.2658,14.2669,14.2684,14.2699,14.2713,14.2731,14.2752,14.2772,14.2791,14.2807,14.2823,14.284,14.2861,14.2876,14.2894,14.2918,14.294,14.2963,14.2984,14.3004,14.3024,14.3052,14.3069,14.3097,14.3117,14.314,14.3166,14.3186,14.321,14.3237,14.3262,14.3287,14.3306,14.3332,14.3363,14.339,14.3419,14.3449,14.3473,14.3499,14.3523,14.3551,14.3584,14.3608,14.363,14.3647,14.368,14.3707,14.373,14.3752,14.3778,14.3805,14.3835,14.3859,14.3891,14.3918,14.395,14.398,14.4004,14.4028,14.4046,14.4067,14.4081,14.4098,14.4121,14.4142,14.416,14.418,14.4199,14.4222,14.423,14.4239,14.4251,14.4259,14.4256,14.4259,14.4257,14.4257,14.426,14.4261,14.4258,14.4259,14.4266,14.4262,14.4263,14.4266,14.4267,14.4271,14.4275,14.428,14.4287,14.43,14.4314,14.4335,14.435,14.4366,14.4383,14.4391,14.4406,14.4421,14.4437,14.4451,14.4472,14.448,14.449,14.4501,14.4506,14.4515,14.4523,14.4525,14.4529,14.4536,14.4545,14.4554,14.4561,14.4586,14.4601,14.4611,14.4624,14.4646,14.4663,14.4677,14.4698,14.4712,14.4735,14.475,14.4765,14.4776,14.4797,14.4814,14.4834,14.4855,14.4871,14.4891,14.4909,14.4921,14.494,14.4955,14.4963,14.4992,14.5001,14.5015,14.503,14.5044,14.5061,14.508,14.5096,14.5111,14.5132,14.5151,14.517,14.5184,14.5197,14.5203,14.5224,14.5233,14.5252,14.5276,14.529,14.5303,14.5316,14.5328,14.535,14.5373,14.5386,14.541,14.5432,14.5455,14.5466,14.5484,14.5503,14.5521,14.5543,14.5575,14.5585,14.5599,14.5615,14.5624,14.5638,14.5655,14.5677,14.5691,14.5704,14.5715,14.5728,14.5741,14.5753,14.5766,14.5784,14.5795,14.5811,14.5831,14.5845,14.5864,14.588,14.5892,14.5904,14.5917,14.5934,14.5948,14.5964,14.5987,14.6006,14.6024,14.6037,14.6043,14.6052,14.606,14.6057,14.606,14.606,14.606,14.6061,14.6061,14.606,14.6059,14.6058,14.6057,14.6057,14.6059,14.6059,14.6056,14.6053,14.6054,14.6051,14.6047,14.6046,14.6046,14.6045,14.6046,14.6043,14.6042,14.6038,14.6042,14.6042,14.604,14.6038,14.6036,14.6038,14.604,14.6044,14.6044,14.605,14.6053,14.6057,14.606,14.6067,14.6071,14.6072,14.6072,14.6079,14.6089,14.6093,14.6099,14.6117,14.6122,14.613,14.6142,14.615,14.6161,14.6168,14.6176,14.6182,14.6195,14.6201,14.6217,14.6231,14.6252,14.6258,14.6273,14.6281,14.6302,14.6314,14.6322,14.6332,14.6341,14.6352,14.6363,14.6371,14.6382,14.6422,14.6441,14.646,14.6487,14.6513,14.6536,14.6564,14.659,14.6619,14.6646,14.6669,14.6682,14.67,14.6717,14.6737,14.6744,14.6777,14.6802,14.6831,14.6849,14.6865,14.6898,14.6921,14.6945,14.6964,14.6997,14.7019,14.704,14.7061,14.7091,14.7125,14.7149,14.7175,14.7193,14.7207,14.7228,14.7258,14.7282,14.7312,14.7339,14.7363,14.7395,14.7408,14.7429,14.7452,14.747,14.749,14.7519,14.7543,14.7558,14.758,14.7616,14.764,14.7661,14.7682,14.7699,14.773,14.7749,14.7768,14.7777,14.7793,14.7925,14.7936,14.7946,14.7961,14.7978,14.7991,14.8005,14.8027,14.8044,14.8057,14.8073,14.812,14.8148,14.8173,14.8187,14.8206,14.822,14.823,14.8247,14.8264,14.8281,14.83,14.8319,14.8329,14.8344,14.8363,14.838,14.8398,14.8411,14.8449,14.8468,14.8479,14.8491,14.8512,14.8529,14.8553,14.8591,14.8612,14.8632,14.8651,14.8664,14.8676,14.8688,14.8701,14.8718,14.8729,14.8742,14.8759,14.8772,14.8784,14.8796,14.8805,14.8822,14.8844,14.8858,14.8865,14.8874,14.8887,14.8905,14.8922,14.8938,14.895,14.8963,14.8988,14.9022,14.905,14.9053,14.9068,14.9082,14.9093,14.9108,14.9128,14.9143,14.9155,14.9166,14.918,14.9199,14.9219,14.9235,14.9253,14.9273,14.9292,14.9335,14.9362,14.9402,14.9451,14.948,14.9497,14.9506,14.9532,14.9553,14.9566,14.9574,14.9586,14.96,14.9618,14.9626,14.9642,14.9657,14.9673,14.9688,14.9704,14.9708,14.9724,14.9737,14.9757,14.9777,14.9796,14.9811,14.9824,14.9842,14.9862,14.9883,14.9897,14.9918,14.9945,14.9955,14.9967,14.997,14.9976,14.9978,14.9982,14.9984,14.9993,15.0008,15.0026,15.0036,15.0061,15.0084,15.0109,15.0125,15.0136,15.0143,15.0148,15.0153,15.0162,15.0176,15.0194,15.0213,15.022,15.0222,15.0235,15.0238,15.0243,15.0243,15.0243,15.0244,15.0245,15.0246,15.0246,15.0244,15.0242,15.0239,15.0236,15.0234,15.023,15.0228,15.0226,15.0222,15.0218,15.0221,15.0219,15.0212,15.021,15.021,15.0207,15.0199,15.0192,15.0189,15.018,15.0171,15.0165,15.0154,15.0144,15.0137,15.0127,15.0119,15.0107,15.0094,15.0086,15.0067,15.0057,15.0048,15.0039,15.0027,15.0018,15.001,15.0009,14.9982,14.9975,14.9969,14.9963,14.9957,14.9954,14.9951,14.9947,14.9938,14.9934,14.9925,14.9921,14.9913,14.9902,14.9903,14.9898,14.9888,14.9879,14.987,14.986,14.985,14.9844,14.9836,14.982,14.9806,14.9799,14.9783,14.9771,14.9761,14.9751,14.9732,14.9722,14.971,14.9694,14.9677,14.9668,14.965,14.9633,14.962,14.9598,14.9575,14.956,14.9546,14.953,14.9528,14.9525,14.9519,14.952,14.9527,14.9535,14.9534,14.9544,14.9558,14.9562,14.9574,14.9579,14.9589,14.9592,14.9603,14.961,14.962,14.9625,14.9628,14.9629,14.9648,14.9662,14.9682,14.9695,14.9712,14.9732,14.9752,14.9767,14.9784,14.9803,14.9827,14.9844,14.9865,14.9887,14.9908,14.9934,14.9953,15.0006,15.0015,15.0058,15.0078,15.0092,15.0115,15.0118,15.0161,15.0182,15.0219,15.024,15.0273,15.0306,15.0337,15.037,15.0417,15.0454,15.0491,15.0561,15.0617,15.0663,15.0711,15.0766,15.0805,15.0863,15.0912,15.0934,15.0954,15.0971,15.1003,15.1014,15.1046,15.1068,15.1079,15.1113,15.1147,15.1162,15.118,15.1194,15.1219,15.1247,15.1278,15.131,15.1345,15.1382,15.1443,15.15,15.1537,15.1563,15.1607,15.163,15.1658,15.1677,15.1692,15.1712,15.1741,15.1771,15.1807,15.1831,15.1855,15.1883,15.1907,15.1941,15.1971,15.2006,15.2049,15.21,15.2149,15.2176,15.2201,15.2232,15.2247,15.2288,15.2313,15.2357,15.2392,15.2425,15.2467,15.2498,15.2536,15.257,15.2596,15.2618,15.2646,15.2672,15.269,15.2696,15.2703,15.2717,15.2722,15.2725,15.273,15.2728,15.2725,15.2723,15.2706,15.2691,15.2662,15.2634,15.2596,15.2548,15.2521,15.2509,15.2488,15.249,15.25,15.2513,15.2537,15.2561,15.2579,15.2609,15.2636,15.266,15.2689,15.2715,15.2743,15.2758,15.2769,15.2793,15.282,15.2847,15.2878,15.2925,15.2957,15.2994,15.3029,15.3056,15.3088,15.3105,15.3136,15.3167,15.3196,15.3213,15.3235,15.3242,15.3239,15.3233,15.3213,15.3192,15.3174,15.3164,15.3171,15.3188,15.3212,15.3233,15.3245,15.3255,15.3257,15.3252,15.3258,15.3276,15.3311,15.3339,15.3368,15.3406,15.3431,15.3455,15.3482,15.3511,15.3536,15.356,15.3577,15.3591,15.3591,15.3582,15.3569,15.356,15.3545,15.3518,15.3499,15.3476,15.3461,15.3434,15.3417,15.3393,15.3378,15.3361,15.3341,15.3332,15.3314,15.3296,15.3281,15.3267,15.325,15.3237,15.3224,15.3208,15.3186,15.3158,15.3147,15.3133,15.3119,15.31,15.3083,15.3066,15.3045,15.3022,15.3001,15.2974,15.2953,15.2935,15.2917,15.29,15.288,15.2868,15.2861,15.2849,15.2844,15.2839,15.2836,15.2828,15.281,15.2795,15.2783,15.2766,15.2743,15.2726,15.2713,15.2703,15.2696,15.2688,15.2686,15.2687,15.2693,15.2662,15.2619,15.2585,15.2524,15.2486,15.2456,15.2389,15.2355,15.2254,15.2224,15.2141,15.2111,15.1987,15.1959,15.1924,15.1886,15.1797,15.1792,15.1783,15.1714,15.1649,15.1573,15.1548,15.154,15.1532,15.1526,15.1496,15.1457,15.1563,15.1621,15.166,15.1692,15.1746,15.1785,15.1824,15.1891,15.1957,15.1985,15.2002,15.2034,15.2086,15.2115,15.2178,15.2226,15.2247,15.226,15.2271,15.2275,15.2283,15.2304,15.2357,15.2367,15.2405,15.2409,15.2426,15.2453,15.2461,15.2486,15.2522,15.2541,15.2561,15.2603,15.2625,15.2629,15.265,15.2667,15.2685,15.2693,15.2722,15.2737,15.2768,15.2796,15.282,15.2835,15.2892,15.2903,15.2904,15.291,15.2917,15.2939,15.2949,15.2948,15.2967]}]],[[{"lng":[-16.5668,-16.5692,-16.5696,-16.5684,-16.5675,-16.5663,-16.5668],"lat":[13.6258,13.6258,13.6263,13.6275,13.6275,13.6262,13.6258]},{"lng":[-16.6396,-16.6417,-16.643,-16.643,-16.6409,-16.6359,-16.6346,-16.6346,-16.6359,-16.6367,-16.6392,-16.6396],"lat":[13.6571,13.6575,13.6562,13.6537,13.6517,13.6525,13.6537,13.6563,13.6567,13.655,13.655,13.6571]},{"lng":[-16.668,-16.6679,-16.6688,-16.6688,-16.6701,-16.6709,-16.6719,-16.6721,-16.6721,-16.6713,-16.6713,-16.6642,-16.6592,-16.6542,-16.6484,-16.6442,-16.6436,-16.6421,-16.6421,-16.6434,-16.645,-16.6484,-16.6509,-16.653,-16.6538,-16.6567,-16.6609,-16.6622,-16.6621,-16.6634,-16.6667,-16.668],"lat":[13.6554,13.6596,13.6604,13.6621,13.6633,13.6634,13.6623,13.662,13.6554,13.6546,13.6529,13.6467,13.645,13.6384,13.6359,13.6358,13.6365,13.6379,13.6404,13.6417,13.6417,13.6391,13.6392,13.6413,13.6437,13.6467,13.6475,13.6487,13.6529,13.6542,13.6542,13.6554]},{"lng":[-16.6546,-16.6559,-16.6567,-16.6608,-16.6621,-16.6576,-16.6546],"lat":[13.678,13.6792,13.6783,13.6783,13.6771,13.6767,13.678]},{"lng":[-16.6346,-16.6346,-16.6313,-16.6313,-16.6321,-16.6355,-16.6355,-16.6355,-16.6355,-16.6359,-16.6371,-16.6371,-16.638,-16.638,-16.6388,-16.6388,-16.6396,-16.6405,-16.6413,-16.6413,-16.6401,-16.6384,-16.6346],"lat":[13.6629,13.6663,13.6695,13.6721,13.6729,13.6779,13.6835,13.6861,13.6912,13.6917,13.6913,13.6837,13.6829,13.6779,13.6771,13.6721,13.6712,13.6654,13.6646,13.6612,13.66,13.66,13.6629]},{"lng":[-16.6759,-16.6722,-16.6734,-16.6768,-16.6792,-16.6808,-16.6813,-16.6784,-16.6759],"lat":[13.6934,13.6937,13.695,13.695,13.6934,13.6933,13.6921,13.6917,13.6934]},{"lng":[-16.6526,-16.6517,-16.6463,-16.6455,-16.6455,-16.6475,-16.6484,-16.6526,-16.6538,-16.653,-16.6538,-16.6549,-16.6588,-16.6588,-16.6605,-16.6605,-16.6584,-16.6526],"lat":[13.7775,13.7775,13.7829,13.7846,13.7855,13.7875,13.7875,13.7839,13.7829,13.7821,13.7804,13.7793,13.7754,13.7746,13.773,13.7721,13.7717,13.7775]},{"lng":[-16.6609,-16.6563,-16.6547,-16.6521,-16.6521,-16.6546,-16.6546,-16.6567,-16.658,-16.658,-16.6588,-16.6596,-16.6621,-16.6621,-16.6638,-16.6638,-16.6646,-16.6634,-16.6609],"lat":[13.7758,13.7804,13.7846,13.7871,13.7887,13.7904,13.7912,13.7934,13.7921,13.7904,13.7895,13.7845,13.7804,13.7795,13.7779,13.7771,13.7763,13.775,13.7758]},{"lng":[-16.6696,-16.6709,-16.6733,-16.6746,-16.6746,-16.6746,-16.6734,-16.6696],"lat":[13.7979,13.8025,13.8017,13.8004,13.7954,13.7945,13.7942,13.7979]},{"lng":[-16.65,-16.6472,-16.643,-16.6421,-16.643,-16.643,-16.645,-16.6463,-16.6555,-16.6555,-16.6563,-16.6525,-16.65],"lat":[13.7925,13.7962,13.8046,13.8055,13.807,13.8096,13.8108,13.8071,13.7979,13.7971,13.7962,13.7917,13.7925]},{"lng":[-16.638,-16.6363,-16.6346,-16.635,-16.6367,-16.6409,-16.6422,-16.6405,-16.6404,-16.6415,-16.6421,-16.6409,-16.638],"lat":[13.8054,13.8088,13.8104,13.8142,13.8134,13.8133,13.812,13.8087,13.8063,13.8041,13.8029,13.8025,13.8054]},{"lng":[-16.6117,-16.6076,-16.6051,-16.5992,-16.5967,-16.5913,-16.5905,-16.5905,-16.5951,-16.6051,-16.6084,-16.6109,-16.6134,-16.6176,-16.6209,-16.6234,-16.6276,-16.6313,-16.6321,-16.6338,-16.6363,-16.6384,-16.6418,-16.6446,-16.6446,-16.6463,-16.6463,-16.6426,-16.64,-16.6392,-16.6351,-16.6343,-16.6292,-16.6284,-16.6267,-16.6242,-16.6209,-16.6151,-16.6142,-16.6117],"lat":[13.8,13.8025,13.8034,13.81,13.8108,13.8171,13.8188,13.8229,13.8283,13.8283,13.8267,13.8267,13.825,13.8242,13.8225,13.8225,13.82,13.8154,13.8121,13.8096,13.8038,13.8017,13.8009,13.797,13.7946,13.7929,13.7888,13.785,13.785,13.7859,13.7867,13.7875,13.7875,13.7883,13.7884,13.7908,13.7925,13.7975,13.7975,13.8]},{"lng":[-16.7346,-16.7363,-16.7363,-16.7371,-16.7371,-16.7371,-16.7384,-16.7405,-16.7405,-16.7396,-16.738,-16.7351,-16.7346],"lat":[13.8413,13.8429,13.8462,13.8471,13.85,13.8554,13.8567,13.8554,13.8504,13.8496,13.8429,13.84,13.8413]},{"lng":[-16.7405,-16.7405,-16.7409,-16.7426,-16.743,-16.743,-16.7425,-16.7417,-16.7405],"lat":[13.9062,13.9071,13.9075,13.9075,13.9071,13.9054,13.905,13.905,13.9062]},{"lng":[-16.7463,-16.7463,-16.7476,-16.7497,-16.7496,-16.7505,-16.7505,-16.7492,-16.7463],"lat":[13.9196,13.9229,13.9242,13.9229,13.9187,13.9179,13.9154,13.915,13.9196]},{"lng":[-15.8028,-15.8034,-15.8104,-15.8161,-15.8201,-15.8234,-15.8264,-15.8305,-15.8342,-15.8385,-15.842,-15.8467,-15.8515,-15.855,-15.8588,-15.862,-15.8646,-15.8651,-15.8644,-15.8633,-15.8626,-15.8623,-15.8626,-15.8635,-15.8686,-15.877,-15.8837,-15.8915,-15.8942,-15.8962,-15.8976,-15.8992,-15.9009,-15.9026,-15.9048,-15.9075,-15.9103,-15.9133,-15.9147,-15.9157,-15.9163,-15.9171,-15.9196,-15.9224,-15.9244,-15.9276,-15.9331,-15.9369,-15.9434,-15.9479,-15.9532,-15.9575,-15.963,-15.9688,-15.972,-15.9744,-15.9777,-15.9811,-15.9854,-15.9873,-15.9884,-15.9891,-15.9899,-15.9904,-15.9911,-15.9916,-15.9922,-15.9929,-15.9936,-15.9944,-15.9954,-15.9962,-15.9968,-15.9973,-15.999,-16.0001,-16.0018,-16.004,-16.0058,-16.0078,-16.0099,-16.0139,-16.0156,-16.0176,-16.0264,-16.0321,-16.0348,-16.0361,-16.0377,-16.0408,-16.0439,-16.0469,-16.0495,-16.0521,-16.0545,-16.0577,-16.0602,-16.0658,-16.0681,-16.0703,-16.0728,-16.0759,-16.0794,-16.0816,-16.0843,-16.0863,-16.0879,-16.0905,-16.0921,-16.0938,-16.0948,-16.0971,-16.0992,-16.1025,-16.1049,-16.1074,-16.1095,-16.1114,-16.1158,-16.1177,-16.1213,-16.1235,-16.1284,-16.1317,-16.1345,-16.1368,-16.1395,-16.1422,-16.147,-16.1486,-16.1522,-16.1544,-16.1577,-16.1605,-16.1638,-16.1674,-16.1705,-16.1735,-16.176,-16.177,-16.1782,-16.1817,-16.1854,-16.188,-16.1904,-16.1932,-16.1945,-16.1966,-16.1988,-16.2014,-16.2036,-16.2068,-16.2101,-16.2123,-16.2153,-16.2183,-16.2206,-16.2237,-16.2253,-16.2275,-16.2294,-16.2308,-16.2311,-16.2316,-16.2322,-16.232,-16.2329,-16.2345,-16.2357,-16.237,-16.2396,-16.2408,-16.2426,-16.2454,-16.2472,-16.2483,-16.2502,-16.253,-16.255,-16.2579,-16.2601,-16.2626,-16.2654,-16.2683,-16.2708,-16.2724,-16.2746,-16.2769,-16.2786,-16.281,-16.2841,-16.2867,-16.2885,-16.2911,-16.2939,-16.2976,-16.2988,-16.3008,-16.3029,-16.3049,-16.3081,-16.3099,-16.313,-16.3188,-16.3198,-16.3225,-16.3254,-16.328,-16.3301,-16.3319,-16.3346,-16.3377,-16.3395,-16.3427,-16.3462,-16.3485,-16.3512,-16.3556,-16.3592,-16.3632,-16.3683,-16.3709,-16.3736,-16.3764,-16.3788,-16.3833,-16.3858,-16.3879,-16.3901,-16.3905,-16.394,-16.3988,-16.4022,-16.406,-16.4084,-16.4107,-16.4138,-16.4169,-16.4198,-16.4232,-16.4263,-16.4291,-16.4314,-16.4334,-16.4353,-16.438,-16.4405,-16.4424,-16.4454,-16.4469,-16.4521,-16.4553,-16.4591,-16.4636,-16.4677,-16.4713,-16.4732,-16.4742,-16.4755,-16.4768,-16.4789,-16.4815,-16.4858,-16.4889,-16.4921,-16.4939,-16.4965,-16.4961,-16.4981,-16.5002,-16.5029,-16.5066,-16.5102,-16.5135,-16.5157,-16.5202,-16.527,-16.5323,-16.5377,-16.5438,-16.5519,-16.5576,-16.5608,-16.5645,-16.5738,-16.5799,-16.5871,-16.5913,-16.5963,-16.6015,-16.607,-16.6143,-16.6211,-16.6288,-16.6349,-16.6402,-16.6444,-16.6479,-16.6509,-16.6522,-16.6542,-16.658,-16.6615,-16.6628,-16.6652,-16.6685,-16.6734,-16.6773,-16.6802,-16.6841,-16.6855,-16.6833,-16.6813,-16.6793,-16.6769,-16.6741,-16.6729,-16.6717,-16.6707,-16.6695,-16.6678,-16.6663,-16.6652,-16.6648,-16.6633,-16.6615,-16.6608,-16.6594,-16.6583,-16.6578,-16.6568,-16.656,-16.6552,-16.6545,-16.6534,-16.6524,-16.6516,-16.6495,-16.6468,-16.6451,-16.6433,-16.6419,-16.6401,-16.6384,-16.6369,-16.636,-16.6346,-16.6331,-16.6319,-16.631,-16.6304,-16.6293,-16.6285,-16.6271,-16.6267,-16.6261,-16.6256,-16.6248,-16.6248,-16.6242,-16.6242,-16.6233,-16.6237,-16.6247,-16.6276,-16.6293,-16.632,-16.6339,-16.6353,-16.6374,-16.6395,-16.6411,-16.6428,-16.6437,-16.6445,-16.6453,-16.6468,-16.6474,-16.6486,-16.6484,-16.6451,-16.6429,-16.6408,-16.6383,-16.6356,-16.6331,-16.6306,-16.6285,-16.6264,-16.6244,-16.6222,-16.6158,-16.6139,-16.613,-16.6109,-16.6093,-16.608,-16.6075,-16.607,-16.6068,-16.6064,-16.6062,-16.6058,-16.6058,-16.6055,-16.6054,-16.6053,-16.6051,-16.6048,-16.6039,-16.6039,-16.6042,-16.6044,-16.605,-16.605,-16.6053,-16.6053,-16.6056,-16.6055,-16.6054,-16.605,-16.6051,-16.6055,-16.6051,-16.6047,-16.6047,-16.6044,-16.6046,-16.6047,-16.6045,-16.6047,-16.6053,-16.6056,-16.6064,-16.6075,-16.6088,-16.6104,-16.6108,-16.6123,-16.6142,-16.6156,-16.6169,-16.6185,-16.6196,-16.6201,-16.6206,-16.6208,-16.6212,-16.6211,-16.6214,-16.6209,-16.6206,-16.6192,-16.6171,-16.6163,-16.615,-16.6134,-16.6118,-16.6102,-16.609,-16.6073,-16.6056,-16.6043,-16.6022,-16.6011,-16.5994,-16.5982,-16.5969,-16.5958,-16.5948,-16.594,-16.5926,-16.5914,-16.5905,-16.5895,-16.5883,-16.587,-16.5862,-16.5849,-16.5845,-16.5829,-16.5819,-16.5812,-16.5802,-16.5791,-16.5776,-16.5768,-16.575,-16.5737,-16.5727,-16.5701,-16.5694,-16.5682,-16.5667,-16.5654,-16.5647,-16.563,-16.562,-16.5608,-16.5591,-16.5568,-16.5562,-16.5546,-16.5533,-16.5516,-16.55,-16.5481,-16.5467,-16.5452,-16.5436,-16.5421,-16.5403,-16.5386,-16.5374,-16.5362,-16.5348,-16.5339,-16.5322,-16.5311,-16.5303,-16.5296,-16.5285,-16.5278,-16.5259,-16.5241,-16.5241,-16.5234,-16.5229,-16.5223,-16.5215,-16.5207,-16.5199,-16.5193,-16.5189,-16.5184,-16.5181,-16.518,-16.5178,-16.5177,-16.5176,-16.5176,-16.5178,-16.5177,-16.518,-16.5186,-16.519,-16.5198,-16.5206,-16.5212,-16.5217,-16.5225,-16.5234,-16.5246,-16.5257,-16.5267,-16.5273,-16.5273,-16.5271,-16.5269,-16.5268,-16.5279,-16.5285,-16.5286,-16.5287,-16.5294,-16.5299,-16.5303,-16.5308,-16.5306,-16.531,-16.5314,-16.5317,-16.5322,-16.533,-16.5335,-16.5346,-16.5349,-16.5348,-16.5358,-16.5361,-16.5361,-16.537,-16.5382,-16.5387,-16.5394,-16.5405,-16.5412,-16.5426,-16.5435,-16.5443,-16.5459,-16.5474,-16.5484,-16.5493,-16.5495,-16.5506,-16.5512,-16.5526,-16.5545,-16.5574,-16.559,-16.5611,-16.5635,-16.5653,-16.5681,-16.5694,-16.5711,-16.5735,-16.5753,-16.5771,-16.5794,-16.5814,-16.5841,-16.5867,-16.589,-16.5906,-16.592,-16.595,-16.5972,-16.5997,-16.6017,-16.6037,-16.6057,-16.6073,-16.6089,-16.6106,-16.6119,-16.614,-16.6155,-16.6175,-16.6202,-16.6217,-16.6234,-16.625,-16.6259,-16.6263,-16.628,-16.6297,-16.6308,-16.6327,-16.6336,-16.6352,-16.6377,-16.6389,-16.6404,-16.6424,-16.6439,-16.6457,-16.6466,-16.6492,-16.6516,-16.6541,-16.6568,-16.6596,-16.6625,-16.6659,-16.6684,-16.671,-16.6733,-16.6755,-16.6777,-16.6795,-16.6811,-16.6832,-16.685,-16.6869,-16.6888,-16.6906,-16.6928,-16.6948,-16.6963,-16.6981,-16.7002,-16.702,-16.7038,-16.7052,-16.7054,-16.7053,-16.7061,-16.7064,-16.7069,-16.7065,-16.7068,-16.7072,-16.7077,-16.7079,-16.7083,-16.7079,-16.7079,-16.7088,-16.7091,-16.7104,-16.7105,-16.7102,-16.71,-16.7098,-16.7103,-16.7102,-16.7097,-16.7092,-16.7089,-16.7089,-16.7079,-16.7078,-16.7079,-16.7074,-16.7075,-16.7071,-16.7073,-16.7078,-16.7085,-16.7088,-16.7096,-16.7109,-16.7124,-16.7138,-16.715,-16.7166,-16.7178,-16.7187,-16.7198,-16.7208,-16.7219,-16.7232,-16.7246,-16.7264,-16.7277,-16.729,-16.7306,-16.7318,-16.7333,-16.7344,-16.7362,-16.7372,-16.739,-16.7409,-16.7425,-16.7437,-16.7452,-16.7463,-16.7472,-16.7475,-16.7484,-16.7488,-16.7488,-16.7486,-16.7474,-16.7476,-16.746,-16.7457,-16.7441,-16.7433,-16.7431,-16.7423,-16.7422,-16.7401,-16.7387,-16.7368,-16.7348,-16.7322,-16.7312,-16.7303,-16.7281,-16.7254,-16.7239,-16.7232,-16.7228,-16.7222,-16.7223,-16.7235,-16.7241,-16.7248,-16.7256,-16.7265,-16.7278,-16.7291,-16.7304,-16.7321,-16.7329,-16.7341,-16.7356,-16.7375,-16.7392,-16.7402,-16.7413,-16.7428,-16.7437,-16.745,-16.7456,-16.7468,-16.7479,-16.7493,-16.7507,-16.7523,-16.7539,-16.7557,-16.7574,-16.7588,-16.7601,-16.7612,-16.7624,-16.7635,-16.7634,-16.763,-16.7617,-16.7597,-16.7574,-16.7565,-16.7559,-16.7551,-16.7535,-16.7533,-16.7527,-16.7526,-16.7533,-16.7537,-16.7548,-16.7559,-16.7573,-16.7585,-16.7593,-16.7615,-16.7625,-16.7633,-16.7641,-16.7655,-16.7677,-16.7693,-16.771,-16.7725,-16.7737,-16.7763,-16.7766,-16.778,-16.7798,-16.7818,-16.7835,-16.7849,-16.7998,-16.7996,-16.7996,-16.7946,-16.7938,-16.793,-16.793,-16.7913,-16.7913,-16.788,-16.7863,-16.7855,-16.7855,-16.7846,-16.7846,-16.7838,-16.7838,-16.7813,-16.778,-16.778,-16.778,-16.7771,-16.7772,-16.7763,-16.7763,-16.7755,-16.7755,-16.7746,-16.7746,-16.7738,-16.7738,-16.7746,-16.7746,-16.7738,-16.7738,-16.7721,-16.773,-16.7713,-16.7713,-16.7705,-16.7705,-16.7696,-16.7688,-16.7688,-16.768,-16.768,-16.7688,-16.7688,-16.768,-16.768,-16.7671,-16.7671,-16.768,-16.768,-16.7671,-16.768,-16.768,-16.7688,-16.768,-16.768,-16.7671,-16.7671,-16.7663,-16.7663,-16.7663,-16.7655,-16.7655,-16.7646,-16.7646,-16.7655,-16.7655,-16.7646,-16.7638,-16.763,-16.763,-16.7613,-16.7613,-16.76,-16.7584,-16.7563,-16.7563,-16.7571,-16.758,-16.7588,-16.7588,-16.758,-16.758,-16.7555,-16.7555,-16.7563,-16.7571,-16.7563,-16.7563,-16.7571,-16.7571,-16.7588,-16.7588,-16.7596,-16.763,-16.7629,-16.7638,-16.763,-16.763,-16.7621,-16.7613,-16.7605,-16.7605,-16.7596,-16.7596,-16.7605,-16.7604,-16.7613,-16.7613,-16.7588,-16.7588,-16.7571,-16.7563,-16.7546,-16.7538,-16.7509,-16.7484,-16.7451,-16.7425,-16.7418,-16.7384,-16.7317,-16.7284,-16.723,-16.7221,-16.7221,-16.723,-16.7238,-16.7238,-16.7238,-16.7255,-16.7226,-16.7192,-16.7176,-16.7155,-16.7155,-16.7146,-16.7125,-16.7109,-16.7105,-16.7092,-16.7059,-16.7025,-16.7,-16.6959,-16.6926,-16.69,-16.6884,-16.688,-16.6905,-16.6942,-16.7009,-16.7047,-16.7059,-16.7109,-16.7142,-16.7167,-16.7205,-16.7213,-16.7213,-16.7205,-16.7184,-16.7159,-16.7125,-16.7092,-16.7084,-16.7059,-16.7034,-16.6976,-16.6834,-16.6801,-16.6796,-16.6805,-16.6776,-16.6767,-16.6792,-16.6834,-16.6867,-16.69,-16.6926,-16.695,-16.7001,-16.7029,-16.7021,-16.7021,-16.7121,-16.7129,-16.7155,-16.7155,-16.7159,-16.7188,-16.7193,-16.7222,-16.7238,-16.7267,-16.7276,-16.7304,-16.7304,-16.7321,-16.743,-16.743,-16.7421,-16.7421,-16.7405,-16.738,-16.738,-16.7346,-16.7346,-16.7355,-16.7396,-16.7405,-16.7438,-16.7425,-16.7367,-16.7359,-16.7251,-16.7242,-16.7184,-16.7176,-16.7258,-16.7317,-16.7363,-16.7363,-16.7355,-16.7355,-16.7338,-16.7338,-16.733,-16.7338,-16.7338,-16.7355,-16.7355,-16.7325,-16.7292,-16.7284,-16.7251,-16.7242,-16.7221,-16.7221,-16.7234,-16.725,-16.728,-16.7271,-16.7271,-16.7271,-16.7288,-16.7288,-16.73,-16.7305,-16.7296,-16.7288,-16.7309,-16.7326,-16.7351,-16.7355,-16.7346,-16.7346,-16.7338,-16.7338,-16.7322,-16.7305,-16.7313,-16.7313,-16.7346,-16.7346,-16.7321,-16.7321,-16.7313,-16.7313,-16.728,-16.7246,-16.7213,-16.7167,-16.7142,-16.7117,-16.7101,-16.7063,-16.7013,-16.7005,-16.6963,-16.6934,-16.6917,-16.6892,-16.6884,-16.6851,-16.6834,-16.6813,-16.6805,-16.6763,-16.6767,-16.68,-16.6825,-16.6843,-16.6859,-16.6888,-16.6855,-16.6855,-16.6838,-16.6834,-16.6826,-16.6801,-16.6784,-16.6771,-16.6763,-16.6771,-16.6772,-16.6788,-16.6759,-16.6734,-16.6717,-16.6684,-16.6679,-16.6688,-16.6688,-16.668,-16.6675,-16.6659,-16.6621,-16.6621,-16.6567,-16.655,-16.6517,-16.6492,-16.6467,-16.6455,-16.6455,-16.648,-16.648,-16.6505,-16.6513,-16.6521,-16.6521,-16.6509,-16.6488,-16.6488,-16.648,-16.6471,-16.6442,-16.6434,-16.6421,-16.6418,-16.6405,-16.64,-16.6384,-16.6367,-16.6359,-16.6342,-16.6309,-16.6275,-16.6226,-16.6201,-16.6142,-16.6134,-16.6117,-16.6092,-16.6059,-16.6034,-16.6009,-16.5959,-16.5951,-16.5934,-16.5909,-16.5876,-16.5851,-16.5825,-16.5809,-16.5767,-16.5742,-16.5692,-16.568,-16.568,-16.5688,-16.5646,-16.5638,-16.5621,-16.5621,-16.5538,-16.553,-16.5509,-16.5488,-16.5488,-16.5472,-16.5471,-16.5413,-16.5421,-16.5421,-16.5446,-16.5522,-16.5521,-16.5492,-16.5475,-16.5467,-16.54,-16.5396,-16.5413,-16.5413,-16.5401,-16.5367,-16.5367,-16.5384,-16.5396,-16.5388,-16.5388,-16.5367,-16.5334,-16.5326,-16.53,-16.5292,-16.5267,-16.5259,-16.5217,-16.52,-16.5184,-16.5146,-16.513,-16.5168,-16.5201,-16.5259,-16.5267,-16.53,-16.5309,-16.5326,-16.5359,-16.5367,-16.5392,-16.5401,-16.5425,-16.5459,-16.5476,-16.5488,-16.5488,-16.5451,-16.5442,-16.5407,-16.5351,-16.5334,-16.5322,-16.5321,-16.5338,-16.5338,-16.5359,-16.5392,-16.5413,-16.5413,-16.5401,-16.5392,-16.5375,-16.5371,-16.5384,-16.54,-16.5409,-16.5433,-16.5467,-16.5476,-16.5534,-16.5542,-16.5567,-16.5575,-16.5609,-16.5617,-16.5651,-16.5659,-16.5676,-16.5692,-16.5725,-16.5746,-16.5746,-16.5717,-16.5705,-16.5767,-16.5817,-16.5851,-16.5867,-16.5938,-16.5942,-16.5984,-16.6013,-16.6013,-16.5996,-16.6009,-16.6017,-16.6051,-16.6067,-16.6113,-16.6155,-16.623,-16.6221,-16.6226,-16.6242,-16.6263,-16.6263,-16.6321,-16.6313,-16.6292,-16.6233,-16.62,-16.615,-16.6142,-16.6109,-16.61,-16.6059,-16.6051,-16.6043,-16.6084,-16.6093,-16.6109,-16.6142,-16.6168,-16.6175,-16.6259,-16.6267,-16.6301,-16.633,-16.6371,-16.6372,-16.633,-16.633,-16.6305,-16.6305,-16.6288,-16.6271,-16.6204,-16.6205,-16.6196,-16.6196,-16.6171,-16.6155,-16.6154,-16.6163,-16.6163,-16.6171,-16.6171,-16.6163,-16.6163,-16.6163,-16.6155,-16.6146,-16.6146,-16.6134,-16.61,-16.6076,-16.6063,-16.6063,-16.6068,-16.6092,-16.6109,-16.613,-16.613,-16.6101,-16.6092,-16.6076,-16.6042,-16.6009,-16.5976,-16.5875,-16.5867,-16.5851,-16.5742,-16.57,-16.5676,-16.5671,-16.5646,-16.5684,-16.5709,-16.5742,-16.5746,-16.5788,-16.5813,-16.5813,-16.5834,-16.5859,-16.5859,-16.5809,-16.5775,-16.5742,-16.5726,-16.5717,-16.5642,-16.5605,-16.5605,-16.5609,-16.5634,-16.5651,-16.5659,-16.5709,-16.5734,-16.575,-16.5771,-16.5763,-16.5763,-16.568,-16.568,-16.5638,-16.5592,-16.5591,-16.5584,-16.5576,-16.5567,-16.5558,-16.5527,-16.5511,-16.5492,-16.5477,-16.5464,-16.5457,-16.5453,-16.5462,-16.5463,-16.5462,-16.5455,-16.545,-16.5447,-16.5448,-16.5449,-16.5425,-16.5248,-16.5073,-16.4992,-16.4895,-16.4717,-16.4539,-16.4361,-16.4289,-16.4186,-16.4184,-16.4009,-16.3831,-16.365,-16.3475,-16.3298,-16.312,-16.2945,-16.2764,-16.2586,-16.2411,-16.2234,-16.2059,-16.1878,-16.17,-16.1523,-16.1348,-16.117,-16.1033,-16.0992,-16.0878,-16.0814,-16.0636,-16.0461,-16.0284,-16.0103,-16,-15.9992,-15.982,-15.9642,-15.9467,-15.9289,-15.9111,-15.8934,-15.8756,-15.8578,-15.8403,-15.8225,-15.8048,-15.787,-15.7692,-15.7514,-15.7339,-15.7159,-15.6981,-15.6806,-15.6628,-15.645,-15.6273,-15.6095,-15.6036,-15.5917,-15.5739,-15.5564,-15.5386,-15.5206,-15.5031,-15.4992,-15.4881,-15.4867,-15.4878,-15.4879,-15.4881,-15.4878,-15.4873,-15.4861,-15.4848,-15.4845,-15.4834,-15.4834,-15.4814,-15.4795,-15.4775,-15.4764,-15.4745,-15.4725,-15.4706,-15.4686,-15.4659,-15.4639,-15.4623,-15.4598,-15.4578,-15.4542,-15.4517,-15.4481,-15.4448,-15.4417,-15.4381,-15.4348,-15.4306,-15.4273,-15.4231,-15.4195,-15.4148,-15.4147,-15.4106,-15.4056,-15.4017,-15.397,-15.3928,-15.3886,-15.3853,-15.3811,-15.3761,-15.372,-15.3678,-15.3636,-15.3589,-15.3542,-15.3492,-15.3442,-15.3389,-15.3331,-15.3275,-15.3206,-15.3139,-15.3084,-15.3031,-15.2986,-15.2931,-15.2886,-15.285,-15.2814,-15.277,-15.2731,-15.2703,-15.2664,-15.2636,-15.2603,-15.2575,-15.2545,-15.2514,-15.2506,-15.2334,-15.2159,-15.1986,-15.1814,-15.1642,-15.1636,-15.1459,-15.1284,-15.1111,-15.0934,-15.0759,-15.0695,-15.0514,-15.0336,-15.0156,-14.9992,-14.9975,-14.9909,-14.9861,-14.9811,-14.9756,-14.97,-14.9642,-14.9586,-14.9523,-14.9453,-14.9389,-14.9323,-14.9259,-14.9195,-14.9136,-14.9081,-14.9023,-14.8973,-14.8923,-14.8873,-14.8823,-14.8778,-14.8734,-14.8684,-14.8631,-14.8595,-14.8556,-14.855,-14.8511,-14.847,-14.8431,-14.8407,-14.8234,-14.8213,-14.82,-14.8174,-14.8158,-14.8139,-14.8125,-14.8107,-14.8093,-14.8082,-14.8065,-14.8047,-14.8023,-14.8016,-14.7993,-14.7961,-14.7906,-14.7845,-14.7792,-14.7757,-14.7706,-14.7643,-14.7593,-14.754,-14.7497,-14.7441,-14.7362,-14.7303,-14.7239,-14.7201,-14.7186,-14.7175,-14.7161,-14.7145,-14.7135,-14.7121,-14.7105,-14.7095,-14.7092,-14.7069,-14.7042,-14.7016,-14.7009,-14.6999,-14.6973,-14.6943,-14.6919,-14.687,-14.682,-14.679,-14.675,-14.6709,-14.6673,-14.6631,-14.6587,-14.655,-14.6516,-14.6489,-14.6452,-14.6422,-14.6379,-14.6338,-14.6309,-14.6251,-14.6225,-14.6186,-14.6123,-14.609,-14.6053,-14.6043,-14.6033,-14.6013,-14.6005,-14.5992,-14.5985,-14.5972,-14.5959,-14.5945,-14.5934,-14.5922,-14.5907,-14.59,-14.5889,-14.5879,-14.5871,-14.5862,-14.5856,-14.5853,-14.5843,-14.5835,-14.5828,-14.5822,-14.5818,-14.5815,-14.5818,-14.5821,-14.5823,-14.5826,-14.5829,-14.5833,-14.5842,-14.5845,-14.5847,-14.5848,-14.5849,-14.5853,-14.5857,-14.5863,-14.5862,-14.5869,-14.5867,-14.5873,-14.5878,-14.5879,-14.5882,-14.5888,-14.5893,-14.5894,-14.5901,-14.5899,-14.5904,-14.5905,-14.5913,-14.5913,-14.5921,-14.5921,-14.5923,-14.5924,-14.5928,-14.5931,-14.5933,-14.5935,-14.593,-14.5926,-14.592,-14.592,-14.5914,-14.5912,-14.5907,-14.5907,-14.5899,-14.5893,-14.5894,-14.5889,-14.5884,-14.588,-14.5876,-14.587,-14.5867,-14.5864,-14.5863,-14.5859,-14.5852,-14.5848,-14.5845,-14.5846,-14.5847,-14.5845,-14.5846,-14.5848,-14.5854,-14.5855,-14.586,-14.5867,-14.5868,-14.5866,-14.5868,-14.5868,-14.5871,-14.5876,-14.5881,-14.5884,-14.5882,-14.5884,-14.5886,-14.5889,-14.5891,-14.5889,-14.5894,-14.5898,-14.5901,-14.5901,-14.5901,-14.5901,-14.59,-14.5903,-14.5908,-14.5909,-14.5915,-14.5912,-14.5918,-14.5919,-14.592,-14.5924,-14.5932,-14.5933,-14.5932,-14.5935,-14.5938,-14.5941,-14.5946,-14.5946,-14.595,-14.5953,-14.5954,-14.5954,-14.5961,-14.5961,-14.5975,-14.6003,-14.6029,-14.6055,-14.6085,-14.611,-14.6131,-14.6161,-14.6201,-14.6235,-14.6267,-14.6292,-14.6304,-14.6317,-14.6361,-14.6376,-14.6393,-14.6403,-14.6388,-14.6373,-14.6358,-14.6343,-14.6328,-14.6314,-14.6299,-14.6284,-14.6277,-14.6269,-14.6251,-14.6249,-14.6326,-14.6428,-14.6485,-14.6528,-14.6571,-14.6577,-14.6581,-14.6587,-14.659,-14.6628,-14.6634,-14.6642,-14.683,-14.6888,-14.6926,-14.7051,-14.7056,-14.7191,-14.7191,-14.7277,-14.7279,-14.7315,-14.7323,-14.7345,-14.7379,-14.7462,-14.7534,-14.7608,-14.7693,-14.7703,-14.7711,-14.7733,-14.7734,-14.7822,-14.7822,-14.7922,-14.796,-14.7981,-14.8002,-14.8027,-14.8028,-14.8036,-14.8111,-14.8203,-14.836,-14.8496,-14.8686,-14.8891,-14.9095,-14.9362,-14.9629,-14.9896,-14.9989,-15.0162,-15.0229,-15.0515,-15.0802,-15.0855,-15.0887,-15.1088,-15.1375,-15.1661,-15.1947,-15.2234,-15.252,-15.2807,-15.3093,-15.3366,-15.3639,-15.3912,-15.4185,-15.4458,-15.4454,-15.4446,-15.4435,-15.443,-15.4423,-15.4419,-15.4416,-15.4408,-15.4405,-15.4402,-15.4395,-15.4392,-15.4378,-15.4377,-15.4371,-15.437,-15.436,-15.4361,-15.4353,-15.4358,-15.4357,-15.4357,-15.4371,-15.4397,-15.4425,-15.4449,-15.4493,-15.4532,-15.4569,-15.4594,-15.4617,-15.4646,-15.4672,-15.4686,-15.471,-15.4731,-15.4759,-15.4786,-15.4817,-15.4844,-15.4864,-15.4884,-15.4904,-15.4923,-15.4953,-15.4973,-15.499,-15.501,-15.503,-15.5046,-15.5054,-15.5078,-15.5103,-15.513,-15.5148,-15.5175,-15.5204,-15.5238,-15.5276,-15.5313,-15.5348,-15.5364,-15.5386,-15.5408,-15.5432,-15.5452,-15.5496,-15.5523,-15.5562,-15.5597,-15.563,-15.5665,-15.5698,-15.5736,-15.5757,-15.5773,-15.5785,-15.5815,-15.5846,-15.5868,-15.5909,-15.5937,-15.5962,-15.5991,-15.6009,-15.603,-15.6064,-15.6078,-15.6099,-15.6122,-15.6148,-15.6185,-15.6218,-15.6251,-15.6272,-15.6307,-15.6322,-15.6351,-15.6379,-15.6408,-15.6434,-15.6457,-15.6484,-15.6509,-15.6544,-15.6565,-15.6587,-15.6608,-15.6624,-15.6649,-15.6666,-15.6689,-15.6711,-15.6726,-15.673,-15.6735,-15.6739,-15.6744,-15.6749,-15.6752,-15.6755,-15.676,-15.6773,-15.6841,-15.6903,-15.692,-15.6876,-15.6751,-15.6621,-15.6525,-15.6406,-15.6354,-15.6356,-15.639,-15.6435,-15.6498,-15.6554,-15.6628,-15.6664,-15.6722,-15.6791,-15.687,-15.6946,-15.7003,-15.7066,-15.7117,-15.7158,-15.7223,-15.7274,-15.7322,-15.7369,-15.7413,-15.7457,-15.7508,-15.7552,-15.7584,-15.7603,-15.7617,-15.7633,-15.7659,-15.7678,-15.7699,-15.7715,-15.772,-15.7727,-15.7729,-15.7737,-15.775,-15.7765,-15.7779,-15.7798,-15.7815,-15.7829,-15.7851,-15.7868,-15.7892,-15.793,-15.7962,-15.7989,-15.8002,-15.8018,-15.8023,-15.8028],"lat":[14.9801,14.9979,14.9993,15.0009,15.0033,15.0071,15.0123,15.0183,15.0234,15.0287,15.0327,15.0389,15.0441,15.0485,15.0511,15.0512,15.0494,15.0455,15.0417,15.0378,15.0334,15.025,15.02,15.0165,15.0109,15.0064,15.0027,14.9964,14.9934,14.9908,14.9887,14.9867,14.9849,14.9835,14.9826,14.9816,14.9802,14.9792,14.9788,14.9786,14.9788,14.9785,14.9776,14.9771,14.9768,14.9766,14.9761,14.9756,14.9749,14.9746,14.9741,14.9736,14.9731,14.9726,14.9722,14.9722,14.9722,14.9721,14.9726,14.9723,14.9721,14.9716,14.9708,14.9701,14.9686,14.967,14.9653,14.9634,14.9618,14.9599,14.9583,14.9564,14.9557,14.9553,14.9553,14.9554,14.9557,14.9555,14.9562,14.9563,14.9559,14.9563,14.9567,14.9568,14.9577,14.9586,14.9589,14.959,14.959,14.9587,14.959,14.9595,14.9593,14.9594,14.9597,14.9599,14.9598,14.9604,14.9605,14.9606,14.961,14.961,14.9611,14.9607,14.9608,14.9602,14.9597,14.9588,14.9585,14.9581,14.958,14.9575,14.9573,14.9566,14.9565,14.9565,14.9564,14.9564,14.957,14.9574,14.9579,14.9585,14.9597,14.9603,14.9609,14.9608,14.9612,14.9617,14.9621,14.9626,14.9631,14.9631,14.9628,14.9626,14.9635,14.9636,14.9632,14.9633,14.963,14.9629,14.9628,14.9625,14.962,14.961,14.9603,14.9592,14.9589,14.9579,14.9574,14.9562,14.9558,14.9544,14.9534,14.9535,14.9527,14.952,14.9519,14.9525,14.9528,14.953,14.9546,14.956,14.9575,14.9598,14.962,14.9633,14.965,14.9668,14.9677,14.9694,14.971,14.9722,14.9732,14.9751,14.9761,14.9771,14.9783,14.9799,14.9806,14.982,14.9836,14.9844,14.985,14.986,14.987,14.9879,14.9888,14.9898,14.9903,14.9902,14.9913,14.9921,14.9925,14.9934,14.9938,14.9947,14.9951,14.9954,14.9957,14.9963,14.9969,14.9975,14.9982,15.0009,15.001,15.0018,15.0027,15.0039,15.0048,15.0057,15.0067,15.0086,15.0094,15.0107,15.0119,15.0127,15.0137,15.0144,15.0154,15.0165,15.0171,15.018,15.0189,15.0192,15.0199,15.0207,15.021,15.021,15.0212,15.0219,15.0221,15.0218,15.0222,15.0226,15.0228,15.023,15.0234,15.0236,15.0239,15.0242,15.0244,15.0246,15.0246,15.0245,15.0244,15.0243,15.0243,15.0243,15.0238,15.0235,15.0222,15.022,15.0213,15.0194,15.0176,15.0162,15.0153,15.0148,15.0143,15.0136,15.0125,15.0109,15.0084,15.0061,15.0036,15.0026,15.0008,14.9993,14.9984,14.9982,14.9978,14.9976,14.997,14.9967,14.9955,14.9945,14.9918,14.9897,14.9883,14.9862,14.9842,14.9824,14.9811,14.9796,14.9777,14.9757,14.9737,14.9724,14.9708,14.9704,14.9688,14.9673,14.9657,14.9642,14.9626,14.9618,14.96,14.9586,14.9574,14.9566,14.9553,14.9532,14.9506,14.9497,14.948,14.9451,14.9402,14.9362,14.9335,14.9292,14.9273,14.9253,14.9235,14.9219,14.9199,14.918,14.9166,14.9155,14.9143,14.9128,14.9108,14.9093,14.9082,14.9068,14.9053,14.905,14.9022,14.8988,14.8963,14.895,14.8938,14.8922,14.8905,14.8887,14.8874,14.8865,14.8858,14.8844,14.8822,14.8805,14.8796,14.8784,14.8772,14.8759,14.8742,14.8729,14.8718,14.8701,14.8688,14.8676,14.8664,14.8651,14.8632,14.8612,14.8591,14.8553,14.8529,14.8512,14.8491,14.8479,14.8468,14.8449,14.8411,14.8398,14.838,14.8363,14.8344,14.8329,14.8319,14.83,14.8281,14.8264,14.8247,14.823,14.822,14.8206,14.8187,14.8173,14.8148,14.812,14.8073,14.8057,14.8044,14.8027,14.8005,14.7991,14.7978,14.7961,14.7946,14.7936,14.7925,14.7793,14.7777,14.7768,14.7749,14.773,14.7699,14.7682,14.7661,14.764,14.7616,14.758,14.7558,14.7543,14.7519,14.749,14.747,14.7452,14.7429,14.7408,14.7395,14.7363,14.7339,14.7312,14.7282,14.7258,14.7228,14.7207,14.7193,14.7175,14.7149,14.7125,14.7091,14.7061,14.704,14.7019,14.6997,14.6964,14.6945,14.6921,14.6898,14.6865,14.6849,14.6831,14.6802,14.6777,14.6744,14.6737,14.6717,14.67,14.6682,14.6669,14.6646,14.6619,14.659,14.6564,14.6536,14.6513,14.6487,14.646,14.6441,14.6422,14.6382,14.6371,14.6363,14.6352,14.6341,14.6332,14.6322,14.6314,14.6302,14.6281,14.6273,14.6258,14.6252,14.6231,14.6217,14.6201,14.6195,14.6182,14.6176,14.6168,14.6161,14.615,14.6142,14.613,14.6122,14.6117,14.6099,14.6093,14.6089,14.6079,14.6072,14.6072,14.6071,14.6067,14.606,14.6057,14.6053,14.605,14.6044,14.6044,14.604,14.6038,14.6036,14.6038,14.604,14.6042,14.6042,14.6038,14.6042,14.6043,14.6046,14.6045,14.6046,14.6046,14.6047,14.6051,14.6054,14.6053,14.6056,14.6059,14.6059,14.6057,14.6057,14.6058,14.6059,14.606,14.6061,14.6061,14.606,14.606,14.606,14.6057,14.606,14.6052,14.6043,14.6037,14.6024,14.6006,14.5987,14.5964,14.5948,14.5934,14.5917,14.5904,14.5892,14.588,14.5864,14.5845,14.5831,14.5811,14.5795,14.5784,14.5766,14.5753,14.5741,14.5728,14.5715,14.5704,14.5691,14.5677,14.5655,14.5638,14.5624,14.5615,14.5599,14.5585,14.5575,14.5543,14.5521,14.5503,14.5484,14.5466,14.5455,14.5432,14.541,14.5386,14.5373,14.535,14.5328,14.5316,14.5303,14.529,14.5276,14.5252,14.5233,14.5224,14.5203,14.5197,14.5184,14.517,14.5151,14.5132,14.5111,14.5096,14.508,14.5061,14.5044,14.503,14.5015,14.5001,14.4992,14.4963,14.4955,14.494,14.4921,14.4909,14.4891,14.4871,14.4855,14.4834,14.4814,14.4797,14.4776,14.4765,14.475,14.4735,14.4712,14.4698,14.4677,14.4663,14.4646,14.4624,14.4611,14.4601,14.4586,14.4561,14.4554,14.4545,14.4536,14.4529,14.4525,14.4523,14.4515,14.4506,14.4501,14.449,14.448,14.4472,14.4451,14.4437,14.4421,14.4406,14.4391,14.4383,14.4366,14.435,14.4335,14.4314,14.43,14.4287,14.428,14.4275,14.4271,14.4267,14.4266,14.4263,14.4262,14.4266,14.4259,14.4258,14.4261,14.426,14.4257,14.4257,14.4259,14.4256,14.4259,14.4251,14.4239,14.423,14.4222,14.4199,14.418,14.416,14.4142,14.4121,14.4098,14.4081,14.4067,14.4046,14.4028,14.4004,14.398,14.395,14.3918,14.3891,14.3859,14.3835,14.3805,14.3778,14.3752,14.373,14.3707,14.368,14.3647,14.363,14.3608,14.3584,14.3551,14.3523,14.3499,14.3473,14.3449,14.3419,14.339,14.3363,14.3332,14.3306,14.3287,14.3262,14.3237,14.321,14.3186,14.3166,14.314,14.3117,14.3097,14.3069,14.3052,14.3024,14.3004,14.2984,14.2963,14.294,14.2918,14.2894,14.2876,14.2861,14.284,14.2823,14.2807,14.2791,14.2772,14.2752,14.2731,14.2713,14.2699,14.2684,14.2669,14.2658,14.2638,14.2623,14.2599,14.2574,14.2552,14.2539,14.2524,14.2504,14.2478,14.2454,14.2432,14.2407,14.2376,14.2348,14.2317,14.2294,14.2273,14.2252,14.2231,14.2222,14.2209,14.2201,14.2189,14.2173,14.2156,14.2144,14.2126,14.2112,14.2111,14.2113,14.2104,14.2091,14.208,14.2064,14.2048,14.2035,14.2011,14.1988,14.1971,14.1952,14.1932,14.1915,14.1895,14.1889,14.1879,14.1864,14.1848,14.1829,14.1823,14.1827,14.1827,14.1825,14.182,14.1806,14.1789,14.1766,14.1754,14.1739,14.1729,14.1715,14.171,14.1708,14.1705,14.1698,14.1691,14.1674,14.1655,14.1635,14.1625,14.1609,14.1594,14.1576,14.1562,14.1546,14.1532,14.1518,14.1502,14.1483,14.1468,14.1458,14.1448,14.1432,14.1415,14.141,14.1403,14.1401,14.1407,14.1412,14.142,14.1424,14.142,14.1412,14.1403,14.1395,14.1385,14.138,14.1379,14.138,14.1387,14.1394,14.1403,14.1414,14.1411,14.1403,14.14,14.1397,14.124,14.1238,14.1221,14.1137,14.1104,14.1096,14.1079,14.1055,14.1037,14.0971,14.0913,14.0904,14.0871,14.0863,14.0746,14.0737,14.0721,14.0696,14.062,14.0587,14.0579,14.0571,14.0546,14.0537,14.0504,14.0496,14.0454,14.0446,14.0379,14.0371,14.0321,14.0312,14.0229,14.0221,14.0204,14.0188,14.0096,14.0071,14.0054,14.0046,13.9996,13.9913,13.9904,13.9863,13.9854,13.9779,13.9771,13.9738,13.9729,13.9671,13.9663,13.9579,13.9571,13.9513,13.9504,13.9496,13.9446,13.9437,13.9421,13.9346,13.9338,13.9213,13.9204,13.918,13.8971,13.8963,13.8896,13.8887,13.8771,13.8763,13.868,13.8671,13.8554,13.8546,13.8496,13.8471,13.8346,13.8333,13.8334,13.8354,13.8471,13.8479,13.8562,13.857,13.8588,13.8596,13.8621,13.8662,13.8688,13.8696,13.8729,13.8738,13.8887,13.8896,13.8921,13.8954,13.8988,13.9004,13.9038,13.9063,13.9071,13.9079,13.9129,13.9138,13.9221,13.9229,13.9354,13.9362,13.9446,13.9454,13.9487,13.9496,13.9579,13.9629,13.9646,13.9662,13.9696,13.9721,13.9762,13.9792,13.98,13.9825,13.9825,13.9817,13.9817,13.9891,13.99,13.9946,13.9963,13.9988,13.9996,14.0004,14.0174,14.0188,14.0204,14.0233,14.0233,14.0242,14.0263,14.0287,14.0304,14.0325,14.0325,14.0313,14.0309,14.0342,14.0359,14.0383,14.0384,14.04,14.0442,14.0442,14.0429,14.0379,14.0342,14.0317,14.0288,14.0267,14.0241,14.0209,14.02,14.0171,14.0154,14.0129,14.0113,14.0092,14.0092,14.0108,14.01,14.0109,14.0117,14.0142,14.0167,14.0167,14.0183,14.0171,14.0154,14.0125,14.0125,14.0062,14.0025,14.0025,14.0041,14.0067,14.0067,14.0042,14.0012,13.9996,13.9987,13.9887,13.9854,13.9829,13.9771,13.9767,13.9796,13.9809,13.978,13.9745,13.9717,13.9717,13.9688,13.9679,13.9654,13.9546,13.9521,13.9513,13.9471,13.9437,13.9404,13.9304,13.9279,13.9271,13.9254,13.9212,13.9179,13.9121,13.9108,13.9109,13.91,13.91,13.9092,13.9092,13.9083,13.9069,13.9059,13.9037,13.9012,13.9004,13.8979,13.8954,13.8871,13.8854,13.8846,13.8829,13.8813,13.8763,13.8734,13.8717,13.8733,13.8692,13.8692,13.8671,13.8663,13.8659,13.8675,13.8671,13.8662,13.8621,13.8571,13.8554,13.8538,13.8525,13.8546,13.8554,13.8604,13.8617,13.8592,13.8584,13.8554,13.8546,13.8504,13.8496,13.8471,13.8454,13.8412,13.8404,13.8379,13.8379,13.8363,13.8312,13.8262,13.8254,13.8238,13.8188,13.8155,13.8088,13.8041,13.8042,13.8025,13.8025,13.7988,13.7904,13.7854,13.7771,13.7742,13.7733,13.7733,13.7742,13.7742,13.775,13.7779,13.7804,13.7846,13.7858,13.7833,13.7825,13.7803,13.7784,13.7787,13.7812,13.7821,13.7838,13.785,13.785,13.7875,13.7875,13.7887,13.7971,13.7979,13.8029,13.8071,13.8075,13.805,13.805,13.8084,13.8079,13.8062,13.8029,13.8021,13.7991,13.7992,13.8029,13.8038,13.81,13.8108,13.8109,13.8125,13.8134,13.8146,13.8196,13.8229,13.8254,13.8288,13.8354,13.8363,13.8379,13.8384,13.8354,13.8321,13.8312,13.8271,13.8242,13.8242,13.8254,13.8284,13.8271,13.8242,13.8233,13.8234,13.8242,13.8242,13.8267,13.8275,13.83,13.8333,13.8325,13.8334,13.8333,13.835,13.8358,13.8375,13.8359,13.8358,13.8367,13.8367,13.8383,13.8392,13.8434,13.8442,13.8433,13.8475,13.8483,13.8534,13.8537,13.8571,13.8579,13.8612,13.8662,13.8679,13.8688,13.8763,13.8788,13.88,13.8829,13.8846,13.8871,13.892,13.8854,13.8846,13.8821,13.8771,13.8696,13.8662,13.8625,13.8625,13.8617,13.8625,13.8654,13.8671,13.8688,13.87,13.87,13.8692,13.8692,13.8679,13.8671,13.8663,13.8642,13.8642,13.865,13.865,13.8658,13.8658,13.8667,13.8683,13.8708,13.8708,13.8746,13.8788,13.845,13.8467,13.8467,13.8459,13.8459,13.845,13.845,13.8433,13.8425,13.8425,13.8417,13.8416,13.84,13.84,13.8388,13.8354,13.8334,13.8325,13.8325,13.8325,13.83,13.8296,13.8279,13.8262,13.8204,13.82,13.8234,13.8237,13.8271,13.8275,13.825,13.825,13.8287,13.83,13.83,13.8309,13.8309,13.8309,13.8317,13.8317,13.8308,13.8309,13.83,13.83,13.8292,13.8292,13.8284,13.8283,13.8267,13.825,13.8221,13.8171,13.8167,13.8146,13.8142,13.8083,13.8067,13.8067,13.7988,13.7975,13.795,13.7929,13.7913,13.7879,13.7875,13.7883,13.7883,13.7875,13.7821,13.7754,13.7687,13.7604,13.76,13.7617,13.7604,13.7555,13.7438,13.732,13.7308,13.7309,13.7325,13.7325,13.7333,13.7334,13.7342,13.7342,13.735,13.7312,13.7308,13.7292,13.73,13.7283,13.7283,13.7275,13.7267,13.7258,13.7258,13.7237,13.7146,13.6996,13.6912,13.6896,13.6846,13.6829,13.6795,13.6779,13.6646,13.6629,13.662,13.6579,13.6521,13.6521,13.6546,13.6554,13.6571,13.6579,13.6613,13.6621,13.664,13.6679,13.6696,13.6696,13.6637,13.6608,13.6609,13.6642,13.6637,13.6562,13.6558,13.6592,13.6592,13.6571,13.6546,13.65,13.65,13.6484,13.6475,13.645,13.6434,13.6433,13.6442,13.6442,13.6558,13.6575,13.6575,13.6596,13.6621,13.6534,13.6508,13.6492,13.6479,13.6437,13.6396,13.6379,13.6358,13.6359,13.6342,13.6342,13.6308,13.6292,13.6292,13.63,13.63,13.6263,13.6229,13.6225,13.625,13.625,13.6242,13.6241,13.6267,13.6267,13.6254,13.6238,13.6204,13.6121,13.6112,13.6071,13.6,13.6,13.5986,13.5974,13.5969,13.597,13.5996,13.6013,13.6031,13.6041,13.6045,13.6016,13.6003,13.6,13.5996,13.5982,13.5954,13.5935,13.5913,13.5902,13.5896,13.5896,13.5899,13.5899,13.59,13.5902,13.5902,13.5902,13.5904,13.5904,13.5904,13.5904,13.5907,13.5907,13.591,13.591,13.591,13.5913,13.5913,13.5915,13.5915,13.5915,13.5918,13.5918,13.5918,13.5921,13.5915,13.5915,13.5915,13.5918,13.5918,13.5918,13.5918,13.5918,13.5921,13.5921,13.5921,13.5918,13.5918,13.5921,13.5921,13.5921,13.5924,13.5915,13.5918,13.5918,13.5918,13.5913,13.5913,13.5913,13.5915,13.5907,13.591,13.591,13.591,13.5904,13.5904,13.5904,13.5907,13.5899,13.5902,13.5902,13.5902,13.5902,13.5896,13.5896,13.5896,13.5896,13.5895,13.589,13.5918,13.5979,13.6,13.6049,13.611,13.6171,13.6229,13.6285,13.6346,13.6402,13.6404,13.6457,13.6507,13.6557,13.6613,13.666,13.671,13.676,13.6807,13.6852,13.6899,13.6949,13.699,13.704,13.7077,13.7118,13.7154,13.719,13.7227,13.7263,13.7299,13.7327,13.7363,13.739,13.7429,13.7457,13.7458,13.7488,13.7515,13.7546,13.7577,13.7604,13.7635,13.7671,13.7699,13.7721,13.7752,13.7779,13.781,13.7832,13.7854,13.7877,13.7902,13.7915,13.7932,13.7949,13.7952,13.7949,13.7935,13.7918,13.789,13.7879,13.7854,13.7821,13.7788,13.7763,13.7729,13.7688,13.7657,13.7615,13.7577,13.7535,13.7496,13.7454,13.7443,13.7538,13.7635,13.7738,13.7835,13.794,13.7946,13.8002,13.8057,13.8118,13.8174,13.8238,13.826,13.819,13.8121,13.8046,13.7986,13.7979,13.7946,13.7954,13.7979,13.7993,13.801,13.8027,13.8043,13.8052,13.8054,13.8049,13.8046,13.804,13.8038,13.8027,13.8013,13.8002,13.7985,13.7965,13.7946,13.7927,13.7902,13.7877,13.7857,13.7838,13.7804,13.7782,13.7779,13.7754,13.7721,13.7693,13.7673,13.7829,13.7849,13.7881,13.7913,13.7935,13.7962,13.7979,13.8004,13.8029,13.8046,13.8064,13.8078,13.8094,13.8095,13.8098,13.8102,13.8107,13.8109,13.8109,13.8112,13.8115,13.8117,13.812,13.8126,13.813,13.814,13.8145,13.8152,13.8162,13.8174,13.8208,13.8232,13.827,13.8304,13.8331,13.8366,13.8398,13.8428,13.8436,13.8488,13.8549,13.8602,13.862,13.8651,13.8714,13.8791,13.8836,13.8858,13.8887,13.8905,13.8913,13.8939,13.8974,13.9,13.9024,13.9051,13.9077,13.9093,13.9111,13.9137,13.9166,13.9188,13.9211,13.9242,13.9263,13.9286,13.9312,13.9335,13.9366,13.9374,13.9382,13.9412,13.9436,13.9473,13.95,13.9522,13.9564,13.9599,13.963,13.9663,13.9698,13.9724,13.9761,13.9787,13.9811,13.9834,13.9855,13.9863,13.9887,13.991,13.994,13.9957,13.9984,14.0001,14.0007,14.0014,14.0022,14.0041,14.0066,14.0079,14.0097,14.0111,14.0127,14.0162,14.0192,14.0221,14.0249,14.0269,14.0307,14.0344,14.039,14.0408,14.0436,14.0462,14.0494,14.0517,14.0537,14.0572,14.0591,14.0617,14.0629,14.0676,14.0695,14.0723,14.0743,14.0773,14.0794,14.0824,14.0837,14.0854,14.0873,14.0892,14.0909,14.0948,14.0959,14.0976,14.0995,14.1021,14.1036,14.1058,14.1076,14.1102,14.1118,14.1143,14.1172,14.1181,14.1214,14.1234,14.1251,14.1278,14.13,14.1325,14.1351,14.1367,14.1378,14.1385,14.1396,14.1416,14.144,14.1465,14.1493,14.1527,14.1551,14.1581,14.161,14.162,14.1638,14.1663,14.1692,14.1721,14.1748,14.1768,14.1784,14.1812,14.1834,14.1859,14.1884,14.1909,14.1928,14.1944,14.1958,14.1971,14.1989,14.2012,14.2026,14.2057,14.2074,14.2091,14.213,14.2154,14.2173,14.2192,14.2213,14.2243,14.2253,14.2282,14.231,14.2348,14.2381,14.2402,14.2431,14.2457,14.2484,14.2514,14.2537,14.2552,14.2581,14.2592,14.2601,14.2623,14.2649,14.2676,14.2699,14.2721,14.274,14.2769,14.2812,14.2832,14.2868,14.2895,14.2906,14.2919,14.2958,14.2974,14.2986,14.301,14.3278,14.3546,14.3813,14.4081,14.4349,14.4617,14.4885,14.5153,14.5282,14.542,14.5623,14.5649,14.5672,14.5703,14.572,14.5733,14.5751,14.5754,14.5755,14.5758,14.5758,14.5761,14.5761,14.5763,14.5793,14.5804,14.5807,14.5817,14.5818,14.5827,14.5827,14.583,14.5831,14.5827,14.5831,14.583,14.5829,14.5829,14.5877,14.5926,14.598,14.5986,14.5992,14.6005,14.6005,14.6057,14.6057,14.612,14.614,14.6151,14.6162,14.6179,14.6181,14.6183,14.6206,14.6233,14.6413,14.6551,14.6551,14.6785,14.7019,14.6851,14.6684,14.6516,14.6457,14.6348,14.6097,14.615,14.6202,14.6212,14.6218,14.6255,14.6307,14.6359,14.6412,14.6464,14.6516,14.6569,14.6621,14.6786,14.695,14.7115,14.728,14.7445,14.7476,14.7532,14.7582,14.763,14.7667,14.7726,14.7766,14.781,14.7842,14.7886,14.7929,14.7966,14.8003,14.8051,14.8078,14.8101,14.8127,14.8155,14.8177,14.8208,14.8222,14.823,14.8241,14.8249,14.8254,14.8259,14.8281,14.8293,14.8308,14.8322,14.8334,14.8346,14.8355,14.8358,14.8363,14.8371,14.8377,14.8384,14.839,14.8395,14.8401,14.8401,14.8398,14.8395,14.8395,14.8393,14.8389,14.8386,14.8377,14.837,14.8366,14.8358,14.8349,14.8339,14.8333,14.8322,14.8314,14.8304,14.8287,14.8272,14.8259,14.8253,14.8245,14.8237,14.8237,14.825,14.8273,14.8292,14.8316,14.8335,14.8352,14.8378,14.8395,14.8417,14.843,14.8441,14.8448,14.8469,14.8482,14.8505,14.8531,14.8551,14.856,14.8577,14.8591,14.8597,14.8611,14.8621,14.8632,14.8634,14.8641,14.8645,14.865,14.8656,14.8659,14.8669,14.8669,14.8671,14.8673,14.8674,14.8674,14.8676,14.8686,14.8688,14.8684,14.8685,14.8687,14.8688,14.869,14.869,14.8692,14.8691,14.8697,14.8695,14.8673,14.866,14.8649,14.8635,14.8626,14.8615,14.8597,14.8572,14.8546,14.8578,14.8622,14.8721,14.8815,14.8859,14.8832,14.8794,14.8794,14.8884,14.8965,14.903,14.9102,14.9151,14.9206,14.9305,14.9421,14.9455,14.9459,14.9416,14.9314,14.9208,14.9132,14.9086,14.9062,14.9023,14.9,14.8987,14.8982,14.8977,14.8984,14.9013,14.9042,14.9085,14.9119,14.915,14.9179,14.9222,14.9249,14.9282,14.9306,14.933,14.9358,14.9391,14.9412,14.9438,14.9459,14.9479,14.9497,14.9521,14.9547,14.9568,14.9588,14.9601,14.9611,14.9619,14.9635,14.9661,14.968,14.9688,14.9801]}]],[[{"lng":[-16.7967,-16.7975,-16.798,-16.798,-16.7959,-16.7951,-16.7946,-16.7946,-16.7967],"lat":[12.4834,12.4834,12.4837,12.4846,12.4867,12.4867,12.4862,12.4854,12.4834]},{"lng":[-16.7959,-16.7947,-16.7946,-16.795,-16.7967,-16.798,-16.798,-16.7976,-16.7959],"lat":[12.4884,12.4896,12.4904,12.4908,12.4908,12.4896,12.4887,12.4884,12.4884]},{"lng":[-16.7922,-16.7921,-16.7925,-16.7942,-16.7955,-16.7942,-16.7926,-16.7922],"lat":[12.492,12.4938,12.4942,12.4942,12.4929,12.4917,12.4917,12.492]},{"lng":[-16.6746,-16.6759,-16.6825,-16.6838,-16.6825,-16.6759,-16.6746],"lat":[12.5038,12.505,12.505,12.5038,12.5025,12.5025,12.5038]},{"lng":[-16.515,-16.5138,-16.5151,-16.5192,-16.5196,-16.5176,-16.515],"lat":[12.5442,12.5471,12.5484,12.5483,12.5462,12.5441,12.5442]},{"lng":[-16.6988,-16.7001,-16.7017,-16.7059,-16.7075,-16.7092,-16.7109,-16.7122,-16.7121,-16.7193,-16.7217,-16.7238,-16.7238,-16.7209,-16.7209,-16.7284,-16.7317,-16.7358,-16.7371,-16.7372,-16.738,-16.738,-16.7388,-16.7388,-16.7396,-16.7396,-16.7371,-16.7355,-16.7346,-16.7267,-16.7242,-16.7208,-16.7168,-16.7151,-16.7142,-16.6951,-16.6942,-16.6909,-16.6901,-16.6817,-16.6788,-16.6788,-16.6796,-16.6796,-16.6805,-16.6805,-16.6796,-16.6796,-16.6788,-16.6788,-16.6817,-16.6859,-16.69,-16.693,-16.693,-16.6905,-16.6905,-16.6934,-16.6938,-16.693,-16.693,-16.6938,-16.6938,-16.6938,-16.6975,-16.6984,-16.7,-16.7034,-16.7059,-16.7071,-16.708,-16.7084,-16.7097,-16.7096,-16.7088,-16.7075,-16.7017,-16.7009,-16.6996,-16.6988,-16.6988],"lat":[12.5579,12.5592,12.5592,12.5567,12.5567,12.5592,12.5592,12.558,12.5571,12.5487,12.5458,12.5462,12.5479,12.5492,12.55,12.55,12.5525,12.5525,12.5512,12.5487,12.5479,12.5429,12.5421,12.5346,12.5337,12.5321,12.5279,12.5263,12.5229,12.515,12.5142,12.5108,12.5092,12.5092,12.5084,12.5084,12.5092,12.5092,12.51,12.51,12.5137,12.5163,12.517,12.5196,12.5204,12.5279,12.5288,12.5312,12.5321,12.5388,12.5408,12.5409,12.5434,12.5429,12.5404,12.5379,12.5371,12.535,12.5363,12.5371,12.5388,12.5396,12.5434,12.5454,12.5458,12.545,12.545,12.5475,12.5475,12.5463,12.5446,12.5408,12.5421,12.5463,12.5479,12.5492,12.5517,12.5534,12.5537,12.5555,12.5579]},{"lng":[-16.3234,-16.3226,-16.3221,-16.3221,-16.3226,-16.3242,-16.3255,-16.3242,-16.3234],"lat":[12.5892,12.5892,12.5896,12.5904,12.5908,12.5909,12.5896,12.5883,12.5892]},{"lng":[-16.3242,-16.3238,-16.3238,-16.323,-16.323,-16.3234,-16.3251,-16.3263,-16.3263,-16.3259,-16.3242],"lat":[12.5917,12.5921,12.5929,12.5937,12.5946,12.595,12.595,12.5938,12.592,12.5917,12.5917]},{"lng":[-16.3276,-16.3263,-16.3263,-16.3276,-16.33,-16.3313,-16.3313,-16.3292,-16.3276],"lat":[12.595,12.5962,12.5988,12.6,12.6,12.5987,12.5971,12.595,12.595]},{"lng":[-16.6967,-16.7009,-16.7021,-16.7021,-16.7021,-16.7009,-16.6993,-16.6984,-16.6925,-16.6896,-16.6896,-16.6918,-16.6934,-16.6967],"lat":[12.6042,12.6058,12.6046,12.5982,12.5954,12.5942,12.5942,12.595,12.595,12.5971,12.5979,12.6,12.6,12.6042]},{"lng":[-16.4663,-16.4663,-16.4684,-16.4705,-16.4705,-16.4684,-16.4663],"lat":[12.6538,12.6562,12.6592,12.6588,12.6554,12.6525,12.6538]},{"lng":[-16.7894,-16.7871,-16.7855,-16.7846,-16.7838,-16.7842,-16.7858,-16.788,-16.7888,-16.7892,-16.7896,-16.7896,-16.7894],"lat":[12.7504,12.7504,12.7538,12.7588,12.7596,12.7609,12.7608,12.7587,12.7554,12.7551,12.7546,12.7504,12.7504]},{"lng":[-12.3659,-12.3711,-12.3767,-12.3811,-12.3861,-12.3906,-12.395,-12.3986,-12.4023,-12.4059,-12.4089,-12.4125,-12.4148,-12.4184,-12.422,-12.4256,-12.43,-12.4336,-12.4375,-12.4411,-12.4448,-12.4478,-12.4481,-12.45,-12.4517,-12.4525,-12.4528,-12.4517,-12.4503,-12.4498,-12.4514,-12.4542,-12.4567,-12.4598,-12.462,-12.4642,-12.4664,-12.4681,-12.4698,-12.4714,-12.4727,-12.4728,-12.4753,-12.4781,-12.482,-12.4856,-12.4906,-12.4956,-12.4992,-12.5014,-12.5064,-12.5123,-12.5173,-12.5223,-12.5275,-12.5317,-12.5367,-12.542,-12.5467,-12.5511,-12.5564,-12.5606,-12.565,-12.5686,-12.5731,-12.5767,-12.5811,-12.5848,-12.5886,-12.5923,-12.5959,-12.5989,-12.6017,-12.6053,-12.6086,-12.6131,-12.6173,-12.6214,-12.6259,-12.6311,-12.6367,-12.6425,-12.6484,-12.6539,-12.6611,-12.6675,-12.6725,-12.6759,-12.6786,-12.682,-12.6853,-12.6856,-12.6861,-12.6641,-12.6444,-12.6247,-12.6237,-12.6235,-12.6228,-12.6219,-12.6211,-12.6213,-12.6211,-12.6212,-12.6209,-12.6213,-12.6205,-12.6192,-12.6177,-12.6163,-12.6155,-12.6146,-12.6141,-12.6141,-12.6149,-12.615,-12.6152,-12.6155,-12.6156,-12.6144,-12.614,-12.6126,-12.6117,-12.6102,-12.611,-12.612,-12.6126,-12.6129,-12.6128,-12.6118,-12.6107,-12.6098,-12.6087,-12.6074,-12.6066,-12.6058,-12.6051,-12.6053,-12.6049,-12.6042,-12.6042,-12.6038,-12.6032,-12.6031,-12.6026,-12.6022,-12.6022,-12.602,-12.6019,-12.6014,-12.6006,-12.6008,-12.6013,-12.6007,-12.6011,-12.6011,-12.6005,-12.6005,-12.6001,-12.6001,-12.6006,-12.6006,-12.6,-12.5992,-12.5987,-12.597,-12.5967,-12.596,-12.5953,-12.5948,-12.594,-12.5937,-12.5928,-12.592,-12.5914,-12.591,-12.5904,-12.5892,-12.5885,-12.5878,-12.587,-12.5866,-12.5859,-12.5942,-12.5954,-12.6049,-12.6144,-12.6147,-12.6239,-12.6334,-12.6429,-12.651,-12.6524,-12.6693,-12.6861,-12.703,-12.7198,-12.7449,-12.7701,-12.7952,-12.8211,-12.847,-12.8729,-12.8988,-12.9247,-12.9506,-12.9537,-12.9745,-12.9989,-13.0174,-13.0363,-13.0659,-13.0669,-13.0956,-13.1252,-13.1548,-13.1845,-13.2106,-13.2141,-13.2438,-13.259,-13.2742,-13.2894,-13.3046,-13.3198,-13.335,-13.3502,-13.3654,-13.3806,-13.3958,-13.4109,-13.4261,-13.4413,-13.4565,-13.4717,-13.4869,-13.5021,-13.5172,-13.5323,-13.5471,-13.577,-13.6069,-13.6368,-13.6667,-13.6966,-13.7162,-13.7265,-13.7564,-13.7834,-13.7864,-13.7868,-13.7874,-13.7886,-13.8163,-13.8194,-13.8462,-13.8717,-13.8972,-13.9227,-13.9482,-13.9645,-13.9737,-13.9948,-13.9989,-14.0016,-14.0253,-14.034,-14.0531,-14.0565,-14.1062,-14.1545,-14.1564,-14.1565,-14.1642,-14.1792,-14.1792,-14.1825,-14.1851,-14.1875,-14.2384,-14.2385,-14.2414,-14.2414,-14.2775,-14.2846,-14.2904,-14.2941,-14.3236,-14.3286,-14.33,-14.3301,-14.3301,-14.3302,-14.3513,-14.3515,-14.3681,-14.3683,-14.3799,-14.3998,-14.408,-14.4166,-14.428,-14.4356,-14.4472,-14.456,-14.4666,-14.4693,-14.4703,-14.4707,-14.4713,-14.4739,-14.4755,-14.4776,-14.4799,-14.4814,-14.4824,-14.4834,-14.4837,-14.4846,-14.4858,-14.4879,-14.4882,-14.4895,-14.4904,-14.4913,-14.4937,-14.4939,-14.494,-14.494,-14.4976,-14.4985,-14.5025,-14.513,-14.519,-14.5247,-14.5316,-14.5413,-14.5484,-14.5566,-14.5604,-14.5658,-14.5674,-14.5717,-14.5733,-14.5738,-14.5821,-14.5894,-14.5993,-14.6053,-14.6076,-14.6271,-14.6329,-14.6514,-14.6577,-14.6571,-14.6528,-14.6485,-14.6428,-14.6326,-14.6249,-14.6251,-14.6269,-14.6277,-14.6284,-14.6299,-14.6314,-14.6328,-14.6343,-14.6358,-14.6373,-14.6388,-14.6403,-14.6393,-14.6376,-14.6361,-14.6317,-14.6304,-14.6292,-14.6267,-14.6235,-14.6201,-14.6161,-14.6131,-14.611,-14.6085,-14.6055,-14.6029,-14.6003,-14.5975,-14.5961,-14.5961,-14.5954,-14.5954,-14.5953,-14.595,-14.5946,-14.5946,-14.5941,-14.5938,-14.5935,-14.5932,-14.5933,-14.5932,-14.5924,-14.592,-14.5919,-14.5918,-14.5912,-14.5915,-14.5909,-14.5908,-14.5903,-14.59,-14.5901,-14.5901,-14.5901,-14.5901,-14.5898,-14.5894,-14.5889,-14.5891,-14.5889,-14.5886,-14.5884,-14.5882,-14.5884,-14.5881,-14.5876,-14.5871,-14.5868,-14.5868,-14.5866,-14.5868,-14.5867,-14.586,-14.5855,-14.5854,-14.5848,-14.5846,-14.5845,-14.5847,-14.5846,-14.5845,-14.5848,-14.5852,-14.5859,-14.5863,-14.5864,-14.5867,-14.587,-14.5876,-14.588,-14.5884,-14.5889,-14.5894,-14.5893,-14.5899,-14.5907,-14.5907,-14.5912,-14.5914,-14.592,-14.592,-14.5926,-14.593,-14.5935,-14.5933,-14.5931,-14.5928,-14.5924,-14.5923,-14.5921,-14.5921,-14.5913,-14.5913,-14.5905,-14.5904,-14.5899,-14.5901,-14.5894,-14.5893,-14.5888,-14.5882,-14.5879,-14.5878,-14.5873,-14.5867,-14.5869,-14.5862,-14.5863,-14.5857,-14.5853,-14.5849,-14.5848,-14.5847,-14.5845,-14.5842,-14.5833,-14.5829,-14.5826,-14.5823,-14.5821,-14.5818,-14.5815,-14.5818,-14.5822,-14.5828,-14.5835,-14.5843,-14.5853,-14.5856,-14.5862,-14.5871,-14.5879,-14.5889,-14.59,-14.5907,-14.5922,-14.5934,-14.5945,-14.5959,-14.5972,-14.5985,-14.5992,-14.6005,-14.6013,-14.6033,-14.6043,-14.6053,-14.609,-14.6123,-14.6186,-14.6225,-14.6251,-14.6309,-14.6338,-14.6379,-14.6422,-14.6452,-14.6489,-14.6516,-14.655,-14.6587,-14.6631,-14.6673,-14.6709,-14.675,-14.679,-14.682,-14.687,-14.6919,-14.6943,-14.6973,-14.6999,-14.7009,-14.7016,-14.7042,-14.7069,-14.7092,-14.7095,-14.7105,-14.7121,-14.7135,-14.7145,-14.7161,-14.7175,-14.7186,-14.7201,-14.7239,-14.7303,-14.7362,-14.7441,-14.7497,-14.754,-14.7593,-14.7643,-14.7706,-14.7757,-14.7792,-14.7845,-14.7906,-14.7961,-14.7993,-14.8016,-14.8023,-14.8047,-14.8065,-14.8082,-14.8093,-14.8107,-14.8125,-14.8139,-14.8158,-14.8174,-14.82,-14.8213,-14.8234,-14.8407,-14.8392,-14.8353,-14.8317,-14.8281,-14.8245,-14.8209,-14.8184,-14.8153,-14.8131,-14.8106,-14.8086,-14.8061,-14.8039,-14.8023,-14.8006,-14.7986,-14.7973,-14.7964,-14.7945,-14.7928,-14.7925,-14.7917,-14.7914,-14.7911,-14.7731,-14.7548,-14.7367,-14.7186,-14.7178,-14.7136,-14.7103,-14.707,-14.7036,-14.7,-14.6959,-14.6925,-14.6884,-14.6842,-14.68,-14.6759,-14.6717,-14.6675,-14.6628,-14.6578,-14.6531,-14.6481,-14.6425,-14.637,-14.6311,-14.625,-14.6186,-14.6123,-14.605,-14.5981,-14.5925,-14.5867,-14.5809,-14.5753,-14.5695,-14.5645,-14.5586,-14.5536,-14.5492,-14.5445,-14.5392,-14.5348,-14.5306,-14.5267,-14.5225,-14.5184,-14.5153,-14.5109,-14.5073,-14.5042,-14.5006,-14.4992,-14.4967,-14.4945,-14.4917,-14.4892,-14.487,-14.4845,-14.4823,-14.4806,-14.4798,-14.4781,-14.477,-14.4761,-14.475,-14.4742,-14.4734,-14.4723,-14.472,-14.4661,-14.4606,-14.455,-14.4506,-14.4463,-14.4456,-14.4425,-14.44,-14.4373,-14.4342,-14.4306,-14.427,-14.4237,-14.4234,-14.4181,-14.4145,-14.41,-14.4056,-14.4014,-14.3975,-14.3934,-14.3892,-14.3848,-14.3803,-14.3759,-14.3714,-14.3667,-14.3623,-14.357,-14.3523,-14.347,-14.342,-14.337,-14.3311,-14.3267,-14.3217,-14.3161,-14.3106,-14.3048,-14.2992,-14.2936,-14.2886,-14.2839,-14.2789,-14.2725,-14.2664,-14.2606,-14.2542,-14.2473,-14.24,-14.2361,-14.2325,-14.23,-14.2273,-14.2239,-14.2206,-14.2205,-14.2203,-14.217,-14.2128,-14.2086,-14.2045,-14.1998,-14.1948,-14.1906,-14.1856,-14.18,-14.1753,-14.1689,-14.1634,-14.1575,-14.1514,-14.1456,-14.1386,-14.1323,-14.1259,-14.1231,-14.1225,-14.1198,-14.1164,-14.1175,-14.115,-14.1114,-14.1067,-14.1017,-14.097,-14.0923,-14.087,-14.0817,-14.0761,-14.0698,-14.0639,-14.0567,-14.0506,-14.0448,-14.0398,-14.0345,-14.0298,-14.0256,-14.0214,-14.0181,-14.0139,-14.0098,-14.0056,-14.0014,-13.9992,-13.9964,-13.9923,-13.9867,-13.9798,-13.9725,-13.9675,-13.9623,-13.9581,-13.9539,-13.95,-13.9436,-13.9373,-13.93,-13.9236,-13.9181,-13.9131,-13.9086,-13.9042,-13.8998,-13.8961,-13.8925,-13.8889,-13.8853,-13.8823,-13.88,-13.8775,-13.8761,-13.8745,-13.8728,-13.8698,-13.8686,-13.8678,-13.8589,-13.8501,-13.8499,-13.8492,-13.8403,-13.8334,-13.8267,-13.8211,-13.817,-13.8135,-13.8114,-13.8053,-13.8006,-13.7986,-13.7978,-13.7978,-13.7978,-13.8006,-13.8048,-13.8136,-13.8198,-13.8267,-13.8342,-13.8417,-13.8492,-13.8575,-13.8664,-13.8781,-13.8911,-13.8986,-13.9056,-13.9175,-13.93,-13.9425,-13.9539,-13.9631,-13.972,-13.9814,-13.992,-13.9992,-14.0042,-14.0173,-14.0295,-14.0406,-14.0481,-14.0584,-14.0686,-14.0823,-14.0948,-14.107,-14.115,-14.1198,-14.1253,-14.1331,-14.1417,-14.1509,-14.1603,-14.1714,-14.1817,-14.1858,-14.1934,-14.207,-14.2139,-14.2261,-14.2364,-14.2473,-14.257,-14.2673,-14.2789,-14.2911,-14.3036,-14.3106,-14.3242,-14.337,-14.3503,-14.3617,-14.372,-14.3831,-14.3928,-14.4,-14.4092,-14.4178,-14.427,-14.4336,-14.4406,-14.4459,-14.4517,-14.4611,-14.4748,-14.4881,-14.4992,-14.5003,-14.5111,-14.5203,-14.5306,-14.5381,-14.5456,-14.547,-14.5525,-14.5586,-14.5653,-14.5736,-14.5825,-14.592,-14.6025,-14.6148,-14.6264,-14.64,-14.647,-14.6598,-14.6723,-14.6834,-14.6942,-14.7045,-14.7139,-14.7228,-14.7311,-14.7381,-14.7436,-14.7489,-14.7545,-14.7592,-14.7642,-14.775,-14.7886,-14.8011,-14.8128,-14.8228,-14.8325,-14.8409,-14.847,-14.855,-14.8684,-14.8756,-14.8889,-14.9017,-14.9142,-14.9245,-14.9334,-14.9436,-14.9539,-14.9628,-14.9717,-14.98,-14.9889,-14.997,-14.9992,-14.9997,-15.0061,-15.0142,-15.0223,-15.0311,-15.0403,-15.0486,-15.0553,-15.0628,-15.0703,-15.0778,-15.0848,-15.0909,-15.0978,-15.1039,-15.1109,-15.1217,-15.1339,-15.145,-15.1548,-15.16,-15.165,-15.1703,-15.1742,-15.1767,-15.1839,-15.1923,-15.1998,-15.2017,-15.2025,-15.203,-15.203,-15.2034,-15.2034,-15.2048,-15.2039,-15.2067,-15.2081,-15.212,-15.217,-15.2239,-15.23,-15.2381,-15.2464,-15.2548,-15.2636,-15.2739,-15.2842,-15.295,-15.3075,-15.3198,-15.3334,-15.3377,-15.3459,-15.3595,-15.3717,-15.3834,-15.3936,-15.4056,-15.4173,-15.4286,-15.4398,-15.45,-15.4617,-15.4725,-15.4828,-15.4923,-15.4992,-15.5048,-15.5159,-15.5242,-15.5294,-15.5303,-15.5386,-15.5473,-15.5564,-15.567,-15.5803,-15.5939,-15.6009,-15.6084,-15.6206,-15.6323,-15.6448,-15.655,-15.6659,-15.6761,-15.6842,-15.6925,-15.7017,-15.712,-15.7239,-15.7364,-15.7503,-15.757,-15.7639,-15.7714,-15.7784,-15.7906,-15.8028,-15.8086,-15.8078,-15.8078,-15.8078,-15.8078,-15.8086,-15.8086,-15.8086,-15.8086,-15.8086,-15.8086,-15.8089,-15.8261,-15.8439,-15.8617,-15.8795,-15.8943,-15.8948,-15.8973,-15.9153,-15.9331,-15.9503,-15.9681,-15.9859,-15.9981,-15.9992,-16.0161,-16.0339,-16.0509,-16.0686,-16.0867,-16.1045,-16.1223,-16.14,-16.1573,-16.175,-16.1928,-16.2106,-16.2284,-16.2461,-16.2642,-16.2811,-16.2989,-16.3167,-16.3348,-16.3525,-16.3703,-16.3873,-16.4053,-16.4231,-16.4409,-16.4586,-16.4764,-16.4942,-16.4992,-16.5114,-16.5292,-16.547,-16.5648,-16.5828,-16.5836,-16.6006,-16.6175,-16.6353,-16.6534,-16.6711,-16.6889,-16.6998,-16.7048,-16.7075,-16.7075,-16.7142,-16.7203,-16.7286,-16.7328,-16.7306,-16.7324,-16.7315,-16.7285,-16.7285,-16.7291,-16.7302,-16.7309,-16.7315,-16.7326,-16.7342,-16.7359,-16.7373,-16.738,-16.7382,-16.7392,-16.74,-16.7412,-16.7446,-16.7459,-16.7461,-16.7471,-16.7469,-16.7469,-16.7465,-16.7454,-16.745,-16.7447,-16.7446,-16.743,-16.743,-16.7421,-16.7421,-16.7429,-16.743,-16.7421,-16.7421,-16.743,-16.7438,-16.7447,-16.7446,-16.7455,-16.7455,-16.7463,-16.7463,-16.7471,-16.7472,-16.748,-16.7488,-16.7488,-16.7504,-16.7513,-16.753,-16.753,-16.7555,-16.7547,-16.7554,-16.7555,-16.7563,-16.7563,-16.7571,-16.7571,-16.758,-16.758,-16.7588,-16.7588,-16.7596,-16.7605,-16.7613,-16.7613,-16.7621,-16.7621,-16.763,-16.763,-16.7663,-16.7671,-16.7705,-16.7713,-16.773,-16.7746,-16.7755,-16.7755,-16.7771,-16.7779,-16.7796,-16.7796,-16.7805,-16.7804,-16.7813,-16.7813,-16.7821,-16.783,-16.7837,-16.7838,-16.7846,-16.7846,-16.7855,-16.7855,-16.7863,-16.7863,-16.7871,-16.7871,-16.788,-16.7871,-16.788,-16.788,-16.7871,-16.7871,-16.7855,-16.7855,-16.7834,-16.7817,-16.7805,-16.7813,-16.7822,-16.7821,-16.7813,-16.7813,-16.7788,-16.7788,-16.778,-16.778,-16.7788,-16.7788,-16.778,-16.778,-16.7771,-16.7772,-16.7763,-16.7755,-16.7746,-16.7738,-16.7705,-16.7688,-16.768,-16.7671,-16.7671,-16.7663,-16.7663,-16.7655,-16.7663,-16.7663,-16.768,-16.768,-16.7688,-16.7721,-16.7738,-16.7738,-16.7709,-16.77,-16.7684,-16.7667,-16.7609,-16.7542,-16.7534,-16.75,-16.7484,-16.7425,-16.7384,-16.7317,-16.7313,-16.7358,-16.7417,-16.7434,-16.7442,-16.7467,-16.7484,-16.7501,-16.755,-16.7575,-16.7592,-16.7626,-16.7634,-16.7651,-16.7684,-16.7705,-16.7705,-16.7684,-16.7659,-16.7655,-16.7701,-16.773,-16.773,-16.7746,-16.7746,-16.7734,-16.7721,-16.7721,-16.7751,-16.7763,-16.7763,-16.7771,-16.7771,-16.7779,-16.7813,-16.783,-16.783,-16.7846,-16.7846,-16.7838,-16.7821,-16.7804,-16.7796,-16.7796,-16.7805,-16.7813,-16.7838,-16.7838,-16.7821,-16.7809,-16.7776,-16.7763,-16.7763,-16.7771,-16.7758,-16.7717,-16.7709,-16.7646,-16.7638,-16.7647,-16.7621,-16.7613,-16.7575,-16.7551,-16.7518,-16.7488,-16.7488,-16.7496,-16.7492,-16.7443,-16.7434,-16.7417,-16.7326,-16.7292,-16.7276,-16.7238,-16.7238,-16.7271,-16.7288,-16.7284,-16.7259,-16.7234,-16.7171,-16.7155,-16.7138,-16.7138,-16.7096,-16.7096,-16.7126,-16.7159,-16.7176,-16.7192,-16.7301,-16.7359,-16.7434,-16.7488,-16.7488,-16.7505,-16.7513,-16.7521,-16.7521,-16.7513,-16.7534,-16.7546,-16.7529,-16.753,-16.7534,-16.7542,-16.7567,-16.7605,-16.7621,-16.7638,-16.7646,-16.7663,-16.7663,-16.768,-16.768,-16.7696,-16.7705,-16.7705,-16.773,-16.7738,-16.7738,-16.7746,-16.7746,-16.7755,-16.7763,-16.7796,-16.7797,-16.775,-16.7717,-16.7701,-16.7676,-16.7642,-16.7625,-16.7601,-16.7584,-16.7571,-16.7621,-16.7634,-16.765,-16.7671,-16.7675,-16.7688,-16.7692,-16.7709,-16.773,-16.773,-16.7721,-16.773,-16.7713,-16.7696,-16.7705,-16.7696,-16.7696,-16.7688,-16.7688,-16.768,-16.7671,-16.7663,-16.7655,-16.7646,-16.7646,-16.7646,-16.7634,-16.7631,-16.7605,-16.7605,-16.7588,-16.7596,-16.7613,-16.7613,-16.7629,-16.7655,-16.7655,-16.7647,-16.7646,-16.7646,-16.763,-16.763,-16.7584,-16.7534,-16.7526,-16.745,-16.74,-16.738,-16.738,-16.7396,-16.7413,-16.7409,-16.7396,-16.7396,-16.7384,-16.7334,-16.7284,-16.7217,-16.7209,-16.7176,-16.7134,-16.7109,-16.7101,-16.705,-16.7042,-16.7017,-16.6992,-16.6967,-16.6951,-16.693,-16.693,-16.6942,-16.7009,-16.7017,-16.705,-16.7055,-16.7038,-16.7038,-16.7046,-16.7055,-16.7063,-16.7071,-16.7067,-16.7055,-16.7043,-16.7026,-16.7,-16.6992,-16.6925,-16.6859,-16.6784,-16.6709,-16.6683,-16.6671,-16.6663,-16.6646,-16.6646,-16.6675,-16.67,-16.6709,-16.6742,-16.6747,-16.6726,-16.6709,-16.6684,-16.6638,-16.6634,-16.6617,-16.6592,-16.6584,-16.6525,-16.6459,-16.6451,-16.64,-16.6392,-16.6329,-16.6329,-16.6359,-16.6368,-16.6371,-16.6359,-16.6342,-16.6317,-16.6305,-16.6305,-16.6313,-16.6313,-16.6338,-16.6338,-16.6346,-16.6355,-16.638,-16.638,-16.6388,-16.6388,-16.6396,-16.6396,-16.6394,-16.6371,-16.6351,-16.633,-16.6313,-16.6305,-16.6305,-16.6288,-16.6288,-16.6275,-16.6242,-16.6238,-16.6246,-16.6284,-16.6296,-16.6297,-16.6288,-16.628,-16.6226,-16.6209,-16.6167,-16.6159,-16.6117,-16.6109,-16.6083,-16.6075,-16.6042,-16.5988,-16.5988,-16.6013,-16.6009,-16.5975,-16.5942,-16.5917,-16.59,-16.588,-16.588,-16.5855,-16.5855,-16.5813,-16.5809,-16.5788,-16.5788,-16.5734,-16.5709,-16.5692,-16.565,-16.5588,-16.5588,-16.5576,-16.5572,-16.5563,-16.5563,-16.5555,-16.5546,-16.553,-16.5529,-16.5567,-16.5588,-16.5588,-16.5609,-16.5655,-16.5663,-16.5676,-16.5696,-16.5697,-16.5634,-16.5609,-16.5567,-16.5534,-16.5451,-16.5442,-16.5438,-16.5409,-16.5276,-16.5234,-16.5225,-16.52,-16.5184,-16.5146,-16.5146,-16.5075,-16.5063,-16.5063,-16.5013,-16.5013,-16.4976,-16.4934,-16.4921,-16.4922,-16.493,-16.4938,-16.4955,-16.4955,-16.4938,-16.4921,-16.4921,-16.493,-16.4917,-16.4876,-16.4847,-16.4855,-16.4871,-16.4871,-16.488,-16.488,-16.4901,-16.4922,-16.4905,-16.4905,-16.493,-16.493,-16.4905,-16.4904,-16.4892,-16.4876,-16.4846,-16.4846,-16.4863,-16.4859,-16.4775,-16.4746,-16.4746,-16.4759,-16.478,-16.4788,-16.4804,-16.4809,-16.4825,-16.485,-16.4855,-16.4846,-16.4834,-16.4817,-16.4792,-16.4788,-16.4788,-16.4805,-16.4805,-16.4826,-16.4851,-16.4868,-16.488,-16.488,-16.4817,-16.4776,-16.4768,-16.475,-16.473,-16.473,-16.4717,-16.47,-16.4675,-16.465,-16.4646,-16.4663,-16.468,-16.4671,-16.4651,-16.4642,-16.4621,-16.4613,-16.4592,-16.4584,-16.4571,-16.4596,-16.4596,-16.4638,-16.4638,-16.4659,-16.4701,-16.473,-16.473,-16.4738,-16.4738,-16.4775,-16.4792,-16.4805,-16.4805,-16.4796,-16.4797,-16.478,-16.478,-16.4805,-16.4813,-16.4834,-16.4884,-16.4905,-16.4905,-16.4884,-16.4842,-16.4813,-16.4796,-16.4776,-16.4717,-16.4642,-16.4592,-16.4567,-16.4542,-16.4525,-16.4513,-16.4522,-16.4538,-16.4538,-16.4546,-16.4542,-16.453,-16.4438,-16.4438,-16.4409,-16.4401,-16.4384,-16.4375,-16.4346,-16.4346,-16.4338,-16.4338,-16.4317,-16.4284,-16.4255,-16.4255,-16.4242,-16.4234,-16.4205,-16.4209,-16.4217,-16.423,-16.423,-16.428,-16.4271,-16.425,-16.4159,-16.415,-16.4134,-16.4126,-16.4059,-16.4017,-16.4009,-16.3934,-16.3918,-16.3892,-16.3792,-16.3776,-16.3767,-16.3717,-16.3709,-16.3651,-16.36,-16.3555,-16.3522,-16.3521,-16.3509,-16.3484,-16.3471,-16.3479,-16.3488,-16.3505,-16.353,-16.3521,-16.3521,-16.3505,-16.3496,-16.3497,-16.3467,-16.3434,-16.34,-16.3392,-16.335,-16.3292,-16.3284,-16.3267,-16.3259,-16.3234,-16.3225,-16.3134,-16.3125,-16.3068,-16.3042,-16.2992,-16.2959,-16.2942,-16.2909,-16.29,-16.2851,-16.2817,-16.2809,-16.275,-16.2684,-16.2663,-16.2663,-16.2688,-16.2688,-16.2676,-16.2642,-16.2617,-16.2584,-16.2576,-16.2542,-16.2509,-16.25,-16.2446,-16.2446,-16.2421,-16.2421,-16.2375,-16.2326,-16.2317,-16.2293,-16.2197,-16.2197,-16.2226,-16.2284,-16.2367,-16.24,-16.2417,-16.2467,-16.2484,-16.2501,-16.2542,-16.2592,-16.2675,-16.2684,-16.2709,-16.2751,-16.2801,-16.2826,-16.2851,-16.2867,-16.2909,-16.2918,-16.2938,-16.2955,-16.2955,-16.2976,-16.3009,-16.3017,-16.3051,-16.3067,-16.3159,-16.3176,-16.3192,-16.3225,-16.3233,-16.3259,-16.3284,-16.3309,-16.3326,-16.3359,-16.3421,-16.3429,-16.3463,-16.3463,-16.3471,-16.3471,-16.348,-16.348,-16.3505,-16.3567,-16.3609,-16.3642,-16.365,-16.3717,-16.3725,-16.3767,-16.3775,-16.3826,-16.3842,-16.3884,-16.3896,-16.3876,-16.3871,-16.3888,-16.3888,-16.3909,-16.3926,-16.3934,-16.3967,-16.3984,-16.4009,-16.4033,-16.4042,-16.4109,-16.4167,-16.4184,-16.4209,-16.425,-16.43,-16.4317,-16.4367,-16.4409,-16.4434,-16.4442,-16.4475,-16.4496,-16.4496,-16.4521,-16.4542,-16.4546,-16.4538,-16.4538,-16.4496,-16.4505,-16.4521,-16.4517,-16.4492,-16.4475,-16.4446,-16.4438,-16.4451,-16.4476,-16.4492,-16.4517,-16.4526,-16.453,-16.4555,-16.4555,-16.4584,-16.4625,-16.4634,-16.4709,-16.4717,-16.4809,-16.4817,-16.4842,-16.485,-16.4876,-16.4884,-16.4934,-16.4942,-16.4976,-16.5026,-16.5108,-16.5117,-16.5142,-16.5167,-16.5209,-16.5217,-16.5234,-16.5276,-16.5325,-16.5338,-16.5338,-16.533,-16.5321,-16.5309,-16.5293,-16.5246,-16.5259,-16.5301,-16.5305,-16.5296,-16.5296,-16.5288,-16.5288,-16.5242,-16.5184,-16.5175,-16.5158,-16.5117,-16.5076,-16.5067,-16.5038,-16.5105,-16.5113,-16.5134,-16.515,-16.5184,-16.52,-16.5217,-16.5247,-16.5251,-16.5263,-16.5263,-16.5254,-16.5255,-16.523,-16.523,-16.5251,-16.5275,-16.5301,-16.5317,-16.5321,-16.5321,-16.533,-16.533,-16.5346,-16.5347,-16.5355,-16.5355,-16.5371,-16.5371,-16.5396,-16.5396,-16.5413,-16.5413,-16.5384,-16.5371,-16.5363,-16.5347,-16.5346,-16.5372,-16.5372,-16.5492,-16.55,-16.5559,-16.5592,-16.5634,-16.5655,-16.5659,-16.5676,-16.5688,-16.5692,-16.5726,-16.575,-16.5775,-16.5784,-16.5817,-16.5851,-16.5867,-16.5913,-16.5922,-16.5921,-16.5913,-16.595,-16.5967,-16.6025,-16.6097,-16.6105,-16.6105,-16.6113,-16.6113,-16.6126,-16.615,-16.618,-16.6213,-16.6238,-16.6263,-16.6263,-16.6284,-16.633,-16.6329,-16.6397,-16.6397,-16.6413,-16.6413,-16.643,-16.6434,-16.6451,-16.6534,-16.6567,-16.6584,-16.6626,-16.6659,-16.6667,-16.67,-16.6755,-16.6755,-16.6788,-16.6788,-16.6771,-16.6771,-16.6738,-16.6738,-16.6746,-16.6755,-16.6763,-16.6763,-16.6705,-16.6696,-16.6663,-16.6613,-16.6646,-16.6634,-16.66,-16.6596,-16.6642,-16.6646,-16.6705,-16.6705,-16.6671,-16.6663,-16.665,-16.6647,-16.6667,-16.6701,-16.6718,-16.6776,-16.6834,-16.6847,-16.6792,-16.675,-16.6725,-16.6713,-16.6713,-16.6696,-16.6696,-16.668,-16.668,-16.6663,-16.6663,-16.6654,-16.6667,-16.6717,-16.6742,-16.6855,-16.6855,-16.6876,-16.6892,-16.6951,-16.6963,-16.6963,-16.6955,-16.6955,-16.6926,-16.6909,-16.6884,-16.6834,-16.6805,-16.6796,-16.6755,-16.6755,-16.6746,-16.6747,-16.6755,-16.6792,-16.6859,-16.69,-16.6917,-16.6921,-16.6901,-16.6875,-16.6863,-16.6917,-16.695,-16.6976,-16.7051,-16.7083,-16.7167,-16.7217,-16.7259,-16.7309,-16.7338,-16.7363,-16.7363,-16.7388,-16.7396,-16.7405,-16.7405,-16.7413,-16.7446,-16.7438,-16.7438,-16.7446,-16.7455,-16.748,-16.748,-16.7488,-16.75,-16.7517,-16.7534,-16.755,-16.7584,-16.7625,-16.768,-16.768,-16.7719,-16.778,-16.7805,-16.7805,-16.7822,-16.783,-16.7838,-16.7846,-16.7863,-16.7863,-16.788,-16.788,-16.7921,-16.7921,-16.793,-16.793,-16.7921,-16.7921,-16.7913,-16.7913,-16.7921,-16.7921,-16.7913,-16.7905,-16.7896,-16.7896,-16.7888,-16.7888,-16.7871,-16.7871,-16.7855,-16.7855,-16.7846,-16.7838,-16.783,-16.783,-16.7813,-16.7813,-16.7805,-16.7796,-16.7755,-16.7755,-16.7742,-16.7725,-16.7688,-16.768,-16.7679,-16.7667,-16.7634,-16.7597,-16.7584,-16.7521,-16.7513,-16.7513,-16.7501,-16.7459,-16.7442,-16.738,-16.7355,-16.7291,-16.7198,-16.7156,-16.7153,-16.7143,-16.7109,-16.7106,-16.685,-16.6692,-16.6692,-16.6692,-16.6678,-16.6675,-16.6662,-16.6655,-16.6639,-16.6656,-16.6665,-16.6666,-16.6665,-16.6445,-16.6378,-16.632,-16.6264,-16.6211,-16.617,-16.6106,-16.6036,-16.5967,-16.5903,-16.5831,-16.5767,-16.5706,-16.565,-16.5586,-16.5517,-16.5461,-16.5417,-16.5364,-16.5314,-16.5264,-16.5214,-16.5142,-16.5089,-16.5034,-16.4992,-16.4978,-16.4923,-16.4859,-16.48,-16.4748,-16.4692,-16.4628,-16.4578,-16.4517,-16.4453,-16.4384,-16.432,-16.4261,-16.4214,-16.4161,-16.4111,-16.4056,-16.3984,-16.3936,-16.3911,-16.3906,-16.3809,-16.3636,-16.3464,-16.3289,-16.312,-16.2945,-16.2773,-16.2603,-16.2428,-16.2328,-16.2259,-16.2084,-16.2011,-16.1834,-16.1653,-16.1634,-16.1592,-16.1542,-16.1503,-16.1461,-16.1406,-16.1356,-16.1309,-16.1261,-16.1206,-16.1142,-16.1078,-16.1014,-16.095,-16.0881,-16.0811,-16.0748,-16.0684,-16.062,-16.0564,-16.0506,-16.0448,-16.0392,-16.0342,-16.0289,-16.0234,-16.0184,-16.0142,-16.0098,-16.0053,-15.9998,-15.9992,-15.9945,-15.9895,-15.9845,-15.9795,-15.9745,-15.97,-15.965,-15.9614,-15.9564,-15.9392,-15.9217,-15.9039,-15.8878,-15.8864,-15.8792,-15.8614,-15.8439,-15.8261,-15.8086,-15.7906,-15.7728,-15.7553,-15.7375,-15.72,-15.702,-15.6845,-15.6831,-15.6656,-15.6484,-15.6311,-15.6139,-15.5967,-15.5792,-15.562,-15.5445,-15.5319,-15.5318,-15.5273,-15.51,-15.4992,-15.4925,-15.4753,-15.4578,-15.4406,-15.4236,-15.4186,-15.4061,-15.3945,-15.3775,-15.3606,-15.3587,-15.3434,-15.3336,-15.3164,-15.2986,-15.2814,-15.2684,-15.2514,-15.2342,-15.2173,-15.2117,-15.1939,-15.1761,-15.1586,-15.1409,-15.1364,-15.1231,-15.1056,-15.0878,-15.0703,-15.0525,-15.0348,-15.017,-14.9995,-14.9992,-14.9931,-14.9756,-14.9578,-14.9423,-14.94,-14.9225,-14.9048,-14.8873,-14.8736,-14.8692,-14.8517,-14.8342,-14.8164,-14.7986,-14.7809,-14.7634,-14.7459,-14.7281,-14.7103,-14.6925,-14.675,-14.6575,-14.6395,-14.622,-14.6045,-14.5867,-14.5689,-14.5511,-14.5336,-14.5161,-14.4992,-14.4981,-14.4806,-14.4631,-14.4453,-14.4275,-14.4098,-14.3923,-14.3745,-14.357,-14.3411,-14.3409,-14.3392,-14.3214,-14.3039,-14.2864,-14.2684,-14.2509,-14.2331,-14.2156,-14.1978,-14.18,-14.1625,-14.1589,-14.145,-14.1273,-14.1095,-14.0917,-14.0742,-14.0567,-14.0386,-14.0211,-14.0036,-13.9992,-13.9978,-13.98,-13.9625,-13.9448,-13.9273,-13.9095,-13.8906,-13.8567,-13.8389,-13.8211,-13.8034,-13.7859,-13.7684,-13.7503,-13.7328,-13.7153,-13.7123,-13.6975,-13.6798,-13.6623,-13.6445,-13.627,-13.6258,-13.6198,-13.602,-13.5845,-13.5667,-13.5517,-13.5489,-13.5314,-13.5136,-13.4992,-13.4961,-13.4781,-13.4606,-13.4468,-13.4428,-13.4253,-13.4075,-13.3898,-13.372,-13.3596,-13.3559,-13.3523,-13.345,-13.3406,-13.3361,-13.3278,-13.3225,-13.3186,-13.3148,-13.3092,-13.3069,-13.3036,-13.2984,-13.2906,-13.2861,-13.2817,-13.277,-13.2723,-13.2681,-13.2631,-13.2567,-13.2542,-13.2523,-13.2484,-13.2423,-13.235,-13.23,-13.2242,-13.2195,-13.2159,-13.2125,-13.2056,-13.2006,-13.197,-13.1925,-13.1898,-13.1859,-13.1817,-13.1759,-13.1734,-13.1692,-13.1636,-13.1578,-13.1509,-13.1453,-13.1414,-13.1386,-13.1353,-13.1317,-13.1292,-13.1261,-13.1253,-13.1225,-13.1175,-13.1117,-13.1053,-13.1006,-13.0948,-13.0884,-13.082,-13.0764,-13.07,-13.065,-13.0603,-13.0573,-13.0536,-13.0506,-13.0478,-13.0467,-13.0461,-13.045,-13.0428,-13.0411,-13.0409,-13.0423,-13.0434,-13.0461,-13.0481,-13.0506,-13.0525,-13.0548,-13.055,-13.057,-13.0606,-13.0639,-13.0673,-13.0686,-13.0684,-13.0678,-13.0656,-13.0653,-13.065,-13.0634,-13.0622,-13.0611,-13.0589,-13.0559,-13.0528,-13.0486,-13.0445,-13.04,-13.037,-13.03,-13.0245,-13.0181,-13.012,-13.0048,-12.9992,-12.9984,-12.9925,-12.9878,-12.9834,-12.9781,-12.9734,-12.9678,-12.9628,-12.9586,-12.9545,-12.9525,-12.952,-12.9509,-12.9503,-12.9501,-12.9492,-12.9486,-12.9475,-12.9461,-12.945,-12.9423,-12.9389,-12.9348,-12.9306,-12.9256,-12.9209,-12.9145,-12.9075,-12.9003,-12.8953,-12.8906,-12.887,-12.8836,-12.8809,-12.8778,-12.8742,-12.87,-12.8656,-12.8606,-12.855,-12.85,-12.845,-12.8392,-12.8328,-12.8292,-12.8298,-12.8319,-12.832,-12.8331,-12.8381,-12.8414,-12.8428,-12.8431,-12.8403,-12.8348,-12.8289,-12.8253,-12.8225,-12.8181,-12.8125,-12.8067,-12.7998,-12.7925,-12.7861,-12.7806,-12.7742,-12.7684,-12.7628,-12.7634,-12.7645,-12.7659,-12.7656,-12.7603,-12.7534,-12.747,-12.74,-12.7331,-12.7261,-12.7189,-12.7125,-12.7056,-12.6984,-12.692,-12.685,-12.6786,-12.6717,-12.6648,-12.6584,-12.6514,-12.645,-12.6378,-12.6323,-12.6286,-12.625,-12.622,-12.6192,-12.6173,-12.615,-12.6131,-12.6106,-12.6075,-12.6031,-12.597,-12.5911,-12.5875,-12.5873,-12.5895,-12.5906,-12.587,-12.5842,-12.582,-12.5781,-12.5717,-12.5675,-12.5634,-12.5592,-12.5548,-12.5498,-12.5442,-12.5392,-12.5345,-12.5303,-12.5253,-12.5203,-12.5156,-12.5106,-12.5073,-12.5017,-12.4992,-12.4961,-12.4903,-12.4842,-12.4784,-12.4725,-12.4678,-12.462,-12.4561,-12.4506,-12.445,-12.4398,-12.4356,-12.432,-12.427,-12.4228,-12.4164,-12.41,-12.4042,-12.3986,-12.3942,-12.3906,-12.3873,-12.3839,-12.3817,-12.3795,-12.3773,-12.375,-12.3728,-12.3698,-12.3661,-12.3617,-12.3575,-12.3548,-12.3523,-12.35,-12.347,-12.3448,-12.3434,-12.3386,-12.3339,-12.3289,-12.3248,-12.32,-12.315,-12.3103,-12.3048,-12.3009,-12.2998,-12.2948,-12.29,-12.2845,-12.2786,-12.2731,-12.2675,-12.2625,-12.2584,-12.2536,-12.2489,-12.2425,-12.2356,-12.2284,-12.2214,-12.2142,-12.2078,-12.2017,-12.1953,-12.1911,-12.1856,-12.1798,-12.175,-12.1695,-12.1645,-12.1595,-12.155,-12.1509,-12.1467,-12.1425,-12.1389,-12.1348,-12.1309,-12.1273,-12.1248,-12.1211,-12.1178,-12.1142,-12.1111,-12.107,-12.1034,-12.0978,-12.0909,-12.0845,-12.0773,-12.0703,-12.0631,-12.0567,-12.0511,-12.0461,-12.042,-12.0375,-12.0331,-12.0275,-12.0211,-12.0148,-12.0103,-12.0056,-12.0006,-11.9992,-11.9948,-11.9878,-11.9814,-11.975,-11.9681,-11.9625,-11.9584,-11.9542,-11.9498,-11.9453,-11.9411,-11.937,-11.932,-11.9278,-11.9228,-11.9189,-11.9148,-11.9106,-11.9073,-11.9042,-11.902,-11.8995,-11.8973,-11.8945,-11.8895,-11.8836,-11.8781,-11.8725,-11.8673,-11.8631,-11.8589,-11.8545,-11.8503,-11.845,-11.8409,-11.8348,-11.8289,-11.8234,-11.8178,-11.8123,-11.805,-11.7981,-11.7917,-11.7861,-11.7809,-11.7775,-11.7723,-11.7681,-11.7631,-11.7567,-11.7498,-11.7425,-11.7356,-11.73,-11.7242,-11.72,-11.7159,-11.7111,-11.7048,-11.6992,-11.6928,-11.6878,-11.6836,-11.6798,-11.6748,-11.6706,-11.6664,-11.6623,-11.6581,-11.6539,-11.6498,-11.6456,-11.6423,-11.6373,-11.6317,-11.6273,-11.6231,-11.6189,-11.6125,-11.6061,-11.6011,-11.5956,-11.59,-11.5839,-11.5775,-11.5711,-11.5648,-11.5586,-11.5514,-11.545,-11.5381,-11.5317,-11.5253,-11.5189,-11.512,-11.5056,-11.4992,-11.4986,-11.4923,-11.4853,-11.4789,-11.4725,-11.4661,-11.4603,-11.4548,-11.4492,-11.4436,-11.4378,-11.432,-11.4264,-11.4214,-11.4173,-11.4136,-11.4098,-11.4056,-11.4006,-11.395,-11.3898,-11.3842,-11.3792,-11.3731,-11.3723,-11.3714,-11.3725,-11.3734,-11.3742,-11.3745,-11.3739,-11.3728,-11.37,-11.3664,-11.3636,-11.3603,-11.3575,-11.3611,-11.3614,-11.3622,-11.3625,-11.3631,-11.3653,-11.3671,-11.3672,-11.3698,-11.3742,-11.3741,-11.3739,-11.3753,-11.3759,-11.3786,-11.3842,-11.3864,-11.3859,-11.3861,-11.3861,-11.3884,-11.3917,-11.3942,-11.395,-11.3967,-11.3995,-11.4036,-11.4081,-11.4123,-11.4123,-11.4106,-11.4106,-11.4136,-11.4189,-11.4225,-11.4256,-11.4284,-11.4298,-11.4309,-11.4331,-11.4373,-11.4436,-11.4478,-11.4509,-11.4489,-11.4439,-11.4392,-11.4336,-11.4286,-11.4236,-11.4195,-11.4161,-11.4156,-11.415,-11.4175,-11.42,-11.4234,-11.4256,-11.4284,-11.4292,-11.4273,-11.4253,-11.4217,-11.4216,-11.4206,-11.4195,-11.42,-11.4211,-11.4217,-11.4234,-11.4261,-11.4314,-11.4384,-11.4448,-11.4492,-11.4506,-11.4506,-11.4495,-11.4475,-11.4461,-11.4425,-11.4384,-11.4314,-11.4278,-11.4275,-11.4284,-11.4278,-11.4264,-11.4209,-11.4153,-11.4089,-11.4034,-11.3975,-11.3914,-11.3864,-11.3823,-11.3823,-11.3859,-11.3875,-11.387,-11.3886,-11.3923,-11.3964,-11.4009,-11.4009,-11.3959,-11.3925,-11.3906,-11.3892,-11.3895,-11.3911,-11.392,-11.3923,-11.3914,-11.3911,-11.3925,-11.3942,-11.3998,-11.4048,-11.4111,-11.4156,-11.417,-11.4156,-11.4131,-11.41,-11.4081,-11.4095,-11.412,-11.4136,-11.4123,-11.4103,-11.4061,-11.4034,-11.4009,-11.3986,-11.3989,-11.4011,-11.4048,-11.4084,-11.4114,-11.4109,-11.4067,-11.4017,-11.3967,-11.3923,-11.3867,-11.3803,-11.3753,-11.3711,-11.3673,-11.3667,-11.3673,-11.3684,-11.3692,-11.37,-11.3711,-11.3717,-11.3734,-11.375,-11.3773,-11.3809,-11.3856,-11.3878,-11.3884,-11.3881,-11.3892,-11.3914,-11.3939,-11.3992,-11.4034,-11.4084,-11.4139,-11.4189,-11.4225,-11.4206,-11.4175,-11.4164,-11.4175,-11.4203,-11.4217,-11.4214,-11.4185,-11.4184,-11.4175,-11.4173,-11.4131,-11.4081,-11.4056,-11.4075,-11.4125,-11.4184,-11.4214,-11.4228,-11.4236,-11.4253,-11.4289,-11.4313,-11.4331,-11.4403,-11.4425,-11.4428,-11.44,-11.4367,-11.4336,-11.4339,-11.437,-11.4414,-11.4478,-11.4548,-11.4598,-11.4656,-11.4703,-11.4739,-11.477,-11.4806,-11.485,-11.4881,-11.4914,-11.4959,-11.4992,-11.5028,-11.507,-11.5128,-11.5159,-11.5164,-11.5181,-11.5203,-11.5234,-11.527,-11.5292,-11.5306,-11.5286,-11.5253,-11.5234,-11.5223,-11.5214,-11.5217,-11.5217,-11.5248,-11.5289,-11.5353,-11.542,-11.5461,-11.5506,-11.5534,-11.555,-11.5567,-11.5573,-11.5561,-11.5542,-11.5514,-11.5486,-11.5461,-11.5434,-11.5406,-11.537,-11.5359,-11.5381,-11.5409,-11.5445,-11.547,-11.5492,-11.5514,-11.5545,-11.5595,-11.565,-11.5686,-11.5717,-11.5748,-11.5784,-11.5828,-11.5864,-11.587,-11.5881,-11.5895,-11.5925,-11.5953,-11.597,-11.597,-11.5967,-11.597,-11.597,-11.6011,-11.6039,-11.6073,-11.6131,-11.6189,-11.6245,-11.6289,-11.6311,-11.6298,-11.6292,-11.6292,-11.6317,-11.6359,-11.6423,-11.6473,-11.6509,-11.6542,-11.6606,-11.667,-11.6725,-11.6789,-11.6853,-11.6903,-11.6939,-11.6961,-11.6984,-11.7006,-11.7042,-11.7075,-11.7123,-11.7186,-11.725,-11.732,-11.7384,-11.742,-11.7434,-11.7439,-11.7459,-11.7478,-11.7503,-11.7531,-11.7559,-11.7581,-11.7592,-11.7598,-11.7586,-11.7592,-11.7634,-11.7689,-11.7756,-11.7806,-11.7839,-11.7864,-11.7886,-11.7911,-11.7939,-11.7973,-11.8014,-11.8056,-11.8106,-11.8156,-11.8192,-11.8223,-11.8239,-11.8261,-11.8275,-11.8278,-11.8289,-11.8289,-11.8325,-11.8384,-11.8431,-11.8489,-11.8539,-11.8584,-11.8611,-11.8642,-11.8673,-11.8714,-11.8786,-11.8842,-11.8873,-11.8873,-11.8873,-11.8853,-11.8834,-11.8811,-11.8792,-11.8789,-11.8784,-11.8764,-11.8731,-11.8695,-11.8667,-11.8653,-11.8656,-11.8673,-11.8673,-11.8692,-11.8711,-11.8742,-11.8786,-11.885,-11.892,-11.8975,-11.902,-11.9056,-11.9073,-11.9089,-11.9103,-11.9116,-11.9136,-11.9175,-11.9228,-11.9264,-11.93,-11.9342,-11.9375,-11.9409,-11.9459,-11.9517,-11.9586,-11.9639,-11.9681,-11.9711,-11.9734,-11.9734,-11.9761,-11.9792,-11.982,-11.9859,-11.99,-11.9931,-11.9961,-11.9989,-11.9992,-12.002,-12.0036,-12.005,-12.0067,-12.0098,-12.0125,-12.0167,-12.0211,-12.0256,-12.03,-12.0336,-12.035,-12.0367,-12.0348,-12.0334,-12.0336,-12.0359,-12.0389,-12.0425,-12.0461,-12.0498,-12.0534,-12.0564,-12.0586,-12.0609,-12.0617,-12.0611,-12.0614,-12.0636,-12.0673,-12.0717,-12.0739,-12.0753,-12.0758,-12.081,-12.0811,-12.0814,-12.0784,-12.0734,-12.0692,-12.0645,-12.0595,-12.0553,-12.0506,-12.0453,-12.0406,-12.0364,-12.0323,-12.0273,-12.0231,-12.0184,-12.0142,-12.0092,-12.005,-12.0009,-11.9992,-11.9967,-11.9925,-11.9889,-11.9848,-11.9814,-11.9781,-11.9739,-11.9703,-11.967,-11.9636,-11.9609,-11.9578,-11.9559,-11.9539,-11.952,-11.9509,-11.9498,-11.9484,-11.947,-11.9464,-11.9453,-11.9439,-11.9428,-11.9417,-11.9417,-11.9411,-11.9414,-11.9409,-11.9411,-11.9411,-11.9414,-11.9436,-11.9467,-11.9503,-11.9545,-11.9581,-11.962,-11.9656,-11.97,-11.9736,-11.9781,-11.9823,-11.9859,-11.99,-11.9939,-11.9975,-11.9992,-12.0003,-12.0042,-12.007,-12.0086,-12.0089,-12.0089,-12.0092,-12.0095,-12.01,-12.0103,-12.0125,-12.0142,-12.0136,-12.0123,-12.0098,-12.007,-12.005,-12.0053,-12.0061,-12.0056,-12.0036,-12.0017,-11.9992,-11.9975,-11.9939,-11.9903,-11.987,-11.985,-11.9845,-11.9856,-11.9878,-11.99,-11.9909,-11.99,-11.9881,-11.9861,-11.9842,-11.9817,-11.9798,-11.9792,-11.9786,-11.9786,-11.9798,-11.9814,-11.9836,-11.9859,-11.9881,-11.9898,-11.9911,-11.9928,-11.9959,-11.9986,-11.9992,-12.0009,-12.0034,-12.0056,-12.0078,-12.0095,-12.0117,-12.0134,-12.0148,-12.017,-12.0186,-12.0209,-12.0239,-12.027,-12.03,-12.0336,-12.04,-12.0456,-12.05,-12.057,-12.0628,-12.0686,-12.0728,-12.0778,-12.0823,-12.0867,-12.0909,-12.0945,-12.0975,-12.0989,-12.102,-12.105,-12.1086,-12.1123,-12.1117,-12.1061,-12.1003,-12.0961,-12.097,-12.0986,-12.0959,-12.0917,-12.0898,-12.0911,-12.0942,-12.0978,-12.1031,-12.1092,-12.115,-12.1214,-12.1245,-12.1267,-12.1295,-12.1353,-12.1411,-12.1453,-12.15,-12.1545,-12.1581,-12.1611,-12.1648,-12.1695,-12.1761,-12.1811,-12.1856,-12.1898,-12.1928,-12.2056,-12.2045,-12.2042,-12.2041,-12.2039,-12.2031,-12.2031,-12.2039,-12.204,-12.2042,-12.2036,-12.2017,-12.2003,-12.1998,-12.2014,-12.2045,-12.2061,-12.2067,-12.2081,-12.2084,-12.2086,-12.2117,-12.2159,-12.2198,-12.2225,-12.2242,-12.2249,-12.2249,-12.2253,-12.2245,-12.222,-12.2192,-12.2173,-12.2159,-12.2156,-12.2161,-12.2181,-12.2159,-12.2103,-12.2053,-12.1981,-12.1925,-12.1884,-12.1861,-12.185,-12.1867,-12.1878,-12.1878,-12.1856,-12.1823,-12.1789,-12.1784,-12.1778,-12.175,-12.1717,-12.1667,-12.1603,-12.1561,-12.1528,-12.1498,-12.1478,-12.1467,-12.147,-12.1486,-12.1514,-12.1559,-12.1609,-12.1678,-12.1736,-12.1778,-12.1803,-12.1814,-12.177,-12.1714,-12.1678,-12.1653,-12.1634,-12.1642,-12.1706,-12.1739,-12.1781,-12.1811,-12.182,-12.1873,-12.1928,-12.1984,-12.205,-12.2092,-12.2117,-12.2139,-12.2148,-12.2136,-12.2123,-12.2153,-12.22,-12.2253,-12.2311,-12.2361,-12.2398,-12.2434,-12.2442,-12.2442,-12.2441,-12.2439,-12.2439,-12.2475,-12.252,-12.2567,-12.2611,-12.2656,-12.2706,-12.2745,-12.2784,-12.2828,-12.2864,-12.2903,-12.2945,-12.2984,-12.302,-12.3056,-12.3092,-12.3125,-12.3162,-12.3189,-12.322,-12.3256,-12.33,-12.3348,-12.3392,-12.3445,-12.3503,-12.3559,-12.3603,-12.3659],"lat":[14.8396,14.8415,14.8427,14.8454,14.8474,14.8499,14.8527,14.856,14.8593,14.8627,14.8665,14.8699,14.8746,14.8779,14.8813,14.8849,14.8874,14.8907,14.894,14.8974,14.9007,14.9046,14.9052,14.9096,14.9149,14.9213,14.9279,14.9335,14.939,14.9454,14.9507,14.9549,14.9596,14.9638,14.9685,14.9732,14.9779,14.9835,14.9888,14.9943,14.9996,14.9999,15.0046,15.0085,15.0118,15.0152,15.0171,15.019,15.0199,15.0204,15.0221,15.0235,15.0254,15.0274,15.0293,15.0318,15.0338,15.0357,15.0377,15.0402,15.0421,15.0446,15.0474,15.0507,15.0532,15.0565,15.0593,15.0627,15.066,15.0693,15.0727,15.0768,15.0807,15.084,15.0882,15.0907,15.0935,15.096,15.0985,15.1004,15.1018,15.1029,15.104,15.1054,15.1052,15.1043,15.1021,15.0985,15.0943,15.0907,15.0877,15.0871,15.087,15.0543,15.0273,15.0004,14.9973,14.9955,14.9933,14.9916,14.9905,14.9882,14.9872,14.9861,14.985,14.9837,14.9829,14.9824,14.9822,14.9813,14.9808,14.98,14.9794,14.9783,14.9775,14.9766,14.9759,14.9751,14.9742,14.973,14.9726,14.9724,14.9719,14.9703,14.969,14.9682,14.9675,14.9666,14.9656,14.965,14.9647,14.964,14.9634,14.9632,14.9628,14.9623,14.9608,14.9593,14.958,14.9568,14.956,14.9552,14.953,14.9521,14.9512,14.9502,14.9492,14.9478,14.9466,14.9454,14.944,14.943,14.9412,14.9401,14.9386,14.9375,14.9363,14.9356,14.9341,14.9333,14.9321,14.9308,14.9294,14.9279,14.9263,14.9245,14.9229,14.922,14.9206,14.9194,14.9179,14.9169,14.9151,14.9139,14.9129,14.912,14.911,14.91,14.9091,14.9079,14.9071,14.9057,14.9052,14.8798,14.8761,14.8474,14.8188,14.818,14.7902,14.7616,14.7329,14.7087,14.7043,14.6816,14.6589,14.6362,14.6135,14.6025,14.5915,14.5804,14.5767,14.573,14.5694,14.5657,14.562,14.5583,14.5592,14.5654,14.5726,14.5764,14.5803,14.5847,14.5849,14.5892,14.5936,14.5981,14.6025,14.6064,14.6069,14.6114,14.5999,14.5884,14.577,14.5655,14.554,14.5425,14.5311,14.5196,14.5081,14.4967,14.4852,14.4737,14.4622,14.4507,14.4393,14.4278,14.4163,14.4048,14.4052,14.4055,14.4062,14.4069,14.4076,14.4083,14.4089,14.4094,14.4096,14.4103,14.4109,14.4109,14.4109,14.411,14.411,14.4116,14.4117,14.4123,14.4246,14.437,14.4493,14.4617,14.4696,14.4741,14.484,14.4859,14.4864,14.4911,14.493,14.4973,14.4984,14.5145,14.5304,14.531,14.531,14.5335,14.5385,14.5385,14.5395,14.5404,14.5412,14.5578,14.5578,14.5588,14.5588,14.5709,14.573,14.5747,14.5758,14.5845,14.586,14.5864,14.5864,14.5864,14.5864,14.5838,14.5838,14.5823,14.5823,14.5812,14.5801,14.58,14.5799,14.5791,14.5802,14.5807,14.579,14.5787,14.5787,14.5779,14.5776,14.577,14.5743,14.5724,14.5702,14.5675,14.5665,14.5658,14.5655,14.5654,14.5651,14.5647,14.564,14.5639,14.5634,14.563,14.5626,14.562,14.562,14.5619,14.5619,14.5621,14.5621,14.5623,14.5627,14.5632,14.5628,14.5623,14.562,14.5604,14.561,14.561,14.5627,14.5638,14.5667,14.5678,14.568,14.5717,14.5742,14.5745,14.5746,14.5747,14.5738,14.5735,14.5745,14.5754,14.5751,14.5733,14.572,14.5703,14.5672,14.5649,14.5623,14.542,14.5282,14.5153,14.4885,14.4617,14.4349,14.4081,14.3813,14.3546,14.3278,14.301,14.2986,14.2974,14.2958,14.2919,14.2906,14.2895,14.2868,14.2832,14.2812,14.2769,14.274,14.2721,14.2699,14.2676,14.2649,14.2623,14.2601,14.2592,14.2581,14.2552,14.2537,14.2514,14.2484,14.2457,14.2431,14.2402,14.2381,14.2348,14.231,14.2282,14.2253,14.2243,14.2213,14.2192,14.2173,14.2154,14.213,14.2091,14.2074,14.2057,14.2026,14.2012,14.1989,14.1971,14.1958,14.1944,14.1928,14.1909,14.1884,14.1859,14.1834,14.1812,14.1784,14.1768,14.1748,14.1721,14.1692,14.1663,14.1638,14.162,14.161,14.1581,14.1551,14.1527,14.1493,14.1465,14.144,14.1416,14.1396,14.1385,14.1378,14.1367,14.1351,14.1325,14.13,14.1278,14.1251,14.1234,14.1214,14.1181,14.1172,14.1143,14.1118,14.1102,14.1076,14.1058,14.1036,14.1021,14.0995,14.0976,14.0959,14.0948,14.0909,14.0892,14.0873,14.0854,14.0837,14.0824,14.0794,14.0773,14.0743,14.0723,14.0695,14.0676,14.0629,14.0617,14.0591,14.0572,14.0537,14.0517,14.0494,14.0462,14.0436,14.0408,14.039,14.0344,14.0307,14.0269,14.0249,14.0221,14.0192,14.0162,14.0127,14.0111,14.0097,14.0079,14.0066,14.0041,14.0022,14.0014,14.0007,14.0001,13.9984,13.9957,13.994,13.991,13.9887,13.9863,13.9855,13.9834,13.9811,13.9787,13.9761,13.9724,13.9698,13.9663,13.963,13.9599,13.9564,13.9522,13.95,13.9473,13.9436,13.9412,13.9382,13.9374,13.9366,13.9335,13.9312,13.9286,13.9263,13.9242,13.9211,13.9188,13.9166,13.9137,13.9111,13.9093,13.9077,13.9051,13.9024,13.9,13.8974,13.8939,13.8913,13.8905,13.8887,13.8858,13.8836,13.8791,13.8714,13.8651,13.862,13.8602,13.8549,13.8488,13.8436,13.8428,13.8398,13.8366,13.8331,13.8304,13.827,13.8232,13.8208,13.8174,13.8162,13.8152,13.8145,13.814,13.813,13.8126,13.812,13.8117,13.8115,13.8112,13.8109,13.8109,13.8107,13.8102,13.8098,13.8095,13.8094,13.8078,13.8064,13.8046,13.8029,13.8004,13.7979,13.7962,13.7935,13.7913,13.7881,13.7849,13.7829,13.7673,13.766,13.7629,13.7596,13.7563,13.7529,13.7496,13.7449,13.7407,13.736,13.7313,13.7265,13.7218,13.7171,13.7118,13.7063,13.7007,13.6954,13.6893,13.6838,13.6782,13.6713,13.6652,13.6585,13.6529,13.6432,13.6329,13.6227,13.6121,13.6115,13.6143,13.6179,13.6215,13.6252,13.6288,13.6318,13.6352,13.6382,13.641,13.644,13.6468,13.6499,13.6529,13.6552,13.6574,13.6596,13.6618,13.6635,13.6649,13.6665,13.6674,13.6685,13.6693,13.6696,13.6699,13.6688,13.6677,13.6665,13.6652,13.664,13.6624,13.661,13.6593,13.6565,13.6546,13.6529,13.6502,13.6477,13.6443,13.6418,13.639,13.6365,13.634,13.6307,13.6279,13.6246,13.6234,13.6213,13.6165,13.6127,13.6079,13.6032,13.5985,13.5935,13.5882,13.5821,13.5765,13.5704,13.5643,13.5579,13.5518,13.5457,13.5396,13.534,13.5329,13.5318,13.5304,13.5279,13.5263,13.526,13.5221,13.5174,13.5132,13.5093,13.506,13.5027,13.4997,13.4993,13.4974,13.494,13.4915,13.489,13.4863,13.4829,13.4804,13.4779,13.4752,13.4727,13.4699,13.4674,13.4654,13.4629,13.461,13.4593,13.4574,13.4554,13.4535,13.4538,13.456,13.4582,13.4599,13.4615,13.4629,13.4646,13.4663,13.4685,13.4707,13.4729,13.4738,13.4749,13.4763,13.4774,13.4777,13.4779,13.4807,13.4843,13.4885,13.4927,13.4963,13.4996,13.4997,13.4999,13.5035,13.5063,13.5093,13.5121,13.5143,13.5168,13.5196,13.5218,13.5235,13.5257,13.5265,13.5282,13.5299,13.5307,13.5324,13.5327,13.5335,13.5343,13.5352,13.5352,13.5379,13.5415,13.5396,13.5424,13.546,13.5488,13.551,13.554,13.5563,13.5585,13.5602,13.5618,13.5627,13.5629,13.5632,13.5627,13.5615,13.5596,13.5577,13.5557,13.5574,13.5602,13.5638,13.5665,13.5696,13.5727,13.5754,13.5764,13.5777,13.5804,13.5821,13.5824,13.5827,13.5807,13.5788,13.5763,13.5738,13.5704,13.5699,13.5693,13.5696,13.5693,13.5679,13.5663,13.5635,13.561,13.5582,13.5552,13.5518,13.5485,13.5452,13.541,13.5363,13.5315,13.526,13.5207,13.5152,13.511,13.5091,13.5077,13.5035,13.4997,13.4996,13.4993,13.494,13.4874,13.4802,13.4715,13.4621,13.4579,13.4554,13.4474,13.4382,13.426,13.4193,13.4127,13.406,13.3888,13.3793,13.3674,13.3593,13.3521,13.346,13.3393,13.3327,13.3274,13.3221,13.3188,13.3174,13.3174,13.3174,13.3188,13.3193,13.3215,13.3207,13.316,13.3107,13.306,13.3021,13.301,13.3002,13.2988,13.3002,13.3007,13.294,13.2902,13.2868,13.2854,13.2846,13.2868,13.2807,13.2713,13.2635,13.2568,13.2513,13.246,13.2421,13.2379,13.234,13.2331,13.2313,13.2307,13.2302,13.2321,13.2346,13.2379,13.2427,13.246,13.2435,13.2421,13.2393,13.2393,13.2379,13.2368,13.2374,13.2402,13.2427,13.246,13.2507,13.256,13.2613,13.2668,13.2721,13.2788,13.286,13.2949,13.3027,13.306,13.3049,13.3049,13.3066,13.3068,13.3093,13.3135,13.3215,13.3282,13.334,13.3354,13.3415,13.3488,13.356,13.3621,13.3588,13.354,13.3502,13.3474,13.3449,13.344,13.3435,13.3449,13.3468,13.3493,13.3521,13.356,13.3602,13.3654,13.3715,13.3788,13.3868,13.3954,13.404,13.414,13.4227,13.4254,13.424,13.4254,13.4274,13.4315,13.4349,13.4407,13.4488,13.4549,13.4549,13.4549,13.4549,13.456,13.4574,13.4602,13.464,13.4682,13.4715,13.4768,13.4821,13.4882,13.4927,13.4979,13.4993,13.4996,13.5035,13.5093,13.5146,13.5193,13.5246,13.5307,13.5374,13.544,13.5507,13.5574,13.5646,13.5721,13.5799,13.5879,13.5954,13.5935,13.5907,13.5874,13.5827,13.5746,13.5654,13.5574,13.5524,13.5493,13.5435,13.5374,13.5307,13.5188,13.5068,13.4998,13.4996,13.494,13.4807,13.4688,13.4615,13.4507,13.4388,13.4288,13.4202,13.4135,13.4054,13.3993,13.3935,13.3882,13.3827,13.3782,13.374,13.3707,13.3688,13.3668,13.3654,13.3656,13.366,13.366,13.3682,13.3693,13.3715,13.3735,13.3749,13.3774,13.3802,13.3835,13.3854,13.3888,13.3921,13.396,13.3968,13.3974,13.394,13.3888,13.3819,13.3807,13.3749,13.3693,13.364,13.3607,13.3593,13.3582,13.3574,13.3574,13.3582,13.3607,13.3621,13.3649,13.3682,13.3715,13.3674,13.3615,13.356,13.3521,13.3493,13.3482,13.3468,13.346,13.346,13.3454,13.3454,13.3468,13.3474,13.3307,13.3135,13.2968,13.2793,13.2621,13.2449,13.2446,13.2274,13.2102,13.1927,13.176,13.1593,13.1593,13.1602,13.1593,13.1602,13.1602,13.1602,13.1602,13.1607,13.1613,13.1613,13.1621,13.1621,13.1627,13.1627,13.1627,13.1627,13.1627,13.1621,13.1627,13.1627,13.1627,13.1627,13.1635,13.1635,13.1635,13.1627,13.1627,13.1635,13.1635,13.1635,13.1635,13.164,13.1635,13.1635,13.1635,13.164,13.164,13.164,13.1635,13.1635,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.164,13.1646,13.164,13.164,13.164,13.164,13.1607,13.1513,13.1407,13.1307,13.1232,13.1152,13.1093,13.1027,13.0921,13.0899,13.0883,13.0813,13.0806,13.0788,13.0774,13.0754,13.0738,13.0732,13.0727,13.0727,13.0726,13.073,13.0731,13.0735,13.0735,13.0734,13.0701,13.068,13.0657,13.0648,13.0636,13.0635,13.0617,13.0559,13.0527,13.0521,13.0429,13.0396,13.0346,13.0338,13.0213,13.0205,12.9996,12.9987,12.9871,12.9862,12.9796,12.9788,12.9746,12.9738,12.9712,12.9704,12.9662,12.9654,12.9496,12.9488,12.9471,12.9429,12.9404,12.9371,12.9354,12.9338,12.9287,12.9279,12.9271,12.9221,12.9212,12.9162,12.9154,12.9129,12.9121,12.9096,12.9087,12.9062,12.9054,12.9004,12.8996,12.8971,12.8962,12.8871,12.8862,12.8821,12.8754,12.8712,12.8646,12.8579,12.8546,12.8455,12.8446,12.8413,12.8379,12.8329,12.8296,12.8196,12.8188,12.8137,12.8129,12.8096,12.8088,12.8005,12.7997,12.7996,12.7927,12.7921,12.7912,12.7871,12.7862,12.7837,12.7829,12.7779,12.777,12.7763,12.7755,12.7721,12.7713,12.7696,12.7671,12.7654,12.7642,12.7642,12.7654,12.7688,12.7696,12.7738,12.7746,12.7762,12.7804,12.7829,12.7838,12.7896,12.7904,12.7938,12.7946,12.7996,12.8004,12.8054,12.8062,12.8129,12.8138,12.8221,12.8271,12.8296,12.8363,12.8371,12.8388,12.8387,12.8371,12.8363,12.8354,12.8313,12.8287,12.8263,12.8254,12.8121,12.8088,12.8021,12.7992,12.7992,12.8009,12.8009,12.805,12.805,12.8058,12.8059,12.8067,12.8116,12.8142,12.8142,12.8129,12.8117,12.8084,12.8066,12.8067,12.805,12.805,12.8033,12.8025,12.8009,12.8,12.8,12.7992,12.7992,12.7975,12.7954,12.7912,12.7892,12.7892,12.7879,12.7859,12.783,12.7821,12.7804,12.7771,12.7758,12.7754,12.7746,12.7742,12.7729,12.7629,12.7621,12.7588,12.7571,12.7537,12.7504,12.7487,12.7471,12.7429,12.7421,12.7388,12.7379,12.7362,12.7296,12.7287,12.7154,12.7096,12.7079,12.7046,12.7034,12.7034,12.7046,12.7079,12.7096,12.71,12.7142,12.7142,12.7204,12.7229,12.7238,12.7254,12.7279,12.7316,12.7325,12.735,12.7387,12.7404,12.7413,12.7458,12.7409,12.7409,12.7425,12.7425,12.74,12.7392,12.7396,12.7404,12.7437,12.7471,12.7492,12.7467,12.7458,12.7387,12.7354,12.7346,12.7338,12.7296,12.7279,12.7275,12.7292,12.7291,12.7309,12.7358,12.7392,12.7391,12.7338,12.7329,12.7313,12.7279,12.7271,12.7246,12.7237,12.7233,12.7221,12.7187,12.7154,12.715,12.7167,12.7175,12.7146,12.7104,12.7095,12.7021,12.6987,12.6971,12.6946,12.6929,12.6904,12.6888,12.6829,12.677,12.6754,12.672,12.6712,12.6646,12.6638,12.6612,12.6571,12.6529,12.6509,12.6508,12.6517,12.655,12.6558,12.6584,12.6583,12.6591,12.6579,12.6538,12.6508,12.6525,12.6513,12.6475,12.6479,12.6492,12.6492,12.6471,12.6454,12.6446,12.6379,12.6363,12.6329,12.6321,12.6313,12.6287,12.6279,12.6246,12.6238,12.6187,12.6179,12.6121,12.6112,12.5964,12.5938,12.5917,12.592,12.5954,12.5979,12.5979,12.5937,12.5913,12.5895,12.5879,12.5854,12.5771,12.5764,12.5763,12.5746,12.5729,12.5712,12.5667,12.5642,12.5634,12.5633,12.5675,12.5687,12.5704,12.5738,12.5754,12.58,12.5796,12.5754,12.5734,12.5733,12.5792,12.5783,12.5791,12.5792,12.5825,12.5825,12.5833,12.5842,12.585,12.585,12.5867,12.5884,12.5884,12.5904,12.5912,12.5925,12.5925,12.5933,12.5934,12.5938,12.5954,12.5996,12.6004,12.6045,12.6054,12.6121,12.6125,12.6113,12.6084,12.6075,12.6075,12.6067,12.6042,12.5975,12.5975,12.6033,12.6042,12.6054,12.6079,12.6096,12.6113,12.6142,12.6142,12.6133,12.6134,12.6137,12.6158,12.6167,12.6167,12.6129,12.6109,12.6109,12.6125,12.6125,12.6183,12.6217,12.6217,12.6267,12.6267,12.6321,12.6329,12.6358,12.6358,12.6371,12.6375,12.6358,12.6358,12.637,12.6446,12.6454,12.6496,12.6546,12.6571,12.6579,12.6621,12.6679,12.6713,12.6721,12.6821,12.6829,12.6938,12.6969,12.6988,12.6993,12.6976,12.6937,12.693,12.6904,12.6887,12.6862,12.685,12.685,12.6846,12.6829,12.6825,12.6812,12.6696,12.6687,12.6646,12.6583,12.6575,12.6575,12.6567,12.6559,12.655,12.655,12.6559,12.6567,12.6621,12.6629,12.6654,12.6675,12.6642,12.6625,12.6625,12.6634,12.6654,12.6662,12.6696,12.6738,12.6796,12.6809,12.6821,12.6829,12.6892,12.69,12.6917,12.6925,12.697,12.7012,12.7033,12.6971,12.6963,12.6929,12.6921,12.6879,12.6854,12.6838,12.6834,12.6863,12.6913,12.6917,12.6871,12.6845,12.6842,12.6821,12.6763,12.6691,12.6684,12.665,12.6641,12.6583,12.6584,12.6571,12.655,12.6467,12.6416,12.6417,12.6392,12.6392,12.6354,12.6346,12.6275,12.6271,12.6262,12.6213,12.6204,12.6166,12.6158,12.6171,12.6213,12.6221,12.6271,12.6304,12.6395,12.643,12.6446,12.6488,12.6496,12.6516,12.6533,12.6563,12.6613,12.6629,12.6754,12.6762,12.6779,12.6809,12.6812,12.6838,12.6896,12.6921,12.6938,12.6954,12.7046,12.7059,12.7059,12.7079,12.7096,12.7121,12.715,12.7158,12.7137,12.7113,12.7092,12.7087,12.7046,12.7038,12.7016,12.7025,12.7025,12.6995,12.6987,12.6942,12.6942,12.6959,12.6954,12.6937,12.6921,12.6896,12.6867,12.6866,12.6858,12.6837,12.6813,12.6758,12.6734,12.6725,12.6725,12.6705,12.6671,12.6658,12.6658,12.6684,12.6709,12.6696,12.6679,12.6638,12.6596,12.6567,12.6567,12.6587,12.6613,12.6634,12.6634,12.6621,12.6588,12.658,12.6529,12.6521,12.6509,12.6508,12.6546,12.6604,12.6613,12.6629,12.6675,12.6675,12.6663,12.6637,12.6629,12.6613,12.6596,12.6546,12.6504,12.6438,12.6417,12.64,12.6379,12.6321,12.63,12.63,12.6263,12.6221,12.62,12.6167,12.6167,12.6134,12.6125,12.6109,12.6109,12.6121,12.6171,12.6188,12.6204,12.6212,12.6234,12.6204,12.6104,12.6096,12.6067,12.6067,12.6042,12.6042,12.6013,12.5996,12.5987,12.5971,12.5942,12.5942,12.5987,12.6012,12.6017,12.6,12.5987,12.5975,12.5975,12.5963,12.5912,12.5855,12.5821,12.58,12.58,12.5792,12.5792,12.5784,12.5783,12.575,12.575,12.5675,12.5667,12.5667,12.5617,12.5617,12.5608,12.56,12.5592,12.5591,12.5617,12.5671,12.5738,12.5796,12.5809,12.5817,12.5837,12.5913,12.5921,12.5962,12.5996,12.6012,12.6046,12.6046,12.6029,12.6004,12.5983,12.5984,12.6,12.6,12.6042,12.6058,12.6067,12.6067,12.6075,12.6075,12.6067,12.6067,12.6075,12.6075,12.6058,12.605,12.6025,12.6025,12.5992,12.5992,12.595,12.5933,12.5925,12.5925,12.5983,12.5979,12.5971,12.5946,12.5912,12.5892,12.5884,12.59,12.59,12.5909,12.5917,12.5942,12.5942,12.5996,12.6004,12.6029,12.6038,12.6083,12.6117,12.6125,12.6108,12.6104,12.6096,12.6067,12.6033,12.595,12.5933,12.5933,12.5908,12.5892,12.5892,12.5867,12.5842,12.585,12.5859,12.5859,12.5875,12.5875,12.5859,12.5858,12.5875,12.5883,12.5892,12.5904,12.5938,12.5954,12.5975,12.5983,12.5992,12.5992,12.5984,12.5883,12.5875,12.5867,12.5867,12.5858,12.5859,12.5833,12.5834,12.585,12.585,12.5771,12.5746,12.5679,12.5662,12.5654,12.5621,12.5612,12.5596,12.5546,12.5483,12.5459,12.545,12.5442,12.5442,12.545,12.545,12.5442,12.5442,12.5459,12.5459,12.5471,12.5475,12.5487,12.5513,12.5529,12.555,12.555,12.5567,12.5584,12.5583,12.5609,12.5616,12.5625,12.5625,12.5592,12.5592,12.5625,12.5633,12.5659,12.5659,12.5684,12.5717,12.5717,12.5725,12.5725,12.5704,12.5679,12.5638,12.5625,12.5646,12.5654,12.5671,12.5721,12.5746,12.5762,12.5767,12.575,12.575,12.5771,12.5788,12.5792,12.5809,12.5809,12.5833,12.5833,12.5846,12.5862,12.5871,12.59,12.59,12.5892,12.5892,12.59,12.59,12.5909,12.5908,12.5917,12.5917,12.5925,12.5925,12.5934,12.5942,12.5966,12.5967,12.5959,12.5958,12.5942,12.5942,12.5934,12.5934,12.59,12.5875,12.5862,12.5829,12.5821,12.5771,12.5758,12.575,12.5746,12.5734,12.5734,12.5704,12.5696,12.5671,12.5662,12.5562,12.5517,12.5517,12.5525,12.5525,12.5558,12.5559,12.555,12.5546,12.5487,12.5454,12.5434,12.5425,12.5425,12.545,12.545,12.5404,12.5383,12.5404,12.5421,12.5429,12.5446,12.5462,12.5479,12.5492,12.55,12.5525,12.5517,12.5521,12.5563,12.5571,12.5587,12.5604,12.5621,12.5629,12.5663,12.5687,12.5704,12.5754,12.5913,12.5929,12.5938,12.5942,12.5963,12.6004,12.6029,12.6054,12.6096,12.6104,12.6225,12.6225,12.6284,12.63,12.63,12.6279,12.6267,12.6267,12.6279,12.6333,12.6342,12.6358,12.6358,12.6367,12.6367,12.635,12.635,12.6304,12.6288,12.6262,12.6246,12.6209,12.62,12.62,12.6129,12.6113,12.6079,12.6071,12.6037,12.6025,12.6017,12.5987,12.5921,12.588,12.5854,12.5845,12.5833,12.5788,12.5779,12.5713,12.5704,12.5687,12.5679,12.567,12.5659,12.565,12.565,12.5667,12.5667,12.5617,12.5617,12.5608,12.5591,12.5546,12.5538,12.5504,12.5479,12.5454,12.5438,12.5405,12.5371,12.5363,12.5237,12.5229,12.5162,12.5096,12.5054,12.4996,12.4996,12.4963,12.495,12.495,12.4938,12.4925,12.4913,12.4854,12.4846,12.4813,12.4787,12.4784,12.4771,12.4767,12.4792,12.4792,12.4733,12.4733,12.4712,12.4683,12.4642,12.4642,12.4629,12.4613,12.4596,12.4579,12.4563,12.4546,12.4529,12.4512,12.4504,12.4492,12.4542,12.4542,12.4654,12.4662,12.4683,12.4683,12.4625,12.4629,12.4671,12.4679,12.4696,12.4725,12.4725,12.4742,12.4767,12.4796,12.4821,12.4871,12.4904,12.4913,12.4962,12.4979,12.5,12.5,12.4975,12.4975,12.4979,12.5,12.5008,12.5038,12.5042,12.5025,12.5009,12.5008,12.5025,12.5058,12.5058,12.5092,12.5092,12.5137,12.5162,12.5171,12.5196,12.5221,12.5229,12.5271,12.5287,12.5295,12.5312,12.5345,12.5354,12.5379,12.5412,12.5454,12.5471,12.5484,12.5484,12.5467,12.5467,12.5434,12.5417,12.5354,12.5346,12.5303,12.5238,12.5188,12.5171,12.5154,12.5104,12.5096,12.5054,12.5046,12.5021,12.4987,12.4971,12.4887,12.4846,12.4837,12.4813,12.4804,12.4737,12.4729,12.4654,12.4646,12.4571,12.4562,12.452,12.4513,12.4487,12.4479,12.4454,12.443,12.4421,12.4404,12.4371,12.4362,12.4296,12.4287,12.4263,12.4246,12.4229,12.422,12.4163,12.4088,12.4062,12.405,12.405,12.4021,12.4004,12.3979,12.3967,12.3967,12.3929,12.3909,12.3871,12.3854,12.382,12.3809,12.3808,12.38,12.3712,12.3663,12.3562,12.3458,12.3395,12.3391,12.3388,12.3354,12.3352,12.3595,12.3595,12.3612,12.3616,12.3654,12.3659,12.3675,12.3679,12.3662,12.3638,12.3608,12.36,12.3595,12.3595,12.3582,12.3571,12.356,12.3543,12.3515,12.3513,12.3515,12.3518,12.3529,12.3532,12.354,12.3552,12.3568,12.3577,12.3579,12.3571,12.3543,12.3527,12.3507,12.3488,12.3471,12.3474,12.349,12.3507,12.352,12.3524,12.354,12.3549,12.3565,12.3582,12.3599,12.3607,12.3629,12.364,12.3649,12.3654,12.3649,12.3638,12.3621,12.3602,12.3582,12.3571,12.3577,12.3599,12.364,12.3688,12.3743,12.3827,12.3915,12.3999,12.409,12.4174,12.4265,12.4349,12.4438,12.4487,12.4521,12.4604,12.4574,12.4513,12.4452,12.4468,12.4496,12.4518,12.4549,12.4579,12.4593,12.4618,12.464,12.4663,12.4679,12.4688,12.4699,12.471,12.4718,12.4721,12.4724,12.4721,12.4715,12.4713,12.4702,12.469,12.4677,12.4668,12.4649,12.4629,12.4618,12.4602,12.4574,12.4549,12.4524,12.4513,12.4511,12.4493,12.4474,12.4457,12.4438,12.4418,12.4393,12.4374,12.4357,12.4365,12.4393,12.4421,12.4449,12.4477,12.4479,12.4468,12.4454,12.4435,12.4421,12.4402,12.4388,12.4377,12.4357,12.4343,12.4324,12.431,12.4296,12.4304,12.4388,12.4463,12.454,12.4624,12.4699,12.4782,12.486,12.4935,12.4996,12.4996,12.5018,12.5093,12.5145,12.5177,12.5254,12.5329,12.5413,12.5488,12.551,12.5565,12.5674,12.5827,12.5977,12.5992,12.6121,12.6154,12.6202,12.6252,12.6293,12.6407,12.6552,12.6696,12.6843,12.6838,12.6838,12.6838,12.6832,12.6832,12.6832,12.6832,12.6832,12.6832,12.6827,12.6827,12.6827,12.6827,12.6827,12.6826,12.6824,12.6824,12.6824,12.6824,12.6824,12.6824,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.6815,12.681,12.681,12.681,12.681,12.681,12.681,12.681,12.681,12.6802,12.6802,12.6802,12.6802,12.6802,12.6802,12.6802,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6793,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6785,12.6782,12.6777,12.6777,12.6774,12.6774,12.6774,12.6774,12.6774,12.6774,12.6765,12.6765,12.6765,12.6765,12.6763,12.6763,12.6763,12.6763,12.6758,12.6757,12.6757,12.6763,12.6763,12.6763,12.6763,12.6771,12.676,12.6765,12.6765,12.6765,12.6763,12.6763,12.6763,12.6768,12.6768,12.6768,12.6768,12.6765,12.6765,12.6765,12.6771,12.677,12.6765,12.676,12.6752,12.6743,12.6736,12.6735,12.6727,12.6721,12.6714,12.6713,12.6704,12.6696,12.669,12.6688,12.6682,12.6674,12.6665,12.6657,12.6651,12.6649,12.6615,12.6549,12.6524,12.6496,12.6445,12.6427,12.6393,12.6407,12.6493,12.6499,12.6507,12.6529,12.6546,12.6521,12.6493,12.6515,12.6538,12.6568,12.659,12.6585,12.6538,12.649,12.6457,12.6452,12.6454,12.6435,12.6424,12.6446,12.6482,12.6515,12.6518,12.6502,12.6465,12.644,12.6399,12.6368,12.634,12.6343,12.6385,12.6413,12.6429,12.6432,12.6435,12.6435,12.6402,12.6363,12.6396,12.6432,12.6474,12.6435,12.6427,12.6402,12.6382,12.6371,12.6365,12.6346,12.6335,12.6329,12.6338,12.6354,12.6363,12.6385,12.6393,12.6354,12.6321,12.6279,12.624,12.6177,12.6115,12.6054,12.6007,12.5952,12.5896,12.584,12.5785,12.5743,12.5696,12.5652,12.5604,12.5554,12.5493,12.5443,12.5407,12.5371,12.5338,12.5296,12.5282,12.5246,12.5199,12.5143,12.5077,12.5021,12.4997,12.4974,12.4927,12.4885,12.4846,12.4818,12.4793,12.4768,12.4727,12.4729,12.4743,12.4754,12.4763,12.4765,12.4761,12.476,12.4746,12.4729,12.4702,12.4682,12.4663,12.4679,12.4702,12.4732,12.476,12.4807,12.4871,12.4927,12.4988,12.4996,12.5043,12.5107,12.5163,12.5218,12.5271,12.5315,12.5352,12.5379,12.5407,12.5429,12.5452,12.5463,12.5463,12.5465,12.5446,12.5427,12.5393,12.5354,12.5313,12.5274,12.524,12.5213,12.5188,12.5168,12.5157,12.5138,12.5118,12.5107,12.5102,12.5068,12.5018,12.4997,12.4996,12.4985,12.4963,12.4927,12.4871,12.4835,12.4796,12.4782,12.4771,12.4738,12.4699,12.4671,12.466,12.4646,12.4649,12.4652,12.4646,12.4635,12.4629,12.4618,12.4604,12.4543,12.4488,12.4432,12.4363,12.4343,12.4346,12.4354,12.4357,12.436,12.4363,12.4365,12.4374,12.4377,12.4377,12.4371,12.4374,12.4371,12.4371,12.4374,12.4368,12.4371,12.4379,12.4382,12.4371,12.4338,12.4304,12.4263,12.4224,12.4168,12.4121,12.4074,12.4027,12.3985,12.396,12.3954,12.394,12.391,12.3854,12.3804,12.3749,12.3715,12.3677,12.3632,12.3593,12.3588,12.3618,12.3646,12.3677,12.3699,12.3721,12.3735,12.3757,12.3779,12.381,12.3829,12.3852,12.3874,12.3896,12.3932,12.3949,12.3955,12.3963,12.3965,12.396,12.3949,12.3935,12.3915,12.3904,12.3893,12.3879,12.3868,12.3849,12.3824,12.379,12.3771,12.3743,12.374,12.3735,12.3721,12.371,12.3682,12.3649,12.3615,12.3577,12.3529,12.3482,12.3435,12.3385,12.3338,12.3299,12.3265,12.3238,12.3213,12.3171,12.3124,12.3077,12.3035,12.3013,12.3024,12.3046,12.3068,12.309,12.3118,12.314,12.3163,12.3185,12.3199,12.3216,12.3221,12.3243,12.3265,12.3282,12.3296,12.3313,12.3327,12.3349,12.3379,12.3402,12.3424,12.3432,12.3435,12.3435,12.3438,12.344,12.3449,12.3443,12.3452,12.3468,12.3482,12.3499,12.3518,12.3535,12.3557,12.3579,12.3602,12.3629,12.366,12.3688,12.3724,12.3752,12.3788,12.3824,12.3865,12.3902,12.3943,12.3977,12.4021,12.4049,12.4085,12.4099,12.4102,12.411,12.4113,12.4115,12.4115,12.411,12.4099,12.4079,12.4054,12.4027,12.4002,12.3988,12.3982,12.3979,12.3952,12.3932,12.3913,12.391,12.3902,12.3904,12.3899,12.3907,12.391,12.3924,12.3952,12.3982,12.401,12.4032,12.406,12.409,12.4113,12.414,12.4163,12.419,12.4221,12.4249,12.4285,12.4271,12.4224,12.4177,12.4129,12.4088,12.4068,12.4057,12.4043,12.4032,12.4013,12.3985,12.396,12.3932,12.3907,12.3888,12.386,12.3871,12.3885,12.3902,12.3915,12.3932,12.3932,12.3935,12.3929,12.3918,12.3899,12.3865,12.3846,12.3818,12.3799,12.3793,12.3796,12.3799,12.3802,12.3815,12.3832,12.386,12.3888,12.391,12.3918,12.3907,12.3902,12.3924,12.3952,12.3979,12.4002,12.4032,12.406,12.4088,12.4118,12.4146,12.4174,12.4204,12.4238,12.426,12.4263,12.4235,12.421,12.4182,12.419,12.4199,12.4221,12.4238,12.4252,12.426,12.4268,12.4277,12.4288,12.4296,12.4296,12.4304,12.4307,12.4315,12.4324,12.434,12.4343,12.4349,12.4351,12.4352,12.436,12.4363,12.4357,12.4352,12.4346,12.4335,12.4321,12.431,12.4296,12.4285,12.4274,12.426,12.424,12.4215,12.4179,12.4146,12.4121,12.4102,12.409,12.4071,12.4057,12.4038,12.4046,12.4074,12.4102,12.4163,12.4227,12.4288,12.4357,12.4418,12.4474,12.4515,12.4552,12.4593,12.4629,12.4752,12.4763,12.4815,12.4856,12.4877,12.4938,12.4985,12.4996,12.4997,12.5013,12.501,12.4997,12.4949,12.4921,12.4888,12.4843,12.4829,12.4877,12.494,12.4996,12.5007,12.5054,12.509,12.5138,12.5199,12.5254,12.5293,12.5321,12.5346,12.5318,12.5263,12.5207,12.5154,12.5118,12.5124,12.5157,12.5199,12.5238,12.5293,12.5354,12.5402,12.5429,12.5435,12.546,12.5502,12.5538,12.5557,12.5565,12.5582,12.5604,12.5627,12.5654,12.569,12.5752,12.5813,12.5863,12.5902,12.5943,12.599,12.6029,12.6093,12.614,12.619,12.6232,12.6237,12.6288,12.6329,12.6377,12.644,12.6502,12.6557,12.6596,12.6615,12.6615,12.6618,12.6646,12.6674,12.6715,12.6757,12.6804,12.686,12.6896,12.691,12.6913,12.6949,12.701,12.7071,12.7135,12.719,12.7204,12.7193,12.7188,12.7188,12.7204,12.7213,12.7235,12.7263,12.7318,12.7352,12.7407,12.7468,12.7524,12.7557,12.7582,12.761,12.7665,12.7688,12.7721,12.7771,12.7827,12.7893,12.7949,12.8013,12.8079,12.8143,12.8204,12.826,12.8315,12.8313,12.829,12.8282,12.831,12.8363,12.8418,12.8446,12.8488,12.8538,12.8565,12.8613,12.8668,12.8724,12.8771,12.8802,12.8843,12.8885,12.8932,12.9002,12.9049,12.9082,12.9118,12.9157,12.9221,12.9249,12.9257,12.9238,12.921,12.9199,12.9207,12.9229,12.9257,12.9285,12.9349,12.941,12.9471,12.9532,12.9596,12.9657,12.9718,12.9774,12.9829,12.9877,12.991,12.9929,12.9879,12.9818,12.9749,12.9693,12.9646,12.9602,12.9582,12.9607,12.9627,12.9613,12.959,12.9624,12.9671,12.9715,12.9768,12.9832,12.9871,12.9927,12.9977,12.9996,12.9997,13.0003,13.0004,13.0032,13.0054,13.0096,13.0143,13.0163,13.0177,13.0215,13.0271,13.0332,13.0352,13.0388,13.0402,13.0413,13.0479,13.0527,13.0582,13.0624,13.066,13.0702,13.0757,13.0799,13.0824,13.0829,13.0829,13.0849,13.086,13.0879,13.0913,13.0954,13.0988,13.1013,13.1054,13.1088,13.1113,13.1113,13.1113,13.1082,13.1096,13.1138,13.1199,13.1252,13.1302,13.134,13.1374,13.1421,13.1477,13.1527,13.156,13.161,13.1665,13.1699,13.174,13.181,13.1852,13.1877,13.1882,13.1888,13.1913,13.194,13.1982,13.2035,13.209,13.2152,13.2207,13.2257,13.2299,13.234,13.2382,13.2424,13.2468,13.2502,13.2557,13.2604,13.2646,13.2679,13.2727,13.2777,13.2824,13.2863,13.2882,13.2896,13.2929,13.2968,13.301,13.3043,13.3071,13.3104,13.3165,13.3227,13.3282,13.3321,13.3363,13.3418,13.3474,13.3521,13.359,13.3643,13.3657,13.3615,13.3579,13.3565,13.3577,13.359,13.3615,13.3663,13.3718,13.3782,13.3849,13.3899,13.3924,13.3929,13.3907,13.3871,13.3838,13.384,13.3846,13.3832,13.3824,13.3829,13.3849,13.3882,13.3929,13.3977,13.4024,13.4057,13.4099,13.4118,13.4124,13.4113,13.4113,13.4104,13.4068,13.4013,13.3952,13.3902,13.3852,13.381,13.3768,13.3727,13.3677,13.3624,13.356,13.3499,13.3438,13.3407,13.3393,13.3399,13.3377,13.334,13.3299,13.3249,13.3202,13.3157,13.3124,13.3093,13.3065,13.3057,13.3077,13.311,13.3152,13.3204,13.3252,13.3307,13.3377,13.3438,13.3493,13.3527,13.354,13.3543,13.3557,13.3577,13.3602,13.3643,13.3682,13.3724,13.3752,13.3749,13.3746,13.3788,13.3822,13.3857,13.3904,13.3954,13.4004,13.4052,13.4113,13.4177,13.4224,13.426,13.4296,13.4338,13.4393,13.4449,13.4504,13.4571,13.4627,13.4674,13.4715,13.474,13.4746,13.4743,13.4757,13.4782,13.4815,13.4871,13.4927,13.4979,13.4996,13.5021,13.5049,13.5065,13.5102,13.5135,13.516,13.5202,13.5235,13.5254,13.5252,13.5249,13.5268,13.5296,13.5335,13.5371,13.5377,13.5418,13.5457,13.5499,13.5532,13.556,13.5599,13.564,13.5679,13.5683,13.5721,13.5774,13.5829,13.5885,13.5927,13.5965,13.5993,13.6018,13.6046,13.6071,13.6104,13.616,13.6213,13.6268,13.6327,13.6379,13.6427,13.6468,13.6502,13.6535,13.6568,13.6602,13.6643,13.669,13.6738,13.6785,13.6849,13.6904,13.6952,13.6985,13.701,13.7029,13.7043,13.7046,13.7077,13.7077,13.7079,13.709,13.7113,13.7143,13.7165,13.7188,13.7215,13.7238,13.726,13.7282,13.731,13.7338,13.736,13.739,13.7413,13.744,13.7463,13.749,13.7521,13.7532,13.7549,13.7577,13.7613,13.764,13.7677,13.7713,13.774,13.7777,13.7813,13.7849,13.789,13.7932,13.7979,13.8029,13.8079,13.8132,13.8188,13.8243,13.8299,13.8363,13.8418,13.8474,13.8529,13.8585,13.864,13.8702,13.8771,13.8832,13.8902,13.8971,13.904,13.9088,13.9127,13.916,13.9188,13.9221,13.9254,13.9288,13.9313,13.9349,13.9374,13.9399,13.9432,13.946,13.9493,13.9527,13.9552,13.9568,13.9602,13.964,13.9696,13.9752,13.9821,13.9888,13.9957,13.9996,14.0018,14.0068,14.0121,14.0185,14.024,14.0282,14.0324,14.0374,14.044,14.0504,14.0565,14.0613,14.0663,14.0679,14.069,14.0727,14.0763,14.0799,14.0846,14.091,14.0971,14.1018,14.1065,14.1113,14.1177,14.1224,14.1274,14.1324,14.1365,14.1413,14.1477,14.1538,14.1607,14.1668,14.1724,14.1771,14.1818,14.1865,14.1921,14.1977,14.2029,14.2071,14.211,14.2122,14.216,14.2207,14.2254,14.2302,14.2357,14.2404,14.246,14.2513,14.256,14.2615,14.2663,14.2704,14.2743,14.2785,14.2818,14.281,14.2821,14.2849,14.2846,14.2843,14.2857,14.2871,14.2888,14.2915,14.294,14.2968,14.3002,14.304,14.3096,14.3138,14.3177,14.321,14.3243,14.3279,14.3279,14.3268,14.3296,14.336,14.3413,14.3457,14.3485,14.3535,14.3588,14.3629,14.3663,14.3682,14.3671,14.3657,14.3663,14.3704,14.3752,14.379,14.3802,14.3788,14.376,14.3738,14.3763,14.3796,14.3838,14.3871,14.3877,14.3868,14.3888,14.3913,14.3938,14.3965,14.3949,14.3996,14.4004,14.4015,14.4027,14.406,14.4079,14.4143,14.4171,14.421,14.426,14.4307,14.4363,14.4427,14.4479,14.4521,14.4555,14.4568,14.4624,14.4693,14.476,14.4802,14.4829,14.486,14.4902,14.4957,14.4996,14.4997,14.5018,14.5079,14.5124,14.5165,14.5213,14.5268,14.5332,14.5393,14.5449,14.5496,14.5513,14.5507,14.5507,14.5524,14.5552,14.5602,14.5657,14.5713,14.5774,14.5829,14.5877,14.5913,14.5949,14.601,14.6071,14.6115,14.6149,14.6171,14.6182,14.621,14.6246,14.6288,14.6335,14.639,14.646,14.6515,14.6554,14.6582,14.6588,14.6585,14.6596,14.6624,14.6671,14.6732,14.6763,14.6777,14.6813,14.6854,14.6904,14.6965,14.6957,14.6921,14.6921,14.696,14.7021,14.704,14.704,14.7024,14.7029,14.7057,14.7104,14.7152,14.7213,14.7268,14.7324,14.7365,14.7385,14.7402,14.7415,14.7435,14.7468,14.7502,14.7563,14.7582,14.7605,14.7624,14.764,14.7652,14.7677,14.7696,14.7724,14.7749,14.7768,14.7802,14.7829,14.7854,14.7888,14.7921,14.7949,14.7982,14.8015,14.8049,14.8082,14.8121,14.8154,14.8196,14.8238,14.8271,14.8296,14.8315,14.8343,14.836,14.8374,14.8385,14.8385,14.8396]}]]],null,null,{"interactive":true,"className":"","stroke":false,"color":"#03F","weight":5,"opacity":0.5,"fill":true,"fillColor":["#5AC864","#306A8E","#B5DE2B","#FDE725","#87D549","#48186A","#4DC36B","#CFE11C","#29AF7F","#440154","#55C667","#A7DB35"],"fillOpacity":1,"smoothFactor":0.3,"noClip":false},null,null,["Nord (&gt;2010): 94","Ouest: 55","Centre (&gt;2010): 114","Sud (&gt;2010): 132","Nord (&gt;2010): 104","Ouest: 38","Centre (&gt;2010): 91","Sud (&gt;2010): 120","Nord (&gt;2010): 81","Ouest: 35","Centre (&gt;2010): 93","Sud (&gt;2010): 111"],{"interactive":false,"permanent":false,"direction":"auto","opacity":1,"offset":[0,0],"textsize":"10px","textOnly":false,"className":"","sticky":true},null]},{"method":"addLegend","args":[{"colors":["#440154 , #482475 9.7018880335799%, #38588C 27.047764681159%, #26828E 44.3936413287381%, #25AC82 61.7395179763173%, #74D055 79.0853946238964%, #E7E419 96.4312712714756%, #FDE725 "],"labels":["40","50","63","79","100","126"],"na_color":null,"na_label":"NA","opacity":1,"position":"topright","type":"numeric","title":"Age specific fertility rate: 15-19","extra":{"p_1":0.097018880335799,"p_n":0.964312712714756},"layerId":null,"className":"info legend","group":null}]}],"limits":{"lat":[12.3013,16.6927],"lng":[-17.5455,-11.3575]}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->
---
title: "Anemia prevalence: an `rdhs` example"
author: "OJ Watson, Jeff Eaton"
date: "2018-09-24"
output: 
  rmarkdown::html_vignette:
    keep_md: TRUE
vignette: >
  %\VignetteIndexEntry{Anemia prevalence among women: an `rdhs` example}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Anemia is a common cause of fatigue, and women of childbearing age
are at particularly high risk for anemia. The package `rdhs` can be used
to compare estimates of the prevalence of any anemia among women from
Demographic and Health Surveys (DHS) conducted in Armenia, Cambodia,
and Lesotho.

# Setup

Load the `rdhs` package and other useful packages for analysing data.


```r
## devtools::install_github("ropensci/rdhs")
library(rdhs)
library(data.table)
library(ggplot2)
library(survey)
library(haven)
```

# Using calculated indicators from STATcompiler

Anemia prevalence among women is reported as a core indicator through the DHS STATcompiler (https://www.statcompiler.com/).
These indicators can be accessed directly from R via the DHS API with the function `dhs_data()`.

Query the API for a list of all StatCompiler indicators, and then search the indicators for those
that have `"anemia"` in the indicator name. API calls return `data.frame` objects, so if you prefer
to use `data.table` objects then convert afterwards, or we can set this up within our config using
`set_rdhs_config`.


```r
library(rdhs)
set_rdhs_config(data_frame = "data.table::as.data.table")

indicators <- dhs_indicators()
tail(indicators[grepl("anemia", Label), .(IndicatorId, ShortName, Label)])
```

```
##      IndicatorId                 ShortName                       Label
## 1: CN_ANMC_C_SEV Severe anemia (<7.0 g/dl) Children with severe anemia
## 2: AN_ANEM_W_ANY                      Any        Women with any anemia
## 3: AN_ANEM_W_MLD                     Mild       Women with mild anemia
## 4: AN_ANEM_W_MOD                 Moderate   Women with moderate anemia
## 5: AN_ANEM_W_SEV                   Severe     Women with severe anemia
## 6: AN_ANEM_M_ANY                Any anemia         Men with any anemia
```

The indicator ID `"AN_ANEM_W_ANY"` reports the percentage of women with any anemia.
The function `dhs_data()` will query the indicator dataset for the value of this indicator
for our three countries of interest. First, use `dhs_countries()` to query the
list of DHS countries to identify the DHS country code for each country.


```r
countries <- dhs_countries()
dhscc <- countries[CountryName %in% c("Armenia", "Cambodia", "Lesotho"), DHS_CountryCode]
dhscc
```

```
## [1] "AM" "KH" "LS"
```

Now query the indicators dataset for the women with any anemia indicator for these three countries.


```r
statcomp <- dhs_data(indicatorIds = "AN_ANEM_W_ANY", countryIds = dhscc)
statcomp[,.(Indicator, CountryName, SurveyYear, Value, DenominatorWeighted)]
```

```
##                 Indicator CountryName SurveyYear Value DenominatorWeighted
##  1: Women with any anemia     Armenia       2000  12.4                6137
##  2: Women with any anemia     Armenia       2005  24.6                6080
##  3: Women with any anemia     Armenia       2016  13.4                5769
##  4: Women with any anemia    Cambodia       2000  58.8                3634
##  5: Women with any anemia    Cambodia       2005  46.7                8219
##  6: Women with any anemia    Cambodia       2010  44.4                9229
##  7: Women with any anemia    Cambodia       2014  45.4               11286
##  8: Women with any anemia     Lesotho       2004  32.9                3008
##  9: Women with any anemia     Lesotho       2009  26.3                3839
## 10: Women with any anemia     Lesotho       2014  27.3                3297
```

```r
ggplot(statcomp, aes(SurveyYear, Value, col=CountryName)) +
  geom_point() + geom_line()
```

<img src="https://github.com/ropensci/rdhs/raw/development/vignettes/anemia_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

# Analyse DHS microdata

## Identify surveys that include anemia testing

The DHS API provides the facility to filter surveys according to particular characteristics.
We first query the list of survey characteristics and identify the `SurveyCharacteristicID`
that indicates the survey included anemia testing. The first command below queries the API
for the full list of survey characteristics, and the second uses `grepl()` to search
`SurveyCharacteristicName`s that include the word 'anemia'.


```r
surveychar <- dhs_survey_characteristics()
surveychar[grepl("anemia", SurveyCharacteristicName, ignore.case=TRUE)]
```

```
##    SurveyCharacteristicID SurveyCharacteristicName
## 1:                     15         Anemia questions
## 2:                     41           Anemia testing
```

The `SurveyCharacteristicID = 41` indicates that the survey included anemia testing. Next we
query the API to identify the surveys that have this characteristic and were conducted
in our countries of interest.


```r
surveys <- dhs_surveys(surveyCharacteristicIds = 41, countryIds = dhscc)
surveys[,.(SurveyId, CountryName, SurveyYear, NumberOfWomen, SurveyNum, FieldworkEnd)]
```

```
##      SurveyId CountryName SurveyYear NumberOfWomen SurveyNum FieldworkEnd
##  1: AM2000DHS     Armenia       2000          6430       203   2000-12-01
##  2: AM2005DHS     Armenia       2005          6566       262   2005-12-01
##  3: AM2016DHS     Armenia       2016          6116       492   2016-04-01
##  4: KH2000DHS    Cambodia       2000         15351       140   2000-07-01
##  5: KH2005DHS    Cambodia       2005         16823       257   2006-03-01
##  6: KH2010DHS    Cambodia       2010         18754       310   2011-01-01
##  7: KH2014DHS    Cambodia       2014         17578       464   2014-12-01
##  8: LS2004DHS     Lesotho       2004          7095       256   2005-01-01
##  9: LS2009DHS     Lesotho       2009          7624       317   2010-01-01
## 10: LS2014DHS     Lesotho       2014          6621       462   2014-12-01
```

Finally, query the API identify the individual recode (IR) survey datasets for each of these surveys


```r
datasets <- dhs_datasets(surveyIds = surveys$SurveyId, fileType = "IR", fileFormat="flat")
datasets[, .(SurveyId, SurveyNum, FileDateLastModified, FileName)]
```

```
##      SurveyId SurveyNum        FileDateLastModified     FileName
##  1: AM2000DHS       203   October, 05 2006 14:22:40 AMIR42FL.ZIP
##  2: AM2005DHS       262  February, 02 2010 10:38:12 AMIR54FL.zip
##  3: AM2016DHS       492 September, 21 2017 16:10:15 AMIR71FL.ZIP
##  4: KH2000DHS       140   October, 08 2007 12:31:53 KHIR42FL.zip
##  5: KH2005DHS       257   October, 18 2011 13:53:19 KHIR51FL.zip
##  6: KH2010DHS       310   October, 26 2011 11:11:07 KHIR61FL.ZIP
##  7: KH2014DHS       464      July, 28 2017 10:58:10 KHIR73FL.ZIP
##  8: LS2004DHS       256      July, 31 2007 13:14:31 LSIR41FL.ZIP
##  9: LS2009DHS       317  November, 10 2015 10:51:05 LSIR61FL.ZIP
## 10: LS2014DHS       462      June, 14 2016 11:35:19 LSIR71FL.ZIP
```

## Download datasets

To download datasets we need to first log in to our DHS account, by providing our credentials and setting up our configuration using `set_rdhs_config()`. This will require providing as arguments your `email` and `project` for which you want to download datasets from. You will then be prompted for your password. You can also specify a directory for datasets and API calls to be cached to using `cache_path`. In order to comply with CRAN, this function will also ask you for your permission to write to files outside your temporary directory, and you must type out the filename for the `config_path` - "rdhs.json". (See [introduction vignette](https://docs.ropensci.org/rdhs/articles/introduction.html) for specific format for config, or `?set_rdhs_config`).


```r
## set up your credentials
set_rdhs_config(email = "jeffrey.eaton@imperial.ac.uk",
                project = "Joint estimation of HIV epidemic trends and adult mortality")
```

After this the function `get_datasets()` returns a list of file paths where the desired datasets are saved in the cache. The first time a dataset is accessed, `rdhs` will download the dataset from the DHS program website using the supplied credentials. Subsequently, datasets will be simply be located in the cached repository. 


```r
datasets$path <- unlist(get_datasets(datasets$FileName))
```

```
## Logging into DHS website...
```

```
## Creating Download url list from DHS website...
```

## Identify survey variables

Anemia is defined as having a hemoglobin (Hb) <12.0 g/dL for non-pregnant women
or Hb <11.0 g/dL for currently pregnant women^[https://www.measureevaluation.org/prh/rh_indicators/womens-health/womens-nutrition/percent-of-women-of-reproductive-age-with-anemia].
To calculate anemia prevalence from DHS microdata, we must identify the DHS recode
survey variables for hemoglobin measurement and pregnancy status. This could be
done by consulting the DHS recode manual or the .MAP files accompanying survey
datasets. It is convenient though to do this in R by loading the first
individual recode dataset and searching the metadata for the variable
names corresponding to the hemoglobin measurement and pregnancy status.


```r
head(search_variable_labels(datasets$FileName[10], "hemoglobin")[,1:2])
```

```
##   variable
## 1     v042
## 2    v452c
## 3     v453
## 4     v455
## 5     v456
## 6   hw52_1
##                                                             description
## 1                                     Household selected for hemoglobin
## 2                                   Read consent statement - hemoglobin
## 3                                   Hemoglobin level (g/dl - 1 decimal)
## 4                                    Result of measurement - hemoglobin
## 5 Hemoglobin level adjusted for altitude and smoking (g/dl - 1 decimal)
## 6                                   Read consent statement - hemoglobin
```

Variable `v042` records the household selection for hemoglobin testing.
Variable `v455` reports the outcome of hemoglobin measurement and `v456`
the result of altitude adjusted hemoglobin levels.


```r
ir <- readRDS(datasets$path[10])
table(as_factor(ir$v042))
```

```
## 
## not selected     selected 
##         3203         3418
```

```r
table(as_factor(ir$v455))
```

```
## 
##                          measured                       not present 
##                              3349                                 2 
##                           refused                             other 
##                                35                                 8 
## no measurement found in household                           missing 
##                                 0                                24
```

```r
summary(ir$v456)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    24.0   118.0   130.0   145.6   141.0   999.0    3203
```

Variable `v454` reports the current pregnancy status used for determining the
anemia threshold.


```r
search_variable_labels(datasets$FileName[1], "currently.*pregnant")[,1:2]
```

```
##   variable        description
## 1     v213 Currently pregnant
## 2     v454 Currently pregnant
```

```r
table(as_factor(ir$v454))
```

```
## 
## no/don't know           yes       missing 
##          3276           142             0
```

We also keep a number of other variables related to the survey design and potentially
interesting covariates: country code and phase (`v000`), cluster number (`v001`),
sample weight (`v005`), age (`v012`), region (`v024`), urban/rural residence (`v025`),
and education level (`v106`).


```r
vars <- c("SurveyId", "CountryName", "SurveyYear", "v000", "v001", "v005",
          "v012", "v024", "v025", "v106", "v042", "v454", "v455", "v456")
```

## Extract survey data


```r
datlst <- list()

for(i in 1:nrow(datasets)){

  if(file.exists(datasets$path[i])){
  
  print(paste(i, datasets$SurveyId[i]))
  ir <- readRDS(datasets$path[i])

  ir$SurveyId <- datasets$SurveyId[i]
  ir$CountryName <- datasets$CountryName[i]
  ir$SurveyYear <- datasets$SurveyYear[i]

  datlst[[datasets$SurveyId[i]]] <- ir[vars]
  }
}
```

```
## [1] "1 AM2000DHS"
## [1] "2 AM2005DHS"
## [1] "3 AM2016DHS"
## [1] "4 KH2000DHS"
## [1] "5 KH2005DHS"
## [1] "6 KH2010DHS"
## [1] "7 KH2014DHS"
## [1] "8 LS2004DHS"
## [1] "9 LS2009DHS"
## [1] "10 LS2014DHS"
```

We use `rbind_labelled()` to combine datasets with labelled columns. The argument
`labels` describes to combine variable levels for all datasets for `v024` (region)
while providing a consistent set of value labels to be used for `v454` (currently
pregnant) for all datasets.



```r
dat <- rbind_labelled(datlst,
                      labels = list(v024 = "concatenate",
                                    v454 = c("no/don't know" = 0L,
                                             "yes" = 1L, "missing" = 9L)))
```

```
## Warning in rbind_labelled(datlst, labels = list(v024 = "concatenate", v454 = c(`no/don't know` = 0L, : Some variables have non-matching value labels: v106, v455, v456.
## Inheriting labels from first data frame with labels.
```

```r
sapply(dat, is.labelled)
```

```
##    SurveyId CountryName  SurveyYear        v000        v001        v005 
##       FALSE       FALSE       FALSE       FALSE       FALSE       FALSE 
##        v012        v024        v025        v106        v042        v454 
##       FALSE        TRUE        TRUE        TRUE        TRUE        TRUE 
##        v455        v456     DATASET 
##        TRUE        TRUE       FALSE
```

```r
dat$v456 <- zap_labels(dat$v456)
dat <- as_factor(dat)
```

## Data tabulations

It is a good idea to check basic tabulations of the data, especially by
survey to identify and nuances Exploratory analysis of variables


```r
with(dat, table(SurveyId, v025, useNA="ifany"))
```

```
##            v025
## SurveyId    urban rural
##   AM2000DHS  3545  2885
##   AM2005DHS  4592  1974
##   AM2016DHS  3545  2571
##   KH2000DHS  2627 12724
##   KH2005DHS  4152 12671
##   KH2010DHS  6077 12677
##   KH2014DHS  5667 11911
##   LS2004DHS  1945  5150
##   LS2009DHS  1977  5647
##   LS2014DHS  2202  4419
```

```r
with(dat, table(SurveyId, v106, useNA="ifany"))
```

```
##            v106
## SurveyId    no education primary secondary higher missing
##   AM2000DHS            5      24      5329   1072       0
##   AM2005DHS            7      24      5138   1397       0
##   AM2016DHS            5     406      2580   3125       0
##   KH2000DHS         4849    8182      2276     44       0
##   KH2005DHS         3772    9131      3771    149       0
##   KH2010DHS         3203    8796      6141    614       0
##   KH2014DHS         2233    7826      6535    984       0
##   LS2004DHS          169    4309      2520     97       0
##   LS2009DHS          114    3865      3277    368       0
##   LS2014DHS           81    2665      3354    521       0
```

```r
with(dat, table(SurveyId, v454, useNA="ifany"))
```

```
##            v454
## SurveyId    no/don't know   yes missing  <NA>
##   AM2000DHS          6231   199       0     0
##   AM2005DHS          5967   158     441     0
##   AM2016DHS          5939   177       0     0
##   KH2000DHS          3312   296      62 11681
##   KH2005DHS          7685   501     212  8425
##   KH2010DHS          8906   475       0  9373
##   KH2014DHS         10883   663       0  6032
##   LS2004DHS          2857   203       0  4035
##   LS2009DHS          3740   173     103  3608
##   LS2014DHS          3276   142       0  3203
```

```r
with(dat, table(SurveyId, v455, useNA="ifany"))
```

```
##            v455
## SurveyId    measured not present refused other no measurement found in hh
##   AM2000DHS     6137           5     264    24                          0
##   AM2005DHS     6134           8     294     1                          0
##   AM2016DHS     5807          11     295     0                          0
##   KH2000DHS     3666           0      68     0                          0
##   KH2005DHS     8182           2     185     5                          0
##   KH2010DHS     9225           9     106     0                          0
##   KH2014DHS    11390           8      13     2                          0
##   LS2004DHS     3061          15     377    56                          0
##   LS2009DHS     3896           1      78     5                          0
##   LS2014DHS     3349           2      35     8                          0
##            v455
## SurveyId    missing  <NA>
##   AM2000DHS       0     0
##   AM2005DHS     129     0
##   AM2016DHS       3     0
##   KH2000DHS       3 11614
##   KH2005DHS      24  8425
##   KH2010DHS      41  9373
##   KH2014DHS     133  6032
##   LS2004DHS      29  3557
##   LS2009DHS      36  3608
##   LS2014DHS      24  3203
```

```r
with(dat, table(v042, v454, useNA="ifany"))
```

```
##               v454
## v042           no/don't know   yes missing  <NA>
##   not selected             0     0       0 45778
##   selected             58796  2987     818   579
```

## Calculate anemia prevalence
Create indicator variable for 'any anemia'. The threshold depends on pregnancy status.


```r
dat$v456[dat$v456 == 999] <- NA
with(dat, table(v455, is.na(v456)))
```

```
##                             
## v455                         FALSE  TRUE
##   measured                   60847     0
##   not present                    0    61
##   refused                        0  1715
##   other                          0   101
##   no measurement found in hh     0     0
##   missing                        0   422
```

```r
dat$anemia <- as.integer(dat$v456  < ifelse(dat$v454 == "yes", 110, 120))
dat$anemia_denom <- as.integer(!is.na(dat$anemia))
```

Specify survey design using the `survey` package.


```r
dat$w <- dat$v005/1e6
des <- svydesign(~v001+SurveyId, data=dat, weights=~w)

anemia_prev <- svyby(~anemia, ~SurveyId, des, svyciprop, na.rm=TRUE, vartype="ci")
anemia_denom <- svyby(~anemia_denom, ~SurveyId, des, svytotal, na.rm=TRUE)

anemia_prev <- merge(anemia_prev, anemia_denom[c("SurveyId", "anemia_denom")])
res <- statcomp[,.(SurveyId, CountryName, SurveyYear, Value, DenominatorUnweighted, DenominatorWeighted)][anemia_prev, on="SurveyId"]

res$anemia <- 100*res$anemia
res$ci_l <- 100*res$ci_l
res$ci_u <- 100*res$ci_u
res$anemia_denom0 <- round(res$anemia_denom)
```

The table below compares the prevalence of any anemia calculated from survey microdata
with the estimates from DHS StatCompiler and the weighted denominators for each
calculation. The estimates are identical for most cases. There are some small
differences to be ironed out, which will require looking at the specific countries to check
how their stratification was carried out. (We are hoping to bring this feature in once the DHS
program has compiled how sample strata were constructed for all of their studies).


```r
knitr::kable(res[,.(CountryName, SurveyYear, Value, anemia, ci_l, ci_u, DenominatorWeighted, anemia_denom0)], digits=1)
```



CountryName    SurveyYear   Value   anemia   ci_l   ci_u   DenominatorWeighted   anemia_denom0
------------  -----------  ------  -------  -----  -----  --------------------  --------------
Armenia              2000    12.4     11.7   10.6   13.0                  6137            6137
Armenia              2005    24.6     23.1   21.3   24.9                  6080            6080
Armenia              2016    13.4     13.4   11.8   15.3                  5769            5769
Cambodia             2000    58.8     58.8   56.6   60.9                  3634            3634
Cambodia             2005    46.7     46.7   44.9   48.5                  8219            8219
Cambodia             2010    44.4     44.4   42.8   46.0                  9229            9229
Cambodia             2014    45.4     45.4   44.1   46.7                 11286           11286
Lesotho              2004    32.9     32.7   30.5   35.1                  3008            2789
Lesotho              2009    26.3     25.5   23.8   27.4                  3839            3839
Lesotho              2014    27.3     27.3   25.2   29.4                  3297            3297

```r
ggplot(res, aes(x=SurveyYear, y=anemia, ymin=ci_l, ymax=ci_u,
                col=CountryName, fill=CountryName)) +
  geom_ribbon(alpha=0.4, linetype="blank") + geom_point() + geom_line()
```

<img src="https://github.com/ropensci/rdhs/raw/development/vignettes/anemia_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />


# Regression analysis: relationship between education and anemia

A key use of the survey microdata are to conduct secondary analysis of pooled data
from several surveys, such as regression analysis. Here we investigate the
relationship between anemia prevalence and education level (`v106`) for women using
logistic regression, adjusting for urban/rural (`v025`) and fixed effects for each
survey.


```r
des <- update(des, v106 = relevel(v106, "primary"))
summary(svyglm(anemia ~ SurveyId + v025 + v106, des, family="binomial"))
```

```
## Warning in eval(family$initialize): non-integer #successes in a binomial
## glm!
```

```
## 
## Call:
## svyglm(formula = anemia ~ SurveyId + v025 + v106, des, family = "binomial")
## 
## Survey design:
## update(des, v106 = relevel(v106, "primary"))
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -1.91019    0.06478 -29.489  < 2e-16 ***
## SurveyIdAM2005DHS  0.82908    0.08086  10.253  < 2e-16 ***
## SurveyIdAM2016DHS  0.21583    0.09571   2.255 0.024488 *  
## SurveyIdKH2000DHS  2.14550    0.07503  28.596  < 2e-16 ***
## SurveyIdKH2005DHS  1.68112    0.07260  23.155  < 2e-16 ***
## SurveyIdKH2010DHS  1.61671    0.06961  23.224  < 2e-16 ***
## SurveyIdKH2014DHS  1.66621    0.06406  26.011  < 2e-16 ***
## SurveyIdLS2004DHS  1.13997    0.07960  14.322  < 2e-16 ***
## SurveyIdLS2009DHS  0.82962    0.07756  10.696  < 2e-16 ***
## SurveyIdLS2014DHS  0.93593    0.08164  11.464  < 2e-16 ***
## v025rural          0.11625    0.03220   3.610 0.000332 ***
## v106no education   0.15431    0.03845   4.013 6.75e-05 ***
## v106secondary     -0.11932    0.02787  -4.282 2.16e-05 ***
## v106higher        -0.33508    0.04985  -6.722 4.19e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 0.99451)
## 
## Number of Fisher Scoring iterations: 4
```

The results suggest that anemia prevalence is lower among women with higher education.
---
title: "Country Names and Codes"
author: "OJ Watson"
date: "2018-11-12"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Country Codes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## Overview

The following shows how to return a table of the country names and 2 letter codes. 

We can get this information by querying the API using `dhs_countries`, and
specifying the `returnFields` argument to just return the country name and 
DHS country code. This is a useful reference table when you want to pass in 
country IDs to a number of the API functions, e.g. `dhs_data(countryIds = "BJ")`



```r
library(rdhs)

## what are the countryIds
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
ids
```

```
##     DHS_CountryCode                      CountryName
##  1:              AF                      Afghanistan
##  2:              AL                          Albania
##  3:              AO                           Angola
##  4:              AM                          Armenia
##  5:              AZ                       Azerbaijan
##  6:              BD                       Bangladesh
##  7:              BJ                            Benin
##  8:              BO                          Bolivia
##  9:              BT                         Botswana
## 10:              BR                           Brazil
## 11:              BF                     Burkina Faso
## 12:              BU                          Burundi
## 13:              KH                         Cambodia
## 14:              CM                         Cameroon
## 15:              CV                       Cape Verde
## 16:              CF         Central African Republic
## 17:              TD                             Chad
## 18:              CO                         Colombia
## 19:              KM                          Comoros
## 20:              CG                            Congo
## 21:              CD        Congo Democratic Republic
## 22:              CI                    Cote d'Ivoire
## 23:              DR               Dominican Republic
## 24:              EC                          Ecuador
## 25:              EG                            Egypt
## 26:              ES                      El Salvador
## 27:              EK                Equatorial Guinea
## 28:              ER                          Eritrea
## 29:              ET                         Ethiopia
## 30:              GA                            Gabon
## 31:              GM                           Gambia
## 32:              GH                            Ghana
## 33:              GU                        Guatemala
## 34:              GN                           Guinea
## 35:              GY                           Guyana
## 36:              HT                            Haiti
## 37:              HN                         Honduras
## 38:              IA                            India
## 39:              ID                        Indonesia
## 40:              JO                           Jordan
## 41:              KK                       Kazakhstan
## 42:              KE                            Kenya
## 43:              KY                  Kyrgyz Republic
## 44:              LA Lao People's Democratic Republic
## 45:              LS                          Lesotho
## 46:              LB                          Liberia
## 47:              MD                       Madagascar
## 48:              MW                           Malawi
## 49:              MV                         Maldives
## 50:              ML                             Mali
## 51:              MR                       Mauritania
## 52:              MX                           Mexico
## 53:              MB                          Moldova
## 54:              MA                          Morocco
## 55:              MZ                       Mozambique
## 56:              MM                          Myanmar
## 57:              NM                          Namibia
## 58:              NP                            Nepal
## 59:              NC                        Nicaragua
## 60:              NI                            Niger
## 61:              NG                          Nigeria
## 62:              PK                         Pakistan
## 63:              PG                 Papua New Guinea
## 64:              PY                         Paraguay
## 65:              PE                             Peru
## 66:              PH                      Philippines
## 67:              RW                           Rwanda
## 68:              WS                            Samoa
## 69:              ST            Sao Tome and Principe
## 70:              SN                          Senegal
## 71:              SL                     Sierra Leone
## 72:              ZA                     South Africa
## 73:              LK                        Sri Lanka
## 74:              SD                            Sudan
## 75:              SZ                        Swaziland
## 76:              TJ                       Tajikistan
## 77:              TZ                         Tanzania
## 78:              TH                         Thailand
## 79:              TL                      Timor-Leste
## 80:              TG                             Togo
## 81:              TT              Trinidad and Tobago
## 82:              TN                          Tunisia
## 83:              TR                           Turkey
## 84:              TM                     Turkmenistan
## 85:              UG                           Uganda
## 86:              UA                          Ukraine
## 87:              UZ                       Uzbekistan
## 88:              VN                          Vietnam
## 89:              YE                            Yemen
## 90:              ZM                           Zambia
## 91:              ZW                         Zimbabwe
##     DHS_CountryCode                      CountryName
```
---
title: "Toolkit"
author: "OJ Watson"
date: "2018-07-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Toolkit}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## Overview

The following is no longer needed as the CRAN version of haven is working. This vignette is still here though in case we need a non-CRAN pacakge in the future. 

---

In order to install `rdhs` you will require a working development environment. The following guide details how to do this across different operating systems. Please note this is only temporary while we wait for a development version of haven to make it to CRAN. 
Apologies for the inconvenience

1. **Windows**: Install **[Rtools](https://cran.r-project.org/bin/windows/Rtools/)**. For help on how to install **Rtools** please see the following [guide](https://cran.r-project.org/bin/windows/Rtools/), paying particular attention to the section about adding Rtools to your system `PATH`. 

In order to find out which version of **Rtools** you will need to check which version of R you are running. This can be be found out using the `sessionInfo()` function:

``` r 
sessionInfo()
R version 3.4.4 (2016-06-21)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1
```

2. **Mac**: Install Xcode from the Mac App Store.

3. **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).

Once a working development environment is ready then the following should work for you:

```r
# install.packages("devtools")
# devtools::install_github("ropensci/rdhs")
# library(rdhs)
```
---
title: "How to use rdhs?"
author: "OJ Watson, Jeff Eaton"
date: "2019-01-22"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to use rdhs?}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

`rdhs` is a package for management and analysis of [Demographic and Health Survey (DHS)](https://www.dhsprogram.com/) data. This includes functionality to:

1. Access standard indicator data (i.e. [DHS STATcompiler](https://www.statcompiler.com/)) in R via the [DHS API](https://api.dhsprogram.com/).
1. Identify surveys and datasets relevant to a particular analysis.
1. Download survey datasets from the [DHS website](https://dhsprogram.com/data/available-datasets.cfm).
1. Load datasets and associated metadata into R.
1. Extract variables and combining datasets for pooled multi-survey analyses.

This process is described below and should cover most functionality that will be needed for working with these datasets. 

## 0. Installation

Install rdhs from github with `devtools`:


```r
# install.packages("devtools")
# devtools::install_github("ropensci/rdhs")
library(rdhs)
```

---

 > Before starting the tutorial, if you wish to download survey datasets from the DHS website, you will need to set up an account with the DHS website, which will enable you to request access to the datasets. Instructions on how to do this can be found [here](https://dhsprogram.com/data/Access-Instructions.cfm). The email, password, and project name that were used to create the account will then need to be provided to `rdhs` when attempting to download datasets. You can still interact with the DHS API in section 1-2 below without having an account with the DHS website, however, you will need to create an account if you wish to go through steps 3-5. 

---

## 1. Access standard indicator data via the API

The DHS programme has published an API that gives access to a number of different data sets, which each represent one of the DHS API endpoints (e.g.  https://api.dhsprogram.com/rest/dhs/tags, or https://api.dhsprogram.com/rest/dhs/surveys). These data sets include the standard health indicators that are available within [DHS STATcompiler](https://www.statcompiler.com/) as well as a series of meta data sets that describe the types of surveys that have been conducted as well as which raw dataset files are available from which surveys. Each of these data sets are described within the [DHS API website](https://api.dhsprogram.com/), and there are currently 12 different data sets available from the API. Each of these data sets can be accessed using anyone of `dhs_<>()` functions.  All exported functions within `rdhs` that start *dhs_* interact with a different data set of the [DHS API](https://api.dhsprogram.com/). Their website gives great information about the different search terms and filters that can be used, and we have tried to include all of this within the documentation of each function. Each of these data sets

One of those functions, `dhs_data()`, interacts with the the published set of standard health indicator data calculated by the DHS. This data set contains a set of health indicators that have been sample weighted to give country, subnational estimates that can be further refined by education and wealth brackets. To do this we use the `dhs_data()` function, which we can then either search for specific indicators, or by querying for indicators that have been tagged within specific areas.


```r
## what are the indicators
indicators <- dhs_indicators()
indicators[1,]
```

```
##                                                                                                            Definition
## 1: Age-specific fertility rate for the three years preceding the survey for age group 10-14 expressed per 1,000 women
##    NumberScale IndicatorType MeasurementType IsQuickStat  ShortName
## 1:           0             I            Rate           0 ASFR 10-14
##      IndicatorId    Level1 IndicatorTotalId          Level2 Level3
## 1: FE_FRTR_W_A10 Fertility                  Fertility rates  Women
##         SDRID IndicatorOldId TagIds DenominatorWeightedId
## 1: FEFRTRWA10                                            
##                                 Label IndicatorOrder
## 1: Age specific fertility rate: 10-14       11763005
##                                                                      Denominator
## 1: Per thousand women years exposed in the period 1-36 months prior to interview
##    QuickStatOrder IndicatorSpecial1Id DenominatorUnweightedId
## 1:                                                           
##    IndicatorSpecial2Id
## 1:
```

Each call to the DHS API returns a `data.frame` by default with all the results available by default. 

The DHS has a unique *IndicatorId* for each of the statistics it calculates. The definition and specific string for each indicator is included within the *IndicatorId* and *Definition* variable:


```r
# grab the first 5 alphabetically
indicators[order(indicators$IndicatorId),][1:5,c("IndicatorId", "Definition")]
```

```
##      IndicatorId
## 1: AH_CIGA_M_UNW
## 2: AH_CIGA_W_10P
## 3: AH_CIGA_W_12C
## 4: AH_CIGA_W_35C
## 5: AH_CIGA_W_69C
##                                                             Definition
## 1:                     Number of men who smoke cigarettes (unweighted)
## 2: Percentage of women who smoked 10+ cigarettes in preceding 24 hours
## 3: Percentage of women who smoked 1-2 cigarettes in preceding 24 hours
## 4: Percentage of women who smoked 3-5 cigarettes in preceding 24 hours
## 5: Percentage of women who smoked 6-9 cigarettes in preceding 24 hours
```


Since there are quite a lot of indicators, it might be easier to first query by tags. The DHS tags their indicators by what areas of demography and health they relate to, e.g. anaemia, literacy, malaria parasitaemia are all specific tags. First let's look at what the tags are, by interacting with the `dhs_tags()` function, before grabbing data that related to malaria parasitaemia in the DRC and Tanzania since 2010:


```r
# What are the tags
tags <- dhs_tags()

# Let's say we want to view the tags that relate to malaria
tags[grepl("Malaria", tags$TagName), ]
```

```
##    TagType                   TagName TagID TagOrder
## 1:       0       Malaria Parasitemia    36      540
## 2:       2 Select Malaria Indicators    79     1000
```

```r
# and now let's then grab this data by specifying the countryIds and the survey year starts
data <- dhs_data(tagIds = 36,countryIds = c("CD","TZ"),breakdown="subnational",surveyYearStart = 2010)
data[1,]
```

```
##     DataId                           Indicator  SurveyId IsPreferred Value
## 1: 1945295 Malaria prevalence according to RDT CD2013DHS           1  17.1
##         SDRID Precision        RegionId SurveyYearLabel SurveyType
## 1: MLPMALCRDT         1 CDDHS2013503010         2013-14        DHS
##    SurveyYear IndicatorOrder DHS_CountryCode CILow
## 1:       2013      125136010              CD      
##                  CountryName IndicatorType CharacteristicId
## 1: Congo Democratic Republic             I           503010
##    CharacteristicCategory   IndicatorId CharacteristicOrder
## 1:                 Region ML_PMAL_C_RDT             1503010
##    CharacteristicLabel ByVariableLabel DenominatorUnweighted
## 1:            Kinshasa                                   406
##    DenominatorWeighted CIHigh IsTotal ByVariableId
## 1:                 532              0            0
```


Depending on your analysis this maybe more than enough detail. It is also worth mentioning that this data can also be accessed via [DHS STATcompiler](https://www.statcompiler.com/) if you prefer a click and collect version. However, hopefully one can see that selecting a lot of different indicators for multiple countries and breakdowns should be a lot easier using the `rdhs` API interaction. For example we can very quickly find out the trends in antimalarial use in Africa, and see if perhaps antimalarial prescription has decreased after RDTs were introduced (assumed 2010). 


```r
# Make an api request
resp <- dhs_data(indicatorIds = "ML_FEVT_C_AML", surveyYearStart = 2010,breakdown = "subnational")

# filter it to 12 countries for space
countries  <- c("Angola","Ghana","Kenya","Liberia",
                "Madagascar","Mali","Malawi","Nigeria",
                "Rwanda","Sierra Leone","Senegal","Tanzania")

# and plot the results
library(ggplot2)
ggplot(resp[resp$CountryName %in% countries,],
       aes(x=SurveyYear,y=Value,colour=CountryName)) +
  geom_point() +
  geom_smooth(method = "glm") + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  ylab(resp$Indicator[1]) + 
  facet_wrap(~CountryName,ncol = 6) 
```

<img src="https://github.com/OJWatson/rdhs/raw/development/vignettes/introduction_files/figure-html/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />


If we incorrectly entered a filter query (very possible), `rdhs` will let us know our request was invalid:


```r
# Make an api request
resp <- dhs_data(indicatorIds="ML_FEVT_C_AMasfafasfL",
                 surveyYearStart=202231231306,
                 breakdown="subParTyping")
```

```
## Error in timeout_safe_request(url, timeout, encode = "json"): API Timeout Error: No response after 30 seconds.
##  Either increase timeout using set_rdhs_config(timeout = ...)
##  or check if the API is down by checking:
##  https://api.dhsprogram.com/rest/dhs/dataupdates
```


## 2. Identify surveys relevant for further analysis

You may, however, wish to do more nuanced analysis than the API allows. The following 4 sections detail a very basic example of how to quickly identify, download and extract datasets you are interested in.

Let's say we want to get all DHS survey data from the Democratic Republic of Congo and Tanzania in the last 5 years (since 2013), which covers the use of rapid diagnostic tests (RDTs) for malaria. To begin we'll interact with the DHS API to identify our datasets.

To start our extraction we'll query the *surveyCharacteristics* data set using `dhs_survey_characteristics()` function:


```r
## make a call with no arguments
sc <- dhs_survey_characteristics()
sc[grepl("Malaria", sc$SurveyCharacteristicName), ]
```

```
##    SurveyCharacteristicID SurveyCharacteristicName
## 1:                     96            Malaria - DBS
## 2:                     90     Malaria - Microscopy
## 3:                     89            Malaria - RDT
## 4:                     57          Malaria module 
## 5:                      8 Malaria/bednet questions
```


There are 87 different survey characteristics, with one specific survey characteristic for Malaria RDTs. We'll use this to then find the surveys that include this characteristic. We can also at this point filter for our desired countries and years. The DHS API allows for countries to be filtered using by their *countryIds*, which is one of the arguments in `dhs_surveys()`. To have a look at what each countries countryId is we can use another of the API functions:


```r
## what are the countryIds
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
str(ids)
```

```
## Classes 'data.table' and 'data.frame':	91 obs. of  2 variables:
##  $ DHS_CountryCode: chr  "AF" "AL" "AO" "AM" ...
##  $ CountryName    : chr  "Afghanistan" "Albania" "Angola" "Armenia" ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
# lets find all the surveys that fit our search criteria
survs <- dhs_surveys(surveyCharacteristicIds = 89,
                     countryIds = c("CD","TZ"),
                     surveyType = "DHS",
                     surveyYearStart = 2013)

# and lastly use this to find the datasets we will want to download and let's download the flat files (.dat) datasets (have a look in the dhs_datasets documentation for all argument options, and fileformat abbreviations etc.)
datasets <- dhs_datasets(surveyIds = survs$SurveyId, 
                         fileFormat = "flat")
str(datasets)
```

```
## Classes 'data.table' and 'data.frame':	19 obs. of  13 variables:
##  $ FileFormat          : chr  "Flat ASCII data (.dat)" "Flat ASCII data (.dat)" "Flat ASCII data (.dat)" "Flat ASCII data (.dat)" ...
##  $ FileSize            : int  198561 7030083 3226262 8028957 11426382 4794941 1569680 6595349 63022 996906 ...
##  $ DatasetType         : chr  "HIV Datasets" "Survey Datasets" "Survey Datasets" "Survey Datasets" ...
##  $ SurveyNum           : int  421 421 421 421 421 421 421 421 421 421 ...
##  $ SurveyId            : chr  "CD2013DHS" "CD2013DHS" "CD2013DHS" "CD2013DHS" ...
##  $ FileType            : chr  "HIV Test Results Recode" "Births Recode" "Couples' Recode" "Household Recode" ...
##  $ FileDateLastModified: chr  "November, 14 2014 12:48:34" "November, 17 2014 15:42:54" "November, 17 2014 15:43:04" "September, 19 2016 09:57:20" ...
##  $ SurveyYearLabel     : chr  "2013-14" "2013-14" "2013-14" "2013-14" ...
##  $ SurveyType          : chr  "DHS" "DHS" "DHS" "DHS" ...
##  $ SurveyYear          : int  2013 2013 2013 2013 2013 2013 2013 2013 2013 2013 ...
##  $ DHS_CountryCode     : chr  "CD" "CD" "CD" "CD" ...
##  $ FileName            : chr  "CDAR61FL.ZIP" "CDBR61FL.ZIP" "CDCR61FL.ZIP" "CDHR61FL.ZIP" ...
##  $ CountryName         : chr  "Congo Democratic Republic" "Congo Democratic Republic" "Congo Democratic Republic" "Congo Democratic Republic" ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

Lastly, we recommended to download either the spss (.sav), `fileFormat = "SV"`, or the flat file (.dat), `fileFormat = "FL"` datasets. The flat is quicker, but there are still one or two datasets that don't read correctly, whereas the .sav files are slower to read in but so far no datasets have been found that don't read in correctly.

We can now use this to download our datasets for further analysis. 

## 3. Download survey datasets

We can now go ahead and download our datasets. To be able to download survey datasets from the DHS website, you will need to set up an account with them to enable you to request access to the datasets. Instructions on how to do this can be found [here](https://dhsprogram.com/data/Access-Instructions.cfm). The email, password, and project name that were used to create the account will then need to be provided to `rdhs` when attempting to download datasets. 

Once we have created an account, we need to set up our credentials using the function `set_rdhs_config()`. This will require providing as arguments your `email` and `project` for which you want to download datasets from. You will then be prompted for your password.

You can also specify a directory for datasets and API calls to be cached to using `cache_path`. If you do not provide an argument for `cache_path` you will be prompted to provide permission to `rdhs` to save your datasets and API calls within your user cache directory for your operating system. This is to comply with CRAN's requests for permission to be granted before writing to system files. If you do not grant permission, these will be written within your R temporary directory (as we saw above when we first used one of the functions to query the API). Similarly if you do not also provide an argument for `config_path`, this will be saved within your temp directory unless permission is granted. Your config files will always be called "rdhs.json", so that `rdhs` can find them easily.


```r
## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations")
```

```
## Writing your configuration to:
##    -> /home/oj/.cache/rdhs/rdhs.json
```

```
## Adding /home/oj/.cache/rdhs/rdhs.json to .gitignore
```


Because you may have more than one project set up with the DHS website, you may want to have a separate directory for each set of datasets, and thus you will need to set up a different config file. To do this you need to set up a local config file. This can be achieved by setting the `global` param to `FALSE` (i.e. not global). You will also now need to provide the `config_path` argument, which MUST be **"rdhs.json"**. In order to comply with CRAN, you have to type this in (rather than have it as the default option).


```r
## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations",
                config_path = "rdhs.json",
                cache_path = "project_one",
                global = FALSE)
```

```
## Writing your configuration to:
##    -> rdhs.json
```

You may, however, not have different projects with the DHS website, in which case you may prefer to set up one global config file. If you do not want this to be saved in your user cache directory, you can set `global` to `TRUE` (the default) and this will save it in your R default launch directory. This MUST be **~/.rdhs.json**. (There is not really any difference between saving it at `~/.rdhs.json` vs the user cache directory, but you might want to have it somewhere easy to find etc).


```r
## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations",
                config_path = "~/.rdhs.json",
                global = TRUE)
```

```
## Writing your configuration to:
##    -> ~/.rdhs.json
```


After you have used `set_rdhs_config`, `rdhs` will try and find your config file when you next use `rdhs` in a different R session. It will do so by first looking locally for "rdhs.json", then globally for "~/.rdhs.json", then into your user cache directory, before lastly creating one in your temp directory. This is what was happening when you first used one of the API functions, and as such the config that is created to query the API initially will not be able to download datasets. 

Lastly, if you wish to return a `data.table` from your API requests, rather than a `data.frame` then you can change the default behaviour using the `data_frame` argument. You could also use this to convert them to `tibbles` and so on:


```r
## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations",
                config_path = "~/.rdhs.json",
                data_frame = "data.table::as.data.table",
                global = TRUE)
```

```
## Writing your configuration to:
##    -> ~/.rdhs.json
```


To see what config that is being used by `rdhs` at any point, then use `get_rdhs_config()` to view the config settings.

---

Before we download our datasets, it is worth mentioning that once you have set up your login credentials, your API calls will be cached for your within the cache directory used. This will allow those working remotely or without a good internet connection to be able to still return previous API requests. If you do not, API requests will still be cached within the temp directory, so will be quick to be returned a second time, but they will then be deleted when you start a new R session.


```r
# the first time we call this function, rdhs will make the API request
microbenchmark::microbenchmark(dhs_surveys(surveyYear = 1992),times = 1)
```

```
## Unit: milliseconds
##                            expr      min       lq     mean   median
##  dhs_surveys(surveyYear = 1992) 13.83736 13.83736 13.83736 13.83736
##        uq      max neval
##  13.83736 13.83736     1
```

```r
# with it cached it will be returned much quicker
microbenchmark::microbenchmark(dhs_surveys(surveyYear = 1992), times = 1)
```

```
## Unit: milliseconds
##                            expr      min       lq     mean   median
##  dhs_surveys(surveyYear = 1992) 4.307866 4.307866 4.307866 4.307866
##        uq      max neval
##  4.307866 4.307866     1
```


Now back to our dataset downloads. If we have a look back at our datasets object, we'll see there are 19 datasets listed. However, not all of them will be relevant to our malaria RDT questions. One approach is to head to the DHS website and have a look at the [DHS Recodes](https://dhsprogram.com/publications/publication-dhsg4-dhs-questionnaires-and-manuals.cfm), and look at the recodes that relate to the surveys. The other alternative is to download all the surveys and then query the variables within them. This is what we'll demonstrate here as it also demonstrates more of the package's functionality:

So first we will download all these datasets:


```r
# download datasets
downloads <- get_datasets(datasets$FileName)
```


The function returns a list with a file path to where the downloaded datasets have been saved to. By default the files will download quietly, i.e. no progress is shown. However, if you want to see the progress then you can control this by setting this in your config using the `verbose_download` argument.

## 4. Load datasets and associated metadata into R.

We can now examine what it is we have actually downloaded, by reading in one of these datasets:


```r
# read in our dataset
cdpr <- readRDS(downloads$CDPR61FL)
```

The dataset returned here contains all the survey questions within the dataset.  The dataset is by default stored as a *labelled* class from the [haven package](https://github.com/tidyverse/haven). This class preserves the original semantics and can easily be coerced to factors with `haven::as_factor()`. Special missing values are also preserved. For more info on the *labelled* class have a look at their github.

So if we have a look at what is returned for the variable *hv024*:


```r
head(cdpr$hv024)
```

```
## <Labelled integer>: Province
## [1] 4 4 4 4 4 4
## 
## Labels:
##  value            label
##      1         kinshasa
##      2         bandundu
##      3        bas-congo
##      4         equateur
##      5 kasai-occidental
##      6   kasai-oriental
##      7          katanga
##      8          maniema
##      9        nord-kivu
##     10        orientale
##     11         sud-kivu
```

```r
# and then the dataset
class(cdpr$hv024)
```

```
## [1] "haven_labelled"
```

If we want to get the data dictionary for this dataset, we can use the function `get_variable_labels`, which will return what question each of the variables in our dataset refer to:


```r
# let's look at the variable_names
head(get_variable_labels(cdpr))
```

```
##   variable                                                  description
## 1     hhid                                          Case Identification
## 2    hvidx                                                  Line number
## 3    hv000                                       Country code and phase
## 4    hv001                                               Cluster number
## 5    hv002                                             Household number
## 6    hv003 Respondent's line number (answering Household questionnaire)
```

For many of the survey responses this will give enough information for us to understand what the data is. However, for some questions it may be less clear exactly what the question means and how it may differ to other similar questions. If this is the case, then the DHS website publishes a lot of infrmation about the survey protocols and the surveys. We strongly advise for people to have a look through the [DHS website's documentation about using their datasets for analysis section](https://www.dhsprogram.com/data/Using-Datasets-for-Analysis.cfm), as well as the [recode files](https://www.dhsprogram.com/publications/publication-dhsg4-dhs-questionnaires-and-manuals.cfm) to understand how the surveys are carried out.

---

Above we saw that the default behaviour for the function `get_datasets` was to download the datasets, read them in, and save the resultant data.frame as a .rds object within the cache directory. You can control this behaviour using the `download_option` argument as such:

* `get_datasets(download_option = "zip")` - Just the downloaded zip will be saved
* `get_datasets(download_option = "rds")` - Just the read in rds will be saved
* `get_datasets(download_option = "both")` - The zip is downloaded and saved as well as the read in rds

The other main reason for reading the dataset in straight away as the default option is that `rdhs` will also create a table of all the survey variables and their labels (definitions) and cache them for you, which then allows us to quickly query for particular search terms or survey variables:


```r
# rapid diagnostic test search
questions <- search_variable_labels(datasets$FileName, search_terms = "malaria rapid test")

table(questions$dataset_filename)
```

```
## 
## CDHR61FL CDPR61FL TZHR7AFL TZPR7AFL 
##       24        1       48        1
```


What we see from the questions is that the question "Result of malaria rapid test" appears in a few different datasets. This is because the household member recode datasets (CDPR61SV, TZPR7ASV) stores information about the children in a household, with one row per child, whereas the household recode (CDHR61SV, TZHR7ASV) stores information about the household, and thus flattens the information from each child into different subvariables (hml35$01/02 etc). As such it is easier to extract this information from the household member recodes. 

## 5. Extract variables and combining datasets for pooled multi-survey analyses.

To extract our data we pass our questions object to the function `extract_dhs`, which will create a list with each dataset and its extracted data as a `data.frame`. We also have the option to add any geographic data available, which will download the geographic data files for you and add this data to you resultant extract:


```r
# let's just use the PR files thus
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "FL", fileType = "PR")
downloads <- get_datasets(datasets$FileName)

# and grab the questions from this again along with also questions detailing the province
questions <- search_variable_labels(datasets$FileName, search_terms = c("malaria rapid test"))

# and now extract the data
extract <- extract_dhs(questions, add_geo = FALSE)

# what does our extract look like
str(extract)
```

```
## List of 2
##  $ CDPR61FL:Classes 'dhs_dataset' and 'data.frame':	95949 obs. of  2 variables:
##   ..$ hml35   : 'haven_labelled' int [1:95949] NA NA NA NA NA NA NA 1 0 NA ...
##   .. ..- attr(*, "label")= chr "Result of malaria rapid test"
##   .. ..- attr(*, "labels")= Named int [1:3] 0 1 9
##   .. .. ..- attr(*, "names")= chr [1:3] "negative" "positive" "missing"
##   ..$ SurveyId: chr [1:95949] "CD2013DHS" "CD2013DHS" "CD2013DHS" "CD2013DHS" ...
##  $ TZPR7AFL:Classes 'dhs_dataset' and 'data.frame':	64880 obs. of  2 variables:
##   ..$ hml35   : 'haven_labelled' int [1:64880] NA NA NA NA NA NA NA 0 NA NA ...
##   .. ..- attr(*, "label")= chr "Result of malaria rapid test"
##   .. ..- attr(*, "labels")= Named int [1:3] 0 1 9
##   .. .. ..- attr(*, "names")= chr [1:3] "negative" "positive" "missing"
##   ..$ SurveyId: chr [1:64880] "TZ2015DHS" "TZ2015DHS" "TZ2015DHS" "TZ2015DHS" ...
```


The resultant extract is a list, with a new element for each different dataset that you have extracted. The responses from the dataset are by default stored as a *labelled* class from the [haven package](https://github.com/tidyverse/haven). 

We can also query our datasets for the survey question variables. In the example above the survey variable label was *Result of malaria rapid test* and the variable was *hml35*. So if you knew the survey variables that you wanted (either by looking at the Recode file or by looking through the *variable_names* included in the datasets) then we could search against these. So let's grab the regions using *hv024* using the client function `search_variables()`:


```r
# and grab the questions from this now utilising the survey variables
questions <- search_variables(datasets$FileName, variables = c("hv024","hml35"))

# and now extract the data
extract2 <- extract_dhs(questions, add_geo = FALSE)

# quick check 
head(extract2$CDPR61FL)
```

```
##   hv024 hml35  SurveyId
## 1     4    NA CD2013DHS
## 2     4    NA CD2013DHS
## 3     4    NA CD2013DHS
## 4     4    NA CD2013DHS
## 5     4    NA CD2013DHS
## 6     4    NA CD2013DHS
```

```r
head(extract2$TZPR7AFL)
```

```
##   hv024 hml35  SurveyId
## 1     1    NA TZ2015DHS
## 2     1    NA TZ2015DHS
## 3     1    NA TZ2015DHS
## 4     1    NA TZ2015DHS
## 5     1    NA TZ2015DHS
## 6     1    NA TZ2015DHS
```

```r
# and just to prove that hml35 did actually read in okay (there are just lots of NA)
table(extract2$CDPR61FL$hml35,useNA = "always")
```

```
## 
##     0     1     9  <NA> 
##  5260  2959     8 87722
```


We can now combine our two dataframes for further analysis using the `rdhs` package function `rbind_labelled()`. This function works specifically with our lists of labelled data.frames:


```r
# first let's bind our first extraction, without the hv024
extract_bound <- rbind_labelled(extract)

head(extract_bound)
```

```
##            hml35  SurveyId  DATASET
## CDPR61FL.1    NA CD2013DHS CDPR61FL
## CDPR61FL.2    NA CD2013DHS CDPR61FL
## CDPR61FL.3    NA CD2013DHS CDPR61FL
## CDPR61FL.4    NA CD2013DHS CDPR61FL
## CDPR61FL.5    NA CD2013DHS CDPR61FL
## CDPR61FL.6    NA CD2013DHS CDPR61FL
```

```r
# now let's try our second extraction
extract2_bound <- rbind_labelled(extract2)
```

```
## Warning in rbind_labelled(extract2): Some variables have non-matching value labels: hv024.
## Inheriting labels from first data frame with labels.
```


This hasn't quite done what we might want in the second instance. The *hv024* variable stores the regions for these 2 countries, which will not be the same and thus the labels will be different between the two of them. Without specifying any additional arguments `rbind_labelled()` will simply use the first data.frames labelling as the default, which will mean that some of the Tanzanian provinces will have been encoded as DRC provinces - not good! (This is a similar problem in nature to say trying to add new character strings to a factored data.frame).

There are a few work arounds. Firstly, we can specify a *labels* argument to the function which will detail how we should handle different variables. *labels* is a names list that specifies how to handle each variable. If we simply want to keep all the labels then we us the string "concatenate":


```r
# lets try concatenating the hv024
better_bound <- rbind_labelled(extract2, labels = list("hv024"="concatenate"))

head(better_bound$hv024)
```

```
## <Labelled integer>
## [1] 6 6 6 6 6 6
## 
## Labels:
##  value            label
##      1           arusha
##      2         bandundu
##      3        bas-congo
##      4    dar es salaam
##      5           dodoma
##      6         equateur
##      7            geita
##      8           iringa
##      9           kagera
##     10 kasai-occidental
##     11   kasai-oriental
##     12  kaskazini pemba
##     13 kaskazini unguja
##     14          katanga
##     15           katavi
##     16           kigoma
##     17      kilimanjaro
##     18         kinshasa
##     19     kusini pemba
##     20    kusini unguja
##     21            lindi
##     22          maniema
##     23          manyara
##     24             mara
##     25            mbeya
##     26  mjini magharibi
##     27         morogoro
##     28           mtwara
##     29           mwanza
##     30           njombe
##     31        nord-kivu
##     32        orientale
##     33            pwani
##     34            rukwa
##     35           ruvuma
##     36        shinyanga
##     37           simiyu
##     38          singida
##     39         sud-kivu
##     40           tabora
##     41            tanga
```


We could also specify new labels for a variable. For example, imagine the two datasets encoded their RDT responses differently, with the first one as `c("No","Yes")` and the other as `c("Negative","Positive")`. These would be for our purposes the same response, and so we could either leave it and all our results would use the `c("No","Yes")` labelling. But we may want to use the latter as it's more informative/correct, or we may want to be crystal clear and use `c("NegativeTest","PositiveTest")`. we can do that like this:


```r
# lets try concatenating the hv024 and providing new labels
better_bound <- rbind_labelled(
  extract2,
  labels = list("hv024"="concatenate",
                "hml35"=c("NegativeTest"=0, "PositiveTest"=1))
)

# and our new label
head(better_bound$hml35)
```

```
## <Labelled integer>: Result of malaria rapid test
## [1] NA NA NA NA NA NA
## 
## Labels:
##  value        label
##      0 NegativeTest
##      1 PositiveTest
```


The other option is to not use the labelled class at all. We can control this when we download our datasets, using the argument `reformat=TRUE`. This will ensure that no factors or labels are used and it is just the raw data. When this option is set the object returned by `get_datasets()` no longer has any labelled classes or factors. However, we can still recover the variable table for a dataset using `get_variable_labels()`, which will take any dataset output by `get_datasets()` and return a data.frame describing the survey question variables and definitions.  


```r
# download the datasets with the reformat arguments
downloads <- get_datasets(datasets$FileName, reformat=TRUE)

# grab the questions but specifying the reformat argument
questions <- search_variables(datasets$FileName, variables = c("hv024", "hml35"),
                                     reformat=TRUE)

# and now extract the data
extract3 <- extract_dhs(questions, add_geo = FALSE)

# group our results
bound_no_labels <- rbind_labelled(extract3)

# what does our hv024 look like now
class(bound_no_labels$hv024[1])
```

```
## [1] "character"
```


The *hv024* column is now just characters, which is possibly the best option depending on your downstream analysis/preferences. It's for this reason that the geographic data that is added is never turned into factors or labels.  

Lastly, we can now use our extract dataset to carry out some regression analysis, to investigate the relationship between malaria prevalence and the quality of wall materials. To do this we will need to first grab the sample weights and stratification from the surveys, along with the extra variables and we will then check the RDT prevalence calculated using the raw data versus the API:


```r
# grab the additional variable hv023 and hv024 which have the strata and weights respectively, and hc1 which is the age
questions <- search_variables(datasets$FileName,variables = c("hv005","hv021","hv022","hv023","hv024",
                                                              "hv025","hv214","hml20", "hc1","hml35"))
extraction <- extract_dhs(questions,TRUE)

# now concatenate the provinces as before and remove missing responses
dat <- rbind_labelled(extraction,labels=list("hv024"="concatenate","hv214"="concatenate"))
dat <- dat[-which(dat$hml35==9),] # remove missing responses

# and we are going to compare our extract to the API malaria prevalence by RDT, which is for those between 6 and 59 months
dat <- dat[-which(!dat$hc1>6 & dat$hc1<=60),]

# create a denominator response for hml35
dat$hml35denom <- as.integer(!is.na(dat$hml35))
dat$bricks <- dat$hv214 %in% c(8,18,5,9,10)
dat$net <- as.logical(dat$hml20)

# specify the strata and sample weights
dat$strata <- paste0(dat$hv023,dat$DATASET)
dat$hv005 <- dat$hv005/1e6

# construct a survey design using the survey pacakge
library(survey)

# construct the sample design and calculate the mean and totals 
des <-  survey::svydesign(~CLUSTER+DATASET,data=dat,weight=~hv005)
results <- cbind(survey::svyby(~hml35,by=~DHSREGNA+DATASET, des, survey::svyciprop,na.rm=TRUE),
                 survey::svyby(~hml35denom,by=~DHSREGNA+DATASET, des, survey::svytotal,na.rm=TRUE))
results <- results[order(results$DATASET),]

# grab the same data from the API 
dhs_api_data <- dhs_data(countryIds = c("CD","TZ"),indicatorIds = "ML_PMAL_C_RDT",breakdown = "subnational",surveyYearStart = 2013, surveyYearEnd = 2016)
dhs_api_data <- cbind(dhs_api_data$Value,dhs_api_data$DenominatorWeighted,dhs_api_data$CharacteristicLabel, dhs_api_data$SurveyId)
api <- dhs_api_data[!grepl("\\.\\.",dhs_api_data[,3]),] # remove subregions included in Tanzania
api <- api[order(apply(api[,4:3],1,paste,collapse="")),]

# bind the results and remove duplicate Region Columns
comparison <- cbind(results[,c(1,3,7)],api[])
names(comparison) <- c("Region","Survey_RDT_Prev","Survey_RDT_Denom","API_RDT_Prev","API_RDT_Denom","API_Regions","SurveyID")
head(comparison[,c(1,2,4,3,5,7)])
```

```
##                                     Region Survey_RDT_Prev API_RDT_Prev
## Bandundu.CDPR61FL                 Bandundu       0.2038607         20.2
## Bas-Congo.CDPR61FL               Bas-Congo       0.4761599         47.1
## Equateur.CDPR61FL                 Equateur       0.2732268         27.4
## Kasai-Occidental.CDPR61FL Kasai-Occidental       0.4507341         44.5
## Kasai-Oriental.CDPR61FL     Kasai-Oriental       0.4923113         49.4
## Katanga.CDPR61FL                   Katanga       0.3989890         38.9
##                           Survey_RDT_Denom API_RDT_Denom  SurveyID
## Bandundu.CDPR61FL                1415.6056          1414 CD2013DHS
## Bas-Congo.CDPR61FL                342.6473           347 CD2013DHS
## Equateur.CDPR61FL                1267.4746          1236 CD2013DHS
## Kasai-Occidental.CDPR61FL         604.6050           612 CD2013DHS
## Kasai-Oriental.CDPR61FL           892.6006           894 CD2013DHS
## Katanga.CDPR61FL                  861.2030           844 CD2013DHS
```


It's a little off, with the mean values differing due to maybe the specific cut off they used in terms of which ages were included within between 5 and 69. The variance could also be off due to the specific stratification the DHS Program will have used, as well as potentially how they have grouped the primary sampling units. We are hoping to get this information from the DHS for each survey so we can make this process more streamline for the end user. 

And lastly we will construct a logistic regression to investigate the relationship between a positive malaria RDT and whether the main walls of an individual's house were made of bricks or similar, while adjusting for urban/rural (`hv025`) and fixed effects for each survey.


```r
# contsruct our glm using svyglm and specify quasibinomial to handle the na in hml35
summary(svyglm(hml35 ~ DATASET + hv025 + net + bricks, des, family="quasibinomial"))
```

```
## 
## Call:
## svyglm(formula = hml35 ~ DATASET + hv025 + net + bricks, design = des, 
##     family = "quasibinomial")
## 
## Survey design:
## survey::svydesign(~CLUSTER + DATASET, data = dat, weight = ~hv005)
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     -1.65528    0.23267  -7.114 3.20e-12 ***
## DATASETTZPR7AFL -0.95553    0.12901  -7.406 4.39e-13 ***
## hv025            0.52558    0.13009   4.040 6.03e-05 ***
## netTRUE          0.06529    0.07163   0.911    0.362    
## bricksTRUE      -0.72109    0.13517  -5.335 1.36e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for quasibinomial family taken to be 0.9563804)
## 
## Number of Fisher Scoring iterations: 4
```


What we can see is that a significant negative gradient was associated with walls being made of bricks or similarly good materials in comparison to malaria positivity rates by RDT. What is also interesting is that whether the individual slept under a long lasting insecticidal net (*hml20* that we converted to net) was not significant.  

---

## Summary and further thoughts

Hopefully the above tutorial has shown how the `rdhs` package can facilitate both querying the DHS API and hopefully make downloading and interacting with the raw datasets a smoother, more reproducible process. It is worth bearing in mind though, that creating a harmonised dataset is not always as easy as the example above - a lot of the time survey variables differ across years and surveys, which is hopefully when the `survey_questions` functionality will make it easier to first filter down to those that include the relevant questions before having to decide which survey questions are valid. 

Any suggestions or comments/corrections/errors/ideas please let me know either in the issues or send me an email at "o.watson15@imperial.ac.uk". And if there is any further functionality that you think you would be useful, then also let me know. :)
---
title: "Downloading Shape Files for DHS Surveys"
author: "OJ Watson"
date: "2020-02-11"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Downloading shape files for DHS surveys}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

`rdhs` can be used to download shape files associated with the DHS surveys. For
example, we can link API responses to their shape files and plot the subnational
estimates using the spatial pacakge `sp`. To do this, we use the new function
`download_boundaries` from `rdhs` to download shape files. 


```r
# load our package
library(rdhs)

# make request
d <- dhs_data(countryIds = "SN",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2014,
              breakdown = "subnational")

# get our related spatial data frame object
sp <- download_boundaries(surveyId = d$SurveyId[1])
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: "/tmp/RtmpkjqL1S/file115c2d9e2d9d/shps", layer: "sdr_subnational_boundaries"
## with 4 features
## It has 27 fields
```

```r
# match our values to the regions
m <- d$Value[match(sp$sdr_subnational_boundaries$REG_ID, d$RegionId)]
sp$sdr_subnational_boundaries@data$Value <- m

# Use sp to plot
sp::spplot(sp$sdr_subnational_boundaries, "Value", main = d$Indicator[1])
```

![](https://github.com/OJWatson/rdhs/raw/development/vignettes/boundaries_files/figure-html/normal-1.png)<!-- -->

Or we can use `sf` for our plotting, which offers more user-friendly plotting
options. 


```r
# make request
d <- dhs_data(countryIds = "SN",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2017,
              breakdown = "subnational")

# get our related spatial data frame object
sp <- download_boundaries(surveyId = d$SurveyId[1], method = "sf")
```

```
## Reading layer `sdr_subnational_boundaries' from data source `/tmp/RtmpkjqL1S/file115c5afcf6fd/shps/sdr_subnational_boundaries.dbf' using driver `ESRI Shapefile'
## Simple feature collection with 14 features and 27 fields
## geometry type:  MULTIPOLYGON
## dimension:      XY
## bbox:           xmin: -17.54548 ymin: 12.30127 xmax: -11.35754 ymax: 16.69266
## epsg (SRID):    4326
## proj4string:    +proj=longlat +datum=WGS84 +no_defs
```

```r
# match our values to the regions
m <- d$Value[match(sp$sdr_subnational_boundaries$REG_ID, d$RegionId)]
sp$sdr_subnational_boundaries$Value <- m

# Use ggplot and geom_sf to plot
library(ggplot2)
ggplot(sp$sdr_subnational_boundaries) + 
  geom_sf(aes(fill = Value)) + 
  ggtitle(d$Indicator[1])
```

![](https://github.com/OJWatson/rdhs/raw/development/vignettes/boundaries_files/figure-html/sf-1.png)<!-- -->
---
title: "rdhs client object and internal package design"
author: "OJ Watson, Jeff Eaton"
date: "2018-09-24"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using the rdhs client}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

The following is a similar vignette - "How to use rdhs?". However, this one achieves the same results using the `rdhs` client object directly to download datasets etc, rather than using the easier user interfacing functions. This vignette is included to help demonstrate how the pacakge internals are working and adds more information about how the package was designed, which may be useful for anyone adding extra functionality to `rdhs` in the future. 

If are not bothered by this then please read the main introductory vignette. 

---

`rdhs` is a package for management and analysis of [Demographic and Health Survey (DHS)](https://www.dhsprogram.com/) data. This includes functionality to:

1. Access standard indicator data (i.e. [DHS STATcompiler](https://www.statcompiler.com/)) in R via the [DHS API](https://api.dhsprogram.com/).
1. Identify surveys and datasets relevant to a particular analysis.
1. Download survey datasets from the [DHS website](https://dhsprogram.com/data/available-datasets.cfm).
1. Load datasets and associated metadata into R.
1. Extract variables and combining datasets for pooled multi-survey analyses.

This process is described below and should cover most functionality that will be needed for working with these datasets. 

## 0. Installation

Install rdhs from github with `devtools`:


```r
# install.packages("devtools")
# devtools::install_github("ropensci/rdhs")

library(rdhs)
```


## 1. Access standard indicator data via the API

The DHS API has many *endpoints* that can be accessed using anyone of `dhs_<endpoint>()` functions. All exported functions within `rdhs` that start *dhs_* interact with a different with a different endpoint of the [DHS API](https://api.dhsprogram.com/). Their website gives great information about the different search terms and filters that can be used, and we have tried to include all of this within the documentation of each function.

One of those endpoints, `dhs_data()`, interacts with the the published set of standard health indicator data calculated by the DHS. This endpoint allows use to retrieve a set of health indicators that have been sample weighted to give country, subnational estimates that can be further refined by education and wealth brackets. To do this we use the `dhs_data()` endpoint, which we can then either search for specific indicators, or by querying for indicators that have been tagged within specific areas.


```r
## what are the indicators
indicators <- dhs_indicators()
indicators[1,]
```

```
##                                                                                                           Definition
## 1 Age-specific fertility rate for the three years preceding the survey for age group 10-14 expressed per 1,000 women
##   NumberScale IndicatorType MeasurementType IsQuickStat  ShortName
## 1           0             I            Rate           0 ASFR 10-14
##     IndicatorId    Level1 IndicatorTotalId          Level2 Level3
## 1 FE_FRTR_W_A10 Fertility                  Fertility rates  Women
##        SDRID IndicatorOldId TagIds DenominatorWeightedId
## 1 FEFRTRWA10                                            
##                                Label IndicatorOrder
## 1 Age specific fertility rate: 10-14       11763005
##                                                                     Denominator
## 1 Per thousand women years exposed in the period 1-36 months prior to interview
##   QuickStatOrder IndicatorSpecial1Id DenominatorUnweightedId
## 1                                                           
##   IndicatorSpecial2Id
## 1
```


Each call to a DHS endpoint returns a `data.frame` by default with all the results available by default. 

The DHS has a unique *IndicatorId* for each of the statistics it calculates. The definition and specific string for each indicator is included within the *IndicatorId* and *Definition* variable:


```r
# grab the first 5 alphabetically
indicators[order(indicators$IndicatorId),][1:5,c("IndicatorId", "Definition")]
```

```
##        IndicatorId
## 2031 AH_CIGA_M_UNW
## 2026 AH_CIGA_W_10P
## 2023 AH_CIGA_W_12C
## 2024 AH_CIGA_W_35C
## 2025 AH_CIGA_W_69C
##                                                               Definition
## 2031                     Number of men who smoke cigarettes (unweighted)
## 2026 Percentage of women who smoked 10+ cigarettes in preceding 24 hours
## 2023 Percentage of women who smoked 1-2 cigarettes in preceding 24 hours
## 2024 Percentage of women who smoked 3-5 cigarettes in preceding 24 hours
## 2025 Percentage of women who smoked 6-9 cigarettes in preceding 24 hours
```


Since there are quite a lot of indicators, it might be easier to first query by tags. The DHS tags their indicators by what areas of demography and health they relate to, e.g. anaemia, literacy, malaria parasitaemia are all specific tags. First let's look at what the tags are, by interacting with the `dhs_tags()` endpoint, before grabbing data that related to malaria parasitaemia in the DRC and Tanzania since 2010:


```r
# What are the tags
tags <- dhs_tags()

# Let's say we want to view the tags that relate to malaria
tags[grepl("Malaria", tags$TagName), ]
```

```
##    TagType                   TagName TagID TagOrder
## 31       0       Malaria Parasitemia    36      540
## 43       2 Select Malaria Indicators    79     1000
```

```r
# and now let's then grab this data by specifying the countryIds and the survey year starts
data <- dhs_data(tagIds = 36,countryIds = c("CD","TZ"),breakdown="subnational",surveyYearStart = 2010)
data[1,]
```

```
##   DataId                           Indicator  SurveyId IsPreferred Value
## 1 941297 Malaria prevalence according to RDT CD2013DHS           1  17.1
##        SDRID Precision        RegionId SurveyYearLabel SurveyType
## 1 MLPMALCRDT         1 CDDHS2013503010         2013-14        DHS
##   SurveyYear IndicatorOrder DHS_CountryCode CILow
## 1       2013      125136010              CD      
##                 CountryName IndicatorType CharacteristicId
## 1 Congo Democratic Republic             I           503010
##   CharacteristicCategory   IndicatorId CharacteristicOrder
## 1                 Region ML_PMAL_C_RDT             1503010
##   CharacteristicLabel ByVariableLabel DenominatorUnweighted
## 1            Kinshasa                                   406
##   DenominatorWeighted CIHigh IsTotal ByVariableId
## 1                 532              0            0
```


Depending on your analysis this maybe more than enough detail. It is also worth mentioning that this data can also be accessed via [DHS STATcompiler](https://www.statcompiler.com/) if you prefer a click and collect version. However, hopefully one can see that selecting a lot of different indicators for multiple countries and breakdowns should be a lot easier using the `rdhs` API interaction. For example we can very quickly find out the trends in antimalarial use in Africa, and see if perhaps antimalarial prescription has decreased after RDTs were introduced (assumed 2010). 


```r
# Make an api request
resp <- dhs_data(indicatorIds = "ML_FEVT_C_AML", surveyYearStart = 2010,breakdown = "subnational")

# filter it to 12 countries for space
countries  <- c("Angola","Ghana","Kenya","Liberia",
                "Madagascar","Mali","Malawi","Nigeria",
                "Rwanda","Sierra Leone","Senegal","Tanzania")

# and plot the results
library(ggplot2)
ggplot(resp[resp$CountryName %in% countries,],
       aes(x=SurveyYear,y=Value,colour=CountryName)) +
  geom_point() +
  geom_smooth(method = "glm") + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  ylab(resp$Indicator[1]) + 
  facet_wrap(~CountryName,ncol = 6) 
```

<img src="https://github.com/ropensci/rdhs/raw/development/vignettes/client_files/figure-html/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />


If we incorrectly entered a filter query (very possible), `rdhs` will let us know our request was invalid:


```r
# Make an api request
resp <- dhs_data(indicatorIds="ML_FEVT_C_AMasfafasfL",
                 surveyYearStart=202231231306,
                 breakdown="subParTyping")
```

```
## Error in handle_api_response(resp, TRUE): 
##    -> DHS API Request Failed [500] 
##    -> Error Type: dhs_internal_server_error
```


## 2. Identify surveys relevant for further analysis

You may, however, wish to do more nuanced analysis than the API allows. The following 4 section detail a very basic example of how to quickly identify, download and extract datasets you are interested in.

Let's say we want to get all the survey data from the Democratic Republic of Congo and Tanzania in the last 5 years (since 2013), which covers the use of rapid diagnostic tests (RDTs) for malaria. To begin we'll interact with the DHS API to identify our datasets.

To start our extraction we'll query the *surveyCharacteristics* endpoint using `dhs_surveyCharacteristics()`:


```r
## make a call with no arguments
sc <- dhs_survey_characteristics()
sc[grepl("Malaria", sc$SurveyCharacteristicName), ]
```

```
##    SurveyCharacteristicID SurveyCharacteristicName
## 57                     96            Malaria - DBS
## 58                     90     Malaria - Microscopy
## 59                     89            Malaria - RDT
## 60                     57          Malaria module 
## 61                      8 Malaria/bednet questions
```


There are 87 different survey characteristics, with one specific survey characteristic for Malaria RDTs. We'll use this to then find the surveys that include this characteristic. We can also at this point filter for our desired countries and years. The DHS API allows for countries to be filtered using by their *countryIds*, which is one of the arguments in `dhs_surveys()`. To have a look at what each countries countryId is we can use another of the API endpoints first:


```r
## what are the countryIds
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
str(ids)
```

```
## 'data.frame':	91 obs. of  2 variables:
##  $ DHS_CountryCode: chr  "AF" "AL" "AO" "AM" ...
##  $ CountryName    : chr  "Afghanistan" "Albania" "Angola" "Armenia" ...
```

```r
# lets find all the surveys that fit our search criteria
survs <- dhs_surveys(surveyCharacteristicIds = 89,countryIds = c("CD","TZ"),surveyYearStart = 2013)

# and lastly use this to find the datasets we will want to download and let's download the flat files (.dat) datasets (have a look in the dhs_datasets documentation for all argument options, and fileformat abbreviations etc.)
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "flat")
str(datasets)
```

```
## 'data.frame':	19 obs. of  13 variables:
##  $ FileFormat          : chr  "Flat ASCII data (.dat)" "Flat ASCII data (.dat)" "Flat ASCII data (.dat)" "Flat ASCII data (.dat)" ...
##  $ FileSize            : int  198561 7030083 3226262 8028957 11426382 4794941 1569680 6595349 63022 996906 ...
##  $ DatasetType         : chr  "HIV Datasets" "Survey Datasets" "Survey Datasets" "Survey Datasets" ...
##  $ SurveyNum           : int  421 421 421 421 421 421 421 421 421 421 ...
##  $ SurveyId            : chr  "CD2013DHS" "CD2013DHS" "CD2013DHS" "CD2013DHS" ...
##  $ FileType            : chr  "HIV Test Results Recode" "Births Recode" "Couples' Recode" "Household Recode" ...
##  $ FileDateLastModified: chr  "November, 14 2014 12:48:34" "November, 17 2014 15:42:54" "November, 17 2014 15:43:04" "September, 19 2016 09:57:20" ...
##  $ SurveyYearLabel     : chr  "2013-14" "2013-14" "2013-14" "2013-14" ...
##  $ SurveyType          : chr  "DHS" "DHS" "DHS" "DHS" ...
##  $ SurveyYear          : int  2013 2013 2013 2013 2013 2013 2013 2013 2013 2013 ...
##  $ DHS_CountryCode     : chr  "CD" "CD" "CD" "CD" ...
##  $ FileName            : chr  "CDAR61FL.ZIP" "CDBR61FL.ZIP" "CDCR61FL.ZIP" "CDHR61FL.ZIP" ...
##  $ CountryName         : chr  "Congo Democratic Republic" "Congo Democratic Republic" "Congo Democratic Republic" "Congo Democratic Republic" ...
```

Lastly, we recommended to download either the spss (.sav), `fileFormat = "SV"`, or the flat file (.dat), `fileFormat = "FL"` datasets. The flat is quicker, but there are still one or two datasets that don't read correctly, whereas the .sav files are slower to read in but so far no datasets have been found that don't read in correctly.

We can now use this to download our datasets for further analysis. 

## 3. Download survey datasets

We can now go ahead and download our datasets. To do this we need to first create a `client`. The client is an R6 class (similar to R's built in reference classes and make caching survey and API queries more reproducible) and will be used to log in to your DHS account, download datasets for you, and help query those datasets for the question you are interested in. The client will also cache all of these processes, which really helps increase the reproducibility of your analysis. 

In order to set up our credentials we use the function `set_rdhs_config()`. This will require providing as arguments your `email` and `project` for which you want to download datasets from. You will then be prompted for your password. 

You can also specify a directory for datasets and API calls to be cached to using `cache_path`. If you do not provide an argument for `cache_path` you will be prompted to provide permission to `rdhs` to save your datasets and API calls within your user cache directory for your operating system. This is to comply with CRAN's requests for permission to be granted before writing to system files. If you do not grant permission, these will be written within your R temporary directory (as we saw above when we first used one of the functions to query the API). Similarly if you do not also provide an argument for `config_path`, this will be saved within your temp directory unless permission is granted. Your config files will always be called "rdhs.json", so that `rdhs` can find them easily.


```r
## create a client
## set up your credentials
config <- set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations")

client <- client_dhs(config)
client
```

```
## <client_dhs>
##   Public:
##     available_datasets: function (clear_cache_first = FALSE) 
##     clear_namespace: function (namespace) 
##     dhs_api_request: function (api_endpoint, query = list(), api_key = private$api_key, 
##     extract: function (questions, add_geo = FALSE) 
##     get_cache_date: function () 
##     get_config: function () 
##     get_datasets: function (dataset_filenames, download_option = "rds", reformat = FALSE, 
##     get_downloaded_datasets: function () 
##     get_root: function () 
##     get_variable_labels: function (dataset_filenames = NULL, dataset_paths = NULL, rm_na = FALSE) 
##     initialize: function (config, api_key = NULL, root = rappdirs_rdhs()) 
##     save_client: function () 
##     set_cache_date: function (date) 
##     survey_questions: function (dataset_filenames, search_terms = NULL, essential_terms = NULL, 
##     survey_variables: function (dataset_filenames, variables, essential_variables = NULL, 
##   Private:
##     api_endpoints: data indicators countries surveys surveycharacteristics  ...
##     api_key: ********
##     cache_date: 2018-09-24 16:36:45
##     check_available_datasets: function (filenames) 
##     config: rdhs_config
##     na_s: ^na -|^na-|.*-na$|.* - na$| \{NA\}$|.* NA$|.*NA$
##     package_version: package_version, numeric_version
##     root: C:\Users\Oliver\AppData\Local\Oliver\rdhs\Cache
##     storr: storr, R6
##     url: https://api.dhsprogram.com/rest/dhs/
##     user_declared_root: NULL
```


Before we use our client to download our datasets, it is worth mentioning that the client can be passed as an argument to any of the API functions we have just seen. This will then cache the results for you, so that if you are working remotely or without a good internet connection you can still return your previous API requests. This is what is happening behind the scenes anyway, as `rdhs` will create a client for you when you call an API function (if one does not already exist).

---

Now back to our dataset downloads. If we have a look back at our datasets object, we'll see there are 19 datasets listed. However, not all of them will be relevant to our malaria RDT questions. One approach is to head to the DHS website and have a look at the [DHS Recodes](https://dhsprogram.com/publications/publication-dhsg4-dhs-questionnaires-and-manuals.cfm), and look at the recodes that relate to the surveys. The other alternative is to download all the surveys and then query the variables within them. This is what we'll demonstrate here as it also demonstrates more of the package's functionality:

So first we will download all these datasets:


```r
# download datasets
downloads <- client$get_datasets(datasets$FileName)
```

The function returns a list with a file path to where the downloaded datasets have been saved to. By default the files will download quietly, i.e. no progress is shown. However, if you want to see the progress then you can control this by setting this in your config using the `verbose_download` argument.

## 4. Load datasets and associated metadata into R.

We can now examine what it is we have actually downloaded, by reading in one of these datasets:


```r
# read in our dataset
cdpr <- readRDS(downloads$CDPR61FL)
```


The dataset returned here contains all the survey questions within the dataset, and what their survey variables are:


```r
# let's look at the variable_names
head(get_variable_labels(cdpr))
```

```
##   variable                                                  description
## 1     hhid                                          Case Identification
## 2    hvidx                                                  Line number
## 3    hv000                                       Country code and phase
## 4    hv001                                               Cluster number
## 5    hv002                                             Household number
## 6    hv003 Respondent's line number (answering Household questionnaire)
```

```r
# and then the dataset
class(cdpr$hv024)
```

```
## [1] "labelled"
```


This is the default behaviour for the client function `get_datasets` - it will download the datasets for you, and then by default save them in your client's root directory and then unzip them and read them in for you, and save the resultant data.frame as a .rds object within the client's root directory. You can control this behaviour using the `download_option` argument as such:

* `client$get_datasets(download_option = "zip")` - Just the downloaded zip will be saved
* `client$get_datasets(download_option = "rds")` - Just the read in rds will be saved
* `client$get_datasets(download_option = "both")` - The zip is downloaded and saved as well as the read in rds

The other main reason for reading the dataset in straight away as the default option is that the created table of all the survey variables and their definitions is cached then and there, which then allows us to quickly query for particular search terms or survey variables:


```r
# rapid diagnostic test search
questions <- client$survey_questions(datasets$FileName, search_terms = "malaria rapid test")

table(questions$dataset_filename)
```

```
## 
## CDHR61FL CDPR61FL TZHR7AFL TZPR7AFL 
##       24        1       48        1
```


What we see from the questions is that the question "Result of malaria rapid test" appears in a few different datasets. This is because the household member recode datasets (CDPR61SV, TZPR7HSV) stores information about the children in a household, with one row per child, whereas the household recode (CDHR61SV, TZHR7HSV) stores information about the household, and thus flattens the information from each child into different subvariables (hml35$01/02 etc). As such it is easier to extract this information from the household member recodes. 

## 5. Extract variables and combining datasets for pooled multi-survey analyses.

To extract our data we pass our questions object to the client function `extract`, which will create a list with each dataset and its extracted data as a `data.frame`. We also have the option to add any geographic data available, which will download the geographic data files for you and add this data to you resultant extract:


```r
# let's just use the PR files thus
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "FL", fileType = "PR")
downloads <- client$get_datasets(datasets$FileName)

# and grab the questions from this again along with also questions detailing the province
questions <- client$survey_questions(datasets$FileName, search_terms = c("malaria rapid test"))

# and now extract the data
extract <- client$extract(questions, add_geo = FALSE)

# what does our extract look like
str(extract)
```

```
## List of 2
##  $ CDPR61FL:Classes 'dhs_dataset' and 'data.frame':	95949 obs. of  2 variables:
##   ..$ hml35   : 'labelled' int [1:95949] NA NA NA NA NA NA NA 1 0 NA ...
##   .. ..- attr(*, "label")= chr "Result of malaria rapid test"
##   .. ..- attr(*, "labels")= Named int [1:3] 0 1 9
##   .. .. ..- attr(*, "names")= chr [1:3] "negative" "positive" "missing"
##   ..$ SurveyId: chr [1:95949] "CD2013DHS" "CD2013DHS" "CD2013DHS" "CD2013DHS" ...
##  $ TZPR7AFL:Classes 'dhs_dataset' and 'data.frame':	64880 obs. of  2 variables:
##   ..$ hml35   : 'labelled' int [1:64880] NA NA NA NA NA NA NA 0 NA NA ...
##   .. ..- attr(*, "label")= chr "Result of malaria rapid test"
##   .. ..- attr(*, "labels")= Named int [1:3] 0 1 9
##   .. .. ..- attr(*, "names")= chr [1:3] "negative" "positive" "missing"
##   ..$ SurveyId: chr [1:64880] "TZ2015DHS" "TZ2015DHS" "TZ2015DHS" "TZ2015DHS" ...
```


The resultant extract is a list, with a new element for each different dataset that you have extracted. The responses from the dataset are by default stored as a *labelled* class from the [haven package](https://github.com/tidyverse/haven). This class preserves the original semantics and can easily be coerced to factors with `haven::as_factor()`. Special missing values are also preserved. For more info on the *labelled* class have a look at their github.

We can also query our datasets for the survey question variables. In the example above the survey question was *Result of malaria rapid test* and the variable was *hml35*. So if you knew the survey variables that you wanted (either by looking at the Recode file or by looking through the *variable_names* included in the datasets) then we could search against these. So let's grab the regions using *hv024* using the client function `survey_variables()`:


```r
# and grab the questions from this now utilising the survey variables
questions <- client$survey_variables(datasets$FileName, variables = c("hv024","hml35"))

# and now extract the data
extract2 <- client$extract(questions, add_geo = FALSE)

# quick check 
head(extract2$CDPR61FL)
```

```
##   hv024 hml35  SurveyId
## 1     4    NA CD2013DHS
## 2     4    NA CD2013DHS
## 3     4    NA CD2013DHS
## 4     4    NA CD2013DHS
## 5     4    NA CD2013DHS
## 6     4    NA CD2013DHS
```

```r
head(extract2$TZPR7HFL)
```

```
## NULL
```

```r
# and just to prove that hml35 did actually read in okay (there are just lots of NA)
table(extract2$CDPR61FL$hml35,useNA = "always")
```

```
## 
##     0     1     9  <NA> 
##  5260  2959     8 87722
```


We can now combine our two dataframes for further analysis using the `rdhs` package function `rbind_labelled()`. This function works specifically with our lists of labelled data.frames:


```r
# first let's bind our first extraction, without the hv024
extract_bound <- rbind_labelled(extract)

head(extract_bound)
```

```
##            hml35  SurveyId  DATASET
## CDPR61FL.1    NA CD2013DHS CDPR61FL
## CDPR61FL.2    NA CD2013DHS CDPR61FL
## CDPR61FL.3    NA CD2013DHS CDPR61FL
## CDPR61FL.4    NA CD2013DHS CDPR61FL
## CDPR61FL.5    NA CD2013DHS CDPR61FL
## CDPR61FL.6    NA CD2013DHS CDPR61FL
```

```r
# now let's try our second extraction
extract2_bound <- rbind_labelled(extract2)
```

```
## Warning in rbind_labelled(extract2): Some variables have non-matching value labels: hv024.
## Inheriting labels from first data frame with labels.
```


This hasn't quite done what we might want in the second instance. The *hv024* variable stores the regions for these 2 countries, which will not be the same and thus the labels will be different between the two of them. Without specifying any additional arguments `rbind_labelled()` will simply use the first data.frames labelling as the default, which will mean that some of the Tanzanian provinces will have been encoded as DRC provinces - not good! (This is a similar problem in nature to say trying to add new character strings to a factored data.frame).

There are a few work arounds. Firstly, we can specify a *labels* argument to the function which will detail how we should handle different variables. *labels* is a names list that specifies how to handle each variable. If we simply want to keep all the labels then we us the string "concatenate":


```r
# lets try concatenating the hv024
better_bound <- rbind_labelled(extract2, labels = list("hv024"="concatenate"))

head(better_bound$hv024)
```

```
## <Labelled integer>
## [1] 6 6 6 6 6 6
## 
## Labels:
##  value            label
##      1           arusha
##      2         bandundu
##      3        bas-congo
##      4    dar es salaam
##      5           dodoma
##      6         equateur
##      7            geita
##      8           iringa
##      9           kagera
##     10 kasai-occidental
##     11   kasai-oriental
##     12  kaskazini pemba
##     13 kaskazini unguja
##     14          katanga
##     15           katavi
##     16           kigoma
##     17      kilimanjaro
##     18         kinshasa
##     19     kusini pemba
##     20    kusini unguja
##     21            lindi
##     22          maniema
##     23          manyara
##     24             mara
##     25            mbeya
##     26  mjini magharibi
##     27         morogoro
##     28           mtwara
##     29           mwanza
##     30           njombe
##     31        nord-kivu
##     32        orientale
##     33            pwani
##     34            rukwa
##     35           ruvuma
##     36        shinyanga
##     37           simiyu
##     38          singida
##     39         sud-kivu
##     40           tabora
##     41            tanga
```


We could also specify new labels for a variable. For example, imagine the two datasets encoded their RDT responses differently, with the first one as `c("No","Yes")` and the other as `c("Negative","Positive")`. These would be for our purposes the same response, and so we could either leave it and all our results would use the `c("No","Yes")` labelling. But we may want to use the latter as it's more informative/correct, or we may want to be crystal clear and use `c("NegativeTest","PositiveTest")`. we can do that like this:


```r
# lets try concatenating the hv024 and providing new labels
better_bound <- rbind_labelled(extract2,
                               labels = list("hv024"="concatenate",
                                             "hml35"=c("NegativeTest"=0, "PositiveTest"=1)))

# and our new label
head(better_bound$hml35)
```

```
## <Labelled integer>
## [1] NA NA NA NA NA NA
## 
## Labels:
##  value        label
##      0 NegativeTest
##      1 PositiveTest
```


The other option is to not use the labelled class at all. We can control this when we download our datasets, using the argument `reformat=TRUE`. This will ensure that no factors or labels are used and it is just the raw data. When this option is set the object returned by `client$get_datasets()` no longer has any labelled classes or factors. However, we can still recover the variable table for a dataset using `get_variable_labels()`, which will take any dataset output by `get_datasets()` and return a data.frame describing the survey question variables and definitions.  


```r
# download the datasets with the reformat arguments
downloads <- client$get_datasets(datasets$FileName, reformat=TRUE)

# grab the questions but specifying the reformat argument
questions <- client$survey_variables(datasets$FileName, variables = c("hv024", "hml35"),
                                     reformat=TRUE)

# and now extract the data
extract3 <- client$extract(questions, add_geo = FALSE)

# group our results
bound_no_labels <- rbind_labelled(extract3)

# what does our hv024 look like now
class(bound_no_labels$hv024[1])
```

```
## [1] "character"
```


The *hv024* column is now just characters, which is possibly the best option depending on your downstream analysis/preferences. It's for this reason that the geographic data that is added is never turned into factors or labels.  

Lastly, we can now use our extract dataset to carry out some regression analysis, to investigate the relationship between malaria prevalence and the quality of wall materials. To do this we will need to first grab the sample weights and stratification from the surveys, along with the extra variables and we will then check the RDT prevalence calculated using the raw data versus the API:


```r
# grab the additional variable hv023 and hv024 which have the strata and weights respectively
questions <- client$survey_variables(datasets$FileName,variables = c("hv005","hv021","hv022","hv023","hv024",
                                                              "hv025","hv214","hml20", "hc1","hml35"))
extraction <- client$extract(questions,TRUE)

# now concatenate the provinces as before and remove missing responses
dat <- rbind_labelled(extraction,labels=list("hv024"="concatenate","hv214"="concatenate"))
dat <- dat[-which(dat$hml35==9),] # remove missing responses

# and we are going to compare our extract to the API malaria prevalence by RDT, which is for those between 6 and 59 months
dat <- dat[-which(!dat$hc1>6 & dat$hc1<=60),]

# create a denominator response for hml35
dat$hml35denom <- as.integer(!is.na(dat$hml35))
dat$bricks <- dat$hv214 %in% c(8,18,5,9,10)
dat$net <- as.logical(dat$hml20)

# specify the strata and sample weights
dat$strata <- paste0(dat$hv023,dat$DATASET)
dat$hv005 <- dat$hv005/1e6

# construct a survey design using the survey pacakge
library(survey)

# construct the sample design and calculate the mean and totals 
des <-  survey::svydesign(~CLUSTER+DATASET,data=dat,weight=~hv005)
results <- cbind(survey::svyby(~hml35,by=~DHSREGNA+DATASET, des, survey::svyciprop,na.rm=TRUE),
                 survey::svyby(~hml35denom,by=~DHSREGNA+DATASET, des, survey::svytotal,na.rm=TRUE))
results <- results[order(results$DATASET),]

# grab the same data from the API 
dhs_api_data <- dhs_data(countryIds = c("CD","TZ"),indicatorIds = "ML_PMAL_C_RDT",breakdown = "subnational",surveyYearStart = 2013, surveyYearEnd = 2016)
dhs_api_data <- cbind(dhs_api_data$Value,dhs_api_data$DenominatorWeighted,dhs_api_data$CharacteristicLabel, dhs_api_data$SurveyId)
api <- dhs_api_data[!grepl("\\.\\.",dhs_api_data[,3]),] # remove subregions included in Tanzania
api <- api[order(apply(api[,4:3],1,paste,collapse="")),]

# bind the results and remove duplicate Region Columns
comparison <- cbind(results[,c(1,3,7)],api[])
names(comparison) <- c("Region","Survey_RDT_Prev","Survey_RDT_Denom","API_RDT_Prev","API_RDT_Denom","API_Regions","SurveyID")

head(comparison[,c(1,2,4,3,5,7)])
```

```
##                                     Region Survey_RDT_Prev API_RDT_Prev
## Bandundu.CDPR61FL                 Bandundu       0.2038607         20.2
## Bas-Congo.CDPR61FL               Bas-Congo       0.4761599         47.1
## Equateur.CDPR61FL                 Equateur       0.2732268         27.4
## Kasai-Occidental.CDPR61FL Kasai-Occidental       0.4507341         44.5
## Kasai-Oriental.CDPR61FL     Kasai-Oriental       0.4923113         49.4
## Katanga.CDPR61FL                   Katanga       0.3989890         38.9
##                           Survey_RDT_Denom API_RDT_Denom  SurveyID
## Bandundu.CDPR61FL                1415.6056          1414 CD2013DHS
## Bas-Congo.CDPR61FL                342.6473           347 CD2013DHS
## Equateur.CDPR61FL                1267.4746          1236 CD2013DHS
## Kasai-Occidental.CDPR61FL         604.6050           612 CD2013DHS
## Kasai-Oriental.CDPR61FL           892.6006           894 CD2013DHS
## Katanga.CDPR61FL                  861.2030           844 CD2013DHS
```


It's a little off, which will be due to the specific stratification the DHS Program will have used, as well as potentially how they have grouped the primary sampling units. We are hoping to get this information from the DHS for each survey so we can make this process more streamline again.

And lastly we will construct a logistic regression to investigate the relationship between a positive malaria RDT and whether the main walls of an individual's house were made of bricks or similar, while adjusting for urban/rural (`v025`) and fixed effects for each survey.


```r
# contsruct our glm using svyglm and specify quasibinomial to handle the na in hml35
summary(svyglm(hml35 ~ DATASET + hv025 + net + bricks, des, family="quasibinomial"))
```

```
## 
## Call:
## svyglm(formula = hml35 ~ DATASET + hv025 + net + bricks, des, 
##     family = "quasibinomial")
## 
## Survey design:
## survey::svydesign(~CLUSTER + DATASET, data = dat, weight = ~hv005)
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     -1.65528    0.23267  -7.114 3.20e-12 ***
## DATASETTZPR7AFL -0.95553    0.12901  -7.406 4.39e-13 ***
## hv025            0.52558    0.13009   4.040 6.03e-05 ***
## netTRUE          0.06529    0.07163   0.911    0.362    
## bricksTRUE      -0.72109    0.13517  -5.335 1.36e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for quasibinomial family taken to be 0.9563804)
## 
## Number of Fisher Scoring iterations: 4
```


What we can see is that a significant negative gradient was associated with walls being made of bricks or similarly good materials in comparison to malaria positivity rates by RDT. What is also interesting is that whether the individual slept under a long lasting insecticidal net (*hml20* that we converted to net) was not significant.  

---

## Summary and further thoughts

Hopefully the above tutorial has shown how the `rdhs` package can facilitate both querying the DHS API and hopefully make downloading and interacting with the raw datasets a smoother, more reproducible process. It is worth bearing in mind though, that creating a harmonised dataset is not always as easy as the example above - a lot of the time survey variables differ across years and surveys, which is hopefully when the `survey_questions` functionality will make it easier to first filter down to those that include the relevant questions before having to decide which survey questions are valid. 

Any suggestions or comments/corrections/errors/ideas please let me know either in the issues or send me an email at "o.watson15@imperial.ac.uk". And if there is any further functionality that you think you would be useful, then also let me know. :)
---
title: "Running test suite locally"
author: "OJ Watson"
date: "2018-07-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Testing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## Overview

The following is a guide describing how to be able to run the full test suite locally.

---

## 1. Creating an account with the DHS website

In order to run all the tests in the *rdhs* package you need to have created an account with the DHS website. This requires about 5 mins to fill in their form and then a few days for them to approve your account creation. To start this head to the DHS website and go to the Data tab [here](https://dhsprogram.com/data/dataset_admin/login_main.cfm). A full rationale for this process can be found on their website [here](https://dhsprogram.com/data/Access-Instructions.cfm). In short you provide your contact information and then describe the study you are conducting that requires the DHS survey datasets.

Importantly, when filling in their form describing the description of your study be sure to include some mention of needing geographic datasets as a number of the tests download geographic datasets. After you have filled this in you need to specify which countries and datasets you need. In order for the full test suite to work you must request every country and dataset type possible. After you have done this the DHS programme will take a few days to approve you request and will then email you to let you know at the email address you have provided. 

## 2. Setting up your rdhs config

Having now created your account, you can set up your *rdhs* configuration, which will allow you to download datasets from their website. This will also enable the tests included in "test_downloads.R", "test_extraction.R" and "test_ui.R" to run, as these tests first check for an rdhs config file ("rdhs.json") within the testing directory.

After having cloned/downloaded the repository, you will need to set up your rdhs config file within the testing directory, after which you should be able to run the test suite in full.


```r
pkg_dir <- "C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/rdhs"
setwd(file.path(pkg_dir,"tests/testthat/"))
rdhs::set_rdhs_config(email = "rdhs.tester@gmail.com",
                      project = "Testing Malaria Investigations", 
                      config_path = "rdhs.json", 
                      global = FALSE)

# devtools::test(pkg_dir)
```
---
title: "Interacting with geojson API returns"
author: "OJ Watson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Interacting with geojson API results}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

The API can return geojson objects, which can be very useful for quickly creating maps of the 
DHS indicator data. These are large objects though, and can cause the API to return slowly if
too much data is returned. 

In this demonstration we will be using 2 other packages, `leaflet` and `geojson`.

---

```{r geojson, warnings = FALSE, message = FALSE, fig.width=6, fig.height=6}

# load our package
library(rdhs)

# install other packages
install.packages("geojson")
install.packages("leaflet")

# make request
d <- dhs_data(countryIds = "SN",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2014,
              breakdown = "subnational",
              returnGeometry = TRUE,
              f = "geojson")

# convert to sp
m <- geojsonio::as.json(d)
nc <- geojsonio::geojson_sp(m) 

# plot using leaflet
pal <- leaflet::colorNumeric("viridis", NULL)

leaflet::leaflet(nc[nc$IndicatorId=="FE_FRTR_W_A15",]) %>%
  leaflet::addTiles() %>%
  leaflet::addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal(log10(Value)),
              label = ~paste0(CharacteristicLabel, ": ", formatC(Value, big.mark = ","))) %>%
  leaflet::addLegend(pal = pal, values = ~log10(Value), opacity = 1.0,
            labFormat = leaflet::labelFormat(transform = function(x) round(10^x)),title = ~Indicator[1])


```
---
title: "Anemia prevalence: an `rdhs` example"
author: "OJ Watson, Jeff Eaton"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    keep_md: TRUE
vignette: >
  %\VignetteIndexEntry{Anemia prevalence among women: an `rdhs` example}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Anemia is a common cause of fatigue, and women of childbearing age
are at particularly high risk for anemia. The package `rdhs` can be used
to compare estimates of the prevalence of any anemia among women from
Demographic and Health Surveys (DHS) conducted in Armenia, Cambodia,
and Lesotho.

# Setup

Load the `rdhs` package and other useful packages for analysing data.

```{r, results=FALSE, message = FALSE, warning=FALSE}
## devtools::install_github("ropensci/rdhs")
library(rdhs)
library(data.table)
library(ggplot2)
library(survey)
library(haven)
```

# Using calculated indicators from STATcompiler

Anemia prevalence among women is reported as a core indicator through the DHS STATcompiler (https://www.statcompiler.com/).
These indicators can be accessed directly from R via the DHS API with the function `dhs_data()`.

Query the API for a list of all StatCompiler indicators, and then search the indicators for those
that have `"anemia"` in the indicator name. API calls return `data.frame` objects, so if you prefer
to use `data.table` objects then convert afterwards, or we can set this up within our config using
`set_rdhs_config`.

```{r, R.options=list("rappdir_permission"=TRUE) , message = FALSE}
library(rdhs)
set_rdhs_config(data_frame = "data.table::as.data.table")

indicators <- dhs_indicators()
tail(indicators[grepl("anemia", Label), .(IndicatorId, ShortName, Label)])
```

The indicator ID `"AN_ANEM_W_ANY"` reports the percentage of women with any anemia.
The function `dhs_data()` will query the indicator dataset for the value of this indicator
for our three countries of interest. First, use `dhs_countries()` to query the
list of DHS countries to identify the DHS country code for each country.

```{r }
countries <- dhs_countries()
dhscc <- countries[CountryName %in% c("Armenia", "Cambodia", "Lesotho"), DHS_CountryCode]
dhscc
```

Now query the indicators dataset for the women with any anemia indicator for these three countries.

```{r ,fig.height=2.5, fig.width=4, fig.align="center"}
statcomp <- dhs_data(indicatorIds = "AN_ANEM_W_ANY", countryIds = dhscc)
statcomp[,.(Indicator, CountryName, SurveyYear, Value, DenominatorWeighted)]

ggplot(statcomp, aes(SurveyYear, Value, col=CountryName)) +
  geom_point() + geom_line()
```

# Analyse DHS microdata

## Identify surveys that include anemia testing

The DHS API provides the facility to filter surveys according to particular characteristics.
We first query the list of survey characteristics and identify the `SurveyCharacteristicID`
that indicates the survey included anemia testing. The first command below queries the API
for the full list of survey characteristics, and the second uses `grepl()` to search
`SurveyCharacteristicName`s that include the word 'anemia'.

```{r }
surveychar <- dhs_survey_characteristics()
surveychar[grepl("anemia", SurveyCharacteristicName, ignore.case=TRUE)]
```

The `SurveyCharacteristicID = 41` indicates that the survey included anemia testing. Next we
query the API to identify the surveys that have this characteristic and were conducted
in our countries of interest.

```{r }
surveys <- dhs_surveys(surveyCharacteristicIds = 41, countryIds = dhscc)
surveys[,.(SurveyId, CountryName, SurveyYear, NumberOfWomen, SurveyNum, FieldworkEnd)]
```

Finally, query the API identify the individual recode (IR) survey datasets for each of these surveys

```{r }
datasets <- dhs_datasets(surveyIds = surveys$SurveyId, fileType = "IR", fileFormat="flat")
datasets[, .(SurveyId, SurveyNum, FileDateLastModified, FileName)]
```

## Download datasets

To download datasets we need to first log in to our DHS account, by providing our credentials and setting up our configuration using `set_rdhs_config()`. This will require providing as arguments your `email` and `project` for which you want to download datasets from. You will then be prompted for your password. You can also specify a directory for datasets and API calls to be cached to using `cache_path`. In order to comply with CRAN, this function will also ask you for your permission to write to files outside your temporary directory, and you must type out the filename for the `config_path` - "rdhs.json". (See [introduction vignette](https://docs.ropensci.org/rdhs/articles/introduction.html) for specific format for config, or `?set_rdhs_config`).

```{r , R.options=list("rappdir_permission"=TRUE), message = FALSE }

## set up your credentials
set_rdhs_config(email = "jeffrey.eaton@imperial.ac.uk",
                project = "Joint estimation of HIV epidemic trends and adult mortality")

```

After this the function `get_datasets()` returns a list of file paths where the desired datasets are saved in the cache. The first time a dataset is accessed, `rdhs` will download the dataset from the DHS program website using the supplied credentials. Subsequently, datasets will be simply be located in the cached repository. 

```{r }
datasets$path <- unlist(get_datasets(datasets$FileName))
```

## Identify survey variables

Anemia is defined as having a hemoglobin (Hb) <12.0 g/dL for non-pregnant women
or Hb <11.0 g/dL for currently pregnant women^[https://www.measureevaluation.org/prh/rh_indicators/womens-health/womens-nutrition/percent-of-women-of-reproductive-age-with-anemia].
To calculate anemia prevalence from DHS microdata, we must identify the DHS recode
survey variables for hemoglobin measurement and pregnancy status. This could be
done by consulting the DHS recode manual or the .MAP files accompanying survey
datasets. It is convenient though to do this in R by loading the first
individual recode dataset and searching the metadata for the variable
names corresponding to the hemoglobin measurement and pregnancy status.

```{r }
head(search_variable_labels(datasets$FileName[10], "hemoglobin")[,1:2])
```

Variable `v042` records the household selection for hemoglobin testing.
Variable `v455` reports the outcome of hemoglobin measurement and `v456`
the result of altitude adjusted hemoglobin levels.

```{r }
ir <- readRDS(datasets$path[10])
table(as_factor(ir$v042))
table(as_factor(ir$v455))
summary(ir$v456)
```

Variable `v454` reports the current pregnancy status used for determining the
anemia threshold.

```{r }
search_variable_labels(datasets$FileName[1], "currently.*pregnant")[,1:2]
table(as_factor(ir$v454))
```

We also keep a number of other variables related to the survey design and potentially
interesting covariates: country code and phase (`v000`), cluster number (`v001`),
sample weight (`v005`), age (`v012`), region (`v024`), urban/rural residence (`v025`),
and education level (`v106`).

```{r }
vars <- c("SurveyId", "CountryName", "SurveyYear", "v000", "v001", "v005",
          "v012", "v024", "v025", "v106", "v042", "v454", "v455", "v456")
```

## Extract survey data

```{r }
datlst <- list()

for(i in 1:nrow(datasets)){

  if(file.exists(datasets$path[i])){
  
  print(paste(i, datasets$SurveyId[i]))
  ir <- readRDS(datasets$path[i])

  ir$SurveyId <- datasets$SurveyId[i]
  ir$CountryName <- datasets$CountryName[i]
  ir$SurveyYear <- datasets$SurveyYear[i]

  datlst[[datasets$SurveyId[i]]] <- ir[vars]
  }
}
```

We use `rbind_labelled()` to combine datasets with labelled columns. The argument
`labels` describes to combine variable levels for all datasets for `v024` (region)
while providing a consistent set of value labels to be used for `v454` (currently
pregnant) for all datasets.


```{r }
dat <- rbind_labelled(datlst,
                      labels = list(v024 = "concatenate",
                                    v454 = c("no/don't know" = 0L,
                                             "yes" = 1L, "missing" = 9L)))

sapply(dat, is.labelled)
dat$v456 <- zap_labels(dat$v456)
dat <- as_factor(dat)
```

## Data tabulations

It is a good idea to check basic tabulations of the data, especially by
survey to identify and nuances Exploratory analysis of variables

```{r }
with(dat, table(SurveyId, v025, useNA="ifany"))
with(dat, table(SurveyId, v106, useNA="ifany"))
with(dat, table(SurveyId, v454, useNA="ifany"))
with(dat, table(SurveyId, v455, useNA="ifany"))

with(dat, table(v042, v454, useNA="ifany"))
```

## Calculate anemia prevalence
Create indicator variable for 'any anemia'. The threshold depends on pregnancy status.

```{r }
dat$v456[dat$v456 == 999] <- NA
with(dat, table(v455, is.na(v456)))
dat$anemia <- as.integer(dat$v456  < ifelse(dat$v454 == "yes", 110, 120))
dat$anemia_denom <- as.integer(!is.na(dat$anemia))
```

Specify survey design using the `survey` package.

```{r }
dat$w <- dat$v005/1e6
des <- svydesign(~v001+SurveyId, data=dat, weights=~w)

anemia_prev <- svyby(~anemia, ~SurveyId, des, svyciprop, na.rm=TRUE, vartype="ci")
anemia_denom <- svyby(~anemia_denom, ~SurveyId, des, svytotal, na.rm=TRUE)

anemia_prev <- merge(anemia_prev, anemia_denom[c("SurveyId", "anemia_denom")])
res <- statcomp[,.(SurveyId, CountryName, SurveyYear, Value, DenominatorUnweighted, DenominatorWeighted)][anemia_prev, on="SurveyId"]

res$anemia <- 100*res$anemia
res$ci_l <- 100*res$ci_l
res$ci_u <- 100*res$ci_u
res$anemia_denom0 <- round(res$anemia_denom)
```

The table below compares the prevalence of any anemia calculated from survey microdata
with the estimates from DHS StatCompiler and the weighted denominators for each
calculation. The estimates are identical for most cases. There are some small
differences to be ironed out, which will require looking at the specific countries to check
how their stratification was carried out. (We are hoping to bring this feature in once the DHS
program has compiled how sample strata were constructed for all of their studies).

```{r fig.height=2.5, fig.width=4, fig.align="center"}
knitr::kable(res[,.(CountryName, SurveyYear, Value, anemia, ci_l, ci_u, DenominatorWeighted, anemia_denom0)], digits=1)


ggplot(res, aes(x=SurveyYear, y=anemia, ymin=ci_l, ymax=ci_u,
                col=CountryName, fill=CountryName)) +
  geom_ribbon(alpha=0.4, linetype="blank") + geom_point() + geom_line()
```


# Regression analysis: relationship between education and anemia

A key use of the survey microdata are to conduct secondary analysis of pooled data
from several surveys, such as regression analysis. Here we investigate the
relationship between anemia prevalence and education level (`v106`) for women using
logistic regression, adjusting for urban/rural (`v025`) and fixed effects for each
survey.

```{r }
des <- update(des, v106 = relevel(v106, "primary"))
summary(svyglm(anemia ~ SurveyId + v025 + v106, des, family="binomial"))
```

The results suggest that anemia prevalence is lower among women with higher education.
---
title: "Country Names and Codes"
author: "OJ Watson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Country Codes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## Overview

The following shows how to return a table of the country names and 2 letter codes. 

We can get this information by querying the API using `dhs_countries`, and
specifying the `returnFields` argument to just return the country name and 
DHS country code. This is a useful reference table when you want to pass in 
country IDs to a number of the API functions, e.g. `dhs_data(countryIds = "BJ")`


``` {r set config, message = FALSE}

library(rdhs)

## what are the countryIds
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
ids
```
---
title: "Toolkit"
author: "OJ Watson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Toolkit}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## Overview

The following is no longer needed as the CRAN version of haven is working. This vignette is still here though in case we need a non-CRAN pacakge in the future. 

---

In order to install `rdhs` you will require a working development environment. The following guide details how to do this across different operating systems. Please note this is only temporary while we wait for a development version of haven to make it to CRAN. 
Apologies for the inconvenience

1. **Windows**: Install **[Rtools](https://cran.r-project.org/bin/windows/Rtools/)**. For help on how to install **Rtools** please see the following [guide](https://cran.r-project.org/bin/windows/Rtools/), paying particular attention to the section about adding Rtools to your system `PATH`. 

In order to find out which version of **Rtools** you will need to check which version of R you are running. This can be be found out using the `sessionInfo()` function:

``` r 
sessionInfo()
R version 3.4.4 (2016-06-21)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1
```

2. **Mac**: Install Xcode from the Mac App Store.

3. **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).

Once a working development environment is ready then the following should work for you:

```r
# install.packages("devtools")
# devtools::install_github("ropensci/rdhs")
# library(rdhs)
```
---
title: "How to use rdhs?"
author: "OJ Watson, Jeff Eaton"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to use rdhs?}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

`rdhs` is a package for management and analysis of [Demographic and Health Survey (DHS)](https://www.dhsprogram.com/) data. This includes functionality to:

1. Access standard indicator data (i.e. [DHS STATcompiler](https://www.statcompiler.com/)) in R via the [DHS API](https://api.dhsprogram.com/).
1. Identify surveys and datasets relevant to a particular analysis.
1. Download survey datasets from the [DHS website](https://dhsprogram.com/data/available-datasets.cfm).
1. Load datasets and associated metadata into R.
1. Extract variables and combining datasets for pooled multi-survey analyses.

This process is described below and should cover most functionality that will be needed for working with these datasets. 

## 0. Installation

Install rdhs from github with `devtools`:

```{r gh-installation, warnings = FALSE, message = FALSE}
# install.packages("devtools")
# devtools::install_github("ropensci/rdhs")
library(rdhs)
```

---

 > Before starting the tutorial, if you wish to download survey datasets from the DHS website, you will need to set up an account with the DHS website, which will enable you to request access to the datasets. Instructions on how to do this can be found [here](https://dhsprogram.com/data/Access-Instructions.cfm). The email, password, and project name that were used to create the account will then need to be provided to `rdhs` when attempting to download datasets. You can still interact with the DHS API in section 1-2 below without having an account with the DHS website, however, you will need to create an account if you wish to go through steps 3-5. 

---

## 1. Access standard indicator data via the API

The DHS programme has published an API that gives access to a number of different data sets, which each represent one of the DHS API endpoints (e.g.  https://api.dhsprogram.com/rest/dhs/tags, or https://api.dhsprogram.com/rest/dhs/surveys). These data sets include the standard health indicators that are available within [DHS STATcompiler](https://www.statcompiler.com/) as well as a series of meta data sets that describe the types of surveys that have been conducted as well as which raw dataset files are available from which surveys. Each of these data sets are described within the [DHS API website](https://api.dhsprogram.com/), and there are currently 12 different data sets available from the API. Each of these data sets can be accessed using anyone of `dhs_<>()` functions.  All exported functions within `rdhs` that start *dhs_* interact with a different data set of the [DHS API](https://api.dhsprogram.com/). Their website gives great information about the different search terms and filters that can be used, and we have tried to include all of this within the documentation of each function. Each of these data sets

One of those functions, `dhs_data()`, interacts with the the published set of standard health indicator data calculated by the DHS. This data set contains a set of health indicators that have been sample weighted to give country, subnational estimates that can be further refined by education and wealth brackets. To do this we use the `dhs_data()` function, which we can then either search for specific indicators, or by querying for indicators that have been tagged within specific areas.

```{r inddata}
## what are the indicators
indicators <- dhs_indicators()
indicators[1,]

```

Each call to the DHS API returns a `data.frame` by default with all the results available by default. 

The DHS has a unique *IndicatorId* for each of the statistics it calculates. The definition and specific string for each indicator is included within the *IndicatorId* and *Definition* variable:

```{r inddata defs}
# grab the first 5 alphabetically
indicators[order(indicators$IndicatorId),][1:5,c("IndicatorId", "Definition")]

```


Since there are quite a lot of indicators, it might be easier to first query by tags. The DHS tags their indicators by what areas of demography and health they relate to, e.g. anaemia, literacy, malaria parasitaemia are all specific tags. First let's look at what the tags are, by interacting with the `dhs_tags()` function, before grabbing data that related to malaria parasitaemia in the DRC and Tanzania since 2010:

```{r tags}
# What are the tags
tags <- dhs_tags()

# Let's say we want to view the tags that relate to malaria
tags[grepl("Malaria", tags$TagName), ]

# and now let's then grab this data by specifying the countryIds and the survey year starts
data <- dhs_data(tagIds = 36,countryIds = c("CD","TZ"),breakdown="subnational",surveyYearStart = 2010)
data[1,]
```


Depending on your analysis this maybe more than enough detail. It is also worth mentioning that this data can also be accessed via [DHS STATcompiler](https://www.statcompiler.com/) if you prefer a click and collect version. However, hopefully one can see that selecting a lot of different indicators for multiple countries and breakdowns should be a lot easier using the `rdhs` API interaction. For example we can very quickly find out the trends in antimalarial use in Africa, and see if perhaps antimalarial prescription has decreased after RDTs were introduced (assumed 2010). 

```{r fig.height=4, fig.width=7, fig.align="center", warnings = FALSE}
# Make an api request
resp <- dhs_data(indicatorIds = "ML_FEVT_C_AML", surveyYearStart = 2010,breakdown = "subnational")

# filter it to 12 countries for space
countries  <- c("Angola","Ghana","Kenya","Liberia",
                "Madagascar","Mali","Malawi","Nigeria",
                "Rwanda","Sierra Leone","Senegal","Tanzania")

# and plot the results
library(ggplot2)
ggplot(resp[resp$CountryName %in% countries,],
       aes(x=SurveyYear,y=Value,colour=CountryName)) +
  geom_point() +
  geom_smooth(method = "glm") + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  ylab(resp$Indicator[1]) + 
  facet_wrap(~CountryName,ncol = 6) 

```


If we incorrectly entered a filter query (very possible), `rdhs` will let us know our request was invalid:

```{r API fail, include=TRUE, message = FALSE, error=TRUE, warning = FALSE}
# Make an api request
resp <- dhs_data(indicatorIds="ML_FEVT_C_AMasfafasfL",
                 surveyYearStart=202231231306,
                 breakdown="subParTyping")

```


## 2. Identify surveys relevant for further analysis

You may, however, wish to do more nuanced analysis than the API allows. The following 4 sections detail a very basic example of how to quickly identify, download and extract datasets you are interested in.

Let's say we want to get all DHS survey data from the Democratic Republic of Congo and Tanzania in the last 5 years (since 2013), which covers the use of rapid diagnostic tests (RDTs) for malaria. To begin we'll interact with the DHS API to identify our datasets.

To start our extraction we'll query the *surveyCharacteristics* data set using `dhs_survey_characteristics()` function:

```{r sc}
## make a call with no arguments
sc <- dhs_survey_characteristics()
sc[grepl("Malaria", sc$SurveyCharacteristicName), ]

```


There are 87 different survey characteristics, with one specific survey characteristic for Malaria RDTs. We'll use this to then find the surveys that include this characteristic. We can also at this point filter for our desired countries and years. The DHS API allows for countries to be filtered using by their *countryIds*, which is one of the arguments in `dhs_surveys()`. To have a look at what each countries countryId is we can use another of the API functions:

```{r surv}
## what are the countryIds
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
str(ids)

# lets find all the surveys that fit our search criteria
survs <- dhs_surveys(surveyCharacteristicIds = 89,
                     countryIds = c("CD","TZ"),
                     surveyType = "DHS",
                     surveyYearStart = 2013)

# and lastly use this to find the datasets we will want to download and let's download the flat files (.dat) datasets (have a look in the dhs_datasets documentation for all argument options, and fileformat abbreviations etc.)
datasets <- dhs_datasets(surveyIds = survs$SurveyId, 
                         fileFormat = "flat")
str(datasets)
```

Lastly, we recommended to download either the spss (.sav), `fileFormat = "SV"`, or the flat file (.dat), `fileFormat = "FL"` datasets. The flat is quicker, but there are still one or two datasets that don't read correctly, whereas the .sav files are slower to read in but so far no datasets have been found that don't read in correctly.

We can now use this to download our datasets for further analysis. 

## 3. Download survey datasets

We can now go ahead and download our datasets. To be able to download survey datasets from the DHS website, you will need to set up an account with them to enable you to request access to the datasets. Instructions on how to do this can be found [here](https://dhsprogram.com/data/Access-Instructions.cfm). The email, password, and project name that were used to create the account will then need to be provided to `rdhs` when attempting to download datasets. 

Once we have created an account, we need to set up our credentials using the function `set_rdhs_config()`. This will require providing as arguments your `email` and `project` for which you want to download datasets from. You will then be prompted for your password.

You can also specify a directory for datasets and API calls to be cached to using `cache_path`. If you do not provide an argument for `cache_path` you will be prompted to provide permission to `rdhs` to save your datasets and API calls within your user cache directory for your operating system. This is to comply with CRAN's requests for permission to be granted before writing to system files. If you do not grant permission, these will be written within your R temporary directory (as we saw above when we first used one of the functions to query the API). Similarly if you do not also provide an argument for `config_path`, this will be saved within your temp directory unless permission is granted. Your config files will always be called "rdhs.json", so that `rdhs` can find them easily.

```{r client, R.options=list("rappdir_permission"=TRUE)}

## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations")

```


Because you may have more than one project set up with the DHS website, you may want to have a separate directory for each set of datasets, and thus you will need to set up a different config file. To do this you need to set up a local config file. This can be achieved by setting the `global` param to `FALSE` (i.e. not global). You will also now need to provide the `config_path` argument, which MUST be **"rdhs.json"**. In order to comply with CRAN, you have to type this in (rather than have it as the default option).

```{r client2, R.options=list("rappdir_permission"=TRUE)}
## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations",
                config_path = "rdhs.json",
                cache_path = "project_one",
                global = FALSE)

```

You may, however, not have different projects with the DHS website, in which case you may prefer to set up one global config file. If you do not want this to be saved in your user cache directory, you can set `global` to `TRUE` (the default) and this will save it in your R default launch directory. This MUST be **~/.rdhs.json**. (There is not really any difference between saving it at `~/.rdhs.json` vs the user cache directory, but you might want to have it somewhere easy to find etc).

```{r client3, R.options=list("rappdir_permission"=TRUE)}
## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations",
                config_path = "~/.rdhs.json",
                global = TRUE)

```


After you have used `set_rdhs_config`, `rdhs` will try and find your config file when you next use `rdhs` in a different R session. It will do so by first looking locally for "rdhs.json", then globally for "~/.rdhs.json", then into your user cache directory, before lastly creating one in your temp directory. This is what was happening when you first used one of the API functions, and as such the config that is created to query the API initially will not be able to download datasets. 

Lastly, if you wish to return a `data.table` from your API requests, rather than a `data.frame` then you can change the default behaviour using the `data_frame` argument. You could also use this to convert them to `tibbles` and so on:

```{r client4, R.options=list("rappdir_permission"=TRUE)}
## set up your credentials
set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations",
                config_path = "~/.rdhs.json",
                data_frame = "data.table::as.data.table",
                global = TRUE)

```


To see what config that is being used by `rdhs` at any point, then use `get_rdhs_config()` to view the config settings.

---

Before we download our datasets, it is worth mentioning that once you have set up your login credentials, your API calls will be cached for your within the cache directory used. This will allow those working remotely or without a good internet connection to be able to still return previous API requests. If you do not, API requests will still be cached within the temp directory, so will be quick to be returned a second time, but they will then be deleted when you start a new R session.

```{r client_api_cache}
# the first time we call this function, rdhs will make the API request
microbenchmark::microbenchmark(dhs_surveys(surveyYear = 1992),times = 1)

# with it cached it will be returned much quicker
microbenchmark::microbenchmark(dhs_surveys(surveyYear = 1992), times = 1)

```


Now back to our dataset downloads. If we have a look back at our datasets object, we'll see there are 19 datasets listed. However, not all of them will be relevant to our malaria RDT questions. One approach is to head to the DHS website and have a look at the [DHS Recodes](https://dhsprogram.com/publications/publication-dhsg4-dhs-questionnaires-and-manuals.cfm), and look at the recodes that relate to the surveys. The other alternative is to download all the surveys and then query the variables within them. This is what we'll demonstrate here as it also demonstrates more of the package's functionality:

So first we will download all these datasets:

```{r download, message=FALSE}
# download datasets
downloads <- get_datasets(datasets$FileName)
```


The function returns a list with a file path to where the downloaded datasets have been saved to. By default the files will download quietly, i.e. no progress is shown. However, if you want to see the progress then you can control this by setting this in your config using the `verbose_download` argument.

## 4. Load datasets and associated metadata into R.

We can now examine what it is we have actually downloaded, by reading in one of these datasets:

```{r read a dataset}
# read in our dataset
cdpr <- readRDS(downloads$CDPR61FL)
```

The dataset returned here contains all the survey questions within the dataset.  The dataset is by default stored as a *labelled* class from the [haven package](https://github.com/tidyverse/haven). This class preserves the original semantics and can easily be coerced to factors with `haven::as_factor()`. Special missing values are also preserved. For more info on the *labelled* class have a look at their github.

So if we have a look at what is returned for the variable *hv024*:

```{r haven type}
head(cdpr$hv024)
```
```{r}
# and then the dataset
class(cdpr$hv024)
```

If we want to get the data dictionary for this dataset, we can use the function `get_variable_labels`, which will return what question each of the variables in our dataset refer to:

```{r probe dataset}
# let's look at the variable_names
head(get_variable_labels(cdpr))
```

For many of the survey responses this will give enough information for us to understand what the data is. However, for some questions it may be less clear exactly what the question means and how it may differ to other similar questions. If this is the case, then the DHS website publishes a lot of infrmation about the survey protocols and the surveys. We strongly advise for people to have a look through the [DHS website's documentation about using their datasets for analysis section](https://www.dhsprogram.com/data/Using-Datasets-for-Analysis.cfm), as well as the [recode files](https://www.dhsprogram.com/publications/publication-dhsg4-dhs-questionnaires-and-manuals.cfm) to understand how the surveys are carried out.

---

Above we saw that the default behaviour for the function `get_datasets` was to download the datasets, read them in, and save the resultant data.frame as a .rds object within the cache directory. You can control this behaviour using the `download_option` argument as such:

* `get_datasets(download_option = "zip")` - Just the downloaded zip will be saved
* `get_datasets(download_option = "rds")` - Just the read in rds will be saved
* `get_datasets(download_option = "both")` - The zip is downloaded and saved as well as the read in rds

The other main reason for reading the dataset in straight away as the default option is that `rdhs` will also create a table of all the survey variables and their labels (definitions) and cache them for you, which then allows us to quickly query for particular search terms or survey variables:

```{r questions, message = FALSE}
# rapid diagnostic test search
questions <- search_variable_labels(datasets$FileName, search_terms = "malaria rapid test")

table(questions$dataset_filename)
```


What we see from the questions is that the question "Result of malaria rapid test" appears in a few different datasets. This is because the household member recode datasets (CDPR61SV, TZPR7ASV) stores information about the children in a household, with one row per child, whereas the household recode (CDHR61SV, TZHR7ASV) stores information about the household, and thus flattens the information from each child into different subvariables (hml35$01/02 etc). As such it is easier to extract this information from the household member recodes. 

## 5. Extract variables and combining datasets for pooled multi-survey analyses.

To extract our data we pass our questions object to the function `extract_dhs`, which will create a list with each dataset and its extracted data as a `data.frame`. We also have the option to add any geographic data available, which will download the geographic data files for you and add this data to you resultant extract:

```{r extract_questions, message=FALSE}
# let's just use the PR files thus
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "FL", fileType = "PR")
downloads <- get_datasets(datasets$FileName)

# and grab the questions from this again along with also questions detailing the province
questions <- search_variable_labels(datasets$FileName, search_terms = c("malaria rapid test"))

# and now extract the data
extract <- extract_dhs(questions, add_geo = FALSE)

# what does our extract look like
str(extract)
```


The resultant extract is a list, with a new element for each different dataset that you have extracted. The responses from the dataset are by default stored as a *labelled* class from the [haven package](https://github.com/tidyverse/haven). 

We can also query our datasets for the survey question variables. In the example above the survey variable label was *Result of malaria rapid test* and the variable was *hml35*. So if you knew the survey variables that you wanted (either by looking at the Recode file or by looking through the *variable_names* included in the datasets) then we could search against these. So let's grab the regions using *hv024* using the client function `search_variables()`:

```{r extract_variables, message=FALSE}
# and grab the questions from this now utilising the survey variables
questions <- search_variables(datasets$FileName, variables = c("hv024","hml35"))

# and now extract the data
extract2 <- extract_dhs(questions, add_geo = FALSE)

# quick check 
head(extract2$CDPR61FL)
head(extract2$TZPR7AFL)

# and just to prove that hml35 did actually read in okay (there are just lots of NA)
table(extract2$CDPR61FL$hml35,useNA = "always")
```


We can now combine our two dataframes for further analysis using the `rdhs` package function `rbind_labelled()`. This function works specifically with our lists of labelled data.frames:

```{r rbind_labelled}
# first let's bind our first extraction, without the hv024
extract_bound <- rbind_labelled(extract)

head(extract_bound)

# now let's try our second extraction
extract2_bound <- rbind_labelled(extract2)

```


This hasn't quite done what we might want in the second instance. The *hv024* variable stores the regions for these 2 countries, which will not be the same and thus the labels will be different between the two of them. Without specifying any additional arguments `rbind_labelled()` will simply use the first data.frames labelling as the default, which will mean that some of the Tanzanian provinces will have been encoded as DRC provinces - not good! (This is a similar problem in nature to say trying to add new character strings to a factored data.frame).

There are a few work arounds. Firstly, we can specify a *labels* argument to the function which will detail how we should handle different variables. *labels* is a names list that specifies how to handle each variable. If we simply want to keep all the labels then we us the string "concatenate":

```{r concatenate}
# lets try concatenating the hv024
better_bound <- rbind_labelled(extract2, labels = list("hv024"="concatenate"))

head(better_bound$hv024)

```


We could also specify new labels for a variable. For example, imagine the two datasets encoded their RDT responses differently, with the first one as `c("No","Yes")` and the other as `c("Negative","Positive")`. These would be for our purposes the same response, and so we could either leave it and all our results would use the `c("No","Yes")` labelling. But we may want to use the latter as it's more informative/correct, or we may want to be crystal clear and use `c("NegativeTest","PositiveTest")`. we can do that like this:

```{r concatenate and new label}
# lets try concatenating the hv024 and providing new labels
better_bound <- rbind_labelled(
  extract2,
  labels = list("hv024"="concatenate",
                "hml35"=c("NegativeTest"=0, "PositiveTest"=1))
)

# and our new label
head(better_bound$hml35)
```


The other option is to not use the labelled class at all. We can control this when we download our datasets, using the argument `reformat=TRUE`. This will ensure that no factors or labels are used and it is just the raw data. When this option is set the object returned by `get_datasets()` no longer has any labelled classes or factors. However, we can still recover the variable table for a dataset using `get_variable_labels()`, which will take any dataset output by `get_datasets()` and return a data.frame describing the survey question variables and definitions.  

```{r reformat, message=FALSE}
# download the datasets with the reformat arguments
downloads <- get_datasets(datasets$FileName, reformat=TRUE)

# grab the questions but specifying the reformat argument
questions <- search_variables(datasets$FileName, variables = c("hv024", "hml35"),
                                     reformat=TRUE)

# and now extract the data
extract3 <- extract_dhs(questions, add_geo = FALSE)

# group our results
bound_no_labels <- rbind_labelled(extract3)

# what does our hv024 look like now
class(bound_no_labels$hv024[1])

```


The *hv024* column is now just characters, which is possibly the best option depending on your downstream analysis/preferences. It's for this reason that the geographic data that is added is never turned into factors or labels.  

Lastly, we can now use our extract dataset to carry out some regression analysis, to investigate the relationship between malaria prevalence and the quality of wall materials. To do this we will need to first grab the sample weights and stratification from the surveys, along with the extra variables and we will then check the RDT prevalence calculated using the raw data versus the API:

```{r regression analysis, message=FALSE, warning=FALSE}
# grab the additional variable hv023 and hv024 which have the strata and weights respectively, and hc1 which is the age
questions <- search_variables(datasets$FileName,variables = c("hv005","hv021","hv022","hv023","hv024",
                                                              "hv025","hv214","hml20", "hc1","hml35"))
extraction <- extract_dhs(questions,TRUE)

# now concatenate the provinces as before and remove missing responses
dat <- rbind_labelled(extraction,labels=list("hv024"="concatenate","hv214"="concatenate"))
dat <- dat[-which(dat$hml35==9),] # remove missing responses

# and we are going to compare our extract to the API malaria prevalence by RDT, which is for those between 6 and 59 months
dat <- dat[-which(!dat$hc1>6 & dat$hc1<=60),]

# create a denominator response for hml35
dat$hml35denom <- as.integer(!is.na(dat$hml35))
dat$bricks <- dat$hv214 %in% c(8,18,5,9,10)
dat$net <- as.logical(dat$hml20)

# specify the strata and sample weights
dat$strata <- paste0(dat$hv023,dat$DATASET)
dat$hv005 <- dat$hv005/1e6

# construct a survey design using the survey pacakge
library(survey)

# construct the sample design and calculate the mean and totals 
des <-  survey::svydesign(~CLUSTER+DATASET,data=dat,weight=~hv005)
results <- cbind(survey::svyby(~hml35,by=~DHSREGNA+DATASET, des, survey::svyciprop,na.rm=TRUE),
                 survey::svyby(~hml35denom,by=~DHSREGNA+DATASET, des, survey::svytotal,na.rm=TRUE))
results <- results[order(results$DATASET),]

# grab the same data from the API 
dhs_api_data <- dhs_data(countryIds = c("CD","TZ"),indicatorIds = "ML_PMAL_C_RDT",breakdown = "subnational",surveyYearStart = 2013, surveyYearEnd = 2016)
dhs_api_data <- cbind(dhs_api_data$Value,dhs_api_data$DenominatorWeighted,dhs_api_data$CharacteristicLabel, dhs_api_data$SurveyId)
api <- dhs_api_data[!grepl("\\.\\.",dhs_api_data[,3]),] # remove subregions included in Tanzania
api <- api[order(apply(api[,4:3],1,paste,collapse="")),]

# bind the results and remove duplicate Region Columns
comparison <- cbind(results[,c(1,3,7)],api[])
names(comparison) <- c("Region","Survey_RDT_Prev","Survey_RDT_Denom","API_RDT_Prev","API_RDT_Denom","API_Regions","SurveyID")
head(comparison[,c(1,2,4,3,5,7)])
```


It's a little off, with the mean values differing due to maybe the specific cut off they used in terms of which ages were included within between 5 and 69. The variance could also be off due to the specific stratification the DHS Program will have used, as well as potentially how they have grouped the primary sampling units. We are hoping to get this information from the DHS for each survey so we can make this process more streamline for the end user. 

And lastly we will construct a logistic regression to investigate the relationship between a positive malaria RDT and whether the main walls of an individual's house were made of bricks or similar, while adjusting for urban/rural (`hv025`) and fixed effects for each survey.

```{r altitutde glm}
# contsruct our glm using svyglm and specify quasibinomial to handle the na in hml35
summary(svyglm(hml35 ~ DATASET + hv025 + net + bricks, des, family="quasibinomial"))

```


What we can see is that a significant negative gradient was associated with walls being made of bricks or similarly good materials in comparison to malaria positivity rates by RDT. What is also interesting is that whether the individual slept under a long lasting insecticidal net (*hml20* that we converted to net) was not significant.  

---

## Summary and further thoughts

Hopefully the above tutorial has shown how the `rdhs` package can facilitate both querying the DHS API and hopefully make downloading and interacting with the raw datasets a smoother, more reproducible process. It is worth bearing in mind though, that creating a harmonised dataset is not always as easy as the example above - a lot of the time survey variables differ across years and surveys, which is hopefully when the `survey_questions` functionality will make it easier to first filter down to those that include the relevant questions before having to decide which survey questions are valid. 

Any suggestions or comments/corrections/errors/ideas please let me know either in the issues or send me an email at "o.watson15@imperial.ac.uk". And if there is any further functionality that you think you would be useful, then also let me know. :)
---
title: "Downloading Shape Files for DHS Surveys"
author: "OJ Watson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Downloading shape files for DHS surveys}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

`rdhs` can be used to download shape files associated with the DHS surveys. For
example, we can link API responses to their shape files and plot the subnational
estimates using the spatial pacakge `sp`. To do this, we use the new function
`download_boundaries` from `rdhs` to download shape files. 

```{r normal, warnings = FALSE, message = FALSE, fig.width=6, fig.height=6}

# load our package
library(rdhs)

# make request
d <- dhs_data(countryIds = "SN",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2014,
              breakdown = "subnational")

# get our related spatial data frame object
sp <- download_boundaries(surveyId = d$SurveyId[1])

# match our values to the regions
m <- d$Value[match(sp$sdr_subnational_boundaries$REG_ID, d$RegionId)]
sp$sdr_subnational_boundaries@data$Value <- m

# Use sp to plot
sp::spplot(sp$sdr_subnational_boundaries, "Value", main = d$Indicator[1])

```

Or we can use `sf` for our plotting, which offers more user-friendly plotting
options. 

```{r sf, warnings = FALSE, message = FALSE, fig.width=6, fig.height=6}

# make request
d <- dhs_data(countryIds = "SN",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2017,
              breakdown = "subnational")

# get our related spatial data frame object
sp <- download_boundaries(surveyId = d$SurveyId[1], method = "sf")

# match our values to the regions
m <- d$Value[match(sp$sdr_subnational_boundaries$REG_ID, d$RegionId)]
sp$sdr_subnational_boundaries$Value <- m

# Use ggplot and geom_sf to plot
library(ggplot2)
ggplot(sp$sdr_subnational_boundaries) + 
  geom_sf(aes(fill = Value)) + 
  ggtitle(d$Indicator[1])

```
---
title: "rdhs client object and internal package design"
author: "OJ Watson, Jeff Eaton"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using the rdhs client}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview

The following is a similar vignette - "How to use rdhs?". However, this one achieves the same results using the `rdhs` client object directly to download datasets etc, rather than using the easier user interfacing functions. This vignette is included to help demonstrate how the pacakge internals are working and adds more information about how the package was designed, which may be useful for anyone adding extra functionality to `rdhs` in the future. 

If are not bothered by this then please read the main introductory vignette. 

---

`rdhs` is a package for management and analysis of [Demographic and Health Survey (DHS)](https://www.dhsprogram.com/) data. This includes functionality to:

1. Access standard indicator data (i.e. [DHS STATcompiler](https://www.statcompiler.com/)) in R via the [DHS API](https://api.dhsprogram.com/).
1. Identify surveys and datasets relevant to a particular analysis.
1. Download survey datasets from the [DHS website](https://dhsprogram.com/data/available-datasets.cfm).
1. Load datasets and associated metadata into R.
1. Extract variables and combining datasets for pooled multi-survey analyses.

This process is described below and should cover most functionality that will be needed for working with these datasets. 

## 0. Installation

Install rdhs from github with `devtools`:

```{r gh-installation, warnings = FALSE, message = FALSE}
# install.packages("devtools")
# devtools::install_github("ropensci/rdhs")

library(rdhs)
```


## 1. Access standard indicator data via the API

The DHS API has many *endpoints* that can be accessed using anyone of `dhs_<endpoint>()` functions. All exported functions within `rdhs` that start *dhs_* interact with a different with a different endpoint of the [DHS API](https://api.dhsprogram.com/). Their website gives great information about the different search terms and filters that can be used, and we have tried to include all of this within the documentation of each function.

One of those endpoints, `dhs_data()`, interacts with the the published set of standard health indicator data calculated by the DHS. This endpoint allows use to retrieve a set of health indicators that have been sample weighted to give country, subnational estimates that can be further refined by education and wealth brackets. To do this we use the `dhs_data()` endpoint, which we can then either search for specific indicators, or by querying for indicators that have been tagged within specific areas.

```{r inddata, message = FALSE}
## what are the indicators
indicators <- dhs_indicators()
indicators[1,]

```


Each call to a DHS endpoint returns a `data.frame` by default with all the results available by default. 

The DHS has a unique *IndicatorId* for each of the statistics it calculates. The definition and specific string for each indicator is included within the *IndicatorId* and *Definition* variable:

```{r inddata defs}
# grab the first 5 alphabetically
indicators[order(indicators$IndicatorId),][1:5,c("IndicatorId", "Definition")]

```


Since there are quite a lot of indicators, it might be easier to first query by tags. The DHS tags their indicators by what areas of demography and health they relate to, e.g. anaemia, literacy, malaria parasitaemia are all specific tags. First let's look at what the tags are, by interacting with the `dhs_tags()` endpoint, before grabbing data that related to malaria parasitaemia in the DRC and Tanzania since 2010:

```{r tags}
# What are the tags
tags <- dhs_tags()

# Let's say we want to view the tags that relate to malaria
tags[grepl("Malaria", tags$TagName), ]

# and now let's then grab this data by specifying the countryIds and the survey year starts
data <- dhs_data(tagIds = 36,countryIds = c("CD","TZ"),breakdown="subnational",surveyYearStart = 2010)
data[1,]
```


Depending on your analysis this maybe more than enough detail. It is also worth mentioning that this data can also be accessed via [DHS STATcompiler](https://www.statcompiler.com/) if you prefer a click and collect version. However, hopefully one can see that selecting a lot of different indicators for multiple countries and breakdowns should be a lot easier using the `rdhs` API interaction. For example we can very quickly find out the trends in antimalarial use in Africa, and see if perhaps antimalarial prescription has decreased after RDTs were introduced (assumed 2010). 

```{r fig.height=4, fig.width=7, fig.align="center", warnings = FALSE}
# Make an api request
resp <- dhs_data(indicatorIds = "ML_FEVT_C_AML", surveyYearStart = 2010,breakdown = "subnational")

# filter it to 12 countries for space
countries  <- c("Angola","Ghana","Kenya","Liberia",
                "Madagascar","Mali","Malawi","Nigeria",
                "Rwanda","Sierra Leone","Senegal","Tanzania")

# and plot the results
library(ggplot2)
ggplot(resp[resp$CountryName %in% countries,],
       aes(x=SurveyYear,y=Value,colour=CountryName)) +
  geom_point() +
  geom_smooth(method = "glm") + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  ylab(resp$Indicator[1]) + 
  facet_wrap(~CountryName,ncol = 6) 

```


If we incorrectly entered a filter query (very possible), `rdhs` will let us know our request was invalid:

```{r API fail, include=TRUE, message = FALSE, error=TRUE, warning = FALSE}
# Make an api request
resp <- dhs_data(indicatorIds="ML_FEVT_C_AMasfafasfL",
                 surveyYearStart=202231231306,
                 breakdown="subParTyping")

```


## 2. Identify surveys relevant for further analysis

You may, however, wish to do more nuanced analysis than the API allows. The following 4 section detail a very basic example of how to quickly identify, download and extract datasets you are interested in.

Let's say we want to get all the survey data from the Democratic Republic of Congo and Tanzania in the last 5 years (since 2013), which covers the use of rapid diagnostic tests (RDTs) for malaria. To begin we'll interact with the DHS API to identify our datasets.

To start our extraction we'll query the *surveyCharacteristics* endpoint using `dhs_surveyCharacteristics()`:

```{r sc}
## make a call with no arguments
sc <- dhs_survey_characteristics()
sc[grepl("Malaria", sc$SurveyCharacteristicName), ]

```


There are 87 different survey characteristics, with one specific survey characteristic for Malaria RDTs. We'll use this to then find the surveys that include this characteristic. We can also at this point filter for our desired countries and years. The DHS API allows for countries to be filtered using by their *countryIds*, which is one of the arguments in `dhs_surveys()`. To have a look at what each countries countryId is we can use another of the API endpoints first:

```{r surv}
## what are the countryIds
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
str(ids)

# lets find all the surveys that fit our search criteria
survs <- dhs_surveys(surveyCharacteristicIds = 89,countryIds = c("CD","TZ"),surveyYearStart = 2013)

# and lastly use this to find the datasets we will want to download and let's download the flat files (.dat) datasets (have a look in the dhs_datasets documentation for all argument options, and fileformat abbreviations etc.)
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "flat")
str(datasets)
```

Lastly, we recommended to download either the spss (.sav), `fileFormat = "SV"`, or the flat file (.dat), `fileFormat = "FL"` datasets. The flat is quicker, but there are still one or two datasets that don't read correctly, whereas the .sav files are slower to read in but so far no datasets have been found that don't read in correctly.

We can now use this to download our datasets for further analysis. 

## 3. Download survey datasets

We can now go ahead and download our datasets. To do this we need to first create a `client`. The client is an R6 class (similar to R's built in reference classes and make caching survey and API queries more reproducible) and will be used to log in to your DHS account, download datasets for you, and help query those datasets for the question you are interested in. The client will also cache all of these processes, which really helps increase the reproducibility of your analysis. 

In order to set up our credentials we use the function `set_rdhs_config()`. This will require providing as arguments your `email` and `project` for which you want to download datasets from. You will then be prompted for your password. 

You can also specify a directory for datasets and API calls to be cached to using `cache_path`. If you do not provide an argument for `cache_path` you will be prompted to provide permission to `rdhs` to save your datasets and API calls within your user cache directory for your operating system. This is to comply with CRAN's requests for permission to be granted before writing to system files. If you do not grant permission, these will be written within your R temporary directory (as we saw above when we first used one of the functions to query the API). Similarly if you do not also provide an argument for `config_path`, this will be saved within your temp directory unless permission is granted. Your config files will always be called "rdhs.json", so that `rdhs` can find them easily.

```{r client, R.options=list("rappdir_permission"=TRUE), message = FALSE}
## create a client
## set up your credentials
config <- set_rdhs_config(email = "rdhs.tester@gmail.com",
                project = "Testing Malaria Investigations")

client <- client_dhs(config)
client

```


Before we use our client to download our datasets, it is worth mentioning that the client can be passed as an argument to any of the API functions we have just seen. This will then cache the results for you, so that if you are working remotely or without a good internet connection you can still return your previous API requests. This is what is happening behind the scenes anyway, as `rdhs` will create a client for you when you call an API function (if one does not already exist).

---

Now back to our dataset downloads. If we have a look back at our datasets object, we'll see there are 19 datasets listed. However, not all of them will be relevant to our malaria RDT questions. One approach is to head to the DHS website and have a look at the [DHS Recodes](https://dhsprogram.com/publications/publication-dhsg4-dhs-questionnaires-and-manuals.cfm), and look at the recodes that relate to the surveys. The other alternative is to download all the surveys and then query the variables within them. This is what we'll demonstrate here as it also demonstrates more of the package's functionality:

So first we will download all these datasets:

```{r download, message=FALSE}
# download datasets
downloads <- client$get_datasets(datasets$FileName)
```

The function returns a list with a file path to where the downloaded datasets have been saved to. By default the files will download quietly, i.e. no progress is shown. However, if you want to see the progress then you can control this by setting this in your config using the `verbose_download` argument.

## 4. Load datasets and associated metadata into R.

We can now examine what it is we have actually downloaded, by reading in one of these datasets:

```{r read a dataset}
# read in our dataset
cdpr <- readRDS(downloads$CDPR61FL)
```


The dataset returned here contains all the survey questions within the dataset, and what their survey variables are:

```{r probe dataset}
# let's look at the variable_names
head(get_variable_labels(cdpr))

# and then the dataset
class(cdpr$hv024)
```


This is the default behaviour for the client function `get_datasets` - it will download the datasets for you, and then by default save them in your client's root directory and then unzip them and read them in for you, and save the resultant data.frame as a .rds object within the client's root directory. You can control this behaviour using the `download_option` argument as such:

* `client$get_datasets(download_option = "zip")` - Just the downloaded zip will be saved
* `client$get_datasets(download_option = "rds")` - Just the read in rds will be saved
* `client$get_datasets(download_option = "both")` - The zip is downloaded and saved as well as the read in rds

The other main reason for reading the dataset in straight away as the default option is that the created table of all the survey variables and their definitions is cached then and there, which then allows us to quickly query for particular search terms or survey variables:

```{r questions, message = FALSE}
# rapid diagnostic test search
questions <- client$survey_questions(datasets$FileName, search_terms = "malaria rapid test")

table(questions$dataset_filename)
```


What we see from the questions is that the question "Result of malaria rapid test" appears in a few different datasets. This is because the household member recode datasets (CDPR61SV, TZPR7HSV) stores information about the children in a household, with one row per child, whereas the household recode (CDHR61SV, TZHR7HSV) stores information about the household, and thus flattens the information from each child into different subvariables (hml35$01/02 etc). As such it is easier to extract this information from the household member recodes. 

## 5. Extract variables and combining datasets for pooled multi-survey analyses.

To extract our data we pass our questions object to the client function `extract`, which will create a list with each dataset and its extracted data as a `data.frame`. We also have the option to add any geographic data available, which will download the geographic data files for you and add this data to you resultant extract:

```{r extract_questions, message=FALSE}
# let's just use the PR files thus
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileFormat = "FL", fileType = "PR")
downloads <- client$get_datasets(datasets$FileName)

# and grab the questions from this again along with also questions detailing the province
questions <- client$survey_questions(datasets$FileName, search_terms = c("malaria rapid test"))

# and now extract the data
extract <- client$extract(questions, add_geo = FALSE)

# what does our extract look like
str(extract)
```


The resultant extract is a list, with a new element for each different dataset that you have extracted. The responses from the dataset are by default stored as a *labelled* class from the [haven package](https://github.com/tidyverse/haven). This class preserves the original semantics and can easily be coerced to factors with `haven::as_factor()`. Special missing values are also preserved. For more info on the *labelled* class have a look at their github.

We can also query our datasets for the survey question variables. In the example above the survey question was *Result of malaria rapid test* and the variable was *hml35*. So if you knew the survey variables that you wanted (either by looking at the Recode file or by looking through the *variable_names* included in the datasets) then we could search against these. So let's grab the regions using *hv024* using the client function `survey_variables()`:

```{r extract_variables, message=FALSE}
# and grab the questions from this now utilising the survey variables
questions <- client$survey_variables(datasets$FileName, variables = c("hv024","hml35"))

# and now extract the data
extract2 <- client$extract(questions, add_geo = FALSE)

# quick check 
head(extract2$CDPR61FL)
head(extract2$TZPR7HFL)

# and just to prove that hml35 did actually read in okay (there are just lots of NA)
table(extract2$CDPR61FL$hml35,useNA = "always")
```


We can now combine our two dataframes for further analysis using the `rdhs` package function `rbind_labelled()`. This function works specifically with our lists of labelled data.frames:

```{r rbind_labelled}
# first let's bind our first extraction, without the hv024
extract_bound <- rbind_labelled(extract)

head(extract_bound)

# now let's try our second extraction
extract2_bound <- rbind_labelled(extract2)

```


This hasn't quite done what we might want in the second instance. The *hv024* variable stores the regions for these 2 countries, which will not be the same and thus the labels will be different between the two of them. Without specifying any additional arguments `rbind_labelled()` will simply use the first data.frames labelling as the default, which will mean that some of the Tanzanian provinces will have been encoded as DRC provinces - not good! (This is a similar problem in nature to say trying to add new character strings to a factored data.frame).

There are a few work arounds. Firstly, we can specify a *labels* argument to the function which will detail how we should handle different variables. *labels* is a names list that specifies how to handle each variable. If we simply want to keep all the labels then we us the string "concatenate":

```{r concatenate}
# lets try concatenating the hv024
better_bound <- rbind_labelled(extract2, labels = list("hv024"="concatenate"))

head(better_bound$hv024)

```


We could also specify new labels for a variable. For example, imagine the two datasets encoded their RDT responses differently, with the first one as `c("No","Yes")` and the other as `c("Negative","Positive")`. These would be for our purposes the same response, and so we could either leave it and all our results would use the `c("No","Yes")` labelling. But we may want to use the latter as it's more informative/correct, or we may want to be crystal clear and use `c("NegativeTest","PositiveTest")`. we can do that like this:

```{r concatenate and new label}
# lets try concatenating the hv024 and providing new labels
better_bound <- rbind_labelled(extract2,
                               labels = list("hv024"="concatenate",
                                             "hml35"=c("NegativeTest"=0, "PositiveTest"=1)))

# and our new label
head(better_bound$hml35)
```


The other option is to not use the labelled class at all. We can control this when we download our datasets, using the argument `reformat=TRUE`. This will ensure that no factors or labels are used and it is just the raw data. When this option is set the object returned by `client$get_datasets()` no longer has any labelled classes or factors. However, we can still recover the variable table for a dataset using `get_variable_labels()`, which will take any dataset output by `get_datasets()` and return a data.frame describing the survey question variables and definitions.  

```{r reformat, message=FALSE}
# download the datasets with the reformat arguments
downloads <- client$get_datasets(datasets$FileName, reformat=TRUE)

# grab the questions but specifying the reformat argument
questions <- client$survey_variables(datasets$FileName, variables = c("hv024", "hml35"),
                                     reformat=TRUE)

# and now extract the data
extract3 <- client$extract(questions, add_geo = FALSE)

# group our results
bound_no_labels <- rbind_labelled(extract3)

# what does our hv024 look like now
class(bound_no_labels$hv024[1])

```


The *hv024* column is now just characters, which is possibly the best option depending on your downstream analysis/preferences. It's for this reason that the geographic data that is added is never turned into factors or labels.  

Lastly, we can now use our extract dataset to carry out some regression analysis, to investigate the relationship between malaria prevalence and the quality of wall materials. To do this we will need to first grab the sample weights and stratification from the surveys, along with the extra variables and we will then check the RDT prevalence calculated using the raw data versus the API:

```{r regression analysis, message=FALSE}
# grab the additional variable hv023 and hv024 which have the strata and weights respectively
questions <- client$survey_variables(datasets$FileName,variables = c("hv005","hv021","hv022","hv023","hv024",
                                                              "hv025","hv214","hml20", "hc1","hml35"))
extraction <- client$extract(questions,TRUE)

# now concatenate the provinces as before and remove missing responses
dat <- rbind_labelled(extraction,labels=list("hv024"="concatenate","hv214"="concatenate"))
dat <- dat[-which(dat$hml35==9),] # remove missing responses

# and we are going to compare our extract to the API malaria prevalence by RDT, which is for those between 6 and 59 months
dat <- dat[-which(!dat$hc1>6 & dat$hc1<=60),]

# create a denominator response for hml35
dat$hml35denom <- as.integer(!is.na(dat$hml35))
dat$bricks <- dat$hv214 %in% c(8,18,5,9,10)
dat$net <- as.logical(dat$hml20)

# specify the strata and sample weights
dat$strata <- paste0(dat$hv023,dat$DATASET)
dat$hv005 <- dat$hv005/1e6

# construct a survey design using the survey pacakge
library(survey)

# construct the sample design and calculate the mean and totals 
des <-  survey::svydesign(~CLUSTER+DATASET,data=dat,weight=~hv005)
results <- cbind(survey::svyby(~hml35,by=~DHSREGNA+DATASET, des, survey::svyciprop,na.rm=TRUE),
                 survey::svyby(~hml35denom,by=~DHSREGNA+DATASET, des, survey::svytotal,na.rm=TRUE))
results <- results[order(results$DATASET),]

# grab the same data from the API 
dhs_api_data <- dhs_data(countryIds = c("CD","TZ"),indicatorIds = "ML_PMAL_C_RDT",breakdown = "subnational",surveyYearStart = 2013, surveyYearEnd = 2016)
dhs_api_data <- cbind(dhs_api_data$Value,dhs_api_data$DenominatorWeighted,dhs_api_data$CharacteristicLabel, dhs_api_data$SurveyId)
api <- dhs_api_data[!grepl("\\.\\.",dhs_api_data[,3]),] # remove subregions included in Tanzania
api <- api[order(apply(api[,4:3],1,paste,collapse="")),]

# bind the results and remove duplicate Region Columns
comparison <- cbind(results[,c(1,3,7)],api[])
names(comparison) <- c("Region","Survey_RDT_Prev","Survey_RDT_Denom","API_RDT_Prev","API_RDT_Denom","API_Regions","SurveyID")

head(comparison[,c(1,2,4,3,5,7)])
```


It's a little off, which will be due to the specific stratification the DHS Program will have used, as well as potentially how they have grouped the primary sampling units. We are hoping to get this information from the DHS for each survey so we can make this process more streamline again.

And lastly we will construct a logistic regression to investigate the relationship between a positive malaria RDT and whether the main walls of an individual's house were made of bricks or similar, while adjusting for urban/rural (`v025`) and fixed effects for each survey.

```{r altitutde glm}
# contsruct our glm using svyglm and specify quasibinomial to handle the na in hml35
summary(svyglm(hml35 ~ DATASET + hv025 + net + bricks, des, family="quasibinomial"))

```


What we can see is that a significant negative gradient was associated with walls being made of bricks or similarly good materials in comparison to malaria positivity rates by RDT. What is also interesting is that whether the individual slept under a long lasting insecticidal net (*hml20* that we converted to net) was not significant.  

---

## Summary and further thoughts

Hopefully the above tutorial has shown how the `rdhs` package can facilitate both querying the DHS API and hopefully make downloading and interacting with the raw datasets a smoother, more reproducible process. It is worth bearing in mind though, that creating a harmonised dataset is not always as easy as the example above - a lot of the time survey variables differ across years and surveys, which is hopefully when the `survey_questions` functionality will make it easier to first filter down to those that include the relevant questions before having to decide which survey questions are valid. 

Any suggestions or comments/corrections/errors/ideas please let me know either in the issues or send me an email at "o.watson15@imperial.ac.uk". And if there is any further functionality that you think you would be useful, then also let me know. :)
---
title: "Running test suite locally"
author: "OJ Watson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Toolkit}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## Overview

The following is a guide describing how to be able to run the full test suite locally.

---

## 1. Creating an account with the DHS website

In order to run all the tests in the *rdhs* package you need to have created an account with the DHS website. This requires about 5 mins to fill in their form and then a few days for them to approve your account creation. To start this head to the DHS website and go to the Data tab [here](https://dhsprogram.com/data/dataset_admin/login_main.cfm). A full rationale for this procss can be found on their website [her](https://dhsprogram.com/data/Access-Instructions.cfm). In short you provide your contact information and then describe the stuy you are conducting that requires the DHS survey datasets.

Importantly, when filling in their form describing the description of your study be sure to include some mention of needing geographic datasets as a number of the tests download geographic datasets. After you have filled this in you need to specify which countries and datasets you need. In order for the full test suite to work you must request every country and dataset type possible. After you have done this the DHS programme will take a few days to approve you request and will then email you to let you know at the email address you have provided. 

## 2. Setting up your rdhs config

Having now created your account, you can set up your *rdhs* configuration, which will allow you to download datasets from their website. This will also enable the tests included in "test_downloads.R", "test_extraction.R" and "test_ui.R" to run, as these tests first check for an rdhs confing file ("rdhs.json") within the testing directory.

After having cloned/downloaded the repository, you will need to set up your rdhs config file within the testing directory, after which you should be able to run the test suite in full.

``` {r set config, R.options=list("rappdir_permission"=TRUE), message = FALSE}
pkg_dir <- "C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/rdhs"
setwd(file.path(pkg_dir,"tests/testthat/"))
rdhs::set_rdhs_config(email = "rdhs.tester@gmail.com",
                      project = "Testing Malaria Investigations", 
                      config_path = "rdhs.json", 
                      global = FALSE)

# devtools::test(pkg_dir)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbind_labelled.R
\name{as_factor.labelled}
\alias{as_factor.labelled}
\title{Archived dataset capable as_factor}
\usage{
as_factor.labelled(
  x,
  levels = c("default", "labels", "values", "both"),
  ordered = FALSE,
  ...
)
}
\arguments{
\item{x}{Object to coerce to a factor.}

\item{levels}{How to create the levels of the generated factor:
\itemize{
\item "default": uses labels where available, otherwise the values.
Labels are sorted by value.
\item "both": like "default", but pastes together the level and value
\item "label": use only the labels; unlabelled values become \code{NA}
\item "values: use only the values
}}

\item{ordered}{If \code{TRUE} create an ordered (ordinal) factor, if
\code{FALSE} (the default) create a regular (nominal) factor.}

\item{...}{Other arguments passed down to method.}
}
\description{
Changes in `haven` have meant that `labelled` class are now referred to as
`haven_labelled` classes. If `haven::as_factor` is used on old datasets they
will fail to find the suitable method. rdhs::as_factor.labelled will work on
old archived datasets that have a `labelled` class.
}
\details{
For more details see \code{haven::as_factor}
}
\examples{
\dontrun{
# create a data.frame using the new haven_labelled class
df1 <- data.frame(
area = haven::labelled(c(1L, 2L, 3L), c("reg 1"=1,"reg 2"=2,"reg 3"=3)),
climate = haven::labelled(c(0L, 1L, 1L), c("cold"=0,"hot"=1))
)

# manually change it to the old style
class(df1$area) <- "labelled"
class(df1$climate) <- "labelled"

# with rdhs attached, i.e. library(rdhs), we can now do the following
haven::as_factor(df1$area)

# we can also use this on the data.frame by using the only_labelled argument
haven::as_factor(df1, only_labelled = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{get_datasets}
\alias{get_datasets}
\title{Get Datasets}
\usage{
get_datasets(
  dataset_filenames,
  download_option = "rds",
  reformat = FALSE,
  all_lower = TRUE,
  output_dir_root = NULL,
  clear_cache = FALSE,
  ...
)
}
\arguments{
\item{dataset_filenames}{The desired filenames to be downloaded. These can be
found as one of the returned fields from \code{\link{dhs_datasets}}.
Alternatively you can also pass the desired rows from
code{\link{dhs_datasets}}.}

\item{download_option}{Character specifying whether the dataset should be
just downloaded ("zip"), imported and saved as an .rds object ("rds"), or
both extract and rds ("both"). Conveniently you can just specify any letter
from these options.}

\item{reformat}{Boolean concerning whether to reformat read in datasets by
removing all factors and labels. Default = FALSE.}

\item{all_lower}{Logical indicating whether all value labels should be lower
case. Default to `TRUE`.}

\item{output_dir_root}{Root directory where the datasets will be stored
within. The default will download
datasets to a subfolder of the client root called "datasets"}

\item{clear_cache}{Should your available datasets cache be cleared first.
This will allow newly accessed datasets to be available. Default = `FALSE`}

\item{...}{Any other arguments to be passed to \code{\link{read_dhs_dataset}}}
}
\value{
Depends on the download_option requested, but ultimately it is a file
  path to where the dataset was downloaded to, so that you can interact with
  it accordingly.
}
\description{
Downloads datasets you have access to from the DHS website
}
\details{
Gets datasets from your cache or downloads from the DHS website.
  By providing the filenames, as specified in one of the returned fields from
  \code{\link{dhs_datasets}}, the client will log in for you and download all
  the files you have requested. If any of the requested files are unavailable
  for your log in, these will be flagged up first as a message so you can
  make a note and request them through the DHS website. You also have the
  option to control whether the downloaded zip file is then extracted and
  converted into a more convenient R \code{data.frame}. This converted object
  will then be subsequently saved as a ".rds" object within the client root
  directory datasets folder, which can then be more quickly loaded when
  needed with \code{readRDS}. You also have the option to reformat the
  dataset, which will ensure that the datasets returned are encoded simply
  as character strings, i.e. there are no factors or labels.
}
\examples{
 \dontrun{
# get the model datasets included with the package
model_datasets <- model_datasets

# download one of them
g <- get_datasets(dataset_filenames = model_datasets$FileName[1])
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{search_variable_labels}
\alias{search_variable_labels}
\title{Search Survey Variable Definitions}
\usage{
search_variable_labels(
  dataset_filenames,
  search_terms = NULL,
  essential_terms = NULL,
  regex = NULL,
  ...
)
}
\arguments{
\item{dataset_filenames}{The desired filenames to be downloaded. These can be
found as one of the returned fields from \code{\link{dhs_datasets}}.}

\item{search_terms}{Character vector of search terms. If any of these terms
are found within the survey question definitions, the corresponding survey
variable and definitions will be returned.}

\item{essential_terms}{Character pattern that has to be in the definitions of
survey question definitions. I.e. the function will first find
all survey variable definitions that contain your `search_terms`
(or regex) OR `essential_terms`. It will then remove any questions
that did not contain your `essential_terms`. Default = `NULL`.}

\item{regex}{Regex character pattern for matching. If you want to specify
your regex search pattern, then specify this argument. N.B. If both
`search_terms` and `regex`` are supplied as arguments then regex will
be ignored.}

\item{...}{Any other arguments to be passed to
\code{\link{download_datasets}}}
}
\value{
A \code{data.frame} of the surveys where matches were found
  and then all the resultant codes and descriptions.
}
\description{
Searches across datasets specified for requested survey variable definitions.
This function (or \code{\link{search_variable_labels}}) should be used to
provide the `questions` argument for \code{\link{extract_dhs}}.
}
\details{
Use this function after \code{\link{get_datasets}} to query
  downloaded datasets for what survey questions they asked.
  This function will look for your downloaded and imported survey datasets
  from your cached files, and will download them if not downloaded.
}
\examples{
\dontrun{
# get the model datasets included with the package
model_datasets <- model_datasets

# download two of them
g <- get_datasets(dataset_filenames = model_datasets$FileName[1:2])

# and now seearch within these for survey variable labels of interest
vars <- search_variable_labels(
dataset_filenames = names(g), search_terms = "fever"
)

head(vars)

# if we specify an essential term then no results will be returned from
# a dataset if it does not have any results from the search with this term
search_variable_labels(
dataset_filenames = names(g),
search_terms = "fever",
essential_terms = "primaquine",
)

# we can also use regex queries if we prefer, by passing `regex = TRUE`
vars <- search_variable_labels(
dataset_filenames = names(g), search_terms = "fever|net", regex = TRUE
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{client_cache_date}
\alias{client_cache_date}
\title{Pull last cache date}
\usage{
client_cache_date(root)
}
\arguments{
\item{root}{Character for root path to where client,
caches, surveys etc. will be stored.}
}
\description{
Pull last cache date
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{response_to_json}
\alias{response_to_json}
\title{converts response to json by first converting the response to text}
\usage{
response_to_json(x)
}
\arguments{
\item{x}{A response}
}
\description{
converts response to json by first converting the response to text
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shapes.R
\name{download_boundaries}
\alias{download_boundaries}
\title{DHS Spatial Boundaries}
\usage{
download_boundaries(
  surveyNum = NULL,
  surveyId = NULL,
  countryId = NULL,
  method = "sf",
  quiet_download = FALSE,
  quiet_parse = TRUE
)
}
\arguments{
\item{surveyNum}{Numeric for the survey number to be downloaded. Values for
surveyNum can be found in the datasets or surveys endpoints in the DHS API
that can be accessed using \code{\link{dhs_datasets}} and
\code{\link{dhs_surveys}}. Default is NULL, which will cause the SurveyId
to be used to find the survey.}

\item{surveyId}{Numeric for the survey ID to be downloaded. Values for
surveyId can be found in the datasets or surveys endpoints in the DHS API
that can be accessed using \code{\link{dhs_datasets}} and
\code{\link{dhs_surveys}}. Default is NULL, which will cause the SurveyNum
to be used to find the survey.}

\item{countryId}{2-letter DHS country code for the country of the survey
being downloaded. Default = NULL, which will cause the countrycode to be
looked up from the API.}

\item{method}{Character for how the downloaded shape file is read in.
Default = "sf", which uses \code{sf::st_read}. Currenlty, you can also
specify "rgdal", which reads the file using rgdal::readOGR.
To just return the file paths for the files use method = "zip".}

\item{quiet_download}{Whether to download file quietly. Passed to
[`download_file()`]. Default is `FALSE`.}

\item{quiet_parse}{Whether to read boundaries dataset quietly. Applies to
`method = "sf"`. Default is `TRUE`.}
}
\value{
Returns either the spatial file as a "SpatialPolygonsDataFrame" or
  a vector of the file paths of where the boundary was downloaded to.
}
\description{
Download Spatial Boundaries
}
\details{
Downloads the spatial boundaries from the DHS spatial repository,
  which can be found at \url{https://spatialdata.dhsprogram.com/home/}.
}
\examples{
\dontrun{
 # using the surveyNum
 res <- download_boundaries(surveyNum = 471, countryId = "AF")

 # using the surveyId and no countryID
 res <- download_boundaries(surveyId = "AF2010OTH")

 # using rgdal
 res <- download_boundaries(surveyNum = 471, countryId = "AF", method = "rgdal")
 }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_publications}
\alias{dhs_publications}
\title{API request of DHS Publications}
\usage{
dhs_publications(
  countryIds = NULL,
  selectSurveys = NULL,
  indicatorIds = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  surveyCharacteristicIds = NULL,
  tagIds = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{selectSurveys}{Specify to filter data from the latest survey by
including `selectSurveys=TRUE` in your request. Note: Please use this
parameter in conjunction with countryCode, surveyType, or indicatorIds
for best results.}

\item{indicatorIds}{Specify a comma separated list of indicators ids to
filter by. For a list of indicators use
\code{dhs_indicators(returnFields=c("IndicatorId","Label","Definition"))}}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter
Publications on. surveyYearStart is an inclusive value. Can be used alone
or in conjunction with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Publications
on. surveyYearEnd is an inclusive value. Can be used alone or in
conjunction with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{surveyCharacteristicIds}{Specify a survey characteristic id to filter
publications with countries with the specified survey characteristics.
For a list of survey characteristics use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{tagIds}{Specify a tag id to filter publications with surveys
containing indicators with the specified tag. For a list of tags use
\code{dhs_tags()}}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format is
JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request.
By default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 10 (or less if \code{returnFields} is provided)
  publications with detailed information for each publication. A detailed
  description of all the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/publications/fields}
}
\description{
API request of DHS Publications
}
\examples{

\dontrun{
# A main use for the publications API endpoint is to find which surveys have
# a published report resulting from the conducted survey:

dat <- dhs_publications()

# A complete list of examples for how each argument to the publications
# API endpoint can be provided is given below, which is a
# copy of each of the examples listed in the API at:

# https://api.dhsprogram.com/#/api-publications.cfm


dat <- dhs_publications(countryIds="EG",all_results=FALSE)
dat <- dhs_publications(selectSurveys="latest",all_results=FALSE)
dat <- dhs_publications(indicatorIds="FE_FRTR_W_TFR",all_results=FALSE)
dat <- dhs_publications(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_publications(surveyYear="2010",all_results=FALSE)
dat <- dhs_publications(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_publications(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_publications(surveyType="DHS",all_results=FALSE)
dat <- dhs_publications(surveyCharacteristicIds="32",all_results=FALSE)
dat <- dhs_publications(tagIds=1,all_results=FALSE)
dat <- dhs_publications(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_data_updates}
\alias{dhs_data_updates}
\title{API request of DHS Data Updates}
\usage{
dhs_data_updates(
  lastUpdate = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{lastUpdate}{Specify a date or Unix time to filter the updates by.
Only results for data that have been updated on or after the specified
date will be returned.}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format is
JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request.
By default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 9 (or less if \code{returnFields} is provided)
  indicators or surveys that have been added/updated or removed. A detailed
  description of all the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/dataupdates/fields}
}
\description{
API request of DHS Data Updates
}
\examples{

\dontrun{
# The API endpoint for the data updates available within the DHS
# is a very useful endpoint, which is used a lot within `rdhs`. For example,
# we use it to keep the end user's cache up to date. For example to find all
# updates that have occurred in 2018:

dat <- dhs_data_updates(lastUpdate="20180101")

# A complete list of examples for how each argument to the data updates
# API endpoint can be provided is given below, which is a
# copy of each of the examples listed in the API at:

# https://api.dhsprogram.com/#/api-dataupdates.cfm


dat <- dhs_data_updates(lastUpdate="20150901",all_results=FALSE)
dat <- dhs_data_updates(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/authentication.R
\name{authenticate_dhs}
\alias{authenticate_dhs}
\title{DHS Website Authentication}
\usage{
authenticate_dhs(config)
}
\arguments{
\item{config}{Object of class `rdhs_config` as produced by `read_rdhs_config`
that must contain a valid `email`, `project` and `password`.}
}
\value{
Returns list of length 3:
  \itemize{
      \item{"user_name"}{ your email usually }
      \item{"user_pass"}{ your password you provided }
      \item{"proj_id"}{ your project number }
      }
}
\description{
Authenticate Users for DHS website
}
\details{
If the user has more than one project that contains the first
  30 characters of the provided project they will be prompted to choose
  which project they want. This choice will be saved so they do
  not have to enter it again in this R session.
}
\note{
Credit for some of the function to
  \url{https://github.com/ajdamico/lodown/blob/master/R/dhs.R}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{response_is_json}
\alias{response_is_json}
\title{checks if the response is json or not by looking at the responses headers}
\usage{
response_is_json(x)
}
\arguments{
\item{x}{A response}
}
\description{
checks if the response is json or not by looking at the responses headers
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_datasets.R
\name{read_dhs_dataset}
\alias{read_dhs_dataset}
\title{read in dhs standard file types}
\usage{
read_dhs_dataset(file, dataset, reformat = FALSE, all_lower = TRUE, ...)
}
\arguments{
\item{file}{path to zip file to be read}

\item{dataset}{row from \code{\link{dhs_datasets}} that
corresponds to the file}

\item{reformat}{boolean detailing if datasets should be nicely
reformatted. Default = `FALSE`}

\item{all_lower}{Logical indicating whether all value labels
should be lower case. Default to `TRUE`.}

\item{...}{Extra arguments to be passed to either
\code{\link{read_dhs_dta}} or \code{\link{read_dhs_flat}}}
}
\description{
read in dhs standard file types
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_indicators}
\alias{dhs_indicators}
\title{API request of DHS Indicators}
\usage{
dhs_indicators(
  countryIds = NULL,
  indicatorIds = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  surveyCharacteristicIds = NULL,
  tagIds = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
  by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{indicatorIds}{Specify a comma separated list of indicators ids to
  filter by. For a list of indicators use
\code{dhs_indicators(returnFields=c("IndicatorId","Label","Definition"))}}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
  For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
  "SurveyType","CountryName"))}}

\item{surveyYear}{Specify a survey year to filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Indicators
on. surveyYearStart is an inclusive value. Can be used alone or in
conjunction with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Indicators
on. surveyYearEnd is an inclusive value. Can be used alone or in
conjunction with surveyYearStart.}

\item{surveyType}{Specify a comma separated list of survey years to
filter by.}

\item{surveyCharacteristicIds}{Specify a survey characteristic id to filter
indicators in surveys with the specified survey characteristic. For a list
of survey characteristics use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{tagIds}{Specify a tag id to filter indicators with the specified tag.
For a list of tags use \code{dhs_tags()}}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request. By
default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 18 (or less if \code{returnFields} is provided)
  indicators with attributes for each indicator. A detailed description of
  all the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/indicators/fields}
}
\description{
API request of DHS Indicators
}
\examples{
\dontrun{
# A common use for the indicators data API will be to search for a list of
# health indicators within a given characteristic category, such as anemia
# testing, HIV prevalence, micronutrients etc. For example to return all the
# indicators relating to malaria testing by RDTs:

dat <- dhs_indicators(surveyCharacteristicIds="90")

# A list of the different `surveyCharacteristicIds` can be found
# [here](https://api.dhsprogram.com/rest/dhs/surveycharacteristics?f=html)

# A complete list of examples for how each argument to the indicator API
# endpoint can be provided is given below, which is a copy of each of
# the examples listed in the API at:

# https://api.dhsprogram.com/#/api-indicators.cfm


dat <- dhs_indicators(countryIds="EG",all_results=FALSE)
dat <- dhs_indicators(indicatorIds="FE_FRTR_W_TFR",all_results=FALSE)
dat <- dhs_indicators(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_indicators(surveyYear="2010",all_results=FALSE)
dat <- dhs_indicators(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_indicators(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_indicators(surveyType="DHS",all_results=FALSE)
dat <- dhs_indicators(surveyCharacteristicIds="32",all_results=FALSE)
dat <- dhs_indicators(tagIds="1",all_results=FALSE)
dat <- dhs_indicators(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/authentication.R
\name{available_datasets}
\alias{available_datasets}
\title{Create a data frame of datasets that your log in can download}
\usage{
available_datasets(
  config,
  datasets_api_results = NULL,
  surveys_api_results = NULL
)
}
\arguments{
\item{config}{Object of class `rdhs_config` as produced by `read_rdhs_config`
that must contain a valid `email`, `project` and `password`.}

\item{datasets_api_results}{Data.table for the api results for the datasets
  endpoint. Default = NULL and
generated by default if not declared.}

\item{surveys_api_results}{Data.table for the api results for the surveys
  endpoint. Default = NULL and
generated by default if not declared.}
}
\value{
Returns \code{"data.frame"} of length 14:
\itemize{
      \item{"FileFormat"}
      \item{"FileSize"}
      \item{"DatasetType"}
      \item{"SurveyNum"}
      \item{"SurveyId"}
      \item{"FileType"}
      \item{"FileDateLastModified"}
      \item{"SurveyYearLabel"}
      \item{"SurveyType"}
      \item{"SurveyYear"}
      \item{"DHS_CountryCode"}
      \item{"FileName"}
      \item{"CountryName"}
      \item{"URLS"}
      }
}
\description{
DHS datasets that can be downloaded
}
\note{
Inspiration for function to
  \url{https://github.com/ajdamico/lodown/blob/master/R/dhs.R}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{model_datasets}
\alias{model_datasets}
\title{DHS model datasets}
\format{
A dataframe of 36 observations of 14 variables:

\code{model_datasets}: A dataframe of model datasets
\itemize{
      \item{"FileFormat"}
      \item{"FileSize"}
      \item{"DatasetType"}
      \item{"SurveyNum"}
      \item{"SurveyId"}
      \item{"FileType"}
      \item{"FileDateLastModified"}
      \item{"SurveyYearLabel"}
      \item{"SurveyType"}
      \item{"SurveyYear"}
      \item{"DHS_CountryCode"}
      \item{"FileName"}
      \item{"CountryName"}
      \item{"URLS"}
      }
}
\usage{
data(model_datasets)
}
\description{
The model datasets from the DHS website in a `data.frame` that is analogous
to those returned by `get_available_datasets()`
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdhs.R
\docType{package}
\name{rdhs}
\alias{rdhs}
\alias{rdhs-package}
\title{\pkg{rdhs} DHS database through R}
\description{
Provides a client for (1) querying the DHS API for survey indicators
and metadata, (2) identifying surveys and datasets for analysis,
(3) downloading survey datasets from the DHS website, (4) loading
datasets and associate metadata into R, and (5) extracting variables
and combining datasets for pooled analysis.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/rdhs/}
  \item Report bugs at \url{https://github.com/ropensci/rdhs/issues}
}

}
\author{
\strong{Maintainer}: OJ Watson \email{oj.watson@hotmail.co.uk} (\href{https://orcid.org/0000-0003-2374-0741}{ORCID})

Authors:
\itemize{
  \item Jeff Eaton (\href{https://orcid.org/0000-0001-7728-728X}{ORCID})
}

Other contributors:
\itemize{
  \item Lucy D'Agostino McGowan (\href{https://orcid.org/0000-0001-7297-9359}{ORCID}) [reviewer]
  \item Duncan Gillespie [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dhs_gps_data_format}
\alias{dhs_gps_data_format}
\title{DHS GPS Data Format}
\format{
A dataframe of 20 observations of 3 variables:

\code{dhs_gps_data_format}: A dataframe of GPS data descriptions.
\itemize{
      \item{"Name"}
      \item{"Type"}
      \item{"Description"}
      }
}
\usage{
data(dhs_gps_data_format)
}
\description{
Data frame to describe the data encoded in DHS GPS files
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_ui_updates}
\alias{dhs_ui_updates}
\title{API request of DHS UI Updates}
\usage{
dhs_ui_updates(
  lastUpdate = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{lastUpdate}{Specify a date or Unix time to filter the updates by. Only
results for interfaces that has been updated on or after the specified
date will be returned.}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request. By
default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 3 (or less if \code{returnFields} is provided)
  interfaces that have been added/updated or removed. A detailed description
  of all the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/uiupdates/fields}
}
\description{
API request of DHS UI Updates
}
\examples{

\dontrun{
# The main use for the ui updates API will be to search for the last time
# there was a change to the UI. For example to return all the
# changes since 2018:

dat <- dhs_ui_updates(lastUpdate="20180101")

# A complete list of examples for how each argument to the ui updates API
# endpoint can be provided is given below, which is a copy of each of
# the examples listed in the API at:

# https://api.dhsprogram.com/#/api-uiupdates.cfm


dat <- dhs_ui_updates(lastUpdate="20150901",all_results=FALSE)
dat <- dhs_ui_updates(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_info}
\alias{dhs_info}
\title{API request of DHS Info}
\usage{
dhs_info(
  infoType = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{infoType}{Specify a type of info to obtain the information requested.
Default is version. `infoType="version"`` (default) Provides the version
of the API.
Example: https://api.dhsprogram.com/rest/dhs/info?infoType=version
`infoType="citation"` Provides the citation for the API to include with
your application or data.
Example: https://api.dhsprogram.com/rest/dhs/info?infoType=citation}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request. By
default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 2 (or less if \code{returnFields} is provided)
  fields describing the type of information that was requested and a value
  corresponding to the information requested.
  \url{https://api.dhsprogram.com/rest/dhs/info/fields}
}
\description{
API request of DHS Info
}
\examples{

\dontrun{
# The main use for the info API  will be to confirm the version of the API
# being used to providing the most current citation for the data.

dat <- dhs_info(infoType="version")

# A complete list of examples for how each argument to the info API
# endpoint can be provided is given below, which is a copy of each of
# the examples listed in the API at:

# https://api.dhsprogram.com/#/api-info.cfm


dat <- dhs_info(infoType="version",all_results=FALSE)
dat <- dhs_info(infoType="citation",all_results=FALSE)
dat <- dhs_info(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extraction.R
\name{extraction}
\alias{extraction}
\title{DHS survey questions extracted from datasets}
\usage{
extraction(questions, available_datasets, geo_surveys, add_geo = FALSE)
}
\arguments{
\item{questions}{Output of
\code{R6_client_dhs$public_methods$survey_questions}}

\item{available_datasets}{Datasets that could be available.
Output of \code{R6_client_dhs$public_methods$available_datasets}}

\item{geo_surveys}{Geographic Data Survey file paths.}

\item{add_geo}{Boolean detailing if geographic datasets should be added.}
}
\value{
Returns `data.frame` with variables corresponding to
  the requested variables in the questions object. Will also have
  geographic data related columns if `add_geo=TRUE` is set.
  Lastly a SurveyId variable will also be appended corresponding to
  \code{\link{dhs_datasets}}$SurveyId
}
\description{
Create a list of survey responses extracted using output of
\code{R6_client_dhs$public_methods$survey_questions}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{get_available_datasets}
\alias{get_available_datasets}
\title{Get Available Datasets}
\usage{
get_available_datasets(clear_cache = FALSE)
}
\arguments{
\item{clear_cache}{Boolean detailing if you would like to clear the
cached available datasets first. The default is set to FALSE. This option
is available so that you can make sure your client fetches any new datasets
that you have recently been given access to.}
}
\value{
A \code{data.frame} with 14 variables that detail the surveys you can
  download, their url download links and the country, survey, year etc info
  for that link.
}
\description{
Details the datasets that your login credentials have access to
}
\details{
Searches the DHS website for all the datasets that you can download.
  The results of this function are cached in the client. If you have recently
  requested new datasets from the DHS website then you can specify to clear
  the cache first so that you get the new set of datasets available to you.
  This function is used by \code{\link{get_datasets}} and should thus be
  used with `clear_cache_first = TRUE` before using `get_datasets` if you
  have recently requested new datasets.
}
\examples{

\dontrun{
# grab the datasets
datasets <- get_available_datasets()

# and if we look at the last one it will be the model datasets from DHS
tail(datasets, 1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_datasets}
\alias{dhs_datasets}
\title{API request of DHS Datasets}
\usage{
dhs_datasets(
  countryIds = NULL,
  selectSurveys = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  fileFormat = NULL,
  fileType = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{selectSurveys}{Specify to filter data from the latest survey by
including `selectSurveys=TRUE` in your request. Note: Please use this
parameter in conjunction with countryCode, surveyType, or indicatorIds for
best results.}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Datasets
on. surveyYearStart is an inclusive value. Can be used alone or in
conjunction with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Datasets on.
surveyYearEnd is an inclusive value. Can be used alone or in conjunction
with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{fileFormat}{Specify the file format for the survey. Can use file
format type name (SAS, Stata, SPSS, Flat, Hierarchical) or file format
code. View list of file format codes -
\url{https://dhsprogram.com/data/File-Types-and-Names.cfm}}

\item{fileType}{Specify the type of dataset generated for the survey
(e.g. household, women, men, children, couples, etc.). View list of file
type codes - \url{https://dhsprogram.com/data/File-Types-and-Names.cfm}}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format is
JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request.
By default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 13 (or less if \code{returnFields} is provided)
  datasets with their corresponding details. A detailed description of all
  the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/datasets/fields}
}
\description{
API request of DHS Datasets
}
\examples{

\dontrun{
# The API endpoint for the datasets available within the DHS website
# is a very useful endpoint, which is used a lot within `rdhs`. For example,
# it is used to find the file names and size of the dataset files, as well
# as when they were last modified. This enables us to see which datasets
# have been updated and may thus be out of date. For example to find all
# datasets that have been modified in 2018:

dat <- dhs_datasets()
dates <- rdhs:::mdy_hms(dat$FileDateLastModified)
years <- as.POSIXlt(dates, tz = tz(dates))$year + 1900
modified_in_2018 <- which(years == 2018)

# A complete list of examples for how each argument to the datasets
# API endpoint can be provided is given below, which is a
# copy of each of the examples listed in the API at:

# https://api.dhsprogram.com/#/api-datasets.cfm


dat <- dhs_datasets(countryIds="EG",all_results=FALSE)
dat <- dhs_datasets(selectSurveys="latest",all_results=FALSE)
dat <- dhs_datasets(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_datasets(surveyYear="2010",all_results=FALSE)
dat <- dhs_datasets(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_datasets(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_datasets(surveyType="DHS",all_results=FALSE)
dat <- dhs_datasets(fileFormat="stata",all_results=FALSE)
dat <- dhs_datasets(fileFormat="DT",all_results=FALSE)
dat <- dhs_datasets(fileType="KR",all_results=FALSE)
dat <- dhs_datasets(f="geojson",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/authentication.R
\name{download_datasets}
\alias{download_datasets}
\title{Create a data frame of datasets that your log in can download}
\usage{
download_datasets(
  config,
  desired_dataset,
  download_option = "both",
  reformat = TRUE,
  all_lower = TRUE,
  output_dir_root = NULL,
  ...
)
}
\arguments{
\item{config}{Object of class `rdhs_config` as produced by `read_rdhs_config`
that must contain a valid `email`, `project` and `password`.}

\item{desired_dataset}{Row from \code{available_datasets}}

\item{download_option}{Character dictating how the survey is stored when
downloaded. Must be one of:
\itemize{
    \item{"zip"} - Just the zip. "z", "i", "p" or "zip" will match
    \item{"rds"} - Just the read in and saved rds. "r", "d", "s" or "rdhs"
    will match
    \item{"both"} - Both the rds and extract. "b", "o", "t", "h" or "both"
    will match
    }}

\item{reformat}{Boolean detailing whether dataset rds should be
reformatted for ease of use later. Default = TRUE}

\item{all_lower}{Logical indicating whether all value labels should be
lower case. Default to `TRUE`.}

\item{output_dir_root}{Directory where files are to be downloaded to}

\item{...}{Any other arguments to be passed to
\code{\link{read_dhs_dataset}}}
}
\description{
Download datasets specified using output of \code{available_datasets}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_dhs_flat.R
\name{read_dhs_flat}
\alias{read_dhs_flat}
\title{Read DHS flat file data set}
\usage{
read_dhs_flat(zfile, all_lower = TRUE, meta_source = NULL)
}
\arguments{
\item{zfile}{Path to `.zip` file containing flat file dataset, usually
ending in filename `XXXXXXFL.zip`}

\item{all_lower}{Logical indicating whether all value labels should be
lower case. Default to `TRUE`.}

\item{meta_source}{character string indicating metadata source file for data
dictionary. Default \code{NULL} first tried to use \code{.DCF} and then
{.SPS} if not found.}
}
\value{
A data frame. Value labels for each variable are stored as the
  `labelled` class from `haven`.
}
\description{
This function reads a DHS recode dataset from the zipped flat file dataset.
}
\examples{
mrfl_zip <- tempfile()
download.file("https://dhsprogram.com/data/model_data/dhs/zzmr61fl.zip",
              mrfl_zip,mode="wb")

mr <- rdhs:::read_dhs_flat(mrfl_zip)
attr(mr$mv213, "label")
class(mr$mv213)
head(mr$mv213)
table(mr$mv213)
table(haven::as_factor(mr$mv213))

}
\seealso{
\code{\link[haven]{labelled}}, \code{\link{read_dhs_dta}}.

For more information on the DHS filetypes and contents of distributed
dataset .ZIP files, see
\url{https://dhsprogram.com/data/File-Types-and-Names.cfm#CP_JUMP_10334}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_data}
\alias{dhs_data}
\title{API request of DHS Indicator Data}
\usage{
dhs_data(
  countryIds = NULL,
  indicatorIds = NULL,
  surveyIds = NULL,
  selectSurveys = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  surveyCharacteristicIds = NULL,
  characteristicCategory = NULL,
  characteristicLabel = NULL,
  tagIds = NULL,
  breakdown = NULL,
  returnGeometry = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
  by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{indicatorIds}{Specify a comma separated list of indicator ids to
  filter by. For a list of indicators use
\code{dhs_indicators(returnFields=c("IndicatorId","Label","Definition"))}}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
  For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
  "SurveyType","CountryName"))}}

\item{selectSurveys}{Specify to filter Data from the latest survey by adding
`selectSurveys="latest"` in conjunction with a Country Code and/or Survey
Type. Please Note: Not all indicators are present in the latest surveys.
To filter your API Indicator Data call to return the latest survey data in
which a specific set of indicators is present, add
`selectSurveys="byIndicator"` in conjunction with Indicator IDs, Country
Code, and/or Survey Type instead
of using `selectSurveys="latest"`}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Data on.
surveyYearStart is an inclusive value. Can be used alone or in conjunction
with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Data on.
surveyYearEnd is an inclusive value. Can be used alone or in conjunction
with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{surveyCharacteristicIds}{Specify a survey characteristic id to filter
  data on surveys with the specified survey characteristic. For a list of
  survey characteristics use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
  "SurveyType","CountryName"))}}

\item{characteristicCategory}{Specify a survey characteristic category to
filter data on surveys with the specified survey characteristic category.
This query is case insensitive, but it only recognizes exact phrase
matches. For example, `characteristicCategory="wealth"` will return
results that have a characteristic category of `Wealth` while
`characteristicCategory="wealth quintile"' will return results that have
a characteristic category of `Wealth Quintile`.}

\item{characteristicLabel}{Specify a survey characteristic category to
filter data on surveys with the specified survey characteristic category.
This query is case insensitive, but it only recognizes exact phrase
matches. You can use characteristicLabel on its own or in conjunction with
characteristicCategory.}

\item{tagIds}{Specify a tag id to filter data on indicators with the
specified tag. For a list of tags use \code{dhs_tags()}}

\item{breakdown}{Data can be requested at different levels via the breakdown
parameter. By default national data is returned and provides totals on a
national level. `breakdown="subnational"` data provides values on a
subnational level. `breakdown="background"` provides totals on categorized
basis. Examples are urban/rural, education and wealth level.
`breakdown="all"` provides all the data including disaggregated data.}

\item{returnGeometry}{Coordinates can be requested from the API by including
`returnGeometry=TRUE` in your request. The default for this value is
false.}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format is
JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request.
By default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE}
}
\value{
Returns a \code{data.table} of 27 (or less if \code{returnFields} is provided)
  data for your particular query. Details of properties returned with each
  row of data are provided at
\url{https://api.dhsprogram.com/rest/dhs/data/fields}
}
\description{
API request of DHS Indicator Data
}
\examples{

\dontrun{
# A common use for the indicator data API will be to search for a specific
# health indicator for a given country. For example to return the total
# malaria prevalence according to RDT, given by the indicator ML_PMAL_C_RDT,
# in Senegal since 2010:

dat <- dhs_data(
indicatorIds="ML_PMAL_C_RDT",
countryIds="SN",
surveyYearStart="2006"
)

# A complete list of examples for how each argument to the data api
# endpoint can be provided is given below, which is a copy of each of
# the examples listed in the API at:

# https://api.dhsprogram.com/#/api-data.cfm


dat <- dhs_data(countryIds="EG",all_results=FALSE)
dat <- dhs_data(indicatorIds="FE_FRTR_W_TFR",all_results=FALSE)
dat <- dhs_data(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_data(selectSurveys="latest",all_results=FALSE)
dat <- dhs_data(selectSurveys="byIndicator", indicatorIds="FE_CEBA_W_CH0",
all_results=FALSE)
dat <- dhs_data(surveyYear="2010",all_results=FALSE)
dat <- dhs_data(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_data(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_data(surveyType="DHS",all_results=FALSE)
dat <- dhs_data(surveyCharacteristicIds="32",all_results=FALSE)
dat <- dhs_data(characteristicCategory="wealth quintile",all_results=FALSE)
dat <- dhs_data(breakdown="all", countryIds="AZ", characteristicLabel="6+",
all_results=FALSE)
dat <- dhs_data(tagIds="1",all_results=FALSE)
dat <- dhs_data(breakdown="subnational",all_results=FALSE)
dat <- dhs_data(breakdown="background",all_results=FALSE)
dat <- dhs_data(breakdown="all",all_results=FALSE)
dat <- dhs_data(f="html",all_results=FALSE)
dat <- dhs_data(f="geojson", returnGeometry="true",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_datasets.R
\name{read_zipdata}
\alias{read_zipdata}
\title{Read filetype from a zipped folder based on the file ending}
\usage{
read_zipdata(zfile, pattern = ".dta$", readfn = haven::read_dta, ...)
}
\arguments{
\item{zfile}{Path to `.zip` file containing flat file dataset,
usually ending in filename `XXXXXXFL.zip`}

\item{pattern}{String detailing which filetype is to be read
from within the zip by means of a grep. Default = ".dta$"}

\item{readfn}{Function object to be used for reading in the
identified file within the zip. Default = `haven::read_dta`}

\item{...}{additional arguments to readfn}
}
\description{
Read filetype from a zipped folder based on the file ending
}
\examples{
\dontrun{
# get the model datasets included in the package
model_datasets <- model_datasets

# download just the zip
g <- get_datasets(
dataset_filenames = model_datasets$FileName[1],
download_option = "zip"
)

# and then read from the zip. This function is used internally by rdhs
# when using `get_datasets` with `download_option = .rds` (default)
r <- read_zipdata(
g[[1]], pattern = ".dta"
)

# and we can pass a function to read the file and any other args with ...
r <- read_zipdata(
g[[1]], pattern = ".dta", readfn = haven::read_dta, encoding = "UTF-8"
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_datasets.R
\name{data_and_labels}
\alias{data_and_labels}
\title{Create list of dataset and its variable names}
\usage{
data_and_labels(dataset)
}
\arguments{
\item{dataset}{Any read in dataset created by \code{get_datasets},
either as the file path or after having
been read using \code{readRDS}}
}
\description{
Function to give the former output of get_datasets as it can be nice to have
both the definitions and the dataset attached together
}
\examples{
\dontrun{
# get the model datasets included with the package
model_datasets <- model_datasets

# download one of them
g <- get_datasets(dataset_filenames = model_datasets$FileName[1])
dl <- data_and_labels(g$zzbr62dt)

# now we easily have our survey question labels easily accessible
grep("bed net", dl$variable_names$description, value = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{file_dataset_format}
\alias{file_dataset_format}
\title{Returns what the dataset file ending should be for a given filename}
\usage{
file_dataset_format(file_format)
}
\arguments{
\item{file_format}{FileFormat for a file as taken from the API,
e.g. \code{dhs_datasets(returnFields = "FileFormat")}}
}
\value{
One of "dat","dat","sas7bdat","sav" or "dta"
}
\description{
Returns what the dataset file ending should be for a given filename
}
\examples{
file_format <- "Stata dataset (.dta)"
identical(rdhs:::file_dataset_format(file_format),"dta")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_geometry}
\alias{dhs_geometry}
\title{API request of DHS Geometry}
\usage{
dhs_geometry(
  countryIds = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Geometry
on. surveyYearStart is an inclusive value. Can be used alone or in
conjunction with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Geometry on.
surveyYearEnd is an inclusive value. Can be used alone or in conjunction
with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request.
By default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 7 (or less if \code{returnFields} is provided)
  geometry with their corresponding details. A detailed description of all
  the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/geometry/fields}
}
\description{
API request of DHS Geometry
}
\examples{

\dontrun{
# The geometry API endpoint returns the spatial geometry for countries, and
# can then be used to recreate the spatial polygon for a given country. For
# example to return the coordinates for the Senegal 2010 DHS survey:

dat <- dhs_geometry(surveyIds="SN2010DHS")

# At the moment there is no function to convert the coordinates returned by
# the API but this will be in the next version of rdhs. For those interested
# look at the geojson vignette for an alternative way to produce plots.

# A complete list of examples for how each argument to the geometry
# API endpoint can be provided is given below, which is a
# copy of each of the examples listed in the API at:

# https://api.dhsprogram.com/#/api-geometry.cfm


dat <- dhs_geometry(countryIds="EG",all_results=FALSE)
dat <- dhs_geometry(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_geometry(surveyYear="2010",all_results=FALSE)
dat <- dhs_geometry(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_geometry(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_geometry(surveyType="DHS",all_results=FALSE)
dat <- dhs_geometry(f="geojson",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{search_variables}
\alias{search_variables}
\title{Search Survey Variables}
\usage{
search_variables(dataset_filenames, variables, essential_variables = NULL, ...)
}
\arguments{
\item{dataset_filenames}{The desired filenames to be downloaded.
These can be found as one of the returned fields from
\code{\link{dhs_datasets}}.}

\item{variables}{Character vector of survey variables to be looked up}

\item{essential_variables}{Character vector of variables that need to
present. If any of the codes are not present in that survey,
the survey will not be returned by this function. Default = `NULL`.}

\item{...}{Any other arguments to be passed to
\code{\link{download_datasets}}}
}
\value{
A \code{data.frame} of the surveys where matches were
found and then all the resultant codes and descriptions.
}
\description{
Searches across datasets specified for requested survey variables.
This function (or \code{\link{search_variable_labels}})
should be used to provide the `questions` argument
for \code{\link{extract_dhs}}.
}
\details{
Use this function after \code{\link{get_datasets}} to look up all
  the survey variables that have the required variable.
}
\examples{
\dontrun{
# get the model datasets included with the package
model_datasets <- model_datasets

# download two of them
g <- get_datasets(dataset_filenames = model_datasets$FileName[1:2])

# and now seearch within these for survey variables
search_variables(
dataset_filenames = names(g), variables = c("v002","v102","ml13"),
)

# if we specify an essential variable then that dataset has to have that
# variable or else no variables will be returned for that datasets
search_variables(
dataset_filenames = names(g),
variables = c("v002","v102","ml13"),
essential_variables = "ml13"
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{get_rdhs_config}
\alias{get_rdhs_config}
\title{Get rdhs config}
\usage{
get_rdhs_config()
}
\value{
A \code{data.frame} containing your rdhs config
}
\description{
Gets the rdhs config being used
}
\details{
Returns the config being used by rdhs at the moment. This will
  either be a `data.frame` with class `rdhs_config` or will be NULL if
  this has not been set up yet
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{\link[magrittr:pipe]{\%>\%}} for more details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{collapse_api_responses}
\alias{collapse_api_responses}
\title{collapse API response list}
\usage{
collapse_api_responses(x)
}
\arguments{
\item{x}{List of lists from API to be collapsed}
}
\description{
collapse API response list
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rbind_list_base}
\alias{rbind_list_base}
\title{implementation of data.tables rbindlist}
\usage{
rbind_list_base(x)
}
\arguments{
\item{x}{List of lists to be converted to a data.frame}
}
\description{
implementation of data.tables rbindlist
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbind_labelled.R
\name{rbind_labelled}
\alias{rbind_labelled}
\title{Combine data frames with columns of class `labelled`}
\usage{
rbind_labelled(..., labels = NULL, warn = TRUE)
}
\arguments{
\item{...}{data frames to bind together, potentially with columns of
class "labelled". The first argument can be a list of data frames, similar
to `plyr::rbind.fill`.}

\item{labels}{A named list providing vectors of value labels or describing
how to handle columns of class `labelled`. See details for usage.}

\item{warn}{Logical indicating to warn if combining variables with different
value labels. Defaults to TRUE.}
}
\value{
A data frame.
}
\description{
Combine data frames with columns of class `labelled`
}
\details{
The argument `labels` provides options for how to handle binding of
columns of class `labelled`. Typical use is to provide a named list
with elements for each labelled column. Elements of the list
are either a vector of labels that should be applied to the column
or the character string "concatenated", which indicates that labels
should be concatenated such that all unique labels are distinct
values in the combined vector. This is accomplished by converting
to character strings, binding, and then casting back to labelled.
For labelled columns for which labels are not provided in the `label`
argument, the default behaviour is that the labels from the first
data frame with labels for that column are inherited by the combined
data.

See examples.
}
\examples{
df1 <- data.frame(
area = haven::labelled(c(1L, 2L, 3L), c("reg 1"=1,"reg 2"=2,"reg 3"=3)),
climate = haven::labelled(c(0L, 1L, 1L), c("cold"=0,"hot"=1))
)
df2 <- data.frame(
area    = haven::labelled(c(1L, 2L), c("reg A"=1, "reg B"=2)),
climate = haven::labelled(c(1L, 0L), c("cold"=0, "warm"=1))
)

# Default: all data frames inherit labels from first df. Incorrect if
# "reg 1" and "reg A" are from different countries, for example.
dfA <- rbind_labelled(df1, df2)
haven::as_factor(dfA)

# Concatenate value labels for "area". Regions are coded separately,
# and original integer values are lost (by necessity of more levels now).
# For "climate", codes "1 = hot" and "1 = warm", are coded as the same
# outcome, inheriting "1 = hot" from df1 by default.
dfB <- rbind_labelled(df1, df2, labels=list(area = "concatenate"))
dfB
haven::as_factor(dfB)

# We can specify to code as "1=warm/hot" rather than inheriting "hot".
dfC <- rbind_labelled(df1, df2,
labels=list(area = "concatenate", climate = c("cold"=0, "warm/hot"=1)))

dfC$climate
haven::as_factor(dfC)

# Or use `climate="concatenate"` to code "warm" and "hot" as different.
dfD <- rbind_labelled(df1, df2,
labels=list(area = "concatenate", climate="concatenate"))

dfD
haven::as_factor(dfD)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbind_labelled.R
\name{delabel_df}
\alias{delabel_df}
\title{convert labelled data frame to data frame of just characters}
\usage{
delabel_df(df)
}
\arguments{
\item{df}{data frame to convert labelled elements of. Likely this will be
the output of \code{\link{extract_dhs}}.}
}
\value{
A data frame of de-labelled elements
}
\description{
convert labelled data frame to data frame of just characters
}
\examples{
df1 <- data.frame(
area = haven::labelled(c(1L, 2L, 3L), c("reg 1"=1,"reg 2"=2,"reg 3"=3)),
climate = haven::labelled(c(0L, 1L, 1L), c("cold"=0,"hot"=1))
)

df_char <- delabel_df(df = df1)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{update_rdhs_config}
\alias{update_rdhs_config}
\title{Update your current rdhs config}
\usage{
update_rdhs_config(
  password = FALSE,
  email = NULL,
  project = NULL,
  cache_path = NULL,
  config_path = NULL,
  global = NULL,
  verbose_download = NULL,
  verbose_setup = NULL,
  timeout = NULL,
  data_frame = NULL,
  project_choice = NULL
)
}
\arguments{
\item{password}{Logical for updating your password securely. Default = FALSE}

\item{email}{Character for email used to login to the DHS website.}

\item{project}{Character for the name of the DHS project from which
datasets should be downloaded.}

\item{cache_path}{Character for directory path where datasets and API calls
will be cached. If left bank, a suitable directory will be created within
your user cache directory for your operating system (permission granting).}

\item{config_path}{Character for where the config file should be saved.
For a global configuration, `config_path` must be '~/.rdhs.json'.
For a local configuration, `config_path` must be 'rdhs.json'.
If left bank, the config file will be stored within
your user cache directory for your operating system (permission granting).}

\item{global}{Logical for the config_path to be interpreted as a global
config path or a local one. Default = TRUE.}

\item{verbose_download}{Logical for dataset download progress bars to be
shown. Default = FALSE.}

\item{verbose_setup}{Logical for rdhs setup and messages to be printed.
Default = TRUE.}

\item{timeout}{Numeric for how long in seconds to wait for the DHS API to
respond. Default = 30.}

\item{data_frame}{Function with which to convert API calls into. If left
blank \code{data_frame} objects are returned. Must be passed as a
character. Examples could be:
\code{data.table::as.data.table}
\code{tibble::as.tibble}}

\item{project_choice}{Numeric for project choice. See \code{authenticate_dhs}
for more info.}
}
\description{
\code{update_rdhs_config} allows you to update elements of your
rdhs config, without having to set it completely via \code{set_rdhs_config}.
For each config element, provide the new changes required. To update your
password, set \code{password = TRUE} and you will be asked securely for your
new password.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_datasets.R
\name{get_labels_from_dataset}
\alias{get_labels_from_dataset}
\title{Return variable labels from a dataset}
\usage{
get_labels_from_dataset(data, return_all = TRUE)
}
\arguments{
\item{data}{A \code{data.frame} from which to extract variable labels.}

\item{return_all}{Logical whether to return all variables (\code{TRUE})
or only those with labels.}
}
\value{
A \code{data.frame} consisting of the variable name and labels.
}
\description{
Returns variable labels stored as \code{"label"} attribute.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_surveys}
\alias{dhs_surveys}
\title{API request of DHS Surveys}
\usage{
dhs_surveys(
  countryIds = NULL,
  indicatorIds = NULL,
  selectSurveys = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  surveyStatus = NULL,
  surveyCharacteristicIds = NULL,
  tagIds = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to
filter by. For a list of countries use
 \code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{indicatorIds}{Specify a comma separated list of indicators ids to
filter by. For a list of indicators use
\code{dhs_indicators(returnFields=c("IndicatorId","Label","Definition"))}}

\item{selectSurveys}{Specify to filter data from the latest survey by
including `selectSurveys=TRUE` in your request. Note: Please use this
parameter in conjunction with countryCode, surveyType, or indicatorIds
for best results.}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Surveys on.
surveyYearStart is an inclusive value. Can be used alone or in conjunction
with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Surveys on.
surveyYearEnd is an inclusive value. Can be used alone or in conjunction
with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{surveyStatus}{Every survey is assigned a surveys status and can be
queried based on the surveyStatus parameter. `surveyStatus="available"`
(default) provides a list of all surveys for which the DHS API contains
Indicator Data. `surveyStatus="Completed"` provides a list of all
completed surveys. NOTE: Data may not be available for every completed
survey. `surveyStatus="Ongoing"` provides a list of all ongoing surveys.
`surveyStatus="all"` provides a list of all surveys.}

\item{surveyCharacteristicIds}{Specify a survey characteristic id to filter
surveys with the specified survey characteristic. For a list of survey
characteristics use \code{dhs_surveys(returnFields=c("SurveyId",
"SurveyYearLabel","SurveyType","CountryName"))}}

\item{tagIds}{Specify a tag id to filter surveys containing indicators with
the specified tag. For a list of tags use \code{dhs_tags()}}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request. By
default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 28 (or less if \code{returnFields} is provided)
  surveys with detailed information for each survey. A detailed description
  of all the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/surveys/fields}
}
\description{
API request of DHS Surveys
}
\examples{

\dontrun{
# A common use for the surveys API endpoint is to query which countries
# have conducted surveys since a given year, e.g. since 2010

dat <- dhs_surveys(surveyYearStart="2010")

# Additionally, some countries conduct non DHS surveys, but the data for
# thse is also available within the DHS website/API. To query these:

dat <- dhs_surveys(surveyType="MIS")

# Lastly, you may be interested to know about anything peculiar about a
# particular survey's implementation. This can be found by looking within
# the footnotes variable within the data frame returned. For example, the
# Madagascar 2013 MIS:

dat$Footnotes[dat$SurveyId == "MD2013MIS"]

# A complete list of examples for how each argument to the surveys API
# endpoint can be provided is given below, which is a copy of each of
# the examples listed in the API at:

# https://api.dhsprogram.com/#/api-surveys.cfm


dat <- dhs_surveys(countryIds="EG",all_results=FALSE)
dat <- dhs_surveys(indicatorIds="FE_FRTR_W_TFR",all_results=FALSE)
dat <- dhs_surveys(selectSurveys="latest",all_results=FALSE)
dat <- dhs_surveys(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_surveys(surveyYear="2010",all_results=FALSE)
dat <- dhs_surveys(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_surveys(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_surveys(surveyType="DHS",all_results=FALSE)
dat <- dhs_surveys(surveyStatus="Surveys",all_results=FALSE)
dat <- dhs_surveys(surveyStatus="Completed",all_results=FALSE)
dat <- dhs_surveys(surveyStatus="Ongoing",all_results=FALSE)
dat <- dhs_surveys(surveyStatus="All",all_results=FALSE)
dat <- dhs_surveys(surveyCharacteristicIds="32",all_results=FALSE)
dat <- dhs_surveys(tagIds="1",all_results=FALSE)
dat <- dhs_surveys(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{set_rdhs_config}
\alias{set_rdhs_config}
\title{Set rdhs config}
\usage{
set_rdhs_config(
  email = NULL,
  project = NULL,
  cache_path = NULL,
  config_path = NULL,
  global = TRUE,
  verbose_download = FALSE,
  verbose_setup = TRUE,
  data_frame = NULL,
  timeout = 30,
  password_prompt = FALSE,
  prompt = TRUE
)
}
\arguments{
\item{email}{Character for email used to login to the DHS website.}

\item{project}{Character for the name of the DHS project from which
datasets should be downloaded.}

\item{cache_path}{Character for directory path where datasets and API calls
will be cached. If left bank, a suitable directory will be created within
your user cache directory for your operating system (permission granting).}

\item{config_path}{Character for where the config file should be saved.
For a global configuration, `config_path` must be '~/.rdhs.json'.
For a local configuration, `config_path` must be 'rdhs.json'.
If left bank, the config file will be stored within
your user cache directory for your operating system (permission granting).}

\item{global}{Logical for the config_path to be interpreted as a global
config path or a local one. Default = TRUE.}

\item{verbose_download}{Logical for dataset download progress bars to be
shown. Default = FALSE.}

\item{verbose_setup}{Logical for rdhs setup and messages to be printed.
Default = TRUE.}

\item{data_frame}{Function with which to convert API calls into. If left
blank \code{data_frame} objects are returned. Must be passed as a
character. Examples could be:
\code{data.table::as.data.table}
\code{tibble::as.tibble}}

\item{timeout}{Numeric for how long in seconds to wait for the DHS API to
respond. Default = 30.}

\item{password_prompt}{Logical whether user is asked to type their password,
even if they have previously set it. Default = FALSE. Set to TRUE if you
have mistyped your password when using \code{set_rdhs_config}.}

\item{prompt}{Logical for whether the user should be prompted for
permission to write to files. This should not need be
changed by the user. Default = TRUE.}
}
\value{
Invisibly returns the rdhs config object
}
\description{
Sets the configuration settings for using rdhs.
}
\details{
Setting up a configuration will enable API results to be cached, as
  well as enabling datasets from the DHS website to be downloaded and also
  cached. To enable results to be cached you have to either provide a valid
  `cache_path` argument, or allow rdhs to write to the user cache directory
  for your operating system. To do the later, leave the `cache_path` argument
  blank and you will be explicitly prompted to give permission to `rdhs` to
  save your results in this directory. If you do not then your API calls and
  any downloaded datasets will be saved in the temp directory and deleted
  after your R session closes. To allow `rdhs` to download datasets from the
  DHS website, you have to provide both an `email` and `project` argument.
  You will then be prompted to type in your login password securely.
  Your provided config (email, project, password, cache_path etc) will be
  saved at the location provided by `config_path`. If no argument is provided
  `config_path` will be either set to within your user cache directory if you
  have given permission to do so, otherwise it will be placed within your
  temp directory.

  When creating your config you also have the option to specify whether the
  `config_path` provided should be used as a local configuration or a global
  one. This is controlled using the `global` argument, which by default is
  set equal to `TRUE`. A global config is saved within your R root directory
  (the directory that a new R session will start in). If you set `global` to
  `FALSE` the config file will be saved within the current directory. This
  can be useful if you create a new DHS project for each new piece of work,
  and want to keep the datasets you download for this project separate to
  another. If you want to have your config file saved in a different
  directory, then you must create a file "rdhs.json" first in that directory
  before specifying the full path to it, as well as setting `global` equal to
  `FALSE`.

  As an aside, it is useful for the DHS program to see how the surveys they
  conducted are being used, and thus it is helpful for them if you do create
  a new project for each new piece of work (e.g. a different publication).
  However, we would still recommend setting up a global config and using
  the same `cache_path` for different projects as this will save you time
  downloading the same datasets as you have downloaded before.

  Lastly, you can decide how API calls from the DHS API are formatted by
  providing an argument for `data_frame`. If left blank API calls will be
  returned as `data.frame` objects, however, you could return API calls as
  `data.table` objects using `data.table::as.data.table`.
}
\examples{

\dontrun{
# normal set up we would prvide the email and project, and be prompted for
# the password. (not run as it requires a prompt)
set_rdhs_config(email = "blah@gmail.com", project = "Blahs",
config_path = "rdhs.json", global = FALSE)


# otherwise we can do this by specifying prompt to FALSE
set_rdhs_config(
config_path = "rdhs.json", global = FALSE, prompt = FALSE
)

# you can look at what you have set these to using \code{get_rdhs_config}
config <- get_rdhs_config()
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_datasets.R
\name{factor_format}
\alias{factor_format}
\title{reformat haven and labelled read ins to have no factors or labels}
\usage{
factor_format(res, reformat = FALSE, all_lower = TRUE)
}
\arguments{
\item{res}{dataset to be formatted}

\item{reformat}{Boolean whether to remove all factors and labels and
just return the unfactored data. Default = FALSE}

\item{all_lower}{Logical indicating whether all value labels
should be lower case. Default to `TRUE`.}
}
\value{
list with the formatted dataset and the code descriptions
}
\description{
reformat haven and labelled read ins to have no factors or labels
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{unzip_special}
\alias{unzip_special}
\title{unzip special that catches for 4GB+}
\usage{
unzip_special(
  zipfile,
  files = NULL,
  overwrite = TRUE,
  junkpaths = FALSE,
  exdir = ".",
  unzip = "internal",
  setTimes = FALSE
)
}
\arguments{
\item{zipfile}{The pathname of the zip file: tilde expansion (see
\code{\link{path.expand}} will be performed.)}

\item{files}{A character vector of recorded filepaths to be extracted:
    the default is to extract all files.}

\item{overwrite}{If \code{TRUE}, overwrite existing files (the equivalent
    of \command{unzip -o}), otherwise ignore such files (the equivalent of
    \command{unzip -n}).}

\item{junkpaths}{If \code{TRUE}, use only the basename of the stored
    filepath when extracting.  The equivalent of \command{unzip -j}.}

\item{exdir}{The directory to extract files to (the equivalent of
    \code{unzip -d}).  It will be created if necessary.}

\item{unzip}{The method to be used.  An alternative is to use
    \code{getOption("unzip")}, which on a Unix-alike may be set to the
    path to a \command{unzip} program.}

\item{setTimes}{logical.  For the internal method only, should the
    file times be set based on the times in the zip file?  (NB: this
    applies to included files, not to directories.)}
}
\description{
unzip special that catches for 4GB+
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client.R
\name{client_dhs}
\alias{client_dhs}
\title{Make a dhs client}
\usage{
client_dhs(config = NULL, root = rappdirs_rdhs(), api_key = NULL)
}
\arguments{
\item{config}{config object, as created using \code{read_rdhs_config}}

\item{root}{Character for root directory to where client, caches,
surveys etc. will be stored. Default = \code{rappdirs_rdhs()}}

\item{api_key}{Character for DHS API KEY. Default = NULL}
}
\description{
Make a DHS API client
}
\section{Methods}{


\describe{
\item{\code{dhs_api_request}}{
  Makes a call to the DHS websites API. You can make requests to any of their declared api endpoints (see \code{vignette(rdhs)} for more on these). API queries can be filtered by providing query terms, and you can control how many search results you want returned. The default parameters will return all of the results, and will format it nicely into a data.frame for you.
  N.B. This is easier to now do by using the bespoke functions that are included within the package. These take the form dhs_<endpoint>, e.g. \code{\link{dhs_data}}. These functions can also take your client as an argument that will cache the response for you

  \emph{Usage:}
  \code{dhs_api_request(api_endpoint, query = list(), api_key = private$api_key,
      num_results = 100, just_results = TRUE)}

  \emph{Arguments:}
  \itemize{
    \item{\code{api_endpoint}:   API endpoint. Must be one of the 12 possible endpoints.
    }

    \item{\code{query}:   List of query filters. To see possible query filter terms for each endpoint then head to the DHS api website.
    }

    \item{\code{api_key}:   DHS API key. Default will grab the key provided when the client was created.
    }

    \item{\code{num_results}:   The Number of results to return. Default = "ALL" which will loop through all the api search results pages for you if there are more results than their API will allow you to fetch in one page. If you specify a number this many results will be returned (but probably best to just leave default).
    }

    \item{\code{just_results}:   Boolean whether to return just the results or all the http API response. Default = TRUE (probably best again to leave as this.)
    }
  }

  \emph{Value}:
  Data.frame with search results if just_results=TRUE, otherwise a nested list with all the API responses for each page required.
}
\item{\code{available_datasets}}{
  Searches the DHS website for all the datasets that you can download. The results of this function are cached in the client. If you have recently requested new datasets from the DHS website then you can specify to clear the cache first so that you get the new set of datasets available to you.

  \emph{Usage:}
  \code{available_datasets(clear_cache_first = FALSE)}

  \emph{Arguments:}
  \itemize{
    \item{\code{clear_cache_first}:   Boolean detailing if you would like to clear the cached available datasets first. The default is set to FALSE. This option is available so that you can make sure your client fetches any new datasets that you have recently been given access to.
    }
  }

  \emph{Value}:
  Data.frame object with 14 variables that detail the surveys you can download, their url download links and the country, survey, year etc info for that link.
}
\item{\code{get_datasets}}{
  Gets datasets from your cache or downloads from the DHS website. By providing the filenames, as specified in one of the returned fields from \code{\link{dhs_datasets}}, the client will log in for you and download all the files you have requested. If any of the requested files are unavailable for your log in, these will be flagged up first as a message so you can make a note and request them through the DHS website. You also have the option to control whether the downloaded zip file is then extracted and converted into a more convenient R \code{data.frame}. This converted object will then be subsequently saved as a ".rds" object within the client root directory datasets folder, which can then be more quickly loaded when needed with \code{readRDS}. You also have the option to reformat the dataset, which will ensure that a suitable parser is used to preserve the meta information in your dataset, such as what different survey response codes mean.

  \emph{Usage:}
  \code{get_datasets(dataset_filenames, download_option = "rds", reformat = FALSE,
      all_lower = TRUE, output_dir_root = file.path(private$root,
          "datasets"), clear_cache = FALSE, ...)}

  \emph{Arguments:}
  \itemize{
    \item{\code{dataset_filenames}:   The desired filenames to be downloaded. These can be found as one of the returned fields from \code{\link{dhs_datasets}}. Alternatively you can also pass the desired rows from \code{\link{dhs_datasets}}.
    }

    \item{\code{download_option}:   Character specifying whether the dataset should be just downloaded ("zip"), imported and saved as an .rds object ("rds"), or both extract and rds ("both"). Conveniently you can just specify any letter from these options.
    }

    \item{\code{reformat}:   Boolean concerning whether to reformat read in datasets by removing all factors and labels. Default = FALSE.
    }

    \item{\code{all_lower}:   Logical indicating whether all value labels should be lower case. Default to `TRUE`.
    }

    \item{\code{output_dir_root}:   Root directory where the datasets will be stored within. The default will download datasets to a subfolder of the client root called "datasets"
    }

    \item{\code{clear_cache}:   Should your available datasets cache be cleared first. This will allow newly accessed datasets to be available. Default = `TRUE`
    }

    \item{\code{...}:   Any other arguments to be passed to \code{\link{read_dhs_dataset}}
    }
  }

  \emph{Value}:
  Depends on the download_option requested, but ultimately it is a file path to where the dataset was downloaded to, so that you can interact with it accordingly.
}
\item{\code{survey_questions}}{
  Use this function after download_survey to query downloaded surveys for what questions they asked. This function will look for the downloaded and imported survey datasets from the cache, and will download them if not previously downloaded.

  \emph{Usage:}
  \code{survey_questions(dataset_filenames, search_terms = NULL, essential_terms = NULL,
      regex = NULL, rm_na = TRUE, ...)}

  \emph{Arguments:}
  \itemize{
    \item{\code{dataset_filenames}:   The desired filenames to be downloaded. These can be found as one of the returned fields from \code{\link{dhs_datasets}}.
    }

    \item{\code{search_terms}:   Character vector of search terms. If any of these terms are found within the surveys question descriptions, the corresponding code and description will be returned.
    }

    \item{\code{essential_terms}:   Character pattern that has to be in the description of survey questions. I.e. the function will first find all survey_questions that contain your search terms (or regex) OR essential_terms. It will then remove any questions that did not contain your essential_terms. Default = NULL.
    }

    \item{\code{regex}:   Regex character pattern for matching. If you want to specify your regex search pattern, then specify this argument. N.B. If both search_terms and regex are supplied as arguments then regex will be ignored.
    }

    \item{\code{rm_na}:   Should NAs be removed. Default is `TRUE`
    }

    \item{\code{...}:   Any other arguments to be passed to \code{\link{download_datasets}}
    }
  }

  \emph{Value}:
  Data frame of the surveys where matches were found and then all the resultant codes and descriptions.
}
\item{\code{survey_variables}}{
  Use this function after download_survey to look up all the surveys that have the provided codes.

  \emph{Usage:}
  \code{survey_variables(dataset_filenames, variables, essential_variables = NULL,
      rm_na = TRUE, ...)}

  \emph{Arguments:}
  \itemize{
    \item{\code{dataset_filenames}:   The desired filenames to be downloaded. These can be found as one of the returned fields from \code{\link{dhs_datasets}}.
    }

    \item{\code{variables}:   Character vector of survey variables to be looked up
    }

    \item{\code{essential_variables}:   Character vector of variables that need to present. If any of the codes are not present in that survey, the survey will not be returned by this function. Default = NULL.
    }

    \item{\code{rm_na}:   Should NAs be removed. Default is `TRUE`
    }

    \item{\code{...}:   Any other arguments to be passed to \code{\link{download_datasets}}
    }
  }

  \emph{Value}:
  Data frame of the surveys where matches were found and then all the resultant codes and descriptions.
}
\item{\code{extract}}{
  Function to extract datasets using a set of survey questions as taken from the output from \code{survey_questions}

  \emph{Usage:}
  \code{extract(questions, add_geo = FALSE)}

  \emph{Arguments:}
  \itemize{
    \item{\code{questions}:   Questions to be queried, in the format from \code{survey_questions}
    }

    \item{\code{add_geo}:   Add geographic information to the extract. Default = TRUE
    }
  }
}
\item{\code{get_variable_labels}}{
  Returns information about a dataset's survey variables and definitions.

  \emph{Usage:}
  \code{get_variable_labels(dataset_filenames = NULL, dataset_paths = NULL, rm_na = FALSE)}

  \emph{Arguments:}
  \itemize{
    \item{\code{dataset_filenames}:   Vector of dataset filenames to look up
    }

    \item{\code{dataset_paths}:   Vector of dataset file paths to where datasets have been saved to
    }

    \item{\code{rm_na}:   Should variables and labels with NAs be removed. Default = FALSE
    }
  }

  \emph{Value}:
  Data frame of survey variable names and definitions
}
\item{\code{get_cache_date}}{
  Returns the private member variable cache-date, which is the date the client was last created/validated against the DHS API.

  \emph{Usage:}
  \code{get_cache_date()}

  \emph{Value}:
  POSIXct and POSIXt time
}
\item{\code{get_root}}{
  Returns the file path to the client's root directory

  \emph{Usage:}
  \code{get_root()}

  \emph{Value}:
  Character string file path
}
\item{\code{get_config}}{
  Returns the client's configuration

  \emph{Usage:}
  \code{get_config()}

  \emph{Value}:
  Config data.frame
}
\item{\code{get_downloaded_datasets}}{
  Returns a named list of all downloaded datasets and their file paths

  \emph{Usage:}
  \code{get_downloaded_datasets()}

  \emph{Value}:
  List of dataset names and file paths.
}
\item{\code{set_cache_date}}{
  Sets the private member variable cache-date, which is the date the client was last created/validated against the DHS API. This should never really be needed but is included to demonstrate the cache clearing properties of the client in the vignette.

  \emph{Usage:}
  \code{set_cache_date(date)}

  \emph{Arguments:}
  \itemize{
    \item{\code{date}:   POSIXct and POSIXt time to update cache time to.
    }
  }
}
\item{\code{save_client}}{
  Internally save the client object as an .rds file within the root directory for the client.

  \emph{Usage:}
  \code{save_client()}
}
\item{\code{clear_namespace}}{
  Clear the keys and values associated within a cache context. The dhs client caches a number of different tasks, and places these within specific contexts using the package \code{storr::storr_rds}.

  \emph{Usage:}
  \code{clear_namespace(namespace)}

  \emph{Arguments:}
  \itemize{
    \item{\code{namespace}:   Character string for the namespace to be cleared.
    }
  }
}
}
}

\examples{
\dontrun{
# create an rdhs config file at "rdhs.json"
conf <- set_rdhs_config(
config_path = "rdhs.json",global = FALSE, prompt = FALSE
)
td <- tempdir()
cli <- rdhs::client_dhs(api_key = "DEMO_1234", config = conf, root = td)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{get_variable_labels}
\alias{get_variable_labels}
\title{Get Survey Variable Labels}
\usage{
get_variable_labels(dataset, return_all = TRUE)
}
\arguments{
\item{dataset}{Can be either the file path to a dataset, the dataset as a
`data.frame` or the filenames of datasets. See details for more information}

\item{return_all}{Logical whether to return all variables (\code{TRUE})
or only those with labels.}
}
\value{
A \code{data.frame} consisting of the variable name and labels.
}
\description{
Return variable labels from a dataset
}
\details{
Returns variable names and their labels from a dataset.
  You can pass for the `data` argument any of
  the following:
  \itemize{
      \item The file path to a saved dataset. This would be the direct
      output of \code{\link{get_datasets}}
      \item A read in dataset, i.e. produced by using \code{\link{readRDS}}
      to load a dataset from
      a file path produced by \code{\link{get_datasets}}
      \item Dataset filenames. These can be found as one of the returned
      fields from \code{\link{dhs_datasets}}. If these datasets have not been
      downloaded before this will download them for you.
      }
}
\examples{
\dontrun{
# get the model datasets included with the package
model_datasets <- model_datasets

# download one of them
g <- get_datasets(dataset_filenames = model_datasets$FileName[1])

# we can pass the list of filepaths to the function
head(get_variable_labels(g))

# or we can pass the full dataset
r <- readRDS(g[[1]])
head(get_variable_labels(r))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{get_downloaded_datasets}
\alias{get_downloaded_datasets}
\title{Get Downloaded Datasets}
\usage{
get_downloaded_datasets()
}
\value{
A \code{data.frame} of downloaded datasets
}
\description{
Detail the datasets that you have already downloaded
}
\details{
Returns a \code{data.frame} of the datasets that have been
  downloaded within this client. This could be useful if you are without
  an internet connection and wish to know which saved
  dataset files in your root directory correspond to which dataset
}
\examples{
\dontrun{
# get the model datasets included with the package
model_datasets <- model_datasets

# download one of them
g <- get_datasets(dataset_filenames = model_datasets$FileName[1])

# these will then be stored so that we know what datasets we have downloaded
d <- get_downloaded_datasets()

# which returns a names list of file paths to the datasets
d[1]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_countries}
\alias{dhs_countries}
\title{API request of DHS Countries}
\usage{
dhs_countries(
  countryIds = NULL,
  indicatorIds = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  surveyCharacteristicIds = NULL,
  tagIds = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
  by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{indicatorIds}{Specify a comma separated list of indicators ids to
  filter by. For a list of indicators use
\code{dhs_indicators(returnFields=c("IndicatorId","Label","Definition"))}}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
For a list of surveys use \code{dhs_surveys(returnFields=c("SurveyId",
"SurveyYearLabel","SurveyType","CountryName"))}}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Countries
on. surveyYearStart is an inclusive value. Can be used alone or in
conjunction with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Countries on.
surveyYearEnd is an inclusive value. Can be used alone or in conjunction
with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{surveyCharacteristicIds}{Specify a survey characteristic id to filter
countries in surveys with the specified survey characteristic. For a list
of survey characteristics use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{tagIds}{Specify a tag id to filter countries with surveys containing
indicators with the specified tag. For a list of tags use
\code{dhs_tags()}}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request. By
default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 12 (or less if \code{returnFields} is provided)
  countries with their corresponding details. A detailed description of all
  the attributes returned is provided at
 \url{https://api.dhsprogram.com/rest/dhs/countries/fields}
}
\description{
API request of DHS Countries
}
\examples{

\dontrun{
# A common use for the countries API endpoint is to query which countries
# ask questions about a given topic. For example to find all countries that
# record data on malaria prevalence by RDT:

dat <- dhs_countries(indicatorIds = "ML_PMAL_C_RDT")

# Additionally you may want to know all the countries that have conducted
# MIS (malaria indicator surveys):

dat <- dhs_countries(surveyType="MIS")

# A complete list of examples for how each argument to the countries API
# endpoint can be provided is given below, which is a copy of each of
# the examples listed in the API at:

# https://api.dhsprogram.com/#/api-countries.cfm


dat <- dhs_countries(countryIds="EG",all_results=FALSE)
dat <- dhs_countries(indicatorIds="FE_FRTR_W_TFR",all_results=FALSE)
dat <- dhs_countries(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_countries(surveyYear="2010",all_results=FALSE)
dat <- dhs_countries(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_countries(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_countries(surveyType="DHS",all_results=FALSE)
dat <- dhs_countries(surveyCharacteristicIds="32",all_results=FALSE)
dat <- dhs_countries(tagIds="1",all_results=FALSE)
dat <- dhs_countries(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_dhs_dta.R
\name{read_dhs_dta}
\alias{read_dhs_dta}
\title{Read DHS Stata data set}
\usage{
read_dhs_dta(zfile, mode = "haven", all_lower = TRUE, ...)
}
\arguments{
\item{zfile}{Path to `.zip` file containing Stata dataset, usually ending
in filename `XXXXXXDT.zip`}

\item{mode}{Read mode for Stata `.dta` file. Defaults to "haven", see
'Details' for other options.}

\item{all_lower}{Logical indicating whether all value labels should be lower
case. Default to `TRUE`.}

\item{...}{Other arguments to be passed to \code{\link{read_zipdata}}.
Here this will be arguments to pass to either \code{\link[haven]{read_dta}}
or \code{\link[foreign]{read.dta}} depending on the mode provided}
}
\value{
A data frame. If mode = 'map', value labels for each variable are
  stored as the `labelled` class from `haven`.
}
\description{
This function reads a DHS recode dataset from the zipped Stata dataset.
By default (`mode = "haven"`), it reads in the stata data set using
\code{\link[haven]{read_dta}}
}
\details{
The default `mode="haven"` uses  \code{\link[haven]{read_dta}}
to read in the dataset. We have chosen this option as it is more consistent
with respect to variable labels and descriptions than others.
The other options either use use \code{\link[foreign]{read.dta}}
or they use the `.MAP` dictionary file provided with the DHS Stata datasets
to reconstruct the variable labels and value labels. In this case, value
labels are stored are stored using the the `labelled` class from `haven`.
See `?haven::labelled` for more information. Variable labels are stored in
the "label" attribute of each variable, the same as `haven::read_dta()`.

Currently, `mode="map"` is only implemented for 111
character fixed-width .MAP files, which comprises
the vast majority of recode data files from DHS Phases V,
VI, and VII and some from Phase IV. Parsers
for other .MAP formats will be added in future.

Other available modes read labels from the Stata dataset
with various options available in R:

* `mode="map"` uses the `.MAP` dictionary file provided with the DHS Stata
datasets to reconstruct the variable labels and value labels. In this case,
value labels are stored are stored using the the `labelled` class
from `haven`. See `?haven::labelled` for more information. Variable labels
are stored in the "label" attribute of each variable, the same as
`haven::read_dta()`.

* `mode="haven"`: use `haven::read_dta()` to read dataset.
This option retains the native value codings
with value labels affixed with the 'labelled' class.

* `mode="foreign"`: use `foreign::read.dta()`,
with default options convert.factors=TRUE to add
variable labels. Note that variable labels will
not be added if labels are not present for all
values, but variable labels are available via the "val.labels" attribute.

* `mode="foreignNA"`: use `foreign::read.dta(..., convert.factors=NA)`,
which converts any values without labels to 'NA'. This risks data loss
if labelling is incomplete in Stata datasets.

* `mode="raw"`: use `foreign::read.dta(..., convert.factors=FALSE)`,
which simply loads underlying value coding. Variable labels and value
labels are still available through dataset attributes (see examples).
}
\examples{
mrdt_zip <- tempfile()
download.file("https://dhsprogram.com/data/model_data/dhs/zzmr61dt.zip",
              mrdt_zip, mode="wb")

mr <- rdhs::read_dhs_dta(mrdt_zip,mode="map")
attr(mr$mv213, "label")
class(mr$mv213)
head(mr$mv213)
table(mr$mv213)
table(haven::as_factor(mr$mv213))

## If Stata file codebook is complete, `mode="map"` and `"haven"`
## should be the same.
mr_hav <- rdhs::read_dhs_dta(mrdt_zip, mode="haven")
attr(mr_hav$mv213, "label")
class(mr_hav$mv213)
head(mr_hav$mv213)  # "9=missing" omitted from .dta codebook
table(mr_hav$mv213)
table(haven::as_factor(mr_hav$mv213))

## Parsing codebook when using foreign::read.dta()
# foreign issues with duplicated factors
# Specifying foreignNA can help but often will not as below.
# Thus we would recommend either using mode = "haven" or mode = "raw"
\dontrun{
mr_for <- rdhs::read_dhs_dta(mrdt_zip, mode="foreign")
mr_for <- rdhs::read_dhs_dta(mrdt_zip, mode = "foreignNA")
}
## Don't convert factors
mr_raw <- rdhs::read_dhs_dta(mrdt_zip, mode="raw")
table(mr_raw$mv213)

}
\seealso{
\code{\link[foreign]{read.dta}}, \code{\link[haven]{labelled}},
  \code{\link[haven]{read_dta}}.

For more information on the DHS filetypes and contents of
distributed dataset .ZIP files, see
\url{https://dhsprogram.com/data/File-Types-and-Names.cfm#CP_JUMP_10334}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{last_api_update}
\alias{last_api_update}
\title{Pull last DHS API database update time}
\usage{
last_api_update(timeout = 30)
}
\arguments{
\item{timeout}{Numeric for API timeout. Default = 30}
}
\description{
Pull last DHS API database update time
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_survey_characteristics}
\alias{dhs_survey_characteristics}
\title{API request of DHS Survey Characteristics}
\usage{
dhs_survey_characteristics(
  countryIds = NULL,
  indicatorIds = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{indicatorIds}{Specify a comma separated list of indicators ids to
filter by. For a list of indicators use
\code{dhs_indicators(returnFields=c("IndicatorId","Label","Definition"))}}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Survey
Characteristics on. surveyYearStart is an inclusive value. Can be used
alone or in conjunction with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Survey
Characteristics on. surveyYearEnd is an inclusive value. Can be used alone
or in conjunction with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request. By
default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 2 (or less if \code{returnFields} is provided)
  survey characteristics. A survey can be labelled with one or more of these
  survey characteristics. A description of all the attributes returned is
  provided at
  \url{https://api.dhsprogram.com/rest/dhs/surveycharacteristics/fields}
}
\description{
API request of DHS Survey Characteristics
}
\examples{

\dontrun{
# A good use for the survey characteristics API endpoint is to query what the
# IDs are for each survey characteristic. These are useful for passing as
# arguments to other API endpoints.For example to show all the ids:

dat <- dhs_survey_characteristics()

# Or if your analysis is foucssed on a particular country, and you want to
# see all the characteristics surveyed for e.g. Senegal

dat <- dhs_countries(countryIds="SN")

# A complete list of examples for how each argument to the survey
# characteristics API endpoint can be provided is given below, which is a
# copy of each of the examples listed in the API at:

# https://api.dhsprogram.com/#/api-surveycharacteristics.cfm


dat <- dhs_survey_characteristics(countryIds="EG",all_results=FALSE)
dat <- dhs_survey_characteristics(indicatorIds="FE_FRTR_W_TFR",
all_results=FALSE)
dat <- dhs_survey_characteristics(surveyIds="SN2010DHS,all_results=FALSE")
dat <- dhs_survey_characteristics(surveyYear="2010,all_results=FALSE")
dat <- dhs_survey_characteristics(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_survey_characteristics(surveyYearStart="1991",
surveyYearEnd="2006",all_results=FALSE)
dat <- dhs_survey_characteristics(surveyType="DHS",all_results=FALSE)
dat <- dhs_survey_characteristics(f="html",all_results=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_dhs_dta.R
\name{parse_map}
\alias{parse_map}
\title{Create dictionary from DHS .MAP codebook}
\usage{
parse_map(map, all_lower = TRUE)
}
\arguments{
\item{map}{A character vector containing .MAP file, e.g. from `readLines()`.}

\item{all_lower}{Logical indicating whether all value labels should be
converted to lower case}
}
\value{
A data frame containing metadata, principally variable labels and
  a vector of value labels.
}
\description{
Create dictionary from DHS .MAP codebook
}
\details{
Currently hardcoded for 111 char width .MAP files, which covers the
  vast majority
  of DHS Phase V, VI, and VIII. To be extended in the future and perhaps add other useful options.
}
\examples{
mrdt_zip <- tempfile()
download.file("https://dhsprogram.com/data/model_data/dhs/zzmr61fl.zip",
              mrdt_zip, mode="wb")

map <- rdhs::read_zipdata(mrdt_zip, "\\\\.MAP", readLines)
dct <- rdhs:::parse_map(map)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_dhs_flat.R
\name{parse_meta}
\alias{parse_meta}
\alias{parse_dcf}
\alias{parse_sps}
\alias{parse_do}
\title{Parse fixed-width file metadata}
\usage{
parse_dcf(dcf, all_lower = TRUE)

parse_sps(sps, all_lower = TRUE)

parse_do(do, dct, all_lower = TRUE)
}
\arguments{
\item{dcf}{.DCF file path to parse}

\item{all_lower}{logical indicating whether to convert variable labels to
lower case. Defaults to `TRUE`.}

\item{sps}{.SPS file as character vector (e.g. from readLines / brio::read_lines)}

\item{do}{.DO file as character vector (e.g. from readLines / brio::read_lines)}

\item{dct}{.DCT file as character vector (e.g. from readLines / brio::read_lines)}
}
\value{
data.frame with metadata for parsing fixed-width flat file
}
\description{
Parse dataset metadata
}
\examples{
mrfl_zip <- tempfile()
download.file("https://dhsprogram.com/data/model_data/dhs/zzmr61fl.zip",
              mrfl_zip, mode = "wb")

dcf <- rdhs::read_zipdata(mrfl_zip, "\\\\.DCF", readLines)
dct <- rdhs:::parse_dcf(dcf)

sps <- rdhs::read_zipdata(mrfl_zip, "\\\\.SPS", readLines)
dct <- rdhs:::parse_sps(sps)

do <- rdhs::read_zipdata(mrfl_zip, "\\\\.DO", readLines)
dctin <- rdhs::read_zipdata(mrfl_zip, "\\\\.DCT", readLines)
dct <- rdhs:::parse_do(do, dctin)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui.R
\name{extract_dhs}
\alias{extract_dhs}
\title{Extract Data}
\usage{
extract_dhs(questions, add_geo = FALSE)
}
\arguments{
\item{questions}{Questions to be queried, in the format from
\code{\link{search_variables}} or \code{\link{search_variable_labels}}}

\item{add_geo}{Add geographic information to the extract. Defaut = `TRUE`}
}
\value{
A \code{list} of `data.frames` for each survey data extracted.
}
\description{
Extracts data from your downloaded datasets according to a data.frame of
requested survey variables or survey definitions
}
\details{
Function to extract datasets using a set of survey questions as
  taken from the output from \code{\link{search_variables}}
  or \code{\link{search_variable_labels}}
}
\examples{
\dontrun{
# get the model datasets included with the package
model_datasets <- model_datasets

# download one of them
g <- get_datasets(dataset_filenames = model_datasets$FileName[1])

# create some terms of data me may want to extrac
st <- search_variable_labels(names(g), "bed net")

# and now extract it
ex <- extract_dhs(st)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_endpoints.R
\name{dhs_tags}
\alias{dhs_tags}
\title{API request of DHS Tags}
\usage{
dhs_tags(
  countryIds = NULL,
  indicatorIds = NULL,
  surveyIds = NULL,
  surveyYear = NULL,
  surveyYearStart = NULL,
  surveyYearEnd = NULL,
  surveyType = NULL,
  f = NULL,
  returnFields = NULL,
  perPage = NULL,
  page = NULL,
  client = NULL,
  force = FALSE,
  all_results = TRUE
)
}
\arguments{
\item{countryIds}{Specify a comma separated list of country ids to filter
by. For a list of countries use
\code{dhs_countries(returnFields=c("CountryName","DHS_CountryCode"))}}

\item{indicatorIds}{Specify a comma separated list of indicators ids to
filter by. For a list of indicators use
\code{dhs_indicators(returnFields=c("IndicatorId","Label","Definition"))}}

\item{surveyIds}{Specify a comma separated list of survey ids to filter by.
For a list of surveys use
\code{dhs_surveys(returnFields=c("SurveyId","SurveyYearLabel",
"SurveyType","CountryName"))}}

\item{surveyYear}{Specify a comma separated list of survey years to
filter by.}

\item{surveyYearStart}{Specify a range of Survey Years to filter Tags on.
surveyYearStart is an inclusive value. Can be used alone or in
conjunction with surveyYearEnd.}

\item{surveyYearEnd}{Specify a range of Survey Years to filter Tags on.
surveyYearEnd is an inclusive value. Can be used alone or in conjunction
with surveyYearStart.}

\item{surveyType}{Specify a survey type to filter by.}

\item{f}{You can specify the format of the data returned from the query as
HTML, JSON, PJSON, geoJSON, JSONP, XML or CSV. The default data format
is JSON.}

\item{returnFields}{Specify a list of attributes to be returned.}

\item{perPage}{Specify the number of results to be returned per page. By
default the API will return 100 results.}

\item{page}{Allows specifying a page number to obtain for the API request. By
default the API will return page 1.}

\item{client}{If the API request should be cached, then provide a client
object created by \code{\link{client_dhs}}}

\item{force}{Should we force fetching the API results, and ignore any
cached results we have. Default = FALSE}

\item{all_results}{Boolean for if all results should be returned. If FALSE
then the specified page only will be returned. Default = TRUE.}
}
\value{
Returns a \code{data.table} of 4 (or less if \code{returnFields} is provided)
  tags with detailed information. An indicators can be tagged with one or
  more tags to help identify certain topics an indicator can be identified
  by. A description of the attributes returned is provided at
  \url{https://api.dhsprogram.com/rest/dhs/tags/fields}
}
\description{
API request of DHS Tags
}
\examples{

\dontrun{
# A good use for the tags API endpoint is to query what the
# IDs are for each tag. These are useful for passing as
# arguments to other API endpoints.For example to show all the ids:

dat <- dhs_tags()

# Or if your analysis is foucssed on a particular country, and you want to
# see all the characteristics surveyed for e.g. Senegal

dat <- dhs_tags(countryIds="SN")

# A complete list of examples for how each argument to the survey
# tags API endpoint can be provided is given below, which is a
# copy of each of the examples listed in the API at:

# https://api.dhsprogram.com/#/api-tags.cfm


dat <- dhs_tags(countryIds="EG",all_results=FALSE)
dat <- dhs_tags(indicatorIds="FE_FRTR_W_TFR",all_results=FALSE)
dat <- dhs_tags(surveyIds="SN2010DHS",all_results=FALSE)
dat <- dhs_tags(surveyYear="2010",all_results=FALSE)
dat <- dhs_tags(surveyYearStart="2006",all_results=FALSE)
dat <- dhs_tags(surveyYearStart="1991", surveyYearEnd="2006",
all_results=FALSE)
dat <- dhs_tags(surveyType="DHS",all_results=FALSE)
dat <- dhs_tags(f="html",all_results=FALSE)
}
}
