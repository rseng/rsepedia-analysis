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
