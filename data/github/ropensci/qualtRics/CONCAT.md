---
title: 'qualtRics: retrieve survey data using the Qualtrics API'
authors:
- affiliation: 1
  name: Jasper Ginn
  orcid: 0000-0002-5019-2923
date: "2 February 2018"
output: pdf_document
bibliography: paper.bib
tags:
- R
- survey data
- API
- Qualtrics
affiliations:
- index: 1
  name: Ecole Polytechnique Federale de Lausanne (EPFL)
---

# Summary

Qualtrics [@qualtrics_main] allows users to create and disseminate online surveys. It is used by researchers and other analysts to field responses for the purposes of (academic) research. While users can manually download survey responses from Qualtrics, importing this data into R is cumbersome. The R package `qualtRics` [@ginn_qualtrics_2017] focuses on the retrieval of survey data using the Qualtrics API and aims to reduce the pre-processing steps needed to prepare this data for analysis. Currently, the package is the only package on CRAN that offers such functionality, and is included in the official Qualtrics API documentation [@noauthor_getting_nodate].

The primary goal of the package is to provide a bridge between the Qualtrics user interface (where the survey is designed) and R (where the results can be analyzed) by using as few steps as possible. Users can store their API credentials in a file in an R project root directory that automatically loads when the library is called. This eliminates the need to remember API credentials and prevents the user from accidentally sharing sensitive information if they share their work. The package contains three core functions to retrieve survey data. The first of these functions - `getSurveys()` - retrieves a data frame containing an overview of surveys to which the user has access. Using a unique survey id, the user can download and import their data using `getSurvey()`. This function takes care of requesting, downloading and unpacking the data. It is then imported into R with the `readSurvey()` function. This last function can also be used to import manual data exports and supports both current and legacy data formats.

Apart from the above functionality, the package supports the automatic conversion of single-choice multiple choice questions. Using the rich metadata that Qualtrics provides about surveys, it is possible to automatically convert ordinal data to ordered factors. This functionality will be expanded on an ongoing basis to include other variable types.    

# References
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

# qualtRics

**Authors:** [Julia Silge](https://juliasilge.com/), [Jasper
Ginn](https://jasperhg90.github.io/)<br/> **License:**
[MIT](https://opensource.org/licenses/MIT)

<!-- badges: start -->

[![R-CMD-check](https://github.com/ropensci/qualtRics/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/qualtRics/actions)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/qualtRics)](https://cran.r-project.org/package=qualtRics)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/qualtRics/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/qualtRics?branch=master)
[![rOpenSci](https://badges.ropensci.org/192_status.svg)](https://github.com/ropensci/software-review/issues/192)
[![DOI](https://zenodo.org/badge/70817337.svg)](https://zenodo.org/badge/latestdoi/70817337)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.00690/status.svg)](https://doi.org/10.21105/joss.00690)
[![Downloads](https://cranlogs.r-pkg.org/badges/qualtRics)](https://CRAN.R-project.org/package=qualtRics)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/qualtRics?color=orange)](https://CRAN.R-project.org/package=qualtRics)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

[Qualtrics](https://www.qualtrics.com/) is an online survey and data
collection software platform. Qualtrics is used across many domains in
both academia and industry for online surveys and research. While users
can manually download survey responses from Qualtrics through a browser,
importing this data into R is then cumbersome. The qualtRics R package
implements the retrieval of survey data using the Qualtrics API and aims
to reduce the pre-processing steps needed in analyzing such surveys.
Currently, this package is the only package on CRAN that offers such
functionality, and is included in the official Qualtrics API
documentation.

Note that your institution must support API access and that it must be
enabled for your account. Whoever manages your Qualtrics account can
help you with this. Please refer to the [Qualtrics
documentation](https://api.qualtrics.com/) to find your API token.

The authors of this package are not affiliated with Qualtrics, and
Qualtrics does not offer support for this package. For specific
information about the Qualtrics API, you can refer to the [official
documentation](https://api.qualtrics.com/).

## Installation

This package can be installed from CRAN:

``` r
install.packages("qualtRics")
```

Alternatively, you can install the development version with the
[remotes](https://cran.r-project.org/package=remotes) package (or
alternatively, [devtools](https://cran.r-project.org/package=devtools)):

``` r
install.packages("remotes")
remotes::install_github("ropensci/qualtRics")
```

## Access your Qualtrics data

Currently, the package contains three core functions:

1.  `all_surveys()` fetches a list of all surveys that you own or have
    access to from Qualtrics.
2.  `fetch_survey()` downloads a survey from Qualtrics and loads it
    into R.
3.  `read_survey()` allows you to read CSV files you download manually
    from Qualtrics.

It also contains [multiple helper
functions](https://docs.ropensci.org/qualtRics/reference/index.html#other-helper-functions),
including:

1.  `qualtrics_api_credentials()` stores your API key and base URL in
    environment variables.
2.  `survey_questions()` retrieves a data frame containing questions and
    question IDs for a survey; `extract_colmap()` retrieves a similar
    data frame with more detailed mapping from columns to labels.
3.  `metadata()` retrieves metadata about your survey, such as
    questions, survey flow, number of responses etc.

Note that you can only export surveys that you own, or to which you have
been given administration rights.

## Register your Qualtrics credentials

There are two important credentials you need to authenticate with the
Qualtrics API. These are your **API key** and **datacenter-specific base
URL**. The base URL you pass to the qualtRics package should look like
`yourdatacenterid.qualtrics.com`, without a scheme such as `https://`.
The [Qualtrics API
documentation](https://api.qualtrics.com/instructions/) explains how you
can find your base URL.

You can store your API credentials `QUALTRICS_API_KEY` and
`QUALTRICS_BASE_URL` in your `.Renviron` file for repeated use across
sessions. The qualtRics package has a function to help with this.

``` r
library(qualtRics)

qualtrics_api_credentials(api_key = "<YOUR-QUALTRICS_API_KEY>", 
                          base_url = "<YOUR-QUALTRICS_BASE_URL>",
                          install = TRUE)
```

After you use this function, reload your environment
(`readRenviron("~/.Renviron")`) so you can use the credentials without
restarting R.

## A simple Qualtrics workflow

Once your Qualtrics API credentials are stored, you can see what surveys
are available to you.

``` r
surveys <- all_surveys() 
```

You can then download the data from any of these individual surveys (for
example, perhaps the sixth one) directly into R.

``` r
mysurvey <- fetch_survey(surveyID = surveys$id[6], 
                         verbose = TRUE)
```

See the qualtRics vignette for more details on variable metadata,
automatic conversion of variables, retrieving responses between specific
dates or for specific survey items, and more.

## Related work

-   [Jason Bryer](https://github.com/jbryer/qualtrics) wrote an R
    package to work with the previous version of the Qualtrics API
-   [QualtricsTools](https://github.com/emma-morgan/QualtricsTools)
    creates automatic reports in shiny.
-   [qsurvey](https://github.com/jamesdunham/qsurvey) by James Dunham
    focuses on testing and review of surveys before fielding, and
    analysis of responses afterward.

### Community Guidelines

This project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/qualtRics/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.
Feedback, bug reports (and fixes!), and feature requests are welcome;
file issues or seek support
[here](https://github.com/ropensci/qualtRics/issues).

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# qualtRics (development version)

- Add `fetch_distribution_history()` and `list_distribution_links()` for more handling of distribution data, thanks to @chrisumphlett and @dsen6644 (#221, #239)

- Changed handling of literal `"NA"` text input from users so it is no longer converted to an R `NA` value thanks to @jmobrien (#244)

- Use `httr::RETRY()` instead of `httr::VERB()` in `qualtrics_api_request()` to implement consistent API error-handling across all of the functions in the package. They will be retried up to 3 times if there is any non-4xx error. Thanks to @chrisumphlett (#217)

# qualtRics 3.1.5

- Add `fetch_description()` to download complete survey description metadata from v3 API endpoint (more up-to-date than older `metadata()`) thanks to @jmobrien (#207)
- Warn users about possible incorrect results from API when `breakout_sets` and `label` are both FALSE
- Refactor internal URL creation for API calls thanks to @jmobrien (#225)
- Add `fetch_id()` to return a `surveyID` based on a unique survey name as it appears in the Qualtrics UI thanks to @markjrieke (#230).  

# qualtRics 3.1.4

- Add `fetch_distributions()` to access distribution data for a specific survey thanks to @dsen6644 (#169)
- Handle mailing list embedded data better thanks to @dsen6644 (#175)
- Updated links to API documentation
- Create unique column names for questions using `choiceId` thanks to @lyh970817 (#182)
- Fix bug when `include_questions` only contains one QID thanks to @lyh970817 (#197)
- Generate correct/updated column mapping for survey responses thanks to @jmobrien (#199). These column mappings are available as an attribute on survey results or via the new `extract_colmap()` function. 

# qualtRics 3.1.3

- Update `include_questions` argument to use correct name in API request.
- Build API payloads with jsonlite (#155) thanks to @jmobrien
- Convert tests to webmockr and vcr (#140 and #161) thanks to @shaun-jacks and @dsen6644
- Allow user to specify column types for both `fetch_survey()` and `read_survey()` (#162) thanks to @jntrcs 

# qualtRics 3.1.2

- For empty surveys, return zero row dataframe (#127)
- Remove unnecessary dependency on yaml and deprecate `qualtRicsConfigFile()`, to avoid unexpected behavior
- Deprecate old versions of functions: `getSurveys()`, `getSurveyQuestions()`, `getSurvey()`, `readSurvey()`
- Move to updated version of Qualtrics API (#130) thanks to @jmobrien
- Correctly handle time zone conversions (#137) thanks to @jmobrien
- Add `breakout_sets` parameter thanks to @shaun-jacks
- Fix bug in `infer_data_types()` for answers choices that include HTML
- Deprecate `last_response` argument no longer used by API (#153)

# qualtRics 3.1.1

- Fix bug in `infer_data_types()` to avoid errors with factors/numeric values
- Improvements to documentation, error checking
- Allow user to access column mapping for questions and IDs (#115)
- Deprecate `registerOptions()` to avoid unexpected behavior with options
- Make data import more robust with more condition and error checking, as well as better defaults

# qualtRics 3.1.0

- New maintainer: Julia Silge
- Add all previous contributors to DESCRIPTION as `ctb`
- Declare testthat dependency in DESCRIPTION (reason for previous archiving from CRAN)
- Simpler approach for storing API credentials as environment variables with `qualtrics_api_credentials()` (`registerOptions()` is now soft deprecated with a warning)
- Simplify README (keep all existing detailed workflow documentation in vignette)
- Relicense from GPL-3 to MIT. See [consent from authors here](https://github.com/ropensci/qualtRics/issues/95).
- Improvements to documentation throughout
- Renaming (with warnings on old versions) of key functions for clarity and reduction in confusion, plus improvements:
  - `all_surveys()` (from old version of `getSurveys()`)
  - `survey_questions()` (from old version of `getSurveyQuestions()`)
  - `fetch_survey()` (from old version of `getSurvey()`)
  - `read_survey()` (from old version of `readSurvey()`)


qualtRics 3.0 (2018-02-03)
=========================

### NEW FEATURES

- Added 'metadata' function that allows the user to retrieve detailed metadata about survey.
- User can now convert specific question types automatically. See [this page](https://github.com/JasperHG90/qualtRics#automatic-conversion-of-variables) for more information.

### MINOR IMPROVEMENTS

- Using package [httptest](https://CRAN.R-project.org/package=httptest) for mock API requests so that API calls can be tested. 
- `getSurveys()` and `getSurveyQuestions()` now return a [tibble](https://CRAN.R-project.org/package=tibble)

### BUG FIXES

- Added .onDetach conditions so that environment variables (root url and API key) are removed when package is unloaded. This prevents issues if user decides to load the package again.
- We found that surveys that use new lines in the questions break the readSurvey function.
The problem is, that read.csv (and read.table as well as the readr library implementation) ignore the quote = "\"" option when a skip = 2 or skip = 3 parameter is set. As a result the read function slices off the questions row somewhere in the middle when first importing just the table body using skip.

### DEPRECATED AND DEFUNCT

- convertstandardcolumns deprecated since readr::read_csv does this automagically. It has been changed in config file to 'convertvariables'. 

qualtRics 2.2 (2017-10-27)
=========================

- `readSurvey()` now takes an additional argument, fileEncoding, so that user can import surveys using a specific encoding. 'fileEncoding' can also be passed as optional argument in `getSurvey()`. Added new parameter that reads legacy data format.
- Better argument checking.
- `getSurveyQuestions()` now returns additional information. h/t @phoebewong
- Fixes several bugs and stability issues
- More informative error messages

qualtRics 2.0 (2017-06-16)
=========================

- `registerOptions()` now takes more arguments. User can now set global options. See `qualtRicsConfigFile()` for more information. Same options are now passed through `...` in specific functions.
- Added appveyor testing.
- Added support for a configuration file to store API key and root url in the working directory.
- `registerApiKey()` has been replaced by `registerOptions()`. This function stores both a user's API key and root url. Function also scans for a configuration file `.qualtRics.yml` that contains this information.
- Added a new script called `zzz.R`. When the package is loaded, the .onLoad() function in this file scans the working directory for a `.qualtRics.yml` configuration file so that the user doesn't have to register this information manually.
- Added a new function `qualtRicsConfigFile()` that prints instructions for the user on how to set up a configuration file to the R Console.
- Removed the `root_url` parameter from all functions that required it.
- Dates are now converted without a specific timezone.
- Added a new function `getSurveyQuestions()` that allows the user to download a data frame containing question labels and IDs.
- Added parameter **includedQuestionIds** so user can select questions they want to download. Need to use the QID value from `getSurveyQuestions()`.
- Updated examples and documentation of functions.
- Added the following parameters to `getSurvey()`:
  - **seenUnansweredRecode:**  String. Recode seen but unanswered questions with a string value.
  - **limit:** Integer. Maximum number of responses exported. Defaults to NULL (all responses).
  - **useLocalTime:** Boolean. Use local timezone to determine response date values. 
- `getSurveys()` now retrieves > 100 results.

qualtRics 1.0 (2016-10-13)
=========================

- Added a new function `readSurvey()`. This function is used in the `getSurvey()` function but will also work with surveys downloaded manually from Qualtrics. Standard columns (completed survey/startDate/endDate etc.) are now converted to their proper data types. HT Adrian Brugger & Stefan Borer.
- User can only download surveys in CSV format, no longer in XML or JSON. 
- Added several new parameters to `getSurvey()` function. HT @samuelkaminsky & @eknud
  * *LastResponseId*: If used, only responses that were filled out later than this ID will be downloaded. 
  * *UseLabels*: If TRUE, download will contain character labels. Else, download will contain choice labels.
  * *StartDate*: Only download responses after this date.
  * *EndDate*: Only download responses before this date.
- Survey downloads should be faster now; `getSurvey()` no longer sleeps when checking download status. Also added progress bar.

qualtRics 0.03 [beta] 
=========================

- User can choose to save results directly in a folder through 'save_dir' parameter in `getSurvey()`
- Results can now be requested in .csv, .json or .xml format. The packages `jsonlite` and `XML` are added to 'Suggests' in DESCRIPTION.
- `constructHeader()` is now deprecated and should no longer be used. Instead, users need to call `registerApiKey()`.
- Added a new function `registerApiKey()` which saves the user's API key and header information in the `tempdir()` environment. 

qualtRics 0.02 [beta] 
=========================

- Renamed 'base url' to 'root url' such that it corresponds to Qualtrics documentation.
- The root url no longer requires API-specific endpoints. So e.g. 'https://leidenuniv.eu.qualtrics.com' now works for all functions. The API-specific endpoints are added in the functions itself.
- Institution-specific root url is now required by `getSurvey()`

qualtRics 0.01 [beta] 
=========================

- Added first three functions (`constructHeader`, `getSurvey`, `getSurveyIDs`)
- base_url parameter is now uniform across functions. Parameter is called 'root url' to bring it in line with Qualtrics documentation.
## Release summary

This is the 12th CRAN release of qualtRics (the 6th since it has returned to CRAN since being archived). This release improves user messaging and internal functions, and adds functionality for survey description data (`fetch_description()`) and finding survey IDs (`fetch_id()`).

## Test environments

* local macOS install: release
* macOS 10.15.7 (on GitHub actions): release
* windows server 2019 10.0.17763 (on GitHub actions): release
* ubuntu 20.04 (on GitHub actions): release, devel
* win-builder: release, devel

## R CMD check results

0 errors | 0 warnings | 0 note
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# qualtRics

**Authors:** [Julia Silge](https://juliasilge.com/), [Jasper Ginn](https://jasperhg90.github.io/)<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/qualtRics/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/qualtRics/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/qualtRics)](https://cran.r-project.org/package=qualtRics)
[![Codecov test coverage](https://codecov.io/gh/ropensci/qualtRics/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/qualtRics?branch=master)
[![rOpenSci](https://badges.ropensci.org/192_status.svg)](https://github.com/ropensci/software-review/issues/192)
[![DOI](https://zenodo.org/badge/70817337.svg)](https://zenodo.org/badge/latestdoi/70817337)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.00690/status.svg)](https://doi.org/10.21105/joss.00690)
[![Downloads](https://cranlogs.r-pkg.org/badges/qualtRics)](https://CRAN.R-project.org/package=qualtRics)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/qualtRics?color=orange)](https://CRAN.R-project.org/package=qualtRics)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

[Qualtrics](https://www.qualtrics.com/) is an online survey and data collection software platform. Qualtrics is used across many domains in both academia and industry for online surveys and research. While users can manually download survey responses from Qualtrics through a browser, importing this data into R is then cumbersome. The qualtRics R package implements the retrieval of survey data using the Qualtrics API and aims to reduce the pre-processing steps needed in analyzing such surveys. Currently, this package is the only package on CRAN that offers such functionality, and is included in the official Qualtrics API documentation. 

Note that your institution must support API access and that it must be enabled for your account. Whoever manages your Qualtrics account can help you with this. Please refer to the [Qualtrics documentation](https://api.qualtrics.com/) to find your API token.

The authors of this package are not affiliated with Qualtrics, and Qualtrics does not offer support for this package. For specific information about the Qualtrics API, you can refer to the [official documentation](https://api.qualtrics.com/).


## Installation

This package can be installed from CRAN:

```{r eval=FALSE}
install.packages("qualtRics")
```


Alternatively, you can install the development version with the [remotes](https://cran.r-project.org/package=remotes) package (or alternatively, [devtools](https://cran.r-project.org/package=devtools)):

```{r eval=FALSE}
install.packages("remotes")
remotes::install_github("ropensci/qualtRics")
```


## Access your Qualtrics data

Currently, the package contains three core functions:

1. `all_surveys()` fetches a list of all surveys that you own or have access to from Qualtrics.
2. `fetch_survey()` downloads a survey from Qualtrics and loads it into R.
3. `read_survey()` allows you to read CSV files you download manually from Qualtrics.

It also contains [multiple helper functions](https://docs.ropensci.org/qualtRics/reference/index.html#other-helper-functions), including:

1. `qualtrics_api_credentials()` stores your API key and base URL in environment variables.
2. `survey_questions()` retrieves a data frame containing questions and question IDs for a survey; `extract_colmap()` retrieves a similar data frame with more detailed mapping from columns to labels.
3. `metadata()` retrieves metadata about your survey, such as questions, survey flow, number of responses etc.

Note that you can only export surveys that you own, or to which you have been given administration rights.

## Register your Qualtrics credentials

There are two important credentials you need to authenticate with the Qualtrics API. These are your **API key** and **datacenter-specific base URL**. The base URL you pass to the qualtRics package should look like `yourdatacenterid.qualtrics.com`, without a scheme such as `https://`. The [Qualtrics API documentation](https://api.qualtrics.com/instructions/) explains how you can find your base URL.

You can store your API credentials `QUALTRICS_API_KEY` and `QUALTRICS_BASE_URL` in your `.Renviron` file for repeated use across sessions. The qualtRics package has a function to help with this.

```{r, eval=FALSE}
library(qualtRics)

qualtrics_api_credentials(api_key = "<YOUR-QUALTRICS_API_KEY>", 
                          base_url = "<YOUR-QUALTRICS_BASE_URL>",
                          install = TRUE)
```

After you use this function, reload your environment (`readRenviron("~/.Renviron")`) so you can use the credentials without restarting R.

## A simple Qualtrics workflow

Once your Qualtrics API credentials are stored, you can see what surveys are available to you.

```{r, eval=FALSE}
surveys <- all_surveys() 
```

You can then download the data from any of these individual surveys (for example, perhaps the sixth one) directly into R.

```{r, eval=FALSE}
mysurvey <- fetch_survey(surveyID = surveys$id[6], 
                         verbose = TRUE)
```


See the qualtRics vignette for more details on variable metadata, automatic conversion of variables, retrieving responses between specific dates or for specific survey items, and more.

## Related work

- [Jason Bryer](https://github.com/jbryer/qualtrics) wrote an R package to work with the previous version of the Qualtrics API
- [QualtricsTools](https://github.com/emma-morgan/QualtricsTools) creates automatic reports in shiny.
- [qsurvey](https://github.com/jamesdunham/qsurvey) by James Dunham focuses on testing and review of surveys before fielding, and analysis of responses afterward.


### Community Guidelines

This project is released with a [Contributor Code of Conduct](https://github.com/ropensci/qualtRics/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms. Feedback, bug reports (and fixes!), and feature requests are welcome; file issues or seek support [here](https://github.com/ropensci/qualtRics/issues).


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction to qualtRics"
author: "Julia Silge and Jasper Ginn"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to qualtRics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)
```

[Qualtrics](https://www.qualtrics.com/) is an online survey and data collection software platform. Qualtrics is used across many domains in both academia and industry for online surveys and research. While users can manually download survey responses from Qualtrics through a browser, importing this data into R is then cumbersome. The qualtRics R package implements the retrieval of survey data using the Qualtrics API and aims to reduce the pre-processing steps needed in analyzing such surveys. 

Note that your institution must support API access and that it must be enabled for your account. Whoever manages your Qualtrics account can help you with this. Please refer to the [Qualtrics documentation](https://api.qualtrics.com/instructions/) to find your API token.

The authors and contributors for this R package are not affiliated with Qualtrics and Qualtrics does not offer support for this R package.

## Usage

Currently, the package contains three core functions:

1. `all_surveys()` fetches a list of all surveys that you own or have access to from Qualtrics.
2. `fetch_survey()` downloads a survey from Qualtrics and loads it into R.
3. `read_survey()` allows you to read CSV files you download manually from Qualtrics.

It also contains a number of helper functions, including:

1. `qualtrics_api_credentials()` stores your API key and base url in environment variables.
2. `survey_questions()` retrieves a data frame containing questions and question IDs for a survey; `extract_colmap()` retrieves a similar data frame with more detailed mapping from columns to labels.
3. `metadata()` retrieves metadata about your survey, such as questions, survey flow, number of responses etc.

Note that you can only export surveys that you own, or to which you have been given administration rights.

## Registering your Qualtrics credentials

There are two important credentials you need to authenticate with the Qualtrics API. These are your **API key** and **datacenter-specific base URL**. The base URL you pass to the qualtRics package should look like `yourdatacenterid.qualtrics.com`, without a scheme such as `https://`. The [Qualtrics API documentation](https://api.qualtrics.com/instructions/) explains how you can find your base URL.

You can store your API credentials `QUALTRICS_API_KEY` and `QUALTRICS_BASE_URL` in your `.Renviron` file for repeated use across sessions. The qualtRics package has a function to help with this.

```{r, eval=FALSE}
library(qualtRics)

qualtrics_api_credentials(api_key = "<YOUR-QUALTRICS_API_KEY>", 
                          base_url = "<YOUR-QUALTRICS_BASE_URL>",
                          install = TRUE)
```

After you use this function, reload your environment (`readRenviron("~/.Renviron")`) so you can use the credentials without restarting R.

## A simple Qualtrics workflow

Once your Qualtrics API credentials are stored, you can see what surveys are available to you.

```{r, eval=FALSE}
surveys <- all_surveys() 
```

You can then download the data from any of these individual surveys (for example, perhaps the sixth one) directly into R.

```{r, eval=FALSE}
mysurvey <- fetch_survey(surveyID = surveys$id[6], 
                         verbose = TRUE)
```

## More detailed control

You can add date parameters to only retrieve responses between certain dates:

```{r, eval=FALSE}
mysurvey <- fetch_survey(surveys$id[4],
                         start_date = "2018-10-01",
                         end_date = "2018-10-31",
                         label = FALSE)
```


Note that your date and time settings may not correspond to your own timezone. You can find out how to do this [here](https://www.qualtrics.com/support/survey-platform/managing-your-account/research-core-account-settings/#user-settings). See ["Dates and Times"](https://api.qualtrics.com/instructions/) for more information about how Qualtrics handles dates and times. **Keep in mind that this is important if you plan on using times / dates as cut-off points to filter data**.

You may also reference a response ID; `fetch_survey()` will then download all responses that were submitted after that response:

```{r eval=FALSE}
mysurvey <- fetch_survey(surveys$id[4],
                         last_response = "R_3mmovCIeMllvsER",
                         label = FALSE,
                         verbose = TRUE)
```


You can filter a survey for specific questions:


```{r eval=FALSE}
# what are the questions in a certain survey?
questions <- survey_questions(surveyID = surveys$id[6])

# download that survey, filtering for only certain questions
mysurvey <- fetch_survey(surveyID = surveys$id[6],
                         save_dir = tempdir(),
                         include_questions = c("QID1", "QID2", "QID3"),
                         verbose = TRUE)
```



You can store the results in a specific location if you like:

```{r, eval=FALSE}
mysurvey <- fetch_survey(surveyID = surveys$id[6], 
                         save_dir = "/users/Julia/Desktop/", 
                         verbose = TRUE)
```

Note that surveys that are stored in this way will be saved as an [RDS](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) file rather than e.g. a CSV. Reading an RDS file can be done like so:

```{r, eval=FALSE}
mysurvey <- readRDS(file = "/users/Julia/Desktop/mysurvey.rds")
```


You can read a survey that you downloaded manually from Qualtrics' website via a browser using `read_survey()`:

```{r, eval=FALSE}
mysurvey <- read_survey("/users/Julia/Desktop/mysurvey.csv")
```


To avoid special characters (mainly periods) in header names, `read_survey()` uses question labels as the header names. The question belonging to that label is then added using the [sjlabelled](https://CRAN.R-project.org/package=sjlabelled) package. Qualtrics gives names to these labels automatically, but you can easily change them.

```{r, echo=FALSE, out.width="80%"}
knitr::include_graphics("qualtricsdf.png")
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{getSurveys}
\alias{getSurveys}
\title{Retrieve a data frame of all active surveys on Qualtrics}
\usage{
getSurveys()
}
\description{
This function is deprecated; use \code{\link[=all_surveys]{all_surveys()}}
instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{infer_data_types}
\alias{infer_data_types}
\title{Set proper data types on survey data.}
\usage{
infer_data_types(data, surveyID, verbose = FALSE)
}
\arguments{
\item{data}{Imported Qualtrics survey}

\item{surveyID}{ID of survey}

\item{verbose}{Flag}
}
\description{
Set proper data types on survey data.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_survey.R
\name{read_survey}
\alias{read_survey}
\title{Read a CSV file exported from Qualtrics}
\usage{
read_survey(
  file_name,
  strip_html = TRUE,
  import_id = FALSE,
  time_zone = NULL,
  legacy = FALSE,
  add_column_map = TRUE,
  add_var_labels = TRUE,
  col_types = NULL
)
}
\arguments{
\item{file_name}{String. A CSV data file.}

\item{strip_html}{Logical. If \code{TRUE}, then remove HTML tags. Defaults
to \code{TRUE}.}

\item{import_id}{Logical. If \code{TRUE}, use Qualtrics import IDs instead of
question IDs as column names. Defaults to \code{FALSE}.}

\item{time_zone}{String. A local timezone to determine response date
values. Defaults to \code{NULL} which corresponds to UTC time. See
\href{https://api.qualtrics.com/instructions/}{"Dates and Times"} from Qualtrics
for more information on format.}

\item{legacy}{Logical. If \code{TRUE}, then import "legacy" format CSV files
(as of 2017). Defaults to \code{FALSE}.}

\item{add_column_map}{Logical. If \code{TRUE}, then a column map data frame
will be added as an attribute to the main response data frame.
This column map captures Qualtrics-provided metadata associated with the
response download, such as an item description and internal ID's. Defaults to
\code{TRUE}.}

\item{add_var_labels}{Logical. If \code{TRUE}, then the item description from
each variable (equivalent to the one in the column map) will be added as a
"label" attribute using \code{\link[sjlabelled:set_label]{sjlabelled::set_label()}}. Useful for
reference as well as cross-compatibility with other stats packages (e.g.,
Stata, see documentation in \code{sjlabelled}). Defaults to \code{TRUE}.}

\item{col_types}{Optional. This argument provides a way to manually overwrite
column types that may be incorrectly guessed. Takes a \code{\link[readr:cols]{readr::cols()}}
specification. See example below and \code{\link[readr:cols]{readr::cols()}} for formatting
details. Defaults to \code{NULL}.}
}
\value{
A data frame. Variable labels are stored as attributes. They are not
printed on the console but are visibile in the RStudio viewer.
}
\description{
Reads comma separated CSV files generated by Qualtrics
software. The second line containing the variable labels is imported.
Repetitive introductions to matrix questions are automatically removed.
Variable labels are stored as attributes.
}
\examples{
\dontrun{
# Generic use of read_survey()
df <- read_survey("<YOUR-PATH-TO-CSV-FILE>")
}
# Example using current data format
file <- system.file("extdata", "sample.csv", package = "qualtRics")
df <- read_survey(file)

# Example using legacy data format
file <- system.file("extdata", "sample_legacy.csv", package = "qualtRics")
df <- read_survey(file, legacy = TRUE)

# Example changing column type
file <- system.file("extdata", "sample.csv", package = "qualtRics")
# Force EndDate to be a string
df <- read_survey(file, col_types = readr::cols(EndDate = readr::col_character()))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualtrics_api_credentials.R
\name{qualtrics_api_credentials}
\alias{qualtrics_api_credentials}
\title{Install Qualtrics credentials in your \code{.Renviron} file for repeated use}
\usage{
qualtrics_api_credentials(
  api_key,
  base_url,
  overwrite = FALSE,
  install = FALSE
)
}
\arguments{
\item{api_key}{The API key provided to you from Qualtrics formatted in quotes.
Learn more about Qualtrics API keys at \url{https://api.qualtrics.com/}}

\item{base_url}{The institution-specific base URL for your Qualtrics account,
formatted in quotes. Find your base URL at \url{https://api.qualtrics.com/}}

\item{overwrite}{If TRUE, will overwrite existing Qualtrics
credentials that you already have in your \code{.Renviron} file.}

\item{install}{If TRUE, will install the key in your \code{.Renviron} file
for use in future sessions.  Defaults to FALSE (single session use).}
}
\description{
This function adds your Qualtrics API key and base URL to your
\code{.Renviron} file so it can be called securely without being stored in
your code. After you have installed these two credentials, they can be
called any time with \code{Sys.getenv("QUALTRICS_API_KEY")} or
\code{Sys.getenv("QUALTRICS_BASE_URL")}. If you do not have an
\code{.Renviron} file, the function will create one for you. If you already
have an \code{.Renviron} file, the function will append the key to your
existing file, while making a backup of your original file for disaster
recovery purposes.
}
\examples{

\dontrun{
qualtrics_api_credentials(
  api_key = "<YOUR-QUALTRICS_API_KEY>",
  base_url = "<YOUR-QUALTRICS_BASE_URL>",
  install = TRUE
)
# Reload your environment so you can use the credentials without restarting R
readRenviron("~/.Renviron")
# You can check it with:
Sys.getenv("QUALTRICS_API_KEY")

# If you need to overwrite existing credentials:
qualtrics_api_credentials(
  api_key = "<YOUR-QUALTRICS_API_KEY>",
  base_url = "<YOUR-QUALTRICS_BASE_URL>",
  overwrite = TRUE,
  install = TRUE
)
# Reload your environment to use the credentials
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_surveys.R
\name{all_surveys}
\alias{all_surveys}
\title{Retrieve a data frame of all active surveys on Qualtrics}
\usage{
all_surveys()
}
\description{
Retrieve a data frame of all active surveys on Qualtrics
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of all surveys
surveys <- all_surveys()

# Retrieve a single survey
mysurvey <- fetch_survey(surveyID = surveys$id[6])

mysurvey <- fetch_survey(
  surveyID = surveys$id[6],
  save_dir = tempdir(),
  start_date = "2018-01-01",
  end_date = "2018-01-31",
  limit = 100,
  label = TRUE,
  unanswer_recode = "UNANS",
  verbose = TRUE
)
}

}
\seealso{
See \url{https://api.qualtrics.com/} for documentation on the
Qualtrics API.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{getSurveyQuestions}
\alias{getSurveyQuestions}
\title{Retrieve a data frame containing question IDs and labels}
\usage{
getSurveyQuestions(...)
}
\arguments{
\item{...}{All arguments for \code{survey_questions}}
}
\description{
This function is deprecated; use \code{\link[=survey_questions]{survey_questions()}}
instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{construct_header}
\alias{construct_header}
\title{Construct a header to send to Qualtrics API}
\usage{
construct_header(API_TOKEN)
}
\arguments{
\item{API_TOKEN}{API token. Available in your Qualtrics account (see: \url{https://api.qualtrics.com/})}
}
\description{
Construct a header to send to Qualtrics API
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_mailinglist.R
\name{fetch_mailinglist}
\alias{fetch_mailinglist}
\title{Download a mailing list from Qualtrics}
\usage{
fetch_mailinglist(mailinglistID)
}
\arguments{
\item{mailinglistID}{String. Unique ID for the mailing list you want to
download. Returned as \code{id} by the \link[=all_mailinglists]{all_mailinglists}
function.}
}
\description{
Download a mailing list from Qualtrics
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of all mailing lists
mailinglists <- all_mailinglists()

# Retrieve a single mailing list
mailinglist <- fetch_mailinglist(mailinglists$id[1])
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_distribution_history.R
\name{fetch_distribution_history}
\alias{fetch_distribution_history}
\title{Download distribution history data for a distribution from Qualtrics}
\usage{
fetch_distribution_history(distributionID)
}
\arguments{
\item{distributionID}{String. Unique distribution ID for the distribution history you want to download.}
}
\description{
Download distribution history data for a distribution from Qualtrics
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

surveys <- all_surveys()
distributions <- fetch_distributions(surveys$id[1])
distribution_history <- fetch_distribution_history(distributions$id[1])
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_colmap.R
\name{extract_colmap}
\alias{extract_colmap}
\title{Extract column map from survey data download}
\usage{
extract_colmap(respdata)
}
\arguments{
\item{respdata}{Response data including a column map dataframe as an attribute}
}
\description{
Helper function to extract the column map attached to a response data
download obtained from \code{\link[=fetch_survey]{fetch_survey()}} (using the
default \code{add_column_map = TRUE})
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Retrieve a list of surveys
surveys <- all_surveys()

# Retrieve a single survey
mysurvey <- fetch_survey(surveyID = surveys$id[6])

# Extract column mapping for survey
extract_colmap(mysurvey)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{create_raw_payload}
\alias{create_raw_payload}
\title{Create raw JSON payload to post response exports request}
\usage{
create_raw_payload(
  label = TRUE,
  start_date = NULL,
  end_date = NULL,
  limit = NULL,
  time_zone = NULL,
  unanswer_recode = NULL,
  unanswer_recode_multi = NULL,
  include_display_order = TRUE,
  include_questions = NULL,
  breakout_sets = NULL
)
}
\arguments{
\item{label}{Flag}

\item{start_date}{Flag}

\item{end_date}{Flag}

\item{limit}{Flag}

\item{time_zone}{Flag}

\item{unanswer_recode}{Flag}

\item{unanswer_recode_multi}{Flag}

\item{include_display_order}{Flag}

\item{include_questions}{Flag}

\item{breakout_sets}{Flag}
}
\value{
JSON file with options to send to API
}
\description{
Create raw JSON payload to post response exports request
}
\seealso{
See \code{\link[=all_surveys]{all_surveys()}} for more details on these parameters
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{wrapper_mc}
\alias{wrapper_mc}
\title{Convert multiple choice questions to ordered factors}
\usage{
wrapper_mc(data, question_meta)
}
\arguments{
\item{data}{Imported Qualtrics survey}

\item{question_meta}{Question metadata}
}
\description{
Convert multiple choice questions to ordered factors
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{metadata}
\alias{metadata}
\title{Download metadata for a survey}
\usage{
metadata(surveyID, get = NULL, ...)
}
\arguments{
\item{surveyID}{A string. Unique ID for the survey you want to download.
Returned as "id" by the \link{all_surveys} function.}

\item{get}{A character vector containing any of the following: "metadata",
"questions", "responsecounts", "blocks", "flow", "embedded_data",
or "comments". Will return included elements. By default, the function
returns the "metadata", "questions", and "responsecounts" elements.
See examples below for more information.}

\item{...}{Additional options. User may pass an argument called \code{questions},
a vector containing the names of questions for which you want to
return metadata.}
}
\description{
Using this function, you can retrieve metadata about your survey. This
information includes question metadata (type, options, choices, etc), number
of responses, general metadata, survey flow, etc.
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of surveys
surveys <- all_surveys()

# Get metadata for a survey
md <- metadata(surveyID = surveys$id[6])

# Get metadata with specific elements
md_specific <- metadata(
  surveyID = id,
  get = c("flow")
)

# Get specific question metadata
question_specific <- metadata(
  surveyID = id,
  get = c("questions"),
  questions = c("Q1", "Q2")
)

# Example of a metadata file
file <- system.file("extdata", "metadata.rds", package = "qualtRics")

# Load
metadata_ex <- readRDS(file = file)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_mailinglists.R
\name{all_mailinglists}
\alias{all_mailinglists}
\title{Retrieve a data frame of all mailing lists from Qualtrics}
\usage{
all_mailinglists()
}
\description{
Retrieve a data frame of all mailing lists from Qualtrics
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of all mailing lists
mailinglists <- all_mailinglists()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_description.R
\name{fetch_description}
\alias{fetch_description}
\title{Download complete survey description using the Qualtrics v3 "Get Survey"
API endpoint.}
\usage{
fetch_description(surveyID, elements = NULL, legacy = FALSE, ...)
}
\arguments{
\item{surveyID}{A string. Unique ID for the survey you want to download.
Returned as "id" by the \link{all_surveys} function.}

\item{elements}{A character vector. Lists elements of survey definition to be
maintained.  Possible elements are "metadata", "surveyoptions", "flow",
"blocks", "questions", "responsesets", and/or "scoring" (case-insensitive).
If \code{legacy = TRUE}, then possible elements are "metadata", "questions",
"responsecounts", "blocks", "flow", "embedded_data", and/or "comments".}

\item{legacy}{Logical.  If TRUE, will use older Get Survey API endpoint
via a call to legacy function \link{metadata}.}

\item{...}{Additional options, only used when \code{legacy = TRUE}. User may pass
an argument called \code{questions}, a vector containing the names of questions
for which you want to return metadata.}
}
\value{
A list containing survey description metadata. The contents of the
returned list depend on \code{elements}.
}
\description{
Download complete survey description using the Qualtrics v3 "Get Survey"
API endpoint.
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of surveys
surveys <- all_surveys()

# Get description for a survey
descrip <- fetch_description(surveyID = surveys$id[6])

# Get metadata with specific elements
descrip_specific <- fetch_description(
  surveyID = id,
  elements = c("questions", "flow")
)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_distribution_links.R
\name{list_distribution_links}
\alias{list_distribution_links}
\title{Download distribution links for a distribution from Qualtrics}
\usage{
list_distribution_links(distributionID, surveyID)
}
\arguments{
\item{distributionID}{String. Unique distribution ID for the distribution links you want to download.}

\item{surveyID}{String. Unique ID for the survey you want to download.}
}
\description{
Download distribution links for a distribution from Qualtrics
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

surveys <- all_surveys()
distributions <- fetch_distributions(surveys$id[1])
distribution_links <- list_distribution_links(distributions$id[1], surveyID = surveys$id[1])
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{getSurvey}
\alias{getSurvey}
\title{Download a survey and import it into R}
\usage{
getSurvey(...)
}
\arguments{
\item{...}{All arguments for \code{fetch_survey}}
}
\description{
This function is deprecated; use \code{\link[=fetch_survey]{fetch_survey()}}
instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_distributions.R
\name{fetch_distributions}
\alias{fetch_distributions}
\title{Download distribution data for a survey from Qualtrics}
\usage{
fetch_distributions(surveyID)
}
\arguments{
\item{surveyID}{String. Unique survey ID for the distribution data you want to download.}
}
\description{
Download distribution data for a survey from Qualtrics
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

surveys <- all_surveys()
distributions <- fetch_distributions(surveys$id[1])
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_id.R
\name{fetch_id}
\alias{fetch_id}
\title{Fetch a unique Qualtrics survey ID based on survey name in the Qualtrics UI}
\usage{
fetch_id(.data, survey_name)
}
\arguments{
\item{.data}{Data frame of active surveys created by the function
\code{\link[=all_surveys]{all_surveys()}}.}

\item{survey_name}{Name of the survey as it appears in the Qualtrics UI. Must
be unique to be passed to \code{fetch_id()}.}
}
\description{
Fetch a unique Qualtrics survey ID based on survey name in the Qualtrics UI
}
\details{
Survey names in the Qualtrics platform are not required to be
unique, but the \code{survey_name} argument for this function \emph{must} be unique.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of surveys
surveys <- all_surveys()

# Retrieve surveyID for a unique survey
my_id <- fetch_id(surveys, "Unique Survey Name")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_for_warnings}
\alias{check_for_warnings}
\title{Check if httr GET result contains a warning}
\usage{
check_for_warnings(resp)
}
\arguments{
\item{resp}{object returned by \code{\link[=qualtrics_response_codes]{qualtrics_response_codes()}}}
}
\description{
Check if httr GET result contains a warning
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{generate_url}
\alias{generate_url}
\title{Generate URL for specific API query by type and (if appropriate) ID}
\usage{
generate_url(query, ...)
}
\arguments{
\item{query}{string.  The specific API query desired.  Generally named the
same as associated functions but without underscores, so the request for
\code{fetch_survey()} would be be "fetchsurvey".}

\item{...}{Named elements of URL for specific query desired, such as
\code{surveyID} or \code{mailinglistID}}
}
\value{
Endpoint URL to be passed to querying tools
}
\description{
Generate URL for specific API query by type and (if appropriate) ID
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_survey.R
\name{fetch_survey}
\alias{fetch_survey}
\title{Download a survey and import it into R}
\usage{
fetch_survey(
  surveyID,
  last_response = deprecated(),
  start_date = NULL,
  end_date = NULL,
  unanswer_recode = NULL,
  unanswer_recode_multi = unanswer_recode,
  include_display_order = TRUE,
  limit = NULL,
  include_questions = NULL,
  save_dir = NULL,
  force_request = FALSE,
  verbose = TRUE,
  label = TRUE,
  convert = TRUE,
  import_id = FALSE,
  time_zone = NULL,
  breakout_sets = TRUE,
  add_column_map = TRUE,
  add_var_labels = TRUE,
  col_types = NULL
)
}
\arguments{
\item{surveyID}{String. Unique ID for the survey you want to download.
Returned as \code{id} by the \link[=all_surveys]{all_surveys} function.}

\item{last_response}{Deprecated.}

\item{start_date}{String. Filter to only exports responses recorded after the
specified date. Accepts dates as character strings in format "YYYY-MM-DD".
Defaults to \code{NULL}.}

\item{end_date}{String. Filter to only exports responses recorded before the
specified date. Accepts dates as character strings in format "YYYY-MM-DD".
Defaults to \code{NULL}.}

\item{unanswer_recode}{Integer. Recode seen but unanswered questions with an
integer-like value, such as 999. Defaults to \code{NULL}.}

\item{unanswer_recode_multi}{Integer. Recode seen but unanswered multi-select
questions with an integer-like value, such as 999. Defaults to value for
\code{unaswer_recode}.}

\item{include_display_order}{Display order information (such as for
surveys with randomization).}

\item{limit}{Integer. Maximum number of responses exported. Defaults to
\code{NULL} (all responses).}

\item{include_questions}{Vector of strings (e.g. c('QID1', 'QID2', 'QID3').
Export only specified questions. Defaults to \code{NULL}.}

\item{save_dir}{String. Directory where survey results will be stored.
Defaults to a temporary directory which is cleaned when your R session is
terminated. This argument is useful if you'd like to store survey results.
The downloaded survey will be stored as an RDS file (see
\code{\link[base:readRDS]{base::readRDS()}}).}

\item{force_request}{Logical. fetch_survey() saves each survey in a temporary
directory so that it can quickly be retrieved later. If force_request is
\code{TRUE}, fetch_survey() always downloads the survey from the API instead
of loading it from the temporary directory. Defaults to \code{FALSE}.}

\item{verbose}{Logical. If \code{TRUE}, verbose messages will be printed to
the R console. Defaults to \code{TRUE}.}

\item{label}{Logical. \code{TRUE} to export survey responses as Choice Text
or \code{FALSE} to export survey responses as values.}

\item{convert}{Logical. If \code{TRUE}, then the
\code{\link[=fetch_survey]{fetch_survey()}} function will convert certain question
types (e.g. multiple choice) to proper data type in R. Defaults to \code{TRUE}.}

\item{import_id}{Logical. If \code{TRUE}, use Qualtrics import IDs instead of
question IDs as column names. Will also alter names in the column map, if
used. Defaults to \code{FALSE}.}

\item{time_zone}{String. A local timezone to determine response date
values. Defaults to \code{NULL} which corresponds to UTC time. See
\href{https://api.qualtrics.com/instructions/}{"Dates and Times"} from Qualtrics
for more information on format.}

\item{breakout_sets}{Logical. If \code{TRUE}, then the
\code{\link[=fetch_survey]{fetch_survey()}} function will split multiple
choice question answers into columns. If \code{FALSE}, each multiple choice
question is one column. Defaults to \code{TRUE}.}

\item{add_column_map}{Logical. If \code{TRUE}, then a column map data frame
will be added as an attribute to the main response data frame.
This column map captures Qualtrics-provided metadata associated with the
response download, such as an item description and internal ID's. Defaults to
\code{TRUE}.}

\item{add_var_labels}{Logical. If \code{TRUE}, then the item description from
each variable (equivalent to the one in the column map) will be added as a
"label" attribute using \code{\link[sjlabelled:set_label]{sjlabelled::set_label()}}. Useful for
reference as well as cross-compatibility with other stats packages (e.g.,
Stata, see documentation in \code{sjlabelled}). Defaults to \code{TRUE}.}

\item{col_types}{Optional. This argument provides a way to manually overwrite
column types that may be incorrectly guessed. Takes a \code{\link[readr:cols]{readr::cols()}}
specification. See example below and \code{\link[readr:cols]{readr::cols()}} for formatting
details. Defaults to \code{NULL}. Overwritten by \code{convert = TRUE}.}
}
\description{
Download a Qualtrics survey you own via API and import the survey directly into R.
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of surveys
surveys <- all_surveys()

# Retrieve a single survey
mysurvey <- fetch_survey(surveyID = surveys$id[6])

mysurvey <- fetch_survey(
  surveyID = surveys$id[6],
  save_dir = tempdir(),
  start_date = "2018-01-01",
  end_date = "2018-01-31",
  limit = 100,
  label = TRUE,
  unanswer_recode = 999,
  verbose = TRUE,
  # Manually override EndDate to be a character vector
  col_types = readr::cols(EndDate = readr::col_character())
)

}

}
\seealso{
See \url{https://api.qualtrics.com/} for documentation on
the Qualtrics API.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survey_questions.R
\name{survey_questions}
\alias{survey_questions}
\title{Retrieve a data frame containing question IDs and labels}
\usage{
survey_questions(surveyID)
}
\arguments{
\item{surveyID}{A string. Unique ID for the survey you want to download.
Returned as \code{id} by the \link[=all_surveys]{all_surveys} function.}
}
\description{
Retrieve a data frame containing question IDs and labels
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of surveys
surveys <- all_surveys()

# Retrieve questions for a survey
questions <- survey_questions(surveyID = surveys$id[6])

# Retrieve a single survey, filtering for specific questions
mysurvey <- fetch_survey(
  surveyID = surveys$id[6],
  save_dir = tempdir(),
  include_questions = c("QID1", "QID2", "QID3"),
  verbose = TRUE
)
}

}
\seealso{
See \url{https://api.qualtrics.com/} for documentation on the
Qualtrics API.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_params}
\alias{check_params}
\title{Check if parameters passed to functions are correct}
\usage{
check_params(...)
}
\arguments{
\item{...}{options passed to function}
}
\description{
Check if parameters passed to functions are correct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{readSurvey}
\alias{readSurvey}
\title{Read a CSV file exported from Qualtrics}
\usage{
readSurvey(...)
}
\arguments{
\item{...}{All arguments for \code{\link[=read_survey]{read_survey()}}}
}
\description{
This function is deprecated; use \code{\link[=read_survey]{read_survey()}}
instead.
Reads comma separated CSV files generated by Qualtrics
software. The second line containing the variable labels is imported.
Repetitive introductions to matrix questions are automatically removed.
Variable labels are stored as attributes.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/column_map.R
\name{column_map}
\alias{column_map}
\title{Retrieve a data frame containing survey column mapping}
\usage{
column_map(surveyID)
}
\arguments{
\item{surveyID}{A string. Unique ID for the survey you want to download.
Returned as \code{id} by the \link[=all_surveys]{all_surveys} function.}
}
\description{
Retrieve a data frame containing survey column mapping
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\examples{
\dontrun{
# Register your Qualtrics credentials if you haven't already
qualtrics_api_credentials(
  api_key = "<YOUR-API-KEY>",
  base_url = "<YOUR-BASE-URL>"
)

# Retrieve a list of surveys
surveys <- all_surveys()

# Retrieve column mapping for a survey
mapping <- column_map(surveyID = surveys$id[6])

# Retrieve a single survey, filtering for specific questions
mysurvey <- fetch_survey(
  surveyID = surveys$id[6],
  save_dir = tempdir(),
  include_questions = c("QID1", "QID2", "QID3"),
  verbose = TRUE
)
}

}
\seealso{
See \url{https://api.qualtrics.com/} for documentation on the
Qualtrics API.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{qualtrics_api_request}
\alias{qualtrics_api_request}
\title{Send httr requests to Qualtrics API}
\usage{
qualtrics_api_request(verb = c("GET", "POST"), url = url, body = NULL)
}
\arguments{
\item{verb}{Type of request to be sent (@seealso \code{\link[httr:VERB]{httr::VERB()}})}

\item{url}{Qualtrics endpoint URL created by \code{\link[=generate_url]{generate_url()}} functions}

\item{body}{Options created by \code{\link[=create_raw_payload]{create_raw_payload()}} function}
}
\description{
Send httr requests to Qualtrics API
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{qualtrics_response_codes}
\alias{qualtrics_response_codes}
\title{Checks responses against Qualtrics response codes and returns error message.}
\usage{
qualtrics_response_codes(res, raw = FALSE)
}
\arguments{
\item{res}{Response from httr::GET}

\item{raw}{If TRUE, add 'raw' flag to httr::content() function.}
}
\description{
Checks responses against Qualtrics response codes and returns error message.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{qualtRicsConfigFile}
\alias{qualtRicsConfigFile}
\title{Prints an Example of a QualtRics Configuration File to the Console.}
\usage{
qualtRicsConfigFile(...)
}
\arguments{
\item{...}{All arguments for \code{qualtRicsConfigFile}}
}
\description{
This function is deprecated; use \code{\link[=qualtrics_api_credentials]{qualtrics_api_credentials()}} instead.
}
\examples{
\dontrun{
# Execute this line to get instructions on how to make a .qualtrics.yml config file.
qualtRicsConfigFile()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{download_qualtrics_export}
\alias{download_qualtrics_export}
\title{Download response export}
\usage{
download_qualtrics_export(fetch_url, requestID, verbose = FALSE)
}
\arguments{
\item{fetch_url}{URL provided by Qualtrics API that shows the download percentage completeness}

\item{requestID}{ID}

\item{verbose}{See \code{\link[=fetch_survey]{fetch_survey()}}}
}
\description{
Download response export
}
\details{
If the request to the Qualtrics API made by this function fails, the request
will be retried. If you see these failures on a 500 error (such as a 504
error) be patient while the request is retried; it will typically succeed
on retrying. If you see other types of errors, retrying is unlikely to help.
}
\keyword{internal}
