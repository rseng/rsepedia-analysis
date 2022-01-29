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
