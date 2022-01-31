
<!-- README.md is generated from README.Rmd. Please edit that file -->

# medrxivr <img src="man/figures/hex-medrxivr.png" align="right" width="20%" height="20%" />

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![rOpenSci
Badge](https://badges.ropensci.org/380_status.svg)](https://github.com/ropensci/software-review/issues/380)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02651/status.svg)](https://doi.org/10.21105/joss.02651)
[![CRAN
Downloads.](https://cranlogs.r-pkg.org/badges/grand-total/medrxivr)](https://CRAN.R-project.org/package=medrxivr)
<br> [![R build
status](https://github.com/ropensci/medrxivr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/medrxivr/actions)
[![Travis build
status](https://travis-ci.com/ropensci/medrxivr.svg?branch=master)](https://travis-ci.com/ropensci/medrxivr)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/medrxivr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/medrxivr?branch=master)

<!-- badges: end -->

An increasingly important source of health-related bibliographic content
are preprints - preliminary versions of research articles that have yet
to undergo peer review. The two preprint repositories most relevant to
health-related sciences are [medRxiv](https://www.medrxiv.org/) and
[bioRxiv](https://www.biorxiv.org/), both of which are operated by the
Cold Spring Harbor Laboratory.

The goal of the `medrxivr` R package is two-fold. In the first instance,
it provides programmatic access to the [Cold Spring Harbour Laboratory
(CSHL) API](https://api.biorxiv.org/), allowing users to easily download
medRxiv and bioRxiv preprint metadata (e.g. title, abstract, publication
date, author list, etc) into R. The package also provides access to a
maintained static snapshot of the medRxiv repository (see [Data
sources](#medrxiv-data)). Secondly, `medrxivr` provides functions to
search the downloaded preprint records using regular expressions and
Boolean logic, as well as helper functions that allow users to export
their search results to a .BIB file for easy import to a reference
manager and to download the full-text PDFs of preprints matching their
search criteria.

## Installation

To install the stable version of the package from CRAN:

``` r
install.packages("medrxivr")
library(medrxivr)
```

Alternatively, to install the development version from GitHub, use the
following code:

``` r
install.packages("devtools")
devtools::install_github("ropensci/medrxivr")
library(medrxivr)
```

## Data sources

### medRxiv data

`medrixvr` provides two ways to access medRxiv data:

  - `mx_api_content(server = "medrxiv")` creates a local copy of all
    data available from the medRxiv API at the time the function is run.

<!-- end list -->

``` r
# Get a copy of the database from the live medRxiv API endpoint
preprint_data <- mx_api_content()  
```

  - `mx_snapshot()` provides access to a static snapshot of the medRxiv
    database. The snapshot is created each morning at 6am using
    `mx_api_content()` and is stored as CSV file in the [medrxivr-data
    repository](https://github.com/mcguinlu/medrxivr-data). This method
    does not rely on the API (which can become unavailable during peak
    usage times) and is usually faster (as it reads data from a CSV
    rather than having to re-extract it from the API). Discrepancies
    between the most recent static snapshot and the live database can be
    assessed using `mx_crosscheck()`.

<!-- end list -->

``` r
# Get a copy of the database from the daily snapshot
preprint_data <- mx_snapshot()  
```

The relationship between the two methods for the medRxiv database is
summarised in the figure below:

<img src="vignettes/data_sources.png" width="500px" height="400px" />

### bioRxiv data

Only one data source exists for the bioRxiv repository:

  - `mx_api_content(server = "biorxiv")` creates a local copy of all
    data available from the bioRxiv API endpoint at the time the
    function is run. **Note**: due to it’s size, downloading a complete
    copy of the bioRxiv repository in this manner takes a long time (\~
    1 hour).

<!-- end list -->

``` r
# Get a copy of the database from the live bioRxiv API endpoint
preprint_data <- mx_api_content(server = "biorxiv")
```

## Performing your search

Once you have created a local copy of either the medRxiv or bioRxiv
preprint database, you can pass this object (`preprint_data` in the
examples above) to `mx_search()` to search the preprint records using an
advanced search strategy.

``` r
# Import the medrxiv database
preprint_data <- mx_snapshot()
#> Using medRxiv snapshot - 2021-01-28 09:31

# Perform a simple search
results <- mx_search(data = preprint_data,
                     query ="dementia")
#> Found 192 record(s) matching your search.

# Perform an advanced search
topic1  <- c("dementia","vascular","alzheimer's")  # Combined with Boolean OR
topic2  <- c("lipids","statins","cholesterol")     # Combined with Boolean OR
myquery <- list(topic1, topic2)                    # Combined with Boolean AND

results <- mx_search(data = preprint_data,
                     query = myquery)
#> Found 70 record(s) matching your search.
```

You can also explore which search terms are contributing most to your
search by setting `report = TRUE`:

``` r
results <- mx_search(data = preprint_data,
                     query = myquery,
                     report = TRUE)
#> Found 70 record(s) matching your search.
#> Total topic 1 records: 1078
#> dementia: 192
#> vascular: 917
#> alzheimer's: 0
#> Total topic 2 records: 203
#> lipids: 74
#> statins: 25
#> cholesterol: 136
```

## Further functionality

### Export records identified by your search to a .BIB file

Pass the results of your search above (the `results` object) to the
`mx_export()` to export references for preprints matching your search
results to a .BIB file so that they can be easily imported into a
reference manager (e.g. Zotero, Mendeley).

``` r
mx_export(data = results,
          file = "mx_search_results.bib")
```

### Download PDFs for records returned by your search

Pass the results of your search above (the `results` object) to the
`mx_download()` function to download a copy of the PDF for each record
found by your search.

``` r
mx_download(results,        # Object returned by mx_search(), above
            "pdf/",         # Directory to save PDFs to 
            create = TRUE)  # Create the directory if it doesn't exist
```

## Accessing the raw API data

By default, the `mx_api_*()` functions clean the data returned by the
API for use with other `medrxivr` functions.

To access the raw data returned by the API, the `clean` argument should
set to `FALSE`:

``` r
mx_api_content(to_date = "2019-07-01", clean = FALSE)
```

See [this
article](https://docs.ropensci.org/medrxivr/articles/medrxiv-api.html#accessing-the-raw-api-data)
for more details.

## Detailed guidance

Detailed guidance, including advice on how to design complex search
strategies, is available on the [`medrxivr`
website.](https://docs.ropensci.org/medrxivr/)

## Linked repositories

See here for the [code used to take the daily
snapshot](https://github.com/mcguinlu/medrxivr-data) and [the code that
powers the `medrxivr` web
app](https://github.com/mcguinlu/medrxivr-app).

## Other tools/packages for working with medRxiv/bioRxiv data

The focus of `medrxivr` is on providing tools to allow users to import
and then search medRxiv and bioRxiv data. Below are a list of
complementary packages that provide distinct but related functionality
when working with medRxiv and bioRxiv data:

  - [`rbiorxiv`](https://github.com/nicholasmfraser/rbiorxiv) by
    [Nicholas Fraser](https://github.com/nicholasmfraser) provides
    access to the same medRxiv and bioRxiv *content* data as `medrxivr`,
    but also provides access to the *usage* data (e.g. downloads per
    month) that the Cold Spring Harbour Laboratory API offers. This is
    useful if you wish to explore, for example, [how the number of PDF
    downloads from bioRxiv has grown over
    time.](https://github.com/nicholasmfraser/rbiorxiv#pdf-downloads-over-time)

## Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Disclaimer

This package and the data it accesses/returns are provided “as is”, with
no guarantee of accuracy. Please be sure to check the accuracy of the
data yourself (and do let me know if you find an issue so I can fix it
for everyone\!)
# medrxivr (development version)

# medrxivr 0.0.5

Major changes:

* Improved error handling to address a common bug that causes extraction from the API to fail. The "total number" of records element of the API metadata is frequently artificially inflated. This leads to an overestimation of the number of pages of records, which in turn caused the extraction function to fail at the very end when `mx_api_content()` encounters an empty page. This error has been changed to informative messaging about the expected (as per the metadata) and actual (`nrows()` of returned dataset) number of retrievable records.
* New functionality added to `mx_search()` allows users to view the number of "hits" (records returned) for each individual element of the search. An extra parameter called `report` has been added, which gives the user the option to switch this functionality on or off. The default value for this parameter is set to FALSE. This functionality was added by [James O'Hare](https://github.com/jamesohare1) in response to [Issue #13.](https://github.com/ropensci/medrxivr/issues/13)
* Users can now pass a vector of terms to the `NOT` parameter rather than a single exclusion term.
* New functionality to allow for user-friendly search operators, including wildcards ("randomi*ation" will now find "randomisation" and "randomization") and the NEAR operator ("systematic NEAR1 review" will find "systematic review" and "systematic _<any-other-word>_ review")
* A new argument, `auto_caps`, in the `mx_search()` function to allow for automatic capitalisation of the first character in each search term (e.g. with `auto_caps = TRUE`, "dementia" will be automatically converted to "[Dd]ementia" which will find "**d**ementia" and also "**D**ementia"). This replaces the recommendation that users capitalise the first character themselves using square brackets. However, if user defined alternative are already in place for the first character of the search term, then these are left untouched.
* A helper function, `mx_caps()`, allows users to wrap search terms to find all possible combinations of upper- and lower-case letters in that term. For example, `mx_caps("ncov")` converts the term to "[Nn][Cc][Oo][Vv]" which will find "NCOV", "Ncov", "nCOV", "nCoV", etc.

# medrxivr 0.0.4

Major changes:

* Fixed error which occurred when downloading the whole bioRxiv database. This was caused by any record above 100000 being presented in scientific notation (e.g. 1e+05), which meant the API returned an invalid response.
* Change tests to fix runtime regardless of future growth of the repositories

# medrxivr 0.0.3

Version created for submission to JOSS and CRAN, and onboarded to rOpenSci following peer-review. 

Major changes:

* `mx_snapshot()` now takes a `commit` argument, allowing you to specify exactly which snapshot of the database you would like to use. Details on the commit keys needed are [here](https://github.com/mcguinlu/medrxivr-data/commits/master/snapshot.csv). In addition, the process of taking the snapshot is now managed by GitHub actions, meaning it should be a lot more robust/regular/
* Importing the snapshot to R is now significantly faster, as `vroom::vroom()` is used in place of `read.csv()`
* All functions that return a data frame now return ungrouped tibbles.
* The  to/from date arguments for both `mx_search()` and `mx_api_content()` have been standardized to snake case and now expect the same "YYYY-MM-DD" character format.
* A progress indicator has been added to `mx_api_content()` provide useful information when downloading from the API.
* Some refactoring of code has taken place to reduce duplication of code chunks and to make future maintenance easier.

Minor changes:

* `mx_crosscheck()` no longer uses web-scraping when providing the number of 
* Documentation has been updated to reflect the changes, and some additional sections added to the vignettes. This includes removing references to older versions of the functions names (e.g. `mx_raw()`).
* Additional test have been written, and the overall test coverage has been increased. Some lines (handling exceptional rare errors that can't be mocked) have been marked as `#nocov`.
* \dontrun had been replaced with \donttest in all examples across the package. 
* All examples for mx_download() and mx_export() now use tempfile() and tempdir(), so as not to modify the users home filespace when running the examples.




# medrxivr 0.0.2

Major changes:  

* Following the release of the [medRxiv API](https://api.biorxiv.org/), the way the snapshot of the medRxiv site is taken has changed, resulting in a more accurate snapshot of the entire repository being taken daily (as opposed to just new articles being captured, as was previously the case). This has introduced some breaking changes (e.g. in the `fields` argument, "subject" has become "category", and "link" has become "doi"), but will result in better long-term stability of the package.
* Two new functions, `mx_api_content()` and `mx_api_doi()`, have been added to allow users to interact with the medRxiv API endpoints directly. A new vignette documenting these functions has been added. 
* The API has also allowed for improved data collection. The "authors" variable searched/returned now contains all authors of the paper as opposed to just the first one. Several additional fields are now returned, including corresponding author's institution, preprint license, and the DOI of the published, peer-reviewed version of preprint (if available).
* A companion app was launched, which allows you to build the search strategy using a user-friendly interface and then export the code needed to run it directly from R. 
* You can now define the field(s) you wish to search. By default, the Title, Abstract, First Author, Subject, and Link (which includes the DOI) fields are searched. 
* There is no longer a limit on the number of distinct topics you can search for (previously it was 5).
* The output of `mx_search()` has been cleaned to make it more useful to future end-users. Of note, some of the columns names have changed, and the "pdf_name" and "extraction_date" variables are no longer returned.


# medrxivr 0.0.1

* Added a `NEWS.md` file to track changes to the package.
## Required submission to retain package on CRAN
This submission was prompted by an email from CRAN, requiring that tests fail 
gracefully when the internet resource the package is based on is unavailable.

## Test environments

* local windows R installation, R 4.0.2
* windows-latest (via GitHub actions), (release)
* macOS-latest (via GitHub actions), (release)
* ubuntu-20.04 (via GitHub actions), (release, devel)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## CRAN Checks 

CRAN checks are currently failing on a single platform 
(r-devel-windows-ix86+x86_64), but the changes contained in this release will 
address this.

## Downstream dependencies

There are currently no downstream dependencies for this package.
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

Please note that the medrxivr project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

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
---
title: 'medrxivr: Accessing and searching medRxiv and bioRxiv preprint data in R'
tags:
  - R
  - systematic review
  - evidence synthesis
  - bibliographic database
authors:
  - name: Luke A McGuinness
    orcid: 0000-0003-0872-7098
    affiliation: "1, 2" 
  - name: Lena Schmidt
    orcid: 0000-0003-0709-8226
    affiliation: 1
affiliations:
  - name: Department of Population Health Science, University of Bristol
    index: 1
  - name: MRC Intergrative Epidemiology Unit, University of Bristol
    index: 2
date: 11 August 2020
bibliography: paper.bib
---

# Summary

An increasingly important source of health-related bibliographic content are preprints: preliminary versions of research articles that have yet to undergo peer review. The two preprint repositories most relevant to health-related sciences are [medRxiv](https://www.medrxiv.org) and [bioRxiv](https://www.biorxiv.org), both of which are operated by the Cold Spring Harbor Laboratory, a not-for-profit research and educational institution [@rawlinson2019].

The goal of the `medrxivr` R package is two-fold. In the first instance, it provides programmatic access to the Cold Spring Harbour Laboratory (CSHL) API, allowing users to download medRxiv and bioRxiv preprint metadata (e.g., title, abstract, author list.) This functionality will be of interest to anyone who wishes to import medRxiv and/or bioRxiv preprint metadata into R, for example to explore the distribution of preprints by subject area or by publication year. Examples of this type of usage have already been reported [e.g., by @Brierley].

In the second instance, the package provides functions that allow users to search the downloaded preprint metadata for relevant preprints using complex search strings, including functionality such as search term truncation, Boolean operators (AND, OR, NOT), and term proximity. Helper functions are provided that allow users to export the results of their search to a .bib file for import into a reference manager (e.g., Zotero) and to download the full-text PDFs of preprints matching their search. This aspect of the package will be more relevant to systematic reviewers, health librarians and others performing literature searches, allowing them to perform and document transparent and reproducible searches in these important evidence sources.

# Acknowledgements

We acknowledge funding  from NIHR (LAM through NIHR Doctoral Research Fellowship (DRF-2018-11-ST2-048), and LS through NIHR Systematic Reviews Fellowship (RM-SR-2017-09-028)). LAM is a member of the MRC Integrative Epidemiology Unit at the University of Bristol. The views expressed in this article are those of the authors and do not necessarily represent those of the NHS, the NIHR, MRC, or the Department of Health and Social Care.

# References
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = FALSE,
  message = FALSE,
  warning = FALSE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)

library(medrxivr)
```

# medrxivr <img src="man/figures/hex-medrxivr.png" align="right" width="20%" height="20%" />


<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![rOpenSci Badge](https://badges.ropensci.org/380_status.svg)](https://github.com/ropensci/software-review/issues/380)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02651/status.svg)](https://doi.org/10.21105/joss.02651)
[![CRAN Downloads.](https://cranlogs.r-pkg.org/badges/grand-total/medrxivr)](https://CRAN.R-project.org/package=medrxivr)
<br>
[![R build status](https://github.com/ropensci/medrxivr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/medrxivr/actions)
[![Travis build status](https://travis-ci.com/ropensci/medrxivr.svg?branch=master)](https://travis-ci.com/ropensci/medrxivr)
[![Codecov test coverage](https://codecov.io/gh/ropensci/medrxivr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/medrxivr?branch=master)

<!-- badges: end -->

An increasingly important source of health-related bibliographic content are preprints - preliminary versions of research articles that have yet to undergo peer review. The two preprint repositories most relevant to health-related sciences are [medRxiv](https://www.medrxiv.org/) and [bioRxiv](https://www.biorxiv.org/), both of which are operated by the Cold Spring Harbor Laboratory.

The goal of the `medrxivr` R package is two-fold. In the first instance, it provides programmatic access to the [Cold Spring Harbour Laboratory (CSHL) API](https://api.biorxiv.org/), allowing users to easily download medRxiv and bioRxiv preprint metadata (e.g. title, abstract, publication date, author list, etc) into R. The package also provides access to a maintained static snapshot of the medRxiv repository (see [Data sources](#medrxiv-data)). Secondly, `medrxivr` provides functions to search the downloaded preprint records using regular expressions and Boolean logic, as well as helper functions that allow users to export their search results to a .BIB file for easy import to a reference manager and to download the full-text PDFs of preprints matching their search criteria.

## Installation


To install the stable version of the package from CRAN:

``` {r}
install.packages("medrxivr")
library(medrxivr)
```

Alternatively, to install the development version from GitHub, use the following code:

``` {r}
install.packages("devtools")
devtools::install_github("ropensci/medrxivr")
library(medrxivr)
```

## Data sources

### medRxiv data

`medrixvr` provides two ways to access medRxiv data:  

  - `mx_api_content(server = "medrxiv")` creates a local copy of all data available from the medRxiv API at the time the function is run.
  
``` {r}
# Get a copy of the database from the live medRxiv API endpoint
preprint_data <- mx_api_content()  
```
  
  - `mx_snapshot()` provides access to a static snapshot of the medRxiv database. The snapshot is created each morning at 6am using `mx_api_content()` and is stored as CSV file in the [medrxivr-data repository](https://github.com/mcguinlu/medrxivr-data). This method does not rely on the API (which can become unavailable during peak usage times) and is usually faster (as it reads data from a CSV rather than having to re-extract it from the API). Discrepancies between the most recent static snapshot and the live database can be assessed using `mx_crosscheck()`.
  
``` {r}
# Get a copy of the database from the daily snapshot
preprint_data <- mx_snapshot()  
```
  
The relationship between the two methods for the medRxiv database is summarised in the figure below:
  
``` {r eval = TRUE, echo = FALSE, out.width = "500px", out.height = "400px"}

knitr::include_graphics("vignettes/data_sources.png")
```
  
### bioRxiv data

Only one data source exists for the bioRxiv repository: 

  - `mx_api_content(server = "biorxiv")` creates a local copy of all data available from the bioRxiv API endpoint at the time the function is run. __Note__: due to it's size, downloading a complete copy of the bioRxiv repository in this manner takes a long time (~ 1 hour). 

``` {r}
# Get a copy of the database from the live bioRxiv API endpoint
preprint_data <- mx_api_content(server = "biorxiv")
```

## Performing your search

Once you have created a local copy of either the medRxiv or bioRxiv preprint database, you can pass this object (`preprint_data` in the examples above) to `mx_search()` to search the preprint records using an advanced search strategy.

``` {r, eval = TRUE, message = TRUE}
# Import the medrxiv database
preprint_data <- mx_snapshot()

# Perform a simple search
results <- mx_search(data = preprint_data,
                     query ="dementia")

# Perform an advanced search
topic1  <- c("dementia","vascular","alzheimer's")  # Combined with Boolean OR
topic2  <- c("lipids","statins","cholesterol")     # Combined with Boolean OR
myquery <- list(topic1, topic2)                    # Combined with Boolean AND

results <- mx_search(data = preprint_data,
                     query = myquery)

```

You can also explore which search terms are contributing most to your search by setting `report = TRUE`:

```{r, eval = TRUE, message = TRUE}
results <- mx_search(data = preprint_data,
                     query = myquery,
                     report = TRUE)
```

## Further functionality

### Export records identified by your search to a .BIB file

Pass the results of your search above (the `results` object) to the `mx_export()` to export references for preprints matching your search results to a .BIB file so that they can be easily imported into a reference manager (e.g. Zotero, Mendeley).

```{r, eval = FALSE}
mx_export(data = results,
          file = "mx_search_results.bib")

```

### Download PDFs for records returned by your search

Pass the results of your search above (the `results` object) to the `mx_download()` function to download a copy of the PDF for each record found by your search.

```{r}
mx_download(results,        # Object returned by mx_search(), above
            "pdf/",         # Directory to save PDFs to 
            create = TRUE)  # Create the directory if it doesn't exist

```

## Accessing the raw API data

By default, the `mx_api_*()` functions clean the data returned by the API for use with other `medrxivr` functions. 

To access the raw data returned by the API, the `clean` argument should  set to `FALSE`:

``` {r}
mx_api_content(to_date = "2019-07-01", clean = FALSE)

```

See [this article](https://docs.ropensci.org/medrxivr/articles/medrxiv-api.html#accessing-the-raw-api-data) for more details.

## Detailed guidance

Detailed guidance, including advice on how to design complex search strategies, is available on the [`medrxivr` website.](https://docs.ropensci.org/medrxivr/)

## Linked repositories

See here for the [code used to take the daily snapshot](https://github.com/mcguinlu/medrxivr-data) and [the code that powers the `medrxivr` web app](https://github.com/mcguinlu/medrxivr-app).

## Other tools/packages for working with medRxiv/bioRxiv data

The focus of `medrxivr` is on providing tools to allow users to import and then search medRxiv and bioRxiv data. Below are a list of complementary packages that provide distinct but related functionality when working with medRxiv and bioRxiv data:

* [`rbiorxiv`](https://github.com/nicholasmfraser/rbiorxiv) by [Nicholas Fraser](https://github.com/nicholasmfraser) provides access to the same medRxiv and bioRxiv _content_ data as `medrxivr`, but also provides access to the _usage_ data (e.g. downloads per month) that the Cold Spring Harbour Laboratory API offers. This is useful if you wish to explore, for example, [how the number of PDF downloads from bioRxiv has grown over time.](https://github.com/nicholasmfraser/rbiorxiv#pdf-downloads-over-time)

## Code of conduct

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By
contributing to this project, you agree to abide by its terms.

## Disclaimer

This package and the data it accesses/returns are provided "as is", with no guarantee of accuracy. Please be sure to check the accuracy of the data yourself (and do let me know if you find an issue so I can fix it for everyone!)

---
title: "Interacting with the Cold Spring Harbour Laboratory API"
author: "Luke A McGuinness"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    toc: yes
    toc_depth: 4
vignette: >
  %\VignetteIndexEntry{medrxiv-api}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r,echo=FALSE, message=FALSE, warning=FALSE}
# Delete when done
library(medrxivr)
library(dplyr)

knitr::opts_chunk$set(
  collapse = TRUE,
  eval = TRUE,
  warning = FALSE,
  comment = "#>"
)

```

## Background

The [Cold Spring Harbour Laboratory API](https://api.biorxiv.org/) provides a direct interface to the medRxiv and bioRxiv databases. **However, the API does not allow you to perform searches,** instead providing two endpoints that return either all content between two specified dates or all information held on a particular DOI. 

`medrxivr` provides two convenience functions for importing the data provided by these endpoints in R: `mx_api_content()` and `mx_api_doi()`. The results of either function can then be passed to `mx_search()` for searching.

### By date range (`mx_api_content()`)

The format of this endpoint is https://api.biorxiv.org/details/[server]/[interval]/[cursor] where 'interval' must be two YYYY-MM-DD dates separated by '/'. Where metadata for multiple papers is returned, results are paginated with 100 papers served in a call. The 'cursor' value can be used to iterate through the result. 

`mx_api_content()` automatically moves through the pages for you, capturing all records returned by the endpoint and returning them as an R object. For instance, https://api.biorxiv.org/details/medrxiv/2020-01-01/2020-01-31/0 will output 100 results (if that many remain) within the date range of 2020-01-01 to 2020-01-31 beginning from result 1. To import this into R as a dataframe:

``` {r}
medrxiv_data <- mx_api_content(from_date = "2020-01-01", 
                               to_date = "2020-01-05")


biorxiv_data <- mx_api_content(server = "biorxiv",
                               from_date = "2020-01-01", 
                               to_date = "2020-01-05")

```

### By DOI (`mx_api_doi()`)

https://api.biorxiv.org/details/[server]/[DOI] returns detail for a single manuscript. For instance, https://api.biorxiv.org/details/medrxiv/10.1101/2020.02.25.20021568 will output metadata for the medRxiv paper with DOI 10.1101/2020.02.25.20021568. To import the results from this endpoint into R as a dataframe:

``` {r}
mx_api_doi(doi = "10.1101/2020.02.25.20021568")

```

## Accessing the raw API data

Both functions contain a `clean` argument with is set to `TRUE` by default. This is to ensure that the datasets returned by the `mx_api_*()` functions can immediately be passed to `mx_search()`. However, there may be occasions where this is not required, and so setting this argument to `FALSE` will return the raw data provided by the API endpoints. For example:

``` {r}
mx_api_content(to_date = "2019-07-01", clean = FALSE)
```
---
title: "Get started"
author: "Luke A McGuinness"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    toc: yes
    toc_depth: 4
vignette: >
  %\VignetteIndexEntry{medrxivr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r,echo=FALSE, message=FALSE,warning=FALSE}
# Delete when done
library(medrxivr)
library(dplyr)

knitr::opts_chunk$set(
  collapse = TRUE,
  eval = FALSE,
  comment = "#>"
)

```

An increasingly important source of health-related bibliographic content are preprints - preliminary versions of research articles that have yet to undergo peer review. The two preprint repositories most relevant to health-related sciences are [medRxiv](https://www.medrxiv.org/) and [bioRxiv](https://www.biorxiv.org/), both of which are operated by the Cold Spring Harbor Laboratory.

The goal of the `medrxivr` R package is two-fold. In the first instance, it provides programmatic access to the [Cold Spring Harbour Laboratory (CSHL) API](https://api.biorxiv.org/), allowing users to easily download medRxiv and bioRxiv preprint metadata (e.g. title, abstract, publication date, author list, etc) into R. The package also provides access to a maintained static snapshot of the medRxiv repository (see [Data sources](#medrxiv-data)). Secondly, `medrxivr` provides functions to search the downloaded preprint records using regular expressions and Boolean logic, as well as helper functions that allow users to export their search results to a .BIB file for easy import to a reference manager and to download the full-text PDFs of preprints matching their search criteria.

## Installation

You can install the development version of this package using:

``` {r}
devtools::install_github("mcguinlu/medrxivr")
library(medrxivr)
```

## Data sources

### medRxiv data

`medrixvr` provides two ways to access medRxiv data:  

  - `mx_api_content(server = "medrxiv")` creates a local copy of all data available from the medRxiv API at the time the function is run.
  
``` {r}
# Get a copy of the database from the live medRxiv API endpoint
preprint_data <- mx_api_content()  
```
  
  - `mx_snapshot()` provides access to a static snapshot of the medRxiv database. The snapshot is created each morning at 6am using `mx_api_content()` and is stored as CSV file in the [medrxivr-data repository](https://github.com/mcguinlu/medrxivr-data). This method does not rely on the API (which can become unavailable during peak usage times) and is usually faster (as it reads data from a CSV rather than having to re-extract it from the API). Discrepancies between the most recent static snapshot and the live database can be assessed using `mx_crosscheck()`.
  
``` {r}
# Get a copy of the database from the daily snapshot
preprint_data <- mx_snapshot()  
```
  
The relationship between the two methods for the medRxiv database is summarised in the figure below:
  
``` {r eval = TRUE, echo = FALSE, out.width = "500px", out.height = "400px"}

knitr::include_graphics("data_sources.png")
```
  
### bioRxiv data

Only one data source exists for the bioRxiv repository: 

  - `mx_api_content(server = "biorxiv")` creates a local copy of all data available from the bioRxiv API endpoint at the time the function is run. __Note__: due to it's size, downloading a complete copy of the bioRxiv repository in this manner takes a long time (~ 1 hour). 

``` {r}
# Get a copy of the database from the live bioRxiv API endpoint
preprint_data <- mx_api_content(server = "biorxiv")
```

## Performing your search

Once you have created a local copy of either the medRxiv or bioRxiv preprint database, you can pass this object (`preprint_data` in the examples above) to `mx_search()` to search the preprint records using an advanced search strategy.

``` {r}

# Perform a simple search
results <- mx_search(data = preprint_data,
                     query ="dementia")

# Perform an advanced search
topic1  <- c("dementia","vascular","alzheimer's")  # Combined with Boolean OR
topic2  <- c("lipids","statins","cholesterol")     # Combined with Boolean OR
myquery <- list(topic1, topic2)                    # Combined with Boolean AND

results <- mx_search(data = preprint_data,
                     query = myquery)

```

## Dataset description

The dataset (in this case, `results`) returned by the search function above contains 14 variables: 

```{r, eval = TRUE, echo = FALSE}

mx_variables <-
  data.frame(
    Variable = c(
         "ID"      ,
         "title"   ,
         "abstract",
         "authors" ,
         "date"    ,
         "category",
         "doi"     ,
         "version" ,
         "author_corresponding",
         "author_corresponding_institution",
         "link_page",
         "link_pdf" ,
         "license"  ,
         "published"
    ),
    Description = c(
      "Unique identifier",
      "Preprint title",
      "Preprint abstract",
      "Author list in the format 'LastName, InitalOfFirstName.' (e.g. McGuinness, L.). Authors are seperated by a semi-colon.",
      "Date the preprint was posted, in the format YYYYMMDD.",
      "On submission, medRxiv asks authors to classify their preprint into one of a set number of subject categories.",
      "Preprint Digital Object Identifier.",
      "Preprint version number. As authors can update their preprint at any time, this indicates which version of a given preprint the record refers to.", 
      "Corresponding authors name.",
      "Corresponding author's institution.",
      "Link to preprint webpage. The \"?versioned=TRUE\" is required, as otherwise, the URL will resolve to the most recent version of the article (assuming there is >1 version available).",
      "Link to preprint PDF. This is used by `mx_download()` to download a copy of the PDF for that preprint.",
      "Preprint license",
      "If the preprint was subsequently published in a peer-reviewed journal, this variable contains the DOI of the published version."
    )
  )


knitr::kable(mx_variables, format = "html") %>%
  kableExtra::kable_styling(full_width = F) %>%
  kableExtra::column_spec(1, bold = T, border_right = T) %>%
  kableExtra::column_spec(2, width = "30em")
```

## Export records identified by your search to a .BIB file

`medrxivr` provides a helper function to export your search results to a .BIB file so that they can be easily imported into a reference manager (e.g. Zotero, Mendeley)

```{r, eval = FALSE}

mx_export(data = mx_results,
          file = tempfile(fileext = ".bib"))

```

## Download PDFs for records identified by your search

Pass the results of your search above (`results`) to the `mx_download()` function to download a copy of the PDF for each record. 

```{r, eval = FALSE}

mx_download(results,        # Object returned by mx_search
            tempdir(),      # Temporary directory to save PDFs to 
            create = TRUE)  # Create the directory if it doesn't exist

```


## Further guidance

Please see the *[medrxivr website](https://docs.ropensci.org/medrxivr/index.html)* vignette for extended guidance on developing search strategies and for detailed instructions on interacting with the Cold Springs Harbour API for medRxiv and bioRxiv.
---
title: "Building complex search strategies"
author: "Luke A McGuinness"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    toc: yes
    toc_depth: 4
vignette: >
  %\VignetteIndexEntry{building-complex-search-strategies}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(medrxivr)
library(dplyr)

knitr::opts_chunk$set(
  collapse = TRUE,
  eval = TRUE,
  warning = FALSE,
  comment = "#>"
)
```


## Building your search with Boolean operators

First load the `medrxivr` package:

```{r setup}
library(medrxivr)

```

To find records that contain any of many terms, pass the terms as a vector to the `mx_search()` function, as in the code chunk below. Query terms can include regular expression syntax - see the [section at the end of this document](#regex) on common regular expression that may be useful when searching.

``` {r}
myquery <- c("dementia","vascular","alzheimer's") # Combined with Boolean OR

mx_results <- mx_search(data = mx_snapshot(),     # Use daily snapshot for data
                        query = myquery)

```

To find records relevant to more than one topic domain, create a vector for each topic (note: there is no upper limit on the number of topics your can have) and combine these vectors into a list which is then passed to the `mx_search()` function:

``` {r}
topic1  <- c("dementia","vascular","alzheimer's")  # Combined with Boolean OR
topic2  <- c("lipids","statins","cholesterol")     # Combined with Boolean OR
myquery <- list(topic1, topic2)                    # Combined with Boolean AND

mx_results <- mx_search(data = mx_snapshot(),
                        query = myquery)

```

## Additional filters and options

### Limit search by field

By default, a range of fields (title, abstract, first author, subject, link (which contains DOI)) are searched, but you can limit the search to a subset of these using the `fields` argument:

```{r}

# Limit search to title/abstract
mx_results <- mx_search(data = mx_snapshot(),
                        query = "dementia",
                        fields = c("title","abstract"))

# Search by DOI
mx_results <- mx_search(data = mx_snapshot(),
                        query = "10.1101/2020.01.30.20019836",
                        fields = "link")

```

### Exclude records containing certain terms

Often it is useful to be able to exclude records that contain a certain term that is not relevant to your search. For example, in the search below, we are looking for records related to "dementia" alone by excluding those that mention "mild cognitive impairment":

```{r}
mx_results <- mx_search(data = mx_snapshot(),
                        query = "dementia",
                        NOT = "[Mm]ild cognitive impairment")
```

### Limit by date posted

You can define either/both of the earliest and latest date you wish to include records from. Note: the search is inclusive of both dates specified:

```{r}
mx_results <- mx_search(data = mx_snapshot(),
                        query = "dementia",
                        from_date = "2020-01-01",      # 1st Jan 2020
                        to_date = "2020-01-08")        # 8th Jan 2020
```

### Return multiple versions of a record

_medRxiv_ allows authors to upload a new version of their preprint as often as they like. By default, `medrxivr` only returns the most recent version of the preprint, but if you are interested in exploring how a record changed over time, you can retrieve all versions of the preprint by setting `deduplicate = FALSE` 

```{r}
mx_results <- mx_search(data = mx_snapshot(),
                        query = "10.1101/2020.01.30.20019836",
                        fields = "link",
                        deduplicate = FALSE)
```

## Useful syntax for the systematic reviewer {#regex}

### Capitalisation

__Example regex:__ `[Dd]ementia`  
__Description:__ The search is case sensitive, so this syntax allows you to find both <b>D</b>ementia and <b>d</b>ementia using a single term, rather than having to enter them separately. However, setting the `autocaps` argument of `mx_search()` to `TRUE` will automatically search for both capitalised and uncapitalised versions of your search terms (e.g. with `auto_caps = TRUE` you just need to search for "dementia" to find both <b>D</b>ementia and <b>d</b>ementia - behind the scenes, "dementia" is converted to "[Dd]ementia".

### Wildcard

__Example regex:__ `randomi*ation`  
__Description:__ The wildcard operator "*" defines any single alphanumeric character - in this case, the term will find both randomi<b>s</b>ation and randomi<b>z</b>ation. 


### NEAR

__Example regex:__ `systematic NEAR4 review`  
__Description:__ The "NEAR4" operator defines that up to 4 words can be between <b>systematic</b> and <b>review</b> and the search will still find it. To change how far apart the terms are allowed to be, simply change the number following NEAR (e.g. to find terms that are only one word apart, the syntax would be `systematic NEAR1 review`). **Please note that the search is directional, in that the example term here will find "systematic methods for the review", but will not find "the review was systematic".**

### Word limits

__Example regex:__ `\\bNCOV\\b`  
__Description:__ Sometimes it is useful to be able to define the start and end of terms. For example, if you were searching for NCOV-19, simply using `ncov` as your search term would also return records containing u<b>ncov</b>ered. Using `\\b` allows you to define where the term beings and ends, thus excluding false positive matches.

### Example using these regexes

To find records that contain "Mendelian" within 4 words of "randomisation" (with varying capitalisation of "Mendelian" and UK/US spellings of "randomisation"), the following syntax is correct:

``` {r}
mx_results <- mx_search(data = mx_snapshot(),
                        query = "mendelian NEAR4 randomi*ation", 
                        auto_caps = TRUE)

```

### Regex tester

To check whether your search term will find what you expect it to, there is a useful [regex tester](https://spannbaueradam.shinyapps.io/r_regex_tester/), designed by [Adam Spannbauer](https://adamspannbauer.github.io/2018/01/16/r-regex-tester-shiny-app/).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_syntax.R
\name{mx_caps}
\alias{mx_caps}
\title{Search term wrapper that allows for different capitalisation of term}
\usage{
mx_caps(x)
}
\arguments{
\item{x}{Search term to be formatted}
}
\value{
The input string is return, but with each non-space character
  repeated in lower- and upper-case, and enclosed in square brackets. For
  example, mx_caps("ncov") returns "[Nn][Cc][Oo][Vv]"
}
\description{
Inspired by the varying capitalisation of "NCOV" during the
  coronavirus pandemic (e.g. ncov, nCoV, NCOV, nCOV), this function allows
  for all possible configurations of lower- and upper-case letters in your
  search term.
}
\examples{
\donttest{

query <- c("coronavirus", mx_caps("ncov"))

mx_search(mx_snapshot("6c4056d2cccd6031d92ee4269b1785c6ec4d555b"), query)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{api_to_df}
\alias{api_to_df}
\title{Convert API data to data frame}
\usage{
api_to_df(url)
}
\arguments{
\item{url}{API endpoint from which to extract and format data}
}
\value{
Raw API data in a dataframe
}
\description{
Convert API data to data frame
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_api.R
\name{mx_api_content}
\alias{mx_api_content}
\title{Access medRxiv/bioRxiv data via the Cold Spring Harbour Laboratory API}
\usage{
mx_api_content(
  from_date = "2013-01-01",
  to_date = as.character(Sys.Date()),
  clean = TRUE,
  server = "medrxiv",
  include_info = FALSE
)
}
\arguments{
\item{from_date}{Earliest date of interest, written as "YYYY-MM-DD". Defaults
to 1st Jan 2013 ("2013-01-01"), ~6 months prior to earliest preprint
registration date.}

\item{to_date}{Latest date of interest, written as "YYYY-MM-DD". Defaults to
current date.}

\item{clean}{Logical, defaulting to TRUE, indicating whether to clean the
data returned by the API. If TRUE, variables containing absolute paths to
the preprints web-page ("link_page") and PDF ("link_pdf") are generated
from the "server", "DOI", and "version" variables returned by the API. The
"title", "abstract" and "authors" variables are converted to title case.
Finally, the "type" and "server" variables are dropped.}

\item{server}{Specify the server you wish to use: "medrxiv" (default) or
"biorxiv"}

\item{include_info}{Logical, indicating whether to include variables
containing information returned by the API (e.g. API status, cursor number,
total count of papers, etc). Default is FALSE.}
}
\value{
Dataframe with 1 record per row
}
\description{
Provides programmatic access to all preprints available through
  the Cold Spring Harbour Laboratory API, which serves both the medRxiv and
  bioRxiv preprint repositories.
}
\examples{
if(interactive()){
mx_data <- mx_api_content(from_date = "2020-01-01",
to_date = "2020-01-07")
}
}
\seealso{
Other data-source: 
\code{\link{mx_api_doi}()},
\code{\link{mx_snapshot}()}
}
\concept{data-source}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{skip_if_api_message}
\alias{skip_if_api_message}
\title{Skips API tests if API isn't working correctly}
\usage{
skip_if_api_message()
}
\description{
Skips API tests if API isn't working correctly
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_export.R
\name{mx_export}
\alias{mx_export}
\title{Export references for preprints returning by a search to a .bib file}
\usage{
mx_export(data, file = "medrxiv_export.bib")
}
\arguments{
\item{data}{Dataframe returned by mx_search() or mx_api_*() functions}

\item{file}{File location to save to. Must have the .bib file extension}
}
\value{
Exports a formatted .BIB file, for import into a reference manager
}
\description{
Export references for preprints returning by a search to a .bib file
}
\examples{
\donttest{
mx_results <- mx_search(mx_snapshot(), query = "brain")
mx_export(mx_results, tempfile(fileext = ".bib"))
}

}
\seealso{
Other helper: 
\code{\link{mx_crosscheck}()},
\code{\link{mx_download}()}
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medrxivr.R
\docType{package}
\name{medrxivr}
\alias{medrxivr}
\title{medrxivr: Accessing medRxiv and bioRxiv preprint data from R}
\description{
The medrxivr package enables users to access data on preprints in the medRxiv
and bioRxiv preprints repositories, both of which are run by the Cold Spring
Harbour Laboratory.It also provides functions to search the preprint data,
export it to a .bib file, and download the PDFs associated with specified
records.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_search.R
\name{mx_search}
\alias{mx_search}
\title{Search preprint data}
\usage{
mx_search(
  data = NULL,
  query = NULL,
  fields = c("title", "abstract", "authors", "category", "doi"),
  from_date = NULL,
  to_date = NULL,
  auto_caps = FALSE,
  NOT = "",
  deduplicate = TRUE,
  report = FALSE
)
}
\arguments{
\item{data}{The preprint dataset that is to be searched, created either using
mx_api_content() or mx_snapshot()}

\item{query}{Character string, vector or list}

\item{fields}{Fields of the database to search - default is Title, Abstract,
Authors, Category, and DOI.}

\item{from_date}{Defines earliest date of interest. Written in the format
"YYYY-MM-DD". Note, records published on the date specified will also be
returned.}

\item{to_date}{Defines latest date of interest. Written in the format
"YYYY-MM-DD". Note, records published on the date specified will also be
returned.}

\item{auto_caps}{As the search is case sensitive, this logical specifies
whether the search should automatically allow for differing capitalisation
of search terms. For example, when TRUE, a search for "dementia" would find
both "dementia" but also "Dementia". Note, that if your term is multi-word
(e.g. "systematic review"), only the first word is automatically
capitalised (e.g your search will find both "systematic review" and
"Systematic review" but won't find "Systematic Review". Note that this
option will format terms in the query and NOT arguments (if applicable).}

\item{NOT}{Vector of regular expressions to exclude from the search. Default
is "".}

\item{deduplicate}{Logical. Only return the most recent version of a record.
Default is TRUE.}

\item{report}{Logical. Run mx_reporter. Default is FALSE.}
}
\description{
Search preprint data
}
\examples{
\donttest{
# Using the daily snapshot
mx_results <- mx_search(data = mx_snapshot(), query = "dementia")
}
}
\seealso{
Other main: 
\code{\link{mx_reporter}()},
\code{\link{print_full_results}()},
\code{\link{run_search}()}
}
\concept{main}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_crosscheck.R
\name{mx_crosscheck}
\alias{mx_crosscheck}
\title{Check how up-to-date the maintained medRxiv snapshot is}
\usage{
mx_crosscheck()
}
\description{
Provides information on how up-to-date the maintained medRxiv
  snapshot provided by `mx_snapshot()` is by checking whether there have been
  any records added to, or updated in, the medRxiv repository since the last
  snapshot was taken.
}
\examples{
\donttest{
mx_crosscheck()
}
}
\seealso{
Other helper: 
\code{\link{mx_download}()},
\code{\link{mx_export}()}
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{api_link}
\alias{api_link}
\title{Create link for API}
\usage{
api_link(...)
}
\arguments{
\item{...}{Arguments to specify the path to the API endpoint}
}
\value{
Formatted link to API endpoint
}
\description{
Create link for API
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_info.R
\name{mx_info}
\alias{mx_info}
\title{Provide information on the medRxiv snapshot used to perform the search}
\usage{
mx_info(commit = "master")
}
\arguments{
\item{commit}{Commit hash for the snapshot, taken from
https://github.com/mcguinlu/medrxivr-data. Defaults to "master", which will
return info on the most recent snapshot.}
}
\value{
Message with snapshot details
}
\description{
Provide information on the medRxiv snapshot used to perform the search
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_search.R
\name{print_full_results}
\alias{print_full_results}
\title{Search for terms in the dataset}
\usage{
print_full_results(num_results, deduplicate)
}
\arguments{
\item{num_results}{number of searched terms returned}

\item{deduplicate}{Logical. Only return the most recent version of a record.
Default is TRUE.}
}
\description{
Search for terms in the dataset
}
\seealso{
Other main: 
\code{\link{mx_reporter}()},
\code{\link{mx_search}()},
\code{\link{run_search}()}
}
\concept{main}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_snapshot.R
\name{mx_snapshot}
\alias{mx_snapshot}
\title{Access a static snapshot of the medRxiv repository}
\usage{
mx_snapshot(commit = "master")
}
\arguments{
\item{commit}{Commit hash for the snapshot, taken from
https://github.com/mcguinlu/medrxivr-data. Allows for reproducible
searching by specifying the exact snapshot used to perform the searches.
Defaults to "master", which will return the most recent snapshot.}
}
\value{
Formatted dataframe
}
\description{
[Available for medRxiv only] Rather than downloading a copy of
  the medRxiv database from the API, which can become unavailable at peak
  usage times, this allows users to import a maintained static snapshot of
  the medRxiv repository.
}
\examples{
\donttest{
mx_data <- mx_snapshot()
}

}
\seealso{
Other data-source: 
\code{\link{mx_api_content}()},
\code{\link{mx_api_doi}()}
}
\concept{data-source}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_syntax.R
\name{fix_caps}
\alias{fix_caps}
\title{Allow for capitalisation of search terms}
\usage{
fix_caps(x)
}
\arguments{
\item{x}{Search query to be formatted. Note, any search term already
containing a square bracket will not be reformatted to preserve
user-defined regexes.}
}
\value{
The same list or vector search terms, but with proper regular
  expression syntax to allow for capitalisation of the first letter of each
  term.
}
\description{
Allow for capitalisation of search terms
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_search.R
\name{mx_reporter}
\alias{mx_reporter}
\title{Search and print output for individual search items}
\usage{
mx_reporter(mx_data, num_results, query, fields, deduplicate, NOT)
}
\arguments{
\item{mx_data}{The mx_dataset filtered for the date limits}

\item{num_results}{The number of results returned by the overall search}

\item{query}{Character string, vector or list}

\item{fields}{Fields of the database to search - default is Title, Abstract,
Authors, Category, and DOI.}

\item{deduplicate}{Logical. Only return the most recent version of a record.
Default is TRUE.}

\item{NOT}{Vector of regular expressions to exclude from the search. Default
is "".}
}
\description{
Search and print output for individual search items
}
\seealso{
Other main: 
\code{\link{mx_search}()},
\code{\link{print_full_results}()},
\code{\link{run_search}()}
}
\concept{main}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_syntax.R
\name{fix_wildcard}
\alias{fix_wildcard}
\title{Replace user-friendly 'wildcard' operator with appropriate regex syntax}
\usage{
fix_wildcard(x)
}
\arguments{
\item{x}{Search query to be reformatted}
}
\description{
Replace user-friendly 'wildcard' operator with appropriate regex syntax
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_download.R
\name{mx_download}
\alias{mx_download}
\title{Download PDF's of preprints returned by a search}
\usage{
mx_download(
  mx_results,
  directory,
  create = TRUE,
  name = c("ID", "DOI"),
  print_update = 10
)
}
\arguments{
\item{mx_results}{Vector containing the links to the medRxiv PDFs}

\item{directory}{The location you want to download the PDF's to}

\item{create}{TRUE or FALSE. If TRUE, creates the directory if it doesn't
exist}

\item{name}{How to name the downloaded PDF. By default, both the ID number of
the record and the DOI are used.}

\item{print_update}{How frequently to print an update}
}
\description{
Download PDF's of all the papers in your search results
}
\examples{
\donttest{
mx_results <- mx_search(mx_snapshot(), query = "10.1101/2020.02.25.20021568")
mx_download(mx_results, directory=tempdir())
}
}
\seealso{
Other helper: 
\code{\link{mx_crosscheck}()},
\code{\link{mx_export}()}
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_api.R
\name{mx_api_doi}
\alias{mx_api_doi}
\title{Access data on a single medRxiv/bioRxiv record via the Cold Spring Harbour
Laboratory API}
\usage{
mx_api_doi(doi, server = "medrxiv", clean = TRUE)
}
\arguments{
\item{doi}{Digital object identifier of the preprint you wish to retrieve
data on.}

\item{server}{Specify the server you wish to use: "medrxiv" (default) or
"biorxiv"}

\item{clean}{Logical, defaulting to TRUE, indicating whether to clean the
data returned by the API. If TRUE, variables containing absolute paths to
the preprints web-page ("link_page") and PDF ("link_pdf") are generated
from the "server", "DOI", and "version" variables returned by the API. The
"title", "abstract" and "authors" variables are converted to title case.
Finally, the "type" and "server" variables are dropped.}
}
\value{
Dataframe containing details on the preprint identified by the DOI.
}
\description{
Provides programmatic access to data on a single preprint
  identified by a unique Digital Object Identifier (DOI).
}
\examples{
if(interactive()){
mx_data <- mx_api_doi("10.1101/2020.02.25.20021568")
}
}
\seealso{
Other data-source: 
\code{\link{mx_api_content}()},
\code{\link{mx_snapshot}()}
}
\concept{data-source}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{clean_api_df}
\alias{clean_api_df}
\title{Helper script to clean data from API to make it compatible with mx_search()}
\usage{
clean_api_df(df)
}
\arguments{
\item{df}{Raw dataframe from API}
}
\value{
Cleaned dataframe
}
\description{
Helper script to clean data from API to make it compatible with mx_search()
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_syntax.R
\name{fix_near}
\alias{fix_near}
\title{Replace user-friendly 'NEAR' operator with appropriate regex syntax}
\usage{
fix_near(x)
}
\arguments{
\item{x}{Search query to be reformatted}
}
\description{
Replace user-friendly 'NEAR' operator with appropriate regex syntax
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{internet_check}
\alias{internet_check}
\title{Checks whether the user has internet, and returns a helpful message it not.}
\usage{
internet_check()
}
\value{
Informative error if not connected to the internet
}
\description{
Checks whether the user has internet, and returns a helpful message it not.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mx_search.R
\name{run_search}
\alias{run_search}
\title{Search for terms in the dataset}
\usage{
run_search(mx_data, query, fields, deduplicate, NOT = "")
}
\arguments{
\item{mx_data}{The mx_dataset filtered for the date limits}

\item{query}{Character string, vector or list}

\item{fields}{Fields of the database to search - default is Title, Abstract,
Authors, Category, and DOI.}

\item{deduplicate}{Logical. Only return the most recent version of a record.
Default is TRUE.}

\item{NOT}{Vector of regular expressions to exclude from the search. Default
is NULL.}
}
\description{
Search for terms in the dataset
}
\seealso{
Other main: 
\code{\link{mx_reporter}()},
\code{\link{mx_search}()},
\code{\link{print_full_results}()}
}
\concept{main}
