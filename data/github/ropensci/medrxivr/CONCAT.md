
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
