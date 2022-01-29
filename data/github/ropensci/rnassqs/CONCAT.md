---
title: "rnassqs: An `R` package to access agricultural data via the USDA 
  National Agricultural Statistics Service (USDA-NASS) 'Quick Stats' API"
authors:
  - affiliation: 1
    email: potter.nicholas@gmail.com
    name: Nicholas A. Potter
    orcid: 0000-0002-3410-3732
date: "09 June 2019"
output: pdf_document
bibliography: paper.bib
tags:
  - R
  - API
  - reproducibility
  - agriculture
  - economics
affiliations:
  - index: 1
    name: Washington State University
---

# Summary

The [rnassqs](https://github.com/ropensci/rnassqs) `R` package
provides a simple interface for accessing the United States Department of
Agriculture National Agricultural Statistics Service (USDA-NASS)
'[Quick Stats](https://quickstats.nass.usda.gov/)' API. The
core functionality allows the user to query agricultural data from 'Quick Stats'
in a reproducible and automated way. The primary benefit of `rnassqs` is
that users need not download data through repeated use of the Quick Stats
point-and-click interface, which reduces the chance of errors, eliminates the
need for repeated manual downloads of new data over time or space, and allows
for automated updates of web applications that rely on new releases of USDA-NASS
data over time.

`rnassqs` manages API authentication by setting a system environmental variable
for the duration of the `R` session. Convenience functions facilitate querying
common data. Users can also use `rnassqs` to query the list of data parameters
and available values for a given parameter (for example, to see the commodities
available in a particular county and year). The query requests data as a JSON
object and parses that data into a `data.frame` object. 


# About USDA-NASS data and 'Quick Stats'

'Quick Stats' is a web interface to access data produced by USDA-NASS. The data 
comes primarily from the Census of Agriculture, but also includes
data from USDA-NASS surveys on a wide range of topics. The Census of Agriculture
is conducted every five years in years ending in '2' and '7'. The earliest year
available on Quick Stats is 1997. Surveys have different collection periods,
but most are collected annually. Some specific data such as the average value of
agricultural land and buildings by state is reported as early as 1850.

Aggregate data from the census and surveys is released primarily at the 
national, state, and county level, though some data may be released for 
congressional districts, watersheds, and zip codes. It includes a range of data
classified under five sectors: Animals & Products, Crops, Demographics, 
Economics, and Environmental. Examples of data available in these sectors 
include counts of farms, farm operators, acres of cropland, farm sales, farm 
expenses, and crop yields to name a few.


# Benefits of `rnassqs` over 'Quick Stats'

'Quick Stats' provides a number of selection fields in which the user can select
categories of data. Each selection causes the options in other fields to update
to reflect available options available options based on other selections that
the user has made. This makes the selection of multiple different variables an
at times frustrating process. In addition, data requests are limited to 50,000
records. If a user wants to access more records they must manually subset their
data request. For example, requesting crop yields by county for all counties
and census years for 1997 to 2017 requires either downloading each state's counties
separately or, where possible, downloading all counties for each year separately.
The Quick Stats interface works well for quick access or a single use. However,
there are several cases in which the 'Quick Stats' interface is not ideal:

- Requests for a combination of measures, years, and 
  geographies that exceed 50,000 records.
- Requests for newly released data that are identical to previous requests.
- Requests that are reproducible.

`rnassqs` addresses each of these issues by making the 'Quick Stats' API 
accessible with `R` code. This allows the user to loop over a series of requests
to address the first issue, to execute (perhaps automated on a schedule) a data 
request repeatedly to access new data with the same query to address the 
second, and to make code available that allows others to reproducibly access 
the same data to address the third.

For example, there are currently there are currently 332,125 records of crop 
yields in all U.S. counties from 2000 to 2018. Accessing this data through 
'Quick Stats' would require manually selecting either a set of years or a set 
of states to reduce each request to less than 50,000 records and then 
aggregating that data. With `rnassqs` this can be done with:


```r
# Access yields for all counties and all crops
params <- list(sector_desc = "CROPS",
               group_desc = c("FIELD CROPS", "FRUIT & TREE NUTS", 
                              "HORTICULTURE", "VEGETABLES"),
               statisticcat_desc = "YIELD", 
               agg_level_desc = "COUNTY")

# Get all years from 2000 to 2018 in a list of data.frames
data_list <- lapply(2000:2018, function(yr) { 
  params$year <- yr
  rnassqs::nassqs(params, url_only = TRUE)
})

# Aggregate the list of data.frames into a single data.frame
d <- do.call("rbind", data_list)
```

This results in significant time savings, increases the reproducibility of
the data, and allows for easy updating of the request when a new year is made
available.


# Alternatives to `rnassqs`

USDA-NASS also provides FTP access to text data files^[Available at: 
[ftp://ftp.nass.usda.gov/quickstats/](ftp://ftp.nass.usda.gov/quickstats/)]. 
By accessing the data via FTP users can avoid using the selection interface of
'Quick Stats' and avoid limitations on the number of records per request, but do 
not resolve issues of automated repeated requests or of making data requests 
reproducible. Other `R` packages have also been released since the development
of `rnassqs` that provide access to Quick Stats data, though somewhat differently,
include `usdarnass` and `tidyUSDA`, which uses US Census maps to create maps
of Quick Stats data.


# Acknowledgements

Thank you to Jonathan Adams, Julia Piaskowski, and Joseph Stachelek for code
contributions. A huge thanks to reviewers Adam H Sparks and Neal Richardson
for their thoughtful feedback to improve the package and documentation.
### Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment include:

- Using welcoming and inclusive language
- Being respectful of differing viewpoints and experiences
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

- The use of sexualized language or imagery and unwelcome sexual attention or advances
- Trolling, insulting/derogatory comments, and personal or political attacks
- Public or private harassment
- Publishing others’ private information, such as a physical or electronic address, without explicit permission
- Other conduct which could reasonably be considered inappropriate in a professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

### Scope

This Code of Conduct applies within all project spaces, and it also applies when an individual is representing the project or its community in public spaces. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at [potter.nicholas@gmail.com]. All complaints will be reviewed and investigated and will result in a response that is deemed necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project’s leadership.

### Attribution

This Code of Conduct is adapted from the Contributor Covenant, version 1.4, available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see https://www.contributor-covenant.org/faq<!-- README.md is generated from README.Rmd. Please edit that file -->

<table class="table">

<thead>

<tr class="header">

<th align="left">

rnassqs

</th>

<th align="left">

Usage

</th>

<th align="left">

Release

</th>

<th align="left">

Development

</th>

</tr>

</thead>

<tbody>

<tr class="odd">

<td rowspan="5">

<a href="https://docs.ropensci.org/rnassqs"><img src="man/figures/logo.png" alt="rnassqs" align="right" height="139"></a>

<p style="font-size:xx-small;">

(Wheat image from
<a href="https://www.flickr.com/photos/53018729@N00/2669034542">here</a>.)

</p>

</td>

<td align="left">

<a href="http://choosealicense.com/licenses/mit/"><img src="https://img.shields.io/github/license/mashape/apistatus.svg" alt="License"></a>

</td>

<td align="left">

<a href="https://cran.r-project.org/package=rnassqs"><img src="http://www.r-pkg.org/badges/version-last-release/rnassqs" alt="CRAN"></a>

</td>

<td align="left">

<a href="https://github.com/ropensci/rnassqs/commits/master"><img src="https://img.shields.io/badge/last%20change-2020--02--13-brightgreen.svg" alt="Last Change"></a>

</td>

</tr>

<tr class="even">

<td align="left">

<a href="https://CRAN.R-project.org/package=rnassqs"><img src="https://cranlogs.r-pkg.org/badges/rnassqs" alt="downloads"></a>

</td>

<td align="left">

<a href="https://zenodo.org/badge/latestdoi/37335585"><img src="https://zenodo.org/badge/37335585.svg" alt="Zenodo"></a>

</td>

<td align="left">

<a href="https://travis-ci.org/ropensci/rnassqs"><img src="https://travis-ci.org/ropensci/rnassqs.svg?branch=master" alt="Build Status"></a>

</td>

</tr>

<tr class="odd">

<td align="left">

</td>

<td align="left">

<a href="https://github.com/ropensci/onboarding/issues/298" alt="rOpensci reviewed!"><img src="https://badges.ropensci.org/298_status.svg"></a>

</td>

<td align="left">

<a href="https://codecov.io/gh/ropensci/rnassqs"><img src="https://codecov.io/gh/ropensci/rnassqs/branch/master/graph/badge.svg" alt="Coverage Status"></a>

</td>

</tr>

<tr class="even">

<td align="left">

</td>

<td align="left">

<a href="https://orcid.org/0000-0002-3410-3732"><img src="https://img.shields.io/badge/ORCiD-0000--0002--3410--3732-green.svg" alt="ORCID"></a>

</td>

<td align="left">

<a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg" alt="Project Status: Active – The project has reached a stable, usable state and is being actively developed." /></a>

</td>

</tr>

<tr class="even">

<td align="left">

</td>

<td align="left">

<a style="border-width:0" href="https://doi.org/10.21105/joss.01880">
<img src="https://joss.theoj.org/papers/10.21105/joss.01880/status.svg" alt="DOI:10.21105/joss.01880" >
</a>

</td>

<td align="left">

<a href="https://www.tidyverse.org/lifecycle/#maturing"><img src="https://img.shields.io/badge/lifecycle-maturing-blue.svg" alt="Project Status: Maturing." /></a>

</td>

</tr>

</tbody>

</table>

<br>

**As required by the NASS Terms of Use: This product uses the NASS API
but is not endorsed or certified by NASS.**

## rnassqs (R NASS Quick Stats)

`rnassqs` allows users to access the USDA’s National Agricultural
Statistics Service (NASS) ‘Quick Stats’ data through their API. It is
simple and easy to use, and provides some functions to help navigate the
bewildering complexity of some Quick Stats data.

For docs and code examples, visit the package web page here:
<https://docs.ropensci.org/rnassqs/>.

## Installing

Install the package via `devtools` or CRAN:

``` r
    # Via devtools
    library(devtools)
    install_github('ropensci/rnassqs')
    
    # Via CRAN
    install.packages("rnassqs")
```

## API Key

To use the NASS Quick Stats API you need an [API
key](http://quickstats.nass.usda.gov/api). The API key should in general
not be included in scripts. One way of making the key available without
defining it in a script is by setting it in your `.Renviron` file, which
is usually located in your home directory. If you are an `rstudio` user,
you can use `usethis::edit_r_environ()` to open your `.Renviron` file
and add a line that looks like:

``` r
    NASSQS_TOKEN="<your api key here>"
```

Alternatively, you can set it explicitly in the console with
`nassqs_auth(key = <your api key>)`. This will set the environmental
variable NASSQS\_TOKEN, which is used to access the API. You can also
set this directly with `Sys.setenv("NASSQS_TOKEN" = <your api key>)`.

## Usage

See the examples in [inst/examples](inst/examples) for quick recipes to
download data.

The primary function is `nassqs()`, with which you can make any query of
variables. For example, to mirror the request that is on the [NASS API
documentation](http://quickstats.nass.usda.gov/api), you can use:

``` r
    library(rnassqs)
    
    # You must set your api key before requesting data
    nassqs_auth(key = <your api key>)
    
    # Parameters to query on and data call
    params <- list(commodity_desc = "CORN", year__GE = 2012, state_alpha = "VA")
    d <- nassqs(params)
```

Parameters **do not** need to be capitalized, and also do not need to be
in a list format. The following works just as well:

``` r
    d <- nassqs(commodity_desc = "corn", year__GE = 2012, state_alpha = "va")
```

You can request data for multiple values of the same parameter by using
a simple list as follows:

``` r
    params <- list(commodity_desc = "CORN", year__GE = 2012, state_alpha = c("VA", "WA"))
    d <- nassqs(params)
```

NASS does not allow GET requests that pull more than 50,000 records in
one request. The function will inform you if you try to do that. It will
also inform you if you’ve requested a set of parameters for which there
are no records.

Other useful functions include:

``` r
    # returns a set of unnique values for the parameter "STATISTICCAT_DESC"
    nassqs_param_values("statisticcat_desc")
    
    # returns a count of the number of records for a given query
    nassqs_record_count(params=params)
    
    # Get yields specifically
    # Equivalent to including "'statisticat_desc' = 'YIELD'" in your parameter list. 
    nassqs_yields(params)
    
    # Get acres specifically
    # Equivalent to including all "AREA" values in statisticcat_desc
    nassqs_acres(params)
    
    # Specifies just "AREA HARVESTED" values of statisticcat_desc
    nassqs_acres(params, area = "AREA HARVESTED")
```

### Handling inequalities and operators other than “=”

The NASS API handles other operators by modifying the variable name. The
API can accept the following modifications:

  - \_\_LE: \<=
  - \_\_LT: \<
  - \_\_GT: \>
  - \_\_GE: \>=
  - \_\_LIKE: like
  - \_\_NOT\_LIKE: not like
  - \_\_NE: not equal

For example, to request corn yields in Virginia and Pennsylvania for all
years since 2000, you would use something like:

``` r
    params <- list(commodity_desc = "CORN", 
                  year__GE = 2000, 
                  state_alpha = c("VA", "PA"), 
                  statisticcat_desc = "YIELD")
    df <- nassqs(params) #returns data as a data frame.
```

See the
[vignette](https://docs.ropensci.org/rnassqs/articles/rnassqs.html) for
more examples and details on usage.

## Contributing

Contributions are more than welcome, and there are several ways to
contribute:

  - Examples: More examples are always helpful. If you use `rnassqs` to
    query data from ‘Quick Stats’ and would like to contribute your
    query, consider submitting a pull request adding your query as a
    file in
    [inst/examples/](https://github.com/ropensci/rnassqs/tree/master/inst/examples).
  - File an issue: If there is functionality you’d like to see added or
    something that is confusing, consider [creating an
    issue](https://github.com/ropensci/rnassqs/issues/new). The best
    issue contains an example of the problem or feature. Consider the
    excellent package [reprex](https://github.com/tidyverse/reprex) in
    creating a reproducible example.
  - Contributing documentation: Clarifying and expanding the
    documentation is always appreciated, especially if you find an area
    that is lacking and would like to improve it. `rnassqs` uses
    roxygen2, which means the documentation is at the top of each
    function definition. Please submit any improvements as a pull
    request.
  - Contributing code: if you see something that needs improving and
    you’d like to make the changes, contributed code is very welcome.
    Begin by filing a new issue to discuss the proposed change, and then
    submit a pull request to address the issue. `rnassqs` follows the
    style outlined in Hadley Wickham’s [R
    Packages](http://r-pkgs.had.co.nz/style.html). Following this style
    makes the pull request and review go more smoothly.

## Alternatives

In June 2019 the `usdarnass` package was released on
[CRAN](https://cran.r-project.org/package=usdarnass) and is also
available to install via [github](https://github.com/rdinter/usdarnass).
`usdarnass` has similar functionality to this package.

NASS also provides a daily tarred and gzipped file of their entire
dataset. At the time of writing it is approaching 1 GB. You can download
that file via their [FTP site](ftp://ftp.nass.usda.gov/quickstats).

The FTP link also contains builds for: NASS census (every 5 years ending
with 2 and 7), or data for one of their specific sectors (CROPS,
ECONOMICS, ANIMALS & PRODUCTS). At the time of this writing, specific
files for the ENVIRONMENTAL and DEMOGRAPHICS sectors are not available.

### Acknowledgements

Thank you to rOpensci reviewers Adam Sparks and Neal Richardson and
editor Lincoln Mullen, for their fantastic feedback and assistance. User
feedback and use case contributions have been a huge help to make
`rnassqs` more accessible and user-friendly. More use cases or feature
requests are always welcome\!

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# rnassqs 0.6.0

- `nassqs_record_count()` now validates parameters. 
- add defined parameters for common parameters in function call (thanks rdinter!).
- add option to see valid parameter values given a query of other values in `get_param_values()`.

# rnassqs 0.5.0

- approval for rOpensci inclusion!
- additional testing to improve code coveral by @nealrichardson
- small changes for rOpensci review process
- switch to rOpensci repository
- Change in syntax to allow for lower case query parameter values
- Change in syntax to allow for specifying each parameter as a separate function argument rather than as a single list (in addition to specifying a single list)
- Create package website with pkgdown
- Standardize code style in package code, examples, and vignette
- Simplify authentication
- Expanded test coverage with use of httptest::with_mock_api()
- Better clarification in documentation and documentation examples
- Improved README and vignette

# rnassqs 0.4.0.9000

- Development version

# rnassqs 0.4.0

- Add automated unit tests that work locally and others that work on CRAN.
- Improve documentation for core functions.
- Add parsing for CSV formatted data.
- Improve authentication.
- Simplify function calls to eliminate redundant calls.
- Add working examples and tests.
- fix name error in the function `nassqs_params_values` to `nassqs_param_values`

# rnassqs 0.3.0

- Prepare package for CRAN submission.
- Vignettes and README.md are up to date with respect to current functions.
- Fix tests.
- Minor spelling fixes contributed by Julia Piaskowski <@jpiaskowski>
- Remove test code that couldn't be run due to API needing authentication.
### Contributing

Contributions are more than welcome. Please begin by reading our [code of conduct](CONDUCT.md). There are several ways to contribute:

- Examples: More examples are always helpful. If you use `rnassqs` to query data from 'Quick Stats' and would like to contribute your query, consider submitting a pull request adding your query as a file in [inst/examples/](https://github.com/ropensci/rnassqs/tree/master/inst/examples).
- File an issue: If there is functionality you'd like to see added or something that is confusing, consider [creating an issue](https://github.com/ropensci/rnassqs/issues/new). The best issue contains an example of the problem or feature. Consider the excellent package [reprex](https://github.com/tidyverse/reprex) in creating a reproducible example.
- Contributing documentation: Clarifying and expanding the documentation is always appreciated, especially if you find an area that is lacking and would like to improve it. `rnassqs` uses roxygen2, which means the documentation is at the top of each function definition. Please submit any improvements as a pull request.
- Contributing code: if you see something that needs improving and you'd like to make the changes, contributed code is very welcome. Begin by filing a new issue to discuss the proposed change, and then submit a pull request to address the issue. `rnassqs` follows the style outlined in Hadley Wickham's [R Packages](http://r-pkgs.had.co.nz/style.html). Following this style makes the pull request and review go more smoothly.
## Test environments
* local Arch Linux Install, kernel 4.20.3-arch1-1-ARCH, R 3.6.1
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.3
* win-builder (devel and release)
* r-hub builder

## R CMD check results
There were no ERRORs, no WARNINGs, and no NOTEs

## Downstream dependencies
There are no downstream dependencies
