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
---
output:
  md_document:
    variant: gfm
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


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
<p style="font-size:xx-small;">(Wheat image from <a href="https://www.flickr.com/photos/53018729@N00/2669034542">here</a>.)</p>
</td>
<td align="left">
<a href="http://choosealicense.com/licenses/mit/"><img src="https://img.shields.io/github/license/mashape/apistatus.svg" alt="License"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/package=rnassqs"><img src="http://www.r-pkg.org/badges/version-last-release/rnassqs" alt="CRAN"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/rnassqs/commits/master"><img src="https://img.shields.io/badge/last%20change-`r gsub('-', '--', Sys.Date())`-brightgreen.svg" alt="Last Change"></a>
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


__As required by the NASS Terms of Use: This product uses the NASS API but is not endorsed or certified by NASS.__


## rnassqs (R NASS Quick Stats)

`rnassqs` allows users to access the USDA's National Agricultural Statistics Service (NASS) 'Quick Stats' data through their API. It is simple and easy to use, and provides some functions to help navigate the bewildering complexity of some Quick Stats data.

For docs and code examples, visit the package web page here: [https://docs.ropensci.org/rnassqs/](https://docs.ropensci.org/rnassqs/).

## Installing

Install the package via `devtools` or CRAN:

```{r eval=FALSE}
    # Via devtools
    library(devtools)
    install_github('ropensci/rnassqs')
    
    # Via CRAN
    install.packages("rnassqs")
```

## API Key

To use the NASS Quick Stats API you need an [API key](http://quickstats.nass.usda.gov/api). The API key should in general not be included in scripts. One way of making the key available without defining it in a script is by setting it in your `.Renviron` file, which is usually located in your home directory. If you are an `rstudio` user, you can use `usethis::edit_r_environ()` to open your `.Renviron` file and add a line that looks like:

```{r eval=FALSE}
    NASSQS_TOKEN="<your api key here>"
```

Alternatively, you can set it explicitly in the console with `nassqs_auth(key = <your api key>)`. This will set the environmental variable NASSQS_TOKEN, which is used to access the API. You can also set this directly with `Sys.setenv("NASSQS_TOKEN" = <your api key>)`. 

## Usage

See the examples in [inst/examples](inst/examples) for quick recipes to download data.

The primary function is `nassqs()`, with which you can make any query of variables. 
For example, to mirror the request that is on the [NASS API documentation](http://quickstats.nass.usda.gov/api), you can use:

```{r eval=FALSE}
    library(rnassqs)
    
    # You must set your api key before requesting data
    nassqs_auth(key = <your api key>)
    
    # Parameters to query on and data call
    params <- list(commodity_desc = "CORN", year__GE = 2012, state_alpha = "VA")
    d <- nassqs(params)
```

Parameters __do not__ need to be capitalized, and also do not need to be in a list format. The following works just as well:

```{r eval=FALSE}
    d <- nassqs(commodity_desc = "corn", year__GE = 2012, state_alpha = "va")
```

You can request data for multiple values of the same parameter by using a simple list as follows:
    
```{r eval=FALSE}
    params <- list(commodity_desc = "CORN", year__GE = 2012, state_alpha = c("VA", "WA"))
    d <- nassqs(params)
```

NASS does not allow GET requests that pull more than 50,000 records in one request. The function will inform you if you try to do that. It will also inform you if you've requested a set of parameters for which there are no records.

Other useful functions include:

```{r eval=FALSE}
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


### Handling inequalities and operators other than "="
The NASS API handles other operators by modifying the variable name. The API can accept the following modifications:

* __LE: <= 
* __LT: < 
* __GT: > 
* __GE: >= 
* __LIKE: like 
* __NOT_LIKE: not like 
* __NE: not equal 

For example, to request corn yields in Virginia and Pennsylvania for all years since 2000, you would use something like:

```{r eval=FALSE}
    params <- list(commodity_desc = "CORN", 
                  year__GE = 2000, 
                  state_alpha = c("VA", "PA"), 
                  statisticcat_desc = "YIELD")
    df <- nassqs(params) #returns data as a data frame.
```

See the [vignette](https://docs.ropensci.org/rnassqs/articles/rnassqs.html) for more examples and details on usage.

## Contributing

Contributions are more than welcome, and there are several ways to contribute:

- Examples: More examples are always helpful. If you use `rnassqs` to query data from 'Quick Stats' and would like to contribute your query, consider submitting a pull request adding your query as a file in [inst/examples/](https://github.com/ropensci/rnassqs/tree/master/inst/examples).
- File an issue: If there is functionality you'd like to see added or something that is confusing, consider [creating an issue](https://github.com/ropensci/rnassqs/issues/new). The best issue contains an example of the problem or feature. Consider the excellent package [reprex](https://github.com/tidyverse/reprex) in creating a reproducible example.
- Contributing documentation: Clarifying and expanding the documentation is always appreciated, especially if you find an area that is lacking and would like to improve it. `rnassqs` uses roxygen2, which means the documentation is at the top of each function definition. Please submit any improvements as a pull request.
- Contributing code: if you see something that needs improving and you'd like to make the changes, contributed code is very welcome. Begin by filing a new issue to discuss the proposed change, and then submit a pull request to address the issue. `rnassqs` follows the style outlined in Hadley Wickham's [R Packages](http://r-pkgs.had.co.nz/style.html). Following this style makes the pull request and review go more smoothly.


## Alternatives

In June 2019 the `usdarnass` package was released on [CRAN](https://cran.r-project.org/package=usdarnass) and is also available to install via [github](https://github.com/rdinter/usdarnass). `usdarnass` has similar functionality to this package.

NASS also provides a daily tarred and gzipped file of their entire dataset. At the time of writing it is approaching 1 GB. You can download that file via their [FTP site](ftp://ftp.nass.usda.gov/quickstats).

The FTP link also contains builds for: NASS census (every 5 years ending with 2 and 7), or data for one of their specific sectors (CROPS, ECONOMICS, ANIMALS & PRODUCTS). At the time of this writing, specific files for the ENVIRONMENTAL and DEMOGRAPHICS sectors are not available.




### Acknowledgements

Thank you to rOpensci reviewers Adam Sparks and Neal Richardson and editor Lincoln Mullen, for their fantastic feedback and assistance. User feedback and use case contributions have been a huge help to make `rnassqs` more accessible and user-friendly. More use cases or feature requests are always welcome!

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```


---
title: "Using rnassqs"
author: "Nicholas A Potter"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using rnassqs}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

`rnassqs` is a package to access the QuickStats API from national agricultural statistics service (NASS) at the USDA. There are at least two good reasons to do this:

1. **Reproducibility**. downloading the data via an R script creates a trail that you can revisit later to see exactly what you downloaded. It also makes it much easier for people seeking to replicate your results to ensure they have the same data that you do.

2. **DRY**. Don't repeat yourself. Downloading data via API makes it easier to download new data as it is released, and to fetch multiple variables, geographies, or time frames without having to manually click through the QuickStats tool for each data request.

In the beginning it can be more confusing, and potentially take more time, but as you become familiar with the variables and calls of the `rnassqs` package and the QuickStats database, you'll be able to quickly and easily download new data.

## API Information

The USDA-NASS Quick Stats API has a graphic interface here: [https://quickstats.nass.usda.gov](https://quickstats.nass.usda.gov). 
Information on the query parameters is found at [https://quickstats.nass.usda.gov/api#param_define](https://quickstats.nass.usda.gov/api#param_define).

## A quick example

First, obtain an API key from the 'Quick Stats' service: [https://quickstats.nass.usda.gov/api](https://quickstats.nass.usda.gov/api). Then we can make a query. Here we request the number of farm operators by operation acreage in Oregon in 2012.

    library(rnassqs)
    
    # Specify the query parameters
    params <- list(
      commodity_desc = "OPERATORS",
      domain_desc = "AREA OPERATED"
      agg_level_desc = "STATE",
      state_alpha = "OR",
      year = 2012
    )
    
    # Check that our record request is under the 50,000 limit
    nassqs_record_count(params)
    
    # Get the data
    d <- nassqs(params)

Parameters need not be specified in a list and need not be capitalized. The following is equivalent

    # Get the data specifying each parameter as a separate argument to the 
    # function `rnassqs`
    d <- nassqs(commodity_desc = "operators", 
                domain_desc = "area operated",
                agg_level_desc = "state",
                state_alpha = "or",
                year = 2012)

## Convenience functions

A growing list of convenience functions makes querying simpler. For example, you can retrieve yields and acres with

    # Set parameters
    params <- list(
      commodity_desc = "APPLES",
      domaincat_desc = "NOT SPECIFIED"
      agg_level_desc = "STATE",
      state_alpha = "OR",
      year = 2012
    )
    
    # Yields and Acres
    yields <- nassqs_yields(params)
    acres <- nassqs_acres(params)



## Detailed usage


### Step 1: Authentication

the QuickStats API requires authentication. You can get an API Key [here](https://quickstats.nass.usda.gov/api). Once you have a key, you can use it in any of the following ways:


#### Add it to your .Renviron file

In your home directory create or edit the `.Renviron` file, and add `NASSQS_TOKEN = <your api key>` to the file. `R` sessions will have the variable set automatically, and `rnassqs` will detect this when querying data. If you use Rstudio, you can also use `usethis::edit_r_environ` to open your `.Renviron` file and add the key. This will create a new system environmental variable when you start a new `R` session. You can also set the environmental variable directly with `Sys.setenv(NASSQS_TOKEN = <your api key>`.


#### Put it in a file

You can add a file to your project directory and ignore it via `.gitignore` if you're using github. The advantage of this method is that you don't have to think about the API key for the rest of the project, but you have to repeat this process for every new project, and you risk forgetting to add it to `.gitignore`. Once the api key is in a file, you can use it like this:

    # Load the api key
    api_key <- readLines("<file name with api key>")
    nassqs_auth(key = api_key)


#### Add it interactively

If you don't want to add the API key to a file or store it in your `.Renviron`, you can enter it in the console in a session. This is less easy because you have to enter (or copy-paste) the key each time you begin an `R` session. In addition, you won't be able to automate running your script, since it will stop and ask you to provide an api key.

    # Checks if the api key is set and prints it. 
    # If it is not set, asks the user to set the value in the console.
    nassqs_auth()


### Step 2: Building Queries

The QuickStats API offers a bewildering array of fields on which to query. `rnassqs` tries to help navigate query building with some functions that return parameter names and valid values for those parameters. `nassqs_params()` provides the parameter names, which at the time of this writing are

```{r}
library(rnassqs)

# returns a list of fields that you can query
nassqs_params()
```

Including parameter names in `nassqs_params` will return a description of the parameter(s) in question:

```{r}
nassqs_params("agg_level_desc", "source_desc")
```

Documentation on all of the parameters is available at [https://quickstats.nass.usda.gov/api#param_define](https://quickstats.nass.usda.gov/api#param_define).

A list of the valid values for a given field is available via `nassqs_param_values(param = <parameter name>)`. For example, 

    nassqs_param_values(param = 'source_desc')

returns a list of valid values for the `source_desc` parameter. 

Building a query often involves some trial and error. One way of developing the query is to use the [QuickStats web interface](https://quickstats.nass.usda.gov/). This is often the fastest method and provides quick feedback on the subset of values for a given query. Alternatively, you can query values for each field as above and iteratively build your query. The query in the end takes the form of a list of parameters that looks like

    params <- list(commodity_desc = "CORN", year__GE = 2012, state_alpha = "VA")

#### Querying a range of values

Most queries will probably be for specific values such as `year = 2012`, but you may also want to query ranges of values. For those queries, append one of the following to the field you'd like to modify:

* __LE: less than or equal
* __LT: less than
* __GE: greater than or equal
* __GT: greater than
* __LIKE: like 
* __NOT_LIKE: not like 
* __NE: not equal

In the above parameter list, `year__GE` is the `year` field with the `__GE` modifier attached to it. The returned data includes all records with year greater than or equal to 2012.

#### Querying multiple values

Multiple values can be queried at once by including them in a simple list with `c()`. For example, if you'd like data from both Washington and Oregon, you can write `state_alpha = c('WA', 'OR')`.

#### Query limits

The API only returns queries that return 50,000 or less records, so it's a good idea to check that before running a query. Do do so, you can use `nassqs_record_count()`. Combined with an assert from the `assertthat` package, you can ensure that your queries are valid before attempting to access the data:

    # Check that the number of returned records will be less than 50000
    params <- list(commodity_desc = "CORN", year__GE = 2012, state_alpha = "VA")
    records <- nassqs_record_count(params)
    assertthat::assert_that(as.integer(records$count) <= 50000)


### Step 3: Running Queries

Once you've built a query, running it is easy:

    # Run a query given a set of parameters and an API key
    nassqs(params = params, key = api_key)


### Step 4. Putting it all together

Putting all of the above together, we have a script that looks like:

    library(rnassqs)
    library(assertthat) #for checking the size of the query

    # Check for the API key. This prints the key if it is set, or asks for it
    # if the session is interactive
    nassqs_auth()
    
    # Get a list of available fields
    parameters <- nassq_params()
    
    # Get valid values for 'commodity_desc'
    nassqs_param_values(param = 'source_desc')
    
    # Set a list of parameters to query on
    params <- list(commodity_desc = "CORN", year__GE = 2012, state_alpha = "VA")
    
    # Check that the number of returned records will be less than 50000
    records <- nassqs_record_count(params)
    assert_that(as.integer(records$count) <= 50000)
    
    # Run a query given a set of parameters and an API key
    d <- nassqs(params = params, key = api_key)
    
    # Run the same query but parse into a data.frame separately
    raw <- nassqs_GET(params = params, key = api_key)
    parsed <- nassqs_parse(raw, as = 'data.frame')


## Lists of parameters and dealing with large queries

The ability of `rnassqs` to iterate over lists of parameters is especially helpful. In some cases you may wish to collect many different sets of data, and in others your queries may be larger than the API restriction of 50,000 records. In both cases iterating over a list of parameters is helpful.

### Iterating to reduce individual query size

Generally the best way to deal with large queries is to make multiple queries subset by year if possible, and by geography if not. Some care is needed if subsetting by geography. Due to suppression of data, the _sum of all counties in a state will not necessarily equal the state value_. Moreover, some data is collected only at specific geographies. It is best to start by iterating over years, so that if you want say all county cash rents on irrigated land for every year since they became available in 2008, you can iterate by doing the following:

    # Define the list of parameters to use repeatedly
    param_list <- list(
      sector_desc = "ECONOMICS",
      commodity_desc = "RENT",
      prodn_practice_desc = "IRRIGATED",
      class_desc = "CASH, CROPLAND",
      agg_level_desc = "COUNTY",
      domaincat_desc = "NOT SPECIFIED")
    
    # Iterate through each year to get data  
    data_list <- lapply(2008:2017, function(yr) {
      params <- param_list
      params[['year']] <- yr
      nassqs(params)
    })
    
    # Using dplyr to bind the data list
    library(dplyr)
    df <- rbind_list(data_list)
    
    # Using data.table to bind the data list
    library(data.table)
    dt <- rbindlist(data_list)

Subsetting by geography works similarly, looping over the geography variable (usually `state_alpha` or `county_code` or the like) in lapply.

### Iterating over lists of parameters

Similar to above, at times it is helpful to make multiple queries and bind the data into a single `data.frame`. For example, you may want to collect the many different categories of acres for every Agricultural Census since 1997, which you can do with something like

    # First define a base parameter list to modify for each new query
    base_params <- list(
      source_desc = "CENSUS",
      sector_desc = "ECONOMICS",
      commodity_desc = "AG LAND",
      agg_level_desc = "COUNTY",
      unit_desc = "ACRES",
      statisticcat_desc = "AREA",
      domain_desc = "TOTAL",
      domaincat_desc = "NOT SPECIFIED",
      year_GE = 1997
    )
    
    # List of parameters that vary for each query
    param_list <- list(
      ag_land_other = list(
        class_desc = "(EXCL CROPLAND & PASTURELAND & WOODLAND)"), 
      ag_land_irr = list(
        prodn_practice_desc = "IRRIGATED",
        class_desc = "ALL CLASSES"),
      ag_woodland = list(
        class_desc = "WOODLAND"),
      ag_pastureland = list(
        class_desc = "PASTURELAND, (EXCL CROPLAND & WOODLAND)"),
      ag_cropland = list(
        class_desc = "CROPLAND"),
      ag_cropland_excl_harvested = list(
        class_desc = "CROPLAND, (EXCL HARVESTED & PASTURED)"),
      ag_cropland_harvested = list(
        class_desc = "CROPLAND, HARVESTED",
        prodn_practice_desc = "ALL PRODUCTION PRACTICES"),
      ag_cropland_harvested_irr = list(
        class_desc = "CROPLAND, HARVESTED",
        prodn_practice_desc = "IRRIGATED")
      )
      
      # Iterate through different variable queries
      data_list <- lapply(param_list, function(var_params) {
        # Create the new parameter list and append the query items that vary
        # by query
        params <- base_params
        for(n in names(var_params)) { 
          params[[n]] <- var_params[[n]]
        }
        nassqs(params)
      })
      
      # Then rbind_list() or rbindlist() as above


## Under the hood

`nassqs` is a wrapper around the `nassqs_GET` function, which uses `httr::GET` to make an HTTP GET request to the Quick Stats API. If you need to access the underlying request object generated by the GET call, you can use `nassqs_GET` to return the request object. The `rnassqs` package also has a `nassqs_parse` function that will process a request object into a data.frame, list, or raw text. `nassqs` does handles both together, but you can replicate that functionality with low-level functions as follows:

    # Make a HTTP GET request and parse into a data.frame with separate
    # function calls. The below is equivalent to 
    # 'nassqs(params, key = api_key)'
    request <- nassqs_GET(params = params, key = api_key)
    parsed <- nassqs_parse(request, as = 'data.frame')
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request.R
\name{nassqs}
\alias{nassqs}
\title{Get data and return a data frame}
\usage{
nassqs(
  ...,
  source_desc = NULL,
  sector_desc = NULL,
  group_desc = NULL,
  commodity_desc = NULL,
  short_desc = NULL,
  domain_desc = NULL,
  domaincat_desc = NULL,
  agg_level_desc = NULL,
  statisticcat_desc = NULL,
  state_name = NULL,
  asd_desc = NULL,
  county_name = NULL,
  region_desc = NULL,
  zip_5 = NULL,
  watershed_desc = NULL,
  year = NULL,
  freq_desc = NULL,
  reference_period_desc = NULL,
  as = c("data.frame", "text", "list")
)
}
\arguments{
\item{...}{either a named list of parameters or a series of parameters to
form the query}

\item{source_desc}{"Program" - Source of data ("CENSUS" or "SURVEY"). Census
program includes the Census of Ag as well as follow up projects. Survey
program includes national, state, and county surveys.}

\item{sector_desc}{"Sector" - Five high level, broad categories useful to
narrow down choices. ("ANIMALS & PRODUCTS", "CROPS", "DEMOGRAPHICS",
"ECONOMICS", or "ENVIRONMENTAL")}

\item{group_desc}{"Group" - Subsets within sector (e.g., under sector_desc =
"CROPS", the groups are "FIELD CROPS", "FRUIT & TREE NUTS", "HORTICULTURE",
and "VEGETABLES").}

\item{commodity_desc}{"Commodity" - The primary subject of interest (e.g.,
"CORN", "CATTLE", "LABOR", "TRACTORS", "OPERATORS").}

\item{short_desc}{"Data Item" - A concatenation of six columns:
commodity_desc, class_desc, prodn_practice_desc, util_practice_desc,
statisticcat_desc, and unit_desc.}

\item{domain_desc}{"Domain" - Generally another characteristic of operations
that produce a particular commodity (e.g., "ECONOMIC CLASS", "AREA
OPERATED", "NAICS CLASSIFICATION", "SALES"). For chemical usage data, the
domain describes the type of chemical applied to the commodity. The
domain_desc = "TOTAL" will have no further breakouts; i.e., the data value
pertains completely to the short_desc.}

\item{domaincat_desc}{"Domain Category" - Categories or partitions within a
domain (e.g., under domain_desc = "SALES", domain categories include $1,000
TO $9,999, $10,000 TO $19,999, etc).}

\item{agg_level_desc}{"Geographic Level" - Aggregation level or geographic
granularity of the data. ("AGRICULTURAL DISTRICT", "COUNTY",
"INTERNATIONAL", "NATIONAL", "REGION : MULTI-STATE", "REGION : SUB-STATE",
"STATE", "WATERSHED", or "ZIP CODE")}

\item{statisticcat_desc}{"Category" - The aspect of a commodity being
measured (e.g., "AREA HARVESTED", "PRICE RECEIVED", "INVENTORY", "SALES").}

\item{state_name}{"State" - State full name.}

\item{asd_desc}{"Ag District" - Ag statistics district name.}

\item{county_name}{"County" - County name.}

\item{region_desc}{"Region" - NASS defined geographic entities not readily
defined by other standard geographic levels. A region can be a less than a
state (SUB-STATE) or a group of states (MULTI-STATE), and may be specific
to a commodity.}

\item{zip_5}{"Zip Code" - US Postal Service 5-digit zip code.}

\item{watershed_desc}{"Watershed" - Name assigned to the HUC.}

\item{year}{"Year" - The numeric year of the data and can be either a
character or numeric vector. Conditional values are also possible, for
example a character vector of ">=1999" of "1999<=" will give years greater
than or equal to 1999. Right now the logical values can either be
greater/less than or equal to with the logical at either the beginning or
end of a string with the year.}

\item{freq_desc}{"Period Type" - Length of time covered ("ANNUAL", "SEASON",
"MONTHLY", "WEEKLY", "POINT IN TIME"). "MONTHLY" often covers more than one
month. "POINT IN TIME" is as of a particular day.}

\item{reference_period_desc}{"Period" - The specific time frame, within a
freq_desc.}

\item{as}{whether to return a data.frame, list, or text string
\code{\link[=nassqs_GET]{nassqs_GET()}}}
}
\value{
a data frame, list, or text string of requested data.
}
\description{
The primary function in the \code{rnassqs} package, \code{nassqs} makes a HTTP GET
request to the USDA-NASS Quick Stats API and returns the data parsed as a
data.frame, plain text, or list. Various other functions make use of \code{nassqs}
to make specific queries. For a data request the Quick Stats API returns
JSON that when parsed to a data.frame contains 39 columns and a varying
number of rows depending on the query. Unfortunately there is not a way to
restrict the number of columns.
}
\examples{
\donttest{
  # Get corn yields in Virginia in 2012
  params <- list(commodity_desc = "CORN",
                 year = 2012,
                 agg_level_desc = "COUNTY",
                 state_alpha = "VA",
                 statisticcat_desc = "YIELD")
  yields <- nassqs(params)
  head(yields)
}
}
\seealso{
\code{\link[=nassqs_GET]{nassqs_GET()}}, \code{\link[=nassqs_yields]{nassqs_yields()}}, \code{\link[=nassqs_acres]{nassqs_acres()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nassqs.R
\docType{package}
\name{rnassqs-package}
\alias{rnassqs-package}
\alias{rnassqs}
\title{rnassqs-package: Access the NASS 'Quick Stats' API}
\description{
rnassqs is a wrapper for the United States Department of Agriculture's
National Agricultural Statistical Service (NASS) 'Quick Stats' API to enable
getting NASS 'Quick Stats' data directly from \R.  Based on the httr API
package guide.
}
\details{
The functions in this package facilitate getting data from NASS 'Quick Stats'.
It handles the API key checking and storage, authorization, and fetching of
data.
}
\references{
\url{http://quickstats.nass.usda.gov}
}
\seealso{
\url{http://quickstats.nass.usda.gov/api}
}
\author{
Nicholas Potter
}
\keyword{package}
\keyword{rnassqs-package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request.R
\name{nassqs_check}
\alias{nassqs_check}
\title{Check the response.}
\usage{
nassqs_check(response)
}
\arguments{
\item{response}{a \code{\link[httr:GET]{httr::GET()}} request result returned from the API.}
}
\value{
nothing if check is passed, or an informative error if not passed.
}
\description{
Check that the response is valid, i.e. that it doesn't exceed 50,000 records
and that all the parameter values are valid. This is used to ensure that
the query is valid before querying to reduce wait times before receiving an
error.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request.R
\name{nassqs_GET}
\alias{nassqs_GET}
\title{Issue a GET request to the NASS 'Quick Stats' API}
\usage{
nassqs_GET(..., api_path = c("api_GET", "get_param_values", "get_counts"))
}
\arguments{
\item{...}{either a named list of parameters or a series of parameters to
use in the query}

\item{api_path}{the API path that determines the type of request being made}
}
\value{
a \code{\link[httr:GET]{httr::GET()}} response object
}
\description{
This is the workhorse of the package that provides the core request
functionality to the NASS 'Quick Stats' API:
\url{https://quickstats.nass.usda.gov/api}.
In most cases \code{\link[=nassqs]{nassqs()}} or other high-level functions should be used.
\code{nassqs_GET()} uses \code{\link[httr:GET]{httr::GET()}} to make a HTTP GET request, which returns a
request object which must then be parsed to a data.frame, list, or other \code{R}
object. Higher-level functions will do that parsing automatically. However,
if you need access to the request object directly, \code{nassqs_GET()} provides
that.
}
\examples{
\donttest{
  # Yields for corn in 2012 in Washington
  params <- list(commodity_desc = "CORN",
                 year = 2012,
                 agg_level_desc = "STATE",
                 state_alpha = "WA",
                 statisticcat_desc = "YIELD")

  # Returns a request object that must be parsed either manually or
  # by using nassqs_parse()
  response <- nassqs_GET(params)
  yields <- nassqs_parse(response)
  head(yields)

  # Get the number of records that would be returned for a given request
  # Equivalent to 'nassqs_record_count(params)'
  response <- nassqs_GET(params, api_path = "get_counts")
  records <- nassqs_parse(response)
  records

  # Get the list of allowable values for the parameters 'statisticcat_desc'
  # Equivalent to 'nassqs_param_values("statisticcat_desc")'
  req <- nassqs_GET(list(param = "statisticcat_desc"),
                    api_path = "get_param_values")
  statisticcat_desc_values <- nassqs_parse(req, as = "list")
  head(statisticcat_desc_values)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{nassqs_record_count}
\alias{nassqs_record_count}
\title{Get a count of number of records for given parameters.}
\usage{
nassqs_record_count(...)
}
\arguments{
\item{...}{either a named list of parameters or a series of parameters to
form the query}
}
\value{
integer that is the number of records that are returned from the
API in response to the query
}
\description{
Returns the number of records that fit a set of parameters. Useful if your
current parameter set returns more than the 50,000 record limit.
}
\examples{
\donttest{
  # Check the number of records returned for corn in 1995, Washington state
  params <- list(
    commodity_desc = "CORN",
    year = "2005",
    agg_level_desc = "STATE",
    state_name = "WASHINGTON"
  )
  
  records <- nassqs_record_count(params) 
  records  # returns 17
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{nassqs_acres}
\alias{nassqs_acres}
\title{Get NASS Area given a set of parameters.}
\usage{
nassqs_acres(
  ...,
  area = c("AREA", "AREA PLANTED", "AREA BEARING", "AREA BEARING & NON-BEARING",
    "AREA GROWN", "AREA HARVESTED", "AREA IRRIGATED", "AREA NON-BEARING", "AREA PLANTED",
    "AREA PLANTED, NET")
)
}
\arguments{
\item{...}{either a named list of parameters or a series of parameters to
form the query}

\item{area}{the type of area to return. Default is all types.}
}
\value{
a data.frame of acres data
}
\description{
Get NASS Area given a set of parameters.
}
\examples{
\donttest{
  # Get Area bearing for Apples in Washington, 2012.
  params <- list(
    commodity_desc = "APPLES",
    year = "2012",
    state_name = "WASHINGTON",
    agg_level_desc = "STATE"
  )
  area <- nassqs_acres(params, area = "AREA BEARING")
  head(area)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/params.R
\name{nassqs_param_values}
\alias{nassqs_param_values}
\title{Get all values for a specific parameter.}
\usage{
nassqs_param_values(param, ...)
}
\arguments{
\item{param}{the name of a NASS quickstats parameter}

\item{...}{additional parameters for which to filter the valid responses.}
}
\value{
a list containing all valid values for that parameter
}
\description{
Returns a list of all possible values for a given parameter. Including
additional parameters will restrict the list of valid values to those for
data meeting the additional parameter restrictions. However, this is only
possible by requesting the entire dataset and then filtering for unique
values. It is recommended to make the query as small as possible if
including additional parameters
}
\examples{
\donttest{
  # See all values available for the statisticcat_desc field. Values may not
  # be available in the context of other parameters you set, for example
  # a given state may not have any 'YIELD' in blueberries if they don't grow
  # blueberries in that state.
  # Requires an API key:
  
  nassqs_param_values("source_desc")

  # Valid values for a parameter given a specific set of additional
  # parameters
  nassqs_param_values("commodity_desc", state_fips_code = "53", county_code = "077", year = 2017, group_desc = "EXPENSES")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request.R
\name{nassqs_parse}
\alias{nassqs_parse}
\title{Parse a response object from \code{nassqs_GET()}.}
\usage{
nassqs_parse(req, as = c("data.frame", "list", "text"), ...)
}
\arguments{
\item{req}{the GET response from \code{\link[=nassqs_GET]{nassqs_GET()}}}

\item{as}{whether to return a data.frame, list, or text string}

\item{...}{additional parameters passed to \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}} or
\code{\link[utils:read.csv]{utils::read.csv()}}}
}
\value{
a data frame, list, or text string of the content from the response.
}
\description{
Returns a data frame, list, or text string. If a data.frame, all columns
except \code{year} strings because the 'Quick Stats' data returns suppressed data
as '(D)', '(Z)', or other character indicators which mean different things.
Converting the value to a numerical results in NA, which loses that
information.
}
\examples{
\donttest{
  # Set parameters and make the request
  params <- list(commodity_desc = "CORN",
                 year = 2012,
                 agg_level_desc = "STATE",
                 state_alpha = "WA",
                 statisticcat_desc = "YIELD")
  response <- nassqs_GET(params)

  # Parse the response to a data frame
  corn <- nassqs_parse(response, as = "data.frame")
  head(corn)

  # Parse the response into a raw character string.
  corn_text<- nassqs_parse(response, as = "text")
  head(corn_text)

  # Get a list of parameter values and parse as a list
  response <- nassqs_GET(list(param = "statisticcat_desc"),
                    api_path = "get_param_values")
  statisticcat_desc_values <- nassqs_parse(response, as = "list")
  head(statisticcat_desc_values)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{nassqs_auth}
\alias{nassqs_auth}
\title{Get/Set the environmental variable NASSQS_TOKEN to the API key}
\usage{
nassqs_auth(key)
}
\arguments{
\item{key}{the API key (obtained from \url{https://quickstats.nass.usda.gov/api})}
}
\description{
If the API key is provided, sets the environmental variable. You can set
your API key in four ways:
}
\details{
\enumerate{
\item directly or as a variable from your \code{R}
program: \verb{nassqs_auth(key = "<your api key>"}
\item by setting \code{NASSQS_TOKEN} in your \code{R} environment file (you'll never have
to enter it again).
\item by entering it into the console when asked (it will be stored for the
rest of the session.)
}
}
\examples{
# Set the API key
nassqs_auth(key = "<your api key>")
Sys.getenv("NASSQS_TOKEN")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/params.R
\name{nassqs_params}
\alias{nassqs_params}
\title{Return list of NASS QS parameters.}
\usage{
nassqs_params(...)
}
\arguments{
\item{...}{a parameter, series of parameters, or a list of parameters that
you would like a description of. If missing, a list of all available
parameters is returned.}
}
\value{
a list of all available parameters or a description of a subset
}
\description{
Contains a simple hard-coded list of all available parameters. If no
parameter name is provided, returns a list of all parameters. More
information can be found in the API documentation on parameters found at
\url{https://quickstats.nass.usda.gov/api#param_define}.
}
\examples{
# Get a list of all available parameters
nassqs_params()

# Get information about specific parameters
nassqs_params("source_desc", "group_desc")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/params.R
\name{nassqs_fields}
\alias{nassqs_fields}
\title{Deprecated: Return list of NASS QS parameters.}
\usage{
nassqs_fields(...)
}
\arguments{
\item{...}{a parameter, series of parameters, or a list of parameters that
you would like a description of. If missing, a list of all available
parameters is returned.}
}
\description{
Deprecated. Use \code{\link[=nassqs_params]{nassqs_params()}} instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{nassqs_yields}
\alias{nassqs_yields}
\title{Get yield records for a specified crop.}
\usage{
nassqs_yields(...)
}
\arguments{
\item{...}{either a named list of parameters or a series of parameters to
form the query}
}
\value{
a data.frame of yields data
}
\description{
Returns yields for other specified parameters. This function is intended to
simplify common requests.
}
\examples{
\donttest{
  # Get yields for wheat in 2012, all geographies
  params <- list(
    commodity_desc = "WHEAT", 
    year = "2012", 
    agg_level_desc = "STATE",
    state_alpha = "WA")
    
  yields <- nassqs_yields(params)
  head(yields)
}
}
