---
title: 'SymbiotaR2: An R Package for Accessing Symbiota2 Data'
authors: 
- affiliation: 1
  name: Austin Koontz
  orcid: 0000-0002-6103-5894
- affiliation: 3
  name: Benjamin Brandt
- affiliation: 4
  name: Curtis Dyreson
- affiliation: "1,2"
  name: William D Pearse
  orcid: 0000-0002-6241-3164
affiliations:
- index: 1
  name: Department of Biology & Ecology Center, Utah State University, Logan, Utah,
    USA
- index: 2
  name: Department of Life Sciences, Imperial College London, Silwood Park Campus, 
    Buckhurst Rd., Ascot, Berkshire, SL5 7PY UK
- index: 3
  name: Northern Arizona University, Arizona, USA
- index: 4
  name: Department of Computer Science, Utah State University, Logan, Utah, USA
date: "10/08/2020"
output: html_document
bibliography: paper.bib
tags:
  - R
  - Symbiota
  - specimen-records
  - biodiversity
---
# Summary

`SymbiotaR2` is an R [@R2018] package for easily accessing and
handling specimen-based Symbiota2 data within an R environment,
allowing anyone to download digitized biological collection data from
any established database. Symbiota2 is the updated version of Symbiota
[@Symbiota2014], a widely used software platform that grants access to
data from >750 museums and herbaria worldwide. Through a complete
refactoring of the Symbiota code, the structure of Symbiota2 places an
emphasis on modularity and accessibility.

# Statement of need

The release of `SymbiotaR2` is motivated by the ongoing development 
of Symbiota2. Several R packages already exist for accessing 
data from standard Symbiota portals; because Symbiota2 
is a complete refactoring of the original framework, 
a new package is required to interface with new Symbiota2 portals.

The goal of this package is to allow users to access digitized 
biological specimen data hosted via Symbiota2 portals quickly and 
efficiently. This data could include geographic, taxonomic, or genetic 
information tied to a recorded specimen, as well as pertinent collection 
information, images, and publications. By making this data accessible in R, 
`SymbiotaR2` is part of the broader Symbiota effort of creating a 
collaborative environment to share biodiversity data more widely. 

# Demonstration

The code below provides an example of a user accessing two different attributes
from the same entry within a Symbiota2  database. First, data from a particular 
entry (here, `id=28`) is pulled into R using the "Occurrences" API endpoint. 
Then, different attributes associated with this entry (`reproductiveCondition`,
`decimalLatitude`, and `decimalLongitude`) are printed. The specification of the 
`url` argument (at the beginning of the example) allows for users to reference 
the web address corresponding to any Symbiota2 portal of interest.
```{R}
# Specify web address from which to access Symbiota2 API
> url <- "http://demo.portal.address/api/"

# Pull occurrence information associated with database entry (id=28)
> test.Occ <- Occurrences(id=28, url=url)

# Find reproductive condition for entry of interest
> print(test.Occ$reproductiveCondition)
[1]   "bud"
# Retrieve and concatenate coordinates for entry of interest
> test.coord <- c(test.Occ$decimalLatitude, test.Occ$decimalLongitude)
> print(test.coord)
[1]   50.70 -103.65
```
From here, these R objects be manipulated and used for any downstream analysis. 
The library provides built in commands for all of the default API endpoints 
included for any given Symbiota2 portal. These endpoints include access to 
several different resources: trait data, taxonomic classification, 
collector information, and much more. 

# References
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/SymbiotaR2)](https://cran.r-project.org/package=SymbiotaR2)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://api.travis-ci.org/ropensci/Symbiota2.svg)](https://travis-ci.org/ropensci/SymbiotaR2)
[![codecov](https://codecov.io/gh/ropensci/SymbiotaR2/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/SymbiotaR2)
[![DOI](https://zenodo.org/badge/190439935.svg)](https://zenodo.org/badge/latestdoi/190439935)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02917/status.svg)](https://doi.org/10.21105/joss.02917)
<!-- badges: end -->

# SymbiotaR2

Austin Koontz, Benjamin Brandt, Curtis Dyreson, and William D. Pearse

## Overview

Designed to assist students, taxonomists, educators, and anyone
working with biological specimens,
[Symbiota](https://symbiota.org/docs/) is an open-source content
management system designed to integrate virtual biodiversity
databases. Over 750 natural history collections use Symbiota, which
has motivated the development of the Symbiota2 update to the Symbiota
platform. This package, _SymbiotaR2_, allows users to access and
download specimen- and observation-based data from a published
Symbiota2 portal. If you're looking for something similar for the
original Symbiota, take a look at
[rSymbiota](https://github.com/FranzKrah/rSymbiota).

More information about Symbiota2 can be found on the
[website](https://symbiota.org/docs/symbiota2-project/). The GitHub
page for Symbiota2 can be found
[here](https://github.com/Symbiota2/Symbiota2); if you want to set up
a new Symbiota2 portal, please follow the instructions on its
[documentation
site](https://symbiota2.github.io/Symbiota2/setup/installation.html).
Finally, a review of the original Symbiota platform is offered in
[Gries et al., 2014](https://bdj.pensoft.net/articles.php?id=1114).

## Example workflow for SymbiotaR2

In general, there are four steps for using SymbiotaR2:

1. Determine the URL of the Symbiota2 portal you wish to access data
   from; its "API endpoint" is probably its web address with "api"
   appended to it (see below). Remember that your particular portal
   may not have enabled data download.
2. Load the SymbiotaR2 package (see [Installation](#inst) below for
   install instructions).
3. Find the function corresponding to the kind of data you wish to
pull from the Symbiota2 portal (e.g., `Coordinates` for co-ordinate
data).
    - Functions are named after the resources they download, and are
  grouped according to the relevant API call.
    - Note that, because each Symbiota2 portal owner can load their
  own plugins into the API, it's possible that not every API endpoint
  will be covered.
    - You can find a full list by typing `library(help=SymbiotaR2)`
4. Call the function, specifying the Symbiota2 portal using the `url`
argument (see [Example](#ex) and [Portal Specification](#portspec)
below).

The R code below provides an example of how to install a `Coordinates` 
resource from an example Symbiota2 portal

```{R}
# Step 1 - Find the URL
myURL <- "http://imaginary-symbiota2-portal.com/api"
# Step 2 - Install the package (only required the first time you're using it)
install.packages("SymbiotaR2") 
# Step 3 - Load the package
library(SymbiotaR2)
# Step 4 - Choose the kind of data you want
library(help=SymbiotaR2)
# Step 5 - Download your data
myCoordinates <- Coordinates(url = myURL)
```

From here, the `myCoordinates` object can be used as desired.

## <a name="portspec"></a>Setting up a default portal for download

Calling `SymbiotaR2_setup` will specify a default portal URL to use for 
all subsequent SymbiotaR2 calls made within your R session.
Specifing a different `url` argument will let you refer to
a portal besides the default. The code below offers an example:

```{R}
SymbiotaR2_setup("http://imaginary-symbiota-portal.com/api")
Coordinates() # Download from http://imaginary-symbiota-portal.com/api
Coordinates("http://another-imaginary-portal.com/api") # Download from a different portal
```

If the `append` argument in your `SymbiotaR2_setup` call is set to 
`TRUE`, then the specified `url` argument will be saved as the default
to your `.Rprofile`, allowing it to be used everytime you start up `R`.

## <a name="inst"></a>Installation

This package is not currently up on CRAN (as it's being developed), 
but it can be downloaded by calling:

```{R}
library(devtools)
install_github("ropensci/SymbiotaR2")
```

Once it has passed peer review, you will be able to install it by
running:

```{R}
install.packages("SymbiotaR2")
```

Load the package using:

```{R}
library(SymbiotaR2)
```

## (For developers) Unit tests

All of the package functions come with tests, for both pulling a
single Symbiota2 resource (using the `id` argument), or a collection
of resources (using `page`). Tests for each function are contained in
the `tests/testthat` directory. Running these tests requires you have
access to a fully configured Symbiota2 test instance, complete with
demo data, which is both time-consuming to setup and then
time/bandwidth-consuming to run the tests. We therefore release cached
data downloads, generated using `vcr`, for use with this package.

While more information about the `vcr` package can be found on [the
`vcr` page on GitHub](https://github.com/ropensci/vcr), you don't need
to understand how `vcr` works to run the tests for yourself. Instead,
do the following:

1. Build the package as you would normally, with something like `R CMD
   build SymbiotaR2` from the command line.
2. Check the package as you would normally, with something like `R CMD
   check SymbbiotaR2_0.0-1` from the command line.

If you want to add new tests, or new functions that address new API
endpoints (perhaps because you have written a Symbiota2 plugin and
want it to work with this package), do the following:

1. Setup a Symbiota2 instance with the canned example data.
2. If you are adding support for a new API endpoint, make a new file
   in `tests/testthat` for your tests. Otherwise, add to one of the
   existing files.
3. Write your test, following the coding style of the other tests,
   particularly with respect to setting up the `vcr` _cassette_. Note
   that the folder `fixtures` contains the cassettes, and that
   `SymbiotaR2` makes use of the file
   `tests/testthat/helper-SymbiotaR2.R` to setup the automatic
   tests. See point 4 below.
4. When writing/checking your test, set the `url` variable at the top
   of the script to be wherever your test instance is. When committing
   your code to submit a pull request (see point 5), change it to the
   address at the top of the other tests (currently
   `http://a02235015-6.bluezone.usu.edu/api/`).
5. When you are finished, submit a pull request to the `master` branch
   of this repository. Please use the pull request template and follow
   the contributor guidelines.

Here is an example of what a piece of testing code may look like:

```{R}
context("AccessStats")
vcr::use_cassette(name = "AccessStats_id", {
  data <- AccessStats(id = 4, url = url)
})
test_that("AccessStats_id", {
  expect_equal(length(data), 12)
  expect_type(data, "list")
})
```

The `data <- AccessStates(url = url, id = 4)` line is the Symbiota2
call, and the `test_that` block below it contains the test
conditions--here, that the `data` object is a `list` of
length 12. 

## Citation

If you use SymbiotaR2 as a component to publication, we ask that you 
properly cite this R package. Use `citation("SymbiotaR2")` to see how to 
cite SymbiotaR2.

## Contributions

Please check out our [contribution guidelines](https://github.com/ropensci/SymbiotaR2/blob/master/.github/CONTRIBUTING.md).

If you'd like to contribute to the tests in this package (e.g.
to include a test for a new plugin), remember that you'll need an
accessible SymbiotaR2 instance to determine test criteria.

As the underlying Symbiota2 API and the R package are still 
being worked on, we generally recommend holding off
on any package development until these are finalized.
If you're interested in contributing to the
package in the future, though, please [drop one of us an email and
we'll let you know when we're ready](http://pearselab.com/team.html)!

-----

Please note that this package is released with a [Contributor Code
of
Conduct](https://ropensci.org/code-of-conduct/).
By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# Version 0.0-1 (2020-10-01	)

* Package release 

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

Please note that the SymbiotaR2 project is released with a
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
title: "Downloading data from Symbiota2 portals in R"
author: "Austin Koontz & William D. Pearse"
date: "`r Sys.Date()`"
output:
  pdf_document: default
  html_document: default
vignette: >
  %\VignetteEngine{knitr::rmarkdown} 
  %\VignetteIndexEntry{Vignette Title} 
  %\usepackage[utf8]{inputenc}
---

## Overview
[Symbiota](https://symbiota.org/docs/)
is an open-source content management system built for the purpose of 
integrating virtual biodiversity databases. Currently used by over 
700 natural history collections, containing more than 30 million specimens,
Symbiota is an essential tool for digitizing biological specimen data. 
In an effort to expand modularity and accessibility, [Symbiota2](https://symbiota.org/docs/symbiota2-project/)
is an improved, refactored version of the original Symbiota core code 
structure, designed based on user feedback. While packages do exist for 
accessing Symbiota portals (for instance, see the [rSymbiota](https://github.com/FranzKrah/rSymbiota)
package), R users currently cannot access the data offered by Symbiota2. 
Here, we describe **SymbiotaR2**, a package built to address this need 
by allowing users to access Symbiota2 portals in an R environment.

Below, we provide a general workflow for using SymbiotaR2, a description of
the command structure, code for installing the package, and examples of using 
SymbiotaR2 functions. The GitHub page for the Symbiota2 software can be found [here](https://github.com/Symbiota2/Symbiota2),
and instructions for setting up a new Symbiota2 portal can be found on the
[documentation site](https://symbiota2.github.io/Symbiota2/setup/installation.html).
Finally, a review of the original Symbiota platform is offered in 
[Gries et al., 2014](https://bdj.pensoft.net/articles.php?id=1114).

## Workflow
SymbiotaR2 allows R users to download data from specified Symbiota2 portals, 
granting access to thousands of digitized flora and fauna specimen records 
across the United States. It does this by querying endpoints in the 
Symbiota2 API, then downloading a JSON object containing the requested data 
to a temporary directory on the local computer. The JSON object is then 
converted into an R format that is straightforward and easy to use. 
The general argument structure of all SymbiotaR2 functions allows users 
to specify whether they want to pull a single resource or a collection 
of resources from the API.

In general, there are four steps for using SymbiotaR2:

1. Determine the URL of the Symbiota2 portal you wish to access data
   from; its "API endpoint" is probably its web address with "api"
   appended to it (see below). Remember that your particular portal
   may not have enabled data download.
2. Load the SymbiotaR2 package (see [Installation](#inst) below for
   install instructions).
3. Find the function corresponding to the kind of data you wish to
pull from the Symbiota2 portal (e.g., `Coordinates` for co-ordinate
data).
    - Functions are named after the resources they download, and are
  grouped according to the relevant API call.
    - Note that, because each Symbiota2 portal owner can load their
  own plugins into the API, it's possible that not every API endpoint
  will be covered.
    - You can find a full list by typing `library(help=SymbiotaR2)`
4. Call the function, specifying the Symbiota2 portal using the `url`
argument (see [Example](#ex) and [Portal Specification](#portspec)
below).

## <a name="inst"></a>Installation
SymbiotaR2 can be downloaded by calling:

```
library(devtools)
install_github("pearselab/SymbiotaR2")
```

Once it has passed peer review, you will be able to install it by
running:

```
install.packages("SymbiotaR2")
```

Load the package using:

```
library(SymbiotaR2)
```

#### Portal Specification
`SymbiotaR2_setup` will save to your `.Rprofile` a default URL, 
for automatic reference. Specifying a different `url` argument will let you
refer to a portal besides the default. The code below demonstrates this:

```
SymbiotaR2_setup("http://imaginary-symbiota-portal.com/api", append=TRUE)

Coordinates() # Download from http://imaginary-symbiota-portal.com/api
Coordinates("http://another-imaginary-portal.com/api") # Download from a different portal
```

## <a name="ex"></a>Example
SymbiotaR2 consists of commands pulling from the Checklists, 
Collections, Crowdsource, Exsiccati, Glossary, ImageProcessor, Key, Media, 
Occurrence, Reference, Taxa, Traits, and UserRoles API families of the 
specified Symbiota2 portal. Note that because each Symbiota2 portal 
owner can load their own plugins into the API, it's possible that 
not every possible API endpoint from the specified Symbiota2 
instance will be covered. 

Below, we provide an example of pulling a single `Taxa` resource into the R
environment, by specifying an `id` argument in the command call
(using a random, nonexistent URL). Please note that this example won't
work for users (as they need to specify a working Symbiota2 portal
they can access), but is included to demonstrate typical usage:

```
myURL <- "http://imaginary-symbiota2-portal.com/api"
myTaxa <- Taxa(id = 12, url = myURL)
str(myTaxa)
```
```
List of 23
 $ @context             : chr "/api/contexts/Taxa"
 $ @id                  : chr "/api/taxa/12"
 $ @type                : chr "Taxa"
 $ id                   : num 12
 $ rankId               : chr "/api/taxa/ranks/31"
 $ scientificName       : chr "Polygonum bistortoides"
 $ unitIndicator1       : logi NA
 $ unitName1            : chr "Polygonum"
 $ unitIndicator2       : logi NA
 $ unitName2            : chr "bistortoides"
 $ unitIndicator3       : logi NA
 $ unitName3            : logi NA
 $ author               : chr "Pursh"
 $ phylogenySortSequence: logi NA
 $ status               : chr "AZTT-USDA Plants consistant"
 $ source               : logi NA
 $ notes                : logi NA
 $ hybrid               : logi NA
 $ securityStatus       : num 0
 $ modifiedTimestamp    : logi NA
 $ initialTimestamp     : chr "2019-01-11T21:44:39+00:00"
 $ modifiedUserId       : logi NA
 $ taxaAuthorityId      : list()
```

If a collection of resources from the Symbiota2 API needs to come into 
the R environment, then the `page` argument can be specified in place of
`id` to retrieve a list of resources (here, as a `data.frame`):

```
myURL <- "http://imaginary-symbiota2-portal.com/api"
myCoordinates <- Coordinates(page = 1, url = myURL)
str(my.Coordinates)
```
```
'data.frame':	5 obs. of  2 variables:
 $ latitude : num  32.2 32.2 32.2 32.2 32.2
 $ longitude: num  -111 -111 -111 -111 -111
```

If neither an `id` or a `page` argument is provided, the functions are written
to return the list of resources at `page = 1`. Once downloaded, 
these R objects can be taken and manipulated as needed for any
downstream processes.

## Troubleshooting
The code for SymbiotaR2 is structured hierarchically, and includes parameter 
type checking to ensure arguments are provided in the proper format. 
Additionally, all commands include a URL check (`.check.url`), 
which confirms the following:

1. URL provided refers to an accessible website, and
2. that website is a working Symbiota2 portal.

The second step consists of an API call made at the end of the URL check.
If either step fails, the error below will be triggered:

```
Error in .check.url(badURL) : 
  URL http://incorrect-portal-address.com/api cannot be reached; is it a valid Symbiota2 portal API?
```

If this error is received, make sure your portal address is spelled correctly. 
Note that functions are designed such that a forward slash (`/`) at the end 
of the URL is optional. If your URL is correctly spelled, make sure that 
the Symbiota2 portal manager has allowed you access to the portal. 

#### (For developers) Unit tests

All SymbiotaR2 functions come with tests, for both pulling a
single SymbiotaR2 resource (using the `id` argument), or a collection
of resources (using `page`). Tests for each function are contained in
the `tests/testthat` directory. Running these tests requires you have
access to a fully configured SymbiotaR2 test instance, complete with
demo data, which is both time-consuming to setup and then
time/bandwidth-consuming to run the tests. We therefore release cached
data downloads, generated using `vcr`, for use with this package.
Information about the `vcr` package can be found on [the
`vcr` page on GitHub](https://github.com/ropensci/vcr). 

To run the package tests, do the following:

1. Build the package as you would normally, with something like `R CMD
   build SymbiotaR2` from the command line.
2. Check the package as you would normally, with something like `R CMD
   check SymbbiotaR2_0.0-1` from the command line.

If you want to add new tests, or new functions that address new API
endpoints, do the following:

1. Setup a Symbiota2 instance with the canned example data.
2. If you are adding support for a new API endpoint, make a new file
   in `tests/testthat` for your tests. Otherwise, add to one of the
   existing files.
3. Write your test, following the coding style of the other tests,
   particularly with respect to setting up the `vcr` _cassette_. Note
   that the folder `fixtures` contains the cassettes, and that
   `SymbiotaR2` makes use of the file
   `tests/testthat/helper-SymbiotaR2.R` to setup the automatic
   tests. See point 4 below.
4. When writing/checking your test, set the `url` variable at the top
   of the script to be wherever your test instance is. When committing
   your code to submit a pull request (see point 5), change it to the
   address at the top of the other tests (currently
   `http://a02235015-6.bluezone.usu.edu/api/`).
5. When you are finished, submit a pull request to the `master` branch
   of this repository. Please use the pull request template and follow
   the contributor guidelines.

Here is an example of what a piece of testing code may look like:

```
context("AccessStats")
vcr::use_cassette(name = "AccessStats_id", {
  data <- AccessStats(id = 4, url = url)
})
test_that("AccessStats_id", {
  expect_equal(length(data), 12)
  expect_type(data, "list")
})
```

The `data <- AccessStates(url = url, id = 4)` line is the Symbiota2
call, and the `test_that` block below it contains the test
conditions--here, that the `data` object is a `list` of
length 12. 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Key.R
\name{Key}
\alias{Key}
\alias{CharacterHeading}
\alias{Characters}
\alias{CharacterStateImages}
\alias{CharacterStates}
\alias{DescriptionDeletions}
\title{Retrieves Key resources from the Symbiota2 server}
\usage{
CharacterHeading(id, page, url = NULL)

Characters(id, page, url = NULL)

CharacterStateImages(id, page, url = NULL)

CharacterStates(id, page, url = NULL)

DescriptionDeletions(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Key resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a Characters resource (id = 3), from a (nonexistent) dummy portal
object <- Characters(id = 3, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Crowdsource.R
\name{Crowdsource}
\alias{Crowdsource}
\alias{Central}
\alias{Queue}
\title{Retrieves Crowdsource resources from the Symbiota2 server}
\usage{
Central(id, page, url = NULL)

Queue(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Crowdsource resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Acquiring a Queue resource (id = 2), from a (nonexistent) dummy portal
object <- Queue(id = 2, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
\name{pez-internal}
\alias{.Random.seed}
\alias{.eco.evo.clade.regression}
\alias{.eco.evo.regression}
\alias{.eco.null}
\alias{.eco.phy.regression}
\alias{.plot.regression}
\alias{.prepare.regression.output}
\alias{.summary.regression}
\alias{summary.phy.structure}
\title{Internal pez Functions}
\alias{.removeErrors}
\description{
 Internal pez functions
}
\details{


  These are not to be called by the user.}
  \keyword{internal}% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Reference.R
\name{Reference}
\alias{Reference}
\alias{LookupReferenceTypes}
\title{Retrieves Reference resources from the Symbiota2 server}
\usage{
LookupReferenceTypes(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Reference resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a LookupReferenceType resource (id = 1), 
from a (nonexistent) dummy portal
object <- LookupReferenceType(id = 1, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Occurrence.R
\name{Occurrences}
\alias{Occurrences}
\alias{AccessStats}
\alias{Determinations}
\alias{Duplicates}
\alias{EditLocks}
\alias{Edits}
\alias{FullText}
\alias{GuidDeterminations}
\alias{GuidOccurrences}
\alias{LookupChronostratigraphy}
\alias{LookupCounties}
\alias{LookupCountries}
\alias{LookupStateProvinces}
\alias{UploadMappings}
\alias{UploadParameters}
\alias{Verification}
\alias{Associations}
\alias{Comments}
\alias{DatasetLink}
\alias{Datasets}
\alias{Exchange}
\alias{Loans}
\title{Retrieves Occurrence resources from the Symbiota2 server}
\usage{
AccessStats(id, page, url = NULL)

Determinations(id, page, url = NULL)

Duplicates(id, page, url = NULL)

EditLocks(id, page, url = NULL)

Edits(id, page, url = NULL)

FullText(id, page, url = NULL)

GuidDeterminations(id, page, url = NULL)

GuidOccurrences(id, page, url = NULL)

LookupChronostratigraphy(id, page, url = NULL)

LookupCounties(id, page, url = NULL)

LookupCountries(id, page, url = NULL)

LookupStateProvinces(id, page, url = NULL)

UploadMappings(id, page, url = NULL)

UploadParameters(id, page, url = NULL)

Verification(id, page, url = NULL)

Associations(id, page, url = NULL)

Comments(id, page, url = NULL)

DatasetLink(id, page, url = NULL)

Datasets(id, page, url = NULL)

Exchange(id, page, url = NULL)

Loans(id, page, url = NULL)

Occurrences(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Occurrence resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a page of Occurrences, from a (nonexistent) dummy portal
entries <- Occurrence(page = 6, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SymbiotaR2-package.R
\docType{package}
\name{SymbiotaR2}
\alias{SymbiotaR2}
\alias{package-SymbiotaR2}
\alias{SymbiotaR2-package}
\title{Downloading data from Symbiota2 portals into R}
\description{
This package allows users to access and 
download from Symbiota2, a content management system 
for biodiveristy data.
}
\section{About}{

Symbiota2 is the improved, refactored version of Symbiota, 
an open source content management system for biological
specimen data. SymbiotaR2 allows users to access the data
available at Symbiota2 portals. By specifying the URL 
of the relevant portal, and the resource to be downloaded,
users can use SymbiotaR2 to deliver biological specimen-data 
in an R format.
}

\section{Code Structure}{

Package functions are organized by API family, which generally
group the functions by the type of resource they pull from the portal.
Each function can either return an individual resources (through
specifying the `id` argument) or a collection of resources (through
specifying the `page` argument). After providing either the `id` 
or the `page` of resources, and the URL of the relevant portal, 
SymbiotaR2 will return an R object (for `id`, usually a list; for
`page`, usually a data.frame).
}

\section{Portal Specification}{

All SymbiotaR2 commands require a URL that directs to the Symiobta2
portal to download data from. Users need to make sure they are granted
access to a Symbiota2 portal before trying to download data from it. 

The address of a Symbiota2 portal is provided as the `url` string 
argument to each function. To specify a default URL, use the 
`SymbiotaR2_setup` function, which will the default url to your
.Rprofile. 

This package only allows users to access data from existing Symbiota2
portals; to create a new Symbiota2 portal, see the documentation at
https://symbiota2.github.io/Symbiota2/setup/installation.html
}

\examples{
\dontrun{
myURL <- "http://ImaginarySymbiota2Portal.com/api"
myTaxa <- Taxa(id = 12, url = myURL)
str(myTaxa)

myOccurrences <- Occurrence(page = 2, url = myURL)
length(myOccurrences)
}
}
\references{
https://symbiota.org/docs/

Gries, C., Gilbert, E. E., & Franz, N. M. (2014). Symbiota - A virtual platform for creating voucher-based biodiversity information communities. Biodiversity Data Journal, 2, e1114.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImageProcessor.R
\name{ImagePRocessor}
\alias{ImagePRocessor}
\alias{Projects}
\alias{RawLabels}
\title{Retrieves ImageProcessor resources from the Symbiota2 server}
\usage{
Projects(id, page, url = NULL)

RawLabels(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve ImageProcessor resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Acquiring a RawLabels resource (id = 1), from a (nonexistent) dummy portal
object <- RawLabels(id = 1, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Exsiccati.R
\name{Exsiccati}
\alias{Exsiccati}
\alias{Numbers}
\alias{Titles}
\title{Retrieves Exsiccati resources from the Symbiota2 server}
\usage{
Numbers(id, page, url = NULL)

Titles(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Exsiccati resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Acquiring a Titles resource (id = 3), from a (nonexistent) dummy portal
object <- Titles(id = 3, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Media.R
\name{Media}
\alias{Media}
\alias{TagKey}
\title{Retrieves Media resources from the Symbiota2 server}
\usage{
TagKey(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Media resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a page of TagKey resources, from a (nonexistent) dummy portal
object <- TagKey(page = 1, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Checklists.R
\name{Checklists}
\alias{Checklists}
\alias{ChecklistProjects}
\alias{Coordinates}
\alias{TaxaLink}
\title{Retrieves Checklists resources from the Symbiota2 server}
\usage{
ChecklistProjects(id, page, url = NULL)

Coordinates(id, page, url = NULL)

TaxaLink(id, page, url = NULL)

Checklists(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Checklist resources
from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a Coordinates resource (id = 1), from a (nonexistent) dummy portal
object <- Coordinates(id = 1, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Glossary.R
\name{Glossary}
\alias{Glossary}
\alias{TermLink}
\title{Retrieves Glossary resources from the Symbiota2 server}
\usage{
Glossary(id, page, url = NULL)

TermLink(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Glossary resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Acquiring a page of Glossary resources, from a (nonexistent) dummy portal
glossPage <- Glossary(page = 1, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Taxa.R
\name{Taxa}
\alias{Taxa}
\alias{Authorities}
\alias{DescriptionBlock}
\alias{Synonymy}
\title{Retrieves Taxa resources from the Symbiota2 server}
\usage{
Taxa(id, page, url = NULL)

Authorities(id, page, url = NULL)

DescriptionBlock(id, page, url = NULL)

Synonymy(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Taxa resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a page of Taxa resources, from a (nonexistent) dummy portal
object <- Taxa(page = 2, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Miscellaneous.R
\name{Miscellaneous}
\alias{Miscellaneous}
\alias{Configurations}
\alias{LookupLanguages}
\alias{SchemaVersion}
\title{Retrieves Miscellaneous resources from the Symbiota2 server}
\usage{
Configurations(id, page, url = NULL)

LookupLanguages(id, page, url = NULL)

SchemaVersion(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Miscellaneous resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a Configurations resource (id = 4), from a (nonexistent) dummy portal
object <- Configurations(id = 4, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Collection.R
\name{Collection}
\alias{Collection}
\alias{Categories}
\alias{Institutions}
\alias{Stats}
\alias{Collections}
\title{Retrieves Collection resources from the Symbiota2 server}
\usage{
Categories(id, page, url = NULL)

Institutions(id, page, url = NULL)

Stats(id, page, url = NULL)

Collections(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Collection resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling in a page of Institutions, from a (nonexistent) dummy portal
ints <- Institutions(page = 3, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Traits.R
\name{Traits}
\alias{Traits}
\title{Retrieves Traits resources from the Symbiota2 server}
\usage{
Traits(id, page, url = NULL)
}
\arguments{
\item{id}{id value (usually \code{numeric}, but not always) 
used to refer to the specific resource to pull from the database}

\item{page}{\code{numeric} value referring to the page of resources to pull.
If neither an id or a page parameter is provided, 
function will pull the first page of resources (i.e. \code{page=1})}

\item{url}{URL string of the Symbiota2 portal to be connected to. 
A trailing \code{/} will be appended, if it is not given.}
}
\value{
If using \code{id}, the specific resource specified; 
if using page, the \code{page} specified of resources
}
\description{
Functions that retrieve Traits resources from the server previously connected to.
Each function either retrieves an individual resource or a page of resources,
depending on the arguments provided.
}
\note{
To specify a default URL to refer to, see \code{\link[SymbiotaR2:SymbiotaR2_setup]{SymbiotaR2_setup()}}
}
\examples{
\dontrun{
# Pulling a Traits resource (id = 4), from a (nonexistent) dummy portal
object <- Traits(id = 4, url = "http://dummy-portal.com/api/")
}
}
\author{
Austin Koontz
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Defaults.R
\name{SymbiotaR2_setup}
\alias{SymbiotaR2_setup}
\title{Set default URL for Symbiota2 portal download}
\usage{
SymbiotaR2_setup(url, append = FALSE, verbose = TRUE)
}
\arguments{
\item{url}{URL of Symbiota2 portal (a trailing \code{/} will be
appended, if it is not given)}

\item{append}{\code{logical} of whether to attempt to append this
to your \code{.Rprofile} file, making this your default every
time you start up \code{R}}

\item{verbose}{\code{logical} of whether or not to display output}
}
\value{
Invisbly, the URL that has been stored
}
\description{
Sets the \code{SymbiotaR2_url} option for you, optionally, by
appending it to your \code{.Rprofile}. Checks whether you've
specified a valid URL that can be reached, and attempts to pull 
a resource from the API, to confirm that the URL does specify 
a Symbiota2 portal.
}
\examples{
\dontrun{
# An example (that doesn't work because it's not a real portal)
SymbiotaR2_setup("http://nonexistent-portal.com/api/")
# Trying to save a non-existence portal
SymbiotaR2_setup("http://nonexistent-portal.com/api/", TRUE)
}
}
\author{
Will Pearse
}
