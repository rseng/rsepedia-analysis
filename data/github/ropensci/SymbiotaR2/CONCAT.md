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
