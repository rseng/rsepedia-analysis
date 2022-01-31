
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- build with rmarkdown::render("README.Rmd") -->
[infx](https://docs.ropensci.org/infx)
=======================================

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Travis-CI Build Status](https://travis-ci.org/ropensci/infx.svg?branch=master)](https://travis-ci.org/ropensci/infx) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/o6h1088dpuwgvmh4?svg=true)](https://ci.appveyor.com/project/nbenn/infx) [![Coverage status](https://codecov.io/gh/ropensci/infx/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/infx?branch=master) [![rOpenSci status](https://badges.ropensci.org/218_status.svg)](https://github.com/ropensci/onboarding/issues/218)

With discovery of the RNA interference (RNAi) pathway and subsequent work on design, synthesis and delivery of short-interfering RNA molecules (siRNAs), an effective tool was developed for functional genomics. Using such technology, in conjunction with advancements in lab automation, large-scale loss-of-function studies have become experimentally tractable. Fluorescence microscopy imaging is often used as experimental readout, where cellular structures are stained using several fluorophores, yielding large numbers of multi-channel images under varying genetic knockdown conditions. Image processing is typically used to deal with technical artifacts, followed by feature extraction and downstream analysis of the resulting data.

Such experimental set ups can easily generate considerable amounts of data which typically are ingested by specialized data management platforms, such as the Open Biology Information System (openBIS). Screening data collected by [InfectX](http://www.infectx.ch) and [TargetInfectX](https://www.targetinfectx.ch) can be accessed using the presented R client while a browser-based view of the data is available from the [InfectX data browser](http://www.infectx.ch/databrowser). While presented as a tool for InfectC data access, `infx` is not limited to this specific dataset. Any openBIS instance that supports the v1 JSON-RPC API can be accessed using the `infx` client.

Installation
------------

You can install the development version of [infx](https://docs.ropensci.org/infx) from GitHub by running

``` r
source("https://install-github.me/ropensci/infx")
```

Alternatively, if you have the `remotes` package available and are interested in the latest release, you can install from GitHub using `install_github()` as

``` r
# install.packages("remotes")
remotes::install_github("ropensci/infx@*release")
```

InfectX
-------

[InfectX](http://www.infectx.ch) and its successor project [TargetInfectX](https://www.targetinfectx.ch) are large-scale high throughput screening experiments focused on the human infectome of a set of viral and bacterial pathogens. In order to identify host-provided components involved in pathogen entry and host colonization, several RNAi screens were carried out on HeLa cells, using siRNA libraries from vendors including Dharmacon, Quiagen and Ambion. Of the many performed screens, currently the data of kinome-wide screens for five bacterial pathogens (*Bartonella henselae*, *Brucella abortus*, *Listeria monocytogenes*, *Salmonella* typhimurium, and *Shigella flexneri*) and three viruses (Adenovirus, Rhinovirus, and *Vaccinia virus*) is publicly available[1]. Additionally, several genome-wide screens will follow suit in the coming months.

All collected data, including raw imaging data, [CellProfiler](http://cellprofiler.org) derived feature data and infection scoring at single cell resolution, alongside extensive metadata, is hosted by the laboratory information management system [openBIS](https://labnotebook.ch). This R package provides access to the openBIS [JSON-RPC API](https://wiki-bsse.ethz.ch/display/openBISDoc1304/openBIS+JSON+API), enabling listing of data organization objects, searching for and downloading of data sets.

OpenBIS
-------

Only a brief introduction on how to work with openBIS is given here. For more in-depth information on how data is organized in openBIS and how it can be accessed using this package, please refer to the vignette ["Introduction to infx"](https://docs.ropensci.org/infx/articles/infx-intro.html). For an extensive look at what parts of the API are currently implemented and how to extend the package to support further functionality, have a look at the vignettes ["OpenBIS API coverage"](https://docs.ropensci.org/infx/articles/openbis-api.html) and ["JSON object handling"](https://docs.ropensci.org/infx/articles/json-class.html). Documentation of exported functions is available from within the R help system or from [here](https://docs.ropensci.org/infx/reference/index.html).

For every API call, a valid login token is required. Tokens can be created using [`login_openbis()`](https://docs.ropensci.org/infx/reference/login.html) and tested for validity with [`is_token_valid()`](https://docs.ropensci.org/infx/reference/login.html).

``` r
tok <- login_openbis()

is_token_valid(tok)
#> [1] TRUE
```

Using the valid login token, openBIS can now be queried, for example for a list of all projects that are available to the given user, using [`list_projects()`](https://docs.ropensci.org/infx/reference/list_projects.html).

``` r
projects <- list_projects(tok)
print(projects, length = 10L)
#> ┌─█─Project 
#> │ ├─permId = 20130710131815818-2788266 
#> │ ├─spaceCode = INFECTX_PUBLISHED 
#> │ ├─code = ADENO_TEAM 
#> │ ├─description =  
#> │ ├─registrationDetails = █─EntityRegistrationDetails 
#> │ │                       └─... 
#> │ └─id = 39 
#> ├─█─Project 
#> ...
```

Finally, the login token should be destroyed, using [`logout_openbis()`](https://docs.ropensci.org/infx/reference/login.html).

``` r
logout_openbis(tok)
is_token_valid(tok)
#> [1] FALSE
```

While this client has been thoroughly tested with the openBIS instance hosted by InfectX and certain aspects are geared towards high content screening application of openBIS, it is in no way limited to usage with InfectX data. The function [`login_openbis()`](https://docs.ropensci.org/infx/reference/login.html) accepts a `host_url` argument which is stored as `host_url` attribute with the created login token. Any method that issues an API call subsequently uses the login token's `host_url` attribute in order to construct the API endpoint url. As a small example for this functionality, the demo openBIS instance, maintained by the openBIS development team, is queried for available projects.

``` r
tok <- login_openbis(user = "test_observer",
                     pwd = "test_observer",
                     host_url = "https://openbis-eln-lims.ethz.ch")

projects <- list_projects(tok)
print(projects, length = 10L)
#> ┌─█─Project 
#> │ ├─permId = 20150126115738287-33 
#> │ ├─spaceCode = MATERIALS 
#> │ ├─code = BACTERIA 
#> │ ├─description =  
#> │ ├─registrationDetails = █─EntityRegistrationDetails 
#> │ │                       └─... 
#> │ └─id = 3 
#> ├─█─Project 
#> ...

logout_openbis(tok)
```

Acknowledgments
---------------

This work is partially funded by [SystemsX.ch](http://www.systemsx.ch), the Swiss Initiative for Systems Biology via grants 51RT-0\_126008 and 51RTP0\_151029 for the Research and Technology Development (RTD) projects [InfectX](https://infectx.ch) and [TargetInfectX](https://www.targetinfectx.ch) respectively. Further funding is provided by the [Seminar for Statistics](https://www.math.ethz.ch/sfs) at ETH Zurich.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

[1] [*BMC Genomics* 2014 **15**:1162](https://doi.org/10.1186/1471-2164-15-1162)
infx 0.1.0 (2018-03-23)
=========================

### NEW FEATURES

  * ready for rOpenSci onboarding request
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
# Contributing to [infx](https://github.com/ropensci/infx)

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

*  Look at the Travis and AppVeyor build status before and after making changes. The `README` should contain badges for any continuous integration services used by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  Use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the [infx](https://github.com/ropensci/infx) project is released with a [Contributor Code of Conduct](CONDUCT.md). By contributing to this project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
output:
  github_document:
    html_preview: false
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- build with rmarkdown::render("README.Rmd") -->

```{r setup, include = FALSE}
library(infx)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# [infx](https://docs.ropensci.org/infx)

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/infx.svg?branch=master)](https://travis-ci.org/ropensci/infx)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/o6h1088dpuwgvmh4?svg=true)](https://ci.appveyor.com/project/nbenn/infx)
[![Coverage status](https://codecov.io/gh/ropensci/infx/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/infx?branch=master)
[![rOpenSci status](https://badges.ropensci.org/218_status.svg)](https://github.com/ropensci/onboarding/issues/218)

With discovery of the RNA interference (RNAi) pathway and subsequent work on design, synthesis and delivery of short-interfering RNA molecules (siRNAs), an effective tool was developed for functional genomics. Using such technology, in conjunction with advancements in lab automation, large-scale loss-of-function studies have become experimentally tractable. Fluorescence microscopy imaging is often used as experimental readout, where cellular structures are stained using several fluorophores, yielding large numbers of multi-channel images under varying genetic knockdown conditions. Image processing is typically used to deal with technical artifacts, followed by feature extraction and downstream analysis of the resulting data.

Such experimental set ups can easily generate considerable amounts of data which typically are ingested by specialized data management platforms, such as the Open Biology Information System (openBIS). Screening data collected by [InfectX](http://www.infectx.ch) and [TargetInfectX](https://www.targetinfectx.ch) can be accessed using the presented R client while a browser-based view of the data is available from the [InfectX data browser](http://www.infectx.ch/databrowser). While presented as a tool for InfectC data access, `infx` is not limited to this specific dataset. Any openBIS instance that supports the v1 JSON-RPC API can be accessed using the `infx` client.

## Installation

You can install the development version of [infx](https://docs.ropensci.org/infx) from GitHub by running

```{r gh-dev, eval = FALSE}
source("https://install-github.me/ropensci/infx")
```

Alternatively, if you have the `remotes` package available and are interested in the latest release, you can install from GitHub using `install_github()` as

```{r gh-rel, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/infx@*release")
```

## InfectX

[InfectX](http://www.infectx.ch) and its successor project [TargetInfectX](https://www.targetinfectx.ch) are large-scale high throughput screening experiments focused on the human infectome of a set of viral and bacterial pathogens. In order to identify host-provided components involved in pathogen entry and host colonization, several RNAi screens were carried out on HeLa cells, using siRNA libraries from vendors including Dharmacon, Quiagen and Ambion. Of the many performed screens, currently the data of kinome-wide screens for five bacterial pathogens (*Bartonella henselae*, *Brucella abortus*, *Listeria monocytogenes*, *Salmonella* typhimurium, and *Shigella flexneri*) and three viruses (Adenovirus, Rhinovirus, and *Vaccinia virus*) is publicly available^[[*BMC Genomics* 2014 **15**:1162](https://doi.org/10.1186/1471-2164-15-1162)]. Additionally, several genome-wide screens will follow suit in the coming months.

All collected data, including raw imaging data, [CellProfiler](http://cellprofiler.org) derived feature data and infection scoring at single cell resolution, alongside extensive metadata, is hosted by the laboratory information management system [openBIS](https://labnotebook.ch). This R package provides access to the openBIS [JSON-RPC API](https://wiki-bsse.ethz.ch/display/openBISDoc1304/openBIS+JSON+API), enabling listing of data organization objects, searching for and downloading of data sets.

## OpenBIS

Only a brief introduction on how to work with openBIS is given here. For more in-depth information on how data is organized in openBIS and how it can be accessed using this package, please refer to the vignette ["Introduction to infx"](https://docs.ropensci.org/infx/articles/infx-intro.html). For an extensive look at what parts of the API are currently implemented and how to extend the package to support further functionality, have a look at the vignettes ["OpenBIS API coverage"](https://docs.ropensci.org/infx/articles/openbis-api.html) and ["JSON object handling"](https://docs.ropensci.org/infx/articles/json-class.html). Documentation of exported functions is available from within the R help system or from [here](https://docs.ropensci.org/infx/reference/index.html).

For every API call, a valid login token is required. Tokens can be created using [`login_openbis()`](https://docs.ropensci.org/infx/reference/login.html) and tested for validity with [`is_token_valid()`](https://docs.ropensci.org/infx/reference/login.html).

```{r login}
tok <- login_openbis()

is_token_valid(tok)
```

Using the valid login token, openBIS can now be queried, for example for a list of all projects that are available to the given user, using [`list_projects()`](https://docs.ropensci.org/infx/reference/list_projects.html).

```{r projects}
projects <- list_projects(tok)
print(projects, length = 10L)
```

Finally, the login token should be destroyed, using [`logout_openbis()`](https://docs.ropensci.org/infx/reference/login.html).

```{r logout}
logout_openbis(tok)
is_token_valid(tok)
```

While this client has been thoroughly tested with the openBIS instance hosted
by InfectX and certain aspects are geared towards high content screening application of openBIS, it is in no way limited to usage with InfectX data. The function [`login_openbis()`](https://docs.ropensci.org/infx/reference/login.html) accepts a `host_url` argument which is stored as `host_url` attribute with the created login token. Any method that issues an API call subsequently uses the login token's `host_url` attribute in order to construct the API endpoint url. As a small example for this functionality, the demo openBIS instance, maintained by the openBIS development team, is queried for available projects. 

```{r other-openbis}
tok <- login_openbis(user = "test_observer",
                     pwd = "test_observer",
                     host_url = "https://openbis-eln-lims.ethz.ch")

projects <- list_projects(tok)
print(projects, length = 10L)

logout_openbis(tok)
```

## Acknowledgments

This work is partially funded by [SystemsX.ch](http://www.systemsx.ch), the Swiss Initiative for Systems Biology via grants 51RT-0_126008 and 51RTP0_151029 for the Research and Technology Development (RTD) projects [InfectX](https://infectx.ch) and [TargetInfectX](https://www.targetinfectx.ch) respectively. Further funding is provided by the [Seminar for Statistics](https://www.math.ethz.ch/sfs) at ETH Zurich.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "JSON object handling"
author: "Nicolas Bennett"
date: "`r Sys.Date()`"
output:
  html_vignette:
    self_contained: no
vignette: >
  %\VignetteIndexEntry{JSON objects}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(infx)
```

The openBIS JSON-RPC API is powered by the Jackson JSON processor on the server side. As such, object type information is communicated via `@type` fields. Take as an example the class `Dog` defined as

```java
public class Dog {
  private String name;
  private String breed;
  // setters and getters
}
```

An instance of `Dog` can be represented as the following JSON object

```json
{
  "@type": "Dog",
  "name": "Rufus",
  "breed": "english shepherd"
}
```

where the `@type` filed is used by the deserializer to infer what java object should result. Furthermore, on the client side, class information is used for S3 method dispatch. This document will illustrate how typed JSON objects returned by openBIS are converted into S3 classes, how these classes can be manipulated in R and subsequently used for further openBIS queries.

## Creating `json_class` objects

Again a feature of the Jackson serializer, `@id` fields are generated for all objects by an `ObjectIdGenerator` (an `IntSequenceGenerator` to be more specific), which can be used to represent relationships among objects. These object ids are employed whenever the same object is returned multiple times. This happens for example when several objects representing wells on the same plate are returned. The first `WellIdentifier` object will contain a `PlateIdentifier` object and all subsequent wells identify their parent plate via a reference to this `PlateIdentifier` object. A simplified example, demonstrating this structure is shown below.

```{r resp}
response <- list(list(`@type` = "WellIdentifier",
                      `@id` = 1L,
                      row = "A",
                      col = 1L,
                      plate = list(`@type` = "PlateIdentifier",
                                   `@id` = 2L,
                                   barcode = "abcd-123")),
                 list(`@type` = "WellIdentifier",
                      `@id` = 3L,
                      row = "A",
                      col = 2L,
                      plate = 2L))
```

Turning the list `response` into a `json_class` S3 class can be done, using the function `as_json_class()` (or its alias `as.json_class()`). This action is applied recursively, meaning that all sub-lists, containing an `@type` field are turned into `json_class` objects as well.

```{r as-class}
json_class <- as_json_class(response)
str(json_class)
```

In order to resolve object references encoded by `@id` fields, the function `resolve_references()`, which also acts recursively on its input, can be used. This will duplicate all occurrences of referenced objects, introducing a memory penalty, but in turn making objects self-contained. After resolving references, the list of well objects can be subsetted and the second `WellIdentifier` object can be used for a further query on its own.

```{r ref}
json_class <- resolve_references(json_class)
str(json_class)
```

In case the user wants to create a `json_class` object, the constructor `json_class()` is available. It can be used as follows, e.g. for the re-creation of the above `json_class` object.

```{r construct}
construct <- json_class(row = "A", col = 1L,
                        plate = json_class(barcode = "abcd-123",
                                           class = "PlateIdentifier"),
                        class = "WellIdentifier")
identical(construct, json_class[[1L]])
```

The effect of `as_json_class()` can be reversed by `rm_json_class()` or `as_list()` (as well as its alias `as.list()`). The two functions differ in default behavior however. Where `rm_json_class()` removes S3 classes and writes class information into the respective `@type` fields, `as_list()`, by default simple returns its input. This somewhat odd choice is owed to the circumstance that iterating though a `json_class` object with `lapply()` or `sapply()` calls `as.list()` on its first argument.

```{r destruct}
identical(response,
          rm_json_class(json_class))
identical(as.list(json_class, keep_asis = FALSE),
          rm_json_class(json_class))
```

In addition to functions for creating and destroying `json_class` objects, several utility functions for `json_class` are provided as well. `is_json_class()` (and its alias `is.json_class()`) tests whether the object in question is a list inheriting the class attribute `json_class` and has at least one more class attribute in front of `json_class`, which is expected in last place (within the class vector). Similarly, `check_json_class()` can be used to recursively test a list structure to make sure that every node inheriting the `json_class` class attribute is a properly formed `json_class` object. Finally, the two functions `get_subclass()` and `has_subclass()` can be used to extract the sub-class and test whether a given `json_class` object has a specific sub-class.

```{r test}
test <- structure(list(a = structure(list(b = "c"),
                                     class = c("foo", "json_class")),
                       d = structure(c(e = "f"),
                                     class = c("bar", "json_class"))),
                  class = c("foobar", "json_class"))

is_json_class(test)
check_json_class(test)
is_json_class(test$a)
is_json_class(test$d)

has_subclass(test, "foobar")
get_subclass(test)
```

In the above example, `is_json_class(test$d)` returns `FALSE` because it is an object that is not composed as a list with attributes, but as a vector with attributes. This is also the reason why the recursive execution of `is_json_class()` to all contained objects that inherit from `json_class` in `check_json_class(test)` returns `FALSE`.

The two base R generic functions `print()` and `c()` have `json_class`-specific methods implemented. Combining several `json_class` object, using `c()` will yield a `json_vec` object, which is described in the following section. Printing is recursive and recursion depth can be controlled using the argument `depth`. Further printing options are `width`, `length` and a logical flag `fancy` for enabling fancy printing (console output is colored and UTF box characters are used for creating a tree structure instead of ASCII characters).

```{r base-generics}
json_class[[1]]
print(json_class[[1]], depth = 2L)
print(json_class[[1]], depth = 2L, length = 4L)
print(json_class[[1]], depth = 2L, fancy = FALSE)
```

## Using vectors of `json_class` objects

Now that class information of objects fetched from openBIS is available in R, this can be used for method dispatch of S3 generic functions. Assume we have a generic function `list_datasets()`. We could then implement an object-specific method for objects of type `sample` as `list_datasets.sample()` and one for objects of type `experiment` as `list_datasets.experiment()`. Depending on the class of the object passed to `list_datasets()`, datasets for an experiment or for a sample will be listed.

There is an issue with this approach though: listing for example datasets associated with multiple experiments. A straightforward approach could be simply iterating over the list of experiments and for each one, issuing a separate request to openBIS. A more efficient way of doing this would be to query openBIS with a single list of several experiment objects. This, however defeats S3 method dispatch, as the object on which dispatch occurs is no longer `experiment` but `list` (containing several `experiment`s). To work around this, `json_vec` objects are used.

The `json_vec` class wraps around a list of `json_class` and it serves to bring the common sub-class of all child objects to the surface for using method dispatch on. Instantiation of `json_vec` objects is possible in several ways: using the `json_vec()` constructor, by coercing a list structure with `as_json_vec()` or by combining several `json_class` objects using `base::c()`. The reverse action can be performed, using `base::as.list()`. The `json_vec` constructor, as well as functions for coercing to `json_vec` accept a `simplify` argument. When set to `TRUE` (default is `FALSE`), `json_vec` objects of length 1 are simplified to `json_class` objects. 

```{r json-vec}
a <- json_class("a", class = "foo")
b <- json_class("b", class = "foo")

foo_vec <- json_vec(a, b)

str(foo_vec)

identical(foo_vec, as_json_vec(list(a, b)))
identical(foo_vec, c(a, b))

identical(as.list(foo_vec), list(a, b))

identical(a, as_json_vec(list(a), simplify = TRUE))
```

Utility functions available for `json_vec` objects include `has_common_subclass()` and `get_subclass()` to check whether all entries in a list are `json_class` objects with the same sub-class and to extract the common sub-class, as well as `is_json_vec()` (and its alias `is.json_vec()`), which can be used to check whether

  * all child elements are of the same sub-class
  * all child elements are properly formed `json_class` objects
  * the `json_vec` class attribute is in last position of the class vector
  * the remaining class attributes are equal to the common sub-class of the
    children.

```{r vec-utils}
has_common_subclass(list(a, b))

get_subclass(list(a, b))
get_subclass(foo_vec)

is_json_vec(foo_vec)
is_json_vec(list(a, b))
```

Finally, the base R generics for which `json_vec`-specific methods are provided include `print()`, `c()`, as well as accessors and assignment functions.

```{r vec-base}
foo_vec
foo_vec[1]

c(foo_vec[2], foo_vec[1])

class(foo_vec[1])
class(foo_vec[[1]])
```

As shown in the above code block, subsetting a `json_vec` object preserves all class attributes. This is in analogy to base vector objects, where subsetting, for example a character vector again yields a character vector. The single element accessor `[[`, however yields an element of the `json_vec` object, which should not come as surprise given the list nature of `json_vec` objects.
---
title: "Introduction to infx"
author: "Nicolas Bennett"
date: "`r Sys.Date()`"
output:
  html_vignette:
    self_contained: no
vignette: >
  %\VignetteIndexEntry{Introduction to infx}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(infx)
library(tibble)
library(magick)
library(ggplot2)
library(RColorBrewer)
```

OpenBIS (Open Biology Information System) is a laboratory information management system designed for robust data management of large-scale experiments in biological sciences. As storage infrastructure it is therefore well suited for the needs of image-based high throughput screening (HTS) as performed by the InfectX consortium. For data access, JSON-RPC services are provided by openBIS, which can be called from the presented client package `infx`.

This document gives a short introduction to some basic organizational concepts of openBIS with a focus on aspects relevant to HTS and provides some examples of how the `infx` can be used to access various types of data generated by the InfectX experiments. For more information, general openBIS documentation is available [here](https://wiki-bsse.ethz.ch/display/openBISDoc1304/openBIS) and documentation specific to the JSON-RPC API can be accessed from [here](https://wiki-bsse.ethz.ch/display/openBISDoc1304/openBIS+JSON+API). It might help to have a look at the browser-based web GUI available [here](https://infectx.biozentrum.unibas.ch/openbis) alongside this document to help understanding the presented ideas.

## Organizational concepts in openBIS 

An organizational entity central to the openBIS storage logic is an **experiment**. In the context of InfectX, an experiment is a single screen, meaning the combination of a compound library and an experimental condition provided by the presence of a pathogen. For example in the experiment `ADENO-AU-K1`, a kinome-wide siRNA library by Ambion (Silencer Select) was applied in unpooled fashion (3 siRNAs per gene), alongside exposure to the pathogen Adenovirus (cf. the `properties` field of `Experiment` objects as shown below).

Experiments are grouped into **projects** (one per pathogen in the case of InfectX), which in turn are grouped into **spaces** (an unimportant hierarchical level for InfectX). Projects can be listed using `list_projects()` and experiments with `list_experiments()`.

```{r exp-proj}
token <- login_openbis()

projects <- list_projects(token)

print(projects, length = 10L)
length(projects)

adeno_exps <- list_experiments(token, projects[[1L]])

print(adeno_exps, length = 15L)
length(adeno_exps)

str(adeno_exps[[1L]]$properties)
```

In order to access the API, a login token has to be created. Using this token, all available projects are listed using `list_projects()` and all experiments
corresponding to a project are listed with `list_experiments()`. As mentioned, some information on the individual experiments is available in the `properties` entry of `Experiment` objects.

Experiments in high-throughput screening are typically carried out on microtiter plates which lends itself to a natural way of sub-dividing individual experiments. All InfectX screens were performed on 384 well plates, composed of 16 rows (A through P) and 24 columns (1 through 24) and each plate can be uniquely identified by a barcode. The functions for listing plates and wells are `list_plates()` and `list_wells()`, respectively. The following example shows how for a single experiment, all associated plates and for a single plate, all contained wells can be retrieved.


```{r plate-well}
plates <- list_plates(token, adeno_exps[[1L]])

print(plates, length = 15L)
length(plates)

wells <- list_wells(token, plates[[2]])

print(wells, length = 15L, depth = 2L)
length(wells)
```

In terms of openBIS entities, both a plate and a well are considered **samples**. A sample is described as follows by the [openBIS user documentation](https://wiki-bsse.ethz.ch/display/openBISDoc1304/openBIS+User+Documentation#openBISUserDocumentation-Introduction):

> A sample refers to any object that has been observed, measured, or compared to another. It must be uniquely identifiable, which means that any two samples must be distinguishable from one another. Please note that different use cases may use the term "sample" with slightly different meanings, dependent upon the context and utility. ... [T]he term "sample" could [for example] refer to an individual well in a multi-titer plate containing cells of different phenotypes.

The function `list_samples()` retrieves `Sample` objects, which generalize `PlateIdentifier` and `WellIdentifier` objects. As with many other `infx` functions, `list_samples()` is an S3 generic function. If dispatch occurs on an `Experiment` object, the set of plate samples belonging to an experiment is fetched, as with `list_plates()`. Well samples per plate cannot be directly listed as was the case with `list_wells()`. However `list_samples()` dispatched on a set of `WellIdentifier` objects will return the corresponding well samples.

```{r sample}
plate_samp <- list_samples(token, adeno_exps[[1L]])

print(plate_samp, length = 20L)
length(plate_samp)

wells_samp <- list_samples(token, wells[1L:2L])

print(wells_samp, length = 20L)
```

The sample type is encoded in the `sampleTypeCode` field of each `Sample` object and an exhaustive list of available sample types can be shown using `list_sample_types()`.

A further important organizational concept of openBIS is that of a **data set**. On this entity, the [openBIS user documentation](https://wiki-bsse.ethz.ch/display/openBISDoc1304/openBIS+User+Documentation#openBISUserDocumentation-Introduction) notes the following:

> A data set is the computer's representation of a series of sample measurements, or the results of computational processing derived from those measurements. As with samples and experiments, data sets also have specific data set types to better handle searching and analysis needs.

Essentially, a data set represents a collection of files associated with a sample. Furthermore data sets may have (multiple) parent/child relationships among each other, to indicate one data set being derived of another. Retrieving all data sets belonging to a plate can be achieved with calling `list_datasets()` on plate sample objects.

```{r data-set}
data_sets <- list_datasets(token, plate_samp[[2L]])

print(data_sets, length = 30L)
length(data_sets)

unique(get_field(data_sets, "dataSetTypeCode"))

list_datasets(token, wells_samp[[1L]])
```

Several different types of data sets (possibly in multiple versions) are typically associated with a plate. Some of the more interesting data set types are

* `HCS_IMAGE_CONTAINER_RAW`: raw imaging data, 6-9 images per well each available for 3-4 imaging channels
* `HCS_IMAGE_CONTAINER_SEGMENTATION`: image overlays for segmenting images into cells, nuclei, etc.
* `HCS_ANALYSIS_IMAGE_ACQUISITION_METADATA`: microscope image meta data and settings
* `HCS_ANALYSIS_CELL_FEATURES_CC_MAT`: CellProfiler feature data at single cell resolution
* `HCS_ANALYSIS_CELL_CLASSIFICATIONS_MAT`: decision tree-based infection scoring data

The way openBIS is set up for InfectX, data sets are only available on the plate sample level and not per well, as is demonstrated by passing a well sample object to `list_datasets()`, which returns an empty list.

## Searching in openBIS 

Search queries for openBIS are constructed with `search_criteria()` and the search is executed by calling `search_openbis()`. The function `search_criteria()` instantiates a `SearchCriteria` object which consists of a set of match clauses combined with either an `any` or `all` operator. Nesting of `SearchCriteria` objects is possibly by supplying a `SearchCriteria` object as `sub_criteria` argument to a call to `search_criteria()`, in turn creating the enclosing `SearchCriteria` object.

Five different types of match clauses can be constructed:

* `PropertyMatchClause`: A `MatchClause` for checking that a property equals a desired value.
* `AnyPropertyMatchClause`: A `MatchClause` for checking that any of the properties equals a desired value.
* `AnyFieldMatchClause`: A `MatchClause` for checking that any of the properties or attributes equals a desired value.
* `AttributeMatchClause`: A `MatchClause` for checking that an attribute equals a desired value.
* `TimeAttributeMatchClause`: A `MatchClause` for comparing a time attribute to a specified value.

For every match clause, a desired value has to be supplied, as well as a comparison mode which can either be `eq` (equal to), `lte` (less than or equal to) or `gte` (greater than or equal to). Additionally, for a `PropertyMatchClause`, a property code has to be specified (possibilities can be enumerated with `list_property_types()`), for an `AttributeMatchClause`, an attribute^[possible values are `code`, `type`, `perm_id`, `space`, `project`, `project_perm_id`, `metaproject`, `registrator_user_id`, `registrator_first_name`, `registrator_last_name`, `registrator_email`, `modifier_user_id`, `modifier_first_name`, `modifier_last_name` or `modifier_email`] and for a `TimeAttributeMatchClause`, a time attribute (either `registration_date` or `modification_date`).

```{r simple-search}
amb_kin <- search_criteria(property_clause("library", "Ambion"),
                           property_clause("geneset", "Kinome"),
                           operator = "all")

ak_exps <- search_openbis(token, amb_kin,
                          target_object = "experiment")

print(ak_exps, length = 15L)
get_field(ak_exps, "code")
```
In this example, openBIS is queried for all experiments that involve kinome-wide screens with Ambion libraries. First, a `SearchCriteria` object is created containing two property match clauses that both have to be met simultaneously. This `SearchCriteria` object is then passed to `search_openbis()` along with the specification of a target type which can be either `data_set`, `experiment`, `material` or `sample`.

```{r material-search}
mtor_mat <- search_openbis(
  token,
  search_criteria(property_clause("gene_symbol", "MTOR")),
  "material"
)

print(mtor_mat, depth = 2L)

well_refs <- list_references(token, mtor_mat, ak_exps[[1L]])
print(well_refs, length = 15L)
```

A second example for a query, this time for a material object is given above. The search is constructed such that the returned object represents a compound targeting the gene [`MTOR`](https://www.ncbi.nlm.nih.gov/gene/2475). This `MaterialGeneric` object then can be used to list wells on plates, involving this compound using the function `list_references()`. The inverse of this, where for a given plate object all used materials are listed with associated wells, can be achieved using the function `list_material()`.

## Retrieving openBIS data resources

Three different types of data resources are available from openBIS: The most straightforward is `files`. As explained above, each data set contains a set of files, for each of which a download url can be created using `list_download_urls()`. As this openBIS instance is hosting image-based HTS data, a second available data resource is `images`. Raw images can be retrieved as files in a `HCS_IMAGE_CONTAINER_RAW` data set but in addition to that, openBIS can be queried for specific images, instead of the plate-wise access provided by the data set route, and is able to serve transformations of raw images. A final type of data resource is `features`. This data is also available as files in a data set but similar to images is treated specially by openBIS in order to allow fine-grained queries.

### File download

The following example demonstrates how InfectX single cell feature data, calculated by CellProfiler, can be accessed. First a search for data sets of type `HCS_ANALYSIS_CELL_FEATURES_CC_MAT` is carried out. This search is limited to the `ADENO-AU-K1` experiment, using a `search_sub_criteria` object. One of the resulting data sets is then passed to `fetch_files()` together with a regular expression to filter the list of available files (several hundred feature files are typically available for such data sets). The function `read_mat_files()` is passed as `reader` argument and reads the binary Matlab files using `R.matlab::readMat()`.

```{r fetch-counts}
adeno_au_sub <- search_sub_criteria(
  search_criteria(
    property_clause("pathogen", "Adenovirus"),
    property_clause("library", "Ambion"),
    property_clause("geneset", "Kinome"),
    property_clause("replicate", 1L)
  ),
  type = "experiment"
)

adeno_au_mat <- search_criteria(
  attribute_clause("type", "HCS_ANALYSIS_CELL_FEATURES_CC_MAT"),
  sub_criteria = adeno_au_sub
)

cell_ds <- search_openbis(token, adeno_au_mat,
                          target_object = "data_set")

print(cell_ds, length = 30L)
length(cell_ds)

dat <- fetch_files(token, cell_ds[[1L]],
                   file_regex = "Image\\.Count_",
                   reader = read_mat_files)

names(dat) <- sapply(dat, attr, "feature")
dat <- lapply(dat, as.integer)

tibble::as_tibble(lapply(dat, unlist))
```

For each file, `read_mat_files()` will return a list with one entry per imaging site. For this data set, there are 9 imaging sites per well which yields 3456 sites for the entire plate. Additionally, `fetch_files()` returns a list structure per request, containing information on which request corresponds to which data set and file. This is necessary because `fetch_files()` could be called on several data sets at once, each returning multiple files.

As a second example, area features are requested. For the given screen, area measurements are available for the three object types `PeriNuclei`, `Nuclei` and `Cells`. Unlike in the previous example, where a scalar corresponds to each well, here the variables are vector-valued per well. Therefore we need to create a column `Well`, indicating which rows correspond to which wells. Well indices are linearized in row-major fashion with respect to the plate layout.

```{r fetch-areas}
dat <- fetch_files(token, cell_ds[[1L]],
                   file_regex = "AreaShape_Area",
                   reader = read_mat_files)

attributes(dat[[1L]])

names(dat) <- paste0("Area_", sapply(dat, attr, "object"))
well_names <- paste0(rep(LETTERS[1L:16L], each = 24L), rep(1L:24L, 16L))
well_names <- rep(rep(well_names, each = 9L), sapply(dat[[1L]], length))

dat <- lapply(dat, unlist)

tibble::as_tibble(c(list(Well = well_names), lapply(dat, as.integer)))
```

The resulting data matrix $X$ holds all measurements of a plate for the selected features and is structured as

$$X = \begin{bmatrix}
X_{G_1} \\
X_{G_2} \\
... \\
X_{G_m}
\end{bmatrix}
$$

where groups of rows $X_{G_i}$ are $n_i \times p$ matrices holding $p$ features as columns corresponding to $n_i$ single cell measurements under knock-down of gene $G_i$.

Several attributes are set for each requested file. The `read_mat_files()` function extracts object (what type of CellProfiler object the feature was calculated on) and feature information (what kind of CellProfiler measurement was performed) from the read file. In addition, `fetch_files()` stores request information such as dataset and file in order for the user to match responses with requests. 

### Image access

In order to fetch images, again first a search is constructed. Re-using the previous `search_sub_criteria`, the search is targeted at sample objects of type `PLATE`, as image data sets are connected to plates. To find the appropriate data set, the function `list_references()` may be used and since the current target is fetching raw image data, the `type` argument of `list_references()` can be left at default value. The returned `ImageDatasetReference` is then passed to `list_image_metadata()` for some additional information on the image data set, mainly the available channels.

To narrow down the requested set of images, `list_references()` is called again, this time on the `ImageDatasetReference` object and in conjunction with a `WellPosition` object. The returned set of `PlateImageReference` objects precisely specify a single image by containing information on image data set, well position, image tile and imaging channel. The `ImageDatasetReference` corresponding to the tile with index 0 is passed to `fetch_images()`, yielding a single image. As `fetch_images()` can be called on several objects specifying images, each request contributes an entry to the resulting list with some meta data attached as attributes.

```{r fetch-images, fig.align = "center"}
adeno_au_samp <- search_criteria(
  attribute_clause("type", "PLATE"),
  sub_criteria = adeno_au_sub
)

samples <- search_openbis(token, adeno_au_samp,
                          target_object = "sample")

raw_ref <- list_references(token, samples[[2L]])

img_meta <- list_image_metadata(token, raw_ref)
print(img_meta)

well_raw <- list_references(token, raw_ref,
                            wells = well_pos(name = "A2"),
                            channel = img_meta[["channelCodes"]][[1L]])

print(well_raw, depth = 2L, length = 15L)

raw_img <- fetch_images(token, well_raw[[2L]],
                        image_size = json_class(width = 600L,
                                                height = 600L,
                                                class = "ImageSize"))

attributes(raw_img[[1L]])
print(raw_img[[1L]])
```

As further illustration of the capabilities of the openBIS API, the following example combines the previously fetched image with a segmentation mask for cells. The same sample object from above is again passed to `list_references()` but this time an `ImageDatasetReference` object corresponding to an image segmentation dataset is retrieved. Again using `list_image_metadata()`, the available channels are listed and using this information, a request for the image segmentation masks for the desired well and image tile is issued.

```{r fetch-mask, fig.align = "center"}
segm_ref <- list_references(token, samples[[2L]],
                            type = "segmentation")

list_image_metadata(token, segm_ref)

well_segm <- list_references(token, segm_ref,
                             wells = well_pos(name = "A2"),
                             channel = "CELLS__CY5_")

segm_img <- fetch_images(token, well_segm[[2L]],
                         image_size = json_class(width = 600L,
                                                 height = 600L,
                                                 class = "ImageSize"))

cells <- c(raw_img[[1L]], magick::image_transparent(segm_img[[1L]], "black"))
print(magick::image_mosaic(cells))
```

Note that in order to create a usable segmentation mask, the black background of the segmentation data set has to be made transparent using `magick::image_transparent()`. Only after applying this transformation, can the segmentation mask be placed on top of the microscopy image.

### OpenBIS feature data

The third type of openBIS data resources are files treated as feature data sets. As such, these can be queried similarly to images. Again using the plate sample object from the previous search query, `list_references()` with type specification set to `feature` will return `FeatureVectorDatasetReference` objects. Several types of feature data sets are available, each of which may contain several features. In order to list all features contained in a feature data set, the function `list_features()` may be used. Finally, `fetch_features()` will return the requested feature data for the specified plates or wells as `FeatureVectorDataset` object. This in turn contains a list of `FeatureVector` objects, each of which holds feature information for a single well.

```{r feat-data}
feat_ref <- list_references(token, samples[[2L]],
                            type = "feature")
unique(get_field(feat_ref, "dataSetType"))

print(list_features(token, feat_ref[[1L]]), length = 10L)

cell_count <- fetch_features(token, feat_ref[[1L]], "COUNT_CELLS")

print(cell_count)
print(cell_count[["featureVectors"]], depth = 2L, length = 10L)
```

Whenever feature data for an entire plate is fetched it might be more efficient to simply download the associated data set file, especially if several of the contained features are of interest. Such files are `.csv` formatted tables with columns corresponding to features and rows to wells. A convenient aspect of feature data is however that it can be queried per well. Passing the `FeatureVectorDatasetReference` object to `list_references()` alongside a set of `WellPosition` objects will return `FeatureVectorDatasetWellReference` objects which can be used to retrieve subsetted feature information using `fetch_features()`.

For illustration purposes, a heatmap of per-well cell counts, drawn with `ggplot::ggplot()`, is shown below. This is akin to the plate heatmaps that are shown in the openBIS web GUI.

```{r heatmap, echo = FALSE, fig.align = "center"}
wells <- get_field(cell_count[["featureVectors"]], "wellPosition")
heatmap <- tibble(
  WellRow = as.factor(LETTERS[get_field(wells, "wellRow")]),
  WellCol = as.factor(as.integer(get_field(wells, "wellColumn"))),
  CellCount = as.integer(get_field(cell_count[["featureVectors"]],
                                   "values"))
)

ggplot(data = heatmap) +
  geom_tile(aes(x = WellCol, y = WellRow, fill = CellCount)) +
  scale_x_discrete(name = "", position = "top") +
  scale_y_discrete(name = "", limits = rev(levels(heatmap$WellRow))) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank()) +
  ggtitle("Plate heatmap of cell count") +
  scale_fill_gradientn(colours = brewer.pal(9, "OrRd"),
                       name = "Cell count") +
  coord_fixed()
```

For the above example, a feature data set of type `HCS_ANALYSIS_WELL_RESULTS_SUMMARIES` was chosen. Such data sets contain per-well aggregated results from CellProfiler analysis, including object counts, mean area measurements and mean intensity measurements. Other feature data sets hold summarized information on infection scoring, quality metrics such as focus scores, underexposure/overexposure indicators, dynamic range analysis or image acquisition meta data.

```{r logout, echo = FALSE}
logout_openbis(token)
```
---
title: "OpenBIS API coverage"
author: "Nicolas Bennett"
date: "`r Sys.Date()`"
output:
  html_vignette:
    self_contained: no
vignette: >
  %\VignetteIndexEntry{OpenBIS API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(infx)
library(rvest)
library(knitr)
```

```{r api-calls, include = FALSE}
get_openbis_items <- function(x) {
  sections <- tools:::RdTags(x) == "\\section"
  if (sum(sections) > 0L)
    sections[sections] <- sapply(x[sections], function(sec) {
      tolower(tools:::.Rd_get_text(sec[[1L]])) == "openbis"
    })

  if (sum(sections) == 1L) {
    obis <- x[[which(sections)]][[2L]]
    is_itemize <- tools:::RdTags(obis) == "\\itemize"
    unlist(
      lapply(obis[is_itemize], function(y) {
        url_hits <- tools:::RdTags(y) == "\\href"
        if (any(url_hits))
          lapply(y[url_hits], sapply, as.character)
        else
          NULL
      }),
      recursive = FALSE
    )
  } else if (sum(sections) > 1L)
    stop("expecting one or zero openbis sections")
  else
    NULL
}

pkg <- "infx"
rd_db <- tools:::fetchRdDB(file.path(find.package(pkg), "help", pkg))
api_calls <- unlist(lapply(rd_db, get_openbis_items), recursive = FALSE)

urls <- sapply(api_calls, `[[`, 1L)
api_calls <- strsplit(sapply(api_calls, `[[`, 2L), ":")
```

```{r api-info, include = FALSE}
info <- lapply(sort(unique(urls)), function(x) {
  docs <- read_html(x)
  api <- html_nodes(docs,
                    xpath = "/html/body/div[4]/div[2]/ul/li/ul[2]/li/table")

  tab <- apply(html_table(api)[[1]], 1, function(y) {
    y <- y[2]
    fun_name <- sub("\\($", "", regmatches(y, regexpr("^.+?\\(", y)))
    desc <- sub("\\)\n", "", regmatches(y, regexpr("\\)\n.+$", y)))

    if (length(desc) && grepl("^Deprecated", desc))
      return(NULL)

    found <- fun_name %in% sapply(api_calls[urls == x], `[`, 2L)

    list(`Method name` = fun_name,
         Status = ifelse(found, "implemented", "skipped"),
         Description = sub("IDssServiceRpcScreening\\.",
                           "IDssServiceRpcScreening ",
                           sub("ch\\..+\\.dto\\.", "",
                               gsub("\n", "",
                                    sub("^/\\*\\* ", "", desc)))))
  })

  title <- html_text(html_nodes(docs, xpath = "/html/body/div[3]/h2"))
  desc <- html_nodes(docs, xpath = "/html/body/div[4]/div[1]/ul/li/div")
  desc <- strsplit(html_text(desc), "\\.")[[1]]
  desc <- paste0(desc[1], ". More information available [here](", x, ").\n")
  tabl <- do.call(rbind, tab)
  link <- paste0("[", sub("^Interface ", "", title), "](", x, ")")

  list(title = title,
       description = desc,
       table = tabl,
       status = round(mean(tabl[, 2] == "implemented") * 100),
       link = link)
})
```

This vignette serves to document which sections of the openBIS API are implemented by the presented client package and which functionality was omitted. Furthermore, the basic mechanisms for making requests to openBIS are explained such that a user who wishes to add some of the omitted API functions has some information on how to extend `infx`. Some general documentation on the API, offered by the openBIS developers is available [here](https://wiki-bsse.ethz.ch/display/openBISDoc1304/openBIS+JSON+API).

## Creating a REST call

The basic functionality powering all Representational state transfer (REST) calls of `infx` is provided by `do_requests_serial()` and `do_requests_parallel()`. Both functions take roughly the same input and produce identical results, but differ in how curl is called. The former function makes requests in a sequential fashion, iteratively calling `curl::curl_fetch_memory()` and the latter performs asynchronous requests, using `curl::curl_fetch_multi()` to assemble all requests and `curl::multi_run()` to carry them out in asynchronous parallel fashion.

Much of the behavior of `do_requests_serial()` and `do_requests_parallel()` can be customized using three functions passed as arguments `create_handle`, `check` and `finally`. Additionally, `do_requests_*()` takes a vector of `urls` and an optional object `bodies` which is expected to be of the same length as `urls`. Urls can either be passed as a character vector or a list of unevaluated function calls which will be evaluated using `base::eval()` shortly before being used. Both functions support retrying failed requests up to `n_try` times.

For constructing the requests, both `urls` and `bodies` are iterated together. First, the function passed as `create_handle`, which receives as input the current entry of the `bodies` object is used to construct a curl handle using `curl::create_handle()`. Together with the current URL entry, `curl::curl_fetch_*()` is called and the resulting object is passed as first argument to the `check` function, which receives as second argument the current entry of the `bodies` object. The function passed as `check` argument should make sure the request completed successfully and in case of failure return a `simpleError` object, created by `base::simpleError()`. In case of success, it should return the part of the curl response object that is of further interest (most likely the `content` entry). Next the `finally` function is called on the object returned by `check` to do some final processing (e.g. parsing JSON or reading a binary file). If the `check` function signaled a failure and the number of available tries as specified by `n_try` is not used up, the request is made again. If the allowed number of tries is exceeded, a warning is issued.

Asynchronous requests are implemented by first adding `n_con` curl handles to a new multi handle and starting the downloads by calling `curl::multi_run()`. For each successful request, a new handle is added to the pool using the `done` callback function of `curl::multi_add()`. In conjunction with passing urls as unevaluated function calls, this helps with urls that have a limited lifetime in that the URL is only created right before being consumed. Instead of adding all requests at the same time and letting curl handle queuing, only `n_con` requests are handled by curl at any given time.

A very basic use of `do_requests_serial()` is shown in the following example. The default function of the `create_handle` argument creates a new clean curl handle with no options set. The default function of the `check` argument makes sure that the status code is equal to 200 and returns the `content` entry of the list returned by `curl::curl_fetch_memory()`. The default function for `finally` is `base::identity()`. In order to receive a human-readable result instead of a raw vector, a function `process_raw()` is created, which converts the raw vector to a character vector yielding a JSON string.

```{r rest-simple}
urls <- "https://httpbin.org/ip"

pretty_json <- function(x)
  jsonlite::prettify(rawToChar(x))

do_requests_serial(urls, finally = pretty_json)
```

Note that a list of length 1 is returned by `do_requests_serial()`. This is because the `do_requests_*()` functions are vectorized over requests and therefore return a list with as many entries as requests. A slightly more involved example for creating a request using `do_requests_serial()` is as follows.

```{r rest-json}
urls <- rep("https://httpbin.org/post", 2)
json <- list(list(a = "foo"),
             list(a = "bar"))

post_handle <- function(x)
  curl::handle_setheaders(
    curl::new_handle(postfields = charToRaw(jsonlite::toJSON(x))),
    "Content-Type" = "application/json"
  )

process_json <- function(x)
  jsonlite::fromJSON(rawToChar(x))$json

res <- do_requests_serial(urls, json,
                          create_handle = post_handle,
                          finally = process_json)
identical(res, json)
```

In order to customize the curl handle, a `create_handle` function is supplied, which receives for each request the corresponding entry of the object passed as `bodies` argument. The `finally` function in this example parses the returned JSON object into a list. Since for POST requests, `httpbin` mirrors the POST data, the initial `json` list is returned.

## Creating a JSON-RPC request

Building on functionality offered by the `do_requests_*()` functions, `make_requests()` helps with constructing JSON-RPC calls. As per the JSON-RPC 2.0 specification, a request object is a JSON string with the following members^[[**JSON-RPC 2.0 Specification** Rev. 2013-01-04. JSON-RPC Working Group](http://www.jsonrpc.org/specification)]:

* `jsonrpc`: A string specifying the version of the JSON-RPC protocol.
* `method`: A string containing the name of the method to be invoked.
* `params`: A structured value that holds the parameter values to be used during the invocation of the method.
* `id`: An identifier established by the client.

In order to create a list of such request objects, the arguments `methods`, `params`, `ids` and `version`, which map to `method`, `params`, `id` and `jsonrpc`, respectively, are combined and converted to JSON. OpenBIS API endpoint urls are constructed with the helper function `api_url()` which can receive arguments `api_endpoint` and `host_url`, forwarded by `make_requests()`. The argument `n_con` is used to specify the maximal number of simultaneous connections made to the server.

As `make_requests()` is designed to construct several requests at the time, the arguments `methods` and `params` are vectorized. If any of the passed objects is of length 1, it is replicated using `base::rep()` to the required length such that all three objects are of the same length. Care has to be taken with the list passed as `params` such that its length corresponds to the number of requests. This means that if only a single set of parameters is passed, this list has to be wrapped by another list.

As a last argument, a `finally` function can be passed, which defaults to `process_json()`. This function first converts all typed JSON objects to `json_class` objects using `as_json_class()`, resolves object references using `resolve_references()` in order to make all `json_class` objects self-contained and then creates a `json_vec` object using `as_json_vec()`, therefore allowing S3 method dispatch on lists of `json_class` objects. For more information on the specifics of this, please refer to the vignette ["JSON object handling"](json-class.html).

For single requests, a wrapper around `make_requests()` is available as `make_request()`. This function wraps the `params` argument in a list such that a list of length 1 is passed to `make_requests()` and returns the first entry of the list resulting from calling `make_requests()`. The following example shows how `make_request()` is used to implement the API method [`listProjects`](https://svnsis.ethz.ch/doc/openbis/13.04.0/ch/systemsx/cisd/openbis/generic/shared/api/v1/IGeneralInformationService.html#listProjects%28java.lang.String%29), which takes an access token passed as a list as its only argument.

```{r simple-rpc}
token <- login_openbis()

projects <- make_request(api_url(api_endpoint = "gis"),
                         method = "listProjects",
                         params = list(token))
print(projects, length = 10L)

logout_openbis(token)
```

To illustrate the used of `make_requests()`, the API method [`getDownloadUrlForFileForDataSet`](https://svnsis.ethz.ch/doc/openbis/13.04.0/ch/systemsx/cisd/openbis/dss/generic/shared/api/v1/IDssServiceRpcGeneric.html#getDownloadUrlForFileForDataSet%28java.lang.String,%20java.lang.String,%20java.lang.String%29) is implemented. This function generates a download URL for a file in a data set and requires a separate API call for each combination of data set and file path. In order to carry out several of these requests asynchronously, a list has to be passed as `params` argument such that each entry contains a list holding a login token, a dataset code and a file path.

```{r multi-rpc}
token <- login_openbis()

dataset_codes <- c("20120629093035782-603380", "20121011143734361-1359915")
file_path <- c("original/aThresholdedInfectionScoring_bBB01-1I.csv")

params <- lapply(dataset_codes, function(x) list(token, x, file_path))

str(params)

donwload_urls <- make_requests(api_url(api_endpoint = "dsrg"),
                               method = "getDownloadUrlForFileForDataSet",
                               params = params)
str(donwload_urls)

logout_openbis(token)
```

In addition to the JSON-RPC calls enabled by `make_request()`/`make_requests()`, the `do_requests_*()` functions are also used for (asynchronous) file downloads in `fetch_files()`. Slightly different functions are used as `create_handle` and `check` arguments but apart from that the logic remains largely the same.

Any arguments passed as `...` to `make_requests()` will be forwarded to `api_url()` which creates an API endpoint URL and passes this on to the `do_requests_*()` functions. As all functions that issue API calls use `make_requests()`, this makes it possible to not only target the InfectX openBIS instance, but arbitrary openBIS servers that support the v1 JSON-RPC API. The following example serves to illustrate how this mechanism can be used to access the openBIS demo.

```{r third-party, fig.align = "center"}
token <- login_openbis("test_observer", "test_observer",
                       host_url = "https://openbis-eln-lims.ethz.ch")

gel_data <- search_openbis(
  token,
  search_criteria(
    attribute_clause("type", "ELN_PREVIEW"),
    sub_criteria = search_sub_criteria(
      search_criteria(
        attribute_clause(value = "WB_LEXA-ER-B112")
      ),
      type = "sample"
    )
  )
)

gel_img <- fetch_files(token, gel_data,
                       reader = magick::image_read)

logout_openbis(token)

print(gel_img[[1]])
attributes(gel_img[[1]])
```

`api_url()` can accept arguments `api_endpoint`, `host_url` and `full_url` and ignores any further arguments. The first argument offers a selection of hard-coded API endpoints and in combination with the second argument is used to access one of them on the specified openBIS host. In case a URL that is not supported by the `api_endpoint` selector is desired, the complete URL can
be passed as `full_url` which will simply be returned by `api_url()`.

## Summary of API methods

Currently no methods from the API sections

* [IGeneralInformationChangingService](http://svnsis.ethz.ch/doc/openbis/13.04.0/ch/systemsx/cisd/openbis/generic/shared/api/v1/IGeneralInformationChangingService.html)
* [IQueryApiServer](http://svnsis.ethz.ch/doc/openbis/13.04.0/ch/systemsx/cisd/openbis/plugin/query/shared/api/v1/IQueryApiServer.html)
* [IWebInformationService](http://svnsis.ethz.ch/doc/openbis/13.04.0/ch/systemsx/cisd/openbis/generic/shared/api/v1/IWebInformationService.html)

are implemented as they mostly deal with modifying metadata and creating new aggregation reports. The main focus of this package is retrieving data and therefore only methods from the API sections

```{r api-summary, echo = FALSE, results = "asis", tidy = FALSE}
for (i in seq_along(info))
  cat(paste0("* ", info[[i]]$link, ", ", info[[i]]$status), "% implemented\n")
```

are currently available. A more detailed overview of what functionality is implemented in this API client is given in the following sections.

```{r api-tables, echo = FALSE, results = "asis", tidy = FALSE}
for (i in seq_along(info)) {
  cat(paste0("### ", info[[i]]$title, "\n", info[[i]]$description, "\n"))
  print(kable(info[[i]]$table))
  cat("\n")
}
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/json-class.R
\name{json_class}
\alias{json_class}
\alias{as_json_class}
\alias{as.json_class}
\alias{as_json_class.json_class}
\alias{as_json_class.list}
\alias{as_json_class.json_vec}
\alias{as_json_class.default}
\alias{rm_json_class}
\alias{as.list.json_class}
\alias{is_json_class}
\alias{is.json_class}
\alias{check_json_class}
\title{Create and validate JSON class objects}
\usage{
json_class(..., class)

as_json_class(x, ...)

as.json_class(x, ...)

\method{as_json_class}{json_class}(x, ...)

\method{as_json_class}{list}(x, recursive = TRUE, ...)

\method{as_json_class}{json_vec}(x, ...)

\method{as_json_class}{default}(x, force = FALSE, ...)

rm_json_class(x, recursive = TRUE, restore_type = TRUE)

\method{as.list}{json_class}(x, keep_asis = TRUE,
  recursive = !keep_asis, restore_type = !keep_asis, ...)

is_json_class(x)

is.json_class(x)

check_json_class(x, recursive = TRUE)
}
\arguments{
\item{...}{Generic compatibility.}

\item{class}{JSON sub-class name.}

\item{x}{Object to process.}

\item{recursive}{Recursively apply the function.}

\item{force}{Suppress error when casting an object to \code{json_class} that
cannot be converted.}

\item{restore_type}{When removing the \code{json_class} information from an
object, whether to preserve the subclass attribute as \code{@type} filed.}

\item{keep_asis}{Used in \code{as.list()}, if \code{TRUE}, the \code{json_class} object
is returned as-is, if \code{FALSE}, class attributes may be dropped/written to
the list structure into the \code{@type} field.}
}
\value{
JSON objects are represented as S3 classes with two class
attributes: \code{json_class} and an object specific class name (typically to
mirror the openBIS class).
}
\description{
To communicate object type information via JSON to the Jackson-powered
openBIS interface, the \code{@type} field is used. Data received from openBIS is
recursively stripped of \code{@type} fields and the type information is saved as
class attribute. Such objects also have the class \code{json_class} added. The
function \code{as_json_class()} (or its alias \code{as.json_class()}) powers this
recursive conversion of a list with filed \code{@type} into a \code{json_class}
object. The constructor \code{json_class()} is available for non-recursive
instantiation of \code{json_class} objects.
}
\details{
The action of \code{as_json_class()} is reversed by \code{rm_json_class()}. This
removes both the \code{json_class} class attribute and the JSON class attribute
itself, which is subsequently written to a \code{@type} filed. This preserving
of type information can be disabled, by setting the argument \code{restore_type}
to \code{FALSE}. Furthermore, the action can be applied recursively with the
argument \code{recursive}. The function \code{as.list()} can also be used to perform
the above actions, but with default arguments, it does nothing, as
functions such as \code{\link[base:sapply]{base::sapply()}} and \code{\link[base:lapply]{base::lapply()}}, call \code{as.list()}.

JSON class objects have custom sub-setting and printing functions available.
Sub-setting of JSON objects that preserve class and \code{json_class}
attributes. This is useful when objects are created from openBIS results
which are subsequently used in further queries, but the constructors they
are passed to require only a subset of the fetched fields.

The functions \code{is_json_class()} tests whether an object is a proper JSON
class object, meaning that:
\itemize{
\item it is a list
\item it inherits \code{json_class}
\item the last class attribute is \code{json_class}
\item apart from \code{json_class} there exists at least one more class attribute
}

In order to recursively test a \code{json_class} object for being properly
formed, the function \code{check_json_class()} can be used. This recurses through
a list structure and whenever an object inherits from \code{json_class} it is
tested with \code{is_json_class()}.
}
\examples{
lst <- list(`@type` = "foobar",
            a = list(`@type` = "foo", b = "c"),
            d = list(`@type` = "bar", e = "f"))

cls <- as_json_class(lst)
print(cls, depth = 2)

is_json_class(cls)
get_subclass(cls)

# recursive validation of json_class objects with check_json_class()
attr(cls[["d"]], "class") <- "json_class"
is_json_class(cls)
check_json_class(cls)

# as_json_class() is idempotent
identical(as_json_class(lst), as_json_class(as_json_class(lst)))

# rm_json_class() reverses the action of as_json_class()
identical(lst, rm_json_class(as_json_class(lst)))

# json_class objects can be instantiated using the constructor json_class()
identical(as_json_class(lst), 
          json_class(a = json_class(b = "c", class = "foo"),
                     d = json_class(e = "f", class = "bar"),
                     class = "foobar"))

cls <- as_json_class(lst)

# the default of as.list does nothing
identical(cls, as.list(cls))
# this can be disabled, by setting keep_asis to FALSE
identical(lst, as.list(cls, keep_asis = FALSE))
# further options are disabling recursive action
as.list(cls, keep_asis = FALSE, recursive = FALSE)
# and dropping type information
as.list(cls, keep_asis = FALSE, recursive = FALSE, restore_type = FALSE)

}
\seealso{
Other json object handling functions: \code{\link{has_fields.json_class}},
  \code{\link{json_vec}}, \code{\link{print.json_class}}
}
\concept{json object handling functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plate.R
\name{list_plates}
\alias{list_plates}
\alias{list_plates.NULL}
\alias{list_plates.ExperimentIdentifier}
\alias{list_plates.Experiment}
\alias{list_plate_metadata}
\alias{list_plate_metadata.PlateIdentifier}
\alias{list_plate_metadata.Plate}
\alias{list_plate_metadata.Sample}
\alias{plate_id}
\alias{as_plate_id}
\alias{as_plate_id.Plate}
\alias{as_plate_id.Sample}
\alias{as_plate_id.PlateMetadata}
\alias{as_plate_id.PlateIdentifier}
\alias{list_wells}
\alias{list_wells.PlateIdentifier}
\alias{list_wells.Plate}
\alias{list_wells.Sample}
\alias{list_wells.MaterialScreening}
\alias{list_wells.MaterialIdentifierScreening}
\alias{list_wells.MaterialGeneric}
\alias{list_wells.MaterialIdentifierGeneric}
\alias{well_id}
\alias{as_well_id}
\alias{as_well_id.WellMetadata}
\alias{as_well_id.WellIdentifier}
\alias{well_pos}
\title{List plates and wells}
\usage{
list_plates(token, x = NULL, ...)

\method{list_plates}{NULL}(token, x, ...)

\method{list_plates}{ExperimentIdentifier}(token, x, ...)

\method{list_plates}{Experiment}(token, x, ...)

list_plate_metadata(token, x, ...)

\method{list_plate_metadata}{PlateIdentifier}(token, x, ...)

\method{list_plate_metadata}{Plate}(token, x, ...)

\method{list_plate_metadata}{Sample}(token, x, ...)

plate_id(code, space)

as_plate_id(x, ...)

\method{as_plate_id}{Plate}(x, ...)

\method{as_plate_id}{Sample}(x, ...)

\method{as_plate_id}{PlateMetadata}(x, ...)

\method{as_plate_id}{PlateIdentifier}(x, ...)

list_wells(token, x, ...)

\method{list_wells}{PlateIdentifier}(token, x, ...)

\method{list_wells}{Plate}(token, x, ...)

\method{list_wells}{Sample}(token, x, ...)

\method{list_wells}{MaterialScreening}(token, x, experiment = NULL,
  include_datasets = FALSE, ...)

\method{list_wells}{MaterialIdentifierScreening}(token, x,
  experiment = NULL, include_datasets = FALSE, ...)

\method{list_wells}{MaterialGeneric}(token, x, experiment = NULL,
  include_datasets = FALSE, ...)

\method{list_wells}{MaterialIdentifierGeneric}(token, x,
  experiment = NULL, include_datasets = FALSE, ...)

well_id(perm_id, plate, well_pos = NULL, well_code = NULL, ...)

as_well_id(x, ...)

\method{as_well_id}{WellMetadata}(x, ...)

\method{as_well_id}{WellIdentifier}(x, ...)

well_pos(row = NULL, col = NULL, name = NULL)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to limit the number of returned wells or plates.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}

\item{code, space}{Character vectors that together can be used to create
\code{PlateIdentifier} objects.}

\item{experiment}{Additionally, the search can be limited to a single
experiment, specified either as \code{Experiment} or \code{ExperimentIdentifier}.}

\item{include_datasets}{Logical switch indicating whether to also return
the connected image and image analysis data sets.}

\item{perm_id, plate, well_pos}{Character vector, set of plate objects and
set of well position objects, all of the same length or length 1, that
together can be used to create \code{WellIdentifier} objects.}

\item{well_code}{Character vector where each entry is of the form
barcode:well_name, e.g. FOO-BAR-1:A1, FOO-BAR-1:A2, etc.}

\item{row, col}{Character vector of plate row names or numeric vector of
plate row indices and numeric vector of plate column indices, both of the
same length or of length 1.}

\item{name}{Character vector of well name, where each entry is of the form
A1, A2, etc.}
}
\value{
Depending on the number of resulting objects, either a
\code{\link{json_class}} (single object) or a \code{\link{json_vec}} (multiple objects), is
returned. For \code{list_plates()}, the additional class attribute \code{Plate} is
added, while \code{list_plate_metadata()} returns \code{PlateMetadata} and
\code{plate_id()}/\code{as_plate_id()} yield \code{PlateIdentifier} objects. For
individual wells, \code{list_wells()}, \code{well_id()} and \code{as_well_id()} all return
objects with sub-type \code{WellIdentifier}. Finally, \code{WellPosition} objects are
created by \code{well_pos()}.
}
\description{
Plates and wells are special cases of \code{Sample} objects and play an
important organizational role when openBIS is used in HTS experiments. All
InfectX screens were arrayed onto 384-well plates, arranged into 16 rows
(A through P) and 24 columns (1 through 24). \code{Plate} and \code{PlateIdentifier}
objects are used to identify plates while for wells only \code{WellIdentifier}
objects exists for representing individual wells. Additional objects
relevant in this context are \code{PlateMetadata}, which for all associated wells
contain respective \code{WellMetadata} objects and
\code{PlateWellReferenceWithDatasets} objects, each holding a \code{Plate} object and
a \code{WellPosition} object, thereby in a sense also identifying individual
wells.
}
\details{
\code{Plate} objects are listed using \code{list_plates()}, which can either list all
available plates (default or dispatch on \code{NULL}) or restrict the listing to
a set of supplied experiment objects (dispatch on either \code{Experiment} or
\code{ExperimentIdentifier} objects). If multiple experiments are used for
limiting the search (i.e. a \code{json_vec} of experiments), a separate API call
for each object has to be made. \code{PlateMetadata} objects (including all
corresponding \code{WellMetadata} objects) are retrieved by
\code{list_plate_metadata()} which can be dispatched on objects that represent
plates, including \code{Plate}, \code{PlateIdentifier} and \code{Sample} (given that the
sample is of type \code{PLATE}). Finally, \code{PlateIdentifier} can be created
either by calling \code{plate_id()} or though coercion of \code{Plate}, \code{Sample}
(with type \code{PLATE}) or \code{PlateMetadata}  objects using the function
\code{as_plate_id()}. Neither \code{plate_id()} nor \code{as_plate_id()} incur API calls.

Wells can be listed with \code{list_wells()}, which returns \code{WellIdentifier}
objects if dispatch occurs on objects representing plates, including
\code{Plate}, \code{PlateIdentifier} and \code{Sample} (with type \code{PLATE}). In this case
the entire set of well id objects corresponding to the selected plates is
returned and a separate API call is required per plate.

Whenever \code{list_wells()} is dispatched on material objects (any of
\code{MaterialScreening}, \code{MaterialIdentifierScreening}, \code{MaterialGeneric} or
\code{MaterialIdentifierGeneric}), \code{PlateWellReferenceWithDatasets} objects are
returned, representing wells associated with the given material. If multiple
material ids are passed, an API call for each object is issued. The well
search can be limited to an experiment by passing a single \code{Experiment} or
\code{ExperimentIdentifier} object as \code{experiment} argument and image dataset
references as well as feature vector dataset references can be retrieved
as part of the \code{PlateWellReferenceWithDatasets} objects if the logical
switch \code{include_datasets} is set to \code{TRUE}. A separate API call per passed
material object is required.

Instantiation of \code{WellIdentifier} objects can be done either using the
constructor \code{well_id()} or via coercion of \code{WellMetadata} objects by
calling \code{as_well_id()}. Well samples cannot be coerced to well id objects
as they do not contain all fields that are required. A further object type
relevant to this context is that of \code{WellPosition}, encoding the position
of a well within a plate. Such objects can be created using the constructor
\code{well_pos()}.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("sas", "listPlates")}
\item \Sexpr[results=rd]{infx::docs_link("sas", "listPlateWells")}
\item \Sexpr[results=rd]{infx::docs_link("sas", "getPlateMetadataList")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  # search for an experiment, e.g. ADENO-AU-K1
  exp <- search_openbis(tok,
                        search_criteria(
                          property_clause("pathogen", "Adenovirus"),
                          property_clause("library", "Ambion"),
                          property_clause("geneset", "Kinome"),
                          property_clause("replicate", 1L)
                        ),
                        target_object = "experiment")

  # list all plates associated with this experiment
  plates <- list_plates(tok, exp)
  length(plates)
  as_plate_id(plates)

  # for a plate, fetch meta data objects
  meta <- list_plate_metadata(tok, plates[[1L]])
  print(meta, depth = 2L, length = 15L)
  print(meta[["wells"]][[1L]], depth = 2L)

  # search for a sample object corresponding to plate BB01-1I
  samp <- search_openbis(tok,
                         search_criteria(
                           attribute_clause("code",
                                            "/INFECTX_PUBLISHED/BB01-1I")
                         ),
                         target_object = "sample")

  # list all wells for this sample
  wells <- list_wells(tok, samp)
  identical(as_well_id(meta[["wells"]][[1L]]),
            wells[[1L]])

  # search for the material corresponding to MTOR
  mat <- search_openbis(tok,
                        search_criteria(
                          property_clause("gene_symbol", "MTOR")
                        ),
                        target_object = "material")
  # search for associated wells, limited to ADENO-AU-K1
  mtor <- list_wells(tok, mat, exp)
  plates <- get_field(mtor, "experimentPlateIdentifier")
  as_plate_id(plates)
  unique(get_field(mtor, "wellPosition"))

  logout_openbis(tok)
}

}
\seealso{
Other object listing functions: \code{\link{list_datasets}},
  \code{\link{list_experiments}},
  \code{\link{list_material}}, \code{\link{list_projects}},
  \code{\link{list_samples}}
}
\concept{object listing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/json-class.R, R/json-utils.R, R/json-vec.R
\name{has_fields.json_class}
\alias{has_fields.json_class}
\alias{get_field.json_class}
\alias{has_subclass.json_class}
\alias{get_subclass.json_class}
\alias{has_fields}
\alias{has_fields.default}
\alias{get_field}
\alias{has_subclass}
\alias{has_subclass.default}
\alias{get_subclass}
\alias{get_subclass.list}
\alias{remove_null}
\alias{has_fields.json_vec}
\alias{get_field.json_vec}
\alias{has_subclass.json_vec}
\alias{get_subclass.json_vec}
\title{JSON class utilities}
\usage{
\method{has_fields}{json_class}(x, fields, ...)

\method{get_field}{json_class}(x, field, ...)

\method{has_subclass}{json_class}(x, class, ...)

\method{get_subclass}{json_class}(x)

has_fields(x, fields, ...)

\method{has_fields}{default}(x, ...)

get_field(x, field, ...)

has_subclass(x, class, ...)

\method{has_subclass}{default}(x, ...)

get_subclass(x)

\method{get_subclass}{list}(x, ...)

remove_null(x)

\method{has_fields}{json_vec}(x, fields, ...)

\method{get_field}{json_vec}(x, field, ...)

\method{has_subclass}{json_vec}(x, class, ...)

\method{get_subclass}{json_vec}(x, ...)
}
\arguments{
\item{x}{Object to test.}

\item{fields}{Character vector of nonzero length, holding the field names
for which to check.}

\item{...}{Generic compatibility.}

\item{field}{Character vector of length 1, holding the field name to
extract.}

\item{class}{Character vector of nonzero length, holding the class names
to test for.}
}
\value{
Depending on whether a single or a set of multiple objects is
represented, the S3 classes \code{\link{json_class}} or \code{\link{json_vec}} are applied
respectively.
}
\description{
Several utility functions for working with \code{json_class} and \code{json_vec}
objects are provided. This includes \code{has_fields()} for checking whether
certain fields are available in an object, \code{get_field()} to extract
values from an object that correspond to a field with a certain name,
\code{has_subclass()} for testing that an object is of a certain class and
\code{get_subclass()} to extract this class. Finally, NULL fields can be
recursively removed using \code{remove_null()}. More information is available
in the details section.
}
\details{
The generic function \code{has_fields()} tests whether a single \code{json_class}
object contains all of the specified fields or whether each \code{json_class}
object contained in a \code{json_vec} object passes this test. If dispatch
occurs on an object that is neither of class \code{json_class}, nor of class
\code{json_vec}, \code{has_fields()} returns \code{FALSE}. A single field can be extracted
from a \code{json_class} or a \code{json_vec} object, using \code{get_field()}. Iteration
for \code{json_vec} objects happens via \code{\link[base:sapply]{base::sapply()}} so that when possible
the result is simplified.

In order to test whether a \code{json_class} or a \code{json_vec}  object is of a
certain sub-class (can also be a vector of sub-classes), the generic
function \code{has_subclass()} can be used. Dispatch on objects that do not
inherit from either \code{json_class} or \code{json_vec} will return \code{FALSE}. The
sub-class of a \code{json_class} or a \code{json_vec} object can be determined, using
\code{get_subclass}. This will also work if dispatched on a \code{list} of objects if
that list object passes \code{\link[=has_common_subclass]{has_common_subclass()}}.

The function \code{remove_null()} recursively removes all NULL fields from a
nested list structure while preserving \code{json_class} and \code{json_vec} class
attributes. This can be useful when fetching an object form openBIS and
subsequently using this object for a further query: whenever the object
returned by the first API call contains NULL fields, it is safer to remove
all of them, as in some cases this might cause an error in the following
API requests.
}
\examples{
obj_1 <- json_class(a = 1, b = 2, class = "foo")
obj_2 <- json_class(a = 2, b = 4, class = "foo")
obj_3 <- json_class(a = 3, c = 6, class = "foo")

# one or more fields can be tested
has_fields(obj_1, "a")
has_fields(obj_1, c("a", "b"))
# dispatch on json_vec objects is possible as well
has_fields(c(obj_1, obj_2), "a")
has_fields(c(obj_1, obj_2), "b")
has_fields(c(obj_1, obj_3), "b")
has_fields(c(obj_1, obj_3), c("a", "b"))
# other types do not pass the test
has_fields(list(obj_1, obj_3), "a")

get_field(obj_1, "a")
get_field(c(obj_1, obj_3), "a")
\dontrun{
  # the requested field must be available in every instance
  get_field(c(obj_1, obj_3), "b")
  # only a single field may be requested
  get_field(c(obj_1, obj_2), c("a", "b"))
}

obj_4 <- json_class(a = 4, c = 8, class = "bar")

# dispatch on json_class
has_subclass(obj_1, "foo")
# dispatch on json_vec
has_subclass(c(obj_1, obj_2), "foo")
# dispatch on other object types always returns FALSE
has_subclass(list(obj_1, obj_2), "foo")

# dispatch on json_class
get_subclass(obj_1)
# dispatch on json_vec
get_subclass(c(obj_1, obj_2))
# dispatch on list is possible if the list passes has_common_subclass()
get_subclass(list(obj_1, obj_2))
\dontrun{
  get_subclass(list(obj_1, obj_4))
}

tmp <- json_class(a = json_class(b = "c", d = NULL, class = "foo"),
                  e = json_class(f = "g", class = "bar"),
                  h = NULL,
                  class = "foobar")
print(tmp, 2)
print(remove_null(tmp), 2)

}
\seealso{
Other json object handling functions: \code{\link{json_class}},
  \code{\link{json_vec}}, \code{\link{print.json_class}}
}
\concept{json object handling functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.R
\name{list_datasets}
\alias{list_datasets}
\alias{list_datasets.Sample}
\alias{list_datasets.Experiment}
\alias{list_datasets.character}
\alias{list_dataset_ids}
\alias{list_dataset_ids.character}
\alias{list_dataset_ids.DataSet}
\alias{list_references}
\alias{list_references.PlateIdentifier}
\alias{list_references.Plate}
\alias{list_references.PlateMetadata}
\alias{list_references.Sample}
\alias{list_references.MaterialGeneric}
\alias{list_references.MaterialScreening}
\alias{list_references.MaterialIdentifierGeneric}
\alias{list_references.MaterialIdentifierScreening}
\alias{list_references.DatasetIdentifier}
\alias{list_references.DataSet}
\alias{list_references.DatasetReference}
\alias{list_references.FeatureVectorDatasetReference}
\alias{list_references.FeatureVectorDatasetWellReference}
\alias{list_references.ImageDatasetReference}
\alias{list_references.MicroscopyImageReference}
\alias{list_references.PlateImageReference}
\alias{list_dataset_types}
\title{List datasets and dataset reference objects}
\usage{
list_datasets(token, x, ...)

\method{list_datasets}{Sample}(token, x, include = c(NA, "children",
  "parents", "all"), ...)

\method{list_datasets}{Experiment}(token, x, include = c(NA, "children",
  "parents", "all"), ...)

\method{list_datasets}{character}(token, x, include = c(NA, "children",
  "parents", "all"), ...)

list_dataset_ids(token, x, ...)

\method{list_dataset_ids}{character}(token, x, ...)

\method{list_dataset_ids}{DataSet}(token, x, ...)

list_references(token, x, ...)

\method{list_references}{PlateIdentifier}(token, x, type = c("raw",
  "segmentation", "feature"), ...)

\method{list_references}{Plate}(token, x, type = c("raw", "segmentation",
  "feature"), ...)

\method{list_references}{PlateMetadata}(token, x, type = c("raw",
  "segmentation", "feature"), ...)

\method{list_references}{Sample}(token, x, type = c("raw",
  "segmentation", "feature"), ...)

\method{list_references}{MaterialGeneric}(token, x, experiment = NULL,
  ...)

\method{list_references}{MaterialScreening}(token, x, experiment = NULL,
  ...)

\method{list_references}{MaterialIdentifierGeneric}(token, x,
  experiment = NULL, ...)

\method{list_references}{MaterialIdentifierScreening}(token, x,
  experiment = NULL, ...)

\method{list_references}{DatasetIdentifier}(token, x, wells = NULL,
  channels, ...)

\method{list_references}{DataSet}(token, x, wells = NULL, channels, ...)

\method{list_references}{DatasetReference}(token, x, wells = NULL,
  channels, ...)

\method{list_references}{FeatureVectorDatasetReference}(token, x,
  wells = NULL, channels, ...)

\method{list_references}{FeatureVectorDatasetWellReference}(token, x,
  wells = NULL, channels, ...)

\method{list_references}{ImageDatasetReference}(token, x, wells = NULL,
  channels, ...)

\method{list_references}{MicroscopyImageReference}(token, x,
  wells = NULL, channels, ...)

\method{list_references}{PlateImageReference}(token, x, wells = NULL,
  channels, ...)

list_dataset_types(token, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to limit search for datasets/files with.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}

\item{include}{String indicating whether to include parent/child datasets
as well.}

\item{type}{For listing image datasets, it can be specified, whether only
raw image datasets, only segmentation image datasets or any kind of image
datasets (default) are to be listed.}

\item{experiment}{When searching for datasets associated with materials,
the search can be limited to a single experiment.}

\item{wells}{A (set of) \code{WellPosition} object(s) to limit the dataset
listing to.}

\item{channels}{A character vector with imaging channel names to limit the
dataset listing to.}
}
\value{
Depending on the number of resulting objects, either a
\code{\link{json_class}} (single object) or a \code{\link{json_vec}} (multiple objects), is
returned. For the specific sub-class, refer to the \emph{Details} section.
}
\description{
All available datasets for the specified experiment(s), sample(s) or
dataset code(s) are retrieved as \code{DataSet} objects by
\code{list_datasets()}. Each dataset has a type and all realized type
identifiers can be listed using \code{list_dataset_types()}. A more compact
object type, uniquely identifying a \code{DataSet}, is that of a
\code{DatasetIdentifier}. Given either a set of \code{DataSet} objects or a character
vector holding (a) dataset code(s), \code{list_dataset_id()} fetches the
corresponding \code{DatasetIdentifier} objects. Behavior of the function
\code{list_references()}, in particular the returned object type, depends on
input types. For more information, please refer to the details section.
}
\details{
\code{list_datasets()} is an s3 generic function that can be dispatched on
\code{Sample} and \code{Experiment} objects, as well as character vectors containing
dataset codes and it returns sets of \code{DataSet} objects. Additionally it
can be requested that parent or child datasets are to be included as well.

Several classes in addition to \code{DatasetIdentifier} implement the
\code{IDatasetIdentifier} interface, including
\itemize{
\item \code{DatasetReference}
\item \code{FeatureVectorDatasetReference}
\item \code{FeatureVectorDatasetWellReference}
\item \code{ImageDatasetReference}
\item \code{MicroscopyImageReference}
\item \code{PlateImageReference}
}

The return type of \code{list_references()} depends on dispatch object type and
in some cases on additional arguments. If the s3 generic function
\code{list_references()} is dispatched on plate objects (\code{Plate},
\code{PlateIdentifier} or \code{PlateMetadata} or \code{Sample} objects, representing
plates), \code{ImageDatasetReference} objects are returned (except if the type
argument is set to \code{feature}, in which case, if
\code{MaterialIdentifierScreening} objects are used as input,
\code{PlateWellReferenceWithDatasets} objects are returned, which each contain
\code{ImageDatasetReference} and \code{FeatureVectorDatasetReference} objects.

Whenever \code{list_references()} is dispatched on dataset ids or dataset
reference objects, the resulting object type depends on whether a (set of)
\code{WellPosition} object(s) were specified as \code{wells} argument. For its
default value (NULL), a set of \code{MicroscopyImageReference} objects is
returned, while \code{PlateImageReference} objects are returned otherwise.
}
\section{Implementation note}{
 The API function \code{listDataSetsForSample()}
has a parameter \code{areOnlyDirectlyConnectedIncluded}, which is currently
fixed to \code{TRUE}. The  documentation contains the following explanation:

If true, only datasets that are directly connected to the sample are
included, otherwise datasets of child samples are included as well.

This does however not seem to correspond to including child datasets in the
API call to \code{listDataSets()} via its \code{connectionsToGet} argument. As long
as it is not entirely clear how the inclusion of child/parent datasets
differs from setting \code{areOnlyDirectlyConnectedIncluded} to \code{FALSE}, this
option is not exposed to the user.
}

\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listDataSetsForSample")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "listDataSets")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listDataSetsForExperiments")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "getDataSetMetaData")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("sas", "getDatasetIdentifiers")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("sas", "listRawImageDatasets")}
\item \Sexpr[results=rd]{infx::docs_link("sas",
"listSegmentationImageDatasets")}
\item \Sexpr[results=rd]{infx::docs_link("sas", "listFeatureVectorDatasets")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "listImageReferences")}
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "listPlateImageReferences")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listDataSetTypes")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  # search for a sample object corresponding to plate KB2-03-1I
  samp <- search_openbis(tok,
                         search_criteria(
                           attribute_clause("code",
                                            "/INFECTX_PUBLISHED/KB2-03-1I")
                         ),
                         target_object = "sample")

  # list all datasets associated with this plate
  ds <- list_datasets(tok, samp)

  # select a feature dataset, note how the fields "parentCodes" and
  # "childrenCodes" both are not set
  feat_ds <- ds[[grep("FEATURES_CC_MAT",
                      get_field(ds, "dataSetTypeCode"))]]

  # fetch parent and child datasets and now both the "parentCodes" and
  # "childrenCodes" fields are populated with the corresponding codes
  feat_ds <- list_datasets(tok, get_field(feat_ds, "code"), include = "all")

  # re-using the plate sample from above, an ImageDatasetReference object
  # corresponding to the associated raw imaging dataset is listed
  raw_ref <- list_references(tok, samp)
  # available imaging channels are
  get_field(raw_ref, "properties")[["IMAGE.CHANNEL.LIST"]]

  # a more specific image reference object can be retrieved by passing a
  # well specification to list_references()
  well_ref <- list_references(tok, raw_ref,
                              wells = well_pos(name = "A2"),
                              channel = "DAPI")
  # a reference to 9 images is returned, as there are 3 x 3 imaging tiles
  # per well
  length(well_ref)

  logout_openbis(tok)
}

}
\seealso{
Other object listing functions: \code{\link{list_experiments}},
  \code{\link{list_material}}, \code{\link{list_plates}},
  \code{\link{list_projects}}, \code{\link{list_samples}}
}
\concept{object listing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiment.R
\name{exp_id_str}
\alias{exp_id_str}
\alias{exp_id_str.ExperimentIdentifier}
\alias{exp_id_str.Experiment}
\title{Extract experiment string}
\usage{
exp_id_str(x, ...)

\method{exp_id_str}{ExperimentIdentifier}(x, ...)

\method{exp_id_str}{Experiment}(x, ...)
}
\arguments{
\item{x}{Experiment object(s).}

\item{...}{Generic compatibility.}
}
\description{
Experiments can be uniquely identified by a string made up of space code,
project code and experiment code. This function extracts the relevant fields
from experiment or experiment id objects and returns a vector of experiment
strings.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/json-utils.R
\name{print_json_class}
\alias{print_json_class}
\alias{style}
\title{Helper function for printing JSON objects}
\usage{
print_json_class(x, unnamed_parent = FALSE, cur_depth, max_depth,
  layout = style())

style(fancy = TRUE)
}
\arguments{
\item{x}{The JSON object to print.}

\item{unnamed_parent}{Whether the parent node is named or not (in some
cases, a different box character has to be used if this is true).}

\item{cur_depth}{The current recursion depth.}

\item{max_depth}{The maximum recursion depth.}

\item{layout}{Characters for printing the tree structure and styles to be
applied to the different entities.}

\item{fancy}{Logical switch to enable font styles, colors and UTF box
characters.}
}
\description{
This function powers the \code{json_class} and \code{json_vec} specific methods of the
base generic \code{\link[base:print]{base::print()}}. As it is applied recursively and recursion
depth has to be controllable, the function is aware of both the current
recursion depth (via \code{cur_depth}) and the maximally allowed recursion depth
(via \code{max_depth}). Furthermore the printing style (colored output and UTF
box characters for visualizing the tree structure) can be controlled through
the \code{layout} argument. Under some circumstances, this requires a given node
to know whether the parent node is a named object or not, which is passed
from a parent node to its children through the \code{unnamed_parent} argument.

In order to enable fancy printing (colored output and UTF box characters
for visualizing the tree structure), this function provides the required
styling information. Fancy printing can be disabled by setting the \code{fancy}
argument to \code{FALSE}, which yields ASCII characters for the tree structure
and disables color. This was more or less directly copied from Hadley's
\href{https://git.io/vFMA5}{lobstr} package.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image.R
\name{fetch_images}
\alias{fetch_images}
\alias{fetch_images.DatasetIdentifier}
\alias{fetch_images.DatasetReference}
\alias{fetch_images.ImageDatasetReference}
\alias{fetch_images.MicroscopyImageReference}
\alias{fetch_images.PlateImageReference}
\alias{list_image_metadata}
\alias{list_image_metadata.ExperimentIdentifier}
\alias{list_image_metadata.Experiment}
\alias{list_image_metadata.DatasetIdentifier}
\alias{list_image_metadata.DatasetReference}
\alias{list_image_metadata.ImageDatasetReference}
\alias{list_image_metadata.MicroscopyImageReference}
\alias{list_image_metadata.PlateImageReference}
\title{List image meta data and download images}
\usage{
fetch_images(token, x, ...)

\method{fetch_images}{DatasetIdentifier}(token, x, channels,
  well_positions = NULL, image_size = NULL, thumbnails = FALSE, ...)

\method{fetch_images}{DatasetReference}(token, x, channels,
  well_positions = NULL, image_size = NULL, thumbnails = FALSE, ...)

\method{fetch_images}{ImageDatasetReference}(token, x, channels,
  well_positions = NULL, image_size = NULL, thumbnails = FALSE, ...)

\method{fetch_images}{MicroscopyImageReference}(token, x,
  well_positions = NULL, image_size = NULL, thumbnails = FALSE, ...)

\method{fetch_images}{PlateImageReference}(token, x, image_size = NULL,
  force_png = FALSE, format = NULL, thumbnails = FALSE, ...)

list_image_metadata(token, x, ...)

\method{list_image_metadata}{ExperimentIdentifier}(token, x, ...)

\method{list_image_metadata}{Experiment}(token, x, ...)

\method{list_image_metadata}{DatasetIdentifier}(token, x,
  type = c("metadata", "format"), ...)

\method{list_image_metadata}{DatasetReference}(token, x,
  type = c("metadata", "format"), ...)

\method{list_image_metadata}{ImageDatasetReference}(token, x,
  type = c("metadata", "format"), ...)

\method{list_image_metadata}{MicroscopyImageReference}(token, x,
  type = c("metadata", "format"), ...)

\method{list_image_metadata}{PlateImageReference}(token, x,
  type = c("metadata", "format"), ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to limit the number of returned images}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}

\item{channels}{A character vector of imaging channels}

\item{well_positions}{A (list of) \code{WellPosition} objects. If the object
passed as argument x already contains well position information this can
be NULL.}

\item{image_size}{Either a single \code{ImageSize} object or NULL, in which case
images are returned in full size.}

\item{thumbnails}{Logical switch; if TRUE, thumbnail images are retrieved
in which case the arguments \code{well_positions} and \code{image_size} are expected
to be at their default values.}

\item{force_png}{Logical switch for making sure the returned image is a
png. If NULL or FALSE, the image is returned in the format it is stored.}

\item{format}{If not NULL, a single \code{ImageRepresentationFormat} object.
Cannot be combined with non-default \code{image_size}, \code{force_png} and \code{format}
arguments.}

\item{type}{Switch to specify the type of meta data objects to be returned.}
}
\value{
For \code{list_image_metadata()}, either a \code{\link{json_class}} (single
object) or a \code{\link{json_vec}} (multiple objects), is returned. For the specific
sub-class, refer to the \emph{Details} section. Image data retrieved with
\code{fetch_images()} is read by \code{\link[magick:image_read]{magick::image_read()}} and returned as
(possibly nested) \code{list} of \code{magick-image} objects.
}
\description{
Experiment level image meta data can be listed by passing experiment
representing objects (\code{Experiment} or \code{ExperimentIdentifier}) to
\code{list_image_metadata()} and data set level image meta data can be retrieved
by passing data set identifying objects which can be associated with image
data sets (data set id and data set reference objects). Images themselves
can be retrieved using \code{fetch_images()}. As with meta data listing, this
function can be dispatched on objects referencing or identifying data sets
associated with image data.
}
\details{
Data set level image meta data can be listed by passing objects, which
implement the \code{IDatasetIdentifier} interface and are connected to image
data sets (this rules out feature data set references and leaves
\code{DatasetIdentifier}, \code{DatasetReference}, \code{ImageDatasetReference},
\code{MicroscopyImageReference} and \code{PlateImageReference} objects). Two
different types of meta data objects are returned, depending on the \code{type}
argument: if it is set to \code{metadata} (default), objects of type
\code{ImageDatasetMetadata} and if it is set to \code{format}, objects of type
\code{DatasetImageRepresentationFormats} are returned. For experiment-level
image meta data, \code{ExperimentImageMetadata} objects are returned.

Dispatch of \code{fetch_images()} is available for the same object types as data
set-level image meta data listing: \code{DatasetIdentifier}, \code{DatasetReference},
\code{ImageDatasetReference}, \code{MicroscopyImageReference} and
\code{PlateImageReference}. The highest level of control over which images are
retrieved is achieved with \code{PlateImageReference} objects, which specify an
image data set, a well, a tile and a channel. The returned image format can
be modified by either passing an \code{ImageRepresentationFormat} object as the
\code{format} argument, by passing a single/list of format selection criterion
objects, which will be used to filter the available image representation
format objects or by specifying one or both of the \code{image_size} (expects an
\code{ImageSize} object) and \code{force_png} (logical switch) arguments.

\code{MicroscopyImageReference} objects contain channel information (as well as
tile information, which is not taken into account though). Therefore a
(list of) \code{WellPosition} object(s) has/have to be specified, for which then
all tiles are fetched for the given imaging channel. If the passed list of
\code{MicroscopyImageReference} objects contain instances that only differ in
tile number, redundancies are filtered out. An API call is necessary for
each non-redundant object.

Finally, \code{DatasetIdentifier}, \code{DatasetReference} and \code{ImageDatasetReference}
objects are all handled identically. For each of the specified data sets,
an imaging channel has to be provided and whenever the data set is
associated with an entire plate, a (list of) \code{WellPosition} object(s) as
well. If the data set is associated with a single well, the
\code{well_positions} can be left at its default value (NULL). If several data
sets are passed, an API call is necessary per data set. Possible
redundancies are not filtered.

Images are retrieved as Base64 encoded strings, which are converted to
binary using \code{\link[base64enc:base64decode]{base64enc::base64decode()}} and subsequently read by
\code{\link[magick:image_read]{magick::image_read()}}. Attached to the images is the information, based
on which they were retrieved, including data set object, well positions
(where applicable) and channel (where applicable). This results in a list
with length corresponding to the number of API calls that were necessary.
}
\section{Implementation notes}{

\itemize{
\item For dispatch on \code{PlateImageReference} objects, currently the only options
controlling the returned images are an argument for image size and a flag
for forcing the returned format to png. OpenBIS also supports
pre-defined image transformations to be applied to the images before they
are sent to the requesting party. These transformations can be requested
by a code (options are listed in \code{ImageRepresentationFormat} objects or
in \code{ImageChannel} objects attached to \code{ImageDatasetMetadata} objects).
However, as no such transformations appear to be defined, this is
currently not implemented.
\item When filtering \code{ImageRepresentationFormat} objects associated with a data
set, only \code{SizeCriterion} objects can be used. The remaining criteria
(\code{ColorDepthCriterion}, \code{FileTypeCriterion} and \code{OriginalCriterion}) are
currently disabled as they extend the abstract class
\code{AbstractFormatSelectionCriterion}, which causes an issue with JSON
deserialization.
}
}

\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "loadImagesBase64")}
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "loadThumbnailImagesBase64")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrs",
"loadPhysicalThumbnailsBase64")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("sas", "getExperimentImageMetadata")}
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "listImageMetadata")}
\item \Sexpr[results=rd]{infx::docs_link("dsrs",
"listAvailableImageRepresentationFormats")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  # search for a sample object corresponding to plate KB2-03-1I
  samp <- search_openbis(tok,
                         search_criteria(
                           attribute_clause("code",
                                            "/INFECTX_PUBLISHED/KB2-03-1I")
                         ),
                         target_object = "sample")
  # for the plate sample object, list raw image data set references
  ds_ref <- list_references(tok, samp)

  # the returned image dataset reference can be used to list image meta data
  img_meta <- list_image_metadata(tok, ds_ref)
  channels <- img_meta[["channelCodes"]]

  imgs <- fetch_images(tok, ds_ref,
                       channels = channels[[1L]],
                       well_positions = well_pos(1, 1),
                       image_size = json_class(width = 300, height = 300,
                                               class = "ImageSize"))
  # this yields 9 images, one per tile
  length(imgs[[1L]]) == img_meta[["numberOfTiles"]]
  # and each image is scaled to fit within 300 x 300 pixels
  magick::image_info(imgs[[1L]][[1L]])

  # if not the entire well is of interest, but only certain tiles
  img_ref <- list_references(tok, ds_ref,
                             wells = well_pos(1, 1),
                             channels = channels[[1L]])
  # this yields 9 objects, one reference per tile
  length(img_ref)
  # select a tile, for example the center one
  img <- fetch_images(tok, img_ref[[5L]],
                      image_size = json_class(width = 300, height = 300,
                                              class = "ImageSize"))
  identical(as.raster(img[[1L]]), as.raster(imgs[[1L]][[5L]]))

  logout_openbis(tok)
}

}
\seealso{
Other resource listing/downloading functions: \code{\link{list_download_urls}},
  \code{\link{list_features}}, \code{\link{list_files}}
}
\concept{resource listing/downloading functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project.R
\name{list_projects}
\alias{list_projects}
\title{List projects}
\usage{
list_projects(token, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{...}{Further arguments will be passed to \code{\link[=make_requests]{make_requests()}}.}
}
\value{
Depending on the number of resulting objects, either a
\code{\link{json_class}} (single object) or a \code{\link{json_vec}} (multiple objects), of
sub-type \code{Project} is returned.
}
\description{
A project forms one of the most basic entities in the organizational
hierarchy of openBIS. Each project consists of one or several experiments
and one or more projects are contained in each space, which is the topmost
structure used for grouping experiments. For InfectX, each project
corresponds to a separate pathogen. All registered projects visible to
the current user can be listed by calling the function \code{list_projects()}.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listProjects")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  proj <- list_projects(tok)
  length(proj)
  get_field(proj, "code")

  logout_openbis(tok)
}

}
\seealso{
Other object listing functions: \code{\link{list_datasets}},
  \code{\link{list_experiments}},
  \code{\link{list_material}}, \code{\link{list_plates}},
  \code{\link{list_samples}}
}
\concept{object listing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.R
\name{list_features}
\alias{list_features}
\alias{list_features.FeatureVectorDatasetReference}
\alias{list_features.FeatureVectorDatasetWellReference}
\alias{list_feature_codes}
\alias{list_feature_codes.FeatureVectorDatasetReference}
\alias{list_feature_codes.FeatureVectorDatasetWellReference}
\alias{fetch_features}
\alias{fetch_features.FeatureVectorDatasetReference}
\alias{fetch_features.FeatureVectorDatasetWellReference}
\title{List and download feature data}
\usage{
list_features(token, x, ...)

\method{list_features}{FeatureVectorDatasetReference}(token, x,
  wells = NULL, ...)

\method{list_features}{FeatureVectorDatasetWellReference}(token, x, ...)

list_feature_codes(token, x, ...)

\method{list_feature_codes}{FeatureVectorDatasetReference}(token, x,
  wells = NULL, ...)

\method{list_feature_codes}{FeatureVectorDatasetWellReference}(token, x,
  ...)

fetch_features(token, x, feature_codes = NA, ...)

\method{fetch_features}{FeatureVectorDatasetReference}(token, x,
  feature_codes = NA, wells = NULL, ...)

\method{fetch_features}{FeatureVectorDatasetWellReference}(token, x,
  feature_codes = NA, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to specify the set of feature vector datasets of interest.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_request]{make_request()}}.}

\item{wells}{Set of \code{WellPosition} objects used to limit the returned
feature data.}

\item{feature_codes}{A character vector of feature codes or NA (all
available feature codes).}
}
\value{
\code{list_feature_codes()} returns a character vector of feature codes,
while \code{list_features()} and \code{fetch_features()} return either \code{\link{json_class}}
(single object) or a \code{\link{json_vec}} (multiple objects), dependent on the
number of resulting objects. The sub-type is \code{FeatureInformation} for
\code{list_features()} and \code{FeatureVectorDataset} for \code{fetch_features()}.
}
\description{
In openBIS, features are datasets that are treated differently from
generic datasets. Briefly put, tabular datasets where columns correspond
to features and rows to wells can be marked as feature datasets which makes
it possible to query openBIS for individual feature values for selected
wells instead of having to download the entire table for a plate. The
relevant object types for handling feature data are
\code{FeatureVectorDatasetReference} and \code{FeatureVectorDatasetWellReference},
where the former represents feature data at plate level and the latter at
well level. The function \code{list_features()} can be used to enumerate
available features and \code{fetch_features()} will download feature data.
}
\details{
Listing of features can be performed by calling \code{list_features()} on
\code{FeatureVectorDatasetReference} or \code{FeatureVectorDatasetWellReference}
objects. Plate-level references of feature datasets can for example be
retrieved using \code{\link[=list_references]{list_references()}} and well-level references are created
whenever a \code{wells} argument is supplied to \code{list_features()}, using the
internal function \code{feat_ds_well_ref()}. The returned objects are of type
\code{FeatureInformation} and each contain code, label and description of each
feature. If for different datasets different sets of features are available,
\code{list_features()} provides the union of the features of all datasets.

Similarly, \code{list_feature_codes()} provides the list of all available
features as character vector of feature codes. As for \code{list_features()},
either plate-level or well-level feature dataset reference may be passed
and if a \code{wells} argument is supplied together with a plate-level
reference, the corresponding well-level references are constructed using
the internal function \code{feat_ds_well_ref()}. If for different datasets
different sets of features are available, \code{list_feature_codes()} provides
the union of the features of all datasets.

For a given set of feature vector datasets, \code{fetch_features()} fetches
feature data for the specified feature codes, or for all available features
in case the argument \code{feature_codes} is not specified (or \code{NA}). The
behavior regarding well selection is the same as in \code{list_features()} and
\code{list_feature_codes()}. Either plate-level or well-level dataset references
are passed and whenever plate-level references are passed in combination
with \code{WellPosition} object, the corresponding well-level references are
created and used. If for different datasets different sets of features are
available, the union of the features of all datasets is searched for. The
returned object is of type \code{FeatureVectorDataset}, which for each entry
contains a \code{FeatureVectorDatasetReference} and a set of \code{FeatureVector}(s),
one for each well.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "listAvailableFeatures")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "listAvailableFeatureCodes")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrs", "loadFeatures")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrs",
"loadFeaturesForDatasetWellReferences")}
}
}

\section{Implementation note}{
 Even though there exists a constructor for
\code{FeatureVectorDatasetWellReference} objects, which takes two arguments, one
for the corresponding \code{FeatureVectorDatasetReference} object and one for
a \code{WellPosition} objects, this does not work. Furthermore, class information
cannot be supplied as this will cause an error as well (hence the use of
\code{rm_json_class()}). Why the function \code{loadFeaturesForDatasetWellReferences}
behaves this way is currently unclear.
}

\examples{
\donttest{
  tok <- login_openbis()

  # search for a sample object corresponding to plate KB2-03-1I
  samp <- search_openbis(tok,
                         search_criteria(
                           attribute_clause("code",
                                            "/INFECTX_PUBLISHED/KB2-03-1I")
                         ),
                         target_object = "sample")
  # for the plate sample object, list all feature data sets
  feat_ref <- list_references(tok, samp, type = "feature")

  # several data set types can act as feature data sets
  get_field(feat_ref, "dataSetType")
  feat_ref <- feat_ref[[7L]]

  # for a feature data set, list all features
  feat_info <- list_features(tok, feat_ref)
  feat_info <- feat_info[c(2L, 6L)]

  # for a feature data set, a set of feature codes and a set of wells,
  # retrieve the corresponding feature data
  feats <- fetch_features(tok, feat_ref,
                          feature_codes = get_field(feat_info, "code"),
                          wells = well_pos(1:6, 3:8))

  well_pos <- get_field(feats, "wellPosition")
  values <- get_field(feats, "values")
  tibble::tibble(
    well_row = get_field(well_pos, "wellRow"),
    well_col = get_field(well_pos, "wellColumn"),
    cell_count = as.integer(values[1L, ]),
    cell_area = as.integer(values[2L, ])
  )

  logout_openbis(tok)
}

}
\seealso{
Other resource listing/downloading functions: \code{\link{fetch_images}},
  \code{\link{list_download_urls}},
  \code{\link{list_files}}
}
\concept{resource listing/downloading functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infx.R
\docType{package}
\name{infx}
\alias{infx}
\alias{infx-package}
\title{API access to the InfectX data repository}
\description{
The \href{https://labnotebook.ch}{openBIS} data repository hosted by
\href{http://www.infectx.ch}{InfectX} contains high throughput screening data
from several large-scale gene knockdown experiments. The screens currently
publicly available are RNA interference based, use kinome-wide libraries
from multiple vendors and were carried out on HeLa cells, in presence of
several viral and bacterial pathogens. Further genome-wide screens have been
carried out and their public release is forthcoming. For more information,
please refer to the \href{https://docs.ropensci.org/infx/}{README} or the
\href{../doc/infx-intro.html}{Introduction vignette}.
}
\details{
The provided functionality is not restricted to InfectX data, but applies
to the v1 JSON-RPC based openBIS API in general. Some parts of the API,
geared more towards data curation are currently not supported. For more
information on what API functions are available, have a look at the
\href{../doc/openbis-api.html}{openBIS API vignette}. The basic infrastructure
for creating and executing a request, as well as processing the response, is
exposed and missing functionality can easily be added.

Type information of JSON objects returned from the API is preserved as S3
class attribute and all retrieved JSON objects additionally inherit from the
S3 class \code{json_class}. As such, a \code{foobar} object retrieved from openBIS,
will have two class attributes: \code{foobar} and \code{json_class}. Sets of
\code{json_class} objects that are of the same sub-type can be represented as
\code{json_vec} objects of that sub-type. Several \code{foobar} objects therefore can
be combined into a list structure with S3 classes \code{foobar} and \code{json_vec},
where every entry in turn is an S3 object with types \code{foobar} and
\code{json_class}.\preformatted{examp <- json_vec(
  json_class(a = "foo", class = "foobar"),
  json_class(a = "bar", class = "foobar")
)
str(examp)

#> List of 2
#>  $ :List of 1
#>   ..$ a: chr "foo"
#>   ..- attr(*, "class")= chr [1:2] "foobar" "json_class"
#>  $ :List of 1
#>   ..$ a: chr "bar"
#>   ..- attr(*, "class")= chr [1:2] "foobar" "json_class"
#>  - attr(*, "class")= chr [1:2] "foobar" "json_vec"
}

Such an approach was chosen in order to not only have available generic
function dispatch on individual \code{json_class} objects, but also on sets (or
\emph{vectors}) of \code{json_class} objects. For more information on working with
\code{json_class} and \code{json_vec} objects refer to the
\href{#json-object-handling}{section on JSON objects} and
\href{../doc/json-class.html}{JSON object vignette}.

This documentation makes a distinction between objects in openBIS that exist
mainly for the purpose of organizing/grouping data and objects that
represent actual data resources. The most basic object in the organizational
hierarchy is that of a \code{Project}. Several \code{Experiment} objects may be
associated with a \code{Project} and \code{Sample} objects live in experiments. Given
the HTS-based context of InfectX data, samples typically represent
microtiter plates or individual plate wells. \code{Material} objects describe
agents applied to samples. Many of the InfectX screens are RNA interference-
based and therefore materials may be for example siRNA oligos or targeted
genes. Finally, samples are associated with \code{DataSet} objects that stand for
experimental measurements or data derived thereof.

Any type of data resource available in openBIS can be accessed as files
belonging to data sets. Due to the image-based nature of InfectX screens,
raw experimental data comes in the form of fluorescence microscopy imagery
which consequently constitutes the most basic form of data resource
available. It is therefore no surprise that image data receives special
treatment, allowing for more fine grained access and functionality that
helps with finding specific sub-sets of images. A further data resource
that comes with features similar to those of image data is termed feature
vector data sets. This is mostly tabular data with a single value
corresponding to an imaging site. This is typically used for image
acquisition meta data, summarized image analysis or quality control results.
}
\section{General comments}{

A login token is required for any type of API call. Passing valid login
credentials to \code{\link[=login_openbis]{login_openbis()}} will return a string that can subsequently
be used for making API requests. Login tokens are invalidated by calling
\code{\link[=logout_openbis]{logout_openbis()}} which is performed automatically upon garbage collection
of login tokens returned by \code{\link[=login_openbis]{login_openbis()}} with the \code{auto_disconnect}
switch set to \code{TRUE} (default). Validity of a login token can be checked
with \code{\link[=is_token_valid]{is_token_valid()}}.

All API requests are constructed by \code{\link[=make_requests]{make_requests()}} (or for single
requests by the wrapper function \code{\link[=make_request]{make_request()}}), which helps with putting
together JSON-RPC requests and parses the returned JSON objects by calling
\code{\link[=process_json]{process_json()}}. Processing of JSON involves generation of \code{json_class}
and \code{json_vec} objects using \code{@type} information, as well as resolution of
\code{@id} references. While obviously a feature for reducing data transfer
overhead, this type of data deduplication has the down-side of yielding
objects that are no longer self-contained. If for example plate wells are
listed and each well contains an object referencing the associated plate,
only a single instance of this plate object will be retrieved as part of the
first well object and all subsequent well objects only contain a reference
to this plate object. Sub-setting this list of wells however might yield
well objects with broken references. To circumvent such issues, all
references are resolved by a call to \code{\link[=resolve_references]{resolve_references()}}, initiated by
\code{\link[=process_json]{process_json()}}.

As a side note: while created for and mainly tested with
\href{https://infectx.biozentrum.unibas.ch/openbis}{InfectX} data, all API
methods can be used for accessing other openBIS instances as well.
Functions that issue API calls can all accept a \code{host_url} argument which
is forwarded to \code{\link[=api_url]{api_url()}} in \code{\link[=make_requests]{make_requests()}} in order to create API
endpoint urls. Another publicly available openBIS instance is the
\href{https://openbis-eln-lims.ethz.ch/openbis/}{demo} offered by the openBIS
development team. It can be accessed with both user name and password
\code{test_observer} both via a browser or by passing
\code{https://openbis-eln-lims.ethz.ch} as \code{host_url} to methods which
initiate API calls.

After being assembled by \code{\link[=make_requests]{make_requests()}}, requests are executed by
\code{\link[=do_requests_serial]{do_requests_serial()}} or \code{\link[=do_requests_parallel]{do_requests_parallel()}}, depending on whether
several API calls are constructed at the same time. The argument \code{n_con}
controls the degree of parallelism and if set to \code{1}, forces serial
execution even in cases where several requests are being issued. Failed
requests can be automatically repeated to provide additional stability by
setting the \code{n_try} argument to a value larger than \code{1} (default is \code{2}).
For more information on how to add further functionality using
\code{\link[=make_requests]{make_requests()}} and \code{\link[=do_requests_serial]{do_requests_serial()}}/\code{\link[=do_requests_parallel]{do_requests_parallel()}},
refer to the \href{../doc/openbis-api.html}{openBIS API vignette}.
}

\section{JSON object handling}{

Object structures as returned by openBIS can be instantiated using the
creator \code{\link[=json_class]{json_class()}}. This function takes an arbitrary set of key-value
pairs, followed by a class name and returns a list-based \code{json_class}
object. Existing list-based objects may be coerced to \code{json_class} using
\code{\link[=as_json_class]{as_json_class()}} where \code{@type} list entries are taken to be class types.
The inverse is achieved by calling \code{\link[=rm_json_class]{rm_json_class()}} on a \code{json_class}
object or by calling \code{\link[=as_list]{as_list()}} and passing the \code{keep_asis} argument as
\code{FALSE}. \code{json_class} objects can be validated with \code{\link[=is_json_class]{is_json_class()}} which
is recursively called on any object inheriting from \code{json_class} in
\code{\link[=check_json_class]{check_json_class()}}.

Similarly to \code{json_class} objects, a constructor for \code{json_vec} objects is
provided in the form of \code{\link[=json_vec]{json_vec()}} and existing structures can be coerced
to \code{json_vec} by \code{\link[=as_json_vec]{as_json_vec()}}. The validator function \code{\link[=is_json_vec]{is_json_vec()}}
tests whether an object is a properly formed \code{json_vec} object and the
utility function \code{\link[=has_common_subclass]{has_common_subclass()}} tests whether the passed list
structure consists of \code{json_class} objects of the same sub-type. The inverse
of applying \code{\link[=as_json_vec]{as_json_vec()}} to a list structure is achieved by passing a
\code{json_vec} object to \code{\link[=as_list]{as_list()}}.

Several utility functions are provided that facilitate handling of
\code{json_class} and \code{json_vec} objects. \code{\link[=has_fields]{has_fields()}} tests whether certain
named entries are present in a \code{json_class} object or in each member of a
\code{json_vec}. In order to extract the content of a field, \code{\link[=get_field]{get_field()}} can be
applied to \code{json_class} and \code{json_vec} objects. Analogously,
\code{\link[=has_subclass]{has_subclass()}} and \code{\link[=get_subclass]{get_subclass()}} test for and extract the original JSON
object type from \code{json_class} and \code{json_vec} objects. Finally,
\code{\link[=remove_null]{remove_null()}} recursively removes empty fields (fields containing \code{NULL})
from \code{json_class} and \code{json_vec} objects.

In addition to the mentioned utility functions, several base R generic
functions have \code{json_class} and \code{json_vec} specific methods implemented.
Combining several \code{json_class} objects using \code{\link[base:c]{base::c()}} yields a \code{json_vec}
object, as does repeating objects using \code{\link[base:rep]{base::rep()}}. The same functions
can be applied to \code{json_vec} objects but this only checks for agreement in
sub-type. Custom sum-setting is provided as well, in order to retain class
attributes and replacement functions acting on \code{json_vec} objects make sure
that sub-types remain compatible. Recursive printing of both \code{json_class}
and \code{json_vec} objects is possible by calling \code{\link[base:print]{base::print()}}. Recursion
depth, as well as printing length and width can be controlled via arguments,
as can fancy printing (colors and UTF box characters for visualizing tree
structures).
}

\section{Listing and searching for objects}{

OpenBIS projects can be listed by calling \code{\link[=list_projects]{list_projects()}} and experiments
are enumerated with \code{\link[=list_experiments]{list_experiments()}}. Two objects types are used for
representing experiments: \code{Experiment} and \code{ExperimentIdentifier}.
\code{\link[=as_experiment_id]{as_experiment_id()}} converts a set of \code{Experiment} objects to
\code{ExperimentIdentifier} (requires no API call) and the inverse is possible
by passing a set of \code{ExperimentIdentifier} objects to \code{\link[=list_experiments]{list_experiments()}}
(does require an API call). All available experiments can be listed as
\code{ExperimentIdentifier} objects using \code{\link[=list_experiment_ids]{list_experiment_ids()}} and all
experiments for a set of projects are enumerated by passing \code{Project}
objects to \code{\link[=list_experiments]{list_experiments()}}. Experiments have a type and all realized
types can be listed with \code{\link[=list_experiment_types]{list_experiment_types()}}.

Experiments consist of samples which can be listed by passing a set of
\code{Experiment} or \code{ExperimentIdentifier} objects to \code{\link[=list_samples]{list_samples()}}. Samples
too have a type and all types are retrieved by calling
\code{\link[=list_sample_types]{list_sample_types()}}. Additional object types that are used to represent
samples are plate and well objects, including \code{Plate}, \code{PlateIdentifier},
\code{PlateMetadata}, \code{WellIdentifier} and \code{WellMetadata}, all of which can be
converted to \code{Sample} objects by calling \code{\link[=list_samples]{list_samples()}}. Plate objects
can be listed using \code{\link[=list_plates]{list_plates()}}, which can either return all available
plate objects or plates for a given set of experiments (passed as
\code{Experiment} or \code{ExperimentIdentifier} objects). Plate meta data, which
also contains associated well meta data is retrieved by
\code{\link[=list_plate_metadata]{list_plate_metadata()}} which can act on plate objects (\code{Plate},
\code{PlateIdentifier} or \code{Sample}). Wells of a plate are listed with
\code{\link[=list_wells]{list_wells()}} which too may be dispatched on plate objects. Wells
associated with a material object can be enumerated by passing a set of
\code{MaterialScreening}, \code{MaterialIdentifierScreening}, \code{MaterialGeneric} or
\code{MaterialIdentifierGeneric} to \code{\link[=list_wells]{list_wells()}}.

Data set objects represent the most diverse group of data-organizational
structures. Possible types include \code{DataSet}, \code{DatasetIdentifier},
\code{DatasetReference}, \code{ImageDatasetReference}, \code{MicroscopyImageReference},
\code{PlateImageReference}, \code{FeatureVectorDatasetReference} and
\code{FeatureVectorDatasetWellReference}. Full \code{DataSet} objects are returned by
\code{\link[=list_datasets]{list_datasets()}}, either for a set of plate samples, experiments or data
set codes (passed as character vector). \code{\link[=list_dataset_ids]{list_dataset_ids()}} gives back
\code{DatasetIdentifier} objects, either for a set of \code{DataSet} objects or data
set codes (again passed as character vector). The remaining data set types
are generated by \code{\link[=list_references]{list_references()}}, and return type depends on input
arguments.

Whenever \code{\link[=list_references]{list_references()}} is dispatched on objects identifying a plate
sample (\code{Plate}, \code{PlateIdentifier}, \code{PlateMetadata} or \code{Sample}), a \code{type}
argument is available, which can be any of \code{raw}, \code{segmentation} or
\code{feature}. Depending on \code{type}, \code{ImageDatasetReference} or
\code{FeatureVectorDatasetReference} objects are returned. The former type of
objects represent plate-wise image data sets (either for raw images or
segmentation masks) while the latter type references feature vector data
sets.

Dispatch of \code{\link[=list_references]{list_references()}} is also possible on objects identifying
data sets and again the return type depends on further arguments. If
imaging channels are specified as \code{channels} argument, but not specific
wells are selected, \code{MicroscopyImageReference} objects are retrieved,
representing a plate-wide raw imaging data set per imaging site and imaging
channel. If in addition to imaging channels, wells are specified
(\code{WellPosition} objects, e.g. created by \code{\link[=well_pos]{well_pos()}}, passed as \code{wells}
argument), the return type changes to \code{PlateImageReference}. Such objects
precisely reference an image, by encoding imaging channel, imaging site,
well position and pate-wise imaging data set.

Finally, \code{\link[=list_references]{list_references()}} can be dispatched on material objects,
including \code{MaterialGeneric}, \code{MaterialScreening},
\code{MaterialIdentifierGeneric} and \code{MaterialIdentifierScreening}, in which case
\code{PlateWellReferenceWithDatasets} objects are returned. While themselves
not representing data sets, \code{PlateWellReferenceWithDatasets} contain all
respective \code{ImageDatasetReference} and \code{FeatureVectorDatasetReference}
objects.
}

\section{Search for objects}{

Instead of enumerating objects using the various \code{list_*()} functions,
search queries can be constructed and run against openBIS. A search query
consists of a possibly nested \code{SearchCriteria} object as instantiated by
\code{\link[=search_criteria]{search_criteria()}} and is executed by calling \code{\link[=search_openbis]{search_openbis()}}.
\code{SearchCriteria} objects are composed of a set of match clauses (see
\code{\link[=property_clause]{property_clause()}}, \code{\link[=any_property_clause]{any_property_clause()}}, \code{\link[=any_field_clause]{any_field_clause()}},
\code{\link[=attribute_clause]{attribute_clause()}} and \code{\link[=time_attribute_clause]{time_attribute_clause()}}) which are combined by
an operator (either \code{any} or \code{all}).

Additionally, a single \code{SearchSubCriteria} may be attached to every
\code{SearchCriteria} object which in turn consists of a \code{SearchCriteria} and an
object type to which this search criteria object is applied to. In the call
to \code{\link[=search_openbis]{search_openbis()}} a target type has to be specified as \code{target_object}
argument (default is \code{data_set} and possible alternatives are \code{experiment},
\code{material} as well as \code{sample}) to indicate what object type the search is
targeted at.
}

\section{Downloading data}{

As mentioned earlier, there are three types of data resources that can be
downloaded: files, images and feature vector data. File access is the most
basic method and any type of data (including images and feature data) is
available via this route. Accessing images and feature data using
specialized interfaces however simplifies and makes possibly more specific
data access.

Files can be listed for any object representing a data set as well as for a
character vector of data set codes using \code{\link[=list_files]{list_files()}}. An object type,
specialized for referencing files in a data set is available as
\code{DataSetFileDTO} can also be passed to \code{\link[=list_files]{list_files()}}. This is useful
whenever only a subset of files within a data set, contained in a folder,
are of interest. In any case, \code{\link[=list_files]{list_files()}} returns a set of
\code{FileInfoDssDTO} objects. As no data set information is encoded in
\code{FileInfoDssDTO} objects, \code{\link[=list_files]{list_files()}} saves data set codes as \code{data_set}
attributes with each object. Download of files is done using
\code{\link[=fetch_files]{fetch_files()}}, which requires for every requested file, the data set code
and file path. This information can be passed as separate character vectors,
\code{DataSetFileDTO} objects or \code{FileInfoDssDTO} objects with data set
information passed separately as character vector or as \code{data_set}
attribute with each object. Furthermore data set membership information can
be passed as any type of data set object and if no file paths are
specified, all available files for the given data sets are retrieved.

\code{\link[=fetch_files]{fetch_files()}} internally creates download urls by calling
\code{\link[=list_download_urls]{list_download_urls()}} and uses \code{\link[=do_requests_serial]{do_requests_serial()}} or
\code{\link[=do_requests_parallel]{do_requests_parallel()}} to execute the downloads. Whether downloads are
performed in serial or parallel fashion can be controlled using the \code{n_con}
argument. Additionally a function may be passed to \code{\link[=fetch_files]{fetch_files()}} as
\code{reader} argument which will be called on each downloaded file.

Images are retrieved using \code{\link[=fetch_images]{fetch_images()}}. If dispatch occurs on general
purpose data set objects, including \code{DatasetIdentifier}, \code{DatasetReference}
or \code{ImageDatasetReference}, further arguments for identifying images are
passed as \code{channels} and \code{well_positions}. As \code{MicroscopyImageReference}
objects already contain channel information, only well positions are needed
in order to specify images. Somewhat surprisingly, image tile information
which is also part of \code{MicroscopyImageReference} objects is disregarded and
images are fetched for entire wells. Data sets that are connected to wells
and not plates can be passed to \code{\link[=fetch_images]{fetch_images()}} without additionally
specifying well locations. Images can be scaled down to smaller sizes either
by setting the \code{thumbnails} argument to \code{TRUE} (only possible for data sets
connected to wells instead of plates, as the corresponding API call does
not support selecting wells) or by passing an \code{ImageSize} object as
\code{image_size} argument, in which case returned images will be scaled to fit
within the box specified by the \code{ImageSize} object, while retaining the
original aspect ratio.

\code{PlateImageReference} objects most precisely reference images, as they
contain data set, well location, site location and channel information. If
a set of \code{PlateImageReference} objects is passed to \code{\link[=fetch_images]{fetch_images()}}, image
size can be set using the \code{thumbnails} or \code{image_size} arguments and image
file type can be forced to png using the \code{force_png} switch. Most
fine-grained control over the returned images is achieved by using
\code{ImageRepresentationFormat} objects. Pre-defined format objects can be
retrieved per data set by calling \code{\link[=list_image_metadata]{list_image_metadata()}} with \code{type} set to
\code{format}. General image meta data, such as tile layout and channel
information is returned by \code{\link[=list_image_metadata]{list_image_metadata()}} if the \code{type} argument
is left at default value \code{metadata}.

Two types of objects are central to specifying feature data sets:
\code{FeatureVectorDatasetReference} and \code{FeatureVectorDatasetWellReference}
where the former object type references feature data for an entire plate and
the latter for individual wells on a plate. Both object types may be passed
to \code{\link[=fetch_features]{fetch_features()}} which returns objects of type \code{FeatureVectorDataset}
whenever a full plate is requested and \code{FeatureVectorWithDescription} for
individual wells. Features are selected by passing a character vector of
feature codes as \code{feature_codes} argument, the possible values of which
can be enumerated for a feature vector data set by calling
\code{\link[=list_feature_codes]{list_feature_codes()}} or by extracting the \code{code} entries from
\code{FeatureInformation} objects as retrieved by \code{\link[=list_features]{list_features()}}. In case the
\code{feature_codes} argument is left at default value (\code{NA}), all available
features are returned by \code{\link[=fetch_features]{fetch_features()}}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/url.R
\name{api_url}
\alias{api_url}
\alias{docs_link}
\title{OpenBIS urls}
\usage{
api_url(api_endpoint = c("gis", "gics", "qas", "wis", "dsrg", "sas",
  "dsrs"), host_url = "https://infectx.biozentrum.unibas.ch",
  full_url = NULL, ...)

docs_link(api_endpoint = c("gis", "gics", "qas", "wis", "dsrg", "sas",
  "dsrs"), method_name = NULL, version = "13.04.0")
}
\arguments{
\item{api_endpoint}{Abbreviated name of the API section (e.g. \code{gis} for
IGeneralInformationService).}

\item{host_url}{Host url.}

\item{full_url}{Instead of constructing the API endpoint url, a string can
be passed which will be returned again.}

\item{...}{Further arguments are ignored.}

\item{method_name}{Name of the method for which the link is created.}

\item{version}{OpenBIS version number.}
}
\value{
A character vector of length 1.
}
\description{
The helper function \code{api_url()} is used to create urls to InfectX
openBIS API endpoints and \code{docs_link()} generates latex links to
documentation pages. Both functions support the following endpoints which
are abbreviated as
\itemize{
\item \code{gis}: \Sexpr[results=rd]{infx::docs_link("gis")}
\item \code{gics}: \Sexpr[results=rd]{infx::docs_link("gics")}
\item \code{qas}: \Sexpr[results=rd]{infx::docs_link("qas")}
\item \code{wis}: \Sexpr[results=rd]{infx::docs_link("wis")}
\item \code{dsrg}: \Sexpr[results=rd]{infx::docs_link("dsrg")}
\item \code{sas}: \Sexpr[results=rd]{infx::docs_link("sas")}
\item \code{dsrs}: \Sexpr[results=rd]{infx::docs_link("dsrs")}
}
}
\details{
If for some reason an url is desired from \code{api_url()} that cannot be
constructed by pasting \code{host_url} and one of the hard-coded API endpoints
together, this can be passed as \code{full_url}, which will simply be returned.
}
\examples{
# default endpoint is the GeneralInformationService interface
api_url()
# base url can be customized
api_url(host_url = "https://foobar.xyz")
# ScreeningApiServer interface endpoint
api_url("sas")
# manual url
api_url(full_url = "https://foobar.xyz/openbis/new-api-section-v1.json")

# link to GeneralInformationService interface docs
docs_link()
# add a method name (only to the link text)
docs_link(method_name = "foo_bar")
# link to ScreeningApiServer interface docs
docs_link("sas")
# link to most recent version of docs
docs_link("sas", version = "16.05.6")

}
\seealso{
Other utility functions: \code{\link{login_openbis}},
  \code{\link{make_requests}}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiment.R
\name{list_experiments}
\alias{list_experiments}
\alias{list_experiments.ExperimentIdentifier}
\alias{list_experiments.Project}
\alias{list_experiment_ids}
\alias{list_experiment_types}
\alias{as_experiment_id}
\alias{as_experiment_id.Experiment}
\alias{as_experiment_id.ExperimentIdentifier}
\title{List experiments}
\usage{
list_experiments(token, x, ...)

\method{list_experiments}{ExperimentIdentifier}(token, x, ...)

\method{list_experiments}{Project}(token, x,
  types = list_experiment_types(token), require = c(NA, "DataSets",
  "Samples"), ...)

list_experiment_ids(token, ...)

list_experiment_types(token, ...)

as_experiment_id(x, ...)

\method{as_experiment_id}{Experiment}(x, ...)

\method{as_experiment_id}{ExperimentIdentifier}(x, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to limit the number of returned experiments, e.g. a set of
\code{ExperimentIdentifier} or \code{Project} objects.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}

\item{types}{Either a single or set of \code{ExperimentType} objects}

\item{require}{A switch to require the resulting experiments to contain a
nonzero number of dataset or sample. Default behavior is no requirement.}
}
\value{
Depending on the number of resulting objects, either a
\code{\link{json_class}} (single object) or a \code{\link{json_vec}} (multiple objects), is
returned. Experiments are represented by \code{Experiment}, experiment ids by
\code{ExperimentIdentifier} and experiment types by \code{ExperimentType} objects.
}
\description{
Experiments can be represented by two object classes: \code{Experiment} and
\code{ExperimentIdentifier}. The fields that make up an \code{ExperimentIdentifier}
object are a subset of those required for an \code{Experiment} object. Therefore
an experiment can be turned into an experiment id object without an API
call, using the function \code{as_experiment_id()}. The reverse can be achieved
by calling \code{list_experiments()} on experiment id objects. In general,
experiments and experiment id objects can be listed using
\code{list_experiments()} and \code{list_experiment_ids()}.
}
\details{
By calling \code{list_experiments()} on project objects, all corresponding
experiments are listed. It is possible to limit the search to experiments
that are of a certain type and exclude experiments that either have no
datasets or samples associated. An exhaustive list of realized experiment
types can be retrieved using \code{list_experiment_types()}. If several
experiment types are requested in \code{list_experiments()}, the default is to
iterate over all available types, an API call per experiment type has to be
made.

\code{ExperimentIdentifier} objects present a more compact way of uniquely
representing experiments. All experiments that are available to the current
user can be listed with \code{list_experiment_ids()}. There is no way of limiting
the search for experiments.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listExperiments")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis",
"listExperimentsHavingDataSets")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "listExperimentsHavingSamples")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("sas", "listExperiments")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listExperimentTypes")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  # list all available projects to limit the search for experiments
  proj <- list_projects(tok)

  # list all experiments corresponding to the project with index 1
  exps <- list_experiments(tok, proj[[1L]])

  # convert experiment to experiment ids and back
  exp_ids <- as_experiment_id(exps)
  identical(exps, list_experiments(tok, exp_ids))

  # experiments can also be searched for
  exp <- search_openbis(tok,
                        search_criteria(
                          attribute_clause("code",
                                           get_field(exps[[1L]],
                                                     "identifier"))
                        ),
                        target_object = "experiment")
  identical(exps[[1L]], exp)

  logout_openbis(tok)
}

}
\seealso{
Other object listing functions: \code{\link{list_datasets}},
  \code{\link{list_material}}, \code{\link{list_plates}},
  \code{\link{list_projects}}, \code{\link{list_samples}}
}
\concept{object listing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{search_openbis}
\alias{search_openbis}
\alias{search_criteria}
\alias{search_sub_criteria}
\alias{property_clause}
\alias{any_property_clause}
\alias{any_field_clause}
\alias{attribute_clause}
\alias{time_attribute_clause}
\alias{list_property_types}
\title{Assemble and execute openBIS search queries}
\usage{
search_openbis(token, criteria, target_object = c("data_set",
  "experiment", "material", "sample"), fetch_options = NULL, ...)

search_criteria(..., operator = "all", sub_criteria = NULL)

search_sub_criteria(criteria, type = "sample")

property_clause(property_code, value, ...)

any_property_clause(value, ...)

any_field_clause(value, ...)

attribute_clause(attribute = "code", value, ...)

time_attribute_clause(attribute = "registration", value = Sys.Date(),
  timezone = 0L, ...)

list_property_types(token, with_relations = FALSE, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{criteria}{A single \code{SearchCriteria} object.}

\item{target_object}{The object type the search is targeted at, i.e.
\code{DataSet}s, \code{Experiment}s, etc.}

\item{fetch_options}{If samples are searched for, additional fetch options
can be specified.}

\item{...}{For \code{search_openbis()} passed to \code{\link[=make_request]{make_request()}}, for
\code{search_criteria()} a set of match clause objects, and for match clause
constructors, the comparison mode can be passed as \code{mode} argument, which
may be \code{eq} (==), \code{lte} (<=) or \code{gte} (>=).}

\item{operator}{How to combine search clauses, either \code{all} or \code{any} have
to be fulfilled.}

\item{sub_criteria}{Optionally, one or several \code{SearchSubCriteria} objects
can be used to create a \code{SearchCriteria} object.}

\item{type}{The entity type, a \code{SearchSubCriteria} is applied to.}

\item{property_code}{Code identifying a property to be used for the
comparison.}

\item{value}{The value used in the comparison.}

\item{attribute}{Name of the attribute to be used for the comparison.}

\item{timezone}{A string identifying the timezone of the specified date,
examples include "+1", "-5", "0", etc.}

\item{with_relations}{Logical switch indicating whether relations should
be returned as well.}
}
\value{
Depending on the number of resulting objects, either a
\code{\link{json_class}} (single object) or a \code{\link{json_vec}} (multiple objects), is
returned. The object specific sub-class for objects returned by
\code{search_openbis()} depends on the target object type specified by the
\code{target_object} argument. The remaining functions are used for constructing
search queries and instantiate the required objects: \code{search_criteria()}
returns a \code{SearchCriteria} and \code{search_sub_criteria()} a \code{SearchSubCriteria}
object, while \code{property_clause()} returns a \code{PropertyMatchClause},
\code{any_property_clause()} an \code{AnyPropertyMatchClause}, \code{any_field_clause()}
an \code{AnyFieldMatchClause}, \code{attribute_clause()} an \code{AttributeMatchClause} and
\code{time_attribute_clause()} a \code{TimeAttributeMatchClause} object. Finally,
\code{list_property_types()} returns a list with two slots, one holding a
\code{json_vec} of sub-type \code{ControlledVocabularyPropertyType} and the other
holding a \code{json_vec} of sub-type \code{PropertyType}.
}
\description{
Searching in openBIS presents a powerful alternative to iterative listing
and selecting of objects. As an example, in order to find image data sets
associated with an experiment, instead of first listing all experiments,
selecting the one of interest, then listing all plates of that experiment,
followed by listing image data sets for each of the plates, the requested
data sets can be directly retrieved by constructing a search query with
\code{search_criteria()} and executing the search by calling \code{search_openbis()}.
}
\details{
Searching openBIS can be done by creating a \code{SearchCriteria} object and
passing that to \code{search_openbis()}, alongside specifying what type of
object is being searched for. In case of \code{Sample}s being searched for,
a further argument, \code{fetch_options}, can be specified for controlling the
search, which can contain one or more of the strings
\itemize{
\item \code{ancestors}: Ask for all ancestors.
\item \code{basic}: Samples will have only basic attributes (id, code, type, space
code, experiment identifier, registrator, registration date,
modification date) but no properties.
\item \code{children}: Samples contain also their children samples.
\item \code{contained}: Ask for contained samples.
\item \code{descendants}: Ask for all descendants.
\item \code{metaprojects}: Ask for metaprojects this sample belongs to.
\item \code{parents}: Samples contain also their parent samples.
\item \code{properties}: Samples contain basic attributes and all properties.
}

A \code{SearchCriteria} object can be instantiated using the constructor
\code{search_criteria()}, which takes one or several match clause objects, a
search operator specifying whether to match \code{all} or \code{any} clauses, and
optionally one or several \code{SearchSubCriteria} objects. \code{SearchSubCriteria}
objects in turn can be created with \code{search_sub_criteria()}, which takes a
single \code{SearchCriteria} object alongside a string specifying the entities,
the sub criterion is applied to. Possibilities are
\itemize{
\item \code{data_set_container}
\item \code{data_set_parent}
\item \code{data_set_child}
\item \code{experiment}
\item \code{sample}
\item \code{sample_container}
\item \code{sample_child}
\item \code{sample_parent}
}

\code{SearchCriteria} objects, used for searching openBIS, contain one or
several match clauses. A match clause, broadly speaking, consists of a
desired value, a field to which this value is compared to and a comparison
operator (e.g. equality). Match clauses can be constructed using any of
\code{attribute_clause()}, \code{time_attribute_clause()}, \code{property_clause()},
\code{any_property_clause()}, and \code{any_field_clause()}. Attribute match clauses
have a fixed set of attributes against which the match is performed:
\itemize{
\item time attribute match clauses
\itemize{
\item \code{registration_date}
\item \code{modification_date}
}
\item attribute match clauses
\itemize{
\item \code{code}
\item \code{type}
\item \code{perm_id}
\item \code{space}
\item \code{project}
\item \code{project_perm_id}
\item \code{metaproject}
\item \code{registrator_user_id}
\item \code{registrator_first_name}
\item \code{registrator_last_name}
\item \code{registrator_email}
\item \code{modifier_user_id}
\item \code{modifier_first_name}
\item \code{modifier_last_name}
\item \code{modifier_email}
}
}

In order to determine the possible values that can be supplied to
\code{property_clause()} as \code{property_code}s, \code{list_property_types()} can be
called. This function returns all property types available throughout the
queried openBIS instance. As objects of several types
(\code{ControlledVocabularyPropertyType} and \code{PropertyType}) are returned as
property types by the API, the resulting object is a list with each entry
corresponding to a type and holding a set of objects of the respective type.

The comparison operator (default is equality) can be any of the following
\itemize{
\item \code{equals}
\item \code{less_than_or_equal}, with alias \code{lte}
\item \code{greater_than_or_equal}, with alias \code{gte}
}

All of the option matching is not case-sensitive and is performed with
\code{\link[base:match.arg]{base::match.arg()}} and therefore options may be abbreviated (e.g. \code{eq}
instead of \code{equals}).
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "searchForDataSets")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "searchForExperiments")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "searchForMaterials")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "searchForSamples")}
}
}

\examples{
\donttest{
  tok <- login_openbis()
  
  # search for an experiment, e.g. ADENO-AU-K1
  exp <- search_openbis(tok,
                        search_criteria(
                          property_clause("pathogen", "Adenovirus"),
                          property_clause("library", "Ambion"),
                          property_clause("geneset", "Kinome"),
                          property_clause("replicate", 1L)
                        ),
                        target_object = "experiment")

  # the same can be achieved using the code attribute
  identical(exp,
            search_openbis(tok,
                           search_criteria(
                             attribute_clause(value = "ADENO-AU-K1")
                           ),
                           target_object = "experiment"))

  # or using the perm_id attribute
  identical(exp,
            search_openbis(tok,
                           search_criteria(
                             attribute_clause("perm_id",
                                              "20111223100933426-318017")
                           ),
                           target_object = "experiment"))

  # a search with no matches returns an empty list
  search_openbis(tok,
                 search_criteria(attribute_clause(value = "foo_bar")),
                 target_object = "experiment")

  # search using sub-criteria: all plate samples of experiment ADENO-DU-K1
  sub <- search_sub_criteria(search_criteria(
                               attribute_clause(value = "ADENO-DU-K1")
                             ),
                             type = "experiment")
  all <- search_openbis(
    tok,
    search_criteria(
      attribute_clause("type", "PLATE"),
      sub_criteria = sub
    ),
    target_object = "sample"
  )
  length(as_json_vec(all))

  # now only include ADENO-DU-K1 plate samples registered after Feb 1st 2013
  some <- search_openbis(
    tok,
    search_criteria(
      attribute_clause("type", "PLATE"),
      time_attribute_clause(value = as.Date("2013-02-01"), mode = "gte"),
      sub_criteria = sub
    ),
    target_object = "sample"
  )
  length(as_json_vec(some))

  logout_openbis(tok)
}

}
\concept{object searching functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.R
\name{dataset_code}
\alias{dataset_code}
\alias{dataset_code.DataSet}
\alias{dataset_code.DatasetIdentifier}
\alias{dataset_code.DatasetReference}
\alias{dataset_code.FeatureVectorDatasetReference}
\alias{dataset_code.FeatureVectorDatasetWellReference}
\alias{dataset_code.ImageDatasetReference}
\alias{dataset_code.MicroscopyImageReference}
\alias{dataset_code.PlateImageReference}
\alias{dataset_code.DataSetFileDTO}
\title{Extract dataset code}
\usage{
dataset_code(x, ...)

\method{dataset_code}{DataSet}(x, ...)

\method{dataset_code}{DatasetIdentifier}(x, ...)

\method{dataset_code}{DatasetReference}(x, ...)

\method{dataset_code}{FeatureVectorDatasetReference}(x, ...)

\method{dataset_code}{FeatureVectorDatasetWellReference}(x, ...)

\method{dataset_code}{ImageDatasetReference}(x, ...)

\method{dataset_code}{MicroscopyImageReference}(x, ...)

\method{dataset_code}{PlateImageReference}(x, ...)

\method{dataset_code}{DataSetFileDTO}(x, ...)
}
\arguments{
\item{x}{Dataset object(s).}

\item{...}{Generic compatibility.}
}
\description{
Given a (set of) dataset object(s), for each one extract the dataset code
and return a character vector of codes.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/request.R
\name{make_requests}
\alias{make_requests}
\alias{make_request}
\alias{do_requests_serial}
\alias{do_requests_parallel}
\alias{process_json}
\alias{resolve_references}
\title{Construct and issue JSON-RPC requests}
\usage{
make_requests(urls, methods, params, ids = NULL, version = "2.0",
  n_con = 5L, finally = process_json, ...)

make_request(url, method, params, ...)

do_requests_serial(urls, bodies = vector("list", length(urls)),
  n_try = 2L, create_handle = create_default_handle,
  check = check_default_result, finally = identity, ...)

do_requests_parallel(urls, bodies = vector("list", length(urls)),
  n_con = 5L, n_try = 2L, create_handle = create_default_handle,
  check = check_default_result, finally = identity, ...)

process_json(x)

resolve_references(x)
}
\arguments{
\item{params}{A list structure holding the arguments which, converted to
JSON, will be used to call the supplied method. The \code{@type} entries will be
generated from \code{json_class} attributes.}

\item{ids}{Identifier(s) for the JSON-RPC request (defaults to a random
string of length 7). Can be usually be ignored, as only single JSON-RPC
requests are issued per HTTP request.}

\item{version}{JSON-RPC protocol version to be used (defaults to \code{"2.0"}.}

\item{n_con}{The number of simultaneous connections.}

\item{finally}{A function that is applied to the \code{result} entry of a
successful JSON RPC request.}

\item{...}{Further arguments to \code{make_request()} are passed to
\code{make_requests()} and from \code{make_requests()} to \code{do_requests_serial()} or
\code{do_requests_parallel()}.}

\item{url, urls}{Destination url(s), the request is sent to.}

\item{method, methods}{The API method name(s).}

\item{bodies}{Request bodies: a list where each entry is a list with slots
\code{id}, \code{jsonrpc}, \code{method} and \code{params}.}

\item{n_try}{Number of tries each request is performed in case of failed
requests.}

\item{create_handle}{A function that will receive a single entry of the
\code{bodies} list at the time and should return a curl handle created by
\code{\link[curl:new_handle]{curl::new_handle()}}.}

\item{check}{A function that receives both the result of a request and the
corresponding entry of the \code{bodies} list. Is expected to return NULL in
which case the request is retried or a list with an entry named \code{result}.}

\item{x}{A (possibly nested) list structure for which all \code{@type} fields
are turned into class attributes and \code{@id} fields are recursively removed.}
}
\value{
The return type of \code{make_request()}/\code{make_requests()} and
\code{do_requests_serial()}/\code{do_requests_parallel()} depends on the callback
functions passed as \code{check} and \code{finally} arguments. At default,
\code{make_request()}/\code{make_requests()} return either a \code{\link{json_class}} (single
object) or a \code{\link{json_vec}} (multiple objects), dependent on the number of
resulting objects, while \code{do_requests_serial()}/\code{do_requests_parallel()}
return a list of raw vectors, one per request.
}
\description{
The functions powering all HTTP requests to openBIS are
\code{do_requests_serial()} for sequential calls and \code{do_requests_parallel()}
for asynchronous calls. Constructing requests for the JSON-RPC API of
openBIS is done with the helper function \code{make_requests()} and a wrapper
for single requests is available as \code{make_request()}. This convenience
function calls \code{make_requests()} and returns the first element of the
resulting list of length 1.
}
\details{
Both \code{do_requests_serial()} and \code{do_requests_parallel()} take as \code{urls}
argument a set of urls, either as character vector or list of unevaluated
function calls which will be evaluated using \code{\link[base:eval]{base::eval()}} shortly before
being used (this is used for urls that are only valid for a limited amount
of time). The behavior of \code{do_requests_*()} can be customized with the
three functions passed as arguments \code{create_handle}, \code{check} and \code{finally}
together with the vector (of the same length as \code{urls}) passed as argument
\code{bodies}.

The function passed as \code{create_handle} receives one entry at the time of
the \code{bodies} object and is expected to return a curl handle created by
\code{\link[curl:new_handle]{curl::new_handle()}}. The \code{check} function receives as first argument
the response of a single curl request alongside the corresponding entry of
the \code{bodies} object. This function should check whether the request was
successful or not, e.g. check the HTTP status code, size of the downloaded
file, etc. In case of failure it should return a \code{simpleError} error object,
created by \code{\link[base:simpleError]{base::simpleError()}} and in case of success it should return
the response data, e.g. the \code{content} entry of a curl response. The third
function, \code{finally}, is applied to the object returned by the \code{check}
function (in case of success) and can be used to parse JSON, read a binary
file, etc.

Both \code{do_requests_serial()} and \code{do_requests_parallel()} have the option of
retrying failed requests and the number of allowed retries can be controlled
with the argument \code{n_try}. Furthermore, \code{do_requests_parallel()} offers
control over the number of simultaneous connections using the argument
\code{n_con}. Only \code{n_con} requests are initially added to the curl multi handle
and for each successful one an additional request is added. This
implementation helps with urls that have a limited lifetime.

The function \code{make_requests()} is used to construct JSON-RPC requests. The
arguments \code{methods}, \code{params}, \code{ids} and \code{version} are combined into one
or several request objects according to the JSON-RPC specification which in
turn are passed to \code{do_requests_*()}. The objects passed as \code{methods} and
\code{params} should all be of the same length but in case any are of length 1,
they will be \code{\link[base:rep]{base::rep()}}'ed to the required length. Care has to be taken
that the list passed as \code{params} has the correct degree of nesting. As
\code{make_requests()} iterates over the topmost list level, a single request
should be wrapped in a list such that the topmost list level is of unit
length. The function \code{make_request()} is a wrapper around \code{make_requests()}
that does exactly this.

As part of the \code{process_json()} function, which is the default value passed
as \code{finally} argument in \code{make_requests()}, \code{@type} fields are converted
to \code{json_class} attributes, using \code{\link[=as_json_class]{as_json_class()}}. Additionally, \code{@id}
fields, which may be referenced if an object is used multiple times, are
recursively resolved using \code{\link[=resolve_references]{resolve_references()}} such that each object is
self-contained.
}
\examples{
\donttest{
  tok <- login_openbis()

  # the function list_projects() is implemented as follows
  projects <- make_request(api_url(api_endpoint = "gis"),
                           "listProjects",
                           list(tok))
  print(projects[[1L]])
  # alternatively, using make_requests(), the params argument has to be a
  # list per request and the first entry of the returned list has to be
  # selected
  proj <- make_requests(api_url(api_endpoint = "gis"),
                        "listProjects",
                        list(list(tok)))
  identical(proj[[1L]][[1L]],
            projects[[1]])

  # without using make_request(), one can achieve the same result by
  # calling do_requests_serial() directly
  proj <- do_requests_serial(api_url(api_endpoint = "gis"),
                             list(
                               list(
                                 id = "foobar",
                                 jsonrpc = "2.0",
                                 method = "listProjects",
                                 params = list(tok)
                               )
                             ),
                             create_handle = infx:::create_request_handle,
                             check = infx:::check_request_result,
                             finally = process_json)
  identical(proj[[1L]][[1L]],
            projects[[1]])

  # the do_requests_*() functions can be used for any HTTP request
  req <- do_requests_serial("https://httpbin.org/headers")
  is.raw(req[[1L]])

  # in order to read the returned binary data, a function can be supplied
  # to do_requests_*() as finally argument

  process_json <- function(x) 
    jsonlite::prettify(rawToChar(x))

  req <- do_requests_serial("https://httpbin.org/headers",
                            finally = process_json)
  req[[1L]]

  # to customize the curl handle, a function can be supplied to
  # do_requests_*() as create_handle argument

  post_handle <- function(x)
    curl::handle_setheaders(
      curl::new_handle(postfields = charToRaw(jsonlite::toJSON(x))),
      "Content-Type" = "application/json"
    )
  process_json <- function(x) 
    jsonlite::fromJSON(rawToChar(x), simplifyDataFrame = FALSE)$json

  data <- list(a = "foo",
               b = "bar")

  req <- do_requests_serial("https://httpbin.org/post",
                            list(data),
                            create_handle = post_handle,
                            finally = process_json)

  # httpbin returns POST data, therefore
  identical(data, req[[1L]])

  logout_openbis(tok)
}

}
\seealso{
Other utility functions: \code{\link{api_url}},
  \code{\link{login_openbis}}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.R
\name{list_img_ref}
\alias{list_img_ref}
\alias{list_img_ref.NULL}
\alias{list_img_ref.WellPosition}
\title{List image references}
\usage{
list_img_ref(token, x, wells = NULL, channels, ...)

\method{list_img_ref}{NULL}(token, x, wells, channels, ...)

\method{list_img_ref}{WellPosition}(token, x, wells, channels, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to limit search for datasets/files with.}

\item{wells}{A (set of) \code{WellPosition} object(s) to limit the dataset
listing to.}

\item{channels}{A character vector with imaging channel names to limit the
dataset listing to.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}
}
\description{
Used for double dispatching on the \code{list_datasets()} generic, list image
reference objects either for a specific (set of) \code{WellPosition} object(s)
or for the specified datasets in general.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/json-vec.R
\name{json_vec}
\alias{json_vec}
\alias{as_json_vec}
\alias{as.json_vec}
\alias{as_json_vec.json_vec}
\alias{as_json_vec.json_class}
\alias{as_json_vec.list}
\alias{as_json_vec.default}
\alias{as.list.json_vec}
\alias{is_json_vec}
\alias{is.json_vec}
\alias{has_common_subclass}
\title{Construct and validate JSON object vectors}
\usage{
json_vec(..., .simplify = FALSE)

as_json_vec(x, ...)

as.json_vec(x, ...)

\method{as_json_vec}{json_vec}(x, simplify = FALSE, ...)

\method{as_json_vec}{json_class}(x, simplify = FALSE, ...)

\method{as_json_vec}{list}(x, recursive = TRUE, force = FALSE,
  simplify = FALSE, ...)

\method{as_json_vec}{default}(x, force = FALSE, ...)

\method{as.list}{json_vec}(x, recursive = FALSE, ...)

is_json_vec(x)

is.json_vec(x)

has_common_subclass(x)
}
\arguments{
\item{...}{Individual \code{json_class} objects, or generic compatibility. Might
be passed on to \code{as.list()} for \code{json_class} objects.}

\item{x}{A single/list of \code{json_class} object(s), or other object to coerce}

\item{simplify, .simplify}{Logical switch indicating whether to simplify
\code{json_vec} objects of length 1 to \code{json_class} objects.}

\item{recursive}{Recursively apply the function.}

\item{force}{Suppress error when casting an object to \code{json_vec} that
cannot be converted.}
}
\value{
Multiple \code{\link{json_class}} objects of the same sub-type can be
represented as S3 objects with type \code{json_vec} and the common sub-type as
second class attribute.
}
\description{
In order to allow method dispatch on a set of \code{json_class} objects without
resorting to iterating over the individual set members, vectors of
\code{json_class} objects are wrapped by a \code{json_vec} class. Iterating over
objects is in some cases inefficient because the openBIS API can for some
functions accept lists of objects. Assembling multiple \code{json_class} objects
as a list in R however breaks method dispatch, as the type of this object
is \code{list} instead of the desired \code{json_class} sub-class. A \code{json_vec} object
therefore represents a list of \code{json_class} objects of the same sub-class
and brings this sub-class to the surface of the compound object.
}
\details{
A \code{json_vec} object can be instantiated using the \code{json_vec()} constructor
which takes a list of \code{json_class} objects of the same sub-class. An
existing list of \code{json_class} objects can be coerced to \code{json_vec} using
\code{as_json_vec()}/\code{as.json_vec()} and applying \code{as_list()}/\code{as.list()} to a
\code{json_vec} object reverses the action of \code{as_json_vec()} by removing all
\code{json_vec} related class information.

The function \code{is_json_vec()} and its alias \code{is.json_vec()} can be used to
test whether an object is a proper \code{json_vec} object. This requires that
\itemize{
\item all child elements have to be of the same sub-class
\item all child elements are required to be properly formed \code{json_class}
objects
\item the \code{json_vec} class attribute has to be in last position
\item the remaining class attributes have to be equal to the common sub-class
determined for the children.
}

Testing whether a list structure consists of \code{json_class} objects which are
of the same sub-class can be done with \code{has_common_subclass()}. This always
returns \code{TRUE} if a \code{json_class} object is passed and \code{FALSE} if a non-list
structure is passed.
}
\examples{
a <- json_class(field = "a", class = "foo")
b <- json_class(field = "b", class = "foo")

ab <- json_vec(a, b)

print(ab)

identical(ab, as_json_vec(list(a, b)))
# as_json_vec() is idempotent
identical(as_json_vec(list(a, b)),
          as_json_vec(as_json_vec(list(a, b))))

# a json_class object can be turned into a json_vec of length 1
ab_class <- json_class(foo1 = a, foo2 = b, class = "bar")
length(ab_class)
ab_vec <- as_json_vec(ab_class)
length(ab_vec)
# this can be reversed using as_json_class()
identical(ab_class, as_json_class(ab_vec))
# this might not be desirable in all cases; the argument simplify can be
# used to only create json_vec objects of length greater than 1
identical(as_json_vec(list(a), simplify = TRUE),
          a)

# has_common_subclass() will alway return true for json_class objects
has_common_subclass(a)
# list-based objects are tested
has_common_subclass(list(a, b))
# this includes json_vec objects
has_common_subclass(ab)
# each list entry has to be a json_class object
has_common_subclass(list("a", "b"))
# here sub-classes are "foo" and "bar"
has_common_subclass(list(ab_class, a))

is_json_vec(a)
is_json_vec(list(a, b))
is_json_vec(ab)

}
\seealso{
Other json object handling functions: \code{\link{has_fields.json_class}},
  \code{\link{json_class}}, \code{\link{print.json_class}}
}
\concept{json object handling functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/material.R
\name{list_material}
\alias{list_material}
\alias{list_material.MaterialIdentifierGeneric}
\alias{list_material.MaterialIdentifierScreening}
\alias{list_material.PlateIdentifier}
\alias{list_material.PlateMetadata}
\alias{list_material.Plate}
\alias{list_material.Sample}
\alias{material_id}
\alias{list_material_types}
\alias{as_screening_mat_id}
\alias{as_screening_mat_id.MaterialGeneric}
\alias{as_screening_mat_id.MaterialScreening}
\alias{as_screening_mat_id.MaterialIdentifierGeneric}
\alias{as_screening_mat_id.MaterialIdentifierScreening}
\alias{as_generic_mat_id}
\alias{as_generic_mat_id.MaterialGeneric}
\alias{as_generic_mat_id.MaterialScreening}
\alias{as_generic_mat_id.MaterialIdentifierGeneric}
\alias{as_generic_mat_id.MaterialIdentifierScreening}
\alias{extract_well_material}
\title{List materials}
\usage{
list_material(token, x, ...)

\method{list_material}{MaterialIdentifierGeneric}(token, x, ...)

\method{list_material}{MaterialIdentifierScreening}(token, x, ...)

\method{list_material}{PlateIdentifier}(token, x, material_type = NULL,
  ...)

\method{list_material}{PlateMetadata}(token, x, material_type = NULL,
  ...)

\method{list_material}{Plate}(token, x, material_type = NULL, ...)

\method{list_material}{Sample}(token, x, material_type = NULL, ...)

material_id(code, type = "gene", mode = c("screening", "generic"))

list_material_types(mode = c("screening", "generic"), types = NULL)

as_screening_mat_id(x, ...)

\method{as_screening_mat_id}{MaterialGeneric}(x, ...)

\method{as_screening_mat_id}{MaterialScreening}(x, ...)

\method{as_screening_mat_id}{MaterialIdentifierGeneric}(x, ...)

\method{as_screening_mat_id}{MaterialIdentifierScreening}(x, ...)

as_generic_mat_id(x, ...)

\method{as_generic_mat_id}{MaterialGeneric}(x, ...)

\method{as_generic_mat_id}{MaterialScreening}(x, ...)

\method{as_generic_mat_id}{MaterialIdentifierGeneric}(x, ...)

\method{as_generic_mat_id}{MaterialIdentifierScreening}(x, ...)

extract_well_material(x, row, col)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{A (vector of) \code{MaterialIdentifier} object(s).}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}

\item{material_type}{A \code{MaterialTypeIdentifierScreening} object to restrict
the material listing to a certain type of materials.}

\item{code}{The material code for which an id object is created.}

\item{type}{The material type (possible values depend on mode).}

\item{mode}{Switch between generic and screening material id objects.}

\item{types}{Select one or several material types for which to return the
type id objects. NULL returns all available.}

\item{row}{Either a single integer or a single character specifying a
plate row.}

\item{col}{A single integer specifying a plate column.}
}
\value{
Depending on the number of resulting objects, either a
\code{\link{json_class}} (single object) or a \code{\link{json_vec}} (multiple objects), is
returned. For \code{list_material()} and \code{extract_well_material()}, the
additional class attribute \code{MaterialGeneric} is added and the utility
functions \code{as_generic_mat_id()} and \code{as_screening_mat_id()} return
\code{MaterialIdentifierGeneric} and \code{MaterialIdentifierScreening} objects,
respectively, while \code{material_id()} can return either, depending on the
\code{mode} argument. Finally, \code{list_material_types()} returns
\code{MaterialTypeIdentifierGeneric} or \code{MaterialTypeIdentifierScreening}, again
depending on the \code{mode} argument.
}
\description{
Materials in openBIS can represent a variety of objects. For the InfectX
HTS setup, this is mainly limited to either compounds such as oligos or
small molecule drugs and targeted genes. Three different objects are used to
identify a material: \code{MaterialGeneric}, \code{MaterialIdentifierGeneric} and
\code{MaterialIdentifierScreening}. Converting to id object types can be
achieved with \code{as_generic_mat_id()} and \code{as_screening_mat_id()} while
listing materials as \code{MaterialGeneric} objects is possible with
\code{list_material()}.
}
\details{
Unfortunately in this version of the openBIS JSON-RPC API, there is no
possibility for listing all available materials for a project or
experiment. Methods that return \code{MaterialGeneric} objects include
\code{list_material()}, which can be dispatched on material id objects and
objects representing plates, and \code{\link[=search_openbis]{search_openbis()}} with the target selector
\code{target_object} set to \code{material}. Coercing \code{MaterialGeneric} objects to
material id objects is possible with \code{as_generic_mat_id()} and
\code{as_screening_mat_id()} which do not incur an API call.

Instantiating material id objects is either done manually by calling
\code{material_id()} or by querying openBIS for \code{MaterialGeneric} objects and
converting to \code{MaterialIdentifierGeneric} or \code{MaterialIdentifierScreening}.
A material id object is defined by a material code and a material type.
Available types depend on whether generic or screening material objects are
of interest. For generic material objects, possible ids are
\itemize{
\item compound
\item control
\item esirna
\item gene
\item mirna
\item mirna_inhibitor
\item mirna_mimic
\item pooled_sirna
\item sirna
}

and for screening materials, ids can be
\itemize{
\item compound
\item gene
\item oligo
}

Material type objects can be instantiated by calling
\code{list_material_types()}, where the \code{mode} argument acts as a switch to
choose between generic and screening objects. If only a subset of types
are relevant, the output of \code{list_material_types()} can be limited by
passing a character vector with type names as \code{types} argument. The second
piece of information for constructing material id objects, material codes,
depends on material type. Genes, for example are identified with Entrez
gene ids (e.g. 2475 for MTOR), while for compounds, a manufacturer name is
used (e.g. for Ambion and MTOR, AMBION_S602, AMBION_S603 and AMBION_S604).

Whenever \code{list_material()} is dispatched on a (set of) material id
object(s), a (set of) \code{MaterialGeneric} object(s) is returned. However if
the dispatch occurs on plate objects (\code{Plate}, \code{PlateIdentifier} or
\code{PlateMetadata}), a (set of) \code{PlateWellMaterialMapping} objects is returned.
If \code{material_type} is not specified (i.e. \code{NULL}), the \code{mapping} field in
the returned object will contain \code{NULL} for each well. When passing a set
of \code{MaterialTypeIdentifierScreening} objects, as returned by
\code{list_material_types()}, the \code{mapping} fields will contain material type
information where available. The convenience function
\code{extract_well_material()} can be applied to a \code{PlateWellMaterialMapping}
object and will return the selected \code{MaterialIdentifierScreening} object.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "getMaterialByCodes")}
\item \Sexpr[results=rd]{infx::docs_link("sas", "listPlateMaterialMapping")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  # search for a sample object corresponding to plate KB2-03-1I
  samp <- search_openbis(tok,
                         search_criteria(
                           attribute_clause("code",
                                            "/INFECTX_PUBLISHED/KB2-03-1I")
                         ),
                         target_object = "sample")

  # list all material types
  types <- list_material_types()
  print(types)

  # list all gene targets on plate KB2-03-1I
  mat_map <- list_material(tok, samp, types[[2L]])
  print(mat_map, depth = 5, length = 20L)
  # there are maximally width x height entries arranged in a linear,
  # row-major fashion; missing entries are omitted, but original indices
  # are accessible as original_index attributes
  length(mat_map[["mapping"]])
  attr(mat_map[["mapping"]][[42L]], "original_index")
  # well A24 does not have a gene target as it is a MOCK control well
  extract_well_material(mat_map, "A", 24)
  # well A22 however has a gene target
  a_22 <- extract_well_material(mat_map, "A", 22)
  print(a_22, depth = 2L)

  # search for a material with material code 3480
  igf1r <- search_openbis(tok,
                          search_criteria(attribute_clause("code", 3480)),
                          target_object = "material")

  all.equal(as_screening_mat_id(igf1r), a_22, check.attributes = FALSE)
  identical(igf1r, list_material(tok, a_22))
  identical(igf1r, 
            search_openbis(tok,
                           search_criteria(
                             property_clause("gene_symbol", "IGF1R")
                           ),
                           target_object = "material"))

  # search for an experiment object corresponding to plate KB2-03-1I
  exp <- search_openbis(tok,
                        search_criteria(
                          attribute_clause(
                            "code",
                            samp[["experimentIdentifierOrNull"]]
                          )
                        ),
                        target_object = "experiment")

  # list all wells for the current material within the selected experiment
  wells <- list_wells(tok, a_22, experiment = exp)
  # this yields 3 plates, one of which is KB2-03-1I
  get_field(get_field(wells, "experimentPlateIdentifier"), "plateCode")
  # and the material of interest is in well A22 in each one
  unique(get_field(wells, "wellPosition"))

  logout_openbis(tok)
}

}
\seealso{
Other object listing functions: \code{\link{list_datasets}},
  \code{\link{list_experiments}},
  \code{\link{list_plates}}, \code{\link{list_projects}},
  \code{\link{list_samples}}
}
\concept{object listing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/json-base.R
\name{print.json_class}
\alias{print.json_class}
\alias{[.json_class}
\alias{c.json_class}
\alias{rep.json_class}
\alias{print.json_vec}
\alias{[.json_vec}
\alias{[<-.json_vec}
\alias{[[.json_vec}
\alias{[[<-.json_vec}
\alias{c.json_vec}
\alias{rep.json_vec}
\alias{as_list}
\title{Base generics for JSON objects}
\usage{
\method{print}{json_class}(x, depth = 1L, width = getOption("width"),
  length = 100L, fancy = TRUE, ...)

\method{[}{json_class}(x, i, ...)

\method{c}{json_class}(x, ...)

\method{rep}{json_class}(x, ...)

\method{print}{json_vec}(x, depth = 1L, width = getOption("width"),
  length = 100L, fancy = TRUE, ...)

\method{[}{json_vec}(x, i, ...)

\method{[}{json_vec}(x, i, ...) <- value

\method{[[}{json_vec}(x, i, ...)

\method{[[}{json_vec}(x, i, ...) <- value

\method{c}{json_vec}(x, ...)

\method{rep}{json_vec}(x, ...)

as_list(x, ...)
}
\arguments{
\item{x}{Object to print/combine/subset, etc.}

\item{depth}{The maximum recursion depth for printing.}

\item{width}{Number of columns to maximally print.}

\item{length}{Number of lines to maximally print.}

\item{fancy}{Logical switch to enable font styles, colors and UTF box
characters for printing.}

\item{...}{Generic compatibility.}

\item{i}{Index for sub-setting. See \code{\link[base:[]{base::[()}} and
\code{\link[base:[[]{base::[[()}}.}

\item{value}{New values for replacement. See
\code{\link[base:[<-]{base::[<-()}} and
\code{\link[base:[[<-]{base::[[<-()}}.}
}
\value{
Depending on whether a single or a set of multiple objects is
represented, the S3 classes \code{\link{json_class}} or \code{\link{json_vec}} are applied
respectively.
}
\description{
Available base generic functions for objects that inherit \code{json_class} are
\code{\link[base:print]{base::print()}}, \code{\link[base:[]{base::[()}}, \code{\link[base:c]{base::c()}} and
\code{\link[base:rep]{base::rep()}} and for \code{json_vec} objects, all of the above in addition
to \code{\link[base:[<-]{base::[<-()}} and
\code{\link[base:[[<-]{base::[[<-()}} are implemented. For further
information on how these class-specific functions differ from their base
counterparts, refer to the details section.
}
\details{
Single bracket sub-setting of \code{json_class} objects preserves
class information, such that the resulting object has the same type but only
a subset of fields. Double bracket sub-setting of \code{json_class} objects
removes the enclosing type, as would be expected considering the list nature
of \code{json_class} objects. Combining or repeating \code{json_class} objects yields
\code{json_vec} objects with the same sub-type. Additionally, when combining
\code{json_class} objects using \code{\link[base:c]{base::c()}}, only \code{json_class} objects with the
same subtype as the first argument are allowed as further arguments.

Analogously to sub-setting of \code{json_class} objects, sub-setting a \code{json_vec}
object with \code{\link[base:[]{base::[()}} returns a \code{json_vec} object
with the same sub type as the one used as input, whereas sub-setting a
\code{json_vec} object with \code{\link[base:[[]{base::[[()}} yields the
selected \code{json_class} object. Replacement operators
\code{\link[base:[<-]{base::[<-()}} and
\code{\link[base:[[<-]{base::[[<-()}} mainly ensure that the objects
being inserted are of the correct sub-type, guaranteeing that all
\code{json_class} members of a given \code{json_vec} object are of the same sub-type.
Combining \code{json_vec} objects with \code{\link[base:c]{base::c()}} is possible whenever the
object passed as first argument has the same sub-type as the objects passed
as further arguments, which additionally are required to be \code{json_vec}
objects as well. Repeating a \code{json_vec} object using \code{\link[base:rep]{base::rep()}}, results
in a \code{json_vec} object of the same sub-type.

Printing of both \code{json_class} and \code{json_vec} objects is inspired by the
\code{ast} printing function of Hadley's
\href{https://git.io/vFMA5}{lobstr package} and borrows code from there.
Printing style can either be fancy (colors, UTF box characters. etc.) or
simple (controlled by the \code{fancy} flag) and several options are available
for setting the max printing width/length, as well as a max recursion depth
for nested \code{json_class} objects.
}
\examples{
obj_c <- json_class(a = json_class(b = "c", class = "foo"),
                    d = json_class(e = "f", class = "bar"),
                    class = "foobar")
obj_c
print(obj_c, depth = 2L)
print(obj_c, depth = 2L, length = 4L)
print(obj_c, depth = 2L, fancy = FALSE)

# sub-setting with single brackets preserves class information
obj_c["a"]
# whereas double brackets extract the selected element
obj_c[["a"]]

# vectors of json_class objects are json_vec objects
obj_cc <- rep(obj_c, 2)
identical(obj_cc, c(obj_c, obj_c))

print(obj_cc, depth = 2L, length = 8L)

obj_g <- json_class(a = json_class(b = "g", class = "foo"),
                    d = json_class(e = "h", class = "bar"),
                    class = "foobar")
obj_cg <- c(obj_c, obj_g)

# sub-setting json_vec objects with single brackets yields json_vec objects
class(obj_cg[1L])
# and with double brackets, the selected json_class object is extracted
class(obj_cg[[1L]])
identical(obj_cg[1L], json_vec(obj_cg[[1L]]))

# json_vec objects can also be combined using c
obj_i <- json_class(a = json_class(b = "i", class = "foo"),
                    d = json_class(e = "j", class = "bar"),
                    class = "foobar")

obj_cgi <- c(obj_cg, json_vec(obj_i))
length(obj_cgi)

# and repeated using rep
length(rep(obj_cgi, 2))

# additionally replacement operators are available
obj_cg[[2L]] <- obj_i
obj_cgi[1L:2L] <- obj_cg
identical(obj_cgi, c(obj_c, obj_i, obj_i))

}
\seealso{
Other json object handling functions: \code{\link{has_fields.json_class}},
  \code{\link{json_class}}, \code{\link{json_vec}}
}
\concept{json object handling functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/login.R
\name{login_openbis}
\alias{login_openbis}
\alias{logout_openbis}
\alias{is_token_valid}
\title{Create and destroy a login token}
\usage{
login_openbis(user = "rdgr2014", pwd = "IXPubReview",
  host_url = "https://infectx.biozentrum.unibas.ch",
  auto_disconnect = TRUE, ...)

logout_openbis(token, ...)

is_token_valid(token, ...)
}
\arguments{
\item{user, pwd}{Login credentials for an openBIS instance.}

\item{host_url}{Host url.}

\item{auto_disconnect}{Logical switch for automatically closing the
connection upon garbage collection of the token.}

\item{...}{Further arguments are forwarded to \code{\link[=make_request]{make_request()}}.}

\item{token}{Login token as created by \code{login_openbis()}.}
}
\value{
A login token as returned by \code{login_openbis()} is a character vector
and \code{logout_openbis()} returns \code{NULL} invisibly while \code{is_token_valid()}
returns a logical flag.
}
\description{
Login tokens for openBIS API calls can be created using \code{login_openbis()}.
If the \code{auto_disconnect} option is enabled, the user is automatically
logged out using \code{logout_openbis()} upon garbage collection of the token.
The validity of a token can be checked using \code{is_token_valid()} and a login
token can be manually destroyed by calling \code{logout_openbis()} on the token.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis",
"tryToAuthenticateForAllServices")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "logout")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "isSessionActive")}
}
}

\examples{
\donttest{
  # create a login token
  tok <- login_openbis()
  # check token
  is_token_valid(tok)

  # destroy token
  logout_openbis(tok)
  # token is no longer valid
  is_token_valid(tok)
}

}
\seealso{
Other utility functions: \code{\link{api_url}},
  \code{\link{make_requests}}
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file.R
\name{list_files}
\alias{list_files}
\alias{list_files.character}
\alias{list_files.DataSet}
\alias{list_files.DatasetIdentifier}
\alias{list_files.DatasetReference}
\alias{list_files.FeatureVectorDatasetReference}
\alias{list_files.FeatureVectorDatasetWellReference}
\alias{list_files.ImageDatasetReference}
\alias{list_files.MicroscopyImageReference}
\alias{list_files.PlateImageReference}
\alias{list_files.DataSetFileDTO}
\alias{fetch_files}
\alias{fetch_files.character}
\alias{fetch_files.NULL}
\alias{fetch_files.DataSet}
\alias{fetch_files.DatasetIdentifier}
\alias{fetch_files.DatasetReference}
\alias{fetch_files.FeatureVectorDatasetReference}
\alias{fetch_files.FeatureVectorDatasetWellReference}
\alias{fetch_files.ImageDatasetReference}
\alias{fetch_files.MicroscopyImageReference}
\alias{fetch_files.PlateImageReference}
\alias{fetch_files.DataSetFileDTO}
\alias{fetch_files.FileInfoDssDTO}
\alias{read_mat_files}
\title{List and download files}
\usage{
list_files(token, x, ...)

\method{list_files}{character}(token, x, path = "", recursive = TRUE,
  ...)

\method{list_files}{DataSet}(token, x, path = "", recursive = TRUE,
  ...)

\method{list_files}{DatasetIdentifier}(token, x, path = "",
  recursive = TRUE, ...)

\method{list_files}{DatasetReference}(token, x, path = "",
  recursive = TRUE, ...)

\method{list_files}{FeatureVectorDatasetReference}(token, x, path = "",
  recursive = TRUE, ...)

\method{list_files}{FeatureVectorDatasetWellReference}(token, x,
  path = "", recursive = TRUE, ...)

\method{list_files}{ImageDatasetReference}(token, x, path = "",
  recursive = TRUE, ...)

\method{list_files}{MicroscopyImageReference}(token, x, path = "",
  recursive = TRUE, ...)

\method{list_files}{PlateImageReference}(token, x, path = "",
  recursive = TRUE, ...)

\method{list_files}{DataSetFileDTO}(token, x, ...)

fetch_files(token, x, ...)

\method{fetch_files}{character}(token, x, files = NULL, n_con = 5L,
  reader = identity, ...)

\method{fetch_files}{NULL}(token, x, files, n_con = 5L,
  reader = identity, ...)

\method{fetch_files}{DataSet}(token, x, ...)

\method{fetch_files}{DatasetIdentifier}(token, x, ...)

\method{fetch_files}{DatasetReference}(token, x, ...)

\method{fetch_files}{FeatureVectorDatasetReference}(token, x, ...)

\method{fetch_files}{FeatureVectorDatasetWellReference}(token, x, ...)

\method{fetch_files}{ImageDatasetReference}(token, x, ...)

\method{fetch_files}{MicroscopyImageReference}(token, x, ...)

\method{fetch_files}{PlateImageReference}(token, x, ...)

\method{fetch_files}{DataSetFileDTO}(token, x, ...)

\method{fetch_files}{FileInfoDssDTO}(token, x, data_sets = NULL, ...)

read_mat_files(data)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to limit search for datasets/files with.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}} or \code{\link[=do_requests_serial]{do_requests_serial()}}/\code{\link[=do_requests_parallel]{do_requests_parallel()}}.}

\item{path}{A (vector of) file path(s) to be searched within a dataset.}

\item{recursive}{A (vector of) logicals, indicating whether to list files
recursively.}

\item{files}{Optional set of \code{FileInfoDssDTO} objects. If NULL, all files
corresponding to the specified datasets are assumed. This file list can be
filtered, by passing a regular expression as \code{file_regex} argument via
\code{...}.}

\item{n_con}{The number of simultaneous connections.}

\item{reader}{A function to read the downloaded data. Is forwarded as
finally argument to \code{\link[=do_requests_serial]{do_requests_serial()}}/\code{\link[=do_requests_parallel]{do_requests_parallel()}}.}

\item{data_sets}{Either a single dataset object (anything that has a
\code{dataset_code()} method) or a set of objects of the same length as \code{x}. If
\code{NULL} (default), each \code{FileInfoDssDTO} object passed as \code{x} is expected
to contain a \code{data_set} attribute.}

\item{data}{The data to be read.}
}
\value{
\code{list_files()} either returns a \code{\link{json_class}} or a \code{\link{json_vec}}
object of subtype \code{FileInfoDssDTO}, depending on whether a single or a set
of objects is retrieved. For \code{fetch_files()}, the return type depends on the
callback function passed as \code{reader} argument. At default, a \code{list} is
returned with an entry per file, holding a \code{raw} vector of the file data.
}
\description{
A dataset in openBIS represents a collection of files. The function
\code{list_files()} lists files associated with one or more datasets by
returning a set of \code{FileInfoDssDTO} objects. As this object type does not
contain information on data set association, the data set code is saved
as \code{data_set} attribute with each \code{FileInfoDssDTO} object. Data set files
can be fetched using \code{fetch_files()}, which can either retrieve all
associated files or use file path information, for example from
\code{FileInfoDssDTO} objects to only download a subset of files.
}
\details{
Data sets for \code{list_files()} can be specified as character vector of
dataset codes and therefore all objects for which the internal method
\code{\link[=dataset_code]{dataset_code()}} exists can be used to select datasets. This includes data
set and data set id objects as well as the various flavors of data set
reference objects. In addition to these dataset-representing objects,
dispatch on \code{DataSetFileDTO} objects is possible as well.

File listing can be limited to a certain path within the dataset and the
search can be carried out recursively or non-recursively. In case a set of
objects is passed, the search-tuning arguments \code{path} and \code{recursive} have
to be either of length 1 or of the same length as \code{x}. If dispatch occurs
on \code{DataSetFileDTO} objects, the \code{path} and \code{recursive} arguments are not
needed, as this information is already encoded in the objects passed as \code{x}.
A separate API call is necessary for each of the objects the dispatch
occurs on.

The function \code{fetch_files()} downloads files associated with a dataset.
In order to identify a file, both a data set code and a file path, relative
to the data set root, are required. \code{fetch_files()} can be called in a
variety of ways and internally uses a double dispatch mechanism, first
resolving the data set codes and then calling the non-exported function
\code{fetch_ds_files()} which dispatches on file path objects.

Data set code information can either be communicated using any of the
objects understood by \code{\link[=dataset_code]{dataset_code()}} (including data set, data set id and
data set reference objects) or directly as a character vector, passed as
\code{x} argument. In case data set code information is omitted (passed to \code{x}
as \code{NULL}), the objects encoding file paths have to specify the
corresponding data sets. Furthermore, \code{DataSetFileDTO} objects may be
passed as \code{x} argument to \code{fetch_files()}, which will internally call
\code{fetch_files()} again, setting the argument \code{x} to \code{NULL} and pass the
\code{DataSetFileDTO} objects as files argument. Finally, if \code{FileInfoDssDTO}
are passed to \code{fetch_files()} as \code{x} argument, an optional argument
\code{data_sets} may be specified (it defaults to \code{NULL}) and as above,
\code{fetch_files()} is called again with these two arguments rearranged.

The internal generic function \code{fetch_ds_files()} can be dispatched on
several objects again. When no files are specified (\code{NULL} is passed as
\code{files} argument to \code{fetch_files()}), all available files for the given
data sets are queried. This list can be filtered using the \code{file_regex()}
argument which can be a single regular expression and is applied to file
paths. File paths can be specified as character vector, \code{FileInfoDssDTO} or
\code{DataSetFileDTO} objects. If dispatch occurs on \code{FileInfoDssDTO}, and no
data set code information is available (\code{NULL} passed as \code{x} or \code{data_sets}
argument to \code{fetch_files()}) each \code{FileInfoDssDTO} must contain a \code{data_set}
attribute. Additionally, downloaded files are checked for completeness, as
these objects contain file sizes. If dispatch occurs on \code{DataSetFileDTO}
objects or a character vector, this sanity check is not possible.

Files can only be retrieved after previously having created a corresponding
download url using \code{\link[=list_download_urls]{list_download_urls()}}, as file urls in openBIS have a
limited lifetime and therefore must be used shortly after being created. A
list of \code{call} objects (see \code{\link[base:call]{base::call()}}) is created and passed to either
\code{\link[=do_requests_serial]{do_requests_serial()}} or \code{\link[=do_requests_parallel]{do_requests_parallel()}}. Whether file fetching
is carried out in serial or parallel is controlled by the \code{n_con} argument.
In case a download fails, it is retried again up to the number of times
specified as \code{n_try}. Finally, a function with a single argument can be
passed as the argument \code{done}, which takes the downloaded data as input and
does some processing.

A function for reading the binary data retrieved from openBIS can be
supplied to \code{fetch_files()} as \code{reader} argument. Single cell feature files
as produced by CellProfiler, are stored as Matlab v5.0 \code{.mat} files and
the function \code{read_mat_files()} reads such files using \code{\link[R.matlab:readMat]{R.matlab::readMat()}}
and checks for certain expected attributes and simplifies the read
structure.

The list returned by \code{read_mat_files()} is arranged such that each node
corresponds to a single image and contains a list which is either holding a
single value or a vector of values. For a plate with 16 rows, 24 columns
and 3 x 3 imaging sites this will yield a list of length 3456. Index
linearization is in row-major fashion for both wells and sites.
Furthermore, imaging sites come first such that in this example, the first
three list entries correspond to image row 1 (left to right) of well A1,
the next three entries correspond to row 2 of well A1, images 10 through 12
correspond to row 1 of well A2, etc. Well A2 is located in row 1, column 2
of a plate.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrg", "listFilesForDataSet")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  # search for a cell profiler feature data set from plate KB2-03-1I
  search <- search_criteria(
    attribute_clause("type", "HCS_ANALYSIS_CELL_FEATURES_CC_MAT"),
    sub_criteria = search_sub_criteria(
      search_criteria(attribute_clause("code",
                                       "/INFECTX_PUBLISHED/KB2-03-1I")),
      type = "sample"
    )
  )
  ds <- search_openbis(tok, search)

  # list all files of this data set
  all_files <- list_files(tok, ds)
  length(all_files)

  # select some of the files, e.g. all count features per image
  some_files <- all_files[grepl("Image\\\\.Count_",
                                get_field(all_files, "pathInDataSet"))]
  length(some_files)

  # download the selected files
  data <- fetch_files(tok, some_files)

  # the same can be achieved by passing a file_regex argument to
  # fetch_files(), which internally calls list_files() and filters files
  identical(data, fetch_files(tok, ds, file_regex = "Image\\\\.Count_"))

  # all returned data is raw, the reader argument can be used to supply
  # a function that processes the downloaded data
  sapply(data, class)
  data <- fetch_files(tok, some_files, reader = read_mat_files)
  sapply(data, class)

  logout_openbis(tok)
}

}
\seealso{
Other resource listing/downloading functions: \code{\link{fetch_images}},
  \code{\link{list_download_urls}},
  \code{\link{list_features}}
}
\concept{resource listing/downloading functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.R
\name{feat_ds_well_ref}
\alias{feat_ds_well_ref}
\title{Create well feature dataset reference objects}
\usage{
feat_ds_well_ref(x, wells)
}
\arguments{
\item{x}{A set of \code{FeatureVectorDatasetReference} objects.}

\item{wells}{A set of \code{WellPosition} objects.}
}
\description{
While \code{FeatureVectorDatasetReference} objects represent feature datasets on
a plate level, \code{FeatureVectorDatasetWellReference} reference feature
datasets on well-level. In order to create such well-level representations,
a set of \code{FeatureVectorDatasetReference} and a set of \code{WellPosition} are
combined to \code{FeatureVectorDatasetWellReference} objects where each instance
contains a single object of the inputted sets. If either argument is of
length greater than 1, the other argument has to be of the same length or of
length 1, in which case it will be \code{\link[base:rep]{base::rep()}}eated to the required
length.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.R
\name{list_samples}
\alias{list_samples}
\alias{list_samples.ExperimentIdentifier}
\alias{list_samples.Experiment}
\alias{list_samples.Plate}
\alias{list_samples.PlateIdentifier}
\alias{list_samples.PlateMetadata}
\alias{list_samples.WellIdentifier}
\alias{list_samples.WellMetadata}
\alias{list_sample_types}
\title{List samples and sample types}
\usage{
list_samples(token, x, ...)

\method{list_samples}{ExperimentIdentifier}(token, x, ...)

\method{list_samples}{Experiment}(token, x, ...)

\method{list_samples}{Plate}(token, x, ...)

\method{list_samples}{PlateIdentifier}(token, x, ...)

\method{list_samples}{PlateMetadata}(token, x, ...)

\method{list_samples}{WellIdentifier}(token, x, ...)

\method{list_samples}{WellMetadata}(token, x, ...)

list_sample_types(token, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object to limit the number of returned samples, e.g. a set of
\code{ExperimentIdentifier} or \code{Experiment} objects.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}
}
\value{
Depending on the number of resulting objects, either a
\code{\link{json_class}} (single object) or a \code{\link{json_vec}} (multiple objects), is
returned. Sample objects as returned by \code{list_samples()} additionally
inherit from the \code{Sample} class and sample type objects returned by
\code{list_sample_types()} inherit from \code{SampleType}.
}
\description{
In openBIS, samples can be seen as a generalization of plates and wells (
see \code{\link[=list_plates]{list_plates()}} and \code{\link[=list_wells]{list_wells()}}). In fact, for the HTS focused
openBIS instance of InfectX, wells and plates represent the only types of
samples available. Samples can either be retrieved using \code{\link[=search_openbis]{search_openbis()}}
and setting the argument \code{target_object} to \code{sample} or listed by calling
\code{list_samples()}. Furthermore, all available sample types can be listed
using \code{list_sample_types()}.
}
\details{
\code{list_samples()} can be dispatched on objects identifying experiments
(\code{Experiment} and \code{ExperimentIdentifier}), in which case the associated
plate samples are returned, on objects representing plates (\code{Plate},
\code{PlateIdentifier} and \code{PlateMetadata}) or on objects representing wells
(\code{WellIdentifier} and \code{WellMetadata}). For plates, the corresponding plate
samples and for wells, the corresponding well samples are returned. It is
therefore not possible to list all well samples for a plate. This could
however be achieved by listing all wells of a plate using \code{\link[=list_wells]{list_wells()}}
and calling \code{list_samples()} on the returned set of \code{WellIdentifier}
objects. A separate API call is required for each \code{json_class} object
contained in the \code{json_vec} passed to \code{list_samples()} as \code{x} argument.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listSamplesForExperiment")}
\item \Sexpr[results=rd]{infx::docs_link("sas", "getPlateSample")}
\item \Sexpr[results=rd]{infx::docs_link("sas", "getWellSample")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "listSampleTypes")}
}
}

\examples{
\donttest{
  tok <- login_openbis()

  # search for an experiment, e.g. ADENO-AU-K1
  exp <- search_openbis(tok,
                        search_criteria(
                          property_clause("pathogen", "Adenovirus"),
                          property_clause("library", "Ambion"),
                          property_clause("geneset", "Kinome"),
                          property_clause("replicate", 1L)
                        ),
                        target_object = "experiment")

  # list all plate samples associated with this experiment
  plate_samp <- list_samples(tok, exp)
  length(plate_samp)

  # the same plates are returned using list_plates(), not regarding order
  plates <- list_plates(tok, exp)
  plate_samp <- plate_samp[order(get_field(plate_samp, "code"))]

  identical(as_plate_id(plates), as_plate_id(plate_samp))

  # plates can be converted to samples and back to plates again
  identical(as_plate_id(plates), as_plate_id(list_samples(tok, plates)))

  # the same is not possible for wells: first well ids are listed for a
  # plate
  well_ids <- list_wells(tok, plates[[1L]])
  # for efficiency only the first 5 wells are retained
  well_ids <- well_ids[1:5]
  print(well_ids, length = 10L)
  # the corresponding well samples are fetched
  well_samp <- list_samples(tok, well_ids)
  # from well sample objects it is not possible to directly create well id
  # objects as no plate id information is available
  print(well_samp, length = 20L)
  # with a bit of manual work however it is possible to create well id
  # objects from well samples
  wells <- well_id(get_field(well_samp, "permId"), plates[[1L]],
                   well_code = get_field(well_samp, "code"))
  identical(wells, well_ids)

  logout_openbis(tok)
}

}
\seealso{
Other object listing functions: \code{\link{list_datasets}},
  \code{\link{list_experiments}},
  \code{\link{list_material}}, \code{\link{list_plates}},
  \code{\link{list_projects}}
}
\concept{object listing functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/url.R
\name{list_download_urls}
\alias{list_download_urls}
\alias{list_download_urls.character}
\alias{list_download_urls.DataSet}
\alias{list_download_urls.DatasetIdentifier}
\alias{list_download_urls.DatasetReference}
\alias{list_download_urls.FeatureVectorDatasetReference}
\alias{list_download_urls.FeatureVectorDatasetWellReference}
\alias{list_download_urls.ImageDatasetReference}
\alias{list_download_urls.MicroscopyImageReference}
\alias{list_download_urls.PlateImageReference}
\alias{list_download_urls.DataSetFileDTO}
\alias{list_datastores}
\alias{list_datastore_urls}
\alias{list_datastore_urls.NULL}
\alias{list_datastore_urls.character}
\alias{list_datastore_urls.DataSet}
\alias{list_datastore_urls.DatasetIdentifier}
\alias{list_datastore_urls.DatasetReference}
\alias{list_datastore_urls.FeatureVectorDatasetReference}
\alias{list_datastore_urls.FeatureVectorDatasetWellReference}
\alias{list_datastore_urls.ImageDatasetReference}
\alias{list_datastore_urls.MicroscopyImageReference}
\alias{list_datastore_urls.PlateImageReference}
\title{List data store servers and urls}
\usage{
list_download_urls(token, x, ...)

\method{list_download_urls}{character}(token, x, path, timeout = NA, ...)

\method{list_download_urls}{DataSet}(token, x, path, timeout = NA, ...)

\method{list_download_urls}{DatasetIdentifier}(token, x, path,
  timeout = NA, ...)

\method{list_download_urls}{DatasetReference}(token, x, path,
  timeout = NA, ...)

\method{list_download_urls}{FeatureVectorDatasetReference}(token, x, path,
  timeout = NA, ...)

\method{list_download_urls}{FeatureVectorDatasetWellReference}(token, x,
  path, timeout = NA, ...)

\method{list_download_urls}{ImageDatasetReference}(token, x, path,
  timeout = NA, ...)

\method{list_download_urls}{MicroscopyImageReference}(token, x, path,
  timeout = NA, ...)

\method{list_download_urls}{PlateImageReference}(token, x, path,
  timeout = NA, ...)

\method{list_download_urls}{DataSetFileDTO}(token, x, timeout = NA, ...)

list_datastores(token, ...)

list_datastore_urls(token, x = NULL, ...)

\method{list_datastore_urls}{NULL}(token, x, ...)

\method{list_datastore_urls}{character}(token, x, ...)

\method{list_datastore_urls}{DataSet}(token, x, ...)

\method{list_datastore_urls}{DatasetIdentifier}(token, x, ...)

\method{list_datastore_urls}{DatasetReference}(token, x, ...)

\method{list_datastore_urls}{FeatureVectorDatasetReference}(token, x, ...)

\method{list_datastore_urls}{FeatureVectorDatasetWellReference}(token, x,
  ...)

\method{list_datastore_urls}{ImageDatasetReference}(token, x, ...)

\method{list_datastore_urls}{MicroscopyImageReference}(token, x, ...)

\method{list_datastore_urls}{PlateImageReference}(token, x, ...)
}
\arguments{
\item{token}{Login token as created by \code{login_openbis()}.}

\item{x}{Object representing a (set of) dataset(s), e.g. a vector of dataset
codes, or a set of \code{DataSet}s or \code{DatasetIdentifier}s.}

\item{...}{Generic compatibility. Extra arguments will be passed to
\code{\link[=make_requests]{make_requests()}}.}

\item{path}{A character vector of file paths within datasets.}

\item{timeout}{Time-span (in seconds) for which the file download link
should be valid.}
}
\value{
Both \code{list_download_urls()} and \code{list_datastore_urls()} return
character vectors while \code{list_datastores()} returns either a \code{\link{json_class}}
(single object) or a \code{\link{json_vec}} (multiple objects), dependent on the
number of resulting objects, with sub-type \code{DataStore}.
}
\description{
In order to download files from openBIS, download urls have to be generated
first, which can be done by calling \code{list_download_urls()}. This function
is used in \code{\link[=fetch_files]{fetch_files()}}, which iterates over the selected files, creating
download links and executing the downloads. All data store servers
registered to an openBIS instance are listed by \code{list_datastores()} and data
store server urls per data set can be queried by calling
\code{list_datastore_urls()}.
}
\details{
To specify files for which links are requested by \code{list_download_urls()},
both a data set code and a file path are required. Objects, apart from
character vectors of data set codes, that may be passed to identify the
data set therefore include \code{DataSet}, \code{DatasetIdentifier},
\code{DatasetReference}, \code{FeatureVectorDatasetReference},
\code{FeatureVectorDatasetWellReference}, \code{ImageDatasetReference},
\code{MicroscopyImageReference} and \code{PlateImageReference}. Additionally, dispatch
of \code{list_download_urls()} is possible on \code{DataSetFileDTO} objects which
contain both information on data set and file path of a file. A \code{timeout}
argument may be specified, determining how long (in seconds) the generated
url is valid for. If no specific timeout value is passed the url is valid
for what the openBIS documentation calls "a short time".

\code{list_datastore_urls()} as \code{list_download_urls()} ultimately requires a
character vector of data set codes to make the API call and therefore
dispatch is possible on, in addition to character vector, \code{DataSet},
\code{DatasetIdentifier}, \code{DatasetReference}, \code{FeatureVectorDatasetReference},
\code{FeatureVectorDatasetWellReference}, \code{ImageDatasetReference},
\code{MicroscopyImageReference} and \code{PlateImageReference} objects. Dispatch on
\code{NULL} requests the default data store server url. Datastore sever url
related functionality is uninteresting for the InfectX set-up, as only a
single data store server exists, the url of which can be retrieved by a call
to \code{list_datastores()}.
}
\section{openBIS}{

\itemize{
\item \Sexpr[results=rd]{infx::docs_link("dsrg",
"getDownloadUrlForFileForDataSet")}
\item \Sexpr[results=rd]{infx::docs_link("dsrg",
"getDownloadUrlForFileForDataSetWithTimeout")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis", "listDataStores")}
}


\itemize{
\item \Sexpr[results=rd]{infx::docs_link("gis",
"getDefaultPutDataStoreBaseURL")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "tryGetDataStoreBaseURL")}
\item \Sexpr[results=rd]{infx::docs_link("gis", "getDataStoreBaseURLs")}
}
}

\examples{
\donttest{
  tok <- login_openbis()
  
  # data store server information
  list_datastores(tok)

  # search for a cell profiler feature data set from plate KB2-03-1I
  search <- search_criteria(
    attribute_clause("type", "HCS_ANALYSIS_CELL_FEATURES_CC_MAT"),
    sub_criteria = search_sub_criteria(
      search_criteria(attribute_clause("code",
                                       "/INFECTX_PUBLISHED/KB2-03-1I")),
      type = "sample"
    )
  )
  ds <- search_openbis(tok, search)

  # list all files of this data set
  files <- list_files(tok, ds)
  # extract file paths
  file_paths <- get_field(files, "pathInDataSet")
  # select a file
  file_path <- file_paths[grepl("Count_Cells", file_paths)]

  # generate url
  list_download_urls(tok, ds, file_path)

  # generate url and download file
  dat <- read_mat_files(url(list_download_urls(tok, ds, file_path)[[1L]]))
  attributes(dat)
  str(as.integer(dat))

  # set timeout to 2 sec
  file_url <- list_download_urls(tok, ds, file_path, timeout = 2L)
  tmp <- read_mat_files(url(file_url[[1L]]))

  # let timeout expire
  file_url <- list_download_urls(tok, ds, file_path, timeout = 2L)
  Sys.sleep(4L)
  tmp <- read_mat_files(url(file_url[[1L]]))

  logout_openbis(tok)
}
}
\seealso{
Other resource listing/downloading functions: \code{\link{fetch_images}},
  \code{\link{list_features}}, \code{\link{list_files}}
}
\concept{resource listing/downloading functions}
