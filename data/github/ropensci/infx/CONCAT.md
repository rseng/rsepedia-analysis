
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
