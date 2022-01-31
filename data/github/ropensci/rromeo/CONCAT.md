
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `rromeo` – an R interface for SHERPA/RoMEO API

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis build
status](https://travis-ci.org/ropensci/rromeo.svg?branch=master)](https://travis-ci.org/ropensci/rromeo)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ropensci/rromeo?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/rromeo)
[![codecov](https://codecov.io/gh/ropensci/rromeo/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rromeo)
[![cran
checks](https://cranchecks.info/badges/summary/rromeo)](https://cran.r-project.org/web/checks/check_results_rromeo.html)
[![CRAN-version](https://www.r-pkg.org/badges/version/rromeo)](https://cran.r-project.org/package=rromeo)
[![](https://badges.ropensci.org/285_status.svg)](https://github.com/ropensci/onboarding/issues/285)

`rromeo` is an R client for the [SHERPA/RoMEO
API](http://www.sherpa.ac.uk/romeo/index.php?la=en&fIDnum=&mode=simple).
SHERPA/RoMEO is a database that gives information on editorial policies
of scientific journals regarding the archival of preprint, postprint and
publishers’ manuscripts. `rromeo` is aimed at scientists interested in
archival practices of scientific journals, such as professionals of
[scientometrics](https://en.wikipedia.org/wiki/Scientometrics) but also
at scientist of specific fields interested in the practices of their
fields.

## Install

The latest stable release of `rromeo` is available on CRAN and can be
installed with:

``` r
install.packages("rromeo")
```

You can also install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/rromeo")
```

## API Key

Note that SHERPA/RoMEO lets you run 500 requests per day per IP address,
by [registering for a free API
key](http://www.sherpa.ac.uk/romeo/apiregistry.php) you can bypass this
limit.

`rromeo` can use your registered SHERPA/RoMEO API key; you can either
pass it as a string when querying the data with the argument `key`:

``` r
rr_journal_name("Journal of Geology", key = "Iq83AIL5bss")
```

or you can specify the environment variable `SHERPAROMEO_KEY` in an
`.Rprofile` or in an `.Renviron` file and `rromeo` will automatically
retrieve the API key. See the [specific
vignette](https://docs.ropensci.org/rromeo/articles/setting_up_api_key.html)
to know how to apply and use the API key with `rromeo`.

## Usage

`rromeo` contains functions to retrieve data from the SHERPA/RoMEO API
(for a complete overview please refer to the
[vignette](https://docs.ropensci.org/rromeo/articles/rromeo.html)). The
data is released under the [Creative Commons
Attribution-NonCommercial-ShareAlike 2.5 (CC BY-NC-SA 2.5)
license](https://creativecommons.org/licenses/by-nc-sa/2.5/). A
suggestion of citation is included in `rromeo` via `citation("rromeo")`.

`rromeo` functions are prefixed with `rr_` such as `rr_journal_name()`
that lets you retrieve a journal policy information using the title of a
journal:

``` r
rromeo::rr_journal_name("Journal of Biogeography", qtype = "exact")
#>                     title provided_issn      issn romeocolour preprint
#> 1 Journal of Biogeography          <NA> 0305-0270      yellow      can
#>    postprint    pdf pre_embargo post_embargo pdf_embargo
#> 1 restricted cannot        <NA>    12 months        <NA>
```

the `qtype` argument indicates the type of query to make (`exact` for
exact matching of the title, `contains` for partial matching and `starts
with` to match only the beginning of the title).

You can also retrieve a journal information using its ISSN:

``` r
rromeo::rr_journal_issn("0305-0270")
#>                     title provided_issn      issn romeocolour preprint
#> 1 Journal of Biogeography     0305-0270 0305-0270      yellow      can
#>    postprint    pdf pre_embargo post_embargo pdf_embargo
#> 1 restricted cannot        <NA>    12 months        <NA>
```

`rromeo` also provides a function to retrieve information based on
publisher ID `rr_publisher()`.

SHERPA/RoMEO provides a synthetic “colour” for each journal, the colour
summarizes the editorial policy of a journal:

| RoMEO colour | Archiving policy                                        |
| :----------- | :------------------------------------------------------ |
| `green`      | can archive preprint, postprint and publisher’s version |
| `blue`       | can archive postprint **or** publisher’s version        |
| `yellow`     | can archive preprint                                    |
| `white`      | archiving not formally supported                        |

(Table taken from
<http://www.sherpa.ac.uk/romeo/definitions.php#colours>)

`rromeo` lets you retrieve the policies of all journals of a given
colour using the function `rr_romeo_colour()` (**NOTE:** this function
can be slow as there are many journals to retrieve):

``` r
green_journals = rromeo::rr_romeo_colour("green")
green_journals[8:12,]
#>    romeoid                                                   publisher
#> 8     1128 Association for Information Science and Technology (ASIS&T)
#> 9     1937                                       University of Arizona
#> 10    2951                               Geological Society of America
#> 11    2521                              University of California Press
#> 12    2306                                  Optical Society of America
#>                  alias romeocolour preprint postprint        pdf
#> 8              JASIS&T       green      can       can     cannot
#> 9          Radiocarbon       green      can       can restricted
#> 10           GSA Today       green      can       can        can
#> 11            Collabra       green      can       can        can
#> 12 No Paid Open Access       green      can       can     cannot
```

## Dependency network (Imports only)

<img src="man/figures/README-dependency_network_imports-1.svg" width="100%" />

## Dependency network (Imports and Suggests)

<img src="man/figures/README-dependency_network_full-1.svg" width="100%" />

## Contributing to `rromeo`

We welcome contribution to `rromeo`\! Please read the [contribution
guidelines](https://docs.ropensci.org/rromeo/CONTRIBUTING.html) if you
want to contribute, as well as the below-mentioned Code of Conduct.

## Code of Conduct

Please note that the `rromeo` project is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/).
By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# rromeo 0.1.1

## NEW FEATURES

* add column for provided ISSNs in `rr_journal_*()` functions (fix #54)

## BUG FIXES

* Move to static vignettes to avoid issues with CRAN and others
* Update cassette following `vcr` 0.5.0 update

# rromeo 0.1.0

## NEW FEATURES

* released to CRAN
## Test environments
* local Windows 8.1, R-release
* Ubuntu 16.04 (on travis-ci), R-oldrel
* Ubuntu 16.04 (on travis-ci), R-release
* Ubuntu 16.04 (on travis-ci), R-devel
* Windows Server 2012 R2, R-release (AppVeyor CI)
* win-builder R-oldrel
* win-builder R-release
* win-builder R-devel
* Ubuntu Linux 16.04 LTS, R-release, GCC (rhub)
* Fedora Linux, R-devel, clang, gfortran (rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (rhub)
* Debian Linux, R-patched, GCC (rhub)

## R CMD check results

0 errors | 0 warnings | 0 note
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

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

# Contributing to `rromeo`

We welcome contributions to `rromeo`!
This outlines how to propose a change to rromeo.

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
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the rromeo project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By 
contributing to this project you agree to abide by its terms.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# `rromeo` – an R interface for SHERPA/RoMEO API

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![Travis build status](https://travis-ci.org/ropensci/rromeo.svg?branch=master)](https://travis-ci.org/ropensci/rromeo) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/rromeo?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/rromeo) [![codecov](https://codecov.io/gh/ropensci/rromeo/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rromeo) [![cran checks](https://cranchecks.info/badges/summary/rromeo)](https://cran.r-project.org/web/checks/check_results_rromeo.html) [![CRAN-version]( https://www.r-pkg.org/badges/version/rromeo)](https://cran.r-project.org/package=rromeo) [![](https://badges.ropensci.org/285_status.svg)](https://github.com/ropensci/onboarding/issues/285)

`rromeo` is an R client for the [SHERPA/RoMEO API](http://www.sherpa.ac.uk/romeo/index.php?la=en&fIDnum=&mode=simple). SHERPA/RoMEO is a database that gives information on editorial policies of scientific journals regarding the archival of preprint, postprint and publishers' manuscripts. `rromeo` is aimed at scientists interested in archival practices of scientific journals, such as professionals of [scientometrics](https://en.wikipedia.org/wiki/Scientometrics) but also at scientist of specific fields interested in the practices of their fields.

## Install

The latest stable release of `rromeo` is available on CRAN and can be installed
with:

```{r, eval = FALSE}
install.packages("rromeo")
```

You can also install the development version from GitHub:

```{r}
# install.packages("remotes")
remotes::install_github("ropensci/rromeo")
```

## API Key

Note that SHERPA/RoMEO lets you run 500 requests per day per IP address, by [registering for a free API key](http://www.sherpa.ac.uk/romeo/apiregistry.php) you can bypass this limit.

`rromeo` can use your registered SHERPA/RoMEO API key; you can either pass it as a string when querying the data with the argument `key`:

```{r, eval = FALSE}
rr_journal_name("Journal of Geology", key = "Iq83AIL5bss")
```

or you can specify the environment variable `SHERPAROMEO_KEY` in an `.Rprofile` or in an `.Renviron` file and `rromeo` will automatically retrieve the API key.
See the [specific vignette](https://docs.ropensci.org/rromeo/articles/setting_up_api_key.html) to know how to apply and use the API key with `rromeo`.

## Usage

`rromeo` contains functions to retrieve data from the SHERPA/RoMEO API (for a complete overview please refer to the [vignette](https://docs.ropensci.org/rromeo/articles/rromeo.html)). The data is released under the [Creative Commons Attribution-NonCommercial-ShareAlike 2.5 (CC BY-NC-SA 2.5) license](https://creativecommons.org/licenses/by-nc-sa/2.5/). A suggestion of citation is included in `rromeo` via `citation("rromeo")`.

`rromeo` functions are prefixed with `rr_` such as `rr_journal_name()` that lets you retrieve a journal policy information using the title of a journal:

```{r example_rr_journal_name}
rromeo::rr_journal_name("Journal of Biogeography", qtype = "exact")
```

the `qtype` argument indicates the type of query to make (`exact` for exact matching of the title, `contains` for partial matching and `starts with` to match only the beginning of the title).

You can also retrieve a journal information using its ISSN:

```{r example_rr_journal_issn}
rromeo::rr_journal_issn("0305-0270")
```

`rromeo` also provides a function to retrieve information based on publisher ID `rr_publisher()`.

SHERPA/RoMEO provides a synthetic "colour" for each journal, the colour summarizes the editorial policy of a journal:

  | RoMEO colour | Archiving policy                                        |
  |:-------------|:--------------------------------------------------------|
  | `green`      | can archive preprint, postprint and publisher's version |
  | `blue`       | can archive postprint **or** publisher's version        |
  | `yellow`     | can archive preprint                                    |
  | `white`      | archiving not formally supported                        |
  
(Table taken from <http://www.sherpa.ac.uk/romeo/definitions.php#colours>)

`rromeo` lets you retrieve the policies of all journals of a given colour using the function `rr_romeo_colour()` (**NOTE:** this function can be slow as there are many journals to retrieve):

```{r example_rr_romeo_colour}
green_journals = rromeo::rr_romeo_colour("green")
green_journals[8:12,]
```

## Dependency network (Imports only)

```{r dependency_network_imports, echo=FALSE, message=FALSE, warning=FALSE, dev='svglite'}
source("https://gist.githubusercontent.com/Rekyt/261c2c3040715dcd8c211d7f5d18c9c8/raw/")
plot_dependencies()
```

## Dependency network (Imports and Suggests)

```{r dependency_network_full, echo=FALSE, message=FALSE, warning=FALSE, dev='svglite'}
source("https://gist.githubusercontent.com/Rekyt/261c2c3040715dcd8c211d7f5d18c9c8/raw/")
plot_dependencies(type = "full")
```


## Contributing to `rromeo`

We welcome contribution to `rromeo`!
Please read the [contribution guidelines](https://docs.ropensci.org/rromeo/CONTRIBUTING.html) if you want to contribute, as well as the below-mentioned Code of Conduct.

## Code of Conduct

Please note that the `rromeo` project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/rromeo/CODE_OF_CONDUCT.html). By contributing
to this project, you agree to abide by its terms.

 [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "How to set up and use an API key?"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to set up and use an API key?}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




## What is an API key?

Like many other R packages that are API clients, `rromeo` lets you set up and use an API key.
An API key is a string that you communicate to the server (see [Wikipedia Page](https://en.wikipedia.org/wiki/Application_programming_interface_key)). The key is used by the server to identify you. This helps the server to manage your access, for example, giving access to specific services.


## Why use one?

Many APIs have a limit in the number of queries you can make in a given amount of time. Using an API key, because it identifies you, lets you go over this limit so the server knows it is not under a [Denial of Service (DoS) attack](https://en.wikipedia.org/wiki/Denial-of-service_attack).
SHERPA/RoMEO has a limit of 500 queries per day without an API key. Registering for an API key lets you exceed this limit.

## Registering for an API key for SHERPA/RoMEO

You can register a free API key at <http://www.sherpa.ac.uk/romeo/apiregistry.php>
You have to provide your name, your job title as well as a valid email address to register for an API key.

It is considered good practice to register one API key per project (e.g. software or research project), in order to avoid making too many queries to the API. Do fill out the description of your application to help the people who made SHERPA/RoMEO data freely available make statistics.


## How to set up your API key in `rromeo`

Now that you have registered for an API key, you can use it in `rromeo`. There are several ways to set up your API key in `rromeo`.


### Using `Sys.setenv()` during an interactive session

One way to set up your API key without using the `key` argument is to use **environment variables**. Environment variables are system or user-defined variables that are used to store information during your R session.

`rromeo` searches the environment variable `SHERPAROMEO_KEY` to check if an API key has been defined.
You can access it using the function `Sys.getenv()` with the name of the variable as as a string:

```
Sys.getenv("SHERPAROMEO_KEY")
```

If none has been defined it should return `""`. You can then set the environment variable using the function `Sys.setenv()`:


```r
Sys.setenv(SHERPAROMEO_KEY = "Iq83AIL5bss")
Sys.getenv("SHERPAROMEO_KEY")
#> [1] "Iq83AIL5bss"
```

It is set up for the rest of your R session and `rromeo` will automatically use it when you call the different functions. You can also use the `check_key()` function that should give the same result if your key is well set-up:


```r
library("rromeo")
check_key()
#> [1] "Iq83AIL5bss"
```


### Using the dedicated function `rr_auth()`

`rromeo` provides the function `rr_auth()` that creates a well named environmental variable that can be used in the rest of your session:


```r
rr_auth("Iq83AIL5bs2")
check_key()
#> [1] "Iq83AIL5bs2"
```

This works so that you don't have to specify the name of the environmental variable by hand. Under the hood `rr_auth()` uses the same mechanism as explained in the above mentioned section.


### Setting up your API key in an `.Rprofile` file

Every time R starts it looks for `.Rprofile` files in different locations:

- `R_HOME` the directory in which R is installed,
- `HOME` the user's home directory,
- R's current working directory.

R only loads one `.Rprofile` file per session and thus an `.Rprofile` file at the project-level overrides files in other locations.

The `.Rprofile` file is an R script that is launch each time R starts. Put it at the root of your project and type the following:

```
SHERPAROMEO_KEY = "Iq83AIL5bss"
```

You can then reload your session and check that `rromeo` managed to get your key by using the `check_key()` function:


```r
check_key()
#> [1] "Iq83AIL5bs2"
```

Now `rromeo` can use your SHERPA/RoMEO API key! See the [getting started vignette](rromeo.html) for usage of `rromeo` functions.


### Setting up your API key in an `.Renviron` file

`.Renviron` file follow the same loading rules as `.Rprofile` files, the only difference is that it is is a file whose only purpose is to store environment variables. To use your API key in an `.Renviron` you have to type the following (note the absence of quotes):

```
SHERPAROMEO_KEY=Iq83AIL5bss

```

You can then check that your key has been found using `check_key()`:


```r
check_key()
#> [1] "Iq83AIL5bs2"
```

Now `rromeo` can use your SHERPA/RoMEO API key! See the [getting started vignette](rromeo.html) for usage of `rromeo` functions.


### Using the `key` function argument (**NOT RECOMMENDED**)

**NOTE:** This method is **NOT RECOMMENDED** but still available for testing and development purposes. 

All the functions in `rromeo` that can use an API key have a `key` argument.
To use your API key you have to provide it as a string in the `key` argument of the functions for example:

```
rr_journal_issn("1947-6264", key = "YOUR_API_KEY")
```
However, **we do not recommend** this approach as your API key will be available in your R history and in your R scripts. While your API key should stay secret as it grants unlimited access to the server and can be maliciously used in wrong hands.

Furthermore using this method you have to specify it at each function call while the other methods shown above only needs to be set up once.


---
title: "Overview"
author:
  - "Matthias Grenié"
  - "Hugo Gruson"
date: "2020-03-05"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Query Journal Data

**NOTE:** SHERPA/RoMEO data is released under the [Creative Commons Attribution-NonCommercial-ShareAlike 2.5 (CC BY-NC-SA 2.5) license](https://creativecommons.org/licenses/by-nc-sa/2.5/). A suggestion of citation is included in `rromeo` via `citation("rromeo")`

The SHERPA/RoMEO database contains information on the archival policies of academic journals and publishers. Some journals let you archive a version of a submitted manuscript or book chapter (= preprint) before it's reviewed, some also let you archive the peer-reviewed but not formatted version (= postprint) and some even let you archive the peer-reviewed **and** formatted version of the manuscript (= pdf or publisher's version). The goal of `rromeo` is to make this database accessible through R.

### Query using ISSN

Let's try to access summary data using the ISSN of a journal:


```r
library(rromeo)

rr_journal_issn("1947-6264")
```

```
##                                                       title provided_issn      issn romeocolour preprint
## 1 A Critical Introduction to Media and Communication Theory     1947-6264 1947-6264      yellow      can
##    postprint        pdf pre_embargo post_embargo pdf_embargo
## 1 restricted restricted        <NA>    12 months   12 months
```

From this we see that the archival of the preprint is permitted while some restrictions apply for both postprint and pdf version of the manuscript. The restrictions are visible in the `*_embargo` field: the authors have to wait 12 months after publication before archiving postprint and formatted manuscript.

**Nota Bene**: the `rr_journal_issn()` function can use either the paper edition ISSN of the journal or the ISSN of the electronic version of the journal (e-ISSN or ESSN) which may differ


### Query by Journal Title

We can also search journal by names using `rr_journal_find()` the first argument `name` is the string to look for and the second argument `qtype` defines the way to match. If we want to look at journals that have the word `Biostatistics` in their title we can use the following code:


```r
rr_journal_find("Biostatistics", qtype = "contains")
```

```
## 17 journals match your query terms.
```

```
## Only titles and ISSNs of journals returned. Get more information using `rr_journal_name()`
```

```
##                                                             title provided_issn      issn
## 1                               American Journal of Biostatistics          <NA> 1948-9889
## 2                          Annals of Biometrics and Biostatistics          <NA> 2374-0116
## 4                             Austin Biometrics and Biostatistics          <NA> 2378-9840
## 5                Biometrics & biostatistics international journal          <NA>      <NA>
## 6                                                   Biostatistics          <NA> 1465-4644
## 7                                           Biostatistics -Basel-          <NA> 1527-2486
## 8                Biostatistics, bioinformatics and biomathematics          <NA> 0976-1594
## 9                                Edorium Journal of Biostatistics          <NA>      <NA>
## 10                             Enliven: Biostatistics and Metrics          <NA>      <NA>
## 11                   Epidemiology Biostatistics and Public Health          <NA> 2282-2305
## 12                         International Journal of Biostatistics          <NA> 1557-4679
## 13 International journal of clinical biostatistics and biometrics          <NA>      <NA>
## 14                        Journal of Biometrics and Biostatistics          <NA> 2155-6180
## 15                      Journal of epidemiology and biostatistics          <NA> 1359-5229
## 16                                    JP journal of biostatistics          <NA> 0973-5143
## 17                   Monographs in Epidemiology and Biostatistics          <NA> 0740-0845
```

`rr_journal_find()` will only return the list of titles and ISSN.
You can then use this list to select the exact journal you are looking for (in
this case, it is recommended to use the ISSN).
Alternatively, you may want to get data for all those journals. To achieve this,
you need to use the function `rr_journal_name()` which uses a similar syntax to `rr_journal_find()`:


```r
res <- rr_journal_name("Biostatistics", qtype = "contains")
```

```
## 17 journals match your query terms.
```

```
## Recursively fetching data from each journal. This may take some time...
```

```r
tail(res, 3)
```

```
##                                           title provided_issn      issn romeocolour preprint  postprint
## 14    Journal of epidemiology and biostatistics     1359-5229 1359-5229        <NA>     <NA>       <NA>
## 15                  JP journal of biostatistics     0973-5143 0973-5143        <NA>     <NA>       <NA>
## 16 Monographs in Epidemiology and Biostatistics     0740-0845 0740-0845      yellow      can restricted
##        pdf pre_embargo post_embargo pdf_embargo
## 14    <NA>        <NA>         <NA>        <NA>
## 15    <NA>        <NA>         <NA>        <NA>
## 16 unclear        <NA>    12 months        <NA>
```
The query may run for longer but it gives you all the information on all returned journals.

You can also use other type of name matching with other options for `qtype`, for example when `qtype = "starts"` the query string have to begin the title of the journal:


```r
rr_journal_name("Biostatistics", qtype = "starts")
```

```
## 3 journals match your query terms.
```

```
## Recursively fetching data from each journal. This may take some time...
```

```
##                                              title provided_issn      issn romeocolour preprint
## 1                                    Biostatistics     1465-4644 1465-4644       green      can
## 2                            Biostatistics -Basel-     1527-2486 1527-2486        <NA>     <NA>
## 3 Biostatistics, bioinformatics and biomathematics     0976-1594 0976-1594        <NA>     <NA>
##   postprint    pdf pre_embargo post_embargo pdf_embargo
## 1       can cannot        <NA>         <NA>        <NA>
## 2      <NA>   <NA>        <NA>         <NA>        <NA>
## 3      <NA>   <NA>        <NA>         <NA>        <NA>
```

While using `qtype = exact` the title of the journal should match exactly the used string:


```r
rr_journal_name("Biostatistics", qtype = "exact")
```

```
##           title provided_issn      issn romeocolour preprint postprint    pdf pre_embargo post_embargo
## 1 Biostatistics          <NA> 1465-4644       green      can       can cannot        <NA>         <NA>
##   pdf_embargo
## 1        <NA>
```


## Query Publisher Data

Not only can `rromeo` query journals archival policies but it can do the same for publishers. Indeed the SHERPA/RoMEO database contains much information on publishers' policies, with many ways to retrieve the information.

### Query by Publisher's Name

With `rr_publisher_name()` you can query publishers' information using the name of the publishers. The first argument `name` is the string to match publisher's names the second argument `qtype` gives the type of matching:


```r
rr_publisher_name(name = "Oxford", qtype = "all")
```

```
##   romeoid                                                       publisher alias romeocolour preprint
## 1    2632                                       Oxford Brookes University  <NA>        blue   cannot
## 2    1892                     Oxford Centre for Hebrew and Jewish Studies  <NA>        blue  unclear
## 3     986                       Oxford University Anthropological Society  <NA>       green      can
## 4      55                                         Oxford University Press   OUP      yellow      can
## 5    2072 University of Oxford, Oxford Uehiro Centre for Practical Ethics  <NA>       green      can
##    postprint     pdf
## 1        can     can
## 2        can  cannot
## 3        can     can
## 4 restricted unclear
## 5        can     can
```

When `qtype = "all"` the publishers' names should contain all the words, in any order, included in provided string:

```r
rr_publisher_name(name = "Oxford University", qtype = "all")
```

```
##   romeoid                                                       publisher alias romeocolour preprint
## 1    2632                                       Oxford Brookes University  <NA>        blue   cannot
## 2     986                       Oxford University Anthropological Society  <NA>       green      can
## 3      55                                         Oxford University Press   OUP      yellow      can
## 4    2072 University of Oxford, Oxford Uehiro Centre for Practical Ethics  <NA>       green      can
##    postprint     pdf
## 1        can     can
## 2        can     can
## 3 restricted unclear
## 4        can     can
```

While when `qtype = "exact"` the publishers' names should contain all the words in the same order as provided in the string:

```r
rr_publisher_name(name = "Oxford University", qtype = "exact")
```

```
##   romeoid                                 publisher alias romeocolour preprint  postprint     pdf
## 1     986 Oxford University Anthropological Society  <NA>       green      can        can     can
## 2      55                   Oxford University Press   OUP      yellow      can restricted unclear
```

Finally, when `qtype = "any"` publishers' names can contain any words of the provided string.


### Query by RoMEO's ID

The first column of publishers' policies data frame returned by `rromeo` is named `romeoid` it corresponds to the identifier used by SHERPA/RoMEO to identify publishers in a unique way:


```r
rr_publisher_id(id = 55)
```

```
##   romeoid               publisher alias romeocolour preprint  postprint     pdf
## 1      55 Oxford University Press   OUP      yellow      can restricted unclear
```


### Query by Country or Region

You can also query publishers information based on the country they are in using their two-letters ISO codes (see `?rr_publisher_country` for more information):


```r
rr_publisher_country(country = "IR")
```

```
##    romeoid                                                           publisher                   alias
## 1     1936                              Acta Advances in Agricultural Sciences                    <NA>
## 2     2786                                                   Avicenna Journals Advancements in Science
## 3     1556 Azad University, Mahlkesi Branch, Mechanical Engineering Department                    <NA>
## 4     2741                         Geodynamics Research International Bulletin                    <NA>
## 5     2159       International Journal of Management, Accounting and Economics                    <NA>
## 6     1823                           Islamic Azad University, Najafabad Branch                    <NA>
## 7     1745                            Islamic Azad University, Shahreza Branch                    <NA>
## 8     1726                                                       Islamic Press                    <NA>
## 9     1742                             Jahan Elm Institute of Higher Education                    <NA>
## 10    1414                               Kerman University of Medical Sciences                    <NA>
## 11    1867                                              Oboor Publishing Group                    <NA>
## 12    3078                      Shahid Sadoughi University of Medical Sciences                    <NA>
## 13    2084                                  Shefa Neuroscience Research Center                    <NA>
## 14    1362                                                   Shiraz University                    <NA>
## 15    1753                          Society of Diabetic Nephropathy Prevention                    <NA>
## 16    2922                               Tabriz University of Medical Sciences                    <NA>
## 17    1208                               Tehran University of Medical Sciences                    <NA>
## 18    2845                   University of Maragheh, Department of Mathematics                    <NA>
##    romeocolour preprint  postprint        pdf
## 1         blue   cannot        can        can
## 2        green      can        can        can
## 3         blue   cannot     cannot        can
## 4         blue   cannot     cannot        can
## 5        green      can        can        can
## 6         blue   cannot        can        can
## 7         blue  unclear        can        can
## 8         blue   cannot     cannot        can
## 9         blue  unclear        can        can
## 10        blue   cannot        can        can
## 11       green      can        can        can
## 12      yellow      can restricted restricted
## 13       green      can        can        can
## 14        blue   cannot        can        can
## 15        blue   cannot     cannot        can
## 16        blue   cannot        can        can
## 17       green      can        can        can
## 18       green      can        can        can
```

It is also possible to query publisher's information on a specific region or continent using `rr_publisher_continent()` (see the help page for the list of available regions):


```r
rr_publisher_continent(continent = "Australasia")
```

```
##    romeoid                                                                                publisher
## 1     1514                                                                                ANU Press
## 2     1853                                                                                 ANZAMEMS
## 3      959                                                                               ARRB Group
## 4      736                                                      Association of Occupational Science
## 5      222 Auckland University of Technology, School of Communication Studies, Pacific Media Centre
## 6     2343                                         Australasian Association for Information Systems
## 7     1358                                                             Australasian Medical Journal
## 8      220                     Australasian Society for Computers in Learning in Tertiary Education
## 9       98                                                                Australian Academic Press
## 10     518                                                           Australian Accoustical Society
## 11     638                                               Australian Clearinghouse for Youth Studies
## 12      97                                                          Australian Computer Society Inc
## 13    1813                                                 Australian International Academic Centre
## 14     255                                           Australian Library and Information Association
## 15     248                                                          Australian Mathematical Society
## 16     151                                                     Australian Physiotherapy Association
## 17     161                                                         Australian Psychological Society
## 18     383                                                 Australian Rock Art Research Association
## 19     290                                                      Australian Society of Anaesthetists
## 20    2420                                                                        Bareknuckle Books
## 21     150                          College of Intensive Care Medicine of Australia and New Zealand
## 22     160                                                                         CSIRO Publishing
## 23     553                                                                     e-Content Management
## 24     550                                                                Early Childhood Australia
## 25    2501                                                                   Edith Cowan University
## 26    2453                                                                                     EMAJ
## 27     571                                                       Field Naturalists Club of Victoria
## 28     626                    Griffith University, Griffith Law School, Socio-Legal Research Centre
## 29    1928                                                                           Infinity Press
## 30    1026                                                      Institute of Foresters of Australia
## 31     557                                       International Association for Statistics Education
## 32    3226                                                                     LexisNexis Australia
## 33     163                                                                       Libertas Academica
## 34     745                          Macquarie University, Department of International Communication
## 35     284                                                                           Magnolia Press
## 36    1397                                      Mathematics Education Research Group of Australasia
## 37    1310                                             Melbourne University, Law Review Association
## 38     665                                                                 Monash University ePress
## 39    2471                                                           New Zealand Ecological Society
## 40    2472                                                        New Zealand Institute of Forestry
## 41     598                                                          New Zealand Medical Association
## 42     623                                                          New Zealand Nurses Organisation
## 43     666                                                 New Zealand Society of Animal Production
## 44    2134                                                             Outdoor Council of Australia
## 45     670                                                                            Python Papers
## 46    2755                                                      Queensland University of Technology
## 47    2385                                      Queensland University of Technology, Faculty of Law
## 48     581                                                                          RMIT Publishing
## 49    1767                                       Royal New Zealand College of General Practitioners
## 50     264                                                             Royal Society of New Zealand
## 51    1329                                                                Royal Society of Victoria
## 52    1432                                                                              SETScholars
## 53    3041                                                             Studies in Material Thinking
## 54    3016                                                                        Sydney Law School
## 55     164                                                                          Thomson Reuters
## 56    1946                                                                               UTS ePRESS
## 57    2532                                   Western Australian Institutes for Educational Research
## 58     380                                 World Institute for Engineering and Technology Education
##                                                                           alias romeocolour preprint
## 1                                                                          <NA>       green      can
## 2  Australian and New Zealand Association for Medieval and Early Modern Studies       green      can
## 3                                                                          <NA>        blue   cannot
## 4                                                                          <NA>       white   cannot
## 5                                                                          <NA>        blue   cannot
## 6                                                                          AAIS       green      can
## 7                                                                          <NA>       green      can
## 8                                                                      ascilite       green      can
## 9                                                                          <NA>        blue   cannot
## 10                                                                         <NA>       white   cannot
## 11                                                                         <NA>       green      can
## 12                                                                         <NA>       green      can
## 13                                                                         AIAC       green      can
## 14                                                                         <NA>        blue  unclear
## 15                                                                         <NA>      yellow      can
## 16                                                                         <NA>        blue   cannot
## 17                                                                         <NA>       white  unclear
## 18                                                                         <NA>       white   cannot
## 19                                                                         <NA>       white   cannot
## 20                                                                         <NA>       green      can
## 21                                                                         <NA>       white   cannot
## 22                                                                         <NA>       green      can
## 23                                                                         <NA>       green      can
## 24                                                                         <NA>       white  unclear
## 25                                                                         <NA>        blue   cannot
## 26                                                                         <NA>       green      can
## 27                                                                         <NA>        blue  unclear
## 28                                                                         <NA>        blue   cannot
## 29                                                                         <NA>       green      can
## 30                                                                         <NA>       white   cannot
## 31                                                                         IASE        blue  unclear
## 32                                                                         <NA>        blue  unclear
## 33                                                                         <NA>       green      can
## 34                                                                         <NA>       white  unclear
## 35                                                                         <NA>       white   cannot
## 36                                                                         <NA>        blue  unclear
## 37                                                                         <NA>        blue   cannot
## 38                                                                         <NA>      yellow      can
## 39                                                                         <NA>      yellow      can
## 40                                                                         <NA>      yellow      can
## 41                                                                         <NA>       white   cannot
## 42                                                                         <NA>        blue   cannot
## 43                                                                         <NA>       white  unclear
## 44                                                                         <NA>        blue  unclear
## 45                                                                         <NA>        blue  unclear
## 46                                                                         <NA>       green      can
## 47                                                                         <NA>       green      can
## 48                                                                         <NA>       green      can
## 49                                                                         <NA>        blue   cannot
## 50                                                                         <NA>      yellow      can
## 51                                                                         <NA>       white   cannot
## 52                                                                         <NA>       green      can
## 53                                                                         <NA>        blue  unclear
## 54                                                                         <NA>        blue   cannot
## 55                                                                 Professional        blue   cannot
## 56                                                                         <NA>       green      can
## 57                                                                         <NA>        blue   cannot
## 58                                                                        WIETE       white   cannot
##     postprint        pdf
## 1         can       <NA>
## 2         can        can
## 3         can     cannot
## 4  restricted     cannot
## 5         can        can
## 6         can        can
## 7         can       <NA>
## 8         can        can
## 9         can     cannot
## 10     cannot       <NA>
## 11        can restricted
## 12        can       <NA>
## 13        can        can
## 14        can       <NA>
## 15     cannot restricted
## 16     cannot        can
## 17 restricted     cannot
## 18 restricted       <NA>
## 19 restricted       <NA>
## 20        can        can
## 21     cannot     cannot
## 22        can     cannot
## 23        can     cannot
## 24    unclear       <NA>
## 25        can        can
## 26        can        can
## 27        can        can
## 28        can        can
## 29        can        can
## 30 restricted restricted
## 31        can       <NA>
## 32        can restricted
## 33        can        can
## 34    unclear       <NA>
## 35     cannot     cannot
## 36     cannot        can
## 37        can        can
## 38     cannot       <NA>
## 39 restricted restricted
## 40 restricted restricted
## 41     cannot restricted
## 42     cannot        can
## 43 restricted       <NA>
## 44     cannot        can
## 45        can       <NA>
## 46        can        can
## 47        can        can
## 48        can     cannot
## 49        can        can
## 50 restricted     cannot
## 51     cannot     cannot
## 52        can        can
## 53     cannot        can
## 54     cannot        can
## 55        can       <NA>
## 56        can        can
## 57        can        can
## 58     cannot       <NA>
```

### Query by RoMEO colour

RoMEO assigns a colour depending on the different policies of publishers.

  | RoMEO colour | Archiving policy                                        |
  |:-------------|:--------------------------------------------------------|
  | `green`      | can archive preprint, postprint and publisher's version |
  | `blue`       | can archive postprint **or** publisher's version        |
  | `yellow`     | can archive preprint                                    |
  | `white`      | archiving not formally supported                        |

(Table taken from <http://www.sherpa.ac.uk/romeo/definitions.php#colours>)

You can query journals using this classification with the function `rr_romeo_colour()`:


```r
rr_romeo_colour(romeo_colour = "green")
```

In this example vignette we do not run the query because it can run for quite long as it returns the policies of all publishers of the given colour (you can see the numbers of publishers in each category in the following web page <http://www.sherpa.ac.uk/romeo/statistics.php?la=en&fIDnum=|&mode=simple>).


## Setting up an API key

SHERPA/RoMEO lets you make 500 queries per day per IP address for free. If you get past this limit you will get the following error:
```
You have exceeded the free use limit of 500 requests per day. To go beyond this limit you should register for a free API key available at http://www.sherpa.ac.uk/romeo/apiregistry.php
```
We encourage you to register for a **free** API key at the above-mentioned address.

To provide your API key to `rromeo` you can:
1. provide it as a character string as the `key` arguments of `rromeo` functions as `rr_*(..., key = "my_key_as_a_string")`, we **do not recommend** this method as your API key will be available from your history;
1. you can define the variable `SHERPAROMEO_KEY` in an `.Renviron` file in your working directory, the file should contain the following line
  `SHERPAROMEO_KEY=my_key_without_quotes`;
1. you can also define the variable `SHERPAROMEO_KEY` in an `.Rprofile` file in your working directory, the file should contain the following line
  `SHERPAROMEO_KEY="my_key_with_quotes"`.

## Visualization Example

`rromeo` can be quite useful in bibliometric studies to report on archival policies. For example, we can now have a quick visual overview of the policies of journals in a given field from the results obtained in a query:


```r
library(ggplot2)
theme_set(theme_minimal())

stacked_res <- stack(res[, 4:6])
ggplot(stacked_res, aes(x = ind, fill = values)) +
  geom_bar() +
  labs(x = NULL,
       subtitle = "Archiving Policies of Journals with 'Biostatistics' in their title")
```

![plot of chunk bar_graph](bar_graph-1.png)


```r
ggplot(res, aes(x = "a", fill = romeocolour)) +
  geom_bar() +
  coord_polar("y") +
  labs(x = NULL,
       subtitle = "RoMEO colours of Journals with 'Biostatistics' in their title")
```

![plot of chunk pie_chart](pie_chart-1.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_generic.R
\name{parse_generic}
\alias{parse_generic}
\title{Generic parsing function}
\usage{
parse_generic(api_answer, ...)
}
\arguments{
\item{api_answer}{[\code{httr::response()}]\cr{}
The API answer}

\item{...}{Other options passed to parsing functions}
}
\value{
either results from \code{\link[=parse_journal]{parse_journal()}} or \code{\link[=parse_publisher]{parse_publisher()}}
}
\description{
Generic parsing function
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_romeo_colour.R
\name{rr_romeo_colour}
\alias{rr_romeo_colour}
\alias{rr_romeo_color}
\title{Query publisher by RoMEO colour}
\usage{
rr_romeo_colour(
  romeo_colour = c("green", "blue", "yellow", "white"),
  key = NULL
)
}
\arguments{
\item{romeo_colour}{[\code{character(1)}]\cr{}
in \code{c("green", "blue", "yellow", "white")}
the SHERPA/RoMEO colour to retrieve}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data frame with the following columns:
\itemize{
\item \code{romeoid}     [\code{integer(1)}]\cr{}
the internal index of the publisher in the SHERPA/RoMEO
database
\item \code{publisher}   [\code{character(1)}]\cr{}
the name of the publisher
\item \code{alias}       [\code{character(1)}]\cr{}
if applicable an alternative name of the publisher or the
name of the specific publishing branch
\item \code{romeocolour} [\code{character(1)}]\cr{}
a colour assigned by the database that reflects the default
policies of the publisher
\item \code{preprint}    [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}   [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}         [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
archivable?
}
}
\description{
SHERPA/RoMEO classifies publisher in different colours depending on their
archiving policies.
\itemize{
\item \strong{green} publishers let authors archive preprint and postprint or
publisher's version/PDF,
\item \strong{blue} publishers let authors archive postprint or publisher's
version/PDF,
\item \strong{yellow} publishers let authors archive preprint,
\item \strong{white} publishers do not formally support archival.
}
}
\details{
For more details about the definitions of RoMEO colours check the
\href{http://sherpa.ac.uk/romeo/definitions.php#colours}{FAQ section} of
SHERPA/RoMEO

Note that when using \code{\link[=rr_romeo_colour]{rr_romeo_colour()}} the API returns \strong{all} the
publishers in the selected category, so the results are generally bigger in
size than specific functions like \code{\link[=rr_journal_name]{rr_journal_name()}} or \code{\link[=rr_publisher_id]{rr_publisher_id()}}
}
\examples{
\donttest{
rr_romeo_colour(romeo_colour = "green")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_key}
\alias{check_key}
\title{Check SHERPA/RoMEO API key}
\usage{
check_key(key = NULL)
}
\arguments{
\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
if found the character string of the key, \code{NULL} otherwise
}
\description{
The key can be either specified in various ways see the Details section.
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_publisher_all.R
\name{rr_publisher_all}
\alias{rr_publisher_all}
\title{Get all Publisher Policies}
\usage{
rr_publisher_all(key = NULL)
}
\arguments{
\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data frame with the following columns:
\itemize{
\item \code{romeoid}     [\code{integer(1)}]\cr{}
the internal index of the publisher in the SHERPA/RoMEO
database
\item \code{publisher}   [\code{character(1)}]\cr{}
the name of the publisher
\item \code{alias}       [\code{character(1)}]\cr{}
if applicable an alternative name of the publisher or the
name of the specific publishing branch
\item \code{romeocolour} [\code{character(1)}]\cr{}
a colour assigned by the database that reflects the default
policies of the publisher
\item \code{preprint}    [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}   [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}         [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
archivable?
}
}
\description{
Retrieve all data on publishers policies from SHERPA/RoMEO.
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_journal_find.R
\name{rr_journal_find}
\alias{rr_journal_find}
\title{Find if journals are available in SHERPA/RoMEO}
\usage{
rr_journal_find(name, qtype = c("exact", "contains", "starts"), key = NULL)
}
\arguments{
\item{name}{[\verb{character(1+)}]\cr{}
one or several strings to match the titles of the journals}

\item{qtype}{[\code{character(1)}]\cr{}
in:
* \code{"exact"} full title must be exactly to provided \code{name},
* \code{"contains"} the provided \code{name} must appear anywhere in the
title of the journal,
* \code{"starts"} the provided \code{name} must appear at the start of
title of the journal.}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data frame:
\itemize{
\item \code{title}         [\code{character(1)}]\cr{}
the name of the journal
\item \code{provided_issn} [\code{character(1)}]\cr{}
the ISSN you provided in your query (might differ from the
ISSN returned by the API)
\item \code{issn}          [\code{character(1)}]\cr{}
the ISSN of the journal
}
}
\description{
Find if journals are available in SHERPA/RoMEO
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
\examples{
\donttest{
rr_journal_find(name = "Biostatistics", qtype = "contains")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_publisher_continent.R
\name{rr_publisher_continent}
\alias{rr_publisher_continent}
\title{Get Publisher Policy by Publisher's Continent}
\usage{
rr_publisher_continent(
  continent = c("Africa", "Antarctica", "Asia", "Australasia", "Caribbean",
    "Central America", "Europe", "North America", "Oceania", "South America"),
  key = NULL
)
}
\arguments{
\item{continent}{[\verb{character(1+)}]\cr{}
one or a vector of strings in \code{c("Africa", "Antarctica", "Asia", "Australasia", "Carribean", "Central America", "Europe", "North America", "Oceania", "South America")}\cr{}
the continent name to retrieve}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data frame with the following columns:
\itemize{
\item \code{romeoid}     [\code{integer(1)}]\cr{}
the internal index of the publisher in the SHERPA/RoMEO
database
\item \code{publisher}   [\code{character(1)}]\cr{}
the name of the publisher
\item \code{alias}       [\code{character(1)}]\cr{}
if applicable an alternative name of the publisher or the
name of the specific publishing branch
\item \code{romeocolour} [\code{character(1)}]\cr{}
a colour assigned by the database that reflects the default
policies of the publisher
\item \code{preprint}    [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}   [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}         [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
archivable?
}
}
\description{
Retrieve publisher's policy based on publisher's continent. This function
does not work for unclassified or international publishers.
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
\examples{
\donttest{
rr_publisher_continent(continent = "Caribbean")
rr_publisher_continent(continent = "Central America")
rr_publisher_continent(continent = c("Caribbean", "Central America"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_api_version.R
\name{rr_api_version}
\alias{rr_api_version}
\title{Return SHERPA/RoMEO API version}
\usage{
rr_api_version()
}
\description{
This function queries SHERPA/RoMEO and returns the version of the API.
}
\examples{
\donttest{
rr_api_version()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_publisher_country.R
\name{rr_publisher_country}
\alias{rr_publisher_country}
\title{Get Publisher Policy by Publisher's Country}
\usage{
rr_publisher_country(country, key = NULL)
}
\arguments{
\item{country}{[\verb{character(1+)}]\cr{}
one or a vector of ISO two-letter country code or \code{AA} for
international publisher, \code{ZZ} for publisher of unknown
countries and \verb{__} for publishers without specified country
(case insensitive).}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data frame with the following columns:
\itemize{
\item \code{romeoid}     [\code{integer(1)}]\cr{}
the internal index of the publisher in the SHERPA/RoMEO
database
\item \code{publisher}   [\code{character(1)}]\cr{}
the name of the publisher
\item \code{alias}       [\code{character(1)}]\cr{}
if applicable an alternative name of the publisher or the
name of the specific publishing branch
\item \code{romeocolour} [\code{character(1)}]\cr{}
a colour assigned by the database that reflects the default
policies of the publisher
\item \code{preprint}    [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}   [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}         [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
archivable?
}
}
\description{
Retrieve publisher's policy based on publisher's country. The code should be
the ISO_3166-1_alpha-2 code of the country
\url{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}.
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
\examples{
\donttest{
# Taiwan
rr_publisher_country("TW")
# Egypt
rr_publisher_country("EG")
rr_publisher_country(c("TW", "EG"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validate_country_code}
\alias{validate_country_code}
\title{Validate ISO two-letters country code}
\usage{
validate_country_code(country)
}
\arguments{
\item{country}{[\code{character(1)}]\cr{}
a two-letter country code or \code{AA}, \code{ZZ} or
\verb{__} (special country codes for SHERPA/RoMEO)}
}
\value{
\code{TRUE} if the country code is valid, errors otherwise
}
\description{
If available uses \code{\link[ISOcodes:ISO_3166]{ISOcodes::ISO_3166_1}} to validate country code.
Otherwise assume that the code is valid as long as it is a two-letter code or
\verb{__}. See \code{\link[=rr_publisher_country]{rr_publisher_country()}} for use of country codes.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_journal.R
\name{parse_journal}
\alias{parse_journal}
\title{Parse API answer}
\usage{
parse_journal(xml_source, outcome, hits, type = c("find", "name"), key = NULL)
}
\arguments{
\item{type}{[\code{character(1)} in \code{c("find", "name")}]\cr{}
If \code{type = "find"} returns only \code{title} and \code{issn} columns if
\code{type = "name"} returns full data.frame as specified in Returns
sections.}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data.frame with the following columns:
\itemize{
\item \code{title}         [\code{character(1)}]\cr{}
the name of the journal
\item \code{provided_issn} [\code{character(1)}]\cr{}
the ISSN you provided in your query (might differ from the
ISSN returned by the API)
\item \code{issn}          [\code{character(1)}]\cr{}
the ISSN of the journal
\item \code{romeocolour}   [\code{character(1)}]\cr{}
the SHERPA/RoMEO colour of the journal
\item \code{preprint}      [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}     [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}           [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
\item \code{pre_embargo}   [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the preprint
\item \code{post_embargo}  [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the postprint, if value is \code{"after media"}, it
means that the post-print can be archived after media
embargo has passed
\item \code{pdf_embargo}   [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the publisher's version, if value is \code{"after media"},
it means that the publisher's version can be archived after
media embargo has passed
}
}
\description{
Returns data.frame from parsed xml API answer.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rromeo-package.R
\docType{package}
\name{rromeo-package}
\alias{rromeo}
\alias{rromeo-package}
\title{rromeo: Access Publisher Copyright & Self-Archiving Policies via the 
  'SHERPA/RoMEO' API}
\description{
Fetches information from the 'SHERPA/RoMEO' API
  <http://www.sherpa.ac.uk/romeo/apimanual.php> which indexes policies of
  journal regarding the archival of scientific manuscripts before and/or after
  peer-review as well as formatted manuscripts.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/rromeo/}
  \item \url{https://github.com/ropensci/rromeo}
  \item Report bugs at \url{https://github.com/ropensci/rromeo/issues}
}

}
\author{
\strong{Maintainer}: Matthias Grenié \email{matthias.grenie@gmail.com} (\href{https://orcid.org/0000-0002-4659-7522}{ORCID})

Authors:
\itemize{
  \item Hugo Gruson (\href{https://orcid.org/0000-0002-4094-1476}{ORCID})
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_base_functions.R
\name{rr_GET}
\alias{rr_GET}
\title{rromeo internal GET function}
\usage{
rr_GET(...)
}
\arguments{
\item{...}{additional parameter to \code{\link[httr:GET]{httr::GET}}}
}
\description{
rromeo internal GET function
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_journal_name.R
\name{rr_journal_name}
\alias{rr_journal_name}
\title{Retrieve journals policies by matching title}
\usage{
rr_journal_name(name, qtype = c("exact", "contains", "starts"), key = NULL)
}
\arguments{
\item{name}{[\verb{character(1+)}]\cr{}
one or several strings to match the titles of the journals}

\item{qtype}{[\code{character(1)}]\cr{}
in:
* \code{"exact"} full title must be exactly to provided \code{name},
* \code{"contains"} the provided \code{name} must appear anywhere in the
title of the journal,
* \code{"starts"} the provided \code{name} must appear at the start of
title of the journal.}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data.frame with the following columns:
\itemize{
\item \code{title}         [\code{character(1)}]\cr{}
the name of the journal
\item \code{provided_issn} [\code{character(1)}]\cr{}
the ISSN you provided in your query (might differ from the
ISSN returned by the API)
\item \code{issn}          [\code{character(1)}]\cr{}
the ISSN of the journal
\item \code{romeocolour}   [\code{character(1)}]\cr{}
the SHERPA/RoMEO colour of the journal
\item \code{preprint}      [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}     [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}           [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
\item \code{pre_embargo}   [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the preprint
\item \code{post_embargo}  [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the postprint, if value is \code{"after media"}, it
means that the post-print can be archived after media
embargo has passed
\item \code{pdf_embargo}   [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the publisher's version, if value is \code{"after media"},
it means that the publisher's version can be archived after
media embargo has passed
}
}
\description{
Note that SHERPARoMEO will not return more than 50 journals in a single
query. The function will warn you if you are in this case.
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
\examples{
\donttest{
rr_journal_name(name = "Journal of Geology")
rr_journal_name(name = "Biogeography", qtype = "contains")
# You can also query multiple journals with exact titles in a single call
rr_journal_name(name = c("Journal of Biogeography", "PLoS ONE"),
                qtype = "exact")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_journal_issn.R
\name{rr_journal_issn}
\alias{rr_journal_issn}
\title{Retrieve journal policy using ISSN}
\usage{
rr_journal_issn(issn, key = NULL)
}
\arguments{
\item{issn}{[\verb{character(1+)}]\cr{}
one or a vector of journal(s) ISSN(s) or ESSN(s)}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data.frame with the following columns:
\itemize{
\item \code{title}         [\code{character(1)}]\cr{}
the name of the journal
\item \code{provided_issn} [\code{character(1)}]\cr{}
the ISSN you provided in your query (might differ from the
ISSN returned by the API)
\item \code{issn}          [\code{character(1)}]\cr{}
the ISSN of the journal
\item \code{romeocolour}   [\code{character(1)}]\cr{}
the SHERPA/RoMEO colour of the journal
\item \code{preprint}      [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}     [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}           [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
\item \code{pre_embargo}   [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the preprint
\item \code{post_embargo}  [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the postprint, if value is \code{"after media"}, it
means that the post-print can be archived after media
embargo has passed
\item \code{pdf_embargo}   [\code{character(1)}]\cr{}
if applicable the embargo period before the author(s) can
archive the publisher's version, if value is \code{"after media"},
it means that the publisher's version can be archived after
media embargo has passed
}
}
\description{
Retrieve policy information from the SHERPA/RoMEO API using the ISSN
from the paper edition of the journal or the ISSN of the electronic version
(e-ISSN or ESSN)
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
\examples{
\donttest{
# Query single ISSN
rr_journal_issn(issn = "1947-6264")

# Query multiple ISSN
rr_journal_issn(issn = c("1947-6264", "0030-1299"))

# Query by ESSN
rr_journal_issn("1463-9084")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_auth.R
\name{rr_auth}
\alias{rr_auth}
\title{Store provided API key into Environment Variable}
\usage{
rr_auth(key)
}
\arguments{
\item{key}{[\code{character(1)}]\cr{}
A string giving the API key to save into the environment}
}
\description{
This function stores the provided API key as argument in to an environment
variable \code{SHERPAROMEO_KEY} for further use by other \code{rromeo} functions.
}
\details{
For more information regarding API keys, please refer to dedicated vignette
with the following command
\code{vignette("setting_up_api_key", package = "rromeo")}
}
\examples{
\dontrun{
rr_auth("Iq83AIL5bss")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_publisher_id.R
\name{rr_publisher_id}
\alias{rr_publisher_id}
\title{Get Publisher Policy from Publisher ID}
\usage{
rr_publisher_id(id, key = NULL)
}
\arguments{
\item{id}{[\verb{integer(1+)}]\cr{}
one or a vector of SHERPA/RoMEO publisher's ID}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data frame with the following columns:
\itemize{
\item \code{romeoid}     [\code{integer(1)}]\cr{}
the internal index of the publisher in the SHERPA/RoMEO
database
\item \code{publisher}   [\code{character(1)}]\cr{}
the name of the publisher
\item \code{alias}       [\code{character(1)}]\cr{}
if applicable an alternative name of the publisher or the
name of the specific publishing branch
\item \code{romeocolour} [\code{character(1)}]\cr{}
a colour assigned by the database that reflects the default
policies of the publisher
\item \code{preprint}    [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}   [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}         [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
archivable?
}
}
\description{
Use SHERPA/RoMEO API to retrieve a specific publisher policies on manuscript
archival
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
\examples{
\donttest{
rr_publisher_id(id = 55)
rr_publisher_id(id = c(55, 735))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_publisher.R
\name{parse_publisher}
\alias{parse_publisher}
\title{Parse publisher list}
\usage{
parse_publisher(xml_source, outcome, hits)
}
\arguments{
\item{api_answer}{xml API answer}
}
\value{
Returns a data frame with the following columns:
\itemize{
\item \code{romeoid}     [\code{integer(1)}]\cr{}
the internal index of the publisher in the SHERPA/RoMEO
database
\item \code{publisher}   [\code{character(1)}]\cr{}
the name of the publisher
\item \code{alias}       [\code{character(1)}]\cr{}
if applicable an alternative name of the publisher or the
name of the specific publishing branch
\item \code{romeocolour} [\code{character(1)}]\cr{}
a colour assigned by the database that reflects the default
policies of the publisher
\item \code{preprint}    [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}   [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}         [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
archivable?
}
}
\description{
When the returned content by the SHERPA/RoMEO API is a list of publishers
this function parses the list and returns a structured data.frame with
the default policies of the different publisher.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validate_issn}
\alias{validate_issn}
\title{Checks validity of the ISSN}
\usage{
validate_issn(issn)
}
\arguments{
\item{issn}{[\code{character(1)}]\cr{}
The ISSN of a journal}
}
\value{
\code{NULL} if the ISSN is valid, errors otherwise
}
\description{
International Standard Serial Numbers (ISSNs) have a specific structure and
this function checks that the provided ISSN is valid. For more information
please go to
\url{https://en.wikipedia.org/wiki/International_Standard_Serial_Number}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{parse_embargo}
\alias{parse_embargo}
\title{Parse embargo period from API return}
\usage{
parse_embargo(xml_source, type = c("pre", "post", "pdf"))
}
\arguments{
\item{xml_source}{[\code{xml_document}]\cr{}
a parsed xml document}

\item{type}{[\code{character(1)}]
name of the embargo type must be in \code{"pre"}, \code{"post"}, and
\code{"pdf"}}
}
\value{
the embargo period as a string like \code{"12 months"}
}
\description{
This function provides an easy way to return the embargo period in the
different categories
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_publisher_name.R
\name{rr_publisher_name}
\alias{rr_publisher_name}
\title{Get Publisher Policy by Publisher Name}
\usage{
rr_publisher_name(name, qtype = c("all", "any", "exact"), key = NULL)
}
\arguments{
\item{name}{[\verb{character(1+)}]\cr{}
One or a vector of query string(s) to search publisher name}

\item{qtype}{[\code{character(1)}]\cr{}
in \code{c("all", "any", "exact")} define the type of matching:
\itemize{
\item \code{all} means that all strings in \code{name} must appear in any
order or location
\item \code{any} means that at least one of the strings in \code{name} must
appear
\item \code{exact} means that the \code{name} string must appear in the
publisher's name or its alias.
}}

\item{key}{[\code{character(1)}]\cr{}
a character string containing the API key or \code{NULL}
(see Details section on how to specify it)}
}
\value{
Returns a data frame with the following columns:
\itemize{
\item \code{romeoid}     [\code{integer(1)}]\cr{}
the internal index of the publisher in the SHERPA/RoMEO
database
\item \code{publisher}   [\code{character(1)}]\cr{}
the name of the publisher
\item \code{alias}       [\code{character(1)}]\cr{}
if applicable an alternative name of the publisher or the
name of the specific publishing branch
\item \code{romeocolour} [\code{character(1)}]\cr{}
a colour assigned by the database that reflects the default
policies of the publisher
\item \code{preprint}    [\code{character(1)}]\cr{}
is the preprint (not reviewed) archivable?
\item \code{postprint}   [\code{character(1)}]\cr{}
is the postprint (reviewed but not formatted) archivable?
\item \code{pdf}         [\code{character(1)}]\cr{}
is the publisher's version (reviewed and formatted)
archivable?
}
}
\description{
Use SHERPA/RoMEO API to retrieve a specific publisher policies on manuscript
archival based on matching the name of the publishers.
}
\details{
There are several ways to provide your API key.
The best way to know about them is to refer to the vignette about
"Setting Up Your API key" accessible with the following command:
\code{vignette("setting_up_api_key", package = "rromeo")}.
You can also use \code{\link{rr_auth}} that will use the provided key to store it as
an environmental variable.
}
\examples{
\donttest{
rr_publisher_name(name = "Optical Society", qtype = "all")
rr_publisher_name(name = "Swiss Chemistry", qtype = "any")
rr_publisher_name(name = "Swiss Chemistry", qtype = "exact")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr_base_functions.R
\name{rr_ua}
\alias{rr_ua}
\title{rromeo User Agent}
\usage{
rr_ua()
}
\description{
rromeo User Agent
}
