
# virtuoso <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Travis build
status](https://travis-ci.org/ropensci/virtuoso.svg?branch=master)](https://travis-ci.org/ropensci/virtuoso)
[![Build
status](https://ci.appveyor.com/api/projects/status/x5r18x1cvu6khksd/branch/master?svg=true)](https://ci.appveyor.com/project/cboettig/virtuoso/branch/master)
[![Coverage
status](https://codecov.io/gh/ropensci/virtuoso/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/virtuoso?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/virtuoso)](https://cran.r-project.org/package=virtuoso)
[![Peer
review](http://badges.ropensci.org/271_status.svg)](https://github.com/ropensci/software-review/issues/271)

<!-- README.md is generated from README.Rmd. Please edit that file -->

The goal of virtuoso is to provide an easy interface to Virtuoso RDF
database from R.

## Installation

You can install the development version of virtuoso from GitHub with:

``` r
remotes::install_github("ropensci/virtuoso")
```

## Getting Started

``` r
library(virtuoso)
```

For Mac users, `virtuoso` package includes a utility function to install
and configure a local Virtuoso Open Source instance using Homebrew.
Otherwise, simply install the Virtuoso Open Source edition for your
operating system.

``` r
vos_install()
#> Virtuoso is already installed.
```

We can now start our Virtuoso server from R:

``` r
vos_start()
#> PROCESS 'virtuoso-t', running, pid 14318.
#> Server is now starting up, this may take a few seconds...
#> latest log entry: 21:43:06 Server online at 1111 (pid 14318)
```

Once the server is running, we can connect to the database.

``` r
con <- vos_connect()
```

Our connection is now live, and accepts SPARQL queries directly.

``` r
DBI::dbGetQuery(con, "SPARQL SELECT * WHERE { ?s ?p ?o } LIMIT 4")
#>                                                                              s
#> 1                   http://www.openlinksw.com/virtrdf-data-formats#default-iid
#> 2          http://www.openlinksw.com/virtrdf-data-formats#default-iid-nullable
#> 3          http://www.openlinksw.com/virtrdf-data-formats#default-iid-nonblank
#> 4 http://www.openlinksw.com/virtrdf-data-formats#default-iid-nonblank-nullable
#>                                                 p
#> 1 http://www.w3.org/1999/02/22-rdf-syntax-ns#type
#> 2 http://www.w3.org/1999/02/22-rdf-syntax-ns#type
#> 3 http://www.w3.org/1999/02/22-rdf-syntax-ns#type
#> 4 http://www.w3.org/1999/02/22-rdf-syntax-ns#type
#>                                                         o
#> 1 http://www.openlinksw.com/schemas/virtrdf#QuadMapFormat
#> 2 http://www.openlinksw.com/schemas/virtrdf#QuadMapFormat
#> 3 http://www.openlinksw.com/schemas/virtrdf#QuadMapFormat
#> 4 http://www.openlinksw.com/schemas/virtrdf#QuadMapFormat
```

## DSL

`virtuoso` also provides wrappers around some common queries to make it
easier to work with Virtuoso and RDF.

The bulk loader can be used to quickly import existing sets of triples.

``` r
example <- system.file("extdata", "person.nq", package = "virtuoso")
vos_import(con, example)
```

Can also read in compressed formats as well. Remember to set the pattern
match appropriately. This is convenient because N-Quads compress
particularly well, often by a factor of 20 (or rather, can be
particularly large when uncompressed, owing to the repeated property and
subject URIs).

``` r
ex <- system.file("extdata", "library.nq.gz", package = "virtuoso")
vos_import(con, ex)
```

`vos_import` invisibly returns a table of the loaded files, with error
message and loading times. If a file cannot be imported, an error
message is returned:

``` r
bad_file <- system.file("extdata", "bad_quads.nq", package = "virtuoso")
vos_import(con, bad_file)
#> Error: Error importing: bad_quads.nq 37000 SP029: NQuads RDF loader, line 2: Undefined namespace prefix at ITIS:1000000
```

We can now query the imported data using SPARQL.

``` r
df <- vos_query(con, 
"SELECT ?p ?o 
 WHERE { ?s ?p ?o .
        ?s a <http://schema.org/Person>
       }")
head(df)
#>                                                 p                        o
#> 1 http://www.w3.org/1999/02/22-rdf-syntax-ns#type http://schema.org/Person
#> 2                      http://schema.org/jobTitle                Professor
#> 3                          http://schema.org/name                 Jane Doe
#> 4                     http://schema.org/telephone           (425) 123-4567
#> 5                           http://schema.org/url   http://www.janedoe.com
```

``` r
vos_query(con, 
"SELECT ?p ?o 
 WHERE { ?s ?p ?o .
        ?s a <http://example.org/vocab#Chapter>
       }")
#>                                                 p
#> 1 http://www.w3.org/1999/02/22-rdf-syntax-ns#type
#> 2     http://purl.org/dc/elements/1.1/description
#> 3           http://purl.org/dc/elements/1.1/title
#>                                          o
#> 1         http://example.org/vocab#Chapter
#> 2 An introductory chapter on The Republic.
#> 3                         The Introduction
```

## Server controls

We can control any `virtuoso` server started with `vos_start()` using a
series of helper commands.

``` r
vos_status()
#> latest log entry: 21:43:06 PL LOG: No more files to load. Loader has finished,
#> [1] "sleeping"
```

Advanced usage note: `vos_start()` invisibly returns a `processx` object
which we can pass to other server control functions, or access the
embedded `processx` control methods directly. The `virtuoso` package
also caches this object in an environment so that it can be accessed
directly without having to keep track of an object in the global
environment. Use `vos_process()` to return the `processx` object. For
example:

``` r
library(ps)
p <- vos_process()
ps_is_running(p)
#> [1] TRUE
ps_cpu_times(p)
#>            user          system   children_user children_system 
#>            1.61            0.29            0.00            0.00
ps_suspend(p)
#> NULL
ps_resume(p)
#> NULL
```

## Going further

Please see the package vignettes for more information:

  - [details on Virtuoso Installation &
    configuration](https://docs.ropensci.org/virtuoso/articles/installation.html)
  - [The Data Lake: richer examples of RDF
    use](https://docs.ropensci.org/virtuoso/articles/articles/datalake.html)

-----

Please note that the `virtuoso` R package is released with a
[Contributor Code of
Conduct](https://docs.ropensci.org/virtuoso/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# virtuoso 0.1.8

* Bugfix for CRAN's Solaris configuration

# virtuoso 0.1.7

* `has_virtuoso()` sets PATH on Mac & Windows so that it may evaluate correctly
  without `start_virtuoso()`.  (See issue #35)

# virtuoso 0.1.6

* compatibility with recent changes in rappdirs

# virtuoso 0.1.5

* update documentation links & return to CRAN

# virtuoso 0.1.4

* testing compatibility for Solaris platforms (for CRAN)

# 0.1.3

* tweaks to testing and examples for CRAN
* Revised DESCRIPTION text for CRAN

#  0.1.2

(Version after peer review by rOpenSci)

* See <https://github.com/ropensci/software-review/issues/271>
* Added a `NEWS.md` file to track changes to the package.

# 0.1.1, 2018-12-03

(Version submitted to rOpenSci)
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
Dear CRAN,

This release makes the changes indicated in NEWS.md.

Sincerely,

Carl Boettiger
# Contributing to virtuoso

This outlines how to propose a change to virtuoso. For more detailed
info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

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
*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the virtuoso project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See tidyverse [development contributing guide](https://rstd.io/tidy-contrib)
for further details.
Submitting Author: Carl Boettiger (@cboettig)
Repository: <https://github.com/cboettig/virtuoso>
Version submitted: 0.1.1 (tagged)
Editor: TBD
Reviewer 1: TBD
Reviewer 2: TBD
Archive: TBD
Version accepted: TBD

---

-   Paste the full DESCRIPTION file inside a code block below:

```
Package: virtuoso
Type: Package
Title: R interface to Virtuoso using ODBC
Version: 0.1.1
Authors@R: c(person("Carl", "Boettiger", 
                  email = "cboettig@gmail.com", 
                  role = c("aut", "cre", "cph"),
                  comment = c(ORCID = "0000-0002-1642-628X")),
             person("Bryce", "Mecum", 
                    role = "ctb", 
                    email = "brycemecum@gmail.com",
                    comment = c(ORCID = "0000-0002-0381-3766")))
Description: Virtuoso is a high-performance "universal server," which can act
             as both a relational database (supporting standard SQL queries),
             and an Resource Description Framework (RDF) triplestore, supporting 
             SPARQL queries and semantic reasoning. The virtuoso package R provides
             R users with a DBI-compatible connection to the Virtuoso database. 
             The package also provides helper routines to install, launch, and manage
             a Virtuoso server locally on Mac, Windows and Linux platforms using
             the standard interactive installers from the R command-line.  By 
             automatically handling these setup steps, the package can make Virtuoso
             considerably faster and easier for a most users to deploy in a local
             environment. While this can be used as a normal dplyr backend, Virtuoso 
             excels when used as a RDF triplestore.  Managing the bulk import of triples
             from common serializations with a single intuitive command is another key
             feature of the Virtuoso R package.  Bulk import performance can be tens to
             hundreds of times faster than the comparable imports using existing R tools,
             including rdflib and redland packages.  
License: MIT + file LICENSE
Encoding: UTF-8
LazyData: true
Imports: 
    odbc,
    processx,
    DBI,
    utils,
    ini,
    rappdirs,
    curl,
    fs,
    digest
RoxygenNote: 6.1.1
Suggests: 
    knitr,
    rmarkdown,
    nycflights13,
    testthat,
    covr,
    jsonld,
    rdftools,
    dplyr,
    spelling
VignetteBuilder: knitr
Remotes: cboettig/rdftools
Language: en-US

```


## Scope

- Please indicate which category or categories from our [package fit policies](https://ropensci.github.io/dev_guide/policies.html#package-categories) this package falls under: (Please check an appropriate box below. If you are unsure, we suggest you make a pre-submission inquiry.):

	- [ ] data retrieval
	- [x] data extraction
	- [x] database access
	- [ ] data munging
	- [ ] data deposition
	- [ ] reproducibility
	- [ ] geospatial data
	- [ ] text analysis


- Explain how the and why the package falls under these categories (briefly, 1-2 sentences):

R users confronted with large dump of triples (e.g. nquad, owl, or other file) currently have few ways of reading in this data,
and no performant option that can handle the huge file sizes frequently involved that do not fit into memory. This package provides
a relatively convenient way to import this data into an RDF-capable database and query that data directly from R.  

-   Who is the target audience and what are scientific applications of this package?

Researchers working with RDF / semantic data.

-   Are there other R packages that accomplish the same thing? If so, how does yours differ or meet [our criteria for best-in-category](https://github.com/ropensci/onboarding/blob/master/policies.md#overlap)?

This package overlaps with ropensci package `rdflib` (and thus `redland`, which is `rdflib` uses under the hood.), which primarily provides an in-memory model for working with RDF data, which fails with large triplestores.  `rdflib` & `redland` do have a pluggable backend that can connect to Virtuoso and other databases, but this is not only very complicated to set up (not only does `redland` R package need to be built from source, but so does the redland C library in some cases) but is also much slower.  This package handles the installation easily in a user-friendly and more performant way, and the resulting Virtuoso server can then be used as a backend to `rdflib` (though there is usually little reason to do so since Virtuoso can be called directly though this package already).  

-   If you made a pre-submission enquiry, please paste the link to the corresponding issue, forum post, or other discussion, or @tag the editor you contacted.

## Technical checks

Confirm each of the following by checking the box.  This package:

- [x] does not violate the Terms of Service of any service it interacts with.
- [x] has a CRAN and OSI accepted license.
- [x] contains a README with instructions for installing the development version.
- [x] includes documentation with examples for all functions.
- [x] contains a vignette with examples of its essential functions and uses.
- [x] has a test suite.
- [x] has continuous integration, including reporting of test coverage using services such as Travis CI, Coveralls and/or CodeCov.

## Publication options

- [x] Do you intend for this package to go on CRAN?
- [ ] Do you wish to automatically submit to the [Journal of Open Source Software](http://joss.theoj.org/)? If so:

<details>
 <summary>JOSS Options</summary>

  - [ ] The package has an **obvious research application** according to [JOSS's definition](https://joss.readthedocs.io/en/latest/submitting.html#submission-requirements).
    - [ ] The package contains a `paper.md` matching [JOSS's requirements](https://joss.readthedocs.io/en/latest/submitting.html#what-should-my-paper-contain) with a high-level description in the package root or in `inst/`.
    - [ ] The package is deposited in a long-term repository with the DOI:
    - (*Do not submit your package separately to JOSS*)

</details>

- [ ] Do you wish to submit an Applications Article about your package to [Methods in Ecology and Evolution](http://besjournals.onlinelibrary.wiley.com/hub/journal/10.1111/(ISSN)2041-210X/)? If so:

<details>
<summary>MEE Options</summary>

- [ ] The package is novel and will be of interest to the broad readership of the journal.
- [ ] The manuscript describing the package is no longer than 3000 words.
- [ ] You intend to archive the code for the package in a long-term repository which meets the requirements of the journal (see [MEE's Policy on Publishing Code](http://besjournals.onlinelibrary.wiley.com/hub/journal/10.1111/(ISSN)2041-210X/journal-resources/policy-on-publishing-code.html))
- (*Scope: Do consider MEE's [Aims and Scope](http://besjournals.onlinelibrary.wiley.com/hub/journal/10.1111/(ISSN)2041-210X/aims-and-scope/read-full-aims-and-scope.html) for your manuscript. We make no guarantee that your manuscript will be within MEE scope.*)
- (*Although not required, we strongly recommend having a full manuscript prepared when you submit here.*)
- (*Please do not submit your package separately to Methods in Ecology and Evolution*)
</details>

## Code of conduct

- [x] I agree to abide by [ROpenSci's Code of Conduct](https://github.com/ropensci/onboarding/blob/master/policies.md#code-of-conduct) during the review process and in maintaining my package should it be accepted.
---
output: github_document
---


# virtuoso <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Travis build status](https://travis-ci.org/ropensci/virtuoso.svg?branch=master)](https://travis-ci.org/ropensci/virtuoso)
[![Build status](https://ci.appveyor.com/api/projects/status/x5r18x1cvu6khksd/branch/master?svg=true)](https://ci.appveyor.com/project/cboettig/virtuoso/branch/master)
[![Coverage status](https://codecov.io/gh/ropensci/virtuoso/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/virtuoso?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/virtuoso)](https://cran.r-project.org/package=virtuoso) 
[![Peer review](http://badges.ropensci.org/271_status.svg)](https://github.com/ropensci/software-review/issues/271)

 
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```


The goal of virtuoso is to provide an easy interface to Virtuoso RDF database from R. 

## Installation

You can install the development version of virtuoso from GitHub with:

``` r
remotes::install_github("ropensci/virtuoso")
```

## Getting Started

```{r}
library(virtuoso)
```

For Mac users, `virtuoso` package includes a utility function to install and configure a local Virtuoso Open Source instance using Homebrew.  Otherwise, simply install the Virtuoso Open Source edition for your operating system. 


```{r install}
vos_install()
```

We can now start our Virtuoso server from R:

```{r}
vos_start()
```



 Once the server is running, we can connect to the database.

```{r}
con <- vos_connect()
```

Our connection is now live, and accepts SPARQL queries directly.  

```{r}
DBI::dbGetQuery(con, "SPARQL SELECT * WHERE { ?s ?p ?o } LIMIT 4")
```

## DSL

`virtuoso` also provides wrappers around some common queries to make it easier to work with Virtuoso and RDF.

The bulk loader can be used to quickly import existing sets of triples. 

```{r}
example <- system.file("extdata", "person.nq", package = "virtuoso")
vos_import(con, example)
```

Can also read in compressed formats as well.  Remember to set the pattern match appropriately.  This is convenient because N-Quads compress particularly well, often by a factor of 20 (or rather, can be particularly large when uncompressed, owing to the repeated property and subject URIs).  

```{r}
ex <- system.file("extdata", "library.nq.gz", package = "virtuoso")
vos_import(con, ex)
```

`vos_import` invisibly returns a table of the loaded files, with error message and loading times.  If a file cannot be imported, an error message is returned:

```{r error = TRUE}
bad_file <- system.file("extdata", "bad_quads.nq", package = "virtuoso")
vos_import(con, bad_file)
```

 
We can now query the imported data using SPARQL.  

```{r}
df <- vos_query(con, 
"SELECT ?p ?o 
 WHERE { ?s ?p ?o .
        ?s a <http://schema.org/Person>
       }")
head(df)
```

```{r}
vos_query(con, 
"SELECT ?p ?o 
 WHERE { ?s ?p ?o .
        ?s a <http://example.org/vocab#Chapter>
       }")
```


## Server controls

We can control any `virtuoso` server started with `vos_start()` using a series of helper commands.  

```{r}
vos_status()
```

Advanced usage note: `vos_start()` invisibly returns a `processx` object which we can pass to other server control functions, or access the embedded `processx` control methods directly.  The `virtuoso` package also caches this object in an environment so that it can be accessed directly without having to keep track of an object in the global environment. Use `vos_process()` to return the `processx` object.  For example:

```{r}
library(ps)
p <- vos_process()
ps_is_running(p)
ps_cpu_times(p)
ps_suspend(p)
ps_resume(p)
```

```{r include = FALSE}
vos_kill()
```

## Going further

Please see the package vignettes for more information:

- [details on Virtuoso Installation & configuration](https://docs.ropensci.org/virtuoso/articles/installation.html)
- [The Data Lake: richer examples of RDF use](https://docs.ropensci.org/virtuoso/articles/articles/datalake.html)

---

Please note that the `virtuoso` R package is released with a [Contributor Code of Conduct](https://docs.ropensci.org/virtuoso/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.


```{r include=FALSE}
if(require(codemetar)) codemetar::write_codemeta()
```


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction: Virtuoso Installation and Configuration"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Installation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


[Virtuoso](https://en.wikipedia.org/wiki/Virtuoso_Universal_Server) is a high-performance "universal server" that can act as both a relational database (supporting standard SQL queries) and an RDF triplestore, (supporting SPARQL queries).  
Virtuoso supports communication over the standard ODBC interface, and so R users can potentially connect to Virtuoso merely by installing the server and using the `odbc` R package.  However, installation can present a few gotchas to users unfamiliar with Virtuoso.  This package seeks to streamline the process of installing, managing, and querying a Virtuoso server.  While the package can be also be used merely to provide a standard `DBI` connection to an RDBS, e.g. as a `dplyr` back-end, Virtuoso's popularity and performance is particularly notable with respect to RDF data and SPARQL queries, so most examples focus on those use cases.

## Installation

The `virtuoso` package provides installation helpers for both Mac OSX and Windows users through the function `vos_install()`.  At the time of writing, the Mac OS X installer uses Homebrew to install the Virtuoso Open Source server (similar to the `hugo` installer in RStudio's `blogdown`).  On Windows, `vos_install()` downloads and executes the Windows self-extracting archive (`.exe` file), which will open a standard installation dialog in interactive mode, or be run automatically if not run in an interactive session.  No automated installer is provided for Linux systems; Linux users are encouraged to simply install the appropriate binaries for their distribution (e.g. `apt-get install -y virtuoso-opensource` on Debian/Ubuntu systems). 

## Configuration

Virtuoso Open Source configuration is controlled by a `virtouso.ini` file, which sets, among other things, which directories can be accessed for tasks such as bulk import, as well as performance tweaks such as available memory.  Unfortunately, the Virtuoso server process (`virtuoso-t` application) cannot start without a path to an appropriate config file, and the installers (e.g. on both Windows and Linux) frequently install an example `virtuoso.ini` to a location which can be hard to find and for which users do not have permission to edit directly. Moreover, the file format is not always intuitive to edit.  The `virtuoso` package thus helps locate this file and provides a helper function, `vos_configure()`, to create and modify this configuration file.  Because reasonable defaults are also provided by this function, users should usually not need to call this function manually.  `vos_configure()` is called automatically from `vos_start()` if the path to a `virtuoso.ini` file is not passed to `vos_start()`.  

In addition to configuring Virtuoso's settings through a `virtuoso.ini` file, the other common barrier is setting up the driver for the ODBC connection.  Some installers (Mac, Linux) do not automatically add the appropriate driver to an active `odbcinst.ini` file with a predictable Driver Server Name, which we need to know to initiate the ODBC connection.  An internal helper function handles identifying drivers and establishing the appropriate `odcinst.ini` automatically when necessary. 

## Server management

Lastly, Virtuoso Open Source is often run as a system service, starting when the operating system starts.  This is often undesirable, as the casual laptop user does not want the service running all the time, and can be difficult to control for users unfamiliar with managing such background services on their operating systems.  Instead of this behavior, the `virtuoso` package provides an explicit interface to control the external server. The server only starts when created by `vos_start()`, and ends automatically when the R process ends, or can be killed, paused, or resumed at any time from R (e.g. via `vos_kill()`).  Helper utilities can also query the status and logs of the server from R.  As with most database servers, data persists to disk, at an appropriate location for the OS determined by `rappdirs` package, and a helper utility, `vos_delete_db()` can remove this persistent storage location.  

Users can also connect directly to any existing (local or remote) Virtuoso instance by passing the appropriate information to `vos_connect()`, which can be convenient for queries.  


Note that he Virtuoso back-end provided by the R package `rdflib` can also connect to any Virtuoso server created by the `virtuoso` R package, though queries loading and queries through the `redland` libraries used by `rdflib` will generally be slower than direct calls over ODBC via the `virtuoso` package functions, often dramatically so for large triplestores. 








---
title: "The Data Lake: Schema on Read"
author: "Carl Boettiger"
date: "2020-01-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{datalake}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---






```r
library(virtuoso)
library(dplyr)
library(nycflights13)
library(jsonld)


# needed write_nquads().  Install using: install_github("cboettig/rdftools")
library(rdftools) 
```
Install virtuoso if not already present:


```r
vos_install()
#> Virtuoso is already installed.
```


# Tabular Data

We start up our Virtuoso server, wait for it to come up, and then connect:


```r
vos_start()
#> Virtuoso is already running with pid: 3647
```


```r
con <- vos_connect()
```


We can represent any data as RDF with a little care.  For instance, consider the `nycflights13` data. First, we must represent any primary or foreign keys in any table as URIs, indicated by a prefix, and not by bare strings:


```r
uri_flights <- flights %>% 
  mutate(tailnum = paste0("planes:", tailnum),
         carrier = paste0("airlines:", carrier))
```

We write the `data.frame`s out as nquads.  Recall that each cell of a `data.frame` can be represented as a triple, in which the column is the predicate, the primary key (or row number) the subject, and the cell value the object.  We turn column names and primary keys into URIs using a prefix based on the table name. 


```r
write_nquads(airlines,  "airlines.nq", key = "carrier", prefix = "airlines:")
write_nquads(planes,  "planes.nq", key = "tailnum", prefix = "planes:")
write_nquads(uri_flights,  "flights.nq", prefix = "flights:")
```

We're ready to import all these triples.  This may take a few minutes:



```r
system.time(
  vos_import(con, c("flights.nq", "planes.nq", "airlines.nq"))
)
#>    user  system elapsed 
#>   0.024   0.037 163.698
```


The data from all three tables is now reduced into a single triplestore graph, one triple for each data point. Rather than joining tables, we can write SPARQL query that names the columns we want.




```r
query <- 
'SELECT  ?carrier ?name ?manufacturer ?model ?dep_delay
WHERE {
?flight <flights:tailnum>  ?tailnum .
?flight <flights:carrier>  ?carrier .
?flight <flights:dep_delay>  ?dep_delay .
?tailnum <planes:manufacturer> ?manufacturer .
?tailnum <planes:model> ?model .
?carrier <airlines:name> ?name
}'

system.time(
df <- vos_query(con, query)
)
#>    user  system elapsed 
#>   1.174   0.037   4.393

head(df)
#>       carrier                     name manufacturer     model dep_delay
#> 1 airlines:EV ExpressJet Airlines Inc.      EMBRAER EMB-145XR        -6
#> 2 airlines:EV ExpressJet Airlines Inc.      EMBRAER EMB-145XR        14
#> 3 airlines:EV ExpressJet Airlines Inc.      EMBRAER EMB-145XR         2
#> 4 airlines:EV ExpressJet Airlines Inc.      EMBRAER EMB-145XR        -7
#> 5 airlines:EV ExpressJet Airlines Inc.      EMBRAER EMB-145XR       -10
#> 6 airlines:EV ExpressJet Airlines Inc.      EMBRAER EMB-145XR        -7
```


# List Data




Transform JSON (or list data) into triples.  In this case, we have a large JSON blob (or R list)
containing metadata on all rOpenSci packages:


```r
download.file("https://raw.githubusercontent.com/ropensci/roregistry/gh-pages/raw_cm.json", "raw_cm.json")
nq <- jsonld::jsonld_to_rdf("raw_cm.json") # drops implicit URIs if not base URIs
writeLines(nq, gzfile("ro.nq.gz"))
```



And bulk-import


```r
vos_import(con, "ro.nq.gz")
```

Find all packages where "Carl Boettiger" is an "author", and return:
package name, license, and co-author surnames: 


```r
query <-
"PREFIX schema: <http://schema.org/>
SELECT DISTINCT ?coauthor  ?license ?package 
 WHERE {
 ?s schema:name ?package ;
    schema:author ?author ;
    schema:license ?license ;
    schema:author ?coauth .
 ?author schema:givenName 'Carl' .
 ?author schema:familyName 'Boettiger' .
 ?coauth schema:familyName ?coauthor
}"

vos_query(con, query) %>% distinct() %>%
mutate(license = basename(license), package = basename(package)) # Tidy up URIs into names
#>       coauthor      license                                                                package
#> 1    Boettiger          MIT                               emld: Ecological Metadata as Linked Data
#> 2    Boettiger          MIT      ramlegacy: Download and Read RAM Legacy Stock Assessment Database
#> 3        Gupta          MIT      ramlegacy: Download and Read RAM Legacy Stock Assessment Database
#> 4         Lapp BSD-3-Clause                                               O for the 'NeXML' Format
#> 5          Vos BSD-3-Clause                                               O for the 'NeXML' Format
#> 6    Boettiger BSD-3-Clause                                               O for the 'NeXML' Format
#> 7  Chamberlain BSD-3-Clause                                               O for the 'NeXML' Format
#> 8   Shumelchyk BSD-3-Clause                                               O for the 'NeXML' Format
#> 9    Boettiger          MIT                arkdb: Archive and Unarchive Databases Using Flat Files
#> 10   Boettiger      GPL-3.0                 codemetar: Generate 'CodeMeta' Metadata for R Packages
#> 11      Salmon      GPL-3.0                 codemetar: Generate 'CodeMeta' Metadata for R Packages
#> 12   Boettiger          MIT                 EML: Read and Write Ecological Metadata Language Files
#> 13       Jones          MIT                 EML: Read and Write Ecological Metadata Language Files
#> 14   Boettiger      GPL-3.0                 piggyback: Managing Larger Data on a GitHub Repository
#> 15   Boettiger          MIT                    rdflib: Tools to Manipulate and Query Semantic Data
#> 16   Boettiger          MIT                                  rdryad: Access for Dryad Web Services
#> 17 Chamberlain          MIT                                  rdryad: Access for Dryad Web Services
#> 18         Ram          MIT                                  rdryad: Access for Dryad Web Services
#> 19   Boettiger      CC0-1.0                                rfigshare: An R Interface to 'figshare'
#> 20 Chamberlain      CC0-1.0                                rfigshare: An R Interface to 'figshare'
#> 21        Hart      CC0-1.0                                rfigshare: An R Interface to 'figshare'
#> 22         Ram      CC0-1.0                                rfigshare: An R Interface to 'figshare'
#> 23   Boettiger      CC0-1.0                                   rfishbase: R Interface to 'FishBase'
#> 24 Chamberlain      CC0-1.0                                   rfishbase: R Interface to 'FishBase'
#> 25 Temple Lang      CC0-1.0                                   rfishbase: R Interface to 'FishBase'
#> 26  Wainwright      CC0-1.0                                   rfishbase: R Interface to 'FishBase'
#> 27   Boettiger          MIT                         virtuoso: Interface to 'Virtuoso' using 'ODBC'
#> 28   Boettiger          MIT           datasauce: Create and manipulate Schema.org Dataset metadata
#> 29 Chamberlain          MIT           datasauce: Create and manipulate Schema.org Dataset metadata
#> 30 Chamberlain          MIT                        rcrossref: Client for Various 'CrossRef' 'APIs'
#> 31         Zhu          MIT                        rcrossref: Client for Various 'CrossRef' 'APIs'
#> 32        Jahn          MIT                        rcrossref: Client for Various 'CrossRef' 'APIs'
#> 33   Boettiger          MIT                        rcrossref: Client for Various 'CrossRef' 'APIs'
#> 34         Ram          MIT                        rcrossref: Client for Various 'CrossRef' 'APIs'
#> 35         Ram          MIT      rfisheries: Programmatic Interface to the 'openfisheries.org' API
#> 36   Boettiger          MIT      rfisheries: Programmatic Interface to the 'openfisheries.org' API
#> 37        Dyck          MIT      rfisheries: Programmatic Interface to the 'openfisheries.org' API
#> 38   Boettiger      CC0-1.0          rgpdd: R Interface to the Global Population Dynamics Database
#> 39       Harte      CC0-1.0          rgpdd: R Interface to the Global Population Dynamics Database
#> 40 Chamberlain      CC0-1.0          rgpdd: R Interface to the Global Population Dynamics Database
#> 41         Ram      CC0-1.0          rgpdd: R Interface to the Global Population Dynamics Database
#> 42 Chamberlain          MIT                 rplos: Interface to the Search API for 'PLoS' Journals
#> 43   Boettiger          MIT                 rplos: Interface to the Search API for 'PLoS' Journals
#> 44         Ram          MIT                 rplos: Interface to the Search API for 'PLoS' Journals
#> 45 Chamberlain          MIT                      taxview: Tools for Vizualizing Data Taxonomically
#> 46   Boettiger          MIT                      taxview: Tools for Vizualizing Data Taxonomically
#> 47   Boettiger      CC0-1.0 treebase: Discovery, Access and Manipulation of 'TreeBASE' Phylogenies
#> 48 Temple Lang      CC0-1.0 treebase: Discovery, Access and Manipulation of 'TreeBASE' Phylogenies
```






---
title: "Getting Started"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{getting-started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = !virtuoso:::is_windows()
)
```


The goal of virtuoso is to provide an easy interface to Virtuoso RDF database from R. 


## Installation

You can install the development version of virtuoso from GitHub with:

``` r
remotes::install_github("cboettig/virtuoso")
```

## Getting Started

```{r}
library(virtuoso)
```

For Mac users, `virtuoso` package includes a utility function to install and configure a local Virtuoso Open Source instance using Homebrew.  Otherwise, simply install the Virtuoso Open Source edition for your operating system. 


```{r install}
vos_install()
```

We can now start our Virtuoso server from R:

```{r}
vos_start()
```



 Once the server is running, we can connect to the database.

```{r}
con <- vos_connect()
```

Our connection is now live, and accepts SPARQL queries directly.  

```{r}
DBI::dbGetQuery(con, "SPARQL SELECT * WHERE { ?s ?p ?o } LIMIT 4")
```

## DSL

`virtuoso` also provides wrappers around some common queries to make it easier to work with Virtuoso and RDF.

The bulk loader can be used to quickly import existing sets of triples. 

```{r}
example <- system.file("extdata", "person.nq", package = "virtuoso")
vos_import(con, example)
```

Can also read in compressed formats as well.  Remember to set the pattern match appropriately.  This is convenient because N-Quads compress particularly well, often by a factor of 20 (or rather, can be particularly large when uncompressed, owing to the repeated property and subject URIs).  

```{r}
ex <- system.file("extdata", "library.nq.gz", package = "virtuoso")
vos_import(con, ex)
```

`vos_import` invisibly returns a table of the loaded files, with error message and loading times.  If a file cannot be imported, an error message is returned:

```{r error = TRUE}
bad_file <- system.file("extdata", "bad_quads.nq", package = "virtuoso")
vos_import(con, bad_file)
```

 
We can now query the imported data using SPARQL.  

```{r}
vos_query(con, 
"SELECT ?p ?o 
 WHERE { ?s ?p ?o .
        ?s a <http://schema.org/Person>
       }")
```

```{r}
vos_query(con, 
"SELECT ?p ?o 
 WHERE { ?s ?p ?o .
        ?s a <http://example.org/vocab#Chapter>
       }")
```


## Server controls

We can control any `virtuoso` server started with `vos_start()` using a series of helper commands.  

```{r}
vos_status()
```

Advanced usage note: `vos_start()` invisibly returns a `processx` object which we can pass to other server control functions, or access the embedded `processx` control methods directly.  The `virtuoso` package also caches this object in an environment so that it can be accessed directly without having to keep track of an object in the global environment. Use `vos_process()` to return the `processx` object.  For example:

```{r}
library(ps)
p <- vos_process()
ps_is_running(p)
ps_cpu_times(p)
ps_suspend(p)
ps_resume(p)
```

```{r include = FALSE}
vos_kill()
```

## Going further

Please see the package vignettes for more information:

- [details on Virtuoso Installation & configuration](https://cboettig.github.io/virtuoso/articles/installation.html)
- [The Data Lake: richer examples of RDF use](https://cboettig.github.io/virtuoso/articles/articles/datalake.html)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/virtuoso-package.R
\docType{package}
\name{virtuoso-package}
\alias{virtuoso}
\alias{virtuoso-package}
\title{virtuoso: An R Interface to Virtuoso Using ODBC}
\description{
Virtuoso is a high-performance "universal server," which can act
as both a relational database (supporting standard SQL queries),
and an Resource Description Framework (RDF) triplestore, supporting
SPARQL queries and semantic reasoning. The \code{virtuoso} R package provides
R users with a DBI-compatible connection to the Virtuoso database.
The package also provides helper routines to install, launch, and manage
a Virtuoso server locally on Mac, Windows and Linux platforms using
the standard interactive installers from the R command-line.  By
automatically handling these setup steps, the package can make Virtuoso
considerably faster and easier for a most users to deploy in a local
environment. While this can be used as a normal \code{dplyr} backend, Virtuoso
excels when used as a RDF triplestore.  Managing the bulk import of triples
from common serializations with a single intuitive command is another key
feature of the \code{virtuoso} R package.  Bulk import performance can be tens to
hundreds of times faster than the comparable imports using existing R tools,
including \code{rdflib} and \code{redland} packages.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/virtuoso}
  \item Report bugs at \url{https://github.com/ropensci/virtuoso/issues}
}

}
\author{
\strong{Maintainer}: Carl Boettiger \email{cboettig@gmail.com} (\href{https://orcid.org/0000-0002-1642-628X}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Bryce Mecum \email{brycemecum@gmail.com} (\href{https://orcid.org/0000-0002-0381-3766}{ORCID}) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_process.R
\name{vos_process}
\alias{vos_process}
\title{Return a handle to an existing Virtuoso Process}
\usage{
vos_process(p = NA)
}
\arguments{
\item{p}{a process object, returned by
\code{\link[=vos_process]{vos_process()}} or  \code{\link[=vos_start]{vos_start()}}.
(will be restored from cache if not provided)}
}
\value{
returns the \code{\link[processx:process]{processx::process()}} object cached by \code{\link[=vos_start]{vos_start()}}
to control the external Virtuoso sever process from R.
}
\description{
Generally a user will not need to access this function directly,
though it may be useful for debugging purposes.
}
\examples{
if(has_virtuoso())
vos_process()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_status.R
\name{vos_status}
\alias{vos_status}
\title{Query the server status}
\usage{
vos_status(p = NA, wait = 10)
}
\arguments{
\item{p}{a process object, returned by
\code{\link[=vos_process]{vos_process()}} or  \code{\link[=vos_start]{vos_start()}}.
(will be restored from cache if not provided)}

\item{wait}{number of seconds to wait for server to come online}
}
\value{
a character string indicating the state of the server:
\itemize{
\item "not detected" if no process can be found
\item "dead" process exists but reports that server is not alive.  Server may fail
to come online due to errors in configuration file. see \code{\link[=vos_configure]{vos_configure()}}
\item "running" Server is up and accepting queries.
\item "sleeping" Server is up and accepting queries.
}
}
\description{
Query the server status
}
\details{
Note: Use \code{\link[=vos_log]{vos_log()}} to see the full log
}
\examples{
if(has_virtuoso())
  vos_status()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_install.R
\name{has_virtuoso}
\alias{has_virtuoso}
\title{check for Virtuoso}
\usage{
has_virtuoso()
}
\value{
logical indicating if virtuoso-t binary was found or now.
}
\description{
test if the system has a virtuoso installation on the path
}
\examples{
has_virtuoso()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_query.R
\name{vos_query}
\alias{vos_query}
\title{Run a SPARQL query}
\usage{
vos_query(con, query)
}
\arguments{
\item{con}{a ODBC connection to Virtuoso, from \code{\link[=vos_connect]{vos_connect()}}}

\item{query}{a SPARQL query statement}
}
\value{
a \code{data.frame} containing the results of the query
}
\description{
Run a SPARQL query
}
\details{
SPARQL is a graph query language similar in syntax SQL,
but allows the use of variables to walk through graph nodes.
}
\examples{
vos_status()
\donttest{
if(has_virtuoso()){
vos_start()
con <- vos_connect()

# show first 4 triples in the database
DBI::dbGetQuery(con, "SPARQL SELECT * WHERE { ?s ?p ?o } LIMIT 4")
}
}
}
\references{
\itemize{
\item \url{https://en.wikipedia.org/wiki/SPARQL}
\item \url{https://docs.ropensci.org/rdflib/articles/rdf_intro.html}
}
}
\seealso{
\code{\link[=vos_start]{vos_start()}}, \code{\link[=vos_connect]{vos_connect()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_import.R
\name{vos_import}
\alias{vos_import}
\title{Bulk Import of RDF triples}
\usage{
vos_import(
  con,
  files = NULL,
  wd = ".",
  glob = "*",
  graph = "rdflib",
  n_cores = 1L
)
}
\arguments{
\item{con}{a ODBC connection to Virtuoso, from \code{\link[=vos_connect]{vos_connect()}}}

\item{files}{paths to files to be imported}

\item{wd}{Alternatively, can specify directory and globbing pattern
to import. Note that in this case, wd must be in (or a subdir of)
the \code{AllowedDirs} list of \code{virtuoso.ini} file created by
\code{\link[=vos_configure]{vos_configure()}}. By default, this includes the working directory
where you called \code{\link[=vos_start]{vos_start()}} or \code{\link[=vos_configure]{vos_configure()}}.}

\item{glob}{A wildcard aka globbing pattern (e.g. `"*.nq"``).}

\item{graph}{Name (technically URI) for a graph in the database.
Can leave as default. If a graph is already specified by the
import file (e.g. in nquads), that will be used instead.}

\item{n_cores}{specify the number of available cores for parallel loading.
Particularly useful when importing large numbers of bulk files.}
}
\value{
(Invisibly) returns the status table of the bulk loader,
indicating file loading time or errors.
}
\description{
While triples data can be added one by one over SPARQL queries,
Virtuoso bulk import is by far the fastest way to import large
triplestores in the database.
}
\details{
the bulk importer imports all files matching a pattern
in a given directory.  If given a list of files, these are
temporarily symlinked (or copied on Windows machines) to
the Virtuoso app cache dir in a subdirectory, and the entire
subdirectory is loaded (filtered by the globbing pattern).
If files are not specified, load is called directly on the specified
directory and pattern.  This is particularly useful for loading large
numbers of files.

Note that Virtuoso recommends breaking large files into multiple smaller ones,
which can improve loading time (particularly if using multiple cores.)

Virtuoso Bulk Importer recognizes the following file formats:
\itemize{
\item \code{.grdf}
\item \code{.nq}
\item \code{.owl}
\item \code{.nt}
\item \code{.rdf}
\item \code{.trig}
\item \code{.ttl}
\item \code{.xml}
}

Any of these can optionally be gzipped (with a \code{.gz} extension).
}
\examples{

vos_status()

\donttest{
if(has_virtuoso()){
vos_start()
con <- vos_connect()

example <- system.file("extdata", "person.nq", package = "virtuoso")
vos_import(con, example)
}
}
}
\references{
\url{http://vos.openlinksw.com/owiki/wiki/VOS/VirtBulkRDFLoader}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_start.R
\name{vos_start}
\alias{vos_start}
\title{Start a Virtuoso Server}
\usage{
vos_start(ini = NULL, wait = 30)
}
\arguments{
\item{ini}{path to a virtuoso.ini configuration file. If not
provided, function will attempt to determine the location of the
default configuration file.}

\item{wait}{number of seconds to wait for server to come online}
}
\value{
invisibly returns the \code{\link[processx:process]{processx::process()}} object which can be used
to control the external process from R.  It is not necessary for a user
to store this return object, as \code{\link[=vos_start]{vos_start()}} caches the process object so
it can be automatically accessed by other functions without needing to store
and pass the return object.
}
\description{
This function will attempt to start a virtuoso server
instance that can be managed completely from R.  This allows
the user to easily start, stop, and access server logs and functions
from the R command line.  This server will be automatically shut
down when R exits or restarts, or can be explicitly controlled using
\code{\link[=vos_kill]{vos_kill()}}, \code{\link[=vos_log]{vos_log()}}, and \code{\link[=vos_status]{vos_status()}}.
}
\details{
It can take some time for the server to come up before it is ready to
accept queries.  \code{\link[=vos_start]{vos_start()}} will return as soon as the server is active,
which typically takes about 10 seconds on tested systems. \code{\link[=vos_start]{vos_start()}} monitors
the Virtuoso logs every one second for a maximum time of \code{wait} seconds
(default 30 seconds) to see if the server is ready.  If \code{wait} time is exceeded,
\code{\link[=vos_start]{vos_start()}} will simply return the current server status.  This does not mean
that starting has failed, it may simply need longer before the server is active.
Use \code{\link[=vos_status]{vos_status()}} to continue to monitor the server status manually.

If no \code{virtuoso.ini} configuration file is provided, \code{\link[=vos_start]{vos_start()}} will
automatically attempt to configure one.  For more control over this,
use \code{\link[=vos_configure]{vos_configure()}}, see examples.
}
\examples{
\donttest{

if(has_virtuoso()){
vos_start()
## or with custom config:
vos_start(vos_configure(gigs_ram = 3))

}
}
}
\seealso{
\code{\link[=vos_install]{vos_install()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/virtuoso.R
\name{vos_list_graphs}
\alias{vos_list_graphs}
\title{List graphs}
\usage{
vos_list_graphs(con)
}
\arguments{
\item{con}{a ODBC connection to Virtuoso, from \code{\link[=vos_connect]{vos_connect()}}}
}
\description{
List graphs
}
\examples{
status <- vos_status()
\donttest{
if(has_virtuoso() & is.null(status)){
vos_start()
con <- vos_connect()
vos_list_graphs(con)

}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_connect.R
\name{vos_connect}
\alias{vos_connect}
\title{Connect to a Virtuoso Server over ODBC}
\usage{
vos_connect(
  driver = NULL,
  uid = "dba",
  pwd = "dba",
  host = "localhost",
  port = "1111",
  system_odbcinst = find_odbcinst(),
  local_odbcinst = odbcinst_path()
)
}
\arguments{
\item{driver}{Name of the Driver line in the ODBC configuration}

\item{uid}{User id. Defaults to "dba"}

\item{pwd}{Password. Defaults to "dba"}

\item{host}{IP address of the Virtuoso Server}

\item{port}{Port used by Virtuoso. Defaults to
the Virtuoso standard port, 1111}

\item{system_odbcinst}{Path to the system \code{odbcinst.ini} file. (Does not
require write access.) Default will attempt to find the file for your system.}

\item{local_odbcinst}{Path to the local odbcinst we should use.}
}
\value{
a DBI connection to the Virtuoso database.  This can
be passed to additional virtuoso functions such as \code{\link[=vos_import]{vos_import()}}
or \code{\link[=vos_query]{vos_query()}}, and can also be used as a standard DBI or dplyr
database backend.
}
\description{
Connect to a Virtuoso Server over ODBC
}
\details{
Default parameters are appropriate for the automatic installer
provided by the package and for the default settings typically used by
local Virtuoso installers.  Adjust these only if you are connecting to a
remote virtuoso server that is not controlled from the R package.
}
\examples{
status <- vos_status()
\donttest{
if(has_virtuoso()){
## start up
vos_start()
con <- vos_connect()
}
}
}
\seealso{
\code{\link[=vos_install]{vos_install()}}, \code{\link[=vos_start]{vos_start()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_odbcinst.R
\name{vos_odbcinst}
\alias{vos_odbcinst}
\title{Configure the ODBC Driver for Virtuoso}
\usage{
vos_odbcinst(
  system_odbcinst = find_odbcinst(),
  local_odbcinst = odbcinst_path()
)
}
\arguments{
\item{system_odbcinst}{Path to the system \code{odbcinst.ini} file. (Does not
require write access.) Default will attempt to find the file for your system.}

\item{local_odbcinst}{Path to the local odbcinst we should use.}
}
\value{
the path to the odbcinst file that is created or modified.
}
\description{
ODBC uses an \code{odbcinst.ini} file to point ODBC at the library required
to drive any given database.  This function helps us automatically
locate the driver library on different operating systems and configure
the odbcinst appropriately for each OS.
}
\details{
This function is called automatically by \code{\link[=vos_install]{vos_install()}} and thus
does not usually need to be called by the user.  Users can also manually
configure ODBC as outlined in
\url{https://github.com/r-dbi/odbc#dsn-configuration-files}.
This is merely a convenience function automating that process on most
systems.
}
\examples{
\donttest{
## Configures ODBC and returns silently on success.
vos_odbcinst()

## see where the inst file is located:
inst <- vos_odbcinst()
inst
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_destroy_all.R
\name{vos_destroy_all}
\alias{vos_destroy_all}
\title{Destroy all Virtuoso's directories}
\usage{
vos_destroy_all(force = FALSE)
}
\arguments{
\item{force}{should permissions be changed (if possible) to allow deletion?}
}
\value{
\link{TRUE} if entirely successful in removing all files,
\link{FALSE} otherwise (invisibly).
}
\description{
Provides a clean reset of the system that purges all
data files, config files, cache and log files created
by virtuoso R package. This does not uninstall Virtuoso software
itself, see \code{\link[=vos_uninstall]{vos_uninstall()}} to uninstall.
}
\examples{

\dontshow{
virtuoso:::vos_test_paths()
}
vos_destroy_all()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_kill.R
\name{vos_kill}
\alias{vos_kill}
\title{Stop (kill) the Virtuoso server}
\usage{
vos_kill(p = NA)
}
\arguments{
\item{p}{a process object, returned by
\code{\link[=vos_process]{vos_process()}} or  \code{\link[=vos_start]{vos_start()}}.
(will be restored from cache if not provided)}
}
\description{
Kill ends the process started by \code{\link[=vos_start]{vos_start()}}
}
\details{
vos_kill simply shuts down the local Virtuoso server,
it does not remove any data stored in the database system.
\code{\link[=vos_kill]{vos_kill()}} terminates the process, removing the
process id from the process table.
}
\examples{
\donttest{
if(has_virtuoso()){

  vos_start()
  vos_kill()

  }
}
}
\seealso{
\code{\link[=vos_start]{vos_start()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/virtuoso_paths.R
\name{vos_set_paths}
\alias{vos_set_paths}
\title{set Virtuoso paths}
\usage{
vos_set_paths(
  db_dir = vos_db(),
  config_dir = vos_config(),
  cache_dir = vos_cache(),
  log_dir = vos_logdir(),
  home = virtuoso_home()
)
}
\arguments{
\item{db_dir}{Location of data in the Virtuoso (tables, triplestore)}

\item{config_dir}{Location of configuration files for Virtuoso}

\item{cache_dir}{Location of cache for bulk importing}

\item{log_dir}{Location of Virutoso Server logs}

\item{home}{Location of the Virtuoso installation}
}
\value{
A logical vector, with elements being true
if setting the corresponding variable succeeded
(invisibly).
}
\description{
Set the location of Virtuoso database, configure files,
cache, and logs to your preferred location.  Set home
to the location of your Virtuoso installation.
}
\examples{
if(has_virtuoso())
  vos_set_paths()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_configure.R
\name{vos_configure}
\alias{vos_configure}
\title{Configure Virtuoso Server ini file}
\usage{
vos_configure(
  dirs_allowed = getwd(),
  gigs_ram = 2,
  template = find_virtuoso_ini(),
  db_dir = vos_db()
)
}
\arguments{
\item{dirs_allowed}{Paths (relative or absolute) to directories from which
Virtuoso should have read and write access (e.g. for bulk uploading). Should
be specified as a single comma-separated string.}

\item{gigs_ram}{Indicate approximately the maximum GB of memory Virtuoso can
have access to.  (Used to set NumberOfBuffers & MaxDirtyBuffers in config.)}

\item{template}{Location of an existing virtuoso.ini file which will be used
as a template. By default, \code{vos_configure()} will attempt to locate the
appropriate template for your system.}

\item{db_dir}{location where \code{virtuoso.ini} file should be written.  Other
Virtuoso database log files will also be written here.}
}
\value{
Writes the requested \code{virtuoso.ini} file to the db_dir specified
and returns the path to this file.
}
\description{
Virtuoso Server configuration is determined by a virtuoso.ini file when
server starts. This file includes both system-specific information from
your install (location of server files, addons, etc) and user-configurable
parameters. This helper function provides a way to create and modify an
appropriate \code{virtuoso.ini} file.
}
\examples{
\donttest{
# can take > 5s to test
## configure with typical defaults:
vos_configure()
## Increase or decrease RAM available to virtuoso:
vos_configure(gigs_ram = 1)
}
}
\references{
\url{http://docs.openlinksw.com/virtuoso/dbadm/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_install.R
\name{vos_install}
\alias{vos_install}
\title{Helper method for installing Virtuoso Server}
\usage{
vos_install(ask = is_interactive(), use_brew = FALSE)
}
\arguments{
\item{ask}{Should we ask user for interactive installation?}

\item{use_brew}{Should we use homebrew to install? (MacOS only)}
}
\description{
Installation helper for Mac and Windows machines.  By default,
method will download and launch the official \code{.dmg} or \code{.exe} installer
for your platform, running the standard drag-n-drop installer or
interactive dialog.  Setting \code{ask = FALSE} will allow the installer
to run entirely unsupervised, which is suitable for use in scripts.
Mac users can alternatively opt to install Virtuoso through HomeBrew
by setting \code{use_brew=TRUE}. Linux users should simply install the
\code{virtuoso-opensource} package (e.g. in debian & ubuntu) using the
package manager or by contacting your system administrator.
}
\examples{
\dontshow{ if(has_virtuoso()) }
vos_install()

}
\seealso{
\code{\link[=vos_start]{vos_start()}}, \code{\link[=vos_uninstall]{vos_uninstall()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_log.R
\name{vos_log}
\alias{vos_log}
\title{Query the server logs}
\usage{
vos_log(p = NA, collapse = NULL, just_errors = FALSE)
}
\arguments{
\item{p}{a process object, returned by
\code{\link[=vos_process]{vos_process()}} or  \code{\link[=vos_start]{vos_start()}}.
(will be restored from cache if not provided)}

\item{collapse}{an optional character string to separate the
lines in a single character string.}

\item{just_errors}{logical, default \link{FALSE}. Set to \link{TRUE} to return
just the lines that contain the term "error", which can be useful
in debugging or validating bulk imports.}
}
\value{
Virtuoso logs as a character vector.
}
\description{
Query the server logs
}
\examples{
if(has_virtuoso())
  vos_log()

}
\seealso{
\code{\link[=vos_start]{vos_start()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_uninstall.R
\name{vos_uninstall}
\alias{vos_uninstall}
\title{Uninstall Virtuoso}
\usage{
vos_uninstall()
}
\description{
Automatic uninstaller for Mac OSX and Windows clients.
}
\examples{
\dontrun{
vos_uninstall()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vos_destroy_all.R
\name{vos_delete_db}
\alias{vos_delete_db}
\title{Delete Virtuoso Database}
\usage{
vos_delete_db(ask = is_interactive(), db_dir = vos_db())
}
\arguments{
\item{ask}{ask before deleting?}

\item{db_dir}{location of the directory to delete}
}
\description{
delete the entire Virtuoso database for a fresh start.
}
\examples{

\dontshow{
virtuoso:::vos_test_paths()
}
vos_delete_db()

}
