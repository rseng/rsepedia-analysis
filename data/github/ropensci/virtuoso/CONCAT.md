
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
