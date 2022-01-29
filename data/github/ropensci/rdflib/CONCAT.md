---
title: 'Tools to Manipulate and Query Semantic Data'
tags:
 - linked data
 - rdf
 - sparql
 - semantic
 - json-ld
authors:
 - name: Carl Boettiger
   orcid: 0000-0002-1642-628X
   affiliation: 1
affiliations:
 - name: University of California, Berkeley
   index: 1
date: 2017-12-11
bibliography: paper.bib
---

# Summary

The Resource Description Framework, or RDF [@RDF; @W3C_RDF] is a widely used
data representation model that forms the cornerstone of the 
Semantic Web. RDF represents data as a graph rather than 
the familiar data table or rectangle of relational databases.
The `rdflib` package provides a friendly and concise user interface
for performing common tasks on RDF data, such as reading, writing
and converting between the various serializations of RDF data,
including `rdfxml`, `turtle`, `nquads`, `ntriples`, and `json-ld`;
creating new `rdf` graphs, and performing graph queries using SPARQL [@SPARQL; @W3C_SPARQL].
This package wraps the low level `redland` R package [@redland] which
provides direct bindings to the redland C library.  Additionally,
the package supports the newer and more developer friendly
JSON-LD format through the `jsonld` package [@jsonld; @W3C_jsonld].

# References
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

# rdflib <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Travis-CI Build
Status](https://travis-ci.org/ropensci/rdflib.svg?branch=master)](https://travis-ci.org/ropensci/rdflib)
[![Build
status](https://ci.appveyor.com/api/projects/status/n81e9wsh5bh0xrm6?svg=true)](https://ci.appveyor.com/project/cboettig/rdflib)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/rdflib/master.svg)](https://codecov.io/github/ropensci/rdflib?branch=master)
[![CircleCI](https://app.circleci.com/pipelines/github/ropensci/rdflib.svg?style=svg)](https://app.circleci.com/pipelines/github/ropensci/rdflib "Docker tests")
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rdflib)](https://cran.r-project.org/package=rdflib)
[![](http://badges.ropensci.org/169_status.svg)](https://github.com/ropensci/software-review/issues/169)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/rdflib)](https://CRAN.R-project.org/package=rdflib)
[![DOI](https://zenodo.org/badge/100521776.svg)](https://zenodo.org/badge/latestdoi/100521776)

<!-- README.md is generated from README.Rmd. Please edit that file -->

A friendly and consise user interface for performing common tasks on rdf
data, such as parsing and converting between formats including rdfxml,
turtle, nquads, ntriples, and trig, creating rdf graphs, and performing
SPARQL queries. This package wraps the redland R package which provides
direct bindings to the redland C library. Additionally, the package
supports parsing and serialization of rdf into json-ld through the
json-ld package, which binds the official json-ld javascript API. The
package interface takes inspiration from the Python rdflib library.

## Installation

You can install rdflib from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/rdflib")
```

## Basic use

While not required, `rdflib` is designed to play nicely with `%>%`
pipes, so we will load the `magrittr` package as well:

``` r
library(magrittr)
library(rdflib)
```

Parse a file and serialize into a different format:

``` r
system.file("extdata/dc.rdf", package="redland") %>%
  rdf_parse() %>%
  rdf_serialize("test.nquads", "nquads")
```

Perform SPARQL queries:

``` r
sparql <-
 'PREFIX dc: <http://purl.org/dc/elements/1.1/>
  SELECT ?a ?c
  WHERE { ?a dc:creator ?c . }'

system.file("extdata/dc.rdf", package="redland") %>%
rdf_parse() %>%
rdf_query(sparql)
#> Rows: 1 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): a, c
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> # A tibble: 1 × 2
#>   a                      c           
#>   <chr>                  <chr>       
#> 1 http://www.dajobe.org/ Dave Beckett
```

Initialize graph a new object or add triples statements to an existing
graph:

``` r
x <- rdf()
x <- rdf_add(x, 
    subject="http://www.dajobe.org/",
    predicate="http://purl.org/dc/elements/1.1/language",
    object="en")
x
#> Total of 1 triples, stored in hashes
#> -------------------------------
#> <http://www.dajobe.org/> <http://purl.org/dc/elements/1.1/language> "en" .
```

Change the default display format (`nquads`) for graph objects:

``` r
options(rdf_print_format = "jsonld")
x
#> Total of 1 triples, stored in hashes
#> -------------------------------
#> {
#>   "@id": "http://www.dajobe.org/",
#>   "http://purl.org/dc/elements/1.1/language": "en"
#> }
```

## JSON-LD

We can also work with the JSON-LD format through additional functions
provided in the R package, `jsonld`.

``` r
out <- tempfile()
rdf_serialize(x, out, "jsonld")
rdf_parse(out, format = "jsonld")
#> Total of 1 triples, stored in hashes
#> -------------------------------
#> {
#>   "@id": "http://www.dajobe.org/",
#>   "http://purl.org/dc/elements/1.1/language": "en"
#> }
```

For more information on the JSON-LD RDF API, see
<https://json-ld.org/spec/latest/json-ld-rdf/>.

## Advanced Use

See [articles](https://docs.ropensci.org/rdflib/articles/) from the
documentation for advanced use including applications to large
triplestores, example SPARQL queries, and information about additional
database backends.

------------------------------------------------------------------------

## Citing rdflib

Please also cite the underlying `redland` library when citing `rdflib`

Carl Boettiger. (2018). rdflib: A high level wrapper around the redland
package for common rdf applications (Version 0.1.0). Zenodo.
<https://doi.org/10.5281/zenodo.1098478>

Jones M, Slaughter P, Ooms J, Boettiger C, Chamberlain S (2021).
*redland: RDF Library Bindings in R*. doi: 10.5063/F1VM496B (URL:
<https://doi.org/10.5063/F1VM496B>), R package version 1.0.17-15, \<URL:
<https://github.com/ropensci/redland-bindings/tree/master/R/redland>\>.

[![rofooter](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
# rdflib 0.2.4

* bugfix in write_nquads() for rdf method

# rdflib 0.2.3 2020-01-10

* Drop import of deprecated redland method, getNextResult (#33)

# rdflib 0.2.2 2019-01-15

* Minor patch to fix license file
* Updates documentation with hex


# rdflib 0.2.1 2018-11-25

* Minor patch to make test compatible with breaking change in readr 1.2.0 (#30)

# rdflib 0.2.0 2018-11-13

## New Features

* `rdf()` supports all major storage backends: Virtuoso, SQLite, Postgres, MySQL,
   in addition to existing support for BDB and memory-based storage.
* `length()` method added to report length of triplestore
* `print()` method gains `rdf_max_print()` option and does not print huge triplestores
* `print()` method sumarizes total number of triples and backend

# rdflib 0.1.0 (2018-03-02)

## New Features

* `rdf()` supports BDB backend for disk-based storage for large
   triplestores [#6](https://github.com/ropensci/rdflib/issues/6)
* `rdf_parse()` gains an argument `rdf` to append triples to existing graph
* adds `c()` method to concatenate `rdf` objects
* Performance improvements make it possible to handle triplestores with millions of triples
* Two new vignettes better introduce RDF and package functions.

## Minor Improvements

* `rdf_query` now bypasses the very slow iteration over `getNextResult`
   approach and uses an internal redland function call to access all results
   at once in csv format.
* experimental `as_rdf` method now uses a poor-man's nquad serializer to
  rapidly generate rdf (instead of slowly iterating over `add_rdf`).  

* `rdf_add` argument for `object` can now take all atomic types
   (numeric, integer, string, Date, POSIX, logical) and 
   will automatically declare the appropriate `datatype_uri`
   if the user has not manually specified this. 
* Numerous improvements to documentation from rOpenSci onboarding feedback, see 
  [#9](https://github.com/ropensci/rdflib/issues/9) and 
  [#10](https://github.com/ropensci/rdflib/issues/10) 
* both functions and unit tests are broken out into separate files in
  their respective directories.
* Additional example RDF data added in `extdata`
* `rdf_serialize` passes `...` arguments to serializeToFile (e.g. to set a `baseUri`) 

## Bug Fixes 

* `rdf_free()` will also remove the object from the parent frame, 
  reducing the potential for crashing R by referring to a freed pointer.
* fix encoding with UTF-8 characters (coming from nquads & ntriples)
* `rdf_query()` now coerces data into appropriate type 
   if it recognizes the data URI and can match that 
   to an R type (a few XMLSchema types are recognized,
   otherwise still defaults to character string)
* Memory management: All methods free memory from any 
  temporary objects they initialize, tests free memory.
  (e.g. parsers, serializers, query, statement)
* extend unit tests to cover new features, check UTF-8
* `turtle` parser/serializer fixed

## Deprecated

* `trig` support removed (not working in redland without optional
   libraries and alternative compile configuration)


# rdflib 0.0.3 (2018-01-02)

## Bug Fixes

* add paper.md
* add package level documentation
* set base uri when serializing json-ld to rdf ([#5](https://github.com/ropensci/rdflib/issues/5))


# rdflib 0.0.2 (2018-01-02)

## New Features

* Added a `NEWS.md` file to track changes to the package.
* sparql query returns a data.frame format
* added a vignette
* added pkgdown website for vignette

# rdflib 0.0.1 (2017-12-09)

* Initial prototype


Dear CRAN maintainers,

This release provides a minor update as described in the package NEWS.md

## Test environments

* local OS X install, R 3.6.2
* ubuntu 16.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes


