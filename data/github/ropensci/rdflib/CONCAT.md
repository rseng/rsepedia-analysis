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


---
output: github_document
---

# rdflib <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/rdflib.svg?branch=master)](https://travis-ci.org/ropensci/rdflib) 
[![Build status](https://ci.appveyor.com/api/projects/status/n81e9wsh5bh0xrm6?svg=true)](https://ci.appveyor.com/project/cboettig/rdflib)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/rdflib/master.svg)](https://codecov.io/github/ropensci/rdflib?branch=master)
[![CircleCI](https://app.circleci.com/pipelines/github/ropensci/rdflib.svg?style=svg)](https://app.circleci.com/pipelines/github/ropensci/rdflib "Docker tests")
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rdflib)](https://cran.r-project.org/package=rdflib)
[![](http://badges.ropensci.org/169_status.svg)](https://github.com/ropensci/software-review/issues/169)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/rdflib)](https://CRAN.R-project.org/package=rdflib)
[![DOI](https://zenodo.org/badge/100521776.svg)](https://zenodo.org/badge/latestdoi/100521776)

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```



A friendly and consise user interface for performing common
tasks on rdf data, such as parsing and converting between
formats including rdfxml, turtle, nquads, ntriples,
and trig, creating rdf graphs, and performing SPARQL
queries. This package wraps the redland R package which
provides direct bindings to the redland C library. Additionally,
the package supports parsing and serialization of rdf
into json-ld through the json-ld package, which binds
the official json-ld javascript API. The package
interface takes inspiration from the Python rdflib library.

## Installation

You can install rdflib from GitHub with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/rdflib")
```


## Basic use

While not required, `rdflib` is designed to play nicely with `%>%` pipes, so we will load the `magrittr` package as well:

```{r, message=FALSE}
library(magrittr)
library(rdflib)
```

Parse a file and serialize into a different format:

```{r parse}
system.file("extdata/dc.rdf", package="redland") %>%
  rdf_parse() %>%
  rdf_serialize("test.nquads", "nquads")
```


Perform SPARQL queries:

```{r sparql}
sparql <-
 'PREFIX dc: <http://purl.org/dc/elements/1.1/>
  SELECT ?a ?c
  WHERE { ?a dc:creator ?c . }'

system.file("extdata/dc.rdf", package="redland") %>%
rdf_parse() %>%
rdf_query(sparql)
```

Initialize graph a new object or add triples statements to an existing graph:

```{r}
x <- rdf()
x <- rdf_add(x, 
    subject="http://www.dajobe.org/",
    predicate="http://purl.org/dc/elements/1.1/language",
    object="en")
x
```

Change the default display format (`nquads`) for graph objects:

```{r}
options(rdf_print_format = "jsonld")
x
```


## JSON-LD

We can also work with the JSON-LD format through additional functions provided in the 
R package, `jsonld`. 

```{r}
out <- tempfile()
rdf_serialize(x, out, "jsonld")
rdf_parse(out, format = "jsonld")
```

For more information on the JSON-LD RDF API, see <https://json-ld.org/spec/latest/json-ld-rdf/>.

```{r include=FALSE}
unlink("test.nquads")
unlink(out)
rdf_free(x)
```


## Advanced Use

See [articles](https://docs.ropensci.org/rdflib/articles/) from the documentation for advanced use including applications to large triplestores, example SPARQL queries, and information about additional database backends.  


----

## Citing rdflib


Please also cite the underlying `redland` library when citing `rdflib`

```{r results="asis", warning=FALSE, echo=FALSE}
print(citation("rdflib"), "textVersion")
```

```{r results="asis", warning=FALSE, echo=FALSE}
print(citation("redland"), "text")
```

```{r include=FALSE}
codemeta::write_codemeta()
```


[![rofooter](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
---
title: "A Tidyverse Lover's Introduction to RDF in R"
author:
  - name: Carl Boettiger
    affiliation: UC Berkeley
    address:
    - ESPM Department, University of California,
    - 130 Mulford Hall Berkeley, CA 94720-3114, USA
    - ORCiD \href{https://orcid.org/0000-0002-1642-628X}{0000-0002-1642-628X}
    email:  cboettig@berkeley.edu
abstract: >
  The Resource Description Framework, or RDF, offers a powerful concept of data abstraction and representation which forms the core of the semantic web.  RDF triplestores seek to improve upon the model of relational databases by addressing core challenges such as universal compatibility and extensibility of the underlying data schema, and provide the SPARQL query language to provide graph queries.  However, while R users have long enjoyed rich libraries providing bindings to various SQL-based relational databases such as Postgres and MariaDB, support for RDF-based data has been more limited.  This article reviews several packages that provide support for reading, creating, and querying RDF data from R.  More fundamentally, this article also seeks to introduce users to the potential utility of an RDF-based approach and related concepts.
  
output:
  rticles::rjournal_article:
    includes:
      in_header: preamble.tex
---

## Introduction

Hadley Wickham's classic article, "Tidy Data" \citep{tidydata}, introduced many R users for the first time to well established concepts of relational database design such as Cobb's Third Normal Form, which he dubbed "tidy data," giving rise to the suite of packages now known as the tidyverse \citep{Ross2017}. 


In the world of data science, RDF is a bit of an ugly duckling.  Like XML and Java, only without the massive-adoption-that-refuses-to-die part.  In fact RDF is most frequently expressed in XML, and RDF tools are written in Java, which help give RDF has the aesthetics of *steampunk*, of some technology for some futuristic Semantic Web[^1] in a toolset that feels about as lightweight and modern as iron dreadnought.

[^1]: "The semantic web is the future of the internet and always will be." - Peter Norvig, Director of Research at Google

But don't let these appearances deceive you. After all, SQL-based approaches have been around for decades longer than RDF, and yet SQL technology continues to improve and spread.  


If you've ever gotten carried away using `tidyr::gather` to make everything into one long table, you may have noticed you can just about always get things down to about three columns, as we see with an obligatory `mtcars` data example for `tidyr::gather`:

```{r message = FALSE, warning=FALSE}
library(rdflib)
library(dplyr)
library(tidyr)
library(tibble)
library(jsonld)
```


```{r}
car_triples <- 
mtcars %>% 
  rownames_to_column("Model") %>% 
  gather(attribute,measurement, -Model)
```

```{r echo=FALSE}
DT::datatable(car_triples)
```

If you like long tables like this, RDF is for you. This layout isn't "Tidy Data," where rows are observations and columns are variables, but it is damn useful sometimes.  This format is very liquid, easy to reshape into other structures -- so much so that `tidyr::gather` was originally known as `melt` in the `reshape2` package. It's also a good way to get started thinking about RDF.

## It's all about the triples

Looking at this table closely, we see that *each row is reduced to the most elementary statement you can make from the data*. A row no longer tells you the measurements (observations) *all* attributes (variables) of a given species (key), instead, you get just one fact per row, `Mazda RX4` gets a `mpg` measurement of `21.0`.  In RDF-world, we think of these three-part statements as something very special, which we call **triples**.  RDF is all about these triples. 

The first column came from the row names in this case, the `Model` of car.  This acts serves as a `key` to index the data.frame, i.e. the **subject** being described. The next column is the variable (also called attribute or property) being measured, (that is, column names, other than the key column(s), from the tidy data), called the property or **predicate** in RDF-speak (slash grammar-school jargon). The third column is the actual value measured, more **object** of the predicate. Call it key-property-value or subject-predicate-object, these are our triples.  We can represent just about any data in fully elementary manner.  


&nbsp;          | &nbsp;  | &nbsp;      | &nbsp;
----------------|---------|-------------|-------------
**RDF**         | subject | predicate   | object         
**JSON**        | object  | property    | value
**spreadsheet** | row id  | column name | cell
**data.frame**  | key     | variable    | measurement
**data.frame**  | key     | attribute   | value  

Table: **Table 1**: The many names for triples.  

Table 1 summarizes the many different names associated with triples.  The first naming convention is the terminology typically associated with RDF.  The second set are terms typically associated with JSON data, while the remaining are all examples in tabular or relational data structures.  


## Subject URIs

Using row names as our subject was intuitive but actually a bit sloppy.  `tidyverse` lovers know that `tidyverse` doesn't like rownames, they aren't tidy and have a way of causing trouble.  Of course, we made rownames into a proper column to use `gather`, but we could have taken this one step further.  In true `tidyverse` fashion, this rownames-column is really just one more variable we can observe, one more attribute of the thing we were describing: say, thing A (Car A) is a `car_model_name` as  `Mazda RX4` and thing A also has `mpg` of `21`. We can accomplish such a greater level of abstraction by keeping the Model as just another variable the row ids themselves as the key (i.e. the *subject*) of our triple:

```{r}
car_triples <- 
mtcars %>% 
  rownames_to_column("Model") %>% 
  rowid_to_column("subject") %>% 
  gather(predicate, object, -subject)
```


```{r echo=FALSE}
DT::datatable(car_triples)
```


This is identical to a `gather` of *all* columns, where we have just made the original row ids an explicit column for reference (diligent reader will recognize we would need this information to reverse the operation and `spread` the data back into it's wide form; without it, our transformation is lossy and irreversible).   Our `subject` column now consists only of simple numeric `id`'s, while we have gained an additional triple for every row in the original data which states `Model` of each `id` number (e.g. `1` is `Model` `Mazda RX4`).  Okay, now you're probably thinking: "wait a minute, `1` is not a very unique or specific key, surely that will cause trouble," and you'd be right. For instance, if we performed the same transformation on the iris data, we get triples in the exact same format, ready to `bind_rows`:



```{r}
iris_triples <- iris %>%
  rowid_to_column("subject") %>%
  gather(key = predicate, value = object, -subject)
```

```{r echo=FALSE}
DT::datatable(iris_triples)
```


but in the `iris` data, `1` corresponds to the first individual Iris flower in the measurement data, and not a Mazda RX4.  If we don't want to get confused, we're going to need to make sure our identifiers are unique: not just kind of unique, but unique in the **World** wide.  And what else is unique world-wide? Yup, you guessed it, we are going to use URLs for our subject identifiers, just like the world wide web.  Think of this as a clever out-sourcing to the whole internet domain registry service.  Here, we'll imagine registering each of these example datasets with a separate **base URL**, so instead of a vague `1` to identify the first observation in the `iris` example data, we'll use the URL `http://example.com/iris#1`, which we can now distinguish from `http://example.com/mtcars#1` (and if you're way ahead of me, yes, we'll have more to say about URI vs URL and the use of blank nodes in just a minute).  For example:

```{r}
iris_triples <- iris %>%
  rowid_to_column("subject") %>%
  mutate(subject = paste0("http://example.com/iris#", subject)) %>%
  gather(key = predicate, value = object, -subject)

```


```{r echo=FALSE}
DT::datatable(iris_triples)
```

## Predicate URIs

A slightly more subtle version of the same problem can arise with our predicates. Different tables may use the same attribute (i.e. originally, a column name of a variable) for different things -- the attribute labeled `cyl` means "number of cylinders" in `mtcars` data.frame, but could mean something very different in different data.  Luckily we've already seen how to make names unique in RDF turn them into URLs.

```{r mesage = FALSE}
iris_triples <- iris %>%
  rowid_to_column("subject") %>%
  mutate(subject = paste0("http://example.com/iris#", subject)) %>%
  gather(key = predicate, value = object, -subject) %>%
  mutate(predicate = paste0("http://example.com/iris#", predicate))

```



```{r echo=FALSE}
DT::datatable(iris_triples)
```


At this point the motivation for the name "Linked Data" is probably becoming painfully obvious.  


## Datatype URIs

One more column to go!  But wait a minute, the `object` column is different, isn't it? These measurements don't suffer from the same ambiguity -- after all, there is no confusion if a car has `4` cylinders and an iris has `4` mm long sepals.  However, a new issue has arisen in the data type (e.g. `string`, `boolean`, `double`, `integer`, `dateTime`, etc).  A close look reveals our `object` column is encoded as a `character` and not `numeric` -- how'd that happen?  `tidyr::gather` has coerced the whole column into character strings because some of the values, that is, the `Species` names in `iris` and the Model names in `mtcars`, are text strings (and it couldn't exactly coerce them into integers).  Perhaps this isn't a big deal -- we can often guess the type of an object just by how it looks (so-called [Duck typing](https://en.wikipedia.org/wiki/Duck_typing), because if it quacks like duck...).  Still, being explicit about data types is a Good Thing, so fortunately there's an explicit way to address this too ... oh no ... not ... yes ... more URLs!  

Luckily we don't have to make up `example.com` URLs this time because there's already a well-established list of data types widely used across the internet that were originally developed for use in XML (I warned you) Schemas, listed in see the [W3C RDF DataTypes](https://www.w3.org/TR/rdf11-concepts/#section-Datatypes).  As the standard shows, familiar types `string`, `double`, `boolean`, `integer`, etc are made explicit using the XML Schema URL: `http://www.w3.org/2001/XMLSchema#`, followed by the type; so an integer would be ``http://www.w3.org/2001/XMLSchema#integer`, a character string `http://www.w3.org/2001/XMLSchema#string` etc.  

Because this case is a little different, the URL is attached directly after the object value, which is set off by quotes, using the symbol `^^` (I dunno, but I think two duck feet), such that `5.1` becomes `"5.1"^^http://www.w3.org/2001/XMLSchema#double`.  Wow[^2].  Most of the time we won't have to worry about the type, because, if it quacks... 

[^2]: Couldn't we just have used another column?  Perhaps, but then it wouldn't be a triple.  More to the point, the datatype modifies `object` alone, not the predicate or subject.  


# Triples in `rdflib`

So far, we have explored the concept of triples using familiar `data.frame` structures, but haven't yet introduced any `rdflib` functions.  Though we've been thinking of RDF data in this explicitly tabular three-column structure, that is really just one potentially convenient representation. Just as the same tabular data can be represented in a `data.frame`, written to disk as a `.csv` file, or stored in a database (like MySQL or PostgreSQL), so it is for RDF to even greater degree.  We have separate abstractions for the information itself compared to how it is represented.  

To take advantage of this abstraction, `rdflib` introduces an `rdf` class object. Depending on how this is initialized, this could utilize storage in memory (the default), on disk, or potentially in an array of different databases, (including relational databases like PostgreSQL and rdf-specific ones like Virtuoso, depending on how the underlying `redland` library is compiled -- a topic beyond our scope here).  Here, we simply initialize an `rdf` object using the default in-memory storage:

```{r}
rdf <- rdf()
```

To add triples to this `rdf` object (often called an RDF Model or RDF Graph), we use the function `rdf_add`, which takes a subject, predicate, and object as arguments, as we have just discussed.  A datatype URI can be inferred from the R type used for the object (e.g. `numeric`, `integer`, `logical`, `character`, etc.)


```{r}
base <- "http://example.com/iris#"

rdf %>% 
  rdf_add(subject = paste0(base, "obs1"), 
          predicate = paste0(base, "Sepal.Length"), 
          object = 5.1)

rdf
```

The result is displayed as a triple discussed above. This is technically an example of the `nquad` notation we will see later.  Note the inferred datatype URI.


## Dialing back the ugly

This `gather` thing started well, but now are data is looking pretty ugly, not to mention cumbersome.  You have some idea why RDF hasn't taken data science by storm, and we haven't even looked at how ugly this gets when you write it in the RDF/XML serialization yet!  On the upside, we've introduced most of the essential concepts that will let us start to work with data as triples.  Before we proceed further, we'll take a quick look at some of the options for expressing triples in different ways, and also introduce some of the different serializations (ways of representing in text) frequently used to express these triples.  

#### Prefixes for URIs

Long URL strings are one of the most obvious ways that what started off looking like a concise, minimal statement got ugly and cumbersome.  Borrowing from the notion of [Namespaces in XML](https://en.wikipedia.org/wiki/XML_namespace), most RDF tools permit custom prefixes to be declared and swapped in for longer URLs.  A prefix is typically a short string[^3] followed by a `:` that is used in place of the shared root URL. For instance, we might use the prefix `iris:Sepal.Length` and `iris:Sepal.Width` where `iris:` is defined to mean `http://example.com/iris#` in our example above.   

 [^3]: Technically I believe it should be a [NCName](https://www.w3.org/TR/xmlschema-2/#NCName), defined by the regexp `[\i-[:]][\c-[:]]*`.  [Essentially](https://stackoverflow.com/questions/1631396), this says it cannot include symbol characters like `:`, `@`, `$`, `%`, `&`, `/`, `+`, `,`, `;`, whitespace characters or different parenthesis. Furthermore an NCName cannot begin with a number, dot or minus character although they can appear later in an NCName.  
 
#### URI vs URL

While I've referred to these things as [URL](https://en.wikipedia.org/wiki/URL)s, (uniform resource locator, aka web address) technically they can be a broader class of things known as [URI](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier)s (uniform resource identifier).  In addition to including anything that is a URL, URIs include things which are not URLs, like `urn:isbn:0-486-27557-4` or `urn:uuid:aac06f69-7ec8-403d-ad84-baa549133dce`, which are URNs: unique resource numbers in some numbering scheme (e.g. book ISBN numbers, or UUIDs), neither of which are URLs but nonetheless enjoy the same globally unique property.

#### Blank nodes

Sometimes we do not need a globally unique identifier, we just want a way to refer to a node (e.g. subject, and sometimes an object) uniquely in our document. This is the role of a [blank node](https://en.wikipedia.org/wiki/Blank_node) (do follow the link for a better overview).  These are frequently denoted with the prefix `_:`, e.g. we could have replaced the sample IDs as `_:1`, `_:2` instead of the URLs such as `http://example.com/iris#1` etc.  Note that RDF operations need not preserve the actual string pattern in a blank ID name, it means the exact same thing if we replace all the `_:1`s with `_:b1` and `_:2` with `_:b2`,  etc.

In `librdf` we can get a blank node by passing an empty string or character string that is not a URI as the subject.  Here we also use a URI that isn't a URL as predicate: 

```{r}
rdf <- rdf()
rdf %>% rdf_add("",   
                "iris:Sepal.Length", 
                object = 5.1)
rdf
```

Note that we get a blank node, `_:`  with a randomly generated string.  

## Triple notation: `nquads` `rdfxml`, `turtle`, and `nquads`

So far we have relied primarily on a three-column tabular format to represent our triples.  We have also seen the default `print` format used for the `rdf` method, known as [N-Quads](https://www.w3.org/TR/n-quads/) above, which displays a bare, space-separated triple, possibly with a datatype URI attached to the object.  The line ends with a dot, which indicates this is part of the same local triplestore (aka RDF graph or RDF Model).  Technically this could be another URI indicating a unique global address for the triplestore in question.   

We can serialize any `rdf` object out to a file in this format with the `rdf_serialize()` function, e.g.

```{r}
rdf_serialize(rdf, "rdf.nq", format = "nquads")
```

 Just as each of these formats can be serialized with `rdf_serialize()`, each can be read by `rdflib` using the function `rdf_parse()`:


```{r}
doc <- system.file("extdata/example.rdf", package="redland")
rdf <- rdf_parse(doc, format = "rdfxml") 
rdf
```


N-Quads are convenient in that each triple is displayed on a unique line, and the format supports the blank node and Datatype URIs in the manner we have just discussed. Other formats are not so concise.  Rather than print to file, we can simply change the default print format used by `rdflib` to explore the textual layout of the other serializations.  Here is one of the most common classical serializations, `RDF/XML` which expresses triples in an XML-based schema:

```{r}
options(rdf_print_format = "rdfxml")
rdf
```

Just looking at this is probably enough to explain why so many alternative serializations were created.  Another popular format, [`turtle`](https://www.w3.org/TR/turtle/), looks more like `nquads`: 


```{r}
options(rdf_print_format = "turtle")
rdf
```

Here, blank nodes are denoted by `[]`.  `turtle` uses indentation to indicate that all three predicates (`creator`, `description`, `title`) are properties of the same subject.  



### JSON-LD

While formats such as `nquads` and `turtle` provide a much cleaner syntax than RDF/XML, they also introduce a custom format rather than building on a familiar standard (like XML) for which users already have a well-developed set of tools and intuition.  After more than a decade of such challenges (RDF specification started 1997, including an the HTML-embedded serialization of [RDFa](https://en.wikipedia.org/wiki/RDFa) in 2004), a more user friendly specification has emerged in the form of JSON-LD (1.0 W3C specification was released in 2014, the 1.1 specification released in February 2018).  JSON-LD uses the familiar *object notation* of JSON, (which is rapidly replacing XML as the ubiquitous data exchange format, and will be more familiar to many readers than the specialized RDF formats or even XML.  Here is our `rdf` data in the JSON-LD serialization:


```{r}
options(rdf_print_format = "jsonld")
rdf
```

In this serialization, our subject corresponds to "the thing in the curly braces," (i.e. the JSON "object") which is identified by the special `@id` property (omitting `@id` corresponds to a blank node).  The predicate-object pairs in the triple are then just JSON key-value pairs within the curly braces of the given object.  We can make this format look even more natural by stripping out the URLs.  While it is possible to use prefixes in place of URLs, it is more natural to pull them out entirely, e.g. by declaring a default vocabulary in the JSON-LD "Context", like so:

```{r}
rdf_serialize(rdf, "example.json", "jsonld") %>% 
  jsonld_compact(context = '{"@vocab": "http://purl.org/dc/elements/1.1/"}')
```

The context of a JSON-LD file can also define datatypes, use multiple namespaces, and permit different names in the JSON keys from that found in the URLs.  While a complete introduction to JSON-LD is beyond our scope, this representation essentially provides a way to map intuitive JSON structures into precise RDF triples. 



## From tables to Graphs

So far we have considered examples where the data could be represented in tabular form. 
We frequently encounter data that cannot be easily represented in such a format.  For instance, consider
the JSON data in this example:

```{r}
ex <- system.file("extdata/person.json", package="rdflib")
cat(readLines(ex), sep = "\n")
#jsonld_compact(ex, "{}")

```


This JSON object for a `Person` has another JSON object nested inside (a `PostalAddress`).  Yet if we look at this data as `nquads`, we see the familiar flat triple structure: 



```{r}
options(rdf_print_format = "nquads")
rdf <- rdf_parse(ex, "jsonld")
rdf
```

So what has happened?  Note that our `address` has been given the blank node URI `_:b0`, which serves both as the object in the `address` line of the `Person` and as the subject of all the properties belonging to the `PostalAddress`.  In JSON-LD, this structure is referred to as being 'flattened': 

```{r}
jsonld_flatten(ex, context = "http://schema.org")
```

Note that our JSON-LD structure now starts with an object called `@graph`.  Unlike our opening examples, this data is not tabular in nature, but rather, is formatted as a nested _graph_.  Such nesting is very natural in JSON, where objects can be arranged in a tree-like structure with a single outer-most set of `{}` indicating a root object.  A graph is just a more generic form of a tree structure, where we are agnostic to the root.  (We could in fact use the `@reverse` property on address to create a root `PostalAddress` that contains the `Person`).  In this way, the notion of data as a `graph` offers a powerful generalization to the notion of tabular data.  The `@graph` above consists of two separate objects: a `PostalAddress` (with `id` of `_:b0`) and a `Person` (with an ORCID id).  This layout acts much like a foreign key in a relational database, or as a list-column in `tidyverse` (e.g. see `tidyr::nest()`).  `rdflib` uses this flattened representation when serializing JSON-LD objects.  Note that JSON-LD provides a rich set of utilities to go back and forth between flattened and nested layouts using `jsonld_frame`.  For instance, we can recover the original structure just by specifying a frame that indicates which type we want as the root:

```{r}
jsonld_flatten(ex) %>%
  jsonld_frame('{"@type": "http://schema.org/Person"}') %>%
  jsonld_compact(context = "http://schema.org")
```

(Recall that compacting just replaces URIs and any type declarations with short names given by the context).  This is somewhat analogous to `join` operations in relational data, or nesting and un-nesting functions in `tidyr`.  However, when working with RDF, the beautiful thing is that the differences between these two representations (nested or flattened) are purely aesthetic.  Both representations have precisely the same semantic meaning, and are thus precisely the same thing in RDF world.  We will never have to orchestrate a join on a foreign key before we can perform desired operations like select and filter on the data.  We don't have to think about how our data is organized, because it is always in the same molten triple format, whatever it is, and however nested it might be.  

Just as we saw `gather` could provide a relatively generic way of transforming a data.frame into RDF triples, JSON-LD defines a relatively simple convention for getting nested data (e.g. lists) into RDF triples.  This convention simply treats JSON `{}` objects as `subjects` (often assigning blank node ids, as we saw with row ids), and key-value pairs (or in R-speak, list names and values) as predicates and objects, respectively.  Any raw JSON file can be treated as JSON-LD, ideally by specifying an appropriate `context`, which serves to map terms into URIs as we saw with data.frames. `JSON-LD` is then already a valid RDF format that we can parse with `rdflib`.  

For instance, here is a simple function for coercing list objects into RDF with a specified context:

```{r}
as_rdf.list <- function(x, context = "http://schema.org"){
  if(length(x) == 1) x <- x[[1]]
  x[["@context"]] <- context
  json <- jsonlite::toJSON(x, pretty = TRUE, auto_unbox = TRUE, force = TRUE)
  rdflib::rdf_parse(json, "jsonld")
}
```

Here we set a default context (http://schema.org), and map a few R terms to corresponding schema terms

```{r}
context <- list("http://schema.org", 
                list(schema = "http://schema.org/",
                     given = "givenName",
                     family = "familyName",
                     title = "name",
                     year = "datePublished",
                     note = "softwareVersion",
                     comment = "identifier",
                     role = "http://www.loc.gov/marc/relators/relaterm.html"))

```


We can now apply our function on arbitrary R `list` objects, such as the `bibentry` class object returned by the `citation()` function:

```{r}
options(rdf_print_format = "nquads") # go back to the default


R <- citation("rdflib")
rdf <- as_rdf.list(R, context)
rdf  
```


## SPARQL: A Graph Query Language

```{r}
#source(system.file("examples/as_rdf.R", package="rdflib"))
source(system.file("examples/tidy_schema.R", package="rdflib"))

## Testing: Digest some data.frames into RDF and extract back
 cars <- mtcars %>% rownames_to_column("Model")
 x1 <- as_rdf(iris, NULL, "iris:")
 x2 <- as_rdf(cars, NULL, "mtcars:")
 rdf <- c(x1,x2)
```


## SPARQL: Getting back to Tidy Tables!


```{r}
sparql <-
  'SELECT  ?Species ?Sepal_Length ?Sepal_Width ?Petal_Length  ?Petal_Width
WHERE {
 ?s <iris:Species>  ?Species .
 ?s <iris:Sepal.Width>  ?Sepal_Width .
 ?s <iris:Sepal.Length>  ?Sepal_Length . 
 ?s <iris:Petal.Length>  ?Petal_Length .
 ?s <iris:Petal.Width>  ?Petal_Width 
}'

iris2 <- rdf_query(rdf, sparql)
```

```{r echo=FALSE}
DT::datatable(iris2)
```


We can automatically create the a SPARQL query that returns "tidy data".  Tidy data has predicates as columns, objects as values, subjects as rows.  

```{r}
sparql <- tidy_schema("Species",  "Sepal.Length", "Sepal.Width", prefix = "iris")

rdf_query(rdf, sparql)
```


## Multi-table queries & non-rectangular data

In the following vignette, "The Data Lake and Schema On Read," we explore more complex examples with multiple tables from a relational database and nested JSON data.  These highlight the power of SPARQL as a *Graph* query language. 


## Summary


## Acknowledgements

This manuscript is based on a vignette of the same name which was written for the `rdflib` package.  

\bibliography{RJreferences}

```{r}
library(nycflights13)
library(dplyr)
library(rdflib)
```



```{r}
source(system.file("examples/as_rdf.R", package="rdflib"))
```


We need to set `rdflib` option to use disk-based rather than in-memory storage,
or it appears that `redland` throws an error (even when the machine has sufficient memory!?)
when importing the 336,776 rows of the `flights` table.  

Note: if BDB is not available (e.g. `berkeley-db` libraries were not found when `redland` was built from source),
then this will fallback on in-memory storage and this vignette will use an abridged version of the flights table.

```{r}
options(rdflib_storage = "BDB")
options(rdflib_storage = "memory")
```


## Tidyverse Style

Operations in `dplyr` on the `nyflights13` dataset are easy to write and fast to execute, (in memory or on disk):  

```{r}
df <- flights %>% 
  left_join(airlines) %>%
  left_join(planes, by="tailnum") %>% 
  select(carrier, name, manufacturer, model) %>% 
  distinct()
head(df)
```



Use a smaller dataset if we do not have a BDB backend: 

```{r}
#if(!rdf_has_bdb()){
flights <- flights %>% 
  filter(distance > 3000) # try smaller dataset
#}
```


Keys, including foreign keys, must be represented as URIs and not literal strings.  

```{r}
as_uri <- function(x, base_uri = "x:") paste0(base_uri, x)

uri_flights <- flights %>% 
  mutate(tailnum = as_uri(tailnum),
         carrier = as_uri(carrier))


```



# RDF Serialization Strategies for large data.frames

We consider a variety of strategies for actually importing `data.frames` into RDF:

- **Via `rdf_add()`**: Iterate over each row/cell with calls to `rdf_add`
- **Via JSON-LD**: Coerce the `data.frame` to JSON (via `jsonlite::toJSON(force=TRUE)`), and parse as JSON-LD
- **Via write.table()**` we `tidyr::gather()` and then hack such that we can call `write.table` on a `data.frame` to get an `nquads` text file

As we'll see, only the third solution has adequate performance here. `rdf_add()` requires an initializer call inside each `redland::addStatement`,
which takes a considerable fraction of a second.  Multiply that by the number of cells in the `data.frame` and things do not scale.

`jsonlite` can convert even the large `data.frame`s into JSON reasonably quickly. 
`jsonld::jsonld_to_rdf()` is then also acceptably fast (despite being Javascript) at converting this to `nquads`,
but unfortunately fails dramatically (i.e. `segfault`) when attempting to serialize the flights data.  (Recall we can only get into redland RDF model from JSON-LD via nquads).
Perhaps that is due to some particular data in `flights` table, but it's not obvious.  
Otherwise, this approach has lots to recommend it.  One nice feature about this approach is that it applies to almost 
any R object (e.g. any list object), though some care should be taken with names and URIs, as always. Another nice
feature is that it handles the basic data types automatically -- JSON already has types for logical, double, integer, 
and string, and these will get automatically encoded with the datatype URIs by the built-in `jsonld_to_rdf` algorithm.


The third approach is something of a poor-man's hack to the second approach.  A single call to `rdf_parse()` results in only
a single call through the redland C API to acually read in all the triples -- so unlike the `rdf_add()` approach, all the work
is being done at the C level -- the amount of R code involved doesn't at all depend on the number of triples. This is still
not nearly as fast as reading in large data.frames with `readr` or even with `read.table()`, but is probably as fast as we can get.
The trick then is to serialize the data.frame into an RDF format as quickly as possible.  We can write large `data.frame`s to text
files rather quickly with good ol `write.table()`, and after all `nquads` looks a lot like a space separated, four-column text file,
modulo a little markup to identify URIs and datatypes.  (`readr::write_delim` might be faster, but it's automatic quoting rules appear 
to be incompatible with the `nquads` use of quotations.)  We're left manually encoding the URI strings and the datatypes onto our 
`data.frame` in advance (which requires more nuiance to handle default data types, blank nodes and missing values than I've currently
implemented), but as a proof of principle here this approach is sufficiently fast, as we will now see.  


```{r}
## generic list-based conversion via JSON-LD 
rdf_planes_from_list <- as_rdf.list(planes)
```


Let's do the smaller tables first.  We declare which column is the `key` (i.e. `subject`), 
and we define a `base_uri` prefix which we use to make sure column names and subjects are treated as URIs.
With tables that have only tens of thousands of cells (triples) this is pretty fast:

```{r}
x1 <- as_rdf(airlines, "carrier", "x:")
x2 <- as_rdf(airports, "faa", "x:") ## a few funny chars, UTF8 issues?
x3 <- as_rdf(planes,  "tailnum", "x:")

x <- c(rdf(), x1,x2,x3)
```


SPARQL queries on the resulting data are also pretty fast:




```{r}
sparql <-
  'SELECT   ?model
WHERE {
 ?tailnum <x:carrier> ?carrier .
 ?tailnum <x:model>  ?model 
}'

out <- rdf_query(x1, sparql)
head(out)
```




Big table via poor-man's `nquads` 165 seconds if this is the full table:

```{r}
system.time(
    x4 <- as_rdf(uri_flights, NULL, "x:")
)
```


The json-ld approach just appears to crash, so we won't run that:

```{r}
## nope, the jsonld method appears to crash R...
#system.time(
#  x4 <- as_rdf.list(na.omit(flights))
#)
```


We can join all of these:  

```{r}
rdf <- c(rdf(), x1,x2,x3,x4)
```



Separate queries: This proves very slow on the full data! Would be much faster if we did not have to iterate over getNextResult but could parse all results as a document.  Hopefully this change is coming to `redland` R library soon!

```{r}
sparql <-
'SELECT  ?tailnum ?dep_delay ?carrier
WHERE { 
  ?flight <x:tailnum>  ?tailnum .
  ?flight <x:dep_delay>  ?dep_delay .
  ?flight <x:carrier>  ?carrier 
}'

system.time(

f1 <- rdf_query(rdf, sparql)
)
```

```{r}
sparql <-
  'SELECT  ?tailnum ?model ?manufacturer
WHERE {
?tailnum <x:manufacturer> ?manufacturer .
?tailnum <x:model> ?model
}'
f2 <- rdf_query(rdf, sparql)

tmp <- inner_join(f1,f2)
```




```{r}
s <- 
  'SELECT  ?carrier ?name ?manufacturer ?model
WHERE {
?flight <x:tailnum>  ?tailnum .
?tailnum <x:manufacturer> ?manufacturer .
?tailnum <x:model> ?model .
?flight <x:carrier>  ?carrier .
?carrier <x:name> ?name
}'

out2 <- rdf_query(rdf, s)
head(out2)
```

---
title: "rdflib Introduction"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    df_print: paged
vignette: >
  %\VignetteIndexEntry{rdflib Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---



`rdflib` is really just a lightweight wrapper around two existing R packages: `redland`, and `jsonld`, which themselves are (less trivial) wrappers around existing libraries (the redland C library, and the JSON-LD javascript implementation) which themselves are instances of a set of W3C Standards for the representation of linked data.  `rdflib` has two key features: a simpler, higher-level interface for the common tasks performed in the `redland` library (user does not have to manage `world`, `model` and `storage` objects by default just to perform standard operations and conversions), and integration between the now popular `json-ld` serialization, which is not part of the `redland` library.


```{r setup, message=FALSE}
library(rdflib)
library(magrittr) # pipes

## JSON toolkit 
library(jsonld)
library(jsonlite)
library(jqr)

## For accessing some remote data sources
library(httr)
library(xml2)
library(readr)

## for some typical data processing
library(dplyr)
library(lubridate)

```



## SPARQL queries on JSON-LD data

One example of this utility is the ability to perform graph queries using SPARQL on JSON-LD datasets.  The SPARQL query language is analgous to other query languages such as SQL, but instead of working on an existing set of tables in a relational database, we can query data from a triplestore.  

This is sometimes called "schema on read", since our query will define the schema (i.e. the structure or skeleton of the `data.frame`) our data is returned in.  Because SPARQL is a graph query language, this makes it easy to construct a query that would be cumbersome in SQL or even other recursive tree queries like `jq`, which both would require some knowledge of how the stored data is organized. 

To illustrate this, consider the following example.  NSF asks me to list all of my co-authors within the past four years as conflicts of interest (COI).  Here is a query that for all papers where I am an author, returns a table of given name, family name and year of publication:

```{r}
ex <- system.file("extdata/vita.json", package="rdflib")
vita <- rdf_parse(ex, "jsonld")

sparql <-
 'PREFIX schema: <http://schema.org/>

  SELECT ?coi_given ?coi_family ?year

  WHERE { 
    ?paper a schema:ScholarlyArticle . 
    ?paper schema:author ?authors .
    ?paper schema:dateCreated ?year . 
    ?authors schema:familyName ?coi_family .
    OPTIONAL { ?authors schema:givenName ?coi_given . }

    FILTER ( ?coi_family != "Boettiger" )
}
'

coi <- vita %>% rdf_query(sparql)
```

```{r echo=FALSE}
DT::datatable(coi)
```


Now we have rectangular data, we can tidy things up a bit more. (In principle I believe we could have done this in SPARQL as well)

```{r}
coi2 <- 
  coi %>% 
  ## join names, year as Date
  mutate(year = as.Date(year), name = paste0(coi_family, ", ", coi_given)) %>% 
  ## interaction in last 4 years
  filter(year > lubridate::today() - lubridate::years(3)) %>% 
  ## For each person, list only most recent date
  group_by(name) %>% 
  summarise(year = max(year)) 
  
  
```

```{r echo=FALSE}
DT::datatable(coi2)
```


Yay, that's wasn't so bad.  Rectangling JSON without a graph query is not as easy.  This can be immensely frustrating to do using basic iteration operations with `purrr`, even with it's powerful syntax.  Rectangling is a bit better with a tree-based query language, like a `jq` query.  The only limitation here is that we have to know just a little about how our data is structured, since there are multiple tree structures that correspond to the same graph (or we could just use JSON-LD frame first, but that adds another step in the puzzle.)

Here's the same extraction on the same data, but with a `jq` query:

```{r}
coi_jq <- 
  readr::read_file(ex) %>%
  jqr::jq(
     '."@reverse".author[]  | 
       { year: .dateCreated, 
         author: .author[] | [.givenName, .familyName]  | join(" ")
       }') %>%
  jqr::combine() %>%
  jsonlite::fromJSON()

```

```{r}
coi3 <- 
  coi_jq %>% 
  mutate(year = as.Date(year)) %>%
  filter(year > lubridate::today() - lubridate::years(3)) %>% 
  filter(!grepl("Boettiger", author)) %>%
  ## For each person, list only most recent date
  group_by(author) %>% 
  summarise(year = max(year))  
```

```{r echo=FALSE}
DT::datatable(coi3)
```




## Turning RDF-XML into more friendly JSON


`rdflib` can also be useful as quick and lossless way to convert parse common data formats (e.g. RDF XML) into something more R-friendly.  In this vignette, we illustrate how this might work in a simple example of some citation data returned in RDF-XML format from CrossRef.

Let's begin by reading in some `RDF/XML` data from CrossRef by querying a DOI and requesting `rdf+xml` MIME type (via Content Negotiation):

```{r eval=FALSE}
xml <- "ex.xml"

"https://doi.org/10.1002/ece3.2314" %>%
  httr::GET(httr::add_headers(Accept="application/rdf+xml")) %>%
  httr::content(as = "parsed", type = "application/xml") %>%
  xml2::write_xml(xml)
```

```{r include=FALSE}
xml<- system.file("extdata/ex.xml", package="rdflib")
```

Our `rdflib` functions perform the simple task of parsing this `rdfxml` file into R (as a `redland` `rdf` class object) and then writing it back out in `jsonld` serialization:


```{r}
rdf_parse(xml, "rdfxml") %>% 
  rdf_serialize("ex.json", "jsonld")
```


and we now have JSON file.  We can clean this file up a bit by replacing the long URIs with short prefixes by "compacting" the file into a specific JSON-LD context. FOAF, OWL, and Dublin Core are all recognized by schema.org, so we need not declare them at all here.  PRISM and BIBO ontologies are not, so we simply declare them as additional prefixes:

```{r}
context <- 
'{ "@context": [
    "http://schema.org",
  {
    "prism": "http://prismstandard.org/namespaces/basic/2.1/",
    "bibo": "http://purl.org/ontology/bibo/"
  }]
}'
json <- jsonld_compact("ex.json", context)

```


```{r include=FALSE}
unlink("ex.xml")
unlink("ex.json")
```
---
title: "A tidyverse lover's intro to RDF"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{rdflib Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---

```{r include=FALSE}
options(rdf_print_format = "nquads")
is_linux <- Sys.info()["sysname"] == "Linux"
knitr::opts_chunk$set(eval=is_linux)
```

In the world of data science, RDF is a bit of an ugly duckling.  Like XML and Java, only without the massive-adoption-that-refuses-to-die part.  In fact RDF is most frequently expressed in XML, and RDF tools are written in Java, which help give RDF has the aesthetics of *steampunk*, of some technology for some futuristic Semantic Web[^1] in a toolset that feels about as lightweight and modern as iron dreadnought.

[^1]: "The semantic web is the future of the internet and always will be." -Peter Norvig, Director of Research at Google

But don't let these appearances deceive you. RDF really is cool. If you've ever gotten carried away using `tidyr::gather` to make everything into one long table, you may have noticed you can just about always get things down to about three columns, as we see with an obligatory `mtcars` data example for `tidyr::gather`:

```{r message = FALSE, warning=FALSE}
library(rdflib)
library(dplyr)
library(tidyr)
library(tibble)
library(jsonld)
```


```{r}
car_triples <- 
mtcars %>% 
  rownames_to_column("Model") %>% 
  gather(attribute,measurement, -Model)
```

```{r echo=FALSE}
DT::datatable(car_triples)
```

If you like long tables like this, RDF is for you. This layout isn't "Tidy Data," where rows are observations and columns are variables, but it is damn useful sometimes.  This format is very liquid, easy to reshape into other structures -- so much so that `tidyr::gather` was originally known as `melt` in the `reshape2` package. It's also a good way to get started thinking about RDF.

## It's all about the triples

Looking at this table closely, we see that *each row is reduced to the most elementary statement you can make from the data*. A row no longer tells you the measurements (observations) *all* attributes (variables) of a given species (key), instead, you get just one fact per row, `Mazda RX4` gets a `mpg` measurement of `21.0`.  In RDF-world, we think of these three-part statements as something very special, which we call **triples**.  RDF is all about these triples. 

The first column came from the row names in this case, the `Model` of car.  This acts serves as a `key` to index the data.frame, i.e. the **subject** being described. The next column is the variable (also called attribute or property) being measured, (that is, column names, other than the key column(s), from the tidy data), called the property or **predicate** in RDF-speak (slash grammar-school jargon). The third column is the actual value measured, more **object** of the predicate. Call it key-property-value or subject-predicate-object, these are our triples.  We can represent just about any data in fully elementary manner.  


&nbsp;          | &nbsp;  | &nbsp;      | &nbsp;
----------------|---------|-------------|-------------
**RDF**         | subject | predicate   | object         
**JSON**        | object  | property    | value
**spreadsheet** | row id  | column name | cell
**data.frame**  | key     | variable    | measurement
**data.frame**  | key     | attribute   | value  

Table: **Table 1**: The many names for triples.  

Table 1 summarizes the many different names associated with triples.  The first naming convention is the terminology typically associated with RDF.  The second set are terms typically associated with JSON data, while the remaining are all examples in tabular or relational data structures.  


## Subject URIs

Using row names as our subject was intuitive but actually a bit sloppy.  `tidyverse` lovers know that `tidyverse` doesn't like rownames, they aren't tidy and have a way of causing trouble.  Of course, we made rownames into a proper column to use `gather`, but we could have taken this one step further.  In true `tidyverse` fashion, this rownames-column is really just one more variable we can observe, one more attribute of the thing we were describing: say, thing A (Car A) is a `car_model_name` as  `Mazda RX4` and thing A also has `mpg` of `21`. We can accomplish such a greater level of abstraction by keeping the Model as just another variable the row ids themselves as the key (i.e. the *subject*) of our triple:

```{r}
car_triples <- 
mtcars %>% 
  rownames_to_column("Model") %>% 
  rowid_to_column("subject") %>% 
  gather(predicate, object, -subject)
```


```{r echo=FALSE}
DT::datatable(car_triples)
```


This is identical to a `gather` of *all* columns, where we have just made the original row ids an explicit column for reference (diligent reader will recognize we would need this information to reverse the operation and `spread` the data back into it's wide form; without it, our transformation is lossy and irreversible).   Our `subject` column now consists only of simple numeric `id`'s, while we have gained an additional triple for every row in the original data which states `Model` of each `id` number (e.g. `1` is `Model` `Mazda RX4`).  Okay, now you're probably thinking: "wait a minute, `1` is not a very unique or specific key, surely that will cause trouble," and you'd be right. For instance, if we performed the same transformation on the iris data, we get triples in the exact same format, ready to `bind_rows`:



```{r}
iris_triples <- iris %>%
  rowid_to_column("subject") %>%
  gather(key = predicate, value = object, -subject)
```

```{r echo=FALSE}
DT::datatable(iris_triples)
```


but in the `iris` data, `1` corresponds to the first individual Iris flower in the measurement data, and not a Mazda RX4.  If we don't want to get confused, we're going to need to make sure our identifiers are unique: not just kind of unique, but unique in the **World** wide.  And what else is unique world-wide? Yup, you guessed it, we are going to use URLs for our subject identifiers, just like the world wide web.  Think of this as a clever out-sourcing to the whole internet domain registry service.  Here, we'll imagine registering each of these example datasets with a separate **base URL**, so instead of a vague `1` to identify the first observation in the `iris` example data, we'll use the URL `http://example.com/iris#1`, which we can now distinguish from `http://example.com/mtcars#1` (and if you're way ahead of me, yes, we'll have more to say about URI vs URL and the use of blank nodes in just a minute).  For example:

```{r}
iris_triples <- iris %>%
  rowid_to_column("subject") %>%
  mutate(subject = paste0("http://example.com/", "iris#", subject)) %>%
  gather(key = predicate, value = object, -subject)

```


```{r echo=FALSE}
DT::datatable(iris_triples)
```

## Predicate URIs

A slightly more subtle version of the same problem can arise with our predicates. Different tables may use the same attribute (i.e. originally, a column name of a variable) for different things -- the attribute labeled `cyl` means "number of cylinders" in `mtcars` data.frame, but could mean something very different in different data.  Luckily we've already seen how to make names unique in RDF turn them into URLs.

```{r mesage = FALSE}
iris_triples <- iris %>%
  rowid_to_column("subject") %>%
  mutate(subject = paste0("http://example.com/", "iris#", subject)) %>%
  gather(key = predicate, value = object, -subject) %>%
  mutate(predicate = paste0("http://example.com/", "iris#", predicate))

```



```{r echo=FALSE}
DT::datatable(iris_triples)
```


At this point the motivation for the name "Linked Data" is probably becoming painfully obvious.  


## Datatype URIs

One more column to go!  But wait a minute, the `object` column is different, isn't it? These measurements don't suffer from the same ambiguity -- after all, there is no confusion if a car has `4` cylinders and an iris has `4` mm long sepals.  However, a new issue has arisen in the data type (e.g. `string`, `boolean`, `double`, `integer`, `dateTime`, etc).  A close look reveals our `object` column is encoded as a `character` and not `numeric` -- how'd that happen?  `tidyr::gather` has coerced the whole column into character strings because some of the values, that is, the `Species` names in `iris` and the Model names in `mtcars`, are text strings (and it couldn't exactly coerce them into integers).  Perhaps this isn't a big deal -- we can often guess the type of an object just by how it looks (so-called [Duck typing](https://en.wikipedia.org/wiki/Duck_typing), because if it quacks like duck...).  Still, being explicit about data types is a Good Thing, so fortunately there's an explicit way to address this too ... oh no ... not ... yes ... more URLs!  

Luckily we don't have to make up `example.com` URLs this time because there's already a well-established list of data types widely used across the internet that were originally developed for use in XML (I warned you) Schemas, listed in see the [W3C RDF DataTypes](https://www.w3.org/TR/rdf11-concepts/#section-Datatypes).  As the standard shows, familiar types `string`, `double`, `boolean`, `integer`, etc are made explicit using the XML Schema URL: `http://www.w3.org/2001/XMLSchema#`, followed by the type; so an integer would be ``http://www.w3.org/2001/XMLSchema#integer`, a character string `http://www.w3.org/2001/XMLSchema#string` etc.  

Because this case is a little different, the URL is attached directly after the object value, which is set off by quotes, using the symbol `^^` (I dunno, but I think two duck feet), such that `5.1` becomes `"5.1"^^http://www.w3.org/2001/XMLSchema#double`.  Wow[^2].  Most of the time we won't have to worry about the type, because, if it quacks... 

[^2]: Couldn't we just have used another column?  Perhaps, but then it wouldn't be a triple.  More to the point, the datatype modifies `object` alone, not the predicate or subject.  


# Triples in `rdflib`

So far, we have explored the concept of triples using familiar `data.frame` structures, but haven't yet introduced any `rdflib` functions.  Though we've been thinking of RDF data in this explicitly tabular three-column structure, that is really just one potentially convenient representation. Just as the same tabular data can be represented in a `data.frame`, written to disk as a `.csv` file, or stored in a database (like MySQL or PostgreSQL), so it is for RDF to even greater degree.  We have separate abstractions for the information itself compared to how it is represented.  

To take advantage of this abstraction, `rdflib` introduces an `rdf` class object. Depending on how this is initialized, this could utilize storage in memory (the default), on disk, or potentially in an array of different databases, (including relational databases like PostgreSQL and rdf-specific ones like Virtuoso, depending on how the underlying `redland` library is compiled -- a topic beyond our scope here).  Here, we simply initialize an `rdf` object using the default in-memory storage:

```{r}
rdf <- rdf()
```

To add triples to this `rdf` object (often called an RDF Model or RDF Graph), we use the function `rdf_add`, which takes a subject, predicate, and object as arguments, as we have just discussed.  A datatype URI can be inferred from the R type used for the object (e.g. `numeric`, `integer`, `logical`, `character`, etc.)


```{r}
base <- paste0("http://example.com/", "iris#")

rdf %>% 
  rdf_add(subject = paste0(base, "obs1"), 
          predicate = paste0(base, "Sepal.Length"), 
          object = 5.1)

rdf
```

The result is displayed as a triple discussed above. This is technically an example of the `nquad` notation we will see later.  Note the inferred datatype URI.


## Dialing back the ugly

This `gather` thing started well, but now are data is looking pretty ugly, not to mention cumbersome.  You have some idea why RDF hasn't taken data science by storm, and we haven't even looked at how ugly this gets when you write it in the RDF/XML serialization yet!  On the upside, we've introduced most of the essential concepts that will let us start to work with data as triples.  Before we proceed further, we'll take a quick look at some of the options for expressing triples in different ways, and also introduce some of the different serializations (ways of representing in text) frequently used to express these triples.  

#### Prefixes for URIs

Long URL strings are one of the most obvious ways that what started off looking like a concise, minimal statement got ugly and cumbersome.  Borrowing from the notion of [Namespaces in XML](https://en.wikipedia.org/wiki/XML_namespace), most RDF tools permit custom prefixes to be declared and swapped in for longer URLs.  A prefix is typically a short string[^3] followed by a `:` that is used in place of the shared root URL. For instance, we might use the prefix `iris:Sepal.Length` and `iris:Sepal.Width` where `iris:` is defined to mean `http://example.com/iris#` in our example above.   

 [^3]: Technically I believe it should be a [NCName](https://www.w3.org/TR/xmlschema-2/#NCName), defined by the regexp `[\i-[:]][\c-[:]]*`.  [Essentially](https://stackoverflow.com/questions/1631396), this says it cannot include symbol characters like `:`, `@`, `$`, `%`, `&`, `/`, `+`, `,`, `;`, whitespace characters or different parenthesis. Furthermore an NCName cannot begin with a number, dot or minus character although they can appear later in an NCName.  
 
#### URI vs URL

While I've referred to these things as [URL](https://en.wikipedia.org/wiki/URL)s, (uniform resource locator, aka web address) technically they can be a broader class of things known as [URI](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier)s (uniform resource identifier).  In addition to including anything that is a URL, URIs include things which are not URLs, like `urn:isbn:0-486-27557-4` or `urn:uuid:aac06f69-7ec8-403d-ad84-baa549133dce`, which are URNs: unique resource numbers in some numbering scheme (e.g. book ISBN numbers, or UUIDs), neither of which are URLs but nonetheless enjoy the same globally unique property.

#### Blank nodes

Sometimes we do not need a globally unique identifier, we just want a way to refer to a node (e.g. subject, and sometimes an object) uniquely in our document. This is the role of a [blank node](https://en.wikipedia.org/wiki/Blank_node) (do follow the link for a better overview).  These are frequently denoted with the prefix `_:`, e.g. we could have replaced the sample IDs as `_:1`, `_:2` instead of the URLs such as `http://example.com/iris#1` etc.  Note that RDF operations need not preserve the actual string pattern in a blank ID name, it means the exact same thing if we replace all the `_:1`s with `_:b1` and `_:2` with `_:b2`,  etc.

In `librdf` we can get a blank node by passing an empty string or character string that is not a URI as the subject.  Here we also use a URI that isn't a URL as predicate: 

```{r}
rdf <- rdf()
rdf %>% rdf_add("",   
                "iris:Sepal.Length", 
                object = 5.1)
rdf
```

Note that we get a blank node, `_:`  with a randomly generated string.  

## Triple notation: `nquads` `rdfxml`, `turtle`, and `nquads`

So far we have relied primarily on a three-column tabular format to represent our triples.  We have also seen the default `print` format used for the `rdf` method, known as [N-Quads](https://www.w3.org/TR/n-quads/) above, which displays a bare, space-separated triple, possibly with a datatype URI attached to the object.  The line ends with a dot, which indicates this is part of the same local triplestore (aka RDF graph or RDF Model).  Technically this could be another URI indicating a unique global address for the triplestore in question.   

We can serialize any `rdf` object out to a file in this format with the `rdf_serialize()` function, e.g.

```{r}
rdf_serialize(rdf, "rdf.nq", format = "nquads")
```

 Just as each of these formats can be serialized with `rdf_serialize()`, each can be read by `rdflib` using the function `rdf_parse()`:


```{r}
doc <- system.file("extdata/example.rdf", package="redland")
rdf <- rdf_parse(doc, format = "rdfxml") 
rdf
```


N-Quads are convenient in that each triple is displayed on a unique line, and the format supports the blank node and Datatype URIs in the manner we have just discussed. Other formats are not so concise.  Rather than print to file, we can simply change the default print format used by `rdflib` to explore the textual layout of the other serializations.  Here is one of the most common classical serializations, `RDF/XML` which expresses triples in an XML-based schema:

```{r}
options(rdf_print_format = "rdfxml")
rdf
```

Just looking at this is probably enough to explain why so many alternative serializations were created.  Another popular format, [`turtle`](https://www.w3.org/TR/turtle/), looks more like `nquads`: 


```{r}
options(rdf_print_format = "turtle")
rdf
```

Here, blank nodes are denoted by `[]`.  `turtle` uses indentation to indicate that all three predicates (`creator`, `description`, `title`) are properties of the same subject.  



### JSON-LD

While formats such as `nquads` and `turtle` provide a much cleaner syntax than RDF/XML, they also introduce a custom format rather than building on a familiar standard (like XML) for which users already have a well-developed set of tools and intuition.  After more than a decade of such challenges (RDF specification started 1997, including an the HTML-embedded serialization of [RDFa](https://en.wikipedia.org/wiki/RDFa) in 2004), a more user friendly specification has emerged in the form of JSON-LD (1.0 W3C specification was released in 2014, the 1.1 specification released in February 2018).  JSON-LD uses the familiar *object notation* of JSON, (which is rapidly replacing XML as the ubiquitous data exchange format, and will be more familiar to many readers than the specialized RDF formats or even XML.  Here is our `rdf` data in the JSON-LD serialization:


```{r}
options(rdf_print_format = "jsonld")
rdf
```

In this serialization, our subject corresponds to "the thing in the curly braces," (i.e. the JSON "object") which is identified by the special `@id` property (omitting `@id` corresponds to a blank node).  The predicate-object pairs in the triple are then just JSON key-value pairs within the curly braces of the given object.  We can make this format look even more natural by stripping out the URLs.  While it is possible to use prefixes in place of URLs, it is more natural to pull them out entirely, e.g. by declaring a default vocabulary in the JSON-LD "Context", like so:

```{r}
rdf_serialize(rdf, "example.json", "jsonld") %>% 
  jsonld_compact(context = '{"@vocab": "http://purl.org/dc/elements/1.1/"}')
```

The context of a JSON-LD file can also define datatypes, use multiple namespaces, and permit different names in the JSON keys from that found in the URLs.  While a complete introduction to JSON-LD is beyond our scope, this representation essentially provides a way to map intuitive JSON structures into precise RDF triples. 



## From tables to Graphs

So far we have considered examples where the data could be represented in tabular form. 
We frequently encounter data that cannot be easily represented in such a format.  For instance, consider
the JSON data in this example:

```{r}
ex <- system.file("extdata/person.json", package="rdflib")
cat(readLines(ex), sep = "\n")
#jsonld_compact(ex, "{}")

```


This JSON object for a `Person` has another JSON object nested inside (a `PostalAddress`).  Yet if we look at this data as `nquads`, we see the familiar flat triple structure: 



```{r}
options(rdf_print_format = "nquads")
rdf <- rdf_parse(ex, "jsonld")
rdf
```

So what has happened?  Note that our `address` has been given the blank node URI `_:b0`, which serves both as the object in the `address` line of the `Person` and as the subject of all the properties belonging to the `PostalAddress`.  In JSON-LD, this structure is referred to as being 'flattened': 

```{r}
jsonld_flatten(ex, context = "https://schema.org/")
```

Note that our JSON-LD structure now starts with an object called `@graph`.  Unlike our opening examples, this data is not tabular in nature, but rather, is formatted as a nested _graph_.  Such nesting is very natural in JSON, where objects can be arranged in a tree-like structure with a single outer-most set of `{}` indicating a root object.  A graph is just a more generic form of a tree structure, where we are agnostic to the root.  (We could in fact use the `@reverse` property on address to create a root `PostalAddress` that contains the `Person`).  In this way, the notion of data as a `graph` offers a powerful generalization to the notion of tabular data.  The `@graph` above consists of two separate objects: a `PostalAddress` (with `id` of `_:b0`) and a `Person` (with an ORCID id).  This layout acts much like a foreign key in a relational database, or as a list-column in `tidyverse` (e.g. see `tidyr::nest()`).  `rdflib` uses this flattened representation when serializing JSON-LD objects.  Note that JSON-LD provides a rich set of utilities to go back and forth between flattened and nested layouts using `jsonld_frame`.  For instance, we can recover the original structure just by specifying a frame that indicates which type we want as the root:

```{r}
jsonld_flatten(ex) %>%
  jsonld_frame('{"@type": "https://schema.org//Person"}') %>%
  jsonld_compact(context = "https://schema.org/")
```

(Recall that compacting just replaces URIs and any type declarations with short names given by the context).  This is somewhat analogous to `join` operations in relational data, or nesting and un-nesting functions in `tidyr`.  However, when working with RDF, the beautiful thing is that the differences between these two representations (nested or flattened) are purely aesthetic.  Both representations have precisely the same semantic meaning, and are thus precisely the same thing in RDF world.  We will never have to orchestrate a join on a foreign key before we can perform desired operations like select and filter on the data.  We don't have to think about how our data is organized, because it is always in the same molten triple format, whatever it is, and however nested it might be.  

Just as we saw `gather` could provide a relatively generic way of transforming a data.frame into RDF triples, JSON-LD defines a relatively simple convention for getting nested data (e.g. lists) into RDF triples.  This convention simply treats JSON `{}` objects as `subjects` (often assigning blank node ids, as we saw with row ids), and key-value pairs (or in R-speak, list names and values) as predicates and objects, respectively.  Any raw JSON file can be treated as JSON-LD, ideally by specifying an appropriate `context`, which serves to map terms into URIs as we saw with data.frames. `JSON-LD` is then already a valid RDF format that we can parse with `rdflib`.  

For instance, here is a simple function for coercing list objects into RDF with a specified context:

```{r}
as_rdf.list <- function(x, context = "https://schema.org/"){
  if(length(x) == 1) x <- x[[1]]
  x[["@context"]] <- context
  json <- jsonlite::toJSON(x, pretty = TRUE, auto_unbox = TRUE, force = TRUE)
  rdflib::rdf_parse(json, "jsonld")
}
```

Here we set a default context (https://schema.org/), and map a few R terms to corresponding schema terms

```{r}
context <- list("https://schema.org/", 
                list(schema = "https://schema.org//",
                     given = "givenName",
                     family = "familyName",
                     title = "name",
                     year = "datePublished",
                     note = "softwareVersion",
                     comment = "identifier",
                     role = "http://www.loc.gov/marc/relators/relaterm.html"))

```


We can now apply our function on arbitrary R `list` objects, such as the `bibentry` class object returned by the `citation()` function:

```{r}
options(rdf_print_format = "nquads") # go back to the default


R <- citation("rdflib")
rdf <- as_rdf.list(R, context)
rdf  
```


## SPARQL: A Graph Query Language

So far, we have spent a lot of words describing how to transform data into RDF, and not much actually _doing anything_ cool with said data. 


_Still working on writing this section_

```{r}
#source(system.file("examples/as_rdf.R", package="rdflib"))
source(system.file("examples/tidy_schema.R", package="rdflib"))

## Testing: Digest some data.frames into RDF and extract back
 cars <- mtcars %>% rownames_to_column("Model")
 x1 <- as_rdf(iris, NULL, "iris:")
 x2 <- as_rdf(cars, NULL, "mtcars:")
 rdf <- c(x1,x2)
```


## SPARQL: Getting back to Tidy Tables!


```{r}
sparql <-
  'SELECT  ?Species ?Sepal_Length ?Sepal_Width ?Petal_Length  ?Petal_Width
WHERE {
 ?s <iris:Species>  ?Species .
 ?s <iris:Sepal.Width>  ?Sepal_Width .
 ?s <iris:Sepal.Length>  ?Sepal_Length . 
 ?s <iris:Petal.Length>  ?Petal_Length .
 ?s <iris:Petal.Width>  ?Petal_Width 
}'

iris2 <- rdf_query(rdf, sparql)
```

```{r echo=FALSE}
DT::datatable(iris2)
```


We can automatically create the a SPARQL query that returns "tidy data".  Tidy data has predicates as columns, objects as values, subjects as rows.  

```{r}
sparql <- tidy_schema("Species",  "Sepal.Length", "Sepal.Width", prefix = "iris")

rdf_query(rdf, sparql)
```


```{r include=FALSE}
unlink("rdf.nq")
unlink("example.json")
```
---
title: The Data Lake and Schema On Read
author: Carl Boettiger
date: "2020-01-09"
output: 
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{The Data Lake and Schema On Read}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---





A provocative recent analogy (e.g. [Archer, 2017](https://www.w3.org/blog/2017/04/dive-into-the-semantic-data-lake/)) for thinking about RDF is that of the Data Lake.  Whereas adding new data to a traditional relational database frequently involves laborious wrangling into the an existing rigid data *schema* of tables and columns, RDF can enable a far more simple approach of just 'tossing everything into the data lake.' While this may sound like a recipe for a terrible mess down the road, RDF gives the promise of *schema-on-read*: rather than dictating the shape of your data when first adding it to the database (i.e. schema-on-write), a SPARQL query auto-magically reaches into the lake and retrieves all the relevant data, returning it in exactly the schema you ask for; that is, as a nice single table containing only the requested columns.  

This can serve as a very effective means of data integration (provided a reasonably consistent and diligent use of URIs in identifying subjects and properties (predicates)), since just about any data can be added to the lake without worrying about whether it comes in a schema that matches the existing architecture of the database.  It is this flexibility not to have to define your database schema at the start that is the primary strength of the RDF approach. 

A key advantage of this approach is that it is equally as easy to extract your desired data.frame from non-rectangular data as it is from tables or relational database systems. While the previous vignette focused on simple examples including data from single `data.frames` or small JSON files, this vignette showcases the Data Lake approach on a more complex relational database organized over several tables (`nyflights13` data, used to teach `joins` and other relational data in [Grolemund & Wickham, 2017](http://r4ds.had.co.nz)), and also in a large non-tabular file returned from the GitHub API (e.g. as used in [Bryan & Wickham, 2017 ](https://dcl-2017-04.github.io/curriculum/rectangling.html) to teach data rectangling.) 

Load the necessary libraries to get started: 


```r
## Data 
library(nycflights13)
library(repurrrsive)
```

```
## Error in library(repurrrsive): there is no package called 'repurrrsive'
```

```r
## for comparison approaches
library(dplyr)
library(purrr)
library(jsonlite)

## Our focal package:
library(rdflib)
## experimental functions for rdflib package
source(system.file("examples/tidy_schema.R", package="rdflib"))
```

Configure RDF storage to use the BDB backend for on-disk storage.


```r
rdf <- rdf(storage = "BDB", new_db = TRUE)
```

# Relational data


## The `tidyverse` approach

`tidyverse` operations are incredibly effective for working with relational data.  These `dplyr` on the `nyflights13` dataset are easy to write and fast to execute:


```r
df <- flights %>% 
  left_join(airlines) %>%
  left_join(planes, by="tailnum") %>% 
  select(carrier, name, manufacturer, model) %>% 
  distinct()
head(df)
```

```
## # A tibble: 6 x 4
##   carrier name                   manufacturer     model    
##   <chr>   <chr>                  <chr>            <chr>    
## 1 UA      United Air Lines Inc.  BOEING           737-824  
## 2 AA      American Airlines Inc. BOEING           757-223  
## 3 B6      JetBlue Airways        AIRBUS           A320-232 
## 4 DL      Delta Air Lines Inc.   BOEING           757-232  
## 5 UA      United Air Lines Inc.  BOEING           737-924ER
## 6 B6      JetBlue Airways        AIRBUS INDUSTRIE A320-232
```
 
Still, joins are often a challenge in data preparation.  Tabular formats can often be sloppy about what is a key column and what is a literal value, and also whether a column with the same name in different tables means the same thing in both.  Both of these things pose challenges for later use when joining data.  RDF representation encourages greater discipline through the use of URIs (though we'll run a bit roughshod over that with the caviler use of `x:` here.)  

This example uses only data that is already part of the relational database.  Adding additional information to the database is frequently more tricky, and can result in a rapidly expanding number of tables that can become difficult to work across.
 

## RDF approach

Okay, now let's dump the `nyflights13` into the data lake. Foreign keys in any table must be represented as URIs and not literal strings: 


```r
uri_flights <- flights %>% 
  mutate(tailnum = paste0("planes:", tailnum),
         carrier = paste0("airlines:", carrier))
```

We write the `data.frame`s out as nquads.  Recall that each cell of a `data.frame` can be represented as a triple, in which the column is the predicate, the primary key (or row number) the subject, and the cell value the object.  We turn column names and primary keys into URIs using a prefix based on the table name.  (Note that `rdflib` does this conversion by merely munging cells and calling `write.table`, it is not a standard `redland` library transform).


```r
system.time({
  write_nquads(airlines,  "airlines.nq", key = "carrier", prefix = "airlines:")
  write_nquads(planes,  "planes.nq", key = "tailnum", prefix = "planes:")
  write_nquads(uri_flights,  "flights.nq", prefix = "flights:")
})
```

```
##    user  system elapsed 
##  16.807   2.561  19.128
```

We can now read these into our RDF data lake:


```r
system.time({
  read_nquads("airlines.nq", rdf = rdf)
  read_nquads("flights.nq", rdf = rdf)
  read_nquads("planes.nq", rdf = rdf)

})
```

```
##    user  system elapsed 
##  79.267  51.865  70.052
```

Note that flights does not have a natural key (somewhat surprisingly, `flight` number is not a unique key for this table, as the same flight number is reused on the same route at different times.)  So, we will treat each row as a unique anonymous key by setting the key to `NULL`.

## Schema on read

We simply define the columns we want and we immediately get back the desired `data.frame`:



```r
s <- 
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
df <- rdf_query(rdf, s)
)
```

```
##    user  system elapsed 
##  16.078   1.164  13.734
```

```r
head(df)
```

```
## # A tibble: 6 x 5
##   carrier     name                   manufacturer model     dep_delay
##   <chr>       <chr>                  <chr>        <chr>     <chr>    
## 1 airlines:UA United Air Lines Inc.  BOEING       737-824   2        
## 2 airlines:UA United Air Lines Inc.  BOEING       737-824   4        
## 3 airlines:AA American Airlines Inc. BOEING       757-223   2        
## 4 airlines:B6 JetBlue Airways        AIRBUS       A320-232  -1       
## 5 airlines:DL Delta Air Lines Inc.   BOEING       757-232   -6       
## 6 airlines:UA United Air Lines Inc.  BOEING       737-924ER -4
```



Note that in place of joins, we give more semantically meaningful statements about the data:
e.g. `manufacturer` is a property of a `tailnum` (corresponding to a particular physical aircraft), not of a `flight` number.  Departure delay `dep_delay` is a property of a flight, not of an aircraft (`tailnum`).  

This is reminiscent of  the way in which these data are organized in the relational database tables to begin with: we find `deb_delay` in the `flights` table and `manufacturer` in the `planes` table. Good relational design encourages this, but to work with the data the user is often left having to do the required joins, which also creates tables where these semantics are less clear.  




# Non-tabular data

## The `tidyverse` approach

We start with data from [Bryan & Wickham, 2017 ](https://dcl-2017-04.github.io/curriculum/rectangling.html) lesson on data rectangling:


```r
f <- system.file("extdata/gh_repos.json", package="repurrrsive")
gh_data <- jsonlite::read_json(f)
```

```
## Error in readBin(structure(4L, class = c("file", "connection"), conn_id = <pointer: 0x263>), : can only read from a binary connection
```


The original lesson illustrates the power and reasonably concise syntax of the `tidyverse` package, `purrr`, to iterate over this complex structure to extract the necessary data.  In this approach, nesting of the data is largely a nuisance to overcome rather than an asset to the data analyst:


```r
gh_flat <- gh_data %>% purrr::flatten()  # abandon nested structure and hope we didn't need it
```

```
## Error in eval(lhs, parent, parent): object 'gh_data' not found
```

```r
gh_tibble <- tibble(
  name =     gh_flat %>% map_chr("name"),
  issues =   gh_flat %>% map_int("open_issues_count"),
  wiki =     gh_flat %>% map_lgl("has_wiki"),
  homepage = gh_flat %>% map_chr("homepage", .default = ""),
  owner =    gh_flat %>% map_chr(c("owner", "login"))
)
```

```
## Error in eval(lhs, parent, parent): object 'gh_flat' not found
```

```r
gh_tibble %>% arrange(name) %>% head()
```

```
## Error in eval(lhs, parent, parent): object 'gh_tibble' not found
```

## RDF on non-tabular data

The RDF approach merely treats JSON as JSON-LD within a given vocabulary.  In this context, nesting implicitly provides important semantic information about relationships between the data, which are captured in the RDF triples. Here, we import the JSON data as RDF (and add it to our existing triplestore just for fun)


```r
gh_rdf <- as_rdf(gh_data, rdf = rdf, prefix = "gh:")
```

```
## Error in as_rdf(gh_data, rdf = rdf, prefix = "gh:"): object 'gh_data' not found
```





And we can query it back out of the lake just by selecting the columns of interest. 


```r
s <- 
  'SELECT ?name ?issues ?wiki ?homepage ?owner
WHERE {
?repo <gh:homepage>  ?homepage .
?repo <gh:has_wiki> ?wiki .
?repo <gh:open_issues_count> ?issues .
?repo <gh:name> ?name .
?repo <gh:owner> ?owner_id .
?owner_id <gh:login>  ?owner 
}'

system.time(
rdf_tibble <- rdf_query(rdf, s)
)
```

```
##    user  system elapsed 
##   7.992   0.287   8.279
```

```r
head(rdf_tibble)
```

```
## # A tibble: 0 x 5
## # … with 5 variables: name <chr>, issues <chr>, wiki <chr>, homepage <chr>, owner <chr>
```

## Going further: A `dplyr`-style syntax for SPARQL?

`dplyr` provides a reasonably intuitive, powerful, and concise interface for many common SQL commands.  Indeed, `dplyr` function calls are literally serialized to the corresponding SQL query when working with a relational database back-end (via `dbplyr`).  Can a similar API be developed for SPARQL queries?

SPARQL syntax is obviously inspired by SQL syntax, and includes many of the same operations found in SQL and `dplyr` (e.g. SELECT, WHERE, FILTER, DISTINCT) as well as other more RDF specific queries. 

If we are willing to make some assumptions about the most common queries, we can start to make some simplifying functions.  For instance, in the above patterns, the variables being returned are always objects from the triples, with columns being named using the corresponding predicate.  Assuming some additional convention to define the prefixes and indicate graph traversal (i.e. nested values), we could have constructed the above query from a call:

```r
tidy_schema(name, open_issues_count, has_wiki, homepage owner.login)
```

Though to generalize to arbitrary labels and predicates this might need to support something like:

```r
tidy_schema(list = list(name = "gh:name", 
                        issues = "gh:open_issues_count", 
                        wiki = "gh:has_wiki", 
                        homepage = "gh:homepage",
                        owner = list("gh:owner", login = "gh:login")))
```




---
title: "Working with Database Storage Backends in `rdflib`"
author: "Carl Boettiger"
date: "2020-01-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




By default, `rdflib` stores triples in memory, which can be read in from and written out to files on disk using a variety of serializations (`rdfxml`, `json-ld`, `ntriples`, `nquads` or `turtle`).  This is analgous to reading in a `data.frame` object from a `csv`, `tsv`, or other text-delimited file, which requires the full data object to be read into memory before a user can manipulate it.  This approach runs into limitations with very large datasets, which may exceed the memory available.  Just as R packages such as `dplyr` offer the ability to perform manipulations with very large data through a connection to a database backend such as SQLite, MySQL or PostgreSQL, `rdflib` can rely on these and other databases to provide a disk-based backend. 

### Installation

Unfortunately, at this time, support for these storage devices is not included in the prebuild Mac and Windows binaries of the `redland` R package.  Users wanting to take advantage of disk-based storage must instead build and install the `redland` R package from source; and in some cases, also build the `redland` C libraries from source.  This package provides a [Dockerfile](https://github.com/ropensci/rdflib/blob/master/inst/docker/Dockerfile) containing a portable recipe for building the library with support for the 5 main backend storage devices: SQLite, MySQL, PostgreSQL, Virtuoso, and Berkeley DB.  This vignette documents the installation and use of these devices.  

Also see the [official documentation](http://librdf.org/docs/api/redland-storage-modules.html) for the redland C library discussing all supported storage devices.  


# Benchmarks

Here we show examples of reading in a modest set of 190,000 triples from an `nquads` file and executing a simple SPARQL query using each the five backends supported by `rdflib`.  First we load the libraries and prepare an example file by triple-ifying 10,000 rows of the `flights` data.frame:


```r
library(rdflib)
library(nycflights13)
```


```r
example <- flights[1e4,]

system.time(
write_nquads(example, "example.nq", prefix = "flights")
)
#>    user  system elapsed 
#>   0.009   0.001   0.008

system.time(
write_nquads(flights, "flights.nq", prefix = "flights")
)
#>    user  system elapsed 
#>  15.107   0.474  15.632
```

## In Memory

Because the dataset is small enough to easily fit in memory, the default in-memory option is an obvious choice with excellent performance.  Note that this option will not be possible with larger triplestores (e.g. with millions or more triples).  Our testing has found that even on machines with 100GB+ of RAM, the redland in-memory backend is not able to take advantage of that additional memory and disk-based backends are required. 


```r
triplestore <- rdf()
```


```r

system.time(
  read_nquads("example.nq", rdf = triplestore) # smaller set
)
#>    user  system elapsed 
#>   0.002   0.000   0.002
```


```r
query <- 
  'SELECT  ?carrier ?flight ?origin ?dep_delay
WHERE {
?flight <flights:carrier>  ?carrier .
?flight <flights:dep_delay>  ?dep_delay .
?flight <flights:origin> ?origin
}'
```


```r
system.time(
df <- rdf_query(triplestore, query)
)
#>    user  system elapsed 
#>   0.011   0.008   0.019
```


```r
rdf_free(triplestore)
```


## BDB

The [Berkeley DB](https://en.wikipedia.org/wiki/Berkeley_DB) is a simple key-value store 

> It is the most mature and primary persistent store and suitable for large models, tested in the 2-3 million range.

Berkeley DB is a simple disk-based storage option.  Install both the redland libraries and the berkeley-db (e.g. `bd-dev` on Debian/Ubuntu) libraries, and then install `redland` from source.  BDB is relatively fast for data too large for memory.



```r
triplestore <- rdf(storage="BDB", new_db = TRUE)
system.time(
  read_nquads("flights.nq", rdf = triplestore)
)
#>    user  system elapsed 
#>  55.985  15.686  71.787
```



```r
system.time(
df <- rdf_query(triplestore, query)
)
#>    user  system elapsed 
#>  15.001   0.216  15.222
```


```r
rdf_free(triplestore)
## Becuase BDB is just a hash table, redland needs three separate files:
unlink("rdflib-po2s.db")
unlink("rdflib-so2p.db")
unlink("rdflib-sp2o.db")
```

## Virtuoso

Unlike the other backends which use general purpose relational databases or key-value stores, [Virtuoso](https://en.wikipedia.org/wiki/Virtuoso_Universal_Server) dedicated to RDF-based data.  Virtuoso is a popular open source database for RDF with a rich set of built-in interfaces and features, but we can also interact with it directly through the `redland` bindings just like we do any other backend in `rdflib`.  Virtuoso setup may be slightly more involved for individuals unfamiliar with it, but will probably provide the best performance in the case of very large triplestores.  The example below shows a new database setup, but `rdflib` can also connect to any existing Virtuoso database. 


```r
triplestore <- rdf(storage = "virtuoso", 
                   user = "dba", 
                   password = "dba", 
                   dsn = "Local Virtuoso"
                   )

system.time(
  read_nquads("flights.nq", rdf = triplestore)
)
#>     user   system  elapsed 
#>  111.500   84.856 3336.717
```



```r
system.time(
df <- rdf_query(triplestore, query)
)
#>    user  system elapsed 
#>  17.899   9.056 133.143
df <- rdf_query(triplestore, "SELECT ?s ?p ?o WHERE{ ?s ?p ?o }")
#> Warning: 1391187 parsing failures.
#>     row col expected     actual         file
#> 1338850   o a double _:r1011167 literal data
#> 1338851   o a double _:r1011168 literal data
#> 1338852   o a double _:r1011169 literal data
#> 1338853   o a double _:r1011170 literal data
#> 1338854   o a double _:r1012106 literal data
#> ....... ... ........ .......... ............
#> See problems(...) for more details.
df
#> # A tibble: 6,398,744 x 3
#>    s          p                o
#>    <chr>      <chr>        <dbl>
#>  1 flights:1  flights:year  2013
#>  2 flights:2  flights:year  2013
#>  3 flights:3  flights:year  2013
#>  4 flights:4  flights:year  2013
#>  5 flights:5  flights:year  2013
#>  6 flights:6  flights:year  2013
#>  7 flights:7  flights:year  2013
#>  8 flights:8  flights:year  2013
#>  9 flights:9  flights:year  2013
#> 10 flights:10 flights:year  2013
#> # … with 6,398,734 more rows
```



```r
rdf_free(triplestore)
```

Or on remote virtuoso:



```r
triplestore <- rdf(storage = "virtuoso", 
                   user = "dba", 
                   password = "dba", 
                   host = "virtuoso:1111"
                   )

system.time(
  read_nquads("flights.nq", rdf = triplestore)
)
#>       user     system    elapsed 
#> 129710.695   2036.842 149941.651
```




## POSTGRES

Postgres and MySQL are ubiquitious relational databases.  This backend requires the `redland` binaries are built from source with this support enabled, which is not the case for pre-built Mac or Linux binaries.   


```r
triplestore <- rdf(storage = "postgres", 
                   host = "postgres", 
                   user = "postgres", 
                   password = "rdflib", 
                   new_db = TRUE)
#> Warning in rdf_storage(storage, world, host, port, user, password,
#> database, : postgres driver not found. Falling back on in-memory storage
system.time(
  read_nquads("flights.nq", rdf = triplestore)
)
#>    user  system elapsed 
#>  31.803  20.715  52.586
```




```r
system.time(
df <- rdf_query(triplestore, query)
)
#>    user  system elapsed 
#>   5.449   0.032   5.483
```


```r
rdf_free(triplestore)
```

## MySQL


```r
triplestore <- rdf(storage = "mysql", 
                   host = "mariadb", 
                   user = "root", 
                   password = "rdflib", 
                   database = "mysql",
                   new_db = TRUE
                  )
#> Warning in rdf_storage(storage, world, host, port, user, password,
#> database, : mysql driver not found. Falling back on in-memory storage
  read_nquads("flights.nq", rdf = triplestore)
#> Total of 1996496 triples, stored in hashes
#> -------------------------------
#>  
#>  (preview supressed for performance)
```



```r
system.time(df <- rdf_query(triplestore, query))
#>    user  system elapsed 
#>   5.510   0.004   5.515
```


```r
rdf_free(triplestore)
```




## SQLite


[SQLite](https://en.wikipedia.org/wiki/SQLite) is relatively easy to set up, but appears to have rather poor overall performance. Requires SQLite development libraries installed (should work 'out-of-the-box' with Mac binaries for `redland` package).  



```r
triplestore <- rdf(storage="sqlite", new_db = TRUE, name="rdflib.sqlite")
#> Warning in rdf_storage(storage, world, host, port, user, password,
#> database, : sqlite driver not found. Falling back on in-memory storage

system.time(
  read_nquads("flights.nq", rdf = triplestore)
)
#>    user  system elapsed 
#>  31.607  20.219  51.857
```



```r
system.time(
df <- rdf_query(triplestore, query)
)
#>    user  system elapsed 
#>   5.559   0.028   5.588
```


```r
rdf_free(triplestore)
```

# Building redland with full backend support

Getting full support for the above backend databases through the `redland` R package is not trivial.  The `redland` R package provides bindings to the redland libraries.  Unfortunately, commonly available binary versions of those libraries, such as `librdf0-dev` on Debian, `redland` on Mac OSX brew, and the statically linked versions for Mac and Windows shipping in the R package, do not build those libraries with the optional support for all backends. (NOTE: it is the C library itself which must be built from source with these options, not just the R package source).  Consequently, users must build `librdf` from the original sources, <https://github.com/dajobe/librdf> with all backend linking libraries available, and then also build the `redland` R package from source, to be able to access all of these backends. On a Debian or Ubuntu system this looks like the following:

```bash
apt-get update && apt-get -yq install \
libxml2-dev \
libcurl4-openssl-dev \
libssl-dev \
git \
automake \
libtool \
gtk-doc-tools \
bison \
flex \
libgmp-dev  \
libmhash-dev \
libgcrypt20-dev \
libpcre3-dev \
libv8-dev \
libjq-dev \
libpq-dev \
libdb-dev \
libsqlite3-dev \
libmariadbclient-dev \
librdf-storage-virtuoso \
virtuoso-server \
unixodbc-dev
```

Now we can build `raptor` (parsers), `rasqal` (sparql queries) and `rdflib` from source:

```
git clone git://github.com/dajobe/raptor.git && \
cd raptor && \
./autogen.sh && \
make && \
make install && \
cd .. && \
git clone git://github.com/dajobe/rasqal.git && \
cd rasqal && \
./autogen.sh && \
make && \
make install && \
cd .. && \
git clone git://github.com/dajobe/librdf.git && \
cd librdf && \
./autogen.sh && \
make && \
make install
```

A See the [Dockerfile](https://github.com/ropensci/rdflib/tree/master/inst/docker/Dockerfile) in `inst/docker` for an example of this, or simply use the [Rocker-based](https://rocker-project.org) image `ropensci/rdflib`.  

# Testing

`rdflib` uses Circle-CI to test database backends using `docker-compose`.  `docker-compose` pulls down dedicated docker containers for `postgres` and `mariadb`, along with the `ropensci/rdflib` container, which includes a version of `redland` compiled with support for all five major backend storage systems.  See the [Dockerfile](https://github.com/ropensci/rdflib/tree/master/inst/docker/Dockerfile) in `inst/docker` and associated [docker-compose.yml]((https://github.com/ropensci/rdflib/tree/master/docker-compose.yml) used in testing for an example of this configuration.  You can also pull the Docker image `ropensci/rdflib` from Docker Hub for testing with all these libraries installed.     





The badge below (also on the package README) indicates that these dedicated tests are passing.

[![CircleCI](https://circleci.com/gh/ropensci/rdflib.svg?style=svg)](https://circleci.com/gh/ropensci/rdflib)
---
title: "RDF for Tidyverse Lovers"
subtitle: ""
author: "Carl Boettiger"
date: "2018/03/12"
output:
  xaringan::moon_reader:
    lib_dir: libs
    chakra: libs/remark-latest.min.js
    nature:
      highlightStyle: github
      highlightLines: true
      countIncrementalSlides: false
---



---
class: left, top, inverse
background-image: url(img/uglyduckling.jpg)

# RDF: Ugly Duckling


```{r setup, include=FALSE, message = FALSE}
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(comment = NA)
library(dplyr)
library(tidyr)
library(rdflib)
library(jsonlite)
library(tibble)
options(max.print = 50)

mtcars <- mtcars %>% rownames_to_column("Model")

source(system.file("examples/tidy_schema.R", package="rdflib"))
```


```{r}
cat(readLines(system.file("extdata/ex2.xml", package="rdflib")), sep= "\n")
```

---
class: center, middle, inverse

# The Semantic Web is the Future of the Internet...

---
class: center, middle, inverse

# ... and always will be.

 -- Peter Norvig,  
    Director of Research,  
    Google Inc.


---
class: center, top, inverse
background-image: url(img/steampunk.jpg)

# RDF as Steampunk?


---
class: center, top, inverse
background-image: url(img/tetris.jpg)

# All Data are Tabular

---
class: center, top, inverse
background-image: url(img/tetris-lose.jpg)

# All Data are Definitely Not Tabular



---
class: center, middle, inverse
background-image: url(img/factory-farm.jpg)


# Factory Farm Data

---
class: center, middle, inverse
background-image: url(img/organic-farm.png)


# Or Handcrafted, Organic Data?

---
class: center, middle, inverse

# Heterogeneous Data is Hard



---
class: center, top, inverse
background-image: url(img/field-notes.jpg)

# Heterogeneous Data in Ecology 


---
class: center, top, inverse
background-image: url(img/neon.png)

# Heterogeneous Data in Ecology 


---
class: center, top, inverse
background-image: url(img/integration.png)


# Ecological Metadata Language


---
class: center, middle
background-image: url(img/codemeta.png)

# CodeMeta




---
class: center, top, inverse
background-image: url(img/no-data-lake.jpg)

# The Data Lake

---
class: center, top, inverse
background-image: url(img/data-lake.jpg)

# The Data Lake


---
class: center, middle, inverse

# From Schema on Write

# To Schema on Read

---
class: center, middle, inverse


# All Data Really is Tabular


---
class: center, middle, inverse

# `tidyr::gather()` all the things!


---
class: left, top, inverse

# `tidyr::gather()` all the things!

```{r}
mtcars %>% 
  rowid_to_column("id") %>% 
  gather(property, value, -id)
```


---
class: center, middle, inverse

# Atomizing your data



---
class: left, middle

# Row, Column, Cell

```{r echo=FALSE}
knitr::kable(head(mtcars, 20), "html")
```



---
class: left, middle, inverse

# Object, Property, Value

```{r echo=FALSE}
toJSON(mtcars, pretty = TRUE)
```

---
class: left, middle, inverse

# Subject, Predicate, Object

```{r echo=FALSE}
rdf_ex <- as_rdf(mtcars, prefix = "mtcars:")
rdf_ex
```

---
class: left, top, inverse

# Triples

---
class: center, middle, inverse
background-image: url(img/no-data-lake.jpg)

# Into the Lake: Data Frames


```{r}
triplestore <- rdf()

as_rdf(mtcars, triplestore, "mtcars:")
as_rdf(iris, triplestore, "iris:")

```






---
class: left, middle, inverse
background-image: url(img/no-data-lake.jpg)

# Into the Lake: Lists

Example JSON data returned from the [GitHub API](https://api.github.com/users/cboettig/events)

```{r echo=FALSE}
github.json <- system.file("extdata/github.json", package="rdflib")
cat(readLines(github.json, n = 20), sep="\n")
```




---
class: left, middle, inverse
background-image: url(img/no-data-lake.jpg)

# Into the Lake: Lists

```{r include=FALSE}
events <- read_json(github.json)
```
```{r eval = FALSE}
events <- read_json("https://api.github.com/users/cboettig/events")
```
```{r}
as_rdf(events, triplestore, "gh:")
```


---
class: left, middle

# Schema on read: SPARQL

```{r, results="hide"}
rdf_query(triplestore,
'SELECT  ?Model ?mpg ?cyl ?disp  ?hp
WHERE {
 ?s <mtcars:Model>  ?Model ;
    <mtcars:mpg>  ?mpg ;
    <mtcars:cyl>  ?cyl ; 
    <mtcars:disp>  ?disp ;
    <mtcars:hp>  ?hp 
}')

```

---
class: left, middle

# Schema on read: SPARQL

```{r, echo=FALSE}
rdf_query(triplestore,
'SELECT  ?Model ?mpg ?cyl ?disp  ?hp
WHERE {
 ?s <mtcars:Model>  ?Model ;
    <mtcars:mpg>  ?mpg ;
    <mtcars:cyl>  ?cyl ; 
    <mtcars:disp>  ?disp ;
    <mtcars:hp>  ?hp 
}')

```



---
class: left, middle

# Data Rectangling


```{r, results="hide"}
rdf_query(triplestore, 
'SELECT ?type ?user ?repo ?when
WHERE {
?s <gh:type> ?type ;
   <gh:created_at> ?when ;
   <gh:repo> ?repo_id ;
   <gh:actor> ?actor .
?actor <gh:login> ?user .
?repo_id <gh:name> ?repo
}')
```


---
class: left, middle

# Data Rectangling


```{r, echo=FALSE}
rdf_query(triplestore, 
'SELECT ?type ?user ?repo ?when
WHERE {
?s <gh:type> ?type ;
   <gh:created_at> ?when ;
   <gh:repo> ?r ;
   <gh:actor> ?actor .
?r <gh:name> ?repo .
?actor <gh:login> ?user .
}')
```

---
class: left, middle

# Data Rectangling: Graph Queries


```{r }
df <- rdf_query(triplestore, 
'SELECT DISTINCT ?property ?value
WHERE {
?s <gh:url> "https://api.github.com/repos/cboettig/noise-phenomena" .
?parent ?p ?s .
?parent ?property ?value
}')

```

---
class: left, middle

# Data Rectangling: Graph Queries

```{r echo=FALSE}
df
```


---
class: center, middle, inverse

# Potential Issues

- Potential column name or property collisions
- Dealing with data types (numeric, logical, dates, etc)
- Potential row name or object id collisions

---
class: center, middle, inverse

# U say URL
# I say IRI 


---
class: left, middle, inverse

# Internationalized Resource Identifiers


- `https://example.com`
- `https://schema.org/givenName`
- `isbn:978-0-387-98140-6`
- `urn:uuid:0aae8482-93b9-4b22-879e-aa71af0d3fd1`


---
class: center, middle, inverse

# Unique variable/column names

- `https://schema.org/givenName`
- `https://schema.org/programmingLanguage`
- `https://schema.org/softwareRequirements`



---
class: center, middle, inverse

# Data types

- `http://www.w3.org/2001/XMLSchema#decimal`
- `http://www.w3.org/2001/XMLSchema#dateTime`


---
class: center, middle, inverse

# Subject IRIs

- `https://example.com`
- `isbn:978-0-387-98140-6`
- `doi:10.1007/978-0-387-98141-3`
- `_:` Blank nodes

---
class: center, middle, inverse

# Object Types and Resource Nodes

- `https://schema.org/SoftwareSoureCode`
- `doi:10.1007/978-0-387-98141-3`


---
class: center, middle, inverse

# Practical Issues

- `SQL` -> `dplyr`
- `SPARQL` -> `???`




---
class: center, middle, inverse


# Explore & Contribute

- <https://github.com/ropensci/rdflib>
- <https://cran.r-project.org/package=rdflib>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf_has_bdb.R
\name{rdf_has_bdb}
\alias{rdf_has_bdb}
\title{Check for BDB support}
\usage{
rdf_has_bdb()
}
\value{
TRUE if BDB support is detected, false otherwise
}
\description{
Detect whether Berkeley Database for disk-based storage of RDF graphs
is available.  Disk-based storage requires redland package
to be installed from source with support for the Berkeley DB
(libdb-dev on Ubuntu, berkeley-db on homebrew), otherwise \code{rdf()} will
fall back to in-memory storage with a warning.
}
\examples{
rdf_has_bdb()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_rdf.R
\name{as_rdf}
\alias{as_rdf}
\title{Coerce an object into RDF}
\usage{
as_rdf(
  x,
  rdf = NULL,
  prefix = NULL,
  base = getOption("rdf_base_uri", "localhost://"),
  context = NULL,
  key_column = NULL
)
}
\arguments{
\item{x}{an object to coerce into RDF (list, list-like, or data.frame)}

\item{rdf}{An existing rdf object, (by default a new object will be initialized)}

\item{prefix}{A default vocabulary (URI prefix) to assume for all predicates}

\item{base}{A base URI to assume for blank subject nodes}

\item{context}{a named list mapping any string to a URI}

\item{key_column}{name of a column which should be
treated as the primary key in a table. must be unique}
}
\description{
Coerce an object into RDF
}
\examples{
as_rdf(mtcars)
as_rdf(list(repo = "rdflib", owner = list("id", "ropensci")))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf_add.R
\name{rdf_add}
\alias{rdf_add}
\title{Add RDF Triples}
\usage{
rdf_add(
  rdf,
  subject,
  predicate,
  object,
  subjectType = as.character(NA),
  objectType = as.character(NA),
  datatype_uri = as.character(NA)
)
}
\arguments{
\item{rdf}{an rdf object}

\item{subject}{character string containing the subject}

\item{predicate}{character string containing the predicate}

\item{object}{character string containing the object}

\item{subjectType}{the Node type of the subject, i.e. "uri", "blank"}

\item{objectType}{the Node type of the object, i.e. "literal", "uri", "blank"}

\item{datatype_uri}{the datatype URI to associate with a object literal value}
}
\value{
Silently returns the updated RDF graph (rdf object).
Since the rdf object simply contains external pointers
to the model object in C code, note that the input object is modified
directly, so you need not assign the output of rdf_add() to anything.
}
\description{
add a triple (subject, predicate, object) to the RDF graph
}
\details{
\code{rdf_add()} will automatically 'duck type' nodes (if looks like a duck...).
That is, strings that look like URIs will be declared as URIs. (See
\href{https://en.wikipedia.org/wiki/Uniform_Resource_Identifier}{URI}).
Predicate should always be a URI (e.g. URL or  a \code{prefix:string}),
cannot be blank or literal.  Subjects that look like strings will be
treated as \href{https://en.wikipedia.org/wiki/Blank_node}{Blank Nodes} (i.e.
will be prefixed with \verb{_:}).  An empty subject, \code{""}, will create a
blank node with random name.  Objects that look like URIs will be
typed as resource nodes, otherwise as literals.  An empty object \code{""}
will be treated as blank node.  Set \code{subjectType} or \code{objectType}
explicitly to override this behavior, e.g. to treat an object URI
as a literal string.  NAs are also treated as blank nodes in subject
or object  See examples for details.
}
\examples{
rdf <- rdf()
rdf_add(rdf, 
    subject="http://www.dajobe.org/",
    predicate="http://purl.org/dc/elements/1.1/language",
    object="en")
    
## non-URI string in subject indicates a blank subject
## (prefixes to "_:b0")
rdf_add(rdf, "b0", "http://schema.org/jobTitle", "Professor") 

## identically a blank subject.  
## Note rdf is unchanged when we add the same triple twice.
rdf_add(rdf, "b0", "http://schema.org/jobTitle", "Professor", 
        subjectType = "blank") 
        
## blank node with empty string creates a default blank node id
rdf_add(rdf, "", "http://schema.org/jobTitle", "Professor")   
                    

## Subject and Object both recognized as URI resources:
rdf_add(rdf, 
        "https://orcid.org/0000-0002-1642-628X",
        "http://schema.org/homepage", 
        "http://carlboettiger.info")  

 ## Force object to be literal, not URI resource        
rdf_add(rdf, 
        "https://orcid.org/0000-0002-1642-628X",
        "http://schema.org/homepage", 
        "http://carlboettiger.info",
        objectType = "literal")  
        

}
\references{
\url{https://en.wikipedia.org/wiki/Uniform_Resource_Identifier}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf_free.R
\name{rdf_free}
\alias{rdf_free}
\title{Free Memory Associated with RDF object}
\usage{
rdf_free(rdf, rm = TRUE)
}
\arguments{
\item{rdf}{an rdf object}

\item{rm}{logical, default TRUE. Remove pointer from parent.frame()?
Usually a good idea since referring to a pointer after it has been
removed can crash R.}
}
\description{
Free Memory Associated with RDF object
}
\details{
Free all pointers associated with an rdf object.
Frees memory associated with the storage, world, and model
objects.
}
\examples{
rdf <- rdf()
rdf_free(rdf)
rm(rdf)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf.R
\docType{package}
\name{rdflib-package}
\alias{rdflib}
\alias{rdflib-package}
\title{rdflib: Tools to Manipulate and Query Semantic Data}
\description{
The Resource Description Framework, or RDF is a widely used
data representation model that forms the cornerstone of the
Semantic Web. 'RDF' represents data as a graph rather than
the familiar data table or rectangle of relational databases.
}
\details{
It has three main goals:

\itemize{
\item Easily read, write, and convert between all major RDF serialization formats
\item Support SPARQL queries to extract data from an RDF graph into a data.frame
\item Support JSON-LD format as a first-class citizen in RDF manipulations
}

For more information, see the Wikipedia pages for RDF, SPARQL, and JSON-LD:

\itemize{
\item \url{https://en.wikipedia.org/wiki/Resource_Description_Framework}
\item \url{https://en.wikipedia.org/wiki/SPARQL}
\item \url{https://en.wikipedia.org/wiki/JSON-LD}
}

To learn more about rdflib, start with the vignettes:
\code{browseVignettes(package = "rdflib")}

Configurations via \code{options()}

\code{rdf_print_format}:
\itemize{
\item NULL or "nquads" (default)
\item any valid serializer name: e.g. "rdfxml", "jsonld", "turtle",  "ntriples"
}

\code{rdf_base_uri}:
\itemize{
\item Default base URI to use (when serializing JSON-LD only at this time)
default is "localhost://"
}

\code{rdf_max_print}:
\itemize{
\item maximum number of lines to print from rdf, default 10
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/rdflib}
  \item Report bugs at \url{https://github.com/ropensci/rdflib/issues}
}

}
\author{
\strong{Maintainer}: Carl Boettiger \email{cboettig@gmail.com} (\href{https://orcid.org/0000-0002-1642-628X}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Bryce Mecum (\href{https://orcid.org/0000-0002-0381-3766}{ORCID}) [reviewer]
  \item Anna Krystalli (\href{https://orcid.org/0000-0002-2378-4915}{ORCID}) [reviewer]
  \item Viktor Senderov \email{vsenderov@gmail.com} (\href{https://orcid.org/0000-0003-3340-5963}{ORCID}) [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf_parse.R
\name{rdf_parse}
\alias{rdf_parse}
\title{Parse RDF Files}
\usage{
rdf_parse(
  doc,
  format = c("guess", "rdfxml", "nquads", "ntriples", "turtle", "jsonld"),
  rdf = NULL,
  base = getOption("rdf_base_uri", "localhost://"),
  ...
)
}
\arguments{
\item{doc}{path, URL, or literal string of the rdf document to parse}

\item{format}{rdf serialization format of the doc,
one of "rdfxml", "nquads", "ntriples", "turtle"
or "jsonld". If not provided, will try to guess based
on file extension and fall back on rdfxml.}

\item{rdf}{an existing rdf triplestore to extend with triples from
the parsed file.  Default will create a new rdf object.}

\item{base}{the base URI to assume for any relative URIs (blank nodes)}

\item{...}{additional parameters (not implemented)}
}
\value{
an rdf object, containing the redland world
and model objects
}
\description{
Parse RDF Files
}
\examples{
doc <- system.file("extdata", "dc.rdf", package="redland")
rdf <- rdf_parse(doc)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_nquads.R
\name{write_nquads}
\alias{write_nquads}
\title{write object out as nquads}
\usage{
write_nquads(x, file, ...)
}
\arguments{
\item{x}{an object that can be represented as nquads}

\item{file}{output filename}

\item{...}{additional parameters, see examples}
}
\description{
write object out as nquads
}
\examples{
tmp <- tempfile(fileext = ".nq")
library(datasets)

## convert data.frame to nquads
write_nquads(iris, tmp)
rdf <- read_nquads(tmp)

## or starting a native rdf object
write_nquads(rdf, tempfile(fileext = ".nq"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf.R
\name{rdf}
\alias{rdf}
\title{Initialize an \code{rdf} Object}
\usage{
rdf(
  storage = c("memory", "BDB", "sqlite", "postgres", "mysql", "virtuoso"),
  host = NULL,
  port = NULL,
  user = NULL,
  password = NULL,
  database = NULL,
  charset = NULL,
  dir = NULL,
  dsn = "Local Virtuoso",
  name = "rdflib",
  new_db = FALSE,
  fallback = TRUE
)
}
\arguments{
\item{storage}{Storage backend to use; see details}

\item{host}{host address for mysql, postgres, or virtuoso storage}

\item{port}{port for mysql (mysql storage defaults to mysql standard port, 3306)
or postgres (postgres storage defaults to postgres standard port, 4321)}

\item{user}{user name for postgres, mysql, or virtuoso}

\item{password}{password for postgres, mysql, or virtuoso}

\item{database}{name of the database to be created/used}

\item{charset}{charset for virtuoso database, if desired}

\item{dir}{directory of where to write sqlite or berkeley database.}

\item{dsn}{Virtuoso dsn, either "Local Virtuoso" or "Remote Virtuoso"}

\item{name}{name for the storage object created. Default is usually fine.}

\item{new_db}{logical, default FALSE. Create new database or connect to existing?}

\item{fallback}{logical, default TRUE. If requested storage system cannot initialize,
should \code{rdf()} fall back on memory (default) or throw an error (fallback=FALSE)?}
}
\value{
an rdf object
}
\description{
Initialize an \code{rdf} Object
}
\details{
an rdf Object is a list of class 'rdf', consisting of
three pointers to external C objects managed by the redland library.
These are the \code{world} object: basically a top-level pointer for
all RDF models, and a \code{model} object: a collection of RDF statements,
and a \code{storage} object, indicating how these statements are stored.

\code{rdflib} defaults to an in-memory hash-based storage structure.
which should be best for most use cases. For very large triplestores,
disk-based storage will be necessary. Enabling external storage devices
will require additional libraries and custom compiling. See the storage
vignette for details.
}
\examples{
x <- rdf()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf_serialize.R
\name{rdf_serialize}
\alias{rdf_serialize}
\title{Serialize an RDF Document}
\usage{
rdf_serialize(
  rdf,
  doc = NULL,
  format = c("guess", "rdfxml", "nquads", "ntriples", "turtle", "jsonld"),
  namespace = NULL,
  prefix = names(namespace),
  base = getOption("rdf_base_uri", "localhost://"),
  ...
)
}
\arguments{
\item{rdf}{an existing rdf triplestore to extend with triples from
the parsed file.  Default will create a new rdf object.}

\item{doc}{file path to write out to. If null, will write to character.}

\item{format}{rdf serialization format of the doc,
one of "rdfxml", "nquads", "ntriples", "turtle"
or "jsonld". If not provided, will try to guess based
on file extension and fall back on rdfxml.}

\item{namespace}{a named character containing the prefix to namespace bindings. \code{names(namespace)} are the prefixes, whereas \code{namespace} are the namespaces}

\item{prefix}{(optional) for backward compatibility. See \code{namespace}. It contains the matching prefixes to the namespaces in \code{namespace} and is set automatically if you provide \code{namespace} as a named character vector.}

\item{base}{the base URI to assume for any relative URIs (blank nodes)}

\item{...}{additional arguments to \code{redland::serializeToFile}}
}
\value{
rdf_serialize returns the output file path \code{doc} invisibly.
This makes it easier to use rdf_serialize in pipe chains with
\code{\link{rdf_parse}}.
}
\description{
Serialize an RDF Document
}
\examples{
infile <- system.file("extdata", "dc.rdf", package="redland")
out <- tempfile("file", fileext = ".rdf")

some_rdf <- rdf_parse(infile)
rdf_add(some_rdf,
    subject = "http://www.dajobe.org/dave-beckett",
    predicate = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type",
    object = "http://xmlns.com/foaf/0.1/Person")
rdf_serialize(some_rdf, out)

## With a namespace
rdf_serialize(some_rdf,
          out,
          format = "turtle",
          namespace = c(dc = "http://purl.org/dc/elements/1.1/",
          foaf = "http://xmlns.com/foaf/0.1/")
          )

readLines(out)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf_methods.R
\name{c.rdf}
\alias{c.rdf}
\title{Concatenate rdf Objects.}
\usage{
\method{c}{rdf}(...)
}
\arguments{
\item{...}{objects to be concatenated}
}
\description{
All subsequent rdf objects will be appended to the first rdf object
Note: this does not free memory from any of the individual rdf objects
Note: It is generally better to avoid the use of this function by passing
an existing rdf object to and rdf_parse or rdf_add objects. Multiple active
rdf objects can cause problems when using disk-based storage backends.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdf_query.R
\name{rdf_query}
\alias{rdf_query}
\title{Perform a SPARQL Query}
\usage{
rdf_query(rdf, query, data.frame = TRUE, ...)
}
\arguments{
\item{rdf}{an rdf object (e.g. from \code{\link{rdf_parse}})}

\item{query}{a SPARQL query, as text string}

\item{data.frame}{logical, should the results be returned as a data.frame?}

\item{...}{additional arguments to a redland initialize-Query}
}
\value{
a data.frame of all query results (default.)  Columns will
be named according to variable names in the SPARQL query. Returned
object values will be coerced to match the corresponding R type
to any associated datatype URI, if provided. If a column would
result in mixed classes (e.g. strings and numerics), all types
in the column will be coerced to character strings. If \code{data.frame}
is false, results will be returned as a list with each element
typed by its data URI.
}
\description{
Perform a SPARQL Query
}
\examples{
doc <- system.file("extdata", "dc.rdf", package="redland")

sparql <-
'PREFIX dc: <http://purl.org/dc/elements/1.1/>
 SELECT ?a ?c
 WHERE { ?a dc:creator ?c . }'

rdf <- rdf_parse(doc)
rdf_query(rdf, sparql)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_nquads.R
\name{read_nquads}
\alias{read_nquads}
\title{read an nquads file}
\usage{
read_nquads(file, ...)
}
\arguments{
\item{file}{path to nquads file}

\item{...}{additional arguments to \code{\link[=rdf_parse]{rdf_parse()}}}
}
\value{
an rdf object.  See \code{\link[=rdf_parse]{rdf_parse()}}
}
\description{
read an nquads file
}
\examples{
tmp <- tempfile(fileext = ".nq")
library(datasets)
write_nquads(iris, tmp)
read_nquads(tmp)

}
