---
title: 'Ecological Metadata as Linked Data'
tags:
 - linked data
 - EML
 - JSON-LD
 - RDF
 - Ecology
 
authors:
 - name: Carl Boettiger
   orcid: 0000-0002-1642-628X
   affiliation: 1
affiliations:
 - name: University of California, Berkeley
   index: 1
date: 2019-01-25
bibliography: paper.bib
---

The Ecological Metadata Language, (EML), is a widely used metadata
standard in the ecological and environmental sciences \[@Michener1997;
@Jones2006; @Leinfelder2010\]. Efforts such as the Long Term Ecological
Research (LTER) and the National Center for Ecological Analysis and
Synthesis (NCEAS) have used and driven the development of EML over the
past two decades, and now major efforts such as the NSF-sponsored
National Ecological Observatory Network (NEON), the DataONE Network, are
making an unprecedented amount of ecological data available with rich,
machine-readable metadata descriptions using EML \[@NEON; @Keller2008;
@Overpeck2011\]. EML defines a strict XML schema that ensures metadata
are machine readable and inter-operable through an XML-schema validation
process.

While this provides a predictable but extensible format, it also poses a
significant technical barrier to many researchers seeking to analyze
existing EML documents or create their own. The creation of EML in the
late 1990s also pre-dated the advent of more developer-friendly and
user-friendly technologies such as the popular JSON serialization, as
well as the advent of semantic or linked data model and it’s application
in ecoinformatics \[@Michener2012\]. More recently still, the creation
of the W3C JSON-LD \[@jsonld-w3c\] has brought combined the
developer-friendly JSON format with the powerful semantics of the
Resource Description Format (RDF) popular in informatics \[@rdf-w3c\].
This package seeks to bring both of these advantages to bear on the
rich-but-clunky XML structure used by EML. By automatically mapping EML
into JSON-LD, we are able to realize the benefits easier creation and
manipulation of EML structures as JSON objects or R list objects, while
relying on the power of JSON-LD semantics (in particular, context,
compaction and framing, \[@jsonld\]) to transform that data to and from
EML-schema-compliant XML. This approach opens numerous research
applications, ranging from the simple parsing or creation of EML files
as R list objects to complex queries and semantic reasoning over large
numbers of EML documents with semantic SPARQL queries \[@sparql-w3c\].
The package vignette provides scripted examples of each of these.

``` r
library(emld)
library(jsonlite)
library(magrittr) # for pipes
library(jqr)      # for JQ examples only
library(rdflib)   # for RDf examples only
```

## Reading EML

The `EML` package can get particularly cumbersome when it comes to
extracting and manipulating existing metadata in highly nested EML
files. The `emld` approach can leverage a rich array of tools for
reading, extracting, and manipulating existing EML files.

We can parse a simple example and manipulate is as a familiar list
object (S3 object):

``` r
f <- system.file("extdata/example.xml", package="emld")
eml <- as_emld(f)

cat(eml$dataset$title)
```

    ## Data from Cedar Creek LTER on productivity and species richness
    ##   for use in a workshop titled "An Analysis of the Relationship between
    ##   Productivity and Diversity using Experimental Results from the Long-Term
    ##   Ecological Research Network" held at NCEAS in September 1996.

## Writing EML

Because `emld` objects are just nested lists, we can create EML just by
writing
lists:

``` r
me <- list(individualName = list(givenName = "Carl", surName = "Boettiger"))

eml <- list(dataset = list(
              title = "dataset title",
              contact = me,
              creator = me),
              system = "doi",
              packageId = "10.xxx")

as_xml(eml, "ex.xml")
testthat::expect_true(eml_validate("ex.xml") )
```

Note that we don’t have to worry about the order of the elements here,
`as_xml` will re-order if necessary to validate. (For instance, in valid
EML the `creator` becomes listed before `contact`. Of course this is a
very low-level interface that does not help the user know what an EML
looks like. Creating EML from scratch without knowledge of the schema is
a job for the `EML` package and beyond the scope of the lightweight
`emld`.

# Working with EML as JSON-LD

For many applications, it is useful to merely treat EML as a list
object, as seen above, allowing the R user to leverage a standard tools
and intuition in working with these files. However, `emld` also opens
the door to new possible directions by thinking of EML data in terms of
a JSON-LD serialization rather than an XML serialization. First, owing
to it’s comparative simplicity and native data typing (e.g. of
Boolean/string/numeric data), JSON is often easier for many developers
to work with than EML’s native XML format.

## As JSON: Query with JQ

For example, JSON can be queried with with JQ, a [simple and powerful
query language](https://stedolan.github.io/jq/manual/) that also gives
us a lot of flexibility over the return structure of our results. JQ
syntax is both intuitive and well documented, and often easier than the
typical munging of JSON/list data using `purrr`. Here’s an example query
that turns EML to JSON and then extracts the north and south bounding
coordinates:

``` r
hf205 <- system.file("extdata/hf205.xml", package="emld")

as_emld(hf205) %>% 
  as_json() %>% 
  jq('.dataset.coverage.geographicCoverage.boundingCoordinates | 
       { northLat: .northBoundingCoordinate, 
         southLat: .southBoundingCoordinate }') %>%
  fromJSON()
```

    ## $northLat
    ## [1] "+42.55"
    ## 
    ## $southLat
    ## [1] "+42.42"

Nice features of JQ include the ability to do recursive descent (common
to XPATH but not possible in `purrr`) and specify the shape and names of
the return object (e.g. as a list with elements named `northLat` and
`southLat` in this case.)

## As semantic data: SPARQL queries

Another side-effect of the JSON-LD representation is that we can treat
EML as *semantic* data. This can provide a way to integrate EML records
with other data sources, and means we can query the EML using semantic
SPARQL queries. One nice thing about SPARQL queries is that, in contrast
to XPATH, JQ, or other graph queries, SPARQL always returns a
`data.frame` – a particularly convenient and familiar format for R
users.

First, we render the EML XML document into JSON-LD file:

``` r
f <- system.file("extdata/hf205.xml", package="emld")
as_emld(f) %>%
  as_json("hf205.json")
```

We can now construct a SPARQL query. SPARQL queries look like SQL
queries in that we name the columns we want with a `SELECT` command.
Unlike SQL, SPARQL allows us to walk the graph by treating these names
as variables, indicated by prefixing a `?` to the variable name. We then
use a WHERE block to define how these variables relate to each other. In
this case, we ask for the genus and species name and bounding box found
in the EML file.

``` r
sparql <-
  'PREFIX eml: <https://eml.ecoinformatics.org/eml-2.2.0/>

  SELECT ?genus ?species ?northLat ?southLat ?eastLong ?westLong 

  WHERE { 
    ?y eml:taxonRankName "genus" .
    ?y eml:taxonRankValue ?genus .
    ?y eml:taxonomicClassification ?s .
    ?s eml:taxonRankName "species" .
    ?s eml:taxonRankValue ?species .
    ?x eml:northBoundingCoordinate ?northLat .
    ?x eml:southBoundingCoordinate ?southLat .
    ?x eml:eastBoundingCoordinate ?eastLong .
    ?x eml:westBoundingCoordinate ?westLong .
  }
'
```

We can now use the `rdflib` library to execute this SPARQL query on the
EML document and display the resulting `data.frame`:

``` r
rdf <- rdf_parse("hf205.json", "jsonld")
df <- rdf_query(rdf, sparql)
df
```

    ## # A tibble: 1 x 6
    ##   genus      species  northLat southLat eastLong westLong
    ##   <chr>      <chr>       <dbl>    <dbl>    <dbl>    <dbl>
    ## 1 Sarracenia purpurea     42.6     42.4    -72.1    -72.3

# References

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build
Status](https://travis-ci.org/ropensci/emld.svg?branch=master)](https://travis-ci.org/ropensci/emld)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/cboettig/emld?branch=master&svg=true)](https://ci.appveyor.com/project/cboettig/emld)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/emld/master.svg)](https://codecov.io/github/ropensci/emld?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/emld)](https://cran.r-project.org/package=emld)
[![](https://badges.ropensci.org/269_status.svg)](https://github.com/ropensci/software-review/issues/269)
[![DOI](https://zenodo.org/badge/108223439.svg)](https://zenodo.org/badge/latestdoi/108223439)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01276/status.svg)](https://doi.org/10.21105/joss.01276)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# emld

The goal of emld is to provide a way to work with EML metadata in the
JSON-LD format. At it’s heart, the package is simply a way to translate
an EML XML document into JSON-LD and be able to reverse this so that any
semantically equivalent JSON-LD file can be serialized into EML-schema
valid XML. The package has only three core functions:

  - `as_emld()` Convert EML’s `xml` files (or the `json` version created
    by this package) into a native R object (an S3 class called `emld`,
    essentially just a `list`).
  - `as_xml()` Convert the native R format, `emld`, back into XML-schema
    valid EML.
  - `as_json()` Convert the native R format, `emld`, into `json`(LD).

## Installation

You can install emld from github with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/emld")
```

## Motivation

In contrast to the existing [EML
package](https://docs.ropensci.org/EML/), this package aims to a very
light-weight implementation that seeks to provide both an intuitive data
format and make maximum use of existing technology to work with that
format. In particular, this package emphasizes tools for working with
linked data through the JSON-LD format. This package is not meant to
replace `EML`, as it does not support the more complex operations found
in that package. Rather, it provides a minimalist but powerful way of
working with EML documents that can be used by itself or as a backend
for those complex operations. Version 2.0 of the EML R package uses
`emld` under the hood.

Note that the JSON-LD format is considerably less rigid than the EML
schema. This means that there are many valid, semantically equivalent
representations on the JSON-LD side that must all map into the same or
nearly the same XML format. At the extreme end, the JSON-LD format can
be serialized into RDF, where everything is flat set of triples
(e.g. essentially a tabular representation), which we can query
directly with semantic tools like SPARQL, and also automatically coerce
back into the rigid nesting and ordering structure required by EML. This
ability to “flatten” EML files can be particularly convenient for
applications consuming and parsing large numbers of EML files. This
package may also make it easier for other developers to build on the
EML, since the S3/list and JSON formats used here have proven more
appealing to many R developers than S4 and XML serializations.

``` r
library(emld)
library(jsonlite)
library(magrittr) # for pipes
library(jqr)      # for JQ examples only
library(rdflib)   # for RDf examples only
```

## Reading EML

The `EML` package can get particularly cumbersome when it comes to
extracting and manipulating existing metadata in highly nested EML
files. The `emld` approach can leverage a rich array of tools for
reading, extracting, and manipulating existing EML files.

We can parse a simple example and manipulate is as a familiar list
object (S3 object):

``` r
f <- system.file("extdata/example.xml", package="emld")
eml <- as_emld(f)
eml$dataset$title
#> [1] "Data from Cedar Creek LTER on productivity and species richness\n  for use in a workshop titled \"An Analysis of the Relationship between\n  Productivity and Diversity using Experimental Results from the Long-Term\n  Ecological Research Network\" held at NCEAS in September 1996."
```

## Writing EML

Because `emld` objects are just nested lists, we can create EML just by
writing lists:

``` r

me <- list(individualName = list(givenName = "Carl", surName = "Boettiger"))

eml <- list(dataset = list(
              title = "dataset title",
              contact = me,
              creator = me),
              system = "doi",
              packageId = "10.xxx")

ex.xml <- tempfile("ex", fileext = ".xml") # use your preferred file path

as_xml(eml, ex.xml)
#> NULL
eml_validate(ex.xml)
#> [1] TRUE
#> attr(,"errors")
#> character(0)
```

Note that we don’t have to worry about the order of the elements here,
`as_xml` will re-order if necessary to validate. (For instance, in valid
EML the `creator` becomes listed before `contact`.) Of course this is a
very low-level interface that does not help the user know what an EML
looks like. Creating EML from scratch without knowledge of the schema is
a job for the `EML` package and beyond the scope of the lightweight
`emld`.

# Working with EML as JSON-LD

For many applications, it is useful to merely treat EML as a list
object, as seen above, allowing the R user to leverage a standard tools
and intuition in working with these files. However, `emld` also opens
the door to new possible directions by thinking of EML data in terms of
a JSON-LD serialization rather than an XML serialization. First, owing
to it’s comparative simplicity and native data typing (e.g. of
Boolean/string/numeric data), JSON is often easier for many developers
to work with than EML’s native XML format.

## As JSON: Query with JQ

For example, JSON can be queried with with JQ, a [simple and powerful
query language](https://stedolan.github.io/jq/manual/) that also gives
us a lot of flexibility over the return structure of our results. JQ
syntax is both intuitive and well documented, and often easier than the
typical munging of JSON/list data using `purrr`. Here’s an example query
that turns EML to JSON and then extracts the north and south bounding
coordinates:

``` r
hf205 <- system.file("extdata/hf205.xml", package="emld")

as_emld(hf205) %>% 
  as_json() %>% 
  jq('.dataset.coverage.geographicCoverage.boundingCoordinates | 
       { northLat: .northBoundingCoordinate, 
         southLat: .southBoundingCoordinate }') %>%
  fromJSON()
#> $northLat
#> [1] "+42.55"
#> 
#> $southLat
#> [1] "+42.42"
```

Nice features of JQ include the ability to do recursive descent (common
to XPATH but not possible in `purrr`) and specify the shape of the
return object. Some prototype examples of how we can use this to
translate between EML and <https://schema.org/Dataset> representations
of the same metadata can be found in
<https://github.com/ropensci/emld/blob/master/notebook/jq_maps.md>

## As semantic data: SPARQL queries

Another side-effect of the JSON-LD representation is that we can treat
EML as “semantic” data. This can provide a way to integrate EML records
with other data sources, and means we can query the EML using semantic
SPARQL queries. One nice thing about SPARQL queries is that, in contrast
to XPATH, JQ, or other graph queries, SPARQL always returns a
`data.frame` which is a particularly convenient format. SPARQL queries
look like SQL queries in that we name the columns we want with a
`SELECT` command. Unlike SQL, these names act as variables. We then use
a WHERE block to define how these variables relate to each other.

``` r
f <- system.file("extdata/hf205.xml", package="emld")
hf205.json <- tempfile("hf205", fileext = ".json") # Use your preferred filepath

as_emld(f) %>%
  as_json(hf205.json)

prefix <- paste0("PREFIX eml: <eml://ecoinformatics.org/", eml_version(), "/>\n")
sparql <- paste0(prefix, '

  SELECT ?genus ?species ?northLat ?southLat ?eastLong ?westLong 

  WHERE { 
    ?y eml:taxonRankName "genus" .
    ?y eml:taxonRankValue ?genus .
    ?y eml:taxonomicClassification ?s .
    ?s eml:taxonRankName "species" .
    ?s eml:taxonRankValue ?species .
    ?x eml:northBoundingCoordinate ?northLat .
    ?x eml:southBoundingCoordinate ?southLat .
    ?x eml:eastBoundingCoordinate ?eastLong .
    ?x eml:westBoundingCoordinate ?westLong .
  }
')
  
rdf <- rdf_parse(hf205.json, "jsonld")
df <- rdf_query(rdf, sparql)
df
#> # A tibble: 0 x 0
```

-----

Please note that the `emld` project is released with a [Contributor Code
of Conduct](https://docs.ropensci.org/emld/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# emld 0.5.1

- don't build vignette on machines that don't have packages listed in Suggests
  (these packages must only be used conditionally)

# emld 0.5.0

User-facing changes:

- Rewrote logic around how `emld` works with `schemaLocation`:
  - The `schemaLocation` argument on `as_xml` used to fill in a default value that pointed to a local copy of `eml.xsd`. Now, it automatically fills in a web-resolvable location to make validating documents easier out of the box. 
  
  See [#44](https://github.com/ropensci/emld/issues/44) & [#45](https://github.com/ropensci/emld/issues/45).

- Re-worked the logic for handling `schemaLocation` when serializing to XML:
  - When `schemaLocation` is present on the `emld` object, it is used verbatim (no change in API here) regardless of the `schemaLocation` argument of `as_xml`.
  - When `schemaLocation` is absent on the `emld` object:
    - `as_xml(..., schemaLocation = TRUE)` causes a value to be guessed.
    - `as_xml(..., schemaLocation = FALSE)` explicitly prevents a value from being filled in when serialized.
    - `as_xml(..., schemaLocation = "Some value)` explicitly sets the provided value.

  See [#44](https://github.com/ropensci/emld/issues/44) & [#45](https://github.com/ropensci/emld/issues/45).

- `emld::eml_validate` no longer depends on `schemaLocation` to determine the correct XSD to use during schema validation and now uses two helpers (See below) to find the correct schema file. See [#52](https://github.com/ropensci/emld/issues/44) & [#45](https://github.com/ropensci/emld/issues/53).
- `emld::eml_version` now allows specifying the version without the `eml-` prefix, like `eml_version("2.1.1"), and will throw a warning when it gets output that doesn't 'look right rather than silently failing.
- Fixed a bug where the EML 2.1.1 units dictionary was being used for EML 2.2.0 docs which would cause spurious validation errors. See [#56](https://github.com/ropensci/emld/issues/56).

Developer-facing (non-exported) changes:

- Added two new helper methods:
  1. `find_real_root_name(doc : xml_document) : list(prefix : character, name: character)` which returns the namespace prefix and the local name of the root element on an `xml_document`.
  2. `guess_root_schema(doc : xml_document) : list(module : character, version : character, namespace : character)` which returns the module, schema version, and namespace URI of the root element on an `xml_document`.
- `schemaLocation` is now ignored during roundtrip testing because of the new (above) behavior of `emld` with respect to `schemaLocation`.
- Roundtrip testing can now handle documents that are supposed to be invalid but still roundtripped. Specify intentionally invalid files by adding "invalid" (case insensitive) to the filename in `inst/tests`.

Other changes:

- Minor tweaks to the README. Thanks @jeanetteclark

# emld 0.4.0

- Fixed serialization bug for `references` attributes [#48](https://github.com/ropensci/emld/issues/48)
- Fixed validation bug: packageId is now used as an identifier for checking uniqueness, and no more errors for annotation elements in additionalMetadata. [#49](https://github.com/ropensci/emld/pull/49)
- Fixed validation bug: XPath was referencing the element rather than the attribute `references`. [#47](https://github.com/ropensci/emld/pull/47)

# emld 0.3.0

- Updated package to support version 2.2.0 of EML. [#40](https://github.com/ropensci/emld/pull/40). See the [EML website](https://eml.ecoinformatics.org/whats-new-in-eml-2-2-0.html) for more information on the 2.2.0 release.
- Fixed a minor XML serialization issue with `TextType` nodes where extra whitespace was being added. [#37](https://github.com/ropensci/emld/pull/37).
- Relaxed `eml_validate`'s behavior when validating custom units. [#35](https://github.com/ropensci/emld/pull/35).

# emld 0.2.0

* Implemented changes requested by rOpenSci review, as detailed in 
  [#30](https://github.com/ropensci/emld/pull/30)

# emld 0.1.1

* Version submitted to rOpenSci review
* Added a `NEWS.md` file to track changes to the package.
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


Contributing Guidelines
=======================


Repository structure
--------------------

This repository is structured as a standard R package
following the conventions outlined in the [Writing R
extensions](http://cran.r-project.org/doc/manuals/R-exts.html) manual.
A few additional files are provided that are not part of the built
R package and are listed in `.Rbuildignore`, such as `.travis.yml`,
which is used for continuous testing and integration.


Code
----

All code for this package is found in `R/`, (except compiled source
code, if used, which is in `/src`).  All functions should be thoroughly
documented with `roxygen2` notation; see Documentation. Code should
conform to our [Style guide](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md)

Testing
-------

Any new feature or bug-fix should include a unit-test demonstrating the
change.  Unit tests follow the `testthat` framework with files in
`tests/testthat`.  Please make sure that the testing suite passes
before issuing a pull request.  This can be done by running `check()`
from the `devtools` package, which will also check for consistent
documentation, etc.


This package uses the [travis](https://github.com/craigcitro/r-travis)
continuous testing mechanism for R to ensure that the test suite is run
on each push to Github.  An icon at the top of the README.md indicates
whether or not the tests are currently passing. 

This package also uses
[codecov.io](https://codecov.io/) to
measure test coverage.  While not all code can be covered by automated 
tests (in particular, functions involving user prompts), try to avoid
decreasing coverage by writing unit tests for any contributed code. 
Codecov.io will flag PRs that decrease coverage. 


Documentation
-------------

All of the function documentation is generated automatically.
Please do not edit any of the documentation files in `man/`
or the `NAMESPACE`.  Instead, construct the appropriate
[roxygen2](https://github.com/klutometis/roxygen) documentation in the
function files in `R/` themselves.  The documentation is then generated
by running the `document()` function from the `devtools` package.  Please
consult the [Advanced R programming](http://adv-r.had.co.nz/) guide if
this workflow is unfamiliar to you.  Note that functions should include
examples in the documentation. Please use `\dontrun` for examples that
take more than a few seconds to execute or require an internet connection.

Likewise, the README.md file in the base directory should not be edited
directly.  This file is created automatically from code that runs the
examples shown, helping to ensure that they are functioning as advertised
and consistent with the package README vignette.  Instead, edit the
`README.Rmd` source file in `manuscripts` and run `make` to build
the README.

General Development Goals & Guidelines
---------------------------------------

1. Not having too many high-level functions, 
2. Using sensible defaults, (driven by use cases).
3. Docs should point advanced users to the lower-level API when they need special cases.  
4. Maintain a consistent user-facing API.   
Dear CRAN,

This release provides the requested documentation fix in use of `...`, and addresses several 
unrelated bugs, as documented in NEWS.

Cheers,

Carl Boettiger

Mapping between Schema.Org and EML via JQ
=========================================

We will use J-Query as XSLT-like stylesheets. These two examples shipping in this package are currently a work in progress:

``` r
library(jqr)
```

``` r
eml_to_schema <- readr::read_file(system.file("jq/eml_to_schema.jq", package="emld"))
schema_to_eml <-  readr::read_file(system.file("jq/schema_to_eml.jq", package="emld"))
```

Let's map a more complete EML document into schema.org:

``` r
eml <- readr::read_file("https://raw.githubusercontent.com/cboettig/emld/master/inst/extdata/hf205.json")
jq(eml, eml_to_schema)
```

    ## {
    ##     "id": "HF205",
    ##     "type": null,
    ##     "temporalCoverage": "2012-06-01/2013-12-31",
    ##     "spatialCoverage": {
    ##         "description": "Harvard Forest Greenhouse, Tom Swamp Tract (Harvard Forest)",
    ##         "geo": {
    ##             "box": "+42.42 -72.29 +42.55 -72.10"
    ##         }
    ##     }
    ## }

Convert a dataset marked up in <http://schema.org/Dataset> terms into EML

``` r
schema <- readr::read_file("../inst/extdata/schema-org-dataset.json")

jq(schema, schema_to_eml)
```

    ## {
    ##     "@id": null,
    ##     "@type": "Dataset",
    ##     "coverage": {
    ##         "temporalCoverage": {
    ##             "rangeOfDates": {
    ##                 "beginDate": "1950-01-01",
    ##                 "endDate": "2013-12-18"
    ##             }
    ##         },
    ##         "geographicCoverage": {
    ##             "geographicDescription": null,
    ##             "boundingCoordinates": {
    ##                 "westBoundingCoordinate": "-65.0",
    ##                 "eastBoundingCoordinate": "172.0",
    ##                 "northBoundingCoordinate": "72.0",
    ##                 "southBoundingCoordinate": "18.0"
    ##             }
    ##         }
    ##     }
    ## }
R Notebook
================

``` r
library(jsonld)
library(jsonlite)
library(tidyverse)
library(printr)
```

Switching Contexts: Expansion & Compaction
------------------------------------------

-   Good: homonyms (same word has different meanings), synonyms (different word has ~ same meaning)
-   Fails: Mapping isn't 1:1

### A simple example

``` r
ex <- 
 '{ 
      "date": "1993-10-02",
      "site": "N654",
      "scientificName": "Picea rubens",
      "area": 2,
      "count": 26
 }'
```

No need to alter your own data files/source data files. Simply define a context that maps those terms into the reference vocabulary. Data can use whatever fields it likes.

``` r
map <- 
'{ 
  "@context": {
    "@vocab": "http://example.com/",
    "scientificName": "http://example.com/species"
  }
}'
target <- '{"@context": {"@vocab": "http://example.com/"}}'
add_context(ex, map) %>% 
  jsonld_compact(target)
```

    ## {
    ##   "@context": {
    ##     "@vocab": "http://example.com/"
    ##   },
    ##   "area": 2,
    ##   "count": 26,
    ##   "date": "1993-10-02",
    ##   "site": "N654",
    ##   "species": "Picea rubens"
    ## }

### Using tabular schema

``` r
df <- 
tribble(
  ~Date,       ~Site, ~Species, ~Area, ~Count,
  "10/1/1993", "N654", "PIRU",   2.0,     26,
  "10/3/1994", "N654", "PIRU",   2.0,     29,
  "10/1/1993", "N654", "BEPA",   1.0,     3
)

map <- 
'{ 
  "@context": {
    "@vocab": "http://example.com/",
    "Species": "http://example.com/scientificName"
  }
}'
target <- '{"@context": {"@vocab": "http://example.com/"}}'

df %>% 
  toJSON() %>% 
  add_context(map) %>% 
  jsonld_compact(target) %>% 
  drop_context() %>%
  fromJSON()
```

|     Area|     Count| Date         | Site    | scientificName      |
|--------:|---------:|:-------------|:--------|:--------------------|
|        2|        26| 10/1/1993    | N654    | PIRU                |
|        2|        29| 10/3/1994    | N654    | PIRU                |
|        1|         3| 10/1/1993    | N654    | BEPA                |
|  Wow doe|  s that s| eem like a r | ound-ab | out way to just do: |

``` r
df %>% rename(scientificName = Species)
```

| Date      | Site | scientificName |  Area|  Count|
|:----------|:-----|:---------------|-----:|------:|
| 10/1/1993 | N654 | PIRU           |     2|     26|
| 10/3/1994 | N654 | PIRU           |     2|     29|
| 10/1/1993 | N654 | BEPA           |     1|      3|

and it is. But a few things are worth noting:

-   Robust to what columns are actually present. Instead of a default vocabulary we can map all terms explicitly.
-   Terms not expanded by the map (that is, not part of our context) will be dropped. (`Count` in the example below)
-   Terms not in our target context remain as external keys, and are not compacted.
-   Typesafe transforms. Terms (keys/column names) with the wrong type will not compact if the types do not match the target transform. This way they are not actually lost, but cannot be confused.

``` r
map <- 
'{ 
  "@context": {
    "ex": "http://example.com/",
    "Species": "ex:scientificName",
    "Site": "ex:Site",
    "Date": "ex:Date",
    "Area": "ex:Area"
  }
}'

target <- '{
"@context": {
  "ex": "http://example.com/",
  "scientificName": "ex:scientificName",
  "Site": "ex:Site",
  "Count": "ex:count",
  "Date": {"@id": "ex:Date", "@type": "ex:Date"}
}}'

df %>% 
  toJSON() %>% 
  add_context(map) %>% 
  jsonld_compact(target) %>% 
  drop_context() %>%
  fromJSON() 
```

|  ex:Area| ex:Date   | Site | scientificName |
|--------:|:----------|:-----|:---------------|
|        2| 10/1/1993 | N654 | PIRU           |
|        2| 10/3/1994 | N654 | PIRU           |
|        1| 10/1/1993 | N654 | BEPA           |

``` r
map <- '{ "@context": { "@vocab": "https://data.berkeley.edu/" }}'

target <- '{
"@context": {
  "ucb": "https://data.berkeley.edu/",
  "hav": "http://data.harvard.edu/",
  "scientificName": "ucb:Species",
  "Site": "ucb:Site",
  "Count": "ucb:Count",
  "Date": "ucb:Date",
  "Area": "ucb:Area",
  "Density": "hav:density"
}}'

df %>% 
  toJSON() %>% 
  add_context(map) %>% 
  jsonld_compact(target) %>% 
  drop_context() %>%
  fromJSON() -> ucsb_df

ucsb_df %>% mutate(density = Count / Area)
```

|  Area|  Count| Date      | Site | scientificName |  density|
|-----:|------:|:----------|:-----|:---------------|--------:|
|     2|     26| 10/1/1993 | N654 | PIRU           |     13.0|
|     2|     29| 10/3/1994 | N654 | PIRU           |     14.5|
|     1|      3| 10/1/1993 | N654 | BEPA           |      3.0|

------------------------------------------------------------------------

Integration: Flattening shares id
---------------------------------

``` r
context <- '{   "@context": { "@vocab": "http://example.com/" } }'

ex <- '{
  "@context": { "@vocab": "http://example.com/" },
  "datasets": [
    {
      "@id": "http://global_unique_id",
      "@type": "dataset",
      "image": "http://image_uri",
      "name": "stuff"
    },
  
    {
      "@type": "dataset",
      "@id": "http://global_unique_id",
      "name": "Species name"
    }
  ]
}'

jsonld_flatten(ex, context)
```

    ## {
    ##   "@context": {
    ##     "@vocab": "http://example.com/"
    ##   },
    ##   "@graph": [
    ##     {
    ##       "@id": "_:b0",
    ##       "datasets": {
    ##         "@id": "http://global_unique_id"
    ##       }
    ##     },
    ##     {
    ##       "@id": "http://global_unique_id",
    ##       "@type": "dataset",
    ##       "image": "http://image_uri",
    ##       "name": [
    ##         "stuff",
    ##         "Species name"
    ##       ]
    ##     }
    ##   ]
    ## }

------------------------------------------------------------------------

hdf5: column data, with units, metadata


> In the `Motivations` section of `README.md` the purpose and method of "extending EML with other semantic vocabularies"  isn't clear (to me). Consider adding a section such as "Extending EML" with an example, if this is within the scope of the package.

Good point, I've removed this now since it is merely confusing.  What I had in mind here would be contingent on EML 2.2.0's support for arbitrary RDF annotations in EML.  In principle, a user could just add another term, say, `schema:editor` onto an `eml:dataset` object, and `emld` would translate it into the appropriate (and potentially less intuitive to construct) EML `annotation` element.  But this isn't implemented yet (see issue #1), and would probably be more trouble than help.  

> The help text for the `as_xml` `Description` contains only the text `as_xml` and does not describe what the function does. Also, the return type is not specified. The same is true for both `as_emld` and `as_json`, `template`.

fixed. 


> There is no documentation index is available for the package.

not actually sure to what this refers?  complete pdf & pkgdown-based versions of the package manual / docs should now be in place.

> There is no documentation for `?emld'` or `?emld-package`.

Thanks, fixed!

> This may be outside the scope of the package, but no functions, strategies or examples are presented in the documentation for simple editing of EML data, i.e. (reading, simple edit, write out). This would be very useful, but if this functionality presents to much overlap with the `EML` package, pleas disregard this comment.

Right, this is really the scope of the `EML` package.  As the README shows, it is possible to a list-based approach to create a simple (schema-valid) EML file, or make a minor edit, but this is unlikely to scale well without richer functions in EML. 

#### Functionality

> The `as_json`, `as_emld`, `as_xml` functions have a clear purpose and work as expected.

> In addition to testing using the provided code samples, the following checks for a complex EML document were performed for a couple of EML documents such as https://goa.nceas.ucsb.edu/#view/urn:uuid:3249ada0-afe3-4dd6-875e-0f7928a4c171:

> - verified that `as_json` produced valid JSON-LD
> - verified that the following commands produced the original EML file (round trip):

```
x <- as_emld('metadata.xml') # using
as_xml(x, "new.xml")
```

yay, thanks for checking this!

#### R Source

The following potentially unresolved items exist, as show from a quick scan of the source code:

```
$ grep FIXME *
as_emld.R:    ## FIXME technically this assumes only our context
as_xml.R:## FIXME drop NAs too?
eml_validate.R:    ##  FIXME shouldn't have to write to tempfile,
eml_validate.R:  # FIXME technically must be unique inside parent system only...
emld-methods.R:## FIXME: Print method should drop context
```

Good call, these have now been resolved and the notes removed.  (some had already been resolved, or weren't quite accurate)

#### Misc

devtools::check() reports non-standard file `LICENSE.md`

fixed (added to `.Rbuildignore`)
The primary goal of this project is to determine
experimentally the amount of lead time required to prevent a state
change. To achieve this goal, we will (1) experimentally induce state
changes in a natural aquatic ecosystem - the Sarracenia microecosystem;
(2) use proteomic analysis to identify potential indicators of states
and state changes; and (3) test whether we can forestall state changes
by experimentally intervening in the system. This work uses state-of-the
art molecular tools to identify early warning indicators in the field
of aerobic to anaerobic state changes driven by nutrient enrichment
in an aquatic ecosystem. The study tests two general hypotheses: (1)
proteomic biomarkers can function as reliable indicators of impending
state changes and may give early warning before increasing variances
and statistical flickering of monitored variables; and (2) well-timed
intervention based on proteomic biomarkers can avert future state changes
in ecological systems.
## General Protocols

Field methods. All experiments will be carried out in the greenhouse at Harvard Forest. We have developed an instrumentation system that allows us to collect continuous dissolved [O2] measurements: dedicated micro-probes (DO-166MT; Lazar Research Laboratories: http://www.shelfscientific.com/) connected to multiplexers and data loggers (AM16/32B multiplexer, CR-1000 datalogger and control system [Campbell Scientific: http://www.cambellsci.com]). The initial ecosystem composition in all experimental plants will be standardized by seeding each pitcher with a 10-ml inoculum of liquid collected from pitchers growing at Tom Swamp.  In all experiments, prey will be supplied to pitchers as standardized aliquots of dried and finely ground bald-faced hornets (Dolichovespula maculata; Hymenoptera: Vespidae), which we collect in quantity throughout New England. Both hornets and ants (the latter are the dominant prey of S. purpurea) are hymenoptera, and have nearly identical C:N ratios (hornets: 3.97; common bog-dwelling ants [Tapinoma sessile and Myrmica lobifrons]: 3.37), but on average hornets have greater than 100 times the dry mass of these ants, and are easier to collect and process as a standardized food source. Additions of prey, either as large "pulses" or chronic "presses" are analogous to the enrichment and eutrophication that occur in aquatic "green" food webs in which phytoplankton abundance is boosted through addition of limiting nutrients. In "brown" food webs such as the Sarracenia microecosystem, detritus - not primary production - is at the base of the web, and our treatments boost this material as would happen through increases in arthropod prey capture78 or through nitrogen-enriched precipitation.

Proteomic analysis. Proteomic profiles of microbial communities are determined after separating the microbial fraction from the pitcher fluid, prey, and other detritus. The microbial "pellet" is subjected to SDS-PAGE (sodium dodecyl sulfate polyacrylamide gel) electrophoresis; bands are cut out and digested in-gel with trypsin. Tryptic peptides are subjected to LC-MS/MS (liquid chromatography tandem mass spectrometry) for peptide and protein identification. Absolute abundance of peptides and proteins are quantified using AQUA (Absolute QUAntification) analysis109.
      
## Specific Experiments

Experiment #1. Effects of nutrient enrichment on state changes and [O2] profiles. This experiment alters nutrient enrichment rates to characterize the [O2] profile and the transition to the anaerobic state. The experimental design is a one-way layout with 5 treatment groups: one control (no enrichment) and 4 enrichment levels (0.125, 0.25, 0.5, 1.0 mg prey added ml-1 d-1). One plant is assigned to each treatment group, and the entire set is replicated 6 times over successive weeks. [O2] is monitored continuously for 4 days to characterize state changes and tipping points under different enrichment rates. This experiment tracks a continuous [O2] profile but does not include proteomic analysis. The purpose of Experiment #1 is to identify an enrichment rate E that generates a long pre-tipping period before transition time T to the anaerobic state. This enrichment rate will be used in Experiments #2 - #4.

Experiment #2. Identification of early intervention time and characterization of aerobic and anaerobic proteomes. This experiment will use the single enrichment rate E determined from Experiment #1 and impose different intervention times I at which nutrient enrichment will be terminated. Thus, this experiment will identify the latest time I* at which it is possible to intervene and stop the transition to the anaerobic state by halting enrichment. The [O2] profile will again be monitored continuously over 10 days to measure the state of the system. From Experiment #1, the transition time T to the anaerobic state with no intervention will be known. We will use one control group (no prey addition) and ten levels of intervention time (all with the same enrichment rate E) as a proportion of T (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0). Six plants will be assigned randomly to each of the 11 treatments in a randomized one-way layout and [O2] profiles will be monitored continuously. In addition to the [O2] profiles, we will also characterize the protein profiles of aerobic and anaerobic pitchers in all 11 treatment groups at the end of the experiment.

After the plants are harvested, we will create proteomic profiles of the predominantly bacterial portion (centrifuged pellet) of the pitcher fluid from each plant, as described in General Protocols. Thus, 66 separate pellet-fraction samples will be analyzed by SDS-PAGE. After examining the SDS-PAGE profiles, approximately ten proteins that show dynamic patterns consistent with state change and five that do not change will be cut from the gel, subjected to in-gel tryptic digestion and a portion of the tryptic peptides will be analyzed by LC-MS/MS. Using these data, we will choose three identified peptides from each protein for peptide synthesis such that each synthesized peptide contains stable isotope labels (AQUA peptides) to distinguish them by mass from the native peptides. We will then quantify all 45 of the native peptides from the original samples using a known amount of each AQUA peptide spiked into the tryptic digest. The AQUA analysis of proteins that do not show changes will be used for normalization between samples. These data will be used to independently identify the current state of the system and forecast the time-to-transition.

We will use Sequest searches for initial identification of peptides; relevant scores including Xcorr and Delta Cn values will be given for each peptide. Other peptides will be identified by de novo sequencing using PepNovo; all PepNovo scores will likewise be given including any N- or C-terminal gaps. Mass error in ppm will be reported for each precursor ion. We will use standard multivariate analysis to search for distinctive proteomic profiles114 that characterize aerobic and anaerobic ecosystems, as well as ecosystems that developed with and without inputs of photosynthetic O2 and plant metabolites.
        
Experiment #3. Identification of diagnostic proteins. Using Experiments #1 and 2, we will have identified an enrichment rate E with a long pre-tipping period and an intervention time I* before which mitigation and termination of enrichment will prevent eutrophication. Experiment #3 will characterize the mean and variance of the protein profile before and after I*. We are especially interested in identifying proteins that increase rapidly in abundance (or variance) well before the onset of flickering in [O2] and before the transition time T from the aerobic to the anaerobic state.

A cohort of 100 plants all will be fed at rate E (determined from Experiment #1), with intervention time I* determined from Experiment #2, although no intervention will be used in this "press" experiment so that we can contrast proteins before and after the state change. At seven times before I* and three times after I*, we will harvest 10 randomly chosen plants. At each prescribed harvest time, we will measure [O2] and collect samples from each plant for proteomic screening using both SDS-PAGE and AQUA analysis. This experiment will identify proteins that rise quickly in abundance during the pre- I* period and can be used as early indicators of a future tipping point. Because different plants will be harvested at each time period, this is a one-way ANOVA design, with pre-and post- I* a priori contrasts. A randomization test will be used to determine whether variances in protein expression differ through time. During these analyses we will use the data from the AQUA peptides and from known amounts of protein standards, such as bovine serum albumin, to approximate the amount of protein in a given coomassie-stained SDS-PAGE gel band. The reason for doing this is to provide a fast "real-time" assay based just on expression in the SDS-PAGE. This rapid assay will be used in Experiment #4.

Experiment # 4. Proof-of-application. This experiment will provide a benchmark test of our methods and their ability to correctly identify tipping points. A cohort of 100 plants will each be fitted with [O2] probes and started on the enrichment regime. Two times per day, we will collect 3 plants each, pool their contents, and conduct a rapid screen in the lab with SDS-PAGE for the diagnostic proteins that were identified in Experiment #3. We will use the protein expression in the gel to delineate an "early" and a "late" mitigation strategy. As soon as diagnostic proteins measured in the SDS-gels are at abundances that signal we are at 0.5×I* - approximately one-half of the way to the latest intervention time - we will randomly select one third of the remaining plants for mitigation and termination of enrichment (the "early" mitigation strategy). We will continue to harvest plants from the remainder of the cohort and monitor proteins. As soon as diagnostic proteins signal we are at 0.75 times I*, we will randomly select one half of the remaining plants for mitigation and termination of enrichment (the "late" mitigation strategy). The remaining plants (approximately one sixth to one third of the original cohort) will continue to be enriched. We will monitor [O2] in all 3 groups (no-mitigation control, early mitigation, late mitigation) until all plants reach a new [O2] equilibrium. If the protein markers are successful, the proportion of food webs that remain aerobic will be significantly higher in the two mitigation treatments than in the no-mitigation control.
