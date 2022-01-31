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
---
output: github_document
---

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/emld.svg?branch=master)](https://travis-ci.org/ropensci/emld) 
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/cboettig/emld?branch=master&svg=true)](https://ci.appveyor.com/project/cboettig/emld)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/emld/master.svg)](https://codecov.io/github/ropensci/emld?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/emld)](https://cran.r-project.org/package=emld)
[![](https://badges.ropensci.org/269_status.svg)](https://github.com/ropensci/software-review/issues/269)
[![DOI](https://zenodo.org/badge/108223439.svg)](https://zenodo.org/badge/latestdoi/108223439)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01276/status.svg)](https://doi.org/10.21105/joss.01276)

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# emld

The goal of emld is to provide a way to work with EML metadata in the JSON-LD format. At it's heart, the package is simply a way to translate an EML XML document into JSON-LD and be able to reverse this so that any semantically equivalent JSON-LD file can be serialized into EML-schema valid XML.  The package has only three core functions:

- `as_emld()` Convert EML's `xml` files (or the `json` version created by this package) into a native R object (an S3 class called `emld`, essentially just a `list`).
- `as_xml()` Convert the native R format, `emld`, back into XML-schema valid EML.
- `as_json()` Convert the native R format, `emld`, into `json`(LD).


## Installation

You can install emld from github with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/emld")
```



## Motivation


In contrast to the existing [EML package](https://docs.ropensci.org/EML/), this package aims to a very light-weight implementation that seeks to provide both an intuitive data format and make maximum use of existing technology to work with that format.  In particular, this package emphasizes tools for working with linked data through the JSON-LD format.  This package is not meant to replace `EML`, as it does not support the more complex operations found in that package.  Rather, it provides a minimalist but powerful way of working with EML documents that can be used by itself or as a backend for those complex operations.  Version 2.0 of the EML R package uses `emld` under the hood.  

Note that the JSON-LD format is considerably less rigid than the EML schema.  This means that there are many valid, semantically equivalent representations on the JSON-LD side that must all map into the same or nearly the same XML format.  At the extreme end, the JSON-LD format can be serialized into RDF, where everything is flat set of triples (e.g. essentially a tabular representation), which we can query directly with semantic tools like SPARQL, and also automatically coerce back into the rigid nesting and ordering structure required by EML.  This ability to "flatten" EML files can be particularly convenient for applications consuming and parsing large numbers of EML files.  This package may also make it easier for other developers to build on the EML, since the S3/list and JSON formats used here have proven more appealing to many R developers than S4 and XML serializations.


```{r}
library(emld)
library(jsonlite)
library(magrittr) # for pipes
library(jqr)      # for JQ examples only
library(rdflib)   # for RDf examples only

```


## Reading EML

The `EML` package can get particularly cumbersome when it comes to extracting and manipulating existing metadata in highly nested EML files.  The `emld` approach can leverage a rich array of tools for reading, extracting, and manipulating existing EML files. 

We can parse a simple example and manipulate is as a familiar list object (S3 object):

```{r}
f <- system.file("extdata/example.xml", package="emld")
eml <- as_emld(f)
eml$dataset$title
```


## Writing EML

Because `emld` objects are just nested lists, we can create EML just by writing lists: 

```{r}

me <- list(individualName = list(givenName = "Carl", surName = "Boettiger"))

eml <- list(dataset = list(
              title = "dataset title",
              contact = me,
              creator = me),
              system = "doi",
              packageId = "10.xxx")

ex.xml <- tempfile("ex", fileext = ".xml") # use your preferred file path

as_xml(eml, ex.xml)
eml_validate(ex.xml)
```

Note that we don't have to worry about the order of the elements here, `as_xml` will re-order if necessary to validate. (For instance, in valid EML the `creator` becomes listed before `contact`.)   Of course this is a very low-level interface that does not help the user know what an EML looks like. Creating EML from scratch without knowledge of the schema is a job for the `EML` package and beyond the scope of the lightweight `emld`.  


# Working with EML as JSON-LD 

For many applications, it is useful to merely treat EML as a list object, as seen above, allowing the R user to leverage a standard tools and intuition in working with these files.  However, `emld` also opens the door to new possible directions by thinking of EML data in terms of a JSON-LD serialization rather than an XML serialization.  First, owing to it's comparative simplicity and native data typing (e.g. of Boolean/string/numeric data), JSON is often easier for many developers to work with than EML's native XML format.  


## As JSON: Query with JQ 

For example, JSON can be queried with with JQ, a [simple and powerful query language](https://stedolan.github.io/jq/manual/) that also gives us a lot of flexibility over the return structure of our results.  JQ syntax is both intuitive and well documented, and often easier than the typical munging of JSON/list data using `purrr`.  Here's an example query that turns EML to JSON and then extracts the north and south bounding coordinates:


```{r}
hf205 <- system.file("extdata/hf205.xml", package="emld")

as_emld(hf205) %>% 
  as_json() %>% 
  jq('.dataset.coverage.geographicCoverage.boundingCoordinates | 
       { northLat: .northBoundingCoordinate, 
         southLat: .southBoundingCoordinate }') %>%
  fromJSON()
```

Nice features of JQ include the ability to do recursive descent (common to XPATH but not possible in `purrr`) and specify the shape of the return object.  Some prototype examples of how we can use this to translate between EML and <https://schema.org/Dataset> representations of the same metadata can be found in <https://github.com/ropensci/emld/blob/master/notebook/jq_maps.md>



## As semantic data: SPARQL queries


Another side-effect of the JSON-LD representation is that we can treat EML as "semantic" data.  This can provide a way to integrate EML records with other data sources, and means we can query the EML using semantic SPARQL queries.  One nice thing about SPARQL queries is that, in contrast to XPATH, JQ, or other graph queries, SPARQL always returns a `data.frame` which is a particularly convenient format. SPARQL queries look like SQL queries in that we name the columns we want with a `SELECT` command.  Unlike SQL, these names act as variables.  We then use a WHERE block to define how these variables relate to each other.  


```{r}
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
```






----

Please note that the `emld` project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/emld/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
output: github_document
---



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


The Ecological Metadata Language, (EML), is a widely used metadata standard in the ecological and environmental sciences [@Michener1997; @Jones2006; @Leinfelder2010]. Efforts such as the Long Term Ecological Research (LTER) and the National Center for Ecological Analysis and Synthesis (NCEAS) have used and driven the development of EML over the past two decades, and now major efforts such as the NSF-sponsored National Ecological Observatory Network (NEON), the DataONE Network, are making an unprecedented amount of ecological data available with rich, machine-readable metadata descriptions using EML  [@NEON; @Keller2008; @Overpeck2011]. EML defines a strict XML schema that ensures metadata are machine readable and inter-operable through an XML-schema validation process.  

While this provides a predictable but extensible format, it also poses a significant technical barrier to many researchers seeking to analyze existing EML documents or create their own. The creation of EML in the late 1990s also pre-dated the advent of more developer-friendly and user-friendly technologies such as the popular JSON serialization, as well as the advent of semantic or linked data model and it's application in ecoinformatics [@Michener2012].  More recently still, the creation of the W3C JSON-LD [@jsonld-w3c] has brought combined the developer-friendly JSON format with the powerful semantics of the Resource Description Format (RDF) popular in informatics [@rdf-w3c].  This package seeks to bring both of these advantages to bear on the rich-but-clunky XML structure used by EML.  By automatically mapping EML into JSON-LD, we are able to realize the benefits easier creation and manipulation of EML structures as JSON objects or R list objects, while relying on the power of JSON-LD semantics (in particular, context, compaction and framing, [@jsonld]) to transform that data to and from EML-schema-compliant XML.  This approach opens numerous research applications, ranging from the simple parsing or creation of EML files as R list objects to complex queries and semantic reasoning over large numbers of EML documents with semantic SPARQL queries [@sparql-w3c].  The package vignette provides scripted examples of each of these. 



```{r}
library(emld)
library(jsonlite)
library(magrittr) # for pipes
library(jqr)      # for JQ examples only
library(rdflib)   # for RDf examples only

```


## Reading EML

The `EML` package can get particularly cumbersome when it comes to extracting and manipulating existing metadata in highly nested EML files.  The `emld` approach can leverage a rich array of tools for reading, extracting, and manipulating existing EML files. 

We can parse a simple example and manipulate is as a familiar list object (S3 object):

```{r}
f <- system.file("extdata/example.xml", package="emld")
eml <- as_emld(f)

cat(eml$dataset$title)
```


## Writing EML

Because `emld` objects are just nested lists, we can create EML just by writing lists: 

```{r}

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

Note that we don't have to worry about the order of the elements here, `as_xml` will re-order if necessary to validate. (For instance, in valid EML the `creator` becomes listed before `contact`. Of course this is a very low-level interface that does not help the user know what an EML looks like. Creating EML from scratch without knowledge of the schema is a job for the `EML` package and beyond the scope of the lightweight `emld`.  


# Working with EML as JSON-LD 

For many applications, it is useful to merely treat EML as a list object, as seen above, allowing the R user to leverage a standard tools and intuition in working with these files.  However, `emld` also opens the door to new possible directions by thinking of EML data in terms of a JSON-LD serialization rather than an XML serialization.  First, owing to it's comparative simplicity and native data typing (e.g. of Boolean/string/numeric data), JSON is often easier for many developers to work with than EML's native XML format.  


## As JSON: Query with JQ 

For example, JSON can be queried with with JQ, a [simple and powerful query language](https://stedolan.github.io/jq/manual/) that also gives us a lot of flexibility over the return structure of our results.  JQ syntax is both intuitive and well documented, and often easier than the typical munging of JSON/list data using `purrr`.  Here's an example query that turns EML to JSON and then extracts the north and south bounding coordinates:


```{r}
hf205 <- system.file("extdata/hf205.xml", package="emld")

as_emld(hf205) %>% 
  as_json() %>% 
  jq('.dataset.coverage.geographicCoverage.boundingCoordinates | 
       { northLat: .northBoundingCoordinate, 
         southLat: .southBoundingCoordinate }') %>%
  fromJSON()
```

Nice features of JQ include the ability to do recursive descent (common to XPATH but not possible in `purrr`) and specify the shape and names of the return object (e.g. as a list with elements named `northLat` and `southLat` in this case.)



## As semantic data: SPARQL queries


Another side-effect of the JSON-LD representation is that we can treat EML as *semantic* data.  This can provide a way to integrate EML records with other data sources, and means we can query the EML using semantic SPARQL queries.  One nice thing about SPARQL queries is that, in contrast to XPATH, JQ, or other graph queries, SPARQL always returns a `data.frame` -- a particularly convenient and familiar format for R users. 

First, we render the EML XML document into JSON-LD file:

```{r}
f <- system.file("extdata/hf205.xml", package="emld")
as_emld(f) %>%
  as_json("hf205.json")
```

We can now construct a SPARQL query. SPARQL queries look like SQL queries in that we name the columns we want with a `SELECT` command. Unlike SQL, SPARQL allows us to walk the graph by treating these names as variables, indicated by prefixing a `?` to the variable name.  We then use a WHERE block to define how these variables relate to each other.  In this case, we ask for the genus and species name and bounding box found in the EML file.  

```{r}
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

We can now use the `rdflib` library to execute this SPARQL query on the EML document and display the resulting `data.frame`:

```{r}
rdf <- rdf_parse("hf205.json", "jsonld")
df <- rdf_query(rdf, sparql)
df
```






```{r include=FALSE}
unlink("eml.xml")
unlink("test.xml")
unlink("ex.xml")
unlink("hf205.json")
```



# References
---
---


```{r}
library(tidyverse)
library(glue)
library(xml2)
library(jsonlite)
xsd_files <- list.files("../inst/xsd/eml-2.2.0", full.names = TRUE)
```

Find arbitrary element

```{r}
found <- map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:complexType[@name='ConferenceProceedings']")})

names(found) <- xsd_files
compact(found)
```


```{r}
found <- map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:complexContent/xs:restriction")})

names(found) <- xsd_files
compact(found)
```

# Explore

Root contains: 
 - `attribute`, `attributeGroup`, `complexType`, `element`, `group`, `simpleType`
 - `annotation`, `import` 
 
`element`: contains nothing (but has a type declaration), or has no type declaration and 
contains one of:
  - simpleType
  - complexType
 (plus every element has an annotation child element)

`complexType`: contains
  - `annotation`
  - `attribute`
  - `group`
  - `choice`
  - `sequence`
  - `complexContent`
  - `simpleContent`

`choice`, `sequence`, `group` all contain any number of:
- `choice`, `sequence`, `group`, `element`  

(note that `sequence` also appears in `restriction` & `extension`)
(groups only contain sequences)

- `attribute` occassionally uses `simpleType`.  `attributeGroup` has no child elements(?)
- Both `simpleContent` and `complexContent` contains one of either `xs:extension` or `xs:restriction`
- `simpleType` contains one of `restriction`, `union`, or `list`



```{r}
map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:element") %>% map(xml_children) %>% map(xml_name) })
map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:extension") %>% map(xml_children) %>% map(xml_name) }) %>% unlist() %>% table()
map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:sequence/xs:any") }) %>% compact() -> a
```



```{r}
map(xsd_files, function(xsd){
    read_xml(xsd) %>% xml_root() %>% xml_children() %>% map_chr(xml_name)
}) %>% unlist() %>% table()

terms <-
  map(xsd_files, function(xsd){
    read_xml(xsd) %>% xml_find_all("//xs:element") %>% xml_children() %>% map_chr(xml_name)
  }) %>% unlist()
types <- terms[!terms=="annotation"]
table(types)


map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:attribute") %>% xml_children() %>% map_chr(xml_name)
}) %>% unlist() %>% table()

map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:complexType") %>% xml_children() %>% map_chr(xml_name)
}) %>% unlist() %>% table()
map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:simpleType") %>% xml_children() %>% map_chr(xml_name)
}) %>% unlist() %>% table()
map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:choice") %>% xml_children() %>% map_chr(xml_name)
}) %>% unlist() %>% table()
map(xsd_files, function(xsd){
  read_xml(xsd) %>% xml_find_all("//xs:sequence") %>% xml_children() %>% map_chr(xml_name)
}) %>% unlist() %>% table()
```
---
output: github_document
---

# Mapping between Schema.Org and EML via JQ 



We will use J-Query as XSLT-like stylesheets. These two examples shipping in this package are currently a work in progress:

```{r}
library(jqr)
```


```{r}
eml_to_schema <- readr::read_file(system.file("jq/eml_to_schema.jq", package="emld"))
schema_to_eml <-  readr::read_file(system.file("jq/schema_to_eml.jq", package="emld"))

```


Let's map a more complete EML document into schema.org:

```{r}
eml <- readr::read_file("https://raw.githubusercontent.com/cboettig/emld/master/inst/extdata/hf205.json")
jq(eml, eml_to_schema)
```



Convert a dataset marked up in <http://schema.org/Dataset> terms into EML

```{r}
schema <- readr::read_file("../inst/extdata/schema-org-dataset.json")

jq(schema, schema_to_eml)
```


---
title: "R Notebook"
output: github_document
---


```{r message=FALSE}
library(jsonld)
library(jsonlite)
library(tidyverse)
library(printr)

```


```{r helper_fns, echo=FALSE}
is.url <- function(x){ 
  if(!is.character(x))
    return(FALSE)
  grepl("^http[s]?://.*", x)
}

add_context <- function(json, context = NULL){
  if(is.character(json)){ ## does not handle case of local files
    json <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    if(length(json) > 1)
      json <- list("@graph" = json)
  }
  if(is.url(context)){
     json[["@context"]] <- context[["@context"]]
  } else { 
    if(is.character(context))
      context <- jsonlite::fromJSON(context, simplifyVector = FALSE)
    json[["@context"]] <- context[["@context"]]
  }
  jsonlite::toJSON(json, pretty=TRUE, auto_unbox = TRUE)
}

drop_context <- function(json){
  if(is.character(json)){
    json <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  }
  json[["@context"]] <- NULL
  if(!is.null(json[["@graph"]]))
    json <- json[["@graph"]]
  jsonlite::toJSON(json, pretty=TRUE, auto_unbox = TRUE)

}

```

## Switching Contexts: Expansion & Compaction

- Good: homonyms (same word has different meanings), synonyms (different word has ~ same meaning)
- Fails: Mapping isn't 1:1

### A simple example

```{r}
ex <- 
 '{ 
      "date": "1993-10-02",
      "site": "N654",
      "scientificName": "Picea rubens",
      "area": 2,
      "count": 26
 }'
```

No need to alter your own data files/source data files.  Simply define a context that maps those terms into the reference vocabulary.  Data can use whatever fields it likes.


```{r}
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

### Using tabular schema

```{r}
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
Wow does that seem like a round-about way to just do:

```{r}
df %>% rename(scientificName = Species)
```

and it is.  But a few things are worth noting:

- Robust to what columns are actually present.  Instead of a default vocabulary we can map all terms explicitly.  
  - Terms not expanded by the map (that is, not part of our context) will be dropped. (`Count` in the example below)
  - Terms not in our target context remain as external keys, and are not compacted.  
- Typesafe transforms. Terms (keys/column names) with the wrong type will not compact if the types do not match the target transform.  This way they are not actually lost, but cannot be confused.  

```{r}
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


```{r}
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

----------




## Integration: Flattening shares id

```{r}
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


--------

hdf5: column data, with units, metadata
---
title: "emld Spec"
author: "Carl Boettiger"
date: "1/13/2018"
output: github_document
---

This document outlines the conventions used in mapping EML XML into JSON-LD.

- XML node names become JSON keys
- XML node contents (either text/data or node names) become the corresponding JSON value to the key
- Repeated elements become list-valued:

```xml
<url>http://example.com/1</url>
<url>http://example.com/2</url>
<url>http://example.com/3</url>
```

becomes:

```json
"url": ["http://example.com/1", "http://example.com/2", "http://example.com/3"]
```

- XML node attribute names are converted into JSON keys prefixed with `#`, unless the attribute name is `id`, in which case it uses the special JSON-LD prefix, `@`. For instance: 


```xml
<dataset id="abc123", system="knb"><title>The Title<title></dataset>
```

becomes

```json
"dataset": {
  "@id": "abc123", 
  "#system": "knb",
  "title": "The Title"
}
  
```


- If a node contains text/data and an attribute, the node name is repeated to denote the value as well:

```xml
<url function="download">http://example.com</url>
```

becomes:

```json
"url": {
  "#function"="download"
  "url": "http://example.com"
}
```

- Contents of `para` and `section` (the two choices for components of any `TextType` content) are encoded as literal character strings and not parsed.  

- references are expanded:

- Semantic annotations (EML 2.2) are expanded into native JSON-LD



---
title: "emld Tutorial"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{emld tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = (require("jqr") && require("rdflib"))
)
```


# emld

The goal of emld is to provide a way to work with EML metadata in the JSON-LD format. At it's heart, the package is simply a way to translate an EML XML document into JSON-LD and be able to reverse this so that any semantically equivalent JSON-LD file can be serialized into EML-schema valid XML.  The package has only three core functions:

- `as_emld()` Convert EML's `xml` files (or the `json` version created by this package) into a native R object (an S3 class called `emld`, essentially just a `list`).
- `as_xml()` Convert the native R format, `emld`, back into XML-schema valid EML.
- `as_json()` Convert the native R format, `emld`, into `json`(LD).


## Installation

You can install emld from github with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/emld")
```



## Motivation


In contrast to the existing [EML package](https://docs.ropensci.org/EML/), this package aims to a very light-weight implementation that seeks to provide both an intuitive data format and make maximum use of existing technology to work with that format.  In particular, this package emphasizes tools for working with linked data through the JSON-LD format.  This package is not meant to replace `EML`, as it does not support the more complex operations found in that package.  Rather, it provides a minimalist but powerful way of working with EML documents that can be used by itself or as a backend for those complex operations.  The next release of the EML R package will use `emld` under the hood.  

Note that the JSON-LD format is considerably less rigid than the EML schema.  This means that there are many valid, semantically equivalent representations on the JSON-LD side that must all map into the same or nearly the same XML format.  At the extreme end, the JSON-LD format can be serialized into RDF, where everything is flat set of triples (e.g. essentially a tabular representation), which we can query directly with semantic tools like SPARQL, and also automatically coerce back into the rigid nesting and ordering structure required by EML.  This ability to "flatten" EML files can be particularly convenient for applications consuming and parsing large numbers of EML files.  This package may also make it easier for other developers to build on the EML, since the S3/list and JSON formats used here have proven more appealing to many R developers than S4 and XML serializations.

JSON-LD also makes it easier to extend EML with other existing semantic vocabularies.  The standard JSON-LD operations (e.g. framing, compaction) make it easy for developers to specify desired data structures, filter unnecessary terms and provide defaults for needed ones, or even define custom property names, rather than working with the often cumbersome prefixes or URIs of linked data. 



```{r}
library(emld)
library(jsonlite)
library(magrittr)
```


## Reading EML

The `EML` package can get particularly cumbersome when it comes to extracting and manipulating existing metadata in highly nested EML files.  The `emld` approach can leverage a rich array of tools for reading, extracting, and manipulating existing EML files. 

We can parse a simple example and manipulate is as a familiar list object (S3 object):

```{r}
f <- system.file("extdata/example.xml", package="emld")
eml <- as_emld(f)
eml$dataset$title
```


## Writing EML

Because `emld` objects are just nested lists, we can create EML just by writing lists: 

```{r}

me <- list(individualName = list(givenName = "Carl", surName = "Boettiger"))

eml <- list(dataset = list(
              title = "dataset title",
              contact = me,
              creator = me),
              system = "doi",
              packageId = "10.xxx")

ex.xml <- tempfile("ex", fileext = ".xml") # use your preferred file path

as_xml(eml, ex.xml)
testthat::expect_true(eml_validate(ex.xml) )
```

Note that we don't have to worry about the order of the elements here, `as_xml` will re-order if necessary to validate. (For instance, in valid EML the `creator` becomes listed before `contact`.)   Of course this is a very low-level interface that does not help the user know what an EML looks like. Creating EML from scratch without knowledge of the schema is a job for the `EML` package and beyond the scope of the lightweight `emld`.  


# Working with EML as JSON-LD 

For many applications, it is useful to merely treat EML as a list object, as seen above, allowing the R user to leverage a standard tools and intuition in working with these files.  However, `emld` also opens the door to new possible directions by thinking of EML data in terms of a JSON-LD serialization rather than an XML serialization.  First, owing to it's comparative simplicity and native data typing (e.g. of Boolean/string/numeric data), JSON is often easier for many developers to work with than EML's native XML format.  


## As JSON: Query with JQ 

For example, JSON can be queried with with JQ, a [simple and powerful query language](https://stedolan.github.io/jq/manual/) that also gives us a lot of flexibility over the return structure of our results.  JQ syntax is both intuitive and well documented, and often easier than the typical munging of JSON/list data using `purrr`.  Here's an example query that turns EML to JSON and then extracts the north and south bounding coordinates:


```{r}

if(require(jqr) && require(magrittr)){
  
hf205 <- system.file("extdata/hf205.xml", package="emld")

as_emld(hf205) %>% 
  as_json() %>% 
  jq('.dataset.coverage.geographicCoverage.boundingCoordinates | 
       { northLat: .northBoundingCoordinate, 
         southLat: .southBoundingCoordinate }') %>%
  fromJSON()

}
```

Nice features of JQ include the ability to do recursive descent (common to XPATH but not possible in `purrr`) and specify the shape of the return object.  Some prototype examples of how we can use this to translate between EML and <https://schema.org/Dataset> representations of the same metadata can be found in <https://github.com/ropensci/emld/blob/master/notebook/jq_maps.md>



## As semantic data: SPARQL queries


Another side-effect of the JSON-LD representation is that we can treat EML as "semantic" data.  This can provide a way to integrate EML records with other data sources, and means we can query the EML using semantic SPARQL queries.  One nice thing about SPARQL queries is that, in contrast to XPATH, JQ, or other graph queries, SPARQL always returns a `data.frame` which is a particularly convenient format. SPARQL queries look like SQL queries in that we name the columns we want with a `SELECT` command.  Unlike SQL, these names are act as variables.  We then use a WHERE block to define how these variables relate to each other.  



```{r}
if(require(rdflib) && require(magrittr)){
  
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

}
```




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_validate.R
\name{eml_validate}
\alias{eml_validate}
\title{eml_validate}
\usage{
eml_validate(eml, encoding = "UTF-8", schema = NULL)
}
\arguments{
\item{eml}{file path, xml_document,}

\item{encoding}{optional encoding for files, default UTF-8.}

\item{schema}{path to schema}
}
\value{
Whether the document is valid (logical)
}
\description{
eml_validate processes an EML document using the XSD schema for the
appropriate version of EML and determines if the document is schema-valid
as defined by the XSD specification
}
\examples{
\donttest{

 f <- system.file("extdata", "example.xml", package = "emld")

 ## validate file directly from disk:
 eml_validate(f)

 ## validate an eml object:
 eml <- as_emld(f)
 eml_validate(eml)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_version.R
\name{eml_version}
\alias{eml_version}
\title{Set or check the EML version default}
\usage{
eml_version(version = getOption("emld_db", "eml-2.2.0"))
}
\arguments{
\item{version}{EML version, currently either eml-2.2.0 (current version), or
eml-2.1.1. The 'eml-' prefix can be omitted.}
}
\value{
returns the EML version string. As a side-effect, sets the
requested version as the default version by setting the \code{emld_db}
variable in \code{\link[=options]{options()}}.
}
\description{
Set or check the EML version default
}
\examples{
eml_version()
eml_version("2.1.1")
eml_version("eml-2.1.1")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emld.R
\docType{package}
\name{emld-package}
\alias{emld}
\alias{emld-package}
\title{emld: Ecological Metadata as Linked Data}
\description{
The goal of emld is to provide a way to work with EML metadata
in the JSON-LD format. At it's heart, the package is simply a
way to translate an EML XML document into JSON-LD and be able
to reverse this so that any semantically equivalent JSON-LD
file can be serialized into EML-schema valid XML.
}
\details{
The package has only three core functions:
\itemize{
\item \code{\link[=as_emld]{as_emld()}} Convert EML's \code{xml} files (or the \code{json} version created
by this package) into a native R object (an S3 class called \code{emld},
essentially just a \code{list}).
\item \code{\link[=as_xml]{as_xml()}} Convert the native R format, \code{emld}, back into
XML-schema valid EML.
\item \code{\link[=as_json]{as_json()}} Convert the native R format, \code{emld}, into \code{json}(LD).
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/emld/}
  \item \url{https://github.com/ropensci/emld}
  \item Report bugs at \url{https://github.com/ropensci/emld/issues}
}

}
\author{
\strong{Maintainer}: Carl Boettiger \email{cboettig@gmail.com} (\href{https://orcid.org/0000-0002-1642-628X}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Matthew B. Jones \email{jones@nceas.ucsb.edu} (\href{https://orcid.org/0000-0003-0077-4738}{ORCID}) [copyright holder]
  \item Bryce Mecum \email{mecum@nceas.ucsb.edu} (\href{https://orcid.org/0000-0002-0381-3766}{ORCID}) [copyright holder]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_version.R
\name{eml_ns}
\alias{eml_ns}
\title{Get the XML namespace for a version of EML}
\usage{
eml_ns(version = eml_version())
}
\arguments{
\item{version}{EML version, currently either eml-2.2.0 (current version) or
eml-2.1.1. Defaults to current version.}
}
\value{
returns the full XML namespace URI for the specified version of the
schema
}
\description{
Utility function for use when filling in \code{xmlns}, \code{schemaLocation}, or
\code{vocab} in various representations of EML. This is a little more future-proof
than keeping a dictionary for each version since this won't break on the next
release.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_json.R
\name{as_json}
\alias{as_json}
\title{Coerce an emld object into JSON}
\usage{
as_json(x, file = NULL)
}
\arguments{
\item{x}{an emld object}

\item{file}{optional path to write out to file.
Otherwise, defaults to NULL and will return a json object.}
}
\value{
a json object. Or if a file path is provided, the metadata
is written out in JSON file and the function returns \code{NULL} invisibly.
}
\description{
Coerce an emld object into JSON
}
\details{
Note: since emld list object maintains a 1:1 correspondence with JSON,
following the conventions of jsonlite, this function is basically trivial. The
only purpose is to default to auto_unbox = TRUE in serializing lists to JSON.
}
\examples{
f <- system.file("extdata/example.xml", package = "emld")
emld <- as_emld(f)
json <- as_json(emld)
## can also write a json file to disk:
json_file <- tempfile()
as_json(emld, json_file)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_validate.R
\name{guess_root_schema}
\alias{guess_root_schema}
\title{Find the root schema module and version}
\usage{
guess_root_schema(doc)
}
\arguments{
\item{doc}{An \code{xml_document}}
}
\value{
If found, a list with names 'version', 'module', and `namespace. If
not found, throws an error.
}
\description{
Find the root schema module and version
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_xml.R
\name{as_xml}
\alias{as_xml}
\title{Coerce an emld object into XML (EML's standard format)}
\usage{
as_xml(x, file = NULL, root = "eml", ns = "eml", schemaLocation = TRUE)
}
\arguments{
\item{x}{an emld object}

\item{file}{optional path to write out to file.
Otherwise, defaults to NULL and will return an xml_document object.}

\item{root}{name for the root node; default to 'eml'}

\item{ns}{namespace abbreviation on root node, default 'eml'}

\item{schemaLocation}{If not explicitly set on \code{x}, automatically set
\code{xsi:schemaLocation} based upon the root namespace (\code{TRUE}, default), do not
set a \code{xsi:schemaLocation} (\code{FALSE}), or set a specific \code{xsi:schemaLocation}
value (\code{"Your value here..."}). See Examples.}
}
\value{
a xml_document object. Or if a file path is provided, the metadata
is written out in XML file and the function returns \code{NULL} invisibly.
}
\description{
Coerce an emld object into XML (EML's standard format)
}
\details{
Unlike as_json, this function cannot rely on the existing
convention of serializing a list to xml, eg, as defined by xml2::as_xml_document()
Instead, this relies on a modified version, as_eml_document.  In addition
further steps must be taken when working with JSON-LD to deal with
different possible framings and namespaces from the JSON-LD context
element. Thus this \code{as_xml} function is particular to EML and \code{emld}
objects alone.
}
\examples{
f <- system.file("extdata/example.xml", package = "emld")
emld <- as_emld(f)
xml <- as_xml(emld)

## can also write directly to a file:
xml_file <- tempfile()
as_xml(emld, xml_file)

## if you don't want the `xsi:schemaLocation` attribute set
as_xml(emld, schemaLocation = FALSE)

## or if you want to set your own value
as_xml(emld, schemaLocation = "https://eml.ecoinformatics.org/eml-2.2.0
http://example.com/eml-2.2.0/eml.xsd")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_emld.R
\name{as_emld}
\alias{as_emld}
\title{Coerce an EML file or object into an emld object.}
\usage{
as_emld(x, from = c("guess", "xml", "json", "list"))
}
\arguments{
\item{x}{path to an EML file}

\item{from}{explicit type for the input format. By default, will
attempt to guess the format, but it always safer to specify the
input format. This is essential for literal text strings or raw
vectors where the type cannot be guessed by the R object class
or file extension of the input.}
}
\value{
an emld object
}
\description{
Coerce an EML file or object into an emld object.
}
\examples{
 hf205 <- system.file("extdata/hf205.xml", package="emld")
 as_emld(hf205)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_validate.R
\name{find_real_root_name}
\alias{find_real_root_name}
\title{Get the real \code{QName} for the root element, including its prefix}
\usage{
find_real_root_name(doc)
}
\arguments{
\item{doc}{An \code{xml_document}}
}
\value{
A \code{list} with elements \code{prefix} and \code{name}. \code{prefix} will be \code{NULL}
if the element has no namespace prefix but \code{name} will always be a
\code{character}.
}
\description{
Note that if a default namespace is used, the prefix will be \code{d1}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/template.R
\name{template}
\alias{template}
\title{Create a template for an EML object}
\usage{
template(object)
}
\arguments{
\item{object}{the name of an eml object to create}
}
\value{
a list with elements named according to the properties of the object.
This can be coerced into EML, see vignettes. NULL-valued elements (~)
can take a data entry directly, while empty list()-valued elements ({})
indicate properties that take other eml objects as values.
}
\description{
Create a template for an EML object
}
\details{
Note: while this function can be called in recursions, doing so may be a bad idea.
}
\examples{
template("creator")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_version.R
\name{guess_schema_location}
\alias{guess_schema_location}
\title{Guess an appropriate \code{schemaLocation} value for a given version of the schema}
\usage{
guess_schema_location(version = eml_version())
}
\arguments{
\item{version}{Optional. Override the version of the schema. Defaults to the
current version returned by \code{eml_version}. See \code{eml_version} for information
on how to change the current version.}
}
\value{
Returns a string suitable as a value for \code{schemaLocation} or \code{NULL}
if a value wasn't found.
}
\description{
This is a simple helper to make filling in the \code{schemaLocation} attribute
on documents this package creates. Supports EML 2.1.1 and newer.
}
\examples{
\dontrun{
# Get an appropriate schemaLocation value for the current version fo EML
guess_schema_location()

# Get an appropriate value for EML 2.1.1
guess_schema_location("eml-2.1.1")
}
}
