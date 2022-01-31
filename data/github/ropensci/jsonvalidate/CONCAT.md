# jsonvalidate

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/ropensci/jsonvalidate.svg?branch=master)](https://travis-ci.org/ropensci/jsonvalidate)
[![codecov.io](https://codecov.io/github/ropensci/jsonvalidate/coverage.svg?branch=master)](https://codecov.io/github/ropensci/jsonvalidate?branch=master)
[![](http://www.r-pkg.org/badges/version/jsonvalidate)](https://cran.r-project.org/package=jsonvalidate)
<!-- badges: end -->


Validate JSON against a schema using [`is-my-json-valid`](https://github.com/mafintosh/is-my-json-valid) or [`ajv`](https://github.com/ajv-validator/ajv).  This package is a thin wrapper around these node libraries, using the [V8](https://cran.r-project.org/package=V8) package.

## Usage

Directly validate `json` against `schema`

```r
jsonvalidate::json_validate(json, schema)
```

or create a validator for multiple uses

```r
validate <- jsonvalidate::json_validator(schema)
validate(json)
validate(json2) # etc
```

See the [package vignette](https://docs.ropensci.org/jsonvalidate/articles/jsonvalidate.html) for complete examples.

## Installation

Install from CRAN with

```r
install.packages("jsonvalidate")
```

Alternatively, the current development version can be installed from GitHub with

```r
devtools::install_github("ropensci/jsonvalidate")
```

## License

MIT + file LICENSE © [Rich FitzJohn](https://github.com/richfitz).

 Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/jsonvalidate/blob/master/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
# jsonvalidate 1.4.0

* Support for safely serialising objects to json, guided by the schema, with new function `json_serialise`
* New object `json_schema` for construction of reusable validation and serialisation functions

# jsonvalidate 1.3.2

* Always uses ES5 version of Ajv, which allows use in both current and "legacy" V8 (#51)

# jsonvalidate 1.3.0

* Upgrade to ajv version 8.5.0
* Add arg `strict` to `json_validate` and `json_validator` to allow evaluating schema in strict mode for ajv only. This is off (`FALSE`) by default to use permissive behaviour detailed in JSON schema

# jsonvalidate 1.2.3

* Schemas can use references to other files with JSON pointers i.e. schemas can reference parts of other files e.g. `definitions.json#/definitions/hello`
* JSON can be validated against a subschema (#18, #19, @AliciaSchep)
* Validation with `error = TRUE` now returns `TRUE` (not `NULL)` on success
* Schemas can span multiple files, being included via `"$ref": "filename.json"` - supported with the ajv engine only (#20, #21, @r-ash).
* Validation can be performed against a fraction of the input data (#25)

# jsonvalidate 1.1.0

* Add support for JSON schema draft 06 and 07 using the [`ajv`](https://github.com/ajv-validator/ajv) node library.  This must be used by passing the `engine` argument to `json_validate` and `json_validator` at present (#2, #11, #15, #16, #17, @karawoo & @ijlyttle)

# jsonvalidate 1.0.1

* Initial CRAN release
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or imagery, derogatory comments or personal attacks, trolling, public or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant (http://contributor-covenant.org), version 1.0.0, available at http://contributor-covenant.org/version/1/0/0/
# docker/es5 support

This dockerfile exists to make it easier to test that things still work in an es5 environment, as that is still fairly common (notably Solaris, but also ubuntu 18.04 LTS and RHEL 7).

From the root directory of jsonvalidate source, run

```
docker build --tag richfitz/jsonvalidate:es5 docker
```

to build the image; this should not take that long. It installs the current r-release along with the CRAN versions of jsonvalidate and testthat (which ensures all core dependencies are present).

Once setup, you can bring up a container with:

```
docker run --rm -it -v $PWD:/src:ro richfitz/jsonvalidate:es5 bash
```

which mounts the current directory read-only into the container at `/src`.  That version of the source (rather than the CRAN one installed in the base image) can be installed with `R CMD INSTALL /src`

To run the whole test suite run:

```
Rscript -e 'testthat::test_local("/src")'
```

More simply, to just confirm that the bundle is valid you can do

```
Rscript -e 'V8::new_context()$source("/src/inst/bundle.js")'
```

which will error if the bundle is invalid.

To do a full reverse dependencies check with old libv8, you can bring up R in this container:

```
docker run --rm -it -v $PWD:/src:ro richfitz/jsonvalidate:es5 bash
```

Then install revdepcheck itself

```
install.packages("remotes")
remotes::install_github("r-lib/revdepcheck", upgrade = TRUE)
```

Additional packages that we need, for some reason these did not get installed automatically though we'd have expected them to (see [this issue](https://github.com/r-lib/revdepcheck/issues/209))

```
install.packages(c(
  "cinterpolate",
  "deSolve",
  "devtools",
  "golem",
  "inTextSummaryTable",
  "patientProfilesVis",
  "reticulate",
  "rlist",
  "shiny",
  "shinyBS",
  "shinydashboard",
  "shinyjs",
  "tableschema.r",
  "xml2"))
```

At this point you will need to cycle the R session because the package DB will be corrupted by all the installations.

Finally we can run the reverse dependency check:

```
unlink("/tmp/src", recursive = TRUE)
file.copy("/src", "/tmp", recursive = TRUE)
revdepcheck::revdep_check("/tmp/src", num_workers = 4)
```
---
title: "Introduction to jsonvalidate"
author: "Rich FitzJohn"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to jsonvalidate}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo = FALSE, results = "hide"}
knitr::opts_chunk$set(error = FALSE)
```

This package wraps
[is-my-json-valid](https://github.com/mafintosh/is-my-json-valid)
using [V8](https://cran.r-project.org/package=V8) to do JSON schema
validation in R.

You need a JSON schema file; see
[json-schema.org](http://json-schema.org) for details on writing
these.  Often someone else has done the hard work of writing one
for you, and you can just check that the JSON you are producing or
consuming conforms to the schema.

The examples below come from the [JSON schema
website](http://json-schema.org/learn/getting-started-step-by-step.html)

They describe a JSON based product catalogue, where each product
has an id, a name, a price, and an optional set of tags.  A JSON
representation of a product is:

```json
{
    "id": 1,
    "name": "A green door",
    "price": 12.50,
    "tags": ["home", "green"]
}
```

The schema that they derive looks like this:

```json
{
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Product",
    "description": "A product from Acme's catalog",
    "type": "object",
    "properties": {
        "id": {
            "description": "The unique identifier for a product",
            "type": "integer"
        },
        "name": {
            "description": "Name of the product",
            "type": "string"
        },
        "price": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": true
        },
        "tags": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "minItems": 1,
            "uniqueItems": true
        }
    },
    "required": ["id", "name", "price"]
}
```

This ensures the types of all fields, enforces presence of `id`,
`name` and `price`, checks that the price is not negative and
checks that if present `tags` is a unique list of strings.

There are two ways of passing the schema in to R; as a string or as
a filename.  If you have a large schema loading as a file will
generally be easiest!  Here's a string representing the schema
(watch out for escaping quotes):

```{r}
schema <- '{
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Product",
    "description": "A product from Acme\'s catalog",
    "type": "object",
    "properties": {
        "id": {
            "description": "The unique identifier for a product",
            "type": "integer"
        },
        "name": {
            "description": "Name of the product",
            "type": "string"
        },
        "price": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": true
        },
        "tags": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "minItems": 1,
            "uniqueItems": true
        }
    },
    "required": ["id", "name", "price"]
}'
```

Create a schema object, which can be used to validate a schema:

```{r}
obj <- jsonvalidate::json_schema$new(schema)
```

If we'd saved the json to a file, this would work too:

```{r}
path <- tempfile()
writeLines(schema, path)
obj <- jsonvalidate::json_schema$new(path)
```

```{r include = FALSE}
file.remove(path)
```

The returned object is a function that takes as its first argument
a json string, or a filename of a json file.  The empty list will
fail validation because it does not contain any of the required fields:

```{r}
obj$validate("{}")
```

To get more information on why the validation fails, add `verbose = TRUE`:

```{r}
obj$validate("{}", verbose = TRUE)
```

The attribute "errors" is a data.frame and is present only when the
json fails validation.  The error messages come straight from
`ajv` and they may not always be that informative.

Alternatively, to throw an error if the json does not validate, add
`error = TRUE` to the call:

```{r error = TRUE}
obj$validate("{}", error = TRUE)
```

The JSON from the opening example works:

```{r}
obj$validate('{
    "id": 1,
    "name": "A green door",
    "price": 12.50,
    "tags": ["home", "green"]
}')
```

But if we tried to enter a negative price it would fail:

```{r}
obj$validate('{
    "id": 1,
    "name": "A green door",
    "price": -1,
    "tags": ["home", "green"]
}', verbose = TRUE)
```

...or duplicate tags:

```{r}
obj$validate('{
    "id": 1,
    "name": "A green door",
    "price": 12.50,
    "tags": ["home", "home"]
}', verbose = TRUE)
```

or just basically everything wrong:
```{r}
obj$validate('{
    "id": "identifier",
    "name": 1,
    "price": -1,
    "tags": ["home", "home", 1]
}', verbose = TRUE)
```

The names comes from within the `ajv` source, and may be annoying to work with programmatically.

There is also a simple interface where you take the schema and the
json at the same time:

```{r}
json <- '{
    "id": 1,
    "name": "A green door",
    "price": 12.50,
    "tags": ["home", "green"]
}'
jsonvalidate::json_validate(json, schema, engine = "ajv")
```

However, this will be much slower than building the schema object once and using it repeatedly.

Prior to 1.4.0, the recommended way of building a reusable validator object was to use `jsonvalidate::json_validator`; this is still supported but note that it has different defaults to `jsonvalidate::json_schema` (using imjv for backward compatibility).

```{r}
v <- jsonvalidate::json_validator(schema, engine = "ajv")
v(json)
```

While we do not intend on removing this old interface, new code should prefer both `jsonvalidate::json_schema` and the `ajv` engine.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.R
\name{json_validator}
\alias{json_validator}
\title{Create a json validator}
\usage{
json_validator(schema, engine = "imjv", reference = NULL, strict = FALSE)
}
\arguments{
\item{schema}{Contents of the json schema, or a filename
containing a schema.}

\item{engine}{Specify the validation engine to use.  Options are
"imjv" (the default; which uses "is-my-json-valid") and "ajv"
(Another JSON Schema Validator).  The latter supports more
recent json schema features.}

\item{reference}{Reference within schema to use for validating against a
sub-schema instead of the full schema passed in. For example
if the schema has a 'definitions' list including a definition for a
'Hello' object, one could pass "#/definitions/Hello" and the validator
would check that the json is a valid "Hello" object. Only available if
\code{engine = "ajv"}.}

\item{strict}{Set whether the schema should be parsed strictly or not.
If in strict mode schemas will error to "prevent any unexpected
behaviours or silently ignored mistakes in user schema". For example
it will error if encounters unknown formats or unknown keywords. See
https://ajv.js.org/strict-mode.html for details. Only available in
\code{engine = "ajv"}.}
}
\value{
A function that can be used to validate a
schema. Additionally, the function has two attributes assigned:
\code{v8} which is the javascript context (used internally) and
\code{engine}, which contains the name of the engine used.
}
\description{
Create a validator that can validate multiple json files.
}
\section{Validation Engines}{


We support two different json validation engines, \code{imjv}
("is-my-json-valid") and \code{ajv} ("Another JSON
Validator"). \code{imjv} was the original validator included in
the package and remains the default for reasons of backward
compatibility. However, users are encouraged to migrate to
\code{ajv} as with it we support many more features, including
nested schemas that span multiple files, meta schema versions
later than draft-04, validating using a subschema, and
validating a subset of an input data object.

If your schema uses these features we will print a message to
screen indicating that you should update when running
interactively. We do not use a warning here as this will be
disruptive to users. You can disable the message by setting the
option \code{jsonvalidate.no_note_imjv} to \code{TRUE}. Consider using
\code{\link[withr:with_options]{withr::with_options()}} (or simply \code{\link[=suppressMessages]{suppressMessages()}}) to
scope this option if you want to quieten it within code you do
not control.  Alternatively, setting the option
\code{jsonvalidate.no_note_imjv} to \code{FALSE} will print the message
even noninteractively.

Updating the engine should be simply a case of adding \verb{\{engine = "ajv"} to your \code{json_validator} or \code{json_validate}
calls, but you may see some issues when doing so.
\itemize{
\item Your json now fails validation: We've seen this where schemas
spanned several files and are silently ignored. By including
these, your data may now fail validation and you will need to
either fix the data or the schema.
\item Your code depended on the exact payload returned by \code{imjv}: If
you are inspecting the error result and checking numbers of
errors, or even the columns used to describe the errors, you
will likely need to update your code to accommodate the slightly
different format of \code{ajv}
\item Your schema is simply invalid: If you reference an invalid
metaschema for example, jsonvalidate will fail
}
}

\section{Using multiple files}{


Multiple files are supported.  You can have a schema that
references a file \code{child.json} using \code{{"$ref": "child.json"}} -
in this case if \code{child.json} includes an \code{id} or \verb{$id} element
it will be silently dropped and the filename used to reference
the schema will be used as the schema id.

The support is currently quite limited - it will not (yet) read
sub-child schemas relative to child schema \verb{$id} url, and
does not support reading from URLs (only local files are
supported).
}

\examples{
# A simple schema example:
schema <- '{
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Product",
    "description": "A product from Acme\'s catalog",
    "type": "object",
    "properties": {
        "id": {
            "description": "The unique identifier for a product",
            "type": "integer"
        },
        "name": {
            "description": "Name of the product",
            "type": "string"
        },
        "price": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": true
        },
        "tags": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "minItems": 1,
            "uniqueItems": true
        }
    },
    "required": ["id", "name", "price"]
}'

# Create a validator function
v <- jsonvalidate::json_validator(schema)

# Test if some (invalid) json conforms to the schema
v("{}", verbose = TRUE)

# Test if some (valid) json conforms to the schema
v('{
    "id": 1,
    "name": "A green door",
    "price": 12.50,
    "tags": ["home", "green"]
}')

# Using features from draft-06 or draft-07 requires the ajv engine:
schema <- "{
  '$schema': 'http://json-schema.org/draft-06/schema#',
  'type': 'object',
  'properties': {
    'a': {
      'const': 'foo'
    }
  }
}"

# Create the validator
v <- jsonvalidate::json_validator(schema, engine = "ajv")

# This confirms to the schema
v('{"a": "foo"}')

# But this does not
v('{"a": "bar"}')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/schema.R
\name{json_schema}
\alias{json_schema}
\title{Interact with JSON schemas}
\description{
Interact with JSON schemas, using them to validate
json strings or serialise objects to JSON safely.

This interface supercedes \link{json_schema} and changes
some default arguments.  While the old interface is not going
away any time soon, users are encouraged to switch to this
interface, which is what we will develop in the future.
}
\examples{
# This is the schema from ?json_validator
schema <- '{
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Product",
    "description": "A product from Acme\'s catalog",
    "type": "object",
    "properties": {
        "id": {
            "description": "The unique identifier for a product",
            "type": "integer"
        },
        "name": {
            "description": "Name of the product",
            "type": "string"
        },
        "price": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": true
        },
        "tags": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "minItems": 1,
            "uniqueItems": true
        }
    },
    "required": ["id", "name", "price"]
}'

# We're going to use a validator object below
v <- jsonvalidate::json_validator(schema, "ajv")

# And this is some data that we might generate in R that we want to
# serialise using that schema
x <- list(id = 1, name = "apple", price = 0.50, tags = "fruit")

# If we serialise to json, then 'id', 'name' and "price' end up a
# length 1-arrays
jsonlite::toJSON(x)

# ...and that fails validation
v(jsonlite::toJSON(x))

# If we auto-unbox then 'fruit' ends up as a string and not an array,
# also failing validation:
jsonlite::toJSON(x, auto_unbox = TRUE)
v(jsonlite::toJSON(x, auto_unbox = TRUE))

# Using json_serialise we can guide the serialisation process using
# the schema:
jsonvalidate::json_serialise(x, schema)

# ...and this way we do pass validation:
v(jsonvalidate::json_serialise(x, schema))

# It is typically much more efficient to construct a json_schema
# object first and do both operations with it:
obj <- jsonvalidate::json_schema$new(schema)
json <- obj$serialise(x)
obj$validate(json)
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{schema}}{The parsed schema, cannot be rebound}

\item{\code{engine}}{The name of the schema validation engine}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{json_schema$new()}}
\item \href{#method-validate}{\code{json_schema$validate()}}
\item \href{#method-serialise}{\code{json_schema$serialise()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{json_schema} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{json_schema$new(schema, engine = "ajv", reference = NULL, strict = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{schema}}{Contents of the json schema, or a filename
containing a schema.}

\item{\code{engine}}{Specify the validation engine to use.  Options are
"ajv" (the default; "Another JSON Schema Validator") or "imjv"
("is-my-json-valid", the default everywhere in versions prior
to 1.4.0, and the default for \link{json_validator}.
\emph{Use of \code{ajv} is strongly recommended for all new code}.}

\item{\code{reference}}{Reference within schema to use for validating
against a sub-schema instead of the full schema passed in.
For example if the schema has a 'definitions' list including a
definition for a 'Hello' object, one could pass
"#/definitions/Hello" and the validator would check that the json
is a valid "Hello" object. Only available if \code{engine = "ajv"}.}

\item{\code{strict}}{Set whether the schema should be parsed strictly or not.
If in strict mode schemas will error to "prevent any unexpected
behaviours or silently ignored mistakes in user schema". For example
it will error if encounters unknown formats or unknown keywords. See
https://ajv.js.org/strict-mode.html for details. Only available in
\code{engine = "ajv"} and silently ignored for "imjv".
Validate a json string against a schema.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validate"></a>}}
\if{latex}{\out{\hypertarget{method-validate}{}}}
\subsection{Method \code{validate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{json_schema$validate(
  json,
  verbose = FALSE,
  greedy = FALSE,
  error = FALSE,
  query = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{json}}{Contents of a json object, or a filename containing
one.}

\item{\code{verbose}}{Be verbose?  If \code{TRUE}, then an attribute
"errors" will list validation failures as a data.frame}

\item{\code{greedy}}{Continue after the first error?}

\item{\code{error}}{Throw an error on parse failure?  If \code{TRUE},
then the function returns \code{NULL} on success (i.e., call
only for the side-effect of an error on failure, like
\code{stopifnot}).}

\item{\code{query}}{A string indicating a component of the data to
validate the schema against.  Eventually this may support full
\href{https://www.npmjs.com/package/jsonpath}{jsonpath} syntax, but
for now this must be the name of an element within \code{json}.  See
the examples for more details.
Serialise an R object to JSON with unboxing guided by the schema.
See \link{json_serialise} for details on the problem and
the algorithm.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialise"></a>}}
\if{latex}{\out{\hypertarget{method-serialise}{}}}
\subsection{Method \code{serialise()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{json_schema$serialise(object)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{object}}{An R object to serialise}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serialise.R
\name{json_serialise}
\alias{json_serialise}
\title{Safe JSON serialisation}
\usage{
json_serialise(
  object,
  schema,
  engine = "ajv",
  reference = NULL,
  strict = FALSE
)
}
\arguments{
\item{object}{An object to be serialised}

\item{schema}{A schema (string or path to a string, suitable to be
passed through to \link{json_validator} or a validator
object itself.}

\item{engine}{The engine to use. Only ajv is supported, and trying
to use \code{imjv} will throw an error.}

\item{reference}{Reference within schema to use for validating against a
sub-schema instead of the full schema passed in. For example
if the schema has a 'definitions' list including a definition for a
'Hello' object, one could pass "#/definitions/Hello" and the validator
would check that the json is a valid "Hello" object. Only available if
\code{engine = "ajv"}.}

\item{strict}{Set whether the schema should be parsed strictly or not.
If in strict mode schemas will error to "prevent any unexpected
behaviours or silently ignored mistakes in user schema". For example
it will error if encounters unknown formats or unknown keywords. See
https://ajv.js.org/strict-mode.html for details. Only available in
\code{engine = "ajv"}.}
}
\value{
A string, representing \code{object} in JSON format. As for
\code{jsonlite::toJSON} we set the class attribute to be \code{json} to
mark it as serialised json.
}
\description{
Safe serialisation of json with unboxing guided by the schema.
}
\details{
When using \link[jsonlite:fromJSON]{jsonlite::toJSON} we are forced to deal with the
differences between R's types and those available in JSON. In
particular:
\itemize{
\item R has no scalar types so it is not clear if \code{1} should be
serialised as a number or a vector of length 1; jsonlite
provides support for "automatically unboxing" such values
(assuming that length-1 vectors are scalars) or never unboxing
them unless asked to using \link[jsonlite:unbox]{jsonlite::unbox}
\item JSON has no date/time values and there are many possible string
representations.
\item JSON has no \link{data.frame} or \link{matrix} type and there are several
ways of representing these in JSON, all equally valid (e.g., row-wise,
column-wise or as an array of objects).
\item The handling of \code{NULL} and missing values (\code{NA}, \code{NaN}) are different
\item We need to chose the number of digits to write numbers out at,
balancing precision and storage.
}

These issues are somewhat lessened when we have a schema because
we know what our target type looks like.  This function attempts
to use the schema to guide serialsation of json safely.  Currently
it only supports detecting the appropriate treatment of length-1
vectors, but we will expand functionality over time.

For a user, this function provides an argument-free replacement
for \code{jsonlite::toJSON}, accepting an R object and returning a
string with the JSON representation of the object. Internally the
algorithm is:
\enumerate{
\item serialise the object with \link[jsonlite:fromJSON]{jsonlite::toJSON}, with
\code{auto_unbox = FALSE} so that length-1 vectors are serialised as a
length-1 arrays.
\item operating entirely within JavaScript, deserialise the object
with \code{JSON.parse}, traverse the object and its schema
simultaneously looking for length-1 arrays where the schema
says there should be scalar value and unboxing these, and
re-serialise with \code{JSON.stringify}
}

There are several limitations to our current approach, and not all
unboxable values will be found - at the moment we know that
schemas contained within a \code{oneOf} block (or similar) will not be
recursed into.
}
\section{}{
 Warning:

Direct use of this function will be slow!  If you are going to
serialise more than one or two objects with a single schema, you
should use the \code{serialise} method of a
\link{json_schema} object which you create once and pass around.
}

\examples{
# This is the schema from ?json_validator
schema <- '{
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Product",
    "description": "A product from Acme\'s catalog",
    "type": "object",
    "properties": {
        "id": {
            "description": "The unique identifier for a product",
            "type": "integer"
        },
        "name": {
            "description": "Name of the product",
            "type": "string"
        },
        "price": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": true
        },
        "tags": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "minItems": 1,
            "uniqueItems": true
        }
    },
    "required": ["id", "name", "price"]
}'

# We're going to use a validator object below
v <- jsonvalidate::json_validator(schema, "ajv")

# And this is some data that we might generate in R that we want to
# serialise using that schema
x <- list(id = 1, name = "apple", price = 0.50, tags = "fruit")

# If we serialise to json, then 'id', 'name' and "price' end up a
# length 1-arrays
jsonlite::toJSON(x)

# ...and that fails validation
v(jsonlite::toJSON(x))

# If we auto-unbox then 'fruit' ends up as a string and not an array,
# also failing validation:
jsonlite::toJSON(x, auto_unbox = TRUE)
v(jsonlite::toJSON(x, auto_unbox = TRUE))

# Using json_serialise we can guide the serialisation process using
# the schema:
jsonvalidate::json_serialise(x, schema)

# ...and this way we do pass validation:
v(jsonvalidate::json_serialise(x, schema))

# It is typically much more efficient to construct a json_schema
# object first and do both operations with it:
obj <- jsonvalidate::json_schema$new(schema)
json <- obj$serialise(x)
obj$validate(json)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.R
\name{json_validate}
\alias{json_validate}
\title{Validate a json file}
\usage{
json_validate(
  json,
  schema,
  verbose = FALSE,
  greedy = FALSE,
  error = FALSE,
  engine = "imjv",
  reference = NULL,
  query = NULL,
  strict = FALSE
)
}
\arguments{
\item{json}{Contents of a json object, or a filename containing
one.}

\item{schema}{Contents of the json schema, or a filename
containing a schema.}

\item{verbose}{Be verbose?  If \code{TRUE}, then an attribute
"errors" will list validation failures as a data.frame}

\item{greedy}{Continue after the first error?}

\item{error}{Throw an error on parse failure?  If \code{TRUE},
then the function returns \code{NULL} on success (i.e., call
only for the side-effect of an error on failure, like
\code{stopifnot}).}

\item{engine}{Specify the validation engine to use.  Options are
"imjv" (the default; which uses "is-my-json-valid") and "ajv"
(Another JSON Schema Validator).  The latter supports more
recent json schema features.}

\item{reference}{Reference within schema to use for validating against a
sub-schema instead of the full schema passed in. For example
if the schema has a 'definitions' list including a definition for a
'Hello' object, one could pass "#/definitions/Hello" and the validator
would check that the json is a valid "Hello" object. Only available if
\code{engine = "ajv"}.}

\item{query}{A string indicating a component of the data to
validate the schema against.  Eventually this may support full
\href{https://www.npmjs.com/package/jsonpath}{jsonpath} syntax, but
for now this must be the name of an element within \code{json}.  See
the examples for more details.}

\item{strict}{Set whether the schema should be parsed strictly or not.
If in strict mode schemas will error to "prevent any unexpected
behaviours or silently ignored mistakes in user schema". For example
it will error if encounters unknown formats or unknown keywords. See
https://ajv.js.org/strict-mode.html for details. Only available in
\code{engine = "ajv"}.}
}
\description{
Validate a single json against a schema.  This is a convenience
wrapper around \code{json_validator(schema)(json)} or
\code{json_schema$new(schema, engine = "ajv")$validate(json)}.  See
\code{\link[=json_validator]{json_validator()}} for further details.
}
\examples{
# A simple schema example:
schema <- '{
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Product",
    "description": "A product from Acme\'s catalog",
    "type": "object",
    "properties": {
        "id": {
            "description": "The unique identifier for a product",
            "type": "integer"
        },
        "name": {
            "description": "Name of the product",
            "type": "string"
        },
        "price": {
            "type": "number",
            "minimum": 0,
            "exclusiveMinimum": true
        },
        "tags": {
            "type": "array",
            "items": {
                "type": "string"
            },
            "minItems": 1,
            "uniqueItems": true
        }
    },
    "required": ["id", "name", "price"]
}'

# Test if some (invalid) json conforms to the schema
jsonvalidate::json_validate("{}", schema, verbose = TRUE)

# Test if some (valid) json conforms to the schema
json <- '{
    "id": 1,
    "name": "A green door",
    "price": 12.50,
    "tags": ["home", "green"]
}'
jsonvalidate::json_validate(json, schema)

# Test a fraction of a data against a reference into the schema:
jsonvalidate::json_validate(json, schema,
                            query = "tags", reference = "#/properties/tags",
                            engine = "ajv", verbose = TRUE)
}
