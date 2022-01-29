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
