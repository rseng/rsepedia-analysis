<!-- README.md is generated from README.Rmd. Please edit that file and knit -->



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/staypuft/workflows/R-check/badge.svg)](https://github.com/ropensci/staypuft/actions/)

`staypuft` is a port of Python's [marshmallow][] for converting objects to and from R data structures 

Main `Schema` methods:
- `load`: 'deserialize', or validate and deserialize an input R data structure (e.g., list) to an object
- `dump`: 'serialize', or convert any input (e.g., R6 class) to R data structures (e.g., list)
- `load_json`: same as `load`, but accepts JSON
- `dump_json`: same as `dump`, but returns JSON

## Installation


```r
remotes::install_github("ropensci/staypuft")
```


```r
library("staypuft")
```

## example


```r
z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character(),
  num = fields$integer()
)
z
#> <schema: MySchema>
#> fields: name, title, num
x <- list(name = "Jane Doe", title = "Howdy doody", num = 5.5)
z$load(data = x)
#> $name
#> [1] "Jane Doe"
#> 
#> $title
#> [1] "Howdy doody"
#> 
#> $num
#> [1] 5.5
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/staypuft/issues).
* License: MIT
* Get citation information for `staypuft` in R doing `citation(package = 'staypuft')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)


[marshmallow]: https://github.com/marshmallow-code/marshmallow/
[coc]: https://github.com/ropensci/staypuft/blob/master/CODE_OF_CONDUCT.md
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
(https://www.contributor-covenant.org/), version 1.0.0, available at 
https://www.contributor-covenant.org/version/1/0/0/code-of-conduct.html
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

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

# CONTRIBUTING #

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
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html) for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
<!-- README.md is generated from README.Rmd. Please edit that file and knit -->

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/staypuft/workflows/R-check/badge.svg)](https://github.com/ropensci/staypuft/actions/)

`staypuft` is a port of Python's [marshmallow][] for converting objects to and from R data structures 

Main `Schema` methods:
- `load`: 'deserialize', or validate and deserialize an input R data structure (e.g., list) to an object
- `dump`: 'serialize', or convert any input (e.g., R6 class) to R data structures (e.g., list)
- `load_json`: same as `load`, but accepts JSON
- `dump_json`: same as `dump`, but returns JSON

## Installation

```{r eval=FALSE}
remotes::install_github("ropensci/staypuft")
```

```{r}
library("staypuft")
```

## example

```{r}
z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character(),
  num = fields$integer()
)
z
x <- list(name = "Jane Doe", title = "Howdy doody", num = 5.5)
z$load(data = x)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/staypuft/issues).
* License: MIT
* Get citation information for `staypuft` in R doing `citation(package = 'staypuft')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)


[marshmallow]: https://github.com/marshmallow-code/marshmallow/
[coc]: https://github.com/ropensci/staypuft/blob/master/CODE_OF_CONDUCT.md
---
title: "Introduction to staypuft"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to staypuft}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Installation

```{r eval=FALSE}
remotes::install_github("ropensci/staypuft")
```

```{r}
library("staypuft")
```

## hello world

```{r}
z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character(),
  num = fields$integer()
)
z
x <- list(name = "Jane Doe", title = "Howdy doody", num = 5.5)
z$load(data = x)
z$load(data = x, as_df = TRUE)
z$load_json(jsonlite::toJSON(x, auto_unbox=TRUE))
```

strict mode for integer

```{r error=TRUE}
z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character(),
  num = fields$integer(strict = TRUE)
)
z$fields$num
x <- list(name = "Jane Doe", title = "Howdy doody", num = 5.5)
z$load(data = x)
```

another example

```{r error=TRUE}
z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character(),
  num = fields$integer(),
  uuid = fields$uuid(),
  date = fields$date(),
  foo = fields$boolean()
)
x <- list(name = "Jane Doe", title = "Howdy doody", num = 5.5, 
    uuid = "9a5f6bba-4101-48e9-a7e3-b5ac456a04b5", date = "2020/06/16",
    foo = TRUE)
z$load(data = x)

# invalid uuid
x$uuid <- "foo-bar"
z$load(data = x)

# invalid boolean
x$uuid <- "9a5f6bba-4101-48e9-a7e3-b5ac456a04b5"
x$foo <- "bar"
z$load(data = x)
```
---
title: "Use cases"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Use cases}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Validate JSON

```{r}
library(staypuft)
```

An example `package.json` for a Python package

```{r}
json <- '{
  "name": "dunderscore",
  "version": "1.2.3",
  "description": "The Pythonic JavaScript toolkit",
  "devDependencies": {
    "pest": "^23.4.1"
  },
  "main": "index.js",
  "scripts": {
    "test": "pest"
  },
  "license": "MIT"
}'
```

Define a schema

```{r}
check_json <- Schema$new("CheckJSON",
  name = fields$character(required=TRUE),
  version = fields$character(required=TRUE),
  description = fields$character(required=TRUE),
  main = fields$character(),
  homepage = fields$url(),
  scripts = fields$named_list(keys=fields$character(), values=fields$character()),
  license = fields$character(required=TRUE),
  dependencies = fields$named_list(keys=fields$character(), values=fields$character()),
  dev_dependencies = fields$named_list(
    keys=fields$character(),
    values=fields$character(),
    data_key="devDependencies",
  ),
  unknown = "include"
)
```

Check the `json` against the schema. In this case it's valid, so we get the list output

```{r}
check_json$load_json(json)
```

In the above case, the JSON is valid. What if we give it bad JSON? Here, a required field (`license`) is missing:


```{r}
json_bad <- '{
  "name": "dunderscore",
  "version": "1.2.3",
  "description": "The Pythonic JavaScript toolkit",
  "devDependencies": {
    "pest": "^23.4.1"
  },
  "main": "index.js",
  "scripts": {
    "test": "pest"
  }
}'
```

Check it

```{r error=TRUE}
check_json$load_json(json_bad)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Any}
\alias{Any}
\title{Any}
\description{
A field that applies no formatting
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Any}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-serialize_}{\code{Any$serialize_()}}
\item \href{#method-deserialize_}{\code{Any$deserialize_()}}
\item \href{#method-clone}{\code{Any$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="initialize">}\href{../../staypuft/html/Field.html#method-initialize}{\code{staypuft::Field$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Any$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Any$deserialize_(value, attr = NULL, data = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Any$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_field_schema.R
\name{FieldABC}
\alias{FieldABC}
\title{FieldABC}
\description{
Abstract base class from which all Field classes inherit
}
\examples{
\dontrun{
x <- FieldABC$new()
x
# x$serialize()
# x$deserialize()
# x$serialize_()
# x$deserialize_()
}
}
\keyword{internal}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{parent}}{xxx}

\item{\code{name}}{xxx}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-serialize}{\code{FieldABC$serialize()}}
\item \href{#method-deserialize}{\code{FieldABC$deserialize()}}
\item \href{#method-serialize_}{\code{FieldABC$serialize_()}}
\item \href{#method-deserialize_}{\code{FieldABC$deserialize_()}}
\item \href{#method-clone}{\code{FieldABC$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize"></a>}}
\if{latex}{\out{\hypertarget{method-serialize}{}}}
\subsection{Method \code{serialize()}}{
serialize fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FieldABC$serialize()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize}{}}}
\subsection{Method \code{deserialize()}}{
deserialize fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FieldABC$deserialize()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
serialize_ fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FieldABC$serialize_()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
deserialize_ fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FieldABC$deserialize_()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{FieldABC$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Schema.R
\name{Schema}
\alias{Schema}
\title{Schema}
\description{
Base schema class with which to define custom schemas
}
\examples{
z <- Schema$new("FooBar",
  name = fields$character(),
  title = fields$character()
)
z
z$fields
names(z$fields)

x <- list(name = "Jane Doe", title = "Howdy doody")
x
z$dump(x)
z$dump_json(x)
z$dump_json(x, auto_unbox = TRUE)


z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character()
)
z
x <- list(name = "Jane Doe", title = "Howdy doody")
z$load(x)
z$load_json(jsonlite::toJSON(x, auto_unbox=TRUE))

# unknown field
# x <- list(name = "Jane Doe", my_title = "Howdy doody")
# z$load(x)
# z$load_json(jsonlite::toJSON(x, auto_unbox=TRUE))

# as data.frame
z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character()
)
x <- list(name = "Jane Doe", title = "hello world")
z$load(x, as_df = TRUE)

# list of lists
z <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character()
)
x <- list(
  list(name = "Jane Doe", title = "hello world"),
  list(name = "Alice Water", title = "bye mars")
)
z$load(x, many = TRUE)
# z$load(x, many = FALSE)

# data.frame's
x <- data.frame(name = "jill", title = "boss", stringsAsFactors = FALSE)
x2 <- data.frame(name = c("jill", "jane"), title = c("boss", "ceo"),
  stringsAsFactors = FALSE)
x2 <- data.frame(name = c("jill", "jane"), title = c("boss", "ceo"),
  stringsAsFactors = FALSE)
z <- Schema$new("FooBar",
  name = fields$character(),
  title = fields$character()
)
z$load_df(x)
z$load_df(x2)
z$load_df(x2, many = TRUE, simplifyVector = FALSE)

# nested
artist_schema <- Schema$new("ArtistSchema",
  name = fields$character(),
  role = fields$character(),
  instrument = fields$character()
)
album_schema <- Schema$new("AlbumSchema",
  title = fields$character(),
  release_date = fields$date(),
  artist = fields$nested(artist_schema)
)
artist_schema
album_schema
bowie <- list(name="David Bowie", role="lead", instrument="voice")
album <- list(title="Hunky Dory", release_date="12-17-1971", artist=bowie)
album_schema$dump(album)
album_schema$load(album)
## many
albums <- list(
  list(title="Hunky Dory", release_date="12-17-1971", artist=bowie),
  list(title="Mars and Venus", release_date="03-05-1974", artist=bowie)
)
album_schema$dump(albums, many=TRUE)
album_schema$load(albums, many=TRUE)
## bad
album$artist <- list(stuff = "things")
if (interactive()) album_schema$load(album)

# Deserialize/load and create object with post_load
z <- Schema$new("ArtistSchema",
  name = fields$character(),
  role = fields$character(),
  instrument = fields$character(),
  post_load = {
    function(x) structure(x, class = "Artist", attr = "hello")
  }
)
z$post_load
w <- list(name="David Bowie", role="lead", instrument="voice")
z$load(w)
print.Artist <- function(x) {
  cat("Artist\n")
  cat(sprintf("  name: \%s\n", x$name))
  cat(sprintf("  role: \%s\n", x$role))
  cat(sprintf("  instrument: \%s\n", x$instrument))
}
z$load(w)

# from json
json <- jsonlite::toJSON(w)
z$load_json(json)
## many
ww <- list(
  list(name="David Bowie", role="lead", instrument="voice"),
  list(name="Michael Jackson", role="lead", instrument="voice")
)
json <- jsonlite::toJSON(ww)
z$load_json(json, simplifyVector = FALSE, many = TRUE)
}
\section{Super class}{
\code{\link[staypuft:SchemaABC]{staypuft::SchemaABC}} -> \code{Schema}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{schema_name}}{the schema name}

\item{\code{fields}}{field names}

\item{\code{post_load}}{field names}

\item{\code{many}}{xxxx}

\item{\code{only}}{xxxx}

\item{\code{exclude}}{xxxx}

\item{\code{ordered}}{xxxx}

\item{\code{load_only}}{xxxx}

\item{\code{dump_only}}{xxxx}

\item{\code{partial}}{xxxx}

\item{\code{unknown}}{xxxx}

\item{\code{context}}{xxxx}

\item{\code{opts}}{field names}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Schema$new()}}
\item \href{#method-print}{\code{Schema$print()}}
\item \href{#method-dump}{\code{Schema$dump()}}
\item \href{#method-dump_json}{\code{Schema$dump_json()}}
\item \href{#method-load}{\code{Schema$load()}}
\item \href{#method-load_json}{\code{Schema$load_json()}}
\item \href{#method-load_df}{\code{Schema$load_df()}}
\item \href{#method-clone}{\code{Schema$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Schema} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$new(
  schema_name,
  ...,
  post_load = NULL,
  only = NULL,
  exclude = NULL,
  many = FALSE,
  context = NULL,
  load_only = NULL,
  dump_only = NULL,
  partial = FALSE,
  unknown = "raise"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{schema_name}}{(character) the schema name}

\item{\code{...}}{additional arguments, passed to \code{fields}}

\item{\code{only}}{Whitelist of the declared fields to select when
instantiating the Schema. If None, all fields are used. Nested fields
can be represented with dot delimiters.}

\item{\code{exclude}}{Blacklist of the declared fields to exclude
when instantiating the Schema. If a field appears in both \code{only} and
\code{exclude}, it is not used. Nested fields can be represented with dot
delimiters.}

\item{\code{many}}{Should be set to \code{True} if \code{obj} is a collection
so that the object will be serialized to a list.}

\item{\code{context}}{Optional context passed to :class:\code{fields.Method} and
:class:\code{fields.Function} fields.}

\item{\code{load_only}}{Fields to skip during serialization (write-only fields)}

\item{\code{dump_only}}{Fields to skip during deserialization (read-only fields)}

\item{\code{partial}}{Whether to ignore missing fields and not require
any fields declared. Propagates down to \code{Nested} fields as well. If
its value is an iterable, only missing fields listed in that iterable
will be ignored. Use dot delimiters to specify nested fields.}

\item{\code{unknown}}{Whether to exclude, include, or raise an error for unknown
fields in the data. Use \code{EXCLUDE}, \code{INCLUDE} or \code{RAISE}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for \code{Schema} objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dump"></a>}}
\if{latex}{\out{\hypertarget{method-dump}{}}}
\subsection{Method \code{dump()}}{
Convert various objects to a list
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$dump(x, many = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{input}

\item{\code{many}}{(logical) Should be set to \code{TRUE} if \code{obj} is a list of
lists. default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dump_json"></a>}}
\if{latex}{\out{\hypertarget{method-dump_json}{}}}
\subsection{Method \code{dump_json()}}{
Same as \code{dump()}, but returns JSON
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$dump_json(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{input}

\item{\code{...}}{additional params passed to \code{\link[jsonlite:toJSON]{jsonlite::toJSON()}}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
JSON (character)
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-load"></a>}}
\if{latex}{\out{\hypertarget{method-load}{}}}
\subsection{Method \code{load()}}{
Load data
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$load(data, many = FALSE, partial = FALSE, unknown = NULL, as_df = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{a named list}

\item{\code{many}}{(logical) Should be set to \code{TRUE} if \code{obj} is a list of
lists. default: \code{FALSE}}

\item{\code{partial}}{(logical) not implemented yet}

\item{\code{unknown}}{(character) one or "raise", "exclude", or "include".
default: "raise"}

\item{\code{as_df}}{(logical) convert to tibble? default: \code{FALSE}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
xxxx
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-load_json"></a>}}
\if{latex}{\out{\hypertarget{method-load_json}{}}}
\subsection{Method \code{load_json()}}{
Same as \code{load()}, but takes JSON as input
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$load_json(data, many = FALSE, partial = FALSE, unknown = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{a named list}

\item{\code{many}}{(logical) Should be set to \code{TRUE} if \code{obj} is a list of
lists. default: \code{FALSE}}

\item{\code{partial}}{(logical) not implemented yet}

\item{\code{unknown}}{(character) one or "raise", "exclude", or "include".
default: "raise"}

\item{\code{...}}{additional params passed to \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-load_df"></a>}}
\if{latex}{\out{\hypertarget{method-load_df}{}}}
\subsection{Method \code{load_df()}}{
Same as \code{load()}, but takes a data.frame as input
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$load_df(data, many = FALSE, partial = FALSE, unknown = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{a data.frame}

\item{\code{many}}{(logical) Should be set to \code{TRUE} if \code{obj} is a list of
lists. default: \code{FALSE}}

\item{\code{partial}}{(logical) not implemented yet}

\item{\code{unknown}}{(character) one or "raise", "exclude", or "include".
default: "raise"}

\item{\code{...}}{additional params passed to \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a list
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Schema$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_field_schema.R
\name{SchemaABC}
\alias{SchemaABC}
\title{SchemaABC}
\description{
Abstract base class from which all Schemas inherit
}
\examples{
\dontrun{
x <- SchemaABC$new()
x
# x$dump()
# x$dump_json()
# x$load()
# x$load_json()
}
}
\keyword{internal}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-dump}{\code{SchemaABC$dump()}}
\item \href{#method-dump_json}{\code{SchemaABC$dump_json()}}
\item \href{#method-load}{\code{SchemaABC$load()}}
\item \href{#method-load_json}{\code{SchemaABC$load_json()}}
\item \href{#method-clone}{\code{SchemaABC$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dump"></a>}}
\if{latex}{\out{\hypertarget{method-dump}{}}}
\subsection{Method \code{dump()}}{
dump fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SchemaABC$dump()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dump_json"></a>}}
\if{latex}{\out{\hypertarget{method-dump_json}{}}}
\subsection{Method \code{dump_json()}}{
dump_json fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SchemaABC$dump_json()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-load"></a>}}
\if{latex}{\out{\hypertarget{method-load}{}}}
\subsection{Method \code{load()}}{
load fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SchemaABC$load()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-load_json"></a>}}
\if{latex}{\out{\hypertarget{method-load_json}{}}}
\subsection{Method \code{load_json()}}{
load_json fun: to be replaced by child class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SchemaABC$load_json()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SchemaABC$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/staypuft-package.R
\docType{package}
\name{staypuft-package}
\alias{staypuft-package}
\alias{staypuft}
\title{staypuft}
\description{
Convert Complex Objects to and from R Data Structures
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Url}
\alias{Url}
\title{Url}
\description{
A validated URL field. Validation occurs during both
serialization and deserialization
}
\examples{
aschema <- Schema$new("aSchema",
  url = fields$url()
)
aschema
aschema$load(list(url = "https://ropensci.org/")) # good
if (interactive()) aschema$load(list(url = 6)) # bad

sch <- Schema$new("anotherschema",
  url = fields$url(schemes = c("https", "ftps"))
)
if (interactive()) sch$load(list(url = "http://google.com")) # bad
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Url}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Url$new()}}
\item \href{#method-validate_url}{\code{Url$validate_url()}}
\item \href{#method-clone}{\code{Url$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize_">}\href{../../staypuft/html/Field.html#method-deserialize_}{\code{staypuft::Field$deserialize_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize_">}\href{../../staypuft/html/Field.html#method-serialize_}{\code{staypuft::Field$serialize_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Url object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Url$new(..., relative = FALSE, schemes = NULL, require_tld = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{relative}}{Whether to allow relative URLs. NOT WORKING YET}

\item{\code{schemes}}{Valid schemes. By default \code{http}, \code{https}, \code{ftp}, and
\code{ftps} are allowed}

\item{\code{require_tld}}{Whether to reject non-FQDN hostnames. NOT WORKING YET}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validate_url"></a>}}
\if{latex}{\out{\hypertarget{method-validate_url}{}}}
\subsection{Method \code{validate_url()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Url$validate_url(relative = FALSE, schemes = NULL, require_tld = TRUE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Url$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Email}
\alias{Email}
\title{Email}
\description{
A validated email field. Validation occurs during both
serialization and deserialization
}
\examples{
z <- Schema$new("email", email = fields$email())
z
z$load(list(email = "blueberries@yahoo.com")) # good
if (interactive()) z$load(list(email = 'foobar')) # bad
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Email}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Email$new()}}
\item \href{#method-validate_email}{\code{Email$validate_email()}}
\item \href{#method-clone}{\code{Email$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize_">}\href{../../staypuft/html/Field.html#method-deserialize_}{\code{staypuft::Field$deserialize_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize_">}\href{../../staypuft/html/Field.html#method-serialize_}{\code{staypuft::Field$serialize_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Email object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Email$new(...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validate_email"></a>}}
\if{latex}{\out{\hypertarget{method-validate_email}{}}}
\subsection{Method \code{validate_email()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Email$validate_email()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Email$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Nested}
\alias{Nested}
\title{Nested}
\description{
Nest a Schema inside a field
}
\examples{
artist_schema <- Schema$new("ArtistSchema",
  name = fields$character()
)
x <- Nested$new(artist_schema)
x
x$nested
x$deserialize_(value = list(name = 6)) # good
if (interactive()) x$deserialize_(value = list(foobar = 6)) # bad
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Nested}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Nested$new()}}
\item \href{#method-schema}{\code{Nested$schema()}}
\item \href{#method-serialize_}{\code{Nested$serialize_()}}
\item \href{#method-test_list}{\code{Nested$test_list()}}
\item \href{#method-load_}{\code{Nested$load_()}}
\item \href{#method-deserialize_}{\code{Nested$deserialize_()}}
\item \href{#method-clone}{\code{Nested$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Nested object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Nested$new(
  nested,
  default = miss_ing,
  exclude = NULL,
  only = NULL,
  many = FALSE,
  unknown = "raise",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nested}}{The Schema class or class name (character)
to nest, or "self" to nest the \code{Schema} within itself}

\item{\code{exclude}}{A list or tuple of fields to exclude}

\item{\code{only}}{A list or tuple of fields to marshal. If \code{NULL}, all fields
are marshalled. This parameter takes precedence over \code{exclude}.}

\item{\code{many}}{Whether the field is a collection of objects.}

\item{\code{unknown}}{Whether to exclude, include, or raise an
error for unknown fields in the data. Use "raise", "exclude",
or "include"}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-schema"></a>}}
\if{latex}{\out{\hypertarget{method-schema}{}}}
\subsection{Method \code{schema()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Nested$schema()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Nested$serialize_(value, nested_obj, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-test_list"></a>}}
\if{latex}{\out{\hypertarget{method-test_list}{}}}
\subsection{Method \code{test_list()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Nested$test_list(self, value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-load_"></a>}}
\if{latex}{\out{\hypertarget{method-load_}{}}}
\subsection{Method \code{load_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Nested$load_(value, data = NULL, partial = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Nested$deserialize_(value, attr = NULL, data = NULL, partial = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Nested$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Number}
\alias{Number}
\title{Number}
\description{
A Number field
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Number}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Number$new()}}
\item \href{#method-format_num}{\code{Number$format_num()}}
\item \href{#method-validated}{\code{Number$validated()}}
\item \href{#method-to_string}{\code{Number$to_string()}}
\item \href{#method-serialize_}{\code{Number$serialize_()}}
\item \href{#method-deserialize_}{\code{Number$deserialize_()}}
\item \href{#method-clone}{\code{Number$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Integer object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Number$new(..., as_string = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{as_string}}{If \code{TRUE}, serialize to a string instead of
a \code{numeric} type}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-format_num"></a>}}
\if{latex}{\out{\hypertarget{method-format_num}{}}}
\subsection{Method \code{format_num()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Number$format_num(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validated"></a>}}
\if{latex}{\out{\hypertarget{method-validated}{}}}
\subsection{Method \code{validated()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Number$validated(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_string"></a>}}
\if{latex}{\out{\hypertarget{method-to_string}{}}}
\subsection{Method \code{to_string()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Number$to_string(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Number$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Number$deserialize_(value, attr = NULL, data = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Number$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{NamedList}
\alias{NamedList}
\title{NamedList}
\description{
A class for lists with key-value pairs - aka in R: named lists
}
\examples{
x <- Schema$new("foo",
  title = fields$character(),
  age = fields$integer(strict = TRUE),
  name = fields$named_list(keys=fields$character(), values=fields$number())
)
x
# good
x$load(list(
  title = "barry", 
  age = 3L, 
  name = list(foo = 3.4, iff = 5)))
# bad
if (interactive()) {
x$load(list(
  title = "barry", 
  age = 3L, 
  name = list(foo = 3.4, iff = "a")))

x$load(list(
  title = "barry", 
  age = 3L, 
  name = list("bar", "else")))
}
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{NamedList}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{NamedList$new()}}
\item \href{#method-serialize_}{\code{NamedList$serialize_()}}
\item \href{#method-deserialize_}{\code{NamedList$deserialize_()}}
\item \href{#method-clone}{\code{NamedList$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NamedList$new(keys, values, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NamedList$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NamedList$deserialize_(value, attr = NULL, data = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NamedList$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{XDate}
\alias{XDate}
\title{Date}
\description{
A formatted date string
}
\note{
e.g., value: 2014-12-22T03:12:58.019077+00:00
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Date}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{XDate$new()}}
\item \href{#method-format_date}{\code{XDate$format_date()}}
\item \href{#method-to_string}{\code{XDate$to_string()}}
\item \href{#method-serialize_}{\code{XDate$serialize_()}}
\item \href{#method-validated}{\code{XDate$validated()}}
\item \href{#method-deserialize_}{\code{XDate$deserialize_()}}
\item \href{#method-clone}{\code{XDate$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Date object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{XDate$new(format = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{@param}}{format Either "rfc" (for RFC822) or "iso" (for ISO8601)}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-format_date"></a>}}
\if{latex}{\out{\hypertarget{method-format_date}{}}}
\subsection{Method \code{format_date()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{XDate$format_date(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_string"></a>}}
\if{latex}{\out{\hypertarget{method-to_string}{}}}
\subsection{Method \code{to_string()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{XDate$to_string(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{XDate$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validated"></a>}}
\if{latex}{\out{\hypertarget{method-validated}{}}}
\subsection{Method \code{validated()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{XDate$validated(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{XDate$deserialize_(value, attr = NULL, data = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{XDate$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{UUID}
\alias{UUID}
\title{UUID}
\description{
A UUID field
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{\link[staypuft:Character]{staypuft::Character}} -> \code{UUID}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-validated}{\code{UUID$validated()}}
\item \href{#method-serialize_}{\code{UUID$serialize_()}}
\item \href{#method-deserialize_}{\code{UUID$deserialize_()}}
\item \href{#method-clone}{\code{UUID$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="initialize">}\href{../../staypuft/html/Field.html#method-initialize}{\code{staypuft::Field$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validated"></a>}}
\if{latex}{\out{\hypertarget{method-validated}{}}}
\subsection{Method \code{validated()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UUID$validated(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UUID$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UUID$deserialize_(value, attr = NULL, data = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{UUID$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Character}
\alias{Character}
\title{Character}
\description{
A string field
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Character}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-serialize_}{\code{Character$serialize_()}}
\item \href{#method-deserialize_}{\code{Character$deserialize_()}}
\item \href{#method-clone}{\code{Character$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="initialize">}\href{../../staypuft/html/Field.html#method-initialize}{\code{staypuft::Field$initialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Character$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Character$deserialize_(value, attr = NULL, data = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Character$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fields.R
\docType{data}
\name{fields}
\alias{fields}
\title{staypuft fields}
\format{
An object of class \code{list} of length 13.
}
\usage{
fields
}
\description{
staypuft fields
}
\section{Types of fields}{

\itemize{
\item \link{Missing}
\item \link{Character}
\item \link{UUID}
\item \link{Number}
\item \link{Integer}
\item \link{Boolean}
\item \link{Any}
\item \link{Date}
\item \link{Nested}
\item \link{Url}
\item \link{Email}
\item \link{NamedList}
\item more coming soon
}
}

\section{Usage}{

You can use any supported fields in the \code{fields} object.

You can call the field like \code{fields$character()},
or pass supported arguments like \code{fields$character(data_key = "foobar")}
}

\examples{
fields
fields$character()
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_registry.R
\name{ClassRegistry}
\alias{ClassRegistry}
\title{ClassRegistry}
\description{
Base schema class with which to define custom schemas
}
\examples{
x <- ClassRegistry$new()
x
x$registry
z <- Schema$new("FooBar",
  name = fields$character(),
  title = fields$character()
)
w <- Schema$new("MySchema",
  name = fields$character(),
  title = fields$character()
)
x
x$registry
x$register("FooBar", z)
x$register("MySchema", w)
x
x$registry

x$get_class("FooBar")
x$get_class("MySchema")
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{registry}}{(list) list of classes registered}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{ClassRegistry$print()}}
\item \href{#method-register}{\code{ClassRegistry$register()}}
\item \href{#method-get_class}{\code{ClassRegistry$get_class()}}
\item \href{#method-clone}{\code{ClassRegistry$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClassRegistry$print(x)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-register"></a>}}
\if{latex}{\out{\hypertarget{method-register}{}}}
\subsection{Method \code{register()}}{
Add a class to the registry of serializer classes
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClassRegistry$register(class_name, cls)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{class_name}}{(character) the class name}

\item{\code{cls}}{(Schema) the class object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
nothing, registers internally
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_class"></a>}}
\if{latex}{\out{\hypertarget{method-get_class}{}}}
\subsection{Method \code{get_class()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClassRegistry$get_class(class_name)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClassRegistry$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Missing}
\alias{Missing}
\title{Missing}
\description{
Missing class
}
\keyword{internal}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-bool}{\code{Missing$bool()}}
\item \href{#method-print}{\code{Missing$print()}}
\item \href{#method-clone}{\code{Missing$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-bool"></a>}}
\if{latex}{\out{\hypertarget{method-bool}{}}}
\subsection{Method \code{bool()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Missing$bool()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Missing$print(x, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Missing$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Boolean}
\alias{Boolean}
\title{Boolean}
\description{
A boolean field
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{Boolean}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Boolean$new()}}
\item \href{#method-serialize_}{\code{Boolean$serialize_()}}
\item \href{#method-deserialize_}{\code{Boolean$deserialize_()}}
\item \href{#method-clone}{\code{Boolean$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Boolean object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Boolean$new(..., truthy = NULL, falsy = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{truthy}}{Values that will (de)serialize to \code{TRUE}. If an
empty set, any non-falsy value will deserialize to \code{TRUE}. If \code{NULL},
xx will be used.}

\item{\code{falsy}}{Values that will (de)serialize to \code{FALSE}. If \code{NULL},
xx will be used.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Boolean$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Boolean$deserialize_(value, attr = NULL, data = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Boolean$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Field.R
\name{Field}
\alias{Field}
\title{Field}
\description{
Basic field from which other fields should extend
}
\details{
It applies no formatting by default, and should only be used
in cases where data does not need to be formatted before being
serialized or deserialized. On error, the name of the field will be
returned.
}
\examples{
x <- fields$field()
x
x$error_messages

z <- fields$character()
z
z$error_messages
z$serialize(attr = "foo", obj = list(foo = "bar"))
z$deserialize("foo")
z$deserialize(fields$missing())

x <- Schema$new("x", cow = fields$character(data_key = "stuff"))
x
x$fields$cow$data_key
if (interactive()) x$load(list(cow = 5))
x$load(list(stuff = 5))
}
\section{Super class}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{Field}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{class_name}}{(character) xxx}

\item{\code{CHECK_ATTRIBUTE}}{(logical) default: \code{TRUE}}

\item{\code{creation_index}}{(integer) xxx}

\item{\code{default}}{a class, default: \code{Missing}}

\item{\code{attribute}}{(character) xxx}

\item{\code{data_key}}{(character) xxx}

\item{\code{validate}}{xxx}

\item{\code{required}}{(logical) xxx}

\item{\code{allow_none}}{(logical) xxx}

\item{\code{load_only}}{(logical) xxx}

\item{\code{dump_only}}{(logical) xxx}

\item{\code{missing}}{(logical) xxx}

\item{\code{metadata}}{Extra arguments to be stored as metadata.}

\item{\code{error_messages}}{(list) xxx}

\item{\code{validators}}{(list) xxx}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Field$new()}}
\item \href{#method-print}{\code{Field$print()}}
\item \href{#method-get_value}{\code{Field$get_value()}}
\item \href{#method-validate_}{\code{Field$validate_()}}
\item \href{#method-fail}{\code{Field$fail()}}
\item \href{#method-validate_missing_}{\code{Field$validate_missing_()}}
\item \href{#method-serialize}{\code{Field$serialize()}}
\item \href{#method-deserialize}{\code{Field$deserialize()}}
\item \href{#method-bind_to_schema}{\code{Field$bind_to_schema()}}
\item \href{#method-serialize_}{\code{Field$serialize_()}}
\item \href{#method-deserialize_}{\code{Field$deserialize_()}}
\item \href{#method-context}{\code{Field$context()}}
\item \href{#method-root}{\code{Field$root()}}
\item \href{#method-clone}{\code{Field$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Field object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$new(
  default = miss_ing,
  attribute = NULL,
  data_key = NULL,
  validate = NULL,
  required = FALSE,
  allow_none = NULL,
  load_only = FALSE,
  dump_only = FALSE,
  missing = miss_ing,
  error_messages = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{default}}{If set, this value will be used during serialization if
the input value is missing. If not set, the field will be excluded from
the serialized output if the input value is missing. May be a value or
a callable.}

\item{\code{attribute}}{The name of the key to get the value from when
deserializing. If \code{None}, assumes the key has the same name as the
field.}

\item{\code{data_key}}{The name of the key to get the value from when
deserializing. If \code{None}, assumes the key has the same name as the field.}

\item{\code{validate}}{Validator or collection of validators that
are called during deserialization. Validator takes a field's input
value as its only parameter and returns a boolean. If it returns \code{FALSE},
an \code{ValidationError} is raised.}

\item{\code{required}}{Raise a \code{ValidationError} if the field value
is not supplied during deserialization.}

\item{\code{allow_none}}{Set this to \code{TRUE} if \code{None} should be considered a
valid value during validation/deserialization. If \code{missing=NULL}
and \code{allow_none} is unset, will default to \code{TRUE}. Otherwise, the
default is \code{FALSE}.}

\item{\code{load_only}}{If \code{TRUE} skip this field during serialization,
otherwise its value will be present in the serialized data.}

\item{\code{dump_only}}{If \code{TRUE} skip this field during deserialization,
otherwise its value will be present in the deserialized object. In the
context of an HTTP API, this effectively marks the field as "read-only".}

\item{\code{missing}}{Default deserialization value for the field if the field
is not found in the input data. May be a value or a callable.}

\item{\code{error_messages}}{Overrides for \code{Field.default_error_messages}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Field} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for Field objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$print(x, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_value"></a>}}
\if{latex}{\out{\hypertarget{method-get_value}{}}}
\subsection{Method \code{get_value()}}{
Return the value for a given key from an object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$get_value(obj, attr, accessor = NULL, default = miss_ing)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{obj}}{The object to get the value from}

\item{\code{attr}}{The attribute/key in \code{obj} to get the value from.}

\item{\code{accessor}}{(callback) A callable used to retrieve the value of \code{attr}}

\item{\code{default}}{If set, this value will be used during serialization if
the input value is missing. If not set, the field will be excluded from
the serialized output if the input value is missing. May be a value or
a callable.
from the object \code{obj}. Defaults to \code{marshmallow.utils.get_value}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validate_"></a>}}
\if{latex}{\out{\hypertarget{method-validate_}{}}}
\subsection{Method \code{validate_()}}{
Perform validation on \code{value}. Raise a \code{ValidationError}
if validation does not succeed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$validate_(value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{a value}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fail"></a>}}
\if{latex}{\out{\hypertarget{method-fail}{}}}
\subsection{Method \code{fail()}}{
A helper method that simply raises a
\code{ValidationError}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$fail(key)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{key}}{a key}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-validate_missing_"></a>}}
\if{latex}{\out{\hypertarget{method-validate_missing_}{}}}
\subsection{Method \code{validate_missing_()}}{
Validate missing values. Raise a \code{ValidationError}
if \code{value} should be considered missing.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$validate_missing_(value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{a value}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize"></a>}}
\if{latex}{\out{\hypertarget{method-serialize}{}}}
\subsection{Method \code{serialize()}}{
Pulls the value for the given key from the object,
applies the field's formatting and returns the result.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$serialize(attr, obj, accessor = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{attr}}{(character) The attribute or key to get from the object.}

\item{\code{obj}}{(character) The object to pull the key from.}

\item{\code{accessor}}{(callback) Function used to pull values from \code{obj}.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
raise ValidationError: In case of formatting problem
}

\subsection{Returns}{
xxxx
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize}{}}}
\subsection{Method \code{deserialize()}}{
Deserialize \code{value}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$deserialize(value, attr = NULL, data = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{The value to be deserialized.}

\item{\code{attr}}{(character) The attribute/key in \code{data} to be deserialized.}

\item{\code{data}}{(list) The raw input data passed to the \code{Schema.load}.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
raise ValidationError: If an invalid value is passed or if a
required value is missing.
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-bind_to_schema"></a>}}
\if{latex}{\out{\hypertarget{method-bind_to_schema}{}}}
\subsection{Method \code{bind_to_schema()}}{
Update field with values from its parent schema.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$bind_to_schema(field_name, schema)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{field_name}}{(character) Field name set in schema.}

\item{\code{schema}}{Parent schema.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-serialize_"></a>}}
\if{latex}{\out{\hypertarget{method-serialize_}{}}}
\subsection{Method \code{serialize_()}}{
Serializes \code{value} to a basic Python datatype. Noop by
default. Concrete :class:\code{Field} classes should implement this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$serialize_(value, attr = NULL, obj = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{The value to be deserialized.}

\item{\code{attr}}{(character) The attribute/key in \code{data} to be deserialized.}

\item{\code{obj}}{(character) The object to pull the key from.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
raise ValidationError: In case of formatting or validation
failure.
}

\subsection{Returns}{
The serialized value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-deserialize_"></a>}}
\if{latex}{\out{\hypertarget{method-deserialize_}{}}}
\subsection{Method \code{deserialize_()}}{
Deserialize value. Concrete :class:\code{Field} classes should implement this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$deserialize_(value, attr, data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{The value to be deserialized.}

\item{\code{attr}}{(character) The attribute/key in \code{data} to be deserialized.}

\item{\code{data}}{(list) The raw input data passed to the \code{Schema.load}.}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
raise ValidationError: In case of formatting or validation failure.
}

\subsection{Returns}{
The deserialized value
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-context"></a>}}
\if{latex}{\out{\hypertarget{method-context}{}}}
\subsection{Method \code{context()}}{
The context dictionary for the parent \code{Schema}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$context()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-root"></a>}}
\if{latex}{\out{\hypertarget{method-root}{}}}
\subsection{Method \code{root()}}{
Reference to the \code{Schema} that this field belongs
to even if it is buried in a \code{List}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$root()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{None} for unbound fields
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Field$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field-classes.R
\name{Integer}
\alias{Integer}
\title{Integer}
\description{
A Integer field
}
\keyword{internal}
\section{Super classes}{
\code{\link[staypuft:FieldABC]{staypuft::FieldABC}} -> \code{\link[staypuft:Field]{staypuft::Field}} -> \code{\link[staypuft:Number]{staypuft::Number}} -> \code{Integer}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Integer$new()}}
\item \href{#method-format_num}{\code{Integer$format_num()}}
\item \href{#method-clone}{\code{Integer$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="bind_to_schema">}\href{../../staypuft/html/Field.html#method-bind_to_schema}{\code{staypuft::Field$bind_to_schema()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="context">}\href{../../staypuft/html/Field.html#method-context}{\code{staypuft::Field$context()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="deserialize">}\href{../../staypuft/html/Field.html#method-deserialize}{\code{staypuft::Field$deserialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="fail">}\href{../../staypuft/html/Field.html#method-fail}{\code{staypuft::Field$fail()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="get_value">}\href{../../staypuft/html/Field.html#method-get_value}{\code{staypuft::Field$get_value()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="print">}\href{../../staypuft/html/Field.html#method-print}{\code{staypuft::Field$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="root">}\href{../../staypuft/html/Field.html#method-root}{\code{staypuft::Field$root()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="serialize">}\href{../../staypuft/html/Field.html#method-serialize}{\code{staypuft::Field$serialize()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_">}\href{../../staypuft/html/Field.html#method-validate_}{\code{staypuft::Field$validate_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Field" data-id="validate_missing_">}\href{../../staypuft/html/Field.html#method-validate_missing_}{\code{staypuft::Field$validate_missing_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Number" data-id="deserialize_">}\href{../../staypuft/html/Number.html#method-deserialize_}{\code{staypuft::Number$deserialize_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Number" data-id="serialize_">}\href{../../staypuft/html/Number.html#method-serialize_}{\code{staypuft::Number$serialize_()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Number" data-id="to_string">}\href{../../staypuft/html/Number.html#method-to_string}{\code{staypuft::Number$to_string()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="staypuft" data-topic="Number" data-id="validated">}\href{../../staypuft/html/Number.html#method-validated}{\code{staypuft::Number$validated()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Integer object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Integer$new(..., strict = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{strict}}{If \code{TRUE}, only integer types are valid. Otherwise,
any value castable to \code{integer} is valid}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-format_num"></a>}}
\if{latex}{\out{\hypertarget{method-format_num}{}}}
\subsection{Method \code{format_num()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Integer$format_num(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Integer$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
