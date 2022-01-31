handlr
======



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/handlr)](https://cranchecks.info/pkgs/handlr)
[![R-check](https://github.com/ropensci/handlr/workflows/R-check/badge.svg)](https://github.com/ropensci/handlr/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/handlr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/handlr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/handlr)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/handlr)](https://cran.r-project.org/package=handlr)


a tool for converting among citation formats.

heavily influenced by, and code ported from the Ruby gem `bolognese`

supported readers:

- citeproc
- ris
- bibtex
- codemeta
- cff

supported writers:

- citeproc
- ris
- bibtex
- schemaorg
- rdfxml
- codemeta
- cff

not supported yet, but plan to:

- crosscite

## Installation

stable version


```r
install.packages("handlr")
```

dev version


```r
remotes::install_github("ropensci/handlr")
```


```r
library("handlr")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/handlr/issues).
* License: MIT
* Get citation information for `handlr` in R doing `citation(package = 'handlr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
handlr 0.3.0
============

### DEPENDENCIES

* drop `RefManageR` from package Imports as it will likely be archived soon - add package `bibtex` to Suggests for reading/writing bibtex (can't be in Imports because it's Orphaned on CRAN) (#22)

### NEW FEATURES

* handlr gains support for Citation File Format (CFF), "plain text files with human- and machine-readable citation information for software". See https://citation-file-format.github.io/ for more info - new functions: `cff_reader()` and `cff_writer()` and associated changes in `HandlrClient`. Associated with CFF support, handlr gains new Import package `yaml`  (#16)

### MINOR IMPROVEMENTS

* improvements to Citeproc parsing: previously dropped many fields that we didn't support; now including all Citeproc fields that we don't specifically parse into extra fields prefixed with `csl_`  (#20)
* nothing changed, but see discussion of bibtex errors in case you run into them (#9)


handlr 0.2.0
============

### NEW FEATURES

* gains function `handl_to_df()`; converts any `handl` object (output from `HandlClient` or any `*_reader()` functions) to a data.frame for easier downstream data munging; `HandlClient` gains `$as_df()` method which runs `handl_to_df()`; to support this, now importing data.table package (#15) (#19) feature request by @GeraldCNelson

### MINOR IMPROVEMENTS

* now exporting the `print.handl` method. it only affects how a `handl` class object prints in the console, but is useful for making output more brief/concise (#14)
* filled out a lot more details of what a `handl` object contains. see `?handl` for the documentation (#17)


handlr 0.1.0
============

### NEW FEATURES

* Released to CRAN
## Test environments

* local OS X install, R 4.0.3 patched
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 1 reverse dependency. No problems were found. Summary at <https://github.com/ropensci/handlr/blob/master/revdep/README.md> 

---

This version includes dropping RefManageR package, scheduled to be archived soon, along with fixes and some new functions.

Thanks!
Scott Chamberlain
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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/handlr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/handlr.git`
* Make sure to track progress upstream (i.e., on our version of `handlr` at `ropensci/handlr`) by doing `git remote add upstream https://github.com/ropensci/handlr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/handlr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.3 (2020-10-10) |
|os       |macOS Catalina 10.15.7       |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |US/Pacific                   |
|date     |2020-10-14                   |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|handlr  |0.2.0 |0.3.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*handlr
======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/handlr)](https://cranchecks.info/pkgs/handlr)
[![R-check](https://github.com/ropensci/handlr/workflows/R-check/badge.svg)](https://github.com/ropensci/handlr/actions?query=workflow%3AR-check)
[![codecov.io](https://codecov.io/github/ropensci/handlr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/handlr?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/handlr)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/handlr)](https://cran.r-project.org/package=handlr)


a tool for converting among citation formats.

heavily influenced by, and code ported from the Ruby gem `bolognese`

supported readers:

- citeproc
- ris
- bibtex
- codemeta
- cff

supported writers:

- citeproc
- ris
- bibtex
- schemaorg
- rdfxml
- codemeta
- cff

not supported yet, but plan to:

- crosscite

## Installation

stable version

```{r eval=FALSE}
install.packages("handlr")
```

dev version

```{r eval=FALSE}
remotes::install_github("ropensci/handlr")
```

```{r eval=FALSE}
library("handlr")
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/handlr/issues).
* License: MIT
* Get citation information for `handlr` in R doing `citation(package = 'handlr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "handlr"
author: "Scott Chamberlain"
date: "2020-10-14"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to handlr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



**handlr**: convert among citation formats

heavily influenced by, and code ported from <https://github.com/datacite/bolognese>

supported readers:

- [citeproc][]
- [ris][]
- [bibtex][] (requires suggested package `bibtex`)
- [codemeta][]
- [cff][]

supported writers:

- [citeproc][]
- [ris][]
- [bibtex][]
- [schema.org][]
- [rdfxml][] (requires suggested package [jsonld][])
- [codemeta][]
- [cff][]

not supported yet, but plan to:

- crosscite


## Installation

stable version


```r
install.packages("handlr")
# OR
install.packages("handlr", repos = "https://dev.ropensci.org")
```

dev version


```r
remotes::install_github("ropensci/handlr")
```


```r
library("handlr")
```

## All in one

There's a single R6 interface to all readers and writers


```r
z <- system.file("extdata/citeproc.json", package = "handlr")
x <- HandlrClient$new(x = z)
x
#> <handlr> 
#>   doi: 
#>   ext: json
#>   format (guessed): citeproc
#>   path: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/handlr/extdata/citeproc.json
#>   string (abbrev.): none
```

read the file


```r
x$read(format = "citeproc")
```

the parsed content


```r
x$parsed
#> <handl> 
#>   from: citeproc
#>   many: FALSE
#>   count: 1
#>   first 10 
#>     id/doi: https://doi.org/10.5438/4k3m-nyvg
```

write out bibtex


```r
cat(x$write("bibtex"), sep = "\n")
#> @article{https://doi.org/10.5438/4k3m-nyvg,
#>   doi = {10.5438/4k3m-nyvg},
#>   author = {Martin Fenner},
#>   title = {Eating your own Dog Food},
#>   journal = {DataCite Blog},
#>   pages = {},
#>   publisher = {DataCite},
#>   year = {2016},
#> }
```

## Choose your own adventure

Instead of using the `HandlrClient`, you can use the regular functions for each 
reader or writer. They are:

- `citeproc_reader()` / `citeproc_writer()`
- `ris_reader()` / `ris_writer()`
- `bibtex_reader()` / `bibtex_writer()`
- `codemeta_reader()` / `codemeta_writer()`
- `schema_org_writer()`
- `rdf_xml_writer()`

## Convert data to data.frame


```r
z <- system.file('extdata/bib-many.bib', package = "handlr")
res2 <- bibtex_reader(x = z)
handl_to_df(res2)
#>             key                                        id    type bibtex_type
#> 1    Amano_2016 https://doi.org/10.1093%2fbiosci%2fbiw022 article     article
#> 2 Bachelot_2016       https://doi.org/10.1890%2f15-1397.1 article     article
#>     citeproc_type ris_type resource_type_general additional_type
#> 1 article-journal     JOUR                  <NA>  JournalArticle
#> 2 article-journal     JOUR                  <NA>  JournalArticle
#>                     doi                                   b_url
#> 1 10.1093/biosci/biw022 http://dx.doi.org/10.1093/biosci/biw022
#> 2     10.1890/15-1397.1     http://dx.doi.org/10.1890/15-1397.1
#>                                                                                                        title
#> 1                            Spatial Gaps in Global Biodiversity Information and the Role of Citizen Science
#> 2 Long-lasting effects of land use history on soil fungal communities in second-growth tropical rain forests
#>                                                                                           author
#> 1                                                            Amano T., Lamming J., Sutherland W.
#> 2 Bachelot B., Uriarte M., Zimmerman J., Thompson J., Leff J., Asiaii A., Koshner J., McGuire K.
#>                         publisher                         is_part_of
#> 1 Oxford University Press ({OUP})            Periodical;{BioScience}
#> 2                 Wiley-Blackwell Periodical;Ecological Applications
#>   date_published volume first_page last_page    state
#> 1           <NA>     66        393       400 findable
#> 2           <NA>   <NA>       <NA>      <NA> findable
```



[jsonld]: https://github.com/ropensci/jsonld/
[codemeta]: https://codemeta.github.io/
[citeproc]: https://en.wikipedia.org/wiki/CiteProc
[ris]: https://en.wikipedia.org/wiki/RIS_(file_format)
[bibtex]: http://www.bibtex.org/
[schema.org]: https://schema.org/
[rdfxml]: https://en.wikipedia.org/wiki/RDF/XML
[cff]: https://citation-file-format.github.io/
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handlr-package.r
\docType{package}
\name{handlr-package}
\alias{handlr-package}
\alias{handlr}
\title{\strong{Citation format converter}}
\description{
A tool for converting among citation formats
}
\section{supported readers}{

\itemize{
\item citeproc
\item ris
\item bibtex (requires suggested package \code{bibtex})
\item codemeta
\item cff
}
}

\section{supported writers}{

\itemize{
\item citeproc
\item ris
\item bibtex
\item schema.org
\item rdfxml (requires suggested package \code{jsonld})
\item codemeta
\item cff
}
}

\section{links for citation formats}{

\itemize{
\item citeproc: \url{https://en.wikipedia.org/wiki/CiteProc}
\item codemeta: \url{https://codemeta.github.io/}
\item ris: \url{https://en.wikipedia.org/wiki/RIS_(file_format)}
\item bibtex: \url{http://www.bibtex.org/}
\item schema.org: \url{https://schema.org/}
\item rdfxml: \url{https://en.wikipedia.org/wiki/RDF/XML}
\item cff: \url{https://citation-file-format.github.io/}
}
}

\author{
Scott Chamberlain \email{sckott@protonmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bibtex_writer.R
\name{bibtex_writer}
\alias{bibtex_writer}
\title{bibtex writer}
\usage{
bibtex_writer(z, key = NULL)
}
\arguments{
\item{z}{an object of class \code{handl}; see \link{handl} for more}

\item{key}{(character) optional bibtex key to use. if \code{NULL} we
attempt try the following fields in order: \code{key}, \code{identifier},
\code{id}, \code{doi}. if you pass in ouput from \code{\link[=bibtex_reader]{bibtex_reader()}} you're
likely to have a \code{key} field, but otherwise probably not}
}
\value{
an object of class \code{BibEntry}
}
\description{
bibtex writer
}
\examples{
(z <- system.file('extdata/citeproc.json', package = "handlr"))
(tmp <- citeproc_reader(z))
bibtex_writer(z = tmp)
cat(bibtex_writer(z = tmp), sep = "\n")

# give a bibtex key
cat(bibtex_writer(tmp, "foobar89"), sep = "\n")

# many at once
if (requireNamespace("bibtex", quietly=TRUE)) {
(z <- system.file('extdata/bib-many.bib', package = "handlr"))
out <- bibtex_reader(x = z)
bibtex_writer(out)
}
}
\seealso{
Other writers: 
\code{\link{cff_writer}()},
\code{\link{citeproc_writer}()},
\code{\link{codemeta_writer}()},
\code{\link{rdf_xml_writer}()},
\code{\link{ris_writer}()},
\code{\link{schema_org_writer}()}

Other bibtex: 
\code{\link{bibtex_reader}()}
}
\concept{bibtex}
\concept{writers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ris_writer.R
\name{ris_writer}
\alias{ris_writer}
\title{ris writer (Research Information Systems)}
\usage{
ris_writer(z)
}
\arguments{
\item{z}{an object of class \code{handl}; see \link{handl} for more}
}
\value{
text if one RIS citation or list of many
}
\description{
ris writer (Research Information Systems)
}
\examples{
# from a RIS file
z <- system.file('extdata/crossref.ris', package = "handlr")
tmp <- ris_reader(z)
cat(ris_writer(z = tmp))

# peerj
z <- system.file('extdata/peerj.ris', package = "handlr")
tmp <- ris_reader(z)
cat(ris_writer(z = tmp))

# plos
z <- system.file('extdata/plos.ris', package = "handlr")
tmp <- ris_reader(z)
cat(ris_writer(z = tmp))

# elsevier
z <- system.file('extdata/elsevier.ris', package = "handlr")
tmp <- ris_reader(z)
cat(ris_writer(z = tmp))

z <- system.file('extdata/citeproc.json', package = "handlr")
res <- citeproc_reader(z)
cat(ris_writer(z = res))

# many
## combine many RIS in a handl object
z <- system.file('extdata/crossref.ris', package = "handlr")
cr <- ris_reader(z)
z <- system.file('extdata/peerj.ris', package = "handlr")
prj <- ris_reader(z)
c(cr, prj)

# many bibtex to ris via c method
if (requireNamespace("bibtex", quietly=TRUE)) {
a <- system.file('extdata/bibtex.bib', package = "handlr")
b <- system.file('extdata/crossref.bib', package = "handlr")
aa <- bibtex_reader(a)
bb <- bibtex_reader(a)
(res <- c(aa, bb))
cat(ris_writer(res), sep = "\n\n")
}

## manhy Citeproc to RIS
z <- system.file('extdata/citeproc-many.json', package = "handlr")
w <- citeproc_reader(x = z)
ris_writer(w)
cat(ris_writer(w), sep = "\n")
}
\references{
RIS tags https://en.wikipedia.org/wiki/RIS_(file_format)
}
\seealso{
Other writers: 
\code{\link{bibtex_writer}()},
\code{\link{cff_writer}()},
\code{\link{citeproc_writer}()},
\code{\link{codemeta_writer}()},
\code{\link{rdf_xml_writer}()},
\code{\link{schema_org_writer}()}

Other ris: 
\code{\link{ris_reader}()}
}
\concept{ris}
\concept{writers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citeproc_reader.R
\name{citeproc_reader}
\alias{citeproc_reader}
\title{citeproc reader}
\usage{
citeproc_reader(x)
}
\arguments{
\item{x}{(character) a file path or string}
}
\value{
an object of class \code{handl}; see \link{handl} for more
}
\description{
citeproc reader
}
\examples{
# single
z <- system.file('extdata/citeproc.json', package = "handlr")
citeproc_reader(x = z)
w <- system.file('extdata/citeproc2.json', package = "handlr")
citeproc_reader(x = w)

# many
z <- system.file('extdata/citeproc-many.json', package = "handlr")
citeproc_reader(x = z)
}
\seealso{
Other readers: 
\code{\link{bibtex_reader}()},
\code{\link{cff_reader}()},
\code{\link{codemeta_reader}()},
\code{\link{ris_reader}()}

Other citeproc: 
\code{\link{citeproc_writer}()}
}
\concept{citeproc}
\concept{readers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handl.R
\name{handl}
\alias{handl}
\title{handl object}
\description{
handl object
}
\details{
A \code{handl} object is what's returned from the reader functions,
and what is passed to the writer functions. The \code{handl} object is a list,
but using the \code{print.handl} method makes it look something like:\preformatted{<handl>
  from: codemeta
  many: TRUE
  count: 2
  first 10
    id/doi: https://doi.org/10.5063\%2ff1m61h5x
    id/doi: https://doi.org/10.5063\%2ff1m61h5x
}

You can always \code{unclass()} the object to get the list itself.

The \code{handl} object follows \url{https://github.com/datacite/bolognese}, which
uses the Crosscite format as its internal representation. Note that we
don't currently support writing to or reading from Crosscite.

Details on each entry are stored in the named attributes:
\itemize{
\item from: the data type the citations come from
\item many: is there more than 1 citation?
\item count: number of citations
\item finally, some details of the first 10 are printed
}

If you have a \code{handl} object with 1 citation, it is a named list that
you can access with normal key indexing. If the result is length > 1,
the data is an unnamed list of named lists; the top level
list is unnamed, with each list within it being named.

Each named list should have the following components:
\itemize{
\item key: (string) a key for the citation, e.g., in a bibtex file
\item id: (string) an id for the work being referenced, often a DOI
\item type: (string) type of work
\item bibtex_type: (string) bibtex type
\item citeproc_type: (string) citeproc type
\item ris_type: (string) ris type
\item resource_type_general
\item additional_type: (string) additional type
\item doi: (string) DOI
\item b_url: (string) additional URL
\item title: (string) the title of the work
\item author: (list) authors, with each author a named list of
\itemize{
\item type: type, typically "Person"
\item name: full name
\item givenName: given (first) name
\item familyName: family (last) name
}
\item publisher: (string) the publisher name
\item is_part_of: (list) what the work is published in, or part of, a
named list with:
\itemize{
\item type: (string) the type of work
\item title: (string) title of the work, often a journal or edited book
\item issn: (string) the ISSN
}
\item date_published: (string)
\item volume: (string) the volume, if applicable
\item first_page: (string) the first page
\item last_page: (string) the last page
\item description: (string) description of the work, often an abstract
\item license: (string) license of the work, a named list
\item state: (string) the state of the list
\item software_version: (string) software version
}

Citeproc formats may have extra fields that begin with \code{csl_}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/schema_org_writer.R
\name{schema_org_writer}
\alias{schema_org_writer}
\title{Schema org writer}
\usage{
schema_org_writer(z, auto_unbox = TRUE, pretty = TRUE, ...)
}
\arguments{
\item{z}{an object of class \code{handl}; see \link{handl} for more}

\item{auto_unbox}{(logical) automatically "unbox" all atomic
vectors of length 1 (default: \code{TRUE}). passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}

\item{pretty}{(logical) adds indentation whitespace to JSON output
(default: \code{TRUE}), passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}

\item{...}{further params passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}
}
\value{
an object of class \code{json}
}
\description{
Schema org writer
}
\examples{
if (requireNamespace("bibtex", quietly=TRUE)) {
(z <- system.file('extdata/bibtex.bib', package = "handlr"))
(tmp <- bibtex_reader(z))
schema_org_writer(tmp)
schema_org_writer(tmp, pretty = FALSE)
}

# many citeproc to schema 
z <- system.file('extdata/citeproc-many.json', package = "handlr")
w <- citeproc_reader(x = z)
schema_org_writer(w)
schema_org_writer(w, pretty = FALSE)
}
\seealso{
Other writers: 
\code{\link{bibtex_writer}()},
\code{\link{cff_writer}()},
\code{\link{citeproc_writer}()},
\code{\link{codemeta_writer}()},
\code{\link{rdf_xml_writer}()},
\code{\link{ris_writer}()}
}
\concept{schema_org}
\concept{writers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codemeta_reader.R
\name{codemeta_reader}
\alias{codemeta_reader}
\title{codemeta reader}
\usage{
codemeta_reader(x)
}
\arguments{
\item{x}{(character) a file path or string (character or json)}
}
\value{
an object of class \code{handl}; see \link{handl} for more
}
\description{
codemeta reader
}
\examples{
# single
(z <- system.file('extdata/codemeta.json', package = "handlr"))
codemeta_reader(x = z)

# many
(z <- system.file('extdata/codemeta-many.json', package = "handlr"))
codemeta_reader(x = z)
}
\seealso{
Other readers: 
\code{\link{bibtex_reader}()},
\code{\link{cff_reader}()},
\code{\link{citeproc_reader}()},
\code{\link{ris_reader}()}

Other codemeta: 
\code{\link{codemeta_writer}()}
}
\concept{codemeta}
\concept{readers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c-handl.R
\name{c.handl}
\alias{c.handl}
\title{combine many handl objects}
\usage{
\method{c}{handl}(...)
}
\arguments{
\item{...}{one or more objects of class \code{handl}; see \link{handl} for more.
all inputs must be of class \code{handl}. if the first input is not of class
\code{handl}, you will not get back an object of class \code{handl}}
}
\value{
an object of class \code{handl} of length equal to number of
\code{handl} objects passed in
}
\description{
combine many handl objects
}
\examples{
z <- system.file('extdata/crossref.ris', package = "handlr")
cr <- ris_reader(z)
z <- system.file('extdata/peerj.ris', package = "handlr")
prj <- ris_reader(z)
res <- c(cr, prj)
res
invisible(lapply(bibtex_writer(res), cat, sep = "\n\n"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ris_reader.R
\name{ris_reader}
\alias{ris_reader}
\title{ris reader (Research Information Systems)}
\usage{
ris_reader(x)
}
\arguments{
\item{x}{(character) a file path or string}
}
\value{
an object of class \code{handl}; see \link{handl} for more
}
\description{
ris reader (Research Information Systems)
}
\examples{
z <- system.file('extdata/crossref.ris', package = "handlr")
ris_reader(z)

z <- system.file('extdata/peerj.ris', package = "handlr")
ris_reader(z)

z <- system.file('extdata/plos.ris', package = "handlr")
ris_reader(z)

# from a string
z <- system.file('extdata/crossref.ris', package = "handlr")
my_string <- ris_writer(ris_reader(z))
class(my_string)
ris_reader(my_string)

# many
z <- system.file('extdata/multiple-eg.ris', package = "handlr")
ris_reader(z)
}
\references{
RIS tags https://en.wikipedia.org/wiki/RIS_(file_format)
}
\seealso{
Other readers: 
\code{\link{bibtex_reader}()},
\code{\link{cff_reader}()},
\code{\link{citeproc_reader}()},
\code{\link{codemeta_reader}()}

Other ris: 
\code{\link{ris_writer}()}
}
\concept{readers}
\concept{ris}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codemeta_writer.R
\name{codemeta_writer}
\alias{codemeta_writer}
\title{codemeta writer}
\usage{
codemeta_writer(z, auto_unbox = TRUE, pretty = TRUE, ...)
}
\arguments{
\item{z}{an object of class \code{handl}; see \link{handl} for more}

\item{auto_unbox}{(logical) automatically "unbox" all atomic
vectors of length 1 (default: \code{TRUE}). passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}

\item{pretty}{(logical) adds indentation whitespace to JSON output
(default: \code{TRUE}), passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}

\item{...}{further params passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}
}
\value{
an object of class \code{json}
}
\description{
codemeta writer
}
\examples{
if (requireNamespace("bibtex", quietly=TRUE)) {
(x <- system.file('extdata/crossref.bib', package = "handlr"))
(z <- bibtex_reader(x))
codemeta_writer(z)
}

# many citeproc to schema 
z <- system.file('extdata/citeproc-many.json', package = "handlr")
w <- citeproc_reader(x = z)
codemeta_writer(w)
codemeta_writer(w, pretty = FALSE)
}
\seealso{
Other writers: 
\code{\link{bibtex_writer}()},
\code{\link{cff_writer}()},
\code{\link{citeproc_writer}()},
\code{\link{rdf_xml_writer}()},
\code{\link{ris_writer}()},
\code{\link{schema_org_writer}()}

Other codemeta: 
\code{\link{codemeta_reader}()}
}
\concept{codemeta}
\concept{writers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citeproc_writer.R
\name{citeproc_writer}
\alias{citeproc_writer}
\title{citeproc writer}
\usage{
citeproc_writer(z, auto_unbox = TRUE, pretty = TRUE, ...)
}
\arguments{
\item{z}{an object of class \code{handl}; see \link{handl} for more}

\item{auto_unbox}{(logical) automatically "unbox" all atomic
vectors of length 1 (default: \code{TRUE}). passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}

\item{pretty}{(logical) adds indentation whitespace to JSON output
(default: \code{TRUE}), passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}

\item{...}{further params passed to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}
}
\value{
citeproc as JSON
}
\description{
citeproc writer
}
\examples{
z <- system.file('extdata/citeproc.json', package = "handlr")
(tmp <- citeproc_reader(z))
citeproc_writer(z = tmp)
citeproc_writer(z = tmp, pretty = FALSE)
cat(ris_writer(z = tmp))

# many
z <- system.file('extdata/citeproc-many.json', package = "handlr")
w <- citeproc_reader(x = z)
citeproc_writer(w)
}
\seealso{
Other writers: 
\code{\link{bibtex_writer}()},
\code{\link{cff_writer}()},
\code{\link{codemeta_writer}()},
\code{\link{rdf_xml_writer}()},
\code{\link{ris_writer}()},
\code{\link{schema_org_writer}()}

Other citeproc: 
\code{\link{citeproc_reader}()}
}
\concept{citeproc}
\concept{writers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_writer.R
\name{cff_writer}
\alias{cff_writer}
\title{Citation File Format (cff) writer}
\usage{
cff_writer(
  z,
  path = NULL,
  message = "Please cite the following works when using this software."
)
}
\arguments{
\item{z}{an object of class \code{handl}; see \link{handl} for more}

\item{path}{a file path or connection; default: \code{stdout()}}

\item{message}{a message to display.
Defaults to \code{"Please cite the following works when using this software."}}
}
\value{
text if one cff citation or list of many
}
\description{
Citation File Format (cff) writer
}
\details{
uses \code{yaml::write_yaml} to write to yaml format that
CFF uses
}
\section{Converting to CFF from other formats}{

CFF has required fields that can't be missing. This means that
converting from other citation types to CFF will likely require
adding the required CFF fields manually. Adding fields to a
\code{handl} object is easy: it's really just an R list so add
named elements to it. The required CFF fields are:
\itemize{
\item CFF \strong{v1.1.0}:
\itemize{
\item cff-version: add \code{cff_version}
\item message: add \code{message}
\item version: add \code{software_version}
\item title: add \code{title}
\item authors: add \code{author}
\item date-released: add \code{date_published}
}
\item CFF \strong{v1.2.0}:
\itemize{
\item Only fields \code{cff-version}, \code{message}, \code{title} and \code{authors} are
required.
}
}

If \code{cff_version} is not provided, the value by default is "1.2.0".
}

\examples{
(z <- system.file('extdata/citation.cff', package = "handlr"))
res <- cff_reader(x = z)
res
unclass(res)
cff_writer(res)
cat(cff_writer(res))
f <- tempfile()
cff_writer(res, f)
readLines(f)
unlink(f)

# convert from a different citation format
## see "Converting to CFF from other formats" above
z <- system.file('extdata/citeproc.json', package = "handlr")
w <- citeproc_reader(x = z)
# cff_writer(w) # fails unless we add required fields
w$cff_version <- "1.1.0"
w$software_version <- "2.5"
w$title <- "A cool library"
w$date_published <- "2017-12-18"
cff_writer(w)
cat(cff_writer(w))
}
\references{
CFF format:
https://github.com/citation-file-format/citation-file-format
}
\seealso{
Other writers: 
\code{\link{bibtex_writer}()},
\code{\link{citeproc_writer}()},
\code{\link{codemeta_writer}()},
\code{\link{rdf_xml_writer}()},
\code{\link{ris_writer}()},
\code{\link{schema_org_writer}()}

Other cff: 
\code{\link{cff_reader}()}
}
\concept{cff}
\concept{writers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bibtex_reader.R
\name{bibtex_reader}
\alias{bibtex_reader}
\title{bibtex reader}
\usage{
bibtex_reader(x)
}
\arguments{
\item{x}{(character) a file path or a bibtex string}
}
\value{
an object of class \code{handl}; see \link{handl} for more
}
\description{
bibtex reader
}
\note{
requires package \code{bibtex}, an optional package for handlr
}
\examples{
if (requireNamespace("bibtex", quietly=TRUE)) {
(z <- system.file('extdata/crossref.bib', package = "handlr"))
bibtex_reader(x = z)
(z <- system.file('extdata/bibtex.bib', package = "handlr"))
bibtex_reader(x = z)

# many at once 
(z <- system.file('extdata/bib-many.bib', package = "handlr"))
bibtex_reader(x = z)
}
}
\seealso{
Other readers: 
\code{\link{cff_reader}()},
\code{\link{citeproc_reader}()},
\code{\link{codemeta_reader}()},
\code{\link{ris_reader}()}

Other bibtex: 
\code{\link{bibtex_writer}()}
}
\concept{bibtex}
\concept{readers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handl_to_df.R
\name{handl_to_df}
\alias{handl_to_df}
\title{handl to data.frame conversion}
\usage{
handl_to_df(x)
}
\arguments{
\item{x}{an object of class handl}
}
\value{
data.frame with column following \link{handl}, with as many rows
as there are citations
}
\description{
handl to data.frame conversion
}
\note{
requires the Suggested package \code{data.table}
}
\examples{
z <- system.file('extdata/crossref.ris', package = "handlr")
res <- ris_reader(z)
handl_to_df(res)

(x <- HandlrClient$new(x = z))
x$as_df() # empty data.frame
x$read()
x$as_df() # data.frame with citation data

if (requireNamespace("bibtex", quietly=TRUE)) {
(z <- system.file('extdata/bib-many.bib', package = "handlr"))
res2 <- bibtex_reader(x = z)
handl_to_df(res2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a_types.R
\docType{data}
\name{cff_reference_types}
\alias{cff_reference_types}
\title{cff references types}
\format{
An object of class \code{character} of length 47.
}
\usage{
cff_reference_types
}
\description{
cff references types
}
\details{
cff citation format types for references
}
\references{
http://bit.ly/2PRK1Vt
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_reader.R
\name{cff_reader}
\alias{cff_reader}
\title{Citation File Format (cff) reader}
\usage{
cff_reader(x)
}
\arguments{
\item{x}{(character) a file path or a yaml string}
}
\value{
an object of class \code{handl}; see \link{handl} for more
}
\description{
Citation File Format (cff) reader
}
\details{
CFF only supports one citation, so \code{many} will always be
\code{FALSE}.

Required fields:
\itemize{
\item CFF \strong{v1.1.0}: \code{cff-version}, \code{version}, \code{message}, \code{date-released},
\code{title}, \code{authors}.
\item CFF \strong{v1.2.0}: \code{cff-version}, \code{message}, \code{title}, \code{authors}.
}

We'll stop with error if any of these are missing.

You can though have many references in your CFF file
associated with the citation. \code{references} is an optional component in
cff files. If included, we check the following:
\itemize{
\item each reference must have the 3 required fields: type, authors, title
\item type must be in the allowed set, see \link{cff_reference_types}
\item the elements within authors must each be an entity or person object
https://github.com/citation-file-format/citation-file-format#entity-objects
https://github.com/citation-file-format/citation-file-format#person-objects
\item title must be a string
}
}
\examples{
(z <- system.file("extdata/citation.cff", package = "handlr"))
res <- cff_reader(x = z)
res
res$cff_version
res$software_version
res$message
res$id
res$doi
res$title
res$author
res$references

# no references
(z <- system.file("extdata/citation-norefs.cff", package = "handlr"))
out <- cff_reader(x = z)
out
out$references
}
\references{
CFF format:
https://github.com/citation-file-format/citation-file-format
}
\seealso{
Other readers: 
\code{\link{bibtex_reader}()},
\code{\link{citeproc_reader}()},
\code{\link{codemeta_reader}()},
\code{\link{ris_reader}()}

Other cff: 
\code{\link{cff_writer}()}
}
\concept{cff}
\concept{readers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client.R
\name{HandlrClient}
\alias{HandlrClient}
\title{HandlrClient}
\description{
handlr client, read and write to and from all citation formats
}
\details{
The various inputs to the \code{x} parameter are handled in different
ways:
\itemize{
\item file: contents read from file, we grab file extension, and we guess
format based on combination of contents and file extension because
file extensions may belie what's in the file
\item string: string read in, and we guess format based on contents of
the string
\item DOI: we request citeproc-json format from the Crossref API
\item DOI url: we request citeproc-json format from the Crossref API
}
}
\note{
If \verb{$parsed} is \code{NULL} then it's likely \verb{$read()} has not
been run - in which case we attempt to run \verb{$read()} to
populate \verb{$parsed}
}
\examples{
# read() can be run with format specified or not
# if format not given, we attempt to guess the format and then read
z <- system.file('extdata/citeproc.json', package = "handlr")
(x <- HandlrClient$new(x = z))
x$read()
x$read("citeproc")
x$parsed

# you can run read() then write()
# or just run write(), and read() will be run for you if possible
z <- system.file('extdata/citeproc.json', package = "handlr")
(x <- HandlrClient$new(x = z))
cat(x$write("ris"))

# read from a DOI as a url
if (interactive()) {
  (x <- HandlrClient$new('https://doi.org/10.7554/elife.01567'))
  x$parsed
  x$read()
  x$parsed
  x$write('bibtex')
}

# read from a DOI
if (interactive()) {
  (x <- HandlrClient$new('10.7554/elife.01567'))
  x$parsed
  x$read()
  x$write('bibtex')
}

# read in citeproc, write out bibtex
z <- system.file('extdata/citeproc.json', package = "handlr")
(x <- HandlrClient$new(x = z))
x$path
x$ext
x$read("citeproc")
x$parsed
x$write("bibtex")
f <- tempfile(fileext = ".bib")
x$write("bibtex", file = f)
readLines(f)
unlink(f)

# read in ris, write out ris
z <- system.file('extdata/peerj.ris', package = "handlr")
(x <- HandlrClient$new(x = z))
x$path
x$format_guessed
x$read("ris")
x$parsed
x$write("ris")
cat(x$write("ris"))

# read in bibtex, write out ris
(z <- system.file('extdata/bibtex.bib', package = "handlr"))
(x <- HandlrClient$new(x = z))
x$path
x$format_guessed
if (requireNamespace("bibtex", quietly = TRUE)) {
x$read("bibtex")
x$parsed
x$write("ris")
cat(x$write("ris"))
}

# read in bibtex, write out RDF XML
if (requireNamespace("bibtex", quietly = TRUE) && interactive()) {
  (z <- system.file('extdata/bibtex.bib', package = "handlr"))
  (x <- HandlrClient$new(x = z))
  x$path
  x$format_guessed
  x$read("bibtex")
  x$parsed
  x$write("rdfxml")
  cat(x$write("rdfxml"))
}

# codemeta
(z <- system.file('extdata/codemeta.json', package = "handlr"))
(x <- HandlrClient$new(x = z))
x$path
x$format_guessed
x$read("codemeta")
x$parsed
x$write("codemeta")

# cff: Citation File Format
(z <- system.file('extdata/citation.cff', package = "handlr"))
(x <- HandlrClient$new(x = z))
x$path
x$format_guessed
x$read("cff")
x$parsed
x$write("codemeta")

# > 1 citation
z <- system.file('extdata/citeproc-many.json', package = "handlr")
(x <- HandlrClient$new(x = z))
x$parsed
x$read()
x$parsed
## schmea org
x$write("schema_org")
## bibtex
x$write("bibtex")
## bibtex to file
f <- tempfile(fileext=".bib")
x$write("bibtex", f)
readLines(f)
unlink(f)
## to RIS
x$write("ris")
### only one per file, so not combined
files <- replicate(2, tempfile(fileext=".ris"))
x$write("ris", files)
lapply(files, readLines)

# handle strings instead of files
z <- system.file('extdata/citeproc-crossref.json', package = "handlr")
(x <- HandlrClient$new(x = readLines(z)))
x$read("citeproc")
x$parsed
cat(x$write("bibtex"), sep = "\n")
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{path}}{(character) non-empty if file path passed to initialize}

\item{\code{string}}{(character) non-empty if string (non-file) passed to initialize}

\item{\code{parsed}}{after \code{read()} is run, the parsed content}

\item{\code{file}}{(logical) \code{TRUE} if a file passed to initialize, else \code{FALSE}}

\item{\code{ext}}{(character) the file extension}

\item{\code{format_guessed}}{(character) the guessed file format}

\item{\code{doi}}{(character) the DOI, if any found}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-print}{\code{HandlrClient$print()}}
\item \href{#method-new}{\code{HandlrClient$new()}}
\item \href{#method-read}{\code{HandlrClient$read()}}
\item \href{#method-write}{\code{HandlrClient$write()}}
\item \href{#method-as_df}{\code{HandlrClient$as_df()}}
\item \href{#method-clone}{\code{HandlrClient$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for \code{HandlrClient} objects
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HandlrClient$print(x, ...)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{HandlrClient} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HandlrClient$new(x, format = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{(character) a file path (the file must exist), a string
containing contents of the citation, a DOI, or a DOI as a URL.
See Details.}

\item{\code{format}}{(character) one of citeproc, ris, bibtex, codemeta, cff,
or \code{NULL}. If \code{NULL}, we attempt to guess the format, and error if we
can not guess}

\item{\code{...}}{curl options passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{HandlrClient} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-read"></a>}}
\if{latex}{\out{\hypertarget{method-read}{}}}
\subsection{Method \code{read()}}{
read input
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HandlrClient$read(format = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{format}}{(character) one of citeproc, ris, bibtex, codemeta, cff,
or \code{NULL}. If \code{NULL}, we attempt to guess the format, and error if we
can not guess}

\item{\code{...}}{further args to the writer fxn, if any}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-write"></a>}}
\if{latex}{\out{\hypertarget{method-write}{}}}
\subsection{Method \code{write()}}{
write to std out or file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HandlrClient$write(format, file = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{format}}{(character) one of citeproc, ris, bibtex, schema_org,
rdfxml, codemeta, or cff}

\item{\code{file}}{a file path, if NULL to stdout. for \code{format=ris},
number of files must equal number of ris citations}

\item{\code{...}}{further args to the writer fxn, if any}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-as_df"></a>}}
\if{latex}{\out{\hypertarget{method-as_df}{}}}
\subsection{Method \code{as_df()}}{
convert data to a data.frame using \code{\link[=handl_to_df]{handl_to_df()}}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HandlrClient$as_df()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a data.frame
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HandlrClient$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/rdf_xml_writer.R
\name{rdf_xml_writer}
\alias{rdf_xml_writer}
\title{RDF XML writer}
\usage{
rdf_xml_writer(z, ...)
}
\arguments{
\item{z}{an object of class \code{handl}; see \link{handl} for more}

\item{...}{further params passed to \code{\link[jsonld:jsonld]{jsonld::jsonld_to_rdf()}}}
}
\value{
RDF XML
}
\description{
RDF XML writer
}
\details{
package \code{jsonld} required for this writer
}
\examples{
if (require("jsonld") && interactive()) {
  library("jsonld")
  z <- system.file('extdata/citeproc.json', package = "handlr")
  (tmp <- citeproc_reader(z))
 
  if (requireNamespace("bibtex", quietly=TRUE)) {
  (z <- system.file('extdata/bibtex.bib', package = "handlr"))
  (tmp <- bibtex_reader(z))
  rdf_xml_writer(z = tmp)
  cat(rdf_xml_writer(z = tmp))
  }
}
}
\seealso{
Other writers: 
\code{\link{bibtex_writer}()},
\code{\link{cff_writer}()},
\code{\link{citeproc_writer}()},
\code{\link{codemeta_writer}()},
\code{\link{ris_writer}()},
\code{\link{schema_org_writer}()}
}
\concept{rdf-xml}
\concept{writers}
