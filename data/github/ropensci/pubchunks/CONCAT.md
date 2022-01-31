

pubchunks
=========

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/pubchunks)](https://cranchecks.info/pkgs/pubchunks)
[![R-check](https://github.com/ropensci/pubchunks/workflows/R-check/badge.svg)](https://github.com/ropensci/pubchunks/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/pubchunks/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/pubchunks)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/pubchunks)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/pubchunks)](https://cran.r-project.org/package=pubchunks)

## Get chunks of XML articles


## Package API

 - pub_tabularize
 - pub_guess_publisher
 - pub_sections
 - pub_chunks
 - pub_providers

The main workhorse function is `pub_chunks()`. It allows you to pull out sections of articles from many different publishers (see next section below) WITHOUT having to know how to parse/navigate XML. XML has a steep learning curve, and can require quite a bit of Googling to sort out how to get to different parts of an XML document. 

The other main function is `pub_tabularize()` - which takes the output of `pub_chunks()` and coerces into a data.frame for easier downstream processing.

## Supported publishers/sources

- eLife
- PLOS
- Entrez/Pubmed
- Elsevier
- Hindawi
- Pensoft
- PeerJ
- Copernicus
- Frontiers
- F1000 Research

If you know of other publishers or sources that provide XML let us know by [opening an issue](https://github.com/ropensci/pubchunks/issues).

We'll continue adding additional publishers.


## Installation

Stable version


```r
install.packages("pubchunks")
```

Development version from GitHub


```r
remotes::install_github("ropensci/pubchunks")
```

Load library


```r
library('pubchunks')
```

## Working with files


```r
x <- system.file("examples/10_1016_0021_8928_59_90156_x.xml", 
  package = "pubchunks")
```


```r
pub_chunks(x, "abstract")
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: abstract
#>   showing up to first 5: 
#>    abstract (n=1): Abstract
#>                
#>                   This pa ...
pub_chunks(x, "title")
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
pub_chunks(x, "authors")
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: authors
#>   showing up to first 5: 
#>    authors (n=1): Chetaev, D.N
pub_chunks(x, c("title", "refs"))
#> <pub chunks>
#>   from: file
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title, refs
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
#>    refs (n=6): Watson G.N.. 1949. Teoriia besselevykh funktsii. N
```

The output of `pub_chunks()` is a list with an S3 class `pub_chunks` to make 
internal work in the package easier. You can easily see the list structure 
by using `unclass()`.

## Working with the xml already in a string


```r
xml <- paste0(readLines(x), collapse = "")
pub_chunks(xml, "title")
#> <pub chunks>
#>   from: character
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
```

## Working with xml2 class object


```r
xml <- paste0(readLines(x), collapse = "")
xml <- xml2::read_xml(xml)
pub_chunks(xml, "title")
#> <pub chunks>
#>   from: xml_document
#>   publisher/journal: elsevier/Journal of Applied Mathematics and Mechanics
#>   sections: title
#>   showing up to first 5: 
#>    title (n=1): On the driving of a piston with a rigid collar int ...
```

## Working with output of fulltext::ft_get()


```r
install.packages("fulltext")
```


```r
library("fulltext")
x <- fulltext::ft_get('10.1371/journal.pone.0086169')
pub_chunks(fulltext::ft_collect(x), sections="authors")
#> $plos
#> $plos$`10.1371/journal.pone.0086169`
#> <pub chunks>
#>   from: xml_document
#>   publisher/journal: plos/PLoS ONE
#>   sections: authors
#>   showing up to first 5: 
#>    authors (n=4): nested list
#> 
#> 
#> attr(,"ft_data")
#> [1] TRUE
```

## Coerce pub_chunks output into data.frame's


```r
x <- system.file("examples/elife_1.xml", package = "pubchunks")
res <- pub_chunks(x, c("doi", "title", "keywords"))
pub_tabularize(res)
#>                   doi                                          title
#> 1 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 2 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 3 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 4 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 5 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#> 6 10.7554/eLife.03032 MicroRNA-mediated repression of nonsense mRNAs
#>                       keywords .publisher
#> 1                     microRNA      elife
#> 2            nonsense mutation      elife
#> 3 nonsense-mediated mRNA decay      elife
#> 4                          APC      elife
#> 5             intron retention      elife
#> 6  premature termination codon      elife
```

## Get a random XML article


```r
library(rcrossref)
library(dplyr)

res <- cr_works(filter = list(
    full_text_type = "application/xml", 
    license_url="http://creativecommons.org/licenses/by/4.0/"))
links <- bind_rows(res$data$link) %>% filter(content.type == "application/xml")
download.file(links$URL[1], (i <- tempfile(fileext = ".xml")))
pub_chunks(i)
#> <pub chunks>
#>   from: file
#>   publisher/journal: unknown/NA
#>   sections: all
#>   showing up to first 5: 
#>    front (n=0): 
#>    body (n=0): 
#>    back (n=0): 
#>    title (n=0): 
#>    doi (n=0):
download.file(links$URL[13], (j <- tempfile(fileext = ".xml")))
pub_chunks(j)
#> <pub chunks>
#>   from: file
#>   publisher/journal: hindawi/BioMed Research International
#>   sections: all
#>   showing up to first 5: 
#>    front (n=2): nested list
#>    body (n=49): Oxidative stress and Reactive Oxygen Species (ROS)
#>    back (n=4): nested list
#>    title (n=1): Selected Enzyme Inhibitory Effects of Euphorbia ch ...
#>    doi (n=1): 10.1155/2018/1219367
download.file(links$URL[20], (k <- tempfile(fileext = ".xml")))
pub_chunks(k)
#> <pub chunks>
#>   from: file
#>   publisher/journal: hindawi/Case Reports in Pathology
#>   sections: all
#>   showing up to first 5: 
#>    front (n=2): nested list
#>    body (n=16): Bonnetti et al. first noted in 1992 an unusual cel
#>    back (n=3): nested list
#>    title (n=1): An Inguinal Perivascular Epithelioid Cell Tumor Me ...
#>    doi (n=1): 10.1155/2018/5749421
```




## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pubchunks/issues).
* License: MIT
* Get citation information for `pubchunks`: `citation(package = 'pubchunks')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
pubchunks 0.4.0
===============

### NEW FEATURES

* `pub_chunks()` gains new parameter `extract` to determine the final step of extracting text from XML format articles - use either `as.character()` or `xml2::xml_text()` (#11)

### MINOR IMPROVEMENTS

* fix to `sections="aff"` in `pub_chunks()` - keep all rows in merging data from different affilitation metadata parts to not lose any data (#7)


pubchunks 0.3.0
===============

### MINOR IMPROVEMENTS

* improvements to `pub_chunks()` with `sections="refs"` (fetching references from an article). all reference extraction was simply extracting the entire reference, which lost all white space. now there are custom reference parsers for the various types of reference formats. there still are outstanding issues, so do please point out where ref. extraction is not correct (#2)
* fix pub_guess_publisher examples error (#10)

pubchunks 0.2.2
===============

### MINOR IMPROVEMENTS

* fixed failing example (#9)

pubchunks 0.2.0
===============

### MINOR IMPROVEMENTS

* most section options in `pub_chunks()` now have defaults for extracting the section, and return NULL/empty list when not found (#3) (#4)
* improvements to `print.pub_chunks` so that the printed object contains more information (publisher/journal title) and more accurate ('character' used to include xml as character string and file paths, but are separated now). in addition, we state that the first 5 sections are printed so the user knows there could be more (#8)
* fix `pub_tabularize()` to accept list outputs from `pub_chunks()` (#5)

pubchunks 0.1.0
===============

### NEW FEATURES

* released to CRAN
## Test environments

* local OS X install, R 4.0.3
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There were no problems with the 1 reverse dependency.

--------

This version introduced a new parameter in a function and fixes an XML parsing issue.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/pubchunks/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/pubchunks.git`
* Make sure to track progress upstream (i.e., on our version of `pubchunks` at `ropensci/pubchunks`) by doing `git remote add upstream https://github.com/ropensci/pubchunks.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/pubchunks`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  cache.path = "inst/cache/"
)
```

pubchunks
=========

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/pubchunks)](https://cranchecks.info/pkgs/pubchunks)
[![R-check](https://github.com/ropensci/pubchunks/workflows/R-check/badge.svg)](https://github.com/ropensci/pubchunks/actions?query=workflow%3AR-check)
[![codecov](https://codecov.io/gh/ropensci/pubchunks/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/pubchunks)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/pubchunks)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/pubchunks)](https://cran.r-project.org/package=pubchunks)

## Get chunks of XML articles


## Package API

```{r echo=FALSE, comment=NA, results='asis'}
cat(paste(" -", paste(getNamespaceExports("pubchunks"), collapse = "\n - ")))
```

The main workhorse function is `pub_chunks()`. It allows you to pull out sections of articles from many different publishers (see next section below) WITHOUT having to know how to parse/navigate XML. XML has a steep learning curve, and can require quite a bit of Googling to sort out how to get to different parts of an XML document. 

The other main function is `pub_tabularize()` - which takes the output of `pub_chunks()` and coerces into a data.frame for easier downstream processing.

## Supported publishers/sources

- eLife
- PLOS
- Entrez/Pubmed
- Elsevier
- Hindawi
- Pensoft
- PeerJ
- Copernicus
- Frontiers
- F1000 Research

If you know of other publishers or sources that provide XML let us know by [opening an issue](https://github.com/ropensci/pubchunks/issues).

We'll continue adding additional publishers.


## Installation

Stable version

```{r eval=FALSE}
install.packages("pubchunks")
```

Development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/pubchunks")
```

Load library

```{r}
library('pubchunks')
```

## Working with files

```{r}
x <- system.file("examples/10_1016_0021_8928_59_90156_x.xml", 
  package = "pubchunks")
```

```{r}
pub_chunks(x, "abstract")
pub_chunks(x, "title")
pub_chunks(x, "authors")
pub_chunks(x, c("title", "refs"))
```

The output of `pub_chunks()` is a list with an S3 class `pub_chunks` to make 
internal work in the package easier. You can easily see the list structure 
by using `unclass()`.

## Working with the xml already in a string

```{r}
xml <- paste0(readLines(x), collapse = "")
pub_chunks(xml, "title")
```

## Working with xml2 class object

```{r}
xml <- paste0(readLines(x), collapse = "")
xml <- xml2::read_xml(xml)
pub_chunks(xml, "title")
```

## Working with output of fulltext::ft_get()

```{r eval=FALSE}
install.packages("fulltext")
```

```{r}
library("fulltext")
x <- fulltext::ft_get('10.1371/journal.pone.0086169')
pub_chunks(fulltext::ft_collect(x), sections="authors")
```

## Coerce pub_chunks output into data.frame's

```{r}
x <- system.file("examples/elife_1.xml", package = "pubchunks")
res <- pub_chunks(x, c("doi", "title", "keywords"))
pub_tabularize(res)
```

## Get a random XML article

```{r cache=TRUE}
library(rcrossref)
library(dplyr)

res <- cr_works(filter = list(
    full_text_type = "application/xml", 
    license_url="http://creativecommons.org/licenses/by/4.0/"))
links <- bind_rows(res$data$link) %>% filter(content.type == "application/xml")
download.file(links$URL[1], (i <- tempfile(fileext = ".xml")))
pub_chunks(i)
download.file(links$URL[13], (j <- tempfile(fileext = ".xml")))
pub_chunks(j)
download.file(links$URL[20], (k <- tempfile(fileext = ".xml")))
pub_chunks(k)
```

```{r echo=FALSE}
unlink(i)
unlink(j)
unlink(k)
```


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/pubchunks/issues).
* License: MIT
* Get citation information for `pubchunks`: `citation(package = 'pubchunks')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pub_providers.R
\name{pub_providers}
\alias{pub_providers}
\title{Providers}
\usage{
pub_providers()
}
\value{
character vector
}
\description{
The possible providers to select from
}
\examples{
pub_providers()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabularize.R
\name{pub_tabularize}
\alias{pub_tabularize}
\title{Tabularize chunks output}
\usage{
pub_tabularize(x, bind = FALSE)
}
\arguments{
\item{x}{the output of \code{\link[=pub_chunks]{pub_chunks()}}}

\item{bind}{(logical) whether to bind list of data.frames or not.
ignored unless \code{list} input to \code{x}. default: \code{FALSE}}
}
\value{
a data.frame or list
}
\description{
Tabularize chunks output
}
\examples{
\dontrun{
# one at a time
## example 1, a file path
x <- system.file("examples/elife_1.xml", package = "pubchunks")
(res <- pub_chunks(x, c("doi", "title", "keywords")))
pub_tabularize(res)

## example 2, a file path
y <- system.file("examples/frontiers_1.xml", package = "pubchunks")
(res2 <- pub_chunks(y, c("doi", "title", "keywords")))
pub_tabularize(res2)

# > 1, a list of file paths
x <- system.file("examples/elife_1.xml", package = "pubchunks")
y <- system.file("examples/frontiers_1.xml", package = "pubchunks")
(res <- pub_chunks(list(x, y), c("doi", "title", "keywords")))
pub_tabularize(res)
pub_tabularize(res, bind = TRUE)

# using output of fulltext::ft_get()
if (requireNamespace("fulltext", quietly = TRUE)) {
  dois <- c('10.1371/journal.pone.0086169', '10.1371/journal.pone.0155491', 
    '10.7554/eLife.03032')
  x <- fulltext::ft_get(dois)
  (tmp <- pub_chunks(fulltext::ft_collect(x), sections=c("doi","title")))
  pub_tabularize(tmp)
  pub_tabularize(tmp, bind = TRUE)
}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chunks.R
\name{pub_chunks}
\alias{pub_chunks}
\title{Extract chunks of data from articles}
\usage{
pub_chunks(x, sections = "all", provider = NULL, extract = "xml_text")
}
\arguments{
\item{x}{one of the following:
\itemize{
\item file path for an XML file
\item a character string of XML, a list (of file paths, or XML in a character
string, or \code{xml_document} objects)
\item or an object of class \code{fulltext::ft_data}, the output from a call to
\code{fulltext::ft_get()}
}}

\item{sections}{(character) What elements to get, can be one or more in
a vector or list. See \code{\link[=pub_sections]{pub_sections()}} for options. optional. Default is
to get all sections. See Details.}

\item{provider}{(character) a single publisher name. see
\code{\link[=pub_providers]{pub_providers()}} for options. required. If you select the wrong provider
for the XML you have you may or may not get what you need :). By default
this is \code{NULL} and we use \code{\link[=pub_guess_publisher]{pub_guess_publisher()}} to guess the
publisher; we may get it wrong. You can override our guessing by passing
in a name.}

\item{extract}{(character) one of 'xml_text' (default) or 'as.character'.
The final step of extracting each part of an article is converting to
a character string. By default, we'll use \code{\link[xml2:xml_text]{xml2::xml_text()}}, but if you
prefer you can use \code{\link[=as.character]{as.character()}} which. The latter can be useful if
the chunk being extracted has html tags in it that you do not want
removed.}
}
\value{
A list, named by the section selected. sections not found or
not in accepted list return \code{NULL} or zero length list. A ".publisher"
list element gets attached to each list output, even when no
data is found. When \code{fulltext::ft_get} output is passed in here, the
list is named by the publisher, then within each publisher is a list
of articles named by their identifiers (e.g. DOIs).
}
\description{
\code{pub_chunks} makes it easy to extract sections of an article.
You can extract just authors across all articles, or all references
sections, or the complete text of each article. Then you can pass the
output downstream for visualization and analysis.
}
\details{
Options for the \code{sections} parameter:
\itemize{
\item front - Publisher, journal and article metadata elements
\item body - Body of the article
\item back - Back of the article, acknowledgments, author contributions,
references
\item title - Article title
\item doi - Article DOI
\item categories - Publisher's categories, if any
\item authors - Authors
\item aff - Affiliation (includes author names)
\item keywords - Keywords
\item abstract - Article abstract
\item executive_summary - Article executive summary
\item refs - References
\item refs_dois - References DOIs - if available
\item publisher - Publisher name
\item journal_meta - Journal metadata
\item article_meta - Article metadata
\item acknowledgments - Acknowledgments
\item permissions - Article permissions
\item history - Dates, recieved, published, accepted, etc.
}
}
\examples{
# a file path to an XML file
x <- system.file("examples/elsevier_1.xml", package = "pubchunks")
pub_chunks(x, "title")
pub_chunks(x, "authors")
pub_chunks(x, "acknowledgments")
pub_chunks(x, "refs")
pub_chunks(x, c("title", "refs"))

\dontrun{
# works the same with the xml already in a string
xml <- paste0(readLines(x), collapse = "")
pub_chunks(xml, "title")

# also works if you've already read in the XML (with xml2 pkg)
xml <- paste0(readLines(x), collapse = "")
xml <- xml2::read_xml(xml)
pub_chunks(xml, "title")

# Hindawi
x <- system.file("examples/hindawi_1.xml", package = "pubchunks")
pub_chunks(x, "abstract")$abstract
pub_chunks(x, "abstract", extract="as.character")$abstract
pub_chunks(x, "authors")
pub_chunks(x, "aff")
pub_chunks(x, "title")
pub_chunks(x, "refs")$refs
pub_chunks(x, c("abstract", "title", "authors", "refs"))

# Pensoft
x <- system.file("examples/pensoft_1.xml", package = "pubchunks")
pub_chunks(x, "abstract")
pub_chunks(x, "aff")
pub_chunks(x, "title")
pub_chunks(x, "refs")$refs
pub_chunks(x, c("abstract", "title", "authors", "refs"))

# Peerj
x <- system.file("examples/peerj_1.xml", package = "pubchunks")
pub_chunks(x, "abstract")
pub_chunks(x, "authors")
pub_chunks(x, "aff")
pub_chunks(x, "title")
pub_chunks(x, "refs")$refs
pub_chunks(x, c("abstract", "title", "authors", "refs"))

# Frontiers
x <- system.file("examples/frontiers_1.xml", package = "pubchunks")
pub_chunks(x, "authors")
pub_chunks(x, "aff")
pub_chunks(x, "refs")$refs
pub_chunks(x, c("doi", "abstract", "title", "authors", "refs", "abstract"))

# eLife
x <- system.file("examples/elife_1.xml", package = "pubchunks")
pub_chunks(x, "authors")
pub_chunks(x, "aff")
pub_chunks(x, "refs")$refs
pub_chunks(x, c("doi", "title", "authors", "refs"))

# f1000research
x <- system.file("examples/f1000research_3.xml", package = "pubchunks")
pub_chunks(x, "title")
pub_chunks(x, "aff")
pub_chunks(x, "refs")$refs
pub_chunks(x, c("doi", "title", "authors", "keywords", "refs"))

# Copernicus
x <- system.file("examples/copernicus_1.xml", package = "pubchunks")
pub_chunks(x, c("doi", "abstract", "title", "authors", "refs"))
pub_chunks(x, "aff")
pub_chunks(x, "refs")$refs

# MDPI
x <- system.file("examples/mdpi_1.xml", package = "pubchunks")
x <- system.file("examples/mdpi_2.xml", package = "pubchunks")
pub_chunks(x, "title")
pub_chunks(x, "aff")
pub_chunks(x, "refs")$refs
vv <- pub_chunks(x, c("doi", "title", "authors", "keywords", "refs", 
  "abstract", "categories"))
vv$doi
vv$title
vv$authors
vv$keywords
vv$refs
vv$abstract
vv$categories

# Many inputs at once
x <- system.file("examples/frontiers_1.xml", package = "pubchunks")
y <- system.file("examples/elife_1.xml", package = "pubchunks")
z <- system.file("examples/f1000research_1.xml", package = "pubchunks")
pub_chunks(list(x, y, z), c("doi", "title", "authors", "refs"))

# non-XML files/content are xxx?
# pub_chunks('foo bar')

# Pubmed brief XML files (abstract only)
x <- system.file("examples/pubmed_brief_1.xml", package = "pubchunks")
pub_chunks(x, "title")

# Pubmed full XML files
x <- system.file("examples/pubmed_full_1.xml", package = "pubchunks")
pub_chunks(x, "title")

# using output of fulltext::ft_get()
if (requireNamespace("fulltext", quietly = TRUE)) {
  library("fulltext")

  # single
  x <- fulltext::ft_get('10.7554/eLife.03032')
  pub_chunks(fulltext::ft_collect(x), sections="authors")

  # many
  dois <- c('10.1371/journal.pone.0086169', '10.1371/journal.pone.0155491', 
    '10.7554/eLife.03032')
  x <- fulltext::ft_get(dois)
  pub_chunks(fulltext::ft_collect(x), sections="authors")

  # as.ft_data() function
  x <- ft_collect(as.ft_data())
  names(x)
  x$cached
  pub_chunks(x, "title")
  pub_chunks(x, "title") \%>\% pub_tabularize()
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pubchunks-package.R
\docType{package}
\name{pubchunks-package}
\alias{pubchunks-package}
\alias{pubchunks}
\title{Get chunks of XML articles}
\description{
Pass XML in various forms, and get out specific sections
of articles. It makes not knowing what XPath/CSS selectors are totally
fine.
}
\author{
Scott Chamberlain
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pub_guess_publisher.R
\name{pub_guess_publisher}
\alias{pub_guess_publisher}
\title{Guess the publisher from an XML document}
\usage{
pub_guess_publisher(x)
}
\arguments{
\item{x}{an XML file, a character string of XML, or a
\code{xml_document} object (as from \code{xml2::read_xml})}
}
\value{
a list, with two named character strings, one for
\code{full_name} and the other a \code{short_name}
}
\description{
Guess the publisher from an XML document
}
\examples{
\dontrun{
(x <- system.file("examples/pensoft_1.xml", package = "pubchunks"))
pub_guess_publisher(x)

(x <- system.file("examples/copernicus_2.xml", package = "pubchunks"))
pub_guess_publisher(x)

(x <- system.file("examples/peerj_1.xml", package = "pubchunks"))
pub_guess_publisher(x)

(x <- system.file("examples/hindawi_1.xml", package = "pubchunks"))
pub_guess_publisher(x)

(x <- system.file("examples/frontiers_1.xml", package = "pubchunks"))
pub_guess_publisher(x)

(x <- system.file("examples/elife_1.xml", package = "pubchunks"))
pub_guess_publisher(x)

(x <- system.file("examples/elsevier_1.xml", package = "pubchunks"))
pub_guess_publisher(x)

x <- system.file("examples/f1000research_1.xml", package = "pubchunks")
pub_guess_publisher(x)

x <- system.file("examples/plos_1.xml", package = "pubchunks")
pub_guess_publisher(x)

x <- system.file("examples/mdpi_1.xml", package = "pubchunks")
pub_guess_publisher(x)

x <- system.file("examples/pubmed_brief_1.xml", package = "pubchunks")
pub_guess_publisher(x)

x <- system.file("examples/pubmed_full_1.xml", package = "pubchunks")
pub_guess_publisher(x)

x <- system.file("examples/pubmed_full_2.xml", package = "pubchunks")
pub_guess_publisher(x)

x <- system.file("examples/pubmed_full_3.xml", package = "pubchunks")
pub_guess_publisher(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pub_sections.R
\name{pub_sections}
\alias{pub_sections}
\title{Sections}
\usage{
pub_sections()
}
\value{
character vector
}
\description{
The possible sections of an XML article that are supported
for retrieval
}
\examples{
pub_sections()
}
