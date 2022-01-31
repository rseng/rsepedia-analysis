
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bib2df)](https://cran.r-project.org/package=bib2df) [![Travis-CI Build Status](https://travis-ci.org/ropensci/bib2df.svg?branch=master)](https://travis-ci.org/ropensci/bib2df) [![Build status](https://ci.appveyor.com/api/projects/status/6k3q7272ddnjh20o?svg=true)](https://ci.appveyor.com/project/ottlngr/bib2df) [![](http://cranlogs.r-pkg.org/badges/bib2df)](http://cran.rstudio.com/web/packages/bib2df/index.html) [![codecov](https://codecov.io/gh/ropensci/bib2df/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/bib2df) [![](https://badges.ropensci.org/124_status.svg)](https://github.com/ropensci/onboarding/issues/124)

`bib2df` - Parse a BibTeX file to a tibble
------------------------------------------

Everyone writing reports and articles with LaTeX has probably used BibTeX before. BibTeX is the de facto standard for reference management and grounds its functionality on a list of references stored in local text file. Depending on the reference type, several fields are necessary to define a reference properly. An exemplary BibTeX entry looks as follows:

    @Article{Binmore2008,
      Title = {Do Conventions Need to Be Common Knowledge?},
      Author = {Binmore, Ken},
      Journal = {Topoi},
      Year = {2008},
      Number = {1},
      Pages = {17--27},
      Volume = {27}
    }

Parse the BibTeX file to a tibble
---------------------------------

The BibTeX format is not convenient for any kind of analysis or visualization. Many R applications require a `data.frame` (or `tibble`) and `bib2df` offers a straightforward framework to parse a BibTeX file to a `tibble`.

``` r
library(bib2df)

path <- system.file("extdata", "LiteratureOnCommonKnowledgeInGameTheory.bib", package = "bib2df")

df <- bib2df(path)
df
#> # A tibble: 37 x 27
#>    CATEGORY BIBTEXKEY ADDRESS ANNOTE AUTHOR BOOKTITLE CHAPTER CROSSREF
#>    <chr>    <chr>     <chr>   <chr>  <list> <chr>     <chr>   <chr>   
#>  1 ARTICLE  Arrow1986 <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#>  2 ARTICLE  AumannBr… <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#>  3 ARTICLE  Aumann19… <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#>  4 INCOLLE… Bacharac… Cambri… <NA>   <chr … Knowledg… 17      <NA>    
#>  5 ARTICLE  Basu1988  <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#>  6 ARTICLE  Bernheim… <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#>  7 ARTICLE  Bicchier… <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#>  8 ARTICLE  Binmore2… <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#>  9 ARTICLE  Brandenb… <NA>    <NA>   <chr … <NA>      <NA>    <NA>    
#> 10 INCOLLE… Brandenb… New Yo… <NA>   <chr … The Econ… 3       <NA>    
#> # … with 27 more rows, and 19 more variables: EDITION <chr>,
#> #   EDITOR <list>, HOWPUBLISHED <chr>, INSTITUTION <chr>, JOURNAL <chr>,
#> #   KEY <chr>, MONTH <chr>, NOTE <chr>, NUMBER <chr>, ORGANIZATION <chr>,
#> #   PAGES <chr>, PUBLISHER <chr>, SCHOOL <chr>, SERIES <chr>, TITLE <chr>,
#> #   TYPE <chr>, VOLUME <chr>, YEAR <dbl>, DOI <chr>
```

The `df2bib()` function makes it possible to write this `tibble` back to disk, enabling programmatic manipulation of a .bib file.

Installation
------------

The latest version of `bib2df` can be installed from GitHub using `devtools::install_github()`:

    devtools::install_github("ropensci/bib2df")

Version 1.1.1 is now available on **CRAN**:

    install.packages("bib2df")

Community Guidelines
--------------------

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

------------------------------------------------------------------------

[![ropensci\_footer](./ropensci_footer.png)](https://ropensci.org)
# bib2df 1.1.0

* This version fixes failing unit tests on CRAN.
* This version will be published on CRAN.

# bib2df 1.0.1

* This package versions fixes some issues:
    * `bib2df()` can now read `.bib` files with a single entry
    * `bib2df()` now allows `@` symbol within titles
* This version will be published on CRAN since version 1.0.0 had beed removed from CRAN due to check errors.

# bib2df 1.0.0

* This package version went successfully through the `rOpenSci` onboarding process.
* The package now support AppVeyor builds.
* This version contains aditional unit tests.

# bib2df 0.2.2

* Changes due to the ropensci review process:
    * Improved behavior on malformed .bib files and entries
    * Improved behavior on type conversion
    * `df2bib()` is now able to append to an existing file
    * Minor changes in documentation and vignette

# bib2df 0.2.1

* Changes due to the ropensci review process:
    * Improved formatting
    * Removed exports of internal functions
    * Improved vignette

# bib2df 0.2

* Standardized return value for `bib2df()` (#5, @leeper)
* Added functionality to write `data.frame` back to BibTeX file (`df2bib()`) (#5, @leeper)
* Adds a testthat-based tests for both `bib2df()` and `df2bib()` (#5, @leeper)
* Added functionality of the `humaniformat` package to split up names (#4)
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
## Test environments

* local Ubuntu 18.04, R 3.4.4, via
* Travis-CI
* AppVeyor
* R-hub builder

## R CMD check results

Status: 0 errors | 0 warnings | 1 notes

```
checking CRAN incoming feasibility ... NOTE
  
  Maintainer: 'Philipp Ottolinger <philipp@ottolinger.de>'
  
  New submission
  
  Package was archived on CRAN
```

## Resubmission

Resubmission due to failing unit tests: https://cran-archive.r-project.org/web/checks/2019-04-02_check_results_bib2df.html
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bib2df)](https://cran.r-project.org/package=bib2df) [![Travis-CI Build Status](https://travis-ci.org/ropensci/bib2df.svg?branch=master)](https://travis-ci.org/ropensci/bib2df) [![Build status](https://ci.appveyor.com/api/projects/status/6k3q7272ddnjh20o?svg=true)](https://ci.appveyor.com/project/ottlngr/bib2df) [![](http://cranlogs.r-pkg.org/badges/bib2df)](http://cran.rstudio.com/web/packages/bib2df/index.html) [![codecov](https://codecov.io/gh/ropensci/bib2df/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/bib2df) [![](https://badges.ropensci.org/124_status.svg)](https://github.com/ropensci/onboarding/issues/124)





## `bib2df` - Parse a BibTeX file to a tibble

Everyone writing reports and articles with LaTeX has probably used BibTeX before. BibTeX is the de facto standard for reference management and grounds its functionality on a list of references stored in local text file. Depending on the reference type, several fields are necessary to define a reference properly. An exemplary BibTeX entry looks as follows:

```
@Article{Binmore2008,
  Title = {Do Conventions Need to Be Common Knowledge?},
  Author = {Binmore, Ken},
  Journal = {Topoi},
  Year = {2008},
  Number = {1},
  Pages = {17--27},
  Volume = {27}
}
```

## Parse the BibTeX file to a tibble

The BibTeX format is not convenient for any kind of analysis or visualization. Many R applications require a `data.frame` (or `tibble`) and `bib2df` offers a straightforward framework to parse a BibTeX file to a `tibble`.


```{r, warning=FALSE}
library(bib2df)

path <- system.file("extdata", "LiteratureOnCommonKnowledgeInGameTheory.bib", package = "bib2df")

df <- bib2df(path)
df
```

The `df2bib()` function makes it possible to write this `tibble` back to disk, enabling programmatic manipulation of a .bib file.

## Installation

The latest version of `bib2df` can be installed from GitHub using `devtools::install_github()`:

```
devtools::install_github("ropensci/bib2df")
```

Version 1.1.1 is now available on **CRAN**:

```
install.packages("bib2df")
```

## Community Guidelines

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

------------------------------------------
[![ropensci_footer](./ropensci_footer.png)](https://ropensci.org)
---
title: "bib2df - Parse a BibTeX file to a data.frame"
author: "Philipp Ottolinger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{bib2df}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE
)
```
 
## BibTeX

BibTeX is typically used together with LaTeX to manage references. The BibTeX file format is simple as it follows rather simple but strict rules and represents a reference's core data as a list of partly mandatory fields. 

The resulting BibTeX file can tell much about the work you use it for, may it be an academic paper, a dissertation or any other report that at least partly appoints to the work of others. The average age of the referenced works might tell if one addresses to a current topic or if one digs into the history of a certain field. Does one cite many works of just a few authors or occurs every author at most once? The BibTeX file is definitely able to answer these questions.

## Why bib2df?

As mentioned above, BibTeX represents the entries as lists of named fields, some kind of similar to the `JSON` format. If you want to gain insights from your BibTeX file using R, you will have to transform the data to fit into a more usual data structure. Such a data structure, speaking of R,  is the `tibble`. `bib2df` does exactly this: It takes a BibTeX file and parses it into a `tibble` so you can work with your bibliographic data just the way you do with other data.

Given this `tibble` you can manipulate entries in a familiar way and write the updated references back to a valid BibTeX file. 

## How to use

### Parse the BibTeX file

To parse a BibTeX file to a `tibble` you may want to use the function `bib2df()`. The first argument of `bib2df()` is the path to the file you want to parse.

```{r, eval = FALSE}
install.packages("bib2df")
```

```{r}
library(bib2df)

path <- system.file("extdata", "LiteratureOnCommonKnowledgeInGameTheory.bib", package = "bib2df")

df <- bib2df(path)
df
```

`bib2df()` returns a `tibble` with each row representing one entry of the initial BibTeX file while the columns hold the data originally stored in the named fields. If a field was not present in a particular entry, the respective column gets the value `NA`. As some works can be the work of multiple authors or editors, these fields are converted to a list to avoid having the names of multiple persons concatenated to a single character string:

```{r}
head(df$AUTHOR)
```

The second argument of `bib2df()` is `separate_names` and calls, if `TRUE`, the functionality of the `humaniformat` [^1] package to automatically split persons' names into pieces:

```{r}
df <- bib2df(path, separate_names = TRUE)
head(df$AUTHOR)
```

### Parsing multiple BibTex files

Multiple BibTeX files can be parsed using `lapply()`. The paths to the BibTeX files must be stored in a vector. Using this vector to call `bib2df()` within `lapply()` results in a list of `tibble`, which can be bound, e.g. using `bind_rows()`:

```{r}
bib1 <- system.file("extdata", "LiteratureOnCommonKnowledgeInGameTheory.bib", package = "bib2df")
bib2 <- system.file("extdata", "r.bib", package = "bib2df")
paths <- c(bib1, bib2)

x <- lapply(paths, bib2df)
class(x)
head(x)

res <- dplyr::bind_rows(x)
class(res)
head(res)
```


## Analyze and visualize your references

Since the BibTeX entries are now converted to rows and columns in a `tibble`, one can start to analyze and visualize the data with common tools like `ggplot2`, `dplyr` and `tidyr`.

For example, one can ask which journal is cited most among the references

```{r, message=FALSE}
library(dplyr)
library(ggplot2)
library(tidyr)

df %>%
  filter(!is.na(JOURNAL)) %>%
  group_by(JOURNAL) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  slice(1:3)
```

or what the median age of the cited works is:

```{r}
df %>%
  mutate(age = 2017 - YEAR) %>%
  summarize(m = median(age))
```

Also plotting is possible:

```{r, fig.height = 5, fig.width = 7}
df %>% 
  select(YEAR, AUTHOR) %>% 
  unnest() %>% 
  ggplot() + 
  aes(x = YEAR, y = reorder(full_name, desc(YEAR))) + 
  geom_point()
```


## Manipulate your references

Since all the BibTeX entries are represented by rows in a `tibble`, manipulating the BibTeX entries is now easily possible.

One of the authors of the 10th reference in our file does not have his full first name:

```{r}
df$AUTHOR[[10]]
```

The 'E.' in 'E. Dekel' is for Eddie, so lets change the value of that field:

```{r}
df$AUTHOR[[10]]$first_name[2] <- "Eddie"
df$AUTHOR[[10]]$full_name[2] <- "Eddie Dekel"

df$AUTHOR[[10]]
```


## Write back to BibTeX

Especially when single values of the parsed BibTeX file were changed it is useful to write the parsed `tibble` back to a valid BibTeX file one can use in combination with LaTeX. Just like `bib2df()` parses a BibTeX file, `df2bib()` writes a BibTeX file:

```{r,eval=FALSE}
newFile <- tempfile()
df2bib(df, file = newFile)
```

The just written BibTeX file of course contains the values, we just changed in the `tibble`:

```
@Incollection{BrandenburgerDekel1989,
  Address = {New York},
  Author = {Brandenburger, Adam and Dekel, Eddie},
  Booktitle = {The Economics of Missing Markets, Information and Games},
  Chapter = {3},
  Pages = {46 - 61},
  Publisher = {Oxford University Press},
  Title = {The Role of Common Knowledge Assumptions in Game Theory},
  Year = {1989}
}
```

To append BibTeX entries to an existing file use `append = TRUE` within `df2bib()`. 

[^1]: <https://CRAN.R-project.org/package=humaniformat>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bib2df-package.R
\docType{package}
\name{bib2df-package}
\alias{bib2df-package}
\title{bib2df: Parse a BibTeX file to \code{data.frame}}
\description{
This package provides functions to parse and write BibTeX files. BibTeX files
can be parsed to \code{data.frame} to make them accessible with popular tools like
\code{dplyr}, \code{tidyr}, \code{ggplot2} and many more.
}
\details{
BibTeX entries represented in a \code{data.frame} can be altered in a familiar way and
written back to a valid BibTeX file.

To learn more about \code{bib2df}, start with the vignettes:
 \code{browseVignettes(package = "bib2df")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bib2df.R
\name{bib2df}
\alias{bib2df}
\title{Parse a BibTeX file to a \code{tibble}}
\usage{
bib2df(file, separate_names = FALSE)
}
\arguments{
\item{file}{character, path or URL to a .bib file.}

\item{separate_names}{logical, should authors' and editors' names be separated into first and given name?}
}
\value{
A \code{tibble}.
}
\description{
The BibTeX file is read, parsed, tidied and written to a \code{tibble}
}
\details{
For simplicity \code{bib2df()} unifies the reading, parsing and tidying of a BibTeX file while being aware of a standardized output format, different BibTeX styles and missing values in the BibTeX file.

When \code{separate_names = TRUE}, the respective columns contain a \code{data.frame} for each row. When \code{FALSE}, the respective columns contain character strings.
}
\examples{
# Read from .bib file:
path <- system.file("extdata", "bib2df_testfile_3.bib", package = "bib2df")
bib <- bib2df(path)
str(bib)

# Read from .bib file and separate authors' and editors' names:
bib <- bib2df(path, separate_names = TRUE)
str(bib)
}
\seealso{
\code{\link{df2bib}}
}
\author{
Philipp Ottolinger
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df2bib.R
\name{df2bib}
\alias{df2bib}
\title{Export a BibTeX \code{tibble} to a .bib file}
\usage{
df2bib(x, file = "", append = FALSE)
}
\arguments{
\item{x}{\code{tibble}, in the format as returned by \code{\link{bib2df}}.}

\item{file}{character, file path to write the .bib file. An empty character string writes to \code{stdout} (default).}

\item{append}{logical, if \code{TRUE} the \code{tibble} will be appended to an existing file.}
}
\value{
\code{file} as a character string, invisibly.
}
\description{
The BibTeX \code{tibble} is written to a .bib file
}
\examples{
# Read from .bib file:
path <- system.file("extdata", "bib2df_testfile_3.bib", package = "bib2df")
bib <- bib2df(path)

# Write to .bib file:
# bibFile <- tempfile()
# df2bib(bib, bibFile)

# Use `append = TRUE` to add lines to an existing .bib file:
# df2bib(bib, bibFile, append = TRUE)
}
\references{
\url{http://www.bibtex.org/Format/}
}
\seealso{
\code{\link{bib2df}}
}
\author{
Thomas J. Leeper
}
