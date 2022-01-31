namext
======



[![R-check](https://github.com/ropensci/namext/workflows/R-check/badge.svg)](https://github.com/ropensci/namext/actions/)


`namext` - Extract scientific names from text using gnfinder

## Install gnfinder

gnfinder is required to use this package.

Installation instructions can be found at the [gnfinder repo](https://github.com/gnames/gnfinder). 

## Install namext


```r
remotes::install_github("ropensci/namext")
```


```r
library("namext")
```

## extract names


```r
x <- system.file("examples/BoegerEtal2015.pdf", package="namext")
name_extract(x)
#> $metadata
#> # A tibble: 1 x 9
#>   date  gnfinderVersion withBayes tokensAround language detectLanguage
#>   <chr> <chr>           <lgl>            <int> <chr>    <lgl>         
#> 1 2020… v0.11.1-2-g551… TRUE                 0 eng      FALSE         
#> # … with 3 more variables: totalWords <int>, totalCandidates <int>,
#> #   totalNames <int>
#> 
#> $names
#> # A tibble: 331 x 8
#>    cardinality verbatim   name      odds start   end annotationNomen… annotation
#>          <int> <chr>      <chr>    <dbl> <int> <int> <chr>            <chr>     
#>  1           1 (Teleoste… Teleos…   145.    84    95 NO_ANNOT         ""        
#>  2           1 Sciaenida… Sciaen… 38075.    96   107 NO_ANNOT         ""        
#>  3           1 MARIANA    Mariana  2171.   188   195 NO_ANNOT         ""        
#>  4           1 (Teleoste… Teleos…   145.   398   409 NO_ANNOT         ""        
#>  5           1 Sciaenida… Sciaen… 38075.   410   421 NO_ANNOT         ""        
#>  6           1 Pachyurin… Pachyu… 10451.   673   684 NO_ANNOT         ""        
#>  7           1 (Sciaenid… Sciaen… 38075.   685   697 NO_ANNOT         ""        
#>  8           1 Sciaenida… Sciaen… 38075.  1047  1058 NO_ANNOT         ""        
#>  9           1 Haemulidae Haemul… 38075.  1059  1069 NO_ANNOT         ""        
#> 10           1 Polypteri… Polypt… 86951.  1074  1086 NO_ANNOT         ""        
#> # … with 321 more rows
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/namext/issues)
* License: MIT
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

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

Please note that the {{{ package }}} project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
namext
======

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

[![R-check](https://github.com/ropensci/namext/workflows/R-check/badge.svg)](https://github.com/ropensci/namext/actions/)


`namext` - Extract scientific names from text using gnfinder

## Install gnfinder

gnfinder is required to use this package.

Installation instructions can be found at the [gnfinder repo](https://github.com/gnames/gnfinder). 

## Install namext

```{r eval=FALSE}
remotes::install_github("ropensci/namext")
```

```{r}
library("namext")
```

## extract names

```{r}
x <- system.file("examples/BoegerEtal2015.pdf", package="namext")
name_extract(x)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/namext/issues)
* License: MIT
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_extract.R
\name{name_extract}
\alias{name_extract}
\title{name_extract}
\usage{
name_extract(
  path,
  verify = FALSE,
  language = NULL,
  no_bayes = FALSE,
  odds_details = FALSE,
  data_sources = NULL,
  tokens = NULL
)
}
\arguments{
\item{path}{(character) path to a file}

\item{verify}{(logical) verify names? default: \code{FALSE}}

\item{language}{(logical) text's language. by default it's
automatically detected. default: \code{NULL}}

\item{no_bayes}{(logical) do not run Bayes algorithms. default: \code{FALSE}}

\item{odds_details}{(logical) show details of odds calculation.
default: \code{FALSE}}

\item{data_sources}{(numeric vector) IDs of data sources to display for
matches. default: \code{NULL}. e.g., \code{c(1, 11, 179)}}

\item{tokens}{(logical) xxx. default: \code{NULL}}
}
\value{
list of two data.frames:
\itemize{
\item meta: metadata
\item names: names and their parts, varies based on function parameters
}
}
\description{
extract names using gnfinder
}
\examples{
x <- system.file("examples/BoegerEtal2015.pdf", package="namext")
x
name_extract(x)
name_extract(x, verify = TRUE)
name_extract(x, language = "eng")
name_extract(x, no_bayes = TRUE)
name_extract(x, odds_details = TRUE)
name_extract(x, data_sources = c(4, 12))
name_extract(x, tokens = 5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/namext-package.R
\docType{package}
\name{namext-package}
\alias{namext-package}
\alias{namext}
\title{namext}
\description{
Extract scientific names from text using gnfinder
}
\keyword{package}
