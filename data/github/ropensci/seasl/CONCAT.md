seasl
=====



[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
![R build status](https://github.com/ropensci/seasl/workflows/R-CMD-check/badge.svg)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/seasl?color=C9A115)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/seasl)](https://cran.r-project.org/package=seasl)
<!-- [![Build Status](https://travis-ci.org/ropensci/seasl.svg?branch=master)](https://travis-ci.org/ropensci/seasl) -->

`seasl` is an R client for exploring CSL styles.

This package is inspired by the Ruby gem `csl`: https://github.com/inukshuk/csl-ruby

The Citation Style Language 1.0.1 specification: http://docs.citationstyles.org/en/1.0.1/specification.html

Package API:

 - `csl_locales`
 - `as.location`
 - `csl_styles`
 - `csl_locale_exists`
 - `csl_cache`
 - `csl_fetch_styles`
 - `csl_style_find`
 - `csl_style_xml`
 - `csl_style_load`
 - `csl_style_exists`
 - `csl_locale_load`
 - `csl_fetch_locales`

## Install


```r
install.packages("seasl")
```

or


```r
remotes::install_github("ropensci/seasl")
```


```r
library("seasl")
```

## Download styles and locales

First, you may want to download style and locale files. `csl_fetch_styles()` and `csl_fetch_locales()`
download the files to your machine. See `?csl_cache` for caching information, including
how to change the cache location.

Styles retrieved from the Github repo at https://github.com/citation-style-language/styles-distribution


```r
csl_fetch_styles()
#>
#> Done! Files put in /Users/sckott/Library/Caches/R/seasl/styles
```

Locales retrieved from the Github repo at https://github.com/citation-style-language/locales


```r
csl_fetch_locales()
#>
#> Done! Files put in /Users/sckott/Library/Caches/R/seasl/locales
```

## File paths to CSL styles and locales

calling `csl_styles` without inputs gives all styles, with separate lists for 
dependent and independent styles


```r
csl_styles()
#> $independent
#>    [1] "academy-of-management-review"                                                                                        
#>    [2] "accident-analysis-and-prevention"                                                                                    
#>    [3] "aci-materials-journal"                                                                                               
#>    [4] "acm-sig-proceedings-long-author-list"                                                                                
#>    [5] "acm-sig-proceedings"                                                                                                 
#>    [6] "acm-sigchi-proceedings-extended-abstract-format"                                                                     
#>    [7] "acm-sigchi-proceedings"                                                                                              
#>    [8] "acm-siggraph"                                                                                                        
#>    [9] "acme-an-international-journal-for-critical-geographies"                                                              
...
```

calling `csl_styles` with an input gives the path to that style, if found


```r
csl_styles("apa")
#> [1] "/Users/sckott/Library/Caches/R/seasl/styles/apa.csl"
csl_styles("archeologie-medievale")
#> [1] "/Users/sckott/Library/Caches/R/seasl/styles/archeologie-medievale.csl"
```

Same patterns go for locales (note that there are far fewer locales than styles)


```r
# just locale names
csl_locales()
#>  [1] "locales-af-ZA" "locales-ar"    "locales-bg-BG" "locales-ca-AD"
#>  [5] "locales-cs-CZ" "locales-cy-GB" "locales-da-DK" "locales-de-AT"
#>  [9] "locales-de-CH" "locales-de-DE" "locales-el-GR" "locales-en-GB"
#> [13] "locales-en-US" "locales-es-CL" "locales-es-ES" "locales-es-MX"
#> [17] "locales-et-EE" "locales-eu"    "locales-fa-IR" "locales-fi-FI"
#> [21] "locales-fr-CA" "locales-fr-FR" "locales-he-IL" "locales-hr-HR"
#> [25] "locales-hu-HU" "locales-id-ID" "locales-is-IS" "locales-it-IT"
#> [29] "locales-ja-JP" "locales-km-KH" "locales-ko-KR" "locales-la"   
#> [33] "locales-lt-LT" "locales-lv-LV" "locales-mn-MN" "locales-nb-NO"
#> [37] "locales-nl-NL" "locales-nn-NO" "locales-pl-PL" "locales-pt-BR"
#> [41] "locales-pt-PT" "locales-ro-RO" "locales-ru-RU" "locales-sk-SK"
#> [45] "locales-sl-SI" "locales-sr-RS" "locales-sv-SE" "locales-th-TH"
#> [49] "locales-tr-TR" "locales-uk-UA" "locales-vi-VN" "locales-zh-CN"
#> [53] "locales-zh-TW"
```


```r
# when locale given, you get the full path
csl_locales("fr-FR")
#> [1] "/Users/sckott/Library/Caches/R/seasl/locales/locales-fr-FR.xml"
```

Alternatively, you can try to find a style by using `csl_style_find()`


```r
# single match
csl_style_find(x = "American Journal of Epidemiology")
#> [1] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-epidemiology.csl"
```


```r
# many matches
csl_style_find(x = "American Journal")
#>  [1] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-agricultural-economics.csl"                                     
#>  [2] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-archaeology.csl"                                                
#>  [3] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-botany.csl"                                                     
#>  [4] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-climate-change.csl"                                             
#>  [5] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-clinical-pathology.csl"                                         
#>  [6] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-enology-and-viticulture.csl"                                    
#>  [7] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-epidemiology.csl"                                               
#>  [8] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-health-behavior.csl"                                            
#>  [9] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-hypertension.csl"                                               
#> [10] "/Users/sckott/Library/Caches/R/seasl/styles/american-journal-of-medical-genetics.csl"                                           
...
```

## Load CSL style from a URL


```r
jps <- csl_style_load('http://www.zotero.org/styles/american-journal-of-political-science')
```

## Query style information


```r
jps$info
#> $title
#> [1] "American Journal of Political Science"
#> 
#> $`title-short`
#> [1] "AJPS"
#> 
#> $id
#> [1] "http://www.zotero.org/styles/american-journal-of-political-science"
#> 
#> $contributor
...
```


```r
jps$info$title
#> [1] "American Journal of Political Science"
```


```r
jps$macros
#> [[1]]
#> [[1]]$name
#> [1] "editor"
#> 
#> [[1]][[2]]
#> [[1]][[2]]$names
#> [[1]][[2]]$names$variable
#> [1] "editor"
#> 
#> [[1]][[2]]$names$delimiter
...
```


```r
jps$citation
#> $sort
#> $sort$key
#> list()
#> attr(,"macro")
#> [1] "author-short"
#> 
#> $sort$key
#> list()
#> attr(,"macro")
#> [1] "year-date"
...
```


```r
jps$bibliography
#> $attributes
#> $attributes$`hanging-indent`
#> [1] "true"
#> 
#> $attributes$`et-al-min`
#> [1] "4"
#> 
#> $attributes$`et-al-use-first`
#> [1] "1"
#> 
...
```

## Get raw XML


```r
csl_style_xml('http://www.zotero.org/styles/american-journal-of-political-science')
#> {xml_document}
#> <style class="in-text" version="1.0" demote-non-dropping-particle="sort-only" default-locale="en-US" xmlns="http://purl.org/net/xbiblio/csl">
#>  [1] <info>\n  <title>American Journal of Political Science</title>\n  <title ...
#>  [2] <macro name="editor">\n  <names variable="editor" delimiter=", ">\n    < ...
#>  [3] <macro name="author">\n  <names variable="author">\n    <name name-as-so ...
#>  [4] <macro name="author-short">\n  <names variable="author">\n    <name form ...
#>  [5] <macro name="access">\n  <choose>\n    <if type="legal_case" match="none ...
#>  [6] <macro name="title">\n  <choose>\n    <if type="bill book graphic legal_ ...
#>  [7] <macro name="legal_case">\n  <group prefix=" " delimiter=" ">\n    <text ...
#>  [8] <macro name="publisher">\n  <choose>\n    <if type="thesis" match="none" ...
...
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/seasl/issues).
* License: MIT
* Citation: execute `citation(package = 'seasl')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/seasl/blob/master/CODE_OF_CONDUCT.md
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
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
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that seasl is released with a [Contributor Code of Conduct][coc]. By contributing to this project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html) for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 

[coc]: https://github.com/ropensci/seasl/blob/master/CODE_OF_CONDUCT.md
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
seasl
=====

```{r echo=FALSE}
library("knitr")
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
![R build status](https://github.com/ropensci/seasl/workflows/R-CMD-check/badge.svg)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/seasl?color=C9A115)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/seasl)](https://cran.r-project.org/package=seasl)
<!-- [![Build Status](https://travis-ci.org/ropensci/seasl.svg?branch=master)](https://travis-ci.org/ropensci/seasl) -->

`seasl` is an R client for exploring CSL styles.

This package is inspired by the Ruby gem `csl`: https://github.com/inukshuk/csl-ruby

The Citation Style Language 1.0.1 specification: http://docs.citationstyles.org/en/1.0.1/specification.html

Package API:

```{r echo=FALSE, comment=NA, results='asis'}
# cat(paste(" -", paste(sort(getNamespaceExports("seasl")), collapse = "\n - ")))
cat(paste(" -", paste(sprintf("`%s`", getNamespaceExports("seasl")), collapse = "\n - ")))
```

## Install

```{r eval=FALSE}
install.packages("seasl")
```

or

```{r eval=FALSE}
remotes::install_github("ropensci/seasl")
```

```{r}
library("seasl")
```

## Download styles and locales

First, you may want to download style and locale files. `csl_fetch_styles()` and `csl_fetch_locales()`
download the files to your machine. See `?csl_cache` for caching information, including
how to change the cache location.

Styles retrieved from the Github repo at https://github.com/citation-style-language/styles-distribution

```{r eval=FALSE}
csl_fetch_styles()
#>
#> Done! Files put in /Users/sckott/Library/Caches/R/seasl/styles
```

Locales retrieved from the Github repo at https://github.com/citation-style-language/locales

```{r eval=FALSE}
csl_fetch_locales()
#>
#> Done! Files put in /Users/sckott/Library/Caches/R/seasl/locales
```

## File paths to CSL styles and locales

calling `csl_styles` without inputs gives all styles, with separate lists for 
dependent and independent styles

```{r output.lines = 1:10}
csl_styles()
```

calling `csl_styles` with an input gives the path to that style, if found

```{r}
csl_styles("apa")
csl_styles("archeologie-medievale")
```

Same patterns go for locales (note that there are far fewer locales than styles)

```{r}
# just locale names
csl_locales()
```

```{r}
# when locale given, you get the full path
csl_locales("fr-FR")
```

Alternatively, you can try to find a style by using `csl_style_find()`

```{r}
# single match
csl_style_find(x = "American Journal of Epidemiology")
```

```{r output.lines = 1:10}
# many matches
csl_style_find(x = "American Journal")
```

## Load CSL style from a URL

```{r}
jps <- csl_style_load('http://www.zotero.org/styles/american-journal-of-political-science')
```

## Query style information

```{r output.lines = 1:10}
jps$info
```

```{r}
jps$info$title
```

```{r output.lines = 1:10}
jps$macros
```

```{r output.lines = 1:10}
jps$citation
```

```{r output.lines = 1:10}
jps$bibliography
```

## Get raw XML

```{r output.lines=1:10}
csl_style_xml('http://www.zotero.org/styles/american-journal-of-political-science')
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/seasl/issues).
* License: MIT
* Citation: execute `citation(package = 'seasl')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/seasl/blob/master/CODE_OF_CONDUCT.md
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_locales.R
\name{csl_fetch_locales}
\alias{csl_fetch_locales}
\title{Get CSL locales from the web}
\usage{
csl_fetch_locales(...)
}
\arguments{
\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
path (character) to the files (invisibly)
}
\description{
Pulls down all CSL locale files into a directory on your machine,
returning the path to that directory
}
\details{
Files are stored in
\code{paste0(rappdirs::user_cache_dir(), "/R/seasl/styles")}. See \link{csl_cache}
for more
}
\examples{
\dontrun{
csl_fetch_locales()
}
}
\references{
https://github.com/citation-style-language/locales
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caching.R
\name{csl_cache}
\alias{csl_cache}
\title{Caching}
\description{
Manage cached \code{seasl} files with \pkg{hoardr}
}
\details{
The dafault cache directory is
\code{paste0(rappdirs::user_cache_dir(), "/R/seasl")}, but you can set
your own path using \code{cache_path_set()}

\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}
\section{Useful user functions}{

\itemize{
\item \code{csl_cache$cache_path_get()} get cache path
\item \code{csl_cache$cache_path_set()} set cache path
\item \code{csl_cache$list()} returns a character vector of full
path file names
\item \code{csl_cache$files()} returns file objects with metadata
\item \code{csl_cache$details()} returns files with details
\item \code{csl_cache$delete()} delete specific files
\item \code{csl_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
csl_cache

# list files in cache
csl_cache$list()

# delete certain database files
# csl_cache$delete("file path")
# csl_cache$list()

# delete all files in cache
# csl_cache$delete_all()
# csl_cache$list()

# set a different cache path from the default
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/style_load.R
\name{csl_style_load}
\alias{csl_style_load}
\title{Load a CSL style}
\usage{
csl_style_load(input, ...)
}
\arguments{
\item{input}{URL or local file path}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
named list, including slots for
\itemize{
\item info: basic top level information
\item locale: locale information, if it exists
\item macros: macros, 1 to many
\item citation: the citation format, very messy now as the format is messy,
parsed with \link[xml2:as_list]{xml2::as_list}
\item bibliography: very messy for now, we just run \link[xml2:as_list]{xml2::as_list} on this
element as it's complicated to parse
}
}
\description{
Load a CSL style
}
\details{
This function fetches the style XML document, and parses it into
a more reasonble R list that's easy to navigate. If you want the raw XML,
see \code{\link[=csl_style_xml]{csl_style_xml()}}
}
\examples{
# Load a style from the Zotero style repository
x <- 'http://www.zotero.org/styles/american-journal-of-political-science'
if (crul::ok(x)) {
jps <- csl_style_load(x)

## Query style information
jps$info
jps$locale
jps$macros
jps$citation
jps$bibliography
}

\dontrun{
# fetch styles
csl_fetch_styles()
# Load from a local style file
## just specify the style and we read from the local style files
csl_style_load(input="apa")
csl_style_load("computer-und-recht")
csl_style_load("bulletin-de-correspondance-hellenique")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locales.R
\name{csl_locales}
\alias{csl_locales}
\alias{csl_locale_exists}
\title{List locally stored locales}
\usage{
csl_locales(locale = NULL)

csl_locale_exists(locale)
}
\arguments{
\item{locale}{(character) one locale name}
}
\value{
If \code{locale=NULL}, a list of locales. If \code{locale} is
not \code{NULL}, then a full path to the locale file is returned if the
locale exists.
}
\description{
List locally stored locales
}
\examples{
# setup
csl_cache$cache_path_set("seasl", type = "tempdir")
csl_cache$cache_path_get()

# List all locale files
csl_locales()

# cleanup
csl_cache$delete_all()

\dontrun{
# fetch data first
csl_fetch_locales()

# List all locale files
csl_locales()

# list files
csl_locales("et")
csl_locales("fr-FR")

csl_locale_exists("et")
csl_locale_exists("cc")
csl_locale_exists("fr-FR")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/style_xml.R
\name{csl_style_xml}
\alias{csl_style_xml}
\title{Get XML for a CSL style}
\usage{
csl_style_xml(input, raw = FALSE, ...)
}
\arguments{
\item{input}{(character) URL or local file path. Required.}

\item{raw}{(logical) If \code{FALSE} (default) return parsed XML to class
\code{xml_document}. If \code{TRUE}, get character string of XML.}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
an object of class \code{xml_document}, see \pkg{xml2}
to parse the object
}
\description{
Get XML for a CSL style
}
\details{
This function fetches the style XML document. If you want
parsed data, see \code{\link[=csl_style_load]{csl_style_load()}}.
}
\examples{
ajps <- 'http://zotero.org/styles/american-journal-of-political-science'
if (crul::ok(ajps)) {
csl_style_xml(ajps)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locale_load.R
\name{csl_locale_load}
\alias{csl_locale_load}
\title{Load a CSL locale}
\usage{
csl_locale_load(input, ...)
}
\arguments{
\item{input}{URL or local file path}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Load a CSL locale
}
\details{
This function fetches the style XML document, and parses it into
a more reasonble R list that's easy to navigate. If you want the raw XML,
see \code{locale_xml}
}
\examples{
\dontrun{
# Load a locale from the CSL github repo
de_DE <- 'https://github.com/citation-style-language/locales/raw/master/locales-de-DE.xml'
res <- csl_locale_load(de_DE)
## Query style information
res$info
res$info$translators
res$info$license
res$info$date_updated
res$dates
res$style_options
res$terms

# Load from a local style file
## just specify the style and we read from the local style files
csl_locale_load(input="fr-FR")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_styles.R
\name{csl_fetch_styles}
\alias{csl_fetch_styles}
\title{Get CSL styles from the web}
\usage{
csl_fetch_styles(...)
}
\arguments{
\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
path (character) to the files (invisibly)
}
\description{
Pulls down all CSL style files into a directory on your machine,
returning the path to that directory
}
\details{
Files are stored in
\code{paste0(rappdirs::user_cache_dir(), "/R/seasl/styles")}. See \link{csl_cache}
for more
}
\examples{
\dontrun{
csl_fetch_styles()
}
}
\references{
https://github.com/citation-style-language/styles-distribution
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.location.R
\name{as.location}
\alias{as.location}
\title{Convert a path or URL to a location object.}
\usage{
as.location(x, type = "style")
}
\arguments{
\item{x}{Input}

\item{type}{(character) One of style (default) or locale}
}
\description{
Convert a path or URL to a location object.
}
\examples{
if (length(csl_cache$list()) > 0) {
# Style files
as.location("apa")
as.location("teaching-and-learning-in-nursing")
as.location("regenerative-medicine-research")

# A URL
url <- 'http://zotero.org/styles/american-journal-of-political-science'
as.location(url)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/styles.R
\name{csl_styles}
\alias{csl_styles}
\alias{csl_style_exists}
\title{List locally stored styles}
\usage{
csl_styles(style = NULL)

csl_style_exists(style)
}
\arguments{
\item{style}{(character) Style name}
}
\value{
If \code{style=NULL}, a list of length two, independent and dependent
styles. If \code{style} is not \code{NULL}, then a full path to the style file is
returned if the style exists.
}
\description{
List locally stored styles
}
\examples{
# setup
csl_cache$cache_path_set("seasl", type = "tempdir")
csl_cache$cache_path_get()

# List style files
csl_styles()
csl_styles("apa")
csl_styles("zdm")

csl_style_exists("apa")
csl_style_exists("apaggg")

# cleanup
csl_cache$delete_all()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csl_style_find.R
\name{csl_style_find}
\alias{csl_style_find}
\title{Search for a CSL style}
\usage{
csl_style_find(x)
}
\arguments{
\item{x}{(character) a full or partial journal name}
}
\value{
if no matches \code{NULL}. otherwise, one or more file paths
to the style file on your machine
}
\description{
Search for a CSL style
}
\examples{
# setup
csl_cache$cache_path_set("seasl", type = "tempdir")
csl_cache$mkdir()
dir.create(file.path(csl_cache$cache_path_get(), "styles"))
an <- system.file('inst/examples/acta-naturae.csl', package = 'seasl')
file.copy(an, file.path(csl_cache$cache_path_get(), "styles/acta-naturae.csl"))

# find a style
csl_style_find(x = "Naturae")

# cleanup
csl_cache$delete_all()

\dontrun{
# fetch styles
csl_fetch_styles()

# single match
csl_style_find(x = "American Journal of Epidemiology")

# many matches
csl_style_find(x = "American Journal")
csl_style_find(x = "pediatrics")
csl_style_find(x = "analysis and prevention")

# no matches
csl_style_find(x = "foo bar")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seasl-package.R
\docType{package}
\name{seasl-package}
\alias{seasl-package}
\alias{seasl}
\title{R client for CSL styles}
\description{
R client for Citation Style Language (CSL) styles
}
\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
