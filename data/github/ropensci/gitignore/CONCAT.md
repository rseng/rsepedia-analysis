
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gitignore

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/gitignore)](https://cran.r-project.org/package=gitignore)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/gitignore/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/gitignore?branch=main)
[![DOI](https://zenodo.org/badge/184759416.svg)](https://zenodo.org/badge/latestdoi/184759416)
[![rOpenSci
peer-review](https://badges.ropensci.org/303_status.svg)](https://github.com/ropensci/software-review/issues/303)
[![R-CMD-check](https://github.com/PMassicotte/gitignore/workflows/R-CMD-check/badge.svg)](https://github.com/PMassicotte/gitignore/actions)
<!-- badges: end -->

Based on the definition proposed by
[freecodecamp](https://www.freecodecamp.org/news/gitignore-what-is-it-and-how-to-add-to-repo/):

> The .gitignore file is a text file that tells Git which files or
> folders to ignore in a project. A local .gitignore file is usually
> placed in the root directory of a project. You can also create a
> global .gitignore file and any entries in that file will be ignored in
> all of your Git repositories.

For any project, it is therefore important to have a `.gitignore` file
that is complete and accurate. The package `gitignore` provides a simple
R interface to the
[gitignore.io](https://www.toptal.com/developers/gitignore) API. It can
be used to fetch gitignore templates that can be included into the
`.gitignore` file of you git repository. The `gitignore` R package can
be used with R package, R Studio project or with any `.gitignore` file.
Note that by default, the `usethis` package populates the `.gitignore`
for the R language when you create a R project. However, it is common to
use many different programming languages in a project such as `LaTeX`,
`python`, `matlab`, `julia` and so one. This is where the `gitignore`
package shines as it can be used to programmatically modify the
`.gitignore` file of your project.

## Installation

The CRAN version of `gitignore` can be installed using:

``` r
install.packages("gitignore")
```

The dev version of `gitignore` can be installed from GitHub:

``` r
install.packages("devtools")
devtools::install_github("ropensci/gitignore")
```

## Examples

There are currently two useful functions in the package:

-   `gi_available_templates()` to fetch all supported gitignore
    templates.
-   `gi_fetch_templates()` to fetch one or many gitignore templates.

Show the first 25 templates returned by `gi_available_templates()`.

``` r
library(gitignore)

head(gi_available_templates(), 25)
#>  [1] "1c"                   "1c-bitrix"            "a-frame"             
#>  [4] "actionscript"         "ada"                  "adobe"               
#>  [7] "advancedinstaller"    "adventuregamestudio"  "agda"                
#> [10] "al"                   "alteraquartusii"      "altium"              
#> [13] "amplify"              "android"              "androidstudio"       
#> [16] "angular"              "anjuta"               "ansible"             
#> [19] "ansibletower"         "apachecordova"        "apachehadoop"        
#> [22] "appbuilder"           "appceleratortitanium" "appcode"             
#> [25] "appcode+all"
```

Templates can be fetched using the `gi_fetch_templates()` function.

``` r
gi_fetch_templates("R")

# Created by https://www.toptal.com/developers/gitignore/api/r
# Edit at https://www.toptal.com/developers/gitignore?templates=r

### R ###
# History files
.Rhistory
.Rapp.history

# Session Data files
.RData

# User-specific files
.Ruserdata

# Example code in package build process
*-Ex.R

# Output files from R CMD build
/*.tar.gz

# Output files from R CMD check
/*.Rcheck/

# RStudio files
.Rproj.user/

# produced vignettes
vignettes/*.html
vignettes/*.pdf

# OAuth2 token, see https://github.com/hadley/httr/releases/tag/v0.3
.httr-oauth

# knitr and R markdown default cache directories
*_cache/
/cache/

# Temporary files created by R markdown
*.utf8.md
*.knit.md

# R Environment Variables
.Renviron

# pkgdown site
docs/

# translation temp files
po/*~

### R.Bookdown Stack ###
# R package: bookdown caching files
/*_files/

# End of https://www.toptal.com/developers/gitignore/api/r
```

Multiple templates can be fetched by specifying multiple values:

``` r
gi_fetch_templates(c("java", "c++"))

# Created by https://www.toptal.com/developers/gitignore/api/java,c++
# Edit at https://www.toptal.com/developers/gitignore?templates=java,c++

### C++ ###
# Prerequisites
*.d

# Compiled Object files
*.slo
*.lo
*.o
*.obj

# Precompiled Headers
*.gch
*.pch

# Compiled Dynamic libraries
*.so
*.dylib
*.dll

# Fortran module files
*.mod
*.smod

# Compiled Static libraries
*.lai
*.la
*.a
*.lib

# Executables
*.exe
*.out
*.app

### Java ###
# Compiled class file
*.class

# Log file
*.log

# BlueJ files
*.ctxt

# Mobile Tools for Java (J2ME)
.mtj.tmp/

# Package Files #
*.jar
*.war
*.nar
*.ear
*.zip
*.tar.gz
*.rar

# virtual machine crash logs, see http://www.java.com/en/download/help/error_hotspot.xml
hs_err_pid*

# End of https://www.toptal.com/developers/gitignore/api/java,c++
```

By default, templates are copied into the clipboard. It is also possible
to modify a `.gitignore` file using the `gi_write_gitignore()` function.

``` r
f <- file.path(tempdir(), ".gitignore")
new_lines <- gi_fetch_templates("r")
gi_write_gitignore(fetched_template = new_lines, gitignore_file = f)
```

If `gitignore_file` is not specified, `gitignore` will try to find the
`.gitignore` file of your current project or package.

More examples are provided in the vignette.

``` r
browseVignettes("gitignore")
```

You can also visit the [gitignore
website](https://docs.ropensci.org/gitignore/).

## Code of conduct

Please note that the ‘gitignore’ project is released with a [Contributor
Code of
Conduct](https://docs.ropensci.org/gitignore/CODE_OF_CONDUCT.html). By
[contributing to this
project](https://docs.ropensci.org/gitignore/CONTRIBUTING.html), you
agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# gitignore (development version)

# gitignore 0.1.5

* Using Github Actions for continuous integration.

* Fixing CRAN check results where tests failed when internet connection was not available (#18).

# gitignore 0.1.4

* Change backend from https://www.gitignore.io/ to  https://www.toptal.com/developers/gitignore as the former now redirects to the later (#13 @pat-s).

* Use  `file.path()` instead of `paste0()` to build path. @dpprdan 

# gitignore 0.1.3

* This is a minor update that prevent the use of the clipboard on CRAN Linux systems.

# gitignore 0.1.2

* First release on CRAN.

# gitignore 0.1.1

* First release on rOpenSci.
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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
This patch fixes CRAN check results. Tests are now not performed when internet connection is not available.

## Test environments

* Ubuntu oldrel, release, devel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

## Downstream dependencies

There are currently no reverse dependencies for this package.# Contributing to gitignore

This outlines how to propose a change to gitignore. For more detailed
info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

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
*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the gitignore project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See tidyverse [development contributing guide](https://rstd.io/tidy-contrib)
for further details.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# gitignore

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/gitignore)](https://cran.r-project.org/package=gitignore)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Codecov test coverage](https://codecov.io/gh/ropensci/gitignore/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/gitignore?branch=main)
[![DOI](https://zenodo.org/badge/184759416.svg)](https://zenodo.org/badge/latestdoi/184759416)
[![rOpenSci peer-review](https://badges.ropensci.org/303_status.svg)](https://github.com/ropensci/software-review/issues/303)
[![R-CMD-check](https://github.com/PMassicotte/gitignore/workflows/R-CMD-check/badge.svg)](https://github.com/PMassicotte/gitignore/actions)
<!-- badges: end -->

Based on the definition proposed by [freecodecamp](https://www.freecodecamp.org/news/gitignore-what-is-it-and-how-to-add-to-repo/):

> The .gitignore file is a text file that tells Git which files or folders to ignore in a project. A local .gitignore file is usually placed in the root directory of a project. You can also create a global .gitignore file and any entries in that file will be ignored in all of your Git repositories.

For any project, it is therefore important to have a `.gitignore` file that is complete and accurate. The package `gitignore` provides a simple R interface to the [gitignore.io](https://www.toptal.com/developers/gitignore) API. It can be used to fetch gitignore templates that can be included into the `.gitignore` file of you git repository. The `gitignore` R package can be used with R package, R Studio project or with any `.gitignore` file. Note that by default, the `usethis` package populates the `.gitignore` for the R language when you create a R project. However, it is common to use many different programming languages in a project such as `LaTeX`, `python`, `matlab`, `julia` and so one. This is where the `gitignore` package shines as it can be used to programmatically modify the `.gitignore` file of your project.   

## Installation

The CRAN version of `gitignore` can be installed using:

```{r, eval=FALSE}
install.packages("gitignore")
```

The dev version of `gitignore` can be installed from GitHub:

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/gitignore")
```

## Examples

There are currently two useful functions in the package:

- `gi_available_templates()` to fetch all supported gitignore templates.
- `gi_fetch_templates()` to fetch one or many gitignore templates.

Show the first 25 templates returned by `gi_available_templates()`.

```{r}
library(gitignore)

head(gi_available_templates(), 25)
```

Templates can be fetched using the `gi_fetch_templates()` function.

```{r, results='markup', comment=""}
gi_fetch_templates("R")
```

Multiple templates can be fetched by specifying multiple values:

```{r, results='markup', comment=""}
gi_fetch_templates(c("java", "c++"))
```

By default, templates are copied into the clipboard. It is also possible to modify a `.gitignore` file using the `gi_write_gitignore()` function.

```{r, eval=FALSE}
f <- file.path(tempdir(), ".gitignore")
new_lines <- gi_fetch_templates("r")
gi_write_gitignore(fetched_template = new_lines, gitignore_file = f)
```

If `gitignore_file` is not specified, `gitignore` will try to find the `.gitignore` file of your current project or package.

More examples are provided in the vignette.

```{r, eval=FALSE}
browseVignettes("gitignore")
```

You can also visit the [gitignore website](https://docs.ropensci.org/gitignore/).

## Code of conduct

Please note that the 'gitignore' project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/gitignore/CODE_OF_CONDUCT.html). By [contributing to this project](https://docs.ropensci.org/gitignore/CONTRIBUTING.html), you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Basic idea

Based on the definition proposed by [freecodecamp](https://www.freecodecamp.org/news/gitignore-what-is-it-and-how-to-add-to-repo/):

> The .gitignore file is a text file that tells Git which files or folders to ignore in a project. A local .gitignore file is usually placed in the root directory of a project. You can also create a global .gitignore file and any entries in that file will be ignored in all of your Git repositories.

For any project, it is therefore important to have a `.gitignore` file that is complete and accurate. The package `gitignore` provides a simple R interface to the [gitignore.io](https://www.toptal.com/developers/gitignore) API. It can be used to fetch gitignore templates that can be included into the `.gitignore` file of you git repository. The `gitignore` R package can be used with R package, R Studio project or with any `.gitignore` file. Note that by default, the `usethis` package populates the `.gitignore` for the R language when you create a R project. However, it is common to use many different programming languages in a project such as `LaTeX`, `python`, `matlab`, `julia` and so on. This is where the `gitignore` package shines as it can be used to programmatically modify the `.gitignore` file of your project.

`gitignore` is a simple R package that provide an interface to query [gitignore.io](https://www.toptal.com/developers/gitignore) to fetch gitignore templates that can be included in the .gitignore file. More than 450 templates are currently available. There are actually two functions in the package:

- `gi_available_templates()`: to get a list of all templates available on [gitignore.io](https://www.toptal.com/developers/gitignore).
- `gi_fetch_templates()`: to get one or more template(s). This function can copy the returned template(s) in the clipboard as well as happening the the `.gitignore` file in your project directory.

## Examples

```{r setup}
library(gitignore)
```

Get the list of all available templates from [gitignore.io](https://www.toptal.com/developers/gitignore):

```{r}
head(gi_available_templates(), 50)
```

The function `gi_fetch_templates()` can be used to fetch particular template(s). For example, one could want to get the `java` and `c++` templates as follow:

```{r}
gi_fetch_templates(c("java", "c++"))
```

By default, the template(s) are not copied into the clipboard, this can be changed by setting `copy_to_clipboard = TRUE`:

```{r, eval=FALSE}
gi_fetch_templates(c("java", "c++"), copy_to_clipboard = TRUE)
```

Additionally, you can tell `gi_fetch_templates()` to append automatically the `.gitignore` file in your project by setting `append_gitignore = TRUE`:

```{r, eval=FALSE}
gi_fetch_templates(c("R"), append_gitignore = TRUE)
```

It is also possible to specify the `.gitignore` file to be modified by specifying the `gitignore_file` argument.

```{r, message=TRUE}
f <- file.path(tempdir(), ".gitignore")

file.create(f)

gi_fetch_templates("R", gitignore_file = f, append_gitignore = TRUE)
  
readLines(f)
  
```

---
title: "Introduction"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Basic idea

Based on the definition proposed by [freecodecamp](https://www.freecodecamp.org/news/gitignore-what-is-it-and-how-to-add-to-repo/):

> The .gitignore file is a text file that tells Git which files or folders to ignore in a project. A local .gitignore file is usually placed in the root directory of a project. You can also create a global .gitignore file and any entries in that file will be ignored in all of your Git repositories.

For any project, it is therefore important to have a `.gitignore` file that is complete and accurate. The package `gitignore` provides a simple R interface to the [gitignore.io](https://www.toptal.com/developers/gitignore) API. It can be used to fetch gitignore templates that can be included into the `.gitignore` file of you git repository. The `gitignore` R package can be used with R package, R Studio project or with any `.gitignore` file. Note that by default, the `usethis` package populates the `.gitignore` for the R language when you create a R project. However, it is common to use many different programming languages in a project such as `LaTeX`, `python`, `matlab`, `julia` and so on. This is where the `gitignore` package shines as it can be used to programmatically modify the `.gitignore` file of your project.

`gitignore` is a simple R package that provide an interface to query [gitignore.io](https://www.toptal.com/developers/gitignore) to fetch gitignore templates that can be included in the .gitignore file. More than 450 templates are currently available. There are actually two functions in the package:

- `gi_available_templates()`: to get a list of all templates available on [gitignore.io](https://www.toptal.com/developers/gitignore).
- `gi_fetch_templates()`: to get one or more template(s). This function can copy the returned template(s) in the clipboard as well as happening the the `.gitignore` file in your project directory.

## Examples

```{r setup}
library(gitignore)
```

Get the list of all available templates from [gitignore.io](https://www.toptal.com/developers/gitignore):

```{r}
head(gi_available_templates(), 50)
```

The function `gi_fetch_templates()` can be used to fetch particular template(s). For example, one could want to get the `java` and `c++` templates as follow:

```{r}
gi_fetch_templates(c("java", "c++"))
```

By default, the template(s) are not copied into the clipboard, this can be changed by setting `copy_to_clipboard = TRUE`:

```{r, eval=FALSE}
gi_fetch_templates(c("java", "c++"), copy_to_clipboard = TRUE)
```

Additionally, you can tell `gi_fetch_templates()` to append automatically the `.gitignore` file in your project by setting `append_gitignore = TRUE`:

```{r, eval=FALSE}
gi_fetch_templates(c("R"), append_gitignore = TRUE)
```

It is also possible to specify the `.gitignore` file to be modified by specifying the `gitignore_file` argument.

```{r, message=TRUE}
f <- file.path(tempdir(), ".gitignore")

file.create(f)

gi_fetch_templates("R", gitignore_file = f, append_gitignore = TRUE)
  
readLines(f)
  
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gitignore-package.R
\docType{package}
\name{gitignore-package}
\alias{gitignore}
\alias{gitignore-package}
\title{gitignore: Create Useful .gitignore Files for your Project}
\description{
Simple interface to query gitignore.io to fetch gitignore templates that can be included in the .gitignore file. More than 450 templates are currently available.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/gitignore/}
  \item \url{https://github.com/ropensci/gitignore}
  \item Report bugs at \url{https://github.com/ropensci/gitignore/issues}
}

}
\author{
\strong{Maintainer}: Philippe Massicotte \email{pmassicotte@hotmail.com} (\href{https://orcid.org/0000-0002-5919-4116}{ORCID})

Other contributors:
\itemize{
  \item Amanda Dobbyn \email{amanda.e.dobbyn@gmail.com} [reviewer]
  \item Mauro Lepore \email{maurolepore@gmail.com} (\href{https://orcid.org/0000-0002-1986-7988}{ORCID}) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gi_available_templates.R
\name{gi_available_templates}
\alias{gi_available_templates}
\title{Fetch available templates from gitignore.io}
\usage{
gi_available_templates()
}
\value{
A character with all templates supported by gitignore.io.
}
\description{
This return list of all templates supported by gitignore.io.
}
\details{
The returned list is returned as lower case characters.
}
\examples{
gi_available_templates()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gi_fetch_templates.R
\name{gi_fetch_templates}
\alias{gi_fetch_templates}
\title{Fetch gitignore template(s) from gitignore.io}
\usage{
gi_fetch_templates(
  template_name,
  copy_to_clipboard = FALSE,
  append_gitignore = FALSE,
  gitignore_file = here::here(".gitignore")
)
}
\arguments{
\item{template_name}{A character vector with values included in
\code{\link{gi_available_templates}}.}

\item{copy_to_clipboard}{Logical. Should the returned template(s) be copied
to the clipboard? Otherwise, it will be printed in the console. Default is FALSE.}

\item{append_gitignore}{Logical. Should the .gitignore be modified to include
the returned template(s)?}

\item{gitignore_file}{The path of the .gitignore file to be modified. By
default, it will try to find it in the current package/project using
`here::here(".gitignore")`.}
}
\value{
A character containing gitignore template(s).
}
\description{
Fetch gitignore template(s) from gitignore.io
}
\examples{
# Fetch template for the R language
gi_fetch_templates("R")

# You can combine many templates at once
gi_fetch_templates(c("R", "python", "java"))

# The .gitignore file can be automatically modified with `append_gitignore = TRUE`
gi_fetch_templates(c("R", "python", "java"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gi_write_gitignore.R
\name{gi_write_gitignore}
\alias{gi_write_gitignore}
\title{Append or create a .gitignore file}
\usage{
gi_write_gitignore(fetched_template, gitignore_file = here::here(".gitignore"))
}
\arguments{
\item{fetched_template}{Template(s) returned by `gi_fetch_templates()`.}

\item{gitignore_file}{Path of the .gitignore file to modify.}
}
\value{
TRUE if succeeds to write/append the .gitignore, FALSE otherwise.
}
\description{
Use the returned template(s) to append the existing .gitignore file.
}
\examples{
\dontrun{
f <- file.path(tempdir(), ".gitignore")
new_lines <- gi_fetch_templates("r")
gi_write_gitignore(new_lines, f)

unlink(f)
}
}
