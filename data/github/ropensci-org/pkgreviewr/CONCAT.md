
<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) ![GitHub R package version](https://img.shields.io/github/r-package/v/ropensci-org/pkgreviewr) 
  [![R-CMD-check](https://github.com/ropensci-org/pkgreviewr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-org/pkgreviewr/actions)
[![codecov](https://codecov.io/gh/ropensci-org/pkgreviewr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci-org/pkgreviewr)
 <!-- badges: end -->
 
# pkgreviewr

The goal of pkgreviewr is to facilitate **rOpenSci** reviewers in their package reviews. 

It creates a review project containing populated templates of all the files you'll need to complete your review. It also clones the source code of the package under review to a convenient location, allowing easy checking and testing.

See [Getting started](articles/get_started.html) vignette for more details

## Installation

You can install `pkgreviewr` from GitHub with:


``` r
# install.packages("remotes")
remotes::install_github("ropensci-org/pkgreviewr")
```
<br>

### Git

`pkgreviewr` functions clone the source code of the package under review so require **Git** to be installed. 

### GitHub user configuration

Because rOpenSci reviews are conducted through github repository [`ropensci/software-review`](https://github.com/ropensci/software-review), **`pkgreviewr` uses your GitHub username to prepopulate various fields in the review project files**.

To detect your username correctly, a PAT, Personal Authorisation Token, needs to be set up.
You can use `usethis::browse_github_pat()` to generate a PAT and `usethis::edit_r_environ()` to store it as environment variable `GITHUB_PAT` or `GITHUB_TOKEN` in your .`Renviron` file. For more info, see article on publishing review on GitHub in pkgreviewr documentation.

If you do not have a PAT set up, you will receive an warning and any fields related to your GitHub username will not be correctly populated. However, this shouldn't affect your ability to complete your review. 


### R Notebooks

The package currently also makes use of [**`R Notebooks`**](https://rmarkdown.rstudio.com/r_notebooks.html) (an RMarkdown format) and requires installation of **Rstudio version 1.0** or higher, but we are [considering offering an option to remove the requirement for RStudio](https://github.com/ropenscilabs/pkgreviewr/issues/64).

***

## Review workflow

<br>

#### 1. Create and initialise review project 

```r
library(pkgreviewr)
pkgreview_create(pkg_repo = "ropensci/rdflib", 
                 review_parent = "~/Documents/workflows/rOpenSci/reviews/")
```

The review project directory will contain all the files you'll need to complete the review and will be initialised with git.

```
rdflib-review
├── README.md
├── index.Rmd
├── pkgreview.md
└── rdflib-review.Rproj
```
<br>

#### 2. Perform review

Open `index.Rmd` and work through the review in the notebook. You can make notes either in `index.Rmd` or directly in the `pkgreview.md` response template.

<br>

#### 3. Submit review

Submit your review in the package [`ropensci/software-review`](https://github.com/ropensci/software-review/issues) issue by copying and pasting the completed `pkgreview.md` template.

<br>

#### 4. Publish review*

OPTIONAL: Publish your review on GitHub. See [vignette](articles/publish-review-on-github.html) for further instructions

<br>


***

## `pkgreviewr` for editors 

`pkgreviewr` can now also be used to set up projects for editor checks. See [`pkgreviewr` for editors](articles/editors.html) vignette.


***

Please note that 'pkgreviewr' is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.
# pkgreviewr 0.2.0

* Adjust to new location of templates (thanks @maelle for PR)
* More robust handling of 404 messages
* Use `gh::gh_token` instead of deprecated `usethis::github_token`
* Make testing of created folder contents more robust.

# pkgreviewr 0.1.2

* More robust project handling (thanks @rorynolan for PR)
* Streamlined initialisation messages
* More robust detection of gh username (thanks @benmarwick & @llrs for suggestion) and repository url

# pkgreviewr 0.1.1

* Added ability to specify issue_no


# pkgreviewr 0.1.0

* Added editor templates to enable editor checks and other duties
* git2r updates resolved authentication issues so temporary workaround removed.
* Functionality streamlined more stable

# pkgreviewr 0.0.4

* Added temporary workaround for git2r authentication issues, created by GitHub security protocol change (esp on macOSX).
* `pkgreview_init()` re-instated allowing for the configuration of the review separately form the creation stage.
* More isolated workflow and improved error and exception handling



# pkgreviewr 0.0.3

* Added a `NEWS.md` file to track changes to the package.
* Added package level documentation
* Combined `pkgreview_create()` and `pkgreview_init()`


# pkgreviewr 0.0.2

* Addressed initial installation bugs.
* Added checks for Rstudio and GitHub credentials


# pkgreviewr 0.0.1

* Initial prototype

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
# Contributing to pkgreviewr

This outlines how to propose a change to pkgreviewr. For more detailed
info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, consider opening an issue and
discussing your proposed changes. If you’ve found a
bug, please create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR). 
*  Please **branch off of and make pull requests back to branch `devel`.** 
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

Please note that the pkgreviewr project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

## Package Review

*Please check off boxes as applicable, and elaborate in comments below.  Your review is not limited to these topics, as described in the reviewer guide*

- [ ] As the reviewer I confirm that there are no conflicts of interest for me to review this work (If you are unsure whether you are in conflict, please speak to your editor _before_ starting your review).

#### Documentation

The package includes all the following forms of documentation:

- [ ] **A statement of need** clearly stating problems the software is designed to solve and its target audience in README
- [ ] **Installation instructions:** for the development version of package and any non-standard dependencies in README
- [ ] **Vignette(s)** demonstrating major functionality that runs successfully locally
- [ ] **Function Documentation:** for all exported functions in R help
- [ ] **Examples** for all exported functions in R Help that run successfully locally
- [ ] **Community guidelines** including contribution guidelines in the README or CONTRIBUTING, and DESCRIPTION with `URL`, `BugReports` and `Maintainer` (which may be autogenerated via `Authors@R`).

>##### For packages co-submitting to JOSS
>
>- [ ] The package has an **obvious research application** according to [JOSS's definition](http://joss.theoj.org/about#submission_requirements)
>
>The package contains a `paper.md` matching [JOSS's requirements](http://joss.theoj.org/about#paper_structure) with:
>
>- [ ] **A short summary** describing the high-level functionality of the software
>- [ ] **Authors:**  A list of authors with their affiliations
>- [ ] **A statement of need** clearly stating problems the software is designed to solve and its target audience.
>- [ ] **References:** with DOIs for all those that have one (e.g. papers, datasets, software).

#### Functionality

- [ ] **Installation:** Installation succeeds as documented.
- [ ] **Functionality:** Any functional claims of the software been confirmed.
- [ ] **Performance:** Any performance claims of the software been confirmed.
- [ ] **Automated tests:** Unit tests cover essential functions of the package
   and a reasonable range of inputs and conditions. All tests pass on the local machine.
- [ ] **Packaging guidelines**: The package conforms to the rOpenSci packaging guidelines

#### Final approval (post-review)

- [ ] **The author has responded to my review and made changes to my satisfaction. I recommend approving this package.**

Estimated hours spent reviewing: 

---

### Review Comments
## Package Review

*Please check off boxes as applicable, and elaborate in comments below.  Your review is not limited to these topics, as described in the reviewer guide*

- [ ] As the reviewer I confirm that there are no conflicts of interest for me to review this work (If you are unsure whether you are in conflict, please speak to your editor _before_ starting your review).

#### Documentation

The package includes all the following forms of documentation:

- [ ] **A statement of need** clearly stating problems the software is designed to solve and its target audience in README
- [ ] **Installation instructions:** for the development version of package and any non-standard dependencies in README
- [ ] **Vignette(s)** demonstrating major functionality that runs successfully locally
- [ ] **Function Documentation:** for all exported functions in R help
- [ ] **Examples** for all exported functions in R Help that run successfully locally
- [ ] **Community guidelines** including contribution guidelines in the README or CONTRIBUTING, and DESCRIPTION with `URL`, `BugReports` and `Maintainer` (which may be autogenerated via `Authors@R`).

>##### For packages co-submitting to JOSS
>
>- [ ] The package has an **obvious research application** according to [JOSS's definition](http://joss.theoj.org/about#submission_requirements)
>
>The package contains a `paper.md` matching [JOSS's requirements](http://joss.theoj.org/about#paper_structure) with:
>
>- [ ] **A short summary** describing the high-level functionality of the software
>- [ ] **Authors:**  A list of authors with their affiliations
>- [ ] **A statement of need** clearly stating problems the software is designed to solve and its target audience.
>- [ ] **References:** with DOIs for all those that have one (e.g. papers, datasets, software).

#### Functionality

- [ ] **Installation:** Installation succeeds as documented.
- [ ] **Functionality:** Any functional claims of the software been confirmed.
- [ ] **Performance:** Any performance claims of the software been confirmed.
- [ ] **Automated tests:** Unit tests cover essential functions of the package
   and a reasonable range of inputs and conditions. All tests pass on the local machine.
- [ ] **Packaging guidelines**: The package conforms to the rOpenSci packaging guidelines

#### Final approval (post-review)

- [ ] **The author has responded to my review and made changes to my satisfaction. I recommend approving this package.**

Estimated hours spent reviewing: 

---

### Review Comments
---
output: github_document
---


# `rdflib` - package review repository

##

This repo contains files associated with the **rOpenSci** review of

### **`rdflib`: ropensci/onboarding**  issue [\#169](https://github.com/ropensci/onboarding/issues/169).

<br>


***

## **Reviewer:** [\@annakrystalli](https://github.com/annakrystalli)
### Review Submitted:
**`r cat(sprintf("**Last updated:** %s", Sys.Date()))`**

<br>

### see the review report [here:](https://annakrystalli.github.io/rdflib-review/index.nb.html)

or view the submiited review to rOpenSci [here:](https://github.com/annakrystalli/rdflib-review/blob/master/pkgreview.md)

<br>


## Package info

**Description:**

The Resource Description Framework, or 'RDF' is a widely used
             data representation model that forms the cornerstone of the 
             Semantic Web. 'RDF' represents data as a graph rather than 
             the familiar data table or rectangle of relational databases.
             The 'rdflib' package provides a friendly and concise user interface
             for performing common tasks on 'RDF' data, such as reading, writing
             and converting between the various serializations of 'RDF' data,
             including 'rdfxml', 'turtle', 'nquads', 'ntriples', and 'json-ld';
             creating new 'RDF' graphs, and performing graph queries using 'SPARQL'.
             This package wraps the low level 'redland' R package which
             provides direct bindings to the 'redland' C library.  Additionally,
             the package supports the newer and more developer friendly
             'JSON-LD' format through the 'jsonld' package. The package
             interface takes inspiration from the Python 'rdflib' library.

**Author:** `r person("Carl", "Boettiger", 
                  email = "cboettig@gmail.com", 
                  role = c("aut", "cre", "cph"),
                  comment=c(ORCID = "0000-0002-1642-628X"))`

**repo url:** <https://github.com/cboettig/rdflib>

**website url:** <https://cboettig.github.io/rdflib/>
---
output:
    html_notebook:
        toc: true
        toc_float: true
editor_options:
  chunk_output_type: inline
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)

library(magrittr)
```

# `rdflib` - package review

## **Reviewer:** [\@annakrystalli](https://github.com/annakrystalli)

### Review Submitted:
**`r cat(sprintf("**Last updated:** %s", Sys.Date()))`**

***

<br>

This report contains documents the review of **rOpenSci** submitted package:

### **`rdflib`: ropensci/onboarding**  issue [\#169](https://github.com/ropensci/onboarding/issues/169).

<br>

## Package info

**Description:**

The Resource Description Framework, or 'RDF' is a widely used
             data representation model that forms the cornerstone of the 
             Semantic Web. 'RDF' represents data as a graph rather than 
             the familiar data table or rectangle of relational databases.
             The 'rdflib' package provides a friendly and concise user interface
             for performing common tasks on 'RDF' data, such as reading, writing
             and converting between the various serializations of 'RDF' data,
             including 'rdfxml', 'turtle', 'nquads', 'ntriples', and 'json-ld';
             creating new 'RDF' graphs, and performing graph queries using 'SPARQL'.
             This package wraps the low level 'redland' R package which
             provides direct bindings to the 'redland' C library.  Additionally,
             the package supports the newer and more developer friendly
             'JSON-LD' format through the 'jsonld' package. The package
             interface takes inspiration from the Python 'rdflib' library.

**Author:** `r person("Carl", "Boettiger", 
                  email = "cboettig@gmail.com", 
                  role = c("aut", "cre", "cph"),
                  comment=c(ORCID = "0000-0002-1642-628X"))`

**repo url:** <https://github.com/cboettig/rdflib>

**website url:** <https://cboettig.github.io/rdflib/>

## Review info


#### See [reviewer guidelines](https://github.com/ropensci/onboarding/blob/master/reviewing_guide.md) for further information on the rOpenSci review process.

**key review checks:**

- Does the code comply with **general principles in the [Mozilla reviewing guide](https://mozillascience.github.io/codeReview/review.html)**?
- Does the package **comply with the [ROpenSci packaging guide](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md)**?
- Are there **improvements** that could be made to the **code style?**
- Is there **code duplication** in the package that should be reduced?
- Are there **user interface improvements** that could be made?
- Are there **performance improvements** that could be made?
- Is the **documentation** (installation instructions/vignettes/examples/demos) **clear and sufficient**?

Please be respectful and kind to the authors in your reviews. The rOpenSci [code of conduct](https://github.com/ropensci/onboarding/blob/master/policies.md/#code-of-conduct) is mandatory for everyone involved in our review process.

***

### session info


```{r sessionInfo}
sessionInfo()
```


```{r pkg_dir, echo = F}
pkg_dir <- "/private/var/folders/8p/87cqdx2s34vfvcgh04l6z72w0000gn/T/Rtmp4i9nQx/rdflib-review/../rdflib"
```

## Test installation

### test local `rdflib` install:

```{r test-local}
install(pkg_dir, dependencies = T, build_vignettes = T)
```

```{r github-rm}
remove.packages("rdflib")
```
#### **comments:**

<!-- record comments on local install here -->

***

### test install of `rdflib` from GitHub with:

```{r test-github}
devtools::install_github("cboettig/rdflib", dependencies = T, build_vignettes = T)
```

#### **comments:**

<!-- record comments on github install here -->

***



## Check package integrity

### run checks on `rdflib` source:

```{r check-checks}
devtools::check(pkg_dir)
```
#### **comments:**

<!-- record comments on checks here -->

***

### run tests on `rdflib` source:

```{r check-tests}
devtools::test(pkg_dir)
```
#### **comments:**

<!-- record comments on tests here -->

***


### check `rdflib` for goodpractice:

```{r test-goodpractice}
goodpractice::gp(pkg_dir)
```
#### **comments:**

<!-- record comments on goodpractice here -->

***

## Check package metadata files

### inspect

- #### [README](https://github.com/cboettig/rdflib)
- #### [DESCRIPTION](https://github.com/cboettig/rdflib/blob/master/DESCRIPTION)
- #### [NAMESPACE](https://github.com/cboettig/rdflib/blob/master/NAMESPACE)

### spell check

```{r spell-check}
devtools::spell_check(pkg_dir)
```

#### **comments:**

<!-- record comments on metadata files here -->

***

## Check documentation

online documentation: **<https://cboettig.github.io/rdflib/>**

* Is the documentation (installation instructions/vignettes/examples/demos) clear and sufficient?

### test `rdflib` function help files:

```{r test-help}
help(package = "rdflib")
```

#### **comments:**

<!-- record comments on help files here -->

***

### test `rdflib` vignettes:

```{r test-vignettes}
vignette(package = "rdflib")
```

#### **comments:**

<!-- record comments on vignettes here -->

***

## Test functionality:

- Are there **user interface improvements** that could be made?
- Are there **performance improvements** that could be made?

```{r free-style}
library("rdflib")
```

```{r parse-functions}
exports <-ls("package:rdflib")
exports
```

<!-- experiment with package functions -->

```{r exp-chunk}


```

#### **comments:**

<!-- record comments on rdflib experimentation here -->

***

## Inspect code:

- Does the package **comply with the [ROpenSci packaging guide](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md)**?
    * good [function & variable naming?](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md#funvar)
    * good [dependency management](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md#deps)?
- Are there **improvements** that could be made to the **code style?**
- Is there **code duplication** in the package that should be reduced?

```{r inspect-code}
pkgreviewr::pkgreview_print_source("rdflib")
```
**\* might not be suitable for large packages with many exported functions**


<br>
<br>

#### **comments:**

<!-- record comments on rdflib source code here -->


## Review test suite:

### test coverage

```{r pkg_coverage}
covr::package_coverage(pkg_dir)

```

### inspect [tests](https://github.com/cboettig/rdflib/blob/master/tests/testthat)


#### **comments:**

<!-- record comments on testing suite here -->


***
---
title: "Getting Started"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting Started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

### Example

This is a basic example of **setting up an rOpenSci package review project**:


### 1. create review project

Create the review project, using `pkgreview_create`. The function takes arguments:

* **`pkg_repo`:** the **GitHub repo** details of the **package under review** in the form `username/repo` 
* **`review_parent`:**, the **local directory** in which the **review project (and folder) will be created** and **package source code will be cloned into**.

``` r
library(pkgreviewr)
pkgreview_create(pkg_repo = "ropensci/rdflib", 
                 review_parent = "~/Documents/workflows/rOpenSci/reviews/")
```

The function creates a new review project in the `review_parent` directory following project naming convention `{pkgname}-review` and populates the review templates to create all required documents. 

#### review files

The review project directory will contain all the files you'll need to complete the review and will be initialised with git.

```
rdflib-review
├── README.md
├── index.Rmd
├── pkgreview.md
└── rdflib-review.Rproj
```

#### **`index.Rmd`**

The most important file it creates is the `index.Rmd html_notebook` file. This workbook is prepopulated with all the major steps required to complete the review in an interactive document to perform and record it in. It also extracts useful links, information and parameter values. 

**See example [here](https://github.com/annakrystalli/pkgreviewr/blob/master/inst/examples/example-review-index.Rmd).**

Once rendered to **`index.nb.html`** (`*.nb.html` is the notebook file format), this report can be pushed to GitHub for publication.

#### **`pkgreview.md`** 

Template response form to submit to the package rOpenSci onboarding review issue. 

**See template [here](https://github.com/annakrystalli/pkgreviewr/blob/master/inst/examples/example-pkgreview.md)**.

#### **`README.md`** 

Prepopulated README for the review repo that will present the repo to people navigating to it. 

**See example [here:](https://github.com/annakrystalli/pkgreviewr/blob/master/inst/examples/example-README.md)**.

***





#### clone of package source code

To enable local testing of the package, review creation also clones the review package source code into `review_parent` from the github repository defned in `pkg_repo` . This also makes it available for local review and perhaps even a pull request. Correcting typos in documentation can be a great review contribution, but first you might want to check the contributing guidelines or ask the author if they are open to such pull requests.

The resulting files from a successful review project will look like this: 

```
reviews
├── rdflib
│   ├── DESCRIPTION
│   ├── LICENSE
│   ├── NAMESPACE
│   ├── NEWS.md
│   ├── R
│   │   └── rdf.R
│   ├── README.Rmd
│   ├── README.md
│   ├── appveyor.yml
│   ├── codecov.yml
│   ├── codemeta.json
│   ├── docs
│   │   ├── LICENSE.html
│   │   ├── articles
│   │   │   ├── index.html
│   │   │   ├── rdflib.html
│   │   │   └── rdflib_files
│   │   │       ├── datatables-binding-0.2
│   │   │       │   └── datatables.js
│   │   │       ├── dt-core-1.10.12
│   │   │       │   ├── css
│   │   │       │   │   ├── jquery.dataTables.extra.css
│   │   │       │   │   └── jquery.dataTables.min.css
│   │   │       │   └── js
│   │   │       │       └── jquery.dataTables.min.js
│   │   │       ├── htmlwidgets-0.9
│   │   │       │   └── htmlwidgets.js
│   │   │       └── jquery-1.12.4
│   │   │           ├── LICENSE.txt
│   │   │           └── jquery.min.js
│   │   ├── authors.html
│   │   ├── index.html
│   │   ├── jquery.sticky-kit.min.js
│   │   ├── link.svg
│   │   ├── news
│   │   │   └── index.html
│   │   ├── pkgdown.css
│   │   ├── pkgdown.js
│   │   └── reference
│   │       ├── index.html
│   │       ├── rdf.html
│   │       ├── rdf_add.html
│   │       ├── rdf_parse.html
│   │       ├── rdf_query.html
│   │       ├── rdf_serialize.html
│   │       └── rdflib-package.html
│   ├── inst
│   │   ├── examples
│   │   │   └── rdf_table.R
│   │   └── extdata
│   │       ├── ex.xml
│   │       └── vita.json
│   ├── man
│   │   ├── rdf.Rd
│   │   ├── rdf_add.Rd
│   │   ├── rdf_parse.Rd
│   │   ├── rdf_query.Rd
│   │   ├── rdf_serialize.Rd
│   │   └── rdflib-package.Rd
│   ├── paper.bib
│   ├── paper.md
│   ├── rdflib.Rproj
│   ├── tests
│   │   ├── testthat
│   │   │   └── test-rdf.R
│   │   └── testthat.R
│   └── vignettes
│       └── rdflib.Rmd
└── rdflib-review
    ├── README.md
    ├── index.Rmd
    └── rdflib-review.Rproj

```

<br>

### 2. Perform your review:

Use the index.Rmd notebook to work through the review interactively. The document is designed to guide the process in a logical fashion and bring your attention to relevant aspects and information at different stages of the review. You can make notes and record comments within index.Rmd or directly in the review submission form. 

<br>

### 3. Submit your review:

Currently the workflow is just set up for you to just copy your response from your completed `pkgreview.md` and paste it into the review issue but we're exploring programmatic submission also. Because the response is currently submitted as `.md`, package `reprex` might be useful for inserting reproducible demos of any issues encountered. 

<br>

### 4. Publish your report by pushing to GitHub *

Optional. Have a look at the **Publish pkgreview on GitHub** vignette.

---
title: "pkgreviewr for editors"
author: "Anna Krystalli"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{pkgreviewr for editors}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


`pkgreviewr` now also provides templates for editors. 

```{r, eval = F}
library(pkgreviewr)

pkgreview_create("nldoc/nlrx", 
                 review_parent = "~/Documents/workflows/rOpenSci/editorials/",
                 template = "editor",
                 issue_no = 262)
```

As with initialising a review, the repo is downloaded into `review_parent`. 

It also creates a `pkgreviewr` project containing editor related templates: 

```
.
├── editor.md
├── index.Rmd
├── nlrx-review.Rproj
└── request.Rmd
```

### `index.Rmd`

Similar to a review project, `index.Rmd` is a pre-populated notebook where initial editor checks can be performed 


### `editor.md`

This file contains the editor checks response template

### `request.Rmd`

This contains a prepopulated parametarised `.Rmd` for generating request emails. The parameters it accepts are:

- person's name 
- some friendly banter
- whether JOSS submission is involved

Use:

```{r, eval=FALSE}
render_request()
```

to render the document by first launching a form in the browser to enter parameter values. Alternatively, you can manually edit the YAML header with parameter values and knit.

---
title: "Publish pkgreview on GitHub"
author: "Anna Krystalli"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Publish pkgreview on GitHub}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## 1. Create GitHub repo and add it as a remote to your review project

Currently we've not integrated the creation of an remote GitHub repo but we're exploring this functionality also. In the mean time, you **can create a repo to publish your review in two ways:**

### Manually through GitHub

1. Head over to to your **repositories tab** on [**Github**](https://github.com).

1. Click on **New** to create a **blank** repository.
    - Follow naming convention `"{pkgname}-review"`.
    - Make sure you **don't automatically create any files on GitHub**. These will cause merge conflicts when you try to push to the remote for the first time.
    ![](assets/manual_gh.png)

1. Click on **Clone or download** and **copy the link** displayed

![](assets/clone.png)

1. **Open a terminal** (in Rstudio go to _Tools > Terminal > New Terminal_). In the terminal, ensure you are in the review project.

1. **Add your github repo as a remote** by running the following code in the terminal, **substituting in the link you copied** from GitHub.

    ```
    git remote add origin <the-link-you-copied-from-github>
    git push -u origin master
    ```
  For example, in my case I add my repo to as a remote like so:

    ```
    git remote add origin https://github.com/annakrystalli/rdflib-review.git
    git push -u origin master
    ```
  Follow any authentication steps required

****

### **Programmatically using `usethis::use_github()`**

To use `usethis::use_github`, you'll need to supply a **github personal access token (PAT) token**. The easiest way to set it up for all your r workflows is to store you PAT in a `GITHUB_PAT` system variable in your [.Renviron](https://csgillespie.github.io/efficientR/3-3-r-startup.html#renviron) dotfile.  To do this:

1. **Generate PAT**: use `usethis::browse_github_pat` to launch page to generate a PAT on GitHub.

2. **[Add PAT to your `.Renviron`](https://github.com/jennybc/happy-git-with-r/blob/master/81_github-api-tokens.Rmd) dot file**: Use `usethis::edit_r_environ()` to open your user level `.Renviron`, paste the copied PAT token from github and save it like so:

<img src="assets/renviron.png" width=800>

3. **Create Github repo & add as remote:** Now, while in your review project in Rstudio, run: 

    ```r
    usethis::use_github(protocol = "https")
    ```
  to create a github repository for your review and add it as a remote for your review project. The naming of the github repository is handled automatically.

<div class="alert alert-warning">
  <strong>Warning!</strong> Because of ongoing big changes in dependency `git2r`, you may encounter authentication problems with `usethis::use_github()`. Refer to this <a href="https://community.rstudio.com/t/difficulty-using-usethis-use-github/5579/4"> discussion thread for further details.</a>
</div>

<br>

***

## 2. Commit the review files and push them to github

+ In the **`git` panel** in Rstudio, **select the files you want to share on github**. You can chose to only share `index.nb.html`, the rendered report or include the `index.Rmd`. Also select `README.md` so your repository has an appropriate README from which your review can be accessed.

+ **Commit** the files adding an appropriate commit message

+ **Push** your changes to GitHub

<br>

***

## 3. Enable GitHub Pages

+ In your review GitHub repository click on **Settings**

+ Scroll down to the **GitHub Pages** section and change **Source** location to **master branch**

  <img src="assets/select_gh.png" width=500>

+ **Github Pages is now enabled and your report review [will be published](http://annakrystalli.me/rdflib-review/index.nb.html) at the link displayed:**

  <img src="assets/enabled_gh.png" width=500>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render-templates.R
\name{render_request}
\alias{render_request}
\title{Render request email body}
\usage{
render_request()
}
\value{
renders \code{request.Rmd} using parameters provided.
}
\description{
Launches an interactive input browser tab to complete required parameters:
\itemize{
\item \code{reviewer_first_name}: reviewers first name
\item \code{banter}: character string of custom greeting message
\item \code{JOSS}: logical, whether review includes submission to JOSS
}
}
\examples{
\dontrun{
render_request()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgreviewr-package.R
\docType{package}
\name{pkgreviewr-package}
\alias{pkgreviewr}
\alias{pkgreviewr-package}
\title{pkgreviewr: rOpenSci package review project template}
\description{
Creates files and collects materials necessary to complete an rOpenSci package review. 
    Review files are prepopulated with review package specific metadata. Review package source code is
    also cloned for local testing and inspection.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci-org/pkgreviewr}
  \item Report bugs at \url{https://github.com/ropensci-org/pkgreviewr/issues}
}

}
\author{
\strong{Maintainer}: Anna Krystalli \email{annakrystalli@googlemail.com}

Authors:
\itemize{
  \item Maëlle Salmon \email{maelle.salmon@yahoo.se}
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render-templates.R
\name{pkgreview_index_rmd}
\alias{pkgreview_index_rmd}
\alias{pkgreview_readme_md}
\alias{pkgreview_request}
\title{Create review templates.}
\usage{
pkgreview_index_rmd(pkg_data, template = c("review", "editor"))

pkgreview_readme_md(pkg_data)

pkgreview_request(pkg_data)
}
\arguments{
\item{pkg_data}{package metadata generated by pkgreview_getdata()}

\item{template}{character string, one of \code{review} or \code{editor}.}
}
\description{
Creates skeleton review files:
\itemize{
\item \code{index.Rmd}: \code{html_notebook} to perform and record review in
\item \code{README.md}: prepopulated README for review repo.
}
\code{index.Rmd} will be automatically added to \code{.Rbuildignore}. The resulting templates are populated with default
YAML frontmatter and R fenced code chunks (\code{Rmd}).
}
\examples{
\dontrun{
pkg_data <- pkgreview_getdata(pkg_dir)
pkgreview_index_rmd(pkg_data)
pkgreview_readme_md(pkg_data)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmd-utils.R
\name{pkgreview_print_source}
\alias{pkgreview_print_source}
\title{Print review package function source code}
\usage{
pkgreview_print_source(pkgname)
}
\arguments{
\item{pkgname}{character string. review package name}
}
\value{
prints out function source code for all exported functions.
}
\description{
Print review package function source code
}
\examples{
\dontrun{
library("pkgreviewr")
pkgreview_print_source("pkgreviewr")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgreview.R
\name{pkgreview_getdata}
\alias{pkgreview_getdata}
\title{pkgreview_getdata}
\usage{
pkgreview_getdata(
  pkg_dir = NULL,
  pkg_repo,
  template = c("review", "editor"),
  issue_no = NULL
)
}
\arguments{
\item{pkg_dir}{path to package source directory, cloned from github. Defaults
to the package source code directory in the review parent.}

\item{pkg_repo}{character string of the repo owner and name in the form of
\code{"owner/repo"}.}

\item{template}{character string, one of \code{review} or \code{editor}.}

\item{issue_no}{integer. Issue number of the pkg review in the rOpenSci \href{https://github.com/ropensci/software-review/issues}{\code{software-review} repository}.
If \code{NULL} (default), the issue number is extracted from the rOpenSci \strong{Under Review} badge on the pkg repository README.
Supplying an integer to \code{issue_no} overrides this behaviour and can be useful if a badge has not been added to the README yet.}
}
\value{
a list of package metadata
}
\description{
get package metadata from package source code.
}
\examples{
\dontrun{
# run from within a pkgreviewr project with the package source code in a
sibling directory
pkgreview_getdata("../rdflib")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgreview.R
\name{try_whoami}
\alias{try_whoami}
\title{Try whoami}
\usage{
try_whoami()
}
\value{
a list of whoami token metadata
}
\description{
Try to get whoami info from local  gh token.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgreview.R
\name{pkgreview_init}
\alias{pkgreview_init}
\title{Initialise pkgreview}
\usage{
pkgreview_init(
  pkg_repo,
  review_dir = ".",
  pkg_dir = NULL,
  template = c("review", "editor"),
  issue_no = NULL
)
}
\arguments{
\item{pkg_repo}{character string of the repo owner and name in the form of
\code{"owner/repo"}.}

\item{review_dir}{path to the review directory. Defaults to the working directory.}

\item{pkg_dir}{path to package source directory, cloned from github. Defaults
to the package source code directory in the review parent.}

\item{template}{character string, one of \code{review} or \code{editor}.}

\item{issue_no}{integer. Issue number of the pkg review in the rOpenSci \href{https://github.com/ropensci/software-review/issues}{\code{software-review} repository}.
If \code{NULL} (default), the issue number is extracted from the rOpenSci \strong{Under Review} badge on the pkg repository README.
Supplying an integer to \code{issue_no} overrides this behaviour and can be useful if a badge has not been added to the README yet.}
}
\value{
Initialisation creates pre-populated \code{index.Rmd}, \code{pkgreview.md} and \code{README.md} documents.
To initialise correctly, the function requires that the source code for the
package has been cloned. This might need to be done manually if it failed
during review creation. If setup is correct.
}
\description{
Initialise pkgreview
}
\examples{
\dontrun{
# run from within an uninitialised pkgreviewr project
pkgreview_init(pkg_repo = "ropensci/rdflib")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render-templates.R
\name{use_onboarding_tmpl}
\alias{use_onboarding_tmpl}
\title{Create software review/editor response template}
\usage{
use_onboarding_tmpl(template = c("review", "editor"))
}
\arguments{
\item{template}{character string, one of \code{review} or \code{editor}.}
}
\value{
writes a \verb{\{template\}.md} checklist template file in the project root.
}
\description{
Clone an up to date copy of the specified ropensci software review/editor response template.
}
\examples{
\dontrun{
use_onboarding_tmpl(template = "editor")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkgreview.R
\name{pkgreview_create}
\alias{pkgreview_create}
\title{Create a review project.}
\usage{
pkgreview_create(
  pkg_repo,
  review_parent = ".",
  template = c("review", "editor"),
  issue_no = NULL
)
}
\arguments{
\item{pkg_repo}{character string of the repo owner and name in the form of
\code{"owner/repo"}.}

\item{review_parent}{directory in which to setup review project and source package
source code.}

\item{template}{character string, one of \code{review} or \code{editor}.}

\item{issue_no}{integer. Issue number of the pkg review in the rOpenSci \href{https://github.com/ropensci/software-review/issues}{\code{software-review} repository}.
If \code{NULL} (default), the issue number is extracted from the rOpenSci \strong{Under Review} badge on the pkg repository README.
Supplying an integer to \code{issue_no} overrides this behaviour and can be useful if a badge has not been added to the README yet.}
}
\value{
setup review project with templates
}
\description{
Create and initialise an rOpenSci package review project
}
\examples{
\dontrun{
# for a review project
pkgreview_create(pkg_repo = "ropensci/rdflib", review_parent = "~/Documents/reviews/")
# for editors checks
pkgreview_create(pkg_repo = "ropensci/rdflib", review_parent = "~/Documents/editorials/",
template = "editor")
}
}
