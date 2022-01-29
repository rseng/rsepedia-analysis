
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
