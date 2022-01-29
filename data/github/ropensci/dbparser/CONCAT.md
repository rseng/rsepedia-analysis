
<!-- README.md is generated from README.Rmd. Please edit attributes file -->

# dbparser <img src="man/figures/logo.png" align="right" />

[![Build
Status](https://travis-ci.org/ropensci/dbparser.svg?branch=master)](https://travis-ci.org/ropensci/dbparser)
[![Build
status](https://ci.appveyor.com/api/projects/status/k18sqp55n39f3y5w?svg=true)](https://ci.appveyor.com/project/MohammedFCIS/dbparser)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/dbparser)](https://cran.r-project.org/package=dbparser)
[![codecov](https://codecov.io/gh/ropensci/dbparser/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/dbparser)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/grand-total/dbparser)](https://cran.r-project.org/package=dbparser)
[![Rdoc](https://www.rdocumentation.org/badges/version/dbparser)](https://www.rdocumentation.org/packages/dbparser)
[![CII Best
Practices](https://bestpractices.coreinfrastructure.org/projects/3311/badge)](https://bestpractices.coreinfrastructure.org/projects/3311)
[![](https://badges.ropensci.org/347_status.svg)](https://github.com/ropensci/software-review/issues/347)

## Introduction

[DrugBank](https://www.drugbank.ca/) database is a comprehensive, freely
accessible, online database containing information on drugs and drug
targets. As both a bioinformatics and a cheminformatics resource,
DrugBank combines detailed drug (i.e. chemical, pharmacological and
pharmaceutical) data with comprehensive drug target (i.e. sequence,
structure, and pathway) information. More information about DrugBank can
be found [here](https://www.drugbank.ca/about).

In its raw form, the DrugBank database is a single XML file. Users must
create an [account](https://www.drugbank.ca/public_users/sign_up) with
DrugBank and request permission to
[download](https://www.drugbank.ca/releases/latest) the database. Note
that this may take a couple of days.

The `dbparser` package parses the DrugBank XML database into `R` tibbles
that can be explored and analyzed by the user, check [this
tutorial](https://docs.ropensci.org/dbparser/articles/dbparser.html) for
more details.

Also, the package offers the option to save these tibbles in databases
including **SQL Server DB** and **Maria DB** just by enabling
`save_table` option, check [this
tutorial](https://docs.ropensci.org/dbparser/articles/Database_Saving.html)
for more details.

If you are waiting for access to the DrugBank database, or do not intend
to do a deep dive with the data, you may wish to use the `dbdataset`
[package](https://mohammedfcis.github.io/dbdataset/index.html), which
contains the DrugBank database already parsed into `R` tibbles. Note
that this is a large package that exceeds the limit set by CRAN. It is
only available on GitHub.

`dbparser` is tested against DrugBank versions *5.1.0* through *5.1.6*
successfully. If you find errors with these versions or any other
version please submit an issue
[here](https://github.com/ropensci/dbparser/issues).

## Installation

You can install the released version of dbparser from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dbparser")
```

or you can install the latest updates directly from the repo

``` r
library(devtools)
devtools::install_github("ropensci/dbparser")
```

## Code of Conduct

Please note that the ‘dbparser’ project is released with a [Contributor
Code of
Conduct](https://docs.ropensci.org/dbparser/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.

## Contributing Guide

👍🎉 First off, thanks for taking the time to contribute\! 🎉👍 Please
review our [Contributing
Guide](https://docs.ropensci.org/dbparser/CONTRIBUTING.html).

## Share the love ❤️

Think **dbparser** is useful? Let others discover it, by telling them in
person, via Twitter or a blog post.

Using **dbparser** for a paper you are writing? Consider citing it

``` r
citation("dbparser")
#> 
#> To cite dbparser in publications use:
#> 
#>   Mohammed Ali, Ali Ezzat (). dbparser: DrugBank Database XML Parser. R
#>   package version 1.2.0.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {DrugBank Database XML Parser},
#>     author = {Mohammed Ali and Ali Ezzat},
#>     organization = {Dainanahan},
#>     note = {R package version 1.2.0},
#>     url = {https://CRAN.R-project.org/package=dbparser},
#>   }
```
# dbparser 1.3.0

## Minor Fixes
* Removed RMariaDB dependency (#129)

# dbparser 1.2.0

## UI Changes
* Introduce progress bar in parser functions

## New Parsers
### Collective Parsers
* `drugs`, `cett` and `References` Parsers

### Elements Parsers
* `attachments` parsers for drugs and CETT
* `drug_pharmacology` parser
* Rename `drugs_books` parser to `drugs_textbooks`
* Rename `drug_all` parser to `run_all_parsers`
* renam `drug` parser to `drug_general_information`

## Documentation Update:
* Add returned parsed data structure 
* Explain the returned data functionality as a whole and for each elements
* Point out to related/similar parsers

## Package design
For those who thinking to contribute in `dbparser`, now parsers are implemented
as R6 classes.

## Minor Fixes
* Update database saver functions to accommodate new DrugBank data size.

# dbparser 1.1.2

### Major Changes
* Enhance many memory and performance issues for many parsers.
* Change the drug classification representations to extract more useful
information.
### Minor Changes
* Change some drug tibbles features names
### DEFUNCT
* Size columns in `drugs` main table is no longer exist, will do full 
statistical analysis later using dvminer package.

# dbparser 1.1.1

* Fix column size issue while importing into SQL Server (#91)
* Fix dbparser and upcoming CRAN release of dplyr issues (#92)
* Fix CRAN Notes (#93)
* Fix package documentation and site references


# dbparser 1.1.0
### Major Changes
* Functions have been splitted into 6 categories *DrugBank Database Loading,
Carriers, Targets, Transporters, Drug and common parsers*. All function names
are changed to reflect the function family. The related documentation is also
updated (#66, #75).
* `dbparser` now can cite the package by calling `citation("dbparser")` (#71).
* Adding more user friendly error messages (#76, #81).
* User can now pass `DBI` database connection to parser functions as an 
argument beside *SQLite* and *MariaDB* (#87).

### DEFUNCT
* `open_db`, `open_mdb` and `close_db` functions are no longer supported. 
Creating and maintaining database is completely user responsibility and the 
database connection can be passed to parser functions (#87).

### DOCUMENTATION FIXES
* New tutorials for how to use `dbparser` have been created (#78, #79).
* Contribution guide has been added.
* Code of conduct has been added (#70).
* Enhance function reference documentation to include section for each type (#68).

# dbparser 1.0.4
* Fix save drugs tibbles as csv several issues.
* Update sql database tibbles saver functions.
* Update sql database saver functions documentations.
* Support MariaDB and introduce related functionalities.

# dbparser 1.0.3
* Fix CRAN errors and notes

# dbparser 1.0.2
* Fix zip file location issue
* Replace Secondary and third keys columns from drug framework with *other_keys* column that contains any other keys that might exist in addition to the primary key
* Add **average-mass**, **monoisotopic-mass** and **calculated-properties** parsers.
* Support saving parsed drugs related parsed database as csv

# dbparser 1.0.1
* Fix CRAN Note
* Improve documentation
* Refactor unused functions
* Remove *Count* features from drug data set
* Fix several typos in documentation and code
* Fix consistency issue of CLASS of tibbles Returned by dbparser
* Check if drugbank database exist before parsing
* Add support for *international_brands* and *salts* elements
* Properly rename some features to have clear names
* Reduce datasets size by getting unique rows only
* Support reading zip file containing DrugBank xml database

# dbparser 1.0.0

* Initial release that contains core functionalities
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at mohammed.ali@edu.dsti.institute. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contributing to dbparser

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

First of all, thanks for considering contributing to dbparser! 👍 It's people like you that make it rewarding for us - the project maintainers - to work on dbparser. 😊

dbparser is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/ropensci/dbparser
[issues]: https://github.com/ropensci/dbparser/issues
[new_issue]: https://github.com/ropensci/dbparser/issues/new
[website]: https://docs.ropensci.org/dbparser/
[citation]: https://docs.ropensci.org/dbparser/authors.html
[email]: mohammed.ali@edu.dsti.institute

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think **dbparser** is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using **dbparser** for a paper you are writing? Consider [citing it][citation].

### Ask a question ⁉️

Using our_package and got stuck? Browse the [documentation][website] to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email][email].

### Propose an idea 💡

Have an idea for a new our_package feature? Take a look at the [documentation][website] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue]. While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using our_package and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub][new_issue] so we can fix it. A good bug report makes it easier for us to do so, so please include:

* Your operating system name and version (e.g. Mac OS 10.13.6).
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Improve the documentation 📖

Noticed a typo on the website? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

#### The website

[This website][website] is generated with [`pkgdown`](http://pkgdown.r-lib.org/). That means we don't have to write any html: content is pulled together from documentation in the code, vignettes, [Markdown](https://guides.github.com/features/mastering-markdown/) files, the package `DESCRIPTION` and `_pkgdown.yml` settings. If you know your way around `pkgdown`, you can [propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to improve documentation. If not, [report an issue][new_issue] and we can point you in the right direction.

#### Function documentation

Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you want to improve a function description:

1. Go to `R/` directory in the [code repository][repo].
2. Look for the file with the name of the function.
3. [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the roxygen comments (starting with `#'`).

### Contribute code 📝

Care to fix bugs or implement new functionality for our_package? Awesome! 👏 Have a look at the [issue list][issues] and leave a comment on the things you want to work on. See also the development guidelines below.

## Development guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
5. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors and warnings.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).
## Test environments
* local Windows 10(x86_64-w64-mingw32), R 4.0.0
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.0
* x86_64-w64-mingw32  (on appveyor), R 4.0.0
* Rhub
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
